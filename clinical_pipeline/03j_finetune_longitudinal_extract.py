r"""
03j_finetune_longitudinal_extract.py
====================================
TASK 1 of the T4 longitudinal-MAMBA arm: fine-tune ONE BioClinical-ModernBERT-large model
on ALL clinical visits (contemporaneous CN/MCI/AD label) and extract per-visit embeddings in
ONE consistent space for the visits baseline -> 1 year. No checkpoint is persisted (HPC space
is tight) — embeddings + metrics are exported directly, exactly like 03_encoder_finetune.py.

Why a dedicated script (not a new task in 03_encoder_finetune.py): the training corpus is
multi-row-per-patient (all visits) and the EXTRACTION is per-visit over a month window with
VISCODE tags — neither fits the one-row-per-patient `build_seed_frame`/`write_run_embeddings`
path. We reuse 03_encoder_finetune.py's classes (`ClinicalDataset`, `WeightedTrainer`,
`make_compute_metrics`, `full_metrics`, `load_split`) and _encoder_embed_lib.py's
(`pick_attn_impl`, `pooled_predict`, `make_provenance`) by import — NOT editing them.

Task (fixed)
------------
  label = Label_visit_diag (concurrent dx), 3-class {CN:0, MCI:1, AD:2}, multiclass, full_ft.

Inputs
------
  Per-seed ALL-VISITS splits from 01i_build_longitudinal_visit_splits.py:
    <data_dir>/seed_{N}/{train,val,test}.csv   (default data_dir = verbose/longitudinal_allvisits)
  Each row = one patient-visit, carries Patient_ID, VISCODE_2, visit_months, Generated_Text,
  Label_visit_diag.

Outputs  (out_dir/<model>/T2_long_multiclass/seed_{N}/full_ft/)
-------
  metrics.json                 — config + standard val/test metrics (all visits) + per-VISCODE
                                 (bl/m06/m12) bACC/macro-F1 breakdown on val & test.
  embeddings/embeddings.parquet — ONE ROW PER IN-WINDOW VISIT (VISCODE in {bl,m06,m12}; excludes
                                 screening 'sc' and >m12). Columns: Patient_ID, VISCODE_2,
                                 visit_months, split, seed, label, pred, prob_0/1/2,
                                 emb_pool_0000.. (1024), emb_head_0000.. (1024).
  embeddings/embeddings.npz, embeddings_meta.json — numpy mirror + provenance/QC.
  NO model weights. trainer_ckpts/ is deleted at the end.

  NB: extraction keeps EVERY in-window visit of EVERY patient — no "must have bl+m06+m12"
  filter. Completeness filtering is deferred to the downstream MAMBA step.

Env: conda env `clinical` (torch 2.3.1, transformers >=4.47). CSD3 ampere / dept L4.
Usage
-----
  python clinical_pipeline/03j_finetune_longitudinal_extract.py --seed 0 \
      --data_dir D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified_post_exclusion/verbose/longitudinal_allvisits \
      --out_dir  .../encoder_outputs_no_cdr_post_exclusion_longitudinal
  # Smoke: add --epochs 1 --max_length 256
"""
import argparse
import importlib.util
import json
import os
import shutil
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import balanced_accuracy_score, f1_score
from sklearn.utils.class_weight import compute_class_weight
from transformers import (AutoModelForSequenceClassification, AutoTokenizer,
                          EarlyStoppingCallback, TrainingArguments)

warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def _load_module(name, filename):
    """Import a sibling module whose filename may start with a digit."""
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_ft = _load_module("ft03", "03_encoder_finetune.py")
import _encoder_embed_lib as emb  # no digit prefix -> direct import

ClinicalDataset = _ft.ClinicalDataset
WeightedTrainer = _ft.WeightedTrainer
make_compute_metrics = _ft.make_compute_metrics
full_metrics = _ft.full_metrics
load_split = _ft.load_split

TASK_NAME = "T2_long_multiclass"
TASK_CFG = {
    "label_col":     "Label_visit_diag",
    "num_labels":    3,
    "task_type":     "multiclass",
    "label_map":     {"CN": 0, "MCI": 1, "AD": 2},
    "filter_non_ad": False,
    "description":   "Longitudinal CN/MCI/AD across ALL visits (concurrent dx); single embedding space",
}


# ── In-window per-visit extraction frame ───────────────────────────────────────
def build_visit_frame(data_dir, seed, window_months, text_col="Generated_Text"):
    """Concatenate the seed's train/val/test split rows, tag `split`, map the label,
    and keep ONLY in-window visits (0 <= visit_months <= window). One row per visit.
    Returns a DataFrame with Patient_ID, VISCODE_2, visit_months, split, text, label,
    in_task_cohort (label present)."""
    label_map = TASK_CFG["label_map"]
    frames = []
    for sp in ("train", "val", "test"):
        df = pd.read_csv(Path(data_dir) / f"seed_{seed}" / f"{sp}.csv", low_memory=False)
        df["__split__"] = sp
        frames.append(df)
    full = pd.concat(frames, ignore_index=True)
    full = full[full["visit_months"].between(0, window_months)].copy()

    mapped = full[TASK_CFG["label_col"]].map(label_map)
    out = pd.DataFrame({
        "Patient_ID":   full["Patient_ID"].astype(str).values,
        "VISCODE_2":    full["VISCODE_2"].astype(str).values,
        "visit_months": pd.to_numeric(full["visit_months"], errors="coerce").values,
        "split":        full["__split__"].values,
        "text":         full[text_col].values,
        "label":        pd.to_numeric(mapped, errors="coerce").values,
    })
    out["in_task_cohort"] = out["label"].notna().values
    n_before = len(out)
    out = out[out["text"].notna() & (out["text"].astype(str).str.strip() != "")]
    out = out.reset_index(drop=True)
    out.attrs["n_dropped_no_text"] = int(n_before - len(out))
    return out


def per_viscode_metrics(frame, preds):
    """bACC / macro-F1 / n per (split, VISCODE) over rows with a valid label."""
    df = frame.copy()
    df["pred"] = preds
    df = df[df["label"].notna()]
    res = {}
    for sp in ("val", "test"):
        res[sp] = {}
        sub = df[df["split"] == sp]
        for vc in sorted(sub["VISCODE_2"].unique()):
            s = sub[sub["VISCODE_2"] == vc]
            y, p = s["label"].astype(int).to_numpy(), s["pred"].astype(int).to_numpy()
            res[sp][vc] = {
                "n": int(len(s)),
                "balanced_acc": round(float(balanced_accuracy_score(y, p)), 4) if len(s) else None,
                "macro_f1": round(float(f1_score(y, p, average="macro", zero_division=0)), 4) if len(s) else None,
            }
    return res


# ── Custom per-visit writer (adds VISCODE_2 + visit_months to the schema) ───────
def write_visit_embeddings(run_dir, frame, emb_pool, emb_head, probs, preds, provenance):
    run_dir = Path(run_dir); run_dir.mkdir(parents=True, exist_ok=True)
    N, C = probs.shape
    Dp, Dh = emb_pool.shape[1], emb_head.shape[1]
    pid = frame["Patient_ID"].astype(str).to_numpy()
    vc = frame["VISCODE_2"].astype(str).to_numpy()
    mo = pd.to_numeric(frame["visit_months"], errors="coerce").to_numpy()
    split = frame["split"].astype(str).to_numpy()
    in_cohort = frame["in_task_cohort"].astype(bool).to_numpy()
    label = pd.to_numeric(frame["label"], errors="coerce").astype("float64").to_numpy()

    emb._atomic_npz(
        run_dir / "embeddings.npz",
        patient_id=pid.astype("U"), viscode=vc.astype("U"), visit_months=mo,
        split=split.astype("U"), in_task_cohort=in_cohort, label=label, pred=preds,
        probs=probs, emb_pool=emb_pool, emb_head=emb_head,
        **{f"meta_{k}": np.array(str(provenance.get(k)))
           for k in ("model_slug", "model_id", "task", "strategy", "seed", "created_at")},
    )

    meta = dict(provenance)
    meta.update({
        "n_total": int(N), "n_in_task": int(in_cohort.sum()),
        "n_dropped_no_text": int(frame.attrs.get("n_dropped_no_text", 0)),
        "counts_by_split": {sp: int((split == sp).sum()) for sp in ("train", "val", "test")},
        "counts_by_viscode": {v: int((vc == v).sum()) for v in sorted(set(vc.tolist()))},
        "emb_pool_dim": int(Dp), "emb_head_dim": int(Dh), "num_classes": int(C),
        "qc": {
            "emb_pool_all_finite": bool(np.isfinite(emb_pool).all()),
            "emb_head_all_finite": bool(np.isfinite(emb_head).all()),
            "probs_all_finite": bool(np.isfinite(probs).all()),
            "emb_pool_norm": emb._norm_stats(emb_pool),
            "prob_min": float(probs.min()), "prob_max": float(probs.max()),
        },
    })
    emb._atomic_json(meta, run_dir / "embeddings_meta.json")

    data = {
        "Patient_ID": pid, "VISCODE_2": vc, "visit_months": mo,
        "model_slug": provenance["model_slug"], "task": provenance["task"],
        "strategy": provenance["strategy"], "seed": provenance["seed"],
        "split": split, "in_task_cohort": in_cohort, "label": label, "pred": preds,
    }
    for c in range(C):
        data[f"prob_{c}"] = probs[:, c]
    for j in range(Dp):
        data[f"emb_pool_{j:04d}"] = emb_pool[:, j]
    for j in range(Dh):
        data[f"emb_head_{j:04d}"] = emb_head[:, j]
    try:
        emb._atomic_parquet(pd.DataFrame(data), run_dir / "embeddings.parquet")
    except Exception as e:
        print(f"  [WARN] parquet not written ({e}); npz + meta.json are intact.")
    return run_dir


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--model_id", type=str, default="thomas-sounack/BioClinical-ModernBERT-large")
    p.add_argument("--seed", type=int, required=True, choices=[0, 1, 2])
    p.add_argument("--data_dir", type=str,
                   default=r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion\verbose\longitudinal_allvisits")
    p.add_argument("--out_dir", type=str,
                   default=r"D:\ADNI_BIDS_project\derivatives\encoder_outputs_no_cdr_post_exclusion_longitudinal")
    p.add_argument("--text_col", type=str, default="Generated_Text")
    p.add_argument("--max_length", type=int, default=1024)
    p.add_argument("--extract_window_months", type=int, default=12)
    p.add_argument("--epochs", type=int, default=10)
    p.add_argument("--patience", type=int, default=5)
    p.add_argument("--batch_size", type=int, default=16)
    p.add_argument("--lr", type=float, default=2e-5)
    p.add_argument("--warmup_ratio", type=float, default=0.1)
    p.add_argument("--weight_decay", type=float, default=1e-5)
    p.add_argument("--hf_cache", type=str, default=None)
    return p.parse_args()


def main():
    args = parse_args()
    strategy = "full_ft"
    model_slug = args.model_id.split("/")[-1]
    hf_cache = args.hf_cache or os.environ.get("HF_HOME", None)
    out_dir = Path(args.out_dir) / model_slug / TASK_NAME / f"seed_{args.seed}" / strategy
    out_dir.mkdir(parents=True, exist_ok=True)
    trainer_dir = out_dir / "trainer_ckpts"
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    print(f"\n{'='*65}\n  Model : {args.model_id}\n  Task  : {TASK_CFG['description']}\n"
          f"  Seed  : {args.seed} | full_ft | epochs {args.epochs} | device {device}\n"
          f"  Out   : {out_dir}\n{'='*65}\n")
    torch.manual_seed(args.seed); np.random.seed(args.seed)

    # ── Data (ALL visits; load_split drops NaN-label rows) ─────────────────────
    data_dir = Path(args.data_dir)
    tr_texts, tr_labels = load_split(data_dir / f"seed_{args.seed}" / "train.csv", TASK_CFG, args.text_col)
    va_texts, va_labels = load_split(data_dir / f"seed_{args.seed}" / "val.csv",   TASK_CFG, args.text_col)
    te_texts, te_labels = load_split(data_dir / f"seed_{args.seed}" / "test.csv",  TASK_CFG, args.text_col)
    print(f"  Train {len(tr_labels)}  Val {len(va_labels)}  Test {len(te_labels)} (all-visit rows)")
    uq, ct = np.unique(tr_labels, return_counts=True)
    print(f"  Train class dist: {dict(zip(uq.tolist(), ct.tolist()))}")

    tok = AutoTokenizer.from_pretrained(args.model_id, cache_dir=hf_cache)
    train_ds = ClinicalDataset(tr_texts, tr_labels, tok, args.max_length)
    val_ds   = ClinicalDataset(va_texts, va_labels, tok, args.max_length)
    test_ds  = ClinicalDataset(te_texts, te_labels, tok, args.max_length)

    attn_impl = emb.pick_attn_impl()
    print(f"  Attention impl: {attn_impl}")
    model = AutoModelForSequenceClassification.from_pretrained(
        args.model_id, num_labels=TASK_CFG["num_labels"], cache_dir=hf_cache,
        ignore_mismatched_sizes=True, attn_implementation=attn_impl, reference_compile=False)
    print(f"  Full fine-tune. Trainable params: {sum(p.numel() for p in model.parameters()):,}")

    classes = np.unique(tr_labels)
    cw = compute_class_weight("balanced", classes=classes, y=np.array(tr_labels))
    class_weights = torch.tensor(cw, dtype=torch.float)
    print(f"  Class weights: {dict(zip(classes.tolist(), cw.round(3).tolist()))}")

    use_bf16 = torch.cuda.is_available() and torch.cuda.is_bf16_supported()
    use_fp16 = torch.cuda.is_available() and not use_bf16
    print(f"  Precision: {'bf16' if use_bf16 else 'fp16' if use_fp16 else 'fp32'}")

    training_args = TrainingArguments(
        output_dir=str(trainer_dir), num_train_epochs=args.epochs,
        per_device_train_batch_size=args.batch_size, per_device_eval_batch_size=args.batch_size,
        learning_rate=args.lr, weight_decay=args.weight_decay, warmup_ratio=args.warmup_ratio,
        eval_strategy="epoch", save_strategy="epoch", load_best_model_at_end=True,
        metric_for_best_model="eval_val_loss", greater_is_better=False,
        save_total_limit=1, save_only_model=True, bf16=use_bf16, fp16=use_fp16,
        logging_steps=20, logging_dir=str(out_dir / "tb_logs"), report_to="none",
        seed=args.seed, data_seed=args.seed, remove_unused_columns=False,
        dataloader_num_workers=2, dataloader_pin_memory=True)

    trainer = WeightedTrainer(
        class_weights=class_weights, model=model, args=training_args,
        train_dataset=train_ds, eval_dataset={"val": val_ds, "train": train_ds},
        compute_metrics=make_compute_metrics(TASK_CFG["task_type"]),
        callbacks=[EarlyStoppingCallback(early_stopping_patience=args.patience)])

    print("\n  Training ...")
    trainer.train()

    # ── Standard val/test metrics over ALL visits ──────────────────────────────
    print("\n  Final evaluation (val + test, all visits) ...")
    val_out = trainer.predict(val_ds, metric_key_prefix="val")
    test_out = trainer.predict(test_ds, metric_key_prefix="test")
    val_metrics, _, _ = full_metrics(val_out.label_ids, val_out.predictions, TASK_CFG["task_type"])
    test_metrics, _, _ = full_metrics(test_out.label_ids, test_out.predictions, TASK_CFG["task_type"])
    val_metrics["loss"] = round(float(val_out.metrics.get("val_loss", float("nan"))), 4)
    test_metrics["loss"] = round(float(test_out.metrics.get("test_loss", float("nan"))), 4)
    print(f"  Test (all visits): bACC={test_metrics['balanced_acc']} macroF1={test_metrics['macro_f1']}")

    # ── Per-visit extraction (in-window only) + per-VISCODE metrics ────────────
    print("\n  Extracting per-visit embeddings (in-window) ...")
    frame = build_visit_frame(args.data_dir, args.seed, args.extract_window_months, args.text_col)
    emb_pool, emb_head, probs, preds = emb.pooled_predict(
        model, tok, frame["text"].tolist(), max_length=args.max_length, batch_size=32, device=device)
    pv = per_viscode_metrics(frame, preds)
    print(f"  Extracted {len(frame)} in-window visit rows "
          f"({frame['VISCODE_2'].value_counts().to_dict()})")
    for sp in ("val", "test"):
        print(f"    per-VISCODE {sp} bACC: " +
              ", ".join(f"{vc}={pv[sp][vc]['balanced_acc']}" for vc in pv[sp]))

    prov = emb.make_provenance(model, model_slug, args.model_id, TASK_NAME, strategy,
                               args.seed, args.max_length, TASK_CFG["label_map"],
                               args.text_col, args.data_dir, source="inline")
    write_visit_embeddings(out_dir / "embeddings", frame, emb_pool, emb_head, probs, preds, prov)
    print(f"  Embeddings -> {out_dir / 'embeddings'} (emb_pool={emb_pool.shape[1]})")

    # ── metrics.json ───────────────────────────────────────────────────────────
    config = {
        "model_id": args.model_id, "task": TASK_NAME, "task_description": TASK_CFG["description"],
        "seed": args.seed, "strategy": strategy, "epochs": args.epochs, "lr": args.lr,
        "max_length": args.max_length, "extract_window_months": args.extract_window_months,
        "n_train": len(tr_labels), "n_val": len(va_labels), "n_test": len(te_labels),
        "n_extracted_in_window": int(len(frame)),
        "class_weights": dict(zip(classes.tolist(), cw.round(4).tolist())),
        "best_epoch": int(trainer.state.best_model_checkpoint.split("-")[-1])
                      if trainer.state.best_model_checkpoint else args.epochs,
        "timestamp": datetime.now().isoformat(),
    }
    results = {"config": config, "val_metrics": val_metrics, "test_metrics": test_metrics,
               "per_viscode_metrics": pv}
    with open(out_dir / "metrics.json", "w") as f:
        json.dump(results, f, indent=2)
    pd.DataFrame(trainer.state.log_history).to_csv(out_dir / "training_log.csv", index=False)
    print(f"  Saved -> {out_dir / 'metrics.json'}")

    # ── No weights kept; clean trainer checkpoints ─────────────────────────────
    if trainer_dir.exists():
        shutil.rmtree(trainer_dir, ignore_errors=True)
        print("  Cleaned trainer_ckpts (no weights persisted).")
    print("\n  Done.")


if __name__ == "__main__":
    main()
