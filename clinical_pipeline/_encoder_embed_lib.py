"""
_encoder_embed_lib.py
=====================
Shared helpers to export per-patient clinical-encoder embeddings + predicted
probabilities from a fine-tuned ModernBERT sequence classifier, with full
provenance (model / task / strategy / seed / split) for **leakage-safe**
multimodal fusion.

Used by:
  - 03_encoder_finetune.py            (INLINE: export from the in-memory best model)
  - 03e_build_embeddings_manifest.py  (manifest builder + optional --harvest)

Two embeddings per patient (both dim = hidden_size, float32):
  emb_pool : pooled last hidden state via model.config.classifier_pooling
             ("mean" masked-mean or "cls") — the sentence vector the head consumes.
  emb_head : the classifier's PRE-LOGIT input (post prediction-head), captured via
             a forward-pre-hook on the final classifier Linear, so it does not
             depend on exact submodule names across the 4 models.

Leakage contract
----------------
Every patient is tagged with (seed, split). For fusion at seed s: filter seed==s,
train the fusion head on split=="train", select on "val", report on "test". Those
val/test embeddings are out-of-sample for the seed-s encoder (val only influenced
encoder model-selection), and the same per-seed 80/10/10 split is shared with the
MRI/SNP modalities — so joining on (Patient_ID, seed) keeps every modality's test
patient held out everywhere. `in_task_cohort` lets fusion restrict to a task's real
cohort. `verify_split_integrity` re-checks exported splits against the source CSVs.

Mirrors the npz + provenance/QC-json + atomic-write conventions of
snp_pipeline/27_extract_fm_snp_embeddings.py and the masked mean-pool of
snp_pipeline/nt_v2/03_extract_patient_embeddings.py.
"""
from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

# torch / transformers are imported LAZILY inside load_encoder / pooled_predict only,
# so the consumer side (build_manifest, verify_split_integrity, reading parquet/npz)
# works in a torch-free environment.

ID_COL = "Patient_ID"
TEXT_COL = "Generated_Text"
SPLITS = ("train", "val", "test")


# ── Provenance ────────────────────────────────────────────────────────────────
def git_commit() -> str:
    """Best-effort short git SHA of this repo (empty string if unavailable)."""
    try:
        here = Path(__file__).resolve().parent
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=str(here),
                             capture_output=True, text=True, timeout=5)
        return out.stdout.strip() if out.returncode == 0 else ""
    except Exception:
        return ""


def make_provenance(model, model_slug, model_id, task, strategy, seed, max_length,
                    label_map, text_col, data_dir, source, checkpoint_path=None) -> dict:
    """Provenance dict derived from the loaded model + run identity."""
    try:
        import transformers
        tfx_version = transformers.__version__
    except Exception:
        tfx_version = "unknown"
    cfg = model.config
    try:
        dtype = str(next(model.parameters()).dtype)
    except StopIteration:
        dtype = "unknown"
    return {
        "model_slug": str(model_slug), "model_id": str(model_id),
        "task": str(task), "strategy": str(strategy), "seed": int(seed),
        "num_labels": int(cfg.num_labels), "hidden_size": int(cfg.hidden_size),
        "pooling": str(getattr(cfg, "classifier_pooling", "mean")),
        "max_length": int(max_length),
        "label_map": label_map, "text_col": str(text_col), "data_dir": str(data_dir),
        "source": str(source),
        "checkpoint_path": str(checkpoint_path) if checkpoint_path else None,
        "transformers_version": tfx_version,
        "model_param_dtype": dtype,
        "git_commit": git_commit(),
        "created_at": datetime.now(timezone.utc).isoformat(),
    }


# ── Patient universe for a seed (task-independent 80/10/10) ────────────────────
def build_seed_frame(data_dir, seed, task_cfg, text_col=TEXT_COL) -> pd.DataFrame:
    """Union of the seed's train/val/test split CSVs (the whole ~616-patient cohort,
    since the split is task-independent). Tags each patient with their seed fold
    (`split`), the task `label`, and `in_task_cohort` (had a valid task label AND
    passed the task's filter). Drops rows with empty text. Returns a DataFrame with
    columns: Patient_ID, text, split, label, in_task_cohort.
    `df.attrs['n_dropped_no_text']` records how many rows were dropped."""
    data_dir = Path(data_dir)
    frames = []
    for sp in SPLITS:
        df = pd.read_csv(data_dir / f"seed_{seed}" / f"{sp}.csv", low_memory=False)
        df["__split__"] = sp
        frames.append(df)
    full = pd.concat(frames, ignore_index=True)

    label_col = task_cfg["label_col"]
    if task_cfg.get("label_map") is not None:
        mapped = full[label_col].map(task_cfg["label_map"])
    else:
        mapped = pd.to_numeric(full[label_col], errors="coerce")

    in_cohort = mapped.notna()
    if task_cfg.get("filter_non_ad") and "Label_bl_multi" in full.columns:
        in_cohort = in_cohort & full["Label_bl_multi"].isin(["CN", "MCI"])

    out = pd.DataFrame({
        ID_COL: full[ID_COL].astype(str).values,
        "text": full[text_col].values,
        "split": full["__split__"].values,
        "label": pd.to_numeric(mapped, errors="coerce").values,  # NaN where not in cohort
        "in_task_cohort": in_cohort.values,
    })
    n_before = len(out)
    out = out[out["text"].notna() &
              (out["text"].astype(str).str.strip() != "")].reset_index(drop=True)
    out.attrs["n_dropped_no_text"] = int(n_before - len(out))
    return out


# ── Model loading + embedding extraction ──────────────────────────────────────
def load_encoder(ckpt_path, device=None):
    """Load (model, tokenizer) from a best_checkpoint/ dir (or HF id)."""
    import torch
    from transformers import AutoModelForSequenceClassification, AutoTokenizer
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    tok = AutoTokenizer.from_pretrained(str(ckpt_path))
    # Match the training-time model load: SDPA + reference_compile=False (ModernBERT's
    # torch.compile backward NaNs under torch 2.4.1; harmless to disable for inference).
    model = AutoModelForSequenceClassification.from_pretrained(
        str(ckpt_path), attn_implementation="sdpa", reference_compile=False)
    model.to(device).eval()
    return model, tok


def _find_classifier(model):
    """The final classification Linear (its INPUT is the pre-logit `emb_head`)."""
    import torch
    clf = getattr(model, "classifier", None)
    if isinstance(clf, torch.nn.Linear):
        return clf
    n = getattr(model.config, "num_labels", None)
    last = None
    for m in model.modules():
        if isinstance(m, torch.nn.Linear) and (n is None or m.out_features == n):
            last = m
    return last


def _pool_last_hidden(last_hidden, attention_mask, pooling):
    if pooling == "cls":
        return last_hidden[:, 0]
    m = attention_mask.unsqueeze(-1).to(last_hidden.dtype)          # masked mean
    return (last_hidden * m).sum(1) / m.sum(1).clamp_min(1)


def pooled_predict(model, tokenizer, texts, max_length=1024, batch_size=32, device=None):
    """Run frozen forward over `texts`; return (emb_pool[N,D], emb_head[N,D],
    probs[N,C], preds[N]).

    emb_pool: pooled last hidden state via config.classifier_pooling, captured as
              the INPUT to model.head when available (else recomputed from hidden
              states). emb_head: INPUT to the final classifier Linear (post head),
              captured via forward-pre-hook. Both cast to float32."""
    import torch
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device).eval()
    pooling = str(getattr(model.config, "classifier_pooling", "mean"))

    head_mod = getattr(model, "head", None)        # input == pooled vector
    clf_mod = _find_classifier(model)              # input == pre-logit vector
    captured = {}

    def _pre_hook(key):
        def hook(_module, args):
            if args:
                captured[key] = args[0].detach()
        return hook

    hooks = []
    if isinstance(head_mod, torch.nn.Module):
        hooks.append(head_mod.register_forward_pre_hook(_pre_hook("pool")))
    if clf_mod is not None:
        hooks.append(clf_mod.register_forward_pre_hook(_pre_hook("head")))

    pool_c, head_c, logit_c = [], [], []
    texts = [("" if t is None else str(t)) for t in texts]
    try:
        with torch.no_grad():
            for i in range(0, len(texts), batch_size):
                enc = tokenizer(texts[i:i + batch_size], truncation=True,
                                max_length=max_length, padding=True,
                                return_tensors="pt").to(device)
                captured.clear()
                out = model(**enc)
                logit_c.append(out.logits.float().cpu())

                if "pool" in captured:
                    pool_c.append(captured["pool"].float().cpu())
                else:                                   # fallback: recompute pooled
                    out2 = model(**enc, output_hidden_states=True)
                    pooled = _pool_last_hidden(out2.hidden_states[-1],
                                               enc["attention_mask"], pooling)
                    pool_c.append(pooled.float().cpu())

                head_c.append(captured["head"].float().cpu() if "head" in captured
                              else pool_c[-1].clone())
    finally:
        for h in hooks:
            h.remove()

    emb_pool = torch.cat(pool_c).numpy().astype(np.float32)
    emb_head = torch.cat(head_c).numpy().astype(np.float32)
    logits = torch.cat(logit_c).numpy().astype(np.float32)
    probs = torch.softmax(torch.from_numpy(logits), dim=-1).numpy().astype(np.float32)
    preds = probs.argmax(1).astype(np.int64)
    return emb_pool, emb_head, probs, preds


# ── Atomic writers ────────────────────────────────────────────────────────────
def _atomic_json(obj, path):
    path = Path(path); tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=2, default=str)
    tmp.replace(path)


def _atomic_parquet(df, path):
    path = Path(path); tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_parquet(tmp, index=False)
    tmp.replace(path)


def _atomic_npz(path, **arrays):
    path = Path(path); tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "wb") as f:                       # file handle => no auto ".npz"
        np.savez_compressed(f, **arrays)
    tmp.replace(path)


def _norm_stats(emb):
    nrm = np.linalg.norm(emb, axis=1)
    return {"mean": float(nrm.mean()), "min": float(nrm.min()), "max": float(nrm.max())}


# ── Per-run writer ────────────────────────────────────────────────────────────
def write_run_embeddings(run_embed_dir, frame, emb_pool, emb_head, probs, preds,
                         provenance) -> Path:
    """Write embeddings.parquet + embeddings.npz + embeddings_meta.json (atomic)."""
    run_embed_dir = Path(run_embed_dir)
    run_embed_dir.mkdir(parents=True, exist_ok=True)

    N, C = probs.shape
    Dp, Dh = emb_pool.shape[1], emb_head.shape[1]
    pid = frame[ID_COL].astype(str).to_numpy()
    split = frame["split"].astype(str).to_numpy()
    in_cohort = frame["in_task_cohort"].astype(bool).to_numpy()
    label = pd.to_numeric(frame["label"], errors="coerce").astype("float64").to_numpy()

    # npz FIRST (numpy-only, always works) — compact arrays + scalar metadata.
    scalar_keys = ("model_slug", "model_id", "task", "strategy", "seed", "num_labels",
                   "hidden_size", "pooling", "max_length", "transformers_version",
                   "model_param_dtype", "source", "git_commit", "created_at")
    _atomic_npz(
        run_embed_dir / "embeddings.npz",
        patient_id=pid.astype("U"), split=split.astype("U"),
        in_task_cohort=in_cohort, label=label, pred=preds,
        probs=probs, emb_pool=emb_pool, emb_head=emb_head,
        **{f"meta_{k}": np.array(str(provenance.get(k))) for k in scalar_keys},
    )

    # meta.json (full provenance + QC)
    meta = dict(provenance)
    meta.update({
        "n_total": int(N),
        "n_in_task": int(in_cohort.sum()),
        "n_dropped_no_text": int(frame.attrs.get("n_dropped_no_text", 0)),
        "counts_by_split": {sp: int((split == sp).sum()) for sp in SPLITS},
        "emb_pool_dim": int(Dp), "emb_head_dim": int(Dh), "num_classes": int(C),
        "qc": {
            "emb_pool_all_finite": bool(np.isfinite(emb_pool).all()),
            "emb_head_all_finite": bool(np.isfinite(emb_head).all()),
            "probs_all_finite": bool(np.isfinite(probs).all()),
            "emb_pool_norm": _norm_stats(emb_pool),
            "emb_head_norm": _norm_stats(emb_head),
            "prob_min": float(probs.min()), "prob_max": float(probs.max()),
        },
    })
    _atomic_json(meta, run_embed_dir / "embeddings_meta.json")

    # parquet LAST (needs a parquet engine, e.g. pyarrow). npz + meta already landed,
    # so a missing engine degrades gracefully instead of losing the run.
    data = {
        ID_COL: pid,
        "model_slug": provenance["model_slug"], "model_id": provenance["model_id"],
        "task": provenance["task"], "strategy": provenance["strategy"],
        "seed": provenance["seed"], "split": split,
        "in_task_cohort": in_cohort, "label": label, "pred": preds,
    }
    for c in range(C):
        data[f"prob_{c}"] = probs[:, c]
    for j in range(Dp):
        data[f"emb_pool_{j:04d}"] = emb_pool[:, j]
    for j in range(Dh):
        data[f"emb_head_{j:04d}"] = emb_head[:, j]
    try:
        _atomic_parquet(pd.DataFrame(data), run_embed_dir / "embeddings.parquet")
    except Exception as e:
        print(f"  [WARN] parquet not written ({e}); npz + meta.json are intact. "
              f"Install pyarrow to enable parquet.")
    return run_embed_dir


# ── Manifest + leakage verification ───────────────────────────────────────────
def verify_split_integrity(run_embed_dir, data_dir, seed):
    """Check exported `split` labels match the source seed CSV partition (no leakage
    of fold assignment). Reads the always-present npz (falls back to parquet).
    Returns (ok: bool, message: str)."""
    run_embed_dir = Path(run_embed_dir)
    npz = run_embed_dir / "embeddings.npz"
    if npz.exists():
        z = np.load(npz, allow_pickle=False)
        pids = z["patient_id"].astype(str)
        splits = z["split"].astype(str)
    else:
        df = pd.read_parquet(run_embed_dir / "embeddings.parquet", columns=[ID_COL, "split"])
        pids = df[ID_COL].astype(str).to_numpy()
        splits = df["split"].astype(str).to_numpy()
    truth = {}
    for sp in SPLITS:
        ids = pd.read_csv(Path(data_dir) / f"seed_{seed}" / f"{sp}.csv",
                          usecols=[ID_COL])[ID_COL].astype(str)
        for pid in ids:
            truth[pid] = sp
    bad = [(p, s) for p, s in zip(pids, splits) if truth.get(p) != s]
    dup = int(len(pids) - len(set(pids)))
    if dup:
        return False, f"{dup} duplicate patient rows"
    if bad:
        return False, f"{len(bad)} split mismatches (e.g. {bad[:3]})"
    return True, "ok"


def build_manifest(out_root) -> tuple[Path, pd.DataFrame]:
    """Scan every embeddings/embeddings_meta.json under out_root and write
    embeddings_manifest.csv — the catalogue the fusion step reads."""
    out_root = Path(out_root)
    rows = []
    for meta_path in sorted(out_root.rglob("embeddings/embeddings_meta.json")):
        try:
            meta = json.loads(meta_path.read_text())
        except Exception:
            continue
        run_dir = meta_path.parent
        cby = meta.get("counts_by_split", {})
        rows.append({
            "model_slug": meta.get("model_slug"), "task": meta.get("task"),
            "strategy": meta.get("strategy"), "seed": meta.get("seed"),
            "n_total": meta.get("n_total"),
            "n_train": cby.get("train"), "n_val": cby.get("val"), "n_test": cby.get("test"),
            "n_in_task": meta.get("n_in_task"),
            "emb_pool_dim": meta.get("emb_pool_dim"), "emb_head_dim": meta.get("emb_head_dim"),
            "num_labels": meta.get("num_labels"), "pooling": meta.get("pooling"),
            "transformers_version": meta.get("transformers_version"),
            "model_param_dtype": meta.get("model_param_dtype"),
            "git_commit": meta.get("git_commit"), "source": meta.get("source"),
            "parquet_path": str(run_dir / "embeddings.parquet"),
            "npz_path": str(run_dir / "embeddings.npz"),
            "created_at": meta.get("created_at"),
        })
    man = pd.DataFrame(rows)
    if not man.empty:
        man = man.sort_values(["model_slug", "task", "strategy", "seed"]).reset_index(drop=True)
    out_csv = out_root / "embeddings_manifest.csv"
    man.to_csv(out_csv, index=False)
    return out_csv, man
