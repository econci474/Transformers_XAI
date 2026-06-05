r"""
04_train_classifier_mamba.py  (env: clinical / torch) — runs on L4 or CSD3
==========================================================================
Train the longitudinal-MAMBA 3-class T4-HORIZON CLASSIFIER head (<3y / 3-7y / >=7y) over the backbone/
training/time/align ablation, and export per config x seed the per-patient predicted class probabilities.
This is the OTHER half of the two-head arm (the survival head is 01/02). It is a SEPARATE evaluation on a
DIFFERENT cohort from the survival-derived horizon: the T4 converters on the **T4-aware** leakage-safe
folds (Task 1b: encoder trained out-of-fold w.r.t. the held-out T4 converters; 116/15/15 per seed).

Sweep A only (one head — the classifier; weighted-CE). mamba2 is dropped (its transformers slow-path OOMs
on a single A100, see survival run). Heads/parametrisations are a survival-only concept, so there is no
Sweep B here.

Inputs : encoder_outputs_..._longitudinal_T4split/.../seed_{s}/full_ft/embeddings/embeddings.parquet
         + the longitudinal master (Label_T4 via converter group + years_to_AD bucketing).
Outputs: <out_dir>/<config>/seed_{s}/predictions.parquet (Patient_ID,split,y_true,pred,prob_0..2)
         + meta.json (config, val_ce, n).  NO model weights kept.

Usage:
  conda run -n clinical python integration/T4/mamba_long_clinical/04_train_classifier_mamba.py \
      --seeds 0 1 2 [--only A_default_mamba1_frozen] [--device cuda]
"""
import argparse
import json
import sys
import traceback
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _mamba_seq_lib as lib

EMB_ROOT = Path(r"D:\ADNI_BIDS_project\derivatives\encoder_outputs_no_cdr_post_exclusion_longitudinal_T4split"
                r"\BioClinical-ModernBERT-large\T2_long_multiclass")
MASTER = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion"
              r"\verbose\longitudinal\master_clinical_verbose.csv")
OUT = HERE / "outputs" / "classifier"


def CONFIGS():
    enc = lambda **kw: {**dict(backbone="mamba1", frozen=True, time_concat="append"), **kw}
    return {
        # ── Sweep A (the classifier head; weighted-CE) ──
        "A_default_mamba1_frozen": (enc(), 0.0, 1e-3),
        "A_q2_finetune":           (enc(frozen=False), 0.0, 1e-4),
        "A_q3_time_prepend":       (enc(time_concat="prepend"), 0.0, 1e-3),
        "A_q4_align":              (enc(), 0.5, 1e-3),
        "A_ctrl_meanpool":         (enc(backbone="meanpool"), 0.0, 1e-3),
        "A_ctrl_gru":              (enc(backbone="gru"), 0.0, 1e-3),
    }


def run_config(name, cfg, seed, master_labels, device):
    enc_kw, align_lambda, lr = cfg
    lib.seed_everything(seed)
    emb_parquet = EMB_ROOT / f"seed_{seed}" / "full_ft" / "embeddings" / "embeddings.parquet"
    seq = lib.assemble_sequences(emb_parquet, master_labels, "classifier")

    encoder = lib.make_encoder(emb_dim=seq["emb_dim"], **enc_kw)
    head = lib.make_head("classifier", encoder.d_out, n_classes=3)
    encoder, head, val_ce = lib.train_classifier(
        encoder, head, seq, device=device, lr=lr, align_lambda=align_lambda, n_classes=3)

    import pandas as pd
    idx = np.arange(len(seq["pid"]))
    P = lib.predict_proba(encoder, head, seq, idx, device=device)
    df = pd.DataFrame({"Patient_ID": seq["pid"], "split": seq["split"], "y_true": seq["y"],
                       "pred": P.argmax(1)})
    for k in range(P.shape[1]):
        df[f"prob_{k}"] = P[:, k]
    out_dir = OUT / name / f"seed_{seed}"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_dir / "predictions.parquet", index=False)
    json.dump({"config": name, "seed": seed, "head": "classifier", "encoder": enc_kw,
               "align_lambda": align_lambda, "lr": lr, "val_ce": float(val_ce), "n": int(len(idx))},
              open(out_dir / "meta.json", "w"), indent=2)
    # quick train-side balanced accuracy on test for the log
    from sklearn.metrics import balanced_accuracy_score
    te = seq["split"] == "test"
    bacc = balanced_accuracy_score(seq["y"][te], P[te].argmax(1)) if te.sum() else float("nan")
    print(f"  [{name} seed{seed}] val_ce={val_ce:.4f}  test_bACC={bacc:.4f}  N={len(idx)} -> {out_dir}")


def main():
    global EMB_ROOT, MASTER, OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--only", type=str, nargs="*", default=None, help="subset of config names")
    ap.add_argument("--device", type=str, default=None)
    ap.add_argument("--emb_root", type=str, default=str(EMB_ROOT),
                    help="dir containing seed_*/full_ft/embeddings/embeddings.parquet (T4split)")
    ap.add_argument("--master", type=str, default=str(MASTER))
    ap.add_argument("--out_dir", type=str, default=str(OUT),
                    help="where to write predictions (pass a hpc-work path on CSD3 so nothing lands in repo)")
    args = ap.parse_args()
    EMB_ROOT = Path(args.emb_root); MASTER = Path(args.master); OUT = Path(args.out_dir)
    import torch
    device = args.device or ("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device={device}  emb_root={EMB_ROOT}")

    master_labels = lib.load_master_labels(MASTER)
    cfgs = CONFIGS()
    names = args.only or list(cfgs)
    for name in names:
        for seed in args.seeds:
            try:
                run_config(name, cfgs[name], seed, master_labels, device)
            except Exception:
                print(f"  [FAIL {name} seed{seed}]\n{traceback.format_exc()}")
    print("Done.")


if __name__ == "__main__":
    main()
