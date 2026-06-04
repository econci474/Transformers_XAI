r"""
01_train_survival_mamba.py  (env: clinical / torch) — runs on L4 or CSD3
========================================================================
Train the longitudinal-MAMBA SURVIVAL head over the ablation grid and export, per config × seed, a
per-patient predicted survival curve on a fixed dense time grid (+ event/time/Label_T4). The downstream
02_survival_comparison.py (env survml) scores these via the reused 02l.evaluate() and tables them against
the 02l Cox/Weibull-AFT/RSF tabular baselines.

Sweeps (one-factor-at-a-time; see plan):
  A (backbone/training/time/align around the Weibull-piecewise head):
     default(mamba1,frozen,append,align off), mamba2, finetune, prepend, align-on, +meanpool/+gru controls
  B (survival-head parametrisation around mamba1-frozen):
     weibull_aft, weibull_piecewise(=default), cox, soden(guarded)

Inputs : encoder_outputs_no_cdr_post_exclusion_longitudinal/.../seed_{s}/full_ft/embeddings/embeddings.parquet
         + the longitudinal master (patient-level survival labels).
Outputs: outputs/survival/<config>/seed_{s}/predictions.csv  (Patient_ID,split,event,time,Label_T4,S_000..S_060)
         + meta.json (config, val_nll, grid_times).  NO model weights kept.

Usage:
  conda run -n clinical python integration/T4/mamba_long_clinical/01_train_survival_mamba.py \
      --seeds 0 1 2 [--only A_default] [--device cuda]
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

EMB_ROOT = Path(r"D:\ADNI_BIDS_project\derivatives\encoder_outputs_no_cdr_post_exclusion_longitudinal"
                r"\BioClinical-ModernBERT-large\T2_long_multiclass")
MASTER = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion"
              r"\verbose\longitudinal\master_clinical_verbose.csv")
OUT = HERE / "outputs" / "survival"
GRID_TIMES = np.round(np.linspace(0.0, 15.0, 61), 4)            # dense grid for S export (interp later)

# config = (encoder kwargs, head, align_lambda, lr).  lr: frozen 1e-3 / finetune 1e-4.
def CONFIGS():
    enc = lambda **kw: {**dict(backbone="mamba1", frozen=True, time_concat="append"), **kw}
    return {
        # ── Sweep A (Weibull-piecewise head) ──
        "A_default_mamba1_frozen": (enc(), "weibull_piecewise", 0.0, 1e-3),
        "A_q1_mamba2":             (enc(backbone="mamba2"), "weibull_piecewise", 0.0, 1e-3),
        "A_q2_finetune":           (enc(frozen=False), "weibull_piecewise", 0.0, 1e-4),
        "A_q3_time_prepend":       (enc(time_concat="prepend"), "weibull_piecewise", 0.0, 1e-3),
        "A_q4_align":              (enc(), "weibull_piecewise", 0.5, 1e-3),
        "A_ctrl_meanpool":         (enc(backbone="meanpool"), "weibull_piecewise", 0.0, 1e-3),
        "A_ctrl_gru":              (enc(backbone="gru"), "weibull_piecewise", 0.0, 1e-3),
        # ── Sweep B (head parametrisation; mamba1 frozen append) ──
        "B_weibull_aft":           (enc(), "weibull_aft", 0.0, 1e-3),
        "B_cox":                   (enc(), "cox", 0.0, 1e-3),
        "B_soden":                 (enc(), "soden", 0.0, 1e-3),
        # B_weibull_piecewise == A_default_mamba1_frozen (not duplicated)
    }


def run_config(name, cfg, seed, master_labels, device):
    enc_kw, head_kind, align_lambda, lr = cfg
    lib.seed_everything(seed)
    emb_parquet = EMB_ROOT / f"seed_{seed}" / "full_ft" / "embeddings" / "embeddings.parquet"
    seq = lib.assemble_sequences(emb_parquet, master_labels, "survival")

    encoder = lib.make_encoder(emb_dim=seq["emb_dim"], **enc_kw)
    head = lib.make_head(head_kind, encoder.d_out)               # raises ImportError for soden w/o torchdiffeq
    encoder, head, val_nll = lib.train_survival(
        encoder, head, seq, device=device, lr=lr, align_lambda=align_lambda)

    # export per-patient S on GRID_TIMES (all splits — train needed for IPCW reference downstream)
    import pandas as pd
    idx = np.arange(len(seq["pid"]))
    S = lib.predict_survival_curves(encoder, head, seq, idx, GRID_TIMES, device=device)
    df = pd.DataFrame({"Patient_ID": seq["pid"], "split": seq["split"],
                       "event": seq["event"], "time": seq["time"],
                       "Label_T4": np.where(seq["y"] >= 0, seq["y"], np.nan)})
    for j in range(len(GRID_TIMES)):
        df[f"S_{j:03d}"] = S[:, j]
    out_dir = OUT / name / f"seed_{seed}"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_dir / "predictions.parquet", index=False)
    json.dump({"config": name, "seed": seed, "head": head_kind, "encoder": enc_kw,
               "align_lambda": align_lambda, "lr": lr, "val_nll": float(val_nll),
               "grid_times": GRID_TIMES.tolist(), "n": int(len(idx))},
              open(out_dir / "meta.json", "w"), indent=2)
    print(f"  [{name} seed{seed}] val_nll={val_nll:.4f}  exported N={len(idx)} -> {out_dir}")


def main():
    global EMB_ROOT, MASTER, OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--only", type=str, nargs="*", default=None, help="subset of config names")
    ap.add_argument("--device", type=str, default=None)
    ap.add_argument("--emb_root", type=str, default=str(EMB_ROOT),
                    help="dir containing seed_*/full_ft/embeddings/embeddings.parquet")
    ap.add_argument("--master", type=str, default=str(MASTER))
    ap.add_argument("--out_dir", type=str, default=str(OUT),
                    help="where to write predictions (default repo outputs/ for LOCAL; pass a "
                         "hpc-work path on CSD3 so nothing lands in the repo/script folder)")
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
            except ImportError as e:
                print(f"  [SKIP {name} seed{seed}] {e}")
            except Exception:
                print(f"  [FAIL {name} seed{seed}]\n{traceback.format_exc()}")
    print("Done.")


if __name__ == "__main__":
    main()
