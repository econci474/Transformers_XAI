"""
04_finalize_winner.py
=====================
After a `--metrics_only` wider sweep (e.g. 04_braindino_head_sweep_T2_wide_*),
the losing cells have only metrics.json and the WINNER has no checkpoint. This
script:
  1. selects the HP-winner under --sweep-root = max MEAN val-bACC across seeds
     (same rule as 05b_select_best_hp_per_task.py);
  2. reads that winner's metrics.json `config` to recover the EXACT HPs
     (head / lr / weight_decay / drop_rate / label_smoothing / loss /
      focal_gamma / standardize) — robust vs. parsing the dir-name slug;
  3. re-trains the winner on each --seeds WITHOUT --metrics_only (so it writes
     config + metrics + best_model.pt) into --out-dir (a fresh winner tree, so
     the trainer does not self-skip the sweep's metrics-only cells).

Then extract per-patient probs with:
  python mri_pipeline/05_extract_cached_head_probs.py --arch braindino \
      --tasks T2_multiclass --heads-root <out-dir>

Saves ONLY the best checkpoint's artifacts (config/metrics/best_model.pt) +
embeddings (via the extractor), and the wandb runs are offline-syncable.

Example:
  python mri_pipeline/cached_head_sweep/04_finalize_winner.py \
    --model_name BrainDINO --embed_dim 768 --task T2_multiclass \
    --embeddings_pt .../cached_embeddings/braindino/braindino_pooled.pt \
    --sweep-root .../braindino_outputs/aug_none_T2_wide/BrainDINO_vitb16_frozen_cached \
    --out-dir   .../braindino_outputs/aug_none_T2_wide_winner/BrainDINO_vitb16_frozen_cached \
    --data_dir .../tabular/baseline --matched_labels_csv .../master_*.csv --long all
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

THIS_DIR = Path(__file__).resolve().parent
TRAINER = THIS_DIR / "04_head_finetune_from_embeddings.py"


def select_winner(sweep_root: Path, task: str):
    """Return (hp_dir_name, mean_val_bacc, config_dict) for the HP with the
    highest mean config.best_val_balanced_acc across seeds."""
    rows, cfg_by_hp = [], {}
    for mp in glob.glob(str(sweep_root / task / "seed_*" / "*" / "metrics.json")):
        try:
            cfg = json.load(open(mp)).get("config", {})
        except Exception:
            continue
        hp = os.path.basename(os.path.dirname(mp))
        rows.append((hp, cfg.get("best_val_balanced_acc")))
        cfg_by_hp.setdefault(hp, cfg)             # any seed's config carries the HPs
    if not rows:
        raise FileNotFoundError(f"No metrics.json under {sweep_root/task}")
    df = pd.DataFrame(rows, columns=["hp", "val"]).dropna()
    mean = df.groupby("hp")["val"].mean()
    hp = mean.idxmax()
    return hp, float(mean.max()), cfg_by_hp[hp]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model_name", required=True)
    ap.add_argument("--embeddings_pt", required=True)
    ap.add_argument("--embed_dim", type=int, required=True)
    ap.add_argument("--task", required=True)
    ap.add_argument("--sweep-root", required=True, help="metrics-only sweep tree.")
    ap.add_argument("--out-dir", required=True, help="fresh winner tree to write to.")
    ap.add_argument("--data_dir", required=True)
    ap.add_argument("--matched_labels_csv", required=True)
    ap.add_argument("--long", default="all")
    ap.add_argument("--seeds", default="0,1,2")
    ap.add_argument("--wandb", action="store_true")
    ap.add_argument("--wandb_project", default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    hp, mean_val, cfg = select_winner(Path(args.sweep_root), args.task)
    print(f"[finalize] winner HP = {hp}  (mean val-bACC {mean_val:.4f})")
    print(f"[finalize] config: head={cfg.get('head','mlp')} lr={cfg.get('lr')} "
          f"wd={cfg.get('weight_decay')} drop={cfg.get('drop_rate')} "
          f"ls={cfg.get('label_smoothing')} loss={cfg.get('loss','ce_weighted')} "
          f"focal_gamma={cfg.get('focal_gamma')} standardize={cfg.get('standardize',False)}")

    seeds = [int(s) for s in args.seeds.split(",")]
    for seed in seeds:
        cmd = [sys.executable, str(TRAINER),
               "--model_name", args.model_name,
               "--embeddings_pt", args.embeddings_pt,
               "--embed_dim", str(args.embed_dim),
               "--task", args.task, "--seed", str(seed),
               "--head", str(cfg.get("head", "mlp")),
               "--lr", str(cfg.get("lr", 1e-4)),
               "--weight_decay", str(cfg.get("weight_decay", 1e-5)),
               "--drop_rate", str(cfg.get("drop_rate", 0.2)),
               "--label_smoothing", str(cfg.get("label_smoothing", 0.0)),
               "--loss", str(cfg.get("loss", "ce_weighted")),
               "--long", args.long,
               "--matched_labels_csv", args.matched_labels_csv,
               "--data_dir", args.data_dir,
               "--out_dir", args.out_dir]
        if cfg.get("loss") == "focal" and cfg.get("focal_gamma") is not None:
            cmd += ["--focal_gamma", str(cfg["focal_gamma"])]
        if cfg.get("standardize"):
            cmd += ["--standardize"]
        if args.wandb:
            cmd += ["--wandb"]
            if args.wandb_project:
                cmd += ["--wandb_project", args.wandb_project]
        print("  +", " ".join(cmd))
        if not args.dry_run:
            subprocess.run(cmd, check=True)

    print(f"[finalize] done. Winner re-trained (with best_model.pt) under {args.out_dir}")

    # ── Aggregate the winner's TEST metrics across the re-trained seeds ──
    if not args.dry_run:
        summarize_winner(Path(args.out_dir), args.task)

    print("[finalize] now extract probs:\n"
          f"  python mri_pipeline/05_extract_cached_head_probs.py --arch "
          f"{args.model_name.lower().replace('-','_').replace('vit_mae75','vit_mae')} "
          f"--tasks {args.task} --heads-root {args.out_dir}")


def summarize_winner(out_dir: Path, task: str):
    """Read the winner-tree per-seed metrics.json and report mean±std TEST metrics
    (the headline to compare vs the BrainDINO paper's macro-AUC 0.954). The winner
    tree has exactly one HP dir per seed, so glob it (no slug replication)."""
    import numpy as np
    # test_metrics key aliases across binary/multiclass schemas.
    KEYS = {"test_bACC": ["balanced_acc"],
            "test_macro_AUC": ["auc_roc_ovr", "auc_roc"],
            "test_macro_F1": ["macro_f1", "f1"],
            "test_acc": ["accuracy"]}
    rows, vals = [], {k: [] for k in KEYS}
    for mp in sorted(glob.glob(str(out_dir / task / "seed_*" / "*" / "metrics.json"))):
        d = json.load(open(mp))
        tm = d.get("test_metrics", {}); cfg = d.get("config", {})
        seed = cfg.get("seed")
        row = {"seed": seed, "val_bACC": cfg.get("best_val_balanced_acc"),
               "val_AUC": cfg.get("best_val_auc")}
        for label, aliases in KEYS.items():
            v = next((tm[a] for a in aliases if a in tm), None)
            row[label] = v
            if v is not None:
                vals[label].append(float(v))
        rows.append(row)
    if not rows:
        print("[finalize] no winner metrics.json found to summarise."); return
    summ = pd.DataFrame(rows).sort_values("seed")
    out_csv = out_dir / task / "winner_test_summary.csv"
    summ.to_csv(out_csv, index=False)
    print(f"\n[finalize] WINNER test metrics (n={len(rows)} seeds) — compare macro-AUC vs paper 0.954:")
    for label in KEYS:
        xs = vals[label]
        if xs:
            print(f"    {label:14s} = {np.mean(xs):.4f} ± {np.std(xs):.4f}")
    print(f"[finalize] per-seed table -> {out_csv}")


if __name__ == "__main__":
    main()
