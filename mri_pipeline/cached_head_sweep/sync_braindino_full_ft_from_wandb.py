"""
sync_braindino_full_ft_from_wandb.py
====================================
Synthesise per-run `metrics.json` shims for the BrainDINO full end-to-end
fine-tune sweep (W&B project `braindino_full_ft`). Same trainer as the
LoRA and frozen-aug sweeps -- val_* logged to W&B, test_* only in the
canonical metrics.json on CSD3.

Output layout matches the aggregator's BrainDINO full_ft glob:

    braindino_outputs/ft/aug_<aug>/BrainDINO_vitb16_full_ft/<task>/seed_<n>/metrics.json

Usage:
    python mri_pipeline/cached_head_sweep/sync_braindino_full_ft_from_wandb.py
    python mri_pipeline/cached_head_sweep/sync_braindino_full_ft_from_wandb.py --augment stochastic
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
WANDB_CSV = REPO_ROOT / "mri_pipeline" / "cached_head_sweep" / "wandb_braindino_full_ft.csv"
DERIV_ROOT = Path(r"D:/ADNI_BIDS_project/derivatives/braindino_outputs/ft")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--augment", default=None,
                   help="Restrict to one augment (none|stochastic|plus_original). "
                        "Default: all.")
    p.add_argument("--task", default=None,
                   help="Restrict to a single task. Default: all 5 tasks.")
    return p.parse_args()


def _row_to_metrics(row: pd.Series) -> dict:
    def s(key, default=None):
        v = row.get(f"summary/{key}")
        if pd.isna(v):
            return default
        return v

    val_TP = s("val_TP"); val_FP = s("val_FP")
    val_TN = s("val_TN"); val_FN = s("val_FN")

    return {
        "config": {
            "model_id":              "BrainDINO_vitb16",
            "model_kind":            "full_ft",
            "strategy":              row["config/strategy"],
            "augment":               row["config/augment"],
            "task":                  row["config/task"],
            "seed":                  int(row["config/seed"]),
            "num_labels":            int(row["config/num_labels"]),
            "epochs":                int(row["config/epochs"]),
            "best_epoch":            int(s("epoch")) if s("epoch") is not None else None,
            "best_val_balanced_acc": float(s("best_val_balanced_acc"))
                                     if s("best_val_balanced_acc") is not None else None,
            "lr":                    float(row["config/lr"]),
            "weight_decay":          float(row["config/weight_decay"]),
            "drop_rate":             float(row["config/drop_rate"]),
            "label_smoothing":       float(row["config/label_smoothing"]),
            "patience":              int(row["config/patience"]),
            "n_train":               int(row["config/n_train"]),
            "n_val":                 int(row["config/n_val"]),
            "n_test":                int(row["config/n_test"]),
            "_provenance":           "shim_from_wandb_summary",
        },
        "val_metrics": {
            # checkpoint-selection metric (BEST seen across training, not last
            # epoch -- summary/val_bacc is last; summary/best_val_balanced_acc
            # is best).
            "balanced_acc": float(s("best_val_balanced_acc"))
                            if s("best_val_balanced_acc") is not None else None,
            "loss":         float(s("val_loss")) if s("val_loss") is not None else None,
            "auc_roc":      float(s("val_auc"))  if s("val_auc") is not None else None,
            "f1":           float(s("val_f1"))   if s("val_f1") is not None else None,
            "precision":    float(s("val_prec")) if s("val_prec") is not None else None,
            "recall":       float(s("val_recall")) if s("val_recall") is not None else None,
            "sensitivity":  float(s("val_sens")) if s("val_sens") is not None else None,
            "specificity":  float(s("val_spec")) if s("val_spec") is not None else None,
            "TP":           int(val_TP) if val_TP is not None else None,
            "FP":           int(val_FP) if val_FP is not None else None,
            "TN":           int(val_TN) if val_TN is not None else None,
            "FN":           int(val_FN) if val_FN is not None else None,
        },
        "test_metrics": {},
        "test_metrics_subject": {},
        "test_diagnostics": {},
    }


def main():
    args = parse_args()
    if not WANDB_CSV.exists():
        raise SystemExit(f"W&B CSV not found: {WANDB_CSV}. Run "
                         "06_fetch_wandb_tables.py --projects braindino_full_ft first.")

    df = pd.read_csv(WANDB_CSV)
    df = df[df["state"] == "finished"]
    df = df[df["summary/best_val_balanced_acc"].notna()]
    if args.augment:
        df = df[df["config/augment"] == args.augment]
    if args.task:
        df = df[df["config/task"] == args.task]

    df = df.sort_values("summary/best_val_balanced_acc", ascending=False)
    df = df.drop_duplicates(
        subset=["config/task", "config/seed", "config/augment"], keep="first")

    print(f"Filter   : {'augment=' + args.augment if args.augment else 'all augments'}"
          + (f"  task={args.task}" if args.task else ""))
    print(f"Out tree : {DERIV_ROOT}")
    print(f"Kept {len(df)} unique (task, seed, augment) runs after dedup.")

    n_written = n_skipped = 0
    for _, row in df.iterrows():
        aug = row["config/augment"]
        task = row["config/task"]
        seed = int(row["config/seed"])
        out_dir = DERIV_ROOT / f"aug_{aug}" / "BrainDINO_vitb16_full_ft" / task / f"seed_{seed}"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_file = out_dir / "metrics.json"

        if out_file.exists():
            with open(out_file) as f:
                existing = json.load(f)
            if existing.get("config", {}).get("_provenance") != "shim_from_wandb_summary":
                n_skipped += 1
                continue

        payload = _row_to_metrics(row)
        with open(out_file, "w") as f:
            json.dump(payload, f, indent=2)
        n_written += 1
        val = payload["config"]["best_val_balanced_acc"]
        val_s = f"{val:.3f}" if val is not None else "  -  "
        rel = out_file.relative_to(DERIV_ROOT)
        print(f"  [wrote-shim] {rel}  val={val_s}")

    print(f"\nWrote {n_written} shim(s); skipped {n_skipped} canonical metrics.json")


if __name__ == "__main__":
    main()
