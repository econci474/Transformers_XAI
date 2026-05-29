"""
sync_vit_hi_lr_from_wandb.py
============================
Synthesise per-run `metrics.json` shims for the ViT-MAE75 LLRD-off (γ=1.0)
sweep from W&B project `vit_mae_hi_lr`. The "hi_lr" submit was mis-named:
it didn't change LR, it tested LLRD off (γ=1.0) vs the legacy LLRD on
(γ=0.7). The vit_mae_hi_lr project also covers aug=none and
aug=plus_original variants that the legacy `vit_debugging` project doesn't.

Writes shims at the canonical aggregator-matched paths:

    aug=random:        vit_outputs_hi_lr/ViT_B_mae75/<task>/seed_<n>/<strategy>/metrics.json
    aug=<other>:       vit_outputs_hi_lr/aug_<aug>/ViT_B_mae75/<task>/seed_<n>/<strategy>/metrics.json

(The aggregator's double-nest fallback also picks up
`vit_outputs_hi_lr/vit_outputs_hi_lr/ViT_B_mae75/...` for already-rsynced
files, so this script's shims don't duplicate those.)

The _provenance flag in config marks these as shims; a future rsync of
the canonical metrics.json from CSD3 will overwrite them safely.

Usage:
    python mri_pipeline/cached_head_sweep/sync_vit_hi_lr_from_wandb.py
    python mri_pipeline/cached_head_sweep/sync_vit_hi_lr_from_wandb.py --task T1c_binary
    python mri_pipeline/cached_head_sweep/sync_vit_hi_lr_from_wandb.py --strategy frozen
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
WANDB_CSV = REPO_ROOT / "mri_pipeline" / "cached_head_sweep" / "wandb_vit_mae_hi_lr.csv"
DERIV_BASE = Path(r"D:/ADNI_BIDS_project/derivatives/vit_outputs_hi_lr")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--llrd_gamma", type=float, default=1.0,
                   help="Filter to runs with this llrd_gamma. Default 1.0 "
                        "(the LLRD-off variant that beats LLRD-on on every "
                        "full_ft task).")
    p.add_argument("--strategy", default=None,
                   help="Filter to this strategy (full_ft|frozen). Default: both.")
    p.add_argument("--augment", default=None,
                   help="Filter to this augment (none|random|plus_original). "
                        "Default: all augments present in the CSV.")
    p.add_argument("--task", default=None,
                   help="Restrict to a single task. Default: all 5 tasks.")
    return p.parse_args()


def _row_to_metrics(row: pd.Series) -> dict:
    """Map a flattened W&B run row into the ViT trainer's metrics.json schema."""
    task = row["config/task"]
    is_multi = "multiclass" in task

    def s(key, default=None):
        v = row.get(f"summary/{key}")
        if pd.isna(v):
            return default
        return float(v) if not isinstance(v, (str, bool)) else v

    test_metrics = {
        "accuracy":     s("test/accuracy"),
        "balanced_acc": s("test/balanced_acc"),
        "precision":    s("test/precision_macro" if is_multi else "test/precision"),
        "recall":       s("test/recall_macro"    if is_multi else "test/recall"),
        "sensitivity":  s("test/sensitivity"),
        "specificity":  s("test/specificity"),
        "f1":           s("test/macro_f1"        if is_multi else "test/f1"),
        "auc_roc":      s("test/auc_roc_ovr"     if is_multi else "test/auc_roc"),
        "auc_pr":       s("test/auc_pr_macro"    if is_multi else "test/auc_pr"),
    }
    if is_multi:
        test_metrics["macro_f1"]    = s("test/macro_f1")
        test_metrics["auc_roc_ovr"] = s("test/auc_roc_ovr")

    test_metrics_subject = {
        "accuracy":     s("test_subject/accuracy"),
        "balanced_acc": s("test_subject/balanced_acc"),
        "precision":    s("test_subject/precision_macro" if is_multi
                          else "test_subject/precision"),
        "recall":       s("test_subject/recall_macro"    if is_multi
                          else "test_subject/recall"),
        "sensitivity":  s("test_subject/sensitivity"),
        "specificity":  s("test_subject/specificity"),
        "f1":           s("test_subject/macro_f1"        if is_multi
                          else "test_subject/f1"),
        "auc_roc":      s("test_subject/auc_roc_ovr"     if is_multi
                          else "test_subject/auc_roc"),
        "auc_pr":       s("test_subject/auc_pr_macro"    if is_multi
                          else "test_subject/auc_pr"),
    }

    return {
        "config": {
            "model_id":              "ViT_B_mae75",
            "task":                  task,
            "seed":                  int(row["config/seed"]),
            "strategy":              row["config/strategy"],
            "augment":               row["config/augment"],
            "lr":                    float(row["config/lr"]),
            "llrd_gamma":            float(row["config/llrd_gamma"]),
            "weight_decay":          float(row["config/weight_decay"]),
            "batch_size":            int(row["config/batch_size"]),
            "grad_accum_steps":      int(row["config/grad_accum_steps"]),
            "warmup_epochs":         int(row["config/warmup_epochs"]),
            "patience":              int(row["config/patience"]),
            "epochs":                int(row["config/epochs"]),
            "drop_path_rate":        float(row["config/drop_path_rate"]),
            "attn_dropout":          float(row["config/attn_dropout"]),
            "label_smoothing":       float(row["config/label_smoothing"]),
            "num_labels":            int(row["config/num_labels"]),
            "n_train":               int(row["config/n_train"]),
            "n_val":                 int(row["config/n_val"]),
            "n_test":                int(row["config/n_test"]),
            "best_val_balanced_acc": float(row["summary/best_val_balanced_acc"])
                                     if pd.notna(row.get("summary/best_val_balanced_acc")) else None,
            "_provenance":           "shim_from_wandb_summary",
        },
        "test_metrics":         test_metrics,
        "test_metrics_subject": test_metrics_subject,
        "test_diagnostics":     {},   # CM not in summary (it's a wandb table)
    }


def _out_dir(row: pd.Series) -> Path:
    """Mirror the 04h submit's aug-baked OUT_DIR layout:
        aug=random        -> vit_outputs_hi_lr/ViT_B_mae75/<task>/seed_<n>/<strat>/
        aug=<other>       -> vit_outputs_hi_lr/aug_<aug>/ViT_B_mae75/<task>/seed_<n>/<strat>/
    Matches the aggregator's globs in 05_aggregate_mri_results.py."""
    aug = row["config/augment"]
    task = row["config/task"]
    seed = int(row["config/seed"])
    strat = row["config/strategy"]
    if aug == "random":
        return DERIV_BASE / "ViT_B_mae75" / task / f"seed_{seed}" / strat
    return DERIV_BASE / f"aug_{aug}" / "ViT_B_mae75" / task / f"seed_{seed}" / strat


def main():
    args = parse_args()
    if not WANDB_CSV.exists():
        raise SystemExit(f"W&B CSV not found: {WANDB_CSV}. Run "
                         "06_fetch_wandb_tables.py --projects vit_mae_hi_lr first.")

    df = pd.read_csv(WANDB_CSV)
    df = df[df["state"] == "finished"]
    df = df[df["config/llrd_gamma"] == args.llrd_gamma]
    # Drop runs where best_val_balanced_acc is NaN (incomplete / failed early).
    df = df[df["summary/best_val_balanced_acc"].notna()]
    if args.strategy:
        df = df[df["config/strategy"] == args.strategy]
    if args.augment:
        df = df[df["config/augment"] == args.augment]
    if args.task:
        df = df[df["config/task"] == args.task]

    print(f"Filter   : llrd_gamma={args.llrd_gamma}"
          + (f"  strategy={args.strategy}" if args.strategy else "")
          + (f"  augment={args.augment}" if args.augment else "")
          + (f"  task={args.task}" if args.task else ""))
    print(f"Out base : {DERIV_BASE}")
    print(f"Found {len(df)} finished runs")

    # If multiple runs share the same (task, seed, strategy, augment) -- e.g.
    # a re-submit after a hang -- keep the one with the highest val_bacc.
    df = df.sort_values("summary/best_val_balanced_acc", ascending=False)
    df = df.drop_duplicates(
        subset=["config/task", "config/seed", "config/strategy",
                "config/augment"], keep="first")

    n_written = n_skipped = 0
    for _, row in df.iterrows():
        out_dir = _out_dir(row)
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
        bacc = payload["test_metrics"].get("balanced_acc")
        bacc_s = f"{bacc:.3f}" if bacc is not None else "  -  "
        val = payload["config"]["best_val_balanced_acc"]
        val_s = f"{val:.3f}" if val is not None else "  -  "
        rel = out_file.relative_to(DERIV_BASE)
        print(f"  [wrote-shim] {rel}  val={val_s}  test_bacc={bacc_s}")

    print(f"\nWrote {n_written} shim(s); skipped {n_skipped} canonical metrics.json")


if __name__ == "__main__":
    main()
