"""
sync_agms3d_vanilla_from_wandb.py
=================================
Synthesise per-run `metrics.json` shims for the AG-MS3D-vanilla rescue
sweep, from the W&B summary CSV pulled by
`mri_pipeline/cached_head_sweep/06_fetch_wandb_tables.py`.

The trainer writes metrics.json on the HPC; for now those files are
still on CSD3 and only the W&B summary is locally available. This
script writes a minimal-but-aggregator-compatible metrics.json at the
canonical local path:

    <root>/agms3d_outputs/AGMS3DCNN_vanilla/<task>/seed_<s>/metrics.json

so that `05_aggregate_mri_results.py` / `06_render_cross_model_table.py`
can pick up the rescue results immediately. Rsyncing the real metrics.json
later will overwrite these shims with the canonical files.

Filters out pre-rescue rows (lr != 0.003) so only the rescue config is
materialised on disk.

Usage:
    python mri_pipeline/3d_cnn_vit/sync_agms3d_vanilla_from_wandb.py
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
WANDB_CSV = REPO_ROOT / "mri_pipeline" / "cached_head_sweep" / "wandb_agms3d_vanilla_rescue.csv"
DERIV_ROOT = Path(r"D:/ADNI_BIDS_project/derivatives/agms3d_outputs/AGMS3DCNN_vanilla")

# Only the rescue config (lr=3e-3, batch_size=8) — pre-rescue runs are
# kept in W&B for history but mustn't appear in the comparison table.
RESCUE_LR = 0.003
RESCUE_BS = 8


def _row_to_metrics(row: pd.Series) -> dict:
    """Map a flattened W&B run row into a metrics.json dict matching the
    schema in 05_aggregate_mri_results.py's `read_run`."""
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
        "recall":       s("test/recall_macro" if is_multi else "test/recall"),
        "sensitivity":  s("test/sensitivity"),
        "specificity":  s("test/specificity"),
        "f1":           s("test/macro_f1" if is_multi else "test/f1"),
        "auc_roc":      s("test/auc_roc_ovr" if is_multi else "test/auc_roc"),
        "auc_pr":       s("test/auc_pr_macro" if is_multi else "test/auc_pr"),
    }
    if is_multi:
        test_metrics["macro_f1"]    = s("test/macro_f1")
        test_metrics["auc_roc_ovr"] = s("test/auc_roc_ovr")

    test_metrics_subject = {
        "accuracy":     s("test_subject/accuracy"),
        "balanced_acc": s("test_subject/balanced_acc"),
        "precision":    s("test_subject/precision_macro" if is_multi
                          else "test_subject/precision"),
        "recall":       s("test_subject/recall_macro" if is_multi
                          else "test_subject/recall"),
        "sensitivity":  s("test_subject/sensitivity"),
        "specificity":  s("test_subject/specificity"),
        "f1":           s("test_subject/macro_f1" if is_multi
                          else "test_subject/f1"),
        "auc_roc":      s("test_subject/auc_roc_ovr" if is_multi
                          else "test_subject/auc_roc"),
        "auc_pr":       s("test_subject/auc_pr_macro" if is_multi
                          else "test_subject/auc_pr"),
    }

    return {
        "config": {
            "model_id":    "AGMS3DCNN_vanilla",
            "model_kind":  "vanilla",          # parsed as `variant`
            "task":        task,
            "seed":        int(row["config/seed"]),
            "num_labels":  int(row["config/num_labels"]),
            "epochs":      int(row["config/epochs"]),
            "best_epoch":  int(row["summary/best_epoch"])
                          if pd.notna(row.get("summary/best_epoch")) else None,
            "best_val_balanced_acc": float(row["summary/val_bacc"])
                          if pd.notna(row.get("summary/val_bacc")) else None,
            "lr":          float(row["config/lr"]),
            "batch_size":  int(row["config/batch_size"]),
            "base_filters": int(row["config/base_filters"])
                          if pd.notna(row.get("config/base_filters")) else None,
            "label_smoothing": 0.0,
            "augment":     "-",
            "_provenance": "shim_from_wandb_summary",
        },
        "test_metrics": test_metrics,
        "test_metrics_subject": test_metrics_subject,
        "test_diagnostics": {},   # confusion matrix not in W&B summary;
                                  # _is_degenerate returns False on empty
    }


def main():
    if not WANDB_CSV.exists():
        raise SystemExit(f"W&B CSV not found: {WANDB_CSV}. Run "
                         "06_fetch_wandb_tables.py --projects "
                         "agms3d_vanilla_rescue first.")

    df = pd.read_csv(WANDB_CSV)
    mask = (df["config/lr"] == RESCUE_LR) & (df["config/batch_size"] == RESCUE_BS)
    df = df[mask].sort_values(["config/task", "config/seed"]).reset_index(drop=True)
    print(f"Found {len(df)} rescue runs (lr={RESCUE_LR}, bs={RESCUE_BS}).")

    n_written = n_skipped = 0
    for _, row in df.iterrows():
        task = row["config/task"]
        seed = int(row["config/seed"])
        out_dir = DERIV_ROOT / task / f"seed_{seed}"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_file = out_dir / "metrics.json"
        if out_file.exists():
            # Don't overwrite a real metrics.json (e.g. one rsynced from CSD3).
            with open(out_file) as f:
                existing = json.load(f)
            if existing.get("config", {}).get("_provenance") != "shim_from_wandb_summary":
                print(f"  [skip-real]  {task}/seed_{seed}  (canonical metrics.json present)")
                n_skipped += 1
                continue
        payload = _row_to_metrics(row)
        with open(out_file, "w") as f:
            json.dump(payload, f, indent=2)
        n_written += 1
        print(f"  [wrote-shim] {task}/seed_{seed}  "
              f"val={payload['config']['best_val_balanced_acc']:.3f}  "
              f"test_bacc={payload['test_metrics']['balanced_acc']:.3f}")

    print(f"\nDone. Wrote {n_written} shim(s); skipped {n_skipped} real metrics.json file(s).")
    print(f"Output root: {DERIV_ROOT}")


if __name__ == "__main__":
    main()
