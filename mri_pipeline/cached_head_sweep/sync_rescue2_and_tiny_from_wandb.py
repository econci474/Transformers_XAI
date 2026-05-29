"""
sync_rescue2_and_tiny_from_wandb.py
====================================
Tiny dual-purpose shim writer for two early-stage sweeps that currently
have only 1 run in W&B (T1_binary / seed_0 each):

    agms3d_vanilla_rescue2  -> agms3d_outputs_rescue2/AGMS3DCNN_vanilla_slim/<task>/seed_<n>/metrics.json
    vit_tiny_scratch        -> vit_tiny_baseline/aug_<aug>/ViT_T_scratch/<task>/seed_<n>/<strategy>/metrics.json

Both sweeps collapsed on the seed shown (val_bacc=0.5, all-positives
prediction). Writing the shim surfaces that single-cell evidence in the
cross-model table so the row exists with an explicit "1 seed" annotation,
rather than the row silently missing.

When the user uploads more runs to W&B (or rsyncs canonical metrics.json
from CSD3), re-run this script -- the _provenance flag lets it overwrite
its own previous shims safely without touching real metrics.json files.

Usage:
    python mri_pipeline/cached_head_sweep/sync_rescue2_and_tiny_from_wandb.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
CSV_DIR = REPO_ROOT / "mri_pipeline" / "cached_head_sweep"

AGMS3D_OUT = Path(r"D:/ADNI_BIDS_project/derivatives/agms3d_outputs_rescue2/AGMS3DCNN_vanilla_slim")
TINY_OUT = Path(r"D:/ADNI_BIDS_project/derivatives/vit_tiny_baseline")


def _s(row, key, default=None):
    v = row.get(f"summary/{key}")
    if pd.isna(v):
        return default
    if isinstance(v, (bool, str)):
        return v
    return float(v)


def _agms3d_row_to_metrics(row: pd.Series) -> dict:
    task = row["config/task"]
    is_multi = "multiclass" in task
    test_metrics = {
        "accuracy":     _s(row, "test/accuracy"),
        "balanced_acc": _s(row, "test/balanced_acc"),
        "precision":    _s(row, "test/precision_macro" if is_multi else "test/precision"),
        "recall":       _s(row, "test/recall_macro" if is_multi else "test/recall"),
        "sensitivity":  _s(row, "test/sensitivity"),
        "specificity":  _s(row, "test/specificity"),
        "f1":           _s(row, "test/macro_f1" if is_multi else "test/f1"),
        "auc_roc":      _s(row, "test/auc_roc_ovr" if is_multi else "test/auc_roc"),
        "auc_pr":       _s(row, "test/auc_pr_macro" if is_multi else "test/auc_pr"),
    }
    test_metrics_subject = {
        "accuracy":     _s(row, "test_subject/accuracy"),
        "balanced_acc": _s(row, "test_subject/balanced_acc"),
        "precision":    _s(row, "test_subject/precision_macro" if is_multi
                                else "test_subject/precision"),
        "recall":       _s(row, "test_subject/recall_macro" if is_multi
                                else "test_subject/recall"),
        "sensitivity":  _s(row, "test_subject/sensitivity"),
        "specificity":  _s(row, "test_subject/specificity"),
        "f1":           _s(row, "test_subject/macro_f1" if is_multi
                                else "test_subject/f1"),
        "auc_roc":      _s(row, "test_subject/auc_roc_ovr" if is_multi
                                else "test_subject/auc_roc"),
        "auc_pr":       _s(row, "test_subject/auc_pr_macro" if is_multi
                                else "test_subject/auc_pr"),
    }
    return {
        "config": {
            "model_id":     row.get("config/model", "AGMS3DCNN_vanilla_slim"),
            "model_kind":   "vanilla",       # parsed as `variant`
            "backbone":     row.get("config/backbone"),
            "head":         row.get("config/head"),
            "strong_aug":   bool(row.get("config/strong_aug")),
            "task":         task,
            "seed":         int(row["config/seed"]),
            "num_labels":   int(row["config/num_labels"]),
            "epochs":       int(row["config/epochs"]),
            "best_epoch":   int(_s(row, "best_epoch")) if _s(row, "best_epoch") is not None else None,
            "best_val_balanced_acc": _s(row, "best_val_balanced_acc"),
            "lr":           float(row["config/lr"]),
            "batch_size":   int(row["config/batch_size"]),
            "label_smoothing": float(row["config/label_smoothing"]),
            "augment":      "flips+strong",   # rescue2 implies --strong_aug
            "_provenance":  "shim_from_wandb_summary",
        },
        "test_metrics": test_metrics,
        "test_metrics_subject": test_metrics_subject,
        "test_diagnostics": {},
    }


def _vit_tiny_row_to_metrics(row: pd.Series) -> dict:
    task = row["config/task"]
    is_multi = "multiclass" in task
    test_metrics = {
        "accuracy":     _s(row, "test/accuracy"),
        "balanced_acc": _s(row, "test/balanced_acc"),
        "precision":    _s(row, "test/precision_macro" if is_multi else "test/precision"),
        "recall":       _s(row, "test/recall_macro" if is_multi else "test/recall"),
        "sensitivity":  _s(row, "test/sensitivity"),
        "specificity":  _s(row, "test/specificity"),
        "f1":           _s(row, "test/macro_f1" if is_multi else "test/f1"),
        "auc_roc":      _s(row, "test/auc_roc_ovr" if is_multi else "test/auc_roc"),
        "auc_pr":       _s(row, "test/auc_pr_macro" if is_multi else "test/auc_pr"),
    }
    test_metrics_subject = {
        "accuracy":     _s(row, "test_subject/accuracy"),
        "balanced_acc": _s(row, "test_subject/balanced_acc"),
        "precision":    _s(row, "test_subject/precision_macro" if is_multi
                                else "test_subject/precision"),
        "recall":       _s(row, "test_subject/recall_macro" if is_multi
                                else "test_subject/recall"),
        "sensitivity":  _s(row, "test_subject/sensitivity"),
        "specificity":  _s(row, "test_subject/specificity"),
        "f1":           _s(row, "test_subject/macro_f1" if is_multi
                                else "test_subject/f1"),
        "auc_roc":      _s(row, "test_subject/auc_roc_ovr" if is_multi
                                else "test_subject/auc_roc"),
        "auc_pr":       _s(row, "test_subject/auc_pr_macro" if is_multi
                                else "test_subject/auc_pr"),
    }
    return {
        "config": {
            "model_id":     "ViT_T_scratch",
            "vit_size":     row.get("config/vit_size"),
            "task":         task,
            "seed":         int(row["config/seed"]),
            "strategy":     row["config/strategy"],
            "augment":      row["config/augment"],
            "lr":           float(row["config/lr"]),
            "llrd_gamma":   float(row["config/llrd_gamma"]),
            "weight_decay": float(row["config/weight_decay"]),
            "batch_size":   int(row["config/batch_size"]),
            "grad_accum_steps": int(row["config/grad_accum_steps"]),
            "warmup_epochs":    int(row["config/warmup_epochs"]),
            "patience":     int(row["config/patience"]),
            "epochs":       int(row["config/epochs"]),
            "drop_path_rate":   float(row["config/drop_path_rate"]),
            "attn_dropout":     float(row["config/attn_dropout"]),
            "label_smoothing":  float(row["config/label_smoothing"]),
            "num_labels":   int(row["config/num_labels"]),
            "n_train":      int(row["config/n_train"]),
            "n_val":        int(row["config/n_val"]),
            "n_test":       int(row["config/n_test"]),
            "best_val_balanced_acc": _s(row, "best_val_balanced_acc"),
            "_provenance":  "shim_from_wandb_summary",
        },
        "test_metrics": test_metrics,
        "test_metrics_subject": test_metrics_subject,
        "test_diagnostics": {},
    }


def _maybe_write(out_file: Path, payload: dict, label: str):
    if out_file.exists():
        with open(out_file) as f:
            existing = json.load(f)
        if existing.get("config", {}).get("_provenance") != "shim_from_wandb_summary":
            print(f"  [skip-real]  {label}  (canonical metrics.json present)")
            return False
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file, "w") as f:
        json.dump(payload, f, indent=2)
    val = payload["config"].get("best_val_balanced_acc")
    val_s = f"{val:.3f}" if val is not None else "  -  "
    test_bacc = payload["test_metrics"].get("balanced_acc")
    test_s = f"{test_bacc:.3f}" if test_bacc is not None else "  -  "
    print(f"  [wrote-shim] {label}  val={val_s}  test_bacc={test_s}")
    return True


def main():
    # --- AGMS3D rescue2 -----------------------------------------------------
    csv = CSV_DIR / "wandb_agms3d_vanilla_rescue2.csv"
    if csv.exists():
        df = pd.read_csv(csv)
        df = df[df["state"] == "finished"]
        print(f"\nagms3d_vanilla_rescue2: {len(df)} finished run(s)")
        for _, row in df.iterrows():
            task = row["config/task"]; seed = int(row["config/seed"])
            out_file = AGMS3D_OUT / task / f"seed_{seed}" / "metrics.json"
            _maybe_write(out_file, _agms3d_row_to_metrics(row),
                         f"AGMS3D-r2/{task}/seed_{seed}")
    else:
        print(f"\n[miss] {csv.name} -- run 06_fetch_wandb_tables first.")

    # --- ViT-Tiny scratch ---------------------------------------------------
    csv = CSV_DIR / "wandb_vit_tiny_scratch.csv"
    if csv.exists():
        df = pd.read_csv(csv)
        df = df[df["state"] == "finished"]
        print(f"\nvit_tiny_scratch: {len(df)} finished run(s)")
        for _, row in df.iterrows():
            task = row["config/task"]; seed = int(row["config/seed"])
            aug = row["config/augment"]; strat = row["config/strategy"]
            out_file = (TINY_OUT / f"aug_{aug}" / "ViT_T_scratch"
                        / task / f"seed_{seed}" / strat / "metrics.json")
            _maybe_write(out_file, _vit_tiny_row_to_metrics(row),
                         f"ViT-Tiny/aug_{aug}/{task}/seed_{seed}/{strat}")
    else:
        print(f"\n[miss] {csv.name} -- run 06_fetch_wandb_tables first.")


if __name__ == "__main__":
    main()
