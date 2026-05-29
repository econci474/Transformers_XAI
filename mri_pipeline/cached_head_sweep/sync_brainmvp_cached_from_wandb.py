"""
sync_brainmvp_cached_from_wandb.py
==================================
Synthesise per-run `metrics.json` shims for the cached-head HP sweep projects
(`brainmvp_frozen_cached`, `vit_mae_frozen_cached`, `braindino_frozen_cached`),
from the per-project summary CSV pulled by `06_fetch_wandb_tables.py`. Used
when the canonical metrics.json files are still on CSD3 but you want the
cross-model aggregator to pick the runs up locally.

All three projects share the same trainer (`04_head_finetune_from_embeddings.py`)
so the W&B schema is identical -- only the output tree differs by model.

Important caveat: the cached-head trainer only logs val_* metrics to W&B
(val_bacc, val_auc, val_f1, val_sens/spec, val_TP/FP/TN/FN). It does NOT log
test_metrics to W&B. So the shims have full val coverage but the test_metrics
block is empty -- cells in the TEST cross-model table will read NaN. Once
the canonical metrics.json files are rsynced down from CSD3, this script's
shims should be overwritten (the canonical file has the full test_metrics
block populated; the shim's _provenance flag tells the syncer to overwrite).

Output layout matches each trainer's canonical path so the aggregator globs
pick them up automatically.

Usage:
    # BrainMVP, all tasks
    python mri_pipeline/cached_head_sweep/sync_brainmvp_cached_from_wandb.py
    # BrainMVP, T1c only
    python mri_pipeline/cached_head_sweep/sync_brainmvp_cached_from_wandb.py --task T1c_binary
    # ViT-MAE, all tasks
    python mri_pipeline/cached_head_sweep/sync_brainmvp_cached_from_wandb.py \\
        --project vit_mae_frozen_cached
    # BrainDINO, T1c only
    python mri_pipeline/cached_head_sweep/sync_brainmvp_cached_from_wandb.py \\
        --project braindino_frozen_cached --task T1c_binary
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]

# Per-project metadata: (csv name, local derivatives tree, slug for metrics.json).
# Output tree mirrors the trainer's canonical layout so aggregator globs pick up
# the shims without code change. Slug is what goes into config.model_id.
PROJECTS = {
    "brainmvp_frozen_cached": dict(
        csv="wandb_brainmvp_frozen_cached.csv",
        tree=r"D:/ADNI_BIDS_project/derivatives/brainmvp_debug/aug_none/BrainMVP_uniformer_frozen_cached",
        model_id_fallback="BrainMVP_uniformer",
    ),
    "vit_mae_frozen_cached": dict(
        csv="wandb_vit_mae_frozen_cached.csv",
        tree=r"D:/ADNI_BIDS_project/derivatives/vit_outputs_debug/aug_none/ViT_B_mae75_frozen_cached",
        model_id_fallback="ViT_B_mae75",
    ),
    "braindino_frozen_cached": dict(
        csv="wandb_braindino_frozen_cached.csv",
        # Local dir renamed from aug_none/ to aug_none_hp_tuned/ to disambiguate
        # this cached-head HP sweep from the on-the-fly BrainDINO/frozen/aug_none
        # cells in braindino_outputs/frozen/. Aggregator glob in
        # 05_aggregate_mri_results.py matches the new name.
        tree=r"D:/ADNI_BIDS_project/derivatives/braindino_outputs/aug_none_hp_tuned/BrainDINO_vitb16_frozen_cached",
        model_id_fallback="BrainDINO_vitb16",
    ),
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project", default="brainmvp_frozen_cached",
                   choices=list(PROJECTS.keys()),
                   help="Which cached-head W&B project to materialise. "
                        "Default: brainmvp_frozen_cached.")
    p.add_argument("--task", default=None,
                   help="Restrict to a single task (e.g. T1c_binary). "
                        "Default: all tasks present in the CSV.")
    return p.parse_args()


def _hp_leaf(row: pd.Series) -> str:
    """Mirror the submit script's leaf naming:
        lr$(printf '%.0e' $LR)_d${DROPOUT}_ls${LS}
    Bash's %.0e for 0.001 -> '1e-03'. Python emits '1e-03' too via %.0e.
    """
    lr = row["config/lr"]
    d  = row["config/drop_rate"]
    ls = row["config/label_smoothing"]
    return f"lr{lr:.0e}_d{d}_ls{ls}"


def _row_to_metrics(row: pd.Series, model_id_fallback: str) -> dict:
    """Build a metrics.json dict matching the cached-head trainer schema.
    Only val_* fields are present in W&B summary -- test_metrics + test_diag
    are left empty so the user knows to rsync the canonical file later."""
    def s(key, default=None):
        v = row.get(f"summary/{key}")
        if pd.isna(v):
            return default
        return v

    val_TP = s("val_TP"); val_FP = s("val_FP")
    val_TN = s("val_TN"); val_FN = s("val_FN")
    model_id = (row.get("config/model_name") if "config/model_name" in row.index
                else None) or model_id_fallback

    return {
        "config": {
            "model_id":              model_id,
            "model_kind":            "frozen_cached",
            "strategy":              "frozen_cached",
            "augment":               "none",
            "task":                  row["config/task"],
            "seed":                  int(row["config/seed"]),
            "num_labels":            None,
            "epochs":                int(row["config/epochs"]),
            "best_epoch":            int(s("epoch")) if s("epoch") else None,
            "best_val_balanced_acc": float(s("best_val_bacc"))
                                     if s("best_val_bacc") else None,
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
        # val_* metrics from W&B summary (the cached-head trainer's full set).
        "val_metrics": {
            "balanced_acc": float(s("val_bacc")) if s("val_bacc") else None,
            "loss":         float(s("val_loss")) if s("val_loss") else None,
            "auc_roc":      float(s("val_auc"))  if s("val_auc") else None,
            "f1":           float(s("val_f1"))   if s("val_f1") else None,
            "precision":    float(s("val_prec")) if s("val_prec") else None,
            "recall":       float(s("val_recall")) if s("val_recall") else None,
            "sensitivity":  float(s("val_sens")) if s("val_sens") else None,
            "specificity":  float(s("val_spec")) if s("val_spec") else None,
            "TP":           int(val_TP) if val_TP is not None else None,
            "FP":           int(val_FP) if val_FP is not None else None,
            "TN":           int(val_TN) if val_TN is not None else None,
            "FN":           int(val_FN) if val_FN is not None else None,
        },
        # Empty -- the cached-head trainer DOES compute test metrics but
        # does NOT log them to W&B; they only live in the canonical
        # metrics.json on CSD3. Rsync that down and the test table populates.
        "test_metrics": {},
        "test_metrics_subject": {},
        "test_diagnostics": {},
    }


def main():
    args = parse_args()
    cfg = PROJECTS[args.project]
    wandb_csv = REPO_ROOT / "mri_pipeline" / "cached_head_sweep" / cfg["csv"]
    deriv_root = Path(cfg["tree"])
    model_id_fb = cfg["model_id_fallback"]

    if not wandb_csv.exists():
        raise SystemExit(f"W&B CSV not found: {wandb_csv}. Run "
                         f"06_fetch_wandb_tables.py --projects {args.project} first.")

    df = pd.read_csv(wandb_csv)
    df = df[df["state"] == "finished"]
    if args.task:
        df = df[df["config/task"] == args.task]
        if df.empty:
            raise SystemExit(f"No finished runs for task={args.task}")
    df = df.sort_values(["config/task", "config/seed", "config/lr",
                         "config/drop_rate", "config/label_smoothing"])
    print(f"Project  : {args.project}")
    print(f"Out tree : {deriv_root}")
    print(f"Found {len(df)} finished runs"
          + (f" for task={args.task}" if args.task else " across all tasks"))

    n_written = n_skipped = 0
    for _, row in df.iterrows():
        task = row["config/task"]
        seed = int(row["config/seed"])
        leaf = _hp_leaf(row)
        out_dir = deriv_root / task / f"seed_{seed}" / leaf
        out_dir.mkdir(parents=True, exist_ok=True)
        out_file = out_dir / "metrics.json"

        if out_file.exists():
            with open(out_file) as f:
                existing = json.load(f)
            if existing.get("config", {}).get("_provenance") != "shim_from_wandb_summary":
                n_skipped += 1
                continue

        payload = _row_to_metrics(row, model_id_fb)
        with open(out_file, "w") as f:
            json.dump(payload, f, indent=2)
        n_written += 1

    print(f"Wrote {n_written} shim(s); skipped {n_skipped} canonical metrics.json")
    print()
    print("NOTE: shims contain val_* metrics only (test cells will read '-' in")
    print("the TEST cross-model table). Once you rsync the canonical")
    print("metrics.json from CSD3, the test cells populate automatically.")


if __name__ == "__main__":
    main()
