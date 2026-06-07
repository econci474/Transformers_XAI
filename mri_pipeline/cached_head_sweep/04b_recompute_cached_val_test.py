"""
04b_recompute_cached_val_test.py — local CPU val-AUC/F1 recompute for the
cached-head TABLE WINNERS
=========================================================================
The cached `none†` rows in the cross-model tables were trained before the
cached-head trainer logged val_auc/val_f1, so their metrics.json has no
val AUC/F1. This driver re-evaluates ONLY the displayed winner HP cell per
(model, task, seed) — read from
`mri_pipeline/outputs/cross_model_hp_provenance_table_used.csv` — by calling
`04_head_finetune_from_embeddings.py --val_test`, which loads that cell's own
`best_model.pt`, recomputes val + test, and patches a `val_metrics` block into
the same metrics.json (no retrain, no weight overwrite; >0.02 test guard).

Deterministic (CPU + seeded), so the recomputed test reproduces the published
test (provenance check) and the new val_auc/val_f1 belong to the displayed run.

Run (mri env; CPU; no monai needed):
  KMP_DUPLICATE_LIB_OK=TRUE python mri_pipeline/cached_head_sweep/04b_recompute_cached_val_test.py
  ... --dry-run        # print commands only
  ... --only-missing   # skip cells whose metrics.json already has a val_metrics block
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

MRI = Path(__file__).resolve().parents[1]
TRAINER = MRI / "cached_head_sweep" / "04_head_finetune_from_embeddings.py"
ROOT = Path(r"D:/ADNI_BIDS_project/derivatives")
TABLE_USED = MRI / "outputs" / "cross_model_hp_provenance_table_used.csv"

# Local paths (the HPC paths stored in each config are remapped to these).
MASTER = ROOT / "mri_clinical_matched" / "VISCODE_2_aligned_extended_post_exclusion" / \
    "master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR = ROOT / "clinical" / "no_cdr_stratified_post_exclusion" / "tabular" / "baseline"


def local_embeddings(hpc_path: str) -> Path:
    """Remap `.../cached_embeddings/<sub>/<file>.pt` (HPC) → local D: tree."""
    p = Path(hpc_path)
    return ROOT / "cached_embeddings" / p.parent.name / p.name


def build_cmd(cell_dir: Path, cfg: dict) -> list[str]:
    emb = local_embeddings(cfg["embeddings_pt"])
    cmd = [sys.executable, str(TRAINER), "--val_test",
           "--out_dir", str(cell_dir),
           "--model_name", str(cfg["model_name"]),
           "--embed_dim", str(cfg["embed_dim"]),
           "--embeddings_pt", str(emb),
           "--task", str(cfg["task"]),
           "--seed", str(cfg["seed"]),
           "--drop_rate", str(cfg.get("drop_rate", 0.2)),
           "--lr", str(cfg.get("lr", 1e-4)),
           "--weight_decay", str(cfg.get("weight_decay", 1e-5)),
           "--label_smoothing", str(cfg.get("label_smoothing", 0.0)),
           "--matched_labels_csv", str(MASTER),
           "--data_dir", str(DATA_DIR)]
    # Cohort axis: every cached MRI cell is long_mode=all.
    if cfg.get("long_mode") == "all":
        cmd += ["--long", "all"]
    elif cfg.get("long_mode") == "cutoff" and cfg.get("max_months"):
        cmd += ["--long", str(int(cfg["max_months"] // 12))]
    # Architecture-affecting knobs must match the saved state_dict.
    if cfg.get("head") and cfg["head"] != "mlp":
        cmd += ["--head", str(cfg["head"])]
    if cfg.get("standardize"):
        cmd += ["--standardize"]
    if cfg.get("loss") and cfg["loss"] != "ce_weighted":
        cmd += ["--loss", str(cfg["loss"])]
        if cfg["loss"] == "focal" and cfg.get("focal_gamma") is not None:
            cmd += ["--focal_gamma", str(cfg["focal_gamma"])]
    if cfg.get("select_by") and cfg["select_by"] != "bacc":
        cmd += ["--select_by", str(cfg["select_by"])]
    return cmd


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--table_used", default=str(TABLE_USED))
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--only-missing", action="store_true",
                   help="Skip cells whose metrics.json already has a val_metrics block.")
    a = p.parse_args()

    df = pd.read_csv(a.table_used)
    cached = df[df["model"].str.endswith("-cached")].copy()
    print(f"[plan] {len(cached)} cached winner cells "
          f"({cached['model'].nunique()} models, {cached['task'].nunique()} tasks)")

    ok = skipped = failed = 0
    for _, r in cached.sort_values(["model", "task", "seed"]).iterrows():
        cell_dir = (ROOT / r["source_path"]).parent
        mpath = cell_dir / "metrics.json"
        if not mpath.exists():
            print(f"  [MISS] no metrics.json: {cell_dir}"); failed += 1; continue
        cfg = json.load(open(mpath)).get("config", {})
        if a.only_missing and "val_metrics" in json.load(open(mpath)):
            skipped += 1; continue
        cmd = build_cmd(cell_dir, cfg)
        tag = f"{r['model']}/{r['task']}/seed{r['seed']}"
        if a.dry_run:
            print(f"  [{tag}] " + " ".join(cmd)); ok += 1; continue
        print(f"\n=== {tag} ===")
        rc = subprocess.run(cmd, env=dict(os.environ, KMP_DUPLICATE_LIB_OK="TRUE")).returncode
        if rc == 0:
            ok += 1
        else:
            print(f"  [FAIL rc={rc}] {tag}"); failed += 1

    print(f"\n[done] ok={ok} skipped={skipped} failed={failed}")


if __name__ == "__main__":
    main()
