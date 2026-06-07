"""
04d_rerun_cached_grid_june1.py — re-run the FULL cached-head HP grid on the
CURRENT (June-1) embeddings for the diagnostic tasks
=============================================================================
Decision (user, 2026-06-07): the June-1 frozen-encoder pooled embeddings are
canonical (the May-25 cached numbers predate them). To make the cached `none†`
rows consistent, re-run the ENTIRE 18-HP grid (lr × drop × ls) × 3 seeds for
BrainMVP + ViT-MAE on T1/T1b/T1c/T1d/T2 against the current embeddings — so the
HP-winner selection is internally consistent and every winner cell has a
checkpoint + val_auc/f1.

Deletes the stale flat task-subtree first (the May-25 originals remain in the
double-nested `*/*/aug_none/...` tree as a backup). Standard grid only
(head=mlp, loss=ce_weighted, no standardize), matching the original cached rows.

Run (mri env, CPU): KMP_DUPLICATE_LIB_OK=TRUE python .../04d_rerun_cached_grid_june1.py
  ... --dry-run
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

MRI = Path(__file__).resolve().parents[1]
TRAINER = MRI / "cached_head_sweep" / "04_head_finetune_from_embeddings.py"
ROOT = Path(r"D:/ADNI_BIDS_project/derivatives")
MASTER = ROOT / "mri_clinical_matched" / "VISCODE_2_aligned_extended_post_exclusion" / \
    "master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR = ROOT / "clinical" / "no_cdr_stratified_post_exclusion" / "tabular" / "baseline"

TASKS = ["T1_binary", "T1b_binary", "T1c_binary", "T1d_binary", "T2_multiclass"]
SEEDS = [0, 1, 2]
LRS = ["1e-3", "1e-4", "1e-5"]
DROPS = ["0.1", "0.2", "0.3"]
LSS = ["0.0", "0.1"]

# Outputs go to a CLEARLY-SEPARATE tree so the June-1-embedding runs never mix
# with / overwrite the original May-25 cached runs (which stay in their trees).
JUN01 = ROOT / "cached_embeddings_using_JUN_01"
MODELS = {
    "BrainMVP":   dict(emb=ROOT / "cached_embeddings" / "brainmvp" / "brainmvp_pooled.pt",
                       dim=512, base=JUN01 / "BrainMVP_uniformer_frozen_cached"),
    "ViT-MAE75":  dict(emb=ROOT / "cached_embeddings" / "vit_mae" / "vit_mae_pooled.pt",
                       dim=768, base=JUN01 / "ViT_B_mae75_frozen_cached"),
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--models", nargs="*", default=list(MODELS))
    a = p.parse_args()

    total = len(a.models) * len(TASKS) * len(LRS) * len(DROPS) * len(LSS) * len(SEEDS)
    print(f"[plan] {total} cells "
          f"({len(a.models)} models × {len(TASKS)} tasks × 18 HP × {len(SEEDS)} seeds)")
    ok = fail = 0
    for mname in a.models:
        cfg = MODELS[mname]
        for task in TASKS:
            cell = cfg["base"] / task
            if not a.dry_run and cell.exists():
                shutil.rmtree(cell)            # drop stale (May-25 / mixed) flat cells
            for lr in LRS:
                for drop in DROPS:
                    for ls in LSS:
                        for s in SEEDS:
                            cmd = [sys.executable, str(TRAINER),
                                   "--out_dir", str(cfg["base"]),
                                   "--model_name", mname, "--embed_dim", str(cfg["dim"]),
                                   "--embeddings_pt", str(cfg["emb"]),
                                   "--task", task, "--seed", str(s),
                                   "--lr", lr, "--drop_rate", drop, "--label_smoothing", ls,
                                   "--weight_decay", "1e-5",
                                   "--matched_labels_csv", str(MASTER),
                                   "--data_dir", str(DATA_DIR), "--long", "all"]
                            if a.dry_run:
                                ok += 1
                                if ok <= 2: print("  [dry] " + " ".join(cmd))
                                continue
                            rc = subprocess.run(
                                cmd, env=dict(os.environ, KMP_DUPLICATE_LIB_OK="TRUE")).returncode
                            ok += (rc == 0); fail += (rc != 0)
            print(f"  {mname}/{task}: done ({ok} ok, {fail} fail so far)")
    print(f"[done] ok={ok} fail={fail}")


if __name__ == "__main__":
    main()
