"""
04c_retrain_cached_nockpt.py — fill cached `none†` val AUC/F1 for winner cells
that have NO checkpoint
=============================================================================
Some cached-head HP-winner cells were run `--metrics_only` (no `best_model.pt`),
so `04b` (--val_test recompute) can't load them. The cached-head trainer is
deterministic (CPU + seeded), so we simply RE-TRAIN the winner HP into the same
cell: delete the stale metrics.json (+ last_checkpoint.pt) so the trainer doesn't
skip, run WITHOUT --metrics_only (saves best_model.pt + train_log with val_auc/f1
+ metrics.json best_val_auc/f1), then verify the reproduced best_val_balanced_acc
matches the old one (provenance guard) — if it diverges > tol we flag it.

Winner = max mean best_val_balanced_acc across seeds per (model, task), mirroring
06._filter_cached_to_hp_winners. Only cells whose winner leaf lacks a checkpoint
are re-trained. Runs locally (mri env, CPU).

Usage (mri env):
  KMP_DUPLICATE_LIB_OK=TRUE python mri_pipeline/cached_head_sweep/04c_retrain_cached_nockpt.py
  ... --dry-run     # print what would be retrained
  ... --models BrainMVP-cached ViT-MAE-cached   # restrict
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
PROV = MRI / "outputs" / "cross_model_hp_provenance.csv"
MASTER = ROOT / "mri_clinical_matched" / "VISCODE_2_aligned_extended_post_exclusion" / \
    "master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR = ROOT / "clinical" / "no_cdr_stratified_post_exclusion" / "tabular" / "baseline"
TOL = 0.02


def local_embeddings(hpc_path: str) -> Path:
    p = Path(hpc_path)
    return ROOT / "cached_embeddings" / p.parent.name / p.name


def winners_without_ckpt(prov_csv):
    """Return winner cells (cached) that lack best_model.pt — mirrors the render's
    HP-winner-per-(model,task) selection."""
    df = pd.read_csv(prov_csv)
    c = df[df["model"].str.endswith("-cached")].copy()
    c["hp_leaf"] = c["source_path"].apply(lambda p: os.path.basename(os.path.dirname(p)))
    c["_v"] = pd.to_numeric(c["best_val_balanced_acc"], errors="coerce")
    grp = c.groupby(["model", "task", "hp_leaf"])["_v"].mean().reset_index()
    win = grp.loc[grp.groupby(["model", "task"])["_v"].idxmax()]
    keep = set(zip(win.model, win.task, win.hp_leaf))
    w = c[c.apply(lambda r: (r.model, r.task, r.hp_leaf) in keep, axis=1)]
    rows = []
    for _, r in w.iterrows():
        cell = ROOT / os.path.dirname(r["source_path"])
        if not (cell / "best_model.pt").exists():
            rows.append((r["model"], r["task"], int(r["seed"]), cell))
    return rows


def build_cmd(cfg, base_out):
    emb = local_embeddings(cfg["embeddings_pt"])
    cmd = [sys.executable, str(TRAINER),
           "--out_dir", str(base_out), "--no_resume",
           "--model_name", str(cfg["model_name"]),
           "--embed_dim", str(cfg["embed_dim"]),
           "--embeddings_pt", str(emb),
           "--task", str(cfg["task"]), "--seed", str(cfg["seed"]),
           "--drop_rate", str(cfg.get("drop_rate", 0.2)),
           "--lr", str(cfg.get("lr", 1e-4)),
           "--weight_decay", str(cfg.get("weight_decay", 1e-5)),
           "--label_smoothing", str(cfg.get("label_smoothing", 0.0)),
           "--matched_labels_csv", str(MASTER), "--data_dir", str(DATA_DIR)]
    if cfg.get("long_mode") == "all":
        cmd += ["--long", "all"]
    elif cfg.get("long_mode") == "cutoff" and cfg.get("max_months"):
        cmd += ["--long", str(int(cfg["max_months"] // 12))]
    if cfg.get("head") and cfg["head"] != "mlp":
        cmd += ["--head", str(cfg["head"])]
    if cfg.get("standardize"):
        cmd += ["--standardize"]
    if cfg.get("loss") and cfg["loss"] != "ce_weighted":
        cmd += ["--loss", str(cfg["loss"])]
        if cfg["loss"] == "focal" and cfg.get("focal_gamma") is not None:
            cmd += ["--focal_gamma", str(cfg["focal_gamma"])]
    return cmd


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--prov", default=str(PROV))
    p.add_argument("--models", nargs="*", default=None,
                   help="Restrict to these cached model labels (e.g. BrainMVP-cached).")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()

    cells = winners_without_ckpt(a.prov)
    if a.models:
        cells = [c for c in cells if c[0] in a.models]
    print(f"[plan] {len(cells)} no-checkpoint cached winner cells to re-train")
    ok = warn = fail = 0
    for model, task, seed, cell in sorted(cells, key=lambda x: (x[0], x[1], x[2])):
        mpath = cell / "metrics.json"
        cfg = json.load(open(mpath)).get("config", {}) if mpath.exists() else {}
        old_bacc = cfg.get("best_val_balanced_acc")
        # base_out = the tree root the trainer appends <task>/seed/<hp> to.
        base_out = cell.parent.parent.parent
        tag = f"{model}/{task}/seed{seed}/{cell.name}"
        cmd = build_cmd(cfg, base_out)
        if a.dry_run:
            print(f"  [dry] {tag}\n        " + " ".join(cmd)); ok += 1; continue
        # delete stale outputs so the trainer re-runs (it skips if metrics.json exists)
        for f in ("metrics.json", "last_checkpoint.pt", "best_model.pt"):
            try: (cell / f).unlink()
            except FileNotFoundError: pass
        print(f"\n=== {tag}  (old bACC {old_bacc}) ===")
        rc = subprocess.run(cmd, env=dict(os.environ, KMP_DUPLICATE_LIB_OK="TRUE")).returncode
        if rc != 0:
            print(f"  [FAIL rc={rc}] {tag}"); fail += 1; continue
        new = json.load(open(mpath)).get("config", {}).get("best_val_balanced_acc") \
            if mpath.exists() else None
        if old_bacc is not None and new is not None and abs(float(old_bacc) - float(new)) > TOL:
            print(f"  [WARN] reproduced bACC {new} != old {old_bacc} (>|{TOL}|)"); warn += 1
        else:
            print(f"  [ok] reproduced bACC {new} (old {old_bacc})"); ok += 1
    print(f"\n[done] ok={ok} warn={warn} fail={fail}")


if __name__ == "__main__":
    main()
