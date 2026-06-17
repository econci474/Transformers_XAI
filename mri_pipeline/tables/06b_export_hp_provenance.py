"""
06b_export_hp_provenance.py — per-run HP + provenance table for the thesis
==========================================================================
For every training run that the cross-model tables actually display, dump one
CSV row recording **which exact run/sweep it came from** (provenance) and the
**hyper-parameters** it was trained with. This is the methods-table source for
the thesis: it lets us state, per (model, task, strategy, augment, seed), the
lr / weight-decay / dropout / label-smoothing / LLRD / augment and the on-disk
path the displayed numbers were read from.

Discovery reuses `05.MODEL_TREES` (the single source of truth for which run
feeds each table cell) with the same flat → double-nested fallback as
`06._collect_runs`, so the CSV is one-to-one with the rendered tables — e.g.
ViT-MAE frozen/random resolves to `vit_outputs_debug/.../frozen` while full_ft
+ frozen/{none,plus_original} resolve to `vit_outputs_hi_lr/...`.

Output: mri_pipeline/outputs/cross_model_hp_provenance.csv
Run (no torch needed — pure json/pandas):
  KMP_DUPLICATE_LIB_OK=TRUE python mri_pipeline/06b_export_hp_provenance.py
"""
from __future__ import annotations

import argparse
import glob
import importlib.util
import json
import os
from pathlib import Path

import pandas as pd

MRI = Path(__file__).resolve().parent
DEFAULT_ROOT = r"D:/ADNI_BIDS_project/derivatives"

# Pull MODEL_TREES from 05 so this never drifts from what the tables render.
_spec = importlib.util.spec_from_file_location(
    "_05_aggregate", str(MRI / "05_aggregate_mri_results.py"))
_05 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_05)
MODEL_TREES = _05.MODEL_TREES

# Columns emitted, in order. HP fields are model-dependent; missing ones → "".
HP_FIELDS = [
    "lr", "weight_decay",
    "drop_rate", "drop_path_rate", "attn_dropout",   # dropout variants per arch
    "label_smoothing", "llrd_gamma", "augment", "aug_copies",
    "epochs", "best_epoch", "warmup_epochs", "patience", "batch_size",
    "grad_accum_steps", "effective_batch_size",
]
VAL_FIELDS = ["best_val_balanced_acc", "best_val_auc", "best_val_f1"]


def _test_metric(tm: dict, *keys):
    return next((tm[k] for k in keys if k in tm), None)


def collect(root: str) -> pd.DataFrame:
    rows = []
    for model, rel in MODEL_TREES:
        tree = rel.split("/")[0]
        files = sorted(glob.glob(os.path.join(root, rel)))
        nesting = "flat"
        if not files:                                  # rsync double-nest fallback
            files = sorted(glob.glob(os.path.join(root, tree, rel)))
            nesting = "double-nested" if files else "flat"
        if not files:
            print(f"  [skip] {model:18s} — no metrics.json under {tree}/")
            continue
        for f in files:
            try:
                d = json.load(open(f))
            except Exception as exc:
                print(f"  [WARN] unreadable: {f} ({exc})")
                continue
            cfg = d.get("config", {})
            tm = d.get("test_metrics", {})
            row = {
                "model": model,
                "task": cfg.get("task"),
                "strategy": cfg.get("strategy"),
                "augment": cfg.get("augment"),
                "seed": cfg.get("seed"),
            }
            for k in HP_FIELDS:
                row[k] = cfg.get(k)
            for k in VAL_FIELDS:
                row[k] = cfg.get(k)
            row["test_balanced_acc"] = tm.get("balanced_acc")
            row["test_auc"] = _test_metric(tm, "auc_roc_ovr", "auc_roc")
            row["test_macro_f1"] = _test_metric(tm, "macro_f1", "f1")
            row["pretrained_ckpt"] = cfg.get("pretrained_ckpt") or cfg.get("embeddings_pt")
            row["n_train"], row["n_val"], row["n_test"] = (
                cfg.get("n_train"), cfg.get("n_val"), cfg.get("n_test"))
            # Provenance: sweep tree + path relative to root + on-disk nesting.
            row["source_tree"] = tree
            row["source_nesting"] = nesting
            row["source_path"] = os.path.relpath(f, root).replace(os.sep, "/")
            row["timestamp"] = cfg.get("timestamp")
            rows.append(row)
        print(f"  [ok]   {model:18s} — {len(files)} runs ({tree}/, {nesting})")
    return pd.DataFrame(rows)


def filter_cached_to_hp_winners(df: pd.DataFrame) -> pd.DataFrame:
    """Mirror of `06._filter_cached_to_hp_winners` so the table-used CSV shows
    exactly the cached HP that the cross-model tables display: per (model, task)
    keep only the hp_leaf (`lr<>_d<>_ls<>` dir) with the highest mean
    best_val_balanced_acc across seeds. Non-cached models pass through unchanged
    (each cell is already single-HP)."""
    cached_mask = df["model"].str.endswith("-cached")
    if not cached_mask.any():
        return df
    cached, other = df[cached_mask].copy(), df[~cached_mask].copy()
    cached["hp_leaf"] = cached["source_path"].apply(
        lambda p: os.path.basename(os.path.dirname(p)))
    cached["_val"] = pd.to_numeric(cached["best_val_balanced_acc"], errors="coerce")
    grp = (cached.groupby(["model", "task", "hp_leaf"], dropna=False)["_val"]
                 .mean().reset_index())
    winners = grp.loc[grp.groupby(["model", "task"], dropna=False)["_val"].idxmax()]
    keep = set(zip(winners["model"], winners["task"], winners["hp_leaf"]))
    cached = cached[cached.apply(
        lambda r: (r["model"], r["task"], r["hp_leaf"]) in keep, axis=1)]
    cached = cached.drop(columns=["hp_leaf", "_val"])
    return pd.concat([other, cached], ignore_index=True)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", default=DEFAULT_ROOT)
    p.add_argument("--out", default=str(MRI / "outputs" / "cross_model_hp_provenance.csv"))
    a = p.parse_args()

    df = collect(a.root)
    if df.empty:
        print("[ERROR] no runs collected."); return
    order = ["model", "task", "strategy", "augment", "seed"]
    df = df.sort_values(order, na_position="last").reset_index(drop=True)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, index=False)
    print(f"\n[done] ALL runs: {len(df)} → {a.out}")

    # Table-used subset: cached filtered to the displayed HP winner per cell.
    used = filter_cached_to_hp_winners(df).sort_values(
        order, na_position="last").reset_index(drop=True)
    used_path = a.out.replace(".csv", "_table_used.csv")
    used.to_csv(used_path, index=False)
    print(f"[done] TABLE-USED runs (cached→winner HP): {len(used)} → {used_path}")
    print(f"       {used['model'].nunique()} models, {used['task'].nunique()} tasks, "
          f"{used.groupby(['model','task','strategy','augment']).ngroups} displayed cells")


if __name__ == "__main__":
    main()
