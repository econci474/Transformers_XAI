"""
05b_select_best_hp_per_task.py
==============================
Post-sweep winner-selection for the cached-embedding head HP sweeps
(BrainDINO / BrainMVP / ViT-MAE75). Picks the winning HP combo per
task by mean val_bacc across the 3 seeds, then symlinks (or copies on
Windows) the 3 corresponding `best_model.pt` files into a clean
`_winners/` tree for downstream use (XAI, attention maps, thesis
tables).

WHY this is needed
------------------
Each cell of the head sweep saves its own `best_model.pt` (best
WITHIN that cell's training, by val_bacc). With 270 cells per model
(18 HP combos x 3 seeds x 5 tasks) you get 270 best_model.pt files.
There's no automatic "best per task" -- you'd have to dig manually.

This script does that dig, principled (HPs tuned on val, reported on
test -- no leakage):

  1. Group all 54 cells per task (18 HP x 3 seeds) by HP combo.
  2. Compute mean `best_val_balanced_acc` across the 3 seeds.
  3. Pick the HP combo with the highest mean val_bacc -> "winning HP".
  4. Report the TEST metrics (mean +/- std across 3 seeds) for that
     winning HP -- the headline number for your thesis.
  5. Symlink the 3 `best_model.pt` files for that winning HP to a
     stable path under `_winners/<task>/seed_<n>/best_model.pt`.

Outputs (under `<model_root>/_winners/`)
----------------------------------------
  winners.csv          single row per task: winning HP + test metrics +
                       paths to the 3 winning cells
  top_5_per_task.csv   top-5 HP combos per task by val_bacc (use this
                       for the "median +/- IQR of top-K" treatment if
                       the writeup wants honest variance reporting --
                       see C.8.6 of the plan)
  <task>/seed_<n>/best_model.pt   symlinks to the 3 winning checkpoints

Usage
-----
  python mri_pipeline/cached_head_sweep/05b_select_best_hp_per_task.py \\
      --model_root <model_root>/aug_none/<MODEL>_frozen_cached \\
      --model_name <MODEL>

  # Three model invocations (after each sweep finishes):
  python ... --model_root .../braindino_outputs/aug_none/BrainDINO_vitb16_frozen_cached      --model_name BrainDINO
  python ... --model_root .../brainmvp_debug/aug_none/BrainMVP_uniformer_frozen_cached       --model_name BrainMVP
  python ... --model_root .../vit_outputs_debug/aug_none/ViT_B_mae75_frozen_cached           --model_name ViT-MAE75
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ── Args ──────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--model_root", type=str, required=True,
                   help="Root dir containing <task>/seed_<n>/<hp>/metrics.json "
                        "(e.g. .../BrainDINO_vitb16_frozen_cached/).")
    p.add_argument("--model_name", type=str, required=True,
                   help="Label written into the winners CSV (BrainDINO / "
                        "BrainMVP / ViT-MAE75).")
    p.add_argument("--selection_metric", type=str, default="best_val_balanced_acc",
                   help="Field on metrics.json['config'] to pick the winner "
                        "by (default: best_val_balanced_acc -- val, not test, "
                        "to avoid leakage).")
    p.add_argument("--out_dir", type=str, default=None,
                   help="Where to write winners.csv + symlinks "
                        "(default: <model_root>/_winners/).")
    p.add_argument("--top_k", type=int, default=5,
                   help="How many top-K HP combos per task to record in "
                        "top_k_per_task.csv (default: 5).")
    p.add_argument("--copy", action="store_true",
                   help="Copy best_model.pt files instead of symlinking. Use "
                        "on Windows or any FS where symlinks aren't supported.")
    p.add_argument("--dry-run", action="store_true",
                   help="Print the selection but don't write files / symlinks.")
    return p.parse_args()


# ── Discover + parse metrics.json files ───────────────────────────────────────
_HP_RE = re.compile(
    r"lr(?P<lr>[\w+.-]+)_d(?P<drop_rate>[\d.]+)_ls(?P<label_smoothing>[\d.]+)$"
)


def parse_hp_from_dirname(d: str) -> dict | None:
    """The trainer writes cells into `lr{lr}_d{drop_rate}_ls{label_smoothing}/`.
    Parse the HP back from the dir name as a fallback for runs that for some
    reason lack those fields in metrics.json. Returns None if parse fails."""
    m = _HP_RE.match(d)
    if not m:
        return None
    return {
        "lr":              float(m.group("lr")),
        "drop_rate":       float(m.group("drop_rate")),
        "label_smoothing": float(m.group("label_smoothing")),
    }


def load_runs(model_root: Path) -> pd.DataFrame:
    """Scan all metrics.json under `<model_root>/<task>/seed_<n>/<hp>/`.
    Returns a long DataFrame, one row per cell.
    """
    pattern = str(model_root / "*" / "seed_*" / "*" / "metrics.json")
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"  [ERROR] No metrics.json found under {pattern}", file=sys.stderr)
        sys.exit(1)
    print(f"  Found {len(files)} cells under {model_root}")

    rows = []
    for f in files:
        try:
            with open(f) as fh:
                d = json.load(fh)
        except Exception as exc:
            print(f"  [WARN] unreadable: {f} ({exc})", file=sys.stderr)
            continue
        cfg = d.get("config", {})
        tm  = d.get("test_metrics", {}) or {}

        # Reconstruct HP from dir if missing in config (defensive)
        cell_dir = os.path.dirname(f)
        hp_from_dir = parse_hp_from_dirname(os.path.basename(cell_dir)) or {}

        rows.append({
            "cell_dir":       cell_dir,
            "best_model_pt":  os.path.join(cell_dir, "best_model.pt"),
            "task":           cfg.get("task", "?"),
            "seed":           int(cfg.get("seed", -1)),
            "lr":             cfg.get("lr", hp_from_dir.get("lr")),
            "drop_rate":      cfg.get("drop_rate", hp_from_dir.get("drop_rate")),
            "label_smoothing": cfg.get("label_smoothing",
                                        hp_from_dir.get("label_smoothing")),
            "val_bacc":       cfg.get("best_val_balanced_acc"),
            "test_bacc":      tm.get("balanced_acc"),
            "test_auc":       tm.get("auc_roc", tm.get("auc_roc_ovr")),
            "test_f1":        tm.get("f1", tm.get("macro_f1")),
        })
    df = pd.DataFrame(rows)
    n_missing_pt = int((~df["best_model_pt"].apply(os.path.isfile)).sum())
    if n_missing_pt:
        print(f"  [WARN] {n_missing_pt}/{len(df)} cells have no best_model.pt -- "
              f"likely NaN-aborted before the first epoch saved one.")
    return df


# ── Selection ─────────────────────────────────────────────────────────────────
HP_KEYS = ["lr", "drop_rate", "label_smoothing"]


def select_winners(df: pd.DataFrame, selection_metric_col: str = "val_bacc",
                   top_k: int = 5) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    For each task: group by HP combo, aggregate val_bacc across seeds,
    pick the HP combo with the highest mean val_bacc.

    Returns (winners_df, top_k_df).
    winners_df: one row per task -- the winning HP + its test metrics.
    top_k_df:   top-K HP combos per task by mean val_bacc (long form).
    """
    # Aggregate per (task, hp_combo)
    agg = (
        df.groupby(["task"] + HP_KEYS, dropna=False)
          .agg(
              val_bacc_mean    = (selection_metric_col, "mean"),
              val_bacc_std     = (selection_metric_col, "std"),
              val_bacc_n       = (selection_metric_col, "count"),
              test_bacc_mean   = ("test_bacc", "mean"),
              test_bacc_std    = ("test_bacc", "std"),
              test_auc_mean    = ("test_auc", "mean"),
              test_auc_std     = ("test_auc", "std"),
              test_f1_mean     = ("test_f1", "mean"),
              test_f1_std      = ("test_f1", "std"),
              cell_dirs        = ("cell_dir",
                                  lambda s: ",".join(sorted(s.astype(str)))),
          )
          .reset_index()
    )

    # Per-task top-K (sorted by val_bacc desc, ties by test_bacc desc)
    top_k_chunks = []
    for task, sub in agg.groupby("task"):
        ordered = sub.sort_values(
            ["val_bacc_mean", "test_bacc_mean"],
            ascending=[False, False],
            na_position="last",
        )
        ordered = ordered.assign(rank=range(1, len(ordered) + 1))
        top_k_chunks.append(ordered.head(top_k))
    top_k_df = pd.concat(top_k_chunks, ignore_index=True)

    # Winners = rank-1 per task
    winners_df = top_k_df[top_k_df["rank"] == 1].copy()
    winners_df = winners_df.drop(columns="rank").reset_index(drop=True)
    return winners_df, top_k_df


# ── Symlink (or copy) winning checkpoints ─────────────────────────────────────
def emplace_winning_pts(winners_df: pd.DataFrame, runs_df: pd.DataFrame,
                        out_dir: Path, copy: bool = False,
                        dry_run: bool = False) -> list[dict]:
    """For each winning HP per task, find the 3 cell dirs (one per seed)
    and symlink (or copy) their best_model.pt to
    `<out_dir>/<task>/seed_<n>/best_model.pt`.

    Returns a list of dicts describing each created link/copy."""
    actions = []
    for _, row in winners_df.iterrows():
        task = row["task"]
        hp_mask = (
            (runs_df["task"] == task) &
            (np.isclose(runs_df["lr"], row["lr"])) &
            (np.isclose(runs_df["drop_rate"], row["drop_rate"])) &
            (np.isclose(runs_df["label_smoothing"], row["label_smoothing"]))
        )
        winning_cells = runs_df[hp_mask].sort_values("seed")
        for _, cell in winning_cells.iterrows():
            src = Path(cell["best_model_pt"])
            if not src.is_file():
                continue
            dst_dir = out_dir / task / f"seed_{int(cell['seed'])}"
            dst = dst_dir / "best_model.pt"
            actions.append({
                "task": task, "seed": int(cell["seed"]),
                "src": str(src), "dst": str(dst),
                "method": "copy" if copy else "symlink",
            })
            if dry_run:
                continue
            dst_dir.mkdir(parents=True, exist_ok=True)
            if dst.exists() or dst.is_symlink():
                dst.unlink()
            if copy:
                shutil.copy2(src, dst)
            else:
                # Windows requires admin or Developer Mode for symlinks; if
                # symlink fails, fall back to a hard copy to keep the workflow
                # robust on local mirrors.
                try:
                    os.symlink(src, dst)
                except OSError:
                    shutil.copy2(src, dst)
                    actions[-1]["method"] = "copy_fallback"
    return actions


# ── Pretty printing ───────────────────────────────────────────────────────────
def _fmt(mean, std):
    if mean is None or pd.isna(mean):
        return "  -   "
    if std is None or pd.isna(std):
        return f"{mean:.3f}"
    return f"{mean:.3f}±{std:.3f}"


def print_summary(model_name: str, winners_df: pd.DataFrame):
    print("\n" + "=" * 78)
    print(f"  Winners per task -- {model_name}")
    print("=" * 78)
    print(f"  {'task':18s}  {'lr':>7s}  {'drop':>5s}  {'ls':>5s}  "
          f"{'val_bacc':>13s}  {'test_bacc':>13s}  {'test_auc':>13s}  n")
    print("  " + "-" * 92)
    for _, r in winners_df.iterrows():
        print(f"  {r['task']:18s}  {r['lr']:>7.0e}  {r['drop_rate']:>5.2f}  "
              f"{r['label_smoothing']:>5.2f}  "
              f"{_fmt(r['val_bacc_mean'], r['val_bacc_std']):>13s}  "
              f"{_fmt(r['test_bacc_mean'], r['test_bacc_std']):>13s}  "
              f"{_fmt(r['test_auc_mean'], r['test_auc_std']):>13s}  "
              f"{int(r['val_bacc_n']):d}")
    print("=" * 78)


def main():
    args = parse_args()
    model_root = Path(args.model_root)
    out_dir    = Path(args.out_dir) if args.out_dir else model_root / "_winners"

    if not model_root.is_dir():
        print(f"  [ERROR] --model_root does not exist: {model_root}",
              file=sys.stderr)
        sys.exit(1)

    print(f"  Model: {args.model_name}")
    print(f"  Root : {model_root}")
    print(f"  Out  : {out_dir}")
    print(f"  Sel  : {args.selection_metric}  (top_k={args.top_k}, "
          f"{'copy' if args.copy else 'symlink'}{', dry-run' if args.dry_run else ''})")

    runs_df = load_runs(model_root)

    # The selection metric in the config dict is best_val_balanced_acc; our
    # DataFrame column is val_bacc. Keep the CLI flag forward-looking in case
    # we ever switch to ranking by test_bacc -- but warn if the user does.
    if args.selection_metric != "best_val_balanced_acc":
        print(f"  [WARN] Selecting by {args.selection_metric!r} -- if this is "
              f"a test metric you're cherry-picking; report mean test "
              f"accordingly in the writeup.")
    winners_df, top_k_df = select_winners(runs_df, "val_bacc", top_k=args.top_k)

    # Stamp the model name into both CSVs for cross-model concatenation later
    winners_df.insert(0, "model", args.model_name)
    top_k_df.insert(0, "model", args.model_name)

    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)
        winners_path = out_dir / "winners.csv"
        topk_path    = out_dir / f"top_{args.top_k}_per_task.csv"
        pretty_path  = out_dir / "winners_pretty.csv"
        winners_df.to_csv(winners_path, index=False, float_format="%.4f")
        top_k_df.to_csv(topk_path, index=False, float_format="%.4f")
        # Pretty-formatted headline: one row per task, mean±std strings
        # pre-formatted for direct LaTeX / Excel paste. Raw numerics kept
        # in winners.csv so the user can restyle freely.
        pretty = pd.DataFrame({
            "model":           winners_df["model"],
            "task":            winners_df["task"],
            "winning_hp":      winners_df.apply(
                lambda r: f"lr={r['lr']:.0e} d={r['drop_rate']:.2f} ls={r['label_smoothing']:.2f}",
                axis=1),
            "n_seeds":         winners_df["val_bacc_n"].astype("Int64"),
            "val_bacc":        winners_df.apply(
                lambda r: _fmt(r['val_bacc_mean'], r['val_bacc_std']),
                axis=1),
            "test_bacc":       winners_df.apply(
                lambda r: _fmt(r['test_bacc_mean'], r['test_bacc_std']),
                axis=1),
            "test_auc":        winners_df.apply(
                lambda r: _fmt(r['test_auc_mean'], r['test_auc_std']),
                axis=1),
            "test_f1":         winners_df.apply(
                lambda r: _fmt(r['test_f1_mean'], r['test_f1_std']),
                axis=1),
        })
        pretty.to_csv(pretty_path, index=False)
        print(f"  CSV : {winners_path}")
        print(f"  CSV : {topk_path}")
        print(f"  CSV : {pretty_path}  (mean±std strings, paste-ready)")

    actions = emplace_winning_pts(winners_df, runs_df, out_dir,
                                  copy=args.copy, dry_run=args.dry_run)
    n_link = sum(1 for a in actions if a["method"] == "symlink")
    n_copy = sum(1 for a in actions if a["method"] in ("copy", "copy_fallback"))
    if not args.dry_run:
        print(f"  Checkpoints placed: {n_link} symlinks + {n_copy} copies "
              f"({len(actions)} total)")

    print_summary(args.model_name, winners_df)


if __name__ == "__main__":
    main()
