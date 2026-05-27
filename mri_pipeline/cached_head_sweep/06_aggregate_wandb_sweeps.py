"""
06_aggregate_wandb_sweeps.py
============================
Post-pull aggregator for the 5 wandb_*.csv files produced by
06_fetch_wandb_tables.py. Picks the winning HP combo per (project, task)
by best_val_bacc (mean across seeds), and renders:

  wandb_winners.csv         long-form: one row per (project, task) =
                            winning HP combo + best_val_bacc mean / std / n
  wandb_top5_per_task.csv   long-form: top-5 HP combos per (project, task)
                            for honest variance / "median of top-K" reporting
  wandb_winners.png         publication-style table -- rows = project,
                            cols = task, cell = "bacc ± std (n)" + winning
                            HPs underneath

Selection metric: summary/best_val_bacc (the best val_bacc reached across
training, the same field 05b uses from on-disk metrics.json). Falls back
to summary/val_bacc if best is missing.

Caveat: W&B only has val_bacc for these projects (val_auc/f1/sens/spec etc.
were not logged in earlier runs -- today's trainer change adds them for
future runs). For test metrics, you still need the on-disk metrics.json
files via 05b_select_best_hp_per_task.py.

Usage:
    python mri_pipeline/cached_head_sweep/06_aggregate_wandb_sweeps.py
    python mri_pipeline/cached_head_sweep/06_aggregate_wandb_sweeps.py \
        --csv-dir mri_pipeline/cached_head_sweep \
        --out-dir mri_pipeline/cached_head_sweep \
        --projects vit_mae_frozen_cached braindino_frozen_cached
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_PROJECTS = [
    "vit_mae_frozen_cached",
    "brainmvp_frozen_cached",
    "braindino_frozen_cached",
    "brainmvp_debug",
    "vit_debugging",
]

TASK_ORDER = [
    "T1_binary",
    "T1b_binary",
    "T1c_binary",
    "T1d_binary",
    "T2_multiclass",
]
TASK_LABELS = {
    "T1_binary":     "T1: CN vs MCI+AD",
    "T1b_binary":    "T1b: CN+MCI vs AD",
    "T1c_binary":    "T1c: CN vs AD",
    "T1d_binary":    "T1d: pMCI vs sMCI",
    "T2_multiclass": "T2: CN / MCI / AD",
}

# Config columns we treat as identity/metadata (NOT HPs to group by).
# Anything else under config/ that varies across runs is an HP.
_IDENTITY_COLS = {
    "config/task", "config/seed",
    "config/embeddings_pt", "config/pretrained_ckpt",
    "config/n_train", "config/n_val", "config/n_test",
    "config/n_params", "config/num_labels",
    "config/embed_dim", "config/model", "config/model_name",
    "config/epochs", "config/patience", "config/long_mode",
    "config/max_months", "config/session_policy",
    "config/strategy",         # captured by project name
    "config/augment",          # captured by project name (all cached are aug=none)
    "config/aug_copies",
    "config/warmup_epochs",
    "config/grad_accum_steps",
}


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    here = os.path.dirname(os.path.abspath(__file__))
    p.add_argument("--csv-dir", default=here,
                   help="Directory containing wandb_*.csv files.")
    p.add_argument("--out-dir", default=here,
                   help="Where to write the winners CSV/PNG.")
    p.add_argument("--projects", nargs="+", default=DEFAULT_PROJECTS,
                   help=f"Project names to aggregate. Default: {DEFAULT_PROJECTS}")
    p.add_argument("--metric", default="summary/best_val_bacc",
                   help="Column to select winners by (default best_val_bacc).")
    p.add_argument("--fallback-metric", default="summary/val_bacc",
                   help="Used when --metric is missing for some rows.")
    return p.parse_args()


def _identify_hp_cols(df: pd.DataFrame) -> list[str]:
    """HP columns = config/ columns that vary across runs and aren't in the
    identity denylist. We want a small, comparable set across projects."""
    hp_cols = []
    for c in df.columns:
        if not c.startswith("config/"):
            continue
        if c in _IDENTITY_COLS:
            continue
        if df[c].nunique(dropna=True) <= 1:
            continue
        hp_cols.append(c)
    return sorted(hp_cols)


def _hp_label(row: pd.Series, hp_cols: list[str]) -> str:
    """One-line render of the HP combo for the PNG cell."""
    parts = []
    for c in hp_cols:
        v = row[c]
        short = c.replace("config/", "")
        if isinstance(v, float):
            if np.isnan(v):
                continue
            parts.append(f"{short}={v:g}")
        else:
            parts.append(f"{short}={v}")
    return " ".join(parts)


def _aggregate_project(df: pd.DataFrame, metric: str, fallback: str,
                       project: str):
    """Return (winners_df, top5_df) for one project.

    winners_df : one row per task; the highest-mean HP combo.
    top5_df    : top 5 HP combos per task (long-form).
    """
    # Resolve the metric (fall back to summary/val_bacc per-row).
    if metric in df.columns:
        mvals = pd.to_numeric(df[metric], errors="coerce")
    else:
        mvals = pd.Series(np.nan, index=df.index)
    if fallback in df.columns:
        mvals = mvals.fillna(pd.to_numeric(df[fallback], errors="coerce"))
    df = df.copy()
    df["_metric"] = mvals

    drop_mask = df["_metric"].isna() | df["config/task"].isna() \
                | df["config/seed"].isna()
    if drop_mask.any():
        print(f"    dropping {int(drop_mask.sum())}/{len(df)} rows "
              f"with missing task/seed/metric")
    df = df[~drop_mask].copy()
    if df.empty:
        return pd.DataFrame(), pd.DataFrame()

    hp_cols = _identify_hp_cols(df)
    if not hp_cols:
        # No HPs vary -- everything's the same combo. Group by task only.
        hp_cols = []

    group_cols = ["config/task"] + hp_cols
    agg = (df.groupby(group_cols, dropna=False)["_metric"]
             .agg(["mean", "std", "count"])
             .reset_index()
             .rename(columns={"mean": "metric_mean",
                              "std":  "metric_std",
                              "count": "n_seeds"}))
    # NaN std for n=1 -> 0 (for cosmetics)
    agg["metric_std"] = agg["metric_std"].fillna(0.0)

    winners_rows = []
    top5_rows = []
    for task, sub in agg.groupby("config/task"):
        sub = sub.sort_values("metric_mean", ascending=False)
        # Winner = top row by mean
        win = sub.iloc[0]
        winners_rows.append({
            "project":     project,
            "task":        task,
            "metric_mean": float(win["metric_mean"]),
            "metric_std":  float(win["metric_std"]),
            "n_seeds":     int(win["n_seeds"]),
            "hp_combo":    _hp_label(win, hp_cols),
            **{c.replace("config/", "hp_"): win[c] for c in hp_cols},
        })
        for rank, (_, r) in enumerate(sub.head(5).iterrows(), start=1):
            top5_rows.append({
                "project":     project,
                "task":        task,
                "rank":        rank,
                "metric_mean": float(r["metric_mean"]),
                "metric_std":  float(r["metric_std"]),
                "n_seeds":     int(r["n_seeds"]),
                "hp_combo":    _hp_label(r, hp_cols),
                **{c.replace("config/", "hp_"): r[c] for c in hp_cols},
            })

    return pd.DataFrame(winners_rows), pd.DataFrame(top5_rows)


def _render_winners_png(winners: pd.DataFrame, path: Path, metric_label: str):
    """Publication-style table: rows = project, cols = task, cell = bacc + HPs."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARN] matplotlib not installed -- skipping PNG.")
        return

    projects = winners["project"].drop_duplicates().tolist()
    n_rows = len(projects)
    n_cols = len(TASK_ORDER) + 1

    fig, ax = plt.subplots(
        figsize=(2.2 + 3.4 * len(TASK_ORDER), 0.55 + 0.85 * (n_rows + 1))
    )
    ax.axis("off")

    col_labels = ["Project"] + [TASK_LABELS[t] for t in TASK_ORDER]

    # Build cell text + remember numeric for bold-best.
    cell_text = []
    num_grid = np.full((n_rows, len(TASK_ORDER)), np.nan)
    for i, project in enumerate(projects):
        row = [project]
        sub = winners[winners["project"] == project].set_index("task")
        for j, task in enumerate(TASK_ORDER):
            if task in sub.index:
                r = sub.loc[task]
                txt = (f"{r['metric_mean']:.3f} ± {r['metric_std']:.3f} "
                       f"(n={int(r['n_seeds'])})\n{r['hp_combo']}")
                num_grid[i, j] = r["metric_mean"]
            else:
                txt = "—"
            row.append(txt)
        cell_text.append(row)

    tab = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        loc="center",
        cellLoc="center",
        colLoc="center",
    )
    tab.auto_set_font_size(False)
    tab.set_fontsize(8)
    tab.scale(1.0, 2.0)

    # Bold best per task column
    for j in range(len(TASK_ORDER)):
        col = num_grid[:, j]
        if np.all(np.isnan(col)):
            continue
        i_best = int(np.nanargmax(col))
        tab[i_best + 1, j + 1].set_text_props(weight="bold")

    # Header styling
    for j in range(n_cols):
        tab[0, j].set_facecolor("#ECECEC")
        tab[0, j].set_text_props(weight="bold")
    for i in range(n_rows + 1):
        tab[i, 0].set_text_props(ha="left")
        tab[i, 0].PAD = 0.04

    plt.title(f"Cached-head W&B sweeps -- winners per task "
              f"(by {metric_label}; mean ± std across seeds)",
              pad=14, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  PNG       : {path}")


def main():
    args = parse_args()
    csv_dir = Path(args.csv_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("  06_aggregate_wandb_sweeps -- winners per project x task")
    print(f"  CSV dir  : {csv_dir}")
    print(f"  Out dir  : {out_dir}")
    print(f"  Metric   : {args.metric}  (fallback: {args.fallback_metric})")
    print("=" * 78)

    all_winners = []
    all_top5 = []
    for project in args.projects:
        path = csv_dir / f"wandb_{project}.csv"
        if not path.exists():
            print(f"  [SKIP] {project} -- {path} not found")
            continue
        df = pd.read_csv(path)
        print(f"  [ok]   {project}  -- {len(df)} runs")
        winners, top5 = _aggregate_project(
            df, args.metric, args.fallback_metric, project)
        if winners.empty:
            print(f"        no winnable rows after filtering")
            continue
        all_winners.append(winners)
        all_top5.append(top5)

    if not all_winners:
        sys.exit("No winners aggregated. Check that wandb_<project>.csv exists.")

    winners_df = pd.concat(all_winners, ignore_index=True)
    top5_df = pd.concat(all_top5, ignore_index=True)

    # Stable project ordering for the PNG: input list order, but only kept
    # if any winners landed.
    keep = winners_df["project"].drop_duplicates().tolist()
    project_order = [p for p in args.projects if p in keep]
    winners_df["project"] = pd.Categorical(
        winners_df["project"], categories=project_order, ordered=True)
    winners_df = winners_df.sort_values(["project", "task"]).reset_index(drop=True)
    top5_df["project"] = pd.Categorical(
        top5_df["project"], categories=project_order, ordered=True)
    top5_df = top5_df.sort_values(
        ["project", "task", "rank"]).reset_index(drop=True)

    winners_csv = out_dir / "wandb_winners.csv"
    winners_df.to_csv(winners_csv, index=False, float_format="%.4f")
    print(f"  CSV       : {winners_csv}")

    top5_csv = out_dir / "wandb_top5_per_task.csv"
    top5_df.to_csv(top5_csv, index=False, float_format="%.4f")
    print(f"  CSV       : {top5_csv}")

    metric_label = args.metric.split("/")[-1]
    _render_winners_png(winners_df, out_dir / "wandb_winners.png", metric_label)

    # Stdout summary so the user can read it without opening the PNG
    print("\n" + "=" * 78)
    print("  Per-task ranking summary")
    print("=" * 78)
    for task in TASK_ORDER:
        sub = winners_df[winners_df["task"] == task].sort_values(
            "metric_mean", ascending=False)
        if sub.empty:
            continue
        print(f"\n  {task}")
        print(f"    {'project':35s} {'n':>2s} {'val_bacc':>16s}   HPs")
        print("    " + "-" * 78)
        for _, r in sub.iterrows():
            ba = f"{r['metric_mean']:.3f}±{r['metric_std']:.3f}"
            print(f"    {str(r['project']):35s} {int(r['n_seeds']):>2d} "
                  f"{ba:>16s}   {r['hp_combo']}")
    print()


if __name__ == "__main__":
    main()
