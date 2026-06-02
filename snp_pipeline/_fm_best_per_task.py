"""
_fm_best_per_task.py
====================
Cross-task FM-variant summary helper for the new-task leaderboards
(49a/b/c outputs at outputs/new_tasks/...).

Reads each task's `leaderboard.csv` (49a sCN_vs_pCN, 49b horizon_binary
filtered by task column, 49c horizon_T4_multiclass), parses each row's
`source` column to recover (model_family, attention_mode, attention_source),
then emits:

  outputs/new_tasks/_fm_best_per_task/
    sCN_vs_pCN.tsv              # top-3 FM rows by val metric
    horizon_3y.tsv              # ...
    horizon_5y.tsv              # ...
    horizon_7y.tsv              # ...
    horizon_10y.tsv             # ...
    T4_3class.tsv               # ...
    cross_task_summary.tsv      # wide: variant x task, val_bacc_mean
    cross_task_summary.png      # heatmap of cross_task_summary.tsv

Attention source is recovered from the `source` column's `[attn=...]`
annotation. Wave 1 rows carry `[attn=T1_AD_CN]` or `[attn=fixed_aggregator]`;
Wave 2 rows carry `[attn=<task_name>]`. Wave 2's task-matched attention
is the "task-conditioned" comparison the user asked for; Wave 1 sits
alongside as the Task-1-biased reference.

Usage:
    python snp_pipeline/_fm_best_per_task.py \\
        --ld-config maf1_geno02_r2_0p8                 # optional; defaults to DEFAULT_LD_CONFIG
        --beta-source raw                              # or _filtered / _withdesikan
        --metric val_bacc                              # or val_auc / val_macro_balAcc

Skips a task gracefully if its leaderboard.csv doesn't exist yet.
"""
from __future__ import annotations
import argparse
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _prs_strict_qc_lib import DEFAULT_LD_CONFIG  # noqa: E402

NEW_TASKS_ROOT = Path("D:/ADNI_SNP_Omni2.5M_20140220/outputs/new_tasks")
OUT_ROOT       = NEW_TASKS_ROOT / "_fm_best_per_task"

# (task_label, leaderboard subdir, optional `task` column filter)
TASKS_BINARY = [
    ("sCN_vs_pCN",   "sCN_vs_pCN",       None),
    ("horizon_3y",   "horizon_binary",   "horizon_3y"),
    ("horizon_5y",   "horizon_binary",   "horizon_5y"),
    ("horizon_7y",   "horizon_binary",   "horizon_7y"),
    ("horizon_10y",  "horizon_binary",   "horizon_10y"),
]
TASKS_MULTI = [
    ("T4_3class",    "horizon_T4_multiclass", None),
]

# Order tasks chronologically in the cross-task summary.
TASK_ORDER = ["sCN_vs_pCN", "horizon_3y", "horizon_5y",
               "horizon_7y", "horizon_10y", "T4_3class"]


def _leaderboard_path(subdir: str, ld_config: str, beta_source: str) -> Path:
    if beta_source == "raw":
        return NEW_TASKS_ROOT / subdir / ld_config / "leaderboard.csv"
    return NEW_TASKS_ROOT / f"{subdir}__{beta_source}" / "leaderboard.csv"


_RE_ATTN = re.compile(r"\[attn=([^\]]+)\]")


def _parse_source(s: str) -> dict:
    """Recover (arm, model_family, attn_mode, attn_source) from a leaderboard
    source string. PRS rows return all-None except `arm='PRS'`. FM rows are
    annotated `<MODEL> <mode> chrom_hier <head> [attn=<src>]`."""
    out = {"arm": None, "model_family": None,
            "attn_mode": None, "attn_source": None}
    m = _RE_ATTN.search(s)
    if m is None:
        out["arm"] = "PRS"
        return out
    attn = m.group(1)
    out["attn_source"] = attn
    # Wave 1 attention sources: T1_AD_CN, fixed_aggregator
    # Wave 2 attention sources: task names (sCN_vs_pCN, horizon_5y, T4_3class, ...)
    out["arm"] = "FM_wave1" if attn in ("T1_AD_CN", "fixed_aggregator") else "FM_wave2"
    s_low = s.lower()
    if "bmfm-snp" in s_low or "bmfm_snp" in s_low:
        out["model_family"] = "BMFM-SNP"
    elif "ntv2" in s_low:
        out["model_family"] = "NTv2"
    for mode in ("attn_bias_per_modality", "attn_bias_single", "none"):
        if mode in s_low:
            out["attn_mode"] = mode
            break
    return out


def _read_one(label: str, subdir: str, task_filter: str | None,
                ld_config: str, beta_source: str,
                metric: str) -> pd.DataFrame | None:
    p = _leaderboard_path(subdir, ld_config, beta_source)
    if not p.exists():
        print(f"  [task={label}] missing leaderboard: {p}")
        return None
    df = pd.read_csv(p)
    if task_filter is not None and "task" in df.columns:
        df = df[df["task"].astype(str) == task_filter].copy()
        if df.empty:
            print(f"  [task={label}] no rows for task={task_filter}")
            return None
    df["_task_label"] = label
    parsed = df["source"].astype(str).apply(_parse_source)
    df["arm"]          = [d["arm"] for d in parsed]
    df["model_family"] = [d["model_family"] for d in parsed]
    df["attn_mode"]    = [d["attn_mode"] for d in parsed]
    df["attn_source"]  = [d["attn_source"] for d in parsed]
    if metric not in df.columns:
        # T4 multiclass leaderboard uses 'val_bacc_macro_mean' (macro suffix);
        # binary uses 'val_bacc_mean'. Alias whichever is present so the
        # cross-task pivot can use a single column name.
        for alt in ("val_bacc_macro_mean", "val_macro_bacc_mean",
                     "val_macro_balAcc_mean", "val_bacc_mean"):
            if alt in df.columns:
                df[metric] = df[alt]; break
    return df


def _per_task_best_fm(df: pd.DataFrame, metric: str, top_n: int = 3
                       ) -> pd.DataFrame:
    """Keep only FM rows (wave 1 + 2), sort by `metric` descending, return
    the top `top_n`. Columns kept are minimal + auditable."""
    fm = df[df["arm"].astype(str).str.startswith("FM")].copy()
    if fm.empty: return fm
    fm = fm.sort_values(metric, ascending=False).head(top_n)
    keep = ["source", "covar_mode", "arm", "model_family", "attn_mode",
             "attn_source"]
    # Pull whichever val/test metric cols are present.
    metric_cols = [c for c in fm.columns
                    if c.endswith("_mean") and (c.startswith("val_") or c.startswith("test_"))]
    return fm[keep + metric_cols].reset_index(drop=True)


def _cross_task_summary(all_df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Pivot: rows = (arm, model_family, attn_mode, covar_mode);
    cols = task; values = `metric`. Sorted by mean across tasks (desc)."""
    fm = all_df[all_df["arm"].astype(str).str.startswith("FM")].copy()
    if fm.empty:
        return pd.DataFrame()
    fm["_variant"] = (fm["model_family"].astype(str) + " | "
                       + fm["attn_mode"].astype(str) + " | "
                       + fm["covar_mode"].astype(str) + " | "
                       + fm["arm"].astype(str))
    pv = fm.pivot_table(index="_variant", columns="_task_label",
                          values=metric, aggfunc="mean")
    pv = pv.reindex(columns=[t for t in TASK_ORDER if t in pv.columns])
    pv["_row_mean"] = pv.mean(axis=1, skipna=True)
    pv = pv.sort_values("_row_mean", ascending=False)
    return pv


def _render_heatmap(pv: pd.DataFrame, out_png: Path, metric: str,
                      ld_config: str, beta_source: str) -> None:
    """One-panel heatmap of the cross-task summary. Rows annotated with
    arm tag, cell text = `metric` (3 dp)."""
    if pv.empty:
        print(f"  [heatmap] no FM rows to render")
        return
    data = pv.drop(columns=["_row_mean"], errors="ignore")
    fig_h = max(4.5, 0.30 * len(data) + 2.0)
    fig_w = max(8.0, 1.0 + 1.4 * len(data.columns))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    arr = data.values.astype(float)
    im = ax.imshow(arr, aspect="auto", cmap="viridis",
                     vmin=np.nanmin(arr), vmax=np.nanmax(arr))
    ax.set_xticks(range(len(data.columns)))
    ax.set_xticklabels(data.columns, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(len(data.index)))
    ax.set_yticklabels(data.index, fontsize=7)
    for i in range(len(data.index)):
        for j in range(len(data.columns)):
            v = arr[i, j]
            if np.isnan(v): continue
            txt = f"{v:.3f}"
            color = "white" if v < (np.nanmin(arr) + np.nanmax(arr))/2 else "black"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7, color=color)
    plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label=metric)
    ax.set_title(
        f"FM variant best-by-task summary -- new-task leaderboards\n"
        f"metric = {metric}   LD = {ld_config}   beta = {beta_source}\n"
        f"rows: <model_family> | <attn_mode> | <covar_mode> | <arm>;  "
        f"arm = FM_wave1 (Task-1 attn) | FM_wave2 (task-conditioned attn)",
        fontsize=10, loc="left")
    plt.tight_layout()
    plt.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote: {out_png}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ld-config",   default=DEFAULT_LD_CONFIG)
    ap.add_argument("--beta-source", default="raw")
    ap.add_argument("--metric",      default="val_bacc_mean",
                      help="leaderboard column to rank by (binary: val_bacc_mean; "
                           "T4: val_macro_balAcc_mean — falls back automatically)")
    ap.add_argument("--top-n",       type=int, default=3)
    args = ap.parse_args()

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    print(f"[fm_best_per_task] reading from {NEW_TASKS_ROOT}")
    print(f"  ld_config={args.ld_config}  beta_source={args.beta_source}  metric={args.metric}\n")

    all_frames: list[pd.DataFrame] = []
    metric_binary = args.metric
    # 49c's multiclass leaderboard uses val_bacc_macro_mean (note: macro suffix
    # comes AFTER bacc), not val_macro_balAcc_mean. Auto-fallback if the user
    # passed the binary metric name.
    metric_multi  = args.metric if "macro" in args.metric else "val_bacc_macro_mean"

    for label, subdir, task_filter in TASKS_BINARY:
        df = _read_one(label, subdir, task_filter,
                        args.ld_config, args.beta_source, metric_binary)
        if df is None: continue
        best = _per_task_best_fm(df, metric_binary, args.top_n)
        out = OUT_ROOT / f"{label}.tsv"
        best.to_csv(out, sep="\t", index=False)
        print(f"  wrote: {out}  (n={len(best)} top-FM rows)")
        all_frames.append(df.assign(_metric=df[metric_binary]))

    for label, subdir, task_filter in TASKS_MULTI:
        df = _read_one(label, subdir, task_filter,
                        args.ld_config, args.beta_source, metric_multi)
        if df is None: continue
        best = _per_task_best_fm(df, metric_multi, args.top_n)
        out = OUT_ROOT / f"{label}.tsv"
        best.to_csv(out, sep="\t", index=False)
        print(f"  wrote: {out}  (n={len(best)} top-FM rows)")
        all_frames.append(df.assign(_metric=df[metric_multi]))

    if not all_frames:
        print("[fm_best_per_task] no leaderboards found -- did you run 49a/b/c?")
        return

    all_df = pd.concat(all_frames, ignore_index=True, sort=False)
    pv = _cross_task_summary(all_df, "_metric")
    if not pv.empty:
        tsv = OUT_ROOT / "cross_task_summary.tsv"
        pv.to_csv(tsv, sep="\t")
        print(f"\n  wrote: {tsv}  ({len(pv)} variant rows x {len(pv.columns)-1} tasks)")
        png = OUT_ROOT / "cross_task_summary.png"
        _render_heatmap(pv, png, args.metric, args.ld_config, args.beta_source)


if __name__ == "__main__":
    main()
