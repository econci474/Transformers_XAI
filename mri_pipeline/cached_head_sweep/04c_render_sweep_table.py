"""
04c_render_sweep_table.py — aggregate & render the BrainDINO T2 wide-sweep results
=================================================================================
The wide frozen linear-probe sweep (04_braindino_head_sweep_T2_wide_submit_csd3.sh)
runs 192 HP configs x 3 seeds = 576 cells, each writing metrics.json under
    <sweep-root>/BrainDINO_vitb16_frozen_cached/T2_multiclass/seed_<S>/<hp>/metrics.json
where <hp> encodes the non-default knobs (lr/drop/ls + optional _wd/_hHEAD/_LOSS/_std).

This script reads every metrics.json, aggregates each config across its seeds
(mean +/- std for val-bACC, val-AUC, test-bACC, test-AUC, test-F1, test-acc, with
per-seed columns kept for traceability), ranks by **mean val-bACC** (the agreed
selection criterion), and writes:

  braindino_T2_sweep_3seed.csv      — full per-config table (all configs, all metrics).
  braindino_T2_sweep_summary.png    — styled top-N table + baseline + paper reference rows
                                      (winner row bolded).

Selection = mean val-bACC across seeds; the winner's test metrics are reported for
comparison against the pre-sweep baseline and the paper's macro-AUC 0.954.

Read-only w.r.t. the sweep tree. Run:
  python mri_pipeline/cached_head_sweep/04c_render_sweep_table.py
  python mri_pipeline/cached_head_sweep/04c_render_sweep_table.py \
      --sweep-root /home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/aug_none_T2_wide
"""

from __future__ import annotations

import argparse
import json
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent  # .../Transformers_XAI

DEFAULT_SWEEP_ROOT = Path(
    r"D:/ADNI_BIDS_project/derivatives/braindino_outputs/aug_none_T2_wide"
)
DEFAULT_OUT_DIR = REPO / "mri_pipeline" / "outputs" / "braindino_sweep"
SUBDIR = "BrainDINO_vitb16_frozen_cached/T2_multiclass"

# Pre-sweep BrainDINO T2 frozen cached-head baseline (mlp default head) and the
# published reference, both for the summary table (from cross_model_table_test.csv
# / BrainDINO_2026.pdf macro-OvR AUC, 3-class frozen).
BASELINE = dict(val_bacc=0.602, test_bacc=0.453, test_auc=0.648, test_f1=np.nan)
PAPER_AUC = 0.954

N_SUMMARY = 15  # rows shown in the PNG


def _load(sweep_root: Path) -> pd.DataFrame:
    pat = str(sweep_root / SUBDIR / "seed_*" / "*" / "metrics.json")
    recs = []
    for f in glob.glob(pat):
        with open(f) as fh:
            d = json.load(fh)
        c, t = d["config"], d["test_metrics"]
        recs.append(dict(
            hp=os.path.basename(os.path.dirname(f)),
            seed=int(c["seed"]),
            head=c.get("head"), loss=c.get("loss"), standardize=c.get("standardize"),
            lr=c.get("lr"), weight_decay=c.get("weight_decay"),
            drop_rate=c.get("drop_rate"), label_smoothing=c.get("label_smoothing"),
            best_epoch=c.get("best_epoch"),
            val_bacc=c.get("best_val_balanced_acc"), val_auc=c.get("best_val_auc"),
            test_bacc=t.get("balanced_acc"), test_auc=t.get("auc_roc_ovr"),
            test_f1=t.get("macro_f1"), test_acc=t.get("accuracy"),
        ))
    if not recs:
        raise SystemExit(f"No metrics.json found under {pat}")
    return pd.DataFrame(recs)


def _aggregate(df: pd.DataFrame) -> pd.DataFrame:
    """One row per hp config: mean/std across seeds + per-seed columns."""
    metric_cols = ["val_bacc", "val_auc", "test_bacc", "test_auc", "test_f1", "test_acc"]
    keep_cfg = ["head", "loss", "standardize", "lr", "weight_decay",
                "drop_rate", "label_smoothing"]
    out_rows = []
    for hp, g in df.groupby("hp"):
        g = g.sort_values("seed")
        row = {"hp": hp, "n_seeds": int(g["seed"].nunique())}
        cfg0 = g.iloc[0]
        for k in keep_cfg:
            row[k] = cfg0[k]
        for m in metric_cols:
            vals = g[m].astype(float).values
            row[f"{m}_mean"] = float(np.nanmean(vals)) if len(vals) else np.nan
            row[f"{m}_std"] = float(np.nanstd(vals)) if len(vals) else np.nan
        # per-seed traceability columns
        for s in sorted(g["seed"].unique()):
            gs = g[g.seed == s].iloc[0]
            row[f"val_bacc_s{s}"] = float(gs["val_bacc"])
            row[f"test_auc_s{s}"] = float(gs["test_auc"])
        out_rows.append(row)
    agg = pd.DataFrame(out_rows)
    agg = agg.sort_values("val_bacc_mean", ascending=False).reset_index(drop=True)
    agg.insert(0, "rank", np.arange(1, len(agg) + 1))
    return agg


def _cell(m, s):
    if pd.isna(m):
        return "—"
    return f"{m:.3f} ± {s:.3f}" if pd.notna(s) else f"{m:.3f}"


def render_summary(agg: pd.DataFrame, out_png: Path, n_show: int, n_total: int,
                   n_seeds_mode: int):
    cols = ["#", "HP config (val-bACC ranked)", "val-bACC", "test-bACC",
            "test-AUC", "test-F1"]
    body, kinds = [], []  # kinds: 'data'|'baseline'|'paper' for styling
    for _, r in agg.head(n_show).iterrows():
        body.append([
            f"{int(r['rank'])}", r["hp"],
            _cell(r["val_bacc_mean"], r["val_bacc_std"]),
            _cell(r["test_bacc_mean"], r["test_bacc_std"]),
            _cell(r["test_auc_mean"], r["test_auc_std"]),
            _cell(r["test_f1_mean"], r["test_f1_std"]),
        ])
        kinds.append("data")
    # reference rows
    body.append(["—", "Pre-sweep baseline (mlp default head)",
                 f"{BASELINE['val_bacc']:.3f}", f"{BASELINE['test_bacc']:.3f}",
                 f"{BASELINE['test_auc']:.3f}", "—"])
    kinds.append("baseline")
    body.append(["—", "Paper (BrainDINO, macro-OvR AUC, 3-class frozen)",
                 "—", "—", f"{PAPER_AUC:.3f}", "—"])
    kinds.append("paper")

    fig_h = 0.9 + 0.34 * (len(body) + 1)
    fig, ax = plt.subplots(figsize=(11.5, fig_h))
    ax.axis("off")
    tab = ax.table(cellText=body, colLabels=cols, loc="upper center",
                   cellLoc="center", colLoc="center")
    tab.auto_set_font_size(False)
    tab.set_fontsize(8.5)
    tab.scale(1.0, 1.30)
    n_cols = len(cols)
    # header
    for j in range(n_cols):
        tab[0, j].set_facecolor("#D9E1F2")
        tab[0, j].set_text_props(weight="bold")
    # left-align the config-name column
    for i in range(1, len(body) + 1):
        tab[i, 1].set_text_props(ha="left")
    # row styling
    for i, kind in enumerate(kinds, start=1):
        if kind == "data" and i == 1:        # winner row
            for j in range(n_cols):
                tab[i, j].set_facecolor("#E2EFDA")
                tab[i, j].set_text_props(weight="bold")
        elif kind == "baseline":
            for j in range(n_cols):
                tab[i, j].set_facecolor("#FFF2CC")
        elif kind == "paper":
            for j in range(n_cols):
                tab[i, j].set_facecolor("#FCE4D6")
                tab[i, j].set_text_props(style="italic")

    title = ("BrainDINO T2 (CN/MCI/AD) wide frozen linear-probe sweep — "
             f"top {n_show} of {n_total} configs\n"
             f"selection = mean val-bACC across {n_seeds_mode} seeds "
             "(mean ± std). Winner highlighted; baseline & paper for reference.")
    plt.title(title, fontsize=11, pad=10)
    ax.text(0.0, -0.02, "bACC=balanced accuracy · AUC=macro-OvR ROC-AUC · "
            "F1=macro-F1 (all on TEST unless prefixed val-). _std=feature "
            "standardization, hlinear=linear head, focal=focal loss.",
            transform=ax.transAxes, ha="left", va="top", fontsize=7.5,
            family="monospace")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"  wrote {out_png}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sweep-root", default=str(DEFAULT_SWEEP_ROOT),
                    help="aug_none_T2_wide root (contains BrainDINO_vitb16_frozen_cached/).")
    ap.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    ap.add_argument("--n-show", type=int, default=N_SUMMARY,
                    help="rows in the summary PNG (default 15).")
    args = ap.parse_args()

    sweep_root = Path(args.sweep_root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    raw = _load(sweep_root)
    n_per_seed = raw.groupby("seed")["hp"].nunique().to_dict()
    print(f"Loaded {len(raw)} runs; configs/seed = {n_per_seed}")
    agg = _aggregate(raw)
    n_seeds_mode = int(raw["seed"].nunique())

    # full CSV
    csv_path = out_dir / "braindino_T2_sweep_3seed.csv"
    agg.to_csv(csv_path, index=False, float_format="%.4f")
    print(f"  wrote {csv_path}  ({len(agg)} configs)")

    # summary PNG
    render_summary(agg, out_dir / "braindino_T2_sweep_summary.png",
                   n_show=args.n_show, n_total=len(agg), n_seeds_mode=n_seeds_mode)

    # stdout verdict
    w = agg.iloc[0]
    print("\n" + "=" * 70)
    print(f"WINNER (mean val-bACC, {n_seeds_mode} seeds): {w['hp']}")
    print(f"  val-bACC  {w['val_bacc_mean']:.3f} ± {w['val_bacc_std']:.3f}")
    print(f"  test-bACC {w['test_bacc_mean']:.3f} ± {w['test_bacc_std']:.3f}"
          f"   (baseline {BASELINE['test_bacc']:.3f})")
    print(f"  test-AUC  {w['test_auc_mean']:.3f} ± {w['test_auc_std']:.3f}"
          f"   (baseline {BASELINE['test_auc']:.3f}; paper {PAPER_AUC:.3f})")
    print(f"  test-F1   {w['test_f1_mean']:.3f} ± {w['test_f1_std']:.3f}")
    print("=" * 70)


if __name__ == "__main__":
    main()
