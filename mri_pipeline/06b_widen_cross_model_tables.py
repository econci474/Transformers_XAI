"""
06b_widen_cross_model_tables.py — val CSV + wide (bACC/AUC/F1) cross-model tables
=================================================================================
`06_render_cross_model_table.py` only ever persisted the TEST long-form CSV
(`cross_model_table_test.csv`); the val table was rendered as PNG/TeX only, so
`cross_model_table_val.csv` was never written. This script reads the existing
`cross_model_table_test.csv` (the aggregated snapshot — it already carries test
balanced_acc/auc/f1 + val_bacc per row) and produces, WITHOUT re-reading the
per-model metrics.json trees (which live on CSD3):

  cross_model_table_val.csv          — long-form val view (val_bacc per row).
  cross_model_table_test_wide.{png,csv} — per task: bACC / AUC / F1 sub-columns.
  cross_model_table_val_wide.png     — per task: val bACC (+ val AUC where present;
                                       old runs only tracked val_bacc, so AUC/F1
                                       are '—' until runs that log them land).

Run:  python mri_pipeline/06b_widen_cross_model_tables.py
      python mri_pipeline/06b_widen_cross_model_tables.py --csv <path to test csv>
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DEFAULT_CSV = HERE / "outputs" / "cross_model" / "cross_model_table_test.csv"

TASK_ORDER = ["T1_binary", "T1b_binary", "T1c_binary", "T1d_binary", "T2_multiclass"]
TASK_LABELS = {"T1_binary": "T1 CN|MCI+AD", "T1b_binary": "T1b CN+MCI|AD",
               "T1c_binary": "T1c CN|AD", "T1d_binary": "T1d pMCI|sMCI",
               "T2_multiclass": "T2 CN/MCI/AD"}
# Row ordering (mirrors 06_render_cross_model_table._build_pivot).
MODEL_RANK = {"BrainDINO": 0, "BrainMVP": 2, "ViT-MAE75": 4, "ViT-scratch": 6,
              "AG-MS3D-vanilla": 7, "AG-MS3D-sep": 8, "Spasov-CNN": 9}
VARIANT_RANK = {"full_ft": 0, "lora": 1, "frozen": 2, "vanilla": 3,
                "separable": 4, "baseline": 5, "scratch": 5, "agms3d": 6, "-": 9}
AUG_RANK = {"plus_original": 0, "stochastic": 1, "random": 1, "none": 2, "-": 9}


def _clean_aug(a):
    return str(a).replace("�", "").strip() if pd.notna(a) else "-"


def _sort_key(r):
    a = _clean_aug(r["augment"])
    return (MODEL_RANK.get(r["model"], 99), VARIANT_RANK.get(r["variant"], 99),
            AUG_RANK.get(a, 99), r["model"], str(r["variant"]), a)


def _cell(m, s):
    if pd.isna(m):
        return "—"
    return f"{m:.3f} ± {s:.3f}" if pd.notna(s) else f"{m:.3f}"


def _row_label(r):
    a = _clean_aug(r["augment"])
    base = f"{r['model']} / {r['variant']}"
    return base if a in ("-", "nan") else f"{base} / {a}"


def render_wide(df, metrics, out_png, title, ncol_note=""):
    """metrics = list of (col_in_df_prefix, short_label). Per task, one sub-column
    per metric. df is one row per (model,variant,augment,task)."""
    rows = (df[["model", "variant", "augment"]].drop_duplicates()
            .apply(lambda r: r, axis=1))
    units = sorted({(r.model, r.variant, r.augment) for r in df.itertuples()},
                   key=lambda t: _sort_key({"model": t[0], "variant": t[1], "augment": t[2]}))
    n_metrics = len(metrics)
    # Build body + a parallel numeric matrix for per-column bolding.
    body, numeric = [], []
    for (m, v, a) in units:
        sub = df[(df.model == m) & (df.variant == v) & (df.augment == a)]
        rl = _row_label({"model": m, "variant": v, "augment": a})
        cells, nums = [rl], []
        for task in TASK_ORDER:
            c = sub[sub.task == task]
            for pref, _ in metrics:
                if c.empty or pd.isna(c.iloc[0].get(f"{pref}_mean")):
                    cells.append("—"); nums.append(np.nan)
                else:
                    cells.append(_cell(c.iloc[0][f"{pref}_mean"], c.iloc[0].get(f"{pref}_std")))
                    nums.append(float(c.iloc[0][f"{pref}_mean"]))
        body.append(cells); numeric.append(nums)
    numeric = np.array(numeric, dtype=float)

    sub_labels = ["Model / variant / aug"] + [lbl for _ in TASK_ORDER for _, lbl in metrics]
    n_cols = len(sub_labels)
    fig_w = 2.6 + 0.95 * n_metrics * len(TASK_ORDER)
    fig, ax = plt.subplots(figsize=(fig_w, 0.7 + 0.3 * (len(body) + 2)))
    ax.axis("off")
    tab = ax.table(cellText=body, colLabels=sub_labels, loc="upper center",
                   cellLoc="center", colLoc="center")
    tab.auto_set_font_size(False); tab.set_fontsize(7.5); tab.scale(1.0, 1.25)
    # group header: task names spanning their metric sub-columns (drawn as text)
    for j in range(n_cols):
        tab[0, j].set_facecolor("#ECECEC"); tab[0, j].set_text_props(weight="bold")
    for i in range(1, len(body) + 1):
        tab[i, 0].set_text_props(ha="left")
    # bold best per metric-column (skip the label col 0)
    for j in range(1, n_cols):
        col = numeric[:, j - 1]
        if not np.all(np.isnan(col)):
            bi = int(np.nanargmax(col))
            tab[bi + 1, j].set_text_props(weight="bold")
    # task group labels above the sub-columns
    task_band = " | ".join(f"{TASK_LABELS[t]} [{'/'.join(l for _, l in metrics)}]"
                           for t in TASK_ORDER)
    plt.title(f"{title}\n{task_band}", pad=8, fontsize=10)
    if ncol_note:
        ax.text(0.0, -0.04, ncol_note, transform=ax.transAxes, ha="left", va="top",
                fontsize=7, family="monospace")
    fig.savefig(out_png, dpi=180, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"  wrote {out_png}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=str(DEFAULT_CSV),
                    help="cross_model_table_test.csv (aggregated long-form source).")
    args = ap.parse_args()
    src = Path(args.csv)
    if not src.exists():
        raise SystemExit(f"source not found: {src}\nRun 06_render_cross_model_table.py first.")
    df = pd.read_csv(src)
    df = df[df.task.isin(TASK_ORDER)].copy()
    out = src.parent

    # 1) reproduce the val CSV (long-form val view).
    val_cols = ["model", "variant", "augment", "task", "n_seeds", "n_degenerate",
                "val_bacc_mean", "val_bacc_std"]
    df[val_cols].to_csv(out / "cross_model_table_val.csv", index=False, float_format="%.4f")
    print(f"  wrote {out/'cross_model_table_val.csv'}")

    # 2) wide TEST table: bACC / AUC / F1.
    test_metrics = [("balanced_acc", "bACC"), ("auc", "AUC"), ("f1", "F1")]
    df[["model", "variant", "augment", "task", "n_seeds",
        "balanced_acc_mean", "balanced_acc_std", "auc_mean", "auc_std",
        "f1_mean", "f1_std"]].to_csv(out / "cross_model_table_test_wide.csv",
                                     index=False, float_format="%.4f")
    print(f"  wrote {out/'cross_model_table_test_wide.csv'}")
    render_wide(df, test_metrics, out / "cross_model_table_test_wide.png",
                "MRI models — TEST bACC / AUC / F1 by task (mean ± std across seeds)",
                ncol_note="Best per metric-column in bold. bACC=balanced accuracy, "
                          "AUC=macro-OvR ROC-AUC, F1=macro-F1 (all on TEST).")

    # 3) wide VAL table: only val_bacc tracked historically -> AUC/F1 are '—'
    #    until runs that log them land (the cached-head trainer now logs val_auc).
    val_metrics = [("val_bacc", "bACC"), ("val_auc", "AUC"), ("val_f1", "F1")]
    render_wide(df, val_metrics, out / "cross_model_table_val_wide.png",
                "MRI models — VAL bACC / AUC / F1 by task (mean ± std across seeds)",
                ncol_note="VAL AUC/F1 not recorded for these runs (old trainers tracked only "
                          "val balanced accuracy); they populate from runs that log them "
                          "(cached-head trainer now logs val_auc).")
    print("Done.")


if __name__ == "__main__":
    main()
