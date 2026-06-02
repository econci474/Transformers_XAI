"""
04e_plot_baseline_postexcl_table.py
===================================
Render publication-quality baseline-model tables for the NO-CDR *post-exclusion*
cohort, with VAL and TEST as SEPARATE tables, for either feature set:

  - full    : whole battery  (CSV in outputs/encoder_only_post_exclusions/)
  - totals  : reduced "filtered" set (CSV in outputs/baseline_no_cdr_post_exclusion/,
              suffix _filtered) — one summary score per neurocognitive test +
              demographics/APOE/CSF.

Both produced by 02h_baseline_post_exclusion.py (which now fits on TRAIN and
evaluates VAL + TEST separately, tagging each row with a `Split` column).

The reduced ("filtered") table additionally prints, at the bottom, the exact list
of features used and the total count (read from feature_set*.json written by 02h).

Layout mirrors 02e (the familiar 5-group baseline table): CN vs MCI+AD, CN/MCI/AD,
AD conv. 3/5/10 yr — but rendered once per split.

Outputs (into outputs/baseline_no_cdr_post_exclusion/):
  baseline_model_table_full_{val,test}.png        (if the full CSV exists)
  baseline_model_table_filtered_{val,test}.png     (if the filtered CSV exists)

Run (local Windows, env: snp/clinical):
  python clinical_pipeline/04e_plot_baseline_postexcl_table.py
"""

import json
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"

# ── Paths ──────────────────────────────────────────────────────────────────────
REPO_ROOT = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
OUTPUTS   = REPO_ROOT / "clinical_pipeline" / "outputs"
OUT_DIR   = OUTPUTS / "baseline_no_cdr_post_exclusion"        # all tables land here
OUT_DIR.mkdir(parents=True, exist_ok=True)
SPLIT_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical"
                 r"\no_cdr_stratified_post_exclusion\tabular\baseline")

# feature_set → (comparison CSV, feature-list JSON, output filename stem)
SOURCES = {
    "full": (OUTPUTS / "encoder_only_post_exclusions" / "baseline_model_comparison.csv",
             OUTPUTS / "encoder_only_post_exclusions" / "feature_set.json",
             "baseline_model_table_full"),
    "totals": (OUT_DIR / "baseline_model_comparison_filtered.csv",
               OUT_DIR / "feature_set_filtered.json",
               "baseline_model_table_filtered"),
}

MODELS       = ["LogReg", "SVM", "XGBoost"]
MODEL_LABELS = {"LogReg": "Log. Reg.", "SVM": "SVM", "XGBoost": "XGBoost"}

# 5 core task groups (match 02e / Image #1). `task_sub` must uniquely substring-match
# the `Task` (description) column written by 02h.
GROUPS_TMPL = [
    {"title": "CN vs MCI+AD (binary)", "task": "vs MCI+AD",
     "metrics": ["Accuracy", "BalancedAcc", "AUC-ROC", "Precision", "Recall", "F1"],
     "headers": ["Acc", "Bal.Acc", "AUC", "Prec", "Rec", "F1"]},
    {"title": "CN / MCI / AD (3-class)", "task": "Multi-class",
     "metrics": ["Accuracy", "BalancedAcc", "AUC-ROC (OvR)", "MacroPrecision", "MacroRecall", "MacroF1"],
     "headers": ["Acc", "Bal.Acc", "AUC", "M.Prec", "M.Rec", "M.F1"]},
    {"title": "AD conv. ≤ 3 yrs", "task": "within 3 years",
     "metrics": ["Accuracy", "BalancedAcc", "AUC-ROC", "Precision", "Recall", "F1"],
     "headers": ["Acc", "Bal.Acc", "AUC", "Prec", "Rec", "F1"]},
    {"title": "AD conv. ≤ 5 yrs", "task": "within 5 years",
     "metrics": ["Accuracy", "BalancedAcc", "AUC-ROC", "Precision", "Recall", "F1"],
     "headers": ["Acc", "Bal.Acc", "AUC", "Prec", "Rec", "F1"]},
    {"title": "AD conv. ≤ 10 yrs", "task": "within 10 years",
     "metrics": ["Accuracy", "BalancedAcc", "AUC-ROC", "Precision", "Recall", "F1"],
     "headers": ["Acc", "Bal.Acc", "AUC", "Prec", "Rec", "F1"]},
]


def split_subtitles(split):
    """Per-group N + class breakdown computed from the seed_0 <split> file."""
    s0 = pd.read_csv(SPLIT_DIR / "seed_0" / f"{split}.csv", low_memory=False)
    out = {}
    # T1 CN vs MCI+AD
    t = s0.dropna(subset=["Label_bl_multi"]).copy()
    t["_y"] = t["Label_bl_multi"].map({"CN": 0, "MCI": 1, "AD": 1})
    t = t.dropna(subset=["_y"])
    out["vs MCI+AD"] = f"N({split})={len(t)}, CN={(t['_y']==0).sum()}, MCI+AD={(t['_y']==1).sum()}"
    # T2 multiclass
    t = s0.dropna(subset=["Label_bl_multi"])
    out["Multi-class"] = (f"N({split})={len(t)}, CN={(t['Label_bl_multi']=='CN').sum()}, "
                          f"MCI={(t['Label_bl_multi']=='MCI').sum()}, AD={(t['Label_bl_multi']=='AD').sum()}")
    # T3 conversion horizons
    for yr, lbl in [(3, "Label_3y"), (5, "Label_5y"), (10, "Label_10y")]:
        t = s0[s0["Label_bl_multi"].isin(["CN", "MCI"])].copy()
        t["_y"] = pd.to_numeric(t[lbl], errors="coerce")
        t = t.dropna(subset=["_y"])
        conv = int(t["_y"].sum())
        out[f"within {yr} years"] = (f"N({split})={len(t)}, non-conv={len(t)-conv}, "
                                     f"AD≤{yr}yr={conv}")
    return out


def parse_mean(cell_str):
    s = str(cell_str).strip()
    if s in ("—", "", "nan", "NaN"):
        return -np.inf
    try:
        return float(s.split("±")[0].strip())
    except (ValueError, IndexError):
        return -np.inf


def render(df_split, groups, split, feat_meta, out_stem):
    """Render one table (one split) and save PNG + PDF."""
    def get_cell(task_sub, model, metric):
        row = df_split[df_split["Task"].str.contains(task_sub, na=False, regex=False)
                       & (df_split["Model"] == model)]
        if row.empty or metric not in row.columns:
            return "—"
        val = str(row.iloc[0][metric]).strip()
        return val if val not in ("", "nan", "NaN") else "—"

    # best model per (task, metric)
    best_cells = set()
    for grp in groups:
        for m in grp["metrics"]:
            vals = {mod: parse_mean(get_cell(grp["task"], mod, m)) for mod in MODELS}
            bm = max(vals, key=vals.get)
            if vals[bm] > -np.inf:
                best_cells.add((grp["task"], m, bm))

    COL_W, ROW_H, ROW_H_SUB, MODEL_COL_W = 1.30, 0.36, 0.24, 1.05
    n_data_cols = sum(len(g["metrics"]) for g in groups)
    n_cols = 1 + n_data_cols

    # footnote text (wrapped) — feature list + count
    foot_lines = []
    if feat_meta:
        feats = feat_meta.get("features", [])
        nfeat = feat_meta.get("n_features", len(feats))
        head = f"Features included ({nfeat}): "
        wrapped = textwrap.wrap(head + ", ".join(feats), width=160)
        foot_lines = wrapped
    foot_h = 0.22 * (len(foot_lines) + 1) if foot_lines else 0.0

    fig_w = 0.3 + MODEL_COL_W + n_data_cols * COL_W + 0.3
    fig_h = 0.25 + ROW_H * (3 + len(MODELS)) + ROW_H_SUB + 0.15 + foot_h

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off"); ax.set_xlim(0, fig_w); ax.set_ylim(0, fig_h)

    LEFT = 0.3
    col_widths = [MODEL_COL_W] + [COL_W] * n_data_cols
    col_left = [LEFT]
    for w in col_widths:
        col_left.append(col_left[-1] + w)
    col_cx = [(col_left[i] + col_left[i + 1]) / 2 for i in range(n_cols)]
    RIGHT = col_left[-1]

    TOP = fig_h - 0.15
    y = TOP
    row_tops = []
    row_tops.append(y); y -= ROW_H       # title
    row_tops.append(y); y -= ROW_H       # group names
    row_tops.append(y); y -= ROW_H_SUB   # subtitle
    row_tops.append(y); y -= ROW_H       # sub-headers
    for _ in MODELS:
        row_tops.append(y); y -= ROW_H
    row_tops.append(y)
    bottom = row_tops[4 + len(MODELS)]

    def row_mid(i): return (row_tops[i] + row_tops[i + 1]) / 2

    def hline(yy, lw=1.0, ls="-"):
        ax.plot([LEFT, RIGHT], [yy, yy], color="black", linewidth=lw, linestyle=ls,
                solid_capstyle="butt", zorder=3)

    hline(row_tops[0], 1.5); hline(row_tops[1], 1.5); hline(row_tops[2], 0.8)
    hline(row_tops[3], 0.5, (0, (3, 3))); hline(row_tops[4], 1.2); hline(bottom, 1.5)

    col_cursor = 1
    for grp in groups:
        x = col_left[col_cursor]
        ax.plot([x, x], [bottom, row_tops[1]], color="black", linewidth=0.6,
                linestyle=(0, (3, 3)), zorder=2)
        col_cursor += len(grp["metrics"])

    split_word = {"val": "Validation", "test": "Test"}.get(split, split.title())
    ax.text((LEFT + RIGHT) / 2, row_mid(0),
            f"Clinical Baseline Models: {split_word} Performance at Baseline — "
            f"No CDR, post-exclusion ({feat_meta.get('feature_set','full') if feat_meta else 'full'} "
            f"features; mean ± std, seeds 0/1/2, train-fit / held-out {split})",
            ha="center", va="center", fontsize=8.0, fontweight="bold")

    col_cursor = 1
    for grp in groups:
        n = len(grp["metrics"]); x0 = col_left[col_cursor]; x1 = col_left[col_cursor + n]
        ax.text((x0 + x1) / 2, row_mid(1), grp["title"], ha="center", va="center",
                fontsize=8, fontweight="bold")
        ax.plot([x0 + 0.05, x1 - 0.05], [row_tops[2], row_tops[2]], color="black",
                linewidth=0.7, zorder=4)
        ax.text((x0 + x1) / 2, row_mid(2), grp["subtitle"], ha="center", va="center",
                fontsize=6.5, fontstyle="italic", color="#555555")
        col_cursor += n

    col_cursor = 1
    for grp in groups:
        for h in grp["headers"]:
            ax.text(col_cx[col_cursor], row_mid(3), h, ha="center", va="center",
                    fontsize=7.5, fontstyle="italic")
            col_cursor += 1

    for r_idx, model in enumerate(MODELS):
        row_i = 4 + r_idx
        ax.text(col_left[0] + 0.08, row_mid(row_i), MODEL_LABELS[model], ha="left",
                va="center", fontsize=8, fontweight="bold")
        col_cursor = 1
        for grp in groups:
            for m in grp["metrics"]:
                val = get_cell(grp["task"], model, m)
                is_best = (grp["task"], m, model) in best_cells
                ax.text(col_cx[col_cursor], row_mid(row_i), val, ha="center", va="center",
                        fontsize=7.8, fontweight="bold" if is_best else "normal")
                col_cursor += 1

    ax.add_patch(plt.Rectangle((LEFT, bottom), RIGHT - LEFT, row_tops[0] - bottom,
                               facecolor="none", edgecolor="black", linewidth=1.5, zorder=5))

    # footnote (feature list)
    if foot_lines:
        fy = bottom - 0.18
        for ln in foot_lines:
            ax.text(LEFT, fy, ln, ha="left", va="top", fontsize=6.2, color="#333333")
            fy -= 0.20

    plt.tight_layout(pad=0.1)
    png = OUT_DIR / f"{out_stem}_{split}.png"
    pdf = OUT_DIR / f"{out_stem}_{split}.pdf"
    fig.savefig(png, bbox_inches="tight", dpi=300)
    fig.savefig(pdf, bbox_inches="tight", dpi=300)
    plt.close()
    print(f"  Saved → {png}")


def main():
    any_done = False
    for fs, (csv_path, json_path, stem) in SOURCES.items():
        if not csv_path.exists():
            print(f"[skip] {fs}: no CSV at {csv_path}")
            continue
        df = pd.read_csv(csv_path)
        if "Split" not in df.columns:
            print(f"[skip] {fs}: {csv_path.name} has no 'Split' column "
                  f"(re-run 02h after the val/test update).")
            continue
        feat_meta = json.loads(json_path.read_text()) if json_path.exists() else None
        for split in ["val", "test"]:
            sub = df[df["Split"] == split]
            if sub.empty:
                print(f"[skip] {fs}/{split}: no rows"); continue
            groups = []
            subs = split_subtitles(split)
            for g in GROUPS_TMPL:
                g2 = dict(g); g2["subtitle"] = subs.get(g["task"], "")
                groups.append(g2)
            print(f"[{fs}/{split}]")
            render(sub, groups, split, feat_meta, stem)
            any_done = True
    if not any_done:
        print("Nothing rendered — run 02h_baseline_post_exclusion.py first "
              "(--feature_set full and/or totals).")


if __name__ == "__main__":
    main()
