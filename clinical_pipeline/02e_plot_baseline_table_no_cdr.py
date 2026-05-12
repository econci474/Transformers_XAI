"""
02e_plot_baseline_table_no_cdr.py
==================================
Renders a publication-quality table for the no-CDR baseline models.
Same style as 02b_plot_baseline_table.py but reads from baseline_no_cdr/.

Outputs:
  outputs/baseline_no_cdr/baseline_model_table_no_cdr.pdf
  outputs/baseline_no_cdr/baseline_model_table_no_cdr.png
  outputs/baseline_no_cdr/baseline_model_table_no_cdr.tex

Run:
  python clinical_pipeline/02e_plot_baseline_table_no_cdr.py
"""

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path

matplotlib.rcParams["font.family"] = "DejaVu Serif"

# ── Paths ──────────────────────────────────────────────────────────────────────
REPO_ROOT = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
OUT_DIR   = REPO_ROOT / "clinical_pipeline" / "outputs" / "baseline_no_cdr"
CSV       = OUT_DIR / "baseline_model_comparison.csv"
SPLIT_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified\tabular\baseline")
df        = pd.read_csv(CSV)

# ── Compute test-set class distributions (seed 0) ─────────────────────────────
_s0_test = pd.read_csv(SPLIT_DIR / "seed_0" / "test.csv", low_memory=False)

_t1 = _s0_test.dropna(subset=["Label_bl_multi"]).copy()
_t1["_y"] = _t1["Label_bl_multi"].map({"CN": 0, "MCI": 1, "AD": 1})
_t1 = _t1.dropna(subset=["_y"])
t1_n = len(_t1); t1_cn = int((_t1["_y"] == 0).sum()); t1_mciad = int((_t1["_y"] == 1).sum())

_t2 = _s0_test.dropna(subset=["Label_bl_multi"]).copy()
t2_n = len(_t2)
t2_cn  = int((_t2["Label_bl_multi"] == "CN").sum())
t2_mci = int((_t2["Label_bl_multi"] == "MCI").sum())
t2_ad  = int((_t2["Label_bl_multi"] == "AD").sum())

_t3a = _s0_test[_s0_test["Label_bl_multi"].isin(["CN", "MCI"])].copy()
_t3a["_y"] = pd.to_numeric(_t3a["Label_3y"], errors="coerce")
_t3a = _t3a.dropna(subset=["_y"])
t3a_n = len(_t3a); t3a_conv = int(_t3a["_y"].sum()); t3a_nonconv = t3a_n - t3a_conv

_t3b = _s0_test[_s0_test["Label_bl_multi"].isin(["CN", "MCI"])].copy()
_t3b["_y"] = pd.to_numeric(_t3b["Label_5y"], errors="coerce")
_t3b = _t3b.dropna(subset=["_y"])
t3b_n = len(_t3b); t3b_conv = int(_t3b["_y"].sum()); t3b_nonconv = t3b_n - t3b_conv

_t3c = _s0_test[_s0_test["Label_bl_multi"].isin(["CN", "MCI"])].copy()
_t3c["_y"] = pd.to_numeric(_t3c["Label_10y"], errors="coerce")
_t3c = _t3c.dropna(subset=["_y"])
t3c_n = len(_t3c); t3c_conv = int(_t3c["_y"].sum()); t3c_nonconv = t3c_n - t3c_conv

# ── Table layout ───────────────────────────────────────────────────────────────
GROUPS = [
    {
        "title":   "CN vs MCI+AD (binary)",
        "subtitle": f"N(test)={t1_n}, CN={t1_cn}, MCI+AD={t1_mciad}",
        "task":    "Binary",
        "metrics": ["Accuracy", "BalancedAcc", "AUC-ROC", "Precision", "Recall", "F1"],
        "headers": ["Acc", "Bal.Acc", "AUC", "Prec", "Rec", "F1"],
    },
    {
        "title":   "CN / MCI / AD (3-class)",
        "subtitle": f"N(test)={t2_n}, CN={t2_cn}, MCI={t2_mci}, AD={t2_ad}",
        "task":    "Multi-class",
        "metrics": ["Accuracy", "BalancedAcc", "AUC-ROC (OvR)", "MacroPrecision", "MacroRecall", "MacroF1"],
        "headers": ["Acc", "Bal.Acc", "AUC", "M.Prec", "M.Rec", "M.F1"],
    },
    {
        "title":   "AD conv. \u2264 3 yrs",
        "subtitle": f"N(test)={t3a_n}, non-conv={t3a_nonconv}, AD\u22643yr={t3a_conv}",
        "task":    "3 years",
        "metrics": ["Accuracy", "BalancedAcc", "AUC-ROC", "Precision", "Recall", "F1"],
        "headers": ["Acc", "Bal.Acc", "AUC", "Prec", "Rec", "F1"],
    },
    {
        "title":   "AD conv. \u2264 5 yrs",
        "subtitle": f"N(test)={t3b_n}, non-conv={t3b_nonconv}, AD\u22645yr={t3b_conv}",
        "task":    "5 years",
        "metrics": ["Accuracy", "BalancedAcc", "AUC-ROC", "Precision", "Recall", "F1"],
        "headers": ["Acc", "Bal.Acc", "AUC", "Prec", "Rec", "F1"],
    },
    {
        "title":   "AD conv. \u2264 10 yrs",
        "subtitle": f"N(test)={t3c_n}, non-conv={t3c_nonconv}, AD\u226410yr={t3c_conv}",
        "task":    "10 years",
        "metrics": ["Accuracy", "BalancedAcc", "AUC-ROC", "Precision", "Recall", "F1"],
        "headers": ["Acc", "Bal.Acc", "AUC", "Prec", "Rec", "F1"],
    },
]

MODELS       = ["LogReg", "SVM", "XGBoost"]
MODEL_LABELS = {"LogReg": "Log. Reg.", "SVM": "SVM", "XGBoost": "XGBoost"}

def get_cell(task_sub, model, metric):
    row = df[df["Task"].str.contains(task_sub, na=False) & (df["Model"] == model)]
    if row.empty or metric not in row.columns:
        return "\u2014"
    val = str(row.iloc[0][metric]).strip()
    return val if val not in ("", "nan", "NaN") else "\u2014"

def parse_mean(cell_str):
    """Extract the mean from a 'mean \u00b1 std' string."""
    s = str(cell_str).strip()
    if s in ("\u2014", "", "nan"):
        return -np.inf
    try:
        return float(s.split("\u00b1")[0].strip())
    except (ValueError, IndexError):
        return -np.inf

# Pre-compute which (task, metric) → best model
best_cells = set()  # set of (task_sub, metric, model)
for grp in GROUPS:
    for m in grp["metrics"]:
        vals = {model: parse_mean(get_cell(grp["task"], model, m)) for model in MODELS}
        best_model = max(vals, key=vals.get)
        if vals[best_model] > -np.inf:
            best_cells.add((grp["task"], m, best_model))

# ── Sizing ─────────────────────────────────────────────────────────────────────
COL_W       = 1.30
ROW_H       = 0.36
ROW_H_SUB   = 0.24
MODEL_COL_W = 1.05

n_data_cols = sum(len(g["metrics"]) for g in GROUPS)
n_cols      = 1 + n_data_cols

fig_w = 0.3 + MODEL_COL_W + n_data_cols * COL_W + 0.3
fig_h = 0.25 + ROW_H * (3 + len(MODELS)) + ROW_H_SUB + 0.15

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.axis("off")
ax.set_xlim(0, fig_w)
ax.set_ylim(0, fig_h)

LEFT = 0.3
col_widths = [MODEL_COL_W] + [COL_W] * n_data_cols
col_left   = [LEFT]
for w in col_widths:
    col_left.append(col_left[-1] + w)
col_cx = [(col_left[i] + col_left[i + 1]) / 2 for i in range(n_cols)]
RIGHT = col_left[-1]

TOP = fig_h - 0.15
y_cursor = TOP
row_tops = []
row_tops.append(y_cursor); y_cursor -= ROW_H      # title
row_tops.append(y_cursor); y_cursor -= ROW_H      # group names
row_tops.append(y_cursor); y_cursor -= ROW_H_SUB  # subtitle
row_tops.append(y_cursor); y_cursor -= ROW_H      # sub-headers
for _ in MODELS:
    row_tops.append(y_cursor); y_cursor -= ROW_H
row_tops.append(y_cursor)

def row_mid(i): return (row_tops[i] + row_tops[i + 1]) / 2

def hline(y, lw=1.0, ls="-"):
    ax.plot([LEFT, RIGHT], [y, y], color="black", linewidth=lw, linestyle=ls,
            solid_capstyle="butt", zorder=3)

hline(row_tops[0], lw=1.5)
hline(row_tops[1], lw=1.5)
hline(row_tops[2], lw=0.8)
hline(row_tops[3], lw=0.5, ls=(0, (3, 3)))
hline(row_tops[4], lw=1.2)
hline(row_tops[4 + len(MODELS)], lw=1.5)

# Dotted verticals
col_cursor = 1
for grp in GROUPS:
    x = col_left[col_cursor]
    ax.plot([x, x], [row_tops[4 + len(MODELS)], row_tops[1]],
            color="black", linewidth=0.6, linestyle=(0, (3, 3)), zorder=2)
    col_cursor += len(grp["metrics"])

# Title
ax.text((LEFT + RIGHT) / 2, row_mid(0),
        "Clinical Models: Test Performance at Baseline — No CDR   "
        "(mean \u00b1 std, seeds 0/1/2, stratified 80/10/10 split)",
        ha="center", va="center", fontsize=8.5, fontweight="bold", color="black")

# Group headers
col_cursor = 1
for grp in GROUPS:
    n   = len(grp["metrics"])
    x0  = col_left[col_cursor]
    x1  = col_left[col_cursor + n]
    xm  = (x0 + x1) / 2
    ax.text(xm, row_mid(1), grp["title"],
            ha="center", va="center", fontsize=8, fontweight="bold", color="black")
    pad = 0.05
    ax.plot([x0 + pad, x1 - pad], [row_tops[2], row_tops[2]],
            color="black", linewidth=0.7, zorder=4)
    col_cursor += n

# Subtitle
col_cursor = 1
for grp in GROUPS:
    n   = len(grp["metrics"])
    x0  = col_left[col_cursor]
    x1  = col_left[col_cursor + n]
    xm  = (x0 + x1) / 2
    ax.text(xm, row_mid(2), grp["subtitle"],
            ha="center", va="center", fontsize=6.5,
            fontstyle="italic", color="#555555")
    col_cursor += n

# Sub-headers
col_cursor = 1
for grp in GROUPS:
    for h in grp["headers"]:
        ax.text(col_cx[col_cursor], row_mid(3), h,
                ha="center", va="center", fontsize=7.5,
                fontstyle="italic", color="black")
        col_cursor += 1

# Data rows
for r_idx, model in enumerate(MODELS):
    row_i = 4 + r_idx
    ax.text(col_left[0] + 0.08, row_mid(row_i), MODEL_LABELS[model],
            ha="left", va="center", fontsize=8, fontweight="bold", color="black")
    col_cursor = 1
    for grp in GROUPS:
        for m in grp["metrics"]:
            val = get_cell(grp["task"], model, m)
            is_best = (grp["task"], m, model) in best_cells
            ax.text(col_cx[col_cursor], row_mid(row_i), val,
                    ha="center", va="center", fontsize=7.8, color="black",
                    fontweight="bold" if is_best else "normal")
            col_cursor += 1

# Outer box
rect = plt.Rectangle((LEFT, row_tops[4 + len(MODELS)]),
                      RIGHT - LEFT, row_tops[0] - row_tops[4 + len(MODELS)],
                      facecolor="none", edgecolor="black", linewidth=1.5, zorder=5)
ax.add_patch(rect)

# Save
plt.tight_layout(pad=0.1)
fig.savefig(OUT_DIR / "baseline_model_table_no_cdr.pdf", bbox_inches="tight", dpi=300)
fig.savefig(OUT_DIR / "baseline_model_table_no_cdr.png", bbox_inches="tight", dpi=300)
print(f"Saved PDF \u2192 {OUT_DIR / 'baseline_model_table_no_cdr.pdf'}")
print(f"Saved PNG \u2192 {OUT_DIR / 'baseline_model_table_no_cdr.png'}")
plt.close()

# ── LaTeX booktabs ────────────────────────────────────────────────────────────
def latex_val(s):
    return str(s).replace("\u00b1", r"$\pm$").replace("\u2014", r"---")

col_spec = "l" + "".join(
    r"@{\hspace{3pt}}|@{\hspace{3pt}}" + "c" * len(g["metrics"])
    for g in GROUPS
)
lines = [
    r"\begin{table}[ht]",
    r"\centering",
    r"\caption{Clinical Models (No CDR): Test Performance at Baseline (mean $\pm$ std, 3 seeds $\times$ stratified 80/10/10 split)}",
    r"\label{tab:baseline_models_no_cdr}",
    r"\small\setlength{\tabcolsep}{5pt}",
    r"\begin{tabular}{" + col_spec + "}",
    r"\toprule",
]
grp_row = " & " + " & ".join(
    r"\multicolumn{" + str(len(g["metrics"])) + r"}{c}{\textbf{" + g["title"] + r"}}"
    for g in GROUPS
) + r" \\"
lines.append(grp_row)

sub_row = " & " + " & ".join(
    r"\multicolumn{" + str(len(g["metrics"])) + r"}{c}{\textit{\scriptsize " + g["subtitle"] + r"}}"
    for g in GROUPS
) + r" \\"
lines.append(sub_row)

cursor = 2
cmidrules = ""
for g in GROUPS:
    n = len(g["metrics"])
    cmidrules += r"\cmidrule(lr){" + f"{cursor}-{cursor+n-1}" + "}"
    cursor += n
lines.append(cmidrules)
lines.append("Model & " + " & ".join(
    r"\textit{" + h + r"}" for g in GROUPS for h in g["headers"]
) + r" \\")
lines.append(r"\midrule")
for model in MODELS:
    parts = [r"\textbf{" + MODEL_LABELS[model] + r"}"]
    for grp in GROUPS:
        for m in grp["metrics"]:
            cell = latex_val(get_cell(grp["task"], model, m))
            if (grp["task"], m, model) in best_cells:
                cell = r"\textbf{" + cell + r"}"
            parts.append(cell)
    lines.append(" & ".join(parts) + r" \\")
lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
tex_path = OUT_DIR / "baseline_model_table_no_cdr.tex"
tex_path.write_text("\n".join(lines), encoding="utf-8")
print(f"Saved LaTeX \u2192 {tex_path}")
