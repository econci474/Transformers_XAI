"""
plot_baseline_table.py
======================
Renders a publication-quality table matching the style of the reference image:
  - No shading or fills anywhere
  - All text black
  - Thick horizontal rules: top, below group headers, below sub-headers, bottom
  - Dotted vertical lines between column groups
  - Bold group names, italic sub-headers, bold model labels

Outputs:
  outputs/baseline_model_table.pdf
  outputs/baseline_model_table.png
  outputs/baseline_model_table.tex
"""

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path

matplotlib.rcParams["font.family"] = "DejaVu Serif"

# ── Paths ──────────────────────────────────────────────────────────────────────
OUT_DIR = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\clinical_pipeline\outputs")
CSV     = OUT_DIR / "baseline_model_comparison.csv"
df      = pd.read_csv(CSV)

# ── Table layout ───────────────────────────────────────────────────────────────
GROUPS = [
    {
        "title":   "CN vs MCI+AD (binary)",
        "task":    "Binary",
        "metrics": ["Accuracy", "BalancedAcc", "AUC-ROC", "F1"],
        "headers": ["Acc", "Bal.Acc", "AUC", "F1"],
    },
    {
        "title":   "CN / MCI / AD (3-class)",
        "task":    "Multi-class",
        "metrics": ["Accuracy", "BalancedAcc", "AUC-ROC (OvR)", "MacroF1"],
        "headers": ["Acc", "Bal.Acc", "AUC", "MacroF1"],
    },
    {
        "title":   "AD conversion \u2264 3 yrs",
        "task":    "3 years",
        "metrics": ["Accuracy", "BalancedAcc", "AUC-ROC", "F1"],
        "headers": ["Acc", "Bal.Acc", "AUC", "F1"],
    },
    {
        "title":   "AD conversion \u2264 5 yrs",
        "task":    "5 years",
        "metrics": ["Accuracy", "BalancedAcc", "AUC-ROC", "F1"],
        "headers": ["Acc", "Bal.Acc", "AUC", "F1"],
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

# ── Sizing ─────────────────────────────────────────────────────────────────────
COL_W       = 1.45   # inches per metric column
ROW_H       = 0.36   # inches per row
MODEL_COL_W = 1.05   # inches for the model-name column

n_data_cols = sum(len(g["metrics"]) for g in GROUPS)
n_cols      = 1 + n_data_cols

fig_w = 0.3 + MODEL_COL_W + n_data_cols * COL_W + 0.3
fig_h = 0.25 + ROW_H * (3 + len(MODELS)) + 0.15   # title + group + sub + data rows

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.axis("off")
ax.set_xlim(0, fig_w)
ax.set_ylim(0, fig_h)

# ── Column x-centres and boundaries ───────────────────────────────────────────
LEFT = 0.3
col_widths = [MODEL_COL_W] + [COL_W] * n_data_cols
col_left   = [LEFT]
for w in col_widths:
    col_left.append(col_left[-1] + w)
col_cx = [(col_left[i] + col_left[i + 1]) / 2 for i in range(n_cols)]
RIGHT = col_left[-1]

# ── Row y-midpoints (top to bottom) ───────────────────────────────────────────
TOP = fig_h - 0.15
row_tops = [TOP - i * ROW_H for i in range(3 + len(MODELS) + 1)]
# row_tops[0..N] are the TOP edges; midpoints = (row_tops[i]+row_tops[i+1])/2
def row_mid(i): return (row_tops[i] + row_tops[i + 1]) / 2

# Row indices: 0=title, 1=group names, 2=sub-headers, 3+=data rows

# ── Horizontal rules ───────────────────────────────────────────────────────────
def hline(y, lw=1.0, ls="-"):
    ax.plot([LEFT, RIGHT], [y, y], color="black", linewidth=lw, linestyle=ls,
            solid_capstyle="butt", zorder=3)

hline(row_tops[0], lw=1.5)   # top border
hline(row_tops[1], lw=1.5)   # below title
hline(row_tops[2], lw=0.8)   # below group names
hline(row_tops[3], lw=1.2)   # below sub-headers (midrule)
hline(row_tops[3 + len(MODELS)], lw=1.5)  # bottom border

# ── Dotted vertical lines between column groups ────────────────────────────────
col_cursor = 1
for grp in GROUPS:
    x = col_left[col_cursor]
    ax.plot([x, x], [row_tops[3 + len(MODELS)], row_tops[1]],
            color="black", linewidth=0.6, linestyle=(0, (3, 3)), zorder=2)
    col_cursor += len(grp["metrics"])

# ── Title row (row 0) ─────────────────────────────────────────────────────────
ax.text((LEFT + RIGHT) / 2, row_mid(0),
        "Baseline Clinical Models: Test Performance   "
        "(mean \u00b1 std, seeds 0/1/2, 80/10/10 split)",
        ha="center", va="center", fontsize=8.5, fontweight="bold", color="black")

# ── Group header row (row 1) ──────────────────────────────────────────────────
col_cursor = 1
for grp in GROUPS:
    n   = len(grp["metrics"])
    x0  = col_left[col_cursor]
    x1  = col_left[col_cursor + n]
    xm  = (x0 + x1) / 2
    ax.text(xm, row_mid(1), grp["title"],
            ha="center", va="center", fontsize=8, fontweight="bold", color="black")
    # Short underline beneath the group title (like \cmidrule)
    pad = 0.05
    ax.plot([x0 + pad, x1 - pad], [row_tops[2], row_tops[2]],
            color="black", linewidth=0.7, zorder=4)
    col_cursor += n

# Model column group header
ax.text(col_cx[0], row_mid(1), "", ha="center", va="center")

# ── Sub-header row (row 2) ────────────────────────────────────────────────────
col_cursor = 1
for grp in GROUPS:
    for h in grp["headers"]:
        ax.text(col_cx[col_cursor], row_mid(2), h,
                ha="center", va="center", fontsize=7.5,
                fontstyle="italic", color="black")
        col_cursor += 1

# ── Data rows (rows 3+) ───────────────────────────────────────────────────────
for r_idx, model in enumerate(MODELS):
    row_i = 3 + r_idx
    # Model label (left-aligned, bold)
    ax.text(col_left[0] + 0.08, row_mid(row_i), MODEL_LABELS[model],
            ha="left", va="center", fontsize=8, fontweight="bold", color="black")
    # Data cells
    col_cursor = 1
    for grp in GROUPS:
        for m in grp["metrics"]:
            val = get_cell(grp["task"], model, m)
            ax.text(col_cx[col_cursor], row_mid(row_i), val,
                    ha="center", va="center", fontsize=7.8, color="black")
            col_cursor += 1

# ── Outer box ─────────────────────────────────────────────────────────────────
rect = plt.Rectangle((LEFT, row_tops[3 + len(MODELS)]),
                      RIGHT - LEFT, row_tops[0] - row_tops[3 + len(MODELS)],
                      facecolor="none", edgecolor="black", linewidth=1.5, zorder=5)
ax.add_patch(rect)

# ── Save ──────────────────────────────────────────────────────────────────────
plt.tight_layout(pad=0.1)
OUT_DIR.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_DIR / "baseline_model_table.pdf", bbox_inches="tight", dpi=300)
fig.savefig(OUT_DIR / "baseline_model_table.png", bbox_inches="tight", dpi=300)
print(f"Saved PDF \u2192 {OUT_DIR / 'baseline_model_table.pdf'}")
print(f"Saved PNG \u2192 {OUT_DIR / 'baseline_model_table.png'}")
plt.close()

# ── LaTeX booktabs ────────────────────────────────────────────────────────────
def latex_val(s):
    return str(s).replace("\u00b1", r"$\pm$").replace("\u2014", r"---")

col_spec = "l" + "".join(
    "@{\\hspace{3pt}}|@{\\hspace{3pt}}" + "c" * len(g["metrics"])
    for g in GROUPS
)
lines = [
    r"\begin{table}[ht]",
    r"\centering",
    r"\caption{Baseline Clinical Models: Test Performance (mean $\pm$ std, 3 seeds $\times$ 80/10/10 split)}",
    r"\label{tab:baseline_models}",
    r"\small\setlength{\tabcolsep}{5pt}",
    r"\begin{tabular}{" + col_spec + "}",
    r"\toprule",
]
grp_row = " & " + " & ".join(
    r"\multicolumn{" + str(len(g["metrics"])) + r"}{c}{\textbf{" + g["title"] + r"}}"
    for g in GROUPS
) + r" \\"
lines.append(grp_row)
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
            parts.append(latex_val(get_cell(grp["task"], model, m)))
    lines.append(" & ".join(parts) + r" \\")
lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
tex_path = OUT_DIR / "baseline_model_table.tex"
tex_path.write_text("\n".join(lines), encoding="utf-8")
print(f"Saved LaTeX \u2192 {tex_path}")
