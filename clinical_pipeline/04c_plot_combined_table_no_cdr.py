"""
04c_plot_combined_table_no_cdr.py
==================================
Same as 04b_plot_combined_table.py but uses:
  - No-CDR stratified tabular baselines
  - No-CDR stratified encoder-only results

Outputs go to outputs/encoder-only_no_cdr/

Usage:
  python clinical_pipeline/04c_plot_combined_table_no_cdr.py
"""

import json
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path

matplotlib.rcParams["font.family"] = "DejaVu Serif"

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE      = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\clinical_pipeline")
OUT_DIR   = BASE / "outputs" / "encoder-only_no_cdr"
BL_CSV    = BASE / "outputs" / "baseline_no_cdr" / "baseline_model_comparison.csv"
ENC_DIR   = BASE / "outputs" / "encoder-only_no_cdr"
SPLIT_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified\tabular\baseline")

# ── 1. Load tabular baselines ─────────────────────────────────────────────────
bl_df = pd.read_csv(BL_CSV)

# ── 1b. Compute test-set class distributions (seed 0) ─────────────────────────
_s0_test = pd.read_csv(SPLIT_DIR / "seed_0" / "test.csv", low_memory=False)

_t1 = _s0_test.dropna(subset=["Label_bl_multi"]).copy()
_t1["_y"] = _t1["Label_bl_multi"].map({"CN": 0, "MCI": 1, "AD": 1})
_t1 = _t1.dropna(subset=["_y"])
t1_n = len(_t1); t1_cn = int((_t1["_y"]==0).sum()); t1_mciad = int((_t1["_y"]==1).sum())

_t2 = _s0_test.dropna(subset=["Label_bl_multi"]).copy()
t2_n = len(_t2)
t2_cn = int((_t2["Label_bl_multi"]=="CN").sum())
t2_mci = int((_t2["Label_bl_multi"]=="MCI").sum())
t2_ad = int((_t2["Label_bl_multi"]=="AD").sum())

_t3a = _s0_test[_s0_test["Label_bl_multi"].isin(["CN","MCI"])].copy()
_t3a["_y"] = pd.to_numeric(_t3a["Label_3y"], errors="coerce")
_t3a = _t3a.dropna(subset=["_y"])
t3a_n = len(_t3a); t3a_conv = int(_t3a["_y"].sum()); t3a_nc = t3a_n - t3a_conv

_t3b = _s0_test[_s0_test["Label_bl_multi"].isin(["CN","MCI"])].copy()
_t3b["_y"] = pd.to_numeric(_t3b["Label_5y"], errors="coerce")
_t3b = _t3b.dropna(subset=["_y"])
t3b_n = len(_t3b); t3b_conv = int(_t3b["_y"].sum()); t3b_nc = t3b_n - t3b_conv

_t3c = _s0_test[_s0_test["Label_bl_multi"].isin(["CN","MCI"])].copy()
_t3c["_y"] = pd.to_numeric(_t3c["Label_10y"], errors="coerce")
_t3c = _t3c.dropna(subset=["_y"])
t3c_n = len(_t3c); t3c_conv = int(_t3c["_y"].sum()); t3c_nc = t3c_n - t3c_conv

# ── 2. Load encoder metrics and aggregate across seeds ─────────────────────────
enc_rows = []
n_skipped = 0
for jf in ENC_DIR.rglob("metrics.json"):
    try:
        with open(jf) as f:
            data = json.load(f)
    except (json.JSONDecodeError, ValueError):
        print(f"  [WARN] Skipping corrupt/empty file: {jf}")
        n_skipped += 1
        continue
    cfg = data["config"]
    met = data["test_metrics"]
    enc_rows.append({
        "model":    cfg["model_id"].split("/")[-1],
        "task":     cfg["task"],
        "strategy": cfg["strategy"],
        "seed":     cfg["seed"],
        **met,
    })
if n_skipped:
    print(f"  Skipped {n_skipped} corrupt metrics.json files")

if not enc_rows:
    print("[ERROR] No metrics.json files found in encoder-only_no_cdr/")
    print(f"  Searched: {ENC_DIR}")
    exit(1)

enc_df = pd.DataFrame(enc_rows)
print(f"  Found {len(enc_rows)} encoder metric files")
print(f"  Models: {sorted(enc_df['model'].unique())}")
print(f"  Tasks:  {sorted(enc_df['task'].unique())}")

# Aggregate: mean ± std per (model, task, strategy)
meta = ["model", "task", "strategy", "seed"]
metric_cols = [c for c in enc_df.columns if c not in meta
               and c not in ("loss", "runtime", "samples_per_second",
                             "steps_per_second", "epoch")]

enc_summary = {}
for (model, task, strategy), grp in enc_df.groupby(["model", "task", "strategy"]):
    entry = {}
    for m in metric_cols:
        vals = grp[m].dropna()
        if len(vals) > 0:
            mu = vals.mean()
            sd = vals.std() if len(vals) > 1 else 0.0
            entry[m] = f"{mu:.4f} ± {sd:.4f}"
        else:
            entry[m] = "—"
    enc_summary[(model, task, strategy)] = entry

# Save raw and summary CSVs
enc_df.to_csv(OUT_DIR / "encoder_results_per_run.csv", index=False)
summary_rows = []
for (model, task, strategy), entry in enc_summary.items():
    summary_rows.append({"model": model, "task": task, "strategy": strategy, **entry})
pd.DataFrame(summary_rows).to_csv(OUT_DIR / "encoder_summary.csv", index=False)

# ── Table layout ───────────────────────────────────────────────────────────────
GROUPS = [
    {
        "title":        "CN vs MCI+AD (binary)",
        "subtitle":     f"N(test)={t1_n}, CN={t1_cn}, MCI+AD={t1_mciad}",
        "bl_task":      "Binary",
        "enc_task":     "T1_binary",
        "bl_metrics":   ["Accuracy", "BalancedAcc", "AUC-ROC", "F1"],
        "enc_metrics":  ["accuracy", "balanced_acc", "auc_roc", "f1"],
        "headers":      ["Acc", "Bal.Acc", "AUC", "F1"],
    },
    {
        "title":        "CN / MCI / AD (3-class)",
        "subtitle":     f"N(test)={t2_n}, CN={t2_cn}, MCI={t2_mci}, AD={t2_ad}",
        "bl_task":      "Multi-class",
        "enc_task":     "T2_multiclass",
        "bl_metrics":   ["Accuracy", "BalancedAcc", "AUC-ROC (OvR)", "MacroF1"],
        "enc_metrics":  ["accuracy", "balanced_acc", "auc_roc_ovr", "macro_f1"],
        "headers":      ["Acc", "Bal.Acc", "AUC", "MacroF1"],
    },
    {
        "title":        "AD conversion \u2264 3 yrs",
        "subtitle":     f"N(test)={t3a_n}, non-conv={t3a_nc}, AD\u22643yr={t3a_conv}",
        "bl_task":      "3 years",
        "enc_task":     "T3a_conv3y",
        "bl_metrics":   ["Accuracy", "BalancedAcc", "AUC-ROC", "F1"],
        "enc_metrics":  ["accuracy", "balanced_acc", "auc_roc", "f1"],
        "headers":      ["Acc", "Bal.Acc", "AUC", "F1"],
    },
    {
        "title":        "AD conversion \u2264 5 yrs",
        "subtitle":     f"N(test)={t3b_n}, non-conv={t3b_nc}, AD\u22645yr={t3b_conv}",
        "bl_task":      "5 years",
        "enc_task":     "T3b_conv5y",
        "bl_metrics":   ["Accuracy", "BalancedAcc", "AUC-ROC", "F1"],
        "enc_metrics":  ["accuracy", "balanced_acc", "auc_roc", "f1"],
        "headers":      ["Acc", "Bal.Acc", "AUC", "F1"],
    },
    {
        "title":        "AD conversion \u2264 10 yrs",
        "subtitle":     f"N(test)={t3c_n}, non-conv={t3c_nc}, AD\u226410yr={t3c_conv}",
        "bl_task":      "10 years",
        "enc_task":     "T3c_conv10y",
        "bl_metrics":   ["Accuracy", "BalancedAcc", "AUC-ROC", "F1"],
        "enc_metrics":  ["accuracy", "balanced_acc", "auc_roc", "f1"],
        "headers":      ["Acc", "Bal.Acc", "AUC", "F1"],
    },
]

# ── Model rows ─────────────────────────────────────────────────────────────────
BL_MODELS = [
    ("LogReg",  "Log. Reg."),
    ("SVM",     "SVM"),
    ("XGBoost", "XGBoost"),
]

ENC_MODELS = [
    ("ModernBERT-base",                  "ModernBERT-B"),
    ("ModernBERT-large",                 "ModernBERT-L"),
    ("BioClinical-ModernBERT-base",      "BioClin-MBERT-B"),
    ("BioClinical-ModernBERT-large",     "BioClin-MBERT-L"),
]

ROW_ENTRIES = []

def bl_cell(model_id, grp):
    sub = bl_df[bl_df["Task"].str.contains(grp["bl_task"], na=False)
                & (bl_df["Model"] == model_id)]
    if sub.empty:
        return ["—"] * len(grp["bl_metrics"])
    cells = []
    for m in grp["bl_metrics"]:
        if m not in sub.columns:
            cells.append("—")
            continue
        val = str(sub.iloc[0][m]).strip()
        cells.append(val if val not in ("", "nan", "NaN") else "—")
    return cells

def enc_cell(model_id, strategy, grp):
    key = (model_id, grp["enc_task"], strategy)
    entry = enc_summary.get(key, {})
    cells = []
    for m in grp["enc_metrics"]:
        val = entry.get(m, "—")
        cells.append(val)
    return cells

# Section 1: Tabular baselines
for model_id, label in BL_MODELS:
    ROW_ENTRIES.append(("tabular", label, model_id))

# Section 2: Encoder frozen
for model_id, label in ENC_MODELS:
    ROW_ENTRIES.append(("frozen", label, model_id))

# Section 3: Encoder full fine-tune
for model_id, label in ENC_MODELS:
    ROW_ENTRIES.append(("full_ft", label, model_id))

n_rows = len(ROW_ENTRIES)

# ── Pre-compute all cells and find best per column ─────────────────────────────
def parse_mean(s):
    s = str(s).strip()
    if s in ("", "—", "nan", "NaN"):
        return float("-inf")
    try:
        return float(s.split("±")[0].strip())
    except ValueError:
        return float("-inf")

cell_grid = []
for strategy, label, model_id in ROW_ENTRIES:
    row_cells = []
    for grp in GROUPS:
        if strategy == "tabular":
            row_cells.extend(bl_cell(model_id, grp))
        else:
            row_cells.extend(enc_cell(model_id, strategy, grp))
    cell_grid.append(row_cells)

# Find the best (highest mean) for each column
n_total_cols = sum(len(g["headers"]) for g in GROUPS)
best_mask = [[False] * n_total_cols for _ in range(n_rows)]
for col_idx in range(n_total_cols):
    means = [parse_mean(cell_grid[r][col_idx]) for r in range(n_rows)]
    best_val = max(means)
    if best_val > float("-inf"):
        for r in range(n_rows):
            if means[r] == best_val:
                best_mask[r][col_idx] = True

# ── Sizing ─────────────────────────────────────────────────────────────────────
COL_W       = 1.45
ROW_H       = 0.32
ROW_H_SUB   = 0.22
MODEL_COL_W = 1.45
STRAT_COL_W = 0.90

n_data_cols = sum(len(g["headers"]) for g in GROUPS)
n_cols      = 2 + n_data_cols

fig_w = 0.3 + MODEL_COL_W + STRAT_COL_W + n_data_cols * COL_W + 0.3
fig_h = 0.15 + ROW_H * (3 + n_rows) + ROW_H_SUB + 0.10

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.axis("off")
ax.set_xlim(0, fig_w)
ax.set_ylim(0, fig_h)

# ── Column geometry ────────────────────────────────────────────────────────────
LEFT = 0.3
col_widths = [MODEL_COL_W, STRAT_COL_W] + [COL_W] * n_data_cols
col_left = [LEFT]
for w in col_widths:
    col_left.append(col_left[-1] + w)
col_cx = [(col_left[i] + col_left[i + 1]) / 2 for i in range(len(col_widths))]
RIGHT = col_left[-1]

# ── Row geometry ───────────────────────────────────────────────────────────────
TOP = fig_h - 0.15
y_cur = TOP
row_tops = []
row_tops.append(y_cur); y_cur -= ROW_H
row_tops.append(y_cur); y_cur -= ROW_H
row_tops.append(y_cur); y_cur -= ROW_H_SUB
row_tops.append(y_cur); y_cur -= ROW_H
for _ in range(n_rows):
    row_tops.append(y_cur); y_cur -= ROW_H
row_tops.append(y_cur)

header_rows = 4
total_visual_rows = header_rows + n_rows
def row_mid(i): return (row_tops[i] + row_tops[i + 1]) / 2

# ── Horizontal rules ───────────────────────────────────────────────────────────
def hline(y, lw=1.0, ls="-"):
    ax.plot([LEFT, RIGHT], [y, y], color="black", linewidth=lw, linestyle=ls,
            solid_capstyle="butt", zorder=3)

hline(row_tops[0], lw=1.5)
hline(row_tops[1], lw=1.5)
hline(row_tops[2], lw=0.8)
hline(row_tops[3], lw=0.5, ls=(0, (3, 3)))
hline(row_tops[4], lw=1.2)
hline(row_tops[4 + len(BL_MODELS)], lw=0.6, ls="--")
hline(row_tops[4 + len(BL_MODELS) + len(ENC_MODELS)], lw=0.6, ls="--")
hline(row_tops[total_visual_rows], lw=1.5)

# ── Dotted vertical lines between column groups ────────────────────────────────
col_cursor = 2
for grp in GROUPS:
    x = col_left[col_cursor]
    ax.plot([x, x], [row_tops[total_visual_rows], row_tops[1]],
            color="black", linewidth=0.6, linestyle=(0, (3, 3)), zorder=2)
    col_cursor += len(grp["headers"])

x_strat = col_left[1]
ax.plot([x_strat, x_strat], [row_tops[total_visual_rows], row_tops[3]],
        color="black", linewidth=0.5, linestyle=(0, (3, 3)), zorder=2)

# ── Title row (row 0) ─────────────────────────────────────────────────────────
ax.text((LEFT + RIGHT) / 2, row_mid(0),
        "Clinical Models (No CDR): Test Performance at Baseline   "
        "(mean \u00b1 std, seeds 0/1/2, 80/10/10 stratified split)",
        ha="center", va="center", fontsize=8.5, fontweight="bold", color="black")

# ── Group header row (row 1) ──────────────────────────────────────────────────
col_cursor = 2
for grp in GROUPS:
    n = len(grp["headers"])
    x0 = col_left[col_cursor]
    x1 = col_left[col_cursor + n]
    xm = (x0 + x1) / 2
    ax.text(xm, row_mid(1), grp["title"],
            ha="center", va="center", fontsize=8, fontweight="bold", color="black")
    pad = 0.05
    ax.plot([x0 + pad, x1 - pad], [row_tops[2], row_tops[2]],
            color="black", linewidth=0.7, zorder=4)
    col_cursor += n

# ── Subtitle row (row 2) ─────────────────────────────────────────────────────
col_cursor = 2
for grp in GROUPS:
    n = len(grp["headers"])
    x0 = col_left[col_cursor]
    x1 = col_left[col_cursor + n]
    xm = (x0 + x1) / 2
    ax.text(xm, row_mid(2), grp["subtitle"],
            ha="center", va="center", fontsize=6.5,
            fontstyle="italic", color="#555555")
    col_cursor += n

# ── Sub-header row (row 3) ────────────────────────────────────────────────────
ax.text(col_cx[0], row_mid(3), "", ha="center", va="center")
ax.text(col_cx[1], row_mid(3), "Strategy",
        ha="center", va="center", fontsize=7.5, fontstyle="italic", color="black")
col_cursor = 2
for grp in GROUPS:
    for h in grp["headers"]:
        ax.text(col_cx[col_cursor], row_mid(3), h,
                ha="center", va="center", fontsize=7.5,
                fontstyle="italic", color="black")
        col_cursor += 1

# ── Data rows ──────────────────────────────────────────────────────────────────
STRATEGY_LABELS = {"tabular": "tabular", "frozen": "frozen", "full_ft": "full fine-tune"}

for r_idx, (strategy, label, model_id) in enumerate(ROW_ENTRIES):
    row_i = 4 + r_idx
    ax.text(col_left[0] + 0.06, row_mid(row_i), label,
            ha="left", va="center", fontsize=7.5, fontweight="bold", color="black")
    ax.text(col_cx[1], row_mid(row_i), STRATEGY_LABELS[strategy],
            ha="center", va="center", fontsize=7, color="black")
    col_cursor = 2
    for ci, val in enumerate(cell_grid[r_idx]):
        fw = "bold" if best_mask[r_idx][ci] else "normal"
        ax.text(col_cx[col_cursor], row_mid(row_i), val,
                ha="center", va="center", fontsize=7.3, color="black",
                fontweight=fw)
        col_cursor += 1

# ── Outer box ──────────────────────────────────────────────────────────────────
rect = plt.Rectangle((LEFT, row_tops[total_visual_rows]),
                      RIGHT - LEFT, row_tops[0] - row_tops[total_visual_rows],
                      facecolor="none", edgecolor="black", linewidth=1.5, zorder=5)
ax.add_patch(rect)

# ── Save ───────────────────────────────────────────────────────────────────────
plt.tight_layout(pad=0.1)
OUT_DIR.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_DIR / "combined_model_table.pdf", bbox_inches="tight", dpi=300)
fig.savefig(OUT_DIR / "combined_model_table.png", bbox_inches="tight", dpi=300)
print(f"Saved PDF → {OUT_DIR / 'combined_model_table.pdf'}")
print(f"Saved PNG → {OUT_DIR / 'combined_model_table.png'}")
plt.close()

# ── LaTeX booktabs ─────────────────────────────────────────────────────────────
def latex_val(s):
    return str(s).replace("±", r"$\pm$").replace("—", r"---")

col_spec = "ll" + "".join(
    r"@{\hspace{3pt}}|@{\hspace{3pt}}" + "c" * len(g["headers"])
    for g in GROUPS
)
lines = [
    r"\begin{table}[ht]",
    r"\centering",
    r"\caption{Clinical Models (No CDR): Test Performance at Baseline (mean $\pm$ std, 3 seeds $\times$ 80/10/10 stratified split)}",
    r"\label{tab:combined_models_no_cdr}",
    r"\small\setlength{\tabcolsep}{4pt}",
    r"\begin{tabular}{" + col_spec + "}",
    r"\toprule",
]
grp_row = " & & " + " & ".join(
    r"\multicolumn{" + str(len(g["headers"])) + r"}{c}{\textbf{" + g["title"] + r"}}"
    for g in GROUPS
) + r" \\"
lines.append(grp_row)
cursor = 3
cmidrules = ""
for g in GROUPS:
    n = len(g["headers"])
    cmidrules += r"\cmidrule(lr){" + f"{cursor}-{cursor+n-1}" + "}"
    cursor += n
lines.append(cmidrules)
lines.append("Model & Strategy & " + " & ".join(
    r"\textit{" + h + r"}" for g in GROUPS for h in g["headers"]
) + r" \\")
lines.append(r"\midrule")

prev_strategy = None
for r_idx, (strategy, label, model_id) in enumerate(ROW_ENTRIES):
    if prev_strategy is not None and strategy != prev_strategy:
        lines.append(r"\midrule")
    prev_strategy = strategy
    parts = [r"\textbf{" + label + "}", STRATEGY_LABELS[strategy]]
    for ci, val in enumerate(cell_grid[r_idx]):
        lv = latex_val(val)
        if best_mask[r_idx][ci]:
            lv = r"\textbf{" + lv + "}"
        parts.append(lv)
    lines.append(" & ".join(parts) + r" \\")

lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
tex_path = OUT_DIR / "combined_model_table.tex"
tex_path.write_text("\n".join(lines), encoding="utf-8")
print(f"Saved LaTeX → {tex_path}")
