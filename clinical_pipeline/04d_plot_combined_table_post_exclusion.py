"""
04d_plot_combined_table_post_exclusion.py
=========================================
Combined results table (tabular baselines + encoder LLMs) for the NO-CDR
stratified *POST-EXCLUSION* cohort, expanded to 10 task columns:

  Diagnosis : T1 CN vs MCI+AD | T1b CN+MCI vs AD | T2 CN/MCI/AD (3-class)
  Stratified: T1c sCN vs progressors+AD | T1d pMCI vs sMCI | T1e sCN vs pCN
  Conversion: AD<=3y | AD<=5y | AD<=7y | AD<=10y

Reads:
  - tabular baselines : outputs/encoder_only_post_exclusions/baseline_model_comparison.csv
                        (produced by 02h_baseline_post_exclusion.py)
  - encoder metrics   : outputs/encoder_only_post_exclusions/**/metrics.json
                        (produced by 03_encoder_finetune.py via 03d submit)
Encoder cells render "—" until the HPC runs land; re-run to fill.

Per-model hyperparameters are noted with superscripts a/b/c/d on the model names and
footnotes below the table.

Usage:
  python clinical_pipeline/04d_plot_combined_table_post_exclusion.py
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
OUT_DIR   = BASE / "outputs" / "encoder_only_post_exclusions"
BL_CSV    = OUT_DIR / "baseline_model_comparison.csv"
ENC_DIR   = OUT_DIR
SPLIT_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical"
                 r"\no_cdr_stratified_post_exclusion\tabular\baseline")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── 1. Load tabular baselines (may not exist yet) ─────────────────────────────
if BL_CSV.exists():
    bl_df = pd.read_csv(BL_CSV)
else:
    print(f"  [WARN] {BL_CSV} not found — baseline cells will be '—'. "
          f"Run 02h_baseline_post_exclusion.py first.")
    bl_df = pd.DataFrame(columns=["Task", "Model"])

# ── 1b. Test-set class distributions (seed 0, post-exclusion) ─────────────────
_s0 = pd.read_csv(SPLIT_DIR / "seed_0" / "test.csv", low_memory=False)

def _bl_counts(mapping):
    y = _s0["Label_bl_multi"].map(mapping).dropna()
    return len(y), int((y == 0).sum()), int((y == 1).sum())

def _cg_counts(mapping):
    y = _s0["conversion_group"].map(mapping).dropna()
    return len(y), int((y == 0).sum()), int((y == 1).sum())

def _conv_counts(col):
    sub = _s0[_s0["Label_bl_multi"].isin(["CN", "MCI"])]
    y = pd.to_numeric(sub[col], errors="coerce").dropna()
    return len(y), int((y == 0).sum()), int(y.sum())

t1_n,  t1_cn,  t1_pos  = _bl_counts({"CN": 0, "MCI": 1, "AD": 1})
t1b_n, t1b_neg, t1b_ad = _bl_counts({"CN": 0, "MCI": 0, "AD": 1})
t2_n  = int(_s0["Label_bl_multi"].isin(["CN", "MCI", "AD"]).sum())
t2_cn = int((_s0["Label_bl_multi"] == "CN").sum())
t2_mci = int((_s0["Label_bl_multi"] == "MCI").sum())
t2_ad = int((_s0["Label_bl_multi"] == "AD").sum())
t1c_n, t1c_scn, t1c_pos = _cg_counts({"sCN": 0, "pCN_to_MCI": 1, "pCN_to_AD": 1,
                                      "pMCI": 1, "AD_bl": 1})
t1d_n, t1d_smci, t1d_pmci = _cg_counts({"sMCI": 0, "pMCI": 1})
t1e_n, t1e_scn, t1e_pcn   = _cg_counts({"sCN": 0, "pCN_to_MCI": 1, "pCN_to_AD": 1})
t3a_n, t3a_nc, t3a_conv = _conv_counts("Label_3y")
t3b_n, t3b_nc, t3b_conv = _conv_counts("Label_5y")
t3d_n, t3d_nc, t3d_conv = _conv_counts("Label_7y")
t3c_n, t3c_nc, t3c_conv = _conv_counts("Label_10y")

# ── 2. Load encoder metrics and aggregate across seeds ─────────────────────────
# Non-scalar / non-metric keys excluded from aggregation.
EXCLUDE = {"loss", "runtime", "samples_per_second", "steps_per_second", "epoch",
           "confusion_matrix", "per_class_recall"}

enc_rows = []
n_skipped = 0
for jf in ENC_DIR.rglob("metrics.json"):
    try:
        with open(jf) as f:
            data = json.load(f)
    except (json.JSONDecodeError, ValueError):
        n_skipped += 1
        continue
    cfg = data["config"]
    met = {k: v for k, v in data["test_metrics"].items() if k not in EXCLUDE}
    enc_rows.append({
        "model":    cfg["model_id"].split("/")[-1],
        "task":     cfg["task"],
        "strategy": cfg["strategy"],
        "seed":     cfg["seed"],
        **met,
    })

enc_summary = {}
if enc_rows:
    enc_df = pd.DataFrame(enc_rows)
    print(f"  Found {len(enc_rows)} encoder metric files")
    print(f"  Tasks: {sorted(enc_df['task'].unique())}")
    meta = ["model", "task", "strategy", "seed"]
    metric_cols = [c for c in enc_df.columns if c not in meta and c not in EXCLUDE]
    for (model, task, strategy), grp in enc_df.groupby(["model", "task", "strategy"]):
        entry = {}
        for m in metric_cols:
            vals = pd.to_numeric(grp[m], errors="coerce").dropna()
            if len(vals) > 0:
                sd = vals.std() if len(vals) > 1 else 0.0
                entry[m] = f"{vals.mean():.4f} ± {sd:.4f}"
            else:
                entry[m] = "—"
        enc_summary[(model, task, strategy)] = entry
    enc_df.to_csv(OUT_DIR / "encoder_results_per_run.csv", index=False)
else:
    print(f"  [WARN] No encoder metrics.json under {ENC_DIR} yet — encoder cells will be '—'.")

# ── Table layout: 10 task groups ───────────────────────────────────────────────
# Accuracy omitted — uninformative on the imbalanced cohort; report balanced accuracy.
BINARY_METRICS = (["BalancedAcc", "AUC-ROC", "F1"],
                  ["balanced_acc", "auc_roc", "f1"],
                  ["Bal.Acc", "AUC", "F1"])
MULTI_METRICS  = (["BalancedAcc", "AUC-ROC (OvR)", "MacroF1"],
                  ["balanced_acc", "auc_roc_ovr", "macro_f1"],
                  ["Bal.Acc", "AUC", "MacroF1"])

def grp_binary(title, subtitle, bl_task, enc_task):
    return {"title": title, "subtitle": subtitle, "bl_task": bl_task,
            "enc_task": enc_task, "bl_metrics": BINARY_METRICS[0],
            "enc_metrics": BINARY_METRICS[1], "headers": BINARY_METRICS[2]}

GROUPS = [
    grp_binary("CN vs MCI+AD",  f"N={t1_n}, CN={t1_cn}, MCI+AD={t1_pos}",
               "Binary: CN vs MCI+AD", "T1_binary"),
    grp_binary("CN+MCI vs AD",  f"N={t1b_n}, CN+MCI={t1b_neg}, AD={t1b_ad}",
               "CN+MCI vs AD", "T1b_cnmci_ad"),
    {"title": "CN / MCI / AD (3-class)",
     "subtitle": f"N={t2_n}, CN={t2_cn}, MCI={t2_mci}, AD={t2_ad}",
     "bl_task": "Multi-class: CN / MCI / AD", "enc_task": "T2_multiclass",
     "bl_metrics": MULTI_METRICS[0], "enc_metrics": MULTI_METRICS[1],
     "headers": MULTI_METRICS[2]},
    grp_binary("sCN vs prog.+AD", f"N={t1c_n}, sCN={t1c_scn}, prog+AD={t1c_pos}",
               "sCN vs progressors+AD (excl. sMCI)", "T1c_scn_prog"),
    grp_binary("pMCI vs sMCI",  f"N={t1d_n}, sMCI={t1d_smci}, pMCI={t1d_pmci}",
               "pMCI vs sMCI (baseline MCI)", "T1d_pmci_smci"),
    grp_binary("sCN vs pCN",    f"N={t1e_n}, sCN={t1e_scn}, pCN={t1e_pcn}",
               "sCN vs pCN (baseline CN, to MCI/AD)", "T1e_scn_pcn"),
    grp_binary("AD conversion ≤ 3 yrs", f"N={t3a_n}, non-conv={t3a_nc}, AD≤3yr={t3a_conv}",
               "Conversion to AD within 3 years", "T3a_conv3y"),
    grp_binary("AD conversion ≤ 5 yrs", f"N={t3b_n}, non-conv={t3b_nc}, AD≤5yr={t3b_conv}",
               "Conversion to AD within 5 years", "T3b_conv5y"),
    grp_binary("AD conversion ≤ 7 yrs", f"N={t3d_n}, non-conv={t3d_nc}, AD≤7yr={t3d_conv}",
               "Conversion to AD within 7 years", "T3d_conv7y"),
    grp_binary("AD conversion ≤ 10 yrs", f"N={t3c_n}, non-conv={t3c_nc}, AD≤10yr={t3c_conv}",
               "Conversion to AD within 10 years", "T3c_conv10y"),
]

# ── Model rows + HP superscripts ───────────────────────────────────────────────
BL_MODELS  = [("LogReg", "Log. Reg."), ("SVM", "SVM"), ("XGBoost", "XGBoost")]
ENC_MODELS = [
    ("ModernBERT-base",              "ModernBERT-B"),
    ("ModernBERT-large",             "ModernBERT-L"),
    ("BioClinical-ModernBERT-base",  "BioClin-MBERT-B"),
    ("BioClinical-ModernBERT-large", "BioClin-MBERT-L"),
]
SUPERSCRIPT = {"Log. Reg.": "a", "SVM": "b", "XGBoost": "c",
               "ModernBERT-B": "d", "ModernBERT-L": "d",
               "BioClin-MBERT-B": "d", "BioClin-MBERT-L": "d"}

FOOTNOTES = [
    r"$^{a}$ Logistic Regression: L2 penalty, C=1.0, max_iter=2000, class_weight=balanced.",
    r"$^{b}$ SVM: RBF kernel, C=1.0, gamma=scale, probability=True, class_weight=balanced.",
    r"$^{c}$ XGBoost: 300 trees, max_depth=4, learning_rate=0.05, eval_metric=logloss.   "
    r"All tabular models: median imputation + standardisation, seed=42.",
    r"$^{d}$ Encoders (ModernBERT / BioClinical-ModernBERT, base & large): AdamW, batch=16, "
    r"max_len=1024, weight_decay=1e-5, warmup_ratio=0.1, class-weighted cross-entropy, "
    r"early-stopping patience=5, fp16.",
    r"      Frozen: lr=1e-3, 20 epochs (classifier head only).   "
    r"Full fine-tune: lr=2e-5, 10 epochs.   Seeds 0/1/2; 80/10/10 stratified split.",
]

ROW_ENTRIES = []
for model_id, label in BL_MODELS:
    ROW_ENTRIES.append(("tabular", label, model_id))
for model_id, label in ENC_MODELS:
    ROW_ENTRIES.append(("frozen", label, model_id))
for model_id, label in ENC_MODELS:
    ROW_ENTRIES.append(("full_ft", label, model_id))
n_rows = len(ROW_ENTRIES)


def bl_cell(model_id, grp):
    sub = bl_df[(bl_df["Task"] == grp["bl_task"]) & (bl_df["Model"] == model_id)]
    if sub.empty:
        return ["—"] * len(grp["bl_metrics"])
    cells = []
    for m in grp["bl_metrics"]:
        if m not in sub.columns:
            cells.append("—"); continue
        val = str(sub.iloc[0][m]).strip()
        cells.append(val if val not in ("", "nan", "NaN") else "—")
    return cells


def enc_cell(model_id, strategy, grp):
    entry = enc_summary.get((model_id, grp["enc_task"], strategy), {})
    return [entry.get(m, "—") for m in grp["enc_metrics"]]


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
        row_cells.extend(bl_cell(model_id, grp) if strategy == "tabular"
                         else enc_cell(model_id, strategy, grp))
    cell_grid.append(row_cells)

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
COL_W, ROW_H, ROW_H_SUB, MODEL_COL_W, STRAT_COL_W = 1.40, 0.32, 0.22, 1.55, 0.95
FOOTER_H = 0.16 * len(FOOTNOTES) + 0.25

n_data_cols = sum(len(g["headers"]) for g in GROUPS)
fig_w = 0.3 + MODEL_COL_W + STRAT_COL_W + n_data_cols * COL_W + 0.3
fig_h = 0.15 + ROW_H * (3 + n_rows) + ROW_H_SUB + 0.10 + FOOTER_H

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.axis("off"); ax.set_xlim(0, fig_w); ax.set_ylim(0, fig_h)

LEFT = 0.3
col_widths = [MODEL_COL_W, STRAT_COL_W] + [COL_W] * n_data_cols
col_left = [LEFT]
for w in col_widths:
    col_left.append(col_left[-1] + w)
col_cx = [(col_left[i] + col_left[i + 1]) / 2 for i in range(len(col_widths))]
RIGHT = col_left[-1]

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

total_visual_rows = 4 + n_rows
def row_mid(i): return (row_tops[i] + row_tops[i + 1]) / 2

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

# Dotted vertical separators between task groups
col_cursor = 2
for grp in GROUPS:
    x = col_left[col_cursor]
    ax.plot([x, x], [row_tops[total_visual_rows], row_tops[1]],
            color="black", linewidth=0.6, linestyle=(0, (3, 3)), zorder=2)
    col_cursor += len(grp["headers"])
x_strat = col_left[1]
ax.plot([x_strat, x_strat], [row_tops[total_visual_rows], row_tops[3]],
        color="black", linewidth=0.5, linestyle=(0, (3, 3)), zorder=2)

# Title
ax.text((LEFT + RIGHT) / 2, row_mid(0),
        "Clinical Models (No CDR, post-exclusion): Test Performance at Baseline   "
        "(mean ± std, seeds 0/1/2, 80/10/10 stratified split)",
        ha="center", va="center", fontsize=9, fontweight="bold", color="black")

# Group headers
col_cursor = 2
for grp in GROUPS:
    n = len(grp["headers"])
    x0, x1 = col_left[col_cursor], col_left[col_cursor + n]
    ax.text((x0 + x1) / 2, row_mid(1), grp["title"],
            ha="center", va="center", fontsize=8, fontweight="bold", color="black")
    ax.plot([x0 + 0.05, x1 - 0.05], [row_tops[2], row_tops[2]],
            color="black", linewidth=0.7, zorder=4)
    col_cursor += n

# Subtitles
col_cursor = 2
for grp in GROUPS:
    n = len(grp["headers"])
    x0, x1 = col_left[col_cursor], col_left[col_cursor + n]
    ax.text((x0 + x1) / 2, row_mid(2), grp["subtitle"],
            ha="center", va="center", fontsize=6.0, fontstyle="italic", color="#555555")
    col_cursor += n

# Sub-header row (metric names)
ax.text(col_cx[1], row_mid(3), "Strategy",
        ha="center", va="center", fontsize=7.5, fontstyle="italic", color="black")
col_cursor = 2
for grp in GROUPS:
    for h in grp["headers"]:
        ax.text(col_cx[col_cursor], row_mid(3), h, ha="center", va="center",
                fontsize=7.3, fontstyle="italic", color="black")
        col_cursor += 1

# Data rows
STRATEGY_LABELS = {"tabular": "tabular", "frozen": "frozen", "full_ft": "full fine-tune"}
for r_idx, (strategy, label, model_id) in enumerate(ROW_ENTRIES):
    row_i = 4 + r_idx
    disp = f"{label}$^{{{SUPERSCRIPT[label]}}}$"
    ax.text(col_left[0] + 0.06, row_mid(row_i), disp,
            ha="left", va="center", fontsize=7.5, fontweight="bold", color="black")
    ax.text(col_cx[1], row_mid(row_i), STRATEGY_LABELS[strategy],
            ha="center", va="center", fontsize=7, color="black")
    col_cursor = 2
    for ci, val in enumerate(cell_grid[r_idx]):
        fw = "bold" if best_mask[r_idx][ci] else "normal"
        ax.text(col_cx[col_cursor], row_mid(row_i), val,
                ha="center", va="center", fontsize=7.1, color="black", fontweight=fw)
        col_cursor += 1

# Outer box
rect = plt.Rectangle((LEFT, row_tops[total_visual_rows]),
                     RIGHT - LEFT, row_tops[0] - row_tops[total_visual_rows],
                     facecolor="none", edgecolor="black", linewidth=1.5, zorder=5)
ax.add_patch(rect)

# Footnotes (hyperparameters) below the table
y_fn = row_tops[total_visual_rows] - 0.22
for line in FOOTNOTES:
    ax.text(LEFT, y_fn, line, ha="left", va="top", fontsize=6.3, color="#222222")
    y_fn -= 0.16

plt.tight_layout(pad=0.1)
fig.savefig(OUT_DIR / "combined_model_table.pdf", bbox_inches="tight", dpi=300)
fig.savefig(OUT_DIR / "combined_model_table.png", bbox_inches="tight", dpi=300)
print(f"Saved PDF → {OUT_DIR / 'combined_model_table.pdf'}")
print(f"Saved PNG → {OUT_DIR / 'combined_model_table.png'}")
plt.close()
