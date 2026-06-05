"""
04d_plot_combined_table_post_exclusion.py
=========================================
Combined results table (tabular baselines + encoder LLMs) for the NO-CDR
stratified *POST-EXCLUSION* cohort, expanded to 10 task columns:

  Diagnosis : T1 CN vs MCI+AD | T1b CN+MCI vs AD | T2 CN/MCI/AD (3-class)
  Stratified: T1c sCN vs progressors+AD | T1d pMCI vs sMCI | T1e sCN vs pCN
  Conversion: AD<=3y | AD<=5y | AD<=7y | AD<=10y

VAL and TEST are rendered as SEPARATE tables (both held-out):
  combined_model_table_test.png   (also copied to combined_model_table.png)
  combined_model_table_val.png

Reads:
  - tabular baselines : outputs/encoder_only_post_exclusions/baseline_model_comparison.csv
                        (02h; now carries a `Split` column with val/test rows)
  - encoder metrics   : outputs/encoder_only_post_exclusions/**/metrics.json
                        (03_encoder_finetune.py; stores `val_metrics` + `test_metrics`)
Encoder cells render "—" until the HPC runs land; re-run to fill.

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
# Main-sweep encoder metrics were scp'd here (separate from the repo outputs + the baseline CSV).
ENC_DIR   = Path(r"D:\ADNI_BIDS_project\derivatives\encoder_outputs_no_cdr_post_exclusion")
SPLIT_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical"
                 r"\no_cdr_stratified_post_exclusion\tabular\baseline")
# T4 is an isolated job: its own metrics dir + its own (converter-cohort) split tree.
# (scp'd locally to D:\...\derivatives\ — separate from the main-sweep metrics.)
ENC_DIR_T4   = Path(r"D:\ADNI_BIDS_project\derivatives\encoder_outputs_no_cdr_post_exclusion_T4")
T4_SPLIT_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical"
                    r"\no_cdr_stratified_post_exclusion\tabular\baseline_T4")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── 1. Load tabular baselines (may not exist yet) ─────────────────────────────
if BL_CSV.exists():
    bl_df = pd.read_csv(BL_CSV)
else:
    print(f"  [WARN] {BL_CSV} not found — baseline cells will be '—'. "
          f"Run 02h_baseline_post_exclusion.py first.")
    bl_df = pd.DataFrame(columns=["Task", "Model"])
# Older CSVs (pre val/test refactor) have no Split column → treat everything as test.
if "Split" not in bl_df.columns:
    bl_df["Split"] = "test"
# T4 tabular baseline lives in a separate CSV (its own converter cohort/tree) — append it
# so the T4 column gets LogReg/SVM/XGBoost rows (02m_baseline_T4.py).
BL_CSV_T4 = BASE / "outputs" / "baseline_no_cdr_post_exclusion" / "baseline_model_comparison_T4.csv"
if BL_CSV_T4.exists():
    _t4bl = pd.read_csv(BL_CSV_T4)
    if "Split" not in _t4bl.columns:
        _t4bl["Split"] = "test"
    bl_df = pd.concat([bl_df, _t4bl], ignore_index=True)

# ── 2. Encoder metrics: load BOTH val and test, aggregate across seeds per split ─
EXCLUDE = {"loss", "runtime", "samples_per_second", "steps_per_second", "epoch",
           "confusion_matrix", "per_class_recall"}

enc_rows = {"val": [], "test": []}
n_skipped = 0
_jfiles = list(ENC_DIR.rglob("metrics.json"))
if ENC_DIR_T4.exists():                       # fold in the isolated T4 job's runs
    _jfiles += list(ENC_DIR_T4.rglob("metrics.json"))
for jf in _jfiles:
    try:
        with open(jf) as f:
            data = json.load(f)
    except (json.JSONDecodeError, ValueError):
        n_skipped += 1
        continue
    cfg = data.get("config", {})
    for split in ("val", "test"):
        block = data.get(f"{split}_metrics")
        if not block:
            continue
        met = {k: v for k, v in block.items() if k not in EXCLUDE}
        enc_rows[split].append({
            "model": cfg.get("model_id", "?").split("/")[-1], "task": cfg.get("task"),
            "strategy": cfg.get("strategy"), "seed": cfg.get("seed"), **met})


def _build_enc_summary(rows):
    summary = {}
    if not rows:
        return summary
    edf = pd.DataFrame(rows)
    meta = ["model", "task", "strategy", "seed"]
    mcols = [c for c in edf.columns if c not in meta and c not in EXCLUDE]
    for (model, task, strategy), grp in edf.groupby(["model", "task", "strategy"]):
        entry = {}
        for m in mcols:
            vals = pd.to_numeric(grp[m], errors="coerce").dropna()
            entry[m] = (f"{vals.mean():.4f} ± {(vals.std() if len(vals) > 1 else 0.0):.4f}"
                        if len(vals) > 0 else "—")
        summary[(model, task, strategy)] = entry
    return summary


ENC_SUMMARY = {s: _build_enc_summary(enc_rows[s]) for s in ("val", "test")}
print(f"  Encoder metric files: val={len(enc_rows['val'])}, test={len(enc_rows['test'])}"
      + (f"  ({n_skipped} skipped)" if n_skipped else ""))

# ── Metric column maps ─────────────────────────────────────────────────────────
# Accuracy omitted — uninformative on the imbalanced cohort; report balanced accuracy.
BINARY_METRICS = (["BalancedAcc", "AUC-ROC", "F1"],
                  ["balanced_acc", "auc_roc", "f1"],
                  ["Bal.Acc", "AUC", "F1"])
MULTI_METRICS  = (["BalancedAcc", "AUC-ROC (OvR)", "MacroF1"],
                  ["balanced_acc", "auc_roc_ovr", "macro_f1"],
                  ["Bal.Acc", "AUC", "MacroF1"])


def grp_binary(title, bl_task, enc_task):
    return {"title": title, "bl_task": bl_task, "enc_task": enc_task,
            "bl_metrics": BINARY_METRICS[0], "enc_metrics": BINARY_METRICS[1],
            "headers": BINARY_METRICS[2]}


GROUPS_BASE = [
    grp_binary("CN vs MCI+AD", "Binary: CN vs MCI+AD", "T1_binary"),
    grp_binary("CN+MCI vs AD", "CN+MCI vs AD", "T1b_cnmci_ad"),
    {"title": "CN / MCI / AD (3-class)", "bl_task": "Multi-class: CN / MCI / AD",
     "enc_task": "T2_multiclass", "bl_metrics": MULTI_METRICS[0],
     "enc_metrics": MULTI_METRICS[1], "headers": MULTI_METRICS[2]},
    grp_binary("sCN vs prog.+AD", "sCN vs progressors+AD (excl. sMCI)", "T1c_scn_prog"),
    grp_binary("pMCI vs sMCI", "pMCI vs sMCI (baseline MCI)", "T1d_pmci_smci"),
    grp_binary("sCN vs pCN", "sCN vs pCN (baseline CN, to MCI/AD)", "T1e_scn_pcn"),
    grp_binary("AD conversion ≤ 3 yrs", "Conversion to AD within 3 years", "T3a_conv3y"),
    grp_binary("AD conversion ≤ 5 yrs", "Conversion to AD within 5 years", "T3b_conv5y"),
    grp_binary("AD conversion ≤ 7 yrs", "Conversion to AD within 7 years", "T3d_conv7y"),
    grp_binary("AD conversion ≤ 10 yrs", "Conversion to AD within 10 years", "T3c_conv10y"),
    # T4 — isolated converter-cohort task (3-class horizon); tabular baseline via 02m.
    {"title": "AD horizon (3-class)", "bl_task": "AD horizon (3-class): pMCI + pCN_to_AD",
     "enc_task": "T4_conv_horizon",
     "bl_metrics": MULTI_METRICS[0], "enc_metrics": MULTI_METRICS[1], "headers": MULTI_METRICS[2]},
]


def t4_subtitle(s0):
    """T4 cohort breakdown (3-class horizon over converters) from a baseline_T4 split."""
    if s0 is None or "Label_T4" not in s0.columns:
        return "converters (pMCI + pCN→AD)"
    y = pd.to_numeric(s0["Label_T4"], errors="coerce").dropna()
    return f"N={len(y)}, <3y={(y==0).sum()}, 3–7y={(y==1).sum()}, ≥7y={(y==2).sum()}"


def subtitle_for(grp, s0):
    """N + class breakdown for a group, computed from a split dataframe s0."""
    def bl_counts(mapping):
        y = s0["Label_bl_multi"].map(mapping).dropna()
        return len(y), int((y == 0).sum()), int((y == 1).sum())

    def cg_counts(mapping):
        y = s0["conversion_group"].map(mapping).dropna()
        return len(y), int((y == 0).sum()), int((y == 1).sum())

    def conv_counts(col):
        sub = s0[s0["Label_bl_multi"].isin(["CN", "MCI"])]
        y = pd.to_numeric(sub[col], errors="coerce").dropna()
        return len(y), int((y == 0).sum()), int(y.sum())

    t = grp["enc_task"]
    if t == "T1_binary":
        n, cn, pos = bl_counts({"CN": 0, "MCI": 1, "AD": 1}); return f"N={n}, CN={cn}, MCI+AD={pos}"
    if t == "T1b_cnmci_ad":
        n, neg, ad = bl_counts({"CN": 0, "MCI": 0, "AD": 1}); return f"N={n}, CN+MCI={neg}, AD={ad}"
    if t == "T2_multiclass":
        n = int(s0["Label_bl_multi"].isin(["CN", "MCI", "AD"]).sum())
        return (f"N={n}, CN={(s0['Label_bl_multi']=='CN').sum()}, "
                f"MCI={(s0['Label_bl_multi']=='MCI').sum()}, AD={(s0['Label_bl_multi']=='AD').sum()}")
    if t == "T1c_scn_prog":
        n, scn, pos = cg_counts({"sCN": 0, "pCN_to_MCI": 1, "pCN_to_AD": 1, "pMCI": 1, "AD_bl": 1})
        return f"N={n}, sCN={scn}, prog+AD={pos}"
    if t == "T1d_pmci_smci":
        n, smci, pmci = cg_counts({"sMCI": 0, "pMCI": 1}); return f"N={n}, sMCI={smci}, pMCI={pmci}"
    if t == "T1e_scn_pcn":
        n, scn, pcn = cg_counts({"sCN": 0, "pCN_to_MCI": 1, "pCN_to_AD": 1}); return f"N={n}, sCN={scn}, pCN={pcn}"
    conv_col = {"T3a_conv3y": "Label_3y", "T3b_conv5y": "Label_5y",
                "T3d_conv7y": "Label_7y", "T3c_conv10y": "Label_10y"}[t]
    yr = {"T3a_conv3y": 3, "T3b_conv5y": 5, "T3d_conv7y": 7, "T3c_conv10y": 10}[t]
    n, nc, conv = conv_counts(conv_col); return f"N={n}, non-conv={nc}, AD≤{yr}yr={conv}"


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
    r"All tabular models: median imputation + standardisation (fit on train), seed=42.",
    r"$^{d}$ Encoders (ModernBERT / BioClinical-ModernBERT, base & large): AdamW, batch=16, "
    r"max_len=1024, weight_decay=1e-5, warmup_ratio=0.1, class-weighted cross-entropy, "
    r"early-stopping patience=5, bf16.",
    r"      Frozen: lr=1e-3, 20 epochs (classifier head only).   "
    r"Full fine-tune: lr=2e-5, 10 epochs.   Seeds 0/1/2; 80/10/10 stratified split.",
    r"      AD horizon (T4): 3-class time-to-AD (<3y / 3–7y / ≥7y) over confirmed converters "
    r"(pMCI + pCN→AD), separate split tree; tabular baselines use the full battery joined to that cohort.",
]
VAL_FOOTNOTE = (r"      NB Validation is held-out for both, but the encoders use it for "
                r"early-stopping / checkpoint selection (mildly optimistic); the tabular "
                r"baselines do not. Test is untouched for both.")

ROW_ENTRIES = []
for model_id, label in BL_MODELS:
    ROW_ENTRIES.append(("tabular", label, model_id))
for model_id, label in ENC_MODELS:
    ROW_ENTRIES.append(("frozen", label, model_id))
for model_id, label in ENC_MODELS:
    ROW_ENTRIES.append(("full_ft", label, model_id))
n_rows = len(ROW_ENTRIES)


def parse_mean(s):
    s = str(s).strip()
    if s in ("", "—", "nan", "NaN"):
        return float("-inf")
    try:
        return float(s.split("±")[0].strip())
    except ValueError:
        return float("-inf")


def render_split(split):
    bl_split = bl_df[bl_df["Split"] == split]
    enc_summary = ENC_SUMMARY[split]
    s0 = pd.read_csv(SPLIT_DIR / "seed_0" / f"{split}.csv", low_memory=False)
    s0_t4 = None
    _t4f = T4_SPLIT_DIR / "seed_0" / f"{split}.csv"
    if _t4f.exists():
        s0_t4 = pd.read_csv(_t4f, low_memory=False)

    GROUPS = []
    for g in GROUPS_BASE:
        g2 = dict(g)
        if g["enc_task"] == "T4_conv_horizon":
            g2["subtitle"] = t4_subtitle(s0_t4)
        else:
            g2["subtitle"] = subtitle_for(g, s0)
        GROUPS.append(g2)

    def bl_cell(model_id, grp):
        sub = bl_split[(bl_split["Task"] == grp["bl_task"]) & (bl_split["Model"] == model_id)]
        if sub.empty:
            return ["—"] * len(grp["bl_metrics"])
        out = []
        for m in grp["bl_metrics"]:
            if m not in sub.columns:
                out.append("—"); continue
            v = str(sub.iloc[0][m]).strip()
            out.append(v if v not in ("", "nan", "NaN") else "—")
        return out

    def enc_cell(model_id, strategy, grp):
        entry = enc_summary.get((model_id, grp["enc_task"], strategy), {})
        return [entry.get(m, "—") for m in grp["enc_metrics"]]

    cell_grid = []
    for strategy, label, model_id in ROW_ENTRIES:
        row_cells = []
        for grp in GROUPS:
            row_cells.extend(bl_cell(model_id, grp) if strategy == "tabular"
                             else enc_cell(model_id, strategy, grp))
        cell_grid.append(row_cells)

    # ── CSV backup of the rendered table (mirrors the figure cells 1:1) ───────────
    _strat_lbl = {"tabular": "tabular", "frozen": "frozen", "full_ft": "full fine-tune"}
    _col_names = [f"{g['title']} | {h}" for g in GROUPS for h in g["headers"]]
    _records = []
    for r_idx, (strategy, label, model_id) in enumerate(ROW_ENTRIES):
        rec = {"Model": label, "Strategy": _strat_lbl[strategy]}
        rec.update(dict(zip(_col_names, cell_grid[r_idx])))
        _records.append(rec)
    pd.DataFrame(_records, columns=["Model", "Strategy"] + _col_names).to_csv(
        OUT_DIR / f"combined_model_table_{split}.csv", index=False)
    print(f"  Saved → {OUT_DIR / f'combined_model_table_{split}.csv'}")

    n_total_cols = sum(len(g["headers"]) for g in GROUPS)
    best_mask = [[False] * n_total_cols for _ in range(n_rows)]
    # Bold the best WITHIN each block per column: best tabular baseline AND best encoder,
    # so the baseline winner is highlighted even when an encoder wins overall (and vice-versa).
    baseline_rows = list(range(len(BL_MODELS)))          # 0..2  (LogReg/SVM/XGBoost)
    encoder_rows  = list(range(len(BL_MODELS), n_rows))  # the frozen + full_ft encoder rows
    for col_idx in range(n_total_cols):
        for block in (baseline_rows, encoder_rows):
            finite = [(r, parse_mean(cell_grid[r][col_idx])) for r in block]
            finite = [(r, v) for r, v in finite if v > float("-inf")]
            if not finite:
                continue
            best_val = max(v for _, v in finite)
            for r, v in finite:
                if v == best_val:
                    best_mask[r][col_idx] = True

    # ── Sizing ──────────────────────────────────────────────────────────────
    COL_W, ROW_H, ROW_H_SUB, MODEL_COL_W, STRAT_COL_W = 1.40, 0.32, 0.22, 1.55, 0.95
    footnotes = FOOTNOTES + ([VAL_FOOTNOTE] if split == "val" else [])
    FOOTER_H = 0.16 * len(footnotes) + 0.25
    n_data_cols = n_total_cols
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

    col_cursor = 2
    for grp in GROUPS:
        x = col_left[col_cursor]
        ax.plot([x, x], [row_tops[total_visual_rows], row_tops[1]],
                color="black", linewidth=0.6, linestyle=(0, (3, 3)), zorder=2)
        col_cursor += len(grp["headers"])
    x_strat = col_left[1]
    ax.plot([x_strat, x_strat], [row_tops[total_visual_rows], row_tops[3]],
            color="black", linewidth=0.5, linestyle=(0, (3, 3)), zorder=2)

    split_word = {"val": "Validation", "test": "Test"}[split]
    ax.text((LEFT + RIGHT) / 2, row_mid(0),
            f"Clinical Models (No CDR, post-exclusion): {split_word} Performance at Baseline   "
            f"(mean ± std, seeds 0/1/2, train-fit / held-out {split}, 80/10/10 stratified split)",
            ha="center", va="center", fontsize=9, fontweight="bold", color="black")

    col_cursor = 2
    for grp in GROUPS:
        n = len(grp["headers"]); x0, x1 = col_left[col_cursor], col_left[col_cursor + n]
        ax.text((x0 + x1) / 2, row_mid(1), grp["title"],
                ha="center", va="center", fontsize=8, fontweight="bold", color="black")
        ax.plot([x0 + 0.05, x1 - 0.05], [row_tops[2], row_tops[2]],
                color="black", linewidth=0.7, zorder=4)
        col_cursor += n

    col_cursor = 2
    for grp in GROUPS:
        n = len(grp["headers"]); x0, x1 = col_left[col_cursor], col_left[col_cursor + n]
        ax.text((x0 + x1) / 2, row_mid(2), grp["subtitle"],
                ha="center", va="center", fontsize=6.0, fontstyle="italic", color="#555555")
        col_cursor += n

    ax.text(col_cx[1], row_mid(3), "Strategy",
            ha="center", va="center", fontsize=7.5, fontstyle="italic", color="black")
    col_cursor = 2
    for grp in GROUPS:
        for h in grp["headers"]:
            ax.text(col_cx[col_cursor], row_mid(3), h, ha="center", va="center",
                    fontsize=7.3, fontstyle="italic", color="black")
            col_cursor += 1

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

    ax.add_patch(plt.Rectangle((LEFT, row_tops[total_visual_rows]),
                               RIGHT - LEFT, row_tops[0] - row_tops[total_visual_rows],
                               facecolor="none", edgecolor="black", linewidth=1.5, zorder=5))

    y_fn = row_tops[total_visual_rows] - 0.22
    for line in footnotes:
        ax.text(LEFT, y_fn, line, ha="left", va="top", fontsize=6.3, color="#222222")
        y_fn -= 0.16

    plt.tight_layout(pad=0.1)
    fig.savefig(OUT_DIR / f"combined_model_table_{split}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(OUT_DIR / f"combined_model_table_{split}.png", bbox_inches="tight", dpi=300)
    if split == "test":   # back-compat canonical name
        fig.savefig(OUT_DIR / "combined_model_table.png", bbox_inches="tight", dpi=300)
    plt.close()
    print(f"  Saved → {OUT_DIR / f'combined_model_table_{split}.png'}")


for _split in ("test", "val"):
    print(f"[{_split}]")
    render_split(_split)
