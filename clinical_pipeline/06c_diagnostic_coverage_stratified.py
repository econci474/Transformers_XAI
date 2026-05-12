"""
06c_diagnostic_coverage_stratified.py
======================================
Diagnostic coverage table using the STRATIFIED baseline splits from
01c_build_clinical_csv_stratified.py (baseline_stratified/) instead of the
original random-shuffle splits (baseline/).

Identical rendering logic to 06_diagnostic_coverage_table.py, but reads
from baseline_stratified/ and annotates the output accordingly.

Outputs
-------
  outputs/diagnostic_coverage/diagnostic_coverage_table_stratified.csv
  outputs/diagnostic_coverage/diagnostic_coverage_table_stratified.png
  outputs/diagnostic_coverage/diagnostic_coverage_table_stratified.pdf

Run
---
  python clinical_pipeline/06c_diagnostic_coverage_stratified.py
"""

import re
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"

# ── Paths ──────────────────────────────────────────────────────────────────────
REPO_ROOT   = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
MASTER_CLIN = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\verbose\longitudinal\master_clinical_verbose.csv")
MASTER_MRI  = Path(r"D:\ADNI_BIDS_project\derivatives\mri_clinical_matched\clinical_aligned_intermediates\master_mri_clinical_matched_labels.csv")
SNP_TSV     = Path(r"D:\ADNI_BIDS_project\bids\genotype\subjects_with_snp_and_mri.tsv")
# KEY CHANGE: use stratified splits
SPLIT_DIR   = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\tabular\baseline_stratified")
OUT_DIR     = REPO_ROOT / "clinical_pipeline" / "outputs" / "diagnostic_coverage"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SEEDS  = [0, 1, 2]
SPLITS = ["train", "val", "test"]
MATCHED_STATUSES = ("viscode_exact", "nearest_within_14d")
DX_ORDER = ["CN", "MCI", "AD"]  # display order


# ── Helpers ────────────────────────────────────────────────────────────────────
def to_pid(participant_id: str):
    """'sub-002S0413' -> '002_S_0413'."""
    m = re.match(r"sub-(\d+)S(\d+)", str(participant_id))
    return f"{m.group(1)}_S_{m.group(2)}" if m else None


def session_order(s: str) -> int:
    """Sort key for visit codes: sc < bl < m03 < m06 < ... < m240."""
    s = str(s)
    if s == "sc":  return -1
    if s == "bl":  return 0
    m = re.match(r"^m(\d+)$", s)
    return int(m.group(1)) if m else 9999


# ── Load sources ──────────────────────────────────────────────────────────────
print(f"Loading clinical master    : {MASTER_CLIN}")
clin = pd.read_csv(MASTER_CLIN,
                    usecols=["Patient_ID", "VISCODE_long", "Date",
                             "Label_bl_multi", "Label_visit_diag"],
                    low_memory=False)
clin["Date"] = pd.to_datetime(clin["Date"], errors="coerce")
clin = clin.dropna(subset=["Date"]).copy()

print(f"Loading MRI match table    : {MASTER_MRI}")
mri_match = pd.read_csv(MASTER_MRI, low_memory=False)
mri_match["Patient_ID"] = mri_match["bids_sub"].apply(
    lambda s: f"{str(s)[:3]}_S_{str(s)[4:]}" if str(s)[3:4] == "S" else None)

print(f"Loading SNP+MRI TSV        : {SNP_TSV}")
snp_df  = pd.read_csv(SNP_TSV, sep="\t")
snp_df["Patient_ID"] = snp_df["participant_id"].apply(to_pid)
snp_pids = set(snp_df["Patient_ID"].dropna().unique())

# Restrict to the SNP+MRI cohort
clin      = clin[clin["Patient_ID"].isin(snp_pids)].copy()
mri_match = mri_match[mri_match["Patient_ID"].isin(snp_pids)].copy()
matched_only = mri_match[mri_match["match_status"].isin(MATCHED_STATUSES)].copy()

print(f"  cohort size: {len(snp_pids)} subjects")
print(f"  clinical rows in cohort : {len(clin)}")
print(f"  matched MRI scans       : {len(matched_only)}")

# ── Load STRATIFIED seed splits ───────────────────────────────────────────────
seed_split_map = {}   # {seed: {Patient_ID: "train"/"val"/"test"}}
for seed in SEEDS:
    seed_dir = SPLIT_DIR / f"seed_{seed}"
    mapping = {}
    for split_name in SPLITS:
        csv_path = seed_dir / f"{split_name}.csv"
        df = pd.read_csv(csv_path, usecols=["Patient_ID"], low_memory=False)
        for pid in df["Patient_ID"].unique():
            mapping[pid] = split_name
    seed_split_map[seed] = mapping
    print(f"  seed {seed} (stratified): {len(mapping)} subjects "
          f"(train={sum(v=='train' for v in mapping.values())} "
          f"val={sum(v=='val' for v in mapping.values())} "
          f"test={sum(v=='test' for v in mapping.values())})")

# ── Build the visit list (clinical-aligned sessions) ──────────────────────
sessions = sorted(
    set(clin["VISCODE_long"].dropna().astype(str).unique()),
    key=session_order)
sessions = [s for s in sessions if s != "unknown"]

# ── Per-session diagnostic counts ─────────────────────────────────────────
result_rows = []
for ses in sessions:
    clin_at_ses = clin[clin["VISCODE_long"] == ses].copy()
    clin_pids   = set(clin_at_ses["Patient_ID"].unique())
    clin_N      = len(clin_pids)

    mri_at_ses = matched_only[matched_only["bids_ses"] == ses]
    mri_pids   = set(mri_at_ses["Patient_ID"].dropna().unique())
    mri_N      = len(mri_pids)

    clin_dx = (clin_at_ses
               .dropna(subset=["Label_visit_diag"])
               .drop_duplicates("Patient_ID")
               .set_index("Patient_ID")["Label_visit_diag"]
               .to_dict())

    for seed in SEEDS:
        mapping = seed_split_map[seed]
        row = {"Session": ses, "Seed": seed,
               "Clinical_N": clin_N, "MRI_N": mri_N}

        for modality, pid_set in [("Clinical", clin_pids), ("MRI", mri_pids)]:
            for split in SPLITS:
                pids_in_split = {p for p in pid_set if mapping.get(p) == split}
                dx_counts = {"CN": 0, "MCI": 0, "AD": 0}
                for pid in pids_in_split:
                    dx = clin_dx.get(pid)
                    if dx in dx_counts:
                        dx_counts[dx] += 1
                for dx in DX_ORDER:
                    row[f"{modality}_{split}_{dx}"] = dx_counts[dx]
                row[f"{modality}_{split}_total"] = len(pids_in_split)

        result_rows.append(row)

df_result = pd.DataFrame(result_rows)
df_result.to_csv(OUT_DIR / "diagnostic_coverage_table_stratified.csv", index=False)
print(f"\nSaved CSV → {OUT_DIR / 'diagnostic_coverage_table_stratified.csv'}")
print(df_result.head(9).to_string(index=False))

# ── Trim sessions: keep only those with MRI > 0 (plus always sc, bl) ──────
last_mri_idx = max(
    (sessions.index(s) for s in sessions
     if df_result[df_result["Session"] == s]["MRI_N"].max() > 0),
    default=-1)
display_sessions = [s for i, s in enumerate(sessions)
                    if i <= last_mri_idx or s in ("sc", "bl")]

# ── Render the table ──────────────────────────────────────────────────────
def fmt_dx(total, cn, mci, ad):
    """Format as 'total: CN=X  MCI=Y  AD=Z'."""
    return f"{total}:  CN={cn}  MCI={mci}  AD={ad}"


# Prepare display data
display_data = []
for ses in display_sessions:
    ses_df = df_result[df_result["Session"] == ses]
    clin_N = ses_df.iloc[0]["Clinical_N"]
    mri_N  = ses_df.iloc[0]["MRI_N"]
    display_data.append(("visit_header", ses, clin_N, mri_N))

    for seed in SEEDS:
        seed_row = ses_df[ses_df["Seed"] == seed].iloc[0]
        clin_parts = []
        for split in SPLITS:
            total = int(seed_row[f"Clinical_{split}_total"])
            cn  = seed_row[f"Clinical_{split}_CN"]
            mci = seed_row[f"Clinical_{split}_MCI"]
            ad  = seed_row[f"Clinical_{split}_AD"]
            clin_parts.append(fmt_dx(total, cn, mci, ad))

        mri_parts = []
        for split in SPLITS:
            total = int(seed_row[f"MRI_{split}_total"])
            cn  = seed_row[f"MRI_{split}_CN"]
            mci = seed_row[f"MRI_{split}_MCI"]
            ad  = seed_row[f"MRI_{split}_AD"]
            mri_parts.append(fmt_dx(total, cn, mci, ad))

        display_data.append(("seed_row", seed, clin_parts, mri_parts))

# ── Matplotlib table rendering ────────────────────────────────────────────
n_data_cols = 9
COL_W = [1.10,   # Session
         0.60,   # Clinical N
         2.10,   # Clinical train
         2.10,   # Clinical val
         2.10,   # Clinical test
         1.20,   # MRI N (+ % of Clinical N)
         2.10,   # MRI train
         2.10,   # MRI val
         2.10]   # MRI test

ROW_H_HEADER = 0.30
ROW_H_VISIT  = 0.28
ROW_H_SEED   = 0.24

n_visit_rows = sum(1 for t in display_data if t[0] == "visit_header")
n_seed_rows  = sum(1 for t in display_data if t[0] == "seed_row")

TITLE_H    = 0.70
SUBTITLE_H = 0.50
HROW1_H    = 0.32
HROW2_H    = 0.28
TOP_PAD    = 0.12
BOT_PAD    = 0.12
FOOTNOTE_H = 0.80

LEFT      = 0.25
RIGHT_PAD = 0.25
fig_w = LEFT + sum(COL_W) + RIGHT_PAD
fig_h = (TOP_PAD + TITLE_H + SUBTITLE_H + HROW1_H + HROW2_H
         + ROW_H_VISIT * n_visit_rows
         + ROW_H_SEED  * n_seed_rows
         + FOOTNOTE_H + BOT_PAD)

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.axis("off")
ax.set_xlim(0, fig_w)
ax.set_ylim(0, fig_h)

col_left = [LEFT]
for w in COL_W:
    col_left.append(col_left[-1] + w)
RIGHT = col_left[-1]
col_cx = [(col_left[i] + col_left[i+1]) / 2 for i in range(n_data_cols)]

TOP = fig_h - TOP_PAD
y_cursor = TOP

def hline(y, lw=1.0, ls="-"):
    ax.plot([LEFT, RIGHT], [y, y], color="black", linewidth=lw, linestyle=ls,
            solid_capstyle="butt", zorder=3)

# ── Title ──────────────────────────────────────────────────────────────────
y_title_top = y_cursor
y_cursor -= TITLE_H
ax.text((LEFT + RIGHT) / 2, (y_title_top + y_cursor) / 2,
        "ADNI Diagnostic Coverage by Visit — STRATIFIED Splits\n"
        "(stratified on baseline 3-class diagnosis: CN / MCI / AD)",
        ha="center", va="center", fontsize=11, fontweight="bold",
        color="black", linespacing=1.3)

# ── Subtitle ──────────────────────────────────────────────────────────────
y_sub_top = y_cursor
y_cursor -= SUBTITLE_H
n_subjects = len(snp_pids)
ax.text((LEFT + RIGHT) / 2, (y_sub_top + y_cursor) / 2,
        f"SNP+MRI cohort, n={n_subjects} subjects.\n"
        f"80/10/10 train/val/test splits (seeds 0, 1, 2), "
        f"stratified by Label_bl_multi. "
        f"Diagnosis = per-visit diagnosis (Label_visit_diag).",
        ha="center", va="center", fontsize=9, fontstyle="italic",
        color="black", linespacing=1.4)

hline(y_title_top, lw=1.5)
hline(y_cursor, lw=0.8)

# ── Header row 1 ─────────────────────────────────────────────────────────
y_h1_top = y_cursor
y_cursor -= HROW1_H
y_h1_mid = (y_h1_top + y_cursor) / 2
ax.text(col_cx[0], y_h1_mid, "Session",
        ha="center", va="center", fontsize=10, fontweight="bold", color="black")
clin_cx = (col_left[1] + col_left[5]) / 2
ax.text(clin_cx, y_h1_mid, "Clinical",
        ha="center", va="center", fontsize=10, fontweight="bold",
        fontstyle="italic", color="black")
mri_cx = (col_left[5] + col_left[9]) / 2
ax.text(mri_cx, y_h1_mid, "MRI",
        ha="center", va="center", fontsize=10, fontweight="bold",
        fontstyle="italic", color="black")

# ── Header row 2 ─────────────────────────────────────────────────────────
y_h2_top = y_cursor
y_cursor -= HROW2_H
y_h2_mid = (y_h2_top + y_cursor) / 2
sub_headers = ["N", "train", "val", "test"]
for j, sh in enumerate(sub_headers):
    ax.text(col_cx[1 + j], y_h2_mid, sh,
            ha="center", va="center", fontsize=9, fontstyle="italic", color="black")
    ax.text(col_cx[5 + j], y_h2_mid, sh,
            ha="center", va="center", fontsize=9, fontstyle="italic", color="black")

hline(y_cursor, lw=1.2)

# ── Vertical separators ──────────────────────────────────────────────────
xv_ses = col_left[1]
xv_mid = col_left[5]

# ── Data rows ─────────────────────────────────────────────────────────────
y_data_top = y_cursor

for item in display_data:
    if item[0] == "visit_header":
        _, ses, clin_N, mri_N = item
        y_row_top = y_cursor
        y_cursor -= ROW_H_VISIT
        y_mid = (y_row_top + y_cursor) / 2

        ax.text(col_cx[0], y_mid, ses,
                ha="center", va="center", fontsize=9.5,
                fontweight="bold", color="black")
        ax.text(col_cx[1], y_mid, str(clin_N),
                ha="center", va="center", fontsize=9, color="black")

        pct = (mri_N / clin_N * 100) if clin_N > 0 else 0
        ax.text(col_cx[5], y_mid, f"{mri_N} ({pct:.1f}%)",
                ha="center", va="center", fontsize=9, color="black")

        if ses != display_sessions[0]:
            hline(y_row_top, lw=0.5, ls=(0, (3, 3)))

    elif item[0] == "seed_row":
        _, seed, clin_parts, mri_parts = item
        y_row_top = y_cursor
        y_cursor -= ROW_H_SEED
        y_mid = (y_row_top + y_cursor) / 2

        ax.text(col_cx[0], y_mid, f"Seed {seed}",
                ha="center", va="center", fontsize=7.5,
                fontstyle="italic", color="#444444")

        for j, part in enumerate(clin_parts):
            ax.text(col_cx[2 + j], y_mid, part,
                    ha="center", va="center", fontsize=7, color="#333333")

        for j, part in enumerate(mri_parts):
            ax.text(col_cx[6 + j], y_mid, part,
                    ha="center", va="center", fontsize=7, color="#333333")

BOTTOM = y_cursor

# ── Vertical lines ────────────────────────────────────────────────────────
ax.plot([xv_ses, xv_ses], [BOTTOM, y_data_top],
        color="black", linewidth=0.6, linestyle=(0, (3, 3)), zorder=2)
ax.plot([xv_mid, xv_mid], [BOTTOM, y_h1_top],
        color="black", linewidth=0.8, linestyle="-", zorder=2)

# ── Bottom rule & box ─────────────────────────────────────────────────────
hline(BOTTOM, lw=1.5)
rect = plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT, y_title_top - BOTTOM,
                      facecolor="none", edgecolor="black", linewidth=1.5, zorder=5)
ax.add_patch(rect)

# ── Footnote ──────────────────────────────────────────────────────────────
foot_y = BOTTOM - 0.15
footnote_raw = (
    f"N = total unique subjects with that modality at that visit. "
    f"Train/val/test splits are 80/10/10 STRATIFIED by baseline "
    f"3-class diagnosis (Label_bl_multi: CN / MCI / AD) at the "
    f"Patient_ID level (subject-level), using "
    f"sklearn.model_selection.train_test_split(stratify=...). "
    f"Seeds 0, 1, 2 control the random state. "
    f"Diagnosis counts use per-visit diagnosis (Label_visit_diag)."
)
foot_wrap_chars = int((RIGHT - LEFT) * 17)
footnote_wrapped = "\n".join(textwrap.wrap(footnote_raw, width=foot_wrap_chars))
ax.text(LEFT, foot_y, footnote_wrapped,
        ha="left", va="top", fontsize=7.5, color="black")

plt.tight_layout(pad=0.1)
fig.savefig(OUT_DIR / "diagnostic_coverage_table_stratified.png", bbox_inches="tight", dpi=300)
fig.savefig(OUT_DIR / "diagnostic_coverage_table_stratified.pdf", bbox_inches="tight", dpi=300)
print(f"\nSaved PNG → {OUT_DIR / 'diagnostic_coverage_table_stratified.png'}")
print(f"Saved PDF → {OUT_DIR / 'diagnostic_coverage_table_stratified.pdf'}")
plt.close()
print("Done.")
