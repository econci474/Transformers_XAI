"""
06b_diagnostic_coverage_table_cumul.py
=======================================
Cumulative diagnostic coverage audit per ADNI yearly window, showing:
  - Two modality columns: Clinical, MRI
  - Each with sub-columns: N | train | val | test
  - Rows grouped by cumulative yearly window (≤bl, ≤m12, ≤m24, ≤m36, ...)
  - N = total number of *observations* (patient-visit pairs) within that
    cumulative window (same patient counted at each visit)
  - Within each window, sub-rows per seed (0, 1, 2) showing
    total: CN=X | MCI=Y | AD=Z observation counts

Uses the clinical-aligned master (bids_ses = VISCODE_long for matched scans)
so visit codes align with the modality_coverage table.

Data sources
------------
  master_clinical_verbose.csv : longitudinal clinical visits (Patient_ID,
                                VISCODE_long, Label_bl_multi, Label_visit_diag)
  master_mri_clinical_matched_labels.csv : MRI-to-clinical matching
                                           (clinical-aligned variant)
  subjects_with_snp_and_mri.tsv : cohort membership
  D:\\ADNI_BIDS_project\\derivatives\\clinical\\tabular\\baseline\\seed_{0,1,2}\\
    {train,val,test}.csv : baseline 80/10/10 splits by seed

Outputs
-------
  outputs/diagnostic_coverage/diagnostic_coverage_table_cumul.csv
  outputs/diagnostic_coverage/diagnostic_coverage_table_cumul.png
  outputs/diagnostic_coverage/diagnostic_coverage_table_cumul.pdf

Run
---
  python clinical_pipeline/06b_diagnostic_coverage_table_cumul.py
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
SPLIT_DIR   = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\tabular\baseline")
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


def session_to_months(s: str) -> int:
    """Convert visit code to months: sc -> -1, bl -> 0, m06 -> 6, m12 -> 12, ..."""
    s = str(s)
    if s == "sc":  return -1
    if s == "bl":  return 0
    m_match = re.match(r"^m(\d+)$", s)
    return int(m_match.group(1)) if m_match else 9999


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

# ── Load seed splits (baseline CSVs → Patient_ID → split mapping) ─────────
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
    print(f"  seed {seed}: {len(mapping)} subjects "
          f"(train={sum(v=='train' for v in mapping.values())} "
          f"val={sum(v=='val' for v in mapping.values())} "
          f"test={sum(v=='test' for v in mapping.values())})")

# ── Add month-number column to clinical and MRI data ──────────────────────
clin["_months"] = clin["VISCODE_long"].apply(session_to_months)
matched_only = matched_only.copy()
matched_only["_months"] = matched_only["bids_ses"].apply(session_to_months)

# ── Define cumulative windows (yearly) ────────────────────────────────────
# ≤bl (months ≤ 0, i.e. sc + bl), ≤m12, ≤m24, ≤m36, ...
# Find the max month present in data with MRI
max_month_clin = clin["_months"].max()
max_month_mri  = matched_only["_months"].max()
max_month      = max(max_month_clin, max_month_mri)

# Build yearly cutoffs: 0, 12, 24, 36, ... up to max_month
cutoffs = [0]  # ≤bl
m = 12
while m <= max_month:
    cutoffs.append(m)
    m += 12

# Labels for display
def cutoff_label(c):
    if c == 0:  return "≤bl"
    return f"≤m{c}"

print(f"\nCumulative windows: {[cutoff_label(c) for c in cutoffs]}")

# ── Per-window cumulative observation counts ──────────────────────────────
# N = number of observations (patient-visit pairs), not unique subjects.
# For diagnosis, each observation uses its own Label_visit_diag.

# Pre-build (Patient_ID, VISCODE_long) → diagnosis lookup for MRI rows
clin_dx_lookup = {}
for _, r in clin[["Patient_ID", "VISCODE_long", "Label_visit_diag"]].dropna(
        subset=["Label_visit_diag"]).iterrows():
    clin_dx_lookup[(r["Patient_ID"], r["VISCODE_long"])] = r["Label_visit_diag"]
print(f"  diagnosis lookup built: {len(clin_dx_lookup)} entries")

# Add diagnosis to MRI observations via lookup
matched_only["_dx"] = matched_only.apply(
    lambda r: clin_dx_lookup.get((r["Patient_ID"], r["bids_ses"])), axis=1)

result_rows = []
for cutoff in cutoffs:
    # Clinical observations up to this cutoff
    clin_cumul = clin[clin["_months"] <= cutoff]
    clin_N = len(clin_cumul)  # observation count

    # MRI observations up to this cutoff
    mri_cumul = matched_only[matched_only["_months"] <= cutoff]
    mri_N = len(mri_cumul)  # observation count

    for seed in SEEDS:
        mapping = seed_split_map[seed]
        row = {"Window": cutoff_label(cutoff), "Cutoff_months": cutoff,
               "Seed": seed, "Clinical_N": clin_N, "MRI_N": mri_N}

        # Clinical observations: split by seed, count diagnoses
        clin_cumul_split = clin_cumul["Patient_ID"].map(mapping)
        for split in SPLITS:
            mask = clin_cumul_split == split
            obs_dx = clin_cumul.loc[mask, "Label_visit_diag"]
            dx_vc = obs_dx.value_counts()
            for dx in DX_ORDER:
                row[f"Clinical_{split}_{dx}"] = int(dx_vc.get(dx, 0))
            row[f"Clinical_{split}_total"] = int(mask.sum())

        # MRI observations: split by seed, count diagnoses
        mri_cumul_split = mri_cumul["Patient_ID"].map(mapping)
        for split in SPLITS:
            mask = mri_cumul_split == split
            obs_dx = mri_cumul.loc[mask, "_dx"]
            dx_vc = obs_dx.value_counts()
            for dx in DX_ORDER:
                row[f"MRI_{split}_{dx}"] = int(dx_vc.get(dx, 0))
            row[f"MRI_{split}_total"] = int(mask.sum())

        result_rows.append(row)

df_result = pd.DataFrame(result_rows)
df_result.to_csv(OUT_DIR / "diagnostic_coverage_table_cumul.csv", index=False)
print(f"\nSaved CSV → {OUT_DIR / 'diagnostic_coverage_table_cumul.csv'}")
print(df_result.head(9).to_string(index=False))

# ── Trim windows: keep only those with cumulative MRI > 0 ─────────────────
display_cutoffs = [c for c in cutoffs
                   if df_result[df_result["Cutoff_months"] == c]["MRI_N"].max() > 0
                   or c == 0]

# ── Render the table ──────────────────────────────────────────────────────

def fmt_dx(total, cn, mci, ad):
    """Format as 'total:  CN=X  MCI=Y  AD=Z'."""
    return f"{total}:  CN={cn}  MCI={mci}  AD={ad}"


# Prepare display data
display_data = []  # list of (row_type, values)

for cutoff in display_cutoffs:
    label = cutoff_label(cutoff)
    win_df = df_result[df_result["Cutoff_months"] == cutoff]
    clin_N = win_df.iloc[0]["Clinical_N"]
    mri_N  = win_df.iloc[0]["MRI_N"]

    display_data.append(("window_header", label, clin_N, mri_N))

    for seed in SEEDS:
        seed_row = win_df[win_df["Seed"] == seed].iloc[0]
        # Clinical: train/val/test with total prepended
        clin_parts = []
        for split in SPLITS:
            total = int(seed_row[f"Clinical_{split}_total"])
            cn  = seed_row[f"Clinical_{split}_CN"]
            mci = seed_row[f"Clinical_{split}_MCI"]
            ad  = seed_row[f"Clinical_{split}_AD"]
            clin_parts.append(fmt_dx(total, cn, mci, ad))

        # MRI: train/val/test with total prepended
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
COL_W = [1.10,   # Window
         0.70,   # Clinical N
         2.20,   # Clinical train
         2.20,   # Clinical val
         2.20,   # Clinical test
         1.30,   # MRI N (+ % of Clinical N)
         2.20,   # MRI train
         2.20,   # MRI val
         2.20]   # MRI test

ROW_H_HEADER = 0.30
ROW_H_WINDOW = 0.28
ROW_H_SEED   = 0.24

# Count rows
n_window_rows = sum(1 for t in display_data if t[0] == "window_header")
n_seed_rows   = sum(1 for t in display_data if t[0] == "seed_row")

# Heights
TITLE_H    = 0.70
SUBTITLE_H = 0.55
HROW1_H    = 0.32
HROW2_H    = 0.28
TOP_PAD    = 0.12
BOT_PAD    = 0.12
FOOTNOTE_H = 0.90

LEFT      = 0.25
RIGHT_PAD = 0.25
fig_w = LEFT + sum(COL_W) + RIGHT_PAD
fig_h = (TOP_PAD + TITLE_H + SUBTITLE_H + HROW1_H + HROW2_H
         + ROW_H_WINDOW * n_window_rows
         + ROW_H_SEED   * n_seed_rows
         + FOOTNOTE_H + BOT_PAD)

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.axis("off")
ax.set_xlim(0, fig_w)
ax.set_ylim(0, fig_h)

# Column positions
col_left = [LEFT]
for w in COL_W:
    col_left.append(col_left[-1] + w)
RIGHT = col_left[-1]
col_cx = [(col_left[i] + col_left[i+1]) / 2 for i in range(n_data_cols)]

# Y positions (top to bottom)
TOP = fig_h - TOP_PAD
y_cursor = TOP

def hline(y, lw=1.0, ls="-"):
    ax.plot([LEFT, RIGHT], [y, y], color="black", linewidth=lw, linestyle=ls,
            solid_capstyle="butt", zorder=3)

# ── Title ──────────────────────────────────────────────────────────────────
y_title_top = y_cursor
y_cursor -= TITLE_H
ax.text((LEFT + RIGHT) / 2, (y_title_top + y_cursor) / 2,
        "ADNI Cumulative Diagnostic Coverage by Yearly Window\n"
        "(clinical-aligned: bids_ses = VISCODE_long for matched scans)",
        ha="center", va="center", fontsize=11, fontweight="bold",
        color="black", linespacing=1.3)

# ── Subtitle ──────────────────────────────────────────────────────────────
y_sub_top = y_cursor
y_cursor -= SUBTITLE_H
n_subjects = len(snp_pids)
ax.text((LEFT + RIGHT) / 2, (y_sub_top + y_cursor) / 2,
        f"SNP+MRI cohort, n={n_subjects} subjects.\n"
        f"80/10/10 train/val/test splits (seeds 0, 1, 2). "
        f"N = observations (patient-visit pairs).\n"
        f"Diagnosis = per-visit diagnosis (Label_visit_diag).",
        ha="center", va="center", fontsize=9, fontstyle="italic",
        color="black", linespacing=1.4)

hline(y_title_top, lw=1.5)
hline(y_cursor, lw=0.8)

# ── Header row 1: Window | Clinical (spanning 4) | MRI (spanning 4) ──────
y_h1_top = y_cursor
y_cursor -= HROW1_H
y_h1_mid = (y_h1_top + y_cursor) / 2

ax.text(col_cx[0], y_h1_mid, "Window",
        ha="center", va="center", fontsize=10, fontweight="bold", color="black")

clin_cx = (col_left[1] + col_left[5]) / 2
ax.text(clin_cx, y_h1_mid, "Clinical",
        ha="center", va="center", fontsize=10, fontweight="bold",
        fontstyle="italic", color="black")

mri_cx = (col_left[5] + col_left[9]) / 2
ax.text(mri_cx, y_h1_mid, "MRI",
        ha="center", va="center", fontsize=10, fontweight="bold",
        fontstyle="italic", color="black")

# ── Header row 2: N | train | val | test ─────────────────────────────────
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

first_window = True
for item in display_data:
    if item[0] == "window_header":
        _, label, clin_N, mri_N = item
        y_row_top = y_cursor
        y_cursor -= ROW_H_WINDOW
        y_mid = (y_row_top + y_cursor) / 2

        # Window label bold
        ax.text(col_cx[0], y_mid, label,
                ha="center", va="center", fontsize=9.5,
                fontweight="bold", color="black")

        # Clinical N (observation count)
        ax.text(col_cx[1], y_mid, f"{clin_N:,}",
                ha="center", va="center", fontsize=9, color="black")

        # MRI N (observation count) with percentage of Clinical N
        pct = (mri_N / clin_N * 100) if clin_N > 0 else 0
        ax.text(col_cx[5], y_mid, f"{mri_N:,} ({pct:.1f}%)",
                ha="center", va="center", fontsize=9, color="black")

        # Dotted horizontal line above each window (except the first)
        if not first_window:
            hline(y_row_top, lw=0.5, ls=(0, (3, 3)))
        first_window = False

    elif item[0] == "seed_row":
        _, seed, clin_parts, mri_parts = item
        y_row_top = y_cursor
        y_cursor -= ROW_H_SEED
        y_mid = (y_row_top + y_cursor) / 2

        # Seed label
        ax.text(col_cx[0], y_mid, f"Seed {seed}",
                ha="center", va="center", fontsize=7.5,
                fontstyle="italic", color="#444444")

        # Clinical train/val/test
        for j, part in enumerate(clin_parts):
            ax.text(col_cx[2 + j], y_mid, part,
                    ha="center", va="center", fontsize=7, color="#333333")

        # MRI train/val/test
        for j, part in enumerate(mri_parts):
            ax.text(col_cx[6 + j], y_mid, part,
                    ha="center", va="center", fontsize=7, color="#333333")

BOTTOM = y_cursor

# ── Vertical lines ────────────────────────────────────────────────────────
ax.plot([xv_ses, xv_ses], [BOTTOM, y_data_top],
        color="black", linewidth=0.6, linestyle=(0, (3, 3)), zorder=2)

ax.plot([xv_mid, xv_mid], [BOTTOM, y_h1_top],
        color="black", linewidth=0.8, linestyle="-", zorder=2)

# ── Bottom rule ───────────────────────────────────────────────────────────
hline(BOTTOM, lw=1.5)

# ── Outer box ─────────────────────────────────────────────────────────────
rect = plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT, y_title_top - BOTTOM,
                      facecolor="none", edgecolor="black", linewidth=1.5, zorder=5)
ax.add_patch(rect)

# ── Footnote ──────────────────────────────────────────────────────────────
foot_y = BOTTOM - 0.15
footnote_raw = (
    f"N = total observations (patient-visit pairs) with that modality "
    f"up to and including the window cutoff. "
    f"The same patient contributes one observation per visit attended. "
    f"Train/val/test splits are 80/10/10 by Patient_ID (subject-level), "
    f"with seeds 0, 1, 2 controlling the random shuffle. "
    f"Diagnosis counts use per-visit diagnosis (Label_visit_diag): "
    f"CN = cognitively normal, MCI = mild cognitive impairment, AD = Alzheimer's disease. "
    f"Subjects not in any split (e.g. screening-only) "
    f"are excluded from the seed sub-rows."
)
foot_wrap_chars = int((RIGHT - LEFT) * 17)
footnote_wrapped = "\n".join(textwrap.wrap(footnote_raw, width=foot_wrap_chars))
ax.text(LEFT, foot_y, footnote_wrapped,
        ha="left", va="top", fontsize=7.5, color="black")

plt.tight_layout(pad=0.1)
fig.savefig(OUT_DIR / "diagnostic_coverage_table_cumul.png", bbox_inches="tight", dpi=300)
fig.savefig(OUT_DIR / "diagnostic_coverage_table_cumul.pdf", bbox_inches="tight", dpi=300)
print(f"\nSaved PNG → {OUT_DIR / 'diagnostic_coverage_table_cumul.png'}")
print(f"Saved PDF → {OUT_DIR / 'diagnostic_coverage_table_cumul.pdf'}")
plt.close()
print("Done.")
