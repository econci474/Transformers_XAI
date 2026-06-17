"""
03c_match_mri_to_clinical_viscode2.py
======================================
Match MRI scans to clinical visits using ADNI's native VISCODE2 codes.

Unlike 03b (which matches on the date-derived VISCODE_long), this script
uses the raw ADNI VISCODE2 label from both sides:
  - MRI side:      bids_ses from mri_visit_dates_snp_subjects.csv (same curated
                   spine as 03b — one row per (subject, visit), deduplicated)
  - Clinical side: VISCODE_2 from master_clinical_tabular.csv

Since both are assigned by the same ADNI scheduling system, the join is
direct and unambiguous — no date-based fallback is needed for the
primary match.

Strategy A: direct VISCODE2 join  (bids_sub, bids_ses) ↔ (Patient_ID, VISCODE_2)
Strategy B: nearest-by-date ±14d  (fallback for unmatched)

Inputs
------
- D:\\ADNI_BIDS_project\\bids\\imaging\\mri_visit_dates_snp_subjects.csv
    SPINE: wide MRI dates per subject — one row per subject, one column per
    ADNI VISCODE (bl/m12/m24/...). Same curated spine as 03b.
- D:\\ADNI_BIDS_project\\derivatives\\clinical\\no_cdr_stratified\\tabular\\longitudinal\\master_clinical_tabular.csv
    Clinical master with VISCODE_2 column
- D:\\ADNI_BIDS_project\\bids\\genotype\\subjects_with_snp_and_mri.tsv
    SNP+MRI cohort filter

Outputs (to derivatives/mri_clinical_matched/viscode2_aligned/)
-------
- master_mri_clinical_matched_viscode2.csv
- mri_clinical_match_summary_viscode2.csv / .png / .pdf

Run
---
  python mri_pipeline/03c_match_mri_to_clinical_viscode2.py
"""

import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

matplotlib.rcParams["font.family"] = "DejaVu Serif"

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT  = Path(r"D:\ADNI_BIDS_project")
MRI_DATES_CSV = PROJECT_ROOT / "bids" / "imaging" / "mri_visit_dates_snp_subjects.csv"  # SPINE
MASTER_CLIN   = PROJECT_ROOT / "derivatives" / "clinical" / "no_cdr_stratified" \
                / "tabular" / "longitudinal" / "master_clinical_tabular.csv"
SNP_TSV       = PROJECT_ROOT / "bids" / "genotype" / "subjects_with_snp_and_mri.tsv"

OUT_DIR       = PROJECT_ROOT / "derivatives" / "mri_clinical_matched" / "viscode_2_aligned"
OUT_MASTER    = OUT_DIR / "master_mri_clinical_matched_viscode2.csv"
OUT_SUMMARY_CSV = OUT_DIR / "mri_clinical_match_summary_viscode2.csv"
OUT_SUMMARY_PNG = OUT_DIR / "mri_clinical_match_summary_viscode2.png"
OUT_SUMMARY_PDF = OUT_DIR / "mri_clinical_match_summary_viscode2.pdf"

LABEL_COLS = [
    "Label_bl_multi", "Label_visit_diag", "Label_visit_granular",
    "Label_1y", "Label_2y", "Label_3y", "Label_4y", "Label_5y",
    "Label_6y", "Label_7y", "Label_8y", "Label_9y", "Label_10y",
]

MATCHED_STATUSES = ("viscode2_exact", "nearest_within_14d")


def normalise_subject(sid: str) -> str:
    """002_S_0413 -> 002S0413; sub-002S0413 -> 002S0413."""
    s = str(sid).strip()
    if s.startswith("sub-"):
        s = s[4:]
    return s.replace("_", "")


def viscode_months(ses: str) -> int:
    s = str(ses).strip().lower()
    if s == "sc":  return -1
    if s == "bl":  return 0
    m = re.match(r"^m(\d+)$", s)
    return int(m.group(1)) if m else 9999


def pid_from_bids(bids_sub: str) -> str:
    """002S0413 -> 002_S_0413."""
    s = str(bids_sub).strip()
    if len(s) >= 5 and s[3] == "S":
        return f"{s[:3]}_S_{s[4:]}"
    return s


# ── Load data ─────────────────────────────────────────────────────────────────
print("=" * 70)
print("03c - Match MRI scans to clinical visits via VISCODE2")
print("=" * 70)

OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---- Load: SPINE = mri_visit_dates_snp_subjects.csv (wide -> long) ----------
# Same curated spine as 03b: one cell per (subject, visit), deduplicated.
print(f"\nLoading MRI spine: {MRI_DATES_CSV}")
md = pd.read_csv(MRI_DATES_CSV)
visit_cols = [c for c in md.columns if c != "participant_id"]
md_long = md.melt(id_vars="participant_id", value_vars=visit_cols,
                  var_name="bids_ses", value_name="mri_date_str")
# Drop empty cells ("n/a" or NaN)
md_long = md_long[md_long["mri_date_str"].notna()
                  & (md_long["mri_date_str"].astype(str).str.lower() != "n/a")]
md_long["mri_date"] = pd.to_datetime(md_long["mri_date_str"], errors="coerce")
md_long = md_long.dropna(subset=["mri_date"])
md_long["bids_sub"] = md_long["participant_id"].apply(normalise_subject)
# adni_viscode == bids_ses for this spine (both are ADNI visit codes)
md_long["adni_viscode"] = md_long["bids_ses"]
ses_map = md_long[["bids_sub", "bids_ses", "adni_viscode", "mri_date"]].drop_duplicates(
    ["bids_sub", "bids_ses"]).reset_index(drop=True)
ses_map["Patient_ID"] = ses_map["bids_sub"].apply(pid_from_bids)
print(f"  Spine: {md['participant_id'].nunique()} subjects, "
      f"{len(ses_map):,} (subject, visit) scans")

# Clinical master
print(f"Loading clinical master: {MASTER_CLIN}")
clin = pd.read_csv(MASTER_CLIN, low_memory=False)
clin["clinical_date"] = pd.to_datetime(clin["Date"], errors="coerce")
clin = clin.dropna(subset=["clinical_date"]).copy()
clin["bids_sub"] = clin["Patient_ID"].apply(normalise_subject)
# Dedupe at (bids_sub, VISCODE_2) — keep earliest
n_pre = len(clin)
clin = clin.sort_values("clinical_date").drop_duplicates(
    subset=["bids_sub", "VISCODE_2"], keep="first")
print(f"  {n_pre} -> {len(clin)} after dedupe (subject, VISCODE_2)")

# SNP cohort filter
print(f"Loading SNP+MRI cohort: {SNP_TSV}")
snp_df = pd.read_csv(SNP_TSV, sep="\t")
snp_pids = set(snp_df["participant_id"].apply(normalise_subject))
ses_map = ses_map[ses_map["bids_sub"].isin(snp_pids)].reset_index(drop=True)
clin = clin[clin["bids_sub"].isin(snp_pids)]
print(f"  Cohort: {len(snp_pids)} subjects, {len(ses_map):,} MRI scans, {len(clin)} clinical rows")
print()


# ── Strategy A: direct VISCODE2 join ─────────────────────────────────────────
print("[Strategy A] Direct VISCODE2 join (adni_viscode ↔ VISCODE_2)")
clin_cols = ["bids_sub", "VISCODE_2", "VISCODE_long", "clinical_date"] + \
            [c for c in LABEL_COLS if c in clin.columns]
a = ses_map.merge(
    clin[clin_cols],
    left_on=["bids_sub", "adni_viscode"],
    right_on=["bids_sub", "VISCODE_2"],
    how="left",
)
a["days_diff_A"] = (a["mri_date"] - a["clinical_date"]).abs().dt.days
a["matched_A"]   = a["clinical_date"].notna()
n_a = int(a["matched_A"].sum())
print(f"  matched: {n_a}/{len(a)} ({100 * n_a / len(a):.1f}%)")
if n_a:
    d = a.loc[a["matched_A"], "days_diff_A"]
    print(f"  |days_diff|  median={d.median():.0f}  p95={d.quantile(0.95):.0f}  "
          f"max={d.max():.0f}  n>14d={(d > 14).sum()}")
print()


# ── Strategy B: nearest-by-date ±14d fallback ─────────────────────────────────
print("[Strategy B] Nearest clinical visit by date (±14d)")
sm_b = ses_map.sort_values("mri_date").reset_index(drop=True)
cl_b = clin[clin_cols].sort_values("clinical_date").reset_index(drop=True)
b14 = pd.merge_asof(
    sm_b, cl_b,
    left_on="mri_date", right_on="clinical_date",
    by="bids_sub",
    direction="nearest",
    tolerance=pd.Timedelta(days=14),
)
b14["matched_B14"] = b14["clinical_date"].notna()
n_b14 = int(b14["matched_B14"].sum())
print(f"  matched ±14d: {n_b14}/{len(b14)} ({100 * n_b14 / len(b14):.1f}%)")
print()


# ── Build master ──────────────────────────────────────────────────────────────
print("[Master] Building master CSV")
master = a.copy()
master["match_status"] = np.where(master["matched_A"], "viscode2_exact", "unmatched")

# Fallback: Strategy B for unmatched
unmatched_mask = master["match_status"] == "unmatched"
n_unmatched = int(unmatched_mask.sum())
if n_unmatched:
    sub_unmatched = (master.loc[unmatched_mask, ["bids_sub", "mri_date"]]
                          .sort_values("mri_date").copy())
    sub_unmatched["_idx"] = sub_unmatched.index
    cl_sorted = clin[clin_cols].sort_values("clinical_date")
    recover = pd.merge_asof(
        sub_unmatched.reset_index(drop=True), cl_sorted,
        left_on="mri_date", right_on="clinical_date",
        by="bids_sub", direction="nearest",
        tolerance=pd.Timedelta(days=14),
    )
    recovered = recover[recover["clinical_date"].notna()]
    print(f"  Strategy B recovered: {len(recovered)}/{n_unmatched}")
    if not recovered.empty:
        for _, rec_row in recovered.iterrows():
            orig_idx = sub_unmatched.iloc[rec_row.name]["_idx"] if rec_row.name < len(sub_unmatched) else None
            if orig_idx is not None:
                for col in ["VISCODE_2", "VISCODE_long", "clinical_date"] + \
                           [c for c in LABEL_COLS if c in clin.columns]:
                    if col in recover.columns:
                        master.loc[orig_idx, col] = rec_row[col]
                master.loc[orig_idx, "match_status"] = "nearest_within_14d"

master["days_diff"] = (master["mri_date"] - master["clinical_date"]).abs().dt.days

# Sort and save
master["_sort"] = master["bids_ses"].apply(viscode_months)
master = master.sort_values(["bids_sub", "_sort", "mri_date"]).drop(columns="_sort")

out_cols = ["bids_sub", "bids_ses", "adni_viscode", "mri_date", "Patient_ID",
            "clinical_date", "days_diff", "match_status",
            "VISCODE_2", "VISCODE_long"] + \
           [c for c in LABEL_COLS if c in master.columns]
out_cols = [c for c in out_cols if c in master.columns]
master[out_cols].to_csv(OUT_MASTER, index=False)
print(f"  Saved -> {OUT_MASTER}")
print(f"  match_status counts:")
for k, v in master["match_status"].value_counts().items():
    print(f"    {k:<25s} {v:>6d}")
print()


# ── Per-session summary table ─────────────────────────────────────────────────
print("[Summary] Building per-session summary")
rows = []
for ses, g in master.groupby("adni_viscode"):
    n = len(g)
    n_v2 = int((g["match_status"] == "viscode2_exact").sum())
    n_14d = int((g["match_status"] == "nearest_within_14d").sum())
    n_um  = int((g["match_status"] == "unmatched").sum())
    rows.append({
        "ses": ses, "n_scans": n,
        "n_matched_viscode2": n_v2,
        "n_matched_nearest_14d": n_14d,
        "n_unmatched": n_um,
    })

summary_df = pd.DataFrame(rows)
summary_df["_sort"] = summary_df["ses"].apply(viscode_months)
summary_df = summary_df.sort_values("_sort").drop(columns="_sort").reset_index(drop=True)

# Total row
totals = {
    "ses": "TOTAL",
    "n_scans": int(summary_df["n_scans"].sum()),
    "n_matched_viscode2": int(summary_df["n_matched_viscode2"].sum()),
    "n_matched_nearest_14d": int(summary_df["n_matched_nearest_14d"].sum()),
    "n_unmatched": int(summary_df["n_unmatched"].sum()),
}
summary_df = pd.concat([summary_df, pd.DataFrame([totals])], ignore_index=True)
summary_df.to_csv(OUT_SUMMARY_CSV, index=False)
print(f"  Saved summary CSV -> {OUT_SUMMARY_CSV}")


# ── Render summary table ─────────────────────────────────────────────────────
HEADERS = ["Session", "N scans", "matched\nVISCODE2",
           "matched\n±14d", "unmatched"]
data_cols = ["ses", "n_scans", "n_matched_viscode2",
             "n_matched_nearest_14d", "n_unmatched"]

table_rows = summary_df[data_cols].astype(str).values.tolist()
n_rows = len(table_rows)
n_cols = len(HEADERS)

COL_W = [0.85, 0.85, 1.20, 1.00, 1.00]
ROW_H = 0.34
HDR_H = 0.55

fig_w = 0.30 + sum(COL_W) + 0.30
fig_h = 0.50 + HDR_H + ROW_H * n_rows + 0.20

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.axis("off")
ax.set_xlim(0, fig_w)
ax.set_ylim(0, fig_h)

LEFT = 0.30
col_left = [LEFT]
for w in COL_W:
    col_left.append(col_left[-1] + w)
col_cx = [(col_left[i] + col_left[i + 1]) / 2 for i in range(n_cols)]
RIGHT = col_left[-1]

TOP = fig_h - 0.20
title_y = TOP
hdr_top = TOP - 0.30
hdr_bot = hdr_top - HDR_H
row_tops = [hdr_bot - i * ROW_H for i in range(n_rows + 1)]

def hline(y, lw=1.0):
    ax.plot([LEFT, RIGHT], [y, y], color="black", linewidth=lw,
            solid_capstyle="butt", zorder=3)

hline(hdr_top, lw=1.5)
hline(hdr_bot, lw=1.0)
hline(row_tops[-1], lw=1.5)

# Title
ax.text((LEFT + RIGHT) / 2, title_y + 0.10,
        "MRI–Clinical Match Status by ADNI Session (VISCODE2 alignment)",
        ha="center", va="center", fontsize=10, fontweight="bold")
ax.text((LEFT + RIGHT) / 2, title_y - 0.07,
        f"SNP+MRI cohort (n={len(snp_pids)} subjects, {len(ses_map)} scans)",
        ha="center", va="center", fontsize=8, fontstyle="italic")

# Headers
for j, h in enumerate(HEADERS):
    ax.text(col_cx[j], (hdr_top + hdr_bot) / 2, h,
            ha="center", va="center", fontsize=8, fontstyle="italic")

# Data rows
for i, row in enumerate(table_rows):
    y = (row_tops[i] + row_tops[i + 1]) / 2
    is_total = (row[0] == "TOTAL")
    for j, val in enumerate(row):
        weight = "bold" if j == 0 or is_total else "normal"
        ax.text(col_cx[j], y, val, ha="center", va="center",
                fontsize=8, fontweight=weight)
    # Dotted line above TOTAL row
    if is_total:
        hline(row_tops[i], lw=0.8)

# Outer box
rect = plt.Rectangle((LEFT, row_tops[-1]),
                      RIGHT - LEFT, hdr_top - row_tops[-1],
                      facecolor="none", edgecolor="black", linewidth=1.5, zorder=5)
ax.add_patch(rect)

plt.tight_layout(pad=0.1)
fig.savefig(OUT_SUMMARY_PNG, dpi=300, bbox_inches="tight")
fig.savefig(OUT_SUMMARY_PDF, bbox_inches="tight")
plt.close()
print(f"  Saved summary PNG -> {OUT_SUMMARY_PNG}")
print(f"  Saved summary PDF -> {OUT_SUMMARY_PDF}")
print("\nDone.")
