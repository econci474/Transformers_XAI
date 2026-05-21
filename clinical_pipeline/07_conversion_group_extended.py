r"""
07_conversion_group_extended.py
================================
Three artefacts in one pass:

  1. **Updated `conversion_labels.tsv`** (in place) with new flags:
        CN_to_MCI                    sustained CN -> MCI (not via AD)
        CN_to_AD                     sustained CN -> AD (= pCN>AD)
        AD_bl                        baseline-AD
        AD_final                     last_dx == AD (= AD_bl ∪ ever_AD)
        pMCI                         progressive MCI (= progressive_MCI)
        sMCI                         stable MCI    (= stable_MCI)
        pCN_to_AD                    sustained CN -> AD (CN_to_AD alias)
        pCN_to_MCI                   sustained CN -> MCI (CN_to_MCI alias)
        sCN                          stable CN     (= stable_CN)
        Excluded                     reversions / non-sustained
        conversion_group             single string label per subject

     Existing columns are KEPT (ever_conversion_*, ad_case, …) so older
     scripts continue to work.

  2. **Extended VISCODE2 MRI/clinical mastersheet** at:
        D:/ADNI_BIDS_project/derivatives/mri_clinical_matched/
          VISCODE_2_aligned_extended/
            master_mri_clinical_matched_viscode2_extended.csv
     = `master_mri_clinical_matched_viscode2.csv` joined with the
     per-subject conversion-group info above.
     Adds columns: bl_dx, last_dx, FU_years, years_to_MCI, years_to_AD,
     conversion_group, AD_bl, AD_final, pMCI, sMCI, pCN_to_AD,
     pCN_to_MCI, sCN, CN_to_MCI, CN_to_AD, Excluded.

  3. **Subject-count coverage table** (PNG + CSV) at:
        clinical_pipeline/outputs/diagnostic_coverage/
          conversion_group_coverage_viscode2.{png,pdf,csv}
     Rows: conversion groups (+ TOTAL).  Columns: SNP, MRI (≥1
     matched scan), Clinical (≥1 visit).  Values are SUBJECT counts.

Run
---
  python clinical_pipeline/07_conversion_group_extended.py
"""

import re
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"

# ── Paths ─────────────────────────────────────────────────────────────────
REPO_ROOT     = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
CONV_TSV      = Path(r"D:\ADNI_SNP_Omni2.5M_20140220\conversion_labels\conversion_labels.tsv")
MASTER_MRI    = Path(r"D:\ADNI_BIDS_project\derivatives\mri_clinical_matched\viscode_2_aligned\master_mri_clinical_matched_viscode2.csv")
EXT_MRI_DIR   = Path(r"D:\ADNI_BIDS_project\derivatives\mri_clinical_matched\VISCODE_2_aligned_extended")
EXT_MRI_DIR_POST = Path(r"D:\ADNI_BIDS_project\derivatives\mri_clinical_matched\VISCODE_2_aligned_extended_post_exclusion")
SNP_TSV       = Path(r"D:\ADNI_BIDS_project\bids\genotype\subjects_with_snp_and_mri.tsv")
OUT_DIR       = REPO_ROOT / "clinical_pipeline" / "outputs" / "conversion_coverage"
OUT_DIR.mkdir(parents=True, exist_ok=True)
EXT_MRI_DIR.mkdir(parents=True, exist_ok=True)
EXT_MRI_DIR_POST.mkdir(parents=True, exist_ok=True)

# Allow `from bidsification.exclusions import …` from this script's directory
sys.path.insert(0, str(REPO_ROOT))
from bidsification.exclusions import is_excluded_subject  # noqa: E402

MATCHED_STATUSES = ("viscode2_exact", "nearest_within_14d")


def to_pid(participant_id: str):
    m = re.match(r"sub-(\d+)S(\d+)", str(participant_id))
    return f"{m.group(1)}_S_{m.group(2)}" if m else None


# ── 1. Load conversion_labels.tsv + compute new flags ─────────────────────
print(f"[1/3] Loading {CONV_TSV} …")
cv = pd.read_csv(CONV_TSV, sep="\t")
cv["Patient_ID"] = cv["Patient_ID"].astype(str)
print(f"      {len(cv)} subjects")

def b(x): return x.fillna(0).astype(int).astype(bool)

bl_dx                = cv["bl_dx"].astype(str)
last_dx              = cv["last_dx"].astype(str)
ever_AD              = b(cv["ever_conversion_AD"])         # CN/MCI sustained -> AD
ever_MCI             = b(cv["ever_conversion_MCI"])        # CN sustained -> MCI (may or may not later AD)
progressive_MCI      = b(cv["progressive_MCI"])
stable_MCI           = b(cv["stable_MCI"])
stable_CN            = b(cv["stable_CN"])
reversion_MCI_to_CN  = b(cv["reversion_MCI_to_CN"])

# Group flags (boolean, mutually exclusive partition modulo Excluded)
AD_bl     = (bl_dx == "AD")
pMCI      = (bl_dx == "MCI") & ever_AD                    # MCI -> AD sustained
pCN_to_AD = (bl_dx == "CN")  & ever_AD                    # CN  -> AD sustained
pCN_to_MCI = (bl_dx == "CN") & ever_MCI & (~ever_AD)      # CN  -> MCI sustained, NOT AD
sMCI      = stable_MCI                                    # MCI -> MCI throughout
sCN       = stable_CN                                     # CN  -> CN throughout

# Anything not in the above partition → Excluded (reversions / non-sustained)
in_partition = AD_bl | pMCI | pCN_to_AD | pCN_to_MCI | sMCI | sCN
Excluded     = ~in_partition

# User-requested flags (descriptive; overlap-allowed):
#   CN_to_MCI  = ever_conversion_MCI = baseline-CN reaching MCI ever
#                (INCLUDES the 25 who later went CN→MCI→AD)
#   CN_to_AD   = baseline-CN reaching AD sustained (= pCN_to_AD)
# For the strict-no-AD partition rows use `pCN_to_MCI` (boolean column),
# which is mutually exclusive with `pCN_to_AD` by construction.
CN_to_MCI = ever_MCI                                        # = ever_conversion_MCI (60)
CN_to_AD  = pCN_to_AD.copy()                                # = pCN_to_AD (25)
AD_final  = AD_bl | pMCI | pCN_to_AD                        # last_dx == AD (sustained)

# Single-label partition (string per subject)
def label_one(i):
    if AD_bl.iloc[i]:        return "AD_bl"
    if pMCI.iloc[i]:         return "pMCI"
    if pCN_to_AD.iloc[i]:    return "pCN_to_AD"
    if pCN_to_MCI.iloc[i]:   return "pCN_to_MCI"
    if sMCI.iloc[i]:         return "sMCI"
    if sCN.iloc[i]:          return "sCN"
    return "Excluded"

cv["conversion_group"] = [label_one(i) for i in range(len(cv))]

# Boolean group flags.
# Partition flags (strict, mutually exclusive — sum to N total):
#   AD_bl, pMCI, pCN_to_AD, pCN_to_MCI, sMCI, sCN, Excluded
# Descriptive flags (overlap-allowed):
#   AD_final, CN_to_AD (= pCN_to_AD), CN_to_MCI (= ever_conversion_MCI)
cv["AD_bl"]      = AD_bl.astype(int)
cv["AD_final"]   = AD_final.astype(int)
cv["pMCI"]       = pMCI.astype(int)
cv["sMCI"]       = sMCI.astype(int)
cv["pCN_to_AD"]  = pCN_to_AD.astype(int)
cv["pCN_to_MCI"] = pCN_to_MCI.astype(int)        # strict (35) — partition row
cv["sCN"]        = sCN.astype(int)
cv["CN_to_AD"]   = CN_to_AD.astype(int)          # = pCN_to_AD (25)
cv["CN_to_MCI"]  = CN_to_MCI.astype(int)         # = ever_conversion_MCI (60) — inclusive
cv["Excluded"]   = Excluded.astype(int)

# Years-to-event — re-derived properly from the longitudinal clinical
# master (the legacy `first_conversion_MCI` field is misleading: it carries
# the time of the first FOLLOW-UP MCI visit for baseline-MCI patients,
# which is NOT a conversion from anything).
#
# Strict definitions used here:
#   years_to_MCI = time from baseline visit to FIRST visit with
#                  Label_visit_diag == MCI,   ONLY IF bl_dx == "CN" AND
#                  the patient ever reached MCI.   NaN otherwise.
#                  (For bl_dx == MCI: NaN — already MCI at baseline.)
#                  (For bl_dx == AD : NaN — already AD; skipped MCI.)
#   years_to_AD  = time from baseline visit to FIRST visit with
#                  Label_visit_diag == AD,    ONLY IF bl_dx in {"CN","MCI"}
#                  AND the patient ever reached AD.   NaN otherwise.
#                  (For bl_dx == AD : NaN — already AD at baseline.)
#   FU_years     = time from baseline visit to LAST visit (= followup_years).
MASTER_CLIN = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified\tabular\longitudinal\master_clinical_tabular.csv")
print(f"\n  Re-deriving years_to_MCI / years_to_AD from {MASTER_CLIN.name} …")
clin = pd.read_csv(MASTER_CLIN,
                    usecols=["Patient_ID", "Date", "Label_visit_diag"],
                    low_memory=False)
clin["Date"] = pd.to_datetime(clin["Date"], errors="coerce")
clin = clin.dropna(subset=["Date"]).copy()
clin["Patient_ID"] = clin["Patient_ID"].astype(str)
# Per-subject: baseline (earliest visit) + last visit + first MCI + first AD
per = clin.groupby("Patient_ID").agg(
    _bl_date=("Date", "min"),
    _last_date=("Date", "max"),
).reset_index()
first_mci = (clin[clin["Label_visit_diag"].astype(str).str.upper() == "MCI"]
              .groupby("Patient_ID")["Date"].min()
              .reset_index().rename(columns={"Date": "_first_MCI_date"}))
first_ad  = (clin[clin["Label_visit_diag"].astype(str).str.upper() == "AD"]
              .groupby("Patient_ID")["Date"].min()
              .reset_index().rename(columns={"Date": "_first_AD_date"}))
per = per.merge(first_mci, on="Patient_ID", how="left") \
          .merge(first_ad,  on="Patient_ID", how="left")
per["_FU_years"] = (per["_last_date"] - per["_bl_date"]).dt.days / 365.25
# Restrict to the cohort we're updating
per = per.set_index("Patient_ID")
cv_bl_dx = cv.set_index("Patient_ID")["bl_dx"]
def _y(per_col, allowed_bl):
    out = pd.Series(np.nan, index=cv["Patient_ID"])
    for pid in cv["Patient_ID"]:
        if pid not in per.index:
            continue
        if cv_bl_dx.loc[pid] not in allowed_bl:
            continue
        d = per.loc[pid, per_col]
        if pd.isna(d):
            continue
        out[pid] = (d - per.loc[pid, "_bl_date"]).days / 365.25
    return out.values

cv["FU_years"]     = cv["Patient_ID"].map(per["_FU_years"]).astype(float)
cv["years_to_MCI"] = _y("_first_MCI_date", {"CN"})
cv["years_to_AD"]  = _y("_first_AD_date",  {"CN", "MCI"})

# Sanity print
n_bl_MCI_with_old_yMCI = int(
    ((cv["bl_dx"] == "MCI") & cv["first_conversion_MCI"].notna()).sum())
n_corrected = int((cv["bl_dx"] == "MCI").sum())   # all bl-MCI now have NaN
print(f"    Legacy first_conversion_MCI had {n_bl_MCI_with_old_yMCI} "
      f"non-NaN values for baseline-MCI patients (incorrect).")
print(f"    Corrected years_to_MCI for all {n_corrected} baseline-MCI "
      f"patients to NaN (already MCI at baseline).")
print(f"    Final years_to_MCI: "
      f"NaN={cv['years_to_MCI'].isna().sum()}  "
      f"valid={cv['years_to_MCI'].notna().sum()}")
print(f"    Final years_to_AD : "
      f"NaN={cv['years_to_AD'].isna().sum()}  "
      f"valid={cv['years_to_AD'].notna().sum()}")

# Verify partition counts add up to N
print(f"\n  Partition counts (must sum to {len(cv)}):")
gp_counts = cv["conversion_group"].value_counts()
print(gp_counts.to_string())
print(f"  Σ = {gp_counts.sum()}  (target {len(cv)})  ",
      "✓" if gp_counts.sum() == len(cv) else "MISMATCH")

# Save updated conversion_labels.tsv (in place)
cv.to_csv(CONV_TSV, sep="\t", index=False)
print(f"\n  Updated → {CONV_TSV}")


# ── 2. Build extended VISCODE2 mastersheet ────────────────────────────────
print(f"\n[2/3] Loading VISCODE2 master {MASTER_MRI} …")
mri = pd.read_csv(MASTER_MRI, low_memory=False)
mri["Patient_ID"] = mri["Patient_ID"].astype(str)
print(f"      {len(mri)} scan rows, "
      f"{mri['Patient_ID'].nunique()} unique subjects")

new_cols = ["bl_dx", "last_dx", "FU_years", "years_to_MCI", "years_to_AD",
             "conversion_group", "AD_bl", "AD_final", "pMCI", "sMCI",
             "pCN_to_AD", "pCN_to_MCI", "sCN", "CN_to_MCI", "CN_to_AD",
             "Excluded"]
join = cv[["Patient_ID"] + new_cols]
ext  = mri.merge(join, on="Patient_ID", how="left")

n_unmatched = ext["conversion_group"].isna().sum()
if n_unmatched:
    print(f"  [warn] {n_unmatched} MRI rows have no conversion_labels match "
          "(probably not in SNP+MRI cohort)")

ext_path = EXT_MRI_DIR / "master_mri_clinical_matched_viscode2_extended.csv"
try:
    ext.to_csv(ext_path, index=False)
    print(f"  Saved extended master → {ext_path}  ({len(ext)} rows)")
except PermissionError as e:
    print(f"  [WARN] could not write {ext_path}: {e}")
    print(f"         Likely open in another process — skipping pre-exclusion save.")

# ── Sister save: post-exclusion variant ───────────────────────────────────
# Drop rows for subjects flagged by bidsification.exclusions.is_excluded_subject
# (corrupted-MRI + site-381 + EXCLUDED_PTID_LIST + diagnostic-reversion opt-in).
# This is the canonical input for the post-exclusion regeneration of splits
# and coverage tables.
ext_post = ext[~ext["Patient_ID"].astype(str).apply(
    lambda p: is_excluded_subject(p, include_diagnostic_reversion=True)
)].reset_index(drop=True)
ext_post_path = EXT_MRI_DIR_POST / "master_mri_clinical_matched_viscode2_extended.csv"
try:
    ext_post.to_csv(ext_post_path, index=False)
    n_pre, n_post = len(ext), len(ext_post)
    pids_pre, pids_post = ext["Patient_ID"].nunique(), ext_post["Patient_ID"].nunique()
    print(f"  Saved POST-EXCLUSION master → {ext_post_path}  ({n_post} rows)")
    print(f"    Pre-exclusion : {n_pre} rows / {pids_pre} subjects")
    print(f"    Post-exclusion: {n_post} rows / {pids_post} subjects  "
          f"(dropped {n_pre - n_post} rows, {pids_pre - pids_post} subjects)")
except PermissionError as e:
    print(f"  [WARN] could not write {ext_post_path}: {e}")


# ── 3. Subject-count coverage table ───────────────────────────────────────
print(f"\n[3/3] Building subject-count coverage table by conversion group …")
print(f"      Loading SNP+MRI TSV {SNP_TSV} …")
snp_df = pd.read_csv(SNP_TSV, sep="\t")
snp_df["Patient_ID"] = snp_df["participant_id"].apply(to_pid)
snp_pids = set(snp_df["Patient_ID"].dropna().astype(str).unique())
print(f"      SNP+MRI cohort: {len(snp_pids)} subjects")

mri_matched = mri[mri["match_status"].isin(MATCHED_STATUSES)]
mri_pids    = set(mri_matched["Patient_ID"].astype(str).unique())
print(f"      Subjects with ≥1 matched MRI scan: {len(mri_pids)}")

# All cohort subjects HAVE clinical (conversion_labels comes from clinical)
# so "Clinical" column counts subjects with ≥1 row in conversion_labels.
clin_pids = set(cv["Patient_ID"].astype(str).unique())
print(f"      Subjects with clinical data: {len(clin_pids)}")

# Order of groups in the table
GROUPS = ["AD_bl", "AD_final", "pMCI", "sMCI", "pCN_to_AD", "pCN_to_MCI",
          "sCN", "Excluded"]
GROUP_PRETTY = {
    "AD_bl":      "AD_bl  (baseline AD)",
    "AD_final":   "AD_final  (last visit == AD)",
    "pMCI":       "pMCI  (MCI → AD sustained)",
    "sMCI":       "sMCI  (stable MCI)",
    "pCN_to_AD":  "pCN > AD  (CN → AD sustained)",
    "pCN_to_MCI": "pCN > MCI  (CN → MCI, no AD)",
    "sCN":        "sCN  (stable CN)",
    "Excluded":   "Excluded  (reversions, non-sustained)",
}
# AD_final is an overlap of AD_bl ∪ pMCI ∪ pCN_to_AD (not a partition row);
# we show it for completeness.

rows = []
for g in GROUPS:
    flag = cv.set_index("Patient_ID")[g].astype(int)
    pids_in_g = set(flag[flag == 1].index.astype(str))
    # Restrict to SNP+MRI cohort (all 616 should be in there, but be safe)
    pids_in_g = pids_in_g & snp_pids

    snp_n      = len(pids_in_g & snp_pids)
    mri_n      = len(pids_in_g & mri_pids)
    clin_n     = len(pids_in_g & clin_pids)
    rows.append({"group": g, "label": GROUP_PRETTY[g],
                 "SNP": snp_n, "MRI": mri_n, "Clinical": clin_n})

total_row = {"group": "TOTAL",
              "label": "TOTAL (SNP+MRI cohort)",
              "SNP": len(snp_pids),
              "MRI": len(snp_pids & mri_pids),
              "Clinical": len(snp_pids & clin_pids)}
rows.append(total_row)

table_df = pd.DataFrame(rows)
table_df.to_csv(OUT_DIR / "conversion_group_coverage_viscode2.csv",
                 index=False)
print(f"\n  Saved CSV → "
      f"{OUT_DIR / 'conversion_group_coverage_viscode2.csv'}")
print("\n" + table_df.to_string(index=False))

# Sanity: partition-only groups should sum to TOTAL
partition_groups = ["AD_bl", "pMCI", "sMCI", "pCN_to_AD", "pCN_to_MCI",
                     "sCN", "Excluded"]
psum = sum(int(table_df.loc[table_df["group"] == g, "SNP"].iloc[0])
            for g in partition_groups)
print(f"\n  Partition sanity: Σ partition_groups[SNP] = {psum}  "
      f"(target {total_row['SNP']})  ",
      "✓" if psum == total_row["SNP"] else "MISMATCH")


# ── Render PNG (20b-style) ────────────────────────────────────────────────
n_data_cols = 4   # label, SNP, MRI, Clinical
COL_W = [4.20, 1.30, 1.30, 1.30]

LEFT      = 0.25
RIGHT_PAD = 0.25

TITLE_H    = 0.55
SUBTITLE_H = 0.45
HEADER_H   = 0.42
ROW_H      = 0.36
TOTAL_H    = 0.44
TOP_PAD    = 0.12
BOT_PAD    = 0.12
FOOTNOTE_H = 0.90

fig_w = LEFT + sum(COL_W) + RIGHT_PAD
fig_h = (TOP_PAD + TITLE_H + SUBTITLE_H + HEADER_H
         + ROW_H * (len(rows) - 1)
         + TOTAL_H
         + FOOTNOTE_H + BOT_PAD)

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.axis("off"); ax.set_xlim(0, fig_w); ax.set_ylim(0, fig_h)

col_left = [LEFT]
for w in COL_W:
    col_left.append(col_left[-1] + w)
RIGHT = col_left[-1]
col_cx = [(col_left[i] + col_left[i+1]) / 2 for i in range(n_data_cols)]

TOP = fig_h - TOP_PAD
y = TOP

def hline(yv, lw=1.0, ls="-"):
    ax.plot([LEFT, RIGHT], [yv, yv], color="black",
             linewidth=lw, linestyle=ls,
             solid_capstyle="butt", zorder=3)

# Title
y_title_top = y
y -= TITLE_H
ax.text((LEFT + RIGHT) / 2, (y_title_top + y) / 2,
        "Conversion-group Subject Counts  by Modality  (VISCODE2-aligned)",
        ha="center", va="center", fontsize=11.5, fontweight="bold",
        color="black")
# Subtitle
y_sub_top = y
y -= SUBTITLE_H
ax.text((LEFT + RIGHT) / 2, (y_sub_top + y) / 2,
        f"SNP+MRI cohort, n = {len(snp_pids)} subjects.  "
        f"AD_final is shown for completeness — it is the union "
        f"AD_bl ∪ pMCI ∪ pCN>AD (NOT a partition row).\n"
        f"Cell values are SUBJECT counts (not scans or visits).",
        ha="center", va="center", fontsize=8.7, fontstyle="italic",
        color="black", linespacing=1.4)
hline(y_title_top, lw=1.5)
hline(y, lw=0.8)

# Header
y_hdr_top = y
y -= HEADER_H
y_hdr_mid = (y_hdr_top + y) / 2
headers = ["Conversion group", "SNP", "MRI", "Clinical"]
for i, h in enumerate(headers):
    ax.text(col_cx[i] if i > 0 else col_left[0] + 0.10,
            y_hdr_mid, h, ha=("left" if i == 0 else "center"),
            va="center", fontsize=10, fontweight="bold", color="black")
hline(y, lw=1.2)

# Data rows
for i, row in enumerate(rows):
    is_total = (row["group"] == "TOTAL")
    rh       = TOTAL_H if is_total else ROW_H
    y_row_top = y
    y -= rh
    y_mid = (y_row_top + y) / 2
    fw = "bold" if is_total else "normal"
    fs = 9 if is_total else 8.7

    ax.text(col_left[0] + 0.10, y_mid, row["label"],
            ha="left", va="center", fontsize=fs, fontweight=fw,
            color=("black" if not is_total else "#000000"))
    for k, c in enumerate(["SNP", "MRI", "Clinical"]):
        ax.text(col_cx[k + 1], y_mid, f"{row[c]:,}",
                ha="center", va="center", fontsize=fs,
                fontweight=fw, color="black")

    # Separator above TOTAL + above AD_final (since AD_final is overlap)
    if is_total:
        ax.plot([LEFT, RIGHT], [y_row_top, y_row_top],
                 color="black", linewidth=1.0)
    elif row["group"] == "AD_final":
        hline(y_row_top, lw=0.5, ls=(0, (3, 3)))
    elif row["group"] == "pMCI":
        # back to partition rows after AD_final overlap
        hline(y_row_top, lw=0.5, ls=(0, (3, 3)))

BOTTOM = y
hline(BOTTOM, lw=1.5)

# Outer box
rect = plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT, y_title_top - BOTTOM,
                      facecolor="none", edgecolor="black",
                      linewidth=1.5, zorder=5)
ax.add_patch(rect)

# Footnote
foot_y = BOTTOM - 0.15
footnote_raw = (
    "Partition: AD_bl + pMCI + sMCI + pCN>AD + pCN>MCI + sCN + Excluded = total cohort. "
    "AD_final = AD_bl ∪ pMCI ∪ pCN>AD (overlap, shown above the partition).  "
    "pMCI = baseline-MCI reaching AD sustained (= progressive_MCI).  "
    "pCN>AD = baseline-CN reaching AD sustained.  "
    "pCN>MCI = baseline-CN reaching MCI sustained but NOT AD (= sustained CN→MCI without AD).  "
    "sMCI = baseline-MCI, last visit MCI, never AD.  "
    "sCN = baseline-CN, all visits CN.  "
    "Excluded = reversions (MCI→CN, AD→non-AD) or other non-sustained transitions.  "
    "SNP column = subjects in the SNP+MRI cohort (genotyped at baseline).  "
    "MRI column = subjects with ≥1 matched MRI scan in master_mri_clinical_matched_viscode2.csv "
    "(match_status ∈ {viscode2_exact, nearest_within_14d}).  "
    "Clinical column = subjects with ≥1 visit in conversion_labels.tsv.  "
    "Source: D:/.../mri_clinical_matched/VISCODE_2_aligned_extended/."
)
foot_wrap_chars = int((RIGHT - LEFT) * 13)
footnote_wrapped = "\n".join(textwrap.wrap(footnote_raw, width=foot_wrap_chars))
ax.text(LEFT, foot_y, footnote_wrapped,
         ha="left", va="top", fontsize=7.4, color="black")

plt.tight_layout(pad=0.1)
fig.savefig(OUT_DIR / "conversion_group_coverage_viscode2.png",
             bbox_inches="tight", dpi=300)
fig.savefig(OUT_DIR / "conversion_group_coverage_viscode2.pdf",
             bbox_inches="tight", dpi=300)
print(f"\n  Saved PNG → "
      f"{OUT_DIR / 'conversion_group_coverage_viscode2.png'}")
plt.close()

print("\nDone.")
