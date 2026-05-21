r"""
07b_conversion_group_coverage_extended.py
==========================================
Extended version of `conversion_group_coverage_viscode2.png`: each cell
now shows both the SUBJECT count AND, in parentheses, the total number
of scans (MRI) / visits (Clinical) contributed by those subjects across
all longitudinal timepoints.

Format per cell:
  - SNP column:        N_subjects        (SNP is collected once at baseline)
  - MRI column:        N_subjects (N_scans)
  - Clinical column:   N_subjects (N_visits)

Example: AD_bl has 36 subjects who contribute X MRI scans and Y clinical
visits across all timepoints — cells read "36 (X)" and "36 (Y)".

Data sources:
  - Conversion groups: conversion_labels.tsv (updated by script 07)
  - MRI scan rows:     master_mri_clinical_matched_viscode2.csv,
                       restricted to match_status ∈ {viscode2_exact,
                       nearest_within_14d}
  - Clinical visits:   no_cdr_stratified/tabular/longitudinal/
                       master_clinical_tabular.csv  (each row = 1 visit)

Output:
  Transformers_XAI/clinical_pipeline/outputs/conversion_coverage/
    conversion_group_coverage_viscode2_extension.{png,pdf,csv}
"""
import argparse
import re
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"

REPO_ROOT     = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
sys.path.insert(0, str(REPO_ROOT))
from bidsification.exclusions import is_excluded_subject  # noqa: E402

# ── CLI: --post-exclusion swaps inputs to the sister derivative trees and
#        suffixes all output filenames with `_post_exclusion`.
_ap = argparse.ArgumentParser()
_ap.add_argument("--post-exclusion", action="store_true",
                  help="Apply diagnostic-reversion + corrupted-MRI exclusion, "
                       "read from VISCODE_2_aligned_extended_post_exclusion/, "
                       "and suffix output filenames with _post_exclusion.")
_args = _ap.parse_args()
POST_EXCLUSION = bool(_args.post_exclusion)
SUFFIX         = "_post_exclusion" if POST_EXCLUSION else ""

CONV_TSV      = Path(r"D:\ADNI_SNP_Omni2.5M_20140220\conversion_labels\conversion_labels.tsv")
# In post-exclusion mode, read from the extended_post_exclusion sister master.
if POST_EXCLUSION:
    MASTER_MRI = Path(r"D:\ADNI_BIDS_project\derivatives\mri_clinical_matched"
                       r"\VISCODE_2_aligned_extended_post_exclusion"
                       r"\master_mri_clinical_matched_viscode2_extended.csv")
    MASTER_CLIN = Path(r"D:\ADNI_BIDS_project\derivatives\clinical"
                       r"\no_cdr_stratified_post_exclusion\tabular\longitudinal"
                       r"\master_clinical_tabular.csv")
else:
    MASTER_MRI = Path(r"D:\ADNI_BIDS_project\derivatives\mri_clinical_matched"
                       r"\viscode_2_aligned\master_mri_clinical_matched_viscode2.csv")
    MASTER_CLIN = Path(r"D:\ADNI_BIDS_project\derivatives\clinical"
                       r"\no_cdr_stratified\tabular\longitudinal"
                       r"\master_clinical_tabular.csv")
SNP_TSV       = Path(r"D:\ADNI_BIDS_project\bids\genotype\subjects_with_snp_and_mri.tsv")
OUT_DIR       = REPO_ROOT / "clinical_pipeline" / "outputs" / "conversion_coverage"
OUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"[mode] POST_EXCLUSION = {POST_EXCLUSION}")
print(f"[mode] MASTER_MRI  = {MASTER_MRI}")
print(f"[mode] OUT_SUFFIX  = '{SUFFIX}'")

MATCHED_STATUSES = ("viscode2_exact", "nearest_within_14d")


def to_pid(participant_id: str):
    m = re.match(r"sub-(\d+)S(\d+)", str(participant_id))
    return f"{m.group(1)}_S_{m.group(2)}" if m else None


# ── 1. Load conversion groups ─────────────────────────────────────────────
print(f"Loading conversion labels: {CONV_TSV}")
cv = pd.read_csv(CONV_TSV, sep="\t")
cv["Patient_ID"] = cv["Patient_ID"].astype(str)
print(f"  {len(cv)} subjects (pre-exclusion)")
# Confirm script-07 ran — these columns must exist
required = ["AD_bl", "AD_final", "pMCI", "sMCI", "pCN_to_AD", "pCN_to_MCI",
            "sCN", "Excluded"]
miss = [c for c in required if c not in cv.columns]
if miss:
    raise SystemExit(f"[FATAL] missing flags: {miss} — run script 07 first.")

# Apply exclusion filter when in post-exclusion mode.
if POST_EXCLUSION:
    keep_mask = ~cv["Patient_ID"].apply(
        lambda p: is_excluded_subject(p, include_diagnostic_reversion=True))
    n_drop = int((~keep_mask).sum())
    cv = cv[keep_mask].reset_index(drop=True)
    print(f"  POST_EXCLUSION: dropped {n_drop} subjects → {len(cv)} remain")

# ── 2. SNP+MRI cohort ─────────────────────────────────────────────────────
print(f"Loading SNP+MRI cohort: {SNP_TSV}")
snp_df = pd.read_csv(SNP_TSV, sep="\t")
snp_df["Patient_ID"] = snp_df["participant_id"].apply(to_pid)
snp_pids = set(snp_df["Patient_ID"].dropna().astype(str).unique())
print(f"  SNP+MRI cohort: {len(snp_pids)} subjects (pre-exclusion)")
if POST_EXCLUSION:
    snp_pids = {p for p in snp_pids
                if not is_excluded_subject(p, include_diagnostic_reversion=True)}
    print(f"  POST_EXCLUSION SNP+MRI cohort: {len(snp_pids)} subjects")

# ── 3. Per-subject MRI scan counts (VISCODE2-matched only) ────────────────
# Screening-visit exclusion: by ADNI convention MRI is not collected at
# screening (sc); the MRI master should contain no sc rows. We verify
# this and (defensively) drop any if present.
print(f"Loading MRI VISCODE2 master: {MASTER_MRI}")
mri = pd.read_csv(MASTER_MRI, low_memory=False)
mri["Patient_ID"] = mri["Patient_ID"].astype(str)
mri["_dx"] = mri["Label_visit_diag"].astype(str).str.upper()
if POST_EXCLUSION:
    mri = mri[~mri["Patient_ID"].apply(
        lambda p: is_excluded_subject(p, include_diagnostic_reversion=True))].reset_index(drop=True)
    print(f"  POST_EXCLUSION MRI rows: {len(mri)}")
n_mri_sc = int(((mri["VISCODE_2"].astype(str) == "sc")
                 | (mri["adni_viscode"].astype(str) == "sc")).sum())
print(f"  MRI rows with VISCODE_2/adni_viscode == 'sc' : {n_mri_sc} "
      f"(expected 0; ADNI MRI is not done at screening)")
mri = mri[(mri["VISCODE_2"].astype(str) != "sc")
           & (mri["adni_viscode"].astype(str) != "sc")]
mri_matched = mri[mri["match_status"].isin(MATCHED_STATUSES)]
mri_scans_per_pid = mri_matched.groupby("Patient_ID").size()
# By-label MRI scan counts per PID
mri_AD_per_pid  = mri_matched[mri_matched["_dx"] == "AD"].groupby("Patient_ID").size()
mri_MCI_per_pid = mri_matched[mri_matched["_dx"] == "MCI"].groupby("Patient_ID").size()
mri_CN_per_pid  = mri_matched[mri_matched["_dx"] == "CN"].groupby("Patient_ID").size()
print(f"  matched MRI scans: {len(mri_matched)} across "
      f"{mri_scans_per_pid.size} subjects")
print(f"    of which AD-labelled: {int(mri_AD_per_pid.sum())}, "
      f"MCI-labelled: {int(mri_MCI_per_pid.sum())}, "
      f"CN-labelled: {int(mri_CN_per_pid.sum())}")
mri_pids = set(mri_matched["Patient_ID"].astype(str).unique())

# ── 4. Per-subject clinical visit counts (screening 'sc' excluded) ───────
# User-directed: drop screening visits from clinical counts (they're not
# part of the longitudinal modelling timeline).
print(f"Loading clinical longitudinal master: {MASTER_CLIN}")
clin = pd.read_csv(MASTER_CLIN,
                    usecols=["Patient_ID", "VISCODE_2", "Date",
                             "Label_visit_diag"],
                    low_memory=False)
clin["Patient_ID"] = clin["Patient_ID"].astype(str)
clin = clin.dropna(subset=["Date"])
if POST_EXCLUSION:
    clin = clin[~clin["Patient_ID"].apply(
        lambda p: is_excluded_subject(p, include_diagnostic_reversion=True))].reset_index(drop=True)
    print(f"  POST_EXCLUSION clinical rows: {len(clin)}")
n_clin_sc = int((clin["VISCODE_2"].astype(str) == "sc").sum())
print(f"  Excluding {n_clin_sc} screening (sc) visit rows from clinical")
clin = clin[clin["VISCODE_2"].astype(str) != "sc"]
clin["_dx"] = clin["Label_visit_diag"].astype(str).str.upper()
clin_visits_per_pid = clin.groupby("Patient_ID").size()
clin_AD_per_pid  = clin[clin["_dx"] == "AD"].groupby("Patient_ID").size()
clin_MCI_per_pid = clin[clin["_dx"] == "MCI"].groupby("Patient_ID").size()
clin_CN_per_pid  = clin[clin["_dx"] == "CN"].groupby("Patient_ID").size()
print(f"  clinical visit rows: {len(clin)} across "
      f"{clin_visits_per_pid.size} subjects")
print(f"    of which AD-labelled: {int(clin_AD_per_pid.sum())}, "
      f"MCI-labelled: {int(clin_MCI_per_pid.sum())}, "
      f"CN-labelled: {int(clin_CN_per_pid.sum())}")
clin_pids = set(clin["Patient_ID"].astype(str).unique())

# Defining label per conversion group (for the [n_label, pct%] bracket)
GROUP_LABEL = {
    "AD_bl":      "AD",
    "AD_final":   "AD",
    "pMCI":       "AD",
    "pCN_to_AD":  "AD",
    "sMCI":       "MCI",
    "pCN_to_MCI": "MCI",
    "sCN":        "CN",
    "Excluded":   None,    # mixed → no defining label
}
MRI_BY_LABEL  = {"AD": mri_AD_per_pid,  "MCI": mri_MCI_per_pid,  "CN": mri_CN_per_pid}
CLIN_BY_LABEL = {"AD": clin_AD_per_pid, "MCI": clin_MCI_per_pid, "CN": clin_CN_per_pid}

# ── 5. Per-group counts ───────────────────────────────────────────────────
# Display order — same as the existing coverage table
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

rows = []
for g in GROUPS:
    pids = set(cv.loc[cv[g] == 1, "Patient_ID"].astype(str)) & snp_pids
    pid_list = list(pids)
    # SNP
    snp_n_subj  = len(pids)
    snp_n_rec   = snp_n_subj
    # MRI totals
    mri_n_subj  = len(pids & mri_pids)
    mri_n_rec   = int(mri_scans_per_pid.reindex(pid_list).fillna(0).sum())
    # Clinical totals
    clin_n_subj = len(pids & clin_pids)
    clin_n_rec  = int(clin_visits_per_pid.reindex(pid_list).fillna(0).sum())
    # By defining label
    dl = GROUP_LABEL.get(g)
    if dl is None:
        mri_n_lab = clin_n_lab = -1   # sentinel for "mixed → show '—'"
    else:
        mri_n_lab  = int(MRI_BY_LABEL[dl].reindex(pid_list).fillna(0).sum())
        clin_n_lab = int(CLIN_BY_LABEL[dl].reindex(pid_list).fillna(0).sum())
    rows.append({
        "group": g, "label": GROUP_PRETTY[g], "defining_label": dl,
        "SNP_subj": snp_n_subj, "SNP_rec": snp_n_rec,
        "MRI_subj": mri_n_subj, "MRI_rec": mri_n_rec, "MRI_lab": mri_n_lab,
        "Clinical_subj": clin_n_subj, "Clinical_rec": clin_n_rec, "Clinical_lab": clin_n_lab,
    })

# TOTAL row (over the full SNP+MRI cohort) — defining label is mixed → '—'
all_pids = snp_pids
rows.append({
    "group": "TOTAL", "label": "TOTAL (SNP+MRI cohort)", "defining_label": None,
    "SNP_subj":      len(all_pids),       "SNP_rec":      len(all_pids),
    "MRI_subj":      len(all_pids & mri_pids),
    "MRI_rec":       int(mri_scans_per_pid.reindex(list(all_pids))
                          .fillna(0).sum()),
    "MRI_lab":       -1,
    "Clinical_subj": len(all_pids & clin_pids),
    "Clinical_rec":  int(clin_visits_per_pid.reindex(list(all_pids))
                          .fillna(0).sum()),
    "Clinical_lab":  -1,
})

table_df = pd.DataFrame(rows)
table_df.to_csv(OUT_DIR / f"conversion_group_coverage_viscode2_extension{SUFFIX}.csv",
                 index=False)
print(f"\nSaved CSV → "
      f"{OUT_DIR / ('conversion_group_coverage_viscode2_extension' + SUFFIX + '.csv')}")
print("\n" + table_df.to_string(index=False))

# ── Identify the sMCI subject(s) with 0 matched MRI scans ──────────────
# (User flag: show how many MRI scans + clinical visits they *would*
# contribute, in case they're excluded later.)
sMCI_pids = set(cv.loc[cv["sMCI"] == 1, "Patient_ID"].astype(str)) & snp_pids
sMCI_missing_mri_pids = sMCI_pids - mri_pids
print(f"\nsMCI subjects with 0 matched MRI scans: {len(sMCI_missing_mri_pids)}")
sMCI_missing_info = []
mri_all_per_pid = mri.groupby("Patient_ID").size()
for pid in sorted(sMCI_missing_mri_pids):
    n_clin = int(clin_visits_per_pid.get(pid, 0))
    n_mri_unmatched = int(mri_all_per_pid.get(pid, 0))  # any MRI row (matched or not)
    sMCI_missing_info.append({
        "Patient_ID": pid,
        "clinical_visits": n_clin,
        "mri_rows_total": n_mri_unmatched,   # incl. unmatched
        "mri_scans_matched": 0,
    })
    print(f"  {pid}: clinical_visits={n_clin}  mri_rows(any)={n_mri_unmatched}")

# Save sidecar CSV with the missing-subject detail
miss_df = pd.DataFrame(sMCI_missing_info)
miss_df.to_csv(OUT_DIR / f"conversion_group_coverage_viscode2_extension_sMCI_missing_mri{SUFFIX}.csv",
                index=False)


# ── 6. Render PNG (same house style as the base table) ────────────────────
def fmt(n_subj: int, n_rec: int, show_rec: bool = True,
        n_lab: int = -1, defining_label: str | None = None) -> str:
    """N_subj (N_rec) [N_lab, pct%].  Brackets shown only if defining_label is set."""
    if not show_rec:
        return f"{n_subj:,}"
    base = f"{n_subj:,} ({n_rec:,})"
    if defining_label is None or n_lab < 0:
        return base + "  [—]"
    pct = (100.0 * n_lab / n_rec) if n_rec > 0 else 0.0
    return f"{base}  [{n_lab:,} {defining_label}, {pct:.0f}%]"


n_data_cols = 4
# Numeric cols need to fit "182 (1,430) [1112 AD, 78%]" — needs to be wide.
COL_W = [4.20, 1.40, 3.30, 3.30]
LEFT, RIGHT_PAD = 0.25, 0.25

TITLE_H, HEADER_H = 0.55, 0.40
ROW_H, TOTAL_H = 0.36, 0.44
TOP_PAD, BOT_PAD, FOOTNOTE_H = 0.12, 0.12, 1.40

fig_w = LEFT + sum(COL_W) + RIGHT_PAD
_RIGHT_FOR_WRAP = LEFT + sum(COL_W)

# ── Build subtitle text early so SUBTITLE_H matches the actual wrapped
#    line count. This keeps the subheading text strictly INSIDE the table's
#    bounding rectangle (was overflowing on long footnotes previously).
sub_raw = (
    "Cell format: N_subjects (n_records) [n_label, pct%]. "
    "n_records = total MRI scans (matched VISCODE2) or clinical visits "
    "across all timepoints summed across the subjects in that group. "
    "Square brackets [n_label, pct%] = subset of those records that carry "
    "the group's defining diagnosis label at the visit (Label_visit_diag): "
    "AD for AD_bl/AD_final/pMCI/pCN>AD; MCI for sMCI/pCN>MCI; CN for sCN; "
    "Excluded and TOTAL are mixed → shown as '—'. "
    "Example: AD_bl has 36 subjects contributing 60 MRI scans, of which 60 "
    "are at AD-labelled visits (100%) → 36 (60) [60 AD, 100%]. "
    "SNP is collected once at baseline → N_subj = N_records (single column)."
)
sub_wrap_chars = int((_RIGHT_FOR_WRAP - LEFT) * 16)
sub_wrapped    = "\n".join(textwrap.wrap(sub_raw, width=sub_wrap_chars))
_sub_n_lines   = max(1, sub_wrapped.count("\n") + 1)
# fontsize=8.4 with linespacing=1.4 → ~0.165 in/line in fig coords
SUBTITLE_H = max(0.50, _sub_n_lines * 0.165 + 0.22)

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
y_title_top = y; y -= TITLE_H
ax.text((LEFT + RIGHT) / 2, (y_title_top + y) / 2,
        "Conversion-group Subject + Record Counts by Modality  "
        "(VISCODE2-aligned)",
        ha="center", va="center", fontsize=11.5, fontweight="bold")
# Subtitle (already wrapped above; SUBTITLE_H sized to actual line count
# so the text stays INSIDE the table's bounding rectangle)
y_sub_top = y; y -= SUBTITLE_H
ax.text((LEFT + RIGHT) / 2, (y_sub_top + y) / 2, sub_wrapped,
        ha="center", va="center", fontsize=8.4, fontstyle="italic",
        color="black", linespacing=1.4)
hline(y_title_top, lw=1.5)
hline(y, lw=0.8)

# Header — single-line, simple
y_hdr_top = y; y -= HEADER_H
y_hdr_mid = (y_hdr_top + y) / 2
ax.text(col_left[0] + 0.10, y_hdr_mid, "Conversion group",
        ha="left", va="center", fontsize=10, fontweight="bold", color="black")
ax.text(col_cx[1], y_hdr_mid, "SNP",
        ha="center", va="center", fontsize=10, fontweight="bold", color="black")
ax.text(col_cx[2], y_hdr_mid, "MRI   subj (scans) [at label, %]",
        ha="center", va="center", fontsize=10, fontweight="bold", color="black")
ax.text(col_cx[3], y_hdr_mid, "Clinical   subj (visits) [at label, %]",
        ha="center", va="center", fontsize=10, fontweight="bold", color="black")
hline(y, lw=1.2)

# Data rows
for i, row in enumerate(rows):
    is_total = (row["group"] == "TOTAL")
    rh       = TOTAL_H if is_total else ROW_H
    y_row_top = y; y -= rh
    y_mid = (y_row_top + y) / 2
    fw = "bold" if is_total else "normal"
    fs = 9 if is_total else 8.7

    ax.text(col_left[0] + 0.10, y_mid, row["label"],
            ha="left", va="center", fontsize=fs, fontweight=fw, color="black")
    dl = row.get("defining_label")
    # SNP: subj only (records == subj since SNP is single)
    ax.text(col_cx[1], y_mid, fmt(row["SNP_subj"], row["SNP_rec"], False),
            ha="center", va="center", fontsize=fs, fontweight=fw, color="black")
    # MRI: subj (records) [n_label, pct%].  sMCI gets an asterisk.
    mri_txt = fmt(row["MRI_subj"], row["MRI_rec"], True,
                   n_lab=row["MRI_lab"], defining_label=dl)
    if row["group"] == "sMCI" and sMCI_missing_info:
        mri_txt += " *"
    ax.text(col_cx[2], y_mid, mri_txt,
            ha="center", va="center", fontsize=fs, fontweight=fw, color="black")
    # Clinical: subj (records) [n_label, pct%]
    clin_txt = fmt(row["Clinical_subj"], row["Clinical_rec"], True,
                    n_lab=row["Clinical_lab"], defining_label=dl)
    ax.text(col_cx[3], y_mid, clin_txt,
            ha="center", va="center", fontsize=fs, fontweight=fw, color="black")

    if is_total:
        ax.plot([LEFT, RIGHT], [y_row_top, y_row_top],
                 color="black", linewidth=1.0)
    elif row["group"] == "AD_final":
        hline(y_row_top, lw=0.5, ls=(0, (3, 3)))
    elif row["group"] == "pMCI":
        hline(y_row_top, lw=0.5, ls=(0, (3, 3)))

BOTTOM = y
hline(BOTTOM, lw=1.5)
rect = plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT, y_title_top - BOTTOM,
                      facecolor="none", edgecolor="black",
                      linewidth=1.5, zorder=5)
ax.add_patch(rect)

# Footnote
foot_y = BOTTOM - 0.15
miss_lines = ""
if sMCI_missing_info:
    miss_lines = (
        " * sMCI: " + str(len(sMCI_missing_info))
        + " subject" + ("s" if len(sMCI_missing_info) > 1 else "")
        + " has 0 matched MRI scans and is therefore excluded from the MRI "
        "count (203 sMCI subjects → 202 in MRI column).  Patient_ID = "
        + ", ".join(r["Patient_ID"] for r in sMCI_missing_info)
        + " would otherwise contribute "
        + ", ".join(
            f"{r['mri_rows_total']} MRI rows ({r['mri_rows_total'] - r['mri_scans_matched']} unmatched) "
            f"and {r['clinical_visits']} clinical visits"
            for r in sMCI_missing_info)
        + " — flagged so they can be cleanly excluded downstream "
        "(see sidecar CSV "
        "conversion_group_coverage_viscode2_extension_sMCI_missing_mri.csv).  "
    )
footnote_raw = (
    "Screening visits (VISCODE_2 == 'sc') are EXCLUDED from clinical "
    "counts (they are not part of the longitudinal modelling timeline). "
    "ADNI does not collect MRI at screening, so MRI counts are unaffected "
    f"by this rule ({n_mri_sc} sc-rows found in MRI master — expected 0).  "
    "Partition: AD_bl + pMCI + sMCI + pCN>AD + pCN>MCI + sCN + Excluded = "
    "total cohort. AD_final = AD_bl ∪ pMCI ∪ pCN>AD (overlap, shown above "
    "the partition).  pMCI = baseline-MCI reaching AD sustained.  "
    "pCN>AD = baseline-CN reaching AD sustained.  pCN>MCI = baseline-CN "
    "reaching MCI sustained but NOT AD.  sMCI = baseline-MCI, last visit "
    "MCI, never AD.  sCN = baseline-CN, all visits CN.  Excluded = "
    "reversions / non-sustained.  "
    "SNP column = subjects in the SNP+MRI cohort (1 baseline genotyping "
    "per subject → records = subjects).  "
    "MRI column = subjects with ≥1 matched MRI scan (parens = total "
    "matched scans, summed across the subjects in that group); "
    "match_status ∈ {viscode2_exact, nearest_within_14d}.  "
    "Clinical column = subjects with ≥1 NON-screening clinical visit "
    "(parens = total non-sc visits across all timepoints, summed across "
    "the subjects in that group)." + miss_lines
)
foot_wrap_chars = int((RIGHT - LEFT) * 13)
footnote_wrapped = "\n".join(textwrap.wrap(footnote_raw, width=foot_wrap_chars))
ax.text(LEFT, foot_y, footnote_wrapped,
         ha="left", va="top", fontsize=7.4, color="black")

plt.tight_layout(pad=0.1)
fig.savefig(OUT_DIR / f"conversion_group_coverage_viscode2_extension{SUFFIX}.png",
             bbox_inches="tight", dpi=300)
fig.savefig(OUT_DIR / f"conversion_group_coverage_viscode2_extension{SUFFIX}.pdf",
             bbox_inches="tight", dpi=300)
print(f"\nSaved PNG → "
      f"{OUT_DIR / ('conversion_group_coverage_viscode2_extension' + SUFFIX + '.png')}")
plt.close()
print("Done.")
