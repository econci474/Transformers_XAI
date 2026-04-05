"""
Script 02 — Build QC Selection Table (v4)
==========================================
Selects ONE native-space structural T1w per BIDS session using
FEATURE-based priority (what corrections are present) rather than
LABEL-based priority (MPR vs MT1).

Rationale:
  ADNI1/GO  → "MPR__GradWarp__B1_Correction__N3" naming
  ADNI2/3   → "MT1__GradWarp__N3m" naming
  The label is an artefact of the ADNI processing pipeline version;
  the features present (GradWarp, B1 correction, N3 bias correction)
  are what matter for sMRIprep input quality.

HARD EXCLUSIONS (never included, regardless of QC):
  - *__Mask                           : binary brain masks
  - Spatially_Normalized__Masked*     : pre-normalised + skull-stripped
  - HHP_6_DOF*                        : hippocampal AC-PC registration
  - HarP_135*                         : hippocampus atlas-registered
  - EPI_current*                      : EPI-corrected (not structural T1)
  - SERIES_QUALITY == 4 (unusable)    : flagged unusable by MAYOADIRL QC

FEATURE-BASED PRIORITY (label-agnostic, lower rank = preferred):
  Rank 0 : GradWarp + B1 correction + N3 bias correction (all three)
           [MPR__GradWarp__B1_Correction__N3, MPR-R variant, ±Scaled]
  Rank 1 : GradWarp + N3 (no B1)
           [MPR__GradWarp__N3, MT1__GradWarp__N3m, ±Scaled variants]
  Rank 2 : GradWarp + B1 (no N3)
           [MPR__GradWarp__B1_Correction]
  Rank 3 : N3 only (no GradWarp)
           [MT1__N3m, MPR____N3]
  Rank 4 : GradWarp only (no bias correction)
           [MPR__GradWarp, MT1__GradWarp]
  Rank 5 : Raw / minimal preprocessing
           [MPRAGE, IR-FSPGR, bare MT1/MPR]

Secondary sort (within same rank):
  mayo_qc_quality ASC (1=pass preferred; 4 already excluded)
  ImageUID ASC (tiebreaker)

→ Only ONE image per session is selected as the BIDS run-1.
  The next-best (run-2) is recorded for reference only (not copied by default).

Outputs:
    metadata/scan_selection.csv
"""

import pandas as pd
import numpy as np
import os
import re
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exclusions import is_excluded_subject

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT    = r"D:\ADNI_BIDS_project"
SESSION_MAP_CSV = os.path.join(PROJECT_ROOT, "metadata", "session_map.csv")
OUT_SELECTION   = os.path.join(PROJECT_ROOT, "metadata", "scan_selection.csv")

# ── Hard exclusion patterns (series folder name) ───────────────────────────────
EXCLUDE_SERIES = [
    r"(?i).*__Mask$",
    r"(?i).*__Mask__.*",
    r"(?i).*Spatially_Normalized.*",
    r"(?i).*HHP_6_DOF.*",
    r"(?i).*HarP_135.*",
    r"(?i).*EPI_current.*",
]

def is_excluded_series(name: str) -> bool:
    if not name or pd.isna(name):
        return False
    return any(re.match(p, str(name)) for p in EXCLUDE_SERIES)

# ── Feature-based priority (label-agnostic) ────────────────────────────────────
# Detects which corrections are present in the series folder name.
# Both MPR__ and MT1__ naming conventions are handled by the same rules.
PREF_RULES = [
    # Rank 0 — GradWarp + B1 + N3 (full preprocessing chain)
    (r"(?i).*(GradWarp).*(B1).*(N3).*",           0),
    (r"(?i).*(B1).*(Correction).*(N3).*",         0),
    (r"(?i).*(B1_Correction).*(N3).*",            0),

    # Rank 1 — GradWarp + N3 (no B1); covers MT1__GradWarp__N3m too
    (r"(?i).*(GradWarp).*(N3).*",                  1),
    (r"(?i).*(N3).*(GradWarp).*",                  1),

    # Rank 2 — GradWarp + B1 only (no N3)
    (r"(?i).*(GradWarp).*(B1).*",                  2),
    (r"(?i).*(B1).*(GradWarp).*",                  2),

    # Rank 3 — N3 only (no GradWarp)
    (r"(?i).*N3m.*",                               3),
    (r"(?i).*N3.*",                                3),

    # Rank 4 — GradWarp only
    (r"(?i).*GradWarp.*",                          4),

    # Rank 5 — Raw / standard protocol (no specific preprocessing)
    (r"(?i).*MPRAGE.*",                            5),
    (r"(?i).*IR.FSPGR.*",                          5),
    (r"(?i).*IR.SPGR.*",                           5),
    (r"(?i)^MPR",                                  5),
    (r"(?i)^MT1$",                                 5),
]

def pref_rank(series_name: str) -> int:
    if not series_name or pd.isna(series_name):
        return 99
    s = str(series_name)
    for pattern, rank in PREF_RULES:
        if re.match(pattern, s):
            return rank
    return 99

# ── Load ────────────────────────────────────────────────────────────────────────
print("Loading session_map.csv ...")
sess = pd.read_csv(SESSION_MAP_CSV, low_memory=False)
print(f"  -> {len(sess):,} total images, {sess['SubjectID'].nunique():,} subjects")

# Subject-level exclusions
n0 = len(sess)
sess = sess[~sess["SubjectID"].apply(is_excluded_subject)]
print(f"  -> {n0 - len(sess):,} subject-excluded rows removed")

# Hard series exclusions
sess["excl_series"] = sess["sourcedata_series_name"].apply(is_excluded_series)

# Hard QC exclusion: SERIES_QUALITY == 4 (unusable)
sess["excl_qc4"] = (
    pd.to_numeric(sess.get("mayo_qc_quality", pd.Series(dtype=float)), errors="coerce") == 4
)

sess_ok = sess[~sess["excl_series"] & ~sess["excl_qc4"]].copy()
n_excl_series = sess["excl_series"].sum()
n_excl_qc4   = (~sess["excl_series"] & sess["excl_qc4"]).sum()
print(f"  -> {n_excl_series:,} excluded (series type) + {n_excl_qc4:,} excluded (QC=4 unusable)")
print(f"  -> {len(sess_ok):,} images remain for selection")

# ── Apply feature-based preference rank ───────────────────────────────────────
sess_ok["pref_rank"] = sess_ok["sourcedata_series_name"].apply(pref_rank)

# ImageUID numeric tiebreaker
sess_ok["uid_int"] = pd.to_numeric(
    sess_ok["ImageUID"].astype(str).str.replace("I", "", regex=False),
    errors="coerce"
)

# QC sort: NaN -> 99 (unknown/unavailable); 4 already excluded above
sess_ok["sort_qc"] = pd.to_numeric(
    sess_ok.get("mayo_qc_quality", pd.Series(dtype=float)), errors="coerce"
).fillna(99)

# ── Rank and select ───────────────────────────────────────────────────────────
sess_sorted = sess_ok.sort_values(
    ["bids_sub", "bids_ses", "pref_rank", "sort_qc", "uid_int"],
    ascending=[True, True, True, True, True],
    na_position="last"
)

records = []
for (bids_sub, bids_ses), grp in sess_sorted.groupby(["bids_sub", "bids_ses"], sort=False):
    grp = grp.reset_index(drop=True)
    n = len(grp)
    top = grp.iloc[0]

    has_qc = top["sort_qc"] < 99
    reason = "mayo_qc+feature" if has_qc else ("feature_rank" if n > 1 else "only_image")

    rec = {
        "bids_sub":                   bids_sub,
        "bids_ses":                   bids_ses,
        "SubjectID":                  top["SubjectID"],
        "adni_viscode":               top.get("adni_viscode", ""),
        "n_scans_after_exclusion":    n,
        # Selected
        "ImageUID_selected":          top["ImageUID"],
        "sourcedata_series_selected": top.get("sourcedata_series_name", ""),
        "SeriesDescription_selected": top.get("SeriesDescription", ""),
        "nii_source_selected":        top.get("nii_source_path", ""),
        "pref_rank_selected":         int(top["pref_rank"]),
        "mayo_qc_quality_selected":   top["sort_qc"] if has_qc else None,
        "selection_reason":           reason,
        "Manufacturer":               top.get("Manufacturer", ""),
        "ManufacturerModelName":      top.get("ManufacturerModelName", ""),
        "MagneticFieldStrength":      top.get("MagneticFieldStrength", ""),
        "SoftwareVersion":            top.get("SoftwareVersion", ""),
        "SliceThickness":             top.get("SliceThickness", ""),
        # run-1 (same as selected)
        "ImageUID_run1":              top["ImageUID"],
        "sourcedata_series_run1":     top.get("sourcedata_series_name", ""),
        "nii_source_run1":            top.get("nii_source_path", ""),
        # run-2 (next-best, for reference)
        "ImageUID_run2":              grp.iloc[1]["ImageUID"] if n > 1 else "",
        "sourcedata_series_run2":     grp.iloc[1].get("sourcedata_series_name", "") if n > 1 else "",
        "nii_source_run2":            grp.iloc[1].get("nii_source_path", "") if n > 1 else "",
    }
    records.append(rec)

sel = pd.DataFrame(records)
sel.to_csv(OUT_SELECTION, index=False)
print(f"\nSaved: {OUT_SELECTION}")
print(f"Sessions: {len(sel):,} | Subjects: {sel['SubjectID'].nunique():,}")

print("\nFeature-based rank distribution (selected scan):")
rank_labels = {
    0: "GW + B1 + N3",
    1: "GW + N3",
    2: "GW + B1",
    3: "N3 only",
    4: "GW only",
    5: "Raw",
    99: "Unrecognised",
}
for rank, grp_df in sel.groupby("pref_rank_selected"):
    label = rank_labels.get(rank, "?")
    print(f"  Rank {rank} ({label}): {len(grp_df):,} sessions")

print("\nTop selected series folder names:")
print(sel["sourcedata_series_selected"].value_counts().head(12).to_string())

print("\nSelection reason:")
print(sel["selection_reason"].value_counts().to_string())
print(f"\nWith MAYOADIRL QC score: {(sel['mayo_qc_quality_selected'].notna()).sum():,}")
print(f"Multi-scan sessions (run-2 available): {(sel['n_scans_after_exclusion'] > 1).sum():,}")
