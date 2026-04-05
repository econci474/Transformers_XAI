"""
Script 02 — Build QC Selection Table (v3)
==========================================
For each BIDS (subject, session) selects the preferred T1w image and
documents all alternatives in scan_selection.csv.

EXCLUSION (hard filter — these are never T1w inputs for sMRIprep):
  - *__Mask           : brain mask images (binary)
  - Spatially_Normalized__Masked* : pre-normalised + skull-stripped
  - HHP_6_DOF*        : hippocampal-protocol MPRAGE (AC-PC registered)
  - HarP_135*         : hippocampus atlas-registered version
  - EPI_current*      : EPI distortion-corrected (not T1w)

PREFERENCE HIERARCHY (lower rank = preferred):
  Rank 0 : MPR__GradWarp__B1_Correction__N3 (±R, ±Scaled, ±Scaled_2)
  Rank 1 : MPR__GradWarp__N3 (±R, ±Scaled variants)
  Rank 2 : MT1__GradWarp__N3m
  Rank 3 : MT1__N3m
  Rank 4 : MPR__GradWarp__B1_Correction (no N3)  — (GW+B1 but unbiased)
  Rank 5 : MPR__GradWarp (GW only, no bias correction)
  Rank 6 : MPR____N3 (N3 only, no GW)
  Rank 7 : MPR (raw / minimal preprocessing)
  Rank 99: anything else not explicitly recognised

QC sort (applied before preference rank):
  - mayo_qc_quality (1=pass → 4=unusable) if available (ADNI3 only, ~29 images)
  - Otherwise prefer by rank, then by ImageUID ascending (tiebreaker)

COVERAGE SUMMARY (from session_map, after exclusions):
  After excluding 1,193 images, 6,026 remain across 2,175 sessions.
  Top series: MT1__GradWarp__N3m (2,467), MT1__N3m (1,344),
              MPR__GradWarp__B1_Correction__N3 (245+192+16+6=459),
              MPR__GradWarp__N3 (90+80+63+21+8+7=269)

ADNI4 note: MRIQC has 1,421 ADNI4 T1w records but NONE are in the
downloaded sourcedata — dataset covers ADNI1, ADNI1/GO, ADNI2, ADNI3.

Outputs:
    metadata/scan_selection.csv
"""

import pandas as pd
import numpy as np
import os
import re

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT    = r"D:\ADNI_BIDS_project"
SESSION_MAP_CSV = os.path.join(PROJECT_ROOT, "metadata", "session_map.csv")
OUT_SELECTION   = os.path.join(PROJECT_ROOT, "metadata", "scan_selection.csv")

# ── Hard exclusion patterns ────────────────────────────────────────────────────
# These series should NEVER be used as T1w inputs for sMRIprep.
EXCLUDE_PATTERNS = [
    r"(?i).*__Mask$",                              # brain masks
    r"(?i).*__Mask__.*",                           # mask in middle of name
    r"(?i).*Spatially_Normalized.*",               # pre-normalised + skull-stripped
    r"(?i).*HHP_6_DOF.*",                          # hippocampus AC-PC registered
    r"(?i).*HarP_135.*",                           # hippocampus atlas-registered
    r"(?i).*EPI_current.*",                        # EPI-corrected, not T1w
]

def is_excluded(series_name: str) -> bool:
    if not series_name or pd.isna(series_name):
        return False
    return any(re.match(p, str(series_name)) for p in EXCLUDE_PATTERNS)

# ── Preference ranking ─────────────────────────────────────────────────────────
# Matches on the sourcedata folder name (actual downloaded series).
# Lower rank = more preferred for sMRIprep T1w input.
PREF_RULES = [
    # Rank 0 — GradWarp + B1 correction + N3 bias correction (best chain)
    # Covers: MPR__GradWarp__B1_Correction__N3, MPR-R variant, ±Scaled, ±Scaled_2
    (r"(?i)^MPR-?R?__GradWarp__B1_Correction__N3",   0),

    # Rank 1 — GradWarp + N3 (no B1)
    # Covers: MPR__GradWarp__N3, MPR-R, MPR__GradWarp__N3__Scaled, etc.
    (r"(?i)^MPR-?R?__GradWarp__N3",                   1),
    (r"(?i)^MPR-?R?____N3",                            1),  # no GW, N3 only (MPR____N3)

    # Rank 2 — MT1 with GradWarp + N3m (ADNI standard pipeline)
    (r"(?i)^MT1__GradWarp__N3m",                       2),

    # Rank 3 — MT1 N3m only (no GradWarp)
    (r"(?i)^MT1__N3m",                                 3),

    # Rank 4 — GradWarp + B1, no N3
    (r"(?i)^MPR-?R?__GradWarp__B1_Correction$",        4),

    # Rank 5 — GradWarp only
    (r"(?i)^MPR-?R?__GradWarp$",                       5),

    # Rank 6 — raw MPR (no preprocessing)
    (r"(?i)^MPR",                                       6),
    (r"(?i)^MT1$",                                      6),

    # Rank 7 — Accelerated/standard MPRAGE (ADNI3/4 raw)
    (r"(?i).*MPRAGE.*",                                 7),
    (r"(?i).*IR.FSPGR.*",                               7),
]

def pref_rank(series_name: str) -> int:
    if not series_name or pd.isna(series_name):
        return 99
    for pattern, rank in PREF_RULES:
        if re.match(pattern, str(series_name)):
            return rank
    return 99

# ── Load and filter session map ────────────────────────────────────────────────
print("Loading session_map.csv ...")
sess = pd.read_csv(SESSION_MAP_CSV, low_memory=False)
print(f"  -> {len(sess):,} total images")

# Apply hard exclusions
sess["excluded"] = sess["sourcedata_series_name"].apply(is_excluded)
n_excl = sess["excluded"].sum()
sess_ok = sess[~sess["excluded"]].copy()
print(f"  -> {n_excl:,} images hard-excluded (masks/normalised/hippocampus-registered)")
print(f"  -> {len(sess_ok):,} images remaining for selection")

# Apply preference rank
sess_ok["pref_rank"] = sess_ok["sourcedata_series_name"].apply(pref_rank)

# ImageUID numeric for tiebreaker
sess_ok["uid_int"] = pd.to_numeric(
    sess_ok["ImageUID"].astype(str).str.replace("I", "", regex=False),
    errors="coerce"
)

# QC sort (MAYOADIRL: 1=pass, 4=unusable; NaN=unavailable → sort as 99)
sess_ok["sort_qc"] = pd.to_numeric(
    sess_ok.get("mayo_qc_quality", pd.Series(dtype=float)),
    errors="coerce"
).fillna(99)

# Sort within groups
sess_sorted = sess_ok.sort_values(
    ["bids_sub", "bids_ses", "sort_qc", "pref_rank", "uid_int"],
    ascending=[True, True, True, True, True],
    na_position="last"
)

# ── Build selection table ─────────────────────────────────────────────────────
records = []
for (bids_sub, bids_ses), grp in sess_sorted.groupby(["bids_sub", "bids_ses"], sort=False):
    grp = grp.reset_index(drop=True)
    n = len(grp)
    top = grp.iloc[0]

    has_qc = top["sort_qc"] < 99
    reason = "mayo_qc+pref" if has_qc else ("pref_rank" if n > 1 else "only_image")

    rec = {
        "bids_sub":                   bids_sub,
        "bids_ses":                   bids_ses,
        "SubjectID":                  top["SubjectID"],
        "adni_viscode":               top.get("adni_viscode", ""),
        "n_scans_after_exclusion":    n,
        # Selected (run-1)
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
        # run-2 (second-ranked, if present after exclusion)
        "ImageUID_run2":              grp.iloc[1]["ImageUID"] if n > 1 else "",
        "sourcedata_series_run2":     grp.iloc[1].get("sourcedata_series_name", "") if n > 1 else "",
        "nii_source_run2":            grp.iloc[1].get("nii_source_path", "") if n > 1 else "",
    }
    records.append(rec)

sel = pd.DataFrame(records)
sel.to_csv(OUT_SELECTION, index=False)
print(f"\nSaved: {OUT_SELECTION}")
print(f"Sessions: {len(sel):,} | Subjects: {sel['SubjectID'].nunique():,}")

print("\nPreference rank distribution (selected scan):")
rank_labels = {
    0: "MPR GW+B1+N3", 1: "MPR GW+N3 / MPR N3",
    2: "MT1 GW+N3m",   3: "MT1 N3m",
    4: "MPR GW+B1",    5: "MPR GW only",
    6: "Raw MPR/MT1",  7: "MPRAGE/IR-FSPGR", 99: "Unrecognised"
}
for rank, grp_df in sel.groupby("pref_rank_selected"):
    label = rank_labels.get(rank, "?")
    print(f"  Rank {rank} ({label}): {len(grp_df):,} sessions")

print("\nTop selected series:")
print(sel["sourcedata_series_selected"].value_counts().head(10).to_string())

print("\nSelection reason:")
print(sel["selection_reason"].value_counts().to_string())

print("\nWith MAYOADIRL QC available:", (sel["mayo_qc_quality_selected"].notna()).sum())
print("Multi-scan sessions (run-2 available):", (sel["n_scans_after_exclusion"] > 1).sum())
