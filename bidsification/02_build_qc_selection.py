"""
Script 02 — Build QC Selection Table
======================================
For each (subject, session), rank available T1w scans by:
  1. QC quality score from MAYOADIRL_MRI_IMAGEQC (series_quality: 1=best → 4=fail)
  2. Series preference:
       N3m + GradWarp + B1-corrected  → score 0 (top)
       GradWarp + N3m                 → score 1
       N3m only                       → score 2
       Other accepted                 → score 3
  3. If no QC entry exists, select first by ImageUID (numeric order).

Outputs:
    metadata/scan_selection.csv  — Per-session scan selection table
    
Columns:
    bids_sub, bids_ses, adni_viscode,
    ImageUID_selected, SeriesDescription_selected, selection_reason,
    ImageUID_run1, SeriesDescription_run1,
    ImageUID_run2, SeriesDescription_run2,   (if multiple scans exist)
    qc_quality_selected (1-4 or NaN if no QC),
    series_pref_rank_selected (0-3)
"""

import pandas as pd
import numpy as np
import os
import re

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT = r"D:\ADNI_BIDS_project"
SESSION_MAP_CSV = os.path.join(PROJECT_ROOT, "metadata", "session_map.csv")
MAYOADIRL_CSV = os.path.join(
    PROJECT_ROOT, "sourcedata", "clinical", "MAYOADIRL_MRI_IMAGEQC_12_08_15.csv"
)
OUT_SELECTION = os.path.join(PROJECT_ROOT, "metadata", "scan_selection.csv")

# ── Series preference ranking ──────────────────────────────────────────────────
# Lower rank = more preferred
SERIES_PREF = {
    # GradWarp + B1 + N3 (best) — ADNI-specific preprocessing labels
    r"(?i).*GradWarp.*B1.*N3.*":       0,
    r"(?i).*B1.*GradWarp.*N3.*":       0,
    r"(?i).*N3m.*GradWarp.*":          0,
    r"(?i).*GradWarp.*N3m.*":          0,
    # GradWarp + N3 only
    r"(?i).*GradWarp.*N3.*":           1,
    r"(?i).*N3.*GradWarp.*":           1,
    # N3 corrected only
    r"(?i).*MT1__N3m.*":               2,
    r"(?i).*N3m.*":                    2,
    r"(?i).*N3.*":                     2,
    # GradWarp only
    r"(?i).*GradWarp.*":               3,
    # Accelerated MPRAGE (ADNI3 standard)
    r"(?i).*Accelerated.*MPRAGE.*":    4,
    r"(?i).*MPRAGE.*GRAPPA.*":         4,
    # Standard MPRAGE
    r"(?i).*MPRAGE.*":                 5,
    r"(?i).*MP.RAGE.*":                5,
    r"(?i).*IR.FSPGR.*":               5,
}

def series_pref_rank(desc: str) -> int:
    """Return preference rank for a series description (lower = better)."""
    if pd.isna(desc):
        return 99
    for pattern, rank in SERIES_PREF.items():
        if re.match(pattern, str(desc)):
            return rank
    return 10  # unknown but accepted

# ── Load data ──────────────────────────────────────────────────────────────────
print("Loading session map ...")
sess = pd.read_csv(SESSION_MAP_CSV, low_memory=False)
print(f"  → {len(sess):,} rows, {sess['SubjectID'].nunique():,} subjects")

print("Loading MAYOADIRL QC table ...")
qc = pd.read_csv(MAYOADIRL_CSV, low_memory=False)
# series_quality: 1=pass, 2=acceptable, 3=fail (but selected), 4=unusable, -1=N/A
# loni_image maps to ImageUID
qc_cols = ["PTID", "loni_image", "series_quality", "series_selected",
           "series_type", "series_description", "series_comments",
           "protocol_status", "protocol_comments"]
qc_cols = [c for c in qc_cols if c in qc.columns]
qc_t1 = qc[qc.get("series_type", qc.get("series_type", pd.Series())).isin(["MT1", "MT1_N3m"]) 
            if "series_type" in qc.columns else pd.Series(dtype=bool)].copy()
# If filtering didn't work (series_type mismatch), keep all and merge on image id
if len(qc_t1) == 0:
    qc_t1 = qc.copy()
print(f"  → QC rows available for merge: {len(qc_t1):,}")

# Normalise loni_image → ImageUID
qc_t1["ImageUID_int"] = pd.to_numeric(qc_t1.get("loni_image", pd.Series()), errors="coerce")
sess["ImageUID_int"] = pd.to_numeric(sess["ImageUID"], errors="coerce")

# Merge QC into session map
sess_qc = sess.merge(
    qc_t1[["ImageUID_int", "series_quality", "series_selected", "series_comments"]].drop_duplicates("ImageUID_int"),
    on="ImageUID_int", how="left"
)

# Map series_quality: treat -1 and NaN as 99 (no QC available, neutral)
sess_qc["qc_quality"] = pd.to_numeric(sess_qc.get("series_quality", pd.Series()), errors="coerce")
sess_qc["qc_quality"] = sess_qc["qc_quality"].replace(-1, np.nan)

# Add series preference rank
sess_qc["pref_rank"] = sess_qc["SeriesDescription"].apply(series_pref_rank)

# ── Rank within each (bids_sub, bids_ses) group ────────────────────────────────
# Priority: qc_quality ASC (lower=better), then pref_rank ASC, then ImageUID_int ASC
sess_qc["sort_qc"] = sess_qc["qc_quality"].fillna(99)
sess_qc_sorted = sess_qc.sort_values(
    ["bids_sub", "bids_ses", "sort_qc", "pref_rank", "ImageUID_int"],
    ascending=[True, True, True, True, True],
    na_position="last"
)

records = []
for (bids_sub, bids_ses), grp in sess_qc_sorted.groupby(["bids_sub", "bids_ses"], sort=False):
    grp = grp.reset_index(drop=True)
    n = len(grp)
    
    # Determine selection reason
    has_qc = grp["sort_qc"].min() < 99
    
    rec = {
        "bids_sub":                     bids_sub,
        "bids_ses":                     bids_ses,
        "SubjectID":                    grp.loc[0, "SubjectID"],
        "adni_viscode":                 grp.loc[0, "adni_viscode"] if "adni_viscode" in grp.columns else "",
        "n_scans":                      n,
        # Preferred / selected scan
        "ImageUID_selected":            grp.loc[0, "ImageUID"],
        "SeriesDescription_selected":   grp.loc[0, "SeriesDescription"],
        "nii_source_selected":          grp.loc[0, "nii_source_path"] if "nii_source_path" in grp.columns else "",
        "qc_quality_selected":          grp.loc[0, "qc_quality"] if has_qc else None,
        "pref_rank_selected":           grp.loc[0, "pref_rank"],
        "selection_reason":             "qc_score+series_pref" if has_qc else "first_by_imageuid_no_qc",
        # run-1 (same as selected)
        "ImageUID_run1":                grp.loc[0, "ImageUID"],
        "SeriesDescription_run1":       grp.loc[0, "SeriesDescription"],
        "nii_source_run1":              grp.loc[0, "nii_source_path"] if "nii_source_path" in grp.columns else "",
        # run-2 (only if >1 scan)
        "ImageUID_run2":                grp.loc[1, "ImageUID"] if n > 1 else "",
        "SeriesDescription_run2":       grp.loc[1, "SeriesDescription"] if n > 1 else "",
        "nii_source_run2":              grp.loc[1, "nii_source_path"] if (n > 1 and "nii_source_path" in grp.columns) else "",
    }
    
    # Add QC comment if available
    comments = grp["series_comments"].dropna().tolist()
    rec["qc_comments"] = "; ".join(str(c) for c in comments[:3]) if comments else ""
    
    records.append(rec)

selection_df = pd.DataFrame(records)
selection_df.to_csv(OUT_SELECTION, index=False)
print(f"\nSaved: {OUT_SELECTION}")
print(f"Total sessions: {len(selection_df):,}")
print(f"  With QC data: {(selection_df['qc_quality_selected'].notna()).sum():,}")
print(f"  No QC found:  {(selection_df['qc_quality_selected'].isna()).sum():,}")
print(f"  Multi-scan sessions (run-1 + run-2): {(selection_df['n_scans'] > 1).sum():,}")

# Summary of selection reasons
print("\nSelection reason breakdown:")
print(selection_df["selection_reason"].value_counts().to_string())
