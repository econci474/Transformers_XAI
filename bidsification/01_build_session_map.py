"""
Script 01 — Build Session Map
==============================
Parses MRIQUALITY.csv (which contains VISCODE2 i.e., ADNI visit codes)
and image_series_mapping.csv to create a comprehensive session map.

Outputs:
    metadata/session_map.csv        — Full session-level inventory (one row per image)
    metadata/ses_to_visit_code.csv  — Correspondence table: bids_ses ↔ ADNI VISCODE2
"""

import pandas as pd
import os
import re

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT = r"D:\ADNI_BIDS_project"
MRIQUALITY_CSV = os.path.join(PROJECT_ROOT, "sourcedata", "clinical", "MRIQUALITY.csv")
IMAGE_MAPPING_CSV = os.path.join(PROJECT_ROOT, "metadata", "image_series_mapping.csv")
SOURCEDATA_DIR = os.path.join(PROJECT_ROOT, "sourcedata", "ADNI")
OUT_SESSION_MAP = os.path.join(PROJECT_ROOT, "metadata", "session_map.csv")
OUT_SES_VISIT = os.path.join(PROJECT_ROOT, "metadata", "ses_to_visit_code.csv")

# ── T1w series type codes to keep (from ADNI documentation) ───────────────────
# These correspond to the series_type values in MAYOADIRL / MRIQUALITY
T1W_SERIES_TYPES = {"T1w", "MT1"}          # MRIQUALITY uses 'T1w'
T1W_SERIES_DESC_KEYWORDS = [
    "MPRAGE", "MP-RAGE", "MPRAGE GRAPPA", "IR-FSPGR",
    "Accelerated Sagittal MPRAGE", "Sagittal MPRAGE",
    "MT1", "MT1__N3m", "MT1__GradWarp__N3m",
]

def normalize_subject_id(ptid: str) -> str:
    """Convert e.g. '002_S_0295' → '002S0295' (remove underscores)."""
    return re.sub(r'[_\-]', '', ptid)

def normalize_study_date(date_str: str) -> str:
    """Convert e.g. '2011-06-02' or '2011-06-02 09:30:00' → '20110602'."""
    if pd.isna(date_str):
        return ""
    return str(date_str)[:10].replace("-", "")

print("Loading MRIQUALITY.csv ...")
mri_quality = pd.read_csv(MRIQUALITY_CSV, low_memory=False)
print(f"  → {len(mri_quality):,} rows, columns: {list(mri_quality.columns)}")

# Filter to T1w only
t1w_mask = mri_quality["SeriesType"].isin(T1W_SERIES_TYPES)
mri_t1w = mri_quality[t1w_mask].copy()
print(f"  → T1w rows: {len(mri_t1w):,}")

print("Loading image_series_mapping.csv ...")
img_map = pd.read_csv(IMAGE_MAPPING_CSV, low_memory=False)
print(f"  → {len(img_map):,} rows, columns: {list(img_map.columns)}")
# Normalise column names
img_map.columns = [c.strip() for c in img_map.columns]

# ── Rename key columns in mri_quality for clarity ─────────────────────────────
mri_t1w = mri_t1w.rename(columns={
    "image_id":           "ImageUID",
    "ParticipantID":      "SubjectID",
    "VISCODE2":           "adni_viscode",
    "StudyDate":          "StudyDate",
    "SeriesType":         "SeriesType",
    "SeriesDescription":  "SeriesDescription",
    "ScannerManufacturer":"Manufacturer",
    "ScannerModel":       "ManufacturerModelName",
    "MagneticFieldStrength": "MagneticFieldStrength",
    "StudyInstanceUID":   "StudyInstanceUID",
    "SeriesInstanceUID":  "SeriesInstanceUID",
    "LONIStudy":          "LONIStudy",
    "LONISeries":         "LONISeries",
    "LONIImage":          "LONIImage",
}).copy()

# ── Merge with image_series_mapping (adds SeriesUID and StudyUID) ──────────────
# image_series_mapping columns: ImageUID, SeriesUID, SubjectID, StudyUID
# LONIImage in mri_quality = ImageUID

mri_t1w["ImageUID_int"] = pd.to_numeric(mri_t1w["ImageUID"], errors="coerce")
img_map["ImageUID_int"] = pd.to_numeric(img_map.get("ImageUID", img_map.iloc[:, 0]), errors="coerce")

# Try merge on ImageUID
try:
    merged = mri_t1w.merge(
        img_map[["ImageUID_int", "SeriesUID", "StudyUID"]].drop_duplicates("ImageUID_int"),
        on="ImageUID_int", how="left"
    )
except KeyError:
    print("  WARNING: Could not merge SeriesUID/StudyUID — check image_series_mapping columns.")
    merged = mri_t1w.copy()
    merged["SeriesUID"] = ""
    merged["StudyUID"] = ""

# ── Add BIDS labels ────────────────────────────────────────────────────────────
merged["bids_sub"] = merged["SubjectID"].apply(normalize_subject_id)
merged["bids_ses"] = merged["StudyDate"].apply(normalize_study_date)

# ── Build NIfTI source path (in sourcedata/ADNI/) ─────────────────────────────
def find_nii_path(row):
    """
    Search for the NIfTI file corresponding to this ImageUID in sourcedata/ADNI/.
    Expected structure: sourcedata/ADNI/<SubjectID>/MT1__N3m/<date_dir>/<ImageUID>/<file.nii>
    Returns the path to the .nii file or empty string if not found.
    """
    sub_dir = os.path.join(SOURCEDATA_DIR, str(row["SubjectID"]))
    if not os.path.isdir(sub_dir):
        return ""
    image_uid_str = str(int(row["ImageUID_int"])) if pd.notna(row["ImageUID_int"]) else ""
    if not image_uid_str:
        return ""
    # Walk series subdirectories
    for series_name in os.listdir(sub_dir):
        series_dir = os.path.join(sub_dir, series_name)
        if not os.path.isdir(series_dir):
            continue
        for date_dir in os.listdir(series_dir):
            img_dir = os.path.join(series_dir, date_dir, image_uid_str)
            if os.path.isdir(img_dir):
                niis = [f for f in os.listdir(img_dir) if f.endswith(".nii")]
                if niis:
                    return os.path.join(img_dir, niis[0])
    return ""

print("Locating NIfTI source files (this may take a few minutes) ...")
merged["nii_source_path"] = merged.apply(find_nii_path, axis=1)
found = (merged["nii_source_path"] != "").sum()
print(f"  → {found}/{len(merged)} NIfTI files located in sourcedata/ADNI/")

# ── Save session map ───────────────────────────────────────────────────────────
cols_out = [
    "SubjectID", "ImageUID", "bids_sub", "bids_ses", "adni_viscode",
    "StudyDate", "SeriesType", "SeriesDescription",
    "Manufacturer", "ManufacturerModelName", "MagneticFieldStrength",
    "StudyInstanceUID", "SeriesInstanceUID",
    "SeriesUID", "StudyUID",
    "LONIStudy", "LONISeries", "LONIImage",
    "nii_source_path",
]
cols_out = [c for c in cols_out if c in merged.columns]
merged[cols_out].to_csv(OUT_SESSION_MAP, index=False)
print(f"Saved: {OUT_SESSION_MAP}")

# ── Save visit-code correspondence table ───────────────────────────────────────
ses_visit = (
    merged[["bids_sub", "bids_ses", "adni_viscode", "StudyDate", "SubjectID"]]
    .drop_duplicates(subset=["bids_sub", "bids_ses"])
    .sort_values(["bids_sub", "bids_ses"])
)
ses_visit.to_csv(OUT_SES_VISIT, index=False)
print(f"Saved: {OUT_SES_VISIT}")
print(f"\nDone. {len(merged):,} T1w image records across {merged['SubjectID'].nunique():,} subjects.")
