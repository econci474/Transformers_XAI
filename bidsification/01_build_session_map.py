"""
Script 01 — Build Session Map (REVISED)
========================================
Builds the master session map by joining:

  image_series_mapping.csv (7,219 downloaded images)
    → primary source: ImageUID (I-prefixed), SeriesUID (S-prefixed), SubjectID, StudyUID

  MRIQC_15Feb2026.csv (most up-to-date MRIQUALITY-style table, all ADNI phases)
    → join key: LONISeries == int(SeriesUID[1:])   [strip 'S' prefix]
    → provides: VISCODE2, StudyDate, SeriesDescription, scanner metadata

  MAYOADIRL_MRI_IMAGEQC_12_08_15.csv  (ADNI1/GO/2 QC)
  MAYOADIRL_MRI_QUALITY_ADNI3.csv     (ADNI3 QC, original)
  MAYOADIRL_MRI_QUALITY_ADNI3_15Feb2026.csv  (ADNI3 QC, updated Feb 2026)
    → join key: PTID + series_date (YYYYMMDD) — at session level
    → provides: series_quality (1=pass,2=acceptable,3=fail,4=unusable), series_selected

Join strategy:
  1. img_map is the spine (7,219 rows = all downloaded NIfTIs)
  2. Merge MRIQC on LONISeries (series-level, 1:many reduced to best match)
  3. Merge MAYOADIRL on SubjectID + StudyDate (session-level QC)
  4. Locate NIfTI paths by scanning sourcedata I-dirs (fast pre-index)

Outputs:
    metadata/session_map.csv        — full session/image inventory
    metadata/ses_to_visit_code.csv  — bids_ses <-> ADNI VISCODE2 correspondence
"""

import pandas as pd
import numpy as np
import os
import re
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exclusions import is_excluded_subject

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT   = r"D:\ADNI_BIDS_project"
IMG_MAP_CSV    = os.path.join(PROJECT_ROOT, "metadata", "image_series_mapping.csv")
SOURCEDATA_DIR = os.path.join(PROJECT_ROOT, "sourcedata", "ADNI")
CLINICAL_DIR   = os.path.join(PROJECT_ROOT, "sourcedata", "clinical")

# QC / metadata files (use most up-to-date versions)
MRIQC_CSV      = os.path.join(CLINICAL_DIR, "MRIQC_15Feb2026.csv")   # preferred; fallback below
MRIQUALITY_CSV = os.path.join(CLINICAL_DIR, "MRIQUALITY.csv")        # fallback if MRIQC not present

MAYO1_CSV      = os.path.join(CLINICAL_DIR, "MAYOADIRL_MRI_IMAGEQC_12_08_15.csv")
MAYO3_CSV      = os.path.join(CLINICAL_DIR, "MAYOADIRL_MRI_QUALITY_ADNI3_15Feb2026.csv")
MAYO3_ORIG_CSV = os.path.join(CLINICAL_DIR, "MAYOADIRL_MRI_QUALITY_ADNI3.csv")  # fallback

OUT_SESSION_MAP = os.path.join(PROJECT_ROOT, "metadata", "session_map.csv")
OUT_SES_VISIT   = os.path.join(PROJECT_ROOT, "metadata", "ses_to_visit_code.csv")

# ── Helper functions ───────────────────────────────────────────────────────────
def normalize_subject_id(ptid: str) -> str:
    """'002_S_0295' -> '002S0295' (BIDS sub label, underscores removed)."""
    return re.sub(r"[_\-]", "", str(ptid))

def normalize_date(date_str) -> str:
    """Normalize ADNI study dates to YYYYMMDD string.
    Handles ISO ('2011-06-02'), MRIQC float (20110602.0), and 8-digit strings."""
    if pd.isna(date_str):
        return ""
    s = str(date_str).strip()
    if s.endswith('.0'):          # MRIQC stores dates as float, e.g. 20110602.0
        s = s[:-2]
    if len(s) >= 10 and '-' in s:  # ISO format: 2011-06-02
        return s[:10].replace('-', '')
    return s[:8]                  # already YYYYMMDD

# ── 1. Load image_series_mapping (spine) ──────────────────────────────────────
print("Loading image_series_mapping.csv ...")
img_map = pd.read_csv(IMG_MAP_CSV, low_memory=False)
img_map.columns = [c.strip() for c in img_map.columns]
# Numeric series ID for join: e.g. 'S17197' -> 17197
img_map["series_int"] = pd.to_numeric(
    img_map["SeriesUID"].astype(str).str.replace("S", "", regex=False),
    errors="coerce"
)
# Numeric study ID
img_map["study_int"] = pd.to_numeric(img_map["StudyUID"], errors="coerce")
print(f"  -> {len(img_map):,} rows, {img_map['SubjectID'].nunique():,} subjects")

# Apply subject-level exclusions
n_before = len(img_map)
img_map = img_map[~img_map["SubjectID"].apply(is_excluded_subject)]
print(f"  -> {n_before - len(img_map):,} rows excluded (excluded subject list); {len(img_map):,} remaining")

# ── 2. Load MRIQC (most up-to-date metadata table) ────────────────────────────
mriqc_path = MRIQC_CSV if os.path.isfile(MRIQC_CSV) else MRIQUALITY_CSV
print(f"Loading {os.path.basename(mriqc_path)} ...")
mriqc = pd.read_csv(mriqc_path, low_memory=False)
mriqc_t1 = mriqc[mriqc["SeriesType"] == "T1w"].copy()
print(f"  -> {len(mriqc_t1):,} T1w rows across {mriqc_t1['ParticipantID'].nunique():,} subjects")

# Numeric join key on LONISeries
mriqc_t1["series_int"] = pd.to_numeric(mriqc_t1["LONISeries"], errors="coerce")

# Deduplicate: one MRIQC row per series_int (keep by image_id ascending for reproducibility)
mriqc_t1_dedup = (
    mriqc_t1.sort_values("image_id")
    .drop_duplicates(subset="series_int", keep="first")
)
print(f"  -> {len(mriqc_t1_dedup):,} unique T1w series in MRIQC")

# Merge img_map -> MRIQC on series_int
merged = img_map.merge(
    mriqc_t1_dedup.rename(columns={
        "image_id":            "mriqc_image_id",
        "ParticipantID":       "SubjectID_mriqc",
        "VISCODE2":            "adni_viscode",
        "StudyDate":           "StudyDate",
        "SeriesDescription":   "SeriesDescription",
        "ScannerManufacturer": "Manufacturer",
        "ScannerModel":        "ManufacturerModelName",
        "MagneticFieldStrength": "MagneticFieldStrength",
        "SoftwareVersion":     "SoftwareVersion",
        "Acceleration":        "Acceleration",
        "SliceThickness":      "SliceThickness",
        "AcquisitionType":     "AcquisitionType",
        "StudyInstanceUID":    "StudyInstanceUID",
        "SeriesInstanceUID":   "SeriesInstanceUID",
        "LONIStudy":           "LONIStudy",
        "LONISeries":          "LONISeries",
        "LONIImage":           "LONIImage_mriqc",
    })[
        ["series_int", "mriqc_image_id", "SubjectID_mriqc", "adni_viscode", "StudyDate",
         "SeriesDescription", "Manufacturer", "ManufacturerModelName", "MagneticFieldStrength",
         "SoftwareVersion", "Acceleration", "SliceThickness", "AcquisitionType",
         "StudyInstanceUID", "SeriesInstanceUID", "LONIStudy", "LONISeries", "LONIImage_mriqc"]
    ].drop_duplicates("series_int"),
    on="series_int", how="left"
)
n_matched = merged["adni_viscode"].notna().sum()
print(f"  -> MRIQC join: {n_matched:,}/{len(merged):,} images matched to MRIQC metadata")

# ── 3. Load and concatenate MAYOADIRL QC tables ────────────────────────────────
print("Loading MAYOADIRL QC tables ...")

def load_mayo(path, phase_label):
    """Load a MAYOADIRL QC file, normalise column names to lowercase."""
    df = pd.read_csv(path, low_memory=False)
    df.columns = [c.lower().strip() for c in df.columns]
    df["mayo_phase"] = phase_label
    return df

mayo1 = load_mayo(MAYO1_CSV, "ADNI1/GO/2")
mayo3_path = MAYO3_CSV if os.path.isfile(MAYO3_CSV) else MAYO3_ORIG_CSV
mayo3 = load_mayo(mayo3_path, "ADNI3")
mayo_all = pd.concat([mayo1, mayo3], ignore_index=True)
print(f"  -> Combined MAYOADIRL: {len(mayo_all):,} rows")

# Normalise date column (may be 'series_date' in both)
mayo_all["series_date_norm"] = mayo_all["series_date"].apply(normalize_date)

# Filter to T1w only in MAYOADIRL (series_type codes for T1w)
T1_TYPES = {"T1", "MT1", "MT1w", "T1w"}
if "series_type" in mayo_all.columns:
    mayo_t1 = mayo_all[mayo_all["series_type"].str.upper().isin({t.upper() for t in T1_TYPES})].copy()
    if len(mayo_t1) == 0:
        mayo_t1 = mayo_all  # fallback: no filter if types don't match
else:
    mayo_t1 = mayo_all
print(f"  -> MAYOADIRL T1w rows: {len(mayo_t1):,}")

# Build session-level QC: best quality score per (PTID + series_date)
# series_quality: 1=pass, 2=acceptable, 3=fail, 4=unusable, -1=N/A -> treat as 99
mayo_t1["qc_num"] = pd.to_numeric(mayo_t1["series_quality"], errors="coerce").replace(-1, np.nan)
mayo_t1_sorted = mayo_t1.sort_values("qc_num", ascending=True, na_position="last")
mayo_sess = mayo_t1_sorted.drop_duplicates(subset=["ptid", "series_date_norm"], keep="first")

# Keep only relevant columns
mayo_keep = ["ptid", "series_date_norm", "qc_num", "series_selected", "series_comments",
             "study_overallpass", "series_description", "loni_study", "loni_series", "loni_image",
             "field_strength", "mayo_phase"]
mayo_keep = [c for c in mayo_keep if c in mayo_sess.columns]
mayo_sess = mayo_sess[mayo_keep].rename(columns={
    "ptid":               "SubjectID",
    "series_date_norm":   "StudyDate_norm",
    "qc_num":             "mayo_qc_quality",
    "series_selected":    "mayo_series_selected",
    "series_comments":    "mayo_series_comments",
    "study_overallpass":  "mayo_study_overallpass",
    "loni_study":         "mayo_loni_study",
    "loni_series":        "mayo_loni_series",
    "loni_image":         "mayo_loni_image",
    "field_strength":     "mayo_field_strength",
    "mayo_phase":         "ADNI_phase",
})

# ── 4. Merge MAYOADIRL QC into main table ─────────────────────────────────────
# Need StudyDate_norm on merged first
merged["StudyDate_norm"] = merged["StudyDate"].apply(normalize_date)

merged = merged.merge(
    mayo_sess,
    left_on=["SubjectID", "StudyDate_norm"],
    right_on=["SubjectID", "StudyDate_norm"],
    how="left"
)
n_mayo = merged["mayo_qc_quality"].notna().sum()
print(f"  -> MAYOADIRL join: {n_mayo:,}/{len(merged):,} images matched to QC scores")

# ── 5. MRI3META fallback for images without MRIQC metadata ────────────────────────
# 536 ADNI1 images have old SeriesUIDs not catalogued in MRIQC.
# MRI3META (session-level MRI protocol CRF) provides PTID+VISCODE2+EXAMDATE
# for every 3T MRI session. We recover scan dates from the sourcedata dir structure
# and join via SubjectID + EXAMDATE to fill missing visit info.

_missing_mask = merged["adni_viscode"].isna()
if _missing_mask.sum() > 0:
    print(f"Recovering visit info via MRI3META for {_missing_mask.sum()} images without MRIQC match ...")

    # Load MRI3META (prefer Feb 2026 version)
    _mri3_path = os.path.join(CLINICAL_DIR, "MRI3META_15Feb2026.csv")
    if not os.path.isfile(_mri3_path):
        _mri3_path = os.path.join(CLINICAL_DIR, "MRI3META.csv")
    _mri3 = pd.read_csv(_mri3_path, low_memory=False)
    _mri3["examdate_str"] = _mri3["EXAMDATE"].astype(str).str[:10]
    _mri3_lookup = _mri3[["PTID","VISCODE2","examdate_str"]].drop_duplicates(["PTID","examdate_str"])

    # Recover scan dates from sourcedata timestamp directories
    _missing_idx = merged[_missing_mask].index
    _scan_dates = {}
    for _, row in merged.loc[_missing_idx].iterrows():
        sub_dir = os.path.join(SOURCEDATA_DIR, str(row["SubjectID"]))
        uid_dir = str(row["ImageUID"]) if str(row["ImageUID"]).startswith("I") else "I" + str(row["ImageUID"])
        if not os.path.isdir(sub_dir):
            continue
        for series in os.listdir(sub_dir):
            s_path = os.path.join(sub_dir, series)
            if not os.path.isdir(s_path): continue
            for date_folder in os.listdir(s_path):
                d_path = os.path.join(s_path, date_folder)
                if os.path.isdir(d_path) and uid_dir in os.listdir(d_path):
                    _scan_dates[row["ImageUID"]] = date_folder[:10]
                    break
            if row["ImageUID"] in _scan_dates: break

    merged.loc[_missing_idx, "StudyDate"] = merged.loc[_missing_idx, "ImageUID"].map(_scan_dates)

    # Join MRI3META on SubjectID + scan date
    _tmp = merged.loc[_missing_idx, ["SubjectID","StudyDate"]].copy()
    _tmp["scan_date_str"] = _tmp["StudyDate"].astype(str).str[:10]
    _tmp = _tmp.merge(
        _mri3_lookup.rename(columns={"PTID":"SubjectID","examdate_str":"scan_date_str"}),
        on=["SubjectID","scan_date_str"], how="left"
    )
    _tmp.index = _missing_idx
    merged.loc[_missing_idx, "adni_viscode"] = _tmp["VISCODE2"].values

    _now_matched = merged.loc[_missing_idx, "adni_viscode"].notna().sum()
    print(f"  -> {_now_matched}/{_missing_mask.sum()} images recovered via MRI3META")

# ── 6. Add BIDS labels ─────────────────────────────────────────────────────────
merged["bids_sub"] = merged["SubjectID"].apply(normalize_subject_id)

# BIDS session label: prefer VISCODE2 (de-identified, ADNI-standard);
# fall back to ses-d<YYYYMMDD> for images without MRIQC match.
# ses-d prefix marks date-based fallbacks so they are identifiable.
def make_bids_ses(row) -> str:
    viscode = row.get("adni_viscode")
    if viscode and str(viscode).strip() and str(viscode).strip() != "nan":
        return str(viscode).strip().lower()   # e.g. 'bl', 'm12', 'm24'
    date = row.get("StudyDate")
    norm = normalize_date(date)
    return f"d{norm}" if norm else "unknown"   # 'd' prefix = date-derived fallback

merged["bids_ses"] = merged.apply(make_bids_ses, axis=1)

# ── 6. Build sourcedata NIfTI index ───────────────────────────────────────────
print("Building sourcedata NIfTI index ...")
nii_index = {}        # {uid_int: nii_path}
series_name_index = {}  # {uid_int: series_folder_name}

for sub_name in os.listdir(SOURCEDATA_DIR):
    sub_dir = os.path.join(SOURCEDATA_DIR, sub_name)
    if not os.path.isdir(sub_dir):
        continue
    for series_name in os.listdir(sub_dir):
        series_dir = os.path.join(sub_dir, series_name)
        if not os.path.isdir(series_dir):
            continue
        for date_dir in os.listdir(series_dir):
            date_path = os.path.join(series_dir, date_dir)
            if not os.path.isdir(date_path):
                continue
            for img_dir_name in os.listdir(date_path):
                if not img_dir_name.startswith("I"):
                    continue
                img_dir = os.path.join(date_path, img_dir_name)
                if not os.path.isdir(img_dir):
                    continue
                try:
                    uid_int = int(img_dir_name[1:])
                except ValueError:
                    continue
                niis = [f for f in os.listdir(img_dir)
                        if f.endswith(".nii") or f.endswith(".nii.gz")]
                if niis:
                    nii_index[uid_int] = os.path.join(img_dir, niis[0])
                    series_name_index[uid_int] = series_name

print(f"  -> {len(nii_index):,} NIfTI files indexed (should match img_map rows: {len(img_map):,})")

# Match ImageUID (strip 'I') -> NIfTI path
merged["uid_int"] = pd.to_numeric(
    merged["ImageUID"].astype(str).str.replace("I", "", regex=False),
    errors="coerce"
)
merged["nii_source_path"] = merged["uid_int"].apply(
    lambda x: nii_index.get(int(x), "") if pd.notna(x) else ""
)
merged["sourcedata_series_name"] = merged["uid_int"].apply(
    lambda x: series_name_index.get(int(x), "") if pd.notna(x) else ""
)
found = (merged["nii_source_path"] != "").sum()
print(f"  -> {found:,}/{len(merged):,} NIfTI files located")

# ── 7. Save outputs ────────────────────────────────────────────────────────────
cols_out = [
    "SubjectID", "ImageUID", "bids_sub", "bids_ses", "adni_viscode",
    "StudyDate", "StudyDate_norm", "SeriesUID", "study_int",
    "SeriesDescription", "sourcedata_series_name",
    "Manufacturer", "ManufacturerModelName", "MagneticFieldStrength",
    "SoftwareVersion", "Acceleration", "SliceThickness", "AcquisitionType",
    "StudyInstanceUID", "SeriesInstanceUID",
    "LONIStudy", "LONISeries", "LONIImage_mriqc",
    "mayo_qc_quality", "mayo_series_selected", "mayo_study_overallpass",
    "mayo_series_comments", "mayo_loni_study", "mayo_loni_series", "mayo_loni_image",
    "ADNI_phase", "mriqc_image_id",
    "nii_source_path",
]
cols_out = [c for c in cols_out if c in merged.columns]
merged[cols_out].to_csv(OUT_SESSION_MAP, index=False)
print(f"Saved: {OUT_SESSION_MAP}")

# Visit-code correspondence table
ses_visit = (
    merged[["bids_sub", "bids_ses", "adni_viscode", "StudyDate", "SubjectID"]]
    .drop_duplicates(subset=["bids_sub", "bids_ses"])
    .sort_values(["bids_sub", "bids_ses"])
)
ses_visit.to_csv(OUT_SES_VISIT, index=False)
print(f"Saved: {OUT_SES_VISIT}")

print(f"\nDone.")
print(f"  Total images: {len(merged):,}")
print(f"  Subjects: {merged['SubjectID'].nunique():,}")
print(f"  With MRIQC metadata: {merged['adni_viscode'].notna().sum():,}")
print(f"  With MAYOADIRL QC: {merged['mayo_qc_quality'].notna().sum():,}")
print(f"  With NIfTI path: {found:,}")
print(f"  Unique BIDS sessions: {merged[['bids_sub','bids_ses']].drop_duplicates().shape[0]:,}")
