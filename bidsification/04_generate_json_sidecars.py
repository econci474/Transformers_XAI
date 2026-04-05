"""
Script 04 — Generate JSON Sidecars from MRIQUALITY Metadata
=============================================================
For each NIfTI file already placed in bids/, generates the 
corresponding T1w.json sidecar file using scanner metadata from
MRIQUALITY.csv and the session_map.csv.

BIDS-required/recommended fields populated:
    MagneticFieldStrength, Manufacturer, ManufacturerModelName,
    SoftwareVersions, PulseSequenceType, MRAcquisitionType,
    AcquisitionVoxelSize (SliceThickness), FlipAngle, RepetitionTime,
    ParallelAcquisitionTechnique, MultibandAccelerationFactor

Note: TR, TE, flip angle are NOT in MRIQUALITY; they are DICOM-derived.
      Those fields are left as null and should be filled from DICOMs if 
      available in the future.
"""

import pandas as pd
import os
import json
import glob
import re

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT = r"D:\ADNI_BIDS_project"
SESSION_MAP_CSV = os.path.join(PROJECT_ROOT, "metadata", "session_map.csv")
SELECTION_CSV = os.path.join(PROJECT_ROOT, "metadata", "scan_selection.csv")
BIDS_DIR = os.path.join(PROJECT_ROOT, "bids")

# ── BIDS JSON template for ADNI T1w ───────────────────────────────────────────
def make_sidecar(row_sess, row_sel=None) -> dict:
    """
    Build a BIDS-compliant T1w JSON sidecar from MRIQUALITY metadata.
    
    Parameters
    ----------
    row_sess : pd.Series — row from session_map.csv (scanner metadata)
    row_sel  : pd.Series — row from scan_selection.csv (QC, series info)
    """
    # MagneticFieldStrength
    mfs = row_sess.get("MagneticFieldStrength", None)
    try:
        mfs = float(mfs) if pd.notna(mfs) else None
    except (TypeError, ValueError):
        mfs = None

    # Acceleration
    accel = str(row_sess.get("Acceleration", "")).strip()

    # SeriesDescription -> PulseSequenceType heuristic
    desc = str(row_sess.get("SeriesDescription", "")).upper()
    if "MPRAGE" in desc or "MP-RAGE" in desc:
        pulse_seq = "MPRAGE"
    elif "IR-FSPGR" in desc or "IRFSPGR" in desc:
        pulse_seq = "IR-FSPGR"
    elif "SPGR" in desc or "FSPGR" in desc:
        pulse_seq = "FSPGR"
    else:
        pulse_seq = "UNKNOWN"

    sidecar = {
        "Manufacturer":             row_sess.get("Manufacturer", None),
        "ManufacturerModelName":    row_sess.get("ManufacturerModelName", None),
        "MagneticFieldStrength":    mfs,
        "MRAcquisitionType":        "3D",
        "PulseSequenceType":        pulse_seq,
        "SoftwareVersions":         row_sess.get("SoftwareVersion", None)
                                    if "SoftwareVersion" in row_sess.index else None,
        "SliceThickness":           float(row_sess["SliceThickness"])
                                    if "SliceThickness" in row_sess.index
                                    and pd.notna(row_sess.get("SliceThickness")) else None,
        "RepetitionTime":           None,   # Not in MRIQUALITY; needs DICOM
        "EchoTime":                 None,   # Not in MRIQUALITY; needs DICOM
        "FlipAngle":                None,   # Not in MRIQUALITY; needs DICOM
        "InversionTime":            None,   # Not in MRIQUALITY; needs DICOM
        "ParallelAcquisitionTechnique": accel if accel not in ("", "nan", "None") else None,
        "NumberShots":              None,
        # ADNI-specific provenance
        "ADNISeriesDescription":    row_sess.get("SeriesDescription", None),
        "ADNISeriesType":           row_sess.get("SeriesType", None),
        "ADNIImageUID":             str(int(row_sess["ImageUID_int"]))
                                    if pd.notna(row_sess.get("ImageUID_int")) else None,
    }
    
    # Add QC info if available
    if row_sel is not None:
        sidecar["ADNIQCQuality"] = (
            int(row_sel["qc_quality_selected"])
            if pd.notna(row_sel.get("qc_quality_selected")) else None
        )
        sidecar["ADNISelectionReason"] = row_sel.get("selection_reason", None)

    # Remove None values (optional — keeps JSON cleaner but BIDS allows nulls)
    # Uncomment the next line to strip nulls:
    # sidecar = {k: v for k, v in sidecar.items() if v is not None}
    return sidecar

# ── Load tables ────────────────────────────────────────────────────────────────
print("Loading tables ...")
sess_map = pd.read_csv(SESSION_MAP_CSV, low_memory=False)
sel = pd.read_csv(SELECTION_CSV, low_memory=False)
sess_map["ImageUID_int"] = pd.to_numeric(sess_map.get("ImageUID", pd.Series()), errors="coerce")

# Merge selection into session map for QC info
sess_map = sess_map.merge(
    sel[["bids_sub", "bids_ses", "qc_quality_selected", "selection_reason",
         "ImageUID_selected", "ImageUID_run2"]].rename(
            columns={"ImageUID_selected": "ImageUID_sel", "ImageUID_run2": "ImageUID_run2_sel"}
         ).astype({"ImageUID_sel": str, "ImageUID_run2_sel": str}),
    on=["bids_sub", "bids_ses"], how="left"
)

# Build ImageUID -> row lookup
sess_map["ImageUID_str"] = sess_map["ImageUID_int"].apply(
    lambda x: str(int(x)) if pd.notna(x) else ""
)
img_lookup = sess_map.set_index("ImageUID_str")

# ── Walk BIDS directory and create sidecars ────────────────────────────────────
nii_files = glob.glob(
    os.path.join(BIDS_DIR, "sub-*", "ses-*", "anat", "*_T1w.nii.gz")
)
print(f"Found {len(nii_files):,} T1w NIfTI files in bids/")

created = 0
skipped = 0

for nii_path in sorted(nii_files):
    json_path = nii_path.replace("_T1w.nii.gz", "_T1w.json")
    if os.path.isfile(json_path):
        skipped += 1
        continue
    
    # Parse filename to get sub/ses
    fname = os.path.basename(nii_path)
    m = re.match(r"sub-(\S+?)_ses-(\S+?)(?:_run-\d+)?_T1w\.nii\.gz", fname)
    if not m:
        print(f"  WARNING: Could not parse filename: {fname}")
        skipped += 1
        continue
    
    bids_sub = m.group(1)
    bids_ses = m.group(2)
    is_run2 = "_run-2_" in fname
    
    # Find matching rows in sess_map
    matches = sess_map[
        (sess_map["bids_sub"] == bids_sub) &
        (sess_map["bids_ses"] == bids_ses)
    ]
    
    if len(matches) == 0:
        print(f"  WARNING: No session map entry for sub-{bids_sub} ses-{bids_ses}")
        skipped += 1
        continue
    
    # Pick the right scan (run-1=selected, run-2=second)
    if is_run2:
        run2_uid = str(matches.iloc[0].get("ImageUID_run2_sel", ""))
        row = matches[matches["ImageUID_str"] == run2_uid]
        if len(row) == 0:
            row = matches.iloc[[1]] if len(matches) > 1 else matches.iloc[[0]]
    else:
        sel_uid = str(matches.iloc[0].get("ImageUID_sel", ""))
        row = matches[matches["ImageUID_str"] == sel_uid]
        if len(row) == 0:
            row = matches.iloc[[0]]
    
    row_sess = row.iloc[0]
    
    # Find matching sel row
    sel_row = sel[
        (sel["bids_sub"] == bids_sub) & (sel["bids_ses"] == bids_ses)
    ]
    row_sel = sel_row.iloc[0] if len(sel_row) > 0 else None
    
    sidecar = make_sidecar(row_sess, row_sel)
    
    with open(json_path, "w") as f:
        json.dump(sidecar, f, indent=2, default=str)
    
    created += 1

print(f"\nJSON sidecars created: {created:,}")
print(f"Already existed (skipped): {skipped:,}")
