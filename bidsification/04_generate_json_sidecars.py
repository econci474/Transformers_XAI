"""
Script 04 — Generate T1w JSON Sidecars
=======================================
BIDS requires (or strongly recommends) a JSON sidecar alongside every NIfTI.
This script generates _T1w.json files for all selected scans using acquisition
metadata available from scan_selection.csv (sourced from MRIQC join).

Fields populated:
  MagneticFieldStrength  : always 3 (3T dataset confirmed)
  Manufacturer           : from MRIQC (GE / Siemens / Philips)
  ManufacturerModelName  : scanner model from MRIQC
  SliceThickness         : from MRIQC (mm)
  SeriesDescription      : ADNI series name (e.g. MT1__GradWarp__N3m)

Note: TR, TE, TI, FlipAngle are not available in the ADNI MRIQC table —
these are RECOMMENDED not REQUIRED by BIDS for T1w anatomicals and will
generate INFO-level (not error) warnings in the validator.
"""

import pandas as pd
import os
import json
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exclusions import is_excluded_subject

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT = r"D:\ADNI_BIDS_project"
BIDS_DIR     = os.path.join(PROJECT_ROOT, "bids")
SEL_CSV      = os.path.join(PROJECT_ROOT, "metadata", "scan_selection.csv")

# ── Load scan selection ────────────────────────────────────────────────────────
sel = pd.read_csv(SEL_CSV, dtype={"bids_ses": str}, low_memory=False)
sel = sel[~sel["SubjectID"].apply(is_excluded_subject)]

written = skipped_exists = errors = 0

for _, row in sel.iterrows():
    bids_sub = str(row["bids_sub"])
    bids_ses = str(row["bids_ses"])

    anat_dir = os.path.join(BIDS_DIR, f"sub-{bids_sub}", f"ses-{bids_ses}", "anat")
    nifti    = os.path.join(anat_dir, f"sub-{bids_sub}_ses-{bids_ses}_T1w.nii.gz")
    json_out = nifti.replace(".nii.gz", ".json")

    if not os.path.isfile(nifti):
        continue  # NIfTI not copied

    if os.path.isfile(json_out):
        skipped_exists += 1
        continue

    # Build sidecar — MagneticFieldStrength always 3 (3T confirmed)
    sidecar = {
        "MagneticFieldStrength": 3,
        "Modality": "MRI",
    }
    for field, col in [
        ("Manufacturer",          "Manufacturer"),
        ("ManufacturerModelName", "ManufacturerModelName"),
        ("SliceThickness",        "SliceThickness"),
        ("SeriesDescription",     "SeriesDescription_selected"),
    ]:
        val = row.get(col, None)
        if pd.notna(val) and val != "":
            if field == "SliceThickness":
                try:
                    sidecar[field] = float(val)
                except (ValueError, TypeError):
                    pass
            else:
                sidecar[field] = str(val)

    try:
        os.makedirs(anat_dir, exist_ok=True)
        with open(json_out, "w") as f:
            json.dump(sidecar, f, indent=2)
        written += 1
    except Exception as e:
        print(f"  ERROR: {json_out}: {e}")
        errors += 1

print(f"JSON sidecars written:  {written:,}")
print(f"Already existed:        {skipped_exists:,}")
print(f"Errors:                 {errors}")
