"""
Script 09 — Extract DICOM Tags and Update T1w JSON Sidecars
============================================================
Matches selected processed scans to their original raw DICOM via
SubjectID + acquisition date (date embedded in the zip directory name).

Zip structure: ADNI/<SubjectID>/<SeriesDesc>/<YYYY-MM-DD_HH_MM_SS.0>/<ImageUID>/*.dcm
Processed NIfTI path: .../ADNI/<SubjectID>/<ProcessedSeries>/<YYYY-MM-DD_...>/<ImageUID>/*.nii

Strategy:
  1. Build lookup: (SubjectID, YYYY-MM-DD) -> (zip_path, first_dcm_entry)
     from all 10 zip files
  2. For each selected scan, extract the acquisition date from the
     nii_source_selected path (which encodes date in the folder name)
  3. Read the first DICOM for that (SubjectID, date) and extract tags
  4. Merge extracted tags into the existing _T1w.json sidecar

BIDS recommended tags populated:
  RepetitionTime, EchoTime, InversionTime, FlipAngle,
  ScanningSequence, SequenceVariant, SequenceName,
  PhaseEncodingDirection
"""

import zipfile, io, os, json, logging, re
import pandas as pd
import pydicom
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exclusions import is_excluded_subject

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger(__name__)

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT = r"D:\ADNI_BIDS_project"
BIDS_DIR     = os.path.join(PROJECT_ROOT, "bids")
SEL_CSV      = os.path.join(PROJECT_ROOT, "metadata", "scan_selection.csv")

ZIP_FILES = [
    rf"D:\ADNI_all_3T_original_{i}.zip" for i in range(1, 11)
    if os.path.isfile(rf"D:\ADNI_all_3T_original_{i}.zip")
]
log.info(f"Zip files found: {len(ZIP_FILES)}")

# ── DICOM tags to extract (BIDS recommended for T1w) ──────────────────────────
DICOM_TAGS = {
    "RepetitionTime":         (0x0018, 0x0080),
    "EchoTime":               (0x0018, 0x0081),
    "InversionTime":          (0x0018, 0x0082),
    "FlipAngle":              (0x0018, 0x1314),
    "ScanningSequence":       (0x0018, 0x0020),
    "SequenceVariant":        (0x0018, 0x0021),
    "SequenceName":           (0x0018, 0x0024),
    "PhaseEncodingDirection": (0x0018, 0x1312),
}

# Date pattern in zip folder: YYYY-MM-DD_HH_MM_SS.0
DATE_RE = re.compile(r'(\d{4}-\d{2}-\d{2})_')

def extract_date_from_path(path_str):
    """Extract YYYY-MM-DD from path component like 2011-06-02_07_58_50.0"""
    m = DATE_RE.search(str(path_str))
    return m.group(1) if m else None

# ── Load selected scans ────────────────────────────────────────────────────────
sel = pd.read_csv(SEL_CSV, dtype={"bids_ses": str}, low_memory=False)
sel = sel[~sel["SubjectID"].apply(is_excluded_subject)]

# Build target dict: (SubjectID, YYYY-MM-DD) -> (bids_sub, bids_ses)
target_by_date = {}
for _, row in sel.iterrows():
    nii_src = str(row.get("nii_source_selected", ""))
    acq_date = extract_date_from_path(nii_src)
    if acq_date:
        key = (str(row["SubjectID"]), acq_date)
        target_by_date[key] = (str(row["bids_sub"]), str(row["bids_ses"]))

log.info(f"Target (SubjectID, date) pairs: {len(target_by_date)}")

# ── Phase 1: Build zip index ((SubjectID, date) -> (zip_path, dcm_entry)) ─────
log.info("Phase 1: Indexing zip files...")
date_to_dcm = {}  # (SubjectID, date) -> (zip_path, dcm_entry)

for zip_path in ZIP_FILES:
    log.info(f"  Indexing {os.path.basename(zip_path)} ...")
    with zipfile.ZipFile(zip_path, "r") as z:
        for entry in z.namelist():
            if not entry.endswith(".dcm"):
                continue
            parts = entry.replace("\\", "/").split("/")
            if len(parts) < 5:
                continue
            subject_id = parts[1]   # e.g. 032_S_0187
            date_folder = parts[3]  # e.g. 2006-03-14_10_56_28.0
            acq_date    = extract_date_from_path(date_folder)
            if not acq_date:
                continue
            key = (subject_id, acq_date)
            if key in target_by_date and key not in date_to_dcm:
                date_to_dcm[key] = (zip_path, entry)

log.info(f"DICOMs found for selected scans: {len(date_to_dcm)} / {len(target_by_date)}")

# ── Phase 2: Extract DICOM tags and update JSON sidecars ──────────────────────
log.info("Phase 2: Extracting tags and updating JSON sidecars...")

updated = skipped_exists = missing_json = errors = 0

from collections import defaultdict
by_zip = defaultdict(list)
for key, (zip_path, entry) in date_to_dcm.items():
    by_zip[zip_path].append((key, entry))

for zip_path, items in by_zip.items():
    log.info(f"  Processing {os.path.basename(zip_path)} ({len(items)} sessions)...")
    with zipfile.ZipFile(zip_path, "r") as z:
        for (subject_id, acq_date), entry in items:
            bids_sub, bids_ses = target_by_date[(subject_id, acq_date)]
            json_path = os.path.join(
                BIDS_DIR, f"sub-{bids_sub}", f"ses-{bids_ses}", "anat",
                f"sub-{bids_sub}_ses-{bids_ses}_T1w.json"
            )
            if not os.path.isfile(json_path):
                missing_json += 1
                continue

            with open(json_path) as f:
                sidecar = json.load(f)

            if all(k in sidecar for k in ("RepetitionTime", "ScanningSequence", "SequenceName")):
                skipped_exists += 1
                continue

            try:
                dcm_bytes = z.read(entry)
                ds = pydicom.dcmread(
                    io.BytesIO(dcm_bytes),
                    stop_before_pixels=True,
                    force=True
                )
            except Exception as e:
                errors += 1
                log.warning(f"    DICOM read error ({subject_id}, {acq_date}): {e}")
                continue

            new_fields = {}
            for field_name, tag in DICOM_TAGS.items():
                if field_name in sidecar:
                    continue
                try:
                    val = ds[tag].value
                    if isinstance(val, pydicom.multival.MultiValue):
                        val = [str(v).strip() for v in val] if len(val) > 1 else str(val[0]).strip()
                    elif isinstance(val, pydicom.uid.UID):
                        val = str(val)
                    elif isinstance(val, (int, float)):
                        val = float(val) if isinstance(val, float) else int(val)
                    else:
                        val = str(val).strip()
                    if val not in ("", None, []):
                        new_fields[field_name] = val
                except (KeyError, AttributeError):
                    pass

            if new_fields:
                sidecar.update(new_fields)
                with open(json_path, "w") as f:
                    json.dump(sidecar, f, indent=2)
                updated += 1

not_found = len(target_by_date) - len(date_to_dcm)
print(f"\n=== Summary ===")
print(f"JSON sidecars updated:  {updated:,}")
print(f"Already complete:       {skipped_exists:,}")
print(f"Sessions not in zips:   {not_found:,}")
print(f"Missing JSON sidecars:  {missing_json:,}")
print(f"Errors:                 {errors:,}")

# Show a sample updated sidecar
import glob
sample_jsons = glob.glob(os.path.join(BIDS_DIR, "sub-*", "ses-*", "anat", "*_T1w.json"))
if sample_jsons and updated > 0:
    with open(sample_jsons[0]) as f:
        s = json.load(f)
    print(f"\nSample sidecar ({os.path.relpath(sample_jsons[0], BIDS_DIR)}):")
    print(json.dumps(s, indent=2))
