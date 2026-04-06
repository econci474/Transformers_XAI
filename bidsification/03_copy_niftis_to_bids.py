"""
Script 03 — Copy Selected NIfTIs into BIDS Layout
===================================================
Reads scan_selection.csv and copies + gzip-compresses the SELECTED
(best quality) NIfTI for each session into the BIDS directory.

Design decisions:
  - ONE file per session, NO run label: sub-<sub>_ses-<ses>_T1w.nii.gz
    This is the cleanest input for sMRIprep (no ambiguity about which
    T1w to use). The second-best scan (run-2) is documented in
    scan_selection.csv but NOT copied to bids/ — it remains in sourcedata.
  - Source files are .nii (uncompressed); output is gzip-compressed .nii.gz
  - Skips sessions where nii_source_selected is empty (no downloaded file)
  - Skips if the output file already exists (idempotent)

Output structure:
    bids/sub-<label>/ses-<YYYYMMDD>/anat/sub-<label>_ses-<YYYYMMDD>_T1w.nii.gz
"""

import pandas as pd
import os
import gzip
import shutil
import logging

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT   = r"D:\ADNI_BIDS_project"
SELECTION_CSV  = os.path.join(PROJECT_ROOT, "metadata", "scan_selection.csv")
BIDS_DIR       = os.path.join(PROJECT_ROOT, "bids")
LOG_FILE       = os.path.join(PROJECT_ROOT, "metadata", "03_nifti_copy.log")

# ── Logging ────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE, mode="w"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)

def gzip_copy(src: str, dst: str):
    """Gzip-compress src -> dst. dst should end in .nii.gz"""
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(src, "rb") as f_in, gzip.open(dst, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

# ── Load ────────────────────────────────────────────────────────────────────────
sel = pd.read_csv(SELECTION_CSV, low_memory=False)
log.info(f"Loaded scan_selection.csv: {len(sel):,} sessions")

copied = 0
skipped_no_source = 0
skipped_exists = 0
errors = 0

for _, row in sel.iterrows():
    sub = str(row["bids_sub"])
    ses = str(row["bids_ses"])

    src_path = str(row.get("nii_source_selected", ""))
    if not src_path or src_path == "nan" or not os.path.isfile(src_path):
        skipped_no_source += 1
        log.warning(f"  SKIP (no source): sub-{sub} ses-{ses} | {src_path}")
        continue

    # Single file per session, no run label
    fname = f"sub-{sub}_ses-{ses}_T1w.nii.gz"
    anat_dir = os.path.join(BIDS_DIR, f"sub-{sub}", f"ses-{ses}", "anat")
    dst = os.path.join(anat_dir, fname)

    if os.path.isfile(dst):
        skipped_exists += 1
        continue

    try:
        # If source is already .nii.gz, copy directly; otherwise gzip-compress
        if src_path.endswith(".nii.gz"):
            os.makedirs(anat_dir, exist_ok=True)
            shutil.copy2(src_path, dst)
        else:
            gzip_copy(src_path, dst)
        copied += 1
        log.info(f"  OK: {fname}")
    except Exception as e:
        errors += 1
        log.error(f"  ERROR: {fname}: {e}")

log.info("=" * 60)
log.info(f"Copied:          {copied:,}")
log.info(f"Already existed: {skipped_exists:,}")
log.info(f"No source found: {skipped_no_source:,}")
log.info(f"Errors:          {errors:,}")
log.info(f"Log: {LOG_FILE}")
