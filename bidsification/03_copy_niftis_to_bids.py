"""
Script 03 — NIfTI Renaming + Gzip Compression into BIDS Layout
================================================================
Reads scan_selection.csv and copies/gzips the source NIfTI files
into the BIDS directory structure:
    bids/sub-<label>/ses-<label>/anat/sub-<label>_ses-<label>[_run-N]_T1w.nii.gz

Rules:
  - If n_scans == 1: single file, no run label
  - If n_scans > 1:  run-1 and run-2 labels
  - Source NIfTIs are N3m-corrected (not raw DICOMs)
  - Files are gzip-compressed during copy

Skips files where nii_source_path is empty (not found in sourcedata).
"""

import pandas as pd
import os
import gzip
import shutil
import logging

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT = r"D:\ADNI_BIDS_project"
SELECTION_CSV = os.path.join(PROJECT_ROOT, "metadata", "scan_selection.csv")
BIDS_DIR = os.path.join(PROJECT_ROOT, "bids")
LOG_FILE = os.path.join(PROJECT_ROOT, "metadata", "03_nifti_copy.log")

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
    """Copy src → dst.gz (gzip compress). dst should already end in .nii.gz."""
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(src, "rb") as f_in, gzip.open(dst, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

def bids_anat_dir(bids_root: str, sub: str, ses: str) -> str:
    return os.path.join(bids_root, f"sub-{sub}", f"ses-{ses}", "anat")

# ── Main ───────────────────────────────────────────────────────────────────────
sel = pd.read_csv(SELECTION_CSV, low_memory=False)
log.info(f"Loaded scan_selection.csv: {len(sel):,} sessions")

copied = 0
skipped_no_source = 0
skipped_exists = 0
errors = 0

for _, row in sel.iterrows():
    sub = str(row["bids_sub"])
    ses = str(row["bids_ses"])
    n   = int(row.get("n_scans", 1))
    anat_dir = bids_anat_dir(BIDS_DIR, sub, ses)

    def copy_run(src_path: str, run_label: str):
        """Copy a single NIfTI to BIDS, adding run label if provided."""
        global copied, skipped_no_source, skipped_exists, errors
        if not src_path or not isinstance(src_path, str) or not os.path.isfile(src_path):
            skipped_no_source += 1
            log.warning(f"  SKIP (no source): sub-{sub} ses-{ses} run={run_label} | {src_path}")
            return
        
        fname = f"sub-{sub}_ses-{ses}"
        if run_label:
            fname += f"_{run_label}"
        fname += "_T1w.nii.gz"
        dst = os.path.join(anat_dir, fname)
        
        if os.path.isfile(dst):
            skipped_exists += 1
            log.debug(f"  EXISTS: {dst}")
            return
        
        try:
            gzip_copy(src_path, dst)
            copied += 1
            log.info(f"  OK: {fname}")
        except Exception as e:
            errors += 1
            log.error(f"  ERROR copying {src_path} → {dst}: {e}")

    if n == 1:
        # Single scan: no run label
        copy_run(str(row.get("nii_source_run1", "")), "")
    else:
        # Multiple scans: use run-1, run-2
        copy_run(str(row.get("nii_source_run1", "")), "run-1")
        src2 = str(row.get("nii_source_run2", ""))
        if src2 and src2 != "nan":
            copy_run(src2, "run-2")

log.info("=" * 60)
log.info(f"Copied:           {copied:,}")
log.info(f"Already existed:  {skipped_exists:,}")
log.info(f"No source found:  {skipped_no_source:,}")
log.info(f"Errors:           {errors:,}")
log.info(f"Log: {LOG_FILE}")
