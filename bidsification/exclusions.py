"""
exclusions.py — Shared Subject & Session Exclusion List
=========================================================
Import this in any pipeline script to apply consistent exclusion filters.

Usage:
    from exclusions import EXCLUDED_PTIDS, is_excluded_subject, is_excluded_session

    df = df[~df['SubjectID'].apply(is_excluded_subject)]
    df = df[~df.apply(lambda r: is_excluded_session(r['participant_id'], r['scan_date']), axis=1)]

Categories of exclusion:
    1. Site-381 cohort 10xxx (pattern, ALWAYS excluded)
       — ADNI repository data-quality flag
    2. Corrupted-MRI subjects (set, ALWAYS excluded)
       — per-subject data-quality issue (image unreadable)
    3. Diagnostic-reversion subjects (Excluded==1 in conversion_labels.tsv;
       OPT-IN via include_diagnostic_reversion=True)
       — non-sustained / reversion conversion-group filter, used by the
         post-exclusion cohort regeneration (clinical_pipeline 05b/06c/06d/
         07/07b/01d). Existing callers stay opt-out so their semantics
         don't silently change.
    4. Session-level MALFUNC scans (set, ALWAYS excluded)
       — scanner malfunction recorded in MRI3META.csv
"""

from functools import lru_cache
from pathlib import Path
import re

import pandas as pd

# ── Hard-excluded subject PTIDs ────────────────────────────────────────────────
# These subjects must be excluded from ALL pipeline outputs due to data
# quality concerns flagged by the ADNI data repository.

# Pattern: all PTIDs matching 381_S_10### (site 381, 10xxx IDs)
# 61 unique PTIDs confirmed present in MRIQC clinical tables (none in sourcedata).
EXCLUDED_PTID_PATTERNS = [
    r"^381_S_10\d+",   # Site 381, subject IDs starting with 10
]

# If you need to add individual PTIDs explicitly:
EXCLUDED_PTID_LIST = []  # Add specific PTIDs here if needed beyond pattern

# Corrupted-MRI subjects (data-quality, always excluded).
# Reason: image unreadable. Distinct from MALFUNC (per-session scanner
# malfunction) because the entire subject's imaging is unusable.
CORRUPTED_MRI_PTIDS = {
    "041_S_4629",   # corrupted MRI volumes
}


# ── Diagnostic-reversion subjects (opt-in) ─────────────────────────────────────
# Sourced from conversion_labels.tsv (column Excluded==1), produced by
# clinical_pipeline/07_conversion_group_extended.py. Subjects flagged as
# reverters or non-sustained conversions. Opt-in via kwarg so existing
# callers (e.g. ViT supervised fine-tune) keep their current semantics.
CONVERSION_LABELS_TSV = Path(
    r"D:\ADNI_SNP_Omni2.5M_20140220\conversion_labels\conversion_labels.tsv"
)


@lru_cache(maxsize=1)
def diagnostic_reversion_pids() -> frozenset:
    """Patient_IDs flagged Excluded==1 in conversion_labels.tsv.

    Single source of truth for the post-exclusion cohort. Returns an empty
    frozenset if the TSV is missing (so the module imports cleanly even on
    machines without the D: drive).
    """
    if not CONVERSION_LABELS_TSV.exists():
        return frozenset()
    cv = pd.read_csv(CONVERSION_LABELS_TSV, sep="\t",
                     usecols=["Patient_ID", "Excluded"])
    return frozenset(cv.loc[cv["Excluded"] == 1, "Patient_ID"].astype(str))


def is_excluded_subject(ptid: str,
                         *,
                         include_diagnostic_reversion: bool = False) -> bool:
    """Return True if a PTID should be excluded from the BIDS dataset.

    Args:
        ptid: Patient_ID in canonical form (e.g. '041_S_4629').
        include_diagnostic_reversion: opt-in to also exclude subjects
            flagged Excluded==1 in conversion_labels.tsv (diagnostic
            reversion / non-sustained conversion). Used by the
            post-exclusion regeneration scripts in clinical_pipeline.
            Default False — preserves existing caller semantics.

    Returns:
        True iff the PTID matches ANY of:
            - site-381 10xxx pattern (always)
            - EXCLUDED_PTID_LIST (always)
            - CORRUPTED_MRI_PTIDS (always)
            - diagnostic_reversion_pids()   [only when kwarg is True]
    """
    if not ptid or ptid != ptid:   # NaN check
        return False
    ptid = str(ptid).strip()
    if ptid in EXCLUDED_PTID_LIST:
        return True
    if ptid in CORRUPTED_MRI_PTIDS:
        return True
    if any(re.match(p, ptid) for p in EXCLUDED_PTID_PATTERNS):
        return True
    if include_diagnostic_reversion and ptid in diagnostic_reversion_pids():
        return True
    return False


# Convenience: set of all excluded PTIDs from patterns
# (populated lazily — not pre-computed since pattern-based)
EXCLUDED_PTIDS = set(EXCLUDED_PTID_LIST) | set(CORRUPTED_MRI_PTIDS)

# Summary for logging
EXCLUSION_REASON = {
    r"^381_S_10\d+":       "Data quality concerns flagged by ADNI repository (site 381, cohort 10xxx)",
    "corrupted_mri":       "Corrupted MRI data (image unreadable)",
    "diagnostic_reversion": "Diagnostic reversion / non-sustained conversion (conversion_labels.tsv Excluded==1)",
}

# ── Session-level exclusions ───────────────────────────────────────────────────
# Individual scan sessions excluded due to scanner malfunction (MALFUNC=1 in
# MRI3META.csv). Each entry is (participant_id, scan_date) where scan_date is
# the MRI acquisition date (EXAMDATE in MRI3META), formatted as 'YYYY-MM-DD'.
#
# Source: MRI3META.csv MALFUNC field; verified against our 616-subject cohort.
# MALFUNC=1 means a scanner malfunction was recorded during that exam.
# HAS_QC_ERROR is NOT used (administrative CRF flag, unrelated to image quality).
EXCLUDED_SESSIONS = {
    # (participant_id,  scan_date)       PTID          VISCODE2
    ('sub-023S2068',   '2010-12-01'),  # 023_S_2068    m03
    ('sub-153S2109',   '2011-02-02'),  # 153_S_2109    m03
    ('sub-073S2264',   '2011-02-03'),  # 073_S_2264    scmri
    ('sub-053S0919',   '2013-12-05'),  # 053_S_0919    m84
    ('sub-053S0919',   '2015-11-05'),  # 053_S_0919    m108
    ('sub-137S0994',   '2013-12-03'),  # 137_S_0994    m84
    ('sub-031S4029',   '2011-11-10'),  # 031_S_4029    m06
    ('sub-006S4192',   '2012-05-07'),  # 006_S_4192    m06
    ('sub-129S4371',   '2012-07-17'),  # 129_S_4371    m06
    ('sub-006S4449',   '2013-03-12'),  # 006_S_4449    m12
    ('sub-024S2239',   '2021-10-20'),  # 024_S_2239    m132
}

def is_excluded_session(participant_id: str, scan_date: str) -> bool:
    """Return True if a specific scan session should be excluded.

    Args:
        participant_id: BIDS participant ID, e.g. 'sub-023S2068'
        scan_date:      MRI acquisition date as 'YYYY-MM-DD'

    Returns:
        True if the session has MALFUNC=1 in MRI3META and should be excluded.
    """
    if not participant_id or not scan_date:
        return False
    return (str(participant_id).strip(), str(scan_date).strip()[:10]) in EXCLUDED_SESSIONS
