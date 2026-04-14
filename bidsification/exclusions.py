"""
exclusions.py — Shared Subject & Session Exclusion List
=========================================================
Import this in any pipeline script to apply consistent exclusion filters.

Usage:
    from exclusions import EXCLUDED_PTIDS, is_excluded_subject, is_excluded_session

    df = df[~df['SubjectID'].apply(is_excluded_subject)]
    df = df[~df.apply(lambda r: is_excluded_session(r['participant_id'], r['scan_date']), axis=1)]
"""

import re

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

def is_excluded_subject(ptid: str) -> bool:
    """Return True if a PTID should be excluded from the BIDS dataset."""
    if not ptid or ptid != ptid:   # NaN check
        return False
    ptid = str(ptid).strip()
    if ptid in EXCLUDED_PTID_LIST:
        return True
    return any(re.match(p, ptid) for p in EXCLUDED_PTID_PATTERNS)

# Convenience: set of all excluded PTIDs from patterns
# (populated lazily — not pre-computed since pattern-based)
EXCLUDED_PTIDS = set(EXCLUDED_PTID_LIST)

# Summary for logging
EXCLUSION_REASON = {
    r"^381_S_10\d+": "Data quality concerns flagged by ADNI repository (site 381, cohort 10xxx)",
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
