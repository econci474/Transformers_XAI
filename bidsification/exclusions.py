"""
exclusions.py — Shared Subject Exclusion List
================================================
Import this in any pipeline script to apply consistent exclusion filters.

Usage:
    from exclusions import EXCLUDED_PTIDS, is_excluded_subject

    df = df[~df['SubjectID'].apply(is_excluded_subject)]
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
