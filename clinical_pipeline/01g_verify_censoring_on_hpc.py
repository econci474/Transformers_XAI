"""
01g_verify_censoring_on_hpc.py
==============================
Verify that the censoring-aware Label_Ny patch (01f) landed correctly on
CSD3 after scp -- i.e. the HPC copies of the MRI master CSV and the
no_cdr_stratified baseline/seed splits match the local patched files.

WHY: 01f_apply_censoring_to_Label_Ny.py ran only on local D:. The CSD3
copies were the pre-patch versions (Label_Ny=0 mixed true-negatives with
right-censored subjects). After re-scp, run THIS on CSD3 to confirm the
bytes/labels are the censoring-aware ones the T3a-T3d sweeps will consume.

The EXPECTED_* tables below were captured from the local patched files on
2026-06-02 (master md5 fcc5c47...). The script asserts the HPC files
reproduce them exactly. Any mismatch => the wrong (stale) file is on CSD3,
or the scp did not overwrite.

Run (CSD3, env: mri or any pandas env):
    python clinical_pipeline/01g_verify_censoring_on_hpc.py

Override paths if your layout differs:
    python clinical_pipeline/01g_verify_censoring_on_hpc.py \
        --master /home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv \
        --baseline-dir /home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline

Exit code 0 = all checks pass; 1 = at least one mismatch.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Default CSD3 paths (from the T3abcd submit-script headers)
# --------------------------------------------------------------------------- #
DEFAULT_MASTER = Path(
    "/home/ec474/rds/hpc-work/ADNI_MRI/"
    "master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
)
DEFAULT_BASELINE_DIR = Path(
    "/home/ec474/rds/hpc-work/ADNI_CL/"
    "no_cdr_stratified_post_exclusion/tabular/baseline"
)

# --------------------------------------------------------------------------- #
# Reference values captured from the LOCAL patched files (2026-06-02)
# --------------------------------------------------------------------------- #
# md5 of the local patched master. If scp preserved bytes exactly this will
# match; a CRLF/encoding change can break md5 while the COUNTS still match --
# in that case the md5 line warns but the semantic check is authoritative.
MASTER_MD5 = "fcc5c475d16ef7df2a6804cb4cd51105"
MASTER_ROWS = 1456

# Label_Ny -> (pos, neg, nan) on the master, horizons 1..10.
EXPECTED_MASTER = {
    "Label_1y":  (32, 1420,   4),
    "Label_2y":  (91, 1328,  37),
    "Label_3y": (145, 1129, 182),
    "Label_4y": (189, 1040, 227),
    "Label_5y": (226,  812, 418),
    "Label_6y": (249,  683, 524),
    "Label_7y": (282,  587, 587),
    "Label_8y": (295,  493, 668),
    "Label_9y": (326,  402, 728),
    "Label_10y":(338,  341, 777),
}

# (seed, split) -> {"rows": n, "Label_Ny": (pos, neg, nan)} for the T3 horizons.
EXPECTED_SPLITS = {
    (0, "train"): {"rows": 456, "Label_3y": (48, 323, 85),  "Label_5y": (74, 231, 151), "Label_7y": (88, 173, 195), "Label_10y": (105, 103, 248)},
    (0, "val"):   {"rows": 57,  "Label_3y": (4, 37, 16),    "Label_5y": (8, 22, 27),    "Label_7y": (9, 18, 30),    "Label_10y": (11, 8, 38)},
    (0, "test"):  {"rows": 58,  "Label_3y": (12, 37, 9),    "Label_5y": (13, 26, 19),   "Label_7y": (14, 16, 28),   "Label_10y": (16, 10, 32)},
    (1, "train"): {"rows": 456, "Label_3y": (51, 318, 87),  "Label_5y": (77, 221, 158), "Label_7y": (93, 165, 198), "Label_10y": (113, 96, 247)},
    (1, "val"):   {"rows": 57,  "Label_3y": (4, 40, 13),    "Label_5y": (8, 29, 20),    "Label_7y": (8, 24, 25),    "Label_10y": (9, 16, 32)},
    (1, "test"):  {"rows": 58,  "Label_3y": (9, 39, 10),    "Label_5y": (10, 29, 19),   "Label_7y": (10, 18, 30),   "Label_10y": (10, 9, 39)},
    (2, "train"): {"rows": 456, "Label_3y": (51, 312, 93),  "Label_5y": (79, 222, 155), "Label_7y": (93, 166, 197), "Label_10y": (111, 95, 250)},
    (2, "val"):   {"rows": 57,  "Label_3y": (5, 43, 9),     "Label_5y": (6, 28, 23),    "Label_7y": (8, 18, 31),    "Label_10y": (9, 13, 35)},
    (2, "test"):  {"rows": 58,  "Label_3y": (8, 42, 8),     "Label_5y": (10, 29, 19),   "Label_7y": (10, 23, 25),   "Label_10y": (12, 13, 33)},
}


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def counts(series: pd.Series) -> tuple[int, int, int]:
    s = pd.to_numeric(series, errors="coerce")
    return int((s == 1).sum()), int((s == 0).sum()), int(s.isna().sum())


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------- #
# Checks
# --------------------------------------------------------------------------- #
def check_master(path: Path) -> bool:
    print(f"\n[master] {path}")
    if not path.exists():
        print("  [FAIL] file does not exist on HPC -- scp did not land here.")
        return False

    # Bytes-identical check (informational; counts are authoritative).
    digest = md5_of(path)
    if digest == MASTER_MD5:
        print(f"  [ok]   md5 matches local patched ({digest})")
    else:
        print(f"  [warn] md5 differs (HPC {digest} vs local {MASTER_MD5}).")
        print("         Likely a line-ending/encoding change from transfer; "
              "semantic label check below is authoritative.")

    df = pd.read_csv(path, low_memory=False)
    ok = True
    if len(df) != MASTER_ROWS:
        print(f"  [FAIL] row count {len(df)} != expected {MASTER_ROWS}")
        ok = False
    else:
        print(f"  [ok]   rows = {len(df)}")

    prev_nan = -1
    for col, exp in EXPECTED_MASTER.items():
        if col not in df.columns:
            print(f"  [FAIL] {col} column missing")
            ok = False
            continue
        got = counts(df[col])
        tag = "ok  " if got == exp else "FAIL"
        if got != exp:
            ok = False
        print(f"  [{tag}] {col:9s} pos/neg/nan got={got}  exp={exp}")
        # Monotone censoring sanity (independent of the baked table).
        if got[2] < prev_nan:
            print(f"  [FAIL] {col} censored count {got[2]} < previous horizon "
                  f"{prev_nan} (should be non-decreasing)")
            ok = False
        prev_nan = got[2]

    # The headline patch signature: long horizons MUST have NaNs now.
    n7 = counts(df["Label_7y"])[2] if "Label_7y" in df else 0
    n10 = counts(df["Label_10y"])[2] if "Label_10y" in df else 0
    if n7 == 0 or n10 == 0:
        print("  [FAIL] Label_7y/Label_10y have ZERO NaNs -- this is the "
              "PRE-patch (stale) file. Re-scp the local patched master.")
        ok = False
    return ok


def check_splits(base_dir: Path) -> bool:
    print(f"\n[splits] {base_dir}")
    if not base_dir.exists():
        print("  [FAIL] baseline dir does not exist on HPC.")
        return False
    ok = True
    for (seed, split), exp in EXPECTED_SPLITS.items():
        p = base_dir / f"seed_{seed}" / f"{split}.csv"
        if not p.exists():
            print(f"  [FAIL] missing {p}")
            ok = False
            continue
        df = pd.read_csv(p, low_memory=False)
        row_ok = len(df) == exp["rows"]
        if not row_ok:
            ok = False
        line = [f"seed{seed}/{split:5s} rows={len(df)}({'ok' if row_ok else 'FAIL exp '+str(exp['rows'])})"]
        for col in ("Label_3y", "Label_5y", "Label_7y", "Label_10y"):
            if col not in df.columns:
                line.append(f"{col}:MISSING")
                ok = False
                continue
            got = counts(df[col])
            ecol = exp[col]
            if got != ecol:
                ok = False
                line.append(f"{col}:FAIL got{got}!=exp{ecol}")
            else:
                line.append(f"{col}:ok")
        print("  " + "  ".join(line))
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", type=Path, default=DEFAULT_MASTER)
    ap.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE_DIR)
    args = ap.parse_args()

    print("=" * 64)
    print("  01g -- verify censoring-aware Label_Ny on HPC")
    print("=" * 64)

    m_ok = check_master(args.master)
    s_ok = check_splits(args.baseline_dir)

    print("\n" + "=" * 64)
    if m_ok and s_ok:
        print("  RESULT: PASS -- HPC copies match the local patched files.")
        print("  Safe to run the T3a-T3d sweeps.")
        print("=" * 64)
        return 0
    print("  RESULT: FAIL -- see [FAIL] lines above. Do NOT run T3 sweeps")
    print("  until the correct patched files are in place on CSD3.")
    print("=" * 64)
    return 1


if __name__ == "__main__":
    sys.exit(main())
