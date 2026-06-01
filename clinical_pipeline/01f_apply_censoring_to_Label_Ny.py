"""
01f_apply_censoring_to_Label_Ny.py
====================================
Apply right-censoring to Label_1y..Label_10y in every derivative CSV that
carries them.

Background
----------
`01b_build_clinical_csv.py::compute_prognosis` defines Label_Ny as:
    1 if any AD diagnosis exists in (baseline, baseline + N years]
    0 otherwise (INCLUDING subjects whose last visit was before that
    window closed -- right-censored)

This silently mixes:
    true negatives  (FU >= N y, no AD)
  + censored 0s    (FU <  N y, unknown later state)

For T3c (7y) and especially T3d (10y) the censored share is substantial.

This patch reads `conversion_labels.tsv` (the source of FU_years +
years_to_AD per subject), and for every CSV under the clinical +
mri_clinical_matched derivative trees that has Label_Ny columns,
recomputes them with:

    Label_N = 1            if years_to_AD <= N   (converted in window)
    Label_N = NaN          if (FU_years < N) AND (years_to_AD is NaN)
                              (right-censored, drop downstream)
    Label_N = 0            otherwise            (true negative)

Subjects converting AFTER the horizon (years_to_AD > N) are correct 0s.

Affected directories (recursive walk for *.csv with any Label_<n>y col):
  - D:/ADNI_BIDS_project/derivatives/clinical/...
  - D:/ADNI_BIDS_project/derivatives/mri_clinical_matched/...

For each patched file, writes a sister `.bak` copy on the first run so
nothing is irrecoverable.

Usage:
  python clinical_pipeline/01f_apply_censoring_to_Label_Ny.py --dry-run
  python clinical_pipeline/01f_apply_censoring_to_Label_Ny.py
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
CONV_TSV = Path(
    r"D:/ADNI_SNP_Omni2.5M_20140220/conversion_labels/conversion_labels.tsv"
)
WALK_ROOTS = [
    Path(r"D:/ADNI_BIDS_project/derivatives/clinical"),
    Path(r"D:/ADNI_BIDS_project/derivatives/mri_clinical_matched"),
]

LABEL_COLS = [f"Label_{n}y" for n in range(1, 11)]


# --------------------------------------------------------------------------- #
# Per-subject Label_Ny computation (censored-aware)
# --------------------------------------------------------------------------- #
def compute_label_table(conv: pd.DataFrame) -> pd.DataFrame:
    """Return DataFrame indexed by Patient_ID with columns Label_1y..Label_10y."""
    yta = pd.to_numeric(conv["years_to_AD"], errors="coerce")
    fu = pd.to_numeric(conv["FU_years"], errors="coerce")
    out = pd.DataFrame(index=conv["Patient_ID"].values)
    for n in range(1, 11):
        col = f"Label_{n}y"
        # 1: converted within N
        # 0: not converted AND followed for at least N years
        # NaN: censored (FU < N AND no recorded AD)
        converted = yta <= n
        followed_long_enough = fu >= n
        no_convert_ever = yta.isna()
        is_censored = no_convert_ever & ~followed_long_enough
        val = np.where(
            converted.fillna(False),
            1,
            np.where(is_censored, np.nan, 0),
        )
        out[col] = val
    return out


# --------------------------------------------------------------------------- #
# Patch walker
# --------------------------------------------------------------------------- #
def csv_has_label_cols(path: Path) -> bool:
    """Quick check by reading only the header line."""
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            hdr = fh.readline()
    except OSError:
        return False
    return any(f"Label_{n}y" in hdr for n in range(1, 11))


def patch_one_csv(path: Path, label_table: pd.DataFrame,
                  backup: bool, dry_run: bool) -> tuple[int, int]:
    """Patch Label_Ny columns in one CSV. Returns (n_changed_cells, n_rows)."""
    df = pd.read_csv(path, low_memory=False)
    if "Patient_ID" not in df.columns:
        return 0, 0
    present_label_cols = [c for c in LABEL_COLS if c in df.columns]
    if not present_label_cols:
        return 0, 0
    # Map (ptid -> Label_Ny) row from the corrected table.
    corrected = label_table.reindex(df["Patient_ID"])
    n_changed = 0
    for col in present_label_cols:
        old = pd.to_numeric(df[col], errors="coerce").to_numpy()
        new = corrected[col].to_numpy()
        # Compare ignoring NaN equality
        diff_mask = ~((np.isnan(old) & np.isnan(new)) | (old == new))
        n_changed += int(diff_mask.sum())
        df[col] = new
    if dry_run or n_changed == 0:
        return n_changed, len(df)
    if backup:
        bak = path.with_suffix(path.suffix + ".bak")
        if not bak.exists():
            shutil.copy2(path, bak)
    df.to_csv(path, index=False)
    return n_changed, len(df)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--conv-tsv", type=Path, default=CONV_TSV)
    p.add_argument("--roots", type=Path, nargs="+", default=WALK_ROOTS)
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--no-backup", action="store_true",
                   help="Skip .bak side-copies (NOT recommended).")
    args = p.parse_args()

    if not args.conv_tsv.exists():
        print(f"[ERROR] conversion_labels.tsv missing: {args.conv_tsv}",
              file=sys.stderr)
        return 1

    print(f"[01f] Reading conversion labels: {args.conv_tsv}")
    conv = pd.read_csv(args.conv_tsv, sep="\t")
    print(f"  subjects in conv table: {len(conv)}")

    label_table = compute_label_table(conv)
    label_table.index.name = "Patient_ID"
    print(f"  Censoring-aware label coverage (per horizon):")
    for col in LABEL_COLS:
        ncen = int(label_table[col].isna().sum())
        npos = int((label_table[col] == 1).sum())
        nneg = int((label_table[col] == 0).sum())
        total = len(label_table)
        print(f"    {col}: pos={npos:3d}  neg={nneg:3d}  censored={ncen:3d}  "
              f"(censored rate {100*ncen/total:5.1f}% of {total})")

    print()
    # Walk + patch
    grand_n_files = 0
    grand_n_patched = 0
    grand_n_changed_cells = 0
    grand_n_rows = 0
    for root in args.roots:
        if not root.exists():
            print(f"  [skip] root missing: {root}")
            continue
        print(f"[01f] Walking: {root}")
        for csv_path in sorted(root.rglob("*.csv")):
            grand_n_files += 1
            if not csv_has_label_cols(csv_path):
                continue
            n_changed, n_rows = patch_one_csv(
                csv_path, label_table,
                backup=(not args.no_backup),
                dry_run=args.dry_run,
            )
            if n_changed > 0:
                grand_n_patched += 1
                grand_n_changed_cells += n_changed
                grand_n_rows += n_rows
                rel = csv_path.relative_to(root)
                print(f"  [{ 'DRY' if args.dry_run else 'OK' }] "
                      f"{rel}  changed_cells={n_changed}  rows={n_rows}")

    print()
    print(f"[01f] Summary:")
    print(f"  CSVs scanned         : {grand_n_files}")
    print(f"  CSVs needing changes : {grand_n_patched}")
    print(f"  Total cells changed  : {grand_n_changed_cells}")
    print(f"  Total rows patched   : {grand_n_rows}")
    if args.dry_run:
        print("  [dry-run] No files written.")
    else:
        print("  Done. Original files saved as <path>.bak (re-run is idempotent).")
        print()
        print("  NEXT STEP: scp the patched masters + splits to CSD3 so T3a-T3d sweeps")
        print("  pick up the censored-aware labels. Also commit the updated compute_prognosis")
        print("  in 01b_build_clinical_csv.py so a fresh 01b run produces matching values.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
