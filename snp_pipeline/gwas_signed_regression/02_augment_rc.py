"""
02_augment_rc.py (gwas_signed_regression)
==========================================
Augment the BMFM-DNA-SNP signed regression dataset with reverse complement
sequences.

For each existing row (EA_f / OA_f), generates the reverse complement sequence
with the same z-score (same biological effect, complementary strand).

Input
-----
  D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/bmfm_gwas_signed_regression_without_ukb/
    train.csv, dev.csv, test.csv
    Format: sequence, z_score, seq_id

Output
------
  D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/bmfm_gwas_signed_regression_without_ukb_augmented/
    train.csv, dev.csv, test.csv
    Format: sequence, z_score, seq_id
    seq_id convention:
      rs*_EA_f   — forward strand, effect allele
      rs*_EA_r   — reverse complement, effect allele
      rs*_OA_f   — forward strand, other allele
      rs*_OA_r   — reverse complement, other allele

Augmentation logic
------------------
  Reverse complement preserves the same z-score as the forward sequence,
  because it represents the same SNP on the complementary strand — the
  effect direction relative to the allele does not change.

  RC(ACGT) = ACGT reversed + A↔T, C↔G swap.

Usage
-----
  python 02_augment_rc.py
  python 02_augment_rc.py --check   # dry-run: show counts only
"""

import argparse
import io
import pathlib
import sys

import pandas as pd

# Force UTF-8 stdout on Windows
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR    = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs")
SRC_DIR     = BASE_DIR / "bmfm_gwas_signed_regression_without_ukb"
OUT_DIR     = BASE_DIR / "bmfm_gwas_signed_regression_without_ukb_augmented"
SPLITS      = ["train", "dev", "test"]

# ── CLI ───────────────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser(
    description="Augment BMFM regression dataset with reverse complement sequences.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
_p.add_argument("--src",   default=str(SRC_DIR), help="Source directory (train/dev/test.csv)")
_p.add_argument("--outdir", default=str(OUT_DIR), help="Output directory")
_p.add_argument("--check", action="store_true",
                help="Dry run — print counts only, do not write files")
args = _p.parse_args()

SRC_DIR = pathlib.Path(args.src)
OUT_DIR = pathlib.Path(args.outdir)

# ── Validate ──────────────────────────────────────────────────────────────────
if not SRC_DIR.exists():
    print(f"[ERROR] Source directory not found: {SRC_DIR}")
    sys.exit(1)

if not args.check:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Reverse complement ────────────────────────────────────────────────────────
_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMP)[::-1]


def _rename_seq_id(seq_id: str, suffix: str) -> str:
    """
    Append _f or _r to the base seq_id.

    rs12345_EA  → rs12345_EA_f  (original forward)
    rs12345_OA  → rs12345_OA_f  (original forward)
    """
    return f"{seq_id}_{suffix}"


# ── Process each split ────────────────────────────────────────────────────────
print(f"\nSource : {SRC_DIR}")
print(f"Output : {OUT_DIR}")
print(f"Dry run: {args.check}")
print()

for split in SPLITS:
    src_path = SRC_DIR / f"{split}.csv"
    if not src_path.exists():
        print(f"  [SKIP] {split}.csv not found in source directory")
        continue

    df = pd.read_csv(src_path)
    n_orig = len(df)

    # ── Rename originals: rs*_EA → rs*_EA_f ──────────────────────────────────
    df["seq_id"] = df["seq_id"].apply(lambda s: _rename_seq_id(s, "f"))

    # ── Build reverse complement rows ─────────────────────────────────────────
    rc_rows = df.copy()
    rc_rows["sequence"] = rc_rows["sequence"].apply(reverse_complement)
    rc_rows["seq_id"]   = rc_rows["seq_id"].str.replace("_f$", "_r", regex=True)
    # z_score stays the same — RC represents same effect on complementary strand

    # ── Interleave: EA_f, EA_r, OA_f, OA_r for each SNP ─────────────────────
    augmented = pd.concat([df, rc_rows], ignore_index=True)

    # Sort so EA_f / EA_r / OA_f / OA_r are grouped by rsid
    # (seq_id format: rs*_EA_f, rs*_EA_r, rs*_OA_f, rs*_OA_r)
    augmented = augmented.sort_values("seq_id").reset_index(drop=True)

    n_aug = len(augmented)

    print(f"  {split}.csv : {n_orig:>4} rows → {n_aug:>4} rows after RC augmentation")

    if args.check:
        print(f"    Sample seq_ids: {augmented['seq_id'].head(8).tolist()}")
        continue

    out_path = OUT_DIR / f"{split}.csv"
    augmented[["sequence", "z_score", "seq_id"]].to_csv(out_path, index=False)
    print(f"    Written: {out_path}")

# ── Summary ───────────────────────────────────────────────────────────────────
if not args.check:
    print(f"\nAugmented dataset written to: {OUT_DIR}")
    print("  seq_id convention:")
    print("    rs*_EA_f  forward strand, effect allele   → +z")
    print("    rs*_EA_r  reverse complement, effect allele → +z")
    print("    rs*_OA_f  forward strand, other allele    → -z")
    print("    rs*_OA_r  reverse complement, other allele  → -z")
    print("\n  Each SNP now has 4 rows instead of 2.")
    print("  Total dataset size doubled.")
