"""
prepare_bmfm_dna_snp_inputs.py
===============================
Extract GRCh38 reference sequences around each SNP in the ADNI dataset and
write CSV files ready for BMFM-DNA-SNP inference.

Input format required by bmfm-targets-run dna_predict
------------------------------------------------------
A directory with train.csv / dev.csv / test.csv (or a single all.csv for
inference only), where the first column is 'sequence' (raw ACGT string).

  sequence,seq_id
  ACGTTTACCC...GGTAAG,rs12345_REF
  ACGTTTACCC...GGTAAG,rs12345_ALT

IMPORTANT — sequence length
---------------------------
BMFM-DNA-SNP (ModernBERT) has a maximum context of 8,192 tokens.
With single-nucleotide tokenisation that equals 8,192 bp.
The default flank is 500 bp either side (1,001 bp total).
Using 500 KB flanks would exceed the model's capacity by ~60x.

Dependencies
------------
  pip install pyfaidx
  (pandas, pathlib, etc. are already in the snp conda environment)

Usage
-----
  python prepare_bmfm_dna_snp_inputs.py
  python prepare_bmfm_dna_snp_inputs.py --flank 1000
  python prepare_bmfm_dna_snp_inputs.py --inference-only
  python prepare_bmfm_dna_snp_inputs.py --max-snps 5000   # quick test
"""

import argparse
import io
import pathlib
import sys

import pandas as pd

# Force UTF-8 stdout on Windows
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ── Configuration ─────────────────────────────────────────────────────────────
BASE_DIR  = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220")
BIM_FILE  = BASE_DIR / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr.bim"
FASTA     = BASE_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
OUT_DIR   = BASE_DIR / "bmfm_inputs"

MAX_MODEL_LEN = 8192   # ModernBERT hard limit (tokens ≈ bases)

# ── CLI ───────────────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser(
    description="Prepare BMFM-DNA-SNP inference CSV inputs from ADNI GRCh38 BIM.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
_p.add_argument("--bim",   default=str(BIM_FILE),
                help="Path to GRCh38 PLINK BIM file")
_p.add_argument("--fasta", default=str(FASTA),
                help="Path to GRCh38 reference FASTA (can be .fa.gz if bgzipped)")
_p.add_argument("--outdir", default=str(OUT_DIR),
                help="Output directory for CSV files")
_p.add_argument("--flank", type=int, default=500,
                help="Bases to extract on each side of the SNP (default: 500). "
                     f"Max useful = {MAX_MODEL_LEN // 2}.")
_p.add_argument("--include-alt", action="store_true", dest="include_alt",
                help="Also write an ALT-allele sequence for each SNP.")
_p.add_argument("--inference-only", action="store_true", dest="inference_only",
                help="Write a single all.csv instead of train/dev/test split.")
_p.add_argument("--split", nargs=3, type=float, default=[0.8, 0.1, 0.1],
                metavar=("TRAIN", "DEV", "TEST"),
                help="Train/dev/test proportions (default: 0.8 0.1 0.1)")
_p.add_argument("--max-snps", type=int, default=None, dest="max_snps",
                help="Limit to first N SNPs (for quick testing).")
_p.add_argument("--chunksize", type=int, default=50_000, dest="chunksize",
                help="Write CSVs in chunks of N rows (reduces peak memory).")
args = _p.parse_args()

BIM_FILE  = pathlib.Path(args.bim)
FASTA     = pathlib.Path(args.fasta)
OUT_DIR   = pathlib.Path(args.outdir)
FLANK     = args.flank

# ── Validate ──────────────────────────────────────────────────────────────────
if not BIM_FILE.exists():
    print(f"[ERROR] BIM not found: {BIM_FILE}")
    sys.exit(1)
if not FASTA.exists():
    print(f"[ERROR] FASTA not found: {FASTA}")
    sys.exit(1)
if FLANK * 2 + 1 > MAX_MODEL_LEN:
    print(f"[WARNING] Flank {FLANK} bp => window {FLANK*2+1} bp exceeds "
          f"ModernBERT max context ({MAX_MODEL_LEN} tokens). Capping flank "
          f"at {MAX_MODEL_LEN // 2 - 1} bp.")
    FLANK = MAX_MODEL_LEN // 2 - 1

OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Import pyfaidx (with helpful error) ───────────────────────────────────────
try:
    from pyfaidx import Fasta, FastaIndexingError
except ImportError:
    print("[ERROR] pyfaidx not installed. Run:")
    print("  pip install pyfaidx")
    sys.exit(1)

# ── Load FASTA index ──────────────────────────────────────────────────────────
print(f"\nOpening FASTA: {FASTA}")
print("  (pyfaidx will build a .fai index on first use — this may take a few "
      "minutes for a 3 GB file)")
try:
    fa = Fasta(str(FASTA), rebuild=False)
except FastaIndexingError:
    print("  .fai not found; building index now...")
    fa = Fasta(str(FASTA), rebuild=True)
except Exception as e:
    print(f"[ERROR] Could not open FASTA: {e}")
    print()
    print("  If the file is regular gzip (not bgzipped), pyfaidx cannot index it.")
    print("  Convert with bgzip:")
    print(f"    gunzip -c {FASTA} | bgzip -c > {FASTA.with_suffix('').with_suffix('.bgz.fa.gz')}")
    print("  Then re-run with --fasta pointing to the new file.")
    sys.exit(1)

chrom_names = set(fa.keys())
print(f"  Chromosomes in FASTA: {sorted(chrom_names)[:5]} ...")

def resolve_chrom(chrom_int: str) -> str:
    """Map BIM integer chromosome to FASTA contig name."""
    for candidate in [chrom_int, f"chr{chrom_int}"]:
        if candidate in chrom_names:
            return candidate
    return None

# ── Load BIM ──────────────────────────────────────────────────────────────────
print(f"\nLoading BIM: {BIM_FILE}")
bim = pd.read_csv(BIM_FILE, sep="\t", header=None,
                  names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
                  dtype={"CHR": str, "SNP": str, "CM": str, "BP": int,
                         "A1": str, "A2": str})

if args.max_snps:
    bim = bim.head(args.max_snps)
    print(f"  (subset: first {args.max_snps} SNPs only)")

total = len(bim)
print(f"  {total:,} variants loaded")

# ── Extract sequences ─────────────────────────────────────────────────────────
print(f"\nExtracting {FLANK} bp flanks ({FLANK*2+1} bp windows) ...")
print("  Writing REF sequences" + (" + ALT sequences" if args.include_alt else ""))

records = []
skipped = 0
REPORT_EVERY = 100_000

for i, row in enumerate(bim.itertuples(index=False), 1):
    if i % REPORT_EVERY == 0:
        print(f"  {i:,} / {total:,} ({100*i/total:.1f}%) ...")

    chrom = resolve_chrom(str(row.CHR))
    if chrom is None:
        skipped += 1
        continue

    pos1 = int(row.BP)          # 1-based in BIM
    pos0 = pos1 - 1             # 0-based for slicing
    chrom_len = len(fa[chrom])

    start = max(0, pos0 - FLANK)
    end   = min(chrom_len, pos0 + FLANK + 1)

    seq = str(fa[chrom][start:end]).upper()

    # Position of the SNP within the extracted window
    snp_in_window = pos0 - start

    # REF sequence (reference allele at SNP position)
    ref_allele = str(row.A2).upper()   # A2 is REF in PLINK2 --ref-from-fa output
    ref_seq = list(seq)
    if 0 <= snp_in_window < len(ref_seq):
        ref_seq[snp_in_window] = ref_allele
    ref_seq_str = "".join(ref_seq)

    records.append({
        "sequence": ref_seq_str,
        "seq_id":   f"{row.SNP}_REF",
        "rsid":     row.SNP,
        "chrom":    row.CHR,
        "pos":      pos1,
        "allele":   ref_allele,
        "allele_type": "REF",
    })

    if args.include_alt:
        alt_allele = str(row.A1).upper()   # A1 is ALT in PLINK2 output
        alt_seq = list(seq)
        if 0 <= snp_in_window < len(alt_seq):
            alt_seq[snp_in_window] = alt_allele
        alt_seq_str = "".join(alt_seq)

        records.append({
            "sequence": alt_seq_str,
            "seq_id":   f"{row.SNP}_ALT",
            "rsid":     row.SNP,
            "chrom":    row.CHR,
            "pos":      pos1,
            "allele":   alt_allele,
            "allele_type": "ALT",
        })

print(f"\n  Done. {len(records):,} sequences extracted. Skipped {skipped:,} "
      f"(unresolved chromosomes).")

# ── Build DataFrame ───────────────────────────────────────────────────────────
df = pd.DataFrame(records)

# Save full metadata TSV for downstream use
meta_path = OUT_DIR / "snp_sequence_metadata.tsv"
df.to_csv(meta_path, sep="\t", index=False)
print(f"\n  Metadata TSV: {meta_path}")

# ── Write CSV files ───────────────────────────────────────────────────────────
# BMFM-DNA-SNP expects: sequence column first, then optional labels/ids

def write_csv(subset: pd.DataFrame, path: pathlib.Path):
    """Write BMFM-compatible CSV (sequence, seq_id columns only)."""
    out = subset[["sequence", "seq_id"]].copy()
    out.to_csv(path, index=False)
    print(f"  {path.name}: {len(out):,} rows")

if args.inference_only:
    write_csv(df, OUT_DIR / "all.csv")
    print(f"\n  For inference, run:")
    print(f'    export INPUT_DIRECTORY="{OUT_DIR}"')
    print(f"    bmfm-targets-run -cn dna_predict \\")
    print(f"        input_directory=$INPUT_DIRECTORY \\")
    print(f"        working_dir=/tmp \\")
    print(f"        checkpoint=ibm-research/biomed.dna.snp.modernbert.113m.v1")
else:
    train_frac, dev_frac, test_frac = args.split
    assert abs(sum(args.split) - 1.0) < 1e-6, "Split fractions must sum to 1.0"

    n = len(df)
    n_train = int(n * train_frac)
    n_dev   = int(n * dev_frac)

    # Shuffle with fixed seed for reproducibility
    df_shuffled = df.sample(frac=1, random_state=42).reset_index(drop=True)
    train_df = df_shuffled.iloc[:n_train]
    dev_df   = df_shuffled.iloc[n_train:n_train + n_dev]
    test_df  = df_shuffled.iloc[n_train + n_dev:]

    print(f"\nWriting train/dev/test split ({train_frac:.0%}/{dev_frac:.0%}/{test_frac:.0%}):")
    write_csv(train_df, OUT_DIR / "train.csv")
    write_csv(dev_df,   OUT_DIR / "dev.csv")
    write_csv(test_df,  OUT_DIR / "test.csv")

    # Save split assignments to metadata
    split_col = ["train"] * n_train + ["dev"] * n_dev + ["test"] * (n - n_train - n_dev)
    df_shuffled["split"] = split_col
    df_shuffled.to_csv(meta_path, sep="\t", index=False)

    print(f"\n  For inference/fine-tuning, run:")
    print(f'    export INPUT_DIRECTORY="{OUT_DIR}"')
    print(f"    bmfm-targets-run -cn dna_predict \\")
    print(f"        input_directory=$INPUT_DIRECTORY \\")
    print(f"        working_dir=/tmp \\")
    print(f"        checkpoint=ibm-research/biomed.dna.snp.modernbert.113m.v1")

print(f"\nDone! Output directory: {OUT_DIR}")
print(f"  Window size : {FLANK*2+1} bp ({FLANK} bp each side of SNP)")
print(f"  Total seqs  : {len(df):,}")
print(f"  Metadata    : {meta_path}")
