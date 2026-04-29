"""
prepare_bmfm_gwas_regression_inputs.py
=======================================
Generate BMFM-DNA-SNP fine-tuning CSVs for signed GWAS z-score regression.

For each GW-significant SNP (label=1) with a harmonised z-score:
  EA sequence  →  target = +z   (carrying the risk/effect allele)
  OA sequence  →  target = −z   (carrying the other/protective allele)

Sequence: 500 bp flanking each side (1,001 bp total, SNP centred).

Input
-----
  external_gwas_labels.tsv   (from prepare_bmfm_labels_from_external_gwas.py)
    Required columns: SNP, CHR, BP, REF, EA, OA, z_score, label

  Homo_sapiens.GRCh38.dna.primary_assembly.fa
    GRCh38 reference FASTA (plain .fa or bgzipped .fa.gz)

Output
------
  bmfm_gwas_signed_regression_without_ukb/
    train.csv           — BMFM fine-tuning format (sequence, z_score, seq_id)
    dev.csv
    test.csv
    metadata.tsv        — full record including allele, source, etc.
    finetune_config.yaml
    dataset_summary.txt

BMFM format (docs: https://github.com/BiomedSciAI/biomed-multi-omic/tree/main/run)
------------------------------------------------------------------------------------
  sequence,z_score,seq_id
  ACGTTTACCC...GGTAAG, 5.45, rs61693370_EA
  ACGTTTACCC...GGTAAG,-5.45, rs61693370_OA

Usage
-----
  python prepare_bmfm_gwas_regression_inputs.py
  python prepare_bmfm_gwas_regression_inputs.py --flank 500
  python prepare_bmfm_gwas_regression_inputs.py --max-snps 50   # quick test
"""

import argparse
import io
import pathlib
import random
import sys
import textwrap
from datetime import datetime

import pandas as pd

# Force UTF-8 stdout on Windows
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR   = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220")
LABELS_TSV = BASE_DIR / "bmfm_inputs" / "external_gwas_labels.tsv"
FASTA      = BASE_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
OUT_SUBDIR = "bmfm_gwas_signed_regression_without_ukb"
OUT_DIR    = BASE_DIR / "bmfm_inputs" / OUT_SUBDIR

MAX_MODEL_LEN = 8192   # ModernBERT hard limit

# ── CLI ───────────────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser(
    description="Prepare BMFM-DNA-SNP signed z-score regression fine-tuning inputs.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
_p.add_argument("--labels", default=str(LABELS_TSV),
                help="Path to external_gwas_labels.tsv")
_p.add_argument("--fasta",  default=str(FASTA),
                help="Path to GRCh38 reference FASTA (.fa or bgzipped .fa.gz)")
_p.add_argument("--outdir", default=str(OUT_DIR),
                help="Output directory")
_p.add_argument("--flank", type=int, default=500,
                help="Flanking bases each side of the SNP (default: 500 → 1,001 bp)")
_p.add_argument("--split", nargs=3, type=float, default=[0.8, 0.1, 0.1],
                metavar=("TRAIN", "DEV", "TEST"),
                help="Train/dev/test proportions (must sum to 1; default: 0.8 0.1 0.1)")
_p.add_argument("--seed", type=int, default=42,
                help="Random seed for SNP-level shuffle (default: 42)")
_p.add_argument("--max-snps", type=int, default=None, dest="max_snps",
                help="Cap at first N SNPs after filtering (for quick testing).")
args = _p.parse_args()

LABELS_TSV = pathlib.Path(args.labels)
FASTA      = pathlib.Path(args.fasta)
OUT_DIR    = pathlib.Path(args.outdir)
FLANK      = args.flank

# ── Validate ──────────────────────────────────────────────────────────────────
if not LABELS_TSV.exists():
    print(f"[ERROR] Labels TSV not found: {LABELS_TSV}")
    sys.exit(1)
if not FASTA.exists():
    print(f"[ERROR] FASTA not found: {FASTA}")
    sys.exit(1)
if FLANK * 2 + 1 > MAX_MODEL_LEN:
    print(f"[WARNING] Flank {FLANK} bp => window {FLANK*2+1} bp exceeds "
          f"ModernBERT limit ({MAX_MODEL_LEN}). Capping flank.")
    FLANK = MAX_MODEL_LEN // 2 - 1

OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── pyfaidx ───────────────────────────────────────────────────────────────────
try:
    from pyfaidx import Fasta, FastaIndexingError
except ImportError:
    print("[ERROR] pyfaidx not installed. Run:  pip install pyfaidx")
    sys.exit(1)

# ── Load and filter labels ────────────────────────────────────────────────────
print(f"\nLoading labels: {LABELS_TSV.name}")
labels = pd.read_csv(LABELS_TSV, sep="\t", low_memory=False)
print(f"  Total rows       : {len(labels):,}")

sig = labels[labels["label"] == 1].copy()
print(f"  label=1 (GW-sig) : {len(sig):,}")

# If OA is missing, assume the other allele is the GRCh38 reference (REF = A2 after refcorr)
n_oa_imputed = sig["OA"].isna().sum()
if n_oa_imputed:
    sig.loc[sig["OA"].isna(), "OA"] = sig.loc[sig["OA"].isna(), "REF"]
    print(f"  OA imputed from REF (missing OA → REF) : {n_oa_imputed:,}")

before = len(sig)
dropped_mask = sig["EA"].isna() | sig["OA"].isna() | sig["z_score"].isna()
dropped_snps = sig.loc[dropped_mask, "SNP"].tolist()
sig = sig[~dropped_mask].copy()
if dropped_snps:
    print(f"  Dropped (missing EA / OA / z_score) : {len(dropped_snps):,}")
    for rsid in dropped_snps:
        row = labels.loc[labels["SNP"] == rsid, ["EA", "OA", "z_score"]].iloc[0]
        print(f"    {rsid:<20}  EA={row['EA']}  OA={row['OA']}  z_score={row['z_score']}")
print(f"  Usable SNPs for regression          : {len(sig):,}")

if len(sig) == 0:
    print("[ERROR] No usable SNPs. Run prepare_bmfm_labels_from_external_gwas.py first.")
    sys.exit(1)

if args.max_snps:
    sig = sig.head(args.max_snps).copy()
    print(f"  (subset: first {args.max_snps} SNPs for testing)")

# ── Open FASTA ────────────────────────────────────────────────────────────────
print(f"\nOpening FASTA: {FASTA.name}")
print("  (pyfaidx will build a .fai index on first use — may take a few minutes)")
try:
    fa = Fasta(str(FASTA), rebuild=False)
except FastaIndexingError:
    print("  .fai not found; building index now...")
    fa = Fasta(str(FASTA), rebuild=True)
except Exception as e:
    print(f"[ERROR] Could not open FASTA: {e}")
    sys.exit(1)

chrom_names = set(fa.keys())
print(f"  Chromosomes in FASTA: {sorted(chrom_names)[:5]} ...")


def resolve_chrom(chrom_val: str) -> str:
    """Map BIM integer chromosome to FASTA contig name (handles 'chr' prefix)."""
    for candidate in [str(chrom_val), f"chr{chrom_val}"]:
        if candidate in chrom_names:
            return candidate
    return None


# ── Extract allele-specific sequences ─────────────────────────────────────────
print(f"\nExtracting {FLANK} bp flanks ({FLANK*2+1} bp windows) ...")
print("  EA row → +z   |   OA row → −z")

records = []
skipped = 0
REPORT_EVERY = 50

for i, row in enumerate(sig.itertuples(index=False), 1):
    if i % REPORT_EVERY == 0 or i == len(sig):
        print(f"  {i:,} / {len(sig):,} SNPs processed ...")

    chrom = resolve_chrom(str(row.CHR))
    if chrom is None:
        print(f"  [SKIP] {row.SNP}: chromosome '{row.CHR}' not in FASTA")
        skipped += 1
        continue

    ea = str(row.EA).upper().strip()
    oa = str(row.OA).upper().strip()
    z  = float(row.z_score)

    if not ea or ea in {"NA", "NAN", "NONE", "<NA>"} \
            or not oa or oa in {"NA", "NAN", "NONE", "<NA>"}:
        skipped += 1
        continue

    pos1      = int(row.BP)          # 1-based GRCh38 coordinate from BIM
    pos0      = pos1 - 1             # 0-based for slicing
    chrom_len = len(fa[chrom])

    start = max(0, pos0 - FLANK)
    end   = min(chrom_len, pos0 + FLANK + 1)
    seq   = str(fa[chrom][start:end]).upper()

    snp_in_win = pos0 - start        # index of the SNP within the extracted window

    def _make_seq(allele: str) -> str:
        s = list(seq)
        if 0 <= snp_in_win < len(s):
            s[snp_in_win] = allele.upper()
        return "".join(s)

    ea_seq = _make_seq(ea)
    oa_seq = _make_seq(oa)

    base = {
        "rsid"       : row.SNP,
        "chrom"      : row.CHR,
        "pos_hg38"   : pos1,
        "REF"        : row.REF,
        "source"     : row.source,
        "z_gwas"     : z,
    }

    # Effect-allele sequence → target +z
    records.append({**base,
        "sequence"   : ea_seq,
        "z_score"    : z,
        "seq_id"     : f"{row.SNP}_EA",
        "allele"     : ea,
        "allele_type": "EA",
    })

    # Other-allele sequence → target −z
    records.append({**base,
        "sequence"   : oa_seq,
        "z_score"    : -z,
        "seq_id"     : f"{row.SNP}_OA",
        "allele"     : oa,
        "allele_type": "OA",
    })

n_snp_pairs = len(records) // 2
print(f"\n  Sequences generated : {len(records):,}  ({n_snp_pairs:,} EA+OA pairs)")
print(f"  Skipped             : {skipped:,}")

if not records:
    print("[ERROR] No sequences generated.")
    sys.exit(1)

# ── DataFrame + metadata ──────────────────────────────────────────────────────
df = pd.DataFrame(records)

meta_path = OUT_DIR / "metadata.tsv"
df.to_csv(meta_path, sep="\t", index=False)
print(f"  Metadata TSV        : {meta_path}")

# ── SNP-level train/dev/test split ────────────────────────────────────────────
# Shuffle at the SNP level so EA and OA pairs always end up in the same split.
train_frac, dev_frac, test_frac = args.split
assert abs(sum(args.split) - 1.0) < 1e-6, "Split fractions must sum to 1.0"

snp_list = df["rsid"].unique().tolist()
random.seed(args.seed)
random.shuffle(snp_list)

n_snps  = len(snp_list)
n_train = int(n_snps * train_frac)
n_dev   = int(n_snps * dev_frac)

train_snps = set(snp_list[:n_train])
dev_snps   = set(snp_list[n_train:n_train + n_dev])
test_snps  = set(snp_list[n_train + n_dev:])

train_df = df[df["rsid"].isin(train_snps)].copy()
dev_df   = df[df["rsid"].isin(dev_snps)].copy()
test_df  = df[df["rsid"].isin(test_snps)].copy()

# BMFM-DNA-SNP fine-tuning CSV: sequence first, then label column(s), then seq_id
CSV_COLS = ["sequence", "z_score", "seq_id"]


def write_split(subset: pd.DataFrame, path: pathlib.Path):
    subset[CSV_COLS].to_csv(path, index=False)
    print(f"  {path.name:<12}: {len(subset):,} rows  "
          f"({subset['rsid'].nunique():,} SNP pairs)")


print(f"\nWriting train/dev/test split "
      f"({train_frac:.0%} / {dev_frac:.0%} / {test_frac:.0%}):")
write_split(train_df, OUT_DIR / "train.csv")
write_split(dev_df,   OUT_DIR / "dev.csv")
write_split(test_df,  OUT_DIR / "test.csv")

# ── Fine-tuning YAML ──────────────────────────────────────────────────────────
yaml_path = OUT_DIR / "finetune_config.yaml"
yaml_content = textwrap.dedent(f"""\
    # BMFM-DNA-SNP fine-tuning config — signed GWAS z-score regression
    # Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}
    #
    # The z_score column encodes signed effect size (harmonised to BIM A1/ALT):
    #   EA sequence  →  +z  (effect allele, higher AD risk)
    #   OA sequence  →  −z  (other allele)
    #
    # Usage:
    #   export INPUT_DIRECTORY="{OUT_DIR}"
    #   bmfm-targets-run -cn dna_finetune_train_and_test_config \\
    #       input_directory=$INPUT_DIRECTORY \\
    #       output_directory=/tmp/bmfm_gwas_regression \\
    #       checkpoint=ibm-research/biomed.dna.snp.modernbert.113m.v1

    label_columns:
      - _target_: bmfm_targets.config.LabelColumnInfo
        label_column_name: "z_score"
        is_regression_label: true

    data_module:
      defaults: dna_base_seq_cls
      max_length: {FLANK * 2 + 1}
      dataset_kwargs:
        processed_data_source: ${{input_directory}}

    trainer:
      learning_rate: 1.0e-5
      losses:
        - name: mse
          label_column_name: z_score
""")
yaml_path.write_text(yaml_content, encoding="utf-8")
print(f"\n  Fine-tuning YAML : {yaml_path}")

# ── Summary ───────────────────────────────────────────────────────────────────
SEP = "=" * 65
summary = "\n".join([
    SEP,
    "  BMFM Signed Regression Dataset — GWAS z-score",
    f"  Generated : {datetime.now().strftime('%Y-%m-%d %H:%M')}",
    SEP,
    "",
    f"  GW-sig SNPs with z-score  : {n_snp_pairs:,}",
    f"  Sequences generated       : {len(records):,}  (EA + OA per SNP)",
    f"  Skipped (FASTA/allele)    : {skipped:,}",
    f"  Window size               : {FLANK*2+1} bp  ({FLANK} bp each side of SNP)",
    "",
    f"  Split ({train_frac:.0%} / {dev_frac:.0%} / {test_frac:.0%}):",
    f"    train : {len(train_df):,} rows  ({train_df['rsid'].nunique():,} SNP pairs)",
    f"    dev   : {len(dev_df):,} rows  ({dev_df['rsid'].nunique():,} SNP pairs)",
    f"    test  : {len(test_df):,} rows  ({test_df['rsid'].nunique():,} SNP pairs)",
    "",
    f"  Output : {OUT_DIR}",
    "",
    "  To fine-tune (from Linux/HPC):",
    f'    export INPUT_DIRECTORY="{OUT_DIR}"',
    "    bmfm-targets-run -cn dna_finetune_train_and_test_config \\",
    "        input_directory=$INPUT_DIRECTORY \\",
    "        output_directory=/tmp/bmfm_gwas_regression \\",
    "        checkpoint=ibm-research/biomed.dna.snp.modernbert.113m.v1",
    SEP,
])
print("\n" + summary)

summary_path = OUT_DIR / "dataset_summary.txt"
summary_path.write_text(summary, encoding="utf-8")
print(f"Summary written : {summary_path}")
