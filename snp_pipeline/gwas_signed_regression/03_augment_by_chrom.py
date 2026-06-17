"""
03_augment_by_chrom.py (gwas_signed_regression)
================================================
Generate chromosome-level multi-SNP sequences (up to ~8,192 bp) for
BMFM-DNA-SNP signed regression fine-tuning.

Rationale
---------
Sequence context matters for DNA language models. Rather than 1,001 bp windows
centred on each SNP independently (08c), this script groups all GW-sig SNPs on
the same chromosome into contiguous 8 kb windows, preserving their true genomic
spacing. Each SNP within a window is the "target" SNP in turn: its allele is
substituted (EA or OA), while all other positions in the sequence remain the
GRCh38 reference. This exposes the model to multi-SNP context within a single
forward pass.

Algorithm
---------
  For each chromosome:
    1. Sort GW-sig SNPs by position.
    2. Greedily pack consecutive SNPs into windows ≤ MAX_LEN bp.
       (A new window starts when adding the next SNP would exceed MAX_LEN.)
    3. For each window (N SNPs, positions p1 … pN):
       - Extract reference sequence: [p1 − left_flank … pN + right_flank]
         where flanks are set so total length = MAX_LEN (split symmetrically).
       - For each target SNP k (k = 0 … N−1):
           EA row: substitute position k → EA, all other SNP positions = REF
           OA row: substitute position k → OA, all other SNP positions = REF
           z_score = +z_k  /  −z_k
           seq_id  = chrC_w{win}_rs*_EA_f  /  _OA_f

Train/dev/test split
--------------------
  Inherited from the original 08c dataset:
  each target SNP's seq_id is assigned to the same split as in train/dev/test.csv.

Output
------
  bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom/
    forward/
      train.csv, dev.csv, test.csv   (forward strand only)
    forward_and_reverse/
      train.csv, dev.csv, test.csv   (forward + reverse complement)

Usage
-----
  python 03_augment_by_chrom.py
  python 03_augment_by_chrom.py --max-len 4096  # smaller windows
  python 03_augment_by_chrom.py --min-flank 20  # left-flank window
  python 03_augment_by_chrom.py --check          # dry run
"""

import argparse
import io
import pathlib
import sys
from collections import defaultdict

import pandas as pd

# Force UTF-8 stdout on Windows
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR   = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220")
LABELS_TSV = BASE_DIR / "bmfm_inputs" / "external_gwas_labels.tsv"
FASTA      = BASE_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
SRC_DIR    = BASE_DIR / "bmfm_inputs" / "bmfm_gwas_signed_regression_without_ukb"
OUT_DIR    = BASE_DIR / "bmfm_inputs" / "bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom"

SPLITS     = ["train", "dev", "test"]

# ── CLI ───────────────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser(
    description="Multi-SNP chromosome-level 8 kb context augmentation for BMFM regression.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
_p.add_argument("--labels",  default=str(LABELS_TSV))
_p.add_argument("--fasta",   default=str(FASTA))
_p.add_argument("--src",     default=str(SRC_DIR),
                help="Original 08c forward dataset (for split assignment)")
_p.add_argument("--outdir",  default=str(OUT_DIR))
_p.add_argument("--max-len",   type=int, default=8192, dest="max_len",
                help="Maximum sequence length in bp (default: 8192)")
_p.add_argument("--min-flank", type=int, default=50,  dest="min_flank",
                help="Minimum bp flanking each end of the window (default: 20)")
_p.add_argument("--check",   action="store_true",
                help="Dry run — print window statistics, do not write files")
args = _p.parse_args()

LABELS_TSV = pathlib.Path(args.labels)
FASTA      = pathlib.Path(args.fasta)
SRC_DIR    = pathlib.Path(args.src)
OUT_DIR    = pathlib.Path(args.outdir)
MAX_LEN   = args.max_len
MIN_FLANK = args.min_flank

if 2 * MIN_FLANK >= MAX_LEN:
    print(f"[ERROR] min_flank ({MIN_FLANK}) * 2 >= max_len ({MAX_LEN})")
    sys.exit(1)

MAX_INNER_SPAN = MAX_LEN - 2 * MIN_FLANK  # max bp between first and last SNP in a window

# ── Validate ──────────────────────────────────────────────────────────────────
for p, name in [(LABELS_TSV, "Labels TSV"), (FASTA, "FASTA"), (SRC_DIR, "Source dir")]:
    if not p.exists():
        print(f"[ERROR] {name} not found: {p}")
        sys.exit(1)

FWD_DIR = OUT_DIR / "forward"
RC_DIR  = OUT_DIR / "forward_and_reverse"
if not args.check:
    FWD_DIR.mkdir(parents=True, exist_ok=True)
    RC_DIR.mkdir(parents=True, exist_ok=True)

# ── Helpers ───────────────────────────────────────────────────────────────────
_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _make_seq(ref_seq: str, snp_positions_in_window: list, target_idx: int,
              allele: str) -> str:
    """
    Substitute target SNP with `allele`; all other positions are REF (already
    present in ref_seq from FASTA). Returns the modified sequence string.

    snp_positions_in_window : list of 0-based offsets within ref_seq for each SNP
    target_idx              : index into that list for the focal SNP
    """
    s = list(ref_seq)
    tgt_offset = snp_positions_in_window[target_idx]
    if 0 <= tgt_offset < len(s):
        s[tgt_offset] = allele.upper()
    return "".join(s)


def _make_seq_all(ref_seq: str, snp_offsets: list, alleles: list) -> str:
    """
    Substitute ALL SNPs in the window simultaneously with their respective alleles.
    snp_offsets : list of 0-based offsets within ref_seq
    alleles     : list of allele strings, one per SNP (same order as snp_offsets)
    """
    s = list(ref_seq)
    for offset, allele in zip(snp_offsets, alleles):
        if 0 <= offset < len(s):
            s[offset] = allele.upper()
    return "".join(s)


# ── Load split assignments from original 08c dataset ─────────────────────────
print(f"\nLoading split assignments from: {SRC_DIR}")
snp_to_split = {}   # rsid → 'train' | 'dev' | 'test'
for split in SPLITS:
    path = SRC_DIR / f"{split}.csv"
    if not path.exists():
        print(f"  [WARN] {split}.csv not found in source directory — skipping")
        continue
    df = pd.read_csv(path)
    for seq_id in df["seq_id"]:
        # seq_id format from 08c: rs*_EA  or rs*_OA
        rsid = seq_id.rsplit("_", 1)[0]   # strip _EA / _OA suffix
        snp_to_split[rsid] = split
print(f"  SNPs with split assignment: {len(snp_to_split):,}")

# ── Load GW-sig SNPs ──────────────────────────────────────────────────────────
print(f"\nLoading labels: {LABELS_TSV.name}")
labels = pd.read_csv(LABELS_TSV, sep="\t", low_memory=False)
sig = labels[labels["label"] == 1].copy()

# OA imputation from REF (mirrors 08c logic)
oa_missing = sig["OA"].isna()
if oa_missing.any():
    sig.loc[oa_missing, "OA"] = sig.loc[oa_missing, "REF"]

sig = sig.dropna(subset=["EA", "OA", "z_score"]).copy()
sig["CHR"] = sig["CHR"].astype(str)
sig = sig.sort_values(["CHR", "BP"]).reset_index(drop=True)
print(f"  GW-sig SNPs with z_score : {len(sig):,}")
print(f"  Chromosomes              : {sorted(sig['CHR'].unique())}")

# ── Open FASTA ────────────────────────────────────────────────────────────────
try:
    from pyfaidx import Fasta, FastaIndexingError
except ImportError:
    print("[ERROR] pyfaidx not installed. Run:  pip install pyfaidx")
    sys.exit(1)

print(f"\nOpening FASTA: {FASTA.name}")
try:
    fa = Fasta(str(FASTA), rebuild=False)
except FastaIndexingError:
    print("  Building .fai index (first use — may take a few minutes) ...")
    fa = Fasta(str(FASTA), rebuild=True)

chrom_names = set(fa.keys())

def resolve_chrom(chrom_val: str) -> str:
    for c in [str(chrom_val), f"chr{chrom_val}"]:
        if c in chrom_names:
            return c
    return None


# ── Window packing ────────────────────────────────────────────────────────────
def pack_windows(snp_list, max_inner_span: int) -> list:
    """
    Greedy packing of SNPs into non-overlapping windows.
    A new window starts when the span from the first SNP to the next SNP
    would exceed max_inner_span (= MAX_LEN - 2 * MIN_FLANK).
    This guarantees at least MIN_FLANK bp of flank on each side.
    """
    windows = []
    current = [snp_list[0]]
    for snp in snp_list[1:]:
        span = snp.BP - current[0].BP
        if span >= max_inner_span:     # would leave < MIN_FLANK on each side
            windows.append(current)
            current = [snp]
        else:
            current.append(snp)
    windows.append(current)
    return windows


# ── Build records ─────────────────────────────────────────────────────────────
print(f"\nBuilding multi-SNP {MAX_LEN} bp windows (min flank: {MIN_FLANK} bp) ...")

records_by_split = defaultdict(list)   # split → list of dicts
win_stats = defaultdict(int)           # chr → n_windows
skipped = 0

for chrom, grp in sig.groupby("CHR"):
    chrom_fa = resolve_chrom(chrom)
    if chrom_fa is None:
        print(f"  [SKIP] CHR {chrom} not in FASTA")
        skipped += grp.shape[0]
        continue

    snp_rows = list(grp.itertuples(index=False))
    windows  = pack_windows(snp_rows, MAX_INNER_SPAN)
    win_stats[chrom] = len(windows)

    chrom_len = len(fa[chrom_fa])

    for win_idx, window in enumerate(windows):
        n_snps = len(window)
        first_bp = window[0].BP       # 1-based
        last_bp  = window[-1].BP      # 1-based

        # Compute flanks ensuring MIN_FLANK on each side
        inner_span  = last_bp - first_bp           # bp between first and last SNP
        available   = MAX_LEN - inner_span - 2 * MIN_FLANK  # extra flank to distribute
        left_flank  = MIN_FLANK + available // 2
        right_flank = MIN_FLANK + (available - available // 2)

        win_start0 = max(0, first_bp - 1 - left_flank)   # 0-based
        win_end0   = min(chrom_len, last_bp + right_flank) # exclusive

        ref_seq = str(fa[chrom_fa][win_start0:win_end0]).upper()
        actual_len = len(ref_seq)

        # 0-based offsets of each SNP within ref_seq
        snp_offsets = [(row.BP - 1 - win_start0) for row in window]

        if args.check:
            # Only gather stats (2 rows per window: allEA + allOA)
            splits_in_win = [snp_to_split.get(row.SNP, "train") for row in window]
            split = max(set(splits_in_win), key=splits_in_win.count)  # majority vote
            records_by_split[split].append(1)
            records_by_split[split].append(1)  # EA + OA
            continue

        # ── All SNPs in window substituted simultaneously ──────────────────────
        eas      = [str(row.EA).upper() for row in window]
        oas      = [str(row.OA).upper() for row in window]
        z_total  = sum(float(row.z_score) for row in window)   # additive z-score
        rsids_str = "|".join(row.SNP for row in window)
        pos_str   = "|".join(str(row.BP) for row in window)
        source_str = "|".join(str(row.source) for row in window)

        # Split assignment: majority vote across SNPs in window
        splits_in_win = [snp_to_split.get(row.SNP, "train") for row in window]
        split = max(set(splits_in_win), key=splits_in_win.count)

        win_label = f"chr{chrom}_w{win_idx:03d}"

        # ── Per-SNP metadata (suffixed by _snp1, _snp2, ...) ──────────────────
        for alleles, sign, allele_tag in [(eas, z_total, "allEA"),
                                          (oas, -z_total, "allOA")]:
            seq = _make_seq_all(ref_seq, snp_offsets, alleles)

            snp_type = allele_tag.replace("all", "")   # "EA" or "OA"

            record = {
                "sequence"        : seq,
                "z_score"         : sign,
                "seq_id"          : f"{win_label}_{allele_tag}_f",
                "window"          : win_label,
                "chrom"           : chrom,
                "n_snps_in_window": n_snps,
                "win_len_bp"      : actual_len,
            }

            for k, row in enumerate(window):
                suf = f"_snp{k + 1}"
                raw_z = float(row.z_score)
                record[f"rsid{suf}"]           = row.SNP
                record[f"position_hg38{suf}"]  = row.BP
                record[f"allele_type{suf}"]     = snp_type
                record[f"z_score{suf}"]         = raw_z if snp_type == "EA" else -raw_z

            records_by_split[split].append(record)

# ── Report windows ────────────────────────────────────────────────────────────
total_windows = sum(win_stats.values())
print(f"\n  Windows created (≤{MAX_LEN} bp each):")
for chrom in sorted(win_stats, key=lambda x: (len(x), x)):
    print(f"    CHR {chrom:>3}: {win_stats[chrom]:>3} windows")
print(f"    Total    : {total_windows:>3} windows")

if args.check:
    print("\n  Dry-run row counts (approximate):")
    for split in SPLITS:
        print(f"    {split}: {len(records_by_split.get(split, [])):>5} rows")
    print("\n(Dry run — no files written)")
    sys.exit(0)

# ── Write forward strand CSVs ─────────────────────────────────────────────────
print(f"\nWriting forward-strand CSVs → {FWD_DIR}")

# Fixed prefix columns; per-SNP columns appended dynamically
_PREFIX_COLS = ["sequence", "z_score", "seq_id",
                "window", "chrom", "n_snps_in_window", "win_len_bp"]

def _ordered_cols(df: pd.DataFrame) -> list:
    """Return prefix columns first, then per-SNP columns sorted naturally."""
    snp_cols = sorted(
        [c for c in df.columns if c not in _PREFIX_COLS],
        key=lambda c: (c.rsplit("_snp", 1)[-1].zfill(4),  # snp number
                       c.rsplit("_snp", 1)[0])             # field name
    )
    return _PREFIX_COLS + snp_cols

for split in SPLITS:
    rows = records_by_split.get(split, [])
    if not rows:
        continue
    df = pd.DataFrame(rows)
    out = FWD_DIR / f"{split}.csv"
    df[_ordered_cols(df)].to_csv(out, index=False)
    print(f"  {split}.csv : {len(df):,} rows  ({df['n_snps_in_window'].sum():,} SNP-slots, "
          f"{df['window'].nunique():,} windows)")


# ── Write forward + reverse complement CSVs ───────────────────────────────────
print(f"\nWriting forward+reverse CSVs → {RC_DIR}")

for split in SPLITS:
    rows = records_by_split.get(split, [])
    if not rows:
        continue
    fwd = pd.DataFrame(rows)

    rc = fwd.copy()
    rc["sequence"] = rc["sequence"].apply(reverse_complement)
    rc["seq_id"]   = rc["seq_id"].str.replace("_f$", "_r", regex=True)

    combined = pd.concat([fwd, rc], ignore_index=True).sort_values("seq_id")
    out = RC_DIR / f"{split}.csv"
    combined[_ordered_cols(combined)].to_csv(out, index=False)
    print(f"  {split}.csv : {len(combined):,} rows  ({combined['seq_id'].nunique():,} unique)")

# ── Final summary ─────────────────────────────────────────────────────────────
print(f"\n{'='*65}")
print(f"  Multi-SNP context augmentation complete")
print(f"  Max window length : {MAX_LEN} bp")
print(f"  Total windows     : {total_windows}")
print(f"  Output:")
print(f"    {FWD_DIR}  (forward only)")
print(f"    {RC_DIR}  (forward + RC)")
print(f"{'='*65}")
