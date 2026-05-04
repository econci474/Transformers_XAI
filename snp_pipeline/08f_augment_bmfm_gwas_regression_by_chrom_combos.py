"""
08f_augment_bmfm_gwas_regression_by_chrom_combos.py
=====================================================
Combinatorial EA/OA augmentation of chromosome-level multi-SNP windows.

Extends 08e by generating ALL 2^N EA/OA allele combinations for each window
that contains N SNPs, rather than only the two extreme cases (all-EA / all-OA).

Dataset size (estimated from 08e windows): ~7,400 forward rows / ~14,800 fwd+RC.
This is small enough that no combinatorial cap is needed.

For a window with 3 SNPs the generated sequences are:

  mask (SNP1, SNP2, SNP3)   allele_tag          z_score
  (OA, OA, OA)               allOA_f             −(z1+z2+z3)
  (EA, OA, OA)               EA_OA_OA_f          +z1−z2−z3
  (OA, EA, OA)               OA_EA_OA_f          −z1+z2−z3
  (OA, OA, EA)               OA_OA_EA_f          −z1−z2+z3
  (EA, EA, OA)               EA_EA_OA_f          +z1+z2−z3
  (EA, OA, EA)               EA_OA_EA_f          +z1−z2+z3
  (OA, EA, EA)               OA_EA_EA_f          −z1+z2+z3
  (EA, EA, EA)               allEA_f             +(z1+z2+z3)

  Each sequence is also written as its reverse complement (_r suffix).

Train/dev/test split
--------------------
Split is at the WINDOW level (not SNP level) to avoid data leakage.
All sequences from the same window (all 2^N EA/OA combinations, forward
and reverse) are always assigned to the same split.
Windows are randomly shuffled (seed=42) then assigned 70/15/15.

Output
------
  bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom_combos/
    forward/
      train.csv, dev.csv, test.csv
    forward_and_reverse/
      train.csv, dev.csv, test.csv

Usage
-----
  python 08f_augment_bmfm_gwas_regression_by_chrom_combos.py
  python 08f_augment_bmfm_gwas_regression_by_chrom_combos.py --check
  python 08f_augment_bmfm_gwas_regression_by_chrom_combos.py --split 0.7 0.15 0.15 --seed 42
"""

import argparse
import io
import itertools
import pathlib
import random as _random
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
OUT_DIR    = BASE_DIR / "bmfm_inputs" / "bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom_combos"

SPLITS = ["train", "dev", "test"]

# ── CLI ───────────────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser(
    description="Combinatorial EA/OA augmentation for BMFM chromosome-level windows.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
_p.add_argument("--labels",         default=str(LABELS_TSV))
_p.add_argument("--fasta",          default=str(FASTA))
_p.add_argument("--outdir",         default=str(OUT_DIR))
_p.add_argument("--max-len",        type=int, default=8192, dest="max_len",
                help="Maximum sequence length in bp (default: 8192)")
_p.add_argument("--min-flank",      type=int, default=50,   dest="min_flank",
                help="Minimum bp flanking each end of the window (default: 50)")
_p.add_argument("--max-combo-snps", type=int, default=999,  dest="max_combo_snps",
                help="Cap on number of SNPs varied per window (default: unlimited). "
                     "For our dataset all windows have ≤10 SNPs so this has no effect.")
_p.add_argument("--split",          type=float, nargs=3, default=[0.80, 0.10, 0.10],
                metavar=("TRAIN", "DEV", "TEST"),
                help="Train/dev/test fractions at WINDOW level (default: 0.80 0.10 0.10)")
_p.add_argument("--seed",           type=int, default=42,
                help="Random seed for window shuffle (default: 42)")
_p.add_argument("--check",          action="store_true",
                help="Dry run — print statistics, do not write files")
args = _p.parse_args()

LABELS_TSV     = pathlib.Path(args.labels)
FASTA          = pathlib.Path(args.fasta)
OUT_DIR        = pathlib.Path(args.outdir)
MAX_LEN        = args.max_len
MIN_FLANK      = args.min_flank
MAX_COMBO_SNPS = args.max_combo_snps
TRAIN_FRAC, DEV_FRAC, TEST_FRAC = args.split
SEED           = args.seed
assert abs(TRAIN_FRAC + DEV_FRAC + TEST_FRAC - 1.0) < 1e-6, "Split fractions must sum to 1."

if 2 * MIN_FLANK >= MAX_LEN:
    print(f"[ERROR] min_flank ({MIN_FLANK}) * 2 >= max_len ({MAX_LEN})")
    sys.exit(1)

MAX_INNER_SPAN = MAX_LEN - 2 * MIN_FLANK

# ── Validate inputs ───────────────────────────────────────────────────────────
for p, name in [(LABELS_TSV, "Labels TSV"), (FASTA, "FASTA")]:
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


def _make_seq_all(ref_seq: str, snp_offsets: list, alleles: list) -> str:
    """Substitute all SNPs simultaneously with their respective alleles."""
    s = list(ref_seq)
    for offset, allele in zip(snp_offsets, alleles):
        if 0 <= offset < len(s):
            s[offset] = allele.upper()
    return "".join(s)


def _allele_tag(mask: tuple) -> str:
    """
    Encode the EA/OA mask as a string for seq_id.
      (1,1,1) → 'allEA'   (0,0,0) → 'allOA'   (1,0,1) → 'EA_OA_EA'
    """
    if all(m == 1 for m in mask):
        return "allEA"
    if all(m == 0 for m in mask):
        return "allOA"
    return "_".join("EA" if m else "OA" for m in mask)


# ── Load GW-sig SNPs ──────────────────────────────────────────────────────────
print(f"\nLoading labels: {LABELS_TSV.name}")
labels = pd.read_csv(LABELS_TSV, sep="\t", low_memory=False)
sig    = labels[labels["label"] == 1].copy()

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
    print("  Building .fai index (first use) ...")
    fa = Fasta(str(FASTA), rebuild=True)

chrom_names = set(fa.keys())

def resolve_chrom(chrom_val: str) -> str:
    for c in [str(chrom_val), f"chr{chrom_val}"]:
        if c in chrom_names:
            return c
    return None


# ── Window packing ────────────────────────────────────────────────────────────
def pack_windows(snp_list, max_inner_span: int) -> list:
    windows = []
    current = [snp_list[0]]
    for snp in snp_list[1:]:
        span = snp.BP - current[0].BP
        if span >= max_inner_span:
            windows.append(current)
            current = [snp]
        else:
            current.append(snp)
    windows.append(current)
    return windows


# ── Build records ─────────────────────────────────────────────────────────────
print(f"\nBuilding combinatorial EA/OA sequences (all 2^N per window)")
print(f"  Max sequence length : {MAX_LEN} bp")
print(f"  Min flank           : {MIN_FLANK} bp")
if MAX_COMBO_SNPS < 999:
    print(f"  Combo prefix cap    : {MAX_COMBO_SNPS} SNPs (SNPs beyond this fixed to OA)")
else:
    print(f"  Combo prefix cap    : none (full 2^N for every window)")
print(f"  Split strategy      : window-level random ({TRAIN_FRAC:.0%}/{DEV_FRAC:.0%}/{TEST_FRAC:.0%}, seed={SEED})")

# all_records: flat list — window-level split assigned after the loop
all_records = []
win_meta    = []   # (win_label, n_combos) — used in check mode
win_stats   = defaultdict(int)
skipped     = 0

for chrom, grp in sig.groupby("CHR"):
    chrom_fa = resolve_chrom(chrom)
    if chrom_fa is None:
        print(f"  [SKIP] CHR {chrom} not in FASTA")
        skipped += grp.shape[0]
        continue

    snp_rows  = list(grp.itertuples(index=False))
    windows   = pack_windows(snp_rows, MAX_INNER_SPAN)
    win_stats[chrom] = len(windows)
    chrom_len = len(fa[chrom_fa])

    for win_idx, window in enumerate(windows):
        n_snps   = len(window)
        first_bp = window[0].BP
        last_bp  = window[-1].BP

        inner_span  = last_bp - first_bp
        available   = MAX_LEN - inner_span - 2 * MIN_FLANK
        left_flank  = MIN_FLANK + available // 2
        right_flank = MIN_FLANK + (available - available // 2)

        win_start0 = max(0, first_bp - 1 - left_flank)
        win_end0   = min(chrom_len, last_bp + right_flank)

        win_label = f"chr{chrom}_w{win_idx:03d}"
        combo_n   = min(n_snps, MAX_COMBO_SNPS)
        fixed_n   = n_snps - combo_n

        if args.check:
            n_combos = 2 ** combo_n
            win_meta.append((win_label, n_combos))
            continue

        ref_seq     = str(fa[chrom_fa][win_start0:win_end0]).upper()
        actual_len  = len(ref_seq)
        snp_offsets = [(row.BP - 1 - win_start0) for row in window]

        eas = [str(row.EA).upper() for row in window]
        oas = [str(row.OA).upper() for row in window]
        zs  = [float(row.z_score) for row in window]

        for combo_mask in itertools.product([0, 1], repeat=combo_n):
            mask    = combo_mask + tuple([0] * fixed_n)  # variable prefix + OA suffix
            alleles = [eas[k] if mask[k] else oas[k] for k in range(n_snps)]
            z_score = sum(zs[k] if mask[k] else -zs[k] for k in range(n_snps))
            tag     = _allele_tag(mask)
            seq     = _make_seq_all(ref_seq, snp_offsets, alleles)

            record = {
                "sequence"         : seq,
                "z_score"          : z_score,
                "seq_id"           : f"{win_label}_{tag}_f",
                "window"           : win_label,
                "chrom"            : chrom,
                "n_snps_in_window" : n_snps,
                "n_ea_in_window"   : sum(mask),
                "combo_tag"        : tag,
                "combo_n"          : combo_n,
                "win_len_bp"       : actual_len,
            }

            for k, row in enumerate(window):
                suf = f"_snp{k + 1}"
                fixed_note = " (fixed)" if k >= combo_n else ""
                record[f"rsid{suf}"]          = row.SNP
                record[f"position_hg38{suf}"] = row.BP
                record[f"allele_type{suf}"]   = ("EA" if mask[k] else "OA") + fixed_note
                record[f"z_score{suf}"]       = zs[k] if mask[k] else -zs[k]

            all_records.append(record)

# ── Report windows ────────────────────────────────────────────────────────────
total_windows = sum(win_stats.values())
print(f"\n  Windows (≤{MAX_LEN} bp each):")
for chrom in sorted(win_stats, key=lambda x: (len(x), x)):
    print(f"    CHR {chrom:>3}: {win_stats[chrom]:>3} windows")
print(f"    Total    : {total_windows:>3} windows")

# ── Dry-run report ────────────────────────────────────────────────────────────
if args.check:
    win_labels = [w for w, _ in win_meta]
    _random.seed(SEED)
    _random.shuffle(win_labels)
    n_tr = int(len(win_labels) * TRAIN_FRAC)
    n_dv = int(len(win_labels) * DEV_FRAC)
    _wsplit = {}
    for i, w in enumerate(win_labels):
        _wsplit[w] = "train" if i < n_tr else ("dev" if i < n_tr + n_dv else "test")
    _by_split = defaultdict(lambda: [0, 0])  # [n_windows, n_rows]
    for wlabel, n_combos in win_meta:
        sp = _wsplit[wlabel]
        _by_split[sp][0] += 1
        _by_split[sp][1] += n_combos
    print(f"\n  Window-level split (seed={SEED}, {TRAIN_FRAC:.0%}/{DEV_FRAC:.0%}/{TEST_FRAC:.0%}):")
    for sp in SPLITS:
        nw, nr = _by_split[sp]
        print(f"    {sp}: {nw:>4} windows  {nr:>6} rows (fwd)  {nr*2:>6} rows (fwd+RC)")
    print("\n(Dry run — no files written)")
    sys.exit(0)

# ── Window-level random split ─────────────────────────────────────────────────
# All 2^N × 2 (fwd+RC) sequences from the same window go to the same split.
print(f"\nAssigning window-level split "
      f"({TRAIN_FRAC:.0%}/{DEV_FRAC:.0%}/{TEST_FRAC:.0%}, seed={SEED}) ...")

all_win_labels = list(dict.fromkeys(r["window"] for r in all_records))  # unique, ordered
_random.seed(SEED)
_random.shuffle(all_win_labels)

n_tr = int(len(all_win_labels) * TRAIN_FRAC)
n_dv = int(len(all_win_labels) * DEV_FRAC)
win_to_split = {}
for i, w in enumerate(all_win_labels):
    win_to_split[w] = "train" if i < n_tr else ("dev" if i < n_tr + n_dv else "test")

records_by_split = defaultdict(list)
for rec in all_records:
    records_by_split[win_to_split[rec["window"]]].append(rec)

for sp in SPLITS:
    recs = records_by_split.get(sp, [])
    wins = {r["window"] for r in recs}
    print(f"  {sp}: {len(wins):>4} windows  {len(recs):>6} rows (forward)")

# ── Column ordering ───────────────────────────────────────────────────────────
_PREFIX_COLS = ["sequence", "z_score", "seq_id", "window", "chrom",
                "n_snps_in_window", "n_ea_in_window", "combo_tag", "combo_n", "win_len_bp"]

def _ordered_cols(df: pd.DataFrame) -> list:
    snp_cols = sorted(
        [c for c in df.columns if c not in _PREFIX_COLS],
        key=lambda c: (c.rsplit("_snp", 1)[-1].zfill(4), c.rsplit("_snp", 1)[0])
    )
    return _PREFIX_COLS + snp_cols


# ── Write forward-strand CSVs ─────────────────────────────────────────────────
print(f"\nWriting forward-strand CSVs → {FWD_DIR}")

for split in SPLITS:
    rows = records_by_split.get(split, [])
    if not rows:
        continue
    df  = pd.DataFrame(rows)
    out = FWD_DIR / f"{split}.csv"
    df[_ordered_cols(df)].to_csv(out, index=False)
    print(f"  {split}.csv : {len(df):,} rows  "
          f"(n_ea range: {df['n_ea_in_window'].min()}–{df['n_ea_in_window'].max()}, "
          f"{df['window'].nunique()} windows)")


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


# ── Summary ───────────────────────────────────────────────────────────────────
total_fwd = sum(len(records_by_split.get(sp, [])) for sp in SPLITS)
print(f"\n{'='*65}")
print(f"  Combinatorial EA/OA augmentation complete")
print(f"  Max window length       : {MAX_LEN} bp")
print(f"  Total windows           : {total_windows}")
print(f"  Total forward rows      : {total_fwd:,}")
print(f"  Total fwd+RC rows       : {total_fwd*2:,}")
print(f"  Split strategy          : window-level random (seed={SEED})")
print(f"  Output:")
print(f"    {FWD_DIR}")
print(f"    {RC_DIR}")
print(f"{'='*65}")
