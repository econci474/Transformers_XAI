"""
08h_prepare_bmfm_gwas_classification_inputs.py
================================================
Prepare BMFM-DNA-SNP classification dataset with binned z-score labels.

For each label=1 (GW-significant) SNP:
  - Singleton crop (same as 08g) excluding neighbouring GW-sig positions
  - EA row → class = bin(+z);  OA row → class = bin(−z)

For each label=0 (null/non-associated) SNP:
  - Centred window (up to MAX_LEN bp) excluding nearby label=1 positions
  - REF allele → class = 0

Binning scheme
--------------
  σ = std(|z_score|) across label=1 SNPs
  z_bin ∈ {-3, -2, -1, 0, +1, +2, +3}

  For label=1 rows with z_effective (= +z for EA, −z for OA):
    |z_eff| <  σ   → class sign(z_eff) * 1
     σ ≤ |z_eff| < 2σ → class sign(z_eff) * 2
    |z_eff| ≥ 2σ  → class sign(z_eff) * 3

  For label=0 rows:  class = 0

Split
-----
  SNP-level random split (seed=42, default 80/10/10).
  Label=1 singletons use window-level splits (same as 08g) so that SNPs
  within the same window remain in the same split.
  Label=0 singletons are split independently at the SNP level.

Output
------
  bmfm_gwas_classification_without_ukb/
    forward/          train.csv, dev.csv, test.csv
    forward_and_reverse/  train.csv, dev.csv, test.csv
    label_dict.json

Usage
-----
  python 08h_prepare_bmfm_gwas_classification_inputs.py
  python 08h_prepare_bmfm_gwas_classification_inputs.py --check
  python 08h_prepare_bmfm_gwas_classification_inputs.py --max-null 2000
"""

import argparse
import io
import json
import math
import pathlib
import random as _random
import sys
from bisect import bisect_left, bisect_right
from collections import defaultdict

import numpy as np
import pandas as pd

# Force UTF-8 stdout on Windows
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR   = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220")
LABELS_TSV = BASE_DIR / "bmfm_inputs" / "external_gwas_labels.tsv"
FASTA      = BASE_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
OUT_DIR    = BASE_DIR / "bmfm_inputs" / "bmfm_gwas_classification_without_ukb"

SPLITS = ["train", "dev", "test"]

# ── CLI ───────────────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser(
    description="Prepare BMFM classification dataset with binned z-score labels.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
_p.add_argument("--labels",    default=str(LABELS_TSV))
_p.add_argument("--fasta",     default=str(FASTA))
_p.add_argument("--outdir",    default=str(OUT_DIR))
_p.add_argument("--max-len",   type=int, default=8192, dest="max_len",
                help="Maximum sequence length in bp (default: 8192)")
_p.add_argument("--min-flank", type=int, default=50,   dest="min_flank",
                help="Minimum bp flanking each end (default: 50)")
_p.add_argument("--max-null",  type=int, default=None, dest="max_null",
                help="Cap number of label=0 SNPs (default: all)")
_p.add_argument("--split",     type=float, nargs=3, default=[0.80, 0.10, 0.10],
                metavar=("TRAIN", "DEV", "TEST"),
                help="Train/dev/test fractions (default: 0.80 0.10 0.10)")
_p.add_argument("--seed",      type=int, default=42,
                help="Random seed (default: 42)")
_p.add_argument("--check",     action="store_true",
                help="Dry run — print statistics, do not write files")
args = _p.parse_args()

LABELS_TSV     = pathlib.Path(args.labels)
FASTA          = pathlib.Path(args.fasta)
OUT_DIR        = pathlib.Path(args.outdir)
MAX_LEN        = args.max_len
MIN_FLANK      = args.min_flank
TRAIN_FRAC, DEV_FRAC, TEST_FRAC = args.split
SEED           = args.seed
assert abs(TRAIN_FRAC + DEV_FRAC + TEST_FRAC - 1.0) < 1e-6, "Split fractions must sum to 1."

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


def _make_seq_single(ref_seq: str, offset: int, allele: str) -> str:
    """Substitute a single SNP with the given allele."""
    s = list(ref_seq)
    if 0 <= offset < len(s):
        s[offset] = allele.upper()
    return "".join(s)


# ── Window packing (for label=1 SNPs — identical to 08g) ─────────────────────
def pack_windows(snp_list, max_inner_span: int) -> list:
    """Greedy packing: group sorted SNPs into windows of bounded span."""
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


# ── Load labels ───────────────────────────────────────────────────────────────
print(f"\nLoading labels: {LABELS_TSV.name}")
labels = pd.read_csv(LABELS_TSV, sep="\t", low_memory=False)
sig    = labels[labels["label"] == 1].copy()
nulls  = labels[labels["label"] == 0].copy()

# Fill missing OA for label=1
oa_missing = sig["OA"].isna()
if oa_missing.any():
    sig.loc[oa_missing, "OA"] = sig.loc[oa_missing, "REF"]

sig = sig.dropna(subset=["EA", "OA", "z_score"]).copy()
sig["CHR"] = sig["CHR"].astype(str)
sig = sig.sort_values(["CHR", "BP"]).reset_index(drop=True)

nulls["CHR"] = nulls["CHR"].astype(str)
nulls = nulls.sort_values(["CHR", "BP"]).reset_index(drop=True)

if args.max_null and len(nulls) > args.max_null:
    nulls = nulls.sample(n=args.max_null, random_state=SEED).sort_values(
        ["CHR", "BP"]).reset_index(drop=True)
    print(f"  Subsampled label=0 to {len(nulls):,}")

print(f"  label=1 SNPs (GW-sig)  : {len(sig):,}")
print(f"  label=0 SNPs (null)    : {len(nulls):,}")

# ── Compute z-score bins ──────────────────────────────────────────────────────
z_abs = sig["z_score"].abs()
z_sigma = z_abs.std()
z_mu    = z_abs.mean()
print(f"\n  z-score statistics (label=1):")
print(f"    μ(|z|) = {z_mu:.2f}")
print(f"    σ(|z|) = {z_sigma:.2f}")
print(f"    bin edges: ±1σ = ±{z_sigma:.2f},  ±2σ = ±{2*z_sigma:.2f}")

def assign_z_bin(z_effective: float) -> int:
    """Map z_effective → integer class in {-3, -2, -1, +1, +2, +3}."""
    az = abs(z_effective)
    if az < z_sigma:
        level = 1
    elif az < 2 * z_sigma:
        level = 2
    else:
        level = 3
    return level if z_effective >= 0 else -level

# Preview bin distribution for label=1 SNPs
# (each SNP produces EA (+z) and OA (-z) rows)
bin_counts = defaultdict(int)
for z in sig["z_score"]:
    bin_counts[assign_z_bin(+z)] += 1   # EA
    bin_counts[assign_z_bin(-z)] += 1   # OA

print(f"\n  Projected z_bin distribution (label=1, forward only):")
for b in sorted(bin_counts):
    print(f"    class {b:+d}: {bin_counts[b]:>5}")
print(f"    class  0: {len(nulls):>5}  (label=0 nulls)")

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

# ── Build sorted label=1 position lookup (for exclusion in label=0 crops) ────
# sig_pos_by_chrom: { chrom_str → sorted list of 0-based positions }
sig_pos_by_chrom = defaultdict(list)
for row in sig.itertuples(index=False):
    sig_pos_by_chrom[str(row.CHR)].append(row.BP - 1)  # 0-based
for chrom in sig_pos_by_chrom:
    sig_pos_by_chrom[chrom].sort()


# ── Process label=1 SNPs (singleton crops like 08g, no combos) ────────────────
print(f"\nBuilding label=1 singleton sequences ...")
print(f"  Max sequence length : {MAX_LEN} bp")
print(f"  Min flank           : {MIN_FLANK} bp")

all_records    = []
win_labels_l1  = []   # for window-level splitting of label=1
win_stats      = defaultdict(int)
singleton_stats = {"sig_total": 0, "null_total": 0, "min_len": MAX_LEN, "max_len": 0}
skipped        = 0

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
        win_label = f"chr{chrom}_w{win_idx:03d}"
        win_labels_l1.append(win_label)

        # 0-based genomic positions of each SNP in this window
        snp_positions_0 = [row.BP - 1 for row in window]

        for k, row in enumerate(window):
            # Crop boundaries: exclude positions of other SNPs in this window
            left_bound  = (snp_positions_0[k - 1] + 1) if k > 0         else 0
            right_bound = snp_positions_0[k + 1]       if k < n_snps - 1 else chrom_len

            # Centre on focal SNP, cap at MAX_LEN
            snp_0 = snp_positions_0[k]
            avail = right_bound - left_bound
            crop_len = min(MAX_LEN, avail)

            half = crop_len // 2
            crop_start = snp_0 - half
            crop_end   = crop_start + crop_len

            # Shift if crop exceeds allowed bounds
            if crop_start < left_bound:
                crop_start = left_bound
                crop_end   = crop_start + crop_len
            if crop_end > right_bound:
                crop_end   = right_bound
                crop_start = crop_end - crop_len

            # Clamp to chromosome
            crop_start = max(0, crop_start)
            crop_end   = min(chrom_len, crop_end)

            if args.check:
                clen = crop_end - crop_start
                singleton_stats["sig_total"] += 2
                singleton_stats["min_len"] = min(singleton_stats["min_len"], clen)
                singleton_stats["max_len"] = max(singleton_stats["max_len"], clen)
                continue

            crop_seq = str(fa[chrom_fa][crop_start:crop_end]).upper()
            focal_offset = snp_0 - crop_start
            crop_actual_len = len(crop_seq)

            singleton_stats["sig_total"] += 2
            singleton_stats["min_len"] = min(singleton_stats["min_len"], crop_actual_len)
            singleton_stats["max_len"] = max(singleton_stats["max_len"], crop_actual_len)

            ea = str(row.EA).upper()
            oa = str(row.OA).upper()
            z  = float(row.z_score)

            for allele, allele_tag, z_eff in [(ea, "EA", +z), (oa, "OA", -z)]:
                seq = _make_seq_single(crop_seq, focal_offset, allele)
                z_bin = assign_z_bin(z_eff)

                record = {
                    "sequence"    : seq,
                    "z_bin"       : z_bin,
                    "z_score"     : z_eff,
                    "seq_id"      : f"{win_label}_{row.SNP}_{allele_tag}_f",
                    "window"      : win_label,
                    "chrom"       : chrom,
                    "rsid"        : row.SNP,
                    "position"    : row.BP,
                    "allele"      : allele,
                    "allele_type" : allele_tag,
                    "record_type" : "sig_singleton",
                    "win_len_bp"  : crop_actual_len,
                }
                all_records.append(record)

# Report label=1
total_l1_windows = sum(win_stats.values())
print(f"\n  Label=1 windows (≤{MAX_LEN} bp each):")
for chrom in sorted(win_stats, key=lambda x: (len(x), x)):
    print(f"    CHR {chrom:>3}: {win_stats[chrom]:>3} windows")
print(f"    Total    : {total_l1_windows:>3} windows")
print(f"    Rows (fwd): {singleton_stats['sig_total']:,}")

# ── Process label=0 SNPs (null singletons) ────────────────────────────────────
print(f"\nBuilding label=0 singleton sequences ...")

null_snp_ids = []  # for SNP-level splitting

for chrom, grp in nulls.groupby("CHR"):
    chrom_fa = resolve_chrom(chrom)
    if chrom_fa is None:
        skipped += grp.shape[0]
        continue

    chrom_len = len(fa[chrom_fa])

    # Sorted label=1 positions on this chromosome (for exclusion)
    sig_positions = sig_pos_by_chrom.get(chrom, [])

    for row in grp.itertuples(index=False):
        snp_0 = row.BP - 1   # 0-based
        null_id = f"null_{row.SNP}"
        null_snp_ids.append(null_id)

        # Find crop bounds excluding nearby GW-sig SNP positions
        # Left bound: nearest GW-sig position to the left + 1
        idx_left = bisect_left(sig_positions, snp_0) - 1
        left_bound = (sig_positions[idx_left] + 1) if idx_left >= 0 else 0

        # Right bound: nearest GW-sig position to the right (exclusive)
        idx_right = bisect_right(sig_positions, snp_0)
        right_bound = sig_positions[idx_right] if idx_right < len(sig_positions) else chrom_len

        # Centre on focal SNP, cap at MAX_LEN
        avail = right_bound - left_bound
        crop_len = min(MAX_LEN, avail)

        if crop_len < 10:
            # Too close to a GW-sig SNP on both sides; skip
            skipped += 1
            continue

        half = crop_len // 2
        crop_start = snp_0 - half
        crop_end   = crop_start + crop_len

        if crop_start < left_bound:
            crop_start = left_bound
            crop_end   = crop_start + crop_len
        if crop_end > right_bound:
            crop_end   = right_bound
            crop_start = crop_end - crop_len

        crop_start = max(0, crop_start)
        crop_end   = min(chrom_len, crop_end)

        if args.check:
            clen = crop_end - crop_start
            singleton_stats["null_total"] += 1
            singleton_stats["min_len"] = min(singleton_stats["min_len"], clen)
            singleton_stats["max_len"] = max(singleton_stats["max_len"], clen)
            continue

        crop_seq = str(fa[chrom_fa][crop_start:crop_end]).upper()
        focal_offset = snp_0 - crop_start
        crop_actual_len = len(crop_seq)

        singleton_stats["null_total"] += 1
        singleton_stats["min_len"] = min(singleton_stats["min_len"], crop_actual_len)
        singleton_stats["max_len"] = max(singleton_stats["max_len"], crop_actual_len)

        # REF allele at focal position
        ref_allele = str(row.REF).upper() if pd.notna(row.REF) else str(row.A2).upper()
        seq = _make_seq_single(crop_seq, focal_offset, ref_allele)

        record = {
            "sequence"    : seq,
            "z_bin"       : 0,
            "z_score"     : 0.0,
            "seq_id"      : f"{null_id}_REF_f",
            "window"      : null_id,
            "chrom"       : chrom,
            "rsid"        : row.SNP,
            "position"    : row.BP,
            "allele"      : ref_allele,
            "allele_type" : "REF",
            "record_type" : "null_singleton",
            "win_len_bp"  : crop_actual_len,
        }
        all_records.append(record)

print(f"  Label=0 rows (fwd): {singleton_stats['null_total']:,}")
print(f"  Skipped total     : {skipped:,}")
print(f"\n  Singleton crop lengths: "
      f"min={singleton_stats['min_len']:,} bp  "
      f"max={singleton_stats['max_len']:,} bp")

# ── Dry-run report ────────────────────────────────────────────────────────────
if args.check:
    print(f"\n  Projected dataset:")
    n_sig_fwd  = singleton_stats["sig_total"]
    n_null_fwd = singleton_stats["null_total"]
    total_fwd  = n_sig_fwd + n_null_fwd
    print(f"    label=1 (sig)  : {n_sig_fwd:>7,} rows (forward)")
    print(f"    label=0 (null) : {n_null_fwd:>7,} rows (forward)")
    print(f"    total          : {total_fwd:>7,} rows (forward)")
    print(f"    total fwd+RC   : {total_fwd*2:>7,} rows")

    # Show bin distribution
    print(f"\n  z_bin class distribution (forward only):")
    for b in sorted(bin_counts):
        print(f"    class {b:+d}: {bin_counts[b]:>5}")
    print(f"    class  0: {n_null_fwd:>5}")

    print("\n(Dry run — no files written)")
    sys.exit(0)

# ── Split assignment ──────────────────────────────────────────────────────────
# Label=1: window-level split (all SNPs in same window → same split)
# Label=0: SNP-level split (independent)
print(f"\nAssigning splits "
      f"({TRAIN_FRAC:.0%}/{DEV_FRAC:.0%}/{TEST_FRAC:.0%}, seed={SEED}) ...")

# Label=1 window split
unique_l1_wins = list(dict.fromkeys(win_labels_l1))
_random.seed(SEED)
_random.shuffle(unique_l1_wins)
n_tr = int(len(unique_l1_wins) * TRAIN_FRAC)
n_dv = int(len(unique_l1_wins) * DEV_FRAC)
l1_win_to_split = {}
for i, w in enumerate(unique_l1_wins):
    l1_win_to_split[w] = "train" if i < n_tr else ("dev" if i < n_tr + n_dv else "test")

# Label=0 SNP split
unique_null_ids = list(dict.fromkeys(null_snp_ids))
_random.seed(SEED + 1)  # different seed to avoid correlation
_random.shuffle(unique_null_ids)
n_tr_null = int(len(unique_null_ids) * TRAIN_FRAC)
n_dv_null = int(len(unique_null_ids) * DEV_FRAC)
null_to_split = {}
for i, nid in enumerate(unique_null_ids):
    null_to_split[nid] = "train" if i < n_tr_null else ("dev" if i < n_tr_null + n_dv_null else "test")

# Assign splits
records_by_split = defaultdict(list)
for rec in all_records:
    if rec["record_type"] == "sig_singleton":
        sp = l1_win_to_split[rec["window"]]
    else:
        sp = null_to_split[rec["window"]]
    records_by_split[sp].append(rec)

for sp in SPLITS:
    recs = records_by_split.get(sp, [])
    n_sig  = sum(1 for r in recs if r["record_type"] == "sig_singleton")
    n_null = sum(1 for r in recs if r["record_type"] == "null_singleton")
    print(f"  {sp}: {n_sig:>5} sig + {n_null:>6} null = {len(recs):>6} rows (forward)")

# ── Column ordering ───────────────────────────────────────────────────────────
COL_ORDER = ["sequence", "z_bin", "z_score", "seq_id", "window", "chrom",
             "rsid", "position", "allele", "allele_type", "record_type",
             "win_len_bp"]

# ── Write forward-strand CSVs ─────────────────────────────────────────────────
print(f"\nWriting forward-strand CSVs → {FWD_DIR}")

for split in SPLITS:
    rows = records_by_split.get(split, [])
    if not rows:
        continue
    df  = pd.DataFrame(rows)
    out = FWD_DIR / f"{split}.csv"
    df[COL_ORDER].to_csv(out, index=False)
    n_sig  = (df["record_type"] == "sig_singleton").sum()
    n_null = (df["record_type"] == "null_singleton").sum()
    print(f"  {split}.csv : {len(df):,} rows  "
          f"({n_sig:,} sig + {n_null:,} null)")

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
    combined[COL_ORDER].to_csv(out, index=False)
    print(f"  {split}.csv : {len(combined):,} rows  ({combined['seq_id'].nunique():,} unique)")

# ── Write label_dict.json ─────────────────────────────────────────────────────
label_dict = {str(i): i for i in [-3, -2, -1, 0, 1, 2, 3]}
label_dict_path = OUT_DIR / "label_dict.json"
if not args.check:
    with open(label_dict_path, "w") as f:
        json.dump(label_dict, f, indent=2)
    print(f"\nLabel dict → {label_dict_path}")

# ── Summary ───────────────────────────────────────────────────────────────────
total_fwd    = sum(len(records_by_split.get(sp, [])) for sp in SPLITS)
total_sig    = sum(1 for r in all_records if r["record_type"] == "sig_singleton")
total_null   = sum(1 for r in all_records if r["record_type"] == "null_singleton")

# Class distribution
class_dist = defaultdict(int)
for r in all_records:
    class_dist[r["z_bin"]] += 1

print(f"\n{'='*70}")
print(f"  BMFM GWAS Classification dataset prepared")
print(f"  Binning σ (from GW-sig |z|)  : {z_sigma:.2f}")
print(f"  Bin edges                    : ±{z_sigma:.1f}, ±{2*z_sigma:.1f}")
print(f"  Total forward rows           : {total_fwd:,}")
print(f"    ├─ sig singletons (label=1): {total_sig:,}")
print(f"    └─ null singletons (label=0): {total_null:,}")
print(f"  Total fwd+RC rows            : {total_fwd*2:,}")
print(f"  Class distribution (forward):")
for b in sorted(class_dist):
    print(f"    class {b:+d}: {class_dist[b]:>6}")
print(f"  Crop lengths                 : {singleton_stats['min_len']:,}–{singleton_stats['max_len']:,} bp")
print(f"  Split strategy               : window(l=1)/SNP(l=0) random (seed={SEED})")
print(f"  Output:")
print(f"    {FWD_DIR}")
print(f"    {RC_DIR}")
print(f"    {label_dict_path}")
print(f"{'='*70}")
