"""
08g_augment_bmfm_by_chrom_combos_singletons.py
================================================
Combinatorial EA/OA augmentation **plus singleton per-SNP crops**.

This script extends 08f by adding singleton rows for each individual SNP
in a window.  The combinatorial rows are identical to 08f — all 2^N EA/OA
allele combinations for each multi-SNP window.

Singleton rows
--------------
For each SNP k in a window with N SNPs, two singleton rows are produced
(EA and OA).  The sequence is a **crop** of the multi-SNP window centred
on SNP k, trimmed so that no other SNP position from the window appears
in the crop:

  left  boundary : position of SNP k−1 + 1 bp  (or window start if k=0)
  right boundary : position of SNP k+1 − 1 bp  (or window end   if k=N−1)

Within those bounds the crop extends up to MAX_LEN bp, centred on SNP k.
The focal SNP is substituted with EA (z = +z_k) or OA (z = −z_k).
All other positions are GRCh38 reference (no other salient positions are
visible to the model).

This way the model learns:
  • Marginal effect of each individual SNP  (singletons)
  • Combinatorial interaction effects       (combos)

Dataset size (estimated): ~7,400 combo forward rows + ~700 singleton
forward rows ≈ 8,100 forward / ~16,200 fwd+RC.

Train/dev/test split
--------------------
Window-level random split (seed=42, default 80/10/10).
All records from the same window — combinatorial AND singleton, forward
AND reverse — are always assigned to the same split.

Output
------
  bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom_combos_singletons/
    forward/
      train.csv, dev.csv, test.csv
    forward_and_reverse/
      train.csv, dev.csv, test.csv

Usage
-----
  python 08g_augment_bmfm_by_chrom_combos_singletons.py
  python 08g_augment_bmfm_by_chrom_combos_singletons.py --check
  python 08g_augment_bmfm_by_chrom_combos_singletons.py --split 0.7 0.15 0.15 --seed 42
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
OUT_DIR    = BASE_DIR / "bmfm_inputs" / "bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom_combos_singletons"

SPLITS = ["train", "dev", "test"]

# ── CLI ───────────────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser(
    description="Combinatorial EA/OA augmentation + singleton per-SNP crops.",
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
                help="Cap on number of SNPs varied per window (default: unlimited).")
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


def _make_seq_single(ref_seq: str, offset: int, allele: str) -> str:
    """Substitute a single SNP with the given allele."""
    s = list(ref_seq)
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
print(f"\nBuilding combinatorial EA/OA sequences + singleton crops")
print(f"  Max sequence length : {MAX_LEN} bp")
print(f"  Min flank           : {MIN_FLANK} bp")
if MAX_COMBO_SNPS < 999:
    print(f"  Combo prefix cap    : {MAX_COMBO_SNPS} SNPs (SNPs beyond this fixed to OA)")
else:
    print(f"  Combo prefix cap    : none (full 2^N for every window)")
print(f"  Split strategy      : window-level random ({TRAIN_FRAC:.0%}/{DEV_FRAC:.0%}/{TEST_FRAC:.0%}, seed={SEED})")

all_records    = []   # flat list — window-level split assigned after the loop
win_meta       = []   # (win_label, n_combos, n_singletons) — used in check mode
win_stats      = defaultdict(int)
singleton_stats = {"total": 0, "min_len": MAX_LEN, "max_len": 0}
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

        # 0-based genomic positions of each SNP
        snp_positions_0 = [row.BP - 1 for row in window]

        if args.check:
            n_combos = 2 ** combo_n
            n_singletons = 2 * n_snps  # EA + OA per SNP
            win_meta.append((win_label, n_combos, n_singletons))
            singleton_stats["total"] += n_singletons
            # Estimate singleton crop lengths
            for k in range(n_snps):
                left_bound  = (snp_positions_0[k - 1] + 1) if k > 0 else max(0, snp_positions_0[k] - MAX_LEN // 2)
                right_bound = snp_positions_0[k + 1] if k < n_snps - 1 else min(chrom_len, snp_positions_0[k] + MAX_LEN // 2 + 1)
                crop_len = min(MAX_LEN, right_bound - left_bound)
                singleton_stats["min_len"] = min(singleton_stats["min_len"], crop_len)
                singleton_stats["max_len"] = max(singleton_stats["max_len"], crop_len)
            continue

        ref_seq     = str(fa[chrom_fa][win_start0:win_end0]).upper()
        actual_len  = len(ref_seq)
        snp_offsets = [(pos0 - win_start0) for pos0 in snp_positions_0]

        eas = [str(row.EA).upper() for row in window]
        oas = [str(row.OA).upper() for row in window]
        zs  = [float(row.z_score) for row in window]

        # ── Combinatorial rows (identical to 08f) ─────────────────────────────
        for combo_mask in itertools.product([0, 1], repeat=combo_n):
            mask    = combo_mask + tuple([0] * fixed_n)
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
                "record_type"      : "combo",
            }

            for k, row in enumerate(window):
                suf = f"_snp{k + 1}"
                fixed_note = " (fixed)" if k >= combo_n else ""
                record[f"rsid{suf}"]          = row.SNP
                record[f"position_hg38{suf}"] = row.BP
                record[f"allele_type{suf}"]   = ("EA" if mask[k] else "OA") + fixed_note
                record[f"z_score{suf}"]       = zs[k] if mask[k] else -zs[k]

            all_records.append(record)

        # ── Singleton rows ────────────────────────────────────────────────────
        for k, row in enumerate(window):
            # Crop boundaries: exclude positions of other SNPs in this window
            # left_bound / right_bound are 0-based genomic coordinates
            # defining the region that contains NO other salient SNP positions.
            left_bound  = (snp_positions_0[k - 1] + 1) if k > 0         else 0
            right_bound = snp_positions_0[k + 1]       if k < n_snps - 1 else chrom_len
            # right_bound is exclusive (python slice convention)

            # Centre on focal SNP, cap at MAX_LEN
            snp_0 = snp_positions_0[k]
            avail = right_bound - left_bound
            crop_len = min(MAX_LEN, avail)

            half = crop_len // 2
            crop_start = snp_0 - half
            crop_end   = crop_start + crop_len

            # Shift if the centred crop exceeds the allowed bounds
            if crop_start < left_bound:
                crop_start = left_bound
                crop_end   = crop_start + crop_len
            if crop_end > right_bound:
                crop_end   = right_bound
                crop_start = crop_end - crop_len

            # Final clamp to chromosome boundaries
            crop_start = max(0, crop_start)
            crop_end   = min(chrom_len, crop_end)

            crop_seq = str(fa[chrom_fa][crop_start:crop_end]).upper()
            focal_offset = snp_0 - crop_start
            crop_actual_len = len(crop_seq)

            singleton_stats["total"] += 2
            singleton_stats["min_len"] = min(singleton_stats["min_len"], crop_actual_len)
            singleton_stats["max_len"] = max(singleton_stats["max_len"], crop_actual_len)

            for allele, allele_tag, z_sign in [(eas[k], "singleEA", +1),
                                                (oas[k], "singleOA", -1)]:
                seq = _make_seq_single(crop_seq, focal_offset, allele)
                z_score = z_sign * zs[k]

                record = {
                    "sequence"         : seq,
                    "z_score"          : z_score,
                    "seq_id"           : f"{win_label}_{row.SNP}_{allele_tag}_f",
                    "window"           : win_label,
                    "chrom"            : chrom,
                    "n_snps_in_window" : n_snps,
                    "n_ea_in_window"   : 1 if allele_tag == "singleEA" else 0,
                    "combo_tag"        : allele_tag,
                    "combo_n"          : 1,
                    "win_len_bp"       : crop_actual_len,
                    "record_type"      : "singleton",
                }

                # Only the focal SNP gets metadata columns
                suf = f"_snp{k + 1}"
                record[f"rsid{suf}"]          = row.SNP
                record[f"position_hg38{suf}"] = row.BP
                record[f"allele_type{suf}"]   = "EA" if allele_tag == "singleEA" else "OA"
                record[f"z_score{suf}"]       = z_score

                all_records.append(record)

# ── Report windows ────────────────────────────────────────────────────────────
total_windows = sum(win_stats.values())
print(f"\n  Windows (≤{MAX_LEN} bp each):")
for chrom in sorted(win_stats, key=lambda x: (len(x), x)):
    print(f"    CHR {chrom:>3}: {win_stats[chrom]:>3} windows")
print(f"    Total    : {total_windows:>3} windows")

print(f"\n  Singleton crop lengths: "
      f"min={singleton_stats['min_len']:,} bp  "
      f"max={singleton_stats['max_len']:,} bp  "
      f"total singleton rows (fwd)={singleton_stats['total']:,}")

# ── Dry-run report ────────────────────────────────────────────────────────────
if args.check:
    win_labels = [w for w, _, _ in win_meta]
    _random.seed(SEED)
    _random.shuffle(win_labels)
    n_tr = int(len(win_labels) * TRAIN_FRAC)
    n_dv = int(len(win_labels) * DEV_FRAC)
    _wsplit = {}
    for i, w in enumerate(win_labels):
        _wsplit[w] = "train" if i < n_tr else ("dev" if i < n_tr + n_dv else "test")
    _by_split = defaultdict(lambda: [0, 0, 0])  # [n_windows, n_combo_rows, n_singleton_rows]
    for wlabel, n_combos, n_singletons in win_meta:
        sp = _wsplit[wlabel]
        _by_split[sp][0] += 1
        _by_split[sp][1] += n_combos
        _by_split[sp][2] += n_singletons
    print(f"\n  Window-level split (seed={SEED}, {TRAIN_FRAC:.0%}/{DEV_FRAC:.0%}/{TEST_FRAC:.0%}):")
    for sp in SPLITS:
        nw, nc, ns = _by_split[sp]
        total_fwd = nc + ns
        print(f"    {sp}: {nw:>4} windows  "
              f"{nc:>6} combo + {ns:>5} singleton = {total_fwd:>6} rows (fwd)  "
              f"{total_fwd*2:>6} rows (fwd+RC)")
    print("\n(Dry run — no files written)")
    sys.exit(0)

# ── Window-level random split ─────────────────────────────────────────────────
print(f"\nAssigning window-level split "
      f"({TRAIN_FRAC:.0%}/{DEV_FRAC:.0%}/{TEST_FRAC:.0%}, seed={SEED}) ...")

all_win_labels = list(dict.fromkeys(r["window"] for r in all_records))
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
    n_combo = sum(1 for r in recs if r["record_type"] == "combo")
    n_single = sum(1 for r in recs if r["record_type"] == "singleton")
    print(f"  {sp}: {len(wins):>4} windows  "
          f"{n_combo:>6} combo + {n_single:>5} singleton = {len(recs):>6} rows (forward)")

# ── Column ordering ───────────────────────────────────────────────────────────
_PREFIX_COLS = ["sequence", "z_score", "seq_id", "window", "chrom",
                "n_snps_in_window", "n_ea_in_window", "combo_tag",
                "combo_n", "win_len_bp", "record_type"]

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
    n_combo  = (df["record_type"] == "combo").sum()
    n_single = (df["record_type"] == "singleton").sum()
    print(f"  {split}.csv : {len(df):,} rows  "
          f"({n_combo:,} combo + {n_single:,} singleton, "
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
total_fwd    = sum(len(records_by_split.get(sp, [])) for sp in SPLITS)
total_combo  = sum(1 for r in all_records if r["record_type"] == "combo")
total_single = sum(1 for r in all_records if r["record_type"] == "singleton")
print(f"\n{'='*70}")
print(f"  Combinatorial + Singleton augmentation complete")
print(f"  Max window length           : {MAX_LEN} bp")
print(f"  Total windows               : {total_windows}")
print(f"  Total forward rows          : {total_fwd:,}")
print(f"    ├─ combo rows             : {total_combo:,}")
print(f"    └─ singleton rows         : {total_single:,}")
print(f"  Total fwd+RC rows           : {total_fwd*2:,}")
print(f"  Singleton crop lengths      : {singleton_stats['min_len']:,}–{singleton_stats['max_len']:,} bp")
print(f"  Split strategy              : window-level random (seed={SEED})")
print(f"  Output:")
print(f"    {FWD_DIR}")
print(f"    {RC_DIR}")
print(f"{'='*70}")
