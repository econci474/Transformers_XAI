"""
liftover_hg19_to_hg38_wightman_2021.py
========================================
Lift Wightman 2021 AD GWAS summary statistics from GRCh37 → GRCh38 and
annotate each variant with the GRCh38 REF allele from the Ensembl FASTA.

Background
----------
Both Wightman 2021 files store coordinates in GRCh37.  The ADNI BIM used
for BMFM label generation is GRCh38, so positional matching requires a
liftover step.

Important allele note
---------------------
In these summary-stat files the two allele columns are **both ALT alleles**
(the tested/effect allele and the other alternate allele).  Neither is
the GRCh38 reference.  This script therefore fetches REF at each lifted
position from the GRCh38 primary-assembly FASTA and adds a REF_38 column.

Input files  (selected via --with-ukb flag)
-------------------------------------------
  Default (excl UKB + 23&Me):
    PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis (1).txt
    Columns: chromosome  base_pair_location  effect_allele  other_allele
             beta  standard_error  effect_allele_frequency  p_value
             SNPID  N  Neffective  Build

  --with-ukb (incl UKB, excl 23&Me only):
    PGCALZ2sumstatsExcluding23andMe.txt
    Columns: chr  PosGRCh37  testedAllele  otherAllele  z  p  N

Output files
------------
  wightman_without_ukb_GRCh38.tsv   (default)
  wightman_with_ukb_GRCh38.tsv      (--with-ukb)

  Columns in output:
    CHR_37  BP_37  CHR_38  BP_38
    effect_allele  other_allele  REF_38
    [original stat columns: beta/z, p_value, N, ...]
    SNPID_38   (format: CHR:BP:REF:effect_allele)
    lifted     (True/False)

Usage
-----
  pip install pyliftover pyfaidx pandas
  python liftover_hg19_to_hg38_wightman_2021.py
  python liftover_hg19_to_hg38_wightman_2021.py --with-ukb
  python liftover_hg19_to_hg38_wightman_2021.py --all-snps   # lift ALL SNPs, not just GW-sig
  python liftover_hg19_to_hg38_wightman_2021.py --gw-thresh 5e-8
"""

import argparse
import io
import pathlib
import sys
import time
import urllib.request
from datetime import datetime

try:
    import requests as _requests
except ImportError:
    _requests = None

import pandas as pd

if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR  = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220")
GWAS_DIR  = BASE_DIR / "GWAS" / "Wightman_2021"
OUT_DIR   = BASE_DIR / "GWAS" / "Wightman_2021"
LIFT_DIR  = BASE_DIR / "liftover"

CHAIN_URL  = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
CHAIN_FILE = LIFT_DIR / "hg19ToHg38.over.chain.gz"

GRCH38_FA = BASE_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# ADNI hg19 BIM — used to look up rsIDs from Wightman hg19 positions
# (avoids liftover coord-shift problems when matching the ADNI array)
ADNI_BIM_HG19   = BASE_DIR / "SNP_filtered_with_mri_rsid_clean_current.bim"
NCBI_MERGE_ARCH = BASE_DIR / "liftover" / "NCBI" / "RsMergeArch.bcp.gz"

# File definitions
EXCL_UKB_DIR  = GWAS_DIR / "PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis (1).txt"
EXCL_UKB_FILE = EXCL_UKB_DIR / "PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis (1).txt"

INCL_UKB_DIR  = GWAS_DIR / "PGCALZ2sumstatsExcluding23andMe.txt"
INCL_UKB_FILE = INCL_UKB_DIR / "PGCALZ2sumstatsExcluding23andMe.txt"

CANONICAL = set([str(i) for i in range(1, 23)] + ["X", "Y", "MT"])

# ── File format definitions ────────────────────────────────────────────────────
# Each format maps source column names → canonical names used internally.
FORMATS = {
    "excl_ukb": {
        "file":         EXCL_UKB_FILE,
        "sep":          r"\s+",
        "chr_col":      "chromosome",
        "pos_col":      "base_pair_location",
        "ea_col":       "effect_allele",
        "oa_col":       "other_allele",
        "p_col":        "p_value",
        "stat_cols":    ["beta", "standard_error", "effect_allele_frequency",
                         "p_value", "SNPID", "N", "Neffective"],
        "out_stem":     "wightman_without_ukb_GRCh38",
        "label":        "Wightman 2021 (excl UKB + 23&Me)",
    },
    "incl_ukb": {
        "file":         INCL_UKB_FILE,
        "sep":          "\t",
        "chr_col":      "chr",
        "pos_col":      "PosGRCh37",
        "ea_col":       "testedAllele",
        "oa_col":       "otherAllele",
        "p_col":        "p",
        "stat_cols":    ["z", "p", "N"],
        "out_stem":     "wightman_with_ukb_GRCh38",
        "label":        "Wightman 2021 (incl UKB, excl 23&Me)",
    },
}


# ── CLI ───────────────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser(
    description="Liftover Wightman 2021 summary stats GRCh37 → GRCh38.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
_p.add_argument(
    "--with-ukb", action="store_true", dest="with_ukb",
    help="Use the file INCLUDING UK Biobank (excl 23&Me only). "
         "Default: file excluding BOTH UKB and 23&Me.",
)
_p.add_argument(
    "--gw-thresh", type=float, default=5e-8, dest="gw_thresh",
    help="Genome-wide significance threshold for filtering (default: 5e-8).",
)
_p.add_argument(
    "--all-snps", action="store_true", dest="all_snps",
    help="Lift ALL SNPs in the file, not just genome-wide significant ones. "
         "WARNING: ~12M rows; expect ~45 min runtime.",
)
_p.add_argument(
    "--out-dir", default=str(OUT_DIR), dest="out_dir",
    help="Output directory (default: same as input GWAS folder).",
)
_p.add_argument(
    "--fa", default=str(GRCH38_FA), dest="fa",
    help="Path to GRCh38 primary-assembly FASTA (.fa or .fa.gz). "
         "Required for REF_38 annotation.",
)
args = _p.parse_args()

fmt_key = "incl_ukb" if args.with_ukb else "excl_ukb"
FMT     = FORMATS[fmt_key]
OUT_DIR = pathlib.Path(args.out_dir)
OUT_DIR.mkdir(parents=True, exist_ok=True)
LIFT_DIR.mkdir(parents=True, exist_ok=True)
GRCH38_FA = pathlib.Path(args.fa)

SEP = "=" * 65
print(SEP)
print(f"  Wightman 2021 GRCh37 → GRCh38 Liftover")
print(f"  Source  : {FMT['label']}")
print(f"  File    : {FMT['file'].name}")
print(f"  Filter  : {'ALL SNPs' if args.all_snps else f'GW-sig only (p < {args.gw_thresh:.0e})'}")
print(f"  Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
print(SEP)


# ── Step 1: Check input file ───────────────────────────────────────────────────
if not FMT["file"].exists():
    print(f"\n[ERROR] Input file not found: {FMT['file']}")
    print("  Make sure the Wightman files are unzipped into their sub-folders.")
    sys.exit(1)


# ── Step 2: Download chain file ───────────────────────────────────────────────
print("\n[Step 1] Chain file")
if CHAIN_FILE.exists():
    print(f"  Already present ({CHAIN_FILE.stat().st_size / 1e6:.1f} MB): {CHAIN_FILE}")
else:
    print(f"  Downloading hg19→hg38 chain file from UCSC (~20 MB)...")

    def _progress(block, block_size, total):
        pct = min(100, block * block_size * 100 // total)
        print(f"\r  {pct}% ", end="", flush=True)

    urllib.request.urlretrieve(CHAIN_URL, CHAIN_FILE, reporthook=_progress)
    print(f"\n  Saved → {CHAIN_FILE}")


# ── Step 3: Load summary stats ────────────────────────────────────────────────
print(f"\n[Step 2] Load summary statistics")
print(f"  Reading: {FMT['file']}")
print(f"  (This may take a minute for ~12M rows...)")

chr_col = FMT["chr_col"]
pos_col = FMT["pos_col"]
ea_col  = FMT["ea_col"]
oa_col  = FMT["oa_col"]
p_col   = FMT["p_col"]

keep_cols = [chr_col, pos_col, ea_col, oa_col] + [
    c for c in FMT["stat_cols"] if c not in [chr_col, pos_col, ea_col, oa_col]
]

df = pd.read_csv(FMT["file"], sep=FMT["sep"])
df.columns = df.columns.str.strip()   # strip leading/trailing spaces from headers
print(f"  Loaded {len(df):,} SNPs. Columns: {list(df.columns)}")

# Coerce numeric columns
df[p_col] = pd.to_numeric(df[p_col], errors="coerce")
df[pos_col] = pd.to_numeric(df[pos_col], errors="coerce").astype("Int64")
df[chr_col] = df[chr_col].astype(str).str.replace("^chr", "", regex=True).str.strip()

# Filter to GW-significant if requested
if not args.all_snps:
    n_before = len(df)
    df = df[df[p_col] < args.gw_thresh].copy()
    print(f"  Filtered to GW-sig (p < {args.gw_thresh:.0e}): "
          f"{len(df):,} / {n_before:,} SNPs retained")
else:
    print(f"  Keeping all {len(df):,} SNPs (--all-snps)")

df = df.dropna(subset=[chr_col, pos_col]).copy()
df[pos_col] = df[pos_col].astype(int)
print(f"  After dropping missing positions: {len(df):,} SNPs")


# ── Step 4: Run pyliftover ────────────────────────────────────────────────────
print(f"\n[Step 3] Run pyliftover")
try:
    from pyliftover import LiftOver
except ImportError:
    print("\n[ERROR] pyliftover not installed. Run:  pip install pyliftover")
    sys.exit(1)

print("  Loading chain file (first run builds index, ~30 s)...")
lo = LiftOver(str(CHAIN_FILE))
print("  Chain file loaded.")

n = len(df)
new_chrs, new_bps, lifted_flags = [], [], []

print(f"  Lifting {n:,} positions...")
for i, (_, row) in enumerate(df.iterrows()):
    if i % 5_000 == 0 and i > 0:
        print(f"    {i:>7,} / {n:,}  ({100*i/n:.0f}%)", flush=True)

    chrom_in = f"chr{row[chr_col]}"
    pos_in   = int(row[pos_col]) - 1   # 1-based → 0-based

    result = lo.convert_coordinate(chrom_in, pos_in)

    if result is None or len(result) == 0:
        new_chrs.append(None)
        new_bps.append(None)
        lifted_flags.append(False)
    else:
        new_chrom, new_pos, strand, score = result[0]
        chrom_clean = new_chrom.replace("chr", "")
        if chrom_clean not in CANONICAL:
            new_chrs.append(None)
            new_bps.append(None)
            lifted_flags.append(False)
        else:
            new_chrs.append(chrom_clean)
            new_bps.append(int(new_pos) + 1)   # 0-based → 1-based
            lifted_flags.append(True)

print(f"    {n:>7,} / {n:,}  (100%)")

df = df.copy()
df["CHR_38"]  = new_chrs
df["BP_38"]   = new_bps
df["lifted"]  = lifted_flags

n_lifted   = sum(lifted_flags)
n_unmapped = n - n_lifted
print(f"\n  ── Liftover summary ──")
print(f"  Total input         : {n:>8,}")
print(f"  Successfully lifted : {n_lifted:>8,}")
print(f"  Unmapped / dropped  : {n_unmapped:>8,}  ({100*n_unmapped/max(n,1):.2f}%)")


# ── Step 5: Fetch REF allele — pyfaidx (local FASTA) primary, REST API fallback
print(f"\n[Step 4] Annotate REF_38 from GRCh38 FASTA (pyfaidx primary → Ensembl REST fallback)")

df["REF_38"] = None
lifted_df = df[df["lifted"]].copy()
n_fa = len(lifted_df)

if n_fa > 0:
    # ── Tier 1: pyfaidx on local GRCh38 FASTA ────────────────────────────────
    if GRCH38_FA.exists():
        try:
            from pyfaidx import Fasta
            print(f"  Opening FASTA: {GRCH38_FA}")
            fa = Fasta(str(GRCH38_FA))
            fa_keys = set(fa.keys())
            refs_list = []
            print(f"  Fetching REF for {n_fa:,} lifted positions...")
            for i, (_, row) in enumerate(lifted_df.iterrows()):
                if i % 500 == 0 and i > 0:
                    print(f"    {i:>5,} / {n_fa:,}", flush=True)
                chrom = str(row["CHR_38"])
                bp    = int(row["BP_38"])
                ref   = None
                for key in [chrom, f"chr{chrom}"]:
                    if key in fa_keys:
                        try:
                            ref = fa[key][bp - 1:bp].seq.upper()
                        except Exception:
                            pass
                        break
                refs_list.append(ref)
            df.loc[df["lifted"], "REF_38"] = refs_list
            n_found = df["REF_38"].notna().sum()
            print(f"  REF_38 annotated via pyfaidx: {n_found:,} / {n_fa:,}")
        except Exception as e:
            print(f"  [WARN] pyfaidx failed: {e}")
            n_found = 0
    else:
        print(f"  [WARN] GRCh38 FASTA not found: {GRCH38_FA}")
        n_found = 0

    # ── Tier 2: Ensembl REST API for any remaining missing positions ──────────
    n_missing = n_fa - n_found
    if n_missing > 0:
        if _requests is None:
            print(f"  [WARN] 'requests' not installed — skipping REST fallback for {n_missing:,} positions.")
        else:
            print(f"  Ensembl REST fallback for {n_missing:,} remaining positions...")
            ENSEMBL_URL = "https://rest.ensembl.org/sequence/region/human"
            HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}
            missing_rows = df[df["lifted"] & df["REF_38"].isna()].copy()
            positions = [(str(r["CHR_38"]), int(r["BP_38"]), idx)
                         for idx, r in missing_rows.iterrows()]
            BATCH = 50   # Ensembl REST /sequence/region limit is 50 per POST
            for start in range(0, len(positions), BATCH):
                batch = positions[start:start + BATCH]
                regions = [f"{c}:{p}..{p}" for c, p, _ in batch]
                try:
                    resp = _requests.post(ENSEMBL_URL,
                                          json={"regions": regions},
                                          headers=HEADERS, timeout=60)
                    resp.raise_for_status()
                    ref_map = {}
                    for item in resp.json():
                        q = item.get("query", "")
                        s = item.get("seq", "").upper()
                        if q and s:
                            ref_map[f"{q.split(':')[0]}:{q.split(':')[1].split('..')[0]}"] = s
                    for c, p, idx in batch:
                        val = ref_map.get(f"{c}:{p}")
                        if val:
                            df.at[idx, "REF_38"] = val
                except Exception as e:
                    print(f"  [WARN] REST batch {start//BATCH+1} failed: {e}")
                if start + BATCH < len(positions):
                    time.sleep(0.3)
            n_after = df["REF_38"].notna().sum()
            print(f"  After REST fallback: {n_after:,} / {n_fa:,} have REF_38")
else:
    print("  No lifted positions to annotate.")


# ── Step 5b: Look up rsIDs from ADNI hg19 BIM (match by hg19 chr:pos) ────────
print(f"\n[Step 4b] Look up rsIDs via ADNI hg19 BIM")
df["rsID"] = None

if ADNI_BIM_HG19.exists():
    bim19 = pd.read_csv(ADNI_BIM_HG19, sep="\t", header=None,
                        names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
                        dtype={"CHR": str, "BP": int})
    bim19_lookup = {(r.CHR, r.BP): r.SNP for r in bim19.itertuples(index=False)}
    print(f"  ADNI hg19 BIM loaded: {len(bim19_lookup):,} positions")

    rsids = []
    for _, row in df.iterrows():
        if not row.get("lifted"):
            rsids.append(None)
            continue
        # Use hg19 chr:pos for lookup (original Wightman coordinates)
        key = (str(row[chr_col if chr_col in df.columns else "CHR_37"]),
               int(row[pos_col if pos_col in df.columns else "BP_37"]))
        rsids.append(bim19_lookup.get(key))

    df["rsID"] = rsids
    # Tag source for provenance
    df["rsID_source"] = df["rsID"].apply(
        lambda x: "ADNI_hg19_BIM" if pd.notna(x) else None)
    n_rsid = df["rsID"].notna().sum()
    print(f"  rsID matched via hg19 BIM: {n_rsid:,} / {sum(df['lifted']):,} lifted SNPs")
    if n_rsid == 0:
        print("  [NOTE] No hg19 position overlap with ADNI array found.")
        print("         These Wightman lead SNPs may not be directly genotyped.")
else:
    df["rsID"] = None
    df["rsID_source"] = None
    print(f"  [WARN] ADNI hg19 BIM not found: {ADNI_BIM_HG19}")


# ── Step 4c: Ensembl GRCh37 REST fallback (position + strand-aware allele match) ──
need_ens = df[df["rsID"].isna() & df["lifted"]]
print(f"\n[Step 4c] Ensembl GRCh37 fallback for {len(need_ens):,} positions without rsID")

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")
def _allele_ok_ensembl(ea, oa, ens_alleles):
    """True if EA (or its complement) is found in Ensembl allele list."""
    if not ens_alleles:
        return False
    ens_set = {a.upper() for a in ens_alleles if a}
    ea_u = ea.upper()
    return ea_u in ens_set or ea_u.translate(_COMP) in ens_set

if len(need_ens) > 0:
    try:
        import requests, time as _time
        # GRCh37 mirror only supports GET /overlap/region, not batch POST
        ENSEMBL37_BASE = "https://grch37.rest.ensembl.org/overlap/region/human"
        ensembl_rsids = {}   # {df_index: rsID}
        idxs = need_ens.index.tolist()
        n_ok = 0

        for i, idx in enumerate(idxs):
            row   = df.loc[idx]
            chrom = str(row[chr_col] if chr_col in df.columns else row["CHR_37"])
            pos   = int(row[pos_col] if pos_col in df.columns else row["BP_37"])
            ea    = str(row.get(ea_col, "") or "").upper()
            oa    = str(row.get(oa_col, "") or "").upper()
            try:
                resp = requests.get(
                    f"{ENSEMBL37_BASE}/{chrom}:{pos}-{pos}",
                    params={"feature": "variation"},
                    headers={"Accept": "application/json"},
                    timeout=20,
                )
                resp.raise_for_status()
                for v in resp.json():
                    if not v.get("id", "").startswith("rs"):
                        continue
                    if _allele_ok_ensembl(ea, oa, v.get("alleles", [])):
                        ensembl_rsids[idx] = v["id"]
                        n_ok += 1
                        break
            except Exception as e:
                pass   # silent — majority may genuinely have no entry
            # ~14 req/s to stay within Ensembl rate limit
            _time.sleep(0.07)
            if (i + 1) % 500 == 0:
                print(f"    ... {i+1:,} / {len(idxs):,} queried, {n_ok:,} rsIDs found")

        for idx, rsid in ensembl_rsids.items():
            df.at[idx, "rsID"]        = rsid
            df.at[idx, "rsID_source"] = "Ensembl_GRCh37"
        print(f"  Ensembl rsIDs found: {len(ensembl_rsids):,} / {len(need_ens):,}")
    except ImportError:
        print("  [WARN] requests not installed — skipping Ensembl fallback")
else:
    print("  All lifted SNPs already have rsIDs from ADNI hg19 BIM.")


# ── Step 4d: RsMergeArch canonicalization ──────────────────────────────────────
print(f"\n[Step 4d] Canonicalize rsIDs via RsMergeArch")
df["rsID_canonical"] = df["rsID"].copy()
if NCBI_MERGE_ARCH.exists():
    rsid_ints = set()
    for rsid in df["rsID"].dropna():
        num = str(rsid).replace("rs", "").strip()
        if num.isdigit():
            rsid_ints.add(int(num))
    if rsid_ints:
        merge_map = {}   # {low_rsid_int: high_rsid_int}
        import gzip as _gz
        with _gz.open(NCBI_MERGE_ARCH, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = line.split("\t")
                if len(parts) >= 2:
                    try:
                        high, low = int(parts[0]), int(parts[1])
                        if low in rsid_ints:
                            merge_map[low] = high
                    except ValueError:
                        continue
        n_merged = 0
        for idx, rsid in df["rsID"].items():
            if pd.isna(rsid):
                continue
            num = str(rsid).replace("rs", "").strip()
            if num.isdigit():
                low = int(num)
                if low in merge_map:
                    df.at[idx, "rsID_canonical"] = f"rs{merge_map[low]}"
                    n_merged += 1
        print(f"  Deprecated rsIDs updated: {n_merged:,} of {len(rsid_ints):,}")
    else:
        print("  No rsIDs to canonicalize.")
else:
    print(f"  [WARN] RsMergeArch not found: {NCBI_MERGE_ARCH}")
    print("  rsID_canonical = rsID (no merge resolution applied)")


# ── Step 6: Build SNPID_38 and tidy output ────────────────────────────────────
print(f"\n[Step 5] Build output")

# Rename source columns to canonical names
df = df.rename(columns={
    chr_col: "CHR_37",
    pos_col: "BP_37",
    ea_col:  "effect_allele",
    oa_col:  "other_allele",
    p_col:   "p_value",
})

# Build GRCh38 positional SNPID: CHR:BP:REF:EA
def _build_snpid(row):
    if not row.get("lifted"):
        return None
    chrom = row.get("CHR_38")
    bp    = row.get("BP_38")
    ref   = row.get("REF_38") or "?"
    ea    = str(row.get("effect_allele", "?")).upper()
    return f"{chrom}:{int(bp)}:{ref}:{ea}"

df["SNPID_38"] = df.apply(_build_snpid, axis=1)

# Column order
fixed_cols = [
    "CHR_37", "BP_37", "CHR_38", "BP_38",
    "rsID", "rsID_canonical", "rsID_source",
    "effect_allele", "other_allele", "REF_38",
    "SNPID_38", "lifted",
]
# Add remaining stat columns that exist in df
stat_cols_present = [c for c in ["beta", "standard_error",
                                  "effect_allele_frequency",
                                  "p_value", "z", "N", "Neffective", "SNPID"]
                     if c in df.columns]
out_cols = fixed_cols + stat_cols_present

out_df = df[[c for c in out_cols if c in df.columns]].copy()

# Save
out_path = OUT_DIR / f"{FMT['out_stem']}.tsv"
out_df.to_csv(out_path, sep="\t", index=False)

# Trusted TSV: exact ADNI hg19 BIM positional matches only
trusted_df = out_df[out_df["rsID_source"] == "ADNI_hg19_BIM"].copy()
trusted_path = OUT_DIR / f"{FMT['out_stem']}_trusted_adni_hg19.tsv"
trusted_df.to_csv(trusted_path, sep="\t", index=False)
print(f"  Trusted (ADNI hg19 exact): {len(trusted_df):,} SNPs → {trusted_path.name}")

n_out = len(out_df)
n_out_lifted = out_df["lifted"].sum()
print(f"\n  Rows written          : {n_out:,}")
print(f"  Successfully lifted   : {n_out_lifted:,}")
print(f"  With REF_38 annotated : {out_df['REF_38'].notna().sum():,}")

# Sample of lifted rows
if n_out_lifted > 0:
    sample = out_df[out_df["lifted"]].head(5)
    print(f"\n  Sample of lifted GW-sig SNPs:")
    for _, r in sample.iterrows():
        print(f"    {r['CHR_37']}:{r['BP_37']} → "
              f"{r['CHR_38']}:{r['BP_38']}  "
              f"EA={r['effect_allele']}  OA={r['other_allele']}  "
              f"REF={r['REF_38']}  p={r.get('p_value','?'):.2e}")

print(f"\n  Output: {out_path}")
print(f"\nDone!")
