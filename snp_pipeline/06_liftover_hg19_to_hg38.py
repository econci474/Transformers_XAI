"""
liftover_hg19_to_hg38.py
Convert PLINK binary dataset from GRCh37 (hg19) → GRCh38 (hg38).

Inputs
------
  D:/ADNI_SNP_Omni2.5M_20140220/SNP_filtered_with_mri_rsid_clean_current.{bed,bim,fam}

Outputs
-------
  D:/ADNI_SNP_Omni2.5M_20140220/SNP_filtered_with_mri_rsid_clean_current_GRCh38.{bed,bim,fam}
  D:/ADNI_SNP_Omni2.5M_20140220/liftover/unmapped_snps.txt
  D:/ADNI_SNP_Omni2.5M_20140220/liftover/multiallelic_id_rename_grch38.txt  (6 positional ID renames)

Method
------
Uses pyliftover (pip install pyliftover) — a pure-Python implementation of
the UCSC liftOver algorithm that works on Windows, Linux, and Mac with no
compiled dependencies. Uses the same UCSC hg19→hg38 chain file.

Why pyliftover instead of CrossMap?
  CrossMap requires pysam, which wraps htslib C code and cannot be built on
  Windows (causes "FileNotFoundError / '.' is not recognized" build errors).
  pyliftover is pure Python and pip-installs in seconds everywhere.

Usage
------
  pip install pyliftover pandas
  python liftover_hg19_to_hg38.py

  # If PLINK is not on PATH, point to it:
  python liftover_hg19_to_hg38.py --plink D:/ADNI_SNP_Omni2.5M_20140220/plink.exe
"""

import argparse
import io
import pathlib
import shutil
import subprocess
import sys
import urllib.request

# Force UTF-8 output on Windows (avoids cp1252 UnicodeEncodeError)
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import pandas as pd

# ── Paths (resolved from CLI flags at runtime — see __main__) ────────────────
BASE_DIR  = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220")

CHAIN_URL  = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"

# Table of multi-allelic positional IDs produced by resolve_preexisting_dups.py
MULTIALLELIC_TABLE = BASE_DIR / "multiallelic_snps_grch37.tsv"

# GRCh38 reference FASTA for plink2 --ref-from-fa (overridden by --fa in __main__)
GRCH38_FA = BASE_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

# These are set in __main__ after parsing CLI args
INPUT_STEM    = None
OUTPUT_STEM   = None
LIFT_DIR      = None
CHAIN_FILE    = None
UNMAPPED_SNPS = None
BIM_IN        = None
PLINK_EXE     = None

# Canonical chromosomes to keep (drop alt contigs like chr1_random)
CANONICAL = set([str(i) for i in range(1, 23)] + ["X", "Y", "M"])


# ─────────────────────────────────────────────────────────────────────────────
# 1. Download chain file if needed
# ─────────────────────────────────────────────────────────────────────────────
def download_chain():
    if CHAIN_FILE.exists():
        size_mb = CHAIN_FILE.stat().st_size / 1e6
        print(f"  Chain file already present ({size_mb:.1f} MB): {CHAIN_FILE}")
        return
    print(f"  Downloading hg19→hg38 chain file (~20 MB) from UCSC …")
    print(f"  URL: {CHAIN_URL}")

    def _progress(block, block_size, total):
        pct = min(100, block * block_size * 100 // total)
        print(f"\r  {pct}% ", end="", flush=True)

    urllib.request.urlretrieve(CHAIN_URL, CHAIN_FILE, reporthook=_progress)
    print(f"\n  Saved → {CHAIN_FILE}")


# ─────────────────────────────────────────────────────────────────────────────
# 2. Load BIM
# ─────────────────────────────────────────────────────────────────────────────
def load_bim() -> pd.DataFrame:
    print(f"  Reading BIM: {BIM_IN}")
    bim = pd.read_csv(BIM_IN, sep="\t", header=None,
                      names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
    n_before = len(bim)
    bim = bim[bim["BP"] > 0].copy()
    n_dropped = n_before - len(bim)
    if n_dropped:
        print(f"  Dropped {n_dropped:,} SNPs with BP=0 (unplaced)")
    print(f"  Loaded {len(bim):,} SNPs with known positions")
    return bim


# ─────────────────────────────────────────────────────────────────────────────
# 3. Run liftOver via pyliftover (pure Python, Windows-compatible)
# ─────────────────────────────────────────────────────────────────────────────
def run_liftover(bim: pd.DataFrame) -> pd.DataFrame:
    try:
        from pyliftover import LiftOver
    except ImportError:
        print("\n[ERROR] pyliftover not installed.")
        print("  Run:  pip install pyliftover\n")
        sys.exit(1)

    print(f"  Loading chain file into pyliftover (first run builds index, ~30 s) …")
    lo = LiftOver(str(CHAIN_FILE))
    print("  Chain file loaded.")

    # Vectorised-style loop with progress reporting
    n = len(bim)
    new_chrs, new_bps, failed = [], [], []

    print(f"  Lifting {n:,} SNPs …")
    for i, (_, row) in enumerate(bim.iterrows()):
        if i % 50_000 == 0:
            print(f"    {i:>7,} / {n:,}  ({100*i/n:.0f}%)", flush=True)

        # pyliftover uses 0-based coordinates
        chrom_in = f"chr{row['CHR']}"
        pos_in   = int(row["BP"]) - 1   # convert 1-based → 0-based

        result = lo.convert_coordinate(chrom_in, pos_in)

        if result is None or len(result) == 0:
            failed.append(row["SNP"])
            new_chrs.append(None)
            new_bps.append(None)
        else:
            # result[0] = (chrom, pos, strand, score)  — pos is 0-based
            new_chrom, new_pos, strand, score = result[0]
            chrom_clean = new_chrom.replace("chr", "")
            if chrom_clean not in CANONICAL:
                # Mapped to alt contig — treat as failed
                failed.append(row["SNP"])
                new_chrs.append(None)
                new_bps.append(None)
            else:
                new_chrs.append(chrom_clean)
                new_bps.append(int(new_pos) + 1)    # back to 1-based

    print(f"    {n:>7,} / {n:,}  (100%)")

    bim = bim.copy()
    bim["CHR_new"] = new_chrs
    bim["BP_new"]  = new_bps

    # Split into mapped / unmapped
    mapped   = bim[bim["CHR_new"].notna()].copy()
    unmapped = failed

    print(f"\n  ── LiftOver summary ──")
    print(f"  Total input SNPs    : {n:>10,}")
    print(f"  Successfully lifted : {len(mapped):>10,}")
    print(f"  Unmapped / dropped  : {len(unmapped):>10,}  ({100*len(unmapped)/n:.2f}%)")

    with open(UNMAPPED_SNPS, "w") as f:
        f.write("\n".join(unmapped))
    print(f"  Unmapped IDs saved  → {UNMAPPED_SNPS}")

    return mapped, unmapped


# ─────────────────────────────────────────────────────────────────────────────
# 4. Build updated BIM dataframe
# ─────────────────────────────────────────────────────────────────────────────
def build_new_bim(mapped: pd.DataFrame) -> pd.DataFrame:
    new_bim = pd.DataFrame({
        "CHR": mapped["CHR_new"],
        "SNP": mapped["SNP"],
        "CM":  mapped["CM"],
        "BP":  mapped["BP_new"].astype(int),
        "A1":  mapped["A1"],
        "A2":  mapped["A2"],
    })
    return new_bim


# ─────────────────────────────────────────────────────────────────────────────
# 5. Write GRCh38 PLINK files using PLINK --update-map / --update-chr
# ─────────────────────────────────────────────────────────────────────────────
def write_plink_outputs(new_bim: pd.DataFrame, unmapped: list):
    exclude_file    = LIFT_DIR / "unmapped_snps.txt"   # already written above
    update_map_file = LIFT_DIR / "update_map.txt"      # SNP  new_BP
    update_chr_file = LIFT_DIR / "update_chr.txt"      # SNP  new_CHR
    tmp_stem        = LIFT_DIR / "tmp_excluded"

    new_bim[["SNP", "BP"]].to_csv(update_map_file, sep="\t", header=False, index=False)
    new_bim[["SNP", "CHR"]].to_csv(update_chr_file, sep="\t", header=False, index=False)

    if PLINK_EXE is None:
        print("\n[WARNING] PLINK not found — cannot build binary files automatically.")
        print("  Run these commands manually (adjust paths if needed):\n")
        _print_manual_commands(exclude_file, update_map_file, update_chr_file, tmp_stem)
        return

    # Step 1: exclude unmapped SNPs
    print(f"\n  Running PLINK step 1: exclude {len(unmapped):,} unmapped SNPs …")
    _plink([
        "--bfile", str(INPUT_STEM),
        "--exclude", str(exclude_file),
        "--make-bed",
        "--out", str(tmp_stem),
    ])
    print(f"  → {tmp_stem}.bed/bim/fam written")

    # Step 2: update coordinates
    print(f"  Running PLINK step 2: update positions and chromosomes …")
    _plink([
        "--bfile", str(tmp_stem),
        "--update-map", str(update_map_file),
        "--update-chr", str(update_chr_file),
        "--make-bed",
        "--out", str(OUTPUT_STEM),
    ])

    print(f"\n✓ GRCh38 PLINK files written:")
    for ext in [".bed", ".bim", ".fam", ".log"]:
        p = OUTPUT_STEM.with_suffix(ext)
        if p.exists():
            print(f"    {p}  ({p.stat().st_size / 1e6:.1f} MB)")


def _plink(args: list):
    result = subprocess.run([PLINK_EXE] + args, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout[-2000:])
        print(result.stderr[-2000:])
        sys.exit(1)


def _print_manual_commands(exclude_file, update_map_file, update_chr_file, tmp_stem):
    exe = "plink"
    print(f"  {exe} --bfile {INPUT_STEM} \\")
    print(f"        --exclude {exclude_file} \\")
    print(f"        --make-bed \\")
    print(f"        --out {tmp_stem}\n")
    print(f"  {exe} --bfile {tmp_stem} \\")
    print(f"        --update-map {update_map_file} \\")
    print(f"        --update-chr {update_chr_file} \\")
    print(f"        --make-bed \\")
    print(f"        --out {OUTPUT_STEM}")


# ───────────────────────────────────────────────────────────────────────────────
# 6. Rename multi-allelic positional IDs from GRCh37 → GRCh38 format
# ───────────────────────────────────────────────────────────────────────────────
def rename_multiallelic_ids(new_bim: pd.DataFrame):
    """
    After liftover the 6 multi-allelic variants still carry their GRCh37
    positional IDs (e.g. '5:89956851:G:A').  We rename them to the
    corresponding GRCh38 positional ID (e.g. '5:90052501:G:A') using the
    coordinates now in the GRCh38 BIM.
    """
    if not MULTIALLELIC_TABLE.exists():
        print(f"  No multi-allelic table found at {MULTIALLELIC_TABLE} — skipping.")
        return

    multi = pd.read_csv(MULTIALLELIC_TABLE, sep="\t")
    # Columns: old_rsid, chrom, bp_grch37, ref_grch37, alt, new_id_grch37, strand, action

    # Build a lookup: GRCh37 positional ID → (ref, alt)
    grch37_ids = set(multi["new_id_grch37"].dropna())
    if not grch37_ids:
        print("  Multi-allelic table is empty — skipping ID rename.")
        return

    # Find which of these IDs survived liftover (appear in the GRCh38 BIM)
    # The GRCh38 BIM SNP column still has the old GRCh37 positional IDs at this point
    in_bim = new_bim[new_bim["SNP"].isin(grch37_ids)].copy()
    if in_bim.empty:
        print("  None of the multi-allelic positional IDs found in GRCh38 BIM — skipping.")
        return

    # ── Verify REF alleles against GRCh38 Ensembl ────────────────────────────
    # REF may change between assemblies; query the GRCh38 REST endpoint
    # for each lifted position before embedding REF in the new ID.
    import requests, time as _time
    ENSEMBL38_URL = "https://rest.ensembl.org/sequence/region/human"
    HEADERS38 = {"Content-Type": "application/json", "Accept": "application/json"}

    grch38_ref = {}  # "CHR:BP" → confirmed GRCh38 REF base
    positions38 = [f"{int(row['CHR'])}:{int(row['BP'])}..{int(row['BP'])}"
                   for _, row in in_bim.iterrows()]
    try:
        r38 = requests.post(ENSEMBL38_URL, json={"regions": positions38},
                            headers=HEADERS38, timeout=30)
        r38.raise_for_status()
        for item in r38.json():
            query = item.get("query", "")
            seq   = item.get("seq", "").upper()
            if query and ".." in query and seq:
                chrom_q = query.split(":")[0]
                pos_q   = query.split(":")[1].split("..")[0]
                grch38_ref[f"{chrom_q}:{pos_q}"] = seq
        print(f"  GRCh38 REF verified for {len(grch38_ref)}/{len(positions38)} positions.")
    except Exception as e:
        print(f"  [WARN] GRCh38 REF lookup failed ({e}); falling back to GRCh37 REF in IDs.")

    # ── Build rename map: old_grch37_id → new_grch38_id ──────────────────────
    rename_rows = []
    for _, row in in_bim.iterrows():
        old_id  = row["SNP"]           # e.g. '5:89956851:G:A'
        chrom   = str(int(row["CHR"])) # new GRCh38 chr (numeric)
        bp      = int(row["BP"])       # new GRCh38 position
        parts   = old_id.split(":")
        ref_37  = parts[2] if len(parts) >= 4 else "?"
        alt     = parts[3] if len(parts) >= 4 else "?"

        # Use GRCh38 REF if we got it; warn if it differs from GRCh37
        ref_38 = grch38_ref.get(f"{chrom}:{bp}", ref_37)
        if ref_38 != ref_37:
            print(f"  [NOTE] REF changed at {chrom}:{bp}: "
                  f"GRCh37={ref_37} → GRCh38={ref_38}  (updating ID)")

        new_id = f"{chrom}:{bp}:{ref_38}:{alt}" if ref_38 != "?" else f"{chrom}:{bp}"
        rename_rows.append({"old_id": old_id, "new_id": new_id,
                             "ref_grch37": ref_37, "ref_grch38": ref_38,
                             "ref_changed": ref_38 != ref_37})


    rename_df = pd.DataFrame(rename_rows)

    # Write rename file for PLINK --update-name (col1=old, col2=new)
    rename_file = LIFT_DIR / "multiallelic_id_rename_grch38.txt"
    rename_df[["old_id", "new_id"]].to_csv(rename_file, sep="\t", header=False, index=False)
    print(f"  Rename file written: {rename_file}  ({len(rename_df)} variants)")
    for _, r in rename_df.iterrows():
        print(f"    {r['old_id']}  →  {r['new_id']}")

    # Apply with PLINK --update-name if PLINK is available
    if PLINK_EXE is None:
        print("\n  [INFO] Run manually to apply the GRCh38 positional ID renames:")
        print(f'    plink --bfile "{OUTPUT_STEM}" \\')
        print(f'          --update-name "{rename_file}" \\')
        print(f'          --make-bed \\')
        print(f'          --out "{OUTPUT_STEM}"')
        return

    final_stem = OUTPUT_STEM  # overwrite in place
    tmp_rename = LIFT_DIR / "tmp_multiallelic_rename"
    _plink([
        "--bfile", str(OUTPUT_STEM),
        "--update-name", str(rename_file),
        "--make-bed",
        "--out", str(tmp_rename),
    ])
    # Replace OUTPUT_STEM with the renamed version
    for ext in [".bed", ".bim", ".fam", ".log"]:
        src = tmp_rename.with_suffix(ext)
        dst = final_stem.with_suffix(ext)
        if src.exists():
            src.replace(dst)
    print(f"  ✓ GRCh38 positional IDs applied to {final_stem}")


# ─────────────────────────────────────────────────────────────────────────────
# 7. Correct REF/ALT orientation for all variants using GRCh38 FASTA
# ─────────────────────────────────────────────────────────────────────────────
def _ref_from_fa():
    """
    Run plink2 --ref-from-fa to set REF alleles for all GRCh38 variants.
    Skipped gracefully if plink2 is not on PATH or FASTA is missing.

    Why needed:
      Liftover only updates CHR/BP coordinates. After liftover some variants
      may have A1/A2 columns swapped relative to the GRCh38 reference base.
      This step corrects that globally, required for VCF export, imputation,
      and comparison with GRCh38 reference panels.

    Caveats:
      - Palindromic SNPs (A/T or C/G) are ambiguous; plink2 skips them.
      - Requires plink2 (not plink 1.9).
      - FASTA chromosome names must match the BIM (Ensembl numeric naming
        '1, 2, ..., 22, X, Y' is used by both — no conversion needed).
    """
    plink2_exe = shutil.which("plink2")
    if plink2_exe is None:
        print("  plink2 not found on PATH — skipping --ref-from-fa.")
        print("  Run manually after installing plink2:")
        print(f'    plink2 --bfile "{OUTPUT_STEM}" \\')
        print(f'           --fa "{GRCH38_FA}" \\')
        print(f'           --ref-from-fa force \\')
        print(f'           --make-bed \\')
        print(f'           --out "{OUTPUT_STEM}_refcorr"')
        return

    if not GRCH38_FA.exists():
        print(f"  GRCh38 FASTA not found: {GRCH38_FA}")
        print("  Skipping --ref-from-fa. Provide the correct path with --fa.")
        return

    refcorr_stem = OUTPUT_STEM.parent / (OUTPUT_STEM.name + "_refcorr")
    print(f"  FASTA : {GRCH38_FA}")
    print(f"  Input : {OUTPUT_STEM}")
    print(f"  Output: {refcorr_stem}")
    print(f"  NOTE: Palindromic SNPs (A/T, C/G) are skipped — see .log for details.")

    cmd = [
        plink2_exe,
        "--bfile", str(OUTPUT_STEM),
        "--fa", str(GRCH38_FA),
        "--ref-from-fa", "force",
        "--make-bed",
        "--out", str(refcorr_stem),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    for line in (result.stdout + result.stderr).splitlines():
        if any(kw in line for kw in ["variant", "sample", "Warning", "Error",
                                      "ref-from-fa", "palindrom", "written", "time"]):
            print(f"  {line}")

    if result.returncode != 0:
        print(f"  [ERROR] plink2 --ref-from-fa failed. Full output:\n{result.stderr[-2000:]}")
    else:
        print(f"  ✓ REF/ALT corrected → {refcorr_stem}  (this is your final GRCh38 dataset)")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Lift a PLINK binary dataset from GRCh37 (hg19) to GRCh38 (hg38).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # default (LD-pruned dataset):\n"
            "  python liftover_hg19_to_hg38.py\n\n"
            "  # convert a different dataset:\n"
            "  python liftover_hg19_to_hg38.py --input D:/ADNI_SNP_Omni2.5M_20140220/SNP_filtered_with_mri\n\n"
            "  # specify PLINK explicitly:\n"
            "  python liftover_hg19_to_hg38.py --input D:/.../<stem> --plink D:/.../plink.exe"
        ),
    )
    parser.add_argument(
        "--input",
        default=str(BASE_DIR / "SNP_filtered_with_mri_rsid_clean_current"),
        help=(
            "Path to PLINK dataset stem (no extension). "
            "Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--plink",
        default=None,
        help="Path to PLINK executable (auto-detected from PATH or D:/.../plink.exe if omitted).",
    )
    parser.add_argument(
        "--fa",
        default=str(BASE_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"),
        help="Path to GRCh38 reference FASTA (.fa or .fa.gz) for --ref-from-fa. Default: %(default)s",
    )
    args = parser.parse_args()

    # ── Resolve all paths from --input ────────────────────────────────────────
    INPUT_STEM  = pathlib.Path(args.input)
    OUTPUT_STEM = INPUT_STEM.parent / (INPUT_STEM.name + "_GRCh38")
    LIFT_DIR    = INPUT_STEM.parent / "liftover"
    LIFT_DIR.mkdir(parents=True, exist_ok=True)
    CHAIN_FILE    = LIFT_DIR / "hg19ToHg38.over.chain.gz"
    UNMAPPED_SNPS = LIFT_DIR / "unmapped_snps.txt"
    BIM_IN        = INPUT_STEM.with_suffix(".bim")

    PLINK_LOCAL = BASE_DIR / "plink.exe"
    PLINK_EXE   = (args.plink
                   or shutil.which("plink")
                   or (str(PLINK_LOCAL) if PLINK_LOCAL.exists() else None))
    GRCH38_FA   = pathlib.Path(args.fa)

    if not BIM_IN.exists():
        print(f"[ERROR] BIM file not found: {BIM_IN}")
        print(f"  Check that --input points to a valid PLINK stem (without extension).")
        sys.exit(1)

    print("=" * 65)
    print("GRCh37 → GRCh38 LiftOver for PLINK dataset")
    print("Using: pyliftover (pure Python, Windows-compatible)")
    print("=" * 65)

    print("\n[Step 1] Download UCSC hg19→hg38 chain file")
    download_chain()

    print("\n[Step 2] Load BIM file")
    bim = load_bim()

    print("\n[Step 3] Run pyliftover on all SNP positions")
    mapped, unmapped = run_liftover(bim)

    print("\n[Step 4] Build updated BIM")
    new_bim = build_new_bim(mapped)

    print("\n[Step 5] Write GRCh38 PLINK files")
    write_plink_outputs(new_bim, unmapped)

    print("\n[Step 6] Rename multi-allelic positional IDs to GRCh38 coordinates")
    rename_multiallelic_ids(new_bim)

    print("\n[Step 7] Correct REF/ALT orientation using GRCh38 FASTA (plink2 --ref-from-fa)")
    _ref_from_fa()

    print("\nDone!")



