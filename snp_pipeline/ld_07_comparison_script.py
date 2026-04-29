"""
ld_07_comparison_script.py
===========================
Full SNP cleaning pipeline for an LD-pruned (r² ≤ 0.7) comparison dataset.

Pipeline stages
---------------
  0. LD pruning          plink --indep-pairwise 50 5 0.7  →  ld_pruned_07_SNP_filtered_with_mri
  1. rsID assignment     reuse existing update_name.txt   →  ld_pruned_07_..._rsid
  2. Dedup               resolve_preexisting_dups.py      →  ld_pruned_07_..._rsid_deduped.bim
  3. Clean               plink --exclude palindromic list →  ld_pruned_07_..._rsid_clean
  4. Deprecation patch   validate_rsid.py + --update-name →  ld_pruned_07_..._rsid_clean_current
  5. Liftover            liftover_hg19_to_hg38.py         →  ld_pruned_07_..._GRCh38
  6. REF correction      plink2 --ref-from-fa             →  ld_pruned_07_..._GRCh38_refcorr

All output files are prefixed with 'ld_pruned_07_' so they never overwrite
the main pipeline outputs (patch_deprecated_safe.txt, etc.).

Usage
-----
  python ld_07_comparison_script.py

  # Resume from a specific step (0-6):
  python ld_07_comparison_script.py --start-step 3

  # Override PLINK executable path:
  python ld_07_comparison_script.py --plink D:/ADNI_SNP_Omni2.5M_20140220/plink.exe
"""

import argparse
import io
import pathlib
import shutil
import subprocess
import sys

# Force UTF-8 output on Windows (avoids cp1252 UnicodeEncodeError)
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ── Configuration ─────────────────────────────────────────────────────────────
BASE_DIR   = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220")
SCRIPT_DIR = pathlib.Path(__file__).parent   # snp_pipeline/
PFX        = "ld_pruned_07_"                 # prefix for all outputs

# Input dataset (before LD pruning — must exist)
INPUT_STEM = BASE_DIR / "SNP_filtered_with_mri"

# All intermediate and final stems
LD_STEM          = BASE_DIR / f"{PFX}SNP_filtered_with_mri"
RSID_STEM        = BASE_DIR / f"{PFX}SNP_filtered_with_mri_rsid"
CLEAN_STEM       = BASE_DIR / f"{PFX}SNP_filtered_with_mri_rsid_clean"
CURRENT_STEM     = BASE_DIR / f"{PFX}SNP_filtered_with_mri_rsid_clean_current"
GRCH38_STEM      = BASE_DIR / f"{PFX}SNP_filtered_with_mri_rsid_clean_current_GRCh38"
REFCORR_STEM     = BASE_DIR / f"{PFX}SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr"

# Reused files from the main pipeline (manifest-level, not BIM-specific)
UPDATE_NAME_FILE = BASE_DIR / "update_name.txt"          # kgp → rsID map (reuse existing)

# LD pruning parameters  (kb-based window for consistent physical coverage)
LD_WINDOW   = "100kb"  # 100 kilobase physical window
LD_STEP     = 1        # step size in SNP count (1 = maximally thorough)
LD_R2       = 0.7      # r² threshold




# Prefixed output files from resolve_preexisting_dups.py
EXCL_PALINDROMIC = BASE_DIR / f"{PFX}exclude_palindromic_dups.txt"
COLLAPSED_LIST   = BASE_DIR / f"{PFX}collapsed_identical_dups.txt"
DEDUPED_BIM      = BASE_DIR / f"{PFX}SNP_filtered_with_mri_rsid_deduped.bim"

# Prefixed output files from validate_rsid.py
PATCH_SAFE     = BASE_DIR / f"{PFX}patch_deprecated_safe.txt"

# GRCh38 reference FASTA for --ref-from-fa
GRCH38_FA = BASE_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

# ── CLI ───────────────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser(
    description="LD-pruned (r²=0.7) SNP pipeline for comparison analysis.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
_p.add_argument("--start-step", type=int, default=0, metavar="N",
                help="Resume from step N (0=LD prune, 1=rsID, 2=dedup, 3=clean, "
                     "4=deprecation, 5=liftover, 6=ref-from-fa). Default: 0")
_p.add_argument("--plink", default=None,
                help="Path to PLINK 1.9 executable (auto-detected from PATH if omitted).")
_p.add_argument("--plink2", default=None,
                help="Path to PLINK2 executable (auto-detected from PATH if omitted).")
_p.add_argument("--dry-run", action="store_true",
                help="Print commands without executing them.")
args = _p.parse_args()

# ── Resolve executables ───────────────────────────────────────────────────────
PLINK_EXE  = (args.plink  or shutil.which("plink")
              or str(BASE_DIR / "plink.exe")  if (BASE_DIR / "plink.exe").exists() else None)
PLINK2_EXE = (args.plink2 or shutil.which("plink2"))

if PLINK_EXE is None:
    print("[ERROR] plink 1.9 not found on PATH and not at D:/ADNI.../plink.exe")
    print("  Pass the path with --plink <path>")
    sys.exit(1)


# ── Helper functions ──────────────────────────────────────────────────────────
def banner(step: int, title: str):
    print(f"\n{'═'*70}")
    print(f"  Step {step}: {title}")
    print(f"{'═'*70}")


def run(cmd: list):
    """Run a PLINK/shell command, capture output, print last 30 lines."""
    cmd_str = " ".join(str(c) for c in cmd)
    print(f"\n  $ {cmd_str}")
    if args.dry_run:
        return
    result = subprocess.run(cmd, capture_output=True, text=True,
                            encoding="utf-8", errors="replace")
    stdout_text = result.stdout or ""
    stderr_text = result.stderr or ""
    if stdout_text.strip():
        for line in stdout_text.strip().splitlines()[-30:]:
            print(f"    {line}")
    if result.returncode != 0:
        print(f"\n[ERROR] Command failed (exit {result.returncode}):")
        print(stderr_text[-3000:])
        sys.exit(result.returncode)


def run_python(script: pathlib.Path, extra_args: list = None):
    """Run a Python sub-script, streaming output live (no capture)."""
    cmd = [sys.executable, "-u", str(script)] + (extra_args or [])
    cmd_str = " ".join(str(c) for c in cmd)
    print(f"\n  $ {cmd_str}")
    if args.dry_run:
        return
    # Stream output directly to terminal — avoids encoding capture issues
    # and gives real-time progress for long-running steps.
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"\n[ERROR] Script failed (exit {result.returncode})")
        sys.exit(result.returncode)


def skip(step: int) -> bool:
    return step < args.start_step


def check_exists(path: pathlib.Path, label: str):
    if not path.exists():
        print(f"[ERROR] Expected file not found: {path}")
        print(f"  ({label})")
        sys.exit(1)


# ── Step 0: LD pruning ────────────────────────────────────────────────────────
STEP = 0
banner(STEP, f"LD pruning (r² ≤ {LD_R2}, window={LD_WINDOW}, step={LD_STEP})")

if skip(STEP):
    print(f"  Skipping (--start-step {args.start_step})")
else:
    prune_prefix = BASE_DIR / f"{PFX}ld_prune"

    # Phase 1: identify variants to keep
    run([PLINK_EXE,
         "--bfile",            str(INPUT_STEM),
         "--indep-pairwise",   str(LD_WINDOW), str(LD_STEP), str(LD_R2),
         "--out",              str(prune_prefix)])

    prune_in  = prune_prefix.with_suffix(".prune.in")
    check_exists(prune_in, "plink --indep-pairwise .prune.in list")

    # Phase 2: extract pruned variants
    run([PLINK_EXE,
         "--bfile",   str(INPUT_STEM),
         "--extract", str(prune_in),
         "--make-bed",
         "--out",     str(LD_STEM)])

    check_exists(LD_STEM.with_suffix(".bim"), "LD-pruned BIM")
    print(f"\n  ✓ LD-pruned dataset: {LD_STEM}")


# ── Step 1: rsID assignment (reuse existing update_name.txt) ──────────────────
STEP = 1
banner(STEP, "rsID assignment (reuse update_name.txt from main pipeline)")

if skip(STEP):
    print(f"  Skipping (--start-step {args.start_step})")
else:
    check_exists(UPDATE_NAME_FILE,
                 "update_name.txt from previous update_snp_ids.py run")

    run([PLINK_EXE,
         "--bfile",       str(LD_STEM),
         "--update-name", str(UPDATE_NAME_FILE),
         "--make-bed",
         "--out",         str(RSID_STEM)])

    check_exists(RSID_STEM.with_suffix(".bim"), "rsID-updated BIM")
    print(f"\n  ✓ rsIDs applied: {RSID_STEM}")


# ── Step 2: Resolve pre-existing duplicate rsIDs ──────────────────────────────
STEP = 2
banner(STEP, "Resolve pre-existing duplicate rsIDs (resolve_preexisting_dups.py)")

if skip(STEP):
    print(f"  Skipping (--start-step {args.start_step})")
else:
    check_exists(RSID_STEM.with_suffix(".bim"), "rsID BIM for dedup input")

    run_python(SCRIPT_DIR / "resolve_preexisting_dups.py", [
        "--bim",    str(RSID_STEM.with_suffix(".bim")),
        "--outdir", str(BASE_DIR),
        "--prefix", PFX,
    ])

    check_exists(EXCL_PALINDROMIC, "exclude_palindromic_dups.txt")
    check_exists(DEDUPED_BIM,      "deduped BIM (for copying alongside bed/fam)")

    # Copy .bed and .fam (unchanged) to the deduped stem so PLINK can use it
    deduped_stem = DEDUPED_BIM.with_suffix("")
    for ext in [".bed", ".fam"]:
        src = RSID_STEM.with_suffix(ext)
        dst = deduped_stem.with_suffix(ext)
        if not dst.exists() or args.dry_run:
            print(f"  Copying {src.name} → {dst.name}")
            if not args.dry_run:
                shutil.copy2(src, dst)

    print(f"\n  ✓ Dedup outputs in {BASE_DIR} with prefix '{PFX}'")


# ── Step 3: Exclude palindromic reverse-strand probes ────────────────────────
STEP = 3
banner(STEP, "Exclude palindromic reverse-strand probes → _clean")

if skip(STEP):
    print(f"  Skipping (--start-step {args.start_step})")
else:
    deduped_stem = DEDUPED_BIM.with_suffix("")
    check_exists(EXCL_PALINDROMIC, "palindromic exclude list")

    exclude_lists = [str(EXCL_PALINDROMIC)]

    # Also exclude collapsed identical duplicates if present
    if COLLAPSED_LIST.exists():
        n_collapsed = sum(1 for _ in open(COLLAPSED_LIST) if _.strip())
        if n_collapsed > 0:
            exclude_lists.append(str(COLLAPSED_LIST))
            print(f"  Including {n_collapsed} collapsed identical duplicates in exclusion.")

    # PLINK only takes one --exclude file; merge them if there are two
    if len(exclude_lists) == 2:
        merged_excl = BASE_DIR / f"{PFX}exclude_all_dups.txt"
        print(f"  Merging exclude lists → {merged_excl.name}")
        if not args.dry_run:
            with open(merged_excl, "w") as fout:
                for path in exclude_lists:
                    with open(path) as fin:
                        fout.write(fin.read())
        excl_arg = str(merged_excl)
    else:
        excl_arg = exclude_lists[0]

    run([PLINK_EXE,
         "--bfile",   str(deduped_stem),
         "--exclude", excl_arg,
         "--make-bed",
         "--out",     str(CLEAN_STEM)])

    check_exists(CLEAN_STEM.with_suffix(".bim"), "clean BIM")
    print(f"\n  ✓ Palindromic probes excluded: {CLEAN_STEM}")


# ── Step 4: Validate and patch deprecated rsIDs ───────────────────────────────
STEP = 4
banner(STEP, "Validate deprecated rsIDs → _clean_current")

if skip(STEP):
    print(f"  Skipping (--start-step {args.start_step})")
else:
    check_exists(CLEAN_STEM.with_suffix(".bim"), "clean BIM for validate_rsid input")

    run_python(SCRIPT_DIR / "validate_rsid.py", [
        "--bim",    str(CLEAN_STEM.with_suffix(".bim")),
        "--outdir", str(BASE_DIR),
        "--prefix", PFX,
    ])

    check_exists(PATCH_SAFE, "patch_deprecated_safe.txt")

    # Apply the safe renames
    run([PLINK_EXE,
         "--bfile",       str(CLEAN_STEM),
         "--update-name", str(PATCH_SAFE),
         "--make-bed",
         "--out",         str(CURRENT_STEM)])

    check_exists(CURRENT_STEM.with_suffix(".bim"), "current BIM")
    print(f"\n  ✓ Deprecated rsIDs patched: {CURRENT_STEM}")


# ── Step 5: Liftover GRCh37 → GRCh38 ─────────────────────────────────────────
STEP = 5
banner(STEP, "Liftover GRCh37 → GRCh38 (liftover_hg19_to_hg38.py)")

if skip(STEP):
    print(f"  Skipping (--start-step {args.start_step})")
else:
    check_exists(CURRENT_STEM.with_suffix(".bim"), "current BIM for liftover input")

    run_python(SCRIPT_DIR / "liftover_hg19_to_hg38.py", [
        "--input", str(CURRENT_STEM),
        "--fa",    str(GRCH38_FA),
        "--plink", PLINK_EXE,
    ])

    check_exists(GRCH38_STEM.with_suffix(".bim"), "GRCh38 BIM")
    print(f"\n  ✓ Lifted to GRCh38: {GRCH38_STEM}")


# ── Step 6: REF/ALT correction (plink2 --ref-from-fa) ────────────────────────
STEP = 6
banner(STEP, "REF/ALT correction (plink2 --ref-from-fa)")

if skip(STEP):
    print(f"  Skipping (--start-step {args.start_step})")
elif PLINK2_EXE is None:
    print("  plink2 not found on PATH — skipping --ref-from-fa.")
    print("  Run manually:")
    print(f'    plink2 --bfile "{GRCH38_STEM}" \\')
    print(f'           --fa "{GRCH38_FA}" \\')
    print(f'           --ref-from-fa force \\')
    print(f'           --make-bed \\')
    print(f'           --out "{REFCORR_STEM}"')
elif not GRCH38_FA.exists():
    print(f"  GRCh38 FASTA not found: {GRCH38_FA} — skipping --ref-from-fa.")
else:
    check_exists(GRCH38_STEM.with_suffix(".bim"), "GRCh38 BIM for ref correction")

    run([PLINK2_EXE,
         "--bfile",       str(GRCH38_STEM),
         "--fa",          str(GRCH38_FA),
         "--ref-from-fa", "force",
         "--make-bed",
         "--out",         str(REFCORR_STEM)])

    check_exists(REFCORR_STEM.with_suffix(".bim"), "refcorr BIM")
    print(f"\n  ✓ REF/ALT corrected: {REFCORR_STEM}")


# ── Done ──────────────────────────────────────────────────────────────────────
print(f"\n{'═'*70}")
print("  Pipeline complete!")
print(f"{'═'*70}")
print(f"\n  Final dataset: {REFCORR_STEM}")
print(f"  All outputs prefixed with: '{PFX}'")
print(f"\n  Key files:")
print(f"    BIM (GRCh38, ref-corrected): {REFCORR_STEM.with_suffix('.bim')}")
print(f"    Dedup report              : {BASE_DIR / (PFX + 'preexisting_dups_report.tsv')}")
print(f"    Multi-allelic liftover map: {BASE_DIR / (PFX + 'multiallelic_snps_grch37.tsv')}")
print(f"    Deprecated rsID patch     : {PATCH_SAFE}")
