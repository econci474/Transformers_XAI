"""
resolve_preexisting_dups.py
═══════════════════════════
Finds rsIDs that already appear >1 times in the BIM (pre-existing duplicates),
classifies them as palindromic or multi-allelic via Ensembl GRCh37 REF lookup,
and produces a deduplicated BIM ready for the deprecation patch.

Classification:
  PALINDROMIC  — allele pairs are reverse complements → same SNP, different strand
                 → keep first occurrence (rsID retained), write rest to exclude list
  MULTI_ALLELIC — different ALT alleles after REF normalisation
                 → rename each to unique CHR:POS:A1:A2 in the BIM directly
  UNKNOWN_REF  — Ensembl returned no REF → treated as MULTI_ALLELIC (safe fallback)

Outputs:
  preexisting_dups_report.tsv     full classification of every duplicate row
  exclude_palindromic_dups.txt    variant IDs to pass to plink --exclude
  SNP_filtered_with_mri_rsid_deduped.bim   cleaned BIM (multi-allelic renamed in-place)

Pipeline (run BEFORE validate_rsid.py):
  # Step A: exclude palindromic dups, keep first occurrence
  plink --bfile "D:/ADNI_SNP_.../SNP_filtered_with_mri_rsid" \\
        --exclude "D:/ADNI_SNP_.../exclude_palindromic_dups.txt" \\
        --make-bed --out "D:/ADNI_SNP_.../SNP_filtered_with_mri_rsid_deduped"
  # (multi-allelic renames are already baked into the deduped BIM by this script)
"""

import argparse
import io
import pathlib
import sys
import time

# Force UTF-8 output on Windows (avoids cp1252 UnicodeEncodeError for → ═ etc.)
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import pandas as pd
import requests


# ── CLI (all args optional — defaults preserve original behaviour) ──────────────
_p = argparse.ArgumentParser(
    description="Resolve pre-existing duplicate rsIDs in a PLINK BIM.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
_p.add_argument("--bim",    default="D:/ADNI_SNP_Omni2.5M_20140220/SNP_filtered_with_mri_rsid.bim",
                help="Input BIM file (default: original pipeline path)")
_p.add_argument("--outdir", default="D:/ADNI_SNP_Omni2.5M_20140220",
                help="Output directory (default: D:/ADNI_SNP_Omni2.5M_20140220)")
_p.add_argument("--prefix", default="",
                help="Prefix for all output file names, e.g. 'ld_pruned_07_'")
_args = _p.parse_args()

# ── Paths ──────────────────────────────────────────────────────────────────────
BIM_IN   = pathlib.Path(_args.bim)
OUT_DIR  = pathlib.Path(_args.outdir)
PFX      = _args.prefix          # e.g. "ld_pruned_07_"

# Derive deduped BIM name from the input stem (stem already carries any prefix)
_bim_stem    = BIM_IN.stem       # e.g. "SNP_filtered_with_mri_rsid" or "ld_pruned_07_..._rsid"
BIM_OUT      = OUT_DIR / f"{_bim_stem}_deduped.bim"        # no extra PFX — stem has it
REPORT       = OUT_DIR / f"{PFX}preexisting_dups_report.tsv"
EXCL_OUT     = OUT_DIR / f"{PFX}exclude_palindromic_dups.txt"
MULTI_OUT    = OUT_DIR / f"{PFX}multiallelic_snps_grch37.tsv"   # for liftover use
COLLAPSED_OUT = OUT_DIR / f"{PFX}collapsed_identical_dups.txt"  # same REF:ALT after normalisation

ENSEMBL_URL = "https://grch37.rest.ensembl.org/sequence/region/human"
BATCH_SIZE  = 50
API_SLEEP   = 0.15   # seconds between batches (Ensembl rate limit: 15 req/s)

# ── Allele helpers ─────────────────────────────────────────────────────────────
_COMP    = {"A": "T", "T": "A", "C": "G", "G": "C"}
_MISSING = {"0", ".", "N", ""}


def comp(a):
    return _COMP.get(a.upper(), a.upper())


def normalize(ref, a1, a2):
    """
    Normalize a probe's alleles relative to REF.
    Returns (strand '+'/'-', alt_allele) or (None, None) if REF doesn't match.
    """
    ref = ref.upper()
    alleles = {x.upper() for x in [a1, a2] if x.upper() not in _MISSING}
    if not alleles:
        return None, None
    # Forward strand
    if ref in alleles:
        alt = (alleles - {ref}).pop() if len(alleles) > 1 else "."
        return "+", alt
    # Reverse strand — complement probe alleles and check again
    comp_alleles = {comp(a) for a in alleles}
    if ref in comp_alleles:
        alt = (comp_alleles - {ref}).pop() if len(comp_alleles) > 1 else "."
        return "-", alt
    return None, None   # REF doesn't match either strand


# ── Load BIM ──────────────────────────────────────────────────────────────────
print(f"Loading BIM: {BIM_IN.name} ...")
bim = pd.read_csv(BIM_IN, sep="\t", header=None,
                  names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
                  dtype=str)
for col in bim.columns:
    bim[col] = bim[col].str.strip()

print(f"  {len(bim):,} total variants")

counts  = bim["SNP"].value_counts()
dup_ids = counts[counts > 1].index.tolist()
dup_bim = bim[bim["SNP"].isin(dup_ids)].copy()

print(f"  {len(dup_ids):,} rsIDs appear >1 times  ({len(dup_bim):,} rows)")

if dup_bim.empty:
    print("No pre-existing duplicates — nothing to do.")
    sys.exit(0)

# ── Fetch REF alleles from Ensembl GRCh37 ─────────────────────────────────────
positions = dup_bim[["CHR", "BP"]].drop_duplicates()
regions   = [f"{r.CHR}:{r.BP}..{r.BP}" for _, r in positions.iterrows()]
key_to_ref = {}   # "CHR:BP" → REF base

print(f"\nFetching REF for {len(regions):,} positions from Ensembl GRCh37 ...")
n_batches = (len(regions) + BATCH_SIZE - 1) // BATCH_SIZE
_debug_printed = False
for b, start in enumerate(range(0, len(regions), BATCH_SIZE)):
    batch = regions[start:start + BATCH_SIZE]
    try:
        r = requests.post(
            ENSEMBL_URL,
            json={"regions": batch},
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            timeout=30,
        )
        r.raise_for_status()
        data = r.json()
        # Debug: print first item keys on first successful batch so we can verify format
        if not _debug_printed and data:
            print(f"\n  [DEBUG] First item keys: {list(data[0].keys())}")
            print(f"  [DEBUG] First item sample: { {k: str(v)[:60] for k, v in data[0].items()} }")
            _debug_printed = True
        for item in data:
            # Ensembl returns the request string in "query", format: "1:948921..948921"
            # (the "id" field contains the verbose coord-system format, not the request string)
            identifier = item.get("query", "")
            seq        = item.get("seq", "").upper()
            if identifier and ".." in identifier and seq:
                # Parse "CHR:START..END" — START == END for single-base queries
                chrom = identifier.split(":")[0]
                pos   = identifier.split(":")[1].split("..")[0]
                key_to_ref[f"{chrom}:{pos}"] = seq

    except Exception as e:
        print(f"\n  [WARN] batch {b+1}/{n_batches} failed: {e}")
    print(f"  batch {b+1}/{n_batches}", end="\r", flush=True)
    time.sleep(API_SLEEP)
print(f"  Done. REF retrieved for {len(key_to_ref):,}/{len(regions):,} positions.")


# ── Classify each duplicate group ─────────────────────────────────────────────
print("\nClassifying duplicate groups ...")

report_rows        = []
exclude_ids        = []   # (unused legacy — kept for safety)
bim_rename         = {}   # bim row index → new SNP name
collapsed_identical = []  # multi-allelic probes that collapsed to same REF:ALT

for rsid, group in dup_bim.groupby("SNP", sort=False):
    chrom = group.iloc[0]["CHR"]
    bp    = group.iloc[0]["BP"]
    ref   = key_to_ref.get(f"{chrom}:{bp}")

    # Normalise each row
    norm_results = []
    for idx, row in group.iterrows():
        if ref:
            strand, alt = normalize(ref, row["A1"], row["A2"])
        else:
            strand, alt = None, None
        norm_results.append({
            "bim_index": idx,
            "rsid": rsid, "chrom": chrom, "bp": bp,
            "a1": row["A1"], "a2": row["A2"],
            "ref": ref or "?",
            "strand": strand or "?",
            "alt": alt or "?",
        })

    # All normalize to the same REF:ALT  → palindromic duplicates
    unique_alts = {n["alt"] for n in norm_results if n["alt"] not in ("?", ".")}
    all_same    = (len(unique_alts) <= 1) and ref is not None

    if all_same:
        verdict = "PALINDROMIC"
        # Keep the forward-strand probe (strand == "+") — its alleles are in the
        # same orientation as dbSNP, so it is the canonical representation.
        # Fall back to the first row if no "+" strand is found.
        fwd_indices = [i for i, n in enumerate(norm_results) if n["strand"] == "+"]
        keep_idx    = fwd_indices[0] if fwd_indices else 0
        for i, n in enumerate(norm_results):
            is_keeper = (i == keep_idx)
            action = "KEEP (forward strand)" if is_keeper else "EXCLUDE (reverse-strand duplicate)"
            n.update({"verdict": verdict, "action": action,
                       "new_id": rsid if i == 0 else "EXCLUDED"})
    else:
        # Multi-allelic or UNKNOWN_REF: rename ALL probes to CHR:POS:REF:ALT.
        # If two probes at the same position normalise to IDENTICAL REF:ALT
        # (e.g. the same variant submitted twice), keep only the first and
        # exclude the rest — do NOT add a _1/_2 suffix which would imply
        # different variants.
        verdict = "MULTI_ALLELIC" if ref else "UNKNOWN_REF"
        seen_norm: dict[str, int] = {}   # normalized_id → count of collisions seen
        for n in norm_results:
            if ref and n["alt"] != "?":
                base_id = f"{chrom}:{bp}:{ref}:{n['alt']}"
            else:
                base_id = f"{chrom}:{bp}:{n['a1']}:{n['a2']}"
            if base_id in seen_norm:
                # Identical REF:ALT collision — exclude this probe
                seen_norm[base_id] += 1
                dup_id = f"{base_id}__collapse{seen_norm[base_id]}"
                n.update({"verdict": verdict,
                          "action": f"EXCLUDE (collapses to identical {base_id})",
                          "new_id": dup_id})
                bim_rename[n["bim_index"]] = dup_id   # make unique so --exclude works
                collapsed_identical.append(dup_id)
            else:
                seen_norm[base_id] = 0
                n.update({"verdict": verdict,
                          "action": f"RENAME → {base_id}",
                          "new_id": base_id})
                bim_rename[n["bim_index"]] = base_id

    report_rows.extend(norm_results)

# ── Write report ───────────────────────────────────────────────────────────────
report_df = pd.DataFrame(report_rows,
                         columns=["bim_index", "rsid", "chrom", "bp", "a1", "a2",
                                  "ref", "strand", "alt", "new_id", "verdict", "action"])
report_df.to_csv(REPORT, sep="\t", index=False)
print(f"\nReport written: {REPORT.name}  ({len(report_df):,} rows)")

n_palindromic   = (report_df["verdict"] == "PALINDROMIC").sum()
n_multiallelic  = (report_df["verdict"].isin(["MULTI_ALLELIC", "UNKNOWN_REF"])).sum()
print(f"  PALINDROMIC rows   : {n_palindromic:,}")
print(f"  MULTI-ALLELIC rows : {n_multiallelic:,}")

# ── Write exclude list (palindromic redundant rows) ───────────────────────────
# For palindromic dups: all rows share the same rsID, so we can't exclude
# just the extra rows by name. Instead:
# (a) rename ALL rows of each palindromic rsID to unique positional IDs in the BIM
# (b) then exclude everything except the first occurrence
# Simpler: treat palindromic dups the same as multi-allelic — rename each in BIM
# The first gets the "clean" positional ID, extras get _1, _2 suffixes
# Then output exclude list of the suffixed variants.

exclude_ids_final = []
for rsid, group in dup_bim.groupby("SNP", sort=False):
    rows_in_group = report_df[report_df["rsid"] == rsid]
    verdict = rows_in_group.iloc[0]["verdict"]
    if verdict == "PALINDROMIC":
        # Forward-strand probe keeps its rsID — no rename needed.
        # Reverse-strand probe gets a unique positional ID so it can be
        # individually targeted by plink --exclude without touching the keeper.
        fwd_rows = rows_in_group[rows_in_group["strand"] == "+"]
        keep_bim_idx = (
            int(fwd_rows.iloc[0]["bim_index"]) if not fwd_rows.empty
            else int(rows_in_group.iloc[0]["bim_index"])
        )
        suffix_ctr = 0
        for _, row in rows_in_group.iterrows():
            is_keeper = (int(row["bim_index"]) == keep_bim_idx)
            if is_keeper:
                pass   # forward-strand probe: rsID stays, no entry in bim_rename
            else:
                # Reverse-strand probe: give a unique positional name so it
                # can be excluded without also removing the keeper.
                suffix_ctr += 1
                base_id = f"{row['chrom']}:{row['bp']}:{row['ref']}:{row['alt']}"
                if "?" in base_id:
                    base_id = f"{row['chrom']}:{row['bp']}:{row['a1']}:{row['a2']}"
                new_id = f"{base_id}_rev{suffix_ctr}"
                bim_rename[int(row["bim_index"])] = new_id
                exclude_ids_final.append(new_id)


excl_df = pd.DataFrame({"id": exclude_ids_final})
excl_df.to_csv(EXCL_OUT, index=False, header=False)
print(f"\nExclude list written: {EXCL_OUT.name}  ({len(exclude_ids_final):,} IDs)")

# ── Write collapsed identical list ────────────────────────────────────────────
if collapsed_identical:
    pd.DataFrame({"id": collapsed_identical}).to_csv(
        COLLAPSED_OUT, index=False, header=False)
    print(f"Collapsed identical list: {COLLAPSED_OUT.name}  ({len(collapsed_identical):,} IDs)")
else:
    print("No collapsed-identical duplicates found.")

# ── Write multi-allelic table (for liftover) ──────────────────────────────────
multi_rows = report_df[
    report_df["verdict"].isin(["MULTI_ALLELIC", "UNKNOWN_REF"]) &
    ~report_df["action"].str.startswith("EXCLUDE")
].copy()
multi_rows = multi_rows[["rsid", "chrom", "bp", "ref", "alt", "new_id", "strand", "action"]].copy()
multi_rows.columns = ["old_rsid", "chrom", "bp_grch37", "ref_grch37",
                       "alt", "new_id_grch37", "strand", "action"]
multi_rows.to_csv(MULTI_OUT, sep="\t", index=False)
print(f"Multi-allelic table:     {MULTI_OUT.name}  ({len(multi_rows):,} kept variants)")
print(f"  Columns: old_rsid | chrom | bp_grch37 | ref_grch37 | alt | new_id_grch37")
print(f"  Use this table after liftover to rename CHR:POS(hg19) → CHR:POS(hg38).")

# ── Write deduplicated BIM ────────────────────────────────────────────────────
bim_out = bim.copy()
for idx, new_id in bim_rename.items():
    bim_out.at[idx, "SNP"] = new_id
bim_out.to_csv(BIM_OUT, sep="\t", header=False, index=False)
print(f"Deduped BIM written: {BIM_OUT.name}")

# ── Summary ────────────────────────────────────────────────────────────────────
n_palindromic  = (report_df["verdict"] == "PALINDROMIC").sum()
n_multiallelic = report_df["verdict"].isin(["MULTI_ALLELIC", "UNKNOWN_REF"]).sum()
n_collapsed    = len(collapsed_identical)
n_excluded     = len(exclude_ids_final) + n_collapsed
print(f"\n{'═'*58}")
print(f"  Pre-existing duplicate rsIDs  : {len(dup_ids):>6,}")
print(f"  Rows affected                 : {len(dup_bim):>6,}")
print(f"    PALINDROMIC rows            : {n_palindromic:>6,}")
print(f"    MULTI-ALLELIC rows          : {n_multiallelic:>6,}")
print(f"  Rows renamed in BIM           : {len(bim_rename):>6,}")
print(f"  Rows to exclude (palindromic) : {len(exclude_ids_final):>6,}")
print(f"  Rows collapsed (identical)    : {n_collapsed:>6,}")
print(f"{'═'*58}")
print(f"\nNext steps:")
print(f"  1. Copy .bed and .fam alongside the deduped BIM:")
print(f'     cp "{BIM_IN.with_suffix(".bed")}" "{BIM_OUT.with_suffix(".bed")}"')
print(f'     cp "{BIM_IN.with_suffix(".fam")}" "{BIM_OUT.with_suffix(".fam")}"')
print(f"  2. Exclude reverse-strand and collapsed probes:")
all_excl = str(EXCL_OUT)
print(f'     plink --bfile "{BIM_OUT.with_suffix("")}" \\')
if collapsed_identical:
    print(f'           --exclude "{all_excl}" \\'   # palindromic
          f'  # also manually append {COLLAPSED_OUT.name} if needed')
else:
    print(f'           --exclude "{EXCL_OUT}" \\')
print(f'           --make-bed \\')
print(f'           --out "D:/ADNI_SNP_Omni2.5M_20140220/SNP_filtered_with_mri_rsid_clean"')
print(f"  3. Run validate_rsid.py (update BIM_FILE path to _clean)")
print(f"  4. Apply patch_deprecated_safe.txt")
