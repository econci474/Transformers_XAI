"""
validate_rsid.py
Check whether assigned rsIDs are still current in dbSNP using the NCBI
RsMergeArch bulk file — entirely offline after a one-time download.

WHY NOT USE AN API:
  1.57M rsIDs x API round-trip = 45-180 minutes.
  The RsMergeArch file encodes every rsID merge since dbSNP inception
  and can be queried locally in seconds.

DOWNLOAD (one-time, ~100 MB):
  https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/database/organism_data/RsMergeArch.bcp.gz
  Save to MERGE_FILE path below.  No decompression needed.

RsMergeArch column layout (tab-separated, no header):
  0  rsHigh      deprecated rsID (the one that was retired)
  1  rsLow       surviving rsID at time of merge
  2  build_id
  3  orien
  4  create_time
  5  last_updated_time
  6  rsCurrent   FINAL current rsID after all chained merges  <-- we use this
  7  snp3P
  8  comment

Outputs:
  deprecated_rsids.tsv         full table: kgp_id | old_rsid | current_rsid | status
  patch_deprecated_safe.txt    PLINK --update-name: old_rsid  current_rsid
                               (conflicts where target already exists in BIM are excluded)
  patch_deprecated_conflicts.txt  pairs that would create duplicate IDs
                                  (these variants are effectively the same SNP in dbSNP)

Apply patch:
  plink --bfile SNP_filtered_with_mri_rsid \\
        --update-name patch_deprecated_safe.txt \\
        --make-bed --out SNP_filtered_with_mri_rsid_current
"""

import argparse
import gzip
import io
import pathlib
import sys
import time

# Force UTF-8 output on Windows (avoids cp1252 UnicodeEncodeError)
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import pandas as pd

# ── CLI (all args optional — defaults preserve original behaviour) ──────────────
_p = argparse.ArgumentParser(
    description="Validate and patch deprecated rsIDs in a PLINK BIM.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
_p.add_argument("--bim",    default="D:/ADNI_SNP_Omni2.5M_20140220/SNP_filtered_with_mri_rsid_clean.bim",
                help="BIM file to check (default: original pipeline path)")
_p.add_argument("--outdir", default="D:/ADNI_SNP_Omni2.5M_20140220",
                help="Directory for output files")
_p.add_argument("--prefix", default="",
                help="Prefix for output file names, e.g. 'ld_pruned_07_'")
_args = _p.parse_args()

# ── Paths ──────────────────────────────────────────────────────────────────────
MERGE_FILE     = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220/liftover/NCBI/RsMergeArch.bcp.gz")
LOOKUP_TSV     = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220/kgp_to_rs_lookup.tsv")
BIM_FILE       = pathlib.Path(_args.bim)
_outdir        = pathlib.Path(_args.outdir)
_pfx           = _args.prefix
DEPRECATED_TSV = _outdir / f"{_pfx}deprecated_rsids.tsv"
PATCH_SAFE     = _outdir / f"{_pfx}patch_deprecated_safe.txt"
PATCH_CONFLICT = _outdir / f"{_pfx}patch_deprecated_conflicts.txt"

# ── Check files exist ──────────────────────────────────────────────────────────
for p in [MERGE_FILE, LOOKUP_TSV, BIM_FILE]:
    if not p.exists():
        print(f"[ERROR] File not found: {p}")
        if p == MERGE_FILE:
            print("  Download from:")
            print("  https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/database/"
                  "organism_data/RsMergeArch.bcp.gz")
        sys.exit(1)

# ── Step 1: Load RsMergeArch → build deprecated_rs → current_rs dict ──────────
print(f"[1/4] Loading RsMergeArch merge table: {MERGE_FILE.name} ...")
t0 = time.time()

deprecated_to_current = {}   # {rs_high_int: rs_current_int}

opener = gzip.open if MERGE_FILE.suffix == ".gz" else open
with opener(MERGE_FILE, "rt", encoding="utf-8", errors="replace") as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 7:
            continue
        try:
            rs_high    = int(parts[0])
            rs_low     = int(parts[1])
            rs_current = int(parts[6]) if parts[6].strip() not in ("", "0", "\\N") else rs_low
        except (ValueError, IndexError):
            continue
        deprecated_to_current[rs_high] = rs_current

print(f"     Loaded {len(deprecated_to_current):,} merge records  ({time.time()-t0:.1f}s)")

# ── Step 2: Load BIM — collect existing rsIDs ──────────────────────────────────
print(f"\n[2/4] Loading BIM: {BIM_FILE.name} ...")
bim = pd.read_csv(BIM_FILE, sep="\t", header=None,
                  names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
                  dtype=str, usecols=["SNP"])
# Strip whitespace — Windows line endings can add \r to the last field
existing_ids = set(bim["SNP"].str.strip())
print(f"     {len(existing_ids):,} variant IDs in BIM")
print(f"     rs15842 in BIM: {'YES' if 'rs15842' in existing_ids else 'NO'}")

# ── Step 3: Load lookup → find deprecated rsIDs ───────────────────────────────
print(f"\n[3/4] Loading rsID lookup: {LOOKUP_TSV.name} ...")
lookup  = pd.read_csv(LOOKUP_TSV, sep="\t", dtype=str)
rs_rows = lookup[lookup["new_id"].str.match(r"^rs\d+$", na=False)].copy()
print(f"     {len(rs_rows):,} assigned rsIDs to validate")

print("\n[4/4] Checking against merge table ...")
t1 = time.time()

deprecated = []
for _, row in rs_rows.iterrows():
    try:
        rs_num = int(row["new_id"].lstrip("rs"))
    except ValueError:
        continue
    current_num = deprecated_to_current.get(rs_num)
    if current_num is not None and current_num != rs_num:
        deprecated.append({
            "kgp_id":       row["kgp_id"],
            "old_rsid":     row["new_id"],
            "current_rsid": f"rs{current_num}",
            "status":       "merged",
        })

print(f"     Done in {time.time()-t1:.1f}s")

# ── Build results DataFrame ───────────────────────────────────────────────────
df = pd.DataFrame(deprecated) if deprecated else pd.DataFrame(
    columns=["kgp_id", "old_rsid", "current_rsid", "status"])

merged = df[df["status"] == "merged"].copy()

# ── Conflict detection: simulate all renames against a live ID set ────────────
# Instead of enumerating edge-case types, we walk every rename and check
# whether the target already exists in the running ID set.  This handles:
#   - target already in BIM
#   - two patch entries mapping to the same target
#   - chain merges (rs_X→rs15842 AND rs15842→rs_Y in same pass)
# in one unified pass.

live_ids = set(existing_ids)   # all IDs currently in the BIM

safe_rows      = []
conflict_rows  = []

for _, row in merged.iterrows():
    old_rsid = row["old_rsid"].strip()
    new_rsid = row["current_rsid"].strip()

    if new_rsid in live_ids:
        conflict_rows.append(row)
        if new_rsid == "rs15842":
            print(f"  [DEBUG] CONFLICT: {old_rsid} → {new_rsid}  (target already in live_ids)")
    else:
        safe_rows.append(row)
        live_ids.discard(old_rsid)
        live_ids.add(new_rsid)
        if new_rsid == "rs15842" or old_rsid == "rs15842":
            print(f"  [DEBUG] SAFE:     {old_rsid} → {new_rsid}")

safe      = pd.DataFrame(safe_rows)     if safe_rows     else pd.DataFrame(columns=merged.columns)
conflicts = pd.DataFrame(conflict_rows) if conflict_rows else pd.DataFrame(columns=merged.columns)



# ── Outputs ───────────────────────────────────────────────────────────────────
df.to_csv(DEPRECATED_TSV, sep="\t", index=False)
safe[["old_rsid", "current_rsid"]].to_csv(PATCH_SAFE, sep="\t", header=False, index=False)
conflicts[["old_rsid", "current_rsid"]].to_csv(PATCH_CONFLICT, sep="\t", header=False, index=False)

total = len(rs_rows)
print(f"\n── Summary {'─'*40}")
print(f"  rsIDs checked          : {total:>10,}")
print(f"  Deprecated (merged)    : {len(merged):>10,}  ({100*len(merged)/total:.2f}%)")
print(f"    Safe to rename       : {len(safe):>10,}  → patch_deprecated_safe.txt")
print(f"    Conflict (dup IDs)   : {len(conflicts):>10,}  → patch_deprecated_conflicts.txt")
print(f"  Still current          : {total-len(merged):>10,}  ({100*(total-len(merged))/total:.2f}%)")
print(f"{'─'*50}")

if len(conflicts) > 0:
    print(f"\n  NOTE: {len(conflicts):,} renames skipped because the target rsID already exists")
    print(f"  in the BIM as a separate probe. These represent variants that dbSNP has")
    print(f"  merged — both probes measure the same underlying SNP. You may want to")
    print(f"  remove one copy with PLINK --exclude after the rename step.")
    print(f"  See: {PATCH_CONFLICT.name}")

print(f"\nOutputs written to {DEPRECATED_TSV.parent}")
print(f"  {DEPRECATED_TSV.name}")
print(f"  {PATCH_SAFE.name}")
print(f"  {PATCH_CONFLICT.name}")

if len(safe) > 0:
    print(f"\nNext step — apply safe patch:")
    print(f"  plink --bfile {BIM_FILE.with_suffix('').name} \\")
    print(f"        --update-name {PATCH_SAFE} \\")
    print(f"        --make-bed \\")
    print(f"        --out {BIM_FILE.parent / 'SNP_filtered_with_mri_rsid_current'}")
