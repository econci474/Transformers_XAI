"""
inspect_duplicate_rsid.py
═════════════════════════
Inspect the conflict pairs produced by validate_rsid.py to determine whether
each duplicate rsID pair represents:

  (a) Exact duplicate  — same CHR, BP, allele set (safe to drop one)
  (b) Allele mismatch  — same position, different alleles (tri-allelic / strand flip)
  (c) Position clash   — truly different positions sharing an rsID (unusual)

Then runs PLINK2 --rm-dup to report duplicate-position variants in the BIM.

Usage:
    python inspect_duplicate_rsid.py

Outputs printed to terminal + written to inspect_duplicates_report.txt
"""

import pathlib
import subprocess
import sys

import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────────────
BIM_FILE      = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220/SNP_filtered_with_mri_rsid_clean_current.bim")
CONFLICT_FILE = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220/patch_deprecated_conflicts.txt")
BIM_PREFIX    = "D:/ADNI_SNP_Omni2.5M_20140220/SNP_filtered_with_mri_rsid_clean_current"
REPORT_FILE   = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220/inspect_duplicates_report.txt")

for p in [BIM_FILE, CONFLICT_FILE]:
    if not p.exists():
        print(f"[ERROR] File not found: {p}")
        sys.exit(1)

# ── Load BIM ──────────────────────────────────────────────────────────────────
print(f"Loading BIM: {BIM_FILE.name} ...")
bim = pd.read_csv(BIM_FILE, sep="\t", header=None,
                  names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
                  dtype=str)
bim["BP"] = bim["BP"].str.strip()
bim["SNP"] = bim["SNP"].str.strip()
print(f"  {len(bim):,} variants loaded")

# ── Load conflict pairs ───────────────────────────────────────────────────────
print(f"Loading conflict pairs: {CONFLICT_FILE.name} ...")
conflicts = pd.read_csv(CONFLICT_FILE, sep="\t", header=None,
                        names=["old_rsid", "current_rsid"], dtype=str)
conflicts["old_rsid"]     = conflicts["old_rsid"].str.strip()
conflicts["current_rsid"] = conflicts["current_rsid"].str.strip()
print(f"  {len(conflicts):,} pairs\n")

# ── Inspect each pair ─────────────────────────────────────────────────────────
lines = []
exact_dups   = []
mismatches   = []
not_in_bim   = []

sep = "─" * 72

for _, row in conflicts.iterrows():
    old_id = row["old_rsid"]
    new_id = row["current_rsid"]

    old_rows = bim[bim["SNP"] == old_id]
    new_rows = bim[bim["SNP"] == new_id]

    lines.append(f"\n{sep}")
    lines.append(f"Conflict:  {old_id}  →  {new_id}  (old → current rsID)")
    lines.append(f"{sep}")

    if old_rows.empty:
        lines.append(f"  {old_id} : NOT FOUND in BIM (already renamed?)")
    else:
        for _, r in old_rows.iterrows():
            lines.append(f"  {old_id:20s}  CHR={r.CHR:>2}  BP={r.BP:>12}  A1={r.A1}  A2={r.A2}")

    if new_rows.empty:
        lines.append(f"  {new_id} : NOT FOUND in BIM")
    else:
        for _, r in new_rows.iterrows():
            lines.append(f"  {new_id:20s}  CHR={r.CHR:>2}  BP={r.BP:>12}  A1={r.A1}  A2={r.A2}")

    # Classify
    if old_rows.empty or new_rows.empty:
        verdict = "NOT IN BIM — skip"
        not_in_bim.append({"old": old_id, "new": new_id})
    elif len(old_rows) > 1 or len(new_rows) > 1:
        verdict = "MULTI-ROW — complex (investigate manually)"
    else:
        o = old_rows.iloc[0]
        n = new_rows.iloc[0]
        same_pos    = (o["CHR"].strip() == n["CHR"].strip() and
                       o["BP"].strip()  == n["BP"].strip())
        same_alleles = (frozenset({o["A1"].strip(), o["A2"].strip()} - {"0","."}) ==
                        frozenset({n["A1"].strip(), n["A2"].strip()} - {"0","."}))

        if same_pos and same_alleles:
            verdict = "EXACT DUPLICATE — same position & alleles (safe to drop one)"
            exact_dups.append({"old": old_id, "new": new_id})
        elif same_pos and not same_alleles:
            verdict = "ALLELE MISMATCH — same position, different alleles (tri-allelic?)"
            mismatches.append({"old": old_id, "new": new_id, "reason": "allele_mismatch"})
        else:
            verdict = "POSITION CLASH — different positions (unusual — investigate)"
            mismatches.append({"old": old_id, "new": new_id, "reason": "position_clash"})

    lines.append(f"\n  ► Verdict: {verdict}")

lines.append(f"\n{'═'*72}")
lines.append(f"SUMMARY")
lines.append(f"{'═'*72}")
lines.append(f"  Total conflict pairs         : {len(conflicts):>4}")
lines.append(f"  Exact duplicates (drop one)  : {len(exact_dups):>4}")
lines.append(f"  Allele / position mismatches : {len(mismatches):>4}")
lines.append(f"  Old rsID not in BIM          : {len(not_in_bim):>4}")
lines.append(f"{'═'*72}")

if exact_dups:
    lines.append(f"\nFor exact duplicates, PLINK2 --rm-dup is the cleanest fix.")
    lines.append(f"See PLINK2 section below.")
if mismatches:
    lines.append(f"\nWARNING: {len(mismatches)} mismatch(es) need manual review.")
    lines.append(f"These may be tri-allelic sites or mis-annotated rsIDs.")

report = "\n".join(lines)
print(report)
REPORT_FILE.write_text(report, encoding="utf-8")
print(f"\nReport written to: {REPORT_FILE}")

# ── PLINK2 --rm-dup ───────────────────────────────────────────────────────────
print(f"\n{'═'*72}")
print(f"PLINK2  --rm-dup  check")
print(f"{'═'*72}")
print(f"Running: plink2 --bfile {BIM_PREFIX} --rm-dup error")
print(f"(If duplicates exist, this will list them and exit without writing files.)\n")

cmd = ["plink2", "--bfile", BIM_PREFIX, "--rm-dup", "error",
       "--make-bed", "--out", BIM_PREFIX + "_rmdup_check"]
result = subprocess.run(cmd, capture_output=True, text=True)

output = result.stdout + result.stderr
for line in output.splitlines():
    print(f"  {line}")

if result.returncode == 0:
    print("\n  ✓ PLINK2 found no duplicate-position variants.")
else:
    print(f"\n  Duplicates detected (see above).")
    print(f"\n  To REMOVE exact duplicates and keep one copy:")
    print(f'    plink2 --bfile "{BIM_PREFIX}" \\')
    print(f'           --rm-dup force-first \\')
    print(f'           --make-bed \\')
    print(f'           --out "D:/ADNI_SNP_Omni2.5M_20140220/SNP_filtered_with_mri_rsid_rmdup"')
    print(f"\n  To REMOVE mismatched duplicates and keep exact ones:")
    print(f'    plink2 --bfile "{BIM_PREFIX}" \\')
    print(f'           --rm-dup exclude-mismatch \\')
    print(f'           --make-bed \\')
    print(f'           --out "D:/ADNI_SNP_Omni2.5M_20140220/SNP_filtered_with_mri_rsid_rmdup"')
