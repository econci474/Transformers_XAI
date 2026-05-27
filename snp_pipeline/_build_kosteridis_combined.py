"""
_build_kosteridis_combined.py
=============================
Union the two Kosteridis sub-bundles into a single resolution TSV:
  Kosteridis_MTAG_AD       — 114 published lead SNPs (Supp Tab 1), MTAG-AD β.
  Kosteridis_shared_AD_CV  — 33 colocalising AD–CV loci (Supp Tab 3 ∪ Tb 4).
                              17 in Supp1 (MTAG β); 16 sumstats_fallback β.

These come from the same Kosteridis 2024 GWAS (PMID 38826624,
GCST90449060). Treating them as separate "sources" in the leaderboard
double-counts overlapping SNPs and obscures the true Kosteridis signal.

This script writes:
  D:/ADNI_SNP_Omni2.5M_20140220/source_prs/Kosteridis_full_snps_resolution.tsv

Deduplication rule: union by rsID_canonical. When a SNP appears in both
sub-bundles, prefer MTAG_AD β (Supp1 source — that's the canonical
MTAG-AD-boosted β; shared_AD_CV agrees with it for those rsIDs anyway).
For shared_AD_CV-only rsIDs (16 sumstats_fallback entries), keep the
sumstats β from the full GCST90449060 file.

`Kosteridis_novel_AD` (5 SNPs, publication Table 1) is a strict subset
of MTAG_AD by construction (all 5 are in Supp1), so it adds nothing to
the union — we exclude it from the merge.

Usage (env: snp):
  python snp_pipeline/_build_kosteridis_combined.py
  python snp_pipeline/_build_kosteridis_combined.py --dry-run
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

SRC_DIR = Path("D:/ADNI_SNP_Omni2.5M_20140220/source_prs")
OUT_FILE = SRC_DIR / "Kosteridis_full_snps_resolution.tsv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    mtag = pd.read_csv(SRC_DIR / "Kosteridis_MTAG_AD_full_snps_resolution.tsv",
                        sep="\t")
    sh   = pd.read_csv(SRC_DIR / "Kosteridis_shared_AD_CV_full_snps_resolution.tsv",
                        sep="\t")

    print(f"MTAG_AD       : {len(mtag):3d} rows ({(mtag.drop_reason=='ok').sum()} resolved)")
    print(f"shared_AD_CV  : {len(sh):3d} rows ({(sh.drop_reason=='ok').sum()} resolved)")

    # Tag the source-of-β for traceability
    mtag = mtag.assign(_kost_origin="MTAG_AD_Supp1")
    sh   = sh.assign(_kost_origin=sh.get("evidence_source", "shared_AD_CV"))

    # Dedupe priority: MTAG_AD first → for overlapping rsIDs, MTAG_AD wins.
    combined = pd.concat([mtag, sh], ignore_index=True)
    n_before = len(combined)
    combined = combined.drop_duplicates(subset=["rsID_canonical"], keep="first")
    n_after = len(combined)
    print(f"\nUnion before dedup: {n_before}  ->  after: {n_after}  "
          f"(removed {n_before - n_after} overlap rows)")

    # Update the source column so the lib's bookkeeping reads "Kosteridis".
    combined["source"] = "Kosteridis"

    ok = combined[combined["drop_reason"] == "ok"]
    print(f"\nResolved-OK rows in combined: {len(ok)}  "
          f"(with rsID_in_bim: {ok['rsID_in_bim'].notna().sum()})")
    print(f"_kost_origin counts:\n{combined['_kost_origin'].value_counts(dropna=False)}")

    if args.dry_run:
        print(f"\n[DRY] would write {OUT_FILE}")
        return
    combined.to_csv(OUT_FILE, sep="\t", index=False)
    print(f"\nwrote {OUT_FILE}  ({len(combined)} rows)")


if __name__ == "__main__":
    main()
