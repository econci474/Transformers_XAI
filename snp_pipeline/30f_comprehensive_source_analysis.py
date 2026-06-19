r"""
30f_comprehensive_source_analysis.py   (LOCAL, env: snp)
========================================================
After extending the resolution pipeline (30b + 30c) from 4 to 8 GWAS sources
(adding Lambert 2013, DeRojas 2021, Schwanzentruber 2021, Najar 2021),
this script produces the comprehensive cross-source delta analysis:

  - Which SNPs do we already have from the existing 4 sources (B/W/K/D)?
  - Which SNPs are NEW from the 4 added sources?
  - What's the new total number of unique SNPs available from our array?
  - Per-SNP, what's the reason for absence (not genotyped vs filtered out)?

Output → `D:\ADNI_SNP_Omni2.5M_20140220\GWAS_comprehensive\`:
  existing_vs_new_summary.txt   compact text summary (terminal-friendly)
  all_snps_master.tsv           one row per unique published rsID across
                                 all 8 sources; columns: rsID,
                                 in_<source>×8 booleans, on_MAF,
                                 in_unfilt_BIM, filter_reason, gene
  new_snps_from_extension.tsv   SNPs added by Lambert/DeRojas/Schwanzentruber/
                                 Najar that were NOT in B/W/K/D
  delta_on_MAF.tsv              on-MAF SNPs newly available vs the prior set
  delta_recoverable.tsv         on-chip-but-filtered SNPs newly added

Usage:
  conda run -n snp python snp_pipeline/30f_comprehensive_source_analysis.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
RECON = BASE / "source_prs" / "unfiltered_SNP_reconciliation"
# v2 output folder — parallel to v1 (8 sources). The v1 folder
# `GWAS_comprehensive/` remains intact for thesis-reference traceability.
OUT = BASE / "GWAS_comprehensive_v2"

# Pre-extension 4 sources (the original baseline before 2026-05-23).
EXISTING_SOURCES = ("Bellenguez", "Wightman", "Kunkle", "Desikan")
# 2026-05-23 4-source extension.
NEW_SOURCES_V1 = ("Lambert", "DeRojas", "Schwanzentruber", "Najar")
# 2026-05-25 8-source extension.
NEW_SOURCES_V2 = ("Ebenau", "Leonenko", "Vesilievick", "Zhang",
                  "Felsky_MF", "Felsky_IT", "ONeil_NPY", "ONeil_GHR",
                  "Kosteridis_novel_AD", "Kosteridis_shared_AD_CV")
# 2026-05-25 night — Yang 2023 cell-type-partitioned PRS sub-sources.
NEW_SOURCES_V3 = ("Yang_Ex", "Yang_In", "Yang_Ast",
                  "Yang_Mic", "Yang_Oli", "Yang_Opc")
# Audit-only (not in PRS roster but tracked for cross-source overlap).
AUDIT_SOURCES = ("Huang", "ONeil_SST_candidates")
# All PRS sources (22) + audit-only (2) = 24.
PRS_SOURCES = EXISTING_SOURCES + NEW_SOURCES_V1 + NEW_SOURCES_V2 + NEW_SOURCES_V3
ALL_SOURCES = PRS_SOURCES + AUDIT_SOURCES

# Maps pheno_class → coarse binary `pheno_type` for clinical-vs-pathology
# stratification.
PHENO_TYPE = {
    "clinical_AD":                 "clinical",
    "clinical_AD_MTAG":            "clinical",
    "clinical_AAO":                "clinical",
    "pathology_microglia":         "pathology",
    "pathology_amyloid":           "pathology",
    "pathology_plaque":            "pathology",
    "pathology_mixed":             "pathology",
    "pathology_astrocyte":         "pathology",
    "pathology_neuron_excitatory": "pathology",
    "pathology_neuron_inhibitory": "pathology",
    "pathology_oligodendrocyte":   "pathology",
    "pathology_opc":               "pathology",
}


def _load_reconciliation(src: str) -> pd.DataFrame:
    p = RECON / f"{src}_unfiltered_reconciliation.tsv"
    df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    df["source"] = src
    return df


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    # ── Load all 18 source reconciliations (16 PRS + 2 audit) ────────────
    per_src = {}
    for src in ALL_SOURCES:
        try:
            per_src[src] = _load_reconciliation(src)
        except FileNotFoundError as e:
            print(f"  [{src:24s}] MISSING: {e}")
            continue
    for src, df in per_src.items():
        n = len(df)
        audit_tag = " [AUDIT-ONLY]" if src in AUDIT_SOURCES else ""
        pheno = df["pheno_class"].iloc[0] if (
            len(df) and "pheno_class" in df.columns) else "?"
        print(f"  [{src:24s}] {n:4d} rows  pheno_class={pheno}{audit_tag}")

    # ── Build per-rsID master across all sources ──────────────────────────
    # Key by rsID_pub (publication rsID); preserve canonical for cross-ref.
    print("\n[master] building per-rsID master table …")
    rows = []
    for src, df in per_src.items():
        is_audit = src in AUDIT_SOURCES
        for _, r in df.iterrows():
            rs = str(r["rsID_pub"])
            if not rs.startswith("rs"):
                continue
            ph = str(r.get("pheno_class", ""))
            rows.append({
                "rsID": rs,
                "source": src,
                "pheno_class": ph,
                "pheno_type": PHENO_TYPE.get(ph, ""),
                "audit_only": is_audit,
                "rsID_canonical": str(r.get("rsID_canonical", rs)),
                "in_maf_bim": str(r.get("in_maf_bim", "False")) == "True",
                "in_unfilt_bim": str(r.get("in_unfiltered_bim", "False")) == "True",
                "reconciliation_category": str(r.get("reconciliation_category", "")),
                "filter_reason": str(r.get("filter_reason", "")),
                "drop_stage": str(r.get("drop_stage", "")),
                "unfilt_kgp_name": str(r.get("unfilt_kgp_name", "")),
                "kgp_manifest_source": str(r.get("kgp_manifest_source", "")),
                "unfilt_CHR_grch37": str(r.get("unfilt_CHR_grch37", "")),
                "unfilt_BP_grch37": str(r.get("unfilt_BP_grch37", "")),
                "unfilt_A1": str(r.get("unfilt_A1", "")),
                "unfilt_A2": str(r.get("unfilt_A2", "")),
                "F_MISS_812": str(r.get("F_MISS_812", "")),
                "MAF_812": str(r.get("MAF_812", "")),
                "HWE_p_812": str(r.get("HWE_p_812", "")),
                "gene": str(r.get("locus_name", "")),
                "beta_pub": str(r.get("OR_or_beta_pub", "")),
                "se_pub": str(r.get("SE_pub", "")),
            })
    long_df = pd.DataFrame(rows)
    print(f"  long-form rows (rsID × source): {len(long_df)}")

    # Collapse to per-rsID master: union flags + chosen filter_reason
    master_rows = []
    for rs, grp in long_df.groupby("rsID"):
        sources_with_rs = set(grp["source"])
        in_existing = bool(sources_with_rs & set(EXISTING_SOURCES))
        in_v1 = bool(sources_with_rs & set(NEW_SOURCES_V1))
        in_v2 = bool(sources_with_rs & set(NEW_SOURCES_V2))
        in_v3 = bool(sources_with_rs & set(NEW_SOURCES_V3))
        in_audit = bool(sources_with_rs & set(AUDIT_SOURCES))
        is_prs = bool(sources_with_rs & set(PRS_SOURCES))
        on_MAF = grp["in_maf_bim"].any()
        in_unfilt = grp["in_unfilt_bim"].any()
        pheno_classes = sorted({p for p in grp["pheno_class"] if p})
        pheno_types = sorted({p for p in grp["pheno_type"] if p})
        # Pick the row with the most informative status
        # (on_MAF beats filtered_out beats not_on_chip)
        priority = {"on_chip_passed_maf": 0,
                    "on_chip_filtered_out": 1,
                    "on_chip_via_rsID_merge_filtered_out": 1,
                    "on_chip_via_kgp_renamed_filtered_out": 1,
                    "on_chip_via_position_filtered_out": 1,
                    "in_maf_but_not_on_chip_manifest": 2,
                    "not_on_chip": 3,
                    "apoe_haplotype_special": 4}
        grp = grp.copy()
        grp["_prio"] = grp["reconciliation_category"].map(
            lambda c: priority.get(c, 5))
        best = grp.sort_values("_prio").iloc[0]
        master_rows.append({
            "rsID": rs,
            "n_sources": len(sources_with_rs),
            "in_existing_4": in_existing,
            "in_new_v1": in_v1,
            "in_new_v2": in_v2,
            "in_new_v3": in_v3,
            "in_audit": in_audit,
            "is_prs": is_prs,
            "membership": "+".join(sorted(sources_with_rs)),
            "pheno_classes": "|".join(pheno_classes),
            "pheno_types": "|".join(pheno_types),
            "on_MAF": on_MAF,
            "in_unfilt_bim": in_unfilt,
            "reconciliation_category": best["reconciliation_category"],
            "filter_reason": best["filter_reason"],
            "drop_stage": best["drop_stage"],
            "unfilt_kgp_name": best["unfilt_kgp_name"],
            "kgp_manifest_source": best["kgp_manifest_source"],
            "unfilt_CHR_grch37": best["unfilt_CHR_grch37"],
            "unfilt_BP_grch37": best["unfilt_BP_grch37"],
            "unfilt_A1": best["unfilt_A1"],
            "unfilt_A2": best["unfilt_A2"],
            "F_MISS_812": best["F_MISS_812"],
            "MAF_812": best["MAF_812"],
            "HWE_p_812": best["HWE_p_812"],
            "gene": best["gene"],
        })
    master = (pd.DataFrame(master_rows)
              .sort_values(["on_MAF", "n_sources", "rsID"],
                            ascending=[False, False, True]))
    master.to_csv(args.out / "all_snps_master.tsv", sep="\t", index=False)
    print(f"  [out] all_snps_master.tsv  ({len(master)} unique rsIDs)")

    # ── Sub-tables: new SNPs from v2+v3 extension, delta on-MAF, delta recover
    # Baseline v1 = existing(4) + v1(4) = 8 sources. Post-baseline extensions:
    # v2 (10 sources, 2026-05-25 evening) + v3 (6 Yang sources, 2026-05-25 night).
    in_baseline_v1 = master["in_existing_4"] | master["in_new_v1"]
    in_post_baseline = master["in_new_v2"] | master["in_new_v3"]
    new_post = master[(~in_baseline_v1) & in_post_baseline].copy()
    new_post.to_csv(args.out / "new_snps_from_v2_extension.tsv",
                       sep="\t", index=False)
    print(f"  [out] new_snps_from_v2_extension.tsv  ({len(new_post)} rsIDs from v2+v3)")

    delta_maf = new_post[new_post["on_MAF"]].copy()
    delta_maf.to_csv(args.out / "delta_on_MAF.tsv", sep="\t", index=False)
    print(f"  [out] delta_on_MAF.tsv  ({len(delta_maf)} NEW on-MAF SNPs added by v2+v3 extension)")

    delta_recov = new_post[
        new_post["reconciliation_category"].isin([
            "on_chip_filtered_out",
            "on_chip_via_rsID_merge_filtered_out",
            "on_chip_via_kgp_renamed_filtered_out",
            "on_chip_via_position_filtered_out"])].copy()
    delta_recov.to_csv(args.out / "delta_recoverable.tsv",
                       sep="\t", index=False)
    print(f"  [out] delta_recoverable.tsv  ({len(delta_recov)} NEW on-chip-but-filtered)")

    # ── Per-phenoclass stratum sizes (PRS-usable rsIDs only) ─────────────
    # For rsIDs in multiple phenoclasses, prefer the PRS-grade ones over
    # audit-only (clinical_AAO from Huang, pathology_mixed from ONeil-SST).
    _PRIORITY = ["clinical_AD", "clinical_AD_MTAG",
                 "pathology_microglia",      # Felsky + Yang_Mic
                 "pathology_astrocyte",      # Yang_Ast
                 "pathology_neuron_excitatory",
                 "pathology_neuron_inhibitory",
                 "pathology_oligodendrocyte",
                 "pathology_opc",
                 "pathology_plaque", "pathology_amyloid",
                 "clinical_AAO", "pathology_mixed"]
    def _pick_main(s):
        classes = [c for c in str(s).split("|") if c]
        for p in _PRIORITY:
            if p in classes:
                return p
        return classes[0] if classes else ""
    prs_master = master[master["is_prs"]].copy()
    prs_master["main_class"] = prs_master["pheno_classes"].map(_pick_main)
    by_pc = (prs_master.groupby("main_class")
              .agg(unique_rsIDs=("rsID", "nunique"),
                   on_MAF=("on_MAF", "sum"),
                   on_chip=("in_unfilt_bim", "sum"))
              .reset_index())
    by_pc.to_csv(args.out / "per_phenoclass_summary.tsv",
                 sep="\t", index=False)
    print(f"  [out] per_phenoclass_summary.tsv  ({len(by_pc)} phenoclass strata)")

    # ── Per-source × category counts ──────────────────────────────────────
    cat_tab = (long_df.groupby(["source", "reconciliation_category"])
               .size().unstack(fill_value=0))
    cat_tab.to_csv(args.out / "per_source_category_counts.tsv", sep="\t")
    print(f"  [out] per_source_category_counts.tsv")

    # ── Text summary ─────────────────────────────────────────────────────
    n_master = len(master)
    n_master_prs = master["is_prs"].sum()
    n_master_audit_only = ((master["in_audit"]) & ~master["is_prs"]).sum()

    in_baseline_v1 = master["in_existing_4"] | master["in_new_v1"]
    # v1 baseline = 8 sources; v2 = +8 more = 16. We track the v2 delta.
    on_maf_pre = master[in_baseline_v1 & master["on_MAF"]].shape[0]
    on_maf_post = master[master["is_prs"] & master["on_MAF"]].shape[0]

    n_added_to_maf = delta_maf.shape[0]
    n_added_recov = delta_recov.shape[0]
    n_v2_unique = ((master["in_new_v2"] | master["in_new_v3"]) & ~in_baseline_v1).sum()

    lines = [
        "Comprehensive PRS-source analysis v2 — pre/post 2026-05-25 delta",
        "",
        f"Baseline sources (pre 2026-05-23, 4):     {', '.join(EXISTING_SOURCES)}",
        f"v1 extension (2026-05-23, +4):            {', '.join(NEW_SOURCES_V1)}",
        f"v2 extension (2026-05-25, +8 PRS + 2 audit):",
        f"  PRS:    {', '.join(NEW_SOURCES_V2)}",
        f"  audit:  {', '.join(AUDIT_SOURCES)}",
        f"",
        f"PRS sources total = 16; including audit-only = 18.",
        "",
        "=" * 70,
        "## Per-source published-lead counts",
        "",
    ]
    for src in ALL_SOURCES:
        if src not in per_src:
            continue
        df = per_src[src]
        n_pub = len(df)
        n_maf = (df["reconciliation_category"] == "on_chip_passed_maf").sum()
        n_filt = df["reconciliation_category"].isin([
            "on_chip_filtered_out",
            "on_chip_via_rsID_merge_filtered_out",
            "on_chip_via_kgp_renamed_filtered_out",
            "on_chip_via_position_filtered_out"]).sum()
        n_nochip = (df["reconciliation_category"] == "not_on_chip").sum()
        n_hap = (df["reconciliation_category"] == "apoe_haplotype_special").sum()
        if src in AUDIT_SOURCES:
            tag = " (AUDIT)"
        elif src in NEW_SOURCES_V3:
            tag = " (NEW-v3)"
        elif src in NEW_SOURCES_V2:
            tag = " (NEW-v2)"
        elif src in NEW_SOURCES_V1:
            tag = " (NEW-v1)"
        else:
            tag = ""
        pheno = df["pheno_class"].iloc[0] if (
            len(df) and "pheno_class" in df.columns) else ""
        lines.append(
            f"  {src:24s}{tag:10s} [{pheno:18s}]: "
            f"pub={n_pub:>4}  on_MAF={n_maf:>4}  "
            f"filtered_out={n_filt:>3}  not_on_chip={n_nochip:>4}  "
            f"hap={n_hap:>1}")

    lines.extend([
        "",
        "=" * 70,
        "## Unique-rsID master across all 18 sources",
        "",
        f"  total unique published rsIDs (all sources)      : {n_master:>5}",
        f"    of which PRS-usable                           : {n_master_prs:>5}",
        f"    of which audit-only (not in any PRS source)   : {n_master_audit_only:>5}",
        f"  rsIDs in v2 sources only (not in v1 baseline)   : {n_v2_unique:>5}",
        "",
        "=" * 70,
        "## On-MAF availability (the SNPs actually usable for PRS)",
        "",
        f"  pre-extension v1 baseline (8 sources) on-MAF    : {on_maf_pre:>5}",
        f"  post-extension v2 (16 PRS sources) on-MAF       : {on_maf_post:>5}",
        f"  NET ADDITION from v2 extension                  : {on_maf_post - on_maf_pre:>5}",
        "",
        "=" * 70,
        "## Breakdown of v2-extension's contribution",
        "",
        f"  total rsIDs ONLY in the 8 new v2 sources        : {n_v2_unique:>5}",
        f"    → of which on-MAF (immediately usable)        : {n_added_to_maf:>5}",
        f"    → of which on-chip-but-filtered (recoverable) : {n_added_recov:>5}",
        f"    → of which not on chip at all                 : "
        f"{n_v2_unique - n_added_to_maf - n_added_recov:>5}",
        "",
        "=" * 70,
        "## Per-phenoclass stratum (PRS-usable rsIDs only)",
        "",
    ])
    for _, r in by_pc.iterrows():
        lines.append(
            f"  {r['main_class']:24s} : "
            f"unique_rsIDs={int(r['unique_rsIDs']):>4}  "
            f"on_MAF={int(r['on_MAF']):>4}  on_chip={int(r['on_chip']):>4}")
    lines.append("")

    # List the new on-MAF SNPs (the most actionable result)
    if n_added_to_maf > 0:
        lines.extend([
            "=" * 70,
            "## NEW on-MAF SNPs added by v2-extension (immediately usable)",
            "",
        ])
        for _, r in delta_maf.iterrows():
            lines.append(
                f"  {r['rsID']:14s}  {r['membership']:50s}  "
                f"pheno={r['pheno_classes']:<22s}  "
                f"gene={r['gene']:<15s}  chr{r['unfilt_CHR_grch37']}:"
                f"{r['unfilt_BP_grch37']}")
        lines.append("")

    if n_added_recov > 0:
        lines.extend([
            "=" * 70,
            "## NEW on-chip-but-filtered SNPs from v2-extension (recoverable)",
            "",
        ])
        for _, r in delta_recov.iterrows():
            lines.append(
                f"  {r['rsID']:14s}  {r['membership']:40s}  "
                f"pheno={r['pheno_classes']:<22s}  "
                f"filter={r['filter_reason']:<40s}  gene={r['gene']}")
        lines.append("")

    # Categorisation key
    lines.extend([
        "=" * 70,
        "## Reconciliation category key (categories used by 30c)",
        "",
        "  on_chip_passed_maf                       : in MAF BIM AND on Illumina chip",
        "  on_chip_filtered_out                     : direct rsID match; lost to QC",
        "  on_chip_via_rsID_merge_filtered_out      : only via RsMergeArch canonical",
        "  on_chip_via_kgp_renamed_filtered_out     : only via Illumina kgp manifest rename",
        "  on_chip_via_position_filtered_out        : only via GRCh37 position",
        "  in_maf_but_not_on_chip_manifest          : in MAF BIM but absent from unfilt",
        "  not_on_chip                              : intrinsic Illumina manifest gap",
        "  apoe_haplotype_special                   : Desikan APOE haplotype rows",
        "",
        "Reasons for absence — see filter_reason column in all_snps_master.tsv:",
        "  missingness>0.02 (--geno step at QC chain step 2)",
        "  HWE<1e-7         (--hwe step at QC chain step 3)",
        "  MAF<0.01         (--maf step at QC chain step 6; monomorphic = MAF=0)",
        "  rename_or_liftover_or_refcorr  (intermediate step 5; rare)",
        "  not_on_chip        : never genotyped by Illumina Omni2.5M chip",
        "",
        "=" * 70,
        "Files written:",
        f"  - all_snps_master.tsv             ({n_master} rows; per-rsID master)",
        f"  - new_snps_from_v2_extension.tsv  ({n_v2_unique} rows)",
        f"  - delta_on_MAF.tsv                ({n_added_to_maf} rows)",
        f"  - delta_recoverable.tsv           ({n_added_recov} rows)",
        f"  - per_source_category_counts.tsv  ({len(per_src)}×N matrix)",
        f"  - per_phenoclass_summary.tsv      ({len(by_pc)} strata)",
        f"  - existing_vs_new_summary.txt     (this file)",
    ])

    txt = "\n".join(lines) + "\n"
    (args.out / "existing_vs_new_summary.txt").write_text(txt, "utf-8")
    print(f"\n[out] {args.out/'existing_vs_new_summary.txt'}")
    print(f"\n=== headline ===")
    print(f"  Pre-v2-extension on-MAF unique:  {on_maf_pre}")
    print(f"  Post-v2-extension on-MAF unique: {on_maf_post}")
    print(f"  NET ADDITION from v2:            {on_maf_post - on_maf_pre}")
    print(f"  Recoverable from v2:             {n_added_recov}")


if __name__ == "__main__":
    main()
