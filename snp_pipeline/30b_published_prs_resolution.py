r"""
30b_published_prs_resolution.py   (LOCAL, env: snp)
====================================================
S1c (user-directed): produce a per-source **full published SNP list
resolution TSV** documenting **why each published rsID is or isn't
retrievable** on the ADNI Omni2.5M MAF>0.01 array. The full MAF BIM
(1,491,403 SNPs) is the right resolution target (NOT the 430-SNP set
which is Wightman-pre-filter + Bellenguez, with no MAF threshold).

Resolution chain (per rsID):
  1. direct rsID in MAF BIM
  2. RsMergeArch canonicalised rsID in MAF BIM
  3. positional match in MAF BIM at the published CHR:BP (GRCh38)
  4. positional match in MAF BIM at lifted CHR:BP (GRCh37 → GRCh38)
  → drop with reason if none succeed.

Original genome assemblies recorded per source:
  Bellenguez 2022   GRCh38   (GWAS-Catalog accession; modern catalog uses GRCh38)
  Wightman 2021     GRCh37   (paper; positions in `wightman_2021_filtered.tsv`)
  Kunkle 2019       GRCh37   (paper; per user)
  Desikan 2017      GRCh37   (PGS000026 — `genome_build=NR` but BIN1
                              chr2:127,892,810 confirms GRCh37)

Output → `D:\ADNI_SNP_Omni2.5M_20140220\source_prs\`:
  <source>_full_snps_resolution.tsv   per-published-SNP row with drop reasons
  source_prs_resolution_summary.txt   per-source totals + drop histograms

Usage:
  conda run -n snp python snp_pipeline/30b_published_prs_resolution.py
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util
import math
import re
from pathlib import Path

import pandas as pd

# Import the stdlib xlsx helper. Use importlib to dodge the
# numeric-leading filename problem from sibling 30* scripts.
_THIS_DIR = Path(__file__).parent
_HELPER = _THIS_DIR / "_xlsx_lead_extract.py"
_spec_xlsx = importlib.util.spec_from_file_location(
    "_xlsx_lead_extract", _HELPER)
_xlsx = importlib.util.module_from_spec(_spec_xlsx)
_spec_xlsx.loader.exec_module(_xlsx)

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
MAF_BED = BASE / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr_maf001"
RSMERGE = BASE / "liftover" / "NCBI" / "RsMergeArch.bcp.gz"
CHAIN = BASE / "liftover" / "hg19ToHg38.over.chain.gz"
G = BASE / "GWAS"
OUT = BASE / "source_prs"

SOURCES = {
    # ── existing 8 (retro-tagged with pheno_class 2026-05-25) ───────────
    "Bellenguez": {
        "file": G / "Bellenguez_2022" / "gwas-association-downloaded_2026-04-28-pubmedId_35379992.tsv",
        "assembly": "GRCh38",     # GWAS-Catalog modern default
        "kind": "gwas_catalog",
        "pheno_class": "clinical_AD",
        # Enrich the 89 catalog leads with raw β + SE + EA + OA from the
        # full sumstats (GCST90027158, 755 MB, downloaded 2026-05-25 night).
        # The catalog TSV only carries OR (not β) and lacks SE/OA — sumstats
        # gives all four cleanly.
        "sumstats_dir": G / "Bellenguez_2022",
    },
    "Wightman": {
        "file": G / "Wightman_2021" / "wightman_2021_filtered.tsv",
        "assembly": "GRCh37",
        "kind": "wightman_filtered",
        "pheno_class": "clinical_AD",
        # β/SE/EAF enrichment: trusted ADNI-hg19 subset (rsID-keyed, 448 SNPs)
        # primary; PGCALZ2 raw meta-analysis (chr:pos GRCh37, 12.2M variants)
        # as fallback for leads not in the trusted subset. Both files carry
        # the published METAL inverse-variance meta β/SE/EAF directly — no
        # β derivation from Z is performed.
        "wightman_trusted_path": G / "Wightman_2021" / "wightman_without_ukb_GRCh38_trusted_adni_hg19.tsv",
        "wightman_meta_path":    G / "Wightman_2021" / "PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis (1).txt" / "PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis (1).txt",
    },
    "Kunkle": {
        "file": G / "Kunkle_2019" / "Kunkle_2019_gwas-association-downloaded_2026-05-18-accessionId_GCST007511.tsv",
        "assembly": "GRCh37",
        "kind": "gwas_catalog",
        "pheno_class": "clinical_AD",
    },
    "Desikan": {
        "file": G / "Desikan_2019" / "PGS000026.txt",
        "assembly": "GRCh37",
        "kind": "pgs_catalog",
        "pheno_class": "clinical_AD",   # polygenic-hazard, clinical case-control
    },
    "Lambert": {
        "file": G / "Lambert_2013" / "gwas-association-downloaded_2026-05-23-pubmedId_24162737.tsv",
        "assembly": "GRCh37",
        "kind": "gwas_catalog",
        "pheno_class": "clinical_AD",
    },
    "DeRojas": {
        "file": G / "DeRojas_2021" / "PGS000898_hmPOS_GRCh37.txt.gz",
        "assembly": "GRCh37",
        "kind": "pgs_catalog",
        "pheno_class": "clinical_AD",
    },
    # Schwanzentruber upgraded 2026-05-25: catalog (30 leads) → xlsx (34 leads).
    "Schwanzentruber": {
        "file": G / "Schwanzentruber_2021" / "Schwanzentruber_2021_Suppl_41588_2020_776_MOESM3_ESM.xlsx",
        "assembly": "GRCh37",
        "kind": "schwanzentruber_xlsx",
        "pheno_class": "clinical_AD",
    },
    "Najar": {
        "file": G / "Najar_2021" / "PGS000812_hmPOS_GRCh37.txt.gz",
        "assembly": "GRCh37",
        "kind": "pgs_catalog",
        "pheno_class": "clinical_AD",
    },
    # ── new 10 (added 2026-05-25) ──────────────────────────────────────
    "Ebenau": {
        "file": G / "Ebenau_2021" / "PGS001775_hmPOS_GRCh37.txt.gz",
        "assembly": "GRCh37",
        "kind": "pgs_catalog",
        "pheno_class": "clinical_AD",
    },
    "Leonenko": {
        "file": G / "Leonenko_2019" / "PGS000876_hmPOS_GRCh37.txt.gz",
        "assembly": "GRCh37",
        "kind": "pgs_catalog",
        "pheno_class": "clinical_AD",
    },
    "Vesilievick": {
        "file": G / "Vesilievick_2023" / "PGS004898_hmPOS_GRCh37.txt.gz",
        "assembly": "GRCh37",
        "kind": "pgs_catalog",
        "pheno_class": "clinical_AD",
    },
    "Zhang": {
        "file": G / "Zhang_2020" / "PGS000334_hmPOS_GRCh37.txt.gz",
        "assembly": "GRCh37",
        "kind": "pgs_catalog",
        "pheno_class": "clinical_AD",
    },
    "Felsky_MF": {
        "file": G / "Felsky_2019" / "Felsky_2019_41467_2018_8279_MOESM7_ESM.xlsx",
        "assembly": "GRCh37",
        "kind": "felsky_xlsx",
        "region": "MF",
        "bypass_exceptions": {"rs2997325"},
        "pheno_class": "pathology_microglia",
    },
    "Felsky_IT": {
        "file": G / "Felsky_2019" / "Felsky_2019_41467_2018_8279_MOESM7_ESM.xlsx",
        "assembly": "GRCh37",
        "kind": "felsky_xlsx",
        "region": "IT",
        "bypass_exceptions": {"rs183093970"},
        "pheno_class": "pathology_microglia",
    },
    "ONeil_NPY": {
        "file": G / "ONeil_2025" / "ONeil_2025_1-s2.0-S0197458025001241-mmc2.xlsx",
        "assembly": "GRCh37",
        "kind": "oneil_xlsx",
        "sheet": "Table_S9_npy_variants",
        "target_rsids": {"rs3940268"},
        "pheno": "NPB",
        "pheno_class": "pathology_plaque",
    },
    "ONeil_GHR": {
        "file": G / "ONeil_2025" / "ONeil_2025_1-s2.0-S0197458025001241-mmc2.xlsx",
        "assembly": "GRCh37",
        "kind": "oneil_xlsx",
        "sheet": "Table_S10_ghr_variants",
        "target_rsids": {"rs10941583"},
        "pheno": "ABD",
        "pheno_class": "pathology_amyloid",
    },
    # Kosteridis re-scoped 2026-05-25 evening: narrowed from the broad
    # 114-SNP MTAG_AD to publication-grade scope — 5 novel AD loci
    # (Table 1) + 33-SNP shared AD-CV set (Supp Tb 3 ∪ Tb 4). The 16
    # Tb4-only colocalising SNPs draw β/SE/P from the full sumstats
    # GCST90449060.h.tsv.gz (marginal AD GWAS) since they're not in
    # Supp1/Supp3 (which only cover GW-sig AD MTAG signals).
    "Kosteridis_novel_AD": {
        "file": G / "Kosteridis_2024" / "Kosteridis_2024_41467_2024_53452_MOESM3_ESM.xlsx",
        "assembly": "GRCh38",
        "kind": "kosteridis_novel_AD",
        "pheno_class": "clinical_AD_MTAG",
    },
    "Kosteridis_shared_AD_CV": {
        "file": G / "Kosteridis_2024" / "Kosteridis_2024_41467_2024_53452_MOESM3_ESM.xlsx",
        "assembly": "GRCh38",
        "kind": "kosteridis_shared_AD_CV",
        "pheno_class": "clinical_AD_MTAG",
        "sumstats_path": G / "Kosteridis_2024" / "GCST90449060.h.tsv.gz",
    },
}

# ── Yang 2023 cell-type-partitioned PRS sub-sources (added 2026-05-25 night) ─
# Bellenguez 2022 GWAS SNPs (PRset file) filtered to ROSMAP single-cell
# expression markers from Mathys et al. 2019. Lead set = PRset
# `<celltype>==1 AND P<5e-8`. β/SE from raw Bellenguez 2022 sumstats
# (GCST90027158). Six parallel sub-sources, one per cell type.
_YANG_PRSET = G / "Yang_2023" / "PRset_ROSMAP_25OCT2023+direct.snp"
_BELL_SUMSTATS_DIR = G / "Bellenguez_2022"
for _ct, _pheno in [
    ("Ex",  "pathology_neuron_excitatory"),
    ("In",  "pathology_neuron_inhibitory"),
    ("Ast", "pathology_astrocyte"),
    ("Mic", "pathology_microglia"),         # joins Felsky_MF/IT
    ("Oli", "pathology_oligodendrocyte"),
    ("Opc", "pathology_opc"),
]:
    SOURCES[f"Yang_{_ct}"] = {
        "file": _YANG_PRSET,
        "assembly": "GRCh38",
        "kind": "yang_celltype",
        "cell_type": _ct,
        "pheno_class": _pheno,
        "sumstats_dir": _BELL_SUMSTATS_DIR,
    }


def _build_rsmerge_map(needed_rsids: set[str]) -> dict[str, str]:
    """RsMergeArch low→high merge map for the given rsIDs. Returns
    {old_rsID: current_rsID}."""
    out: dict[str, str] = {}
    needles = {r.lstrip("rs").lstrip("0") for r in needed_rsids
               if r.startswith("rs")}
    if not RSMERGE.exists() or not needles:
        return out
    with gzip.open(RSMERGE, "rt") as f:
        for line in f:
            t = line.rstrip("\n").split("\t")
            if len(t) < 3:
                continue
            lo, hi = t[0], t[1]
            if lo in needles:
                out[f"rs{lo}"] = f"rs{hi}"
    return out


def _parse_gwas_catalog(path: Path) -> pd.DataFrame:
    """GWAS-Catalog .tsv → DataFrame[rsID_pub, CHR_pub, BP_pub,
    effect_allele_pub, OR_or_beta_pub, P_pub]. Handles Bellenguez/Kunkle."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    cols = {c.upper().strip(): c for c in df.columns}

    def col(*aliases):
        for a in aliases:
            k = a.upper().strip()
            if k in cols:
                return cols[k]
        return None
    c_rs = col("SNPS")
    c_chr = col("CHR_ID")
    c_bp = col("CHR_POS")
    c_risk = col("STRONGEST SNP-RISK ALLELE")
    c_or = col("OR or BETA", "OR_or_BETA")
    c_p = col("P-VALUE", "P_VALUE")
    keep = []
    for _, r in df.iterrows():
        rs = str(r.get(c_rs, "")).strip()
        if not rs.startswith("rs"):
            continue
        risk = str(r.get(c_risk, "")).strip() if c_risk else ""
        ea = risk.split("-")[-1].upper() if "-" in risk else ""
        if ea in {"", "?"}:
            ea = ""
        keep.append({
            "rsID_pub": rs,
            "CHR_pub": str(r.get(c_chr, "")).strip() if c_chr else "",
            "BP_pub": str(r.get(c_bp, "")).strip() if c_bp else "",
            "effect_allele_pub": ea,
            "OR_or_beta_pub": str(r.get(c_or, "")).strip() if c_or else "",
            "P_pub": str(r.get(c_p, "")).strip() if c_p else "",
        })
    return pd.DataFrame(keep).drop_duplicates("rsID_pub")


def _parse_wightman_filtered(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    return df.rename(columns={
        "Lead_variant": "rsID_pub",
        "Chromosome": "CHR_pub",
        "Position_GRCh37": "BP_pub",
        "A1": "effect_allele_pub",
        "P_value": "P_pub",
    }).assign(OR_or_beta_pub="")[
        ["rsID_pub", "CHR_pub", "BP_pub", "effect_allele_pub",
         "OR_or_beta_pub", "P_pub"]]


def _parse_pgs_catalog(path: Path) -> pd.DataFrame:
    """PGS Catalog .txt or .txt.gz → DataFrame. Skips '#' header lines,
    treats APOE haplotype rows (is_haplotype=TRUE) as special (no rsID
    resolution). PGS files vary in column set — some omit `other_allele`
    (Najar/Zhang) and add `hm_inferOtherAllele` from the hmPOS harmoniser;
    we don't need other_allele in 30b's emitted schema, but we do tolerate
    its absence by using `.get(..., "")` everywhere.

    Multi-allele warn (Vesilievick PGS004898): if `hm_inferOtherAllele`
    contains `/` (e.g. `C/G/T`), the row is retained but flagged with
    `qc_verdict='warn_multi_allele'`.
    """
    rows = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as f:
        in_body = False
        header = None
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            if not in_body:
                header = line.split("\t")
                in_body = True
                continue
            t = line.split("\t")
            if len(t) < len(header):
                t += [""] * (len(header) - len(t))
            r = dict(zip(header, t))
            rs = r.get("rsID", "").strip()
            oa_raw = (r.get("other_allele", "")
                      or r.get("hm_inferOtherAllele", "")).strip()
            qc = "warn_multi_allele" if "/" in oa_raw else ""
            rows.append({
                "rsID_pub": rs,
                "CHR_pub": r.get("chr_name", "").strip(),
                "BP_pub": r.get("chr_position", "").strip(),
                "effect_allele_pub": r.get("effect_allele", "").strip(),
                "other_allele_pub": oa_raw,
                "OR_or_beta_pub": r.get("effect_weight", "").strip(),
                "SE_pub": "",
                "P_pub": "",
                "locus_name": r.get("locus_name", "").strip(),
                "is_haplotype": r.get("is_haplotype", "").strip(),
                "qc_verdict": qc,
            })
    return pd.DataFrame(rows)


def _parse_felsky_xlsx(path: Path, region: str,
                       bypass_exceptions: set | None = None
                       ) -> pd.DataFrame:
    """Felsky 2019 PAM-GWAS xlsx supplementary; filter by region (MF/IT)."""
    df = _xlsx.extract_felsky(path, bypass_exceptions=bypass_exceptions)
    if df.empty:
        return df
    # Bypass-exception rows added with no region — assign target region
    if bypass_exceptions:
        mask_bp = df["published_exception"] & (df["GWAS_Region"] == "")
        df.loc[mask_bp, "GWAS_Region"] = region
    df = df[df["GWAS_Region"] == region].copy()
    df = df.rename(columns={
        "rsID": "rsID_pub",
        "CHR": "CHR_pub",
        "BP": "BP_pub",
        "A1": "effect_allele_pub",
        "A2": "other_allele_pub",
        "Beta": "OR_or_beta_pub",
        "SE": "SE_pub",
        "P": "P_pub",
    })
    return df[[
        "rsID_pub", "CHR_pub", "BP_pub", "effect_allele_pub",
        "other_allele_pub", "OR_or_beta_pub", "SE_pub", "P_pub",
        "published_exception",
    ]]


def _parse_oneil_xlsx(path: Path, sheet: str, target_rsids: set,
                      pheno: str, N: int | None = None) -> pd.DataFrame:
    """ONeil 2025 mmc2 — extract single SNP per source by target rsID."""
    df = _xlsx.extract_oneil(path, sheet=sheet, target_rsids=target_rsids,
                              pheno=pheno, N=N)
    if df.empty:
        return df
    df = df.rename(columns={
        "rsID": "rsID_pub",
        "CHR": "CHR_pub",
        "BP": "BP_pub",
        "Beta": "OR_or_beta_pub",
        "SE": "SE_pub",
        "P": "P_pub",
    })
    df["effect_allele_pub"] = ""
    df["other_allele_pub"] = ""
    df["locus_name"] = df.get("Locus", "")
    return df[[
        "rsID_pub", "CHR_pub", "BP_pub", "effect_allele_pub",
        "other_allele_pub", "OR_or_beta_pub", "SE_pub", "P_pub",
        "locus_name",
    ]]


KOSTERIDIS_NOVEL_AD_RSIDS = (
    "rs7529220",     # HSPG2
    "rs11692604",    # AC019055.1
    "rs73069394",    # ULK4
    "rs77399788",    # KRT18P16
    "rs56365761",    # ACTN4
)


def _lookup_kosteridis_sumstats(sumstats_path: Path,
                                  target_rsids: set[str]
                                  ) -> dict[str, dict]:
    """Stream the gzipped Kosteridis full sumstats (264 MB) once and return
    {rsID: {chr, bp, EA, OA, beta, SE, P, EAF}} for the target rsIDs.
    Used as fallback AD-β source for the Tb4-only colocalising signals
    that don't appear in Supp1/Supp3."""
    import gzip
    out: dict[str, dict] = {}
    remaining = set(target_rsids)
    if not sumstats_path.exists():
        print(f"  [WARN] Kosteridis sumstats not found: {sumstats_path}")
        return out
    print(f"  [sumstats] streaming {sumstats_path.name} for "
          f"{len(remaining)} target rsIDs …")
    with gzip.open(sumstats_path, "rt", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        rs_col = idx.get("rsid", idx.get("rs_id"))
        for line in f:
            if not remaining:
                break
            t = line.rstrip("\n").split("\t")
            if len(t) <= rs_col:
                continue
            rs = t[rs_col]
            if rs in remaining:
                out[rs] = {
                    "CHR": t[idx["chromosome"]],
                    "BP": t[idx["base_pair_location"]],
                    "EA": t[idx["effect_allele"]],
                    "OA": t[idx["other_allele"]],
                    "beta": t[idx["beta"]],
                    "SE": t[idx["standard_error"]],
                    "P": t[idx["p_value"]],
                    "EAF": t[idx.get("effect_allele_frequency",
                                       len(t) - 1)] if idx.get(
                        "effect_allele_frequency") else "",
                }
                remaining.discard(rs)
    print(f"  [sumstats] {len(out)}/{len(target_rsids)} resolved; "
          f"unresolved: {sorted(remaining) if remaining else '(none)'}")
    return out


def _lookup_bellenguez_sumstats(sumstats_path: Path,
                                  target_rsids: set[str]
                                  ) -> dict[str, dict]:
    """Stream the EBI Bellenguez 2022 full sumstats (GCST90027158, 755 MB)
    once and return {rsID: {chr, bp, EA, OA, beta, SE, P, EAF}} for the
    target rsIDs. Used by Yang_<celltype> sub-sources to get raw β/SE for
    cell-type-filtered SNPs that are NOT in Bellenguez's 89-lead Catalog
    TSV.

    EBI GWAS-SSF schema (column names verified):
        variant_id, p_value, chromosome, base_pair_location, effect_allele,
        other_allele, effect_allele_frequency, odds_ratio, ci_lower,
        ci_upper, beta, standard_error, n_cases, n_controls, het_isq,
        het_pvalue, variant_alternate_id
    Note rsID column = `variant_id` (not `rsid` like Kosteridis sumstats)."""
    import gzip
    out: dict[str, dict] = {}
    remaining = set(target_rsids)
    if not sumstats_path.exists():
        print(f"  [WARN] Bellenguez sumstats not found: {sumstats_path}")
        return out
    print(f"  [sumstats] streaming {sumstats_path.name} for "
          f"{len(remaining)} target rsIDs …")
    with gzip.open(sumstats_path, "rt", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        # Auto-detect rsID column — Bellenguez uses `variant_id`, Kosteridis
        # uses `rsid`. Fall back to either.
        rs_col = idx.get("variant_id", idx.get("rsid", idx.get("rs_id")))
        if rs_col is None:
            print(f"  [WARN] no rsID column found in sumstats header: {header[:10]}")
            return out
        for line in f:
            if not remaining:
                break
            t = line.rstrip("\n").split("\t")
            if len(t) <= rs_col:
                continue
            rs = t[rs_col]
            if rs in remaining:
                out[rs] = {
                    "CHR": t[idx["chromosome"]],
                    "BP": t[idx["base_pair_location"]],
                    "EA": t[idx["effect_allele"]],
                    "OA": t[idx["other_allele"]],
                    "beta": t[idx["beta"]],
                    "SE": t[idx["standard_error"]],
                    "P": t[idx["p_value"]],
                    "EAF": (t[idx["effect_allele_frequency"]]
                            if "effect_allele_frequency" in idx else ""),
                }
                remaining.discard(rs)
    print(f"  [sumstats] {len(out)}/{len(target_rsids)} resolved; "
          f"unresolved: {len(remaining)} "
          f"({sorted(remaining)[:5] if remaining else ''})")
    return out


def _parse_yang_celltype(prset_path: Path, cell_type: str,
                          sumstats_dir: Path,
                          p_threshold: float = 5e-8) -> pd.DataFrame:
    """Yang 2023 cell-type-partitioned PRS sub-source.
    Reads the PRset .snp file, filters to <cell_type>==1 AND P<p_threshold,
    then looks up raw β/SE from Bellenguez 2022 sumstats."""
    df = pd.read_csv(prset_path, sep="\t", dtype=str, keep_default_na=False)
    df["P_num"] = pd.to_numeric(df["P"], errors="coerce")
    df[cell_type + "_num"] = pd.to_numeric(df[cell_type],
                                            errors="coerce").fillna(0).astype(int)
    df = df[(df[cell_type + "_num"] == 1)
             & (df["P_num"] < p_threshold)].copy()
    target_rsids = set(df["SNP"])

    # Find the Bellenguez sumstats file by glob
    candidates = list(sumstats_dir.glob("*GCST90027158_buildGRCh38.tsv.gz")) \
        or list(sumstats_dir.glob("*GCST90027158*.tsv.gz")) \
        or list(sumstats_dir.glob("*.tsv.gz"))
    if not candidates:
        raise FileNotFoundError(
            f"Bellenguez sumstats not found in {sumstats_dir} — "
            f"download GCST90027158 from EBI GWAS Catalog.")
    sumstats_path = candidates[0]
    sumstats_hits = _lookup_bellenguez_sumstats(sumstats_path, target_rsids)

    rows = []
    for _, r in df.iterrows():
        rs = r["SNP"]
        h = sumstats_hits.get(rs, {})
        rows.append({
            "rsID_pub": rs,
            "CHR_pub": r["CHR"],
            "BP_pub": r["BP"],
            "effect_allele_pub": h.get("EA", ""),
            "other_allele_pub": h.get("OA", ""),
            "OR_or_beta_pub": h.get("beta", ""),
            "SE_pub": h.get("SE", ""),
            "P_pub": r["P"],   # PRset P == Bellenguez P (same source)
            "locus_name": "",
            "evidence_source": ("Bellenguez_sumstats" if rs in sumstats_hits
                                  else "sumstats_not_resolved"),
            "cell_type": cell_type,
        })
    return pd.DataFrame(rows)


def _parse_kosteridis_novel_AD(path: Path) -> pd.DataFrame:
    """Kosteridis 2024 publication Table 1 — 5 novel AD loci.
    Pulls β/SE/P from Supp Data 1 (all 5 are in there, verified)."""
    s1 = _xlsx.extract_kosteridis_supp1(path)
    df = s1[s1["rsID"].isin(KOSTERIDIS_NOVEL_AD_RSIDS)].copy()
    df["evidence_source"] = "Supp1_MTAG"
    df["colocalisation_PP"] = ""
    df["colocalisation_trait"] = ""
    df = df.rename(columns={
        "rsID": "rsID_pub",
        "CHR": "CHR_pub",
        "BP": "BP_pub",
        "EA": "effect_allele_pub",
        "OA": "other_allele_pub",
        "Beta": "OR_or_beta_pub",
        "SE": "SE_pub",
        "P": "P_pub",
        "Gene": "locus_name",
    })
    return df[[
        "rsID_pub", "CHR_pub", "BP_pub", "effect_allele_pub",
        "other_allele_pub", "OR_or_beta_pub", "SE_pub", "P_pub",
        "locus_name", "MTAG_origins", "evidence_source",
        "colocalisation_PP", "colocalisation_trait",
    ]]


def _parse_kosteridis_shared_AD_CV(path: Path,
                                     sumstats_path: Path | None = None
                                     ) -> pd.DataFrame:
    """Kosteridis 2024 shared AD-CV set = Supp Tb 3 (15 CV_GWS=yes) ∪
    Supp Tb 4 (21 colocalising) = 33 unique rsIDs.
    For rsIDs in Supp1 (17/33): use AD MTAG β/SE/P from Supp1.
    For rsIDs ONLY in Tb4 (16/33): use marginal AD β/SE/P from the full
    sumstats GCST90449060.h.tsv.gz as fallback (tagged
    evidence_source='sumstats_fallback').
    Each row also carries colocalisation_PP + colocalisation_trait from
    Tb4 if present, else empty."""
    s1 = _xlsx.extract_kosteridis_supp1(path)
    s1_by_rs = {r["rsID"]: r for _, r in s1.iterrows()}
    s3 = _xlsx.extract_kosteridis_supp3(path, supp1_df=s1)
    s4 = _xlsx.extract_kosteridis_supp4(path)
    union_rs = set(s3["rsID"]) | set(s4["rsID"])
    s4_by_rs = {r["rsID"]: r for _, r in s4.iterrows()}

    # Build a lookup for Tb3 origin tracking
    tb3_rsids = set(s3["rsID"])
    tb4_rsids = set(s4["rsID"])
    def _origin(rs: str) -> str:
        a, b = rs in tb3_rsids, rs in tb4_rsids
        if a and b: return "Tb3+Tb4"
        if a:       return "Tb3_only"
        if b:       return "Tb4_only"
        return ""

    # SNPs that need sumstats fallback (in union but not in Supp1)
    in_s1 = set(s1_by_rs)
    fallback_targets = union_rs - in_s1
    sumstats_hits = (_lookup_kosteridis_sumstats(sumstats_path, fallback_targets)
                     if sumstats_path else {})

    rows = []
    for rs in sorted(union_rs):
        s4row = s4_by_rs.get(rs, {})
        if rs in in_s1:
            r = s1_by_rs[rs]
            rows.append({
                "rsID_pub": rs,
                "CHR_pub": r["CHR"],
                "BP_pub": r["BP"],
                "effect_allele_pub": r["EA"],
                "other_allele_pub": r["OA"],
                "OR_or_beta_pub": r["Beta"],
                "SE_pub": r["SE"],
                "P_pub": r["P"],
                "locus_name": r["Gene"],
                "MTAG_origins": r.get("MTAG_origins", ""),
                "evidence_source": "Supp1_MTAG",
                "colocalisation_PP": s4row.get("PP", "") if rs in tb4_rsids else "",
                "colocalisation_trait": s4row.get("Traits", "") if rs in tb4_rsids else "",
                "tb_origin": _origin(rs),
            })
        elif rs in sumstats_hits:
            h = sumstats_hits[rs]
            rows.append({
                "rsID_pub": rs,
                "CHR_pub": h["CHR"],
                "BP_pub": h["BP"],
                "effect_allele_pub": h["EA"],
                "other_allele_pub": h["OA"],
                "OR_or_beta_pub": h["beta"],
                "SE_pub": h["SE"],
                "P_pub": h["P"],
                "locus_name": s4row.get("Gene", ""),
                "MTAG_origins": "",
                "evidence_source": "sumstats_fallback",
                "colocalisation_PP": s4row.get("PP", ""),
                "colocalisation_trait": s4row.get("Traits", ""),
                "tb_origin": _origin(rs),
            })
        else:
            # Last-resort: no β anywhere — present row with NA β
            rows.append({
                "rsID_pub": rs,
                "CHR_pub": s4row.get("CHR", ""),
                "BP_pub": s4row.get("BP", ""),
                "effect_allele_pub": "",
                "other_allele_pub": "",
                "OR_or_beta_pub": float("nan"),
                "SE_pub": float("nan"),
                "P_pub": float("nan"),
                "locus_name": s4row.get("Gene", ""),
                "MTAG_origins": "",
                "evidence_source": "no_AD_beta",
                "colocalisation_PP": s4row.get("PP", ""),
                "colocalisation_trait": s4row.get("Traits", ""),
                "tb_origin": _origin(rs),
            })
    df = pd.DataFrame(rows)
    return df[[
        "rsID_pub", "CHR_pub", "BP_pub", "effect_allele_pub",
        "other_allele_pub", "OR_or_beta_pub", "SE_pub", "P_pub",
        "locus_name", "MTAG_origins", "evidence_source",
        "colocalisation_PP", "colocalisation_trait", "tb_origin",
    ]]


def _parse_schwanzentruber_xlsx(path: Path) -> pd.DataFrame:
    """Schwanzentruber 2021 supplementary `1-AD loci` sheet — 37 leads."""
    df = _xlsx.extract_schwanzentruber_supp1(path)
    if df.empty:
        return df
    df = df.rename(columns={
        "rsID": "rsID_pub",
        "CHR": "CHR_pub",
        "BP": "BP_pub",
        "EA": "effect_allele_pub",
        "OA": "other_allele_pub",
        "Beta": "OR_or_beta_pub",
        "SE": "SE_pub",
        "P": "P_pub",
        "Nearest_gene": "locus_name",
    })
    return df[[
        "rsID_pub", "CHR_pub", "BP_pub", "effect_allele_pub",
        "other_allele_pub", "OR_or_beta_pub", "SE_pub", "P_pub",
        "locus_name",
    ]]


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    # ── Load resolution targets ──────────────────────────────────────────
    print(f"[bim] loading MAF BIM …")
    bim = pd.read_csv(MAF_BED.with_suffix(".bim"), sep="\t", header=None,
                       names=["CHR", "rsID", "cM", "BP", "A1", "A2"],
                       dtype=str)
    bim["BP"] = pd.to_numeric(bim["BP"], errors="coerce").astype("Int64")
    bim_rs = set(bim["rsID"])
    bim_by_rs = bim.set_index("rsID")
    bim_by_pos = bim.set_index(
        bim["CHR"].astype(str) + ":" + bim["BP"].astype(str))
    print(f"  MAF BIM: {len(bim)} SNPs × 616 patients")

    # liftover
    lo = None
    try:
        from pyliftover import LiftOver
        if CHAIN.exists():
            lo = LiftOver(str(CHAIN))
            print(f"[liftover] hg19→hg38 chain ready")
    except ImportError:
        pass

    def _lift_grch37_to_grch38(chrom: str, bp_str: str) -> int | None:
        if lo is None or not chrom or not bp_str:
            return None
        try:
            bp = int(bp_str) - 1                          # liftover wants 0-based
        except ValueError:
            return None
        hits = lo.convert_coordinate(f"chr{chrom}", bp)
        if not hits:
            return None
        return int(hits[0][1]) + 1                       # 1-based back

    # ── Process each source ──────────────────────────────────────────────
    # Kosteridis Supp3 depends on Supp1's SE values — parse Supp1 first
    # and cache the DF so Supp3 can join.
    _supp1_cache: dict[str, pd.DataFrame] = {}
    all_summary = []
    for src, spec in SOURCES.items():
        path = spec["file"]
        assembly = spec["assembly"]
        kind = spec["kind"]
        pheno_class = spec.get("pheno_class", "")
        print(f"\n[{src}] parsing {path.name} (assembly: {assembly}, "
              f"pheno_class: {pheno_class})")
        if kind == "gwas_catalog":
            df = _parse_gwas_catalog(path)
        elif kind == "wightman_filtered":
            df = _parse_wightman_filtered(path)
        elif kind == "pgs_catalog":
            df = _parse_pgs_catalog(path)
        elif kind == "felsky_xlsx":
            df = _parse_felsky_xlsx(path, region=spec["region"],
                                     bypass_exceptions=spec.get("bypass_exceptions"))
        elif kind == "oneil_xlsx":
            df = _parse_oneil_xlsx(path, sheet=spec["sheet"],
                                    target_rsids=spec["target_rsids"],
                                    pheno=spec["pheno"], N=spec.get("N"))
        elif kind == "kosteridis_novel_AD":
            df = _parse_kosteridis_novel_AD(path)
        elif kind == "kosteridis_shared_AD_CV":
            df = _parse_kosteridis_shared_AD_CV(
                path, sumstats_path=spec.get("sumstats_path"))
        elif kind == "yang_celltype":
            df = _parse_yang_celltype(path,
                                       cell_type=spec["cell_type"],
                                       sumstats_dir=spec["sumstats_dir"])
        elif kind == "schwanzentruber_xlsx":
            df = _parse_schwanzentruber_xlsx(path)
        else:
            print(f"  [WARN] unknown kind: {kind} — skipping {src}")
            continue

        # Optional: enrich OR_or_beta_pub / SE_pub / EA / OA from full
        # sumstats for sources where the parser couldn't get them. For
        # Bellenguez, the catalog gives OR with no SE — sumstats has clean
        # β + SE per rsID.
        sumstats_dir = spec.get("sumstats_dir")
        if sumstats_dir is not None and kind != "yang_celltype":
            candidates = (list(sumstats_dir.glob("*GCST90027158_buildGRCh38.tsv.gz"))
                           or list(sumstats_dir.glob("*GCST*.tsv.gz")))
            if candidates and "rsID_pub" in df.columns:
                rsids = set(df["rsID_pub"])
                hits = _lookup_bellenguez_sumstats(candidates[0], rsids)
                if hits:
                    def _enrich(row):
                        h = hits.get(row["rsID_pub"])
                        if not h:
                            return row
                        # Preserve catalog values in *_catalog columns for
                        # cross-check; then take EA + OA + β + SE from
                        # sumstats as one coherent triple (mixing sources
                        # risks EA/OA orientation mismatch with β sign).
                        row["OR_or_beta_pub_catalog"] = row.get(
                            "OR_or_beta_pub", "")
                        row["effect_allele_pub_catalog"] = row.get(
                            "effect_allele_pub", "")
                        row["effect_allele_pub"] = h["EA"]
                        row["other_allele_pub"] = h["OA"]
                        row["OR_or_beta_pub"] = h["beta"]
                        row["SE_pub"] = h["SE"]
                        return row
                    df = df.apply(_enrich, axis=1)
                    print(f"  [enrich] {len(hits)}/{len(rsids)} {src} rsIDs "
                          f"populated with β/SE/EA/OA from sumstats")

        # Wightman-specific enrichment: trusted ADNI-hg19 subset (rsID) → fallback
        # to PGCALZ2 raw meta (chr:pos GRCh37). Both have published METAL β/SE/EAF.
        if src == "Wightman" and spec.get("wightman_trusted_path"):
            tpath = spec["wightman_trusted_path"]
            mpath = spec.get("wightman_meta_path")
            tdf = pd.read_csv(tpath, sep="\t") if tpath.exists() else None
            t_by_rsid = ({} if tdf is None else
                         {str(r["rsID"]): r for _, r in tdf.iterrows()})
            # Build chr:pos lookup of leads not in trusted, for meta-fallback
            need_pos = set()
            if "rsID_pub" in df.columns and "BP_pub" in df.columns:
                for _, r in df.iterrows():
                    if str(r["rsID_pub"]) not in t_by_rsid:
                        try:
                            need_pos.add((int(r["CHR_pub"]), int(r["BP_pub"])))
                        except (ValueError, TypeError):
                            pass
            m_by_pos = {}
            if mpath and mpath.exists() and need_pos:
                for ch in pd.read_csv(mpath, sep=r"\s+", engine="python",
                                       chunksize=300000):
                    sub = ch[ch.set_index(["chromosome", "base_pair_location"])
                               .index.isin(need_pos)]
                    for _, r in sub.iterrows():
                        m_by_pos[(int(r["chromosome"]),
                                    int(r["base_pair_location"]))] = r
                    if len(m_by_pos) >= len(need_pos):
                        break

            n_t = n_m = 0

            def _enrich_w(row):
                nonlocal n_t, n_m
                rs = str(row["rsID_pub"])
                if rs in t_by_rsid:
                    h = t_by_rsid[rs]
                    row["effect_allele_pub"] = str(h["effect_allele"]).upper()
                    row["other_allele_pub"] = str(h["other_allele"]).upper()
                    row["OR_or_beta_pub"] = h["beta"]
                    row["SE_pub"] = h["standard_error"]
                    n_t += 1
                    return row
                try:
                    key = (int(row["CHR_pub"]), int(row["BP_pub"]))
                except (ValueError, TypeError):
                    return row
                if key in m_by_pos:
                    h = m_by_pos[key]
                    row["effect_allele_pub"] = str(h["effect_allele"]).upper()
                    row["other_allele_pub"] = str(h["other_allele"]).upper()
                    row["OR_or_beta_pub"] = h["beta"]
                    row["SE_pub"] = h["standard_error"]
                    n_m += 1
                return row

            # Ensure columns exist
            for col in ("other_allele_pub", "SE_pub"):
                if col not in df.columns:
                    df[col] = ""
            df = df.apply(_enrich_w, axis=1)
            print(f"  [enrich-W] Wightman β: {n_t} from trusted file, "
                  f"{n_m} from PGCALZ2 raw meta (chr:pos), "
                  f"{len(df) - n_t - n_m} unfilled")

        n_pub = len(df)
        print(f"  published entries: {n_pub}")

        # RsMergeArch lookup for ALL rsIDs in this source
        rsmap = _build_rsmerge_map(set(df["rsID_pub"]))
        df["rsID_canonical"] = df["rsID_pub"].map(
            lambda r: rsmap.get(r, r) if isinstance(r, str) else r)
        if rsmap:
            print(f"  RsMergeArch hits: {len(rsmap)}")

        rows = []
        for _, r in df.iterrows():
            rs_pub = r["rsID_pub"]
            rs_can = r["rsID_canonical"]
            chrom = str(r["CHR_pub"]).strip()
            bp_str = str(r["BP_pub"]).strip()
            ea = str(r.get("effect_allele_pub", "")).strip()
            method, rs_match, chr38, bp38, a1, a2 = (None, "", "", "",
                                                       "", "")
            drop = "ok"

            # APOE haplotype rows (Desikan) — handle separately
            if (kind == "pgs_catalog"
                    and str(r.get("is_haplotype", "")).upper() == "TRUE"):
                drop = "apoe_haplotype_special"
                method = "haplotype"
            elif not rs_pub.startswith("rs"):
                drop = "no_rsID_in_source"
            elif rs_pub in bim_rs:
                method, rs_match = "rsID", rs_pub
            elif rs_can in bim_rs:
                method, rs_match = "rsID_merged", rs_can
            elif chrom and bp_str:
                # build-specific positional matching
                if assembly == "GRCh38":
                    key = (f"{chrom}:{bp_str}"
                            if bp_str.lstrip("-").isdigit() else None)
                    if key and key in bim_by_pos.index:
                        row = bim_by_pos.loc[key]
                        if isinstance(row, pd.DataFrame):
                            row = row.iloc[0]
                        method = "position_grch38"
                        rs_match = str(row.name) if hasattr(row, "name") else ""
                # liftover GRCh37→GRCh38 if not GRCh38 source OR if GRCh38
                # match failed
                if method is None and assembly in ("GRCh37", "NR"):
                    bp_lifted = _lift_grch37_to_grch38(chrom, bp_str)
                    if bp_lifted is not None:
                        key = f"{chrom}:{bp_lifted}"
                        if key in bim_by_pos.index:
                            row = bim_by_pos.loc[key]
                            if isinstance(row, pd.DataFrame):
                                row = row.iloc[0]
                            method = "position_lifted_grch37_to_38"
                            rs_match = (str(row.name)
                                          if hasattr(row, "name") else "")
                            chr38, bp38 = chrom, str(bp_lifted)

            if method is not None and rs_match and rs_match in bim_by_rs.index:
                br = bim_by_rs.loc[rs_match]
                if isinstance(br, pd.DataFrame):
                    br = br.iloc[0]
                chr38 = chr38 or str(br["CHR"])
                bp38 = bp38 or str(br["BP"])
                a1 = str(br["A1"]).upper()
                a2 = str(br["A2"]).upper()
            elif drop == "ok":
                drop = "not_on_array_or_below_maf"
                method = "none"

            rows.append({
                "source": src,
                "pheno_class": pheno_class,
                "rsID_pub": rs_pub,
                "rsID_canonical": rs_can,
                "original_assembly": assembly,
                "CHR_pub": chrom,
                "BP_pub": bp_str,
                "effect_allele_pub": ea,
                "other_allele_pub": (r.get("other_allele_pub", "")
                                       if "other_allele_pub" in r else ""),
                "OR_or_beta_pub": r.get("OR_or_beta_pub", ""),
                "SE_pub": (r.get("SE_pub", "") if "SE_pub" in r else ""),
                "P_pub": r.get("P_pub", ""),
                "locus_name": r.get("locus_name", "")
                                if "locus_name" in r else "",
                "is_haplotype": r.get("is_haplotype", "")
                                if "is_haplotype" in r else "",
                "MTAG_origins": (r.get("MTAG_origins", "")
                                   if "MTAG_origins" in r else ""),
                "evidence_source": (r.get("evidence_source", "")
                                      if "evidence_source" in r else ""),
                "colocalisation_PP": (r.get("colocalisation_PP", "")
                                        if "colocalisation_PP" in r else ""),
                "colocalisation_trait": (r.get("colocalisation_trait", "")
                                           if "colocalisation_trait" in r else ""),
                "tb_origin": (r.get("tb_origin", "")
                                if "tb_origin" in r else ""),
                "cell_type": (r.get("cell_type", "")
                                if "cell_type" in r else ""),
                "qc_verdict": (r.get("qc_verdict", "")
                                 if "qc_verdict" in r else ""),
                "published_exception": (r.get("published_exception", "")
                                         if "published_exception" in r else ""),
                "rsID_in_bim": rs_match,
                "CHR_GRCh38": chr38,
                "BP_GRCh38": bp38,
                "A1_in_bim": a1,
                "A2_in_bim": a2,
                "resolved_method": method,
                "drop_reason": drop,
            })

        out_df = pd.DataFrame(rows)
        out_path = args.out / f"{src}_full_snps_resolution.tsv"
        out_df.to_csv(out_path, sep="\t", index=False)
        n_ok = int((out_df["drop_reason"] == "ok").sum())
        n_hap = int((out_df["drop_reason"] == "apoe_haplotype_special").sum())
        n_drop = int(((out_df["drop_reason"] != "ok") &
                      (out_df["drop_reason"] != "apoe_haplotype_special")).sum())
        print(f"  resolved: {n_ok}, APOE-haplotype: {n_hap}, dropped: {n_drop}")
        print(f"  drop reasons: "
              f"{dict(out_df[out_df['drop_reason'].isin(['not_on_array_or_below_maf', 'no_rsID_in_source'])]['drop_reason'].value_counts())}")
        print(f"  resolved methods: "
              f"{dict(out_df[out_df['drop_reason'] == 'ok']['resolved_method'].value_counts())}")
        all_summary.append({
            "source": src, "assembly": assembly, "n_published": n_pub,
            "n_resolved": n_ok, "n_haplotype": n_hap, "n_dropped": n_drop,
            "by_method": dict(out_df[out_df["drop_reason"] == "ok"]
                              ["resolved_method"].value_counts()),
        })
        print(f"  [out] {out_path}")

    # ── Summary text ─────────────────────────────────────────────────────
    lines = [
        "Published-PRS SNP resolution summary",
        f"Resolution target: MAF>0.01 BIM ({len(bim)} SNPs × 616 patients)",
        "Resolution chain: direct rsID → RsMergeArch canonical → positional",
        "  (GRCh38 if source assembly is GRCh38; else liftover GRCh37→GRCh38)",
        "",
    ]
    for s in all_summary:
        lines.append(f"## {s['source']}  (assembly: {s['assembly']})")
        lines.append(f"  published entries  : {s['n_published']}")
        lines.append(f"  resolved (in MAF)  : {s['n_resolved']}")
        if s["n_haplotype"]:
            lines.append(f"  APOE-haplotype     : {s['n_haplotype']}")
        lines.append(f"  dropped            : {s['n_dropped']}  "
                     f"(reason: not_on_array_or_below_maf)")
        if s["by_method"]:
            lines.append("  resolved by method :")
            for m, n in s["by_method"].items():
                lines.append(f"    {m:32s} {n}")
        lines.append("")
    (args.out / "source_prs_resolution_summary.txt").write_text(
        "\n".join(lines), "utf-8")
    print(f"\n[out] {args.out/'source_prs_resolution_summary.txt'}")


if __name__ == "__main__":
    main()
