"""
_prs_strict_qc_lib.py
=====================
Shared helpers for strict-QC + LD-pruned PRS analyses (scripts 45/46/47).

Given a source's per-SNP resolution TSV and the LD-pruned strict-QC dosage
matrix, this module:
  - filters to resolved SNPs (`drop_reason == 'ok'`) intersected with the
    LD-pruned strict-QC SNP set,
  - orientation-aligns published β to the in-BIM A1 (flip sign if needed),
  - converts OR → log(OR) where the source reports OR not beta,
  - drops palindromic / mismatched SNPs that can't be safely flipped,
  - returns (1) a per-source SNP table and (2) per-patient PRS series.

Per the user's existing pipeline ([[feedback-strict-qc-ld-pruning-first]]):
raw β only; no shrinkage. LD pruning was applied upstream (snp_pipeline/44).
"""
from __future__ import annotations
import os
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd

ROOT = Path("D:/ADNI_SNP_Omni2.5M_20140220")
SRC_DIR = ROOT / "source_prs"
QC_DIR  = ROOT / "GWAS_comprehensive_v2" / "QC_strict"
DEFAULT_LD_CONFIG = "ld_1000kb_r2_0.8"


def get_ldpruned_paths(ld_config: str = DEFAULT_LD_CONFIG):
    """Return (bim_path, dosage_path) for the given LD config subdir."""
    d = QC_DIR / ld_config
    return (d / "recover_all_pool_strictQC_LDpruned.bim",
            d / "recover_all_pool_strictQC_LDpruned_dosage.tsv")

# Sources known to report OR (positive, log-OR-needs-conversion).
# Detected by inspect on 2026-05-27. Others report beta directly.
OR_SOURCES = {"Kunkle", "Lambert", "ONeil_GHR"}

# Order matches the per-source resolution TSV files in source_prs/.
# Note: "Kosteridis" alone has no resolution TSV — the main Kosteridis 2024
# sumstats are split into MTAG_AD / shared_AD_CV / novel_AD per the
# project_prs_source_extension memory.
ALL_SOURCES = [
    "Bellenguez", "Wightman", "Kunkle", "Schwanzentruber",
    "Lambert", "DeRojas", "Desikan", "Ebenau", "Najar", "Leonenko",
    "Vesilievick", "Zhang",
    "Felsky_MF", "Felsky_IT",
    "ONeil_NPY", "ONeil_GHR",
    "Kosteridis_novel_AD", "Kosteridis_shared_AD_CV", "Kosteridis_MTAG_AD",
]

# Per-source effective GWAS N (from paper abstracts; used to break ties when
# multiple sources contribute the same rsID in the dedup combined PRS).
# Largest-N source wins. See feedback_multi_method_prs_comparison / project_prs_source_extension.
N_GWAS = {
    "Bellenguez":               487511,
    "Wightman":                 472868,
    "Schwanzentruber":          446092,
    "Kosteridis_MTAG_AD":       436000,
    "Kosteridis_shared_AD_CV":  436000,
    "Kosteridis_novel_AD":      436000,
    "Desikan":                   70000,   # IGAP-combined PGS bundle
    "Kunkle":                    63926,
    "Lambert":                   54162,
    "DeRojas":                   25000,   # SPIGAPUK2 effective N
    "Najar":                      8000,   # Swedish cohort
    "Vesilievick":                7000,
    "Ebenau":                     6000,
    "Leonenko":                   6000,
    "Zhang":                      5000,
    "Felsky_IT":                  2000,
    "Felsky_MF":                  2000,
    "ONeil_GHR":                  1500,
    "ONeil_NPY":                  1500,
}

# Standard ε2-allele extraction from APOE_Alleles column ("3/3", "3/4", "2/4", "4/4", "2/3").
def _apoe2_dosage(alleles_str: str) -> int:
    if not isinstance(alleles_str, str): return 0
    return sum(1 for a in alleles_str.split("/") if a.strip() == "2")


def _complement(base: str) -> str:
    return {"A":"T","T":"A","C":"G","G":"C"}.get(str(base).upper(), str(base).upper())


def _is_palindromic(a1: str, a2: str) -> bool:
    return {a1.upper(), a2.upper()} in ({"A","T"}, {"C","G"})


def load_ldpruned_dosage(ld_config: str = DEFAULT_LD_CONFIG
                            ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return (bim_df, dosage_df). Dosage indexed by PTID, columns = rsIDs."""
    bim_p, dos_p = get_ldpruned_paths(ld_config)
    bim = pd.read_csv(bim_p, sep=r"\s+", header=None,
                       names=["chrom","rsID","cM","bp","A1","A2"], dtype=str)
    dos = pd.read_csv(dos_p, sep="\t", dtype=str)
    dos["PTID"] = dos["PTID"].astype(str)
    dos = dos.set_index("PTID")
    for c in dos.columns:
        dos[c] = pd.to_numeric(dos[c], errors="coerce")
    # cohort-mean impute (used as default; train-fold impute is per-fit)
    dos = dos.apply(lambda c: c.fillna(c.mean()))
    return bim, dos


def load_source_beta(src: str, bim: pd.DataFrame) -> pd.DataFrame:
    """Return per-source SNP table with cols: rsID, beta_A1, A1, A2, drop_reason.

    Logic:
      * Read <SRC>_full_snps_resolution.tsv.
      * Keep rows with drop_reason == 'ok' AND rsID_in_bim is non-null.
      * Intersect rsID_in_bim with the LD-pruned BIM.
      * Resolve sign of beta:
          - If effect_allele_pub == A1_in_bim   →  beta_A1 = +β_pub
          - If effect_allele_pub == A2_in_bim   →  beta_A1 = -β_pub
          - If effect_allele_pub == complement(A1_in_bim) AND not palindromic
                                                →  beta_A1 = +β_pub  (strand-flipped)
          - If effect_allele_pub == complement(A2_in_bim) AND not palindromic
                                                →  beta_A1 = -β_pub
          - Else (palindromic, mismatched)      →  drop
      * If src in OR_SOURCES, convert OR → log(OR) before sign-flip.
    """
    p = SRC_DIR / f"{src}_full_snps_resolution.tsv"
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t")

    # filter to resolved-OK rows
    if "drop_reason" not in df.columns:
        return pd.DataFrame()
    df = df[df["drop_reason"] == "ok"].copy()
    df = df.dropna(subset=["rsID_in_bim"])
    if df.empty:
        return pd.DataFrame()

    # numeric β
    df["beta_pub"] = pd.to_numeric(df["OR_or_beta_pub"], errors="coerce")
    df = df.dropna(subset=["beta_pub"])
    if src in OR_SOURCES:
        df["beta_pub"] = np.log(df["beta_pub"])

    # intersect with LD-pruned BIM
    bim_idx = bim.set_index("rsID")
    keep = df["rsID_in_bim"].isin(bim_idx.index)
    df = df.loc[keep].copy()
    if df.empty:
        return pd.DataFrame()
    df["A1_in_bim"] = df["rsID_in_bim"].map(bim_idx["A1"])
    df["A2_in_bim"] = df["rsID_in_bim"].map(bim_idx["A2"])

    # orient β to A1_in_bim. Handle missing other_allele_pub by only checking
    # effect_allele_pub against A1_in_bim / A2_in_bim (or complements).
    def _orient(row):
        ea = str(row["effect_allele_pub"]).upper()
        if ea in ("NAN", "NA", ""):
            return np.nan
        a1 = str(row["A1_in_bim"]).upper()
        a2 = str(row["A2_in_bim"]).upper()
        if _is_palindromic(a1, a2):
            return np.nan
        # If other_allele_pub is known, sanity-check the allele pair matches
        oa = row.get("other_allele_pub", None)
        oa = str(oa).upper() if pd.notna(oa) else None
        if oa is not None and oa not in ("NAN", "NA", ""):
            if {ea, oa} != {a1, a2} and {ea, oa} != {_complement(a1), _complement(a2)}:
                return np.nan
        # Match ea to A1 / A2 / strand-flip
        if ea == a1:
            return row["beta_pub"]
        if ea == a2:
            return -row["beta_pub"]
        if ea == _complement(a1):
            return row["beta_pub"]
        if ea == _complement(a2):
            return -row["beta_pub"]
        return np.nan
    df["beta_A1"] = df.apply(_orient, axis=1)
    df = df.dropna(subset=["beta_A1"]).copy()
    out = df.rename(columns={"rsID_in_bim": "rsID"})[
        ["rsID", "beta_A1", "A1_in_bim", "A2_in_bim", "beta_pub",
         "effect_allele_pub", "other_allele_pub"]]
    return out.reset_index(drop=True)


def compute_prs(src_snp: pd.DataFrame, dosage: pd.DataFrame) -> pd.Series:
    """Σ_i β_A1_i · dosage_p,i over SNPs present in both src_snp and dosage."""
    cols = [r for r in src_snp["rsID"] if r in dosage.columns]
    if not cols:
        return pd.Series(0.0, index=dosage.index)
    w = src_snp.set_index("rsID").loc[cols, "beta_A1"].astype(float).values
    return dosage[cols].astype(float).mul(w, axis=1).sum(axis=1)


def per_source_prs_table(sources: List[str] | None = None,
                            ld_config: str = DEFAULT_LD_CONFIG,
                            include_dedup: bool = True
                           ) -> Tuple[pd.DataFrame, Dict[str, pd.DataFrame]]:
    """Build per-source PRS for the entire 616-subj cohort.
    Returns:
      df_prs: DataFrame with PTID + one PRS column per source ("PRS_<source>")
              + (if include_dedup) "PRS_prs_all_dedup"
      df_snp: dict source → SNP table (rsID, beta_A1, alleles, etc.)
              + (if include_dedup) "prs_all_dedup" → composed SNP table
    """
    if sources is None:
        sources = ALL_SOURCES
    bim, dos = load_ldpruned_dosage(ld_config)
    cols = {"PTID": list(dos.index)}
    df_snp = {}
    for s in sources:
        st = load_source_beta(s, bim)
        if st.empty:
            cols[f"PRS_{s}"] = [np.nan] * len(dos)
            df_snp[s] = st
            continue
        cols[f"PRS_{s}"] = compute_prs(st, dos).values
        df_snp[s] = st

    if include_dedup:
        dedup_st = build_prs_all_dedup(df_snp)
        if not dedup_st.empty:
            cols["PRS_prs_all_dedup"] = compute_prs(dedup_st, dos).values
        else:
            cols["PRS_prs_all_dedup"] = [np.nan] * len(dos)
        df_snp["prs_all_dedup"] = dedup_st
    return pd.DataFrame(cols), df_snp


def build_prs_all_dedup(df_snp: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Build a single rsID-deduplicated SNP table by taking, for each rsID,
    the β from the source with the largest published GWAS N (per N_GWAS dict).

    Inputs: df_snp = {source_name: SNP_table_for_that_source}
    Output: DataFrame with cols [rsID, beta_A1, A1_in_bim, A2_in_bim, beta_pub,
            effect_allele_pub, other_allele_pub, contributing_source]
    """
    rows = []
    # Sort sources by N_GWAS descending so first-seen wins
    ordered = sorted(df_snp.keys(), key=lambda s: -N_GWAS.get(s, 0))
    seen = set()
    for src in ordered:
        st = df_snp.get(src, pd.DataFrame())
        if st.empty: continue
        for _, r in st.iterrows():
            rs = r["rsID"]
            if rs in seen:
                continue
            seen.add(rs)
            rows.append({**r.to_dict(), "contributing_source": src})
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).reset_index(drop=True)


def get_dedup_dosage_matrix(ld_config: str = DEFAULT_LD_CONFIG
                             ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return (snp_table, dosage_subset) limited to rsIDs that appear in the
    dedup SNP table for this LD config (i.e., at least one source provides β).

    snp_table: same shape as build_prs_all_dedup output.
    dosage_subset: PTID-indexed dosage DataFrame restricted to those rsIDs.
    """
    bim, dos = load_ldpruned_dosage(ld_config)
    df_snp = {}
    for s in ALL_SOURCES:
        st = load_source_beta(s, bim)
        if not st.empty: df_snp[s] = st
    dedup = build_prs_all_dedup(df_snp)
    keep_cols = [r for r in dedup["rsID"] if r in dos.columns]
    return dedup, dos[keep_cols]


def load_subject_covariates(splits_root: Path | None = None) -> pd.DataFrame:
    """Return DataFrame Patient_ID -> {Sex_M, age_at_baseline, APOE4_Dosage,
    APOE2_Dosage}. Reads from the post-exclusion baseline splits' train CSV
    of seed 0 as the canonical source (covariates are subject-level and
    identical across seeds/splits)."""
    if splits_root is None:
        splits_root = Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                            "no_cdr_stratified_post_exclusion/tabular/baseline")
    parts = []
    for sp in ("train", "val", "test"):
        df = pd.read_csv(splits_root / f"seed_0/{sp}.csv", dtype=str)
        parts.append(df)
    df = pd.concat(parts, ignore_index=True).drop_duplicates(subset=["Patient_ID"]).copy()
    df["Patient_ID"] = df["Patient_ID"].astype(str)
    df["Sex_M"] = (df["Sex"] == "Male").astype(int)
    df["APOE4_Dosage"] = pd.to_numeric(df["APOE4_Dosage"], errors="coerce")
    df["APOE2_Dosage"] = df["APOE_Alleles"].apply(_apoe2_dosage)
    yob = pd.to_numeric(df["YOB"], errors="coerce").astype("Int64")
    date = pd.to_datetime(df["Date"], errors="coerce")
    birth = pd.to_datetime(yob.astype("Int64").astype(str) + "-07-01",
                            errors="coerce")
    df["age_at_baseline"] = (date - birth).dt.days / 365.25
    return df[["Patient_ID","Sex_M","APOE4_Dosage","APOE2_Dosage","age_at_baseline"]]


if __name__ == "__main__":
    import sys
    ld_config = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_LD_CONFIG
    df, snps = per_source_prs_table(ld_config=ld_config)
    print(f"\n[LD config = {ld_config}]")
    print(f"Per-source SNP intersect:")
    for s, st in snps.items():
        col = f"PRS_{s}"
        print(f"  {s:30s}  n_SNPs={len(st):3d}  "
              f"PRS_mean={df[col].mean():+.4f}  std={df[col].std():.4f}")
    print(f"\nSubject covariate sample:")
    cov = load_subject_covariates()
    print(cov.head(3).to_string(index=False))
    print(f"  total subjects with covariates: {len(cov)}")
