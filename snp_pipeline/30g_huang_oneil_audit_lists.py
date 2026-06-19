r"""
30g_huang_oneil_audit_lists.py   (LOCAL, env: snp)
==================================================
Emit the two audit-only TSVs that are NOT PRS sources but feed into the
overlap audit (and 30c reconciliation) for cross-source presence checks.

Outputs in `source_prs/`:
  Huang_lead_snps_no_beta.tsv     — 22 AAOS lead variants (Huang 2017,
                                    `is_query_snp == 1`). No β/SE — used
                                    only to check if other sources' PRS
                                    cover the same Huang loci.
  ONeil_SST_candidates.tsv        — ~262 SST-neuron region variants
                                    with P<0.05 in ≥1 of {NPB, ABD,
                                    AD-Status}, sorted by `n_sig_phenos`
                                    desc (3-of-3 first). For user
                                    review before promotion to PRS
                                    source.

Usage:
  conda run -n snp python snp_pipeline/30g_huang_oneil_audit_lists.py
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
G = BASE / "GWAS"
OUT = BASE / "source_prs"
MAF_BED = BASE / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr_maf001"

HUANG_XLSX = G / "Huang_2017" / "Huang_2017_AAO_41593_2017_BFnn4587_MOESM6_ESM.xlsx"
ONEIL_XLSX = G / "ONeil_2025" / "ONeil_2025_1-s2.0-S0197458025001241-mmc2.xlsx"

# Import the stdlib xlsx helper
_HELPER = Path(__file__).with_name("_xlsx_lead_extract.py")
_spec = importlib.util.spec_from_file_location("_xlsx_lead_extract", _HELPER)
_xlsx = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_xlsx)


def _load_bim_rs(bim_path: Path) -> set:
    df = pd.read_csv(bim_path, sep="\t", header=None, usecols=[1],
                     names=["rsID"], dtype=str)
    return set(df["rsID"])


def _resolution_row(rs: str, pheno_class: str, src: str, in_maf: bool,
                    locus: str = "", **extra) -> dict:
    """Minimal 30b-compatible resolution row for audit-only sources."""
    return {
        "source": src,
        "pheno_class": pheno_class,
        "rsID_pub": rs,
        "rsID_canonical": rs,
        "original_assembly": "GRCh37",
        "CHR_pub": extra.get("CHR_pub", ""),
        "BP_pub": extra.get("BP_pub", ""),
        "effect_allele_pub": extra.get("effect_allele_pub", ""),
        "other_allele_pub": extra.get("other_allele_pub", ""),
        "OR_or_beta_pub": extra.get("OR_or_beta_pub", ""),
        "SE_pub": extra.get("SE_pub", ""),
        "P_pub": extra.get("P_pub", ""),
        "locus_name": locus,
        "is_haplotype": "",
        "MTAG_origins": "",
        "qc_verdict": "",
        "published_exception": "",
        "rsID_in_bim": rs if in_maf else "",
        "CHR_GRCh38": "",
        "BP_GRCh38": "",
        "A1_in_bim": "",
        "A2_in_bim": "",
        "resolved_method": "rsID" if in_maf else "none",
        "drop_reason": "ok" if in_maf else "not_on_array_or_below_maf",
        "audit_only": True,
    }


def write_huang(out_dir: Path, maf_rs: set) -> Path:
    df = _xlsx.extract_huang(HUANG_XLSX)
    df["source"] = "Huang"
    df["pheno_class"] = "clinical_AAO"
    df["audit_only"] = True
    cols = ["source", "pheno_class", "rsID", "CHR", "BP_GRCh38",
            "ref", "alt", "audit_only"]
    p = out_dir / "Huang_lead_snps_no_beta.tsv"
    df[cols].to_csv(p, sep="\t", index=False)
    print(f"  [out] {p}  ({len(df)} rows)")

    # 30b-compatible resolution TSV so 30c can process it
    res_rows = [_resolution_row(rs=row["rsID"], pheno_class="clinical_AAO",
                                  src="Huang",
                                  in_maf=row["rsID"] in maf_rs,
                                  CHR_pub=row["CHR"],
                                  BP_pub=row["BP_GRCh38"])
                for _, row in df.iterrows()]
    res = pd.DataFrame(res_rows)
    res_path = out_dir / "Huang_full_snps_resolution.tsv"
    res.to_csv(res_path, sep="\t", index=False)
    print(f"  [out] {res_path}  ({len(res)} rows, schema-compatible)")
    return p


def write_oneil_sst_candidates(out_dir: Path, maf_rs: set,
                                p_threshold: float = 0.05) -> Path:
    df = _xlsx.extract_oneil(ONEIL_XLSX,
                              sheet="Table_S8_sst_variants",
                              p_threshold=p_threshold)
    if df.empty:
        print(f"  [warn] no ONeil SST candidates at P<{p_threshold}")
        return out_dir / "ONeil_SST_candidates.tsv"
    df["source"] = "ONeil_SST_candidates"
    df["pheno_class"] = "pathology_mixed"
    df["audit_only"] = True
    # Order: n_sig_phenos desc, then min-P asc to surface strongest first
    df["_min_p"] = df[["NPB_P", "ABD_P", "AD_P"]].min(axis=1)
    df = df.sort_values(["n_sig_phenos", "_min_p"], ascending=[False, True])
    cols = ["source", "pheno_class", "rsID", "CHR", "BP", "LOC", "MAF",
            "NPB_Beta", "NPB_SE", "NPB_P",
            "ABD_Beta", "ABD_SE", "ABD_P",
            "AD_Beta", "AD_P",
            "n_sig_phenos", "sig_phenos", "Locus", "audit_only"]
    p = out_dir / "ONeil_SST_candidates.tsv"
    df[cols].to_csv(p, sep="\t", index=False)
    n_3 = (df["n_sig_phenos"] == 3).sum()
    n_2 = (df["n_sig_phenos"] == 2).sum()
    n_1 = (df["n_sig_phenos"] == 1).sum()
    print(f"  [out] {p}  ({len(df)} rows; "
          f"3-of-3 sig={n_3}, 2-of-3 sig={n_2}, 1-of-3 sig={n_1})")

    # 30b-compatible resolution TSV
    res_rows = [_resolution_row(rs=row["rsID"], pheno_class="pathology_mixed",
                                 src="ONeil_SST_candidates",
                                 in_maf=row["rsID"] in maf_rs,
                                 locus=str(row.get("Locus") or ""),
                                 CHR_pub=row["CHR"],
                                 BP_pub=row["BP"])
                for _, row in df.iterrows()]
    res = pd.DataFrame(res_rows)
    res_path = out_dir / "ONeil_SST_candidates_full_snps_resolution.tsv"
    res.to_csv(res_path, sep="\t", index=False)
    print(f"  [out] {res_path}  ({len(res)} rows, schema-compatible)")
    return p


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--p-threshold", type=float, default=0.05,
                    help="ONeil SST candidate filter (default 0.05).")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    print(f"[audit] writing audit-only TSVs to {args.out}")
    print(f"[audit] loading MAF BIM rsID set …")
    maf_rs = _load_bim_rs(MAF_BED.with_suffix(".bim"))
    print(f"  MAF BIM: {len(maf_rs):,} rsIDs")
    write_huang(args.out, maf_rs)
    write_oneil_sst_candidates(args.out, maf_rs,
                                p_threshold=args.p_threshold)


if __name__ == "__main__":
    main()
