r"""
25b_export_apoe_control_sets.py   (LOCAL, env: snp)
===================================================
Export two APOE positive-control "SNP sets" into GWAS_filtered/ in the
SAME script-25 schema (<SET>_snps.tsv + <SET>_patient_dosage.tsv +
<SET>_patient_dosage.bim) so the diff-attention pipeline (26 → 27 → 28)
runs on them identically to 9W27B/9W27B5D. The honest test: can the
frozen-FM diff-embedding pipeline detect the strongest genetic AD signal
(APOE) that classical logistic/PRS trivially gets?

  APOE_cluster  the genotyped chr19 APOE LD-block (script-20 `onlyAPOE`
                tier: chr19 44.4–46.5 Mb, minus the isolated rs2072561 and
                the monomorphic/absent ε-SNPs) from the original 430 — real
                array genotype, large β. effect_allele = A1 = ALT
                (dosage-counted), β = harmonised beta_A1.
  APOE_causal   synthetic rs429358 (ε4) + rs7412 (ε2): REF/ALT 1 kb
                sequence comes from the FASTA (script 26, positions below),
                `dosage` = the CLINICAL ε4/ε2 allele count (script-20
                `_covars`, since the Omni2.5M array does not genotype
                rs429358/rs7412), β = Desikan APOE haplotype weights
                (ε4 +1.03, ε2 −0.47). Tests the actual causal variants
                with no array genotype needed.

Reuses (importlib, unmodified) from 20_genetic_baselines.py: build_beta_table,
snps_for_tier, _covars, DOSAGE_TSV/DOSAGE_BIM, SPLITS, E_SNPS_DEFAULT,
APOE_CLUSTER_DEFAULT. Output schema matches 25_export_filtered_prs_sets.py.

Usage:
  conda run -n snp python snp_pipeline/25b_export_apoe_control_sets.py
"""
from __future__ import annotations

import argparse
import importlib.util as _il
from pathlib import Path

import pandas as pd

_S20 = Path(__file__).parent / "20_genetic_baselines.py"
_spec = _il.spec_from_file_location("script20", _S20)
_s20 = _il.module_from_spec(_spec)
_spec.loader.exec_module(_s20)
build_beta_table = _s20.build_beta_table
snps_for_tier = _s20.snps_for_tier
_covars = _s20._covars
DOSAGE_TSV = _s20.DOSAGE_TSV
DOSAGE_BIM = _s20.DOSAGE_BIM
SPLITS = _s20.SPLITS
E_SNPS_DEFAULT = _s20.E_SNPS_DEFAULT
APOE_CLUSTER_DEFAULT = _s20.APOE_CLUSTER_DEFAULT
DATA = _s20.DATA

GWAS_FILTERED = DATA / "GWAS_filtered"
DESIKAN_W = DATA / "genetic_baselines" / "desikan_apoe_weights.tsv"
# Well-established GRCh38 (NC_000019.10) coordinates / alleles for the two
# APOE-defining SNPs. REF = genomic reference; ALT = the ε-risk/ε2 allele
# that beta + dosage refer to (script-25 convention: effect_allele=ALT=A1).
APOE_CAUSAL = {
    "rs429358": {"CHR": "19", "BP": 44908684, "REF": "T", "ALT": "C",
                 "tag": "e4"},   # ε4 determinant; ALT=C is the ε4 allele
    "rs7412":   {"CHR": "19", "BP": 44908822, "REF": "C", "ALT": "T",
                 "tag": "e2"},   # ε2 determinant; ALT=T is the ε2 allele
}
SNP_COLS = ["rsID", "gene", "lead_rsID", "source", "origin", "CHR",
            "BP_GRCh38", "effect_allele", "other_allele", "REF", "ALT",
            "beta_A1"]


def _write_set(name: str, snp_df: pd.DataFrame, dosage: pd.DataFrame,
               bim: pd.DataFrame, readme: str, out_root: Path) -> None:
    d = out_root / name
    d.mkdir(parents=True, exist_ok=True)
    snp_df = snp_df[SNP_COLS].copy()
    snp_df.to_csv(d / f"{name}_snps.tsv", sep="\t", index=False)
    dosage.to_csv(d / f"{name}_patient_dosage.tsv", sep="\t", index=False)
    bim[["rsID", "CHR", "BP", "A1", "A2"]].to_csv(
        d / f"{name}_patient_dosage.bim", sep="\t", index=False)
    (d / "README.txt").write_text(readme, encoding="utf-8")
    print(f"[out] {d}  ({len(snp_df)} SNPs, dosage "
          f"{len(dosage)}×{len(snp_df)})")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=GWAS_FILTERED)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    # ── APOE_cluster (genotyped onlyAPOE tier from the 430) ──────────────
    print("[β] build_beta_table (unfiltered 430) …", flush=True)
    bt = build_beta_table()
    cc = snps_for_tier(bt, list(bt.index), "onlyAPOE",
                       E_SNPS_DEFAULT, APOE_CLUSTER_DEFAULT)
    cc = [s for s in cc if pd.notna(bt.loc[s, "beta_A1"])]
    sub = bt.loc[cc]
    snp_cluster = pd.DataFrame({
        "rsID": cc, "gene": "APOE_region", "lead_rsID": cc,
        "source": sub["source"].values,
        "origin": "onlyAPOE_chr19_LDblock",
        "CHR": sub["CHR"].astype(str).values,
        "BP_GRCh38": sub["BP"].astype(int).values,
        "effect_allele": sub["A1"].astype(str).values,   # = ALT = dosage A1
        "other_allele": sub["A2"].astype(str).values,
        "REF": sub["A2"].astype(str).values,
        "ALT": sub["A1"].astype(str).values,
        "beta_A1": sub["beta_A1"].astype(float).values})

    dos = pd.read_csv(DOSAGE_TSV, sep="\t")
    idc = "PTID" if "PTID" in dos.columns else dos.columns[0]
    dos = dos.rename(columns={idc: "PTID"})
    dos["PTID"] = dos["PTID"].astype(str)
    keep = [c for c in cc if c in dos.columns]
    missing = [c for c in cc if c not in dos.columns]
    if missing:
        print(f"  [warn] {len(missing)} onlyAPOE SNPs not in dosage tsv "
              f"→ dropped: {missing[:5]}")
        snp_cluster = snp_cluster[snp_cluster["rsID"].isin(keep)]
    dos_cluster = dos[["PTID"] + keep]
    bim_all = pd.read_csv(DOSAGE_BIM, sep="\t")
    bim_all["rsID"] = bim_all["rsID"].astype(str)
    bim_cluster = bim_all[bim_all["rsID"].isin(keep)][
        ["rsID", "CHR", "BP", "A1", "A2"]]
    _write_set("APOE_cluster", snp_cluster, dos_cluster, bim_cluster,
               "APOE_cluster — genotyped chr19 APOE LD-block (script-20 "
               "onlyAPOE tier; original 430). effect_allele=A1=ALT=dosage-"
               "counted; beta_A1 harmonised. Genotyped FM positive control "
               "for ad_case vs stable_CN — classical PRS/logistic clears "
               "chance here (OR/SD≈1.6); if the diff-attention pipeline "
               "cannot, the FM embeddings are not allele-informative.\n",
               args.out)

    # ── APOE_causal (synthetic rs429358/rs7412; clinical ε-dosage) ───────
    dw = pd.read_csv(DESIKAN_W, sep="\t", index_col=0)["weight"]
    b_e4 = float(dw.get("APOE_e4", 1.03))
    b_e2 = float(dw.get("APOE_e2", -0.47))
    rows = []
    for rs, m in APOE_CAUSAL.items():
        rows.append({
            "rsID": rs, "gene": "APOE", "lead_rsID": rs,
            "source": "Desikan2019_PGS_APOE",
            "origin": "apoe_causal_synthetic", "CHR": m["CHR"],
            "BP_GRCh38": m["BP"], "effect_allele": m["ALT"],
            "other_allele": m["REF"], "REF": m["REF"], "ALT": m["ALT"],
            "beta_A1": (b_e4 if m["tag"] == "e4" else b_e2)})
    snp_causal = pd.DataFrame(rows)

    # per-patient clinical ε4/ε2 dosage (union of seed_0 splits = all pts)
    root = SPLITS["bl_multi"]
    parts = [pd.read_csv(root / "seed_0" / f"{sp}.csv")
             for sp in ("train", "val", "test")]
    cov = _covars(pd.concat(parts, ignore_index=True)
                  .drop_duplicates("Patient_ID"))
    cov.index = cov.index.astype(str)
    dos_causal = pd.DataFrame({
        "PTID": cov.index,
        "rs429358": pd.to_numeric(cov["APOE4"], errors="coerce").values,
        "rs7412": pd.to_numeric(cov["APOE2"], errors="coerce").values,
    }).drop_duplicates("PTID").reset_index(drop=True)
    bim_causal = pd.DataFrame([
        {"rsID": rs, "CHR": m["CHR"], "BP": m["BP"],
         "A1": m["ALT"], "A2": m["REF"]}     # A1=effect/ALT, A2=other/REF
        for rs, m in APOE_CAUSAL.items()])
    n_e4 = int(dos_causal["rs429358"].notna().sum())
    _write_set("APOE_causal", snp_causal, dos_causal, bim_causal,
               f"APOE_causal — synthetic rs429358(ε4,+{b_e4})/"
               f"rs7412(ε2,{b_e2}). REF/ALT 1 kb sequence built from the "
               "GRCh38 FASTA by script 26 (positions are well-established; "
               "script-26 FASTA-REF QC validates them). dosage = CLINICAL "
               f"ε-allele count from APOE_Alleles ({n_e4}/{len(dos_causal)} "
               "patients with APOE) — the array does not genotype these "
               "SNPs. Tests whether the FM sees the actual causal APOE "
               "variants.\n", args.out)

    print("\n[done] APOE_cluster + APOE_causal exported. Next: "
          "26_build_fm_snp_sequences.py --set APOE_cluster (and APOE_causal),"
          " then 27 (Colab) + 28.")


if __name__ == "__main__":
    main()
