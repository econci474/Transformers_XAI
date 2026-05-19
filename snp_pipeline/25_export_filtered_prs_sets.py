r"""
25_export_filtered_prs_sets.py
==============================
Export the three named filtered GW-significant PRS SNP sets used by the
filtered genetic baselines (script 20, `genetic_baselines_filtered/`
`genetic_baselines_withdesikan/`) and the filtered Cox survival models
(script 23, `cox_survival_filtered/` `cox_withdesikan/`) as
self-contained, training-ready bundles:

  6W27B    29 noAPOE SNPs  = 6 in-array Wightman-2021 GW-sig leads
                              + 27 Bellenguez-2022 leads  (the base β table)
  9W27B    32              = 6W27B + 3 re-extracted Wightman leads
                              (CLNK rs4504245, GRN rs708382, CD33 rs1354106)
  9W27B5D  37              = 9W27B + 5 non-overlapping novel Desikan-2019
                              PGS loci (GPR115/BC043356/ZCWPW1/SLC24A4/TRIP4)

The SNP set + β + dosage are reconstructed by **exactly** the script-23
`covar_matrix` logic (build_beta_table base ∪ enriched novel, restricted to
the noAPOE tier and to genotyped SNPs present in the filtered dosage
matrix) so the export is byte-identical to what the models trained on.

Per set → D:\ADNI_SNP_Omni2.5M_20140220\GWAS_filtered\<SET>\ :
  <SET>_snps.tsv            rsID, gene, lead_rsID, source, origin, CHR,
                            BP_GRCh38, effect_allele, other_allele, REF,
                            ALT, beta_A1                       (per SNP)
  <SET>_patient_dosage.tsv  PTID + one column per set rsID  (additively
                            coded 0/1/2 = ALT/effect-allele dosage)
  <SET>_patient_dosage.bim  plink-style rsID CHR BP A1 A2 (A1 = counted/
                            effect/ALT allele)               (subset)
  README.txt                column + PRS-reconstruction semantics

Effect-allele / dosage convention (same as the models): the dosage counts
the **A1** allele (plink --recode A --keep-allele-order); `beta_A1` is the
published GWAS effect size harmonised to that A1. For the refcorr GRCh38
array A1 = ALT, A2 = REF. The weighted PRS for a set is
  PRS = Σ_i  beta_A1_i × dosage_i        (missing dosage → cohort mean)
APOE-region SNPs are excluded (noAPOE tier) — identical to the models'
`betaPRS_<set>` / `dosage_<set>`.

Usage:  conda run -n snp python snp_pipeline/25_export_filtered_prs_sets.py
"""
from __future__ import annotations

import argparse
import importlib.util as _il
from pathlib import Path

import pandas as pd

# ── Reuse script 20 (digit-leading filename → importlib; same as script 23) ──
_S20 = Path(__file__).parent / "20_genetic_baselines.py"
_spec = _il.spec_from_file_location("script20", _S20)
_s20 = _il.module_from_spec(_spec)
_spec.loader.exec_module(_s20)
build_beta_table = _s20.build_beta_table
snps_for_tier = _s20.snps_for_tier
_load_wightman_allowlist = _s20._load_wightman_allowlist
E_SNPS_DEFAULT = _s20.E_SNPS_DEFAULT
APOE_CLUSTER_DEFAULT = _s20.APOE_CLUSTER_DEFAULT

GBF = Path("D:/ADNI_SNP_Omni2.5M_20140220/genetic_baselines_filtered")
OUT_DEFAULT = Path("D:/ADNI_SNP_Omni2.5M_20140220/GWAS_filtered")
ALLOWLIST = GBF / "wightman_resolved_allowlist.tsv"
ENR_W = GBF / "wfilt_desikan_enriched_weights.tsv"      # all 8 novel (3W+5D)
DOS_TSV = GBF / "wfilt_dosage" / "patient_dosage.tsv"
DOS_BIM = GBF / "wfilt_dosage" / "patient_dosage.bim"
FILT_PRS = GBF / "filtered_prs_snps.tsv"                # provenance, 9W+27B

SNP_COLS = ["rsID", "gene", "lead_rsID", "source", "origin", "CHR",
            "BP_GRCh38", "effect_allele", "other_allele", "REF", "ALT",
            "beta_A1"]


def _provenance() -> dict:
    """rsID → {gene, lead_rsID, source, origin}, merged from
    filtered_prs_snps.tsv (the 9W+27B) and the enriched-weights novel rows
    (3 re-extracted Wightman + 5 novel Desikan)."""
    prov: dict[str, dict] = {}
    fp = pd.read_csv(FILT_PRS, sep="\t")
    for r in fp.itertuples(index=False):
        prov[str(r.rsID)] = {"gene": getattr(r, "gene", ""),
                             "lead_rsID": str(r.rsID),
                             "source": str(getattr(r, "source", "")),
                             "origin": str(getattr(r, "origin", ""))}
    ew = pd.read_csv(ENR_W, sep="\t")
    for r in ew[ew["novel"] == True].itertuples(index=False):   # noqa: E712
        sc = str(getattr(r, "source", ""))
        prov[str(r.rsID)] = {
            "gene": str(getattr(r, "Gene", "") or ""),
            "lead_rsID": str(getattr(r, "lead_rsID", r.rsID) or r.rsID),
            "source": sc,
            "origin": ("desikan_novel" if "desikan" in sc.lower()
                       else "re_extracted")}
    return prov


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    args = ap.parse_args()

    # ── reconstruct comb EXACTLY as 23_cox_survival_models.covar_matrix ──────
    wig = _load_wightman_allowlist(ALLOWLIST)
    print("[β] building harmonised β table …", flush=True)
    bt = build_beta_table(wightman_allowlist=wig)
    ew = pd.read_csv(ENR_W, sep="\t", low_memory=False)
    comb = bt[["beta_A1", "CHR", "BP"]].copy()
    comb["CHR"] = comb["CHR"].astype(str)
    nv = ew[ew["novel"] == True].dropna(subset=["beta_A1"])      # noqa: E712
    for r in nv.itertuples(index=False):
        comb.loc[str(r.rsID)] = {"beta_A1": float(r.beta_A1),
                                 "CHR": str(r.CHR).split(".")[0],
                                 "BP": int(float(r.BP))}
    comb["BP"] = comb["BP"].astype(int)
    src = ew.get("source", pd.Series([""] * len(ew))).astype(str)
    w_novel = [str(r) for r, sc, nvl in zip(ew["rsID"], src, ew["novel"])
               if nvl and "wightman" in sc.lower()]
    d_novel = [str(r) for r, sc, nvl in zip(ew["rsID"], src, ew["novel"])
               if nvl and "desikan" in sc.lower()]
    base = list(bt.index)

    dos = pd.read_csv(DOS_TSV, sep="\t")
    dos["PTID"] = dos["PTID"].astype(str)
    dos_cols = set(dos.columns) - {"PTID"}
    bim = pd.read_csv(DOS_BIM, sep="\t").set_index("rsID")
    prov = _provenance()

    def setcols(rsids):                       # == script-23 _setcols
        ss = snps_for_tier(comb, rsids, "noAPOE",
                            E_SNPS_DEFAULT, APOE_CLUSTER_DEFAULT)
        return [s for s in ss if s in dos_cols]

    sets = {"6W27B": base,
            "9W27B": base + w_novel,
            "9W27B5D": base + w_novel + d_novel}

    args.out.mkdir(parents=True, exist_ok=True)
    summary = []
    for name, rsids in sets.items():
        cc = setcols(rsids)
        d = args.out / name
        d.mkdir(parents=True, exist_ok=True)

        rows = []
        for s in cc:
            a1 = str(bim.loc[s, "A1"]) if s in bim.index else ""
            a2 = str(bim.loc[s, "A2"]) if s in bim.index else ""
            pv = prov.get(s, {"gene": "", "lead_rsID": s,
                              "source": "", "origin": ""})
            rows.append({
                "rsID": s, "gene": pv["gene"], "lead_rsID": pv["lead_rsID"],
                "source": pv["source"], "origin": pv["origin"],
                "CHR": str(comb.loc[s, "CHR"]),
                "BP_GRCh38": int(comb.loc[s, "BP"]),
                "effect_allele": a1, "other_allele": a2,
                "REF": a2, "ALT": a1,
                "beta_A1": float(comb.loc[s, "beta_A1"])})
        snp = (pd.DataFrame(rows, columns=SNP_COLS)
               .sort_values(["CHR", "BP_GRCh38"], key=lambda c: c.map(
                   lambda v: (int(v) if str(v).isdigit() else 99)
                   if c.name == "CHR" else v))
               .reset_index(drop=True))
        snp.to_csv(d / f"{name}_snps.tsv", sep="\t", index=False)

        dsub = dos[["PTID"] + [c for c in dos.columns if c in cc]]
        dsub.to_csv(d / f"{name}_patient_dosage.tsv", sep="\t", index=False)
        (bim.reset_index()[bim.reset_index()["rsID"].isin(cc)]
            [["rsID", "CHR", "BP", "A1", "A2"]]
            .to_csv(d / f"{name}_patient_dosage.bim", sep="\t", index=False))

        nb = int((snp["source"].str.contains("Bellenguez", case=False)).sum())
        nw = int((snp["source"].str.contains("Wightman", case=False)).sum())
        nd = int((snp["source"].str.contains("Desikan", case=False)).sum())
        (d / "README.txt").write_text(
            f"{name} — filtered GW-significant PRS set ({len(cc)} SNPs: "
            f"{nw} Wightman + {nb} Bellenguez"
            + (f" + {nd} Desikan" if nd else "") + ", APOE region excluded "
            "[noAPOE tier]).\n\n"
            f"{name}_snps.tsv columns:\n"
            "  rsID          dbSNP id (GRCh38, RsMergeArch-canonical)\n"
            "  gene          mapped gene (Bellenguez MAPPED_GENE / Desikan "
            "locus_name)\n"
            "  lead_rsID     GWAS lead variant this SNP represents\n"
            "  source        Bellenguez2022 | Wightman2021_lead | "
            "Desikan2019_PGS\n"
            "  origin        bellenguez_in_430 | in_430 | re_extracted | "
            "desikan_novel\n"
            "  CHR,BP_GRCh38 chromosome and GRCh38 base-pair position\n"
            "  effect_allele the allele the GWAS beta refers to AND that the "
            "dosage counts (= A1 = ALT)\n"
            "  other_allele  the non-effect allele (= A2 = REF)\n"
            "  REF/ALT       reference/alternate on the refcorr GRCh38 array "
            "(REF=other_allele, ALT=effect_allele)\n"
            "  beta_A1       published GWAS log-OR effect size on "
            "effect_allele (allele-harmonised)\n\n"
            f"{name}_patient_dosage.tsv: PTID + one column per rsID, "
            "additively coded 0/1/2 = count of effect_allele (ALT). "
            "Missing genotypes are blank here; the models impute them with "
            "the cohort mean of that SNP before scoring.\n\n"
            f"{name}_patient_dosage.bim: plink-style (rsID CHR BP A1 A2); "
            "A1 = counted/effect/ALT allele, A2 = other/REF.\n\n"
            "Weighted PRS reconstruction (identical to the models'\n"
            f"betaPRS_{name}):  PRS = Σ_i beta_A1_i × dosage_i  over the\n"
            f"rsIDs in {name}_snps.tsv (cohort-mean impute missing dosage,\n"
            "then z-score on the training fold). The raw-dosage model\n"
            f"(dosage_{name}) instead uses the per-SNP dosages directly as\n"
            "covariates without the betas.\n", encoding="utf-8")

        print(f"[out] {d}  ({len(cc)} SNPs: {nw}W + {nb}B"
              + (f" + {nd}D" if nd else "") + f"; dosage {len(dsub)}×{len(cc)})")
        summary.append({"set": name, "n_snps": len(cc), "n_wightman": nw,
                        "n_bellenguez": nb, "n_desikan": nd,
                        "n_patients": len(dsub)})

    pd.DataFrame(summary).to_csv(args.out / "sets_summary.tsv", sep="\t",
                                 index=False)
    (args.out / "README.txt").write_text(
        "GWAS_filtered — training-ready filtered GW-significant PRS sets.\n\n"
        "Three nested SNP sets (each its own sub-directory), reconstructed\n"
        "byte-identically to the filtered genetic-baseline (script 20) and\n"
        "Cox survival (script 23) models via the script-23 covar_matrix\n"
        "logic (build_beta_table base ∪ enriched novel → noAPOE tier →\n"
        "genotyped in the filtered dosage matrix):\n\n"
        "  6W27B    29 SNPs  6 Wightman + 27 Bellenguez (base β table)\n"
        "  9W27B    32 SNPs  6W27B + 3 re-extracted Wightman leads\n"
        "  9W27B5D  37 SNPs  9W27B + 5 novel Desikan-2019 PGS loci\n\n"
        "Each <SET>/ has <SET>_snps.tsv (rsID, gene, lead_rsID, source,\n"
        "origin, CHR, BP_GRCh38, effect/other allele, REF/ALT, beta_A1),\n"
        "<SET>_patient_dosage.tsv (PTID + per-rsID 0/1/2 ALT dosage),\n"
        "<SET>_patient_dosage.bim, and a per-set README. sets_summary.tsv\n"
        "tallies the three. Effect allele = A1 = ALT (dosage-counted);\n"
        "PRS = Σ beta_A1 × dosage (cohort-mean impute, train z-score).\n",
        encoding="utf-8")
    print(f"[out] {args.out/'sets_summary.tsv'}  + READMEs")


if __name__ == "__main__":
    main()
