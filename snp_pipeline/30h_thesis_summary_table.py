r"""
30h_thesis_summary_table.py   (LOCAL, env: snp)
==============================================
Render a thesis-grade PNG table summarising the AD-PRS source roster.

One row per source (18 total: 16 PRS-usable + 2 audit-only). Columns:
  study           — source label (e.g. "Bellenguez 2022")
  source_kind     — "GWAS catalog" / "PGS catalog" / "xlsx supplementary"
                    / "GWAS sumstats" / "GWAS sumstats + xlsx"
  pheno_class     — clinical_AD / clinical_AD_MTAG / clinical_AAO /
                    pathology_microglia / pathology_amyloid /
                    pathology_plaque / pathology_mixed
  pub_count       — n published lead-SNPs in the source
  before_QC_count — n on the ADNI Omni2.5M chip BEFORE QC
                    (in_unfiltered_bim==True, any direct-rsID /
                    rsID-merge / Illumina kgp manifest / GRCh37 position)
  after_QC_count  — n on the MAF BIM AFTER QC
                    (in_maf_bim==True, post --geno 0.02 / --hwe 1e-7 /
                    --maf 0.01)

QC threshold footer:
  --geno 0.02 (per-SNP missingness ≤ 2%) | --hwe 1e-7 (HWE p ≥ 1e-7) |
  --maf 0.01 (MAF ≥ 1%) | --mind 0.02 (per-sample missingness) ; on
  GRCh38 after ref-correction via scripts 05a/05b.

Output:
  D:\…\GWAS_comprehensive_v2\thesis_source_summary.png       (300 DPI)
  D:\…\GWAS_comprehensive_v2\thesis_source_summary.tsv

Usage:
  conda run -n snp python snp_pipeline/30h_thesis_summary_table.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
RECON = BASE / "source_prs" / "unfiltered_SNP_reconciliation"
OUT = BASE / "GWAS_comprehensive_v2"

# Source roster (same order as in 30b's SOURCES dict, plus audit)
SOURCE_ORDER = [
    # study_label, source_key, source_kind, citation_year
    ("Bellenguez 2022",        "Bellenguez",        "GWAS catalog",       2022),
    ("Wightman 2021",          "Wightman",          "GWAS sumstats",      2021),
    ("Kunkle 2019",            "Kunkle",            "GWAS catalog",       2019),
    ("Desikan 2019",           "Desikan",           "PGS catalog",        2019),
    ("Lambert 2013 (IGAP)",    "Lambert",           "GWAS catalog",       2013),
    ("DeRojas 2021",           "DeRojas",           "PGS catalog",        2021),
    ("Schwanzentruber 2021",   "Schwanzentruber",   "xlsx supplementary", 2021),
    ("Najar 2021",             "Najar",             "PGS catalog",        2021),
    ("Ebenau 2021",            "Ebenau",            "PGS catalog",        2021),
    ("Leonenko 2019",          "Leonenko",          "PGS catalog",        2019),
    ("Vasiljevic 2023",        "Vesilievick",       "PGS catalog",        2023),
    ("Zhang 2020",             "Zhang",             "PGS catalog",        2020),
    ("Felsky 2019 (MF) *",     "Felsky_MF",         "xlsx supplementary", 2019),
    ("Felsky 2019 (IT) *",     "Felsky_IT",         "xlsx supplementary", 2019),
    ("O'Neil 2025 (NPY) †",    "ONeil_NPY",         "xlsx supplementary", 2025),
    ("O'Neil 2025 (GHR) †",    "ONeil_GHR",         "xlsx supplementary", 2025),
    ("Yang 2023 (Mic) ‡",      "Yang_Mic",          "PRSet + GWAS sumstats", 2023),
    ("Yang 2023 (Ast) ‡",      "Yang_Ast",          "PRSet + GWAS sumstats", 2023),
    ("Yang 2023 (Ex) ‡",       "Yang_Ex",           "PRSet + GWAS sumstats", 2023),
    ("Yang 2023 (In) ‡",       "Yang_In",           "PRSet + GWAS sumstats", 2023),
    ("Yang 2023 (Oli) ‡",      "Yang_Oli",          "PRSet + GWAS sumstats", 2023),
    ("Yang 2023 (Opc) ‡",      "Yang_Opc",          "PRSet + GWAS sumstats", 2023),
    ("Kosteridis 2024 (novel)",       "Kosteridis_novel_AD",     "xlsx supplementary", 2024),
    ("Kosteridis 2024 (shared)",      "Kosteridis_shared_AD_CV", "xlsx + GWAS sumstats", 2024),
    # audit-only
    ("Huang 2017 [audit-only]",        "Huang",                  "xlsx supplementary", 2017),
    ("O'Neil 2025 SST [audit-only]",   "ONeil_SST_candidates",   "xlsx supplementary", 2025),
]

AUDIT_KEYS = {"Huang", "ONeil_SST_candidates"}

PHENO_COLOUR = {
    "clinical_AD":                 "#cfe2ef",  # light blue
    "clinical_AD_MTAG":            "#a8c8e0",  # medium blue
    "clinical_AAO":                "#7fa9c8",  # darker blue
    # pathology palette: orange→red gradient by cell type
    "pathology_microglia":         "#ffd9a8",  # light orange (Felsky + Yang_Mic)
    "pathology_astrocyte":         "#ffc680",  # peach (Yang_Ast)
    "pathology_neuron_excitatory": "#ffb866",  # medium orange (Yang_Ex)
    "pathology_neuron_inhibitory": "#ff9e4d",  # warm orange (Yang_In)
    "pathology_oligodendrocyte":   "#ff8533",  # deeper orange (Yang_Oli)
    "pathology_opc":               "#ff7a1f",  # red-orange (Yang_Opc)
    "pathology_amyloid":           "#ffc080",  # ONeil_GHR
    "pathology_plaque":            "#ffb066",  # ONeil_NPY
    "pathology_mixed":             "#ffa54d",  # ONeil_SST candidates
    "":                            "#f5f5f5",
}


def _load(src: str) -> pd.DataFrame:
    p = RECON / f"{src}_unfiltered_reconciliation.tsv"
    if not p.exists():
        print(f"  [WARN] missing reconciliation TSV: {p}")
        return pd.DataFrame()
    return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)


def _truthy(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin({"true", "1", "yes"})


def build_summary() -> pd.DataFrame:
    rows = []
    for label, key, kind, year in SOURCE_ORDER:
        df = _load(key)
        if df.empty:
            rows.append({"study": label, "source_kind": kind,
                          "pheno_class": "",
                          "pub_count": 0, "before_QC_count": 0,
                          "after_QC_count": 0, "audit_only": key in AUDIT_KEYS})
            continue
        pheno = df["pheno_class"].iloc[0] if "pheno_class" in df.columns else ""
        n_pub = len(df)
        n_chip = int(_truthy(df["in_unfiltered_bim"]).sum())
        n_maf = int(_truthy(df["in_maf_bim"]).sum())
        rows.append({
            "study": label,
            "source_kind": kind,
            "pheno_class": pheno,
            "pub_count": n_pub,
            "before_QC_count": n_chip,
            "after_QC_count": n_maf,
            "audit_only": key in AUDIT_KEYS,
        })
    return pd.DataFrame(rows)


def render_png(df: pd.DataFrame, png_path: Path) -> None:
    # Build the display DataFrame (drop audit_only column from the rendered table)
    disp_cols = ["study", "source_kind", "pheno_class",
                 "pub_count", "before_QC_count", "after_QC_count"]
    disp = df[disp_cols].copy()

    # Figure sizing: tall thin portrait
    n_rows = len(disp)
    fig_h = 0.45 * (n_rows + 4)  # +4 for header + footer
    fig_w = 12.5
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=300)
    ax.axis("off")

    # Table cell text
    cell_text = []
    cell_colours = []
    for _, r in disp.iterrows():
        is_audit = bool(df.loc[r.name, "audit_only"])
        cell_text.append([
            r["study"], r["source_kind"], r["pheno_class"],
            str(r["pub_count"]),
            str(r["before_QC_count"]),
            str(r["after_QC_count"]),
        ])
        base = "#f0f0f0" if is_audit else PHENO_COLOUR.get(r["pheno_class"], "#ffffff")
        cell_colours.append([base] * 6)

    col_labels = ["Study", "Source format", "Phenotype class",
                  "Pub.\nleads", "Before QC\n(on chip)",
                  "After QC\n(on MAF)"]
    col_widths = [0.26, 0.18, 0.20, 0.09, 0.13, 0.14]

    tab = ax.table(
        cellText=cell_text,
        cellColours=cell_colours,
        colLabels=col_labels,
        colWidths=col_widths,
        loc="center",
        cellLoc="left",
    )
    tab.auto_set_font_size(False)
    tab.set_fontsize(9)
    tab.scale(1, 1.3)

    # Header styling
    for j in range(len(col_labels)):
        h = tab[0, j]
        h.set_facecolor("#3a3a3a")
        h.set_text_props(color="white", fontweight="bold",
                          ha="center", va="center")
        h.set_height(0.07)

    # Right-align numeric columns
    for i in range(1, n_rows + 1):
        for j in (3, 4, 5):
            tab[i, j].set_text_props(ha="right", va="center")
        for j in (0, 1, 2):
            tab[i, j].set_text_props(ha="left", va="center")

    # Title
    ax.set_title(
        "AD-PRS source roster — ADNI Omni2.5M cohort, QC-pre/post counts",
        fontsize=14, fontweight="bold", pad=12,
    )

    # Footer with QC thresholds + significance caveats
    footer = (
        "QC thresholds applied (PLINK 1.9): "
        "--mind 0.02 (≤2% per-sample missingness), "
        "--geno 0.02 (≤2% per-SNP missingness), "
        "--hwe 1e-7 (HWE p ≥ 1e-7 in controls), "
        "--maf 0.01 (MAF ≥ 1%).\n"
        "Resolution chain: direct rsID → RsMergeArch canonical → "
        "Illumina kgp probe (b138/b144/b150/b151) → GRCh37 position. "
        "MAF BIM is GRCh38 after liftover + ref-correction.\n"
        "GWAS significance: most sources report leads at the conventional "
        "P < 5×10⁻⁸ threshold (Bellenguez, Wightman, Kunkle, Lambert, "
        "DeRojas, Schwanzentruber, Najar, Zhang, Kosteridis novel + "
        "shared Supp1_MTAG entries) or are curated PGS-Catalog scoring "
        "files (Desikan, Ebenau, Leonenko, Vasiljevic).\n"
        "*  Felsky 2019 (MF + IT) reports leads at the looser suggestive "
        "threshold P < 1×10⁻⁵. The two true GW-significant Felsky leads "
        "(rs2997325 MF P=1.88×10⁻⁸, rs183093970 IT P=5.47×10⁻¹⁰) are "
        "ABSENT from the Omni2.5M chip; resolved Felsky hits are "
        "sub-GW-sig (P ≈ 3-9×10⁻⁶).\n"
        "†  O'Neil 2025 NPY and GHR lead SNPs are targeted gene-region "
        "marginal hits (P=8×10⁻⁵ and 3×10⁻⁴ respectively), not GW-sig. "
        "Huang 2017 leads come from a 22-SNP AAOS LASSO model; no "
        "per-SNP β is reported.\n"
        "‡  Yang 2023 cell-type-partitioned PRS: Bellenguez 2022 GWAS SNPs "
        "filtered to ROSMAP single-nucleus expression markers (Mathys et "
        "al. 2019). Lead set = PRSet `<celltype>==1 AND P<5×10⁻⁸`; β/SE "
        "from raw Bellenguez 2022 sumstats (GCST90027158). Cell types: "
        "Mic (microglia), Ast (astrocytes), Ex/In (excitatory/inhibitory "
        "neurons), Oli (oligodendrocytes), Opc (OPCs).\n"
        "Phenotype classes: clinical_AD (case-control diagnosis) / "
        "clinical_AD_MTAG (MTAG-boosted with cardiovascular co-traits) / "
        "clinical_AAO (age-at-onset) / "
        "pathology_microglia (post-mortem PAM) / "
        "pathology_amyloid (amyloid-β density) / "
        "pathology_plaque (neuritic plaque burden) / "
        "pathology_mixed (audit-only candidate set).\n"
        "[audit-only] sources are tracked for cross-source overlap "
        "but excluded from PRS bundle assembly."
    )
    fig.text(0.04, 0.005, footer, fontsize=8, ha="left", va="bottom",
              wrap=True)

    plt.tight_layout(rect=[0, 0.08, 1, 1])
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [out] {png_path}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    df = build_summary()
    print(f"[summary] {len(df)} sources")
    print(df[["study", "pheno_class", "pub_count",
              "before_QC_count", "after_QC_count"]].to_string(index=False))
    df.to_csv(args.out / "thesis_source_summary.tsv", sep="\t", index=False)
    print(f"\n  [out] {args.out/'thesis_source_summary.tsv'}")
    render_png(df, args.out / "thesis_source_summary.png")


if __name__ == "__main__":
    main()
