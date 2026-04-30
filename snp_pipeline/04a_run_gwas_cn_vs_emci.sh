#!/usr/bin/env bash
# run_gwas_cn_vs_emci.sh
# =======================
# PLINK 1.9 logistic regression: CN (control) vs Early MCI (case)
#
# Prerequisites:
#   1. run_gwas.sh (CN vs AD) has already been run once to confirm setup.
#   2. prepare_pheno_cn_vs_emci.py  → produces phenotype_cn_vs_emci.pheno
#
# Run from: D:/ADNI_SNP_Omni2.5M_20140220/
#   bash "C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\snp_pipeline\run_gwas_cn_vs_emci.sh"

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python}"

PLINK="./plink.exe"
BFILE="SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr"
COVAR="covariates.cov"
OUT_DIR="results/gwas_no_ld_CN_EMCI"
OUT="${OUT_DIR}/gwas_CN_vs_EMCI_GRCh38"

mkdir -p "$OUT_DIR"

# ── Step 1: generate phenotype file ──────────────────────────────────────────
echo "=== Generating CN vs EMCI phenotype file ==="
"$PYTHON" -u "${SCRIPT_DIR}/prepare_pheno_cn_vs_emci.py"

PHENO="phenotype_cn_vs_emci.pheno"

# ── Step 2: logistic regression ───────────────────────────────────────────────
# Phenotype: 1 = CN (control, n=211), 2 = EarlyMCI (case, n=222), -9 = excluded
# Covariates: AGE, SEX, APOE4 (same as CN vs AD run)
echo ""
echo "=== Logistic regression: CN vs Early MCI ==="
"$PLINK" \
    --bfile  "$BFILE"  \
    --pheno  "$PHENO"  \
    --covar  "$COVAR"  \
    --covar-name AGE,SEX,APOE4 \
    --logistic hide-covar \
    --ci 0.95 \
    --out "$OUT"

echo ""
echo "=== Done ==="
echo "Results : ${OUT}.assoc.logistic"
echo "Log     : ${OUT}.log"
echo ""

# ── Step 3: significance report + BMFM labels ────────────────────────────────
echo "=== Running significance analysis ==="
"$PYTHON" -u "${SCRIPT_DIR}/parse_gwas_results.py" \
    --assoc   "${OUT}.assoc.logistic" \
    --out-dir "$OUT_DIR" \
    --alpha 0.05 \
    --fdr-threshold 0.05 \
    --null-p-min 0.5

echo ""
echo "=== Output files ==="
echo "  Significance report : ${OUT_DIR}/gwas_significance_report.txt"
echo "  Significant SNPs    : ${OUT_DIR}/gwas_significant_snps.tsv"
echo "  BMFM labels         : ${OUT_DIR}/gwas_bmfm_labels.tsv"
