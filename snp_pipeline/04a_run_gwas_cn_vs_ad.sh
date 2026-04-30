#!/usr/bin/env bash
# run_gwas.sh
# ===========
# PLINK 1.9 logistic regression: CN (control) vs AD+LMCI (case)
#
# Prerequisites — run these first:
#   1. qc.sh             → produces SNP_filtered_with_mri.bed/.bim/.fam
#   2. prepare_plink_files.py → produces phenotype.pheno, covariates.cov
#
# Run from: D:/ADNI_SNP_Omni2.5M_20140220/
#   bash "C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\snp_pipeline\run_gwas.sh"

set -euo pipefail

PLINK="./plink.exe"
BFILE="SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr"
PHENO="phenotype.pheno"
COVAR="covariates.cov"
OUT_DIR="results/gwas"
OUT="${OUT_DIR}/gwas_CN_vs_AD_GRCh38"

mkdir -p "$OUT_DIR"

# ── Logistic regression ───────────────────────────────────────────────────────
# Phenotype encoding : 1 = CN (control), 2 = AD/LMCI (case), -9 = EMCI (excluded)
# Covariates         : AGE (continuous), SEX (1=M/2=F), APOE4 dosage (0/1/2)
# --logistic         : binary case/control outcome
# hide-covar         : only SNP association rows in output (no covariate rows)
# --ci 0.95          : report 95% confidence intervals on the odds ratio
echo "=== Logistic regression: CN vs AD/LMCI ==="
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

# ── Post-processing: significance report + BMFM labels ───────────────────────
# Requires: python parse_gwas_results.py  (same directory as this script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python}"

echo "=== Running significance analysis ==="
"$PYTHON" -u "${SCRIPT_DIR}/parse_gwas_results.py" \
    --assoc  "${OUT}.assoc.logistic" \
    --out-dir "$OUT_DIR" \
    --alpha 0.05 \
    --fdr-threshold 0.05 \
    --null-p-min 0.5

echo ""
echo "=== Output files ==="
echo "  Significance report : ${OUT_DIR}/gwas_significance_report.txt"
echo "  Significant SNPs    : ${OUT_DIR}/gwas_significant_snps.tsv"
echo "  BMFM labels         : ${OUT_DIR}/gwas_bmfm_labels.tsv"
