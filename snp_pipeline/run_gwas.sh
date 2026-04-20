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
BFILE="SNP_filtered_with_mri_ld_pruned"
PHENO="phenotype.pheno"
COVAR="covariates.cov"
OUT_DIR="results/gwas"
OUT="${OUT_DIR}/gwas_CN_vs_AD"

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
echo "Genome-wide significant hits (p < 5e-8):"
echo "  awk 'NR==1 || \$9 < 5e-8' ${OUT}.assoc.logistic | head -20"
