#!/usr/bin/env bash
# qc.sh
# =====
# Full quality-control pipeline for ADNI Omni2.5M SNP data.
# All steps have already been run. This script is for transparency
# and reproducibility documentation only — do NOT re-run.
#
# Source dataset : WGS_Omni25_BIN_wo_ConsentsIssues.bed/.bim/.fam
# Final QC output: SNP_filtered_hwe.bed/.bim/.fam  (812 subjects)
# MRI-matched    : SNP_filtered_with_mri.bed/.bim/.fam (616 subjects)
#
# Run from: D:/ADNI_SNP_Omni2.5M_20140220/
# Setup: export PATH="/d/ADNI_SNP_Omni2.5M_20140220:$PATH"  # plink.exe on PATH

PLINK="./plink.exe"
BFILE="WGS_Omni25_BIN_wo_ConsentsIssues"

# ══════════════════════════════════════════════════════════════════════════════
# 0. INSPECT INPUT FILES
# ══════════════════════════════════════════════════════════════════════════════

# ── FAM file ──────────────────────────────────────────────────────────────────
# Columns: FID  IID  PAT  MAT  SEX(1=M,2=F,0=?)  PHENO(-9=missing)
#
# Check phenotype distribution (all missing at input — will be provided separately):
#   awk '{print $6}' ${BFILE}.fam | sort -n | uniq -c
#   Result:
#       812 -9
#
# Check sex distribution:
#   awk '{print $5}' ${BFILE}.fam | sort -n | uniq -c
#   Result:
#       448 1   (male)
#       364 2   (female)
#
# Check for subjects with unknown sex (code 0):
#   awk '$5==0 {count++} END{print count+0}' ${BFILE}.fam
#   Result: 0  (all subjects have sex recorded)

# ── BIM file ──────────────────────────────────────────────────────────────────
# Columns: CHR  SNP  cM  BP  A1  A2
# Chr 0  = unplaced/unmapped variants
# Chr 23 = X,  24 = Y,  25 = XY pseudoautosomal,  26 = mitochondrial
#
# Check cM and BP range across all chromosomes:
#   awk 'BEGIN{min=1e99;max=-1e99}{if($3<min)min=$3;if($3>max)max=$3}END{print "cM min="min,"max="max}' ${BFILE}.bim
#   Result: cM min=0  max=280.497
#
#   awk 'BEGIN{min=1e99;max=-1e99}{if($4<min)min=$4;if($4>max)max=$4}END{print "BP min="min,"max="max}' ${BFILE}.bim
#   Result: BP min=0  max=249213967
#
# Count unplaced SNPs (chr 0):
#   awk '$1==0' ${BFILE}.bim | wc -l
#   Result: 7238 unplaced/unmapped SNPs

# ══════════════════════════════════════════════════════════════════════════════
# 1. MINOR ALLELE FREQUENCY (MAF)
# ══════════════════════════════════════════════════════════════════════════════
# Calculate MAF for all SNPs (report only, no filtering):
"$PLINK" --bfile "$BFILE" --freq --out maf

# head maf.frq
#  CHR              SNP   A1   A2          MAF  NCHROBS
#    0        rs2698846    0    G            0     1620
#    0        rs2542903    0    T            0     1600
#    0        rs1059667    0    T            0     1612
#    0       rs35144699    T    A        0.327     1624
#    0       rs35425675    A    T       0.2715     1584
#    0       rs12065999    T    A     0.004353     1608
#
# MAF = 0 means the SNP is monomorphic in this sample (everyone same allele).
# NCHROBS = number of non-missing chromosomes observed for the SNP.
# MAF 0 + high NCHROBS → truly monomorphic; MAF 0 + low NCHROBS → caution.
# These SNPs will be removed by --maf 0.01 in the combined filter step below.

# ══════════════════════════════════════════════════════════════════════════════
# 2. MISSINGNESS
# ══════════════════════════════════════════════════════════════════════════════
mkdir -p results/missingness
"$PLINK" --bfile "$BFILE" --missing --out results/missingness/missing

# Heterozygous haploid genotype calls (.hh file):
#   head missing.hh
#   5   109_S_2237  kgp22761656
#   11  094_S_2216  kgp22761656
#   ...
# These are calls on haploid chromosomes (chrX for males outside PAR) where
# PLINK expects one allele but finds a het call. Most operations treat them
# as missing. Source: PAR coding or sex-coding inconsistency.

# Missingness by individual (.imiss):
#   head missing.imiss
#    FID          IID MISS_PHENO   N_MISS   N_GENO   F_MISS
#      1   067_S_0056          Y     5211  2377294 0.002192
#      2   123_S_0108          Y     6962  2379855 0.002925
#      ...
#
# Count individuals with >2% missingness:
#   awk 'NR>1 && $6>0.02' missing.imiss | wc -l
#   Result: 0  (no subjects exceed 2% missingness threshold)

# Missingness by SNP (.lmiss):
#   head missing.lmiss
#    CHR         SNP   N_MISS   N_GENO   F_MISS
#      0   rs2698846        2      812 0.002463
#      0   rs3876241      365      812   0.4495    ← very high missingness
#      0    rs160179      383      812   0.4717
#      ...
#
# Count SNPs with >5% missingness:
#   awk 'NR>1 && $5>0.05' missing.lmiss | wc -l
#   Result: 20859 SNPs above 5% missingness threshold → removed by --geno 0.02

# ══════════════════════════════════════════════════════════════════════════════
# 3. COMBINED QC FILTER: autosomes + missingness
# ══════════════════════════════════════════════════════════════════════════════
# Remove:
#   - Non-autosomal SNPs (chr 0, X, Y, XY, MT) via --autosome
#   - Individuals with >2% missing genotypes via --mind 0.02
#   - SNPs with >2% missing calls via --geno 0.02
"$PLINK" --bfile "$BFILE" \
    --autosome \
    --mind 0.02 \
    --geno 0.02 \
    --make-bed \
    --out SNP_filtered

# ══════════════════════════════════════════════════════════════════════════════
# 4. HARDY-WEINBERG EQUILIBRIUM (HWE) — pre-filter diagnostic
# ══════════════════════════════════════════════════════════════════════════════
mkdir -p results/hwe
"$PLINK" --bfile "$BFILE" --hardy --out results/hwe/hwe

# Warnings at this stage:
#   Warning: 417855 het. haploid genotypes present (see results/hwe.hh)
#   Warning: Nonmissing nonmale Y chromosome genotype(s) present
#   Total genotyping rate: 0.99744
#
# head results/hwe/hwe.hwe
#  CHR         SNP     TEST   A1   A2              GENO   O(HET)   E(HET)          P
#    0   rs2698846  ALL(NP)    0    G           0/0/810        0        0          1
#    0   rs35425675 ALL(NP)    A    T           0/430/362  0.5429   0.3955   4.33e-36
#    0    rs160179  ALL(NP)    C    G         202/2/225  0.004662   0.4986  1.734e-124
#
# Worst HWE violations (most significant):
#   sort -k9,9g results/hwe/hwe.hwe | head
#   CHR         SNP     TEST   A1  A2              GENO   O(HET)   E(HET)           P
#    22  kgp3432397  ALL(NP)    T   A        280/35/482  0.04391   0.4679  1.446e-168
#     6   rs9259192  ALL(NP)    C   T         197/7/468  0.01042   0.4187  1.215e-162
#     7 kgp13298227  ALL(NP)    T   C         167/6/566 0.008119   0.3542  8.973e-160
#
# Count SNPs below HWE threshold of 1e-7:
#   awk 'NR>1 && $9 < 1e-7 {count++} END{print count+0}' results/hwe/hwe.hwe
#   Result: 5642 (0.24% of all SNPs)

# ── HWE diagnostic on the post-filter dataset ─────────────────────────────────
"$PLINK" --bfile SNP_filtered --hardy --out results/hwe/hwe_filtered

# wc -l results/hwe/hwe_filtered.hwe → 2255863 SNPs remain after combined filter
#
# Count SNPs below HWE threshold of 1e-7 in filtered dataset:
#   awk 'NR>1 && $9 < 1e-7 {count++} END{print count+0}' results/hwe/hwe_filtered.hwe
#   Result: 3389 out of 2255862 = 0.15%

# ══════════════════════════════════════════════════════════════════════════════
# 5. HWE FILTER
# ══════════════════════════════════════════════════════════════════════════════
# Exclude SNPs with HWE exact-test p < 1e-7.
# Threshold chosen to remove clear genotyping errors while retaining
# SNPs under genuine selection (stricter than 1e-6 but not overly aggressive).
"$PLINK" --bfile SNP_filtered \
    --hwe 1e-7 \
    --make-bed \
    --out SNP_filtered_hwe

# Output: SNP_filtered_hwe.bed/.bim/.fam  (812 subjects)

# ══════════════════════════════════════════════════════════════════════════════
# 6. RESTRICT TO MRI-MATCHED SUBJECTS (616)
# ══════════════════════════════════════════════════════════════════════════════
# SNP_filtered_with_mri.fam was produced by prepare_plink_files.py,
# which joins SNP_filtered_hwe.fam against subjects_with_snp_and_mri.tsv.
# 812 total - 196 without MRI = 616 subjects retained.
#
# python "C:\...\snp_pipeline\prepare_plink_files.py"
#
# Produces (in D:/ADNI_SNP_Omni2.5M_20140220/):
#   SNP_filtered_with_mri.fam   (616 subjects)
#   phenotype.pheno             (FID IID PHENO; CN=1, AD/LMCI=2, EMCI=-9)
#   covariates.cov              (FID IID AGE SEX APOE4)
#
# head phenotype.pheno
#   FID     IID           PHENO
#   1       067_S_0056    1
#   3       067_S_0059    1
#   6       131_S_0123    1
#   10      073_S_2264    -9
#   15      023_S_0331    2
#
# head covariates.cov
#   FID     IID           AGE    SEX   APOE4
#   1       067_S_0056    69.6   2     0.0
#   3       067_S_0059    70.9   2     0.0
#   6       131_S_0123    73.3   1     0.0
#   10      073_S_2264    69.3   1     1.0

# Use SNP_filtered_with_mri.fam as a --keep list to subset the binary dataset:
"$PLINK" --bfile SNP_filtered_hwe \
    --keep SNP_filtered_with_mri.fam \
    --make-bed \
    --out SNP_filtered_with_mri

# Output: SNP_filtered_with_mri.bed/.bim/.fam  (616 subjects, same SNP set)
