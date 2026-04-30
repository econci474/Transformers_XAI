"""
prepare_plink_files.py
======================
Build PLINK-compatible phenotype and covariate files by joining:
  - WGS_Omni25_BIN_wo_ConsentsIssues.fam  (in D:/ADNI_SNP_Omni2.5M_20140220/)
  - subjects_with_snp_and_mri.tsv          (in D:/ADNI_BIDS_project/bids/genotype/)

Outputs (written to the same directory as this script, i.e. snp_pipeline/):
  - phenotype.pheno        FID  IID  PHENO          (616 subjects — TSV-matched only)
  - covariates.cov         FID  IID  AGE  SEX  APOE4 (616 subjects — TSV-matched only)

Also writes to D:/ADNI_SNP_Omni2.5M_20140220/:
  - SNP_filtered_hwe_excl.fam   copy of SNP_filtered_hwe.fam with PHENO column
                                 set to -9 for any subject absent from the TSV
                                 (812 total in FAM → 616 kept, 196 excluded)

Sex verification:
  The script will EXIT with an error if any subject is missing sex in either
  the FAM file (col 5 = 0) or the TSV (sex column is null/empty).
  This ensures --allow-no-sex is NOT needed in run_gwas.sh.

PLINK phenotype encoding (binary / case-control):
  1  = control  (CognitivelyNormal)               n=211
  2  = case     (LateMCI + AlzheimersDisease)      n=147+36=183
  -9 = missing  (EarlyMCI, or no match in TSV)    n=222

Run from D:/ADNI_SNP_Omni2.5M_20140220/ as:
  python snp_pipeline/prepare_plink_files.py
"""

import sys
import pathlib
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent          # snp_pipeline/
DATA_DIR   = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220")   # run from here

FAM_FILE   = DATA_DIR / "SNP_filtered_hwe.fam"
TSV_FILE   = pathlib.Path("D:/ADNI_BIDS_project/bids/genotype/subjects_with_snp_and_mri.tsv")
PHENO_OUT  = DATA_DIR / "phenotype.pheno"
COVAR_OUT  = DATA_DIR / "covariates.cov"
FAM_OUT     = DATA_DIR / "SNP_filtered_hwe_excl.fam"    # FAM with non-TSV subjects set to -9
FAM_MRI_OUT = DATA_DIR / "SNP_filtered_with_mri.fam"   # FAM with non-TSV subjects removed (616)

# ── Phenotype label mapping ───────────────────────────────────────────────────
# Verify these against your unique diagnosis_bl values (printed below).
CONTROL_LABELS = {"CognitivelyNormal"}                  # PHENO = 1
CASE_LABELS    = {"LateMCI", "AlzheimersDisease"}        # PHENO = 2
# EarlyMCI → -9 (excluded); provides better case/control balance (183 vs 211)
# CN: 211
# EMCI: 222
# LMCI: 147
# AD: 36
# Total: 616
# ── Load FAM (space-delimited, no header) ─────────────────────────────────────
print(f"Reading FAM : {FAM_FILE}")
fam = pd.read_csv(
    FAM_FILE, sep=r"\s+", header=None,
    names=["FID", "IID", "PAT", "MAT", "SEX_FAM", "PHENO_FAM"]
)
print(f"  {len(fam)} subjects in FAM file")

# ── Load TSV ──────────────────────────────────────────────────────────────────
print(f"Reading TSV : {TSV_FILE}")
tsv = pd.read_csv(TSV_FILE, sep="\t")
print(f"  {len(tsv)} rows in TSV")
print(f"  Columns: {list(tsv.columns)}")

# Show all unique baseline diagnoses so you can verify the label mapping above
print("\n  Unique diagnosis_bl values:")
for d in sorted(tsv["diagnosis_bl"].dropna().unique()):
    print(f"    {d!r}  (n={tsv['diagnosis_bl'].eq(d).sum()})")

# Rename join key to match FAM IID
tsv = tsv.rename(columns={"adni_subject_id": "IID"})

# ── Merge on IID ──────────────────────────────────────────────────────────────
df = fam[["FID", "IID"]].merge(
    tsv[["IID", "diagnosis_bl", "sex", "age", "apoe4_dosage"]],
    on="IID", how="left"
)

n_unmatched = df["diagnosis_bl"].isna().sum()
if n_unmatched > 0:
    print(f"\nWARNING: {n_unmatched} FAM subject(s) had NO match in the TSV — "
          f"they will be set to PHENO=-9 in {FAM_OUT.name}")
else:
    print("\nAll FAM subjects matched successfully in TSV.")

# ── Write exclusion FAM: subjects absent from TSV get PHENO = -9 ───────────────
print(f"\n── FAM exclusion step ────────────────────────────────────────────────")
tsv_iids = set(tsv["IID"])
fam_excl = fam.copy()
# Subjects absent from the TSV → PHENO column set to -9 (PLINK missing)
fam_excl.loc[~fam_excl["IID"].isin(tsv_iids), "PHENO_FAM"] = -9
n_excluded = (fam_excl["PHENO_FAM"] == -9).sum()
n_kept     = len(fam_excl) - n_excluded
print(f"  Subjects in FAM          : {len(fam_excl)}")
print(f"  Subjects in TSV          : {len(tsv)}")
print(f"  Kept (present in TSV)    : {n_kept}")
print(f"  Excluded (PHENO → -9)   : {n_excluded}")
fam_excl.to_csv(FAM_OUT, sep=" ", index=False, header=False)
print(f"  Exclusion FAM written    : {FAM_OUT}")

# Drop excluded subjects entirely → 616-subject FAM
fam_mri = fam[fam["IID"].isin(tsv_iids)].copy()
fam_mri.to_csv(FAM_MRI_OUT, sep=" ", index=False, header=False)
print(f"  MRI-only FAM written     : {FAM_MRI_OUT}  ({len(fam_mri)} subjects)")
print("──────────────────────────────────────────────────────────────────────")

# ── Sex verification ──────────────────────────────────────────────────────────
# 1. FAM sex column (col 5): must be 1 (male) or 2 (female), never 0 (unknown)
print("\n── Sex verification ──────────────────────────────────────────────")
fam_sex_missing = fam[fam["SEX_FAM"] == 0]
if not fam_sex_missing.empty:
    print(f"ERROR: {len(fam_sex_missing)} subject(s) have sex=0 (unknown) in the FAM file:")
    print(fam_sex_missing[["FID", "IID"]].to_string(index=False))
    sys.exit(1)
print(f"  FAM sex  : OK — all {len(fam)} subjects have sex 1 or 2")

# 2. TSV sex column: must not be null
tsv_sex_missing = tsv[tsv["sex"].isna() | (tsv["sex"].astype(str).str.strip() == "")]
if not tsv_sex_missing.empty:
    print(f"ERROR: {len(tsv_sex_missing)} TSV row(s) are missing a sex value:")
    print(tsv_sex_missing[["IID", "sex"]].to_string(index=False))
    sys.exit(1)
print(f"  TSV sex  : OK — all {len(tsv)} TSV subjects have a sex value")

# 3. Cross-check: sex in FAM vs sex in TSV for matched subjects
df_sex_check = df[["FID", "IID"]].merge(fam[["IID", "SEX_FAM"]], on="IID", how="left")
df_sex_check["SEX_TSV"] = df["sex"].apply(
    lambda s: 1 if str(s).strip().upper() == "M" else (2 if str(s).strip().upper() == "F" else -9)
)
df_sex_check = df_sex_check[df_sex_check["SEX_TSV"] != -9]  # only matched subjects
sex_mismatch = df_sex_check[df_sex_check["SEX_FAM"] != df_sex_check["SEX_TSV"]]
if not sex_mismatch.empty:
    print(f"WARNING: {len(sex_mismatch)} subject(s) have mismatched sex between FAM and TSV:")
    print(sex_mismatch.to_string(index=False))
    print("  → The TSV sex values will be used in the covariate file.")
else:
    print(f"  Cross-check: OK — FAM and TSV sex values agree for all matched subjects")
print("─────────────────────────────────────────────────────────────────")

# ── Encode phenotype ──────────────────────────────────────────────────────────
def encode_pheno(diag):
    if pd.isna(diag):
        return -9
    if diag in CASE_LABELS:
        return 2
    if diag in CONTROL_LABELS:
        return 1
    return -9  # MCI variants and any unrecognised labels → missing

# ── Filter df to TSV-matched subjects only (aligns with SNP_filtered_with_mri.fam) ──
# df was built with how="left" so it contains all 812 FAM subjects; restrict to
# the 616 that are present in the TSV before writing phenotype and covariate files.
df = df[df["IID"].isin(tsv_iids)].reset_index(drop=True)
print(f"\nRestricted to {len(df)} TSV-matched subjects for phenotype/covariate outputs.")

df["PHENO"] = df["diagnosis_bl"].apply(encode_pheno)

print("\nPhenotype distribution:")
counts = df["PHENO"].value_counts().sort_index()
labels = {1: "CN / control", 2: "AD / case", -9: "missing / excluded"}
for val, cnt in counts.items():
    print(f"  PHENO={val:>3}  {labels.get(val, '?'):20s}  n={cnt}")

# ── Encode covariates ─────────────────────────────────────────────────────────
# Sex: PLINK convention  1 = male, 2 = female  (matches FAM col 5)
def encode_sex(s):
    if pd.isna(s):
        return -9
    return 1 if str(s).strip().upper() == "M" else 2

df["SEX_COV"] = df["sex"].apply(encode_sex)

# Age: use as-is (continuous); -9 for missing
df["AGE"] = df["age"].where(df["age"].notna(), -9)

# APOE4 dosage: 0, 1, or 2 copies of ε4 allele (continuous additive coding)
df["APOE4"] = df["apoe4_dosage"].where(df["apoe4_dosage"].notna(), -9)

# ── Write phenotype file ───────────────────────────────────────────────────────
pheno_df = df[["FID", "IID", "PHENO"]]
pheno_df.to_csv(PHENO_OUT, sep="\t", index=False, na_rep="-9")
print(f"\nPhenotype file written : {PHENO_OUT}")

# ── Write covariate file ───────────────────────────────────────────────────────
covar_df = df[["FID", "IID", "AGE", "SEX_COV", "APOE4"]].rename(
    columns={"SEX_COV": "SEX"}
)
covar_df.to_csv(COVAR_OUT, sep="\t", index=False, na_rep="-9")
print(f"Covariate file written : {COVAR_OUT}")

print("\nCovariate summary:")
print(covar_df[["AGE", "SEX", "APOE4"]].describe().round(2).to_string())

print("\nDone. Next step: run snp_pipeline/run_gwas.sh from D:/ADNI_SNP_Omni2.5M_20140220/")
