"""
prepare_pheno_cn_vs_emci.py
============================
Generate a PLINK phenotype file for the CN vs Early MCI comparison.

Encoding
--------
  1  = CognitivelyNormal  (control)   n=211
  2  = EarlyMCI           (case)      n=222
  -9 = LateMCI / AlzheimersDisease    (excluded)

Reads the same source TSV used by prepare_plink_files.py.
Writes phenotype_cn_vs_emci.pheno to D:/ADNI_SNP_Omni2.5M_20140220/.

Usage
-----
  python prepare_pheno_cn_vs_emci.py
"""

import io
import pathlib
import sys

import pandas as pd

if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220")
TSV_FILE = pathlib.Path("D:/ADNI_BIDS_project/bids/genotype/subjects_with_snp_and_mri.tsv")
FAM_FILE = DATA_DIR / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr.fam"
PHENO_OUT = DATA_DIR / "phenotype_cn_vs_emci.pheno"

# ── Load ──────────────────────────────────────────────────────────────────────
print(f"Reading FAM : {FAM_FILE}")
fam = pd.read_csv(FAM_FILE, sep=r"\s+", header=None,
                  names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO_FAM"])
print(f"  {len(fam)} subjects")

print(f"Reading TSV : {TSV_FILE}")
tsv = pd.read_csv(TSV_FILE, sep="\t")
tsv = tsv.rename(columns={"adni_subject_id": "IID"})

# ── Merge ─────────────────────────────────────────────────────────────────────
df = fam[["FID", "IID"]].merge(
    tsv[["IID", "diagnosis_bl"]], on="IID", how="left"
)

# ── Encode phenotype ──────────────────────────────────────────────────────────
CONTROL = {"CognitivelyNormal"}   # PHENO = 1
CASE    = {"EarlyMCI"}            # PHENO = 2
# LateMCI + AlzheimersDisease → -9 (excluded from this contrast)

def encode(diag):
    if pd.isna(diag):
        return -9
    if diag in CASE:
        return 2
    if diag in CONTROL:
        return 1
    return -9

df["PHENO"] = df["diagnosis_bl"].apply(encode)

# ── Report ────────────────────────────────────────────────────────────────────
print("\nPhenotype distribution (CN vs EMCI):")
labels = {1: "CognitivelyNormal (control)", 2: "EarlyMCI (case)", -9: "excluded / missing"}
for val, cnt in df["PHENO"].value_counts().sort_index().items():
    print(f"  PHENO={val:>3}  {labels.get(val, '?'):30s}  n={cnt}")

# ── Write ─────────────────────────────────────────────────────────────────────
df[["FID", "IID", "PHENO"]].to_csv(PHENO_OUT, sep="\t", index=False, na_rep="-9")
print(f"\nPhenotype file written: {PHENO_OUT}")
print("Next: run run_gwas_cn_vs_emci.sh from D:/ADNI_SNP_Omni2.5M_20140220/")
