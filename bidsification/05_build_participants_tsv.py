"""
Script 05 — Build participants.tsv and participants.json
==========================================================
Populates bids/participants.tsv with one row per unique subject,
drawing baseline demographics from ADNIMERGE.csv and APOE genotype
from APOERES.csv.

BIDS-required columns: participant_id
BIDS-recommended:      sex, age, handedness (if available)

Additional ADNI-specific columns included:
  diagnosis_bl (DX_bl), apoe4_dosage (APOE4), apoe_genotype,
  education_years, ethnic_category, race_category, marital_status

Output:
    bids/participants.tsv
    bids/participants.json  (JSON sidecar describing each column)
"""

import pandas as pd
import os
import json

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT = r"D:\ADNI_BIDS_project"
ADNIMERGE_CSV = os.path.join(PROJECT_ROOT, "sourcedata", "clinical", "ADNIMERGE.csv")
APOERES_CSV   = os.path.join(PROJECT_ROOT, "sourcedata", "clinical", "APOERES.csv")
SESSION_MAP_CSV = os.path.join(PROJECT_ROOT, "metadata", "session_map.csv")
BIDS_DIR      = os.path.join(PROJECT_ROOT, "bids")
OUT_TSV       = os.path.join(BIDS_DIR, "participants.tsv")
OUT_JSON      = os.path.join(BIDS_DIR, "participants.json")

# ── Load & filter ADNIMERGE to baseline ───────────────────────────────────────
print("Loading ADNIMERGE.csv ...")
merge = pd.read_csv(ADNIMERGE_CSV, low_memory=False)
# Use VISCODE == 'bl' for baseline demographics (age, DX_bl stable at bl)
bl = merge[merge["VISCODE"] == "bl"].copy()
# Deduplicate — keep one baseline row per subject
bl = bl.drop_duplicates(subset="PTID", keep="first")
print(f"  → {len(bl):,} unique subjects at baseline from ADNIMERGE")

# ── Load APOE genotype ─────────────────────────────────────────────────────────
print("Loading APOERES.csv ...")
apoe = pd.read_csv(APOERES_CSV, low_memory=False)
apoe_clean = apoe[["PTID", "APGEN1", "APGEN2"]].drop_duplicates("PTID").copy()
apoe_clean["apoe_genotype"] = apoe_clean["APGEN1"].astype(str) + "/" + apoe_clean["APGEN2"].astype(str)

# ── Load session map to get bids_sub labels and limit to subjects with imaging ─
print("Loading session map ...")
sess_map = pd.read_csv(SESSION_MAP_CSV, low_memory=False)
subjects_with_mri = sess_map[["SubjectID", "bids_sub"]].drop_duplicates("SubjectID")

# ── Merge all data ─────────────────────────────────────────────────────────────
df = subjects_with_mri.merge(bl, left_on="SubjectID", right_on="PTID", how="left")
df = df.merge(apoe_clean[["PTID", "apoe_genotype"]], left_on="SubjectID", right_on="PTID", how="left")

# ── Map ADNI encodings to BIDS-friendly values ─────────────────────────────────
# Sex: PTGENDER in ADNI is 'Male'/'Female'
df["sex"] = df["PTGENDER"].astype(str).str.strip().str.lower().replace(
    {"male": "M", "female": "F", "nan": "n/a"}
)

# Diagnosis at baseline
dx_map = {
    "CN": "CognitivelyNormal",
    "SMC": "SubjectiveMemoryConcern",
    "EMCI": "EarlyMCI",
    "LMCI": "LateMCI",
    "MCI": "MCI",
    "AD": "AlzheimersDisease",
}
df["diagnosis_bl"] = df["DX_bl"].map(dx_map).fillna(df["DX_bl"]).fillna("n/a")

# Education (already in years)
df["education_years"] = pd.to_numeric(df["PTEDUCAT"], errors="coerce").fillna("n/a")

# Ethnic category (ADNI codes: 1=Hisp/Latino, 2=Not, 3=Unknown, 4=NA)
ethnic_map = {1: "HispanicOrLatino", 2: "NotHispanicOrLatino", 3: "Unknown", 4: "n/a"}
df["ethnic_category"] = df["PTETHCAT"].map(ethnic_map).fillna("n/a")

# Race (ADNI codes: 1=Am. Indian, 2=Asian, 3=NativeHawaiian, 4=Black, 5=White, 6=Multi, 7=Unknown)
race_map = {
    1: "AmericanIndianOrAlaskanNative",
    2: "Asian",
    3: "NativeHawaiianOrOtherPacificIslander",
    4: "BlackOrAfricanAmerican",
    5: "White",
    6: "MoreThanOneRace",
    7: "Unknown",
}
df["race_category"] = df["PTRACCAT"].map(race_map).fillna("n/a")

# APOE4 dosage (0/1/2)
df["apoe4_dosage"] = pd.to_numeric(df["APOE4"], errors="coerce").fillna("n/a")

# Age at baseline (ADNI reports to 1 decimal)
df["age"] = pd.to_numeric(df["AGE"], errors="coerce").fillna("n/a")

# Site
df["site"] = df["SITE"].fillna("n/a")

# ── Build output dataframe ─────────────────────────────────────────────────────
out = pd.DataFrame({
    "participant_id":   "sub-" + df["bids_sub"],
    "adni_subject_id":  df["SubjectID"],
    "age":              df["age"],
    "sex":              df["sex"],
    "diagnosis_bl":     df["diagnosis_bl"],
    "education_years":  df["education_years"],
    "ethnic_category":  df["ethnic_category"],
    "race_category":    df["race_category"],
    "apoe4_dosage":     df["apoe4_dosage"],
    "apoe_genotype":    df["apoe_genotype"].fillna("n/a"),
    "site":             df["site"],
})
out = out.sort_values("participant_id").reset_index(drop=True)

# Replace NaN with 'n/a' per BIDS convention
out = out.fillna("n/a")

out.to_csv(OUT_TSV, sep="\t", index=False)
print(f"Saved: {OUT_TSV}  ({len(out):,} subjects)")

# ── JSON sidecar ───────────────────────────────────────────────────────────────
participants_json = {
    "participant_id": {"Description": "Unique participant identifier in BIDS format (sub-<label>)"},
    "adni_subject_id": {"Description": "Original ADNI participant identifier (e.g., 002_S_0295)"},
    "age": {
        "Description": "Age of the participant at baseline visit",
        "Units": "years"
    },
    "sex": {
        "Description": "Biological sex of participant",
        "Levels": {"M": "Male", "F": "Female"}
    },
    "diagnosis_bl": {
        "Description": "Clinical diagnosis group at ADNI baseline visit (VISCODE=bl)",
        "Levels": {
            "CognitivelyNormal": "Cognitively normal control",
            "SubjectiveMemoryConcern": "Subjective memory concern",
            "EarlyMCI": "Early mild cognitive impairment",
            "LateMCI": "Late mild cognitive impairment",
            "MCI": "Mild cognitive impairment (unspecified)",
            "AlzheimersDisease": "Alzheimer's disease dementia",
        }
    },
    "education_years": {
        "Description": "Years of education completed",
        "Units": "years"
    },
    "ethnic_category": {
        "Description": "Ethnicity as reported in ADNI (PTETHCAT)",
        "Levels": {
            "HispanicOrLatino": "Hispanic or Latino",
            "NotHispanicOrLatino": "Not Hispanic or Latino",
            "Unknown": "Unknown or not reported",
            "n/a": "Not available"
        }
    },
    "race_category": {
        "Description": "Self-reported race category (PTRACCAT)",
        "Levels": {
            "AmericanIndianOrAlaskanNative": "American Indian or Alaskan Native",
            "Asian": "Asian",
            "NativeHawaiianOrOtherPacificIslander": "Native Hawaiian or Other Pacific Islander",
            "BlackOrAfricanAmerican": "Black or African American",
            "White": "White",
            "MoreThanOneRace": "More than one race",
            "Unknown": "Unknown or not reported"
        }
    },
    "apoe4_dosage": {
        "Description": "Number of APOE ε4 alleles (0, 1, or 2)",
        "Levels": {"0": "Non-carrier", "1": "Heterozygous carrier", "2": "Homozygous carrier"}
    },
    "apoe_genotype": {
        "Description": "Full APOE genotype (e.g., 3/4), derived from APOERES.csv"
    },
    "site": {
        "Description": "ADNI acquisition site identifier"
    }
}

with open(OUT_JSON, "w") as f:
    json.dump(participants_json, f, indent=2)
print(f"Saved: {OUT_JSON}")

# Summary
print("\nDiagnosis breakdown at baseline:")
print(out["diagnosis_bl"].value_counts().to_string())
print("\nSex breakdown:")
print(out["sex"].value_counts().to_string())
