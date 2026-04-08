"""
Script 06 — Populate phenotype/ Directory
==========================================
Copies selected clinical CSV tables into bids/phenotype/ in BIDS format
and generates JSON sidecar descriptors for each.

Selected tables (covering key clinical domains):
  ADNIMERGE.csv            -> adnimerge.tsv     (core longitudinal: cognitive, biomarkers, volumes, dx)
  CDR.csv                  -> cdr.tsv           (Clinical Dementia Rating)
  ADAS.csv                 -> adas.tsv          (Alzheimer's Disease Assessment Scale)
  MMSE.csv                 -> mmse.tsv          (Mini-Mental State Examination)
  NEUROBAT.csv             -> neurobat.tsv      (Clock Drawing, Trail Making, ANART, Category Fluency)
  ECOGPT.csv               -> ecog_pt.tsv       (Everyday Cognition — self-rated)
  ECOGSP.csv               -> ecog_sp.tsv       (Everyday Cognition — study-partner-rated)
  GDSCALE.csv              -> gds.tsv           (Geriatric Depression Scale)
  FAQ.csv                  -> faq.tsv           (Functional Activities Questionnaire)
  AMNART.csv               -> amnart.tsv        (American National Adult Reading Test)
  DXSUM_PDXCONV_ADNIALL.csv-> dxsum.tsv         (Longitudinal diagnosis + conversion status per visit)
  APOERES.csv              -> apoe.tsv          (APOE genotype)
  ARM.csv                  -> study_arm.tsv     (ADNI arm / protocol assignment)
  BIOMARK.csv              -> biomarkers.tsv    (CSF / plasma biomarkers — if exists)
  BACKMEDS.csv             -> backmeds.tsv      (background medications — if exists)

BIDS phenotype format:
  • TSV: tab-separated, 'participant_id' column, 'session_id' if longitudinal
  • JSON: parallel sidecar with column descriptions

Each TSV maps PTID -> participant_id (BIDS sub- label) using session_map.csv.
VISCODE columns are kept for longitudinal linkage.

Note: The full list of 298 clinical CSVs is intentionally NOT included here —
only the subset most relevant to sMRI / cognitive analyses of ADNI.
"""

import pandas as pd
import os
import json
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exclusions import is_excluded_subject

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT   = r"D:\ADNI_BIDS_project"
CLINICAL_DIR   = os.path.join(PROJECT_ROOT, "sourcedata", "clinical")
SESSION_MAP    = os.path.join(PROJECT_ROOT, "metadata", "session_map.csv")
SES_VISIT      = os.path.join(PROJECT_ROOT, "metadata", "ses_to_visit_code.csv")
PHENO_DIR      = os.path.join(PROJECT_ROOT, "bids", "phenotype")
os.makedirs(PHENO_DIR, exist_ok=True)

# ── Load subject ID mapping ────────────────────────────────────────────────────
print("Loading session map ...")
sess = pd.read_csv(SESSION_MAP, low_memory=False)
# Apply exclusions before building the PTID -> bids_sub lookup
n_before_excl = sess["SubjectID"].nunique()
sess = sess[~sess["SubjectID"].apply(is_excluded_subject)]
print(f"  -> {n_before_excl - sess['SubjectID'].nunique():,} excluded subjects removed from phenotype outputs")
ptid_to_bids = (
    sess[["SubjectID", "bids_sub"]]
    .drop_duplicates("SubjectID")
    .set_index("SubjectID")["bids_sub"]
    .to_dict()
)
print(f"  -> {len(ptid_to_bids):,} subject ID mappings")

# ── BIDS phenotype JSON sidecars ───────────────────────────────────────────────
PHENOTYPE_JSONS = {
    "adnimerge": {
        "_description": "ADNIMERGE longitudinal summary: cognitive tests, biomarkers, brain volumes, diagnosis trajectory.",
        "participant_id": {"Description": "BIDS participant identifier (sub-<label>)"},
        "adni_ptid": {"Description": "Original ADNI participant ID (e.g., 002_S_0295)"},
        "session_id": {"Description": "BIDS session label derived from EXAMDATE (ses-YYYYMMDD)"},
        "VISCODE": {"Description": "ADNI visit code (bl, m06, m12, ...)"},
        "DX": {"Description": "Clinical diagnosis at visit"},
        "DX_bl": {"Description": "Clinical diagnosis at baseline",
                  "Levels": {"CN": "Cognitively Normal", "SMC": "Subjective Memory Concern",
                             "EMCI": "Early MCI", "LMCI": "Late MCI", "AD": "Alzheimer's Disease"}},
        "AGE": {"Description": "Age at baseline", "Units": "years"},
        "CDRSB": {"Description": "CDR Sum of Boxes score"},
        "MMSE": {"Description": "Mini-Mental State Examination total score (0-30)"},
        "ADAS13": {"Description": "ADAS-Cog 13-item total score (higher = worse cognition)"},
        "RAVLT_immediate": {"Description": "Rey Auditory Verbal Learning Test — immediate recall"},
        "FAQ": {"Description": "Functional Activities Questionnaire"},
        "Hippocampus": {"Description": "Hippocampal volume (FreeSurfer)", "Units": "mm3"},
        "Entorhinal": {"Description": "Entorhinal cortex volume (FreeSurfer)", "Units": "mm3"},
        "WholeBrain": {"Description": "Whole brain volume (FreeSurfer)", "Units": "mm3"},
        "ICV": {"Description": "Intracranial volume (FreeSurfer)", "Units": "mm3"},
        "ABETA": {"Description": "CSF Amyloid-beta 1-42 level"},
        "TAU": {"Description": "CSF total tau level"},
        "PTAU": {"Description": "CSF phosphorylated tau-181 level"},
        "AV45": {"Description": "Florbetapir (AV-45) PET SUVR (reference: whole cerebellum)"},
        "APOE4": {"Description": "APOE ε4 allele dosage", "Levels": {"0": "Non-carrier", "1": "1 allele", "2": "2 alleles"}},
    },
    "cdr": {
        "_description": "Clinical Dementia Rating (CDR) assessments across all visits.",
        "participant_id": {"Description": "BIDS participant identifier"},
        "VISCODE": {"Description": "ADNI visit code"},
        "EXAMDATE": {"Description": "Date of examination"},
        "CDGLOBAL": {"Description": "CDR Global score (0=no dementia, 0.5=questionable, 1=mild, 2=moderate, 3=severe)"},
        "CDRSB": {"Description": "CDR Sum of Boxes (0-18)"},
    },
    "adas": {
        "_description": "ADAS-Cog cognitive assessments (11-item and 13-item versions).",
        "participant_id": {"Description": "BIDS participant identifier"},
        "VISCODE": {"Description": "ADNI visit code"},
        "TOTSCORE": {"Description": "ADAS 11-item total score (higher = more impaired, max=70)"},
        "Q4SCORE": {"Description": "ADAS delayed word recall score"},
    },
    "mmse": {
        "_description": "Mini-Mental State Examination scores across all visits.",
        "participant_id": {"Description": "BIDS participant identifier"},
        "VISCODE": {"Description": "ADNI visit code"},
        "MMSCORE": {"Description": "MMSE total score (0-30; lower = more impaired)"},
    },
    "apoe": {
        "_description": "APOE genotyping results (from APOERES.csv).",
        "participant_id": {"Description": "BIDS participant identifier"},
        "APGEN1": {"Description": "APOE allele 1 (2, 3, or 4)"},
        "APGEN2": {"Description": "APOE allele 2 (2, 3, or 4)"},
        "apoe_genotype": {"Description": "Combined APOE genotype (e.g., 3/4)"},
    },
    "study_arm": {
        "_description": "ADNI study protocol arm assignment for each participant.",
        "participant_id": {"Description": "BIDS participant identifier"},
        "COLPROT": {"Description": "Collection protocol (ADNI1, ADNIGO, ADNI2, ADNI3)"},
        "ARM": {"Description": "Study arm / group assignment"},
    },
    "biomarkers": {
        "_description": "Biomarker summary table (if BIOMARK.csv exists).",
        "participant_id": {"Description": "BIDS participant identifier"},
        "VISCODE": {"Description": "ADNI visit code"},
    },
    "neurobat": {
        "_description": "Neuropsychological Battery: AVLT (trial-level), Clock Drawing Test, Trail Making Test (A&B), ANART, Category Fluency.",
        "participant_id": {"Description": "BIDS participant identifier"},
        "VISCODE": {"Description": "ADNI visit code"},
        "AVTOT1":    {"Description": "AVLT Learning Trial 1: total words recalled"},
        "AVTOT2":    {"Description": "AVLT Learning Trial 2: total words recalled"},
        "AVTOT3":    {"Description": "AVLT Learning Trial 3: total words recalled"},
        "AVTOT4":    {"Description": "AVLT Learning Trial 4: total words recalled"},
        "AVTOT5":    {"Description": "AVLT Learning Trial 5: total words recalled"},
        "AVERR1":    {"Description": "AVLT Learning Trial 1: total errors (intrusions + repetitions)"},
        "AVERR2":    {"Description": "AVLT Learning Trial 2: total errors"},
        "AVERR3":    {"Description": "AVLT Learning Trial 3: total errors"},
        "AVERR4":    {"Description": "AVLT Learning Trial 4: total errors"},
        "AVERR5":    {"Description": "AVLT Learning Trial 5: total errors"},
        "AVTOTB":    {"Description": "AVLT List B (interference trial): total words recalled"},
        "AVERRB":    {"Description": "AVLT List B: total errors"},
        "AVTOT6":    {"Description": "AVLT Post-list-B (Trial 6) immediate recall: total words recalled"},
        "AVERR6":    {"Description": "AVLT Post-list-B recall: total errors"},
        "AVENDED":   {"Description": "AVLT: whether the test was ended early (1=yes)"},
        "AVDELBEGAN":{"Description": "AVLT: whether 30-minute delayed recall section began (1=yes)"},
        "AVDEL30MIN":{"Description": "AVLT 30-minute delayed recall: total words recalled"},
        "AVDELERR1": {"Description": "AVLT delayed recall: total intrusion errors"},
        "AVDELTOT":  {"Description": "AVLT delayed recognition: total correct recognitions"},
        "AVDELERR2": {"Description": "AVLT delayed recognition: total false-positive errors"},
        "CLOCKCIRC": {"Description": "Clock Drawing: approximately circular face (1=correct)"},
        "CLOCKSYM":  {"Description": "Clock Drawing: symmetry of number placement (1=correct)"},
        "CLOCKNUM":  {"Description": "Clock Drawing: correctness of numbers (1=correct)"},
        "CLOCKHAND": {"Description": "Clock Drawing: presence of two hands (1=correct)"},
        "CLOCKTIME": {"Description": "Clock Drawing: hands set to ten past eleven (1=correct)"},
        "CLOCKSCOR": {"Description": "Clock Drawing: total score (0-5)"},
        "TRAASCOR":  {"Description": "Trail Making Test Part A: time to complete (seconds)"},
        "TRAAERRCOM":{"Description": "Trail Making Test Part A: errors of commission"},
        "TRAAERROM": {"Description": "Trail Making Test Part A: errors of omission"},
        "TRABSCOR":  {"Description": "Trail Making Test Part B: time to complete (seconds)"},
        "TRABERRCOM":{"Description": "Trail Making Test Part B: errors of commission"},
        "TRABERROM": {"Description": "Trail Making Test Part B: errors of omission"},
        "ANARTERR":  {"Description": "American National Adult Reading Test: total errors"},
        "CATANIMSC": {"Description": "Category Fluency Test (Animals): total correct in 60 seconds"},
    },
    "ecog_pt": {
        "_description": "Everyday Cognition (ECog) — participant self-rated. Subscales: memory, language, visuospatial, planning, organisation, divided attention.",
        "participant_id": {"Description": "BIDS participant identifier"},
        "VISCODE": {"Description": "ADNI visit code"},
        "EcogPtMem":   {"Description": "ECog Memory subscale (self-rated)"},
        "EcogPtLang":  {"Description": "ECog Language subscale (self-rated)"},
        "EcogPtVisspat":{"Description": "ECog Visuospatial subscale (self-rated)"},
        "EcogPtPlan":  {"Description": "ECog Planning subscale (self-rated)"},
        "EcogPtOrgan": {"Description": "ECog Organisation subscale (self-rated)"},
        "EcogPtDivatt":{"Description": "ECog Divided Attention subscale (self-rated)"},
        "EcogPtTotal": {"Description": "ECog Total score (self-rated)"},
    },
    "ecog_sp": {
        "_description": "Everyday Cognition (ECog) — study-partner-rated. Same subscales as ecog_pt.",
        "participant_id": {"Description": "BIDS participant identifier"},
        "VISCODE": {"Description": "ADNI visit code"},
        "EcogSpMem":   {"Description": "ECog Memory subscale (study-partner-rated)"},
        "EcogSpLang":  {"Description": "ECog Language subscale (study-partner-rated)"},
        "EcogSpVisspat":{"Description": "ECog Visuospatial subscale (study-partner-rated)"},
        "EcogSpPlan":  {"Description": "ECog Planning subscale (study-partner-rated)"},
        "EcogSpOrgan": {"Description": "ECog Organisation subscale (study-partner-rated)"},
        "EcogSpDivatt":{"Description": "ECog Divided Attention subscale (study-partner-rated)"},
        "EcogSpTotal": {"Description": "ECog Total score (study-partner-rated)"},
    },
    "gds": {
        "_description": "Geriatric Depression Scale (GDS) — presence and severity of depressive symptoms.",
        "participant_id": {"Description": "BIDS participant identifier"},
        "VISCODE": {"Description": "ADNI visit code"},
        "GDTOTAL": {"Description": "GDS total score (0-15; >=5 suggests depression)"},
    },
    "faq": {
        "_description": "Functional Activities Questionnaire (FAQ) — instrumental activities of daily living.",
        "participant_id": {"Description": "BIDS participant identifier"},
        "VISCODE": {"Description": "ADNI visit code"},
        "FAQTOTAL": {"Description": "FAQ total score (0-30; higher = more dependent)"},
    },
    "amnart": {
        "_description": "American National Adult Reading Test (AMNART) — premorbid IQ estimate.",
        "participant_id": {"Description": "BIDS participant identifier"},
        "VISCODE": {"Description": "ADNI visit code"},
        "AMNARTID": {"Description": "AMNART IQ estimate"},
        "AMNARTERR": {"Description": "AMNART total errors"},
    },
    "dxsum": {
        "_description": "Longitudinal diagnosis and conversion status per visit (DXSUM_PDXCONV_ADNIALL). Use DXCHANGE and DIAGNOSIS to track within-study changes.",
        "participant_id": {"Description": "BIDS participant identifier"},
        "VISCODE": {"Description": "ADNI visit code"},
        "VISCODE2": {"Description": "ADNI visit code (VISCODE2 cumulative format)"},
        "EXAMDATE": {"Description": "Date of clinical assessment"},
        "ORIGPROT": {"Description": "Protocol the subject originally enrolled in (ADNI1/GO/2/3/4)"},
        "COLPROT":  {"Description": "Protocol under which this visit was collected"},
        "DIAGNOSIS": {"Description": "Consensus diagnosis at visit (1=CN, 2=MCI, 3=AD)",
                      "Levels": {"1": "Cognitively Normal", "2": "MCI", "3": "Alzheimer's Disease"}},
        "DXCHANGE": {"Description": "Diagnosis change code relative to previous visit",
                     "Levels": {
                         "1": "Stable: CN", "2": "Stable: MCI", "3": "Stable: AD",
                         "4": "Conversion: CN->MCI", "5": "Conversion: MCI->AD", "6": "Conversion: CN->AD",
                         "7": "Reversion: MCI->CN", "8": "Reversion: AD->MCI", "9": "Reversion: AD->CN"
                     }},
    },
}

# ── Table definitions: (source_csv, output_name, ptid_col, date_col) ──────────
TABLE_DEFS = [
    ("ADNIMERGE.csv",             "adnimerge",  "PTID", "EXAMDATE"),
    ("CDR.csv",                   "cdr",        "PTID", "EXAMDATE"),
    ("ADAS.csv",                  "adas",       "PTID", "EXAMDATE"),
    ("MMSE.csv",                  "mmse",       "PTID", "EXAMDATE"),
    ("NEUROBAT.csv",              "neurobat",   "PTID", "EXAMDATE"),
    ("ECOGPT.csv",                "ecog_pt",    "PTID", "EXAMDATE"),
    ("ECOGSP.csv",                "ecog_sp",    "PTID", "EXAMDATE"),
    ("GDSCALE.csv",               "gds",        "PTID", "EXAMDATE"),
    ("FAQ.csv",                   "faq",        "PTID", "EXAMDATE"),
    ("AMNART.csv",                "amnart",     "PTID", "EXAMDATE"),
    ("DXSUM_PDXCONV_ADNIALL.csv", "dxsum",      "PTID", "EXAMDATE"),
    ("APOERES.csv",               "apoe",       "PTID", None),
    ("ARM.csv",                   "study_arm",  "PTID", None),
    ("BIOMARK.csv",               "biomarkers", "PTID", "EXAMDATE"),  # may not exist
    ("BACKMEDS.csv",              "backmeds",   "PTID", "EXAMDATE"),  # may not exist
]

def normalize_date(date_str: str) -> str:
    if pd.isna(date_str):
        return "n/a"
    return str(date_str)[:10].replace("-", "")

# ── Process each table ─────────────────────────────────────────────────────────
for (csv_name, out_name, ptid_col, date_col) in TABLE_DEFS:
    src = os.path.join(CLINICAL_DIR, csv_name)
    if not os.path.isfile(src):
        print(f"  SKIP (not found): {csv_name}")
        continue
    
    print(f"Processing {csv_name} -> {out_name}.tsv ...")
    df = pd.read_csv(src, low_memory=False)
    
    # Add participant_id column
    if ptid_col not in df.columns:
        print(f"  WARNING: '{ptid_col}' column not found in {csv_name}, skipping.")
        continue
    
    df["participant_id"] = df[ptid_col].map(
        lambda x: f"sub-{ptid_to_bids[x]}" if x in ptid_to_bids else "n/a"
    )
    # Keep only subjects present in the imaging data
    df_filt = df[df["participant_id"] != "n/a"].copy()
    print(f"  -> {len(df_filt):,} rows for {df_filt['participant_id'].nunique():,} imaging subjects")
    
    # Add session_id if date column exists
    if date_col and date_col in df_filt.columns:
        df_filt["session_id"] = "ses-" + df_filt[date_col].apply(normalize_date)
    
    # Move participant_id to first column
    cols = ["participant_id"] + [c for c in df_filt.columns if c != "participant_id"]
    
    # Keep original PTID as adni_ptid for traceability
    if ptid_col in cols:
        cols_renamed = [c if c != ptid_col else "adni_ptid" for c in cols]
        df_filt = df_filt.rename(columns={ptid_col: "adni_ptid"})
        cols = cols_renamed
    
    df_filt = df_filt[cols]
    df_filt = df_filt.fillna("n/a")
    
    out_tsv = os.path.join(PHENO_DIR, f"{out_name}.tsv")
    df_filt.to_csv(out_tsv, sep="\t", index=False)
    print(f"  Saved: {out_tsv}")
    
    # Write JSON sidecar
    out_json_path = os.path.join(PHENO_DIR, f"{out_name}.json")
    pheno_json = PHENOTYPE_JSONS.get(out_name, {})
    # description at top level (remove internal _description key)
    top_desc = pheno_json.pop("_description", f"ADNI {csv_name} data for imaging subjects.")
    pheno_json_out = {"MeasurementToolMetadata": {"Description": top_desc}, **pheno_json}
    
    with open(out_json_path, "w") as f:
        json.dump(pheno_json_out, f, indent=2)
    print(f"  Saved: {out_json_path}")

print("\nDone. phenotype/ directory populated.")
