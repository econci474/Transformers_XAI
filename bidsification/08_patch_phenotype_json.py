"""
Script 08 — Patch Phenotype JSON Sidecars with Full Column Descriptions
========================================================================
Adds ADNI data-dictionary-aligned descriptions for all TSV columns that
currently lack JSON documentation. Descriptions are based on the ADNI
data dictionary (adni.loni.usc.edu), CRF documentation, and ADNIMERGE
codebook.

Run once after script 06. Safe to re-run (only adds missing keys).
"""

import json, os, pandas as pd

PHENO_DIR = r"D:\ADNI_BIDS_project\bids\phenotype"

# ── Shared / administrative columns that appear in multiple tables ─────────────
SHARED = {
    "participant_id":       {"Description": "BIDS participant identifier (sub-<label>)"},
    "adni_ptid":            {"Description": "Original ADNI participant ID (e.g., 002_S_0295)"},
    "session_id":           {"Description": "BIDS session label (ses-<VISCODE2>)"},
    "RID":                  {"Description": "ADNI roster ID (numeric participant identifier)"},
    "ORIGPROT":             {"Description": "Protocol in which the participant originally enrolled (ADNI1 / ADNIGO / ADNI2 / ADNI3 / ADNI4)"},
    "COLPROT":              {"Description": "Protocol under which this visit was collected (ADNI1 / ADNIGO / ADNI2 / ADNI3)"},
    "VISCODE":              {"Description": "ADNI visit code short form (bl, m06, m12, …)"},
    "VISCODE2":             {"Description": "ADNI visit code VISCODE2 — cumulative months from baseline (bl, m06, m12, m24, …)"},
    "VISDATE":              {"Description": "Date of the clinical visit (YYYY-MM-DD)"},
    "EXAMDATE":             {"Description": "Date of examination / assessment (YYYY-MM-DD)"},
    "EXAMDATE_bl":          {"Description": "Date of baseline examination (YYYY-MM-DD)"},
    "PHASE":                {"Description": "ADNI study phase at the time of this assessment"},
    "SITEID":               {"Description": "ADNI site identifier"},
    "SPID":                 {"Description": "Study partner identifier"},
    "SOURCE":               {"Description": "Data source / collection site label"},
    "ID":                   {"Description": "Internal ADNI database record ID"},
    "USERDATE":             {"Description": "Date the record was entered into the ADNI database"},
    "USERDATE2":            {"Description": "Date the record was last updated in the ADNI database"},
    "DD_CRF_VERSION_LABEL": {"Description": "Version label of the CRF (Case Report Form) used for data collection"},
    "LANGUAGE_CODE":        {"Description": "Language code used for test administration (e.g., EN for English)"},
    "HAS_QC_ERROR":         {"Description": "Flag indicating whether a quality control error was flagged for this record (1=yes, 0=no)"},
    "update_stamp":         {"Description": "Timestamp of the most recent database update for this record"},
    "DONE":                 {"Description": "Flag indicating whether the test/assessment was completed (1=yes)"},
    "NDREASON":             {"Description": "Reason assessment was not done / not completed"},
}

# ── Per-table column descriptions ──────────────────────────────────────────────
TABLE_COLS = {

    # ── ADAS ──────────────────────────────────────────────────────────────────
    "adas": {
        "TOTAL13":  {"Description": "ADAS-Cog 13-item total score (higher score = more impaired; max=85)"},
        "VISCODE2": SHARED["VISCODE2"],
        "VISDATE":  SHARED["VISDATE"],
    },

    # ── ADNIMERGE ─────────────────────────────────────────────────────────────
    "adnimerge": {
        "RID":       SHARED["RID"],
        "COLPROT":   SHARED["COLPROT"],
        "ORIGPROT":  SHARED["ORIGPROT"],
        "SITE":      {"Description": "ADNI clinical site code"},
        "EXAMDATE":  SHARED["EXAMDATE"],
        "PTGENDER":  {"Description": "Participant sex as recorded at enrolment (1=Male, 2=Female)"},
        "PTEDUCAT":  {"Description": "Years of education completed"},
        "PTETHCAT":  {"Description": "Ethnicity category (1=Hispanic or Latino, 2=Not Hispanic or Latino, 3=Unknown)"},
        "PTRACCAT":  {"Description": "Race category (1=Am. Indian/Alaskan Native, 2=Asian, 3=Native Hawaiian/Pacific Islander, 4=Black/African American, 5=White, 6=More than one, 7=Unknown)"},
        "PTMARRY":   {"Description": "Marital status (1=Married, 2=Widowed, 3=Divorced, 4=Separated, 5=Never married, 6=Unknown)"},
        "FDG":       {"Description": "FDG-PET mean ROI SUVR (angular/temporal/posterior cingulate; reference: pons)"},
        "PIB":       {"Description": "Pittsburgh Compound-B (PIB) PET amyloid SUVR"},
        "FBB":       {"Description": "Florbetaben (FBB) PET amyloid SUVR (reference: whole cerebellum)"},
        "ADAS11":    {"Description": "ADAS-Cog 11-item total score (higher = more impaired; max=70)"},
        "ADASQ4":    {"Description": "ADAS-Cog Q4: delayed word recall score"},
        "RAVLT_learning":          {"Description": "RAVLT learning score: Trial 5 minus Trial 1 total (learning slope)"},
        "RAVLT_forgetting":        {"Description": "RAVLT forgetting score: Trial 5 minus 30-minute delayed recall total"},
        "RAVLT_perc_forgetting":   {"Description": "RAVLT percent forgetting: forgetting score / Trial 5 × 100"},
        "LDELTOTAL":               {"Description": "Logical Memory II (delayed recall) total score — WMS-R paragraph recall"},
        "DIGITSCOR":               {"Description": "Digit Span total score (forward + backward)"},
        "TRABSCOR":                {"Description": "Trail Making Test Part B: time to complete (seconds)"},
        "MOCA":                    {"Description": "Montreal Cognitive Assessment (MoCA) total score (0-30; <26 = mild cognitive impairment)"},
        "EcogPtMem":               {"Description": "Everyday Cognition (ECog) — participant-rated Memory subscale"},
        "EcogPtLang":              {"Description": "Everyday Cognition (ECog) — participant-rated Language subscale"},
        "EcogPtVisspat":           {"Description": "Everyday Cognition (ECog) — participant-rated Visuospatial Abilities subscale"},
        "EcogPtPlan":              {"Description": "Everyday Cognition (ECog) — participant-rated Planning subscale"},
        "EcogPtOrgan":             {"Description": "Everyday Cognition (ECog) — participant-rated Organisation subscale"},
        "EcogPtDivatt":            {"Description": "Everyday Cognition (ECog) — participant-rated Divided Attention subscale"},
        "EcogPtTotal":             {"Description": "Everyday Cognition (ECog) — participant-rated Total score"},
        "EcogSPMem":               {"Description": "Everyday Cognition (ECog) — study-partner-rated Memory subscale"},
        "EcogSPLang":              {"Description": "Everyday Cognition (ECog) — study-partner-rated Language subscale"},
        "EcogSPVisspat":           {"Description": "Everyday Cognition (ECog) — study-partner-rated Visuospatial Abilities subscale"},
        "EcogSPPlan":              {"Description": "Everyday Cognition (ECog) — study-partner-rated Planning subscale"},
        "EcogSPOrgan":             {"Description": "Everyday Cognition (ECog) — study-partner-rated Organisation subscale"},
        "EcogSPDivatt":            {"Description": "Everyday Cognition (ECog) — study-partner-rated Divided Attention subscale"},
        "EcogSPTotal":             {"Description": "Everyday Cognition (ECog) — study-partner-rated Total score"},
        "FLDSTRENG":               {"Description": "MRI magnetic field strength at the corresponding imaging visit (Tesla)"},
        "FSVERSION":               {"Description": "FreeSurfer version used for structural MRI parcellation"},
        "IMAGEUID":                {"Description": "ADNI image identifier (ImageUID) used for the corresponding MRI"},
        "Ventricles":              {"Description": "Lateral ventricular volume (FreeSurfer; mm³)"},
        "Fusiform":                {"Description": "Fusiform gyrus volume (FreeSurfer; mm³)"},
        "MidTemp":                 {"Description": "Middle temporal gyrus volume (FreeSurfer; mm³)"},
        "mPACCdigit":              {"Description": "Modified Preclinical Alzheimer Cognitive Composite — Digit Symbol Coding subscale"},
        "mPACCtrailsB":            {"Description": "Modified Preclinical Alzheimer Cognitive Composite — Trail Making Test B subscale"},
        # Baseline-value columns (_bl suffix)
        "EXAMDATE_bl":             SHARED["EXAMDATE_bl"],
        "CDRSB_bl":                {"Description": "CDR Sum of Boxes at baseline"},
        "ADAS11_bl":               {"Description": "ADAS-Cog 11-item total score at baseline"},
        "ADAS13_bl":               {"Description": "ADAS-Cog 13-item total score at baseline"},
        "ADASQ4_bl":               {"Description": "ADAS-Cog Q4 delayed word recall score at baseline"},
        "MMSE_bl":                 {"Description": "MMSE total score at baseline"},
        "RAVLT_immediate_bl":      {"Description": "RAVLT immediate recall score at baseline"},
        "RAVLT_learning_bl":       {"Description": "RAVLT learning score at baseline"},
        "RAVLT_forgetting_bl":     {"Description": "RAVLT forgetting score at baseline"},
        "RAVLT_perc_forgetting_bl":{"Description": "RAVLT percent forgetting at baseline"},
        "LDELTOTAL_BL":            {"Description": "Logical Memory II delayed recall score at baseline"},
        "DIGITSCOR_bl":            {"Description": "Digit Span total score at baseline"},
        "TRABSCOR_bl":             {"Description": "Trail Making Test Part B time to complete at baseline (seconds)"},
        "FAQ_bl":                  {"Description": "Functional Activities Questionnaire total score at baseline"},
        "mPACCdigit_bl":           {"Description": "mPACC Digit Symbol Coding subscale at baseline"},
        "mPACCtrailsB_bl":         {"Description": "mPACC Trail Making Test B subscale at baseline"},
        "FLDSTRENG_bl":            {"Description": "MRI field strength at baseline imaging visit (Tesla)"},
        "FSVERSION_bl":            {"Description": "FreeSurfer version used at baseline imaging"},
        "IMAGEUID_bl":             {"Description": "ADNI ImageUID for the baseline MRI scan"},
        "Ventricles_bl":           {"Description": "Lateral ventricular volume at baseline (FreeSurfer; mm³)"},
        "Hippocampus_bl":          {"Description": "Hippocampal volume at baseline (FreeSurfer; mm³)"},
        "WholeBrain_bl":           {"Description": "Whole brain volume at baseline (FreeSurfer; mm³)"},
        "Entorhinal_bl":           {"Description": "Entorhinal cortex volume at baseline (FreeSurfer; mm³)"},
        "Fusiform_bl":             {"Description": "Fusiform gyrus volume at baseline (FreeSurfer; mm³)"},
        "MidTemp_bl":              {"Description": "Middle temporal gyrus volume at baseline (FreeSurfer; mm³)"},
        "ICV_bl":                  {"Description": "Intracranial volume at baseline (FreeSurfer; mm³)"},
        "MOCA_bl":                 {"Description": "MoCA total score at baseline"},
        "EcogPtMem_bl":            {"Description": "ECog participant-rated Memory subscale at baseline"},
        "EcogPtLang_bl":           {"Description": "ECog participant-rated Language subscale at baseline"},
        "EcogPtVisspat_bl":        {"Description": "ECog participant-rated Visuospatial subscale at baseline"},
        "EcogPtPlan_bl":           {"Description": "ECog participant-rated Planning subscale at baseline"},
        "EcogPtOrgan_bl":          {"Description": "ECog participant-rated Organisation subscale at baseline"},
        "EcogPtDivatt_bl":         {"Description": "ECog participant-rated Divided Attention subscale at baseline"},
        "EcogPtTotal_bl":          {"Description": "ECog participant-rated Total score at baseline"},
        "EcogSPMem_bl":            {"Description": "ECog study-partner-rated Memory subscale at baseline"},
        "EcogSPLang_bl":           {"Description": "ECog study-partner-rated Language subscale at baseline"},
        "EcogSPVisspat_bl":        {"Description": "ECog study-partner-rated Visuospatial subscale at baseline"},
        "EcogSPPlan_bl":           {"Description": "ECog study-partner-rated Planning subscale at baseline"},
        "EcogSPOrgan_bl":          {"Description": "ECog study-partner-rated Organisation subscale at baseline"},
        "EcogSPDivatt_bl":         {"Description": "ECog study-partner-rated Divided Attention subscale at baseline"},
        "EcogSPTotal_bl":          {"Description": "ECog study-partner-rated Total score at baseline"},
        "ABETA_bl":                {"Description": "CSF Amyloid-beta 1-42 concentration at baseline (pg/mL)"},
        "TAU_bl":                  {"Description": "CSF total tau concentration at baseline (pg/mL)"},
        "PTAU_bl":                 {"Description": "CSF phosphorylated tau-181 concentration at baseline (pg/mL)"},
        "FDG_bl":                  {"Description": "FDG-PET mean ROI SUVR at baseline"},
        "PIB_bl":                  {"Description": "PIB PET amyloid SUVR at baseline"},
        "AV45_bl":                 {"Description": "Florbetapir (AV-45) PET amyloid SUVR at baseline"},
        "FBB_bl":                  {"Description": "Florbetaben (FBB) PET amyloid SUVR at baseline"},
        "Years_bl":                {"Description": "Time since baseline (years)"},
        "Month_bl":                {"Description": "Time since baseline (months, rounded)"},
        "Month":                   {"Description": "Nominal month visit number derived from VISCODE2"},
        "M":                       {"Description": "Nominal month offset from baseline for this row"},
        "update_stamp":            SHARED["update_stamp"],
    },

    # ── AMNART ────────────────────────────────────────────────────────────────
    "amnart": {
        "VISCODE2":  SHARED["VISCODE2"],
        "DONE":      SHARED["DONE"],
        "NDREASON":  SHARED["NDREASON"],
        "EXAMDATE":  SHARED["EXAMDATE"],
        "SOURCE":    SHARED["SOURCE"],
        "session_id":SHARED["session_id"],
        **{f"AMNART{i}": {"Description": f"American National Adult Reading Test — item {i}: 1=correctly pronounced, 0=error"}
           for i in range(1, 51)},
    },

    # ── APOE ──────────────────────────────────────────────────────────────────
    "apoe": {
        "VISCODE":   SHARED["VISCODE"],
        "GENOTYPE":  {"Description": "Full APOE genotype string (e.g., 3/3, 3/4, 4/4)"},
        "APTESTDT":  {"Description": "Date of APOE blood draw (YYYY-MM-DD)"},
        "APVOLUME":  {"Description": "Volume of blood drawn for APOE genotyping (mL)"},
        "APRECEIVE": {"Description": "Date/time sample received at genotyping laboratory"},
        "APAMBTEMP": {"Description": "Ambient temperature during sample transport (°C)"},
        "APRESAMP":  {"Description": "Flag indicating whether a re-sample was required (1=yes)"},
        "APUSABLE":  {"Description": "Flag indicating whether the sample was usable for genotyping (1=yes, 0=no)"},
    },

    # ── BACKMEDS ──────────────────────────────────────────────────────────────
    "backmeds": {
        "participant_id": SHARED["participant_id"],
        "VISCODE":  SHARED["VISCODE"],
        "VISCODE2": SHARED["VISCODE2"],
        "VISDATE":  SHARED["VISDATE"],
        "KEYMED":   {"Description": "Key medication taken (free-text or coded medication name)"},
    },

    # ── BIOMARKERS ────────────────────────────────────────────────────────────
    "biomarkers": {
        "VISCODE2":     SHARED["VISCODE2"],
        "EXAMDATE":     SHARED["EXAMDATE"],
        "RECNO":        {"Description": "Record number (for studies with multiple biospecimen collections per visit)"},
        "BIBLOOD":      {"Description": "Blood sample collected at this visit (1=yes, 0=no)"},
        "BIURINE":      {"Description": "Urine sample collected at this visit (1=yes, 0=no)"},
        "BLREASON":     {"Description": "Reason blood draw was not completed (if applicable)"},
        "BIFAST":       {"Description": "Participant was fasting before blood draw (1=yes, 0=no)"},
        "BITIME":       {"Description": "Time of blood draw (HH:MM)"},
        "BIREDTIME":    {"Description": "Time red-top tube blood was drawn (HH:MM)"},
        "BIREDAMT":     {"Description": "Amount of blood collected in red-top tube (mL)"},
        "BIREDCENT":    {"Description": "Red-top tube centrifuged at site (1=yes)"},
        "BIREDTRNS":    {"Description": "Red-top tube sample transferred to lab (1=yes)"},
        "BIREDVOL":     {"Description": "Volume of red-top tube aliquot sent to biorepository (mL)"},
        "BIREDFROZ":    {"Description": "Red-top tube sample frozen at site before shipping (1=yes)"},
        "BILAVTIME":    {"Description": "Time lavender-top (EDTA) tube blood was drawn"},
        "BINEEDLE":     {"Description": "Needle size used for blood draw"},
        "BILAVAMT":     {"Description": "Amount of blood collected in lavender-top (EDTA) tube (mL)"},
        "BILAVCENT":    {"Description": "Lavender-top tube centrifuged at site (1=yes)"},
        "BILAVTRNS":    {"Description": "Lavender-top tube sample transferred to lab (1=yes)"},
        "BILAVVOL":     {"Description": "Volume of lavender-top tube aliquot sent to biorepository (mL)"},
        "BILAVFROZ":    {"Description": "Lavender-top tube sample frozen at site (1=yes)"},
        "BIURITIME":    {"Description": "Time urine sample was collected (HH:MM)"},
        "BIURIAMT":     {"Description": "Amount of urine collected (mL)"},
        "BIURITRNS":    {"Description": "Urine sample transferred to lab (1=yes)"},
        "BIURIVOL":     {"Description": "Volume of urine aliquot sent to biorepository (mL)"},
        "BIURIFROZ":    {"Description": "Urine sample frozen at site (1=yes)"},
        "BICSF":        {"Description": "CSF sample collected at this visit (1=yes, 0=no)"},
        "BINONE":       {"Description": "No biospecimen collected at this visit (1=yes)"},
        "REASON":       {"Description": "Reason no biospecimen was collected"},
        "BICSFFAST":    {"Description": "Participant was fasting before lumbar puncture (1=yes)"},
        "BICSFTIME":    {"Description": "Time of lumbar puncture / CSF collection (HH:MM)"},
        "BIMETHOD":     {"Description": "CSF collection method used (e.g., LP with or without fluoroscopy)"},
        "NEEDLESIZE":   {"Description": "Spinal needle gauge used for lumbar puncture"},
        "INTERSPACE":   {"Description": "Intervertebral space used for lumbar puncture (e.g., L3-L4, L4-L5)"},
        "PTPOSITION":   {"Description": "Patient positioning during lumbar puncture (e.g., lateral decubitus, seated)"},
        "COLTUBETYP":   {"Description": "Type of tube used to collect CSF"},
        "SHPTUBETYP":   {"Description": "Type of tube used to ship CSF to biorepository"},
        "TUBEMIN":      {"Description": "Number of tubes of CSF collected"},
        "BICSFAMT":     {"Description": "Total amount of CSF collected (mL)"},
        "BICSFTRNS":    {"Description": "CSF sample transferred to biorepository (1=yes)"},
        "BICSFVOL":     {"Description": "Volume of CSF aliquot sent to biorepository (mL)"},
        "BICSFFROZ":    {"Description": "CSF sample frozen at site before shipping (1=yes)"},
        "BILPPATCH":    {"Description": "Blood patch applied after lumbar puncture (1=yes)"},
        "BILPFLURO":    {"Description": "Fluoroscopy used during lumbar puncture (1=yes)"},
        "BILPSPFILM":   {"Description": "Spot film taken during fluoroscopy-guided LP (1=yes)"},
        "BILPPADATE":   {"Description": "Date blood patch was applied after LP"},
        "BILPFLDATE":   {"Description": "Date fluoroscopy was performed"},
        "BILPSPDATE":   {"Description": "Date spot film was taken"},
        "BILPPADATEYR_DRVD": {"Description": "Year blood patch was applied (derived)"},
        "BILPFLDATEYR_DRVD": {"Description": "Year fluoroscopy was performed (derived)"},
        "BILPSPDATEYR_DRVD": {"Description": "Year spot film was taken (derived)"},
        "BILPOTPROC":   {"Description": "Other procedure performed at LP visit"},
        "BIFEDDATE":    {"Description": "Date samples were shipped to biorepository (FedEx date)"},
        "session_id":   SHARED["session_id"],
    },

    # ── CDR ────────────────────────────────────────────────────────────────────
    "cdr": {
        "VISCODE2":  SHARED["VISCODE2"],
        "VISDATE":   SHARED["VISDATE"],
        "CDSOURCE":  {"Description": "Source of CDR information (1=participant, 2=study partner, 3=both)"},
        "CDVERSION": {"Description": "CDR version used (standard or expanded)"},
        "SPID":      SHARED["SPID"],
        "CDMEMORY":  {"Description": "CDR Memory domain score (0=none, 0.5=questionable, 1=mild, 2=moderate, 3=severe)"},
        "CDORIENT":  {"Description": "CDR Orientation domain score (0–3)"},
        "CDJUDGE":   {"Description": "CDR Judgment and Problem Solving domain score (0–3)"},
        "CDCOMMUN":  {"Description": "CDR Community Affairs domain score (0–3)"},
        "CDHOME":    {"Description": "CDR Home and Hobbies domain score (0–3)"},
        "CDCARE":    {"Description": "CDR Personal Care domain score (0–3)"},
    },

    # ── DXSUM ──────────────────────────────────────────────────────────────────
    "dxsum": {
        "PHASE":     SHARED["PHASE"],
        "DXNORM":    {"Description": "Diagnosis: normal cognition (1=yes, -4=not applicable)"},
        "DXNODEP":   {"Description": "Normal cognition sub-type: no depression (1=yes, -4=N/A)"},
        "DXMCI":     {"Description": "Diagnosis: meets MCI criteria (1=yes, -4=N/A)"},
        "DXMDES":    {"Description": "MCI subtype description (text)"},
        "DXMPTR1":   {"Description": "MCI criterion 1: memory complaint (1=yes)"},
        "DXMPTR2":   {"Description": "MCI criterion 2: objective memory impairment (1=yes)"},
        "DXMPTR3":   {"Description": "MCI criterion 3: intact general cognitive function (1=yes)"},
        "DXMPTR4":   {"Description": "MCI criterion 4: intact activities of daily living (1=yes)"},
        "DXMPTR5":   {"Description": "MCI criterion 5: not demented (1=yes)"},
        "DXMPTR6":   {"Description": "MCI criterion 6: other MCI criterion (1=yes)"},
        "DXMDUE":    {"Description": "MCI: presumed aetiology (1=AD, 2=vascular, 3=other, 4=unknown)"},
        "DXMOTHET":  {"Description": "MCI other aetiology description (text)"},
        "DXDSEV":    {"Description": "Dementia severity rating (1=mild, 2=moderate, 3=severe)"},
        "DXDDUE":    {"Description": "Dementia: presumed aetiology (1=AD, 2=vascular, 3=other)"},
        "DXAD":      {"Description": "Diagnosis: meets criteria for Alzheimer's dementia (1=yes, -4=N/A)"},
        "DXAPP":     {"Description": "AD diagnosis confidence: probable vs possible"},
        "DXAPROB":   {"Description": "AD diagnosis: probable AD (NINCDS-ADRDA criteria; 1=yes)"},
        "DXAPOSS":   {"Description": "AD diagnosis: possible AD (NINCDS-ADRDA criteria; 1=yes)"},
        "DXPARK":    {"Description": "Diagnosis: Parkinson's disease (1=yes, -4=N/A)"},
        "DXPDES":    {"Description": "Parkinson's disease description / comments (text)"},
        "DXPCOG":    {"Description": "Parkinson's disease with cognitive involvement (1=yes)"},
        "DXPATYP":   {"Description": "Atypical parkinsonian syndrome (1=yes)"},
        "DXDEP":     {"Description": "Diagnosis: significant depression present (1=yes, -4=N/A)"},
        "DXOTHDEM":  {"Description": "Diagnosis: other dementia type (1=yes, -4=N/A)"},
        "DXODES":    {"Description": "Other dementia type description (text)"},
        "DXCONFID":  {"Description": "Clinician confidence in the diagnosis (1=uncertain, 2=probable, 3=definite)"},
        "session_id":SHARED["session_id"],
    },

    # ── ECog PT ────────────────────────────────────────────────────────────────
    "ecog_pt": {
        "VISCODE2":  SHARED["VISCODE2"],
        "VISDATE":   SHARED["VISDATE"],
        "CONCERN":   {"Description": "Does the participant have cognitive concerns compared to 10 years ago? (1=yes, 0=no)"},
        "MEMORY1":   {"Description": "ECog Memory item 1: remembering a short list of items"},
        "MEMORY2":   {"Description": "ECog Memory item 2: remembering things that happened recently"},
        "MEMORY3":   {"Description": "ECog Memory item 3: recalling conversations a few days later"},
        "MEMORY4":   {"Description": "ECog Memory item 4: remembering where things are usually kept"},
        "MEMORY5":   {"Description": "ECog Memory item 5: remembering where things have been placed"},
        "MEMORY6":   {"Description": "ECog Memory item 6: remembering to do things"},
        "MEMORY7":   {"Description": "ECog Memory item 7: learning how to use a new device"},
        "MEMORY8":   {"Description": "ECog Memory item 8: remembering things from the past"},
        "LANG1":     {"Description": "ECog Language item 1: forgetting words or names"},
        "LANG2":     {"Description": "ECog Language item 2: referring to objects by description rather than name"},
        "LANG3":     {"Description": "ECog Language item 3: speaking in complete sentences"},
        "LANG4":     {"Description": "ECog Language item 4: communicating thoughts coherently"},
        "LANG5":     {"Description": "ECog Language item 5: following the gist of a conversation"},
        "LANG6":     {"Description": "ECog Language item 6: understanding what others are saying"},
        "LANG7":     {"Description": "ECog Language item 7: reading and understanding books/newspapers"},
        "LANG8":     {"Description": "ECog Language item 8: writing without errors"},
        "LANG9":     {"Description": "ECog Language item 9: spelling"},
        "VISSPAT1":  {"Description": "ECog Visuospatial item 1: finding way in a familiar environment"},
        "VISSPAT2":  {"Description": "ECog Visuospatial item 2: reading a map and giving directions"},
        "VISSPAT3":  {"Description": "ECog Visuospatial item 3: parking without difficulty"},
        "VISSPAT4":  {"Description": "ECog Visuospatial item 4: navigating in a new place"},
        "VISSPAT5":  {"Description": "ECog Visuospatial item 5: judging distances and space"},
        "VISSPAT6":  {"Description": "ECog Visuospatial item 6: reading street signs"},
        "VISSPAT7":  {"Description": "ECog Visuospatial item 7: following instructions involving direction"},
        "VISSPAT8":  {"Description": "ECog Visuospatial item 8: spatial orientation"},
        "PLAN1":     {"Description": "ECog Planning item 1: planning a sequence of stops on an errand"},
        "PLAN2":     {"Description": "ECog Planning item 2: preparing a meal requiring multiple components"},
        "PLAN3":     {"Description": "ECog Planning item 3: developing a new activity or hobby"},
        "PLAN4":     {"Description": "ECog Planning item 4: adapting to a change in plans"},
        "PLAN5":     {"Description": "ECog Planning item 5: thinking ahead and planning ahead"},
        "ORGAN1":    {"Description": "ECog Organisation item 1: keeping living/working space tidy"},
        "ORGAN2":    {"Description": "ECog Organisation item 2: keeping financial records in order"},
        "ORGAN3":    {"Description": "ECog Organisation item 3: keeping important documents organised"},
        "ORGAN4":    {"Description": "ECog Organisation item 4: meeting commitments on time"},
        "ORGAN5":    {"Description": "ECog Organisation item 5: keeping track of different tasks"},
        "ORGAN6":    {"Description": "ECog Organisation item 6: managing multiple tasks simultaneously"},
        "DIVATT1":   {"Description": "ECog Divided Attention item 1: doing two things at once"},
        "DIVATT2":   {"Description": "ECog Divided Attention item 2: attending to one thing while another is occurring"},
        "DIVATT3":   {"Description": "ECog Divided Attention item 3: switching between tasks"},
        "DIVATT4":   {"Description": "ECog Divided Attention item 4: concentrating on a task"},
        "STAFFASST": {"Description": "Assessment was completed with staff assistance (1=yes)"},
        "VALIDITY":  {"Description": "Examiner rating of response validity (1=valid, 2=questionable, 3=invalid)"},
        "SOURCE":    SHARED["SOURCE"],
    },

    # ── ECog SP ────────────────────────────────────────────────────────────────
    "ecog_sp": {
        "VISCODE2":    SHARED["VISCODE2"],
        "VISDATE":     SHARED["VISDATE"],
        "MEMORY1":     {"Description": "ECog Memory item 1 (study-partner-rated): remembering a short list"},
        "MEMORY2":     {"Description": "ECog Memory item 2 (SP): remembering things that happened recently"},
        "MEMORY3":     {"Description": "ECog Memory item 3 (SP): recalling conversations a few days later"},
        "MEMORY4":     {"Description": "ECog Memory item 4 (SP): remembering where things are usually kept"},
        "MEMORY5":     {"Description": "ECog Memory item 5 (SP): remembering where things have been placed"},
        "MEMORY6":     {"Description": "ECog Memory item 6 (SP): remembering to do things"},
        "MEMORY7":     {"Description": "ECog Memory item 7 (SP): learning how to use a new device"},
        "MEMORY8":     {"Description": "ECog Memory item 8 (SP): remembering things from the past"},
        "LANG1":       {"Description": "ECog Language item 1 (SP): forgetting words or names"},
        "LANG2":       {"Description": "ECog Language item 2 (SP): referring to objects by description"},
        "LANG3":       {"Description": "ECog Language item 3 (SP): speaking in complete sentences"},
        "LANG4":       {"Description": "ECog Language item 4 (SP): communicating thoughts coherently"},
        "LANG5":       {"Description": "ECog Language item 5 (SP): following the gist of a conversation"},
        "LANG6":       {"Description": "ECog Language item 6 (SP): understanding what others are saying"},
        "LANG7":       {"Description": "ECog Language item 7 (SP): reading and understanding books"},
        "LANG8":       {"Description": "ECog Language item 8 (SP): writing without errors"},
        "LANG9":       {"Description": "ECog Language item 9 (SP): spelling"},
        "VISSPAT1":    {"Description": "ECog Visuospatial item 1 (SP): finding way in familiar environment"},
        "VISSPAT2":    {"Description": "ECog Visuospatial item 2 (SP): reading a map"},
        "VISSPAT3":    {"Description": "ECog Visuospatial item 3 (SP): parking without difficulty"},
        "VISSPAT4":    {"Description": "ECog Visuospatial item 4 (SP): navigating in a new place"},
        "VISSPAT5":    {"Description": "ECog Visuospatial item 5 (SP): judging distances"},
        "VISSPAT6":    {"Description": "ECog Visuospatial item 6 (SP): reading street signs"},
        "VISSPAT7":    {"Description": "ECog Visuospatial item 7 (SP): following directional instructions"},
        "VISSPAT8":    {"Description": "ECog Visuospatial item 8 (SP): spatial orientation"},
        "PLAN1":       {"Description": "ECog Planning item 1 (SP): planning errands"},
        "PLAN2":       {"Description": "ECog Planning item 2 (SP): preparing a multi-component meal"},
        "PLAN3":       {"Description": "ECog Planning item 3 (SP): developing a new hobby"},
        "PLAN4":       {"Description": "ECog Planning item 4 (SP): adapting to change in plans"},
        "PLAN5":       {"Description": "ECog Planning item 5 (SP): thinking and planning ahead"},
        "ORGAN1":      {"Description": "ECog Organisation item 1 (SP): keeping living space tidy"},
        "ORGAN2":      {"Description": "ECog Organisation item 2 (SP): keeping financial records in order"},
        "ORGAN3":      {"Description": "ECog Organisation item 3 (SP): keeping documents organised"},
        "ORGAN4":      {"Description": "ECog Organisation item 4 (SP): meeting commitments on time"},
        "ORGAN5":      {"Description": "ECog Organisation item 5 (SP): keeping track of tasks"},
        "ORGAN6":      {"Description": "ECog Organisation item 6 (SP): managing multiple tasks"},
        "DIVATT1":     {"Description": "ECog Divided Attention item 1 (SP): doing two things at once"},
        "DIVATT2":     {"Description": "ECog Divided Attention item 2 (SP): attending to one thing while another occurs"},
        "DIVATT3":     {"Description": "ECog Divided Attention item 3 (SP): switching between tasks"},
        "DIVATT4":     {"Description": "ECog Divided Attention item 4 (SP): concentrating on a task"},
        "SOURCE":      SHARED["SOURCE"],
        "EcogSPMem":   {"Description": "ECog Memory subscale total (study-partner-rated)"},
        "EcogSPLang":  {"Description": "ECog Language subscale total (study-partner-rated)"},
        "EcogSPVisspat":{"Description": "ECog Visuospatial subscale total (study-partner-rated)"},
        "EcogSPPlan":  {"Description": "ECog Planning subscale total (study-partner-rated)"},
        "EcogSPOrgan": {"Description": "ECog Organisation subscale total (study-partner-rated)"},
        "EcogSPDivatt":{"Description": "ECog Divided Attention subscale total (study-partner-rated)"},
        "EcogSPTotal": {"Description": "ECog grand total score (study-partner-rated)"},
    },

    # ── FAQ ────────────────────────────────────────────────────────────────────
    "faq": {
        "VISCODE2":  SHARED["VISCODE2"],
        "VISDATE":   SHARED["VISDATE"],
        "SOURCE":    SHARED["SOURCE"],
        "FAQFINAN":  {"Description": "FAQ item 1: Writing cheques, paying bills, balancing chequebook (0=normal, 1=never did/difficulty, 2=needs assistance, 3=dependent)"},
        "FAQFORM":   {"Description": "FAQ item 2: Assembling tax records, business affairs, or other papers (0–3)"},
        "FAQSHOP":   {"Description": "FAQ item 3: Shopping alone for clothes, household necessities, or groceries (0–3)"},
        "FAQGAME":   {"Description": "FAQ item 4: Playing a game of skill such as bridge or chess, working on a hobby (0–3)"},
        "FAQBEVG":   {"Description": "FAQ item 5: Heating water, making a cup of coffee, turning off the stove (0–3)"},
        "FAQMEAL":   {"Description": "FAQ item 6: Preparing a balanced meal (0–3)"},
        "FAQEVENT":  {"Description": "FAQ item 7: Keeping track of current events (0–3)"},
        "FAQTV":     {"Description": "FAQ item 8: Paying attention to and understanding a TV programme, book, or magazine (0–3)"},
        "FAQREM":    {"Description": "FAQ item 9: Remembering appointments, family occasions, holidays, medications (0–3)"},
        "FAQTRAVL":  {"Description": "FAQ item 10: Travelling out of the neighbourhood, driving, or using public transport (0–3)"},
        "SPID":      SHARED["SPID"],
    },

    # ── GDS ────────────────────────────────────────────────────────────────────
    "gds": {
        "VISCODE2":  SHARED["VISCODE2"],
        "VISDATE":   SHARED["VISDATE"],
        "SOURCE":    SHARED["SOURCE"],
        "GDUNABL":   {"Description": "GDS: participant was unable to complete the scale (1=yes)"},
        "GDSATIS":   {"Description": "GDS item 1: Are you basically satisfied with your life? (1=No)"},
        "GDDROP":    {"Description": "GDS item 2: Have you dropped many of your activities and interests? (1=Yes)"},
        "GDEMPTY":   {"Description": "GDS item 3: Do you feel that your life is empty? (1=Yes)"},
        "GDBORED":   {"Description": "GDS item 4: Do you often get bored? (1=Yes)"},
        "GDSPIRIT":  {"Description": "GDS item 5: Are you in good spirits most of the time? (1=No)"},
        "GDAFRAID":  {"Description": "GDS item 6: Are you afraid that something bad is going to happen to you? (1=Yes)"},
        "GDHAPPY":   {"Description": "GDS item 7: Do you feel happy most of the time? (1=No)"},
        "GDHELP":    {"Description": "GDS item 8: Do you often feel helpless? (1=Yes)"},
        "GDHOME":    {"Description": "GDS item 9: Do you prefer to stay at home rather than going out? (1=Yes)"},
        "GDMEMORY":  {"Description": "GDS item 10: Do you feel you have more problems with memory than most? (1=Yes)"},
        "GDALIVE":   {"Description": "GDS item 11: Do you think it is wonderful to be alive now? (1=No)"},
        "GDWORTH":   {"Description": "GDS item 12: Do you feel pretty worthless the way you are now? (1=Yes)"},
        "GDENERGY":  {"Description": "GDS item 13: Do you feel full of energy? (1=No)"},
        "GDHOPE":    {"Description": "GDS item 14: Do you feel that your situation is hopeless? (1=Yes)"},
        "GDBETTER":  {"Description": "GDS item 15: Do you think that most people are better off than you? (1=Yes)"},
    },

    # ── MMSE ────────────────────────────────────────────────────────────────────
    "mmse": {
        "VISCODE2":  SHARED["VISCODE2"],
        "VISDATE":   SHARED["VISDATE"],
        "DONE":      SHARED["DONE"],
        "NDREASON":  SHARED["NDREASON"],
        "SOURCE":    SHARED["SOURCE"],
        "MMDATE":    {"Description": "MMSE Orientation item: correct date (1=correct, 0=incorrect)"},
        "MMYEAR":    {"Description": "MMSE Orientation item: correct year (1=correct, 0=incorrect)"},
        "MMMONTH":   {"Description": "MMSE Orientation item: correct month (1=correct, 0=incorrect)"},
        "MMDAY":     {"Description": "MMSE Orientation item: correct day of week (1=correct, 0=incorrect)"},
        "MMSEASON":  {"Description": "MMSE Orientation item: correct season (1=correct, 0=incorrect)"},
        "MMHOSPIT":  {"Description": "MMSE Orientation item: correct hospital/place name (1=correct, 0=incorrect)"},
        "MMFLOOR":   {"Description": "MMSE Orientation item: correct floor (1=correct, 0=incorrect)"},
        "MMCITY":    {"Description": "MMSE Orientation item: correct city (1=correct, 0=incorrect)"},
        "MMAREA":    {"Description": "MMSE Orientation item: correct county/area (1=correct, 0=incorrect)"},
        "MMSTATE":   {"Description": "MMSE Orientation item: correct state/country (1=correct, 0=incorrect)"},
        "WORDLIST":  {"Description": "MMSE: word list version used for Registration and Recall items"},
        "WORD1":     {"Description": "MMSE Registration: word 1 recalled immediately (1=correct, 0=incorrect)"},
        "WORD2":     {"Description": "MMSE Registration: word 2 recalled immediately (1=correct, 0=incorrect)"},
        "WORD3":     {"Description": "MMSE Registration: word 3 recalled immediately (1=correct, 0=incorrect)"},
        "MMTRIALS":  {"Description": "MMSE Registration: number of trials needed to register all 3 words"},
        "MMD":       {"Description": "MMSE Serial 7s / WORLD backwards: score on 'D' letter (1=correct)"},
        "MML":       {"Description": "MMSE WORLD backwards: score on 'L' letter (1=correct)"},
        "MMR":       {"Description": "MMSE WORLD backwards: score on 'R' letter (1=correct)"},
        "MMO":       {"Description": "MMSE WORLD backwards: score on 'O' letter (1=correct)"},
        "MMW":       {"Description": "MMSE WORLD backwards: score on 'W' letter (1=correct)"},
        "MMLTR1":    {"Description": "MMSE Serial 7s: first subtraction correct (1=correct)"},
        "MMLTR2":    {"Description": "MMSE Serial 7s: second subtraction correct (1=correct)"},
        "MMLTR3":    {"Description": "MMSE Serial 7s: third subtraction correct (1=correct)"},
        "MMLTR4":    {"Description": "MMSE Serial 7s: fourth subtraction correct (1=correct)"},
        "MMLTR5":    {"Description": "MMSE Serial 7s: fifth subtraction correct (1=correct)"},
        "MMLTR6":    {"Description": "MMSE Serial 7s (supplementary): sixth item (1=correct)"},
        "MMLTR7":    {"Description": "MMSE Serial 7s (supplementary): seventh item (1=correct)"},
        "WORLDSCORE":{"Description": "MMSE Attention: WORLD-backwards total score (0–5)"},
        "WORD1DL":   {"Description": "MMSE Recall: word 1 recalled at delayed recall (1=correct, 0=incorrect)"},
        "WORD2DL":   {"Description": "MMSE Recall: word 2 recalled at delayed recall (1=correct, 0=incorrect)"},
        "WORD3DL":   {"Description": "MMSE Recall: word 3 recalled at delayed recall (1=correct, 0=incorrect)"},
        "MMWATCH":   {"Description": "MMSE Naming: name a watch (1=correct, 0=incorrect)"},
        "MMPENCIL":  {"Description": "MMSE Naming: name a pencil (1=correct, 0=incorrect)"},
        "MMREPEAT":  {"Description": "MMSE Repetition: repeat 'No ifs, ands or buts' (1=correct, 0=incorrect)"},
        "MMHAND":    {"Description": "MMSE 3-stage command step 1: take paper in right hand (1=correct)"},
        "MMFOLD":    {"Description": "MMSE 3-stage command step 2: fold paper in half (1=correct)"},
        "MMONFLR":   {"Description": "MMSE 3-stage command step 3: put paper on floor (1=correct)"},
        "MMREAD":    {"Description": "MMSE Reading: read and obey 'Close your eyes' (1=correct)"},
        "MMWRITE":   {"Description": "MMSE Writing: write a sentence spontaneously (1=correct)"},
        "MMDRAW":    {"Description": "MMSE Copying: copy intersecting pentagons (1=correct)"},
    },

    # ── NEUROBAT ────────────────────────────────────────────────────────────────
    "neurobat": {
        "VISCODE2":    SHARED["VISCODE2"],
        "VISDATE":     SHARED["VISDATE"],
        "SOURCE":      SHARED["SOURCE"],
        # Clock Copying Task
        "COPYCIRC":    {"Description": "Clock Copying Task: approximately circular face (1=correct, 0=incorrect)"},
        "COPYSYM":     {"Description": "Clock Copying Task: symmetry of number placement (1=correct, 0=incorrect)"},
        "COPYNUM":     {"Description": "Clock Copying Task: correctness of numbers (1=correct, 0=incorrect)"},
        "COPYHAND":    {"Description": "Clock Copying Task: presence of two hands (1=correct, 0=incorrect)"},
        "COPYTIME":    {"Description": "Clock Copying Task: hands set to correct time — ten past eleven (1=correct, 0=incorrect)"},
        "COPYSCOR":    {"Description": "Clock Copying Task: total score (0–5)"},
        # Logical Memory
        "LMSTORY":     {"Description": "Logical Memory I (immediate recall): story version used (A or B)"},
        "LIMMTOTAL":   {"Description": "Logical Memory I (immediate recall): total correct story units recalled"},
        "LIMMEND":     {"Description": "Logical Memory I: whether testing ended early (1=yes)"},
        # Digit Span
        "DSPANFOR":    {"Description": "Digit Span Forward: longest span length achieved"},
        "DSPANFLTH":   {"Description": "Digit Span Forward: total correct trials across all lengths"},
        "DSPANBAC":    {"Description": "Digit Span Backward: longest span length achieved"},
        "DSPANBLTH":   {"Description": "Digit Span Backward: total correct trials across all lengths"},
        # Category Fluency extras
        "CATANPERS":   {"Description": "Category Fluency (Animals): total perseverations (repeated words)"},
        "CATANINTR":   {"Description": "Category Fluency (Animals): total intrusion errors"},
        "CATVEGESC":   {"Description": "Category Fluency (Vegetables): total correct in 60 seconds"},
        "CATVGPERS":   {"Description": "Category Fluency (Vegetables): total perseverations"},
        "CATVGINTR":   {"Description": "Category Fluency (Vegetables): total intrusion errors"},
        # Digit Symbol
        "DIGITSCOR":   {"Description": "Digit Symbol Coding (WAIS-R): total correct substitutions in 90 seconds"},
        # Logical Memory II (delayed)
        "LDELBEGIN":   {"Description": "Logical Memory II (30-minute delayed recall): whether delayed testing began (1=yes)"},
        "LDELTOTAL":   {"Description": "Logical Memory II (delayed recall): total correct story units recalled"},
        "LDELCUE":     {"Description": "Logical Memory II: number of cues required for recall"},
        # Boston Naming Test
        "BNTND":       {"Description": "Boston Naming Test: number of items not done / skipped"},
        "BNTSPONT":    {"Description": "Boston Naming Test: number of spontaneous correct responses (no cue)"},
        "BNTSTIM":     {"Description": "Boston Naming Test: number of items requiring a stimulus cue"},
        "BNTCSTIM":    {"Description": "Boston Naming Test: number correct after stimulus cue"},
        "BNTPHON":     {"Description": "Boston Naming Test: number of items requiring a phonemic (letter) cue"},
        "BNTCPHON":    {"Description": "Boston Naming Test: number correct after phonemic cue"},
        "BNTTOTAL":    {"Description": "Boston Naming Test: total score (spontaneous correct + correct after cue)"},
        # ANART extras
        "ANARTND":     {"Description": "ANART: number of words not administered"},
        "ANART":       {"Description": "ANART: total number of words administered"},
        # MINT
        "MINTSEMCUE":  {"Description": "Multilingual Naming Test (MINT): number of items requiring semantic cue"},
        "MINTTOTAL":   {"Description": "Multilingual Naming Test (MINT): total correct responses"},
        "MINTUNCUED":  {"Description": "Multilingual Naming Test (MINT): number of correct uncued responses"},
    },

    # ── STUDY_ARM ────────────────────────────────────────────────────────────────
    "study_arm": {
        "PTNO":      {"Description": "Participant number within the site"},
        "TYPE":      {"Description": "Participant type / category within the study arm"},
        "ENROLLED":  {"Description": "Enrolment status (1=enrolled, 0=screen fail)"},
        "RANDDATE":  {"Description": "Date of randomisation or enrolment (YYYY-MM-DD)"},
        "MAPPDATE":  {"Description": "Date the participant was mapped to a study protocol arm (YYYY-MM-DD)"},
    },
}

# ── Apply shared columns to every table ───────────────────────────────────────
for table_name, col_dict in TABLE_COLS.items():
    for shared_col, shared_val in SHARED.items():
        if shared_col not in col_dict:
            col_dict[shared_col] = shared_val

# ── Patch each JSON sidecar ────────────────────────────────────────────────────
patched_tables = 0
total_added = 0

for table_name, new_cols in TABLE_COLS.items():
    json_path = os.path.join(PHENO_DIR, f"{table_name}.json")
    tsv_path  = os.path.join(PHENO_DIR, f"{table_name}.tsv")
    if not os.path.isfile(json_path):
        print(f"SKIP (no JSON): {table_name}")
        continue

    # Load existing JSON
    with open(json_path) as f:
        j = json.load(f)

    # Get actual TSV columns
    tsv_cols = set(pd.read_csv(tsv_path, sep='\t', nrows=1).columns)

    added = 0
    for col, desc in new_cols.items():
        if col in tsv_cols and col not in j:
            j[col] = desc
            added += 1

    # Write back
    with open(json_path, "w") as f:
        json.dump(j, f, indent=2)

    print(f"{table_name}: +{added} descriptions added")
    total_added += added
    patched_tables += 1

print(f"\nDone. {total_added} descriptions added across {patched_tables} tables.")

# ── Final check: remaining undocumented columns ────────────────────────────────
print("\nFinal undocumented column check:")
grand_total = 0
for fname in sorted(os.listdir(PHENO_DIR)):
    if not fname.endswith('.tsv'):
        continue
    t = fname[:-4]
    tsv_cols = set(pd.read_csv(os.path.join(PHENO_DIR, fname), sep='\t', nrows=1).columns)
    jpath = os.path.join(PHENO_DIR, f"{t}.json")
    if not os.path.isfile(jpath):
        continue
    with open(jpath) as f:
        jcols = set(json.load(f).keys()) - {'MeasurementToolMetadata'}
    missing = tsv_cols - jcols
    grand_total += len(missing)
    if missing:
        print(f"  {t}: {len(missing)} remaining — {list(missing)[:5]}")
    else:
        print(f"  {t}: OK")
print(f"\nTotal remaining undocumented: {grand_total}")
