"""
build_clinical_csv.py
Builds a longitudinal master CSV for ADNI subjects with both MRI and SNP data.
Outputs: master CSV, per-visit CSVs, and baseline train/val/test splits.
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
SUBJECTS_TSV   = r"D:\ADNI_BIDS_project\bids\genotype\subjects_with_snp_and_mri.tsv"
CLINICAL_DIR   = r"D:\ADNI_BIDS_project\sourcedata\clinical"
OUTPUT_DIR     = r"D:\ADNI_BIDS_project\derivatives\clinical"
IMAGING_META   = r"D:\ADNI_BIDS_project\bids\imaging\mri_visit_dates_snp_subjects.csv"

Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

# ── Load subject list ──────────────────────────────────────────────────────
subjects_df = pd.read_csv(SUBJECTS_TSV, sep="\t")
# adni_subject_id format: 011_S_0002
valid_ptids = set(subjects_df["adni_subject_id"].str.strip())

# ── Month-label vocabulary from imaging metadata ──────────────────────────
# We don't use imaging dates as the visit backbone — clinical assessments
# happen at more timepoints than MRI scans.  We only read the imaging CSV to
# obtain the ordered list of recognised month labels (bl, m12, m24, …) so
# that per-visit output CSVs use the same naming convention.
imaging_wide = pd.read_csv(IMAGING_META, nrows=0)  # header only
VISIT_COLS = [c for c in imaging_wide.columns
              if c not in ("participant_id",)]
print(f"Month-label vocabulary from imaging metadata: {VISIT_COLS}")

def load(fname, **kwargs):
    path = os.path.join(CLINICAL_DIR, fname)
    df = pd.read_csv(path, low_memory=False, **kwargs)
    # Normalise PTID column name
    for col in ["PTID", "ptid"]:
        if col in df.columns:
            df["PTID"] = df[col].str.strip()
            break
    return df

# ── Demographics (baseline only, keyed by PTID) ───────────────────────────
demo_raw = load("PTDEMOG.csv")
demo_raw["_visdate"] = pd.to_datetime(demo_raw.get("VISDATE", None), errors="coerce")
demo = demo_raw.sort_values("_visdate").drop_duplicates("PTID", keep="first")

# ADNIMERGE -- demographic fallback + per-visit DX labels
adnimerge = load("ADNIMERGE.csv")
adnimerge["PTID"]     = adnimerge["PTID"].astype(str).str.strip()
adnimerge["EXAMDATE"] = pd.to_datetime(adnimerge["EXAMDATE"], errors="coerce")
adnimerge["VISCODE"]  = adnimerge["VISCODE"].astype(str).str.strip().str.lower()
_adni_bl = adnimerge[adnimerge["VISCODE"] == "bl"].drop_duplicates("PTID")
dx_bl_dict   = dict(zip(_adni_bl["PTID"], _adni_bl["DX_bl"]))  # PTID -> CN/EMCI/LMCI/AD
adni_bl_demo = _adni_bl.set_index("PTID")                       # demographic fallback
adnimerge_dx = adnimerge[["PTID","EXAMDATE","DX"]].dropna(subset=["DX"]).copy()

def get_visit_granular(ptid, visit_date):
    """Per-visit granular diagnosis: CN / EMCI / LMCI / MCI / AD."""
    sub = adnimerge_dx[adnimerge_dx["PTID"] == ptid]
    if sub.empty or visit_date is None or pd.isna(visit_date): return np.nan
    idx    = (sub["EXAMDATE"] - visit_date).abs().idxmin()
    dx_val = str(sub.loc[idx, "DX"]).strip()
    if dx_val in ("nan", ""): return np.nan
    if dx_val == "Dementia":  return "AD"
    if dx_val == "MCI":
        bl_val = str(dx_bl_dict.get(ptid, "")).strip()
        return bl_val if bl_val in ("EMCI", "LMCI") else "MCI"
    return dx_val

def gender_str(v):
    if pd.isna(v): return "NaN"
    return str(v)

def dob_year(v):
    """PTDOB is MM/YYYY or PTDOBYY is integer year."""
    if pd.notna(v):
        try:
            return str(int(v))
        except Exception:
            pass
    return "NaN"

def yn(v):
    if pd.isna(v): return "NaN"
    sv = str(v).strip()
    if sv in ("1","Yes","yes","YES"): return "Yes"
    if sv in ("0","No","no","NO"):  return "No"
    return sv

# Build demographics dict; fill NaNs from ADNIMERGE baseline row
demo_dict = {}
for _, row in demo.iterrows():
    pid  = row.get("PTID", "")
    adni = adni_bl_demo.loc[pid] if pid in adni_bl_demo.index else pd.Series(dtype=object)
    def _fb(col, fb_col=None):
        v = row.get(col, np.nan)
        if pd.isna(v) and fb_col is not None:
            try: v = adni[fb_col]
            except KeyError: pass
        return v
    demo_dict[pid] = {
        "sex":       gender_str(_fb("PTGENDER", "PTGENDER")),
        "dob_year":  dob_year(_fb("PTDOBYY")),
        "hand":      _fb("PTHAND"),
        "marital":   _fb("PTMARRY"),
        "education": _fb("PTEDUCAT", "PTEDUCAT"),
        "retired":   yn(_fb("PTNOTRT")),
        "residence": _fb("PTHOME"),
        "language":  _fb("PTPLANG"),
        "ethnicity": _fb("PTRACCAT"),
    }

# ── APOE (raw values for structured formatting) ───────────────────────────
apoe_df = load("APOERES.csv")
apoe_df = apoe_df.drop_duplicates("PTID")

def _parse_apoe(genotype):
    if pd.isna(genotype):
        return {"genotype": "NaN", "e4_count": "NaN", "status": "NaN"}
    g  = str(genotype).strip()
    e4 = g.count("4")
    status = "non-carrier" if e4 == 0 else ("heterozygous" if e4 == 1 else "homozygous")
    return {"genotype": g, "e4_count": e4, "status": status}

apoe_raw_dict = {r["PTID"]: _parse_apoe(r.get("GENOTYPE")) for _, r in apoe_df.iterrows()}

# ── CSF baseline biomarkers (ADNIMERGE _bl columns, carried forward) ──────
_csf_us = ["ABETA_bl", "TAU_bl",  "PTAU_bl"]
_csf_dt = ["ABETA.bl", "TAU.bl",  "PTAU.bl"]
csf_bl_dict = {}
for _ptid in adni_bl_demo.index:
    _row = adni_bl_demo.loc[_ptid]
    def _csf(us, dt):
        for c in [us, dt]:
            if c in _row.index and pd.notna(_row[c]): return _row[c]
        return np.nan
    csf_bl_dict[_ptid] = {
        "abeta": _csf(_csf_us[0], _csf_dt[0]),
        "tau":   _csf(_csf_us[1], _csf_dt[1]),
        "ptau":  _csf(_csf_us[2], _csf_dt[2]),
    }
_n_csf = sum(pd.notna(v["abeta"]) for v in csf_bl_dict.values())
print(f"CSF baseline Abeta present for {_n_csf} of {len(csf_bl_dict)} subjects.")


# ── Diagnosis ─────────────────────────────────────────────────────────────
dx_df = load("DXSUM_PDXCONV_ADNIALL.csv")
dx_df = dx_df[dx_df["PTID"].isin(valid_ptids)].copy()

# ── SMC exclusion ─────────────────────────────────────────────────────────
# SMC (Subjective Memory Concern) is an ADNI3-specific arm.
# Check REGISTRY for ARM column; exclude any PTID whose *only* arm is SMC.
try:
    registry = load("REGISTRY.csv")
    if "ARM" in registry.columns:
        smc_ptids = set(
            registry.groupby("PTID")["ARM"]
            .apply(lambda arms: all(str(a).strip().upper() == "SMC" for a in arms))
            .pipe(lambda s: s[s].index)
        )
        if smc_ptids:
            print(f"Excluding {len(smc_ptids)} SMC-only subjects: {smc_ptids}")
            valid_ptids -= smc_ptids
            dx_df = dx_df[dx_df["PTID"].isin(valid_ptids)]
        else:
            print("No SMC-only subjects found in cohort.")
    else:
        print("WARNING: ARM column not found in REGISTRY.csv — cannot filter SMC.")
except FileNotFoundError:
    print("WARNING: REGISTRY.csv not found — cannot filter SMC automatically. "
          "Please verify cohort manually.")

# ── DXCHANGE → current diagnosis mapping ──────────────────────────────────
# DXCHANGE encodes the *transition*; we derive current dx from it:
#   1, 7, 9  →  CN   (stable NL / MCI→NL / Dementia→NL)
#   2, 4, 8  →  MCI  (stable MCI / NL→MCI / Dementia→MCI)
#   3, 5, 6  →  AD   (stable Dementia / MCI→AD / NL→AD)
# ADNI3 uses DIAGNOSIS (1=CN, 2=MCI, 3=Dementia) where DXCHANGE may be absent.
def map_diagnosis(dxchange, diagnosis=np.nan):
    """Returns (binary_ad, multi_label) collapsing DXCHANGE to current dx.
    Falls back to ADNI3 DIAGNOSIS column when DXCHANGE is missing.
    """
    if pd.notna(dxchange):
        v = float(dxchange)
        if v in (1.0, 7.0, 9.0): return (0, "CN")
        if v in (2.0, 4.0, 8.0): return (0, "MCI")
        if v in (3.0, 5.0, 6.0): return (1, "AD")
    # Fallback: ADNI3 DIAGNOSIS
    if pd.notna(diagnosis):
        v = float(diagnosis)
        if v == 1.0: return (0, "CN")
        if v == 2.0: return (0, "MCI")
        if v == 3.0: return (1, "AD")
    return (np.nan, np.nan)

dx_df["_diag_bin"], dx_df["_diag_multi"] = zip(*dx_df.apply(
    lambda r: map_diagnosis(r.get("DXCHANGE"), r.get("DIAGNOSIS")), axis=1))
dx_df["EXAMDATE"] = pd.to_datetime(dx_df["EXAMDATE"], errors="coerce")
dx_df = dx_df.sort_values("EXAMDATE")

# ── Plasma biomarkers (UPenn) ──────────────────────────────────────────────
plasma = load("UPENN_PLASMA_FUJIREBIO_QUANTERIX.csv")
plasma = plasma[plasma["PTID"].isin(valid_ptids)].copy()
plasma["EXAMDATE"] = pd.to_datetime(plasma.get("EXAMDATE", plasma.get("VISDATE","")), errors="coerce")

# ── pTau181 (UGot) ────────────────────────────────────────────────────────
ptau = pd.read_csv(os.path.join(CLINICAL_DIR, "UGOTPTAU181.csv"), low_memory=False)
ptau["RID"] = ptau["RID"].astype(str).str.strip()
ptau["EXAMDATE"] = pd.to_datetime(ptau["EXAMDATE"], errors="coerce")

# Build RID -> PTID mapping from subjects_df
rid_to_ptid = {}
for _, row in subjects_df.iterrows():
    ptid = row["adni_subject_id"].strip()  # 011_S_0002
    rid  = ptid.split("_S_")[-1].lstrip("0") or "0"
    rid_to_ptid[rid] = ptid
    rid_to_ptid[ptid.split("_S_")[-1]] = ptid  # with leading zeros too

ptau["PTID"] = ptau["RID"].map(lambda r: rid_to_ptid.get(r.lstrip("0"), rid_to_ptid.get(r, np.nan)))
ptau = ptau[ptau["PTID"].notna()]

# ── Cognitive tests ───────────────────────────────────────────────────────
adas   = load("ADAS.csv")
cdr    = load("CDR.csv")
mmse   = load("MMSE.csv")
moca   = load("MOCA.csv")
neuro  = load("NEUROBAT.csv")
faq    = load("FAQ.csv")
gds    = load("GDSCALE.csv")
hach   = load("MODHACH.csv")
not_moca = [adas, cdr, mmse, neuro, faq, gds, hach]
for df in not_moca:
    df["PTID"] = df["PTID"].astype(str).str.strip()
    for dcol in ["VISDATE","EXAMDATE"]:
        if dcol in df.columns:
            df[dcol] = pd.to_datetime(df[dcol], errors="coerce")
# MoCA: date-parse only (score is calculated from items, not the MOCA column)
moca["PTID"] = moca["PTID"].astype(str).str.strip()
if "VISDATE" in moca.columns:
    moca["VISDATE"] = pd.to_datetime(moca["VISDATE"], errors="coerce")

# ── MoCA score calculation from item columns ────────────────────────────
# Items are encoded as 'Correct'/'Incorrect' strings.
# 27 binary items (1 pt each) + Serial 7s (SERIAL1-5, 0-3 pts).
# Immediate recall (IMMT*) = registration only, does NOT score.
_MOCA_BINARY = [
    "TRAILS","CUBE","CLOCKCON","CLOCKNO","CLOCKHAN",
    "LION","RHINO","CAMEL",
    "DIGFOR","DIGBACK","LETTERS",
    "REPEAT1","REPEAT2","FFLUENCY",
    "ABSTRAN","ABSMEAS",
    "DELW1","DELW2","DELW3","DELW4","DELW5",
    "DATE","MONTH","YEAR","DAY","PLACE","CITY",
]
_MOCA_SERIAL = ["SERIAL1","SERIAL2","SERIAL3","SERIAL4","SERIAL5"]

def calc_moca(row):
    """Compute MoCA total (0-30) from item-level Correct/Incorrect columns."""
    if row is None: return np.nan
    score, any_data = 0, False
    for col in _MOCA_BINARY:
        v = str(row.get(col, "")).strip().lower()
        if v in ("correct", "incorrect"):
            any_data = True
            if v == "correct": score += 1
    n_serial = 0
    for col in _MOCA_SERIAL:
        v = str(row.get(col, "")).strip().lower()
        if v in ("correct", "incorrect"):
            any_data = True
            if v == "correct": n_serial += 1
    if   n_serial >= 4: score += 3
    elif n_serial >= 2: score += 2
    elif n_serial == 1: score += 1
    return score if any_data else np.nan

def fmtv(v, decimals=1):
    """Format a numeric value; return 'NaN' string if missing."""
    if pd.isna(v): return "NaN"
    try:
        f = float(v)
        return f"{f:.{decimals}f}"
    except Exception:
        return str(v)

def fmts(v):
    """Format string value."""
    if pd.isna(v): return "NaN"
    return str(v).strip()

def _faq_num(val):
    """Extract numeric score from FAQ value like 'Normal (0)' → 0, or pass through integers."""
    import re as _re
    if pd.isna(val): return np.nan
    s = str(val).strip()
    m = _re.search(r'\((\d+)\)', s)   # match 'Normal (0)' → '0'
    if m: return int(m.group(1))
    try: return int(float(s))           # already a plain number
    except: return np.nan

# ── FAQ item labels ────────────────────────────────────────────────────────
FAQ_LABELS = [
    ("FAQFINAN",  "Writing checks, paying bills, or balancing checkbook"),
    ("FAQFORM",   "Assembling tax records, business affairs, or other papers"),
    ("FAQSHOP",   "Shopping alone for clothes, household necessities, or groceries"),
    ("FAQGAME",   "Playing a game of skill such as bridge or chess, working on a hobby"),
    ("FAQBEVG",   "Heating water, making a cup of coffee, turning off the stove"),
    ("FAQMEAL",   "Preparing a balanced meal"),
    ("FAQEVENT",  "Keeping track of current events"),
    ("FAQTV",     "Paying attention to and understanding a TV program, book, or magazine"),
    ("FAQREM",    "Remembering appointments, family occasions, holidays, medications"),
    ("FAQTRAVL",  "Traveling out of the neighborhood, driving, or arranging to take public transportation"),
]

def closest_row(df, ptid, visit_date, date_col="VISDATE", viscode=None, viscode_col="VISCODE2"):
    """Return the row closest in date to visit_date for the given ptid."""
    sub = df[df["PTID"] == ptid].copy()
    if sub.empty: return None
    if viscode is not None and viscode_col in sub.columns:
        vs = sub[sub[viscode_col] == viscode]
        if not vs.empty: return vs.iloc[0]
    if date_col in sub.columns and visit_date is not None and pd.notna(visit_date):
        sub = sub.dropna(subset=[date_col])
        if sub.empty: return None
        idx = (sub[date_col] - visit_date).abs().idxmin()
        return sub.loc[idx]
    return sub.iloc[0] if not sub.empty else None

def extract_data(ptid, visit_date, viscode, d_row):
    """Extract all data for one visit. Returns a flat dict of raw/formatted values."""
    dem  = demo_dict.get(ptid, {})
    apoe = apoe_raw_dict.get(ptid, {"genotype": "NaN", "e4_count": "NaN", "status": "NaN"})
    csf  = csf_bl_dict.get(ptid, {"abeta": np.nan, "tau": np.nan, "ptau": np.nan})

    # ── Plasma ────────────────────────────────────────────────────────────
    pr    = closest_row(plasma, ptid, visit_date, date_col="EXAMDATE", viscode=None)
    pt181_rows = ptau[ptau["PTID"] == ptid].dropna(subset=["EXAMDATE"])
    if not pt181_rows.empty and visit_date is not None and pd.notna(visit_date):
        idx   = (pt181_rows["EXAMDATE"] - visit_date).abs().idxmin()
        pt181_val = pt181_rows.loc[idx, "PLASMAPTAU181"]
    else:
        pt181_val = np.nan

    # ── GDS ───────────────────────────────────────────────────────────────
    gr      = closest_row(gds,  ptid, visit_date, viscode=viscode)
    gds_val = gr.get("GDTOTAL") if gr is not None else np.nan
    try:    dep_bin = 1 if float(gds_val) >= 5 else 0
    except: dep_bin = np.nan

    # ── Cognitive ─────────────────────────────────────────────────────────
    ar   = closest_row(adas, ptid, visit_date, date_col="VISDATE", viscode=viscode)
    cr   = closest_row(cdr,  ptid, visit_date, date_col="VISDATE", viscode=viscode)
    mr   = closest_row(mmse, ptid, visit_date, date_col="VISDATE", viscode=viscode)
    mocr = closest_row(moca, ptid, visit_date, date_col="VISDATE", viscode=viscode)
    nr   = closest_row(neuro,ptid, visit_date, date_col="VISDATE", viscode=viscode)
    fr   = closest_row(faq,  ptid, visit_date, date_col="VISDATE", viscode=viscode)

    def nget(col): return nr.get(col, np.nan) if nr is not None else np.nan
    def fget(col): return fr.get(col, np.nan) if fr is not None else np.nan

    return {
        # Demographics (static, carried forward)
        "Sex":              fmts(dem.get("sex")),
        "YOB":              fmts(dem.get("dob_year")),
        "Handedness":       fmts(dem.get("hand")),
        "Marital_Status":   fmts(dem.get("marital")),
        "Education_Years":  fmtv(dem.get("education"), 0),
        "Retired":          fmts(dem.get("retired")),
        "Residence":        fmts(dem.get("residence")),
        "Language":         fmts(dem.get("language")),
        "Ethnicity":        fmts(dem.get("ethnicity")),
        # APOE (static)
        "APOE4_Status":     apoe["status"],
        "APOE_Alleles":     apoe["genotype"],
        "APOE4_Dosage":     str(apoe["e4_count"]),
        # CSF (static baseline)
        "CSF_Abeta_bl":     fmtv(csf["abeta"]),
        "CSF_Tau_bl":       fmtv(csf["tau"]),
        "CSF_pTau_bl":      fmtv(csf["ptau"]),
        # Plasma (longitudinal)
        "Plasma_Abeta42":        fmtv(pr.get("AB42_F")      if pr is not None else np.nan),
        "Plasma_Abeta40":        fmtv(pr.get("AB40_F")      if pr is not None else np.nan),
        "Plasma_Abeta42_40":     fmtv(pr.get("AB42_AB40_F") if pr is not None else np.nan),
        "Plasma_pTau217":        fmtv(pr.get("pT217_F")     if pr is not None else np.nan),
        "Plasma_pTau181":        fmtv(pt181_val),
        "Plasma_NfL":            fmtv(pr.get("NfL_Q")       if pr is not None else np.nan),
        "Plasma_GFAP":           fmtv(pr.get("GFAP_Q")      if pr is not None else np.nan),
        # Exclusions / screening
        "GDS_Total":        fmtv(gds_val, 0),
        "Depression_Binary": dep_bin,
        "Hachinski_Score":  fmtv(hach_bl_dict.get(ptid, np.nan), 0),
        # Cognitive tests
        "ADAS_Total13":     fmtv(ar.get("TOTAL13") if ar is not None else np.nan),
        "CDR_Global":       fmtv(cr.get("CDGLOBAL") if cr is not None else np.nan),
        "CDR_SumBoxes":     fmtv(cr.get("CDRSB")    if cr is not None else np.nan),
        "MMSE_Total":       fmtv(mr.get("MMSCORE")  if mr is not None else np.nan, 0),
        "MoCA_Total":       fmtv(calc_moca(mocr), 0),
        # NEUROBAT
        "CDT_Total":        fmtv(nget("CLOCKSCOR")),
        "CategoryFluency_Animals_Total": fmtv(nget("CATANIMSC")),
        "CategoryFluency_Animals_Pers":  fmtv(nget("CATANPERS")),
        "CategoryFluency_Animals_Intr":  fmtv(nget("CATANINTR")),
        "LogicalMemory_I_Total":  fmtv(nget("LIMMTOTAL")),
        "LogicalMemory_II_Total": fmtv(nget("LDELTOTAL")),
        "LogicalMemory_II_Cued":  fmts(nget("LDELCUE")),
        "AVLT_T1":          fmtv(nget("AVTOT1")),  "AVLT_T1_Err": fmtv(nget("AVERR1")),
        "AVLT_T2":          fmtv(nget("AVTOT2")),  "AVLT_T2_Err": fmtv(nget("AVERR2")),
        "AVLT_T3":          fmtv(nget("AVTOT3")),  "AVLT_T3_Err": fmtv(nget("AVERR3")),
        "AVLT_T4":          fmtv(nget("AVTOT4")),  "AVLT_T4_Err": fmtv(nget("AVERR4")),
        "AVLT_T5":          fmtv(nget("AVTOT5")),  "AVLT_T5_Err": fmtv(nget("AVERR5")),
        "AVLT_T6":          fmtv(nget("AVTOT6")),  "AVLT_T6_Err": fmtv(nget("AVERR6")),
        "AVLT_ListB":       fmtv(nget("AVTOTB")),  "AVLT_ListB_Err":  fmtv(nget("AVERRB")),
        "AVLT_Delay30":     fmtv(nget("AVDEL30MIN")), "AVLT_Delay30_Err": fmtv(nget("AVDELERR1")),
        "AVLT_Recognition": fmtv(nget("AVDELTOT")),   "AVLT_Recognition_Err": fmtv(nget("AVDELERR2")),
        "TrailsA_Time":   fmtv(nget("TRAASCOR")),  "TrailsA_ErrCom": fmtv(nget("TRAAERRCOM")), "TrailsA_ErrOm": fmtv(nget("TRAAERROM")),
        "TrailsB_Time":   fmtv(nget("TRABSCOR")),  "TrailsB_ErrCom": fmtv(nget("TRABERRCOM")), "TrailsB_ErrOm": fmtv(nget("TRABERROM")),
        # FAQ — store numeric score only (extracted from strings like 'Normal (0)' → 0)
        "FAQFINAN": _faq_num(fget("FAQFINAN")), "FAQFORM":  _faq_num(fget("FAQFORM")),
        "FAQSHOP":  _faq_num(fget("FAQSHOP")),  "FAQGAME":  _faq_num(fget("FAQGAME")),
        "FAQBEVG":  _faq_num(fget("FAQBEVG")),  "FAQMEAL":  _faq_num(fget("FAQMEAL")),
        "FAQEVENT": _faq_num(fget("FAQEVENT")), "FAQTV":    _faq_num(fget("FAQTV")),
        "FAQREM":   _faq_num(fget("FAQREM")),   "FAQTRAVL": _faq_num(fget("FAQTRAVL")),
        "FAQ_Total": fmtv(fget("FAQTOTAL")),
    }


def format_summarised_text(d):
    """Format a Summarised_Text string from the extract_data dict."""
    lines = []
    lines.append(
        f"Demographics - Sex: {d['Sex']}; YOB: {d['YOB']}; Handedness: {d['Handedness']}; "
        f"Marital status: {d['Marital_Status']}; Education years: {d['Education_Years']}; "
        f"Retired: {d['Retired']}; Residence: {d['Residence']}; "
        f"Primary Language: {d['Language']}; Ethnicity: {d['Ethnicity']}."
    )
    lines.append(
        f"Genotype - APOE4 allele status: {d['APOE4_Status']}; "
        f"alleles: epsilon {d['APOE_Alleles']}; dosage of epsilon 4: {d['APOE4_Dosage']}/2."
    )
    lines.append(
        f"Plasma biomarkers - Aβ42: {d['Plasma_Abeta42']}; Aβ40: {d['Plasma_Abeta40']}; "
        f"Aβ42/40: {d['Plasma_Abeta42_40']}; p-tau217: {d['Plasma_pTau217']}; "
        f"p-tau181: {d['Plasma_pTau181']}; NfL: {d['Plasma_NfL']}; GFAP: {d['Plasma_GFAP']}."
    )
    lines.append(
        f"CSF biomarkers at baseline - Abeta: {d['CSF_Abeta_bl']}; "
        f"tau: {d['CSF_Tau_bl']}; p-tau: {d['CSF_pTau_bl']}."
    )
    lines.append(
        f"Exclusions - Geriatric Depression Scale (GDS) score: {d['GDS_Total']}; "
        f"Depression: {'Yes' if d['Depression_Binary'] == 1 else ('No' if d['Depression_Binary'] == 0 else 'NaN')}; "
        f"Hachinski Ischaemic Score: {d['Hachinski_Score']}."
    )
    lines.append("Cognitive tests -")
    lines.append(
        f"Alzheimer's Disease Assessment Scale-Cog13 total: {d['ADAS_Total13']}; "
        f"Clinical Dementia Rating (CDR) scale global: {d['CDR_Global']}, CDR sum of boxes: {d['CDR_SumBoxes']}."
    )
    lines.append(f"Mini Mental State Examination (MMSE) total: {d['MMSE_Total']}.")
    lines.append(f"Montreal Cognitive Assessment (MoCA) total: {d['MoCA_Total']}.")
    lines.append(
        f"Category Fluency Test Animals total: {d['CategoryFluency_Animals_Total']}, "
        f"Perseverations: {d['CategoryFluency_Animals_Pers']}, "
        f"Intrusions: {d['CategoryFluency_Animals_Intr']}."
    )
    lines.append(
        f"Logical Memory Test Immediate Recall total: {d['LogicalMemory_I_Total']}; "
        f"Delayed Recall total: {d['LogicalMemory_II_Total']}; "
        f"Delayed Recall Cued: {d['LogicalMemory_II_Cued']}."
    )
    lines.append(
        f"Auditory Verbal Learning Test (AVLT) "
        f"Trial 1: Total {d['AVLT_T1']}, Intrusions {d['AVLT_T1_Err']}; "
        f"Trial 2: Total {d['AVLT_T2']}, Intrusions {d['AVLT_T2_Err']}; "
        f"Trial 3: Total {d['AVLT_T3']}, Intrusions {d['AVLT_T3_Err']}; "
        f"Trial 4: Total {d['AVLT_T4']}, Intrusions {d['AVLT_T4_Err']}; "
        f"Trial 5: Total {d['AVLT_T5']}, Intrusions {d['AVLT_T5_Err']}; "
        f"Trial 6: Total {d['AVLT_T6']}, Intrusions {d['AVLT_T6_Err']}; "
        f"List B Total {d['AVLT_ListB']}, Intrusions {d['AVLT_ListB_Err']}; "
        f"30 Minute Delay Total {d['AVLT_Delay30']}, Intrusions {d['AVLT_Delay30_Err']}; "
        f"Recognition Score {d['AVLT_Recognition']}, Intrusions {d['AVLT_Recognition_Err']}."
    )
    lines.append(
        f"Trail Making Test Part A: Time to Complete {d['TrailsA_Time']} seconds, "
        f"Errors of Commission {d['TrailsA_ErrCom']}, Errors of Omission {d['TrailsA_ErrOm']}; "
        f"Part B: Time to complete {d['TrailsB_Time']} seconds, "
        f"Errors of Commission {d['TrailsB_ErrCom']}, Errors of Omission {d['TrailsB_ErrOm']}."
    )
    faq_items = "\n".join(
        f"{label}: {d.get('FAQ' + col.replace('FAQ', ''), 'NaN')}."
        for col, label in FAQ_LABELS
    )
    lines.append(f"Functional Activities Questionnaire (FAQ) -\n{faq_items}\nFAQ total: {d['FAQ_Total']}.")
    return "\n".join(lines)



    # ── Demographics block ─────────────────────────────────────────────
    sex        = fmts(dem.get("sex", np.nan))
    yob        = fmts(dem.get("dob_year", np.nan))
    hand       = fmts(dem.get("hand", np.nan))
    marital    = fmts(dem.get("marital", np.nan))
    edu        = fmtv(dem.get("education", np.nan), 0)
    retired    = fmts(dem.get("retired", np.nan))
    residence  = fmts(dem.get("residence", np.nan))
    language   = fmts(dem.get("language", np.nan))
    ethnicity  = fmts(dem.get("ethnicity", np.nan))

    demo_block = (
        f"The participant's sex is {sex}. "
        f"Their Year of Birth is {yob}. "
        f"Their Handedness is {hand}. "
        f"Their Marital status at baseline is {marital}. "
        f"Their education in years is {edu}. "
        f"Participant Retired? {retired}. "
        f"Type of Participant residence {residence}. "
        f"Participant's Primary Language {language}. "
        f"The participant's Ethnicity is {ethnicity}."
    )

    # ── APOE block ────────────────────────────────────────────────────
    apoe_block = apoe_dict.get(ptid, "NaN")

    # ── Plasma biomarkers ─────────────────────────────────────────────
    pr = closest_row(plasma, ptid, visit_date, date_col="EXAMDATE", viscode=None)
    ab42   = fmtv(pr.get("AB42_F")  if pr is not None else np.nan)
    ab40   = fmtv(pr.get("AB40_F")  if pr is not None else np.nan)
    ratio  = fmtv(pr.get("AB42_AB40_F") if pr is not None else np.nan)
    pt217  = fmtv(pr.get("pT217_F") if pr is not None else np.nan)
    nfl    = fmtv(pr.get("NfL_Q")   if pr is not None else np.nan)
    gfap   = fmtv(pr.get("GFAP_Q")  if pr is not None else np.nan)

    # pTau181 from UGot (matched by RID)
    rid_str = ptid.split("_S_")[-1]
    pt181_rows = ptau[(ptau["PTID"] == ptid)]
    if not pt181_rows.empty and visit_date is not None and pd.notna(visit_date):
        pt181_rows = pt181_rows.dropna(subset=["EXAMDATE"])
        if not pt181_rows.empty:
            idx = (pt181_rows["EXAMDATE"] - visit_date).abs().idxmin()
            pt181 = fmtv(pt181_rows.loc[idx, "PLASMAPTAU181"])
        else:
            pt181 = "NaN"
    else:
        pt181 = "NaN"

    plasma_block = (
        f"Plasma Aβ42 is {ab42}, Aβ40 is {ab40}, Aβ42/40 ratio is {ratio}. "
        f"Plasma p-tau217 is {pt217}, p-tau181 is {pt181}. "
        f"Plasma NfL is {nfl}. "
        f"Plasma GFAP is {gfap}."
    )

    # ── GDS ───────────────────────────────────────────────────────────
    gr = closest_row(gds, ptid, visit_date, viscode=viscode)
    gds_val = gr.get("GDTOTAL") if gr is not None else np.nan
    gds_total = fmtv(gds_val, 0)
    # GDS: score >= 5 indicates depressive symptoms (score < 5 = no depression)
    try:
        gds_float = float(gds_val)
        dep_flag = "Yes" if gds_float >= 5 else "No"
    except (TypeError, ValueError):
        dep_flag = "NaN"
    gds_block = f"Depressive symptoms present (Geriatric Depression Scale)? {dep_flag}. Total GDS score: {gds_total}."

    # ── Hachinski (screening only — carried forward like demographics) ─
    hach_score = fmtv(hach_bl_dict.get(ptid, np.nan), 0)
    hach_block = f"On the Hachinski Ischaemic Score the participant scored: {hach_score}."

    # ── ADAS-Cog13 ────────────────────────────────────────────────────
    ar = closest_row(adas, ptid, visit_date, date_col="VISDATE", viscode=viscode)
    adas_total13 = fmtv(ar.get("TOTAL13") if ar is not None else np.nan)
    adas_block = f"On the Alzheimer's Disease Assessment Scale-Cog13 the participant scored: {adas_total13}."

    # ── CDR ───────────────────────────────────────────────────────────
    cr = closest_row(cdr, ptid, visit_date, date_col="VISDATE", viscode=viscode)
    if cr is not None:
        cdr_global = fmtv(cr.get("CDGLOBAL"))
        cdr_sb     = fmtv(cr.get("CDRSB"))
    else:
        cdr_global = cdr_sb = "NaN"
    cdr_block = f"On the Clinical Dementia Rating scale the participant scored: CDR Global {cdr_global}, CDR Sum of Boxes {cdr_sb}."

    # ── MMSE ──────────────────────────────────────────────────────────
    mr = closest_row(mmse, ptid, visit_date, date_col="VISDATE", viscode=viscode)
    mmse_score = fmtv(mr.get("MMSCORE") if mr is not None else np.nan, 0)
    mmse_block = f"On the Mini Mental State Examination (MMSE) the participant scored: {mmse_score}."

    # ── MoCA ──────────────────────────────────────────────────────────
    mocr = closest_row(moca, ptid, visit_date, date_col="VISDATE", viscode=viscode)
    moca_score = fmtv(calc_moca(mocr), 0)
    moca_block = f"On the Montreal Cognitive Assessment (MoCA) the participant scored: {moca_score}."

    # ── NEUROBAT ──────────────────────────────────────────────────────
    nr = closest_row(neuro, ptid, visit_date, date_col="VISDATE", viscode=viscode)

    def nget(col):
        return nr.get(col, np.nan) if nr is not None else np.nan

    # Clock Drawing
    clock_block = (
        f"On the Clock Drawing Test the participant answered the following questions in this way: "
        f"Approximately circular face {fmts(nget('CLOCKCIRC'))}. "
        f"Symmetry of number placement {fmts(nget('CLOCKSYM'))}. "
        f"Correctness of numbers {fmts(nget('CLOCKNUM'))}. "
        f"Presence of the two hands {fmts(nget('CLOCKHAND'))}. "
        f"Presence of the two hands, set to ten after eleven {fmts(nget('CLOCKTIME'))}. "
        f"Clock Drawing Test Total Score {fmtv(nget('CLOCKSCOR'))}."
    )

    # Category Fluency Animals
    fluency_block = (
        f"On the Category Fluency Test Animals the scores were: "
        f"Total Correct {fmtv(nget('CATANIMSC'))}. "
        f"Perseverations {fmtv(nget('CATANPERS'))}. "
        f"Intrusions {fmtv(nget('CATANINTR'))}."
    )

    # Logical Memory
    lm_block = (
        f"On the Logical Memory Test the participant scored as follows: "
        f"Immediate Recall Total {fmtv(nget('LIMMTOTAL'))}. "
        f"Delayed Recall Total {fmtv(nget('LDELTOTAL'))}. "
        f"Delayed Recall Cued {fmtv(nget('LDELCUE'))}."
    )

    # AVLT
    avlt_block = (
        f"On the Auditory Verbal Learning Test the participant scored as follows in each trial: "
        f"Trial 1 Total {fmtv(nget('AVTOT1'))}. Total Intrusions {fmtv(nget('AVERR1'))}. "
        f"Trial 2 Total {fmtv(nget('AVTOT2'))}. Total Intrusions {fmtv(nget('AVERR2'))}. "
        f"Trial 3 Total {fmtv(nget('AVTOT3'))}. Total Intrusions {fmtv(nget('AVERR3'))}. "
        f"Trial 4 Total {fmtv(nget('AVTOT4'))}. Total Intrusions {fmtv(nget('AVERR4'))}. "
        f"Trial 5 Total {fmtv(nget('AVTOT5'))}. Total Intrusions {fmtv(nget('AVERR5'))}. "
        f"Trial 6 Total {fmtv(nget('AVTOT6'))}. Total Intrusions {fmtv(nget('AVERR6'))}. "
        f"List B Total {fmtv(nget('AVTOTB'))}. Total Intrusions {fmtv(nget('AVERRB'))}. "
        f"30 Minute Delay Total {fmtv(nget('AVDEL30MIN'))}. Total Intrusions {fmtv(nget('AVDELERR1'))}. "
        f"Recognition Score {fmtv(nget('AVDELTOT'))}. Total Intrusions {fmtv(nget('AVDELERR2'))}."
    )

    # Trails A & B
    trails_block = (
        f"Trail Making Test: "
        f"Part A - Time to Complete {fmtv(nget('TRAASCOR'))}. "
        f"Errors of Commission {fmtv(nget('TRAAERRCOM'))}. "
        f"Errors of Omission {fmtv(nget('TRAAERROM'))}. "
        f"Part B - Time to complete {fmtv(nget('TRABSCOR'))}. "
        f"Errors of Commission {fmtv(nget('TRABERRCOM'))}. "
        f"Errors of Omission {fmtv(nget('TRABERROM'))}."
    )

    # ── FAQ ───────────────────────────────────────────────────────────
    fr = closest_row(faq, ptid, visit_date, date_col="VISDATE", viscode=viscode)
    faq_items = ""
    for col, label in FAQ_LABELS:
        val = fmts(fr.get(col, np.nan) if fr is not None else np.nan)
        faq_items += f"{label}. {val}. "
    faq_total = fmtv(fr.get("FAQTOTAL") if fr is not None else np.nan)
    faq_block = (
        f"For the Functional Activities Questionnaire the participant scored as follows for each question: "
        f"{faq_items}Total Score for FAQ is {faq_total}."
    )

    # ── Assemble ──────────────────────────────────────────────────────
    text = "\n\n".join([
        demo_block,
        apoe_block,
        plasma_block,
        gds_block,
        hach_block,
        adas_block,
        cdr_block,
        mmse_block,
        moca_block,
        clock_block,
        fluency_block,
        lm_block,
        avlt_block,
        trails_block,
        faq_block,
    ])
    return text

# ── Prognosis labels (anchored to subject's baseline date) ────────────────
# Labels are computed ONCE at the subject's VISCODE2='bl' date and then
# repeated across all their visit rows.  This makes them static predictors
# rather than recomputing from each visit date.
def compute_prognosis(ptid, years):
    """Return 1 if patient converts to AD within `years` years of their baseline."""
    baseline = bl_date_map.get(ptid)
    if baseline is None: return np.nan
    future_dx = dx_df[
        (dx_df["PTID"] == ptid) &
        (dx_df["EXAMDATE"] > baseline) &
        (dx_df["EXAMDATE"] <= baseline + pd.DateOffset(years=years)) &
        (dx_df["_diag_bin"] == 1)
    ]
    return 1 if not future_dx.empty else 0



# ── Build master rows ──────────────────────────────────────────────────────
# Visit backbone = DXSUM (one row per patient per *clinical* visit).
# The Date column = actual clinical exam date (EXAMDATE from DXSUM).
# This means every visit where cognitive/clinical data was collected gets a
# row, regardless of whether an MRI was performed that day.

dx_df_indexed = dx_df[dx_df["PTID"].isin(valid_ptids)].copy()

# ── Absolute month mapping ────────────────────────────────────────────────────
# ADNI VISCODEs reset between phases (e.g. ADNI1 m24 ≠ ADNI3 Year2 m24).
# Strategy: anchor each subject to their first 'bl' visit date, then compute
# elapsed months from there.  Snap to the nearest standard grid point to
# produce a single unified long-month label per visit.
#

def month_to_label(m):
    if m == 0: return "bl"
    return f"m{m:02d}"   # gives m06, m12, m48, m108, m132, …

def snap_to_grid(elapsed_months, tolerance=4.0):
    """Snap elapsed months to the nearest point in the global MONTH_GRID.
    Returns the grid month if within tolerance, else None."""
    if pd.isna(elapsed_months): return None
    best = min(MONTH_GRID, key=lambda g: abs(elapsed_months - g))
    return best if abs(elapsed_months - best) <= tolerance else None

# Build per-subject absolute baseline date using VISCODE2.
# VISCODE2 == 'bl' is unique per subject across all phases.
# Cross-check: confirm the VISCODE2 'bl' date is the earliest exam date.
bl_date_map = {}
bl_mismatch_count = 0
for ptid, grp in dx_df_indexed.groupby("PTID"):
    grp_sorted = grp.sort_values("EXAMDATE")
    # Primary: VISCODE2 == 'bl'
    v2_bl = grp_sorted[
        grp_sorted["VISCODE2"].astype(str).str.strip().str.lower() == "bl"
    ].dropna(subset=["EXAMDATE"])
    # Earliest exam date across all visits as backup
    earliest_row = grp_sorted.dropna(subset=["EXAMDATE"])
    earliest_date = earliest_row.iloc[0]["EXAMDATE"] if not earliest_row.empty else None
    if not v2_bl.empty:
        v2_bl_date = v2_bl.iloc[0]["EXAMDATE"]
        # Cross-check: earliest date should be within 6 months of VISCODE2 bl
        if earliest_date is not None:
            gap_days = abs((v2_bl_date - earliest_date).days)
            if gap_days > 180:   # >6 months discrepancy
                print(f"  BL date check [{ptid}]: VISCODE2 bl={v2_bl_date.date()} "
                      f"but earliest date={earliest_date.date()} "
                      f"(gap={gap_days}d) — using VISCODE2 bl.")
                bl_mismatch_count += 1
        bl_date_map[ptid] = v2_bl_date
    elif earliest_date is not None:
        # No VISCODE2 bl found — fall back to earliest visit
        print(f"  No VISCODE2 bl for [{ptid}]; using earliest date {earliest_date.date()}.")
        bl_date_map[ptid] = earliest_date

print(f"Baseline dates resolved for {len(bl_date_map)} subjects "
      f"({bl_mismatch_count} with bl/earliest-date mismatch > 6 months).")

# Pre-compute per-subject prognosis now that bl_date_map is available
prognosis_map = {
    ptid: {
        "1y":  compute_prognosis(ptid, 1),
        "2y":  compute_prognosis(ptid, 2),
        "3y":  compute_prognosis(ptid, 3),
        "4y":  compute_prognosis(ptid, 4),
        "5y":  compute_prognosis(ptid, 5),
        "6y":  compute_prognosis(ptid, 6),
        "7y":  compute_prognosis(ptid, 7),
        "8y":  compute_prognosis(ptid, 8),
        "9y":  compute_prognosis(ptid, 9),
        "10y": compute_prognosis(ptid, 10),
    }
    for ptid in valid_ptids
}

# Build MONTH_GRID empirically from VISCODE2 values present in the data.
# This automatically covers whatever timepoints exist, including beyond m132.
import re as _re

def _v2_to_months(v2):
    """Parse a VISCODE2 string to integer months, or None."""
    v = str(v2).strip().lower()
    if v == "bl": return 0
    if v in ("sc", "f"): return None
    m = _re.match(r"^m(\d+)$", v)
    if m: return int(m.group(1))
    return None

_v2_months_in_data = {
    _v2_to_months(v)
    for v in dx_df_indexed["VISCODE2"].dropna().unique()
    if _v2_to_months(v) is not None
}
# Supplement with core early-interval months to ensure full ADNI1/GO coverage
# (some short intervals may not appear in VISCODE2 for this specific cohort).
_core = {0, 3, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78}
MONTH_GRID = sorted(_v2_months_in_data | _core)
print(f"Month grid ({len(MONTH_GRID)} points): {MONTH_GRID}")

def compute_long_viscode(ptid, visit_date, viscode2_raw=None):
    """Compute unified long-month label from elapsed time since the subject's
    VISCODE2=='bl' baseline date. Pure date-based: VISCODE2 is only used to
    locate the baseline (already stored in bl_date_map).

    Returns 'sc' for screening/pre-baseline visits, 'bl' for month 0,
    'mXX' for subsequent months, or a fallback string.
    """
    v2 = str(viscode2_raw).strip().lower() if pd.notna(viscode2_raw) else ""

    # Screening / pre-baseline visits
    if v2 in ("sc", "f"):
        return "sc"

    baseline = bl_date_map.get(ptid)
    if baseline is None or visit_date is None or pd.isna(visit_date):
        return v2 if v2 else "unknown"   # fallback to VISCODE2 label as-is

    elapsed_months = (visit_date - baseline).days / 30.4375

    if elapsed_months < -1:          # pre-baseline visit
        return "sc"

    snapped = snap_to_grid(elapsed_months)
    if snapped is not None:
        return month_to_label(snapped)

    # Beyond the empirical grid: round to nearest month
    return f"m{round(elapsed_months):03d}"


# Pre-load Hachinski baseline score per subject (screening only, carried forward).
# hach has VISDATE for all sc visits; take the first (earliest) record per PTID.
hach_for_bl = hach.sort_values("VISDATE").drop_duplicates("PTID", keep="first")
hach_bl_dict = {
    r["PTID"]: r.get("HMSCORE", np.nan)
    for _, r in hach_for_bl.iterrows()
}
print(f"Hachinski baseline scores loaded for {sum(pd.notna(v) for v in hach_bl_dict.values())} subjects.")

# Baseline diagnosis label per subject (earliest bl/sc visit in DXSUM)
bl_dx = dx_df_indexed[dx_df_indexed["VISCODE2"].isin(["bl", "sc"])]
bl_label_map = {}
for ptid, grp in bl_dx.groupby("PTID"):
    grp = grp.sort_values("EXAMDATE")
    bl_label_map[ptid] = {
        "binary": grp.iloc[0]["_diag_bin"],
        "multi":  grp.iloc[0]["_diag_multi"],
    }

rows_verbose = []
rows_tabular  = []
total = len(dx_df_indexed)
for i, (_, d_row) in enumerate(dx_df_indexed.iterrows()):
    if i % 500 == 0:
        print(f"  Processing row {i}/{total}...", flush=True)
    ptid        = d_row["PTID"]
    viscode2    = d_row.get("VISCODE2", np.nan)
    visit_date  = d_row["EXAMDATE"]
    viscode_long = compute_long_viscode(ptid, visit_date, viscode2)

    # 'f' (failed screen) rows are absent from DXSUM for our cohort; no skip needed.
    # sc rows ARE included as the screening timepoint.

    bl_info = bl_label_map.get(ptid, {})
    prog    = prognosis_map.get(ptid, {})

    data = extract_data(ptid, visit_date, viscode2, d_row)
    text = format_summarised_text(data)

    meta = {
        "Patient_ID":           ptid,
        "VISCODE_long":         viscode_long,
        "VISCODE_2":            viscode2,
        "Date":                 visit_date.date() if pd.notna(visit_date) else np.nan,
        "Label_bl_multi":       bl_info.get("multi",  np.nan),
        "Label_visit_diag":     d_row["_diag_multi"],
        "Label_visit_granular": get_visit_granular(ptid, visit_date),
        "Label_1y":             prog.get("1y",  np.nan),
        "Label_2y":             prog.get("2y",  np.nan),
        "Label_3y":             prog.get("3y",  np.nan),
        "Label_4y":             prog.get("4y",  np.nan),
        "Label_5y":             prog.get("5y",  np.nan),
        "Label_6y":             prog.get("6y",  np.nan),
        "Label_7y":             prog.get("7y",  np.nan),
        "Label_8y":             prog.get("8y",  np.nan),
        "Label_9y":             prog.get("9y",  np.nan),
        "Label_10y":            prog.get("10y", np.nan),
    }
    rows_verbose.append({**meta, "Generated_Text": text})
    rows_tabular.append({**meta, **data})

# ── Build DataFrames ───────────────────────────────────────────────────────
def _sort_master(df):
    def _vn(v):
        v = str(v).strip().lower()
        if v == "sc": return -1
        if v == "bl": return 0
        if v.startswith("m"):
            try: return int(v[1:])
            except ValueError: pass
        return float("inf")
    df = df.copy()
    df["_vsort"] = df["VISCODE_long"].map(_vn)
    return df.sort_values(["Patient_ID", "_vsort"]).drop(columns=["_vsort"]).reset_index(drop=True)

master_verbose = _sort_master(pd.DataFrame(rows_verbose))
master_tabular = _sort_master(pd.DataFrame(rows_tabular))

print(f"Master has {len(master_verbose)} rows, {master_verbose['Patient_ID'].nunique()} unique subjects.")
print("VISCODEs present:", sorted(master_verbose["VISCODE_long"].unique()))

# ── Directory setup ────────────────────────────────────────────────────────
BASE_DIR   = r"D:\ADNI_BIDS_project\derivatives\clinical"
VERBOSE_LI = os.path.join(BASE_DIR, "verbose",  "longitudinal_inputs")
VERBOSE_BV = os.path.join(BASE_DIR, "verbose",  "by_visit")
TABULAR_LI = os.path.join(BASE_DIR, "tabular",  "longitudinal_inputs")
TABULAR_BV = os.path.join(BASE_DIR, "tabular",  "by_visit")
for d in [VERBOSE_LI, VERBOSE_BV, TABULAR_LI, TABULAR_BV]:
    Path(d).mkdir(parents=True, exist_ok=True)

# ── Save master CSVs ───────────────────────────────────────────────────────
master_verbose.to_csv(os.path.join(VERBOSE_LI, "master_clinical.csv"), index=False)
master_tabular.to_csv(os.path.join(TABULAR_LI, "master_clinical.csv"), index=False)
print(f"Saved verbose master → {VERBOSE_LI}")
print(f"Saved tabular master → {TABULAR_LI}")

# ── Save per-visit CSVs ────────────────────────────────────────────────────
all_viscodes = VISIT_COLS + [v for v in sorted(master_verbose["VISCODE_long"].unique())
                             if v not in VISIT_COLS]
for vt in all_viscodes:
    for df, bv_dir in [(master_verbose, VERBOSE_BV), (master_tabular, TABULAR_BV)]:
        sub = df[df["VISCODE_long"] == vt]
        if sub.empty: continue
        sub.to_csv(os.path.join(bv_dir, f"visit_{vt}.csv"), index=False)
    print(f"  Saved visit_{vt}: {len(sub)} rows")

# ── Baseline splits (80/10/10, seeds 0/1/2) ───────────────────────────────
VERBOSE_BL = os.path.join(BASE_DIR, "verbose", "baseline")
TABULAR_BL = os.path.join(BASE_DIR, "tabular", "baseline")
Path(VERBOSE_BL).mkdir(parents=True, exist_ok=True)
Path(TABULAR_BL).mkdir(parents=True, exist_ok=True)

def _make_splits(master, baseline_dir):
    """Save 80/10/10 train/val/test for each of seeds 0, 1, 2."""
    bl_rows = master[master["VISCODE_long"] == "bl"]
    sc_rows = master[master["VISCODE_long"] == "sc"]
    bl_set  = set(bl_rows["Patient_ID"])
    bl_df   = pd.concat([bl_rows,
                         sc_rows[~sc_rows["Patient_ID"].isin(bl_set)]]).drop_duplicates("Patient_ID")
    n       = len(bl_df)
    n_train = int(0.80 * n)
    n_val   = int(0.10 * n)
    for seed in [0, 1, 2]:
        seed_df   = bl_df.sample(frac=1, random_state=seed).reset_index(drop=True)
        seed_dir  = os.path.join(baseline_dir, f"seed_{seed}")
        Path(seed_dir).mkdir(exist_ok=True)
        seed_df.iloc[:n_train].to_csv(os.path.join(seed_dir, "train.csv"), index=False)
        seed_df.iloc[n_train:n_train+n_val].to_csv(os.path.join(seed_dir, "val.csv"), index=False)
        seed_df.iloc[n_train+n_val:].to_csv(os.path.join(seed_dir, "test.csv"), index=False)
        print(f"  seed_{seed} → {seed_dir}  "
              f"Train:{n_train} Val:{n_val} Test:{n-n_train-n_val}")

_make_splits(master_verbose, VERBOSE_BL)
_make_splits(master_tabular, TABULAR_BL)


# ── Save per-subject CSVs ──────────────────────────────────────────────────
# One CSV per patient named {Patient_ID}.csv, rows sorted chronologically.
# Enables leakage-free longitudinal train/val/test splitting at subject level.
VERBOSE_BS = os.path.join(BASE_DIR, "verbose", "by_subject")
TABULAR_BS = os.path.join(BASE_DIR, "tabular", "by_subject")
Path(VERBOSE_BS).mkdir(parents=True, exist_ok=True)
Path(TABULAR_BS).mkdir(parents=True, exist_ok=True)

for ptid in master_verbose["Patient_ID"].unique():
    safe_id = ptid.replace("/", "_")
    master_verbose[master_verbose["Patient_ID"] == ptid].to_csv(
        os.path.join(VERBOSE_BS, f"{safe_id}.csv"), index=False)
    master_tabular[master_tabular["Patient_ID"] == ptid].to_csv(
        os.path.join(TABULAR_BS, f"{safe_id}.csv"), index=False)

n_subj = master_verbose["Patient_ID"].nunique()
print(f"Saved {n_subj} per-subject CSVs → {VERBOSE_BS}")
print(f"Saved {n_subj} per-subject CSVs → {TABULAR_BS}")
print("Done.")
