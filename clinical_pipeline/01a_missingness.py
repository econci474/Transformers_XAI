"""
missingness.py
Generate missingness report for all clinical features across the 616-subject cohort.

Outputs to: D:/ADNI_BIDS_project/derivatives/clinical/missingness/
  - One CSV per feature: rows=patients, columns=visit months, values=exam date (or NaN)
  - master_missingness.csv: features x (% missing, total entries, patients with entries)
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path

SUBJECTS_TSV = r"D:\ADNI_BIDS_project\bids\genotype\subjects_with_snp_and_mri.tsv"
CLINICAL_DIR = r"D:\ADNI_BIDS_project\sourcedata\clinical"
MASTER_CSV   = r"D:\ADNI_BIDS_project\derivatives\clinical\master_clinical.csv"
OUT_DIR      = r"D:\ADNI_BIDS_project\derivatives\clinical\missingness"
Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

subj   = pd.read_csv(SUBJECTS_TSV, sep="\t")
valid  = set(subj["adni_subject_id"].str.strip())
master = pd.read_csv(MASTER_CSV)
all_ptids  = sorted(master["Patient_ID"].unique())

def visit_sort_key(v):
    """Sort visit labels chronologically: bl=0, sc=-1, m06=6, m12=12, ..."""
    v = str(v).strip().lower()
    if v == "bl": return 0
    if v == "sc": return -1
    if v.startswith("m"):
        try: return int(v[1:])
        except ValueError: pass
    return float("inf")

all_visits = sorted(master["VISCODE_long"].dropna().unique(), key=visit_sort_key)

def load(fname):
    df = pd.read_csv(os.path.join(CLINICAL_DIR, fname), low_memory=False)
    if "PTID" in df.columns:
        df["PTID"] = df["PTID"].astype(str).str.strip()
    return df[df["PTID"].isin(valid)] if "PTID" in df.columns else df

# ── Map PTID + exam date -> VISCODE_long using the master CSV ──────────────
# We create a date-to-label lookup per subject from the master.
master["Date"] = pd.to_datetime(master["Date"], errors="coerce")
date_label_map = {}   # (ptid, date_str) -> viscode_long
for _, r in master.iterrows():
    key = (r["Patient_ID"], str(r["Date"]))
    date_label_map[key] = r["VISCODE_long"]

def snap_visit(ptid, exam_date, tol_days=60):
    """Find nearest VISCODE_long for a given ptid + exam_date (within tol_days)."""
    if pd.isna(exam_date): return np.nan
    sub = master[master["Patient_ID"] == ptid].dropna(subset=["Date"])
    if sub.empty: return np.nan
    diffs = (sub["Date"] - exam_date).abs()
    best_idx = diffs.idxmin()
    if diffs[best_idx].days <= tol_days:
        return sub.loc[best_idx, "VISCODE_long"]
    return np.nan

# ── Feature definitions ────────────────────────────────────────────────────
# Each entry: (filename, date_col, val_col)
# val_col: if set, a row is only counted as 'present' when this specific column is non-null.
#           If None, any row for that visit is counted (presence of the test).
FEATURES = {
    # Whole-test presence (any row counts)
    "CDR":                     ("CDR.csv",      "VISDATE",   None),
    "MMSE":                    ("MMSE.csv",     "VISDATE",   None),
    "MOCA":                    ("MOCA.csv",     "VISDATE",   None),
    "ADAS":                    ("ADAS.csv",     "VISDATE",   None),
    "FAQ":                     ("FAQ.csv",      "VISDATE",   None),
    "GDS":                     ("GDSCALE.csv",  "VISDATE",   None),
    "Hachinski":               ("MODHACH.csv",  "VISDATE",   None),
    "Plasma_UPENN":            ("UPENN_PLASMA_FUJIREBIO_QUANTERIX.csv", "EXAMDATE", None),
    "pTau181_UGot":            ("UGOTPTAU181.csv", "EXAMDATE", None),
    # NEUROBAT sub-tests (presence = that specific total score is non-null)
    "CDT_total":               ("NEUROBAT.csv", "VISDATE",   "CLOCKSCOR"),
    "CategoryFluency_Animals": ("NEUROBAT.csv", "VISDATE",   "CATANIMSC"),
    "LogicalMemory_I":         ("NEUROBAT.csv", "VISDATE",   "LIMMTOTAL"),
    "LogicalMemory_II":        ("NEUROBAT.csv", "VISDATE",   "LDELTOTAL"),
    "AVLT":                    ("NEUROBAT.csv", "VISDATE",   "AVTOT1"),
    "TrailsA":                 ("NEUROBAT.csv", "VISDATE",   "TRAASCOR"),
    "TrailsB":                 ("NEUROBAT.csv", "VISDATE",   "TRABSCOR"),
}

summary_rows = []

for feat, (fname, date_col, val_col) in FEATURES.items():
    print(f"Processing {feat}...")
    df = load(fname)
    if date_col not in df.columns:
        print(f"  WARNING: {date_col} not in {fname}, skipping.")
        continue
    df[date_col] = pd.to_datetime(df[date_col], errors="coerce")
    df = df.dropna(subset=[date_col])

    # If a specific score column is required, only keep rows where it is non-null
    if val_col is not None:
        if val_col not in df.columns:
            print(f"  WARNING: val_col '{val_col}' not in {fname}, skipping.")
            continue
        df = df[df[val_col].notna()]

    # pTau181 uses RID not PTID -- map via master
    if feat == "pTau181_UGot":
        rid_map = {}
        for ptid in valid:
            rid = ptid.split("_S_")[-1].lstrip("0") or "0"
            rid_map[rid] = ptid
            rid_map[ptid.split("_S_")[-1]] = ptid
        df["PTID"] = df["RID"].astype(str).str.strip().map(
            lambda r: rid_map.get(r.lstrip("0"), rid_map.get(r, np.nan)))
        df = df.dropna(subset=["PTID"])

    # Build matrix: rows=patients, columns=visit months, values=exam date
    matrix = pd.DataFrame(index=all_ptids, columns=all_visits, dtype=object)
    matrix[:] = np.nan

    for ptid, grp in df.groupby("PTID"):
        if ptid not in valid: continue
        for _, row in grp.iterrows():
            exam_date = row[date_col]
            vlabel    = snap_visit(ptid, exam_date)
            if vlabel is not None and not pd.isna(vlabel):
                matrix.at[ptid, vlabel] = str(exam_date.date())

    # Drop visit columns with no data at all, then sort chronologically
    matrix = matrix.dropna(axis=1, how="all")
    ordered_cols = sorted(matrix.columns, key=visit_sort_key)
    matrix = matrix[ordered_cols]
    out_path = os.path.join(OUT_DIR, f"{feat}.csv")
    matrix.index.name = "Patient_ID"
    matrix.to_csv(out_path)
    print(f"  Saved {feat}.csv  ({matrix.shape})")

    # Summary stats
    total_cells   = matrix.size
    filled_cells  = matrix.notna().sum().sum()
    pct_present   = 100 * filled_cells / total_cells if total_cells else 0
    pct_missing   = 100 - pct_present
    pts_with_data = int(matrix.notna().any(axis=1).sum())
    summary_rows.append({
        "Feature":              feat,
        "Pct_Missing":          round(pct_missing, 1),
        "Pct_Present":          round(pct_present, 1),
        "Total_Entries":        int(filled_cells),
        "Total_Cells":          int(total_cells),
        "Patients_With_Entries": pts_with_data,
        "Total_Patients":       len(all_ptids),
    })

summary = pd.DataFrame(summary_rows).sort_values("Pct_Missing", ascending=False)
summary_path = os.path.join(OUT_DIR, "master_missingness.csv")
summary.to_csv(summary_path, index=False)
print(f"\nMaster missingness summary saved -> {summary_path}")
print(summary.to_string(index=False))
