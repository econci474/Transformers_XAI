"""
01d_build_clinical_csv_stratified_no_cdr.py
============================================
Reads the existing master CSVs produced by 01b, removes CDR columns
(CDR_Global, CDR_SumBoxes) from the tabular data and strips the CDR
sentence from the Generated_Text in the verbose data, then creates
stratified 80/10/10 splits (preserving CN/MCI/AD proportions).

This addresses data leakage from CDR being tightly associated with the
diagnostic label.

Input:
  D:\\ADNI_BIDS_project\\derivatives\\clinical\\verbose\\longitudinal\\master_clinical_verbose.csv
  D:\\ADNI_BIDS_project\\derivatives\\clinical\\tabular\\longitudinal\\master_clinical_tabular.csv

Output (under D:\\ADNI_BIDS_project\\derivatives\\clinical\\no_cdr_stratified\\):
  verbose\\longitudinal\\master_clinical_verbose.csv
  tabular\\longitudinal\\master_clinical_tabular.csv
  verbose\\baseline\\seed_{0,1,2}\\{train,val,test}.csv
  tabular\\baseline\\seed_{0,1,2}\\{train,val,test}.csv

Run:
  python clinical_pipeline/01d_build_clinical_csv_stratified_no_cdr.py
"""

import os
import re
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.model_selection import train_test_split

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = r"D:\ADNI_BIDS_project\derivatives\clinical"
VERBOSE_MASTER = os.path.join(BASE_DIR, "verbose", "longitudinal", "master_clinical_verbose.csv")
TABULAR_MASTER = os.path.join(BASE_DIR, "tabular", "longitudinal", "master_clinical_tabular.csv")

OUT_ROOT    = os.path.join(BASE_DIR, "no_cdr_stratified")
VERBOSE_OUT = os.path.join(OUT_ROOT, "verbose")
TABULAR_OUT = os.path.join(OUT_ROOT, "tabular")

# CDR columns to drop from tabular data
CDR_COLUMNS = ["CDR_Global", "CDR_SumBoxes"]

# Regex to strip the CDR sentence(s) from Generated_Text.
# Two formats exist:
#   Tabular: "...ADAS-Cog13 total: X; Clinical Dementia Rating (CDR) scale global: Y, CDR sum of boxes: Z."
#   Verbose: "On the Clinical Dementia Rating scale the participant scored: CDR Global Y, CDR Sum of Boxes Z."
CDR_TEXT_PATTERNS = [
    # Tabular inline: "; Clinical Dementia Rating..." up to the period
    r";\s*Clinical Dementia Rating \(CDR\) scale global:[^.]+\.",
    # Verbose standalone block (paragraph)
    r"On the Clinical Dementia Rating scale the participant scored:[^\n]+",
    # Catch-all: any sentence containing "CDR" as a standalone sentence
    r"[^.\n]*\bCDR\b[^.\n]*\.\s*",
]

SEEDS = [0, 1, 2]


def strip_cdr_from_text(text: str) -> str:
    """Remove CDR-related sentences from Generated_Text."""
    if pd.isna(text) or not isinstance(text, str):
        return text
    for pattern in CDR_TEXT_PATTERNS:
        text = re.sub(pattern, "", text, flags=re.IGNORECASE)
    # Clean up double/triple newlines left behind
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip()


# ── Load masters ──────────────────────────────────────────────────────────────
print("=" * 60)
print("  BUILD NO-CDR STRATIFIED CSVs")
print("=" * 60)

print(f"\nLoading verbose master : {VERBOSE_MASTER}")
master_verbose = pd.read_csv(VERBOSE_MASTER, low_memory=False)
print(f"  rows: {len(master_verbose)}, subjects: {master_verbose['Patient_ID'].nunique()}")

print(f"Loading tabular master : {TABULAR_MASTER}")
master_tabular = pd.read_csv(TABULAR_MASTER, low_memory=False)
print(f"  rows: {len(master_tabular)}, subjects: {master_tabular['Patient_ID'].nunique()}")

# ── Remove CDR from tabular ───────────────────────────────────────────────────
dropped_cols = [c for c in CDR_COLUMNS if c in master_tabular.columns]
master_tabular = master_tabular.drop(columns=dropped_cols, errors="ignore")
print(f"\nDropped tabular columns: {dropped_cols}")
print(f"  remaining columns: {len(master_tabular.columns)}")

# Also drop CDR from verbose (it has all tabular cols + Generated_Text)
dropped_verbose = [c for c in CDR_COLUMNS if c in master_verbose.columns]
master_verbose = master_verbose.drop(columns=dropped_verbose, errors="ignore")
print(f"Dropped verbose columns: {dropped_verbose}")

# ── Strip CDR from Generated_Text in verbose ──────────────────────────────────
if "Generated_Text" in master_verbose.columns:
    print("Stripping CDR sentences from Generated_Text...")
    before_sample = master_verbose["Generated_Text"].iloc[0][:200] if len(master_verbose) > 0 else ""
    master_verbose["Generated_Text"] = master_verbose["Generated_Text"].apply(strip_cdr_from_text)
    after_sample = master_verbose["Generated_Text"].iloc[0][:200] if len(master_verbose) > 0 else ""
    print(f"  Before: ...{before_sample}...")
    print(f"  After : ...{after_sample}...")

# ── Save master CSVs ──────────────────────────────────────────────────────────
verbose_long_dir = os.path.join(VERBOSE_OUT, "longitudinal")
tabular_long_dir = os.path.join(TABULAR_OUT, "longitudinal")
Path(verbose_long_dir).mkdir(parents=True, exist_ok=True)
Path(tabular_long_dir).mkdir(parents=True, exist_ok=True)

master_verbose.to_csv(os.path.join(verbose_long_dir, "master_clinical_verbose.csv"), index=False)
master_tabular.to_csv(os.path.join(tabular_long_dir, "master_clinical_tabular.csv"), index=False)
print(f"\nSaved verbose master → {verbose_long_dir}")
print(f"Saved tabular master → {tabular_long_dir}")

# ── Save per-visit CSVs ───────────────────────────────────────────────────────
verbose_bv = os.path.join(VERBOSE_OUT, "by_visit")
tabular_bv = os.path.join(TABULAR_OUT, "by_visit")
Path(verbose_bv).mkdir(parents=True, exist_ok=True)
Path(tabular_bv).mkdir(parents=True, exist_ok=True)

for vt in master_verbose["VISCODE_long"].unique():
    for df, bv_dir in [(master_verbose, verbose_bv), (master_tabular, tabular_bv)]:
        sub = df[df["VISCODE_long"] == vt]
        if sub.empty:
            continue
        sub.to_csv(os.path.join(bv_dir, f"visit_{vt}.csv"), index=False)
print(f"Saved per-visit CSVs")

# ── Stratified baseline splits ────────────────────────────────────────────────
def _make_stratified_splits(master, baseline_dir, label):
    """Save 80/10/10 train/val/test with stratification on Label_bl_multi."""
    bl_rows = master[master["VISCODE_long"] == "bl"]
    sc_rows = master[master["VISCODE_long"] == "sc"]
    bl_set  = set(bl_rows["Patient_ID"])
    bl_df   = pd.concat([bl_rows,
                         sc_rows[~sc_rows["Patient_ID"].isin(bl_set)]]).drop_duplicates("Patient_ID")

    n = len(bl_df)
    strat_col = bl_df["Label_bl_multi"].fillna("UNKNOWN")

    print(f"\n  {label}: {n} baseline subjects")
    print(f"  Class distribution: {dict(strat_col.value_counts())}")

    Path(baseline_dir).mkdir(parents=True, exist_ok=True)

    for seed in SEEDS:
        train_df, rest_df = train_test_split(
            bl_df, test_size=0.20, random_state=seed,
            stratify=strat_col)

        rest_strat = rest_df["Label_bl_multi"].fillna("UNKNOWN")
        val_df, test_df = train_test_split(
            rest_df, test_size=0.50, random_state=seed,
            stratify=rest_strat)

        seed_dir = os.path.join(baseline_dir, f"seed_{seed}")
        Path(seed_dir).mkdir(exist_ok=True)
        train_df.to_csv(os.path.join(seed_dir, "train.csv"), index=False)
        val_df.to_csv(os.path.join(seed_dir, "val.csv"),     index=False)
        test_df.to_csv(os.path.join(seed_dir, "test.csv"),   index=False)

        for name, df in [("train", train_df), ("val", val_df), ("test", test_df)]:
            counts = df["Label_bl_multi"].value_counts().to_dict()
            total  = len(df)
            pcts   = {k: f"{v/total*100:.1f}%" for k, v in counts.items()}
            print(f"    seed_{seed} {name:5s}: N={total:3d}  {counts}  ({pcts})")


verbose_bl = os.path.join(VERBOSE_OUT, "baseline")
tabular_bl = os.path.join(TABULAR_OUT, "baseline")

print("\n── Verbose splits ──")
_make_stratified_splits(master_verbose, verbose_bl, "Verbose")

print("\n── Tabular splits ──")
_make_stratified_splits(master_tabular, tabular_bl, "Tabular")

# ── Save per-subject CSVs ─────────────────────────────────────────────────────
verbose_bs = os.path.join(VERBOSE_OUT, "by_subject")
tabular_bs = os.path.join(TABULAR_OUT, "by_subject")
Path(verbose_bs).mkdir(parents=True, exist_ok=True)
Path(tabular_bs).mkdir(parents=True, exist_ok=True)

for ptid in master_verbose["Patient_ID"].unique():
    safe_id = ptid.replace("/", "_")
    master_verbose[master_verbose["Patient_ID"] == ptid].to_csv(
        os.path.join(verbose_bs, f"{safe_id}.csv"), index=False)
    master_tabular[master_tabular["Patient_ID"] == ptid].to_csv(
        os.path.join(tabular_bs, f"{safe_id}.csv"), index=False)

n_subj = master_verbose["Patient_ID"].nunique()
print(f"\nSaved {n_subj} per-subject CSVs")

# ── Verify CDR is gone ────────────────────────────────────────────────────────
print("\n── Verification ──")
test_tab = pd.read_csv(os.path.join(tabular_bl, "seed_0", "test.csv"), low_memory=False)
cdr_remaining = [c for c in test_tab.columns if "CDR" in c.upper()]
print(f"  CDR columns in tabular test: {cdr_remaining if cdr_remaining else 'NONE ✓'}")

test_verb = pd.read_csv(os.path.join(verbose_bl, "seed_0", "test.csv"), low_memory=False)
if "Generated_Text" in test_verb.columns:
    cdr_in_text = test_verb["Generated_Text"].str.contains("CDR", case=False, na=False).sum()
    print(f"  Rows with 'CDR' in Generated_Text: {cdr_in_text} (should be 0)")

print("\nDone.")
