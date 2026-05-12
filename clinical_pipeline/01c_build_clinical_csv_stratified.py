"""
01c_build_clinical_csv_stratified.py
====================================
Re-creates the baseline 80/10/10 splits using *stratified* sampling so that
the CN / MCI / AD proportions from Label_bl_multi are preserved across
train, val, and test sets.

Uses sklearn.model_selection.train_test_split with stratify= for both the
initial 80/20 split (train vs rest) and the subsequent 50/50 split of the
remainder into val and test (each 10% of the full dataset).

Reads the existing master CSVs produced by 01b (verbose and tabular) to
avoid re-building everything from scratch.

Input:
  D:\\ADNI_BIDS_project\\derivatives\\clinical\\verbose\\longitudinal_inputs\\master_clinical.csv
  D:\\ADNI_BIDS_project\\derivatives\\clinical\\tabular\\longitudinal_inputs\\master_clinical.csv

Output:
  D:\\ADNI_BIDS_project\\derivatives\\clinical\\verbose\\baseline_stratified\\seed_{0,1,2}\\{train,val,test}.csv
  D:\\ADNI_BIDS_project\\derivatives\\clinical\\tabular\\baseline_stratified\\seed_{0,1,2}\\{train,val,test}.csv

Run:
  python clinical_pipeline/01c_build_clinical_csv_stratified.py
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.model_selection import train_test_split

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR   = r"D:\ADNI_BIDS_project\derivatives\clinical"
VERBOSE_MASTER = os.path.join(BASE_DIR, "verbose", "longitudinal", "master_clinical_verbose.csv")
TABULAR_MASTER = os.path.join(BASE_DIR, "tabular", "longitudinal", "master_clinical_tabular.csv")

VERBOSE_OUT = os.path.join(BASE_DIR, "verbose", "baseline_stratified")
TABULAR_OUT = os.path.join(BASE_DIR, "tabular", "baseline_stratified")

SEEDS = [0, 1, 2]

# ── Load masters ──────────────────────────────────────────────────────────────
print("Loading master CSVs...")
master_verbose = pd.read_csv(VERBOSE_MASTER, low_memory=False)
master_tabular = pd.read_csv(TABULAR_MASTER, low_memory=False)
print(f"  verbose: {len(master_verbose)} rows, {master_verbose['Patient_ID'].nunique()} subjects")
print(f"  tabular: {len(master_tabular)} rows, {master_tabular['Patient_ID'].nunique()} subjects")


# ── Stratified splitting ──────────────────────────────────────────────────────
def _make_stratified_splits(master, baseline_dir):
    """Save 80/10/10 train/val/test with stratification on Label_bl_multi."""
    # Extract baseline rows (bl, falling back to sc)
    bl_rows = master[master["VISCODE_long"] == "bl"]
    sc_rows = master[master["VISCODE_long"] == "sc"]
    bl_set  = set(bl_rows["Patient_ID"])
    bl_df   = pd.concat([bl_rows,
                         sc_rows[~sc_rows["Patient_ID"].isin(bl_set)]]).drop_duplicates("Patient_ID")

    n = len(bl_df)
    # Stratification column: Label_bl_multi (CN / MCI / AD)
    strat_col = bl_df["Label_bl_multi"].copy()
    # Fill any NaN in stratification column so train_test_split doesn't fail
    strat_col = strat_col.fillna("UNKNOWN")

    print(f"\n  Total baseline subjects: {n}")
    print(f"  Class distribution: {dict(strat_col.value_counts())}")

    Path(baseline_dir).mkdir(parents=True, exist_ok=True)

    for seed in SEEDS:
        # Step 1: split 80% train / 20% remainder
        train_df, rest_df = train_test_split(
            bl_df, test_size=0.20, random_state=seed,
            stratify=strat_col)

        # Step 2: split remainder 50/50 → 10% val + 10% test
        rest_strat = rest_df["Label_bl_multi"].fillna("UNKNOWN")
        val_df, test_df = train_test_split(
            rest_df, test_size=0.50, random_state=seed,
            stratify=rest_strat)

        seed_dir = os.path.join(baseline_dir, f"seed_{seed}")
        Path(seed_dir).mkdir(exist_ok=True)
        train_df.to_csv(os.path.join(seed_dir, "train.csv"), index=False)
        val_df.to_csv(os.path.join(seed_dir, "val.csv"),     index=False)
        test_df.to_csv(os.path.join(seed_dir, "test.csv"),   index=False)

        # Report proportions
        for name, df in [("train", train_df), ("val", val_df), ("test", test_df)]:
            counts = df["Label_bl_multi"].value_counts().to_dict()
            total  = len(df)
            pcts   = {k: f"{v/total*100:.1f}%" for k, v in counts.items()}
            print(f"  seed_{seed} {name:5s}: N={total:3d}  {counts}  ({pcts})")
        print()


print("=" * 60)
print("  STRATIFIED BASELINE SPLITS")
print("=" * 60)

print("\n── Verbose ──")
_make_stratified_splits(master_verbose, VERBOSE_OUT)

print("── Tabular ──")
_make_stratified_splits(master_tabular, TABULAR_OUT)

# ── Quick comparison with original (non-stratified) splits ────────────────────
print("=" * 60)
print("  COMPARISON: Original vs Stratified (seed 0, tabular)")
print("=" * 60)

orig_dir = os.path.join(BASE_DIR, "tabular", "baseline", "seed_0")
strat_dir = os.path.join(BASE_DIR, "tabular", "baseline_stratified", "seed_0")

for name in ["train", "val", "test"]:
    orig  = pd.read_csv(os.path.join(orig_dir, f"{name}.csv"), low_memory=False)
    strat = pd.read_csv(os.path.join(strat_dir, f"{name}.csv"), low_memory=False)

    orig_counts  = orig["Label_bl_multi"].value_counts().sort_index()
    strat_counts = strat["Label_bl_multi"].value_counts().sort_index()

    print(f"\n  {name.upper()}")
    print(f"    Original  (N={len(orig):3d}):  {dict(orig_counts)}")
    print(f"    Stratified(N={len(strat):3d}):  {dict(strat_counts)}")

    # Show proportions
    for label in ["CN", "MCI", "AD"]:
        o_pct = orig_counts.get(label, 0) / len(orig) * 100
        s_pct = strat_counts.get(label, 0) / len(strat) * 100
        print(f"    {label}: {o_pct:.1f}% → {s_pct:.1f}%")

print("\nDone.")
