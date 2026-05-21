"""
01d_build_clinical_csv_stratified_no_cdr.py
============================================
Reads the existing master CSVs produced by 01b, removes CDR columns
(CDR_Global, CDR_SumBoxes) from the tabular data and strips the CDR
sentence from the Generated_Text in the verbose data, then creates
stratified 80/10/10 splits (preserving CN/MCI/AD proportions).

This addresses data leakage from CDR being tightly associated with the
diagnostic label.

The script also writes a sister `_post_exclusion` derivative tree with
diagnostic-reversion + corrupted-MRI subjects filtered out (see
`bidsification.exclusions.is_excluded_subject(..., include_diagnostic_reversion=True)`).
The pre-exclusion tree is preserved unchanged.

Input:
  D:\\ADNI_BIDS_project\\derivatives\\clinical\\verbose\\longitudinal\\master_clinical_verbose.csv
  D:\\ADNI_BIDS_project\\derivatives\\clinical\\tabular\\longitudinal\\master_clinical_tabular.csv

Output trees:
  D:\\ADNI_BIDS_project\\derivatives\\clinical\\no_cdr_stratified\\
                ……_post_exclusion\\
    verbose\\longitudinal\\master_clinical_verbose.csv
    tabular\\longitudinal\\master_clinical_tabular.csv
    verbose\\baseline\\seed_{0,1,2}\\{train,val,test}.csv
    tabular\\baseline\\seed_{0,1,2}\\{train,val,test}.csv
    verbose\\by_visit\\visit_*.csv,  by_subject\\<PID>.csv
    tabular\\by_visit\\visit_*.csv,  by_subject\\<PID>.csv

Run:
  python clinical_pipeline/01d_build_clinical_csv_stratified_no_cdr.py
"""

import os
import re
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.model_selection import train_test_split

# Resolve repo root so we can import bidsification.exclusions
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
from bidsification.exclusions import is_excluded_subject  # noqa: E402

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = r"D:\ADNI_BIDS_project\derivatives\clinical"
VERBOSE_MASTER = os.path.join(BASE_DIR, "verbose", "longitudinal", "master_clinical_verbose.csv")
TABULAR_MASTER = os.path.join(BASE_DIR, "tabular", "longitudinal", "master_clinical_tabular.csv")

# Output roots: pre-exclusion (existing) + post-exclusion (sister)
OUT_ROOT_PRE  = os.path.join(BASE_DIR, "no_cdr_stratified")
OUT_ROOT_POST = os.path.join(BASE_DIR, "no_cdr_stratified_post_exclusion")

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
        for name, df in [("train", train_df), ("val", val_df), ("test", test_df)]:
            try:
                df.to_csv(os.path.join(seed_dir, f"{name}.csv"), index=False)
            except PermissionError as exc:
                print(f"    [WARN] could not write seed_{seed}/{name}.csv: {exc} — skipping.")

        for name, df in [("train", train_df), ("val", val_df), ("test", test_df)]:
            counts = df["Label_bl_multi"].value_counts().to_dict()
            total  = len(df)
            pcts   = {k: f"{v/total*100:.1f}%" for k, v in counts.items()}
            print(f"    seed_{seed} {name:5s}: N={total:3d}  {counts}  ({pcts})")


def _safe_csv(df, path, label=""):
    """to_csv that tolerates locked files (Excel-open etc) — warns and skips."""
    try:
        df.to_csv(path, index=False)
    except PermissionError as exc:
        print(f"    [WARN] could not write {path} ({label}): {exc} — skipping.")


def _emit_tree(out_root, master_verbose, master_tabular, tag):
    """
    Write the full output subtree under out_root, mirroring the existing layout:
        verbose/longitudinal/, baseline/seed_*/, by_visit/, by_subject/
        tabular/   (same)
    `tag` is used purely for log lines (e.g. "PRE", "POST"). Individual file
    writes that fail with PermissionError (e.g. file open in Excel) are
    logged and skipped — the rest of the tree still gets written.
    """
    print(f"\n{'='*60}")
    print(f"  EMIT TREE [{tag}] → {out_root}")
    print(f"{'='*60}")
    print(f"  rows verbose: {len(master_verbose)}, "
          f"subjects: {master_verbose['Patient_ID'].nunique()}")
    print(f"  rows tabular: {len(master_tabular)}, "
          f"subjects: {master_tabular['Patient_ID'].nunique()}")

    verbose_root = os.path.join(out_root, "verbose")
    tabular_root = os.path.join(out_root, "tabular")

    # ── Longitudinal masters ─────────────────────────────────────────────
    verbose_long_dir = os.path.join(verbose_root, "longitudinal")
    tabular_long_dir = os.path.join(tabular_root, "longitudinal")
    Path(verbose_long_dir).mkdir(parents=True, exist_ok=True)
    Path(tabular_long_dir).mkdir(parents=True, exist_ok=True)

    _safe_csv(master_verbose,
              os.path.join(verbose_long_dir, "master_clinical_verbose.csv"),
              "verbose longitudinal master")
    _safe_csv(master_tabular,
              os.path.join(tabular_long_dir, "master_clinical_tabular.csv"),
              "tabular longitudinal master")
    print(f"  longitudinal masters → {verbose_long_dir}, {tabular_long_dir}")

    # ── By-visit CSVs ─────────────────────────────────────────────────────
    verbose_bv = os.path.join(verbose_root, "by_visit")
    tabular_bv = os.path.join(tabular_root, "by_visit")
    Path(verbose_bv).mkdir(parents=True, exist_ok=True)
    Path(tabular_bv).mkdir(parents=True, exist_ok=True)

    for vt in master_verbose["VISCODE_long"].dropna().unique():
        for df, bv_dir, label in [(master_verbose, verbose_bv, "verbose"),
                                   (master_tabular, tabular_bv, "tabular")]:
            sub = df[df["VISCODE_long"] == vt]
            if sub.empty:
                continue
            _safe_csv(sub, os.path.join(bv_dir, f"visit_{vt}.csv"),
                      f"{label} by_visit visit_{vt}")
    print(f"  by_visit CSVs → {verbose_bv}, {tabular_bv}")

    # ── Stratified baseline splits ────────────────────────────────────────
    verbose_bl = os.path.join(verbose_root, "baseline")
    tabular_bl = os.path.join(tabular_root, "baseline")

    print(f"\n  ── Verbose splits ──")
    _make_stratified_splits(master_verbose, verbose_bl, f"Verbose [{tag}]")

    print(f"\n  ── Tabular splits ──")
    _make_stratified_splits(master_tabular, tabular_bl, f"Tabular [{tag}]")

    # ── By-subject CSVs ───────────────────────────────────────────────────
    verbose_bs = os.path.join(verbose_root, "by_subject")
    tabular_bs = os.path.join(tabular_root, "by_subject")
    Path(verbose_bs).mkdir(parents=True, exist_ok=True)
    Path(tabular_bs).mkdir(parents=True, exist_ok=True)

    for ptid in master_verbose["Patient_ID"].unique():
        safe_id = str(ptid).replace("/", "_")
        _safe_csv(master_verbose[master_verbose["Patient_ID"] == ptid],
                  os.path.join(verbose_bs, f"{safe_id}.csv"),
                  f"verbose by_subject {safe_id}")
    for ptid in master_tabular["Patient_ID"].unique():
        safe_id = str(ptid).replace("/", "_")
        _safe_csv(master_tabular[master_tabular["Patient_ID"] == ptid],
                  os.path.join(tabular_bs, f"{safe_id}.csv"),
                  f"tabular by_subject {safe_id}")
    n_subj_v = master_verbose["Patient_ID"].nunique()
    n_subj_t = master_tabular["Patient_ID"].nunique()
    print(f"\n  by_subject CSVs: verbose={n_subj_v}, tabular={n_subj_t}")

    # ── Verification ──────────────────────────────────────────────────────
    print(f"\n  ── Verification [{tag}] ──")
    try:
        test_tab = pd.read_csv(os.path.join(tabular_bl, "seed_0", "test.csv"), low_memory=False)
        cdr_remaining = [c for c in test_tab.columns if "CDR" in c.upper()]
        print(f"    CDR columns in tabular test: {cdr_remaining if cdr_remaining else 'NONE ✓'}")

        test_verb = pd.read_csv(os.path.join(verbose_bl, "seed_0", "test.csv"), low_memory=False)
        if "Generated_Text" in test_verb.columns:
            cdr_in_text = test_verb["Generated_Text"].str.contains("CDR", case=False, na=False).sum()
            print(f"    Rows with 'CDR' in Generated_Text: {cdr_in_text} (should be 0)")
    except FileNotFoundError as exc:
        print(f"    [WARN] verification skipped: {exc}")


# ── Load masters ──────────────────────────────────────────────────────────────
print("=" * 60)
print("  BUILD NO-CDR STRATIFIED CSVs  (pre + post-exclusion)")
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

# ── Pre-exclusion tree ────────────────────────────────────────────────────────
_emit_tree(OUT_ROOT_PRE, master_verbose, master_tabular, tag="PRE")

# ── Post-exclusion tree ───────────────────────────────────────────────────────
# Filter out diagnostic-reversion (Excluded==1) + corrupted-MRI + site-381.
master_verbose["Patient_ID"] = master_verbose["Patient_ID"].astype(str)
master_tabular["Patient_ID"] = master_tabular["Patient_ID"].astype(str)

excluded_mask_verb = master_verbose["Patient_ID"].apply(
    lambda p: is_excluded_subject(p, include_diagnostic_reversion=True))
excluded_mask_tab  = master_tabular["Patient_ID"].apply(
    lambda p: is_excluded_subject(p, include_diagnostic_reversion=True))

n_excl_subj_v = master_verbose.loc[excluded_mask_verb, "Patient_ID"].nunique()
n_excl_subj_t = master_tabular.loc[excluded_mask_tab,  "Patient_ID"].nunique()
print(f"\n[exclusion] dropping {n_excl_subj_v} verbose subjects "
      f"({int(excluded_mask_verb.sum())} rows)")
print(f"[exclusion] dropping {n_excl_subj_t} tabular subjects "
      f"({int(excluded_mask_tab.sum())} rows)")

master_verbose_post = master_verbose[~excluded_mask_verb].reset_index(drop=True)
master_tabular_post = master_tabular[~excluded_mask_tab ].reset_index(drop=True)

_emit_tree(OUT_ROOT_POST, master_verbose_post, master_tabular_post, tag="POST")

print("\nDone.")
