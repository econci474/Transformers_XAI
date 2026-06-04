r"""
01i_build_longitudinal_visit_splits.py
======================================
Build per-seed ALL-VISITS training splits for the single longitudinal T2 encoder
(Task 1 of the T4 longitudinal-MAMBA arm). One BioClinical-ModernBERT-large model is
fine-tuned on EVERY clinical visit (one row per patient-visit) with the CONTEMPORANEOUS
label (`Label_visit_diag`, CN/MCI/AD), giving a single embedding space from which each
visit's embedding can be extracted (see 03j_finetune_longitudinal_extract.py).

This script does NOT re-split. It assigns every visit row of the longitudinal master to
its patient's CANONICAL fold (the leakage-safe, modality-shared 80/10/10 partition), so a
patient's visits never straddle train/test and the held-out test patients stay consistent
with the MRI/SNP arms and 02l.

Inputs
------
  Longitudinal master (one row per patient-visit, has Generated_Text + Label_visit_diag):
    D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion\
        verbose\longitudinal\master_clinical_verbose.csv
  Canonical per-seed patient folds (Patient_ID -> {train,val,test}), shared across modalities:
    D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion\
        tabular\baseline\seed_{0,1,2}\{train,val,test}.csv

Outputs
-------
  D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion\
      verbose\longitudinal_allvisits\seed_{0,1,2}\{train,val,test}.csv
  Each output carries ALL of the master's columns plus `visit_months` (bl=0, m06=6, ...).
  ALL visit rows are kept (incl. the handful with NaN Label_visit_diag): training's
  load_split drops NaN labels automatically, and keeping them lets the extractor still
  emit an embedding for any in-window visit ("extract for ANY in-window visit").

Notes
-----
  * Fold assignment is purely by Patient_ID -> canonical fold. A longitudinal patient absent
    from the canonical split (e.g. excluded upstream) is dropped and reported — this applies
    the same cohort/exclusion decisions the canonical split already encodes.
  * `bidsification/exclusions.py` is applied as a belt-and-braces subject filter.
  * Env: any pandas env (clinical / snp). Local Windows.

Usage
-----
  python clinical_pipeline/01i_build_longitudinal_visit_splits.py [--check]
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# bidsification/exclusions.py — single source of truth for subject exclusions
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bidsification"))
from exclusions import is_excluded_subject  # noqa: E402

DERIV = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion")
LONG_MASTER = DERIV / "verbose" / "longitudinal" / "master_clinical_verbose.csv"
CANON_SPLIT_DIR = DERIV / "tabular" / "baseline"          # Patient_ID -> fold (shared)
OUT_DIR = DERIV / "verbose" / "longitudinal_allvisits"
SEEDS = (0, 1, 2)
SPLITS = ("train", "val", "test")
ID_COL = "Patient_ID"


def viscode_to_months(v) -> float:
    """ADNI VISCODE_2 -> months from baseline. bl=0, sc/screening=-1, m<NN>=NN; else NaN."""
    s = str(v).strip().lower()
    if s in ("bl", "baseline"):
        return 0.0
    if s in ("sc", "scr", "scmri", "screening"):
        return -1.0
    m = re.fullmatch(r"m(\d+)", s)
    return float(m.group(1)) if m else np.nan


def canonical_folds(seed: int) -> dict:
    """Patient_ID -> split for one seed, from the canonical baseline partition."""
    fold = {}
    for sp in SPLITS:
        ids = pd.read_csv(CANON_SPLIT_DIR / f"seed_{seed}" / f"{sp}.csv",
                          usecols=[ID_COL], low_memory=False)[ID_COL].astype(str)
        for pid in ids:
            fold[pid] = sp
    return fold


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="Dry run: report counts, do not write split CSVs.")
    args = ap.parse_args()

    print(f"Reading longitudinal master: {LONG_MASTER}")
    df = pd.read_csv(LONG_MASTER, low_memory=False)
    df[ID_COL] = df[ID_COL].astype(str)
    n0 = len(df)
    print(f"  {n0} visit rows, {df[ID_COL].nunique()} patients")

    df["visit_months"] = df["VISCODE_2"].map(viscode_to_months)

    # belt-and-braces subject exclusion (canonical split already encodes the cohort)
    excl_mask = df[ID_COL].apply(is_excluded_subject)
    if excl_mask.any():
        print(f"  exclusions.py drops {excl_mask.sum()} rows "
              f"({df.loc[excl_mask, ID_COL].nunique()} subjects)")
    df = df[~excl_mask].reset_index(drop=True)

    inwin_total = df["visit_months"].between(0, 12).sum()
    print(f"  in-window (0..12mo) rows available: {inwin_total} "
          f"(viscodes {sorted(df.loc[df['visit_months'].between(0,12),'VISCODE_2'].unique())})")
    print(f"  rows with NaN Label_visit_diag (kept; training drops them): "
          f"{df['Label_visit_diag'].isna().sum()}")

    for seed in SEEDS:
        fold = canonical_folds(seed)
        df["_split"] = df[ID_COL].map(fold)
        unmapped = df["_split"].isna()
        n_unmapped_pt = df.loc[unmapped, ID_COL].nunique()
        seed_df = df[~unmapped].copy()

        # verify patient-disjoint folds
        sets = {sp: set(seed_df.loc[seed_df["_split"] == sp, ID_COL]) for sp in SPLITS}
        overlap = (sets["train"] & sets["val"]) | (sets["train"] & sets["test"]) | \
                  (sets["val"] & sets["test"])
        assert not overlap, f"seed {seed}: patient overlap across folds: {list(overlap)[:5]}"

        print(f"\nseed {seed}: dropped {n_unmapped_pt} patients not in canonical split")
        for sp in SPLITS:
            part = seed_df[seed_df["_split"] == sp].drop(columns="_split")
            iw = part[part["visit_months"].between(0, 12)]
            lab = part["Label_visit_diag"].value_counts().to_dict()
            print(f"  {sp:5s}: {len(part):4d} rows / {part[ID_COL].nunique():3d} pts "
                  f"| in-window {len(iw):4d} ({iw['VISCODE_2'].value_counts().to_dict()}) "
                  f"| dx {lab}")
            if not args.check:
                out = OUT_DIR / f"seed_{seed}" / f"{sp}.csv"
                out.parent.mkdir(parents=True, exist_ok=True)
                part.to_csv(out, index=False)

    df.drop(columns="_split", errors="ignore", inplace=True)
    if args.check:
        print("\n[--check] no files written.")
    else:
        print(f"\nWrote splits -> {OUT_DIR}")


if __name__ == "__main__":
    main()
