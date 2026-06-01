"""
01e_build_T4_labels_and_splits.py
===================================
Build T4_conv_horizon (3-class ordinal: <3y / 3-7y / >7y) labels and
stratified 80/10/10 splits for the converter cohort (pMCI + pCN_to_AD,
146 subjects).

What this script does:

1. Reads the post-exclusion MRI master CSV (the canonical per-scan
   anchor with conversion_group + years_to_AD columns).
2. Filters to `conversion_group in {pMCI, pCN_to_AD}` and drops NaN
   `years_to_AD`. T4 is well-defined ONLY for confirmed converters --
   sMCI/sCN/AD_bl are excluded.
3. Computes `Label_T4` per subject using the user-confirmed 3-bucket
   scheme:
      0 : years_to_AD <  3
      1 : 3 <= years_to_AD <  7
      2 : years_to_AD >= 7
4. Writes per-subject Label_T4 back to the MRI master (in-place column
   addition). Non-T4-eligible subjects get NaN.
5. Generates 3-seed stratified 80/10/10 splits on Label_T4 and writes
   them to a new sister dir:
      derivatives/clinical/no_cdr_stratified_post_exclusion/
        tabular/baseline_T4/seed_{0,1,2}/{train,val,test}.csv

Each split CSV is the bl-row subset of the MRI master columns for the
T4 subjects in that split, plus the new Label_T4 column.

Stratification: on Label_T4 alone (not on conversion_group too) because
the smallest pCN_to_AD bucket has only 1 subject which breaks
per-(bucket, group) stratification. The conversion_group distribution
will fall where it falls; document the per-split breakdown for
transparency.

Distribution (verified):

| Class | n_total | pMCI | pCN_to_AD |
|---|---|---|---|
| 0 (<3y) | 63 | 62 | 1 |
| 1 (3-7y) | 47 | 42 | 5 |
| 2 (>7y) | 36 | 17 | 19 |
| Total | 146 | 121 | 25 |

Inputs:
  MRI master:
    D:/ADNI_BIDS_project/derivatives/mri_clinical_matched/
    VISCODE_2_aligned_extended_post_exclusion/
    master_mri_clinical_matched_viscode2_extended_post_exclusion.csv

Outputs (idempotent):
  - MRI master (in-place addition of Label_T4 column)
  - derivatives/clinical/no_cdr_stratified_post_exclusion/tabular/
      baseline_T4/seed_{0,1,2}/{train,val,test}.csv

Usage:
  python clinical_pipeline/01e_build_T4_labels_and_splits.py
  python clinical_pipeline/01e_build_T4_labels_and_splits.py --dry-run
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
MRI_MASTER = Path(
    r"D:/ADNI_BIDS_project/derivatives/mri_clinical_matched/"
    r"VISCODE_2_aligned_extended_post_exclusion/"
    r"master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
)
SPLIT_OUT_BASE = Path(
    r"D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified_post_exclusion/"
    r"tabular/baseline_T4"
)

SEEDS = [0, 1, 2]

T4_COHORT_GROUPS = ["pMCI", "pCN_to_AD"]


# --------------------------------------------------------------------------- #
# Bucket function
# --------------------------------------------------------------------------- #
def label_T4(years_to_ad: float) -> int | None:
    """User-confirmed 3-bucket scheme: <3 / 3-7 / >=7 years."""
    if years_to_ad is None or (isinstance(years_to_ad, float) and np.isnan(years_to_ad)):
        return None
    try:
        y = float(years_to_ad)
    except (TypeError, ValueError):
        return None
    if y < 3:
        return 0
    if y < 7:
        return 1
    return 2


# --------------------------------------------------------------------------- #
# Steps
# --------------------------------------------------------------------------- #
def compute_label_T4(master: pd.DataFrame) -> pd.DataFrame:
    """Return master with a new Label_T4 column.
    NaN for non-T4 subjects; 0/1/2 for the 146 T4 subjects."""
    out = master.copy()
    out["Label_T4"] = pd.NA
    eligible = out["conversion_group"].isin(T4_COHORT_GROUPS)
    yta = pd.to_numeric(out.loc[eligible, "years_to_AD"], errors="coerce")
    out.loc[eligible, "Label_T4"] = yta.apply(label_T4)
    return out


def build_subject_table(master: pd.DataFrame) -> pd.DataFrame:
    """One row per subject from the master.

    Preference: the row with `adni_viscode == 'bl'` if present, else the
    earliest `mri_date`. We need one well-defined row per subject because
    Label_T4 is subject-level and the trainer's `baseline_anchored`
    session policy filters to bl scans.
    """
    m = master[master["Label_T4"].notna()].copy()
    m["mri_date_dt"] = pd.to_datetime(m["mri_date"], errors="coerce")
    # Sort with bl rows first, then by date.
    m["_bl_priority"] = (m["adni_viscode"] == "bl").astype(int)
    m = m.sort_values(["Patient_ID", "_bl_priority", "mri_date_dt"],
                      ascending=[True, False, True])
    one_per_sub = m.drop_duplicates("Patient_ID", keep="first").copy()
    one_per_sub = one_per_sub.drop(columns=["mri_date_dt", "_bl_priority"])
    one_per_sub["Label_T4"] = one_per_sub["Label_T4"].astype(int)
    return one_per_sub


def make_splits(subject_table: pd.DataFrame, seed: int) -> dict[str, pd.DataFrame]:
    """Stratified 80/10/10 on Label_T4."""
    y = subject_table["Label_T4"]
    # 80/20 first
    train_df, rest_df = train_test_split(
        subject_table, test_size=0.20, random_state=seed,
        stratify=y, shuffle=True,
    )
    # 10/10 from the held-out 20%
    val_df, test_df = train_test_split(
        rest_df, test_size=0.50, random_state=seed,
        stratify=rest_df["Label_T4"], shuffle=True,
    )
    return {"train": train_df.reset_index(drop=True),
            "val":   val_df.reset_index(drop=True),
            "test":  test_df.reset_index(drop=True)}


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--master-csv", type=Path, default=MRI_MASTER)
    p.add_argument("--out-base",   type=Path, default=SPLIT_OUT_BASE)
    p.add_argument("--dry-run",    action="store_true")
    args = p.parse_args()

    print(f"[01e] Reading MRI master: {args.master_csv}")
    if not args.master_csv.exists():
        print(f"[ERROR] missing: {args.master_csv}", file=sys.stderr); return 1

    master = pd.read_csv(args.master_csv, low_memory=False)
    print(f"  master rows: {len(master)}, unique Patient_ID: {master.Patient_ID.nunique()}")

    # Step 1: compute Label_T4 column
    master_with_t4 = compute_label_T4(master)
    eligible = master_with_t4["Label_T4"].notna().sum()
    print(f"  rows with Label_T4 set: {eligible} (of {len(master_with_t4)})")

    # Step 2: build per-subject table (one row per subject)
    subj = build_subject_table(master_with_t4)
    print(f"  unique subjects in T4 cohort: {len(subj)}")
    print(f"  Label_T4 distribution:")
    print(subj["Label_T4"].value_counts().sort_index().to_string())
    print()
    print(f"  conversion_group x Label_T4 cross-tab:")
    print(pd.crosstab(subj["Label_T4"], subj["conversion_group"], margins=True))

    if args.dry_run:
        print("\n[dry-run] not writing outputs.")
        return 0

    # Step 3: write MRI master back IN-PLACE with the new Label_T4 column
    print(f"\n[01e] Writing MRI master back with Label_T4 column: {args.master_csv}")
    master_with_t4.to_csv(args.master_csv, index=False)

    # Step 4: stratified splits per seed
    for seed in SEEDS:
        splits = make_splits(subj, seed=seed)
        seed_dir = args.out_base / f"seed_{seed}"
        seed_dir.mkdir(parents=True, exist_ok=True)
        for split_name, df in splits.items():
            out_path = seed_dir / f"{split_name}.csv"
            df.to_csv(out_path, index=False)
        print(f"\n  seed_{seed}:")
        for split_name, df in splits.items():
            dist = df["Label_T4"].value_counts().sort_index().to_dict()
            group_dist = df["conversion_group"].value_counts().to_dict()
            print(f"    {split_name}: n={len(df)}, "
                  f"Label_T4={dist}, conversion_group={group_dist}")

    print(f"\n[ok] T4 splits written to: {args.out_base}")
    print(f"[ok] MRI master extended with Label_T4 column.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
