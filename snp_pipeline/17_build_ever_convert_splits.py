"""
17_build_ever_convert_splits.py
================================
Generate a "true-CN vs ever-MCI/AD" label set that captures longitudinal
trajectory across the entire ADNI follow-up, rather than just the baseline
diagnosis.

Motivation
----------
The default ``Label_bl_multi`` in the clinical pipeline labels CN-at-baseline
patients as CN even if they later convert to MCI or AD within the study.
For prognostic SNP classification we'd prefer:
    label 1 = patient EVER shows non-CN diagnosis at any follow-up visit
    label 0 = patient remains CN across the entire follow-up

Cohort breakdown (computed once over all 616 patients with both SNP and MRI):
    bl_dx  ever_non_cn      0    1
    CN                    146   66      (66/212 = 31% of baseline-CN convert)
    MCI                     0  368
    AD                      0   36
                          146  470

Output
------
Mirrors the existing no_cdr_stratified split layout, but adds a new column
``Label_ever_convert`` (0/1) to each split CSV. **Patient assignments to
train/val/test are kept identical to the input splits** so that test subjects
remain consistent across SNP / MRI / clinical modalities (per the
cross-modality contract).

Output tree:
    <out_root>/tabular/baseline/seed_{0,1,2}/{train,val,test}.csv

Usage
-----
    python snp_pipeline/17_build_ever_convert_splits.py
    python snp_pipeline/17_build_ever_convert_splits.py \\
        --longitudinal D:/ADNI_BIDS_project/derivatives/clinical/tabular/longitudinal/master_clinical_tabular.csv \\
        --in-splits-root D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified/tabular/baseline \\
        --out-root D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified_ever_convert
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

DEFAULT_LONG = Path(
    r"D:/ADNI_BIDS_project/derivatives/clinical/tabular/longitudinal/master_clinical_tabular.csv"
)
DEFAULT_IN_SPLITS = Path(
    r"D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified/tabular/baseline"
)
DEFAULT_OUT = Path(
    r"D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified_ever_convert"
)


def build_ever_convert(long_csv: Path) -> pd.DataFrame:
    """Per-patient table with the ever_convert label."""
    df = pd.read_csv(long_csv, low_memory=False)
    df["Patient_ID"] = df["Patient_ID"].astype(str)

    per_pt = (
        df.groupby("Patient_ID")
        .agg(
            n_visits=("VISCODE_2", "size"),
            bl_dx=("Label_bl_multi", "first"),
            ever_non_cn=("Label_visit_diag", lambda s: int((s.fillna("CN") != "CN").any())),
            ever_non_cn_granular=(
                "Label_visit_granular", lambda s: int((s.fillna("CN") != "CN").any())
            ),
            last_diag=("Label_visit_diag", lambda s: s.dropna().iloc[-1] if s.dropna().any() else "NA"),
        )
        .reset_index()
    )
    per_pt["Label_ever_convert"] = per_pt["ever_non_cn"].astype(int)
    return per_pt


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--longitudinal", type=Path, default=DEFAULT_LONG)
    p.add_argument("--in-splits-root", type=Path, default=DEFAULT_IN_SPLITS,
                   help="Existing split root (no_cdr_stratified/tabular/baseline)")
    p.add_argument("--out-root", type=Path, default=DEFAULT_OUT,
                   help="Output root (will mirror the in-splits structure)")
    p.add_argument("--check", action="store_true",
                   help="Dry-run: print stats only.")
    args = p.parse_args()

    print(f"Reading longitudinal data: {args.longitudinal}")
    per_pt = build_ever_convert(args.longitudinal)
    print(f"  {len(per_pt)} patients")
    print(f"\nBaseline diagnosis × ever_convert:")
    print(pd.crosstab(per_pt["bl_dx"], per_pt["Label_ever_convert"], margins=True))

    pt_to_label = dict(zip(per_pt["Patient_ID"], per_pt["Label_ever_convert"]))
    pt_to_last  = dict(zip(per_pt["Patient_ID"], per_pt["last_diag"]))

    out_base = args.out_root / "tabular" / "baseline"
    for seed in (0, 1, 2):
        for split in ("train", "val", "test"):
            in_csv = args.in_splits_root / f"seed_{seed}" / f"{split}.csv"
            if not in_csv.exists():
                print(f"  [warn] missing {in_csv}")
                continue
            df = pd.read_csv(in_csv, low_memory=False)
            df["Patient_ID"] = df["Patient_ID"].astype(str)
            df["Label_ever_convert"] = df["Patient_ID"].map(pt_to_label)
            df["Label_last_visit_diag"] = df["Patient_ID"].map(pt_to_last)
            miss = df["Label_ever_convert"].isna().sum()
            if miss:
                print(f"  [warn] seed_{seed}/{split}: {miss} patients with no longitudinal data")
            df["Label_ever_convert"] = df["Label_ever_convert"].fillna(0).astype(int)

            if args.check:
                n_pos = int(df["Label_ever_convert"].sum())
                n_neg = len(df) - n_pos
                print(f"  seed_{seed}/{split}: n={len(df)} pos={n_pos} neg={n_neg} "
                      f"(pos_frac={n_pos/len(df):.3f})")
                continue

            out_csv = out_base / f"seed_{seed}" / f"{split}.csv"
            out_csv.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(out_csv, index=False)

    if not args.check:
        print(f"\nWrote split CSVs under {out_base}")


if __name__ == "__main__":
    main()
