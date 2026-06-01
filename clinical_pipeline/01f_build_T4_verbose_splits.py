"""
01f_build_T4_verbose_splits.py
==============================
Materialise VERBOSE (text) splits for the T4_conv_horizon clinical-encoder task,
REUSING the tabular/baseline_T4 fold assignment from 01e so the clinical T4 task
shares the SAME per-seed folds as the MRI T4 task -- but adding the Generated_Text
the encoder needs.

Why this script exists
----------------------
`tabular/baseline_T4` (built by 01e from the MRI master) is tabular-only and has
NO `Generated_Text`. The clinical `verbose/baseline` (571 subjects, all groups)
DOES have `Generated_Text` (per subject, seed-independent). This script joins that
text onto the baseline_T4 Patient_IDs, KEEPING the baseline_T4 fold, and drops
(and logs) the converters that have no clinical text -- i.e. MRI-matched subjects
that are absent from the clinical text cohort (~24 of 146). Those subjects simply
lack the clinical modality and are handled per-modality by the fusion step.

Inputs
  tabular/baseline_T4/seed_{0,1,2}/{train,val,test}.csv   (Patient_ID, Label_T4, ...; from 01e)
  verbose/baseline/seed_{0,1,2}/{train,val,test}.csv      (Patient_ID, Generated_Text, ...)

Output (mirrors the input fold structure, baseline_T4 columns + Generated_Text)
  verbose/baseline_T4/seed_{0,1,2}/{train,val,test}.csv

Env: clinical or snp (pandas only).
Run:
  python clinical_pipeline/01f_build_T4_verbose_splits.py
  python clinical_pipeline/01f_build_T4_verbose_splits.py --dry-run
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
BASE = Path(
    r"D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified_post_exclusion")
T4_TABULAR = BASE / "tabular" / "baseline_T4"     # input folds (01e)
VERBOSE_BASELINE = BASE / "verbose" / "baseline"  # text source
VERBOSE_T4_OUT = BASE / "verbose" / "baseline_T4"  # output

SEEDS = [0, 1, 2]
SPLITS = ["train", "val", "test"]
ID_COL, TEXT_COL = "Patient_ID", "Generated_Text"


def build_text_map(verbose_baseline: Path) -> pd.Series:
    """Patient_ID -> Generated_Text from the existing verbose/baseline splits.
    Generated_Text is a per-subject baseline summary (seed-independent), so
    concatenating seed_0's three folds and de-duplicating covers every subject."""
    frames = []
    for sp in SPLITS:
        df = pd.read_csv(verbose_baseline / "seed_0" / f"{sp}.csv",
                         usecols=[ID_COL, TEXT_COL], low_memory=False)
        frames.append(df)
    allrows = pd.concat(frames, ignore_index=True)
    allrows[ID_COL] = allrows[ID_COL].astype(str)
    allrows = allrows.dropna(subset=[TEXT_COL]).drop_duplicates(ID_COL, keep="first")
    return allrows.set_index(ID_COL)[TEXT_COL]


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--t4-tabular", type=Path, default=T4_TABULAR)
    p.add_argument("--verbose-baseline", type=Path, default=VERBOSE_BASELINE)
    p.add_argument("--out-base", type=Path, default=VERBOSE_T4_OUT)
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    for d in (args.t4_tabular, args.verbose_baseline):
        if not d.exists():
            print(f"[ERROR] missing input dir: {d}", file=sys.stderr)
            return 1

    print(f"[01f] Text map from: {args.verbose_baseline}")
    text_map = build_text_map(args.verbose_baseline)
    print(f"  Generated_Text available for {len(text_map)} subjects")

    total_kept = total_dropped = 0
    dropped_ids: list[str] = []
    for seed in SEEDS:
        print(f"\n  seed_{seed}:")
        for sp in SPLITS:
            src = args.t4_tabular / f"seed_{seed}" / f"{sp}.csv"
            if not src.exists():
                print(f"[ERROR] missing split: {src}", file=sys.stderr)
                return 1
            df = pd.read_csv(src, low_memory=False)
            df[ID_COL] = df[ID_COL].astype(str)
            df[TEXT_COL] = df[ID_COL].map(text_map)
            has_text = df[TEXT_COL].notna()
            n_drop = int((~has_text).sum())
            if n_drop:
                dropped_ids.extend(df.loc[~has_text, ID_COL].tolist())
            kept = df[has_text].reset_index(drop=True)
            total_kept += len(kept)
            total_dropped += n_drop
            dist = kept["Label_T4"].value_counts().sort_index().to_dict()
            print(f"    {sp:5s}: {len(kept):3d}/{len(df):3d} kept "
                  f"(dropped {n_drop} text-less)  Label_T4={dist}")
            if not args.dry_run:
                out_dir = args.out_base / f"seed_{seed}"
                out_dir.mkdir(parents=True, exist_ok=True)
                kept.to_csv(out_dir / f"{sp}.csv", index=False)

    uniq_dropped = sorted(set(dropped_ids))
    print(f"\n  Kept {total_kept} rows; dropped {total_dropped} text-less rows "
          f"({len(uniq_dropped)} unique subjects with no clinical text).")
    if uniq_dropped:
        print(f"  Dropped (text-less) Patient_IDs: {uniq_dropped}")

    if args.dry_run:
        print("\n[dry-run] not writing outputs.")
        return 0
    print(f"\n[ok] Verbose T4 splits written to: {args.out_base}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
