"""
01h_build_m12_crosssectional_splits.py
======================================
Materialise VERBOSE (text) splits for the **1-year (m12) cross-sectional** clinical-
encoder tasks, REUSING the canonical `verbose/baseline` per-seed fold assignment so the
m12 tasks share the EXACT SAME per-patient train/val/test partition as the baseline tasks
(and therefore the MRI arm, once it adopts the same split) -- but swapping the model
INPUT to the patient's m12-visit clinical summary and the diagnosis LABEL to the
**concurrent** m12 diagnosis (`Label_visit_diag`).

Why this script exists
----------------------
The clinical encoders were trained baseline-only (`verbose/baseline` = 1 row/patient at
VISCODE=bl, label `Label_bl_multi`). For late fusion with MRI -- which is dominated by
m12+ follow-up scans -- we need clinical probabilities at the 1-year visit. This builds a
parallel `verbose/baseline_m12` tree:
  - INPUT  : the m12-visit `Generated_Text` (from verbose/by_visit/visit_m12.csv)
  - LABEL  : `Label_visit_diag` (concurrent diagnosis AT m12; AD count rises vs baseline
             as patients convert) for the diagnosis tasks T2/T1/T1b
  - T1d    : reuses `conversion_group` (carried over from baseline; time-invariant) -- only
             the input timepoint changes.
Patients in the baseline split who have NO m12 visit are DROPPED (and logged); the SAME
per-seed split membership is otherwise preserved (leakage-safe & MRI-alignable).

Inputs
  verbose/baseline/seed_{0,1,2}/{train,val,test}.csv   (Patient_ID, split membership,
                                                        conversion_group, Label_bl_multi, ...)
  verbose/by_visit/visit_m12.csv                       (Patient_ID, Generated_Text [m12],
                                                        Label_visit_diag [concurrent], VISCODE_2)

Output (same fold structure; baseline bookkeeping cols + m12 text + concurrent label)
  verbose/baseline_m12/seed_{0,1,2}/{train,val,test}.csv

Env: clinical or snp (pandas only).
Run:
  python clinical_pipeline/01h_build_m12_crosssectional_splits.py
  python clinical_pipeline/01h_build_m12_crosssectional_splits.py --dry-run
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
VERBOSE_BASELINE = BASE / "verbose" / "baseline"          # fold assignment + bookkeeping
VERBOSE_M12_OUT = BASE / "verbose" / "baseline_m12"       # output
# m12 per-visit text lives in the (pre-exclusion) verbose tree; it is per-(patient,visit)
M12_BY_VISIT = Path(
    r"D:/ADNI_BIDS_project/derivatives/clinical/verbose/by_visit/visit_m12.csv")

SEEDS = [0, 1, 2]
SPLITS = ["train", "val", "test"]
ID_COL = "Patient_ID"
TEXT_COL = "Generated_Text"
CONCURRENT_LABEL = "Label_visit_diag"   # concurrent dx at m12 (used by T2/T1/T1b m12 tasks)
# bookkeeping columns carried over from the baseline split (per-patient, seed-independent)
CARRY_COLS = ["Label_bl_multi", "conversion_group", "years_to_AD", "FU_years"]


def build_m12_map(by_visit_csv: Path) -> pd.DataFrame:
    """Patient_ID -> (m12 Generated_Text, concurrent Label_visit_diag).

    visit_m12.csv groups by the nominal 1-year visit; ~15 patients have both an m06 and
    an m12 ADNI VISCODE_2 row. Prefer the exact `m12` row, else the first available."""
    m = pd.read_csv(by_visit_csv, low_memory=False)
    m[ID_COL] = m[ID_COL].astype(str)
    # rank: exact 'm12' VISCODE_2 first (0), everything else after (1)
    m["_rank"] = (m["VISCODE_2"].astype(str) != "m12").astype(int)
    m = m.sort_values([ID_COL, "_rank"], kind="stable")
    m = m.dropna(subset=[TEXT_COL]).drop_duplicates(ID_COL, keep="first")
    keep = [ID_COL, TEXT_COL, CONCURRENT_LABEL, "VISCODE_2"]
    return m[keep].set_index(ID_COL)


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--verbose-baseline", type=Path, default=VERBOSE_BASELINE)
    p.add_argument("--m12-by-visit", type=Path, default=M12_BY_VISIT)
    p.add_argument("--out-base", type=Path, default=VERBOSE_M12_OUT)
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    if not args.verbose_baseline.exists():
        print(f"[ERROR] missing baseline split dir: {args.verbose_baseline}", file=sys.stderr)
        return 1
    if not args.m12_by_visit.exists():
        print(f"[ERROR] missing m12 by-visit file: {args.m12_by_visit}", file=sys.stderr)
        return 1

    print(f"[01h] m12 text/label map from: {args.m12_by_visit}")
    m12 = build_m12_map(args.m12_by_visit)
    print(f"  m12 visit available for {len(m12)} subjects "
          f"(concurrent {CONCURRENT_LABEL} dist: "
          f"{m12[CONCURRENT_LABEL].value_counts(dropna=False).to_dict()})")

    total_kept = total_dropped = 0
    dropped_ids: list[str] = []
    for seed in SEEDS:
        print(f"\n  seed_{seed}:")
        for sp in SPLITS:
            src = args.verbose_baseline / f"seed_{seed}" / f"{sp}.csv"
            if not src.exists():
                print(f"[ERROR] missing split: {src}", file=sys.stderr)
                return 1
            base = pd.read_csv(src, low_memory=False)
            base[ID_COL] = base[ID_COL].astype(str)

            # carry the per-patient bookkeeping cols that exist in the baseline split
            carry = [c for c in CARRY_COLS if c in base.columns]
            out = base[[ID_COL] + carry].copy()

            # attach m12 input text + concurrent label (drop patients with no m12 visit)
            out = out.join(m12, on=ID_COL)
            has_m12 = out[TEXT_COL].notna()
            n_drop = int((~has_m12).sum())
            if n_drop:
                dropped_ids.extend(out.loc[~has_m12, ID_COL].tolist())
            kept = out[has_m12].copy()
            kept["VISCODE_long"] = "m12"
            kept = kept.reset_index(drop=True)

            total_kept += len(kept)
            total_dropped += n_drop
            dx = kept[CONCURRENT_LABEL].value_counts(dropna=False).sort_index().to_dict()
            print(f"    {sp:5s}: {len(kept):3d}/{len(base):3d} kept "
                  f"(dropped {n_drop} no-m12)  {CONCURRENT_LABEL}={dx}")
            if not args.dry_run:
                out_dir = args.out_base / f"seed_{seed}"
                out_dir.mkdir(parents=True, exist_ok=True)
                kept.to_csv(out_dir / f"{sp}.csv", index=False)

    uniq_dropped = sorted(set(dropped_ids))
    print(f"\n  Kept {total_kept} rows; dropped {total_dropped} no-m12 rows "
          f"({len(uniq_dropped)} unique subjects with no m12 visit).")
    if uniq_dropped:
        print(f"  Dropped (no-m12) Patient_IDs: {uniq_dropped}")

    if args.dry_run:
        print("\n[dry-run] not writing outputs.")
        return 0
    print(f"\n[ok] Verbose m12 cross-sectional splits written to: {args.out_base}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
