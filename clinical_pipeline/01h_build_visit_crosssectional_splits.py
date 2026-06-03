"""
01h_build_visit_crosssectional_splits.py
=========================================
Materialise VERBOSE (text) splits for a SINGLE-VISIT cross-sectional clinical-encoder task
(default m12 = 1 year), REUSING the canonical `verbose/baseline` per-seed fold assignment so
the visit task shares the EXACT SAME per-patient train/val/test partition as the baseline
tasks (and therefore the MRI arm, once it adopts the same split) -- but swapping the model
INPUT to the patient's note at that visit and the diagnosis LABEL to the **concurrent**
diagnosis at that visit (`Label_visit_diag`).

Why this script exists
----------------------
The clinical encoders were trained baseline-only (`verbose/baseline` = 1 row/patient at
VISCODE=bl, label `Label_bl_multi`). For late fusion with MRI -- which is dominated by m12+
follow-up scans -- and for the temporal clinical EN (an elastic net over a patient's visits
up to year 1: bl + m06 + m12), we need clinical probabilities at individual follow-up visits.
This builds a parallel `verbose/baseline_<visit>` tree:
  - INPUT  : the <visit> `Generated_Text` (from verbose/by_visit/visit_<visit>.csv)
  - LABEL  : `Label_visit_diag` (concurrent diagnosis AT that visit) for the diagnosis tasks
             T2/T1/T1b
  - T1d    : reuses `conversion_group` (carried over from baseline; time-invariant) -- only
             the input timepoint changes.
Patients in the baseline split who have NO note at that visit are DROPPED (and logged); the
SAME per-seed split membership is otherwise preserved (leakage-safe & MRI-alignable).

Realistic ADNI clinical visits up to 1 year (canonical 571-patient cohort coverage):
  bl=571, m06=518, m12=529  (m03 is negligible/spurious; sc is pre-baseline screening).

Inputs
  verbose/baseline/seed_{0,1,2}/{train,val,test}.csv   (Patient_ID, split membership,
                                                        conversion_group, Label_bl_multi, ...)
  verbose/by_visit/visit_<visit>.csv                   (Patient_ID, Generated_Text [<visit>],
                                                        Label_visit_diag [concurrent], VISCODE_2)

Output (same fold structure; baseline bookkeeping cols + visit text + concurrent label)
  verbose/baseline_<visit>/seed_{0,1,2}/{train,val,test}.csv

Env: clinical or snp (pandas only).
Run:
  python clinical_pipeline/01h_build_visit_crosssectional_splits.py                 # m12 (default)
  python clinical_pipeline/01h_build_visit_crosssectional_splits.py --visit m06
  python clinical_pipeline/01h_build_visit_crosssectional_splits.py --visit m06 --dry-run
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
# per-visit text lives in the (pre-exclusion) verbose tree; it is per-(patient,visit)
BY_VISIT_DIR = Path(r"D:/ADNI_BIDS_project/derivatives/clinical/verbose/by_visit")

SEEDS = [0, 1, 2]
SPLITS = ["train", "val", "test"]
ID_COL = "Patient_ID"
TEXT_COL = "Generated_Text"
CONCURRENT_LABEL = "Label_visit_diag"   # concurrent dx at the visit (used by T2/T1/T1b tasks)
# bookkeeping columns carried over from the baseline split (per-patient, seed-independent)
CARRY_COLS = ["Label_bl_multi", "conversion_group", "years_to_AD", "FU_years"]


def build_visit_map(by_visit_csv: Path, visit: str) -> pd.DataFrame:
    """Patient_ID -> (<visit> Generated_Text, concurrent Label_visit_diag).

    A visit_<v>.csv groups by the nominal visit; a few patients have an adjacent ADNI
    VISCODE_2 row (e.g. m06 alongside m12). Prefer the exact `<visit>` VISCODE_2 row,
    else the first available."""
    m = pd.read_csv(by_visit_csv, low_memory=False)
    m[ID_COL] = m[ID_COL].astype(str)
    # rank: exact-visit VISCODE_2 first (0), everything else after (1)
    m["_rank"] = (m["VISCODE_2"].astype(str) != visit).astype(int)
    m = m.sort_values([ID_COL, "_rank"], kind="stable")
    m = m.dropna(subset=[TEXT_COL]).drop_duplicates(ID_COL, keep="first")
    keep = [ID_COL, TEXT_COL, CONCURRENT_LABEL, "VISCODE_2"]
    return m[keep].set_index(ID_COL)


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--visit", default="m12",
                   help="ADNI nominal visit to build (default m12). Usable up to 1y: m06, m12.")
    p.add_argument("--verbose-baseline", type=Path, default=VERBOSE_BASELINE)
    p.add_argument("--by-visit-dir", type=Path, default=BY_VISIT_DIR)
    p.add_argument("--out-base", type=Path, default=None,
                   help="Output dir (default: <BASE>/verbose/baseline_<visit>).")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    visit = args.visit
    by_visit_csv = args.by_visit_dir / f"visit_{visit}.csv"
    out_base = args.out_base or (BASE / "verbose" / f"baseline_{visit}")

    if not args.verbose_baseline.exists():
        print(f"[ERROR] missing baseline split dir: {args.verbose_baseline}", file=sys.stderr)
        return 1
    if not by_visit_csv.exists():
        print(f"[ERROR] missing visit by-visit file: {by_visit_csv}", file=sys.stderr)
        return 1

    print(f"[01h] visit={visit}  text/label map from: {by_visit_csv}")
    vmap = build_visit_map(by_visit_csv, visit)
    print(f"  {visit} visit available for {len(vmap)} subjects "
          f"(concurrent {CONCURRENT_LABEL} dist: "
          f"{vmap[CONCURRENT_LABEL].value_counts(dropna=False).to_dict()})")

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

            # attach visit input text + concurrent label (drop patients with no note at this visit)
            out = out.join(vmap, on=ID_COL)
            has_visit = out[TEXT_COL].notna()
            n_drop = int((~has_visit).sum())
            if n_drop:
                dropped_ids.extend(out.loc[~has_visit, ID_COL].tolist())
            kept = out[has_visit].copy()
            kept["VISCODE_long"] = visit
            kept = kept.reset_index(drop=True)

            total_kept += len(kept)
            total_dropped += n_drop
            dx = kept[CONCURRENT_LABEL].value_counts(dropna=False).sort_index().to_dict()
            print(f"    {sp:5s}: {len(kept):3d}/{len(base):3d} kept "
                  f"(dropped {n_drop} no-{visit})  {CONCURRENT_LABEL}={dx}")
            if not args.dry_run:
                out_dir = out_base / f"seed_{seed}"
                out_dir.mkdir(parents=True, exist_ok=True)
                kept.to_csv(out_dir / f"{sp}.csv", index=False)

    uniq_dropped = sorted(set(dropped_ids))
    print(f"\n  Kept {total_kept} rows; dropped {total_dropped} no-{visit} rows "
          f"({len(uniq_dropped)} unique subjects with no {visit} visit).")
    if uniq_dropped:
        print(f"  Dropped (no-{visit}) Patient_IDs: {uniq_dropped}")

    if args.dry_run:
        print("\n[dry-run] not writing outputs.")
        return 0
    print(f"\n[ok] Verbose {visit} cross-sectional splits written to: {out_base}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
