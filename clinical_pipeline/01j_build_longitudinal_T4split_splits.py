r"""
01j_build_longitudinal_T4split_splits.py
========================================
Build per-seed ALL-VISITS longitudinal splits whose fold assignment is **T4-aware**, so a SECOND
longitudinal T2 encoder (re-extraction) yields per-visit embeddings that are LEAKAGE-SAFE for the
downstream T4 horizon CLASSIFIER on the balanced `baseline_T4` split.

The problem 01i did NOT solve for the classifier: the canonical-split encoder was trained on
canonical-train, so a `baseline_T4`-test converter (different stratification) can be in canonical-train →
its embedding is in-sample → leakage if the classifier reports on `baseline_T4`-test.

Fix here: assign every visit row's fold by
  • CONVERTER (conversion_group ∈ {pMCI, pCN_to_AD}, the 146 T4-eligible) → its **`baseline_T4`** fold;
  • non-converter                                                          → its **canonical** fold.
Then every `baseline_T4` val/test converter is OUT-of-fold for the encoder, and the classifier can train
on `baseline_T4`-train / select on val / report on test with NO encoder leakage. Non-converters (never in
the T4 eval) fill the encoder's train/val for representation breadth.

Inputs
------
  Longitudinal master (1 row/patient-visit, has Generated_Text + Label_visit_diag + conversion_group):
    …\no_cdr_stratified_post_exclusion\verbose\longitudinal\master_clinical_verbose.csv
  Converter folds (Label_T4-stratified 80/10/10):  …\tabular\baseline_T4\seed_{0,1,2}\{train,val,test}.csv
  Canonical folds (shared):                        …\tabular\baseline\seed_{0,1,2}\{train,val,test}.csv

Output
------
  …\verbose\longitudinal_allvisits_T4split\seed_{0,1,2}\{train,val,test}.csv  (all visits + visit_months)

Then RE-EXTRACT with the existing 03j_finetune_longitudinal_extract.py pointed at this data_dir + a new
out root (encoder_outputs_no_cdr_post_exclusion_longitudinal_T4split) — see
03j_finetune_longitudinal_T4split_submit_csd3.sh.

Usage:  python clinical_pipeline/01j_build_longitudinal_T4split_splits.py [--check]
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bidsification"))
from exclusions import is_excluded_subject  # noqa: E402

DERIV = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion")
LONG_MASTER = DERIV / "verbose" / "longitudinal" / "master_clinical_verbose.csv"
CANON_SPLIT_DIR = DERIV / "tabular" / "baseline"
T4_SPLIT_DIR = DERIV / "tabular" / "baseline_T4"
OUT_DIR = DERIV / "verbose" / "longitudinal_allvisits_T4split"
SEEDS = (0, 1, 2)
SPLITS = ("train", "val", "test")
ID_COL = "Patient_ID"
CONVERTER_GROUPS = {"pMCI", "pCN_to_AD"}


def viscode_to_months(v) -> float:
    s = str(v).strip().lower()
    if s in ("bl", "baseline"):
        return 0.0
    if s in ("sc", "scr", "scmri", "screening"):
        return -1.0
    m = re.fullmatch(r"m(\d+)", s)
    return float(m.group(1)) if m else np.nan


def folds_from(split_dir: Path, seed: int) -> dict:
    fold = {}
    for sp in SPLITS:
        ids = pd.read_csv(split_dir / f"seed_{seed}" / f"{sp}.csv",
                          usecols=[ID_COL], low_memory=False)[ID_COL].astype(str)
        for pid in ids:
            fold[pid] = sp
    return fold


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="Dry run: report counts, write nothing.")
    args = ap.parse_args()

    print(f"Reading longitudinal master: {LONG_MASTER}")
    df = pd.read_csv(LONG_MASTER, low_memory=False)
    df[ID_COL] = df[ID_COL].astype(str)
    df["visit_months"] = df["VISCODE_2"].map(viscode_to_months)
    df = df[~df[ID_COL].apply(is_excluded_subject)].reset_index(drop=True)
    df["is_converter"] = df["conversion_group"].isin(CONVERTER_GROUPS)
    print(f"  {len(df)} rows | {df[ID_COL].nunique()} patients | "
          f"converters {df.loc[df['is_converter'], ID_COL].nunique()}")

    for seed in SEEDS:
        canon = folds_from(CANON_SPLIT_DIR, seed)
        t4 = folds_from(T4_SPLIT_DIR, seed)

        def assign(pid, is_conv):
            if is_conv and pid in t4:
                return t4[pid]
            return canon.get(pid, np.nan)

        df["_split"] = [assign(p, c) for p, c in zip(df[ID_COL], df["is_converter"])]
        unmapped = df["_split"].isna()
        n_unmapped = df.loc[unmapped, ID_COL].nunique()
        seed_df = df[~unmapped].copy()

        sets = {sp: set(seed_df.loc[seed_df["_split"] == sp, ID_COL]) for sp in SPLITS}
        overlap = (sets["train"] & sets["val"]) | (sets["train"] & sets["test"]) | (sets["val"] & sets["test"])
        assert not overlap, f"seed {seed}: patient overlap across folds: {list(overlap)[:5]}"

        # LEAKAGE CHECK: every baseline_T4 val/test converter must be OUT of encoder-train
        t4_heldout = {p for p, s in t4.items() if s in ("val", "test")}
        enc_train = sets["train"]
        leak = t4_heldout & enc_train
        assert not leak, f"seed {seed}: {len(leak)} baseline_T4 val/test converters leaked into encoder-train: {list(leak)[:5]}"

        conv_by_fold = {sp: len(seed_df.loc[(seed_df['_split'] == sp) & seed_df['is_converter'], ID_COL].unique())
                        for sp in SPLITS}
        print(f"\nseed {seed}: dropped {n_unmapped} unmapped patients | leakage check PASS "
              f"({len(t4_heldout)} T4 held-out converters, 0 in encoder-train)")
        for sp in SPLITS:
            part = seed_df[seed_df["_split"] == sp].drop(columns=["_split", "is_converter"])
            iw = part[part["visit_months"].between(0, 12)]
            print(f"  {sp:5s}: {len(part):4d} rows / {part[ID_COL].nunique():3d} pts "
                  f"| converters {conv_by_fold[sp]:3d} | in-window {len(iw):4d} "
                  f"{iw['VISCODE_2'].value_counts().to_dict()}")
            if not args.check:
                out = OUT_DIR / f"seed_{seed}" / f"{sp}.csv"
                out.parent.mkdir(parents=True, exist_ok=True)
                part.to_csv(out, index=False)

    df.drop(columns=["_split", "is_converter"], errors="ignore", inplace=True)
    print("\n[--check] no files written." if args.check else f"\nWrote splits -> {OUT_DIR}")


if __name__ == "__main__":
    main()
