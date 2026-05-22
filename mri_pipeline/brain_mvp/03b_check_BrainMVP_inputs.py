"""
03b_check_BrainMVP_inputs.py — validate the brainmvp_inputs/ tree
=================================================================
Post-hoc QA for the outputs of 03_prepare_BrainMVP.py. Confirms every
preprocessed volume matches the BrainMVP classification spec before the
fine-tuning sweep consumes it.

Expected per-volume spec (from 03_prepare_BrainMVP.py):
  - path : brainmvp_inputs/sub-{sub}/ses-{ses}/
           sub-{sub}_ses-{ses}_space-BrainMVP128x64_desc-preproc_T1w.nii.gz
  - shape: 128 x 128 x 64           (the model RandSpatialCrops to 96x96x64)
  - dtype: floating point
  - intensities: percentile-clipped [5,95] then rescaled to [0,1]; background 0
                 (small interpolation over/undershoot from the final Resize is
                 expected, so the tolerated range is ~[-0.05, 1.1])
  - orientation: RAS, nominal 1 mm isotropic
  - finite (no NaN / Inf), non-degenerate (real brain content)

Checks, per volume
------------------
  HARD (cause a non-zero exit): unreadable file; wrong shape; any NaN/Inf;
       all-zero volume; ~constant volume; max intensity > 1.5 (not normalised).
  WARN (reported, exit stays 0): min < -0.05 or max in (1.1, 1.5]; zooms not
       ~1 mm; orientation not RAS; non-float dtype; very low brain fraction.

It also cross-checks the file count against brainmvp_manifest.csv and lists any
manifest rows with status failed / missing_input.

Usage (on CSD3, env: mri)
-------------------------
  python mri_pipeline/brain_mvp/03b_check_BrainMVP_inputs.py
  python mri_pipeline/brain_mvp/03b_check_BrainMVP_inputs.py --sample 200
  python mri_pipeline/brain_mvp/03b_check_BrainMVP_inputs.py \\
      --inputs-dir /path/to/brainmvp_inputs --report out.csv

Exit code 0 = all volumes pass the HARD checks; 1 = at least one HARD failure.
"""

from __future__ import annotations

import argparse
import glob
import os
import random
import sys
from collections import Counter

import numpy as np
import nibabel as nib
import pandas as pd

EXPECT_SHAPE = (128, 128, 64)
FNAME_GLOB = "*_space-BrainMVP128x64_desc-preproc_T1w.nii.gz"
DEFAULT_INPUTS = "/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/brainmvp_inputs"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inputs-dir", default=DEFAULT_INPUTS,
                   help="Root of the brainmvp_inputs tree.")
    p.add_argument("--sample", type=int, default=0,
                   help="Check a random N volumes instead of all (0 = all).")
    p.add_argument("--report", default=None,
                   help="CSV path for the per-volume report "
                        "(default: <inputs-dir>/brainmvp_check_report.csv).")
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def check_volume(path: str) -> dict:
    """Validate one volume. Returns a record with `hard` / `warn` issue lists."""
    rec = {"path": path, "shape": "", "dtype": "", "orient": "", "zooms": "",
           "min": "", "max": "", "brain_frac": "", "brain_mean": "",
           "brain_std": "", "hard": [], "warn": []}
    try:
        img = nib.load(path)
        data = np.asanyarray(img.dataobj, dtype=np.float64)
    except Exception as exc:                       # unreadable / corrupt
        rec["hard"].append(f"load error: {exc}")
        return rec

    shape = tuple(int(s) for s in data.shape)
    rec["shape"] = str(shape)
    if shape != EXPECT_SHAPE:
        rec["hard"].append(f"shape {shape} != {EXPECT_SHAPE}")

    rec["dtype"] = str(img.get_data_dtype())
    if not np.issubdtype(img.get_data_dtype(), np.floating):
        rec["warn"].append(f"dtype {rec['dtype']} is not floating point")

    rec["orient"] = "".join(nib.aff2axcodes(img.affine))
    if rec["orient"] != "RAS":
        rec["warn"].append(f"orientation {rec['orient']} != RAS")

    zooms = tuple(round(float(z), 3) for z in img.header.get_zooms()[:3])
    rec["zooms"] = str(zooms)
    if any(abs(z - 1.0) > 0.05 for z in zooms):
        rec["warn"].append(f"zooms {zooms} not ~1 mm")

    n_nan = int(np.isnan(data).sum())
    n_inf = int(np.isinf(data).sum())
    if n_nan or n_inf:
        rec["hard"].append(f"{n_nan} NaN, {n_inf} Inf voxels")
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        rec["hard"].append("no finite voxels")
        return rec

    vmin, vmax = float(finite.min()), float(finite.max())
    rec["min"], rec["max"] = round(vmin, 4), round(vmax, 4)
    if vmax <= 1e-6:
        rec["hard"].append("volume is all-zero")
        return rec
    if vmax > 1.5:
        rec["hard"].append(f"max {vmax:.2f} > 1.5 — not [0,1]-normalised")
    elif vmax > 1.1:
        rec["warn"].append(f"max {vmax:.3f} slightly above 1.1")
    if vmin < -0.1:
        rec["hard"].append(f"min {vmin:.2f} < -0.1 — not [0,1]-normalised")
    elif vmin < -0.05:
        rec["warn"].append(f"min {vmin:.3f} slightly below -0.05")

    brain = finite[finite > 0.05]
    rec["brain_frac"] = round(brain.size / finite.size, 4)
    if rec["brain_frac"] < 0.02:
        rec["warn"].append(f"brain fraction {rec['brain_frac']} very low")
    if brain.size:
        rec["brain_mean"] = round(float(brain.mean()), 4)
        rec["brain_std"] = round(float(brain.std()), 4)
        if rec["brain_std"] < 1e-3:
            rec["hard"].append("brain std ~0 — constant/degenerate volume")
    return rec


def _pct(values, q):
    return round(float(np.percentile(values, q)), 4) if values else float("nan")


def main():
    args = parse_args()
    root = args.inputs_dir
    print("=" * 70)
    print("  03b_check_BrainMVP_inputs — validating brainmvp_inputs/")
    print(f"  Inputs dir : {root}")
    print("=" * 70)

    if not os.path.isdir(root):
        print(f"[ERROR] inputs dir not found: {root}")
        sys.exit(1)

    files = sorted(glob.glob(os.path.join(root, "sub-*", "ses-*", FNAME_GLOB)))
    print(f"  Volumes found on disk : {len(files)}")
    if not files:
        print("[ERROR] no preprocessed volumes found — has 03_prepare_BrainMVP run?")
        sys.exit(1)

    # ── Manifest cross-check ──────────────────────────────────────────────────
    manifest_path = os.path.join(root, "brainmvp_manifest.csv")
    if os.path.isfile(manifest_path):
        mdf = pd.read_csv(manifest_path, dtype=str)
        counts = mdf["status"].value_counts().to_dict() if "status" in mdf else {}
        n_ok = int(counts.get("ok", 0))
        print(f"  Manifest rows         : {len(mdf)}  status={counts}")
        if n_ok != len(files):
            print(f"  [WARN] manifest 'ok'={n_ok} != {len(files)} files on disk")
        bad = mdf[mdf["status"].isin(["failed", "missing_input"])] \
            if "status" in mdf else mdf.iloc[0:0]
        for _, row in bad.iterrows():
            print(f"  [manifest {row['status']}] "
                  f"sub-{row.get('bids_sub','?')}_ses-{row.get('bids_ses','?')}"
                  f"  {str(row.get('error',''))[:90]}")
    else:
        print(f"  [WARN] no brainmvp_manifest.csv at {manifest_path}")

    # ── Per-volume checks ─────────────────────────────────────────────────────
    to_check = files
    if args.sample and args.sample < len(files):
        random.seed(args.seed)
        to_check = sorted(random.sample(files, args.sample))
        print(f"  Sampling {len(to_check)} of {len(files)} volumes.")
    print(f"\n  Checking {len(to_check)} volumes ...")

    records = []
    for i, path in enumerate(to_check, 1):
        records.append(check_volume(path))
        if i % 200 == 0:
            print(f"    {i}/{len(to_check)} ...")

    n_hard = sum(1 for r in records if r["hard"])
    n_warn = sum(1 for r in records if r["warn"] and not r["hard"])
    n_pass = len(records) - n_hard - n_warn

    # ── Cohort summary stats ──────────────────────────────────────────────────
    shapes = Counter(r["shape"] for r in records)
    mins = [r["min"] for r in records if isinstance(r["min"], float)]
    maxs = [r["max"] for r in records if isinstance(r["max"], float)]
    bmean = [r["brain_mean"] for r in records if isinstance(r["brain_mean"], float)]
    bfrac = [r["brain_frac"] for r in records if isinstance(r["brain_frac"], float)]

    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    print(f"  pass: {n_pass}    warn-only: {n_warn}    HARD-FAIL: {n_hard}")
    print(f"  shapes              : {dict(shapes)}")
    if mins:
        print(f"  intensity min  [p1, p50, p99] : "
              f"[{_pct(mins,1)}, {_pct(mins,50)}, {_pct(mins,99)}]")
        print(f"  intensity max  [p1, p50, p99] : "
              f"[{_pct(maxs,1)}, {_pct(maxs,50)}, {_pct(maxs,99)}]")
        print(f"  brain mean     [p1, p50, p99] : "
              f"[{_pct(bmean,1)}, {_pct(bmean,50)}, {_pct(bmean,99)}]")
        print(f"  brain fraction [p1, p50, p99] : "
              f"[{_pct(bfrac,1)}, {_pct(bfrac,50)}, {_pct(bfrac,99)}]")

    # List offenders (hard first, then warn) — capped.
    flagged = [r for r in records if r["hard"]] + \
              [r for r in records if r["warn"] and not r["hard"]]
    if flagged:
        print(f"\n  Flagged volumes ({len(flagged)}; showing up to 30):")
        for r in flagged[:30]:
            tag = "HARD" if r["hard"] else "warn"
            issues = "; ".join(r["hard"] + r["warn"])
            print(f"    [{tag}] {os.path.relpath(r['path'], root)}  ->  {issues}")
        if len(flagged) > 30:
            print(f"    ... and {len(flagged) - 30} more (see the CSV report).")
    else:
        print("\n  All checked volumes pass every check.")

    # ── Per-volume CSV report ─────────────────────────────────────────────────
    report_path = args.report or os.path.join(root, "brainmvp_check_report.csv")
    rep_df = pd.DataFrame([{
        "path": r["path"], "shape": r["shape"], "dtype": r["dtype"],
        "orient": r["orient"], "zooms": r["zooms"], "min": r["min"],
        "max": r["max"], "brain_frac": r["brain_frac"],
        "brain_mean": r["brain_mean"], "brain_std": r["brain_std"],
        "hard_issues": " | ".join(r["hard"]),
        "warn_issues": " | ".join(r["warn"]),
    } for r in records])
    try:
        rep_df.to_csv(report_path, index=False)
        print(f"\n  Per-volume report written: {report_path}")
    except Exception as exc:
        print(f"\n  [WARN] could not write report CSV: {exc}")

    print("=" * 70)
    if n_hard:
        print(f"  RESULT: FAIL — {n_hard} volume(s) failed a HARD check.")
        sys.exit(1)
    print("  RESULT: PASS — all checked volumes match the BrainMVP spec.")
    sys.exit(0)


if __name__ == "__main__":
    main()
