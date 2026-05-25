"""
Script 01b -- BrainDINO preprocessing unit test / checker
=========================================================
Sanity-checks the volumes written by 01_prepare_braindino_inputs.py
against the paper-fidelity recipe (arXiv:2604.27277 Sections 4.1.2 + 4.3).

Designed to FAIL LOUDLY if any one assertion is violated so the user
catches recipe drift before training starts.

Assertions (per sampled volume):

  1. Shape: output volume is exactly (256, 256, 128) in RAS axis order
     (R, A, S=axial). Equivalent to paper's "128 axial x 256 coronal
     x 256 sagittal".
  2. Background preserved: voxels that are 0 in the volume sit OUTSIDE the
     foreground bbox. We can't reload the mask, so we check that the
     foreground fraction is in a plausible range (see assertion 5).
  3. Foreground stats: |mean_nonzero| < 0.05 and |std_nonzero - 1.0| <
     0.05 -- z-score normalisation worked end-to-end.
  4. Percentile clip applied: |min_nonzero| < 4.5 and max_nonzero < 4.5
     (after z-score, the 0.5th-99.5th percentile of a normal distribution
     sits at ~+/- 2.58 sigma; clipped pre-z-score, post-z-score extremes
     shouldn't exceed ~4 sigma; > 4.5 sigma indicates the clip didn't fire).
  5. Brain coverage: non-zero voxel fraction in [0.10, 0.60] of the total
     (8,388,608 voxels). Healthy brain extraction sits ~0.30-0.45.
  6. Axial-slice count is exactly 128 (paper Section 4.3 explicit
     requirement; redundant with assertion 1 but kept as a separate
     named check so the failure message is unambiguous).
  7. No NaN / inf: np.isfinite(volume).all().

Usage:
    python mri_pipeline/brain_dino/01b_check_braindino_inputs.py
    python mri_pipeline/brain_dino/01b_check_braindino_inputs.py --sample 30
    python mri_pipeline/brain_dino/01b_check_braindino_inputs.py --all
    python mri_pipeline/brain_dino/01b_check_braindino_inputs.py \\
        --manifest /path/to/braindino_manifest.csv

Exit code:
    0 -- every sampled volume passes every assertion
    1 -- at least one assertion failed (failing rows printed)
"""

from __future__ import annotations

import argparse
import os
import sys

import nibabel as nib
import numpy as np
import pandas as pd

# Auto-detect local vs HPC manifest location. First existing path wins;
# override with --manifest at runtime.
_MANIFEST_CANDIDATES = [
    "/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/braindino_inputs/braindino_manifest.csv",
    r"D:\ADNI_BIDS_project\derivatives\braindino_inputs\braindino_manifest.csv",
]
DEFAULT_MANIFEST = next(
    (p for p in _MANIFEST_CANDIDATES if os.path.isfile(p)),
    _MANIFEST_CANDIDATES[0],
)

EXPECTED_SHAPE   = (256, 256, 128)            # (R, A, S=axial)
TOTAL_VOXELS     = int(np.prod(EXPECTED_SHAPE))  # 8,388,608
MEAN_TOL         = 0.05                        # |mean_nonzero| <  0.05
STD_TOL          = 0.05                        # |std_nonzero - 1.0| < 0.05
MAX_ABS_Z        = 4.5                         # post-z-score |extremes| < 4.5
BRAIN_FRAC_RANGE = (0.10, 0.60)                # non-zero fraction must sit here
AXIAL_AXIS       = 2                           # in RAS, axis 2 = S = axial


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifest", default=DEFAULT_MANIFEST,
                   help=f"BrainDINO manifest CSV (default: {DEFAULT_MANIFEST}).")
    p.add_argument("--sample", type=int, default=20,
                   help="Number of volumes to sample (default: 20). Use --all "
                        "to check every row in the manifest.")
    p.add_argument("--all", action="store_true",
                   help="Check every volume in the manifest.")
    p.add_argument("--seed", type=int, default=0,
                   help="RNG seed for sampling reproducibility.")
    return p.parse_args()


def check_volume(path: str) -> dict:
    """Run all assertions on one volume. Returns a dict with
    {pass, fails: [str], stats: dict}."""
    result = {"path": path, "pass": True, "fails": [], "stats": {}}
    if not os.path.isfile(path):
        result["pass"] = False
        result["fails"].append(f"file does not exist: {path}")
        return result

    img = nib.load(path)
    vol = img.get_fdata().astype(np.float32)
    result["stats"]["shape"] = vol.shape

    # 1 + 6: shape + axial-slice count
    if vol.shape != EXPECTED_SHAPE:
        result["pass"] = False
        result["fails"].append(
            f"[1] shape {vol.shape} != expected {EXPECTED_SHAPE}")
    if vol.shape[AXIAL_AXIS] != 128:
        result["pass"] = False
        result["fails"].append(
            f"[6] axial dim {vol.shape[AXIAL_AXIS]} != 128")

    # 7: no NaN/inf
    if not np.isfinite(vol).all():
        result["pass"] = False
        result["fails"].append(
            f"[7] {int((~np.isfinite(vol)).sum())} non-finite voxels")

    brain = vol != 0
    n_brain = int(brain.sum())
    result["stats"]["n_brain"] = n_brain
    result["stats"]["brain_frac"] = n_brain / TOTAL_VOXELS

    # 5: brain coverage in plausible range
    lo, hi = BRAIN_FRAC_RANGE
    if not (lo <= result["stats"]["brain_frac"] <= hi):
        result["pass"] = False
        result["fails"].append(
            f"[5] brain fraction {result['stats']['brain_frac']:.3f} "
            f"outside [{lo:.2f}, {hi:.2f}] -- likely botched mask")

    if n_brain == 0:
        result["pass"] = False
        result["fails"].append("[3/4] empty brain -- skipping intensity checks")
        return result

    brain_vals = vol[brain]
    mean_nz = float(brain_vals.mean())
    std_nz  = float(brain_vals.std())
    min_nz  = float(brain_vals.min())
    max_nz  = float(brain_vals.max())
    result["stats"].update({
        "mean_nonzero": round(mean_nz, 5),
        "std_nonzero":  round(std_nz, 5),
        "min_nonzero":  round(min_nz, 4),
        "max_nonzero":  round(max_nz, 4),
    })

    # 3: foreground z-score stats
    if abs(mean_nz) >= MEAN_TOL:
        result["pass"] = False
        result["fails"].append(
            f"[3] |mean_nonzero| = {abs(mean_nz):.4f} >= {MEAN_TOL}")
    if abs(std_nz - 1.0) >= STD_TOL:
        result["pass"] = False
        result["fails"].append(
            f"[3] |std_nonzero - 1| = {abs(std_nz - 1.0):.4f} >= {STD_TOL}")

    # 4: percentile clip applied
    if abs(min_nz) > MAX_ABS_Z:
        result["pass"] = False
        result["fails"].append(
            f"[4] min_nonzero = {min_nz:.3f} > +/-{MAX_ABS_Z} "
            f"-- percentile clip didn't fire")
    if max_nz > MAX_ABS_Z:
        result["pass"] = False
        result["fails"].append(
            f"[4] max_nonzero = {max_nz:.3f} > {MAX_ABS_Z} "
            f"-- percentile clip didn't fire")

    return result


def main():
    args = parse_args()

    if not os.path.isfile(args.manifest):
        print(f"[ERROR] manifest not found: {args.manifest}")
        sys.exit(1)

    df = pd.read_csv(args.manifest, dtype=str)
    df = df[df["status"] == "ok"].reset_index(drop=True)
    if df.empty:
        print(f"[ERROR] manifest has 0 rows with status=='ok' (of "
              f"{len(pd.read_csv(args.manifest, dtype=str))} total)")
        sys.exit(1)

    if args.all or args.sample >= len(df):
        sample = df.copy()
    else:
        sample = df.sample(n=args.sample, random_state=args.seed).reset_index(drop=True)

    print(f"BrainDINO preprocessing checker")
    print(f"  Manifest        : {args.manifest}")
    print(f"  Manifest rows OK: {len(df)}")
    print(f"  Checking        : {len(sample)} volume(s)")
    print(f"  Expected shape  : {EXPECTED_SHAPE} (R, A, S=axial)")
    print(f"  Tolerances      : |mean| < {MEAN_TOL}, "
          f"|std - 1| < {STD_TOL}, |z| < {MAX_ABS_Z}")
    print(f"  Brain fraction  : {BRAIN_FRAC_RANGE}")
    print("-" * 70)

    n_pass = 0
    failures = []
    for _, r in sample.iterrows():
        chk = check_volume(r["output_path"])
        s = chk["stats"]
        if chk["pass"]:
            n_pass += 1
            print(f"  [PASS]  sub-{r['bids_sub']:>10s} ses-{r['bids_ses']:>4s}  "
                  f"shape={s.get('shape')}  brain={s.get('brain_frac', 0):.3f}  "
                  f"mean={s.get('mean_nonzero')}  std={s.get('std_nonzero')}")
        else:
            failures.append((r, chk))
            print(f"  [FAIL]  sub-{r['bids_sub']:>10s} ses-{r['bids_ses']:>4s}  "
                  f"{len(chk['fails'])} issue(s)")
            for f in chk["fails"]:
                print(f"            {f}")
            print(f"          stats: {s}")

    print("-" * 70)
    print(f"  PASSED: {n_pass} / {len(sample)}")
    if failures:
        print(f"  FAILED: {len(failures)}  (see above)")
        sys.exit(1)
    else:
        print("  ALL OK")
        sys.exit(0)


if __name__ == "__main__":
    main()
