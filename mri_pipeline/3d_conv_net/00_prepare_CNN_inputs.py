"""
00_prepare_CNN_inputs.py
========================
Prepare model-ready volumes for the Spasov-style 3D CNN baselines
(`3DCNN.py` / `3DSCNN.py`) from sMRIprep MNI outputs.

WHY A SEPARATE PREPROCESSING FROM THE ViT PIPELINE
--------------------------------------------------
`03_prepare_ViT.py` resamples to 1.75 mm and resizes to a 128^3 cube — the
grid the ViT_recipe_for_AD model expects. The Spasov 3D CNN is a plain
convolutional network: it does NOT need a cube and there is no reason to
downsample. So this script keeps the sMRIprep MNI volume at its **native
1 mm grid** and changes nothing geometric — it only z-score normalises the
intensities. Every subject is already in the same template space
(MNI152NLin2009cAsym, res-1), so every output has the identical shape and a
3D CNN can train on them directly.

Processing chain (per subject, MONAI transforms):
  1. Load sMRIprep MNI-space T1w
  2. Ensure channel-first layout            -> (1, H, W, D)
  3. Reorient to RAS                        -> safety check (MNI is already RAS)
  4. Z-score normalise on non-zero voxels   -> brain voxels only, background = 0

  NO resampling, NO cropping, NO resizing — the volume's geometry is left
  exactly as sMRIprep produced it.

Input / output shape
--------------------
sMRIprep `space-MNI152NLin2009cAsym_res-1` T1w volumes are **193 x 229 x 193**
at 1 mm isotropic, RAS — confirmed uniform across the cohort. The output keeps
that shape. (Spasov 2019 used the 182x218x182 FSL MNI152 grid; the 3D CNN
modules here use a LazyLinear head, so they adapt to the 193x229x193 grid.)

Input (sessionwise sMRIprep, e.g. CSD3):
    .../smriprep/sub-{sub}/ses-{ses}/anat/
        sub-{sub}_ses-{ses}_space-MNI152NLin2009cAsym_res-1_desc-preproc_T1w.nii.gz
Input (subject-level fallback, local non-sessionwise, baseline only):
    .../smriprep/.../sub-{sub}/anat/
        sub-{sub}_space-MNI152NLin2009cAsym_res-1_desc-preproc_T1w.nii.gz

Output (BIDS-correct, always includes ses-*):
    derivatives/cnn_inputs/sub-{sub}/ses-{ses}/
        sub-{sub}_ses-{ses}_space-CNN193_desc-preproc_T1w.nii.gz

Manifest CSV:
    derivatives/cnn_inputs/cnn_manifest.csv

Subject exclusions
------------------
Subjects flagged by `bidsification.exclusions.is_excluded_subject` (the
site-381 cohort and the corrupted-MRI subject) are skipped. The diagnostic-
reversion exclusion is NOT applied here — it is enforced downstream by the
post-exclusion train/val/test splits that `train_3dcnn.py` reads.

Session selection (mutually exclusive):
    --session SESSION   Single session (default 'bl').
    --long N            'all' = every session on disk; integer N = baseline
                        up to ses-m{12*N} inclusive.

Usage:
    python 00_prepare_CNN_inputs.py                       # baseline (ses-bl)
    python 00_prepare_CNN_inputs.py --long all --n-workers 4
    python 00_prepare_CNN_inputs.py --dry-run
    python 00_prepare_CNN_inputs.py --overwrite --sub 002S0413
    python 00_prepare_CNN_inputs.py --smriprep-dir /path --out-root /path
"""

import os
import re
import sys
import glob
import argparse
import datetime
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import nibabel as nib
import pandas as pd

import monai.transforms as mt

# Cross-pipeline subject-level exclusions (single source of truth in the repo).
# This file lives at mri_pipeline/3d_conv_net/ — the repo root is two up.
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(THIS_DIR))
sys.path.insert(0, REPO_ROOT)
from bidsification.exclusions import is_excluded_subject

# ── Paths ──────────────────────────────────────────────────────────────────────
# LOCAL (Windows) paths — used when running on the local machine.
PROJECT_ROOT   = r"D:\ADNI_BIDS_project"
SMRIPREP_DIR   = os.path.join(PROJECT_ROOT, "derivatives", "smriprep",
                              "hpc_first_lex_2026_04_20", "mni_staging")
OUT_ROOT       = os.path.join(PROJECT_ROOT, "derivatives", "cnn_inputs")
MANIFEST_PATH  = os.path.join(OUT_ROOT, "cnn_manifest.csv")

# HPC (Linux) paths — uncomment and comment out the LOCAL block above when
# running directly on the Cambridge HPC.
# PROJECT_ROOT  = "/home/ec474/rds/hpc-work/ADNI_SMRIPREP"
# SMRIPREP_DIR  = os.path.join(PROJECT_ROOT, "derivatives", "smriprep", "smriprep")
# OUT_ROOT      = os.path.join(PROJECT_ROOT, "derivatives", "cnn_inputs")
# MANIFEST_PATH = os.path.join(OUT_ROOT, "cnn_manifest.csv")

# Expected sMRIprep MNI res-1 grid — used only as a sanity check (a volume of a
# different shape is still processed, but flagged in the manifest).
EXPECT_SHAPE = (193, 229, 193)


def build_transform() -> mt.Compose:
    """The 3D-CNN preprocessing pipeline.

    Orientation -> z-score on non-zero voxels. NO Spacing/Crop/Resize: the
    sMRIprep MNI volume is kept at its native 1 mm grid; only the intensities
    are standardised (background stays 0). Exposed as a function so it can be
    unit-tested directly.
    """
    return mt.Compose([
        mt.LoadImaged(keys=["image"]),
        mt.EnsureChannelFirstd(keys=["image"]),
        mt.Orientationd(keys=["image"], axcodes="RAS"),
        mt.NormalizeIntensityd(keys=["image"], nonzero=True),
    ])


# ── BIDS <-> ADNI Patient_ID ──────────────────────────────────────────────────
def bids_sub_to_ptid(bids_sub: str):
    """'002S0413' -> '002_S_0413' (the form bidsification.exclusions expects).
    Returns None if the string does not look like an ADNI subject id."""
    m = re.match(r"^(\d+)S(\d+)$", bids_sub)
    return f"{m.group(1)}_S_{m.group(2)}" if m else None


def session_to_months(ses_label: str):
    """'bl' -> 0, 'm24' -> 24. None for unsupported labels."""
    if ses_label == "bl":
        return 0
    m = re.match(r"^m(\d+)$", ses_label)
    return int(m.group(1)) if m else None


# ── Argument parsing ───────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Prepare 3D-CNN inputs from sMRIprep derivatives.")
parser.add_argument("--dry-run",   action="store_true",
                    help="Print plan without writing files.")
parser.add_argument("--overwrite", action="store_true",
                    help="Reprocess subjects that already have output.")
parser.add_argument("--n-workers", type=int, default=1,
                    help="Number of parallel worker processes (default: 1).")
parser.add_argument("--sub", type=str, default=None,
                    help="Process a single subject only, e.g. '002S0413'.")
parser.add_argument("--smriprep-dir", type=str, default=None,
                    help="Override SMRIPREP_DIR (root with sub-*/.../anat/ MNI T1w).")
parser.add_argument("--out-root", type=str, default=None,
                    help="Override OUT_ROOT (where cnn_inputs/sub-*/ are written).")
session_group = parser.add_mutually_exclusive_group()
session_group.add_argument("--session", type=str, default="bl",
                    help="Single session to process (default 'bl').")
session_group.add_argument("--long", type=str, default=None,
                    help="Longitudinal mode: 'all' for every session, or integer "
                         "N for baseline up to ses-m{12*N} inclusive.")
args = parser.parse_args()

DRY_RUN   = args.dry_run
OVERWRITE = args.overwrite
N_WORKERS = args.n_workers

if args.smriprep_dir:
    SMRIPREP_DIR = args.smriprep_dir
if args.out_root:
    OUT_ROOT = args.out_root
    MANIFEST_PATH = os.path.join(OUT_ROOT, "cnn_manifest.csv")

# Resolve session selection mode
if args.long is None:
    LONG_MODE = None
    SINGLE_SESSION = args.session
    MAX_MONTHS = None
else:
    SINGLE_SESSION = None
    if args.long.lower() == "all" or args.long == "0":
        LONG_MODE = "all"
        MAX_MONTHS = None
    else:
        try:
            n_years = int(args.long)
            if n_years <= 0:
                raise ValueError
        except ValueError:
            parser.error(f"--long must be 'all' or a positive integer (got {args.long!r})")
        LONG_MODE = "cutoff"
        MAX_MONTHS = n_years * 12


# ── Helper: find sMRIprep T1w ─────────────────────────────────────────────────
def find_smriprep_t1w(sub: str, ses: str):
    """Path to the sMRIprep MNI-space preproc T1w for a (subject, session) pair.
    Tries the sessionwise layout sub-XYZ/ses-{ses}/anat/...; for ses == 'bl'
    falls back to the subject-level layout sub-XYZ/anat/... ."""
    sub_label = f"sub-{sub}"
    ses_label = f"ses-{ses}"

    ses_anat_dir = os.path.join(SMRIPREP_DIR, sub_label, ses_label, "anat")
    expected = os.path.join(
        ses_anat_dir,
        f"{sub_label}_{ses_label}_space-MNI152NLin2009cAsym_res-1_desc-preproc_T1w.nii.gz",
    )
    if os.path.isfile(expected):
        return expected
    hits = glob.glob(os.path.join(ses_anat_dir, "*MNI152*desc-preproc_T1w.nii.gz"))
    if hits:
        return hits[0]

    if ses == "bl":
        anat_dir = os.path.join(SMRIPREP_DIR, sub_label, "anat")
        expected = os.path.join(
            anat_dir,
            f"{sub_label}_space-MNI152NLin2009cAsym_res-1_desc-preproc_T1w.nii.gz",
        )
        if os.path.isfile(expected):
            return expected
        hits = glob.glob(os.path.join(anat_dir, "*MNI152*desc-preproc_T1w.nii.gz"))
        if hits:
            return hits[0]

    return None


# ── Core processing function ──────────────────────────────────────────────────
def process_subject(sub: str, ses: str, dry_run: bool = False,
                    overwrite: bool = False) -> dict:
    """Preprocess one (subject, session) pair: z-score the sMRIprep MNI T1w and
    save it (geometry unchanged). Output:
    OUT_ROOT/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_space-CNN193_desc-preproc_T1w.nii.gz
    """
    sub_label = f"sub-{sub}"
    ses_label = f"ses-{ses}"

    record = {
        "bids_sub":       sub,
        "bids_ses":       ses,
        "t1w_path":       "",
        "output_path":    "",
        "input_shape":    "",
        "input_spacing":  "",
        "input_orient":   "",
        "output_shape":   "",
        "brain_voxels":   "",
        "znorm_mean":     "",
        "znorm_std":      "",
        "status":         "unknown",
        "error":          "",
        "timestamp":      datetime.datetime.now().isoformat(),
    }

    out_dir  = os.path.join(OUT_ROOT, sub_label, ses_label)
    out_name = f"{sub_label}_{ses_label}_space-CNN193_desc-preproc_T1w.nii.gz"
    out_path = os.path.join(out_dir, out_name)
    record["output_path"] = out_path

    if os.path.isfile(out_path) and not overwrite:
        record["status"] = "skipped_exists"
        return record

    t1w_path = find_smriprep_t1w(sub, ses)
    record["t1w_path"] = t1w_path or ""
    if not t1w_path:
        record["status"] = "missing_input"
        record["error"]  = f"sMRIprep MNI preproc T1w not found for {sub_label} {ses_label}"
        return record

    # Record input properties
    try:
        src_img = nib.load(t1w_path)
        record["input_shape"]   = str(tuple(src_img.shape))
        record["input_spacing"] = str(tuple(
            round(float(z), 4) for z in src_img.header.get_zooms()[:3]))
        record["input_orient"]  = str(nib.aff2axcodes(src_img.affine))
    except Exception:
        src_img = None

    if dry_run:
        record["status"] = "dry_run"
        return record

    try:
        # ── Preprocess (see build_transform() above) ─────────────────────────
        volume = build_transform()({"image": t1w_path})["image"]
        vol_np = volume.numpy().squeeze()  # (193, 229, 193)
        record["output_shape"] = str(vol_np.shape)

        # Brain stats on the standardised volume.
        brain_mask = np.abs(vol_np) > 1e-6
        n_brain = int(brain_mask.sum())
        record["brain_voxels"] = n_brain
        if n_brain > 0:
            bv = vol_np[brain_mask]
            record["znorm_mean"] = round(float(bv.mean()), 6)
            record["znorm_std"]  = round(float(bv.std()), 6)

        # ── Save ─────────────────────────────────────────────────────────────
        # No geometric transform was applied, so the volume's affine is exactly
        # the source affine — preserve it (do NOT synthesise a new one).
        os.makedirs(out_dir, exist_ok=True)
        affine = src_img.affine if src_img is not None else np.eye(4)
        nib.save(nib.Nifti1Image(vol_np.astype(np.float32), affine=affine),
                 out_path)

        if tuple(vol_np.shape) != EXPECT_SHAPE:
            record["status"] = "ok_unexpected_shape"
            record["error"]  = (f"output shape {vol_np.shape} != expected "
                                 f"{EXPECT_SHAPE}")
        else:
            record["status"] = "ok"

    except Exception:
        record["status"] = "failed"
        record["error"]  = traceback.format_exc(limit=3)

    return record


# ── Discover subject/session pairs ────────────────────────────────────────────
def discover_pairs():
    """Walk the sMRIprep tree for (subject, session) pairs with MNI T1w that
    match the requested mode. Excluded subjects (site-381, corrupted MRI) are
    dropped here via bidsification.exclusions.is_excluded_subject."""
    pairs = []
    n_excluded = 0

    if not os.path.isdir(SMRIPREP_DIR):
        print(f"ERROR: sMRIprep directory not found: {SMRIPREP_DIR}")
        sys.exit(1)

    n_sessions_skipped = 0
    for sub_dir in sorted(os.listdir(SMRIPREP_DIR)):
        # Skip sMRIprep HTML reports (sub-*.html, sub-*_ses-*.html) — they also
        # start with 'sub-' but are files, not directories.
        if not sub_dir.startswith("sub-") or "_" in sub_dir:
            continue
        sub_root = os.path.join(SMRIPREP_DIR, sub_dir)
        if not os.path.isdir(sub_root):
            continue
        sub = sub_dir[4:]

        if args.sub and sub != args.sub:
            continue

        # Subject-level exclusions (site-381 cohort, corrupted MRI).
        ptid = bids_sub_to_ptid(sub)
        if ptid is not None and is_excluded_subject(ptid):
            n_excluded += 1
            continue

        if LONG_MODE is None:
            ses = SINGLE_SESSION
            if find_smriprep_t1w(sub, ses):
                pairs.append((sub, ses))
            continue

        ses_dirs = sorted(d for d in os.listdir(sub_root)
                          if d.startswith("ses-")
                          and os.path.isdir(os.path.join(sub_root, d)))
        for ses_dir in ses_dirs:
            ses = ses_dir[4:]
            if LONG_MODE == "cutoff":
                months = session_to_months(ses)
                if months is None or months > MAX_MONTHS:
                    n_sessions_skipped += 1
                    continue
            if find_smriprep_t1w(sub, ses):
                pairs.append((sub, ses))

    print(f"  Subjects skipped (exclusions): {n_excluded}")
    if LONG_MODE == "cutoff":
        print(f"  Sessions skipped (beyond m{MAX_MONTHS} cutoff / unparseable): "
              f"{n_sessions_skipped}")
    return pairs


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    if LONG_MODE is None:
        mode_str = f"single-session: ses-{SINGLE_SESSION}"
    elif LONG_MODE == "all":
        mode_str = "longitudinal: all sessions"
    else:
        mode_str = f"longitudinal: up to ses-m{MAX_MONTHS}"

    print("=" * 65)
    print("Script 00 — Prepare 3D-CNN inputs from sMRIprep derivatives")
    print("=" * 65)
    print(f"  sMRIprep dir  : {SMRIPREP_DIR}")
    print(f"  Output root   : {OUT_ROOT}")
    print(f"  Mode          : {mode_str}")
    print(f"  Expected shape: {EXPECT_SHAPE} (native 1 mm MNI, z-scored)")
    print(f"  Dry run       : {DRY_RUN}")
    print(f"  Overwrite     : {OVERWRITE}")
    print(f"  Workers       : {N_WORKERS}")
    print()

    pairs = discover_pairs()
    n_subs = len({s for s, _ in pairs})
    print(f"Pairs to process: {len(pairs)} (across {n_subs} unique subjects)")
    print()

    if not pairs:
        print("No (subject, session) pairs found. Check SMRIPREP_DIR / --session / --long.")
        sys.exit(0)

    os.makedirs(OUT_ROOT, exist_ok=True)

    # Load existing manifest to allow resume.
    existing_records = {}
    if os.path.isfile(MANIFEST_PATH):
        old_df = pd.read_csv(MANIFEST_PATH, dtype=str)
        for _, row in old_df.iterrows():
            existing_records[(row.get("bids_sub", ""),
                              row.get("bids_ses", ""))] = row.to_dict()
        print(f"  Resuming: {len(existing_records)} records already in manifest.")

    all_records = dict(existing_records)
    n_ok = n_skip = n_fail = n_missing = n_dry = 0

    def handle_record(rec):
        nonlocal n_ok, n_skip, n_fail, n_missing, n_dry
        all_records[(rec["bids_sub"], rec.get("bids_ses", ""))] = rec
        s = rec["status"]
        if s.startswith("ok"):       n_ok      += 1
        elif s.startswith("sk"):     n_skip    += 1
        elif s == "missing_input":   n_missing += 1
        elif s == "dry_run":         n_dry     += 1
        else:                        n_fail    += 1

        label = f"sub-{rec['bids_sub']}_ses-{rec.get('bids_ses', '')}"
        if s.startswith("ok"):
            tag = "[OK]     " if s == "ok" else "[OK*]    "
            print(f"  {tag}{label}  brain={rec['brain_voxels']}  "
                  f"shape={rec['output_shape']}")
        elif s.startswith("sk"):
            print(f"  [SKIP]    {label}")
        elif s == "missing_input":
            print(f"  [MISSING] {label}")
        elif s == "dry_run":
            print(f"  [DRY-RUN] {label}  <-  {rec['t1w_path']}")
        else:
            print(f"  [FAIL]    {label}  !  {rec['error'][:120]}")

    if N_WORKERS > 1 and not DRY_RUN:
        print(f"Running with {N_WORKERS} parallel workers...")
        with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
            futures = {pool.submit(process_subject, sub, ses, DRY_RUN, OVERWRITE):
                       (sub, ses) for sub, ses in pairs}
            for future in as_completed(futures):
                try:
                    rec = future.result()
                except Exception as exc:
                    sub, ses = futures[future]
                    rec = {"bids_sub": sub, "bids_ses": ses, "status": "failed",
                           "error": str(exc), "t1w_path": "", "output_path": "",
                           "input_shape": "", "input_spacing": "",
                           "input_orient": "", "output_shape": "",
                           "brain_voxels": "", "znorm_mean": "", "znorm_std": "",
                           "timestamp": datetime.datetime.now().isoformat()}
                handle_record(rec)
    else:
        for i, (sub, ses) in enumerate(pairs):
            print(f"[{i+1}/{len(pairs)}] sub-{sub}_ses-{ses} ...", end=" ", flush=True)
            handle_record(process_subject(sub, ses, DRY_RUN, OVERWRITE))

    # ── Save manifest ─────────────────────────────────────────────────────────
    manifest_df = pd.DataFrame(list(all_records.values()))
    manifest_cols = ["bids_sub", "bids_ses", "status", "t1w_path", "output_path",
                     "input_shape", "input_spacing", "input_orient",
                     "output_shape", "brain_voxels", "znorm_mean", "znorm_std",
                     "error", "timestamp"]
    present = [c for c in manifest_cols if c in manifest_df.columns]
    extras  = [c for c in manifest_df.columns if c not in manifest_cols]
    manifest_df = manifest_df[present + extras]
    manifest_df.to_csv(MANIFEST_PATH, index=False)

    # ── Summary ───────────────────────────────────────────────────────────────
    print()
    print("=" * 65)
    print(f"Done.  (sub, ses) pairs processed this run: {len(pairs)}")
    print(f"  OK            : {n_ok}")
    print(f"  Dry-run       : {n_dry}")
    print(f"  Skipped       : {n_skip}")
    print(f"  Missing input : {n_missing}")
    print(f"  Failed        : {n_fail}")
    print(f"  Manifest      : {MANIFEST_PATH}")
    print("=" * 65)

    if n_fail > 0:
        failed = manifest_df[manifest_df["status"] == "failed"][
            ["bids_sub", "bids_ses", "error"]]
        print("\nFailed pairs:")
        print(failed.to_string(index=False))


if __name__ == "__main__":
    main()
