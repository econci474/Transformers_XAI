"""
Script 03 — Prepare ViT Inputs from sMRIprep Outputs
=====================================================
Converts sMRIprep MNI152NLin2009cAsym outputs into the format required by
the ViT_recipe_for_AD model (Kunanbayev et al., MICCAI 2024).

Processing chain (per subject, using MONAI transforms):
  1. Load sMRIprep MNI-space T1w
  2. Ensure channel-first layout          → (1, H, W, D)
  3. Reorient to RAS                      → safety check (MNI is already RAS)
  4. Resample to 1.75 × 1.75 × 1.75 mm   → isotropic voxel spacing
  5. Z-score normalise (nonzero only)     → brain voxels only, background = 0
  6. Crop foreground                      → remove zero-padding around brain
  7. Resize to 128 × 128 × 128           → final output shape

Target specs (from ViT config_vitb_mae.yaml):
  spacing:              [1.75, 1.75, 1.75] mm
  resize:               [128, 128, 128]
  orientation:          RAS
  normalize_non_zero:   True

Input (sessionwise sMRIprep, e.g. CSD3):
    .../smriprep_sessionwise/smriprep/sub-{sub}/ses-{ses}/anat/
        sub-{sub}_ses-{ses}_space-MNI152NLin2009cAsym_res-1_desc-preproc_T1w.nii.gz

Input (subject-level fallback, e.g. local Windows non-sessionwise):
    .../smriprep/.../sub-{sub}/anat/
        sub-{sub}_space-MNI152NLin2009cAsym_res-1_desc-preproc_T1w.nii.gz
    (Used only when --session bl and no ses-bl/ subdirectory exists.)

Output (always BIDS-correct, includes ses-*):
    derivatives/vit_inputs/sub-{sub}/ses-{ses}/
        sub-{sub}_ses-{ses}_space-ViT128_desc-preproc_T1w.nii.gz

Manifest CSV (columns include bids_sub, bids_ses):
    derivatives/vit_inputs/vit_manifest.csv

Session selection (mutually exclusive):
    --session SESSION   Single session to process. Default: 'bl' (baseline).
                        e.g. 'bl', 'm24'.
    --long N            Longitudinal mode. 'all' = every session on disk;
                        N (positive int) = sessions from baseline up to ses-m{12*N}
                        inclusive (e.g. --long 3 = up to m36).

Usage:
    python 03_prepare_ViT.py                           # baseline (ses-bl) — default
    python 03_prepare_ViT.py --session bl --n-workers 4
    python 03_prepare_ViT.py --long 3 --n-workers 4    # baseline through m36
    python 03_prepare_ViT.py --long all --n-workers 4  # all available sessions
    python 03_prepare_ViT.py --dry-run
    python 03_prepare_ViT.py --overwrite
    python 03_prepare_ViT.py --sub 002S0413
    python 03_prepare_ViT.py --smriprep-dir /path/to/smriprep --out-root /path/to/vit_inputs
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

# MONAI transforms for neuroimaging preprocessing
import monai.transforms as mt

# ── Paths ──────────────────────────────────────────────────────────────────────
# LOCAL (Windows) paths — used when running on your local machine
PROJECT_ROOT   = r"D:\ADNI_BIDS_project"
SMRIPREP_DIR   = os.path.join(PROJECT_ROOT, "derivatives", "smriprep",
                              "hpc_first_lex_2026_04_20", "mni_staging")
OUT_ROOT       = os.path.join(PROJECT_ROOT, "derivatives", "vit_inputs")
MANIFEST_PATH  = os.path.join(OUT_ROOT, "vit_manifest.csv")
SNP_TSV        = os.path.join(PROJECT_ROOT, "bids", "genotype",
                              "subjects_with_snp_and_mri.tsv")

# HPC (Linux) paths — uncomment these and comment out the LOCAL block above
# when running this script directly on the Cambridge HPC via SSH.
# PROJECT_ROOT   = "/home/ec474/rds/hpc-work/ADNI_SMRIPREP"
# SMRIPREP_DIR   = os.path.join(PROJECT_ROOT, "derivatives", "smriprep", "smriprep")
# OUT_ROOT       = os.path.join(PROJECT_ROOT, "derivatives", "vit_inputs")
# MANIFEST_PATH  = os.path.join(OUT_ROOT, "vit_manifest.csv")
# SNP_TSV        = os.path.join(PROJECT_ROOT, "bids", "genotype",
#                               "subjects_with_snp_and_mri.tsv")

# ── ViT target dimensions ─────────────────────────────────────────────────────
TARGET_SPACING = (1.75, 1.75, 1.75)   # mm
TARGET_SHAPE   = (128, 128, 128)      # voxels

# ── Argument parsing ───────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Prepare ViT inputs from sMRIprep derivatives.")
parser.add_argument("--dry-run",   action="store_true",
                    help="Print plan without writing files.")
parser.add_argument("--overwrite", action="store_true",
                    help="Reprocess subjects that already have output.")
parser.add_argument("--n-workers", type=int, default=1,
                    help="Number of parallel worker processes (default: 1).")
parser.add_argument("--sub", type=str, default=None,
                    help="Process a single subject only, e.g. '002S0413'.")
parser.add_argument("--smriprep-dir", type=str, default=None,
                    help="Override SMRIPREP_DIR (root with sub-*/anat/ MNI T1w files). "
                         "Lets the same script run on local Windows and on HPC without editing source.")
parser.add_argument("--out-root", type=str, default=None,
                    help="Override OUT_ROOT (where vit_inputs/sub-*/ are written).")
session_group = parser.add_mutually_exclusive_group()
session_group.add_argument("--session", type=str, default="bl",
                    help="Single session to process (e.g. 'bl' for baseline, 'm24' for month-24). "
                         "Default: 'bl'. If sMRIprep is sessionwise, looks under sub-*/ses-{label}/anat/. "
                         "Falls back to subject-level layout sub-*/anat/ when --session bl and no ses- dirs exist.")
session_group.add_argument("--long", type=str, default=None,
                    help="Longitudinal mode: 'all' for every available session per subject, "
                         "or an integer N for sessions from baseline up to ses-m{12*N} inclusive "
                         "(e.g. '3' = up to m36). Mutually exclusive with --session.")
args = parser.parse_args()

DRY_RUN   = args.dry_run
OVERWRITE = args.overwrite
N_WORKERS = args.n_workers

# Apply CLI path overrides (HPC-friendly)
if args.smriprep_dir:
    SMRIPREP_DIR = args.smriprep_dir
if args.out_root:
    OUT_ROOT = args.out_root
    MANIFEST_PATH = os.path.join(OUT_ROOT, "vit_manifest.csv")

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


def session_to_months(ses_label: str):
    """Convert 'bl' -> 0, 'm24' -> 24. Returns None for unsupported labels (date-format etc.)."""
    if ses_label == "bl":
        return 0
    m = re.match(r"^m(\d+)$", ses_label)
    if m:
        return int(m.group(1))
    return None


# ── Helper: find sMRIprep T1w ─────────────────────────────────────────────────
def find_smriprep_t1w(sub: str, ses: str):
    """
    Returns the path to the sMRIprep MNI-space preproc T1w for a (subject, session) pair.
    First tries the sessionwise layout sub-XYZ/ses-{ses}/anat/...; if not found and
    ses == 'bl', falls back to the subject-level layout sub-XYZ/anat/... (so the local
    Windows non-sessionwise sMRIprep tree still works for baseline).
    """
    sub_label = f"sub-{sub}"
    ses_label = f"ses-{ses}"

    # Sessionwise layout (HPC)
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

    # Subject-level fallback (local Windows non-sessionwise) — only valid for baseline
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
    """
    Full preprocessing chain for one (subject, session) pair using MONAI transforms.
    Steps: Orientation → Spacing → NormalizeIntensity → CropForeground → Resize
    Output: OUT_ROOT/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_space-ViT128_desc-preproc_T1w.nii.gz
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
        "output_shape":   str(TARGET_SHAPE),
        "output_spacing": str(TARGET_SPACING),
        "brain_voxels":   "",
        "znorm_mean":     "",
        "znorm_std":      "",
        "status":         "unknown",
        "error":          "",
        "timestamp":      datetime.datetime.now().isoformat(),
    }

    out_dir  = os.path.join(OUT_ROOT, sub_label, ses_label)
    out_name = f"{sub_label}_{ses_label}_space-ViT128_desc-preproc_T1w.nii.gz"
    out_path = os.path.join(out_dir, out_name)
    record["output_path"] = out_path

    # ── Check if already done ────────────────────────────────────────────────
    if os.path.isfile(out_path) and not overwrite:
        record["status"] = "skipped_exists"
        return record

    # ── Locate sMRIprep MNI T1w ──────────────────────────────────────────────
    t1w_path = find_smriprep_t1w(sub, ses)
    record["t1w_path"] = t1w_path or ""

    if not t1w_path:
        record["status"] = "missing_input"
        record["error"]  = (f"sMRIprep MNI preproc T1w not found for {sub_label} {ses_label}")
        return record

    # Record input properties
    try:
        src_img = nib.load(t1w_path)
        record["input_shape"]   = str(src_img.shape)
        record["input_spacing"] = str(tuple(
            round(float(z), 4) for z in src_img.header.get_zooms()[:3]))
        record["input_orient"]  = str(nib.aff2axcodes(src_img.affine))
    except Exception:
        pass

    if dry_run:
        record["status"] = "dry_run"
        return record

    try:
        # ── MONAI transform pipeline ─────────────────────────────────────────
        # This matches the ViT_recipe_for_AD preprocessing exactly:
        #   Orientation → Spacing → NormalizeIntensity → CropForeground → Resize
        transform = mt.Compose([
            mt.LoadImaged(keys=["image"]),
            mt.EnsureChannelFirstd(keys=["image"]),
            mt.Orientationd(keys=["image"], axcodes="RAS"),
            mt.Spacingd(keys=["image"], pixdim=TARGET_SPACING),
            mt.NormalizeIntensityd(keys=["image"], nonzero=True),
            mt.CropForegroundd(keys=["image"], source_key="image"),
            mt.Resized(keys=["image"], spatial_size=TARGET_SHAPE),
        ])

        data_dict = transform({"image": t1w_path})
        volume = data_dict["image"]  # shape: (1, 128, 128, 128)

        # ── Record stats ─────────────────────────────────────────────────────
        vol_np = volume.numpy().squeeze()  # (128, 128, 128)
        brain_mask = np.abs(vol_np) > 1e-6
        n_brain = int(brain_mask.sum())
        record["brain_voxels"] = n_brain

        if n_brain > 0:
            brain_vals = vol_np[brain_mask]
            record["znorm_mean"] = round(float(brain_vals.mean()), 6)
            record["znorm_std"]  = round(float(brain_vals.std()), 6)
        record["output_shape"] = str(vol_np.shape)

        # ── Save as NIfTI ────────────────────────────────────────────────────
        os.makedirs(out_dir, exist_ok=True)

        # Build affine for 128³ @ 1.75mm in RAS
        affine = np.eye(4)
        for i in range(3):
            affine[i, i] = TARGET_SPACING[i]
        # Centre the volume (approximate MNI origin)
        for i in range(3):
            affine[i, 3] = -TARGET_SPACING[i] * TARGET_SHAPE[i] / 2.0

        out_img = nib.Nifti1Image(vol_np, affine=affine)
        nib.save(out_img, out_path)

        record["status"] = "ok"

    except Exception:
        record["status"] = "failed"
        record["error"]  = traceback.format_exc(limit=3)

    return record


# ── Discover subject/session pairs ────────────────────────────────────────────
def discover_pairs():
    """
    Walk the sMRIprep tree to find all (subject, session) pairs with MNI T1w
    that match the requested mode (single-session or longitudinal).
    Returns a sorted list of (sub, ses) tuples.
    """
    pairs = []

    # Load SNP-matched subject list if available
    snp_subjects = None
    if os.path.isfile(SNP_TSV):
        tsv = pd.read_csv(SNP_TSV, sep='\t')
        snp_subjects = set(
            tsv['participant_id'].str.replace('sub-', '', regex=False).tolist()
        )
        print(f"  SNP-matched subjects loaded: {len(snp_subjects)}")
    else:
        print(f"  WARNING: SNP TSV not found at {SNP_TSV}; "
              "processing ALL subjects.")

    if not os.path.isdir(SMRIPREP_DIR):
        print(f"ERROR: sMRIprep directory not found: {SMRIPREP_DIR}")
        sys.exit(1)

    n_sessions_skipped = 0
    for sub_dir in sorted(os.listdir(SMRIPREP_DIR)):
        # Skip sMRIprep's per-subject and per-(subject, session) HTML reports
        # (e.g. 'sub-002S0413.html', 'sub-002S0413_ses-m108.html') which also
        # start with 'sub-' but are files, not directories.
        if not sub_dir.startswith("sub-") or "_" in sub_dir:
            continue
        sub_root = os.path.join(SMRIPREP_DIR, sub_dir)
        if not os.path.isdir(sub_root):
            continue
        sub = sub_dir[4:]

        # Filter to single subject if --sub was given
        if args.sub and sub != args.sub:
            continue

        # Filter to SNP-matched subjects
        if snp_subjects and sub not in snp_subjects:
            continue

        # ── Single-session mode (default: bl) ─────────────────────────────────
        if LONG_MODE is None:
            ses = SINGLE_SESSION
            if find_smriprep_t1w(sub, ses):
                pairs.append((sub, ses))
            continue

        # ── Longitudinal mode: enumerate ses-* under each subject ────────────
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

    if LONG_MODE == "cutoff":
        print(f"  Sessions skipped (beyond m{MAX_MONTHS} cutoff or unparseable): "
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
    print("Script 03 — Prepare ViT inputs from sMRIprep derivatives")
    print("=" * 65)
    print(f"  sMRIprep dir  : {SMRIPREP_DIR}")
    print(f"  Output root   : {OUT_ROOT}")
    print(f"  Mode          : {mode_str}")
    print(f"  Target spacing: {TARGET_SPACING} mm")
    print(f"  Target shape  : {TARGET_SHAPE}")
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

    # ── Load existing manifest to allow resume ────────────────────────────────
    existing_records = {}
    if os.path.isfile(MANIFEST_PATH):
        old_df = pd.read_csv(MANIFEST_PATH, dtype=str)
        for _, row in old_df.iterrows():
            key = (row.get("bids_sub", ""), row.get("bids_ses", ""))
            existing_records[key] = row.to_dict()
        print(f"  Resuming: {len(existing_records)} (sub, ses) records already in manifest.")

    # ── Process pairs ─────────────────────────────────────────────────────────
    all_records = dict(existing_records)
    n_ok = n_skip = n_fail = n_missing = n_dry = 0

    def handle_record(rec):
        nonlocal n_ok, n_skip, n_fail, n_missing, n_dry
        key = (rec["bids_sub"], rec.get("bids_ses", ""))
        all_records[key] = rec
        s = rec["status"]
        if s == "ok":              n_ok      += 1
        elif s.startswith("sk"):   n_skip    += 1
        elif s == "missing_input": n_missing += 1
        elif s == "dry_run":       n_dry     += 1
        else:                      n_fail    += 1

        label = f"sub-{rec['bids_sub']}_ses-{rec.get('bids_ses', '')}"
        if s == "ok":
            print(f"  [OK]      {label}  brain={rec['brain_voxels']}  "
                  f"shape={rec['output_shape']}")
        elif s.startswith("sk"):
            print(f"  [SKIP]    {label}")
        elif s == "missing_input":
            print(f"  [MISSING] {label}")
        elif s == "dry_run":
            print(f"  [DRY-RUN] {label}  ←  {rec['t1w_path']}")
        else:
            print(f"  [FAIL]    {label}  !  {rec['error'][:120]}")

    if N_WORKERS > 1 and not DRY_RUN:
        print(f"Running with {N_WORKERS} parallel workers...")
        with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
            futures = {
                pool.submit(process_subject, sub, ses, DRY_RUN, OVERWRITE): (sub, ses)
                for sub, ses in pairs
            }
            for future in as_completed(futures):
                try:
                    rec = future.result()
                except Exception as exc:
                    sub, ses = futures[future]
                    rec = {
                        "bids_sub": sub, "bids_ses": ses, "status": "failed",
                        "error": str(exc), "t1w_path": "",
                        "output_path": "", "input_shape": "",
                        "input_spacing": "", "input_orient": "",
                        "output_shape": str(TARGET_SHAPE),
                        "output_spacing": str(TARGET_SPACING),
                        "brain_voxels": "", "znorm_mean": "",
                        "znorm_std": "",
                        "timestamp": datetime.datetime.now().isoformat(),
                    }
                handle_record(rec)
    else:
        for i, (sub, ses) in enumerate(pairs):
            print(f"[{i+1}/{len(pairs)}] sub-{sub}_ses-{ses} ...", end=" ",
                  flush=True)
            rec = process_subject(sub, ses, DRY_RUN, OVERWRITE)
            handle_record(rec)

    # ── Save manifest ─────────────────────────────────────────────────────────
    manifest_df = pd.DataFrame(list(all_records.values()))
    manifest_cols = [
        "bids_sub", "bids_ses", "status",
        "t1w_path", "output_path",
        "input_shape", "input_spacing", "input_orient",
        "output_shape", "output_spacing",
        "brain_voxels", "znorm_mean", "znorm_std",
        "error", "timestamp",
    ]
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
