"""
03_prepare_BrainMVP.py — Prepare BrainMVP inputs from sMRIprep outputs
======================================================================
Converts sMRIprep MNI152NLin2009cAsym T1w scans into the format
specified by the BrainMVP paper (CVPR 2025 Highlight, Appendix I
"Classification"):

  1. Load sMRIprep MNI-space T1w
  2. Ensure channel-first layout                → (1, H, W, D)
  3. Reorient to RAS                            → (safety check)
  4. Resample to 1.0 × 1.0 × 1.0 mm            → isometric bilinear
  5. Skull stripping (already done by sMRIprep)  → no-op
  6. Percentile clip [1st, 99th]                → remove outlier voxels
  7. Rescale to [0, 1]                          → linear rescale
  8. Crop foreground                            → remove zero-padding
  9. Resize to 128 × 128 × 64                  → classification base size

Note: The paper further crops to 96 × 96 × 64 during training
(RandSpatialCrop) and inference (CenterSpatialCrop). That cropping
happens in the DataLoader, not here.

BrainMVP paper (Appendix I, Classification):
  "we resize the input image to a fixed size of 128 × 128 × 64 after
   normalizing it ... Subsequently, we randomly crop a fixed region of
   96 × 96 × 64 ... In the inference stage, we crop an area of
   96 × 96 × 64 at the center of the input image."

Input:
    .../smriprep/sub-{sub}/ses-{ses}/anat/
        sub-{sub}_ses-{ses}_space-MNI152NLin2009cAsym_res-1_desc-preproc_T1w.nii.gz

Output:
    derivatives/brainmvp_inputs/sub-{sub}/ses-{ses}/
        sub-{sub}_ses-{ses}_space-BrainMVP128x64_desc-preproc_T1w.nii.gz

Usage:
    python mri_pipeline/brain_mvp/03_prepare_BrainMVP.py --long all --n-workers 4
    python mri_pipeline/brain_mvp/03_prepare_BrainMVP.py --session bl --dry-run
    python mri_pipeline/brain_mvp/03_prepare_BrainMVP.py \\
        --smriprep-dir /path/to/smriprep --out-root /path/to/brainmvp_inputs
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

# ── Paths ──────────────────────────────────────────────────────────────────────
# LOCAL (Windows) paths
PROJECT_ROOT   = r"D:\ADNI_BIDS_project"
SMRIPREP_DIR   = os.path.join(PROJECT_ROOT, "derivatives", "smriprep",
                              "hpc_first_lex_2026_04_20", "mni_staging")
OUT_ROOT       = os.path.join(PROJECT_ROOT, "derivatives", "brainmvp_inputs")
MANIFEST_PATH  = os.path.join(OUT_ROOT, "brainmvp_manifest.csv")
SNP_TSV        = os.path.join(PROJECT_ROOT, "bids", "genotype",
                              "subjects_with_snp_and_mri.tsv")

# ── BrainMVP target dimensions ────────────────────────────────────────────────
TARGET_SPACING = (1.0, 1.0, 1.0)    # mm — isometric 1mm
TARGET_SHAPE   = (128, 128, 64)      # voxels — BrainMVP classification spec

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Prepare BrainMVP inputs (96³ @ 1mm) from sMRIprep derivatives.")
parser.add_argument("--dry-run",   action="store_true")
parser.add_argument("--overwrite", action="store_true")
parser.add_argument("--n-workers", type=int, default=1)
parser.add_argument("--sub", type=str, default=None,
                    help="Process a single subject, e.g. '002S0413'.")
parser.add_argument("--smriprep-dir", type=str, default=None)
parser.add_argument("--out-root", type=str, default=None)
session_group = parser.add_mutually_exclusive_group()
session_group.add_argument("--session", type=str, default="bl")
session_group.add_argument("--long", type=str, default=None)
args = parser.parse_args()

DRY_RUN   = args.dry_run
OVERWRITE = args.overwrite
N_WORKERS = args.n_workers

if args.smriprep_dir:
    SMRIPREP_DIR = args.smriprep_dir
if args.out_root:
    OUT_ROOT = args.out_root
    MANIFEST_PATH = os.path.join(OUT_ROOT, "brainmvp_manifest.csv")

# Resolve session mode
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
        n_years = int(args.long)
        if n_years <= 0:
            parser.error(f"--long must be 'all' or a positive integer (got {args.long!r})")
        LONG_MODE = "cutoff"
        MAX_MONTHS = n_years * 12


def session_to_months(ses_label: str):
    if ses_label == "bl":
        return 0
    m = re.match(r"^m(\d+)$", ses_label)
    return int(m.group(1)) if m else None


# ── Helper: find sMRIprep T1w ─────────────────────────────────────────────────
def find_smriprep_t1w(sub: str, ses: str):
    sub_label = f"sub-{sub}"
    ses_label = f"ses-{ses}"

    ses_anat_dir = os.path.join(SMRIPREP_DIR, sub_label, ses_label, "anat")
    expected = os.path.join(
        ses_anat_dir,
        f"{sub_label}_{ses_label}_space-MNI152NLin2009cAsym_res-1_desc-preproc_T1w.nii.gz")
    if os.path.isfile(expected):
        return expected
    hits = glob.glob(os.path.join(ses_anat_dir, "*MNI152*desc-preproc_T1w.nii.gz"))
    if hits:
        return hits[0]

    if ses == "bl":
        anat_dir = os.path.join(SMRIPREP_DIR, sub_label, "anat")
        expected = os.path.join(
            anat_dir,
            f"{sub_label}_space-MNI152NLin2009cAsym_res-1_desc-preproc_T1w.nii.gz")
        if os.path.isfile(expected):
            return expected
        hits = glob.glob(os.path.join(anat_dir, "*MNI152*desc-preproc_T1w.nii.gz"))
        if hits:
            return hits[0]
    return None


# ── Core processing ──────────────────────────────────────────────────────────
def process_subject(sub: str, ses: str, dry_run: bool = False,
                    overwrite: bool = False) -> dict:
    """
    BrainMVP preprocessing for one (subject, session) pair:
      RAS → 1mm isometric → percentile clip [1st, 99th] → [0,1] → crop → 96³
    """
    sub_label = f"sub-{sub}"
    ses_label = f"ses-{ses}"

    record = {
        "bids_sub": sub, "bids_ses": ses,
        "t1w_path": "", "output_path": "",
        "input_shape": "", "input_spacing": "", "input_orient": "",
        "output_shape": str(TARGET_SHAPE), "output_spacing": str(TARGET_SPACING),
        "brain_voxels": "", "pct_1": "", "pct_99": "",
        "status": "unknown", "error": "",
        "timestamp": datetime.datetime.now().isoformat(),
    }

    out_dir = os.path.join(OUT_ROOT, sub_label, ses_label)
    out_name = f"{sub_label}_{ses_label}_space-BrainMVP128x64_desc-preproc_T1w.nii.gz"
    out_path = os.path.join(out_dir, out_name)
    record["output_path"] = out_path

    if os.path.isfile(out_path) and not overwrite:
        record["status"] = "skipped_exists"
        return record

    t1w_path = find_smriprep_t1w(sub, ses)
    record["t1w_path"] = t1w_path or ""

    if not t1w_path:
        record["status"] = "missing_input"
        record["error"] = f"sMRIprep MNI T1w not found for {sub_label} {ses_label}"
        return record

    try:
        src_img = nib.load(t1w_path)
        record["input_shape"] = str(src_img.shape)
        record["input_spacing"] = str(tuple(
            round(float(z), 4) for z in src_img.header.get_zooms()[:3]))
        record["input_orient"] = str(nib.aff2axcodes(src_img.affine))
    except Exception:
        pass

    if dry_run:
        record["status"] = "dry_run"
        return record

    try:
        # ── Step 1-4: Load → RAS → 1mm isometric ─────────────────────────────
        initial_transform = mt.Compose([
            mt.LoadImaged(keys=["image"]),
            mt.EnsureChannelFirstd(keys=["image"]),
            mt.Orientationd(keys=["image"], axcodes="RAS"),
            mt.Spacingd(keys=["image"], pixdim=TARGET_SPACING, mode="bilinear"),
        ])
        data_dict = initial_transform({"image": t1w_path})
        volume = data_dict["image"]

        # ── Step 5: Percentile clip [1st, 99th] ──────────────────────────────
        vol_np = volume.numpy().squeeze()
        nonzero = vol_np[vol_np > 0]
        if len(nonzero) > 0:
            pct_1 = float(np.percentile(nonzero, 1))
            pct_99 = float(np.percentile(nonzero, 99))
        else:
            pct_1, pct_99 = 0.0, 1.0
        record["pct_1"] = round(pct_1, 4)
        record["pct_99"] = round(pct_99, 4)

        vol_np = np.clip(vol_np, pct_1, pct_99)

        # ── Step 6: Rescale to [0, 1] ─────────────────────────────────────────
        denom = pct_99 - pct_1
        if denom > 0:
            vol_np = (vol_np - pct_1) / denom
        else:
            vol_np = np.zeros_like(vol_np)

        # ── Step 7-8: Crop foreground → Resize to 96³ ────────────────────────
        import torch
        vol_tensor = torch.from_numpy(vol_np).unsqueeze(0).float()  # (1, D, H, W)
        crop_resize = mt.Compose([
            mt.CropForegroundd(keys=["image"], source_key="image"),
            mt.Resized(keys=["image"], spatial_size=TARGET_SHAPE),
        ])
        result = crop_resize({"image": vol_tensor})
        final_vol = result["image"].numpy().squeeze()  # (96, 96, 96)

        # ── Record stats ─────────────────────────────────────────────────────
        brain_mask = final_vol > 1e-6
        record["brain_voxels"] = int(brain_mask.sum())
        record["output_shape"] = str(final_vol.shape)

        # ── Save ──────────────────────────────────────────────────────────────
        os.makedirs(out_dir, exist_ok=True)
        affine = np.eye(4)
        for i in range(3):
            affine[i, i] = TARGET_SPACING[i]
            affine[i, 3] = -TARGET_SPACING[i] * TARGET_SHAPE[i] / 2.0
        out_img = nib.Nifti1Image(final_vol, affine=affine)
        nib.save(out_img, out_path)
        record["status"] = "ok"

    except Exception:
        record["status"] = "failed"
        record["error"] = traceback.format_exc(limit=3)

    return record


# ── Discover subject/session pairs ────────────────────────────────────────────
def discover_pairs():
    pairs = []
    snp_subjects = None
    if os.path.isfile(SNP_TSV):
        tsv = pd.read_csv(SNP_TSV, sep='\t')
        snp_subjects = set(
            tsv['participant_id'].str.replace('sub-', '', regex=False).tolist())
        print(f"  SNP-matched subjects loaded: {len(snp_subjects)}")
    else:
        print(f"  WARNING: SNP TSV not found at {SNP_TSV}; processing ALL subjects.")

    if not os.path.isdir(SMRIPREP_DIR):
        print(f"ERROR: sMRIprep directory not found: {SMRIPREP_DIR}")
        sys.exit(1)

    n_sessions_skipped = 0
    for sub_dir in sorted(os.listdir(SMRIPREP_DIR)):
        if not sub_dir.startswith("sub-") or "_" in sub_dir:
            continue
        sub_root = os.path.join(SMRIPREP_DIR, sub_dir)
        if not os.path.isdir(sub_root):
            continue
        sub = sub_dir[4:]

        if args.sub and sub != args.sub:
            continue
        if snp_subjects and sub not in snp_subjects:
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

    if LONG_MODE == "cutoff":
        print(f"  Sessions skipped (beyond m{MAX_MONTHS}): {n_sessions_skipped}")
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
    print("03_prepare_BrainMVP — BrainMVP inputs from sMRIprep")
    print("=" * 65)
    print(f"  sMRIprep dir  : {SMRIPREP_DIR}")
    print(f"  Output root   : {OUT_ROOT}")
    print(f"  Mode          : {mode_str}")
    print(f"  Target spacing: {TARGET_SPACING} mm (isometric 1mm)")
    print(f"  Target shape  : {TARGET_SHAPE}")
    print(f"  Normalization : Percentile clip [1st, 99th] → [0, 1]")
    print(f"  Dry run       : {DRY_RUN}")
    print(f"  Overwrite     : {OVERWRITE}")
    print(f"  Workers       : {N_WORKERS}")
    print()

    pairs = discover_pairs()
    n_subs = len({s for s, _ in pairs})
    print(f"Pairs to process: {len(pairs)} (across {n_subs} unique subjects)")
    print()

    if not pairs:
        print("No (subject, session) pairs found.")
        sys.exit(0)

    os.makedirs(OUT_ROOT, exist_ok=True)

    existing_records = {}
    if os.path.isfile(MANIFEST_PATH):
        old_df = pd.read_csv(MANIFEST_PATH, dtype=str)
        for _, row in old_df.iterrows():
            key = (row.get("bids_sub", ""), row.get("bids_ses", ""))
            existing_records[key] = row.to_dict()
        print(f"  Resuming: {len(existing_records)} records already in manifest.")

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
                  f"pct=[{rec['pct_1']}, {rec['pct_99']}]")
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
                        "error": str(exc), "t1w_path": "", "output_path": "",
                        "input_shape": "", "input_spacing": "", "input_orient": "",
                        "output_shape": str(TARGET_SHAPE),
                        "output_spacing": str(TARGET_SPACING),
                        "brain_voxels": "", "pct_1": "", "pct_99": "",
                        "timestamp": datetime.datetime.now().isoformat(),
                    }
                handle_record(rec)
    else:
        for i, (sub, ses) in enumerate(pairs):
            print(f"[{i+1}/{len(pairs)}] sub-{sub}_ses-{ses} ...", end=" ",
                  flush=True)
            rec = process_subject(sub, ses, DRY_RUN, OVERWRITE)
            handle_record(rec)

    manifest_df = pd.DataFrame(list(all_records.values()))
    manifest_cols = [
        "bids_sub", "bids_ses", "status",
        "t1w_path", "output_path",
        "input_shape", "input_spacing", "input_orient",
        "output_shape", "output_spacing",
        "brain_voxels", "pct_1", "pct_99",
        "error", "timestamp",
    ]
    present = [c for c in manifest_cols if c in manifest_df.columns]
    extras = [c for c in manifest_df.columns if c not in manifest_cols]
    manifest_df = manifest_df[present + extras]
    manifest_df.to_csv(MANIFEST_PATH, index=False)

    print()
    print("=" * 65)
    print(f"Done.  Pairs processed: {len(pairs)}")
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
