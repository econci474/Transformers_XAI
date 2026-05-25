"""
Script 01 — Prepare BrainDINO Inputs from sMRIprep Outputs
==========================================================
Converts sMRIprep MNI152NLin2009cAsym outputs into the format required by
the BrainDINO 2D-axial ViT downstream pipeline (Wu et al., arXiv:2604.27277,
Section 4.3).

Paper-fidelity finetune recipe (NOT the pretrain recipe — pretrain uses
256³ + Gaussian slice sampling, which is implemented differently and is
not needed for downstream classification):

  > "Each 3D MRI volume was resampled to a spatial resolution of
  > 128×256×256, uniformly sampled into 128 axial slices, processed
  > independently by the pretrained 2D ViT encoder, and aggregated into
  > subject-level representations via mean pooling of slice-wise CLS
  > tokens."          -- Section 4.3

Pre-resampling steps adopt Section 4.1.2's normalisation. The paper's
HD-BET skull-strip is replaced by sMRIprep's explicit brain mask
(`_desc-brain_mask.nii.gz`). **Important deviation from the sibling
ViT/BrainMVP prep scripts**: despite the sidecar JSON's
'SkullStripped: true', the desc-preproc_T1w retains ~5% non-zero
skull residue OUTSIDE the explicit brain mask (verified empirically
on sub-002S0729 ses-m36: 4.78% of voxels are non-zero T1w but lie
outside brain_mask, with 0% gaps inside the mask). We therefore
multiply T1w by the explicit mask before any further processing. An
assertion checks the masks for consistency and writes the diff
fraction to the manifest (info-only; aborts only on catastrophic
disagreement >20%, e.g. shape mismatch or wrong-subject pairing).

Processing chain (per subject, using MONAI transforms):
  1. Load sMRIprep MNI-space preproc T1w + matching brain mask
  2. Mask consistency check (assert (T1w != 0) ~= (mask > 0))
  3. Channel-first  -> (1, H, W, D)
  4. Reorient to RAS -- axis 2 = axial (S = inferior-superior) -- IMPORTANT:
     this convention determines which axis the model slices along
  5. Multiply T1w by mask (strict skull-strip; zeros out the residue)
  6. CropForeground on the masked image (the brain extent)
  7. Resize to (256, 256, 128) in RAS-axis order = 128 axial x 256 cor x
     256 sag (paper's "spatial resolution of 128x256x256" with the axial
     axis = 128). The training script then permutes (C, R, A, S) ->
     (C, S, A, R) = (C, axial, cor, sag) for slice iteration.
  8. Percentile clip nonzero voxels to [0.5, 99.5]
  9. Per-volume z-score using nonzero mean/std
  10. Restore background zeros (so the foreground mask stays intact)

Target specs (from paper Section 4.3 + Section 4.1.2):
  shape (RAS axis order):  (256, 256, 128) = (R, A, S=axial)
  voxel dimensions:        128 axial × 256 cor × 256 sag (paper notation)
  orientation:             RAS
  intensity clip:          [0.5th, 99.5th] percentile on nonzero voxels
  normalisation:           per-volume z-score on nonzero voxels

Input (sessionwise sMRIprep, e.g. CSD3):
    .../smriprep_sessionwise/smriprep/sub-{sub}/ses-{ses}/anat/
        sub-{sub}_ses-{ses}_space-MNI152NLin2009cAsym_res-1_desc-preproc_T1w.nii.gz
        sub-{sub}_ses-{ses}_space-MNI152NLin2009cAsym_res-1_desc-brain_mask.nii.gz

Input (subject-level fallback, e.g. local Windows non-sessionwise):
    .../smriprep/.../sub-{sub}/anat/
        sub-{sub}_space-MNI152NLin2009cAsym_res-1_desc-preproc_T1w.nii.gz
        sub-{sub}_space-MNI152NLin2009cAsym_res-1_desc-brain_mask.nii.gz
    (Used only when --session bl and no ses-bl/ subdirectory exists.)

Output (BIDS-correct, includes ses-*):
    derivatives/braindino_inputs/sub-{sub}/ses-{ses}/
        sub-{sub}_ses-{ses}_space-MNI128x256_desc-braindino_T1w.nii.gz

Manifest CSV:
    derivatives/braindino_inputs/braindino_manifest.csv

Session selection (mutually exclusive):
    --session SESSION   Single session to process. Default: 'bl' (baseline).
    --long N            'all' = every session on disk; N = up to ses-m{12*N}.

Usage:
    python 01_prepare_braindino_inputs.py                       # baseline
    python 01_prepare_braindino_inputs.py --session bl --n-workers 4
    python 01_prepare_braindino_inputs.py --long all --n-workers 4
    python 01_prepare_braindino_inputs.py --dry-run
    python 01_prepare_braindino_inputs.py --overwrite
    python 01_prepare_braindino_inputs.py --sub 002S0413
"""

import os
import re
import sys
import glob
import argparse
import datetime
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

# torch MUST come before numpy/pandas/monai on Windows: otherwise MKL DLL
# collision (the same MKL gets loaded twice and the second load mismatches
# `shm.dll`). On Linux this ordering doesn't matter.
import torch  # noqa: F401
import numpy as np
import nibabel as nib
import pandas as pd
import monai.transforms as mt

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(THIS_DIR))   # mri_pipeline -> repo root
sys.path.insert(0, REPO_ROOT)
from bidsification.exclusions import is_excluded_subject  # noqa: E402

# ── Paths ──────────────────────────────────────────────────────────────────────
# LOCAL (Windows) paths — used when running on your local machine
PROJECT_ROOT   = r"D:\ADNI_BIDS_project"
SMRIPREP_DIR   = os.path.join(PROJECT_ROOT, "derivatives", "smriprep",
                              "smriprep_sessionwise")
OUT_ROOT       = os.path.join(PROJECT_ROOT, "derivatives", "braindino_inputs")
MANIFEST_PATH  = os.path.join(OUT_ROOT, "braindino_manifest.csv")
SNP_TSV        = os.path.join(PROJECT_ROOT, "bids", "genotype",
                              "subjects_with_snp_and_mri.tsv")

# HPC (Linux) paths — uncomment when running on Cambridge HPC. Outputs
# live under ADNI_SMRIPREP/derivatives/ (alongside the existing
# brainmvp_inputs/ and vit_inputs/), NOT under ADNI_MRI/ -- inputs are
# preprocessed-from-sMRIprep, training outputs (model checkpoints,
# wandb, metrics) are what lives in ADNI_MRI/.
# PROJECT_ROOT   = "/home/ec474/rds/hpc-work/ADNI_SMRIPREP"
# SMRIPREP_DIR   = os.path.join(PROJECT_ROOT, "derivatives",
#                               "smriprep_sessionwise", "smriprep")
# OUT_ROOT       = os.path.join(PROJECT_ROOT, "derivatives", "braindino_inputs")
# MANIFEST_PATH  = os.path.join(OUT_ROOT, "braindino_manifest.csv")
# SNP_TSV        = os.path.join(PROJECT_ROOT, "bids", "genotype",
#                               "subjects_with_snp_and_mri.tsv")

# ── BrainDINO target dimensions (paper Section 4.3) ───────────────────────────
# In RAS axis order: (R=256, A=256, S=128 = axial).
# Reported in paper notation: 128 axial x 256 coronal x 256 sagittal.
TARGET_SHAPE_RAS = (256, 256, 128)
PERCENTILE_CLIP  = (0.5, 99.5)


def build_transform() -> mt.Compose:
    """The pre-clip preprocessing pipeline. The intensity clip + z-score
    runs in `process_subject` because MONAI's NormalizeIntensityd does
    not support percentile clipping, and the order matters:
    mask -> bbox -> resize -> clip -> z-score.

    Both T1w and brain_mask are loaded. T1w is multiplied by the mask
    to enforce the strict skull-strip (the desc-preproc_T1w alone
    retains ~5% skull-residue voxels despite the SkullStripped sidecar
    -- see docstring at top of file)."""
    return mt.Compose([
        mt.LoadImaged(keys=["image", "mask"]),
        mt.EnsureChannelFirstd(keys=["image", "mask"]),
        mt.Orientationd(keys=["image", "mask"], axcodes="RAS"),
        # Multiply T1w by mask: enforce the strict skull-strip.
        mt.MaskIntensityd(keys=["image"], mask_key="mask"),
        # Crop to the brain bounding box BEFORE resize so the brain
        # fills the target volume (no wasted voxels on background).
        mt.CropForegroundd(keys=["image", "mask"], source_key="mask",
                           select_fn=lambda x: x > 0),
        # Resize to the paper's target. spatial_size order matches the
        # MONAI tensor (after Orient) = (R, A, S = axial).
        mt.Resized(keys=["image", "mask"], spatial_size=TARGET_SHAPE_RAS,
                   mode=("trilinear", "nearest")),
    ])


def percentile_clip_and_zscore(volume: np.ndarray, mask: np.ndarray):
    """Per-volume percentile clip + nonzero z-score (paper Section 4.1.2).
    Returns (out_volume, stats_dict). The mask is the resized brain
    mask from `build_transform`."""
    brain = mask > 0
    if brain.sum() == 0:
        raise RuntimeError("brain mask is empty after preprocessing")

    brain_vals = volume[brain]
    p_lo, p_hi = np.percentile(brain_vals, PERCENTILE_CLIP)
    clipped = np.where(brain, np.clip(volume, p_lo, p_hi), 0.0)

    # z-score on the CLIPPED nonzero voxels (so post-z-score values
    # respect the percentile bounds).
    brain_vals_clipped = clipped[brain]
    mean = float(brain_vals_clipped.mean())
    std  = float(brain_vals_clipped.std())
    if std < 1e-8:
        raise RuntimeError(f"degenerate z-score: std={std:.3e}")
    z = (clipped - mean) / std
    out = np.where(brain, z, 0.0).astype(np.float32)
    return out, {
        "pct_0_5":   round(float(p_lo), 4),
        "pct_99_5":  round(float(p_hi), 4),
        "znorm_mean_nonzero": round(float(out[brain].mean()), 6),
        "znorm_std_nonzero":  round(float(out[brain].std()), 6),
        "brain_voxels": int(brain.sum()),
    }


# ── Argument parsing ───────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Prepare BrainDINO inputs from sMRIprep derivatives "
                "(paper-fidelity finetune recipe, arXiv:2604.27277 Sec 4.3).")
parser.add_argument("--dry-run",   action="store_true",
                    help="Print plan without writing files.")
parser.add_argument("--overwrite", action="store_true",
                    help="Reprocess subjects that already have output.")
parser.add_argument("--n-workers", type=int, default=1,
                    help="Number of parallel worker processes (default: 1).")
parser.add_argument("--sub", type=str, default=None,
                    help="Process a single subject only, e.g. '002S0413'.")
parser.add_argument("--max-subjects", type=int, default=None,
                    help="Limit to first N subjects (for smoke testing).")
parser.add_argument("--smriprep-dir", type=str, default=None,
                    help="Override SMRIPREP_DIR.")
parser.add_argument("--out-root", type=str, default=None,
                    help="Override OUT_ROOT.")
parser.add_argument("--no-mask-check", action="store_true",
                    help="Skip the T1w/brain_mask consistency check (info-only; "
                         "the pipeline always uses the explicit brain_mask).")
session_group = parser.add_mutually_exclusive_group()
session_group.add_argument("--session", type=str, default="bl",
                    help="Single session label (default: 'bl').")
session_group.add_argument("--long", type=str, default=None,
                    help="Longitudinal mode: 'all' or integer N (= up to m{12N}).")
args = parser.parse_args()

DRY_RUN   = args.dry_run
OVERWRITE = args.overwrite
N_WORKERS = args.n_workers

if args.smriprep_dir:
    SMRIPREP_DIR = args.smriprep_dir
if args.out_root:
    OUT_ROOT = args.out_root
    MANIFEST_PATH = os.path.join(OUT_ROOT, "braindino_manifest.csv")

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
    if ses_label == "bl":
        return 0
    m = re.match(r"^m(\d+)$", ses_label)
    if m:
        return int(m.group(1))
    return None


# ── Safety check: log T1w-vs-mask consistency (info-only) ─────────────────────
CATASTROPHIC_MASK_DIFF = 0.20  # 20% disagreement = abort (wrong subject / shape)


def check_mask_consistency(t1w_path: str, mask_path: str,
                           sub: str, ses: str) -> dict:
    """Compare (T1w != 0) to (mask > 0) and report the diff fraction.

    Empirically (sub-002S0729 ses-m36): ~5% of voxels are non-zero in
    T1w but lie outside the explicit brain_mask -- skull residue that
    survived sMRIprep's skull-strip despite the sidecar's
    'SkullStripped: true'. This is expected; we report the diff to the
    manifest for auditing and only abort on catastrophic disagreement
    (> 20%, which would indicate shape mismatch or wrong-subject
    pairing). The pipeline ALWAYS uses the explicit brain_mask via
    MaskIntensityd downstream -- the diff is informational only."""
    t1 = nib.load(t1w_path).get_fdata()
    mk = nib.load(mask_path).get_fdata()
    if t1.shape != mk.shape:
        raise RuntimeError(
            f"T1w shape {t1.shape} != mask shape {mk.shape} for "
            f"sub-{sub} ses-{ses} -- shape mismatch, aborting.")
    inferred = (t1 != 0)
    explicit = (mk > 0)
    n_voxels = int(inferred.size)
    infer_only = int((inferred & ~explicit).sum())   # skull residue
    mask_only  = int((~inferred & explicit).sum())   # gaps inside mask
    diff_frac  = float((infer_only + mask_only) / n_voxels)
    if diff_frac > CATASTROPHIC_MASK_DIFF:
        raise RuntimeError(
            f"Catastrophic T1w/mask disagreement for sub-{sub} ses-{ses}: "
            f"{diff_frac * 100:.2f}% of voxels differ "
            f"(infer_only={infer_only}, mask_only={mask_only}). Likely "
            f"shape mismatch or wrong-subject pairing.")
    return {
        "mask_check":              "ok",
        "mask_diff_frac":          round(diff_frac, 6),
        "mask_skull_residue_frac": round(infer_only / n_voxels, 6),
        "mask_gaps_inside_frac":   round(mask_only / n_voxels, 6),
    }


# ── Helper: find sMRIprep T1w + brain mask ────────────────────────────────────
def find_smriprep_pair(sub: str, ses: str):
    """Returns (t1w_path, mask_path) or (None, None) if either missing."""
    sub_label = f"sub-{sub}"
    ses_label = f"ses-{ses}"
    suffix_t1 = "_space-MNI152NLin2009cAsym_res-1_desc-preproc_T1w.nii.gz"
    suffix_mk = "_space-MNI152NLin2009cAsym_res-1_desc-brain_mask.nii.gz"

    def _try_dir(d, prefix):
        t1 = os.path.join(d, f"{prefix}{suffix_t1}")
        mk = os.path.join(d, f"{prefix}{suffix_mk}")
        if os.path.isfile(t1) and os.path.isfile(mk):
            return t1, mk
        t1_hits = glob.glob(os.path.join(d, "*MNI152*desc-preproc_T1w.nii.gz"))
        mk_hits = glob.glob(os.path.join(d, "*MNI152*desc-brain_mask.nii.gz"))
        if t1_hits and mk_hits:
            return t1_hits[0], mk_hits[0]
        return None, None

    # Sessionwise layout (HPC)
    ses_anat = os.path.join(SMRIPREP_DIR, sub_label, ses_label, "anat")
    t1, mk = _try_dir(ses_anat, f"{sub_label}_{ses_label}")
    if t1:
        return t1, mk

    # Subject-level fallback (baseline only)
    if ses == "bl":
        anat = os.path.join(SMRIPREP_DIR, sub_label, "anat")
        t1, mk = _try_dir(anat, sub_label)
        if t1:
            return t1, mk

    return None, None


# ── Core processing function ──────────────────────────────────────────────────
def process_subject(sub: str, ses: str, dry_run: bool = False,
                    overwrite: bool = False) -> dict:
    sub_label = f"sub-{sub}"
    ses_label = f"ses-{ses}"

    record = {
        "bids_sub":       sub,
        "bids_ses":       ses,
        "t1w_path":       "",
        "mask_path":      "",
        "output_path":    "",
        "input_shape":    "",
        "input_spacing":  "",
        "input_orient":   "",
        "output_shape":   "",
        "brain_voxels":   "",
        "znorm_mean_nonzero": "",
        "znorm_std_nonzero":  "",
        "pct_0_5":        "",
        "pct_99_5":       "",
        "mask_check":     "",
        "mask_diff_frac": "",
        "mask_skull_residue_frac": "",
        "mask_gaps_inside_frac":   "",
        "status":         "unknown",
        "error":          "",
        "timestamp":      datetime.datetime.now().isoformat(),
    }

    out_dir  = os.path.join(OUT_ROOT, sub_label, ses_label)
    out_name = f"{sub_label}_{ses_label}_space-MNI128x256_desc-braindino_T1w.nii.gz"
    out_path = os.path.join(out_dir, out_name)
    record["output_path"] = out_path

    if os.path.isfile(out_path) and not overwrite:
        record["status"] = "skipped_exists"
        return record

    t1w_path, mask_path = find_smriprep_pair(sub, ses)
    record["t1w_path"]  = t1w_path or ""
    record["mask_path"] = mask_path or ""
    if not t1w_path or not mask_path:
        record["status"] = "missing_input"
        record["error"]  = (f"sMRIprep MNI preproc T1w or brain mask not found "
                            f"for {sub_label} {ses_label}")
        return record

    try:
        src_img = nib.load(t1w_path)
        record["input_shape"]   = str(src_img.shape)
        record["input_spacing"] = str(tuple(round(float(z), 4)
                                            for z in src_img.header.get_zooms()[:3]))
        record["input_orient"]  = str(nib.aff2axcodes(src_img.affine))
    except Exception:
        pass

    if dry_run:
        record["status"] = "dry_run"
        return record

    try:
        # Info-only sanity check: log T1w-vs-mask disagreement to the manifest.
        # Aborts only on catastrophic mismatch (>20%, i.e. shape error or
        # wrong-subject pairing). The pipeline ALWAYS uses the explicit
        # brain_mask downstream regardless of the diff size.
        if not args.no_mask_check:
            record.update(check_mask_consistency(t1w_path, mask_path, sub, ses))

        transform = build_transform()
        data = transform({"image": t1w_path, "mask": mask_path})

        # MONAI tensors are (C, H, W, D); squeeze channel for the numeric op
        vol = data["image"].numpy().squeeze(0).astype(np.float32)
        msk = data["mask"].numpy().squeeze(0).astype(np.float32)

        out, stats = percentile_clip_and_zscore(vol, msk)
        record.update(stats)
        record["output_shape"] = str(out.shape)

        os.makedirs(out_dir, exist_ok=True)
        # Build affine for RAS at 1 mm isotropic with centre at origin.
        # Voxel spacing is implicit -- the model only sees the values, not
        # physical units. Centre at origin for visualisation friendliness.
        affine = np.eye(4)
        for i in range(3):
            affine[i, 3] = -float(out.shape[i]) / 2.0
        nib.save(nib.Nifti1Image(out, affine=affine), out_path)
        record["status"] = "ok"

    except Exception:
        record["status"] = "failed"
        record["error"]  = traceback.format_exc(limit=3)

    return record


# ── Discover subject/session pairs ────────────────────────────────────────────
def discover_pairs():
    pairs = []
    snp_subjects = None
    if os.path.isfile(SNP_TSV):
        tsv = pd.read_csv(SNP_TSV, sep='\t')
        snp_subjects = set(
            tsv['participant_id'].str.replace('sub-', '', regex=False).tolist()
        )
        print(f"  SNP-matched subjects loaded: {len(snp_subjects)}")

    sub_dirs = sorted(glob.glob(os.path.join(SMRIPREP_DIR, "sub-*")))
    if args.sub:
        sub_dirs = [d for d in sub_dirs if d.endswith(f"sub-{args.sub}")]

    for sub_dir in sub_dirs:
        sub = os.path.basename(sub_dir).replace("sub-", "")
        if is_excluded_subject(sub):
            continue
        if snp_subjects is not None and sub not in snp_subjects:
            continue

        if LONG_MODE is None:
            pairs.append((sub, SINGLE_SESSION))
            continue

        # Longitudinal: glob ses-* subdirs
        ses_dirs = sorted(glob.glob(os.path.join(sub_dir, "ses-*")))
        for sd in ses_dirs:
            ses = os.path.basename(sd).replace("ses-", "")
            months = session_to_months(ses)
            if months is None:
                continue
            if MAX_MONTHS is not None and months > MAX_MONTHS:
                continue
            pairs.append((sub, ses))

    if args.max_subjects is not None:
        # Limit to first N unique subjects (keeps per-subject grouping intact).
        seen, kept = set(), []
        for s, e in pairs:
            if len(seen) >= args.max_subjects and s not in seen:
                break
            seen.add(s)
            kept.append((s, e))
        pairs = kept

    return pairs


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    print(f"BrainDINO preprocessing — paper finetune recipe")
    print(f"  SMRIPREP_DIR : {SMRIPREP_DIR}")
    print(f"  OUT_ROOT     : {OUT_ROOT}")
    print(f"  Target shape : {TARGET_SHAPE_RAS} (R, A, S=axial), "
          f"= paper '128 ax x 256 cor x 256 sag'")
    print(f"  Mode         : "
          + ("dry-run" if DRY_RUN else ("overwrite" if OVERWRITE else "skip-if-exists")))
    if LONG_MODE is None:
        print(f"  Session      : {SINGLE_SESSION}")
    else:
        print(f"  Longitudinal : {LONG_MODE}"
              + (f" (cap m{MAX_MONTHS})" if MAX_MONTHS else ""))

    pairs = discover_pairs()
    print(f"  Discovered   : {len(pairs)} (sub, ses) pairs")
    if not pairs:
        print("Nothing to do.")
        return

    records = []
    if N_WORKERS <= 1:
        for sub, ses in pairs:
            r = process_subject(sub, ses, dry_run=DRY_RUN, overwrite=OVERWRITE)
            records.append(r)
            tag = r["status"][:14]
            print(f"  [{tag:>14s}] sub-{sub} ses-{ses}  "
                  f"brain={r['brain_voxels']}  "
                  f"mean={r['znorm_mean_nonzero']}  std={r['znorm_std_nonzero']}")
    else:
        with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
            futures = {pool.submit(process_subject, s, e,
                                   dry_run=DRY_RUN, overwrite=OVERWRITE): (s, e)
                       for s, e in pairs}
            for fut in as_completed(futures):
                r = fut.result()
                records.append(r)
                tag = r["status"][:14]
                print(f"  [{tag:>14s}] sub-{r['bids_sub']} ses-{r['bids_ses']}  "
                      f"brain={r['brain_voxels']}")

    df = pd.DataFrame(records)
    if not DRY_RUN:
        os.makedirs(OUT_ROOT, exist_ok=True)
        # Merge with existing manifest if present (incremental runs).
        if os.path.isfile(MANIFEST_PATH):
            old = pd.read_csv(MANIFEST_PATH, dtype=str)
            key = ["bids_sub", "bids_ses"]
            old = old[~old.set_index(key).index.isin(df.set_index(key).index)]
            df = pd.concat([old, df], ignore_index=True)
        df.to_csv(MANIFEST_PATH, index=False)
        print(f"  Manifest     : {MANIFEST_PATH}  ({len(df)} rows total)")

    print("\nStatus counts:")
    print(df["status"].value_counts().to_string())


if __name__ == "__main__":
    main()
