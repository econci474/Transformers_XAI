"""
Script 13 — Prepare BrainIAC Inputs from sMRIprep Outputs
==========================================================
Converts sMRIprep MNI152NLin2009cAsym outputs into the exact format
required by the BrainIAC model (eugenehp/brainiac on HuggingFace).

Processing chain (per session):
  1. Find sMRIprep MNI-space T1w in derivatives/smriprep/
  2. Resample to 96×96×96      → trilinear zoom (scipy.ndimage, order=1)
  3. Z-score normalise         → on nonzero voxels only
  4. Save as compressed NIfTI  → derivatives/brainiac_inputs/

Note on skull stripping:
  The sMRIprep `desc-preproc_T1w` in MNI space is ALREADY skull-stripped
  (confirmed by FSLeyes inspection). No brain mask application is needed.
  Background voxels are 0; z-score is computed on nonzero voxels only,
  so the background is unaffected.

Input (per session, from sMRIprep):
    .../smriprep/sub-{sub}/[ses-{ses}/]anat/
        sub-{sub}[_ses-{ses}]_space-MNI152NLin2009cAsym_res-1_desc-preproc_T1w.nii.gz

Output (per session):
    derivatives/brainiac_inputs/sub-{sub}/ses-{ses}/
        sub-{sub}_ses-{ses}_space-MNI96_desc-brainiac_T1w.nii.gz

Manifest CSV:
    derivatives/brainiac_inputs/brainiac_manifest.csv
    Columns: bids_sub, bids_ses, t1w_path, output_path,
             input_shape, output_shape, brain_voxels, znorm_mean, znorm_std,
             status, error

Design decisions:
  - Skull stripping: already done by sMRIprep (MNI preproc T1w = brain only).
  - Input space: MNI-registered sMRIprep output is used (not native-space
    INU/N4 images). BrainIAC requires standard MNI space.
  - Resampling: scipy.ndimage.zoom with order=1 (trilinear), matching the
    BrainIAC model card specification of "trilinear" resizing to 96×96×96.
  - Normalisation: z-score on nonzero voxels only (background stays 0),
    matching BrainIAC spec. sMRIprep does NOT z-score; this step is required.

Usage:
    python 13_prepare_brainiac_inputs.py
    python 13_prepare_brainiac_inputs.py --dry-run       # preview, no writes
    python 13_prepare_brainiac_inputs.py --overwrite     # reprocess all
    python 13_prepare_brainiac_inputs.py --n-workers 4   # parallel jobs
    python 13_prepare_brainiac_inputs.py --sub 002S0413  # single subject
"""

import os
import sys
import glob
import argparse
import datetime
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import nibabel as nib
import pandas as pd
from scipy.ndimage import zoom

# ── Paths ──────────────────────────────────────────────────────────────────────
# LOCAL (Windows) paths — used when running on your local machine
PROJECT_ROOT   = r"D:\ADNI_BIDS_project"
SMRIPREP_DIR   = os.path.join(PROJECT_ROOT, "derivatives", "smriprep", "smriprep")
OUT_ROOT       = os.path.join(PROJECT_ROOT, "derivatives", "brainiac_inputs")
MANIFEST_PATH  = os.path.join(OUT_ROOT, "brainiac_manifest.csv")
SNP_TSV        = os.path.join(PROJECT_ROOT, "bids", "genotype", "subjects_with_snp_and_mri.tsv")

# HPC (Linux) paths — uncomment these and comment out the LOCAL block above
# when running this script directly on the Cambridge HPC via SSH.
# PROJECT_ROOT   = "/home/ec474/rds/hpc-work/ADNI_SMRIPREP"
# SMRIPREP_DIR   = "/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/smriprep/smriprep"
# OUT_ROOT       = "/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/brainiac_inputs"
# MANIFEST_PATH  = OUT_ROOT + "/brainiac_manifest.csv"
# SNP_TSV        = "/home/ec474/rds/hpc-work/ADNI_SMRIPREP/bids/genotype/subjects_with_snp_and_mri.tsv"

# ── BrainIAC target dimensions ─────────────────────────────────────────────────
TARGET_SHAPE = (96, 96, 96)

# ── Argument parsing ───────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Prepare BrainIAC inputs from sMRIprep derivatives.")
parser.add_argument("--dry-run",    action="store_true", help="Print plan without writing files.")
parser.add_argument("--overwrite",  action="store_true", help="Reprocess sessions that already have output.")
parser.add_argument("--n-workers",  type=int, default=1,  help="Number of parallel worker processes (default: 1).")
parser.add_argument("--sub",        type=str, default=None, help="Process a single subject only, e.g. '0001234'.")
args = parser.parse_args()

DRY_RUN   = args.dry_run
OVERWRITE = args.overwrite
N_WORKERS = args.n_workers


# ── Helper: find sMRIprep T1w for one subject ─────────────────────────────────
def find_smriprep_t1w(sub: str, ses: str = None):
    """
    Returns the path to the sMRIprep MNI-space preproc T1w for a given subject.

    sMRIprep's LONGITUDINAL pipeline outputs one cross-session anatomical
    template per subject in sub-{sub}/anat/ (NO session directory). Each
    session only gets a small transform file in ses-{ses}/anat/; the full
    MNI-registered T1w is at the subject level.

    Cross-sectional (single-session) runs may also use subject-level anat/.
    The `ses` argument is accepted for interface compatibility but is not
    used to locate the file.
    """
    sub_label  = f"sub-{sub}"
    anat_dir   = os.path.join(SMRIPREP_DIR, sub_label, "anat")

    # Expected filename: no ses- entity (subject-level output)
    expected = os.path.join(
        anat_dir,
        f"{sub_label}_space-MNI152NLin2009cAsym_res-1_desc-preproc_T1w.nii.gz"
    )
    if os.path.isfile(expected):
        return expected

    # Glob fallback (handles minor filename variations)
    hits = glob.glob(os.path.join(anat_dir, "*MNI152*desc-preproc_T1w.nii.gz"))
    if hits:
        return hits[0]

    return None


# ── Core processing function (one session) ─────────────────────────────────────
def process_session(sub: str, ses: str, dry_run: bool = False, overwrite: bool = False) -> dict:
    """
    Full preprocessing chain for one subject/session.
    Steps: resample to 96×96×96 (trilinear) → z-score on nonzero voxels.
    Skull stripping is NOT applied here — the sMRIprep MNI preproc T1w is
    already skull-stripped (confirmed by FSLeyes inspection).
    Returns a record dict for the manifest CSV.
    """
    record = {
        "bids_sub":     sub,
        "bids_ses":     "template",   # sMRIprep longitudinal: one T1w template per subject
        "t1w_path":     "",
        "output_path":  "",
        "input_shape":  "",
        "output_shape": str(TARGET_SHAPE),
        "brain_voxels": "",
        "znorm_mean":   "",
        "znorm_std":    "",
        "status":       "unknown",
        "error":        "",
        "timestamp":    datetime.datetime.now().isoformat(),
    }

    out_dir  = os.path.join(OUT_ROOT, f"sub-{sub}")
    out_name = f"sub-{sub}_space-MNI96_desc-brainiac_T1w.nii.gz"
    out_path = os.path.join(out_dir, out_name)
    record["output_path"] = out_path

    # ── Check if already done ─────────────────────────────────────────────────
    if os.path.isfile(out_path) and not overwrite:
        record["status"] = "skipped_exists"
        return record

    # ── Locate sMRIprep MNI T1w ───────────────────────────────────────────────
    t1w_path = find_smriprep_t1w(sub, ses)
    record["t1w_path"] = t1w_path or ""

    if not t1w_path:
        record["status"] = "missing_input"
        record["error"]  = f"sMRIprep MNI preproc T1w not found for sub-{sub} ses-{ses}"
        return record

    if dry_run:
        record["status"] = "dry_run"
        return record

    try:
        # ── Step 1: Load T1w ──────────────────────────────────────────────────
        # The sMRIprep desc-preproc T1w in MNI space is already:
        #   - Skull-stripped (brain only, background = 0)
        #   - N4 bias field corrected
        #   - Registered to MNI152NLin2009cAsym at 1mm resolution
        t1w_img  = nib.load(t1w_path)
        t1w_data = t1w_img.get_fdata(dtype=np.float32)
        record["input_shape"] = str(t1w_data.shape)

        # ── Step 2: Resample to 96×96×96 (trilinear) ─────────────────────────
        # scipy.ndimage.zoom with order=1 = trilinear interpolation.
        # Zoom factor per axis: 96 / current_size
        zoom_factors = tuple(t / s for t, s in zip(TARGET_SHAPE, t1w_data.shape))
        resampled = zoom(t1w_data, zoom=zoom_factors, order=1, prefilter=False)

        # NOTE: do NOT clip to 0 here. The sMRIprep MNI image already contains
        # small negative values (e.g. min ~-328) at the brain boundary due to
        # registration interpolation artefacts. These are real brain voxels and
        # must be included in z-score computation. Clipping would turn them into
        # fake zeros that look like skull-stripped background.

        # ── Step 3: Z-score normalise on nonzero voxels only ─────────────────
        # Brain mask = |voxel| > tiny threshold. This correctly includes negative
        # brain-boundary voxels while excluding true background (exact 0s from
        # skull stripping). Background is left as 0 after normalisation.
        brain_mask_96 = np.abs(resampled) > 1e-3
        n_brain_voxels = brain_mask_96.sum()

        if n_brain_voxels == 0:
            record["status"] = "failed"
            record["error"]  = "All voxels are zero after resampling — image may be empty."
            return record

        brain_vals = resampled[brain_mask_96]
        mu    = float(brain_vals.mean())
        sigma = float(brain_vals.std())

        if sigma < 1e-6:
            record["status"] = "failed"
            record["error"]  = f"Near-zero std ({sigma:.2e}) — image may be constant."
            return record

        normalised = resampled.copy()
        normalised[brain_mask_96] = (brain_vals - mu) / sigma
        # Background stays 0 — untouched by normalisation

        record["brain_voxels"] = int(n_brain_voxels)
        record["znorm_mean"]   = round(mu,    4)
        record["znorm_std"]    = round(sigma, 4)

        # ── Step 4: Save ──────────────────────────────────────────────────────
        os.makedirs(out_dir, exist_ok=True)

        # Build a new affine scaled to match the new voxel grid.
        # We preserve the origin (centre of MNI space) and adjust voxel sizes.
        old_affine  = t1w_img.affine
        old_voxsz   = np.abs(old_affine[:3, :3]).max(axis=0)   # approx voxel sizes
        new_voxsz   = np.array([s / t for s, t in
                                 zip([t1w_data.shape[i] * old_voxsz[i] for i in range(3)],
                                     TARGET_SHAPE)])
        scale_factor = new_voxsz / old_voxsz

        new_affine = old_affine.copy()
        new_affine[:3, :3] = old_affine[:3, :3] * scale_factor[np.newaxis, :]

        # Inherit header from source image (preserves units=mm, intent codes, etc.)
        # then update shape, voxel sizes and data type for the new 96³ volume.
        out_hdr = t1w_img.header.copy()
        out_hdr.set_data_dtype(np.float32)
        out_hdr.set_data_shape(normalised.shape)
        out_hdr.set_zooms(tuple(new_voxsz))   # pixdim1-3 in mm

        out_img = nib.Nifti1Image(normalised, affine=new_affine, header=out_hdr)
        nib.save(out_img, out_path)

        record["output_shape"] = str(normalised.shape)
        record["status"]       = "ok"

    except Exception as exc:
        record["status"] = "failed"
        record["error"]  = traceback.format_exc(limit=3)

    return record


# ── Discover all subjects ──────────────────────────────────────────────────────
def discover_sessions():
    """
    Walk the sMRIprep output tree to find all subjects with a subject-level
    MNI T1w (sub-xxx/anat/*MNI152*desc-preproc_T1w.nii.gz).

    sMRIprep's longitudinal pipeline outputs ONE T1w per subject (a
    cross-session anatomical template) in sub-xxx/anat/, NOT one per session.
    Each entry is recorded with ses='template'.
    """
    sessions = []   # list of (sub, ses) tuples; ses is always 'template' here

    # Load SNP-matched subject list if available
    snp_subjects = None
    if os.path.isfile(SNP_TSV):
        tsv = pd.read_csv(SNP_TSV, sep='\t')
        snp_subjects = set(
            tsv['participant_id'].str.replace('sub-', '', regex=False).tolist()
        )
        print(f"  SNP-matched subjects loaded: {len(snp_subjects)}")
    else:
        print(f"  WARNING: SNP TSV not found at {SNP_TSV}; processing ALL subjects.")

    if not os.path.isdir(SMRIPREP_DIR):
        print(f"ERROR: sMRIprep directory not found: {SMRIPREP_DIR}")
        sys.exit(1)

    for sub_dir in sorted(os.listdir(SMRIPREP_DIR)):
        if not sub_dir.startswith("sub-"):
            continue
        sub = sub_dir[4:]  # strip "sub-"

        # Filter to single subject if --sub was given
        if args.sub and sub != args.sub:
            continue

        # Filter to SNP-matched subjects
        if snp_subjects and sub not in snp_subjects:
            continue

        # Only include subjects where the subject-level MNI T1w actually exists
        t1w = find_smriprep_t1w(sub)
        if t1w:
            sessions.append((sub, "template"))
        # else: sMRIprep may not have finished for this subject — will be MISSING

    return sessions


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    print("=" * 65)
    print("Script 13 — Prepare BrainIAC inputs from sMRIprep derivatives")
    print("=" * 65)
    print(f"  sMRIprep dir : {SMRIPREP_DIR}")
    print(f"  Output root  : {OUT_ROOT}")
    print(f"  Target shape : {TARGET_SHAPE}")
    print(f"  Dry run      : {DRY_RUN}")
    print(f"  Overwrite    : {OVERWRITE}")
    print(f"  Workers      : {N_WORKERS}")
    print()

    sessions = discover_sessions()
    print(f"Sessions to process: {len(sessions)}")
    print()

    if not sessions:
        print("No sessions found. Check SMRIPREP_DIR path and that sMRIprep has completed.")
        sys.exit(0)

    os.makedirs(OUT_ROOT, exist_ok=True)

    # ── Load existing manifest to allow resume ────────────────────────────────
    existing_records = {}
    if os.path.isfile(MANIFEST_PATH):
        old_df = pd.read_csv(MANIFEST_PATH, dtype=str)
        for _, row in old_df.iterrows():
            key = (row["bids_sub"], row["bids_ses"])
            existing_records[key] = row.to_dict()
        print(f"  Resuming: {len(existing_records)} sessions already in manifest.")

    # ── Process sessions ──────────────────────────────────────────────────────
    all_records = dict(existing_records)   # will be updated as sessions complete
    n_ok = n_skip = n_fail = n_missing = n_dry = 0

    def handle_record(rec):
        nonlocal n_ok, n_skip, n_fail, n_missing, n_dry
        key = (rec["bids_sub"], rec["bids_ses"])
        all_records[key] = rec
        s = rec["status"]
        if s == "ok":              n_ok      += 1
        elif s.startswith("sk"):   n_skip    += 1
        elif s == "missing_input": n_missing += 1
        elif s == "dry_run":       n_dry     += 1
        else:                      n_fail    += 1

        label = f"sub-{rec['bids_sub']} ses-{rec['bids_ses']}"
        if s == "ok":
            print(f"  [OK]      {label}  brain_voxels={rec['brain_voxels']}  "
                  f"znorm μ={rec['znorm_mean']} σ={rec['znorm_std']}")
        elif s.startswith("sk"):
            print(f"  [SKIP]    {label}")
        elif s == "missing_input":
            print(f"  [MISSING] {label}  — sMRIprep MNI T1w not found")
        elif s == "dry_run":
            print(f"  [DRY-RUN] {label}  ←  {rec['t1w_path']}")
        else:
            print(f"  [FAIL]    {label}  !  {rec['error'][:120]}")

    if N_WORKERS > 1 and not DRY_RUN:
        print(f"Running with {N_WORKERS} parallel worker processes...")
        with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
            futures = {
                pool.submit(process_session, sub, ses, DRY_RUN, OVERWRITE): (sub, ses)
                for sub, ses in sessions
            }
            for future in as_completed(futures):
                try:
                    rec = future.result()
                except Exception as exc:
                    sub, ses = futures[future]
                    rec = {
                        "bids_sub": sub, "bids_ses": ses,
                        "status": "failed", "error": str(exc),
                        "t1w_path": "", "output_path": "",
                        "input_shape": "", "output_shape": str(TARGET_SHAPE),
                        "brain_voxels": "", "znorm_mean": "", "znorm_std": "",
                        "timestamp": datetime.datetime.now().isoformat(),
                    }
                handle_record(rec)
    else:
        for i, (sub, ses) in enumerate(sessions):
            print(f"[{i+1}/{len(sessions)}] sub-{sub} ses-{ses} ...", end=" ", flush=True)
            rec = process_session(sub, ses, DRY_RUN, OVERWRITE)
            handle_record(rec)

    # ── Save manifest ─────────────────────────────────────────────────────────
    manifest_df = pd.DataFrame(list(all_records.values()))
    manifest_cols = [
        "bids_sub", "bids_ses", "status",
        "t1w_path", "output_path",
        "input_shape", "output_shape",
        "brain_voxels", "znorm_mean", "znorm_std",
        "error", "timestamp",
    ]
    # Keep only known columns (in correct order), add any extras at end
    present = [c for c in manifest_cols if c in manifest_df.columns]
    extras  = [c for c in manifest_df.columns if c not in manifest_cols]
    manifest_df = manifest_df[present + extras]
    manifest_df.to_csv(MANIFEST_PATH, index=False)

    # ── Summary ───────────────────────────────────────────────────────────────
    print()
    print("=" * 65)
    print(f"Done.  Sessions processed this run: {len(sessions)}")
    print(f"  OK          : {n_ok}")
    print(f"  Dry-run     : {n_dry}")
    print(f"  Skipped     : {n_skip}")
    print(f"  Missing input: {n_missing}  (sMRIprep MNI output not found)")
    print(f"  Failed      : {n_fail}")
    print(f"  Manifest    : {MANIFEST_PATH}")
    print("=" * 65)

    if n_fail > 0:
        failed = manifest_df[manifest_df["status"] == "failed"][["bids_sub", "bids_ses", "error"]]
        print("\nFailed sessions:")
        print(failed.to_string(index=False))


if __name__ == "__main__":
    main()
