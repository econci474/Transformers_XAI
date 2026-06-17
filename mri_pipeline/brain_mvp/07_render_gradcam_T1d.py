"""
07_render_gradcam_T1d.py — overlay T1d Grad-CAM heatmaps on the real MNI T1w
============================================================================
Stage 2 of the T1d BrainMVP Grad-CAM. Consumes the per-scan 96^3 saliency arrays
written by `06_gradcam_T1d.py` and renders heatmap overlays on the subject's
*real, full-resolution* sMRIprep MNI T1w for visual interpretation.

Why back-projection: `06` saved the CAM in the model's own squished input space
(96x96x64, S-axis compressed to 64), and `03_prepare_BrainMVP.py` built that
volume with a synthetic affine and discarded the crop bounds — so the saved
volume cannot be inverted directly. Instead we rebuild a spatial-only, invertible
MONAI chain mirroring 03+06, run it FORWARD on the source T1w (capturing the true
MNI affine + applied_operations), attach the CAM and call `.inverse()` to land
the heatmap on the subject's true `space-MNI152NLin2009cAsym_res-1` grid. Because
every subject's sMRIprep output shares that grid, group means are true-MNI
averages on real anatomy.

Deliberately a plain visual (the reader interprets anatomy directly): no atlas,
no ROI naming. The overlay is semi-transparent (`--alpha`) over the anatomy and
thresholded high (`--thresh`, fraction of map max) so only the strongest
saliency shows.

Outputs (to <gradcam-dir>/figures):
  - per-class group-mean overlay (mean back-projected CAM on the mean T1w) as
    axial / coronal / sagittal montages (PNG + PDF) + the NIfTIs
  - a few confidence-ranked individual exemplar overlays per class

Run in the `mri` env (nilearn + monai). No GPU / no internet needed; runs on the
CSD3 login node where the source T1w, CAMs and nilearn all live.

Usage (CSD3 login node):
  cd /home/ec474/rds/hpc-work/ADNI_MRI/gradcam_T1d
  conda run -n mri python \
    /home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/brain_mvp/07_render_gradcam_T1d.py \
    --gradcam-dir /home/ec474/rds/hpc-work/ADNI_MRI/gradcam_T1d --thresh 0.7 --alpha 0.6
"""
import argparse
import glob
import os
import sys
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

# torch before numpy/pandas (Windows-local DLL-load quirk; matches 05_/06_).
import torch
import numpy as np
import pandas as pd
import nibabel as nib
import monai.transforms as mt

CLASS_NAME = {0: "sMCI", 1: "pMCI"}

# Mirror 03_prepare_BrainMVP.py + 06_gradcam_T1d.py exactly.
TARGET_SPACING = (1.0, 1.0, 1.0)     # mm — 03 isometric resample
TARGET_SHAPE = (128, 128, 64)        # 03 Resize target
BMVP_CROP = (96, 96, 64)             # 06 eval CenterSpatialCrop (model input)

SMRIPREP_DEFAULT = ("/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/"
                    "smriprep_sessionwise/smriprep")


def _clip_rescale(img):
    """Replicate 03's nonzero-percentile clip+rescale so CropForeground sees the
    same intensities (hence the same bounding box) as 03. Intensity-only, so it
    is not inverted; it only affects the forward foreground detection."""
    arr = img.detach().cpu().numpy()
    nz = arr[arr > 0]
    if nz.size:
        p5 = float(np.percentile(nz, 5))
        p95 = float(np.percentile(nz, 95))
    else:
        p5, p95 = 0.0, 1.0
    out = np.clip(arr, p5, p95)
    denom = p95 - p5
    out = (out - p5) / denom if denom > 0 else np.zeros_like(out)
    new = img.clone()
    new.copy_(torch.as_tensor(out, dtype=img.dtype))
    return new


def build_chain() -> mt.Compose:
    """Spatial-only invertible chain: source MNI T1w -> 96^3 model space.
    Forward captures the true affine + applied_operations; inverse maps a CAM
    back onto the source res-1 MNI grid."""
    return mt.Compose([
        mt.LoadImaged(keys=["image"], image_only=False, ensure_channel_first=True),
        mt.Orientationd(keys=["image"], axcodes="RAS"),
        mt.Spacingd(keys=["image"], pixdim=TARGET_SPACING, mode="bilinear"),
        mt.Lambdad(keys=["image"], func=_clip_rescale),  # not invertible (intensity)
        mt.CropForegroundd(keys=["image"], source_key="image"),
        mt.Resized(keys=["image"], spatial_size=TARGET_SHAPE),
        mt.CenterSpatialCropd(keys=["image"], roi_size=BMVP_CROP),
    ], lazy=False)


def find_smriprep_t1w(bids_sub: str, ses: str, root: str):
    """Locate the source sMRIprep MNI T1w (mirrors 03_prepare_BrainMVP.py:122-135)."""
    sub_label = bids_sub if str(bids_sub).startswith("sub-") else f"sub-{bids_sub}"
    ses_label = ses if str(ses).startswith("ses-") else f"ses-{ses}"
    anat = Path(root) / sub_label / ses_label / "anat"
    exp = anat / (f"{sub_label}_{ses_label}_space-MNI152NLin2009cAsym"
                  f"_res-1_desc-preproc_T1w.nii.gz")
    if exp.is_file():
        return str(exp)
    hits = sorted(glob.glob(str(anat / "*MNI152*desc-preproc_T1w.nii.gz")))
    return hits[0] if hits else None


def backproject(cam96: np.ndarray, t1w_path: str, chain: mt.Compose):
    """Run the chain forward on the source T1w, then inverse-map the (96,96,64)
    CAM onto the source grid. Returns (cam_on_source[X,Y,Z], source_affine)."""
    fwd = chain({"image": t1w_path})
    m = fwd["image"]                                   # MetaTensor (1,96,96,64)
    cam_t = torch.as_tensor(cam96[None], dtype=m.dtype)
    cam_mt = m.clone()                                 # carry meta + applied_operations
    cam_mt.copy_(cam_t)
    cam_mt.applied_operations = m.applied_operations
    inv = chain.inverse({"image": cam_mt})
    cam_src = inv["image"]                             # MetaTensor (1,X,Y,Z) on source grid
    arr = cam_src.detach().cpu().numpy()
    arr = arr[0] if arr.ndim == 4 else arr
    return arr.astype(np.float32), np.asarray(cam_src.affine, dtype=np.float64)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gradcam-dir", default="/home/ec474/rds/hpc-work/ADNI_MRI/gradcam_T1d",
                    help="dir with manifest.csv + seed_*/cam_*.npy (from 06)")
    ap.add_argument("--smriprep-dir", default=SMRIPREP_DEFAULT,
                    help="root of sMRIprep MNI T1w (sub-*/ses-*/anat/*_res-1_desc-preproc_T1w.nii.gz)")
    ap.add_argument("--out-dir", default=None, help="default: <gradcam-dir>/figures")
    ap.add_argument("--n-exemplars", type=int, default=3)
    ap.add_argument("--include", nargs="*", default=None,
                    help="Patient_IDs to always render as exemplars regardless of "
                         "ranking / correctness (e.g. --include 005_S_0448)")
    ap.add_argument("--all-scans", action="store_true",
                    help="include misclassified scans in the group mean (default: correct-only)")
    ap.add_argument("--n-cuts", type=int, default=7, help="slices per montage axis")
    ap.add_argument("--thresh", type=float, default=0.7,
                    help="overlay threshold as a fraction of the map max (0 disables)")
    ap.add_argument("--alpha", type=float, default=0.6,
                    help="overlay transparency (0=invisible, 1=opaque) so anatomy shows through")
    ap.add_argument("--sanity", action="store_true",
                    help="back-project an all-ones CAM for the first scan and print Dice vs brain mask")
    args = ap.parse_args()

    gdir = Path(args.gradcam_dir)
    out = Path(args.out_dir) if args.out_dir else gdir / "figures"
    out.mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(gdir / "manifest.csv")
    correct_only = not args.all_scans
    print(f"[07] manifest: {len(man)} scans; classes {man['label'].value_counts().to_dict()}; "
          f"correct_only={correct_only}; thresh={args.thresh} alpha={args.alpha}")

    from nilearn import plotting

    chain = build_chain()

    def resolve_cam(p, seed):
        if os.path.isfile(p):
            return p
        alt = gdir / f"seed_{int(seed)}" / Path(p).name
        return str(alt)

    def _overlay(bg_img, stat_img, mode, title, thr):
        """High-res anatomy + a single semi-transparent hot overlay/colorbar via
        plot_stat_map (it manages the one-colorbar-per-figure rule). `dim=0`
        keeps the anatomy at full brightness; overlay transparency is set via
        `transparency` (nilearn >=0.13) or `alpha` (older); one-sided 0..vmax cbar."""
        kw = dict(bg_img=bg_img, display_mode=mode, cut_coords=args.n_cuts,
                  colorbar=True, threshold=thr, cmap="hot", title=title,
                  black_bg=False, dim=0, symmetric_cbar=False, vmax=1.0)
        for name in ("transparency", "alpha"):   # transparency renamed alpha in 0.13
            try:
                return plotting.plot_stat_map(stat_img, **{name: args.alpha}, **kw)
            except TypeError:
                continue
        return plotting.plot_stat_map(stat_img, **kw)

    def montages(stat_img, bg_img, title, stem, thr):
        for mode, tag in (("z", "axial"), ("y", "coronal"), ("x", "sagittal")):
            disp = _overlay(bg_img, stat_img, mode, f"{title} [{tag}]", thr)
            disp.savefig(out / f"{stem}_{tag}.png", dpi=200)
            if tag == "axial":
                disp.savefig(out / f"{stem}_{tag}.pdf")
            disp.close()

    # ── optional inversion sanity check ──────────────────────────────────────────
    if args.sanity:
        r0 = man.iloc[0]
        t1w0 = find_smriprep_t1w(r0["bids_sub"], r0["adni_viscode"], args.smriprep_dir)
        if t1w0 is None:
            print(f"[sanity] no T1w for {r0['bids_sub']} {r0['adni_viscode']}", file=sys.stderr)
        else:
            ones, _ = backproject(np.ones(BMVP_CROP, dtype=np.float32), t1w0, chain)
            bgd = nib.load(t1w0).get_fdata()
            cam_mask = ones > 0.5
            brain = bgd > 0
            inter = np.logical_and(cam_mask, brain).sum()
            dice = 2 * inter / (cam_mask.sum() + brain.sum() + 1e-9)
            print(f"[sanity] ones-CAM vs brain mask: Dice={dice:.3f} "
                  f"(cam_vox={int(cam_mask.sum())} brain_vox={int(brain.sum())} "
                  f"cam_shape={ones.shape} t1w_shape={bgd.shape})")

    # ── per-class group mean (back-projected CAM + T1w on the shared MNI grid) ────
    cam_sum = {0: None, 1: None}
    t1w_sum = {0: None, 1: None}
    cnt = {0: 0, 1: 0}
    cand = {0: [], 1: []}    # (prob1, pid, vis, cam_path, t1w_path, is_correct)
    include = set(args.include or [])
    miss = 0
    grp_shape = None
    grp_aff = None
    for _, r in man.iterrows():
        lab = int(r["label"])
        pid = r["Patient_ID"]
        is_corr = int(r["correct"]) == 1
        use_for_group = is_corr or not correct_only
        # process the scan if it feeds the group mean OR is a force-included exemplar
        if not (use_for_group or pid in include):
            continue
        cam_p = resolve_cam(r["cam_path"], r["seed"])
        t1w_p = find_smriprep_t1w(r["bids_sub"], r["adni_viscode"], args.smriprep_dir)
        if not (os.path.isfile(cam_p) and t1w_p and os.path.isfile(t1w_p)):
            miss += 1
            continue
        cam96 = np.load(cam_p).astype(np.float32)
        cam_mni, aff = backproject(cam96, t1w_p, chain)
        bg = nib.load(t1w_p).get_fdata().astype(np.float32)
        if grp_shape is None:
            grp_shape, grp_aff = cam_mni.shape, aff
        if cam_mni.shape != grp_shape or bg.shape != grp_shape:
            print(f"  [skip grid mismatch] {r['bids_sub']} {r['adni_viscode']} "
                  f"cam{cam_mni.shape} t1w{bg.shape} vs {grp_shape}", file=sys.stderr)
            miss += 1
            continue
        cand[lab].append((float(r["prob1"]), pid, r["adni_viscode"], cam_p, t1w_p, is_corr))
        if use_for_group:
            cam_sum[lab] = cam_mni if cam_sum[lab] is None else cam_sum[lab] + cam_mni
            t1w_sum[lab] = bg if t1w_sum[lab] is None else t1w_sum[lab] + bg
            cnt[lab] += 1

    print(f"  back-projected: sMCI={cnt[0]} pMCI={cnt[1]} (skipped {miss})")
    if cnt[0] == 0 and cnt[1] == 0:
        print("[ERROR] no scans back-projected — check --gradcam-dir / --smriprep-dir.",
              file=sys.stderr)
        return 2

    for lab in (0, 1):
        if not cnt[lab]:
            continue
        gcam = cam_sum[lab] / cnt[lab]
        ginp = t1w_sum[lab] / cnt[lab]
        rng = float(gcam.max() - gcam.min())
        gcam = (gcam - gcam.min()) / rng if rng > 0 else gcam
        cam_img = nib.Nifti1Image(gcam.astype(np.float32), grp_aff)
        bg_img = nib.Nifti1Image(ginp.astype(np.float32), grp_aff)
        nib.save(cam_img, out / f"gradcam_group_{CLASS_NAME[lab]}.nii.gz")
        nib.save(bg_img, out / f"meaninput_group_{CLASS_NAME[lab]}.nii.gz")
        thr = args.thresh * float(gcam.max()) if args.thresh > 0 else None
        montages(cam_img, bg_img,
                 f"T1d Grad-CAM mean — {CLASS_NAME[lab]} (n={cnt[lab]})",
                 f"gradcam_group_{CLASS_NAME[lab]}", thr)
        print(f"  wrote group {CLASS_NAME[lab]} montages (n={cnt[lab]})")

    # ── exemplars: most-confident correct predictions (+ any --include subjects) ──
    # pMCI panel -> highest p_pMCI; sMCI panel -> lowest p_pMCI (= highest p_sMCI).
    for lab in (0, 1):
        ranked = sorted([c for c in cand[lab] if c[5]], key=lambda t: t[0],
                        reverse=(lab == 1))
        picks = ranked[:args.n_exemplars]
        # force-included subjects, in addition to the top-N (skip dups)
        chosen_ids = {c[1] for c in picks}
        extra = [c for c in cand[lab] if c[1] in include and c[1] not in chosen_ids]
        picks = picks + extra
        if not picks:
            continue
        print(f"  {CLASS_NAME[lab]} exemplars (p_pMCI): "
              f"{['%.2f' % p for p, *_ in picks]}"
              + (f"  [+included: {[c[1] for c in extra]}]" if extra else ""))
        for (p1, pid, vis, cam_p, t1w_p, _corr) in picks:
            cam96 = np.load(cam_p).astype(np.float32)
            cam_mni, aff = backproject(cam96, t1w_p, chain)
            bg = nib.load(t1w_p).get_fdata().astype(np.float32)
            cam_img = nib.Nifti1Image(cam_mni, aff)
            bg_img = nib.Nifti1Image(bg, aff)
            thr = args.thresh * float(cam_mni.max()) if args.thresh > 0 else None
            montages(cam_img, bg_img,
                     f"{CLASS_NAME[lab]} {pid} {vis} (p_pMCI={p1:.2f})",
                     f"exemplar_{CLASS_NAME[lab]}_{pid}_{vis}", thr)
    print(f"  wrote figures to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
