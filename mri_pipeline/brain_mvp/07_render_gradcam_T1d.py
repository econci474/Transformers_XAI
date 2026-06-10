"""
07_render_gradcam_T1d.py — overlay T1d Grad-CAM heatmaps on the model input
===========================================================================
Stage 2 of the T1d BrainMVP Grad-CAM. Consumes the per-scan 96^3 saliency arrays
and the matching 96^3 model-input volumes written by `06_gradcam_T1d.py`, and
renders heatmap overlays for visual interpretation.

Deliberately simple (the reader interprets anatomy directly): no MNI
back-projection, no atlas, no ROI naming. Everything is in the model's own input
space (the RAS 96x96x64 volume the network saw, which `03_prepare_BrainMVP.py`
produced via CropForeground -> Resize, so brains are size-normalised but the
3rd axis is compressed to 64 voxels — slices are the model's view, not true MNI).

Outputs (to <gradcam-dir>/figures):
  - per-class group-mean overlay (mean CAM over correctly-classified scans on the
    mean input) as axial / coronal / sagittal montages (PNG + PDF)
  - a few individual exemplar overlays per class

Run in the `mri` env (nilearn). No GPU / no internet needed; can run on the CSD3
login node, a CPU job, or locally after pulling `gradcam_T1d/` back.

Usage:
  python mri_pipeline/brain_mvp/07_render_gradcam_T1d.py \
      --gradcam-dir /home/ec474/rds/hpc-work/ADNI_MRI/gradcam_T1d
"""
import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import numpy as np
import pandas as pd
import nibabel as nib

CLASS_NAME = {0: "sMCI", 1: "pMCI"}
CROP_SHAPE = (96, 96, 64)

# Canonical RAS 1mm affine for the 96^3 model input (128^3 synthetic affine
# trans (-64,-64,-32), shifted by the +16,+16,0 centre-crop voxel offset).
AFF96 = np.diag([1.0, 1.0, 1.0, 1.0])
AFF96[0, 3], AFF96[1, 3], AFF96[2, 3] = -48.0, -48.0, -32.0


def _img(arr):
    return nib.Nifti1Image(np.asarray(arr, dtype=np.float32), AFF96)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gradcam-dir", default="/home/ec474/rds/hpc-work/ADNI_MRI/gradcam_T1d",
                    help="dir with manifest.csv + seed_*/cam_*.npy + seed_*/input_*.npy (from 06)")
    ap.add_argument("--out-dir", default=None, help="default: <gradcam-dir>/figures")
    ap.add_argument("--n-exemplars", type=int, default=3)
    ap.add_argument("--all-scans", action="store_true",
                    help="include misclassified scans in the group mean (default: correct-only)")
    ap.add_argument("--n-cuts", type=int, default=7, help="slices per montage axis")
    ap.add_argument("--thresh", type=float, default=0.3,
                    help="overlay threshold as a fraction of the map max (0 disables)")
    args = ap.parse_args()

    gdir = Path(args.gradcam_dir)
    out = Path(args.out_dir) if args.out_dir else gdir / "figures"
    out.mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(gdir / "manifest.csv")
    correct_only = not args.all_scans
    print(f"[07] manifest: {len(man)} scans; classes {man['label'].value_counts().to_dict()}; "
          f"correct_only={correct_only}")

    from nilearn import plotting

    def resolve(p, seed, name_prefix):
        if os.path.isfile(p):
            return p
        alt = gdir / f"seed_{int(seed)}" / Path(p).name
        if os.path.isfile(alt):
            return str(alt)
        # last resort: reconstruct name
        return str(alt)

    def montages(stat_img, bg_img, title, stem, thr):
        for mode, tag in (("z", "axial"), ("y", "coronal"), ("x", "sagittal")):
            d = plotting.plot_stat_map(
                stat_img, bg_img=bg_img, display_mode=mode, cut_coords=args.n_cuts,
                colorbar=True, threshold=thr, cmap="hot", title=f"{title} [{tag}]",
                black_bg=False)
            d.savefig(out / f"{stem}_{tag}.png", dpi=200)
            if tag == "axial":
                d.savefig(out / f"{stem}_{tag}.pdf")
            d.close()

    # ── per-class group mean (mean CAM + mean input over selected scans) ──────────
    cam_sum = {0: None, 1: None}
    inp_sum = {0: None, 1: None}
    cnt = {0: 0, 1: 0}
    cand = {0: [], 1: []}   # (prob1, pid, vis, cam_path, inp_path) — confidence-ranked exemplars
    miss = 0
    for _, r in man.iterrows():
        if correct_only and int(r["correct"]) != 1:
            continue
        lab = int(r["label"])
        cam_p = resolve(r["cam_path"], r["seed"], "cam")
        inp_p = resolve(r.get("input_path", ""), r["seed"], "input")
        if not (os.path.isfile(cam_p) and os.path.isfile(inp_p)):
            miss += 1
            continue
        cam = np.load(cam_p).astype(np.float32)
        inp = np.load(inp_p).astype(np.float32)
        cam_sum[lab] = cam if cam_sum[lab] is None else cam_sum[lab] + cam
        inp_sum[lab] = inp if inp_sum[lab] is None else inp_sum[lab] + inp
        cnt[lab] += 1
        cand[lab].append((float(r["prob1"]), r["Patient_ID"], r["adni_viscode"], cam_p, inp_p))

    print(f"  averaged: sMCI={cnt[0]} pMCI={cnt[1]} (skipped {miss})")
    if cnt[0] == 0 and cnt[1] == 0:
        print("[ERROR] no scans found — check --gradcam-dir / that 06 ran.", file=sys.stderr)
        return 2

    for lab in (0, 1):
        if not cnt[lab]:
            continue
        gcam = cam_sum[lab] / cnt[lab]
        ginp = inp_sum[lab] / cnt[lab]
        # renormalise the group-mean CAM to [0,1] for display
        rng = float(gcam.max() - gcam.min())
        gcam = (gcam - gcam.min()) / rng if rng > 0 else gcam
        nib.save(_img(gcam), out / f"gradcam_group_{CLASS_NAME[lab]}.nii.gz")
        nib.save(_img(ginp), out / f"meaninput_group_{CLASS_NAME[lab]}.nii.gz")
        thr = args.thresh * float(gcam.max()) if args.thresh > 0 else None
        montages(_img(gcam), _img(ginp),
                 f"T1d Grad-CAM mean — {CLASS_NAME[lab]} (n={cnt[lab]})",
                 f"gradcam_group_{CLASS_NAME[lab]}", thr)
        print(f"  wrote group {CLASS_NAME[lab]} montages (n={cnt[lab]})")

    # ── exemplars: most-confident correct predictions ────────────────────────────
    # pMCI panel -> highest p_pMCI; sMCI panel -> lowest p_pMCI (= highest p_sMCI).
    for lab in (0, 1):
        ranked = sorted(cand[lab], key=lambda t: t[0], reverse=(lab == 1))
        picks = ranked[:args.n_exemplars]
        print(f"  {CLASS_NAME[lab]} exemplars (p_pMCI): "
              f"{['%.2f' % p for p, *_ in picks]}")
        for (p1, pid, vis, cam_p, inp_p) in picks:
            cam = np.load(cam_p).astype(np.float32)
            inp = np.load(inp_p).astype(np.float32)
            thr = args.thresh * float(cam.max()) if args.thresh > 0 else None
            montages(_img(cam), _img(inp),
                     f"{CLASS_NAME[lab]} {pid} {vis} (p_pMCI={p1:.2f})",
                     f"exemplar_{CLASS_NAME[lab]}_{pid}_{vis}", thr)
    print(f"  wrote figures to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
