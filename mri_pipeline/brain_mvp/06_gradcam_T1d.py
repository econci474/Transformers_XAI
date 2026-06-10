"""
06_gradcam_T1d.py — 3D Grad-CAM for the BrainMVP T1d model (HPC compute stage)
==============================================================================
Computes per-scan 3D Grad-CAM saliency for the strong T1d (sMCI vs pMCI)
BrainMVP UniFormer-Small full_ft / plus_original model, to localise *where in
the brain* the conversion decision is driven.

This is the GPU compute stage and must run on CSD3, because the 128x128x64
BrainMVP input volumes and the full_ft checkpoints live only there. It writes
small per-scan saliency arrays (96^3 float16) + a manifest; the heavy spatial
work (true-MNI back-projection, Harvard-Oxford ROI quantification, nilearn
figures) is done locally afterwards by `07_render_gradcam_T1d.py`.

Grad-CAM recipe (textbook last-feature-map):
  - target layer = the UniFormer encoder's final feature map `features[-1]`,
    shape (B, 512, 4, 6, 6) (post `encoder.norm`, before global-avg-pool).
  - w_k = mean over the (4,6,6) spatial grid of dY_c/dA_k ; CAM = ReLU(sum_k w_k A_k).
  - upsample (4,6,6) -> the encoder's internal permuted grid (64,96,96) then
    transpose back to the input's (96,96,64) axes (forward_features applies
    `x.permute(0,1,4,2,3)`, i.e. (H,W,D)->(D,H,W); we invert it).
  - target class = the scan's TRUE class (so the sMCI group map shows class-0
    evidence and the pMCI group map shows class-1 evidence). Override with
    --target {true,pred,pmci,smci}.

Cohort (matches the reported "MRI at 12 months" T1d reference): per seed, the
rows of brainmvp_embeddings/T1d_binary/aug_plus_original/seed_<n>/embeddings_seed_<n>.csv
with split==test, label in {0,1}, adni_viscode==m12.

Reuses model/transform/checkpoint loaders from 05_extract_brainmvp_full_ft_embeddings.py
(imported via importlib because the filename starts with a digit) so the model
is built and the checkpoint loaded byte-identically to the embedding extraction.

Output:
  <out>/seed_<n>/cam_<PTID>_<vis>.npy      # (96,96,64) float16, normalised [0,1]
  <out>/manifest.csv                       # one row per scan (ids/label/pred/prob/correct/cam_path)

Usage (CSD3, env: mri):
  python mri_pipeline/brain_mvp/06_gradcam_T1d.py --seeds 0 1 2
  python mri_pipeline/brain_mvp/06_gradcam_T1d.py --seeds 0 --limit 1   # smoke test
"""
import argparse
import importlib.util
import os
import sys
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

# Import torch before numpy/pandas (Windows-local DLL-load quirk; matches 05_*).
import torch
import torch.nn.functional as F
import numpy as np
import pandas as pd
import monai.transforms as mt
from monai.data import Dataset as MonaiDataset

# --------------------------------------------------------------------------- #
# Reuse 05's model / transform / checkpoint loaders (filename starts with a digit)
# --------------------------------------------------------------------------- #
_HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location(
    "bmvp05", _HERE / "05_extract_brainmvp_full_ft_embeddings.py")
b05 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(b05)  # self-bootstraps sys.path + imports uniformer_small

BrainMVPClassifier = b05.BrainMVPClassifier
eval_transform = b05.eval_transform
resolve_checkpoint = b05.resolve_checkpoint
resolve_brainmvp_path = b05.resolve_brainmvp_path
patient_id_to_bids_sub = b05.patient_id_to_bids_sub
n_classes_for_task = b05.n_classes_for_task

# --------------------------------------------------------------------------- #
# Defaults (CSD3 paths; mirror 05_*)
# --------------------------------------------------------------------------- #
DEFAULT_INPUTS_DIR = Path("/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/brainmvp_inputs")
DEFAULT_CKPT_ROOT  = Path("/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug")
DEFAULT_EMB_ROOT   = Path("/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_embeddings")
DEFAULT_OUT_DIR    = Path("/home/ec474/rds/hpc-work/ADNI_MRI/gradcam_T1d")

TASK = "T1d_binary"
AUG = "plus_original"
BMVP_CROP = (96, 96, 64)
FEAT_PERM_SHAPE = (64, 96, 96)   # encoder's internal (D,H,W) after permute(0,1,4,2,3)


def load_model(ckpt_root: Path, seed: int, device: torch.device) -> BrainMVPClassifier:
    """Build BrainMVPClassifier + load the full_ft T1d checkpoint (strict — fail loud)."""
    ckpt = resolve_checkpoint(ckpt_root, TASK, AUG, seed)
    if not ckpt.exists():
        raise FileNotFoundError(f"checkpoint not found: {ckpt}")
    model = BrainMVPClassifier(n_classes=n_classes_for_task(TASK), drop_rate=0.1)
    sd = torch.load(ckpt, map_location="cpu", weights_only=False)
    if isinstance(sd, dict) and "net" in sd:
        state = sd["net"]
    elif isinstance(sd, dict) and "model" in sd:
        state = sd["model"]
    elif isinstance(sd, dict) and "state_dict" in sd:
        state = sd["state_dict"]
    else:
        state = sd
    missing, unexpected = model.load_state_dict(state, strict=False)
    if missing:
        raise RuntimeError(
            f"load_state_dict missing {len(missing)} keys (e.g. {missing[:5]}); "
            f"refusing to run on partially-initialised weights (ckpt={ckpt}).")
    if unexpected:
        print(f"  [warn] ignored unexpected keys: {unexpected[:5]}", file=sys.stderr)
    return model.to(device).eval(), str(ckpt)


def gradcam_one(model: BrainMVPClassifier, x: torch.Tensor, target: int):
    """Return (cam_96x96x64 float32 in [0,1], logits np.array) for a single (1,1,96,96,64) input."""
    captured = {}

    def _fwd_hook(_mod, _inp, out):
        A = out[-1]              # (1, 512, 4, 6, 6)
        A.retain_grad()
        captured["A"] = A

    h = model.encoder.register_forward_hook(_fwd_hook)
    try:
        model.zero_grad(set_to_none=True)
        logits = model(x)                      # (1, K)
        score = logits[0, target]
        score.backward()
        A = captured["A"]                      # (1, 512, 4, 6, 6)
        grads = A.grad                         # (1, 512, 4, 6, 6)
        w = grads.mean(dim=(2, 3, 4), keepdim=True)         # (1, 512, 1, 1, 1)
        cam = F.relu((w * A).sum(dim=1))                    # (1, 4, 6, 6)
        cam = F.interpolate(cam.unsqueeze(1), size=FEAT_PERM_SHAPE,
                            mode="trilinear", align_corners=False)  # (1,1,64,96,96)
        cam = cam.squeeze(0).squeeze(0).detach().cpu().numpy()      # (64,96,96)=(D,H,W)
        cam = np.transpose(cam, (1, 2, 0))                          # (96,96,64)=(H,W,D)
        rng = float(cam.max() - cam.min())
        cam = (cam - cam.min()) / rng if rng > 0 else np.zeros_like(cam)
        return cam.astype(np.float32), logits.detach().cpu().numpy().ravel()
    finally:
        h.remove()


def run_seed(seed: int, args, device: torch.device) -> pd.DataFrame:
    model, ckpt_str = load_model(Path(args.ckpt_root), seed, device)
    emb_csv = (Path(args.emb_root) / TASK / f"aug_{AUG}" / f"seed_{seed}" /
               f"embeddings_seed_{seed}.csv")
    df = pd.read_csv(emb_csv, low_memory=False)
    coh = df[(df["split"] == "test") & (df["label"].isin([0, 1])) &
             (df["adni_viscode"] == args.viscode)].copy()
    print(f"[seed {seed}] ckpt={ckpt_str}")
    print(f"  cohort: {len(coh)} scans (split=test, label in {{0,1}}, viscode={args.viscode})")
    if args.limit:
        coh = coh.head(args.limit)

    tfm = eval_transform()
    out_seed = Path(args.out_dir) / f"seed_{seed}"
    out_seed.mkdir(parents=True, exist_ok=True)

    rows, n_match, n_corr = [], 0, 0
    for _, r in coh.iterrows():
        ptid, vis, label = r["Patient_ID"], r["adni_viscode"], int(r["label"])
        csv_pred = int(r["pred"])
        ipath = resolve_brainmvp_path(ptid, vis, Path(args.inputs_dir))
        if not os.path.isfile(ipath):
            print(f"  [skip] missing input: {ipath}", file=sys.stderr)
            continue
        batch = tfm({"image": ipath})
        x = batch["image"].unsqueeze(0).to(device).float()   # (1,1,96,96,64)

        # target class for the CAM
        if args.target == "true":
            target = label
        elif args.target == "pred":
            target = csv_pred
        elif args.target == "pmci":
            target = 1
        else:  # smci
            target = 0

        cam, logits = gradcam_one(model, x, target)
        pred = int(np.argmax(logits))
        prob = torch.softmax(torch.tensor(logits), dim=0).numpy()
        n_match += int(pred == csv_pred)
        n_corr += int(pred == label)

        cam_path = out_seed / f"cam_{ptid}_{vis}.npy"
        np.save(cam_path, cam.astype(np.float16))
        # also save the 96^3 model input (the volume the network saw, RAS) so the
        # renderer can overlay the heatmap on it without re-touching brainmvp_inputs.
        inp = x.squeeze(0).squeeze(0).detach().cpu().numpy().astype(np.float16)  # (96,96,64)
        inp_path = out_seed / f"input_{ptid}_{vis}.npy"
        np.save(inp_path, inp)
        rows.append({
            "Patient_ID": ptid, "bids_sub": patient_id_to_bids_sub(ptid),
            "seed": seed, "adni_viscode": vis, "label": label,
            "pred": pred, "csv_pred": csv_pred, "prob1": float(prob[1]),
            "correct": int(pred == label), "target_class": target,
            "cam_path": str(cam_path), "input_path": str(inp_path),
            "cam_std": float(cam.std()),
        })

    n = len(rows)
    if n:
        # reproduced balanced accuracy over the cohort (sanity vs the table)
        sub = pd.DataFrame(rows)
        sens = sub[sub.label == 1]["correct"].mean() if (sub.label == 1).any() else float("nan")
        spec = sub[sub.label == 0]["correct"].mean() if (sub.label == 0).any() else float("nan")
        bacc = np.nanmean([sens, spec])
        print(f"  done: {n} scans | argmax==csv_pred {n_match}/{n} | "
              f"reproduced bACC={bacc:.3f} (sens={sens:.3f} spec={spec:.3f})")
        if n_match != n:
            print(f"  [WARN] {n - n_match} scans disagree with stored pred — "
                  f"check checkpoint/transform before trusting the CAMs.", file=sys.stderr)
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--viscode", default="m12")
    ap.add_argument("--target", choices=["true", "pred", "pmci", "smci"], default="true")
    ap.add_argument("--inputs-dir", default=str(DEFAULT_INPUTS_DIR))
    ap.add_argument("--ckpt-root", default=str(DEFAULT_CKPT_ROOT))
    ap.add_argument("--emb-root", default=str(DEFAULT_EMB_ROOT))
    ap.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    ap.add_argument("--limit", type=int, default=0, help="cap scans per seed (smoke test)")
    ap.add_argument("--device", default=None)
    args = ap.parse_args()

    device = torch.device(args.device) if args.device else torch.device(
        "cuda" if torch.cuda.is_available() else "cpu")
    print(f"[06_gradcam_T1d] task={TASK} aug={AUG} target={args.target} device={device}")
    print(f"  inputs: {args.inputs_dir}")
    print(f"  ckpts:  {args.ckpt_root}")
    print(f"  out:    {args.out_dir}")

    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    all_rows = []
    for seed in args.seeds:
        all_rows.append(run_seed(seed, args, device))

    manifest = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    mpath = Path(args.out_dir) / "manifest.csv"
    manifest.to_csv(mpath, index=False)
    print(f"\nwrote manifest: {mpath}  ({len(manifest)} scans across seeds {args.seeds})")
    if len(manifest):
        print(f"  by class: {manifest['label'].value_counts().to_dict()}")
        print(f"  correctly classified: {int(manifest['correct'].sum())}/{len(manifest)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
