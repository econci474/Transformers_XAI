"""
05_extract_brainmvp_full_ft_embeddings.py
==========================================
Extract per-scan embeddings + per-class probabilities from a BrainMVP
full_ft checkpoint (specific task / augment / seed) for use in downstream
multimodal modelling.

Runs on CSD3 (or any host with the BrainMVP env + checkpoints +
brainmvp_inputs 128x128x64 volumes). For 4 tasks the user asked about:

  T1_binary    : aug = stochastic, plus_original
  T1b_binary   : aug = stochastic, none
  T1d_binary   : aug = plus_original
  T2_multiclass: aug = stochastic

Total: 6 (task, aug) combinations x 3 seeds = 18 cells. SLURM array submit
script `05_extract_brainmvp_full_ft_submit_csd3.sh` decodes the array ID
into (task, aug, seed).

Pipeline per (task, aug, seed):
  1. Resolve checkpoint best_model.pt from the standard tree
     (the BrainMVP "debug" double-nest path; flag --checkpoint to override).
  2. Build BrainMVPClassifier (uniformer-small + AdaptiveAvgPool3d + Linear).
  3. Load state_dict. Move to GPU.
  4. Iterate scans from the post-exclusion MRI master (1456 rows). For each,
     resolve the 128x128x64 nifti, load + transform (center-crop 96x96x64),
     forward through model. Hook on pool layer captures the 512-d embedding;
     model output is the logits; softmax for probs; argmax for prediction.
  5. Save NPZ + CSV with full identifier set (bids_sub, Patient_ID, RID,
     adni_viscode, VISCODE_2, mri_date) + label + per-seed split + probs.

Output schema (matches `brain_dino/05_extract_embeddings_braindino_frozen.py` in the
sibling top-level mri_pipeline directory, so downstream code reads them the
same way):
  embeddings   : (N, 512) float32
  logits       : (N, K)   float32
  probs        : (N, K)   float32
  pred         : (N,)     int8
  bids_sub     : (N,)     str    sub-XXXSYYYY
  Patient_ID   : (N,)     str    XXX_S_YYYY
  RID          : (N,)     int32  (parsed from Patient_ID)
  adni_viscode : (N,)     str    scan VISCODE_2 (e.g. m12)
  VISCODE_2    : (N,)     str    clinical-side VISCODE_2 (may equal adni_viscode)
  mri_date     : (N,)     str    ISO date
  label        : (N,)     int8   task-specific (-1 if not eligible)
  split        : (N,)     str    train/val/test/'' for the chosen seed
  conversion_group : (N,) str
  model_meta   : dict     task, augment, seed, checkpoint_path, embed_dim, n_classes

Usage on CSD3:
  python mri_pipeline/brain_mvp/05_extract_brainmvp_full_ft_embeddings.py \
      --task T1_binary --augment stochastic --seed 0
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import warnings
from pathlib import Path

# On Windows-local testing only the import-order quirk matters; on CSD3 it's
# harmless. Always import torch before numpy/pandas to match the local script.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import torch  # noqa: E402
import torch.nn as nn  # noqa: E402
import torch.nn.functional as F  # noqa: E402

# CSD3: many DataLoader workers + pin_memory exhausts the open-FD ulimit, which
# surfaces as "received 0 items of ancdata" / "Pin memory thread exited". The
# file_system sharing strategy avoids passing storages as file descriptors.
torch.multiprocessing.set_sharing_strategy("file_system")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import monai.transforms as mt  # noqa: E402
from monai.data import Dataset as MonaiDataset, DataLoader  # noqa: E402

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------- #
# Repo wiring (BrainMVP imports + uniformer model)
# --------------------------------------------------------------------------- #
THIS_DIR = Path(__file__).resolve().parent
MRI_DIR = THIS_DIR.parent
REPO_ROOT = MRI_DIR.parent
sys.path.insert(0, str(MRI_DIR))
sys.path.insert(0, str(REPO_ROOT))

from brain_mvp.uniformer_blocks import uniformer_small  # noqa: E402

# --------------------------------------------------------------------------- #
# Defaults (CSD3 paths; override via CLI for local testing)
# --------------------------------------------------------------------------- #
DEFAULT_INPUTS_DIR = Path(
    "/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/brainmvp_inputs"
)
DEFAULT_CKPT_ROOT = Path(
    "/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug"
)
DEFAULT_MRI_MASTER = Path(
    "/home/ec474/rds/hpc-work/ADNI_MRI/"
    "master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
)
DEFAULT_CLINICAL_SPLITS = Path(
    "/home/ec474/rds/hpc-work/ADNI_CL/"
    "no_cdr_stratified_post_exclusion/tabular/baseline"
)
DEFAULT_OUT_DIR = Path("/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_embeddings")

# BrainMVP crop spec
BMVP_CROP = (96, 96, 64)


# --------------------------------------------------------------------------- #
# Model: same class as 04_supervised_finetuning_BrainMVP.py
# --------------------------------------------------------------------------- #
class BrainMVPClassifier(nn.Module):
    def __init__(self, n_classes: int, drop_rate: float = 0.1):
        super().__init__()
        self.encoder = uniformer_small(in_chans=1)
        self.pool = nn.AdaptiveAvgPool3d(1)
        self.dropout = nn.Dropout(drop_rate)
        self.head = nn.Linear(512, n_classes)

    def forward(self, x):
        features = self.encoder(x)
        x4 = features[-1]
        pooled = self.pool(x4).flatten(1)        # (B, 512)
        return self.head(self.dropout(pooled))


def eval_transform() -> mt.Compose:
    """Same as 04_supervised_finetuning_BrainMVP.py's _bmvp_eval_transform."""
    return mt.Compose([
        mt.LoadImaged(keys=["image"], ensure_channel_first=True, image_only=False),
        mt.EnsureTyped(keys=["image"]),
        mt.Orientationd(keys=["image"], axcodes="RAS"),
        mt.CenterSpatialCropd(keys=["image"], roi_size=BMVP_CROP),
        mt.ToTensord(keys=["image"]),
    ])


# --------------------------------------------------------------------------- #
# Helpers (mirror BrainDINO local script for consistent output schema)
# --------------------------------------------------------------------------- #
def n_classes_for_task(task: str) -> int:
    return 3 if task in ("T2_multiclass", "T4_conv_horizon") else 2


def ptid_to_rid(ptid: str) -> int | None:
    if not isinstance(ptid, str):
        return None
    try:
        return int(ptid.rsplit("_", 1)[-1])
    except (ValueError, IndexError):
        return None


def patient_id_to_bids_sub(ptid: str) -> str:
    return "sub-" + ptid.replace("_", "")


def resolve_brainmvp_path(patient_id: str, viscode: str, inputs_dir: Path) -> str:
    sub = patient_id_to_bids_sub(patient_id)
    ses_label = f"ses-{viscode}"
    return str(inputs_dir / sub / ses_label /
               f"{sub}_{ses_label}_space-BrainMVP128x64_desc-preproc_T1w.nii.gz")


def task_label(task: str, row: pd.Series) -> int | None:
    visit_diag = row.get("Label_visit_diag")
    cg = row.get("conversion_group")
    if task == "T1_binary":
        if visit_diag == "CN":            return 0
        if visit_diag in ("MCI", "AD"):   return 1
        return None
    if task == "T1b_binary":
        if visit_diag == "CN":            return 0
        if visit_diag == "AD":            return 1
        return None
    if task == "T1c_binary":
        if visit_diag == "MCI":           return 0
        if visit_diag == "AD":            return 1
        return None
    if task == "T1d_binary":
        if visit_diag != "MCI":           return None
        if cg == "sMCI":                  return 0
        if cg == "pMCI":                  return 1
        return None
    if task == "T2_multiclass":
        return {"CN": 0, "MCI": 1, "AD": 2}.get(visit_diag, None)
    # conversion horizons (baseline-anchored, subject-level Label_Ny / Label_T4)
    if task in ("T3a_conv3y", "T3b_conv5y", "T3c_conv7y"):
        if row.get("Label_bl_multi") == "AD":   return None     # filter_non_ad=True
        col = {"T3a_conv3y": "Label_3y", "T3b_conv5y": "Label_5y", "T3c_conv7y": "Label_7y"}[task]
        v = row.get(col)
        return int(v) if pd.notna(v) else None
    if task == "T4_conv_horizon":
        v = row.get("Label_T4")
        return int(v) if pd.notna(v) else None
    return None


def load_split_for_seed(splits_root: Path, seed: int) -> dict[str, str]:
    """Patient_ID → train/val/test for the given seed."""
    out: dict[str, str] = {}
    seed_dir = splits_root / f"seed_{seed}"
    for fold in ("train", "val", "test"):
        f = seed_dir / f"{fold}.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f, low_memory=False)
        if "Patient_ID" not in df.columns:
            continue
        for p in df["Patient_ID"].dropna().unique():
            out[p] = fold
    return out


def resolve_checkpoint(ckpt_root: Path, task: str, aug: str, seed: int) -> Path:
    """Standard layout: <ckpt_root>/aug_<aug>/BrainMVP_uniformer/<task>/seed_<n>/full_ft/best_model.pt"""
    return ckpt_root / f"aug_{aug}" / "BrainMVP_uniformer" / task / f"seed_{seed}" / "full_ft" / "best_model.pt"


# --------------------------------------------------------------------------- #
# Inference loop
# --------------------------------------------------------------------------- #
def extract(task: str, aug: str, seed: int, device: torch.device,
            inputs_dir: Path, ckpt_root: Path, master_csv: Path,
            splits_root: Path, out_dir: Path, batch_size: int = 4,
            num_workers: int = 4) -> int:
    # 1. Resolve checkpoint
    ckpt = resolve_checkpoint(ckpt_root, task, aug, seed)
    if not ckpt.exists():
        print(f"[ERROR] Checkpoint not found: {ckpt}", file=sys.stderr)
        return 2
    print(f"  ckpt: {ckpt}")

    # 2. Load model
    K = n_classes_for_task(task)
    model = BrainMVPClassifier(n_classes=K, drop_rate=0.1)
    sd = torch.load(ckpt, map_location="cpu", weights_only=False)
    # 04_supervised_finetuning_BrainMVP.py saves {"net": state_dict, "epoch": ...}.
    # Other/older checkpoints may use "model"/"state_dict" or be a raw state_dict.
    if isinstance(sd, dict) and "net" in sd:
        state = sd["net"]
    elif isinstance(sd, dict) and "model" in sd:
        state = sd["model"]
    elif isinstance(sd, dict) and "state_dict" in sd:
        state = sd["state_dict"]
    else:
        state = sd
    missing, unexpected = model.load_state_dict(state, strict=False)
    # best_model.pt INCLUDES the task-specific head, so a correct checkpoint loads
    # with NO missing keys. Any missing key means we unwrapped the wrong container
    # key and the model would run on RANDOM weights -- fail loud rather than emit
    # meaningless probabilities (this previously slipped through as a warning).
    if missing:
        ck = list(sd.keys()) if isinstance(sd, dict) else type(sd).__name__
        raise RuntimeError(
            f"load_state_dict missing {len(missing)} keys (e.g. {missing[:5]}); "
            f"checkpoint top-level keys were {ck}. Refusing to run on "
            f"partially-initialised weights.")
    if unexpected:
        print(f"  [warn] load_state_dict ignored unexpected keys: {unexpected[:5]}",
              file=sys.stderr)
    model.to(device).eval()

    # 3. Hook pool layer for embedding capture
    embed_buffer: list[torch.Tensor] = []

    def _hook(module, inp, out):
        # out: (B, 512, 1, 1, 1) -> we want flatten(1) -> (B, 512)
        embed_buffer.append(out.flatten(1).detach().cpu())

    h = model.pool.register_forward_hook(_hook)

    # 4. Build full-cohort iteration list from MRI master
    master = pd.read_csv(master_csv, low_memory=False)
    print(f"  master: {master.shape}; unique Patient_ID: {master.Patient_ID.nunique()}")

    items, row_meta, valid = [], [], []
    n_path_missing = 0
    for _, r in master.iterrows():
        ptid = r["Patient_ID"]
        ses = r["adni_viscode"]
        ipath = resolve_brainmvp_path(ptid, ses, inputs_dir)
        exists = os.path.isfile(ipath)
        if not exists:
            n_path_missing += 1
        items.append({"image": ipath, "_idx": len(items)})
        valid.append(exists)
        row_meta.append({
            "bids_sub":         r.get("bids_sub"),
            "Patient_ID":       ptid,
            "RID":              ptid_to_rid(ptid),
            "adni_viscode":     ses,
            "VISCODE_2":        r.get("VISCODE_2"),
            "mri_date":         r.get("mri_date"),
            "conversion_group": r.get("conversion_group"),
            "label_raw":        task_label(task, r),
            "image_path":       ipath,
            "image_exists":     exists,
        })
    N = len(items)
    print(f"  scans: {N}; missing input volume on disk: {n_path_missing}")

    # 5. Iterate dataloader over scans that have files; NaN-fill missing ones
    valid_items = [it for it, v in zip(items, valid) if v]
    if not valid_items:
        print(f"[ERROR] No valid input volumes found for the cohort.", file=sys.stderr)
        return 3
    ds = MonaiDataset(data=valid_items, transform=eval_transform())
    loader = DataLoader(ds, batch_size=batch_size, shuffle=False,
                        num_workers=num_workers, pin_memory=(device.type == "cuda"))

    all_logits = torch.full((N, K), float("nan"))
    n_done = 0
    print(f"  running inference ({device}) ...")
    with torch.no_grad():
        for batch in loader:
            x = batch["image"].to(device, non_blocking=True).float()
            # _idx is a list of ints in the batch
            idxs = batch["_idx"].tolist() if torch.is_tensor(batch["_idx"]) else list(batch["_idx"])
            logits = model(x)                       # (B, K)
            for i_batch, i_global in enumerate(idxs):
                all_logits[i_global] = logits[i_batch].detach().cpu()
            n_done += x.shape[0]
            if n_done % 200 < x.shape[0]:
                print(f"    progress: {n_done}/{len(valid_items)}")

    h.remove()

    # 6. Stitch embeddings back to row order from buffer
    # The hook captured one tensor per BATCH (not per scan), so concat them.
    embeddings_valid = torch.cat(embed_buffer, dim=0)   # (V, 512)
    assert embeddings_valid.shape[0] == len(valid_items), \
        f"embedding count mismatch: {embeddings_valid.shape[0]} vs {len(valid_items)}"

    all_embeddings = torch.full((N, embeddings_valid.shape[1]), float("nan"))
    j = 0
    for i_global, v in enumerate(valid):
        if v:
            all_embeddings[i_global] = embeddings_valid[j]
            j += 1

    # 7. Softmax probs + argmax
    nan_mask = torch.isnan(all_logits).any(dim=1)
    safe_logits = torch.where(nan_mask.unsqueeze(1), torch.zeros_like(all_logits), all_logits)
    probs = F.softmax(safe_logits, dim=1)
    preds = probs.argmax(dim=1)
    probs[nan_mask] = float("nan")
    all_logits[nan_mask] = float("nan")
    preds[nan_mask] = -1

    # 8. Per-seed splits + labels
    ptid_to_split = load_split_for_seed(splits_root, seed)
    splits = [ptid_to_split.get(rm["Patient_ID"], "") for rm in row_meta]
    labels = [rm["label_raw"] if rm["label_raw"] is not None else -1 for rm in row_meta]

    # 9. Save NPZ + CSV
    out_subdir = out_dir / task / f"aug_{aug}" / f"seed_{seed}"
    out_subdir.mkdir(parents=True, exist_ok=True)
    npz_path = out_subdir / f"embeddings_seed_{seed}.npz"
    model_meta = {
        "task":            task,
        "augment":         aug,
        "seed":            seed,
        "model":           "BrainMVP_uniformer_full_ft",
        "embed_dim":       int(embeddings_valid.shape[1]),
        "n_classes":       K,
        "checkpoint_path": str(ckpt),
        "n_scans":         N,
        "n_missing":       int(nan_mask.sum().item()),
        "miss_examples":   [(rm["Patient_ID"], rm["adni_viscode"])
                            for rm in row_meta if not rm["image_exists"]][:10],
    }
    np.savez(
        npz_path,
        embeddings=all_embeddings.numpy().astype("float32"),
        logits=all_logits.numpy().astype("float32"),
        probs=probs.numpy().astype("float32"),
        pred=preds.numpy().astype("int8"),
        bids_sub=np.array([rm["bids_sub"] for rm in row_meta], dtype=object),
        Patient_ID=np.array([rm["Patient_ID"] for rm in row_meta], dtype=object),
        RID=np.array([rm["RID"] if rm["RID"] is not None else -1 for rm in row_meta], dtype="int32"),
        adni_viscode=np.array([rm["adni_viscode"] for rm in row_meta], dtype=object),
        VISCODE_2=np.array([rm["VISCODE_2"] for rm in row_meta], dtype=object),
        mri_date=np.array([str(rm["mri_date"]) for rm in row_meta], dtype=object),
        label=np.array(labels, dtype="int8"),
        split=np.array(splits, dtype=object),
        conversion_group=np.array([rm["conversion_group"] for rm in row_meta], dtype=object),
        model_meta=np.array(json.dumps(model_meta)),
    )
    print(f"  wrote {npz_path}  ({npz_path.stat().st_size/1e6:.2f} MB)")

    # CSV (no embedding columns)
    flat = pd.DataFrame({
        "bids_sub":        [rm["bids_sub"] for rm in row_meta],
        "Patient_ID":      [rm["Patient_ID"] for rm in row_meta],
        "RID":             [rm["RID"] for rm in row_meta],
        "adni_viscode":    [rm["adni_viscode"] for rm in row_meta],
        "VISCODE_2":       [rm["VISCODE_2"] for rm in row_meta],
        "mri_date":        [str(rm["mri_date"]) for rm in row_meta],
        "label":           labels,
        "split":           splits,
        "conversion_group": [rm["conversion_group"] for rm in row_meta],
        "pred":            preds.numpy(),
    })
    for k in range(K):
        flat[f"prob_class_{k}"] = probs[:, k].numpy()
        flat[f"logit_class_{k}"] = all_logits[:, k].numpy()
    flat.to_csv(out_subdir / f"embeddings_seed_{seed}.csv", index=False)

    # 10. Verification
    n_ok = int((~nan_mask).sum().item())
    probs_sums = probs.sum(dim=1)
    ok_sum = float(probs_sums[~nan_mask].mean().item()) if n_ok else float("nan")
    print(f"  verify: ok_rows={n_ok}/{N}  nan_rows={int(nan_mask.sum().item())}  "
          f"emb_dim={all_embeddings.shape[1]}  probs_sum_mean={ok_sum:.6f}")
    valid_lab = np.array(labels) >= 0
    if valid_lab.any():
        labs = np.array(labels)[valid_lab]
        ps = preds.numpy()[valid_lab]
        print(f"  label dist (valid): {pd.Series(labs).value_counts().to_dict()}")
        print(f"  pred==label fraction (all splits mixed): {(labs == ps).mean():.3f}")

    return 0


# --------------------------------------------------------------------------- #
# Main / CLI
# --------------------------------------------------------------------------- #
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--task", required=True,
                   choices=["T1_binary", "T1b_binary", "T1c_binary",
                            "T1d_binary", "T2_multiclass",
                            "T3a_conv3y", "T3b_conv5y", "T3c_conv7y", "T4_conv_horizon"])
    p.add_argument("--augment", required=True,
                   choices=["none", "stochastic", "plus_original"])
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--inputs-dir", type=Path, default=DEFAULT_INPUTS_DIR,
                   help="BrainMVP 128x64 inputs root.")
    p.add_argument("--ckpt-root", type=Path, default=DEFAULT_CKPT_ROOT,
                   help="BrainMVP full_ft checkpoint root "
                        "(contains aug_<aug>/BrainMVP_uniformer/...).")
    p.add_argument("--master-csv", type=Path, default=DEFAULT_MRI_MASTER,
                   help="Post-exclusion MRI master CSV.")
    p.add_argument("--splits-root", type=Path, default=DEFAULT_CLINICAL_SPLITS,
                   help="Per-seed clinical splits root "
                        "(contains seed_0/{train,val,test}.csv etc.).")
    p.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                   help="Output root for per-task NPZ + CSV.")
    p.add_argument("--batch-size", type=int, default=4)
    p.add_argument("--num-workers", type=int, default=4)
    p.add_argument("--device", default=None,
                   help="cuda or cpu; default: cuda if available else cpu.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)

    print(f"[05_brainmvp_full_ft] task={args.task} aug={args.augment} seed={args.seed}")
    print(f"  device: {device}")
    print(f"  inputs: {args.inputs_dir}")
    print(f"  ckpts:  {args.ckpt_root}")
    print(f"  master: {args.master_csv}")
    print(f"  splits: {args.splits_root}")
    print(f"  out:    {args.out_dir}")
    rc = extract(
        task=args.task, aug=args.augment, seed=args.seed,
        device=device, inputs_dir=args.inputs_dir, ckpt_root=args.ckpt_root,
        master_csv=args.master_csv, splits_root=args.splits_root,
        out_dir=args.out_dir, batch_size=args.batch_size,
        num_workers=args.num_workers,
    )
    return rc


if __name__ == "__main__":
    sys.exit(main())
