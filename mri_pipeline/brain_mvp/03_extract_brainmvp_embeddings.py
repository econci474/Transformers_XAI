"""
03_extract_brainmvp_embeddings.py
==================================
One-time extraction of frozen BrainMVP UniFormer-Small pooled stage-4
embeddings for every preprocessed (sub, ses) scan.

Output: paired Parquet + .pt keyed by (bids_sub, bids_ses) with the
row-aligned 512-d pooled feature vector.

Forward (matches BrainMVPClassifier.forward without the head):
    features = encoder(x)              # tuple of stages
    x4       = features[-1]            # (B, 512, D', H', W')
    pooled   = AdaptiveAvgPool3d(1)(x4).flatten(1)  # (B, 512)

Input: (128, 128, 64) NIfTI from brainmvp_inputs/, then CenterSpatialCrop
to (96, 96, 64) per paper Appendix I (classification ROI). Deterministic
-- no augmentation in extract pass.
"""

import argparse
import datetime
import importlib.util
import os
import sys
import warnings
from pathlib import Path

import torch  # noqa: F401  (Windows MKL DLL order)
import numpy as np
import pandas as pd
import torch.nn as nn
import monai.transforms as mt
from monai.data import Dataset as MonaiDataset
from torch.utils.data import DataLoader

THIS_DIR = Path(__file__).resolve().parent
MRI_DIR = THIS_DIR.parent
REPO_ROOT = MRI_DIR.parent
sys.path.insert(0, str(MRI_DIR))
sys.path.insert(0, str(REPO_ROOT))

# Reuse BrainMVP trainer's encoder builder + checkpoint loader
_bmvp_spec = importlib.util.spec_from_file_location(
    "_bmvp_trainer", THIS_DIR / "04_supervised_finetuning_BrainMVP.py")
_bmvp = importlib.util.module_from_spec(_bmvp_spec)
sys.modules["_bmvp_trainer"] = _bmvp
_bmvp_spec.loader.exec_module(_bmvp)

BrainMVPClassifier        = _bmvp.BrainMVPClassifier
load_brainmvp_pretrained  = _bmvp.load_brainmvp_pretrained
BMVP_CROP                 = _bmvp.BMVP_CROP

# Shared ViT-pipeline utilities
_vit_spec = importlib.util.spec_from_file_location(
    "_vit_pipeline", MRI_DIR / "04_supervised_finetuning_ViT.py")
_vit = importlib.util.module_from_spec(_vit_spec)
sys.modules["_vit_pipeline"] = _vit
_vit_spec.loader.exec_module(_vit)

load_matched_labels    = _vit.load_matched_labels
patient_id_to_bids_sub = _vit.patient_id_to_bids_sub
from bidsification.exclusions import is_excluded_subject

warnings.filterwarnings("ignore")

MODEL_NAME = "BrainMVP"
EMBED_DIM  = 512


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pretrained_ckpt",    type=str, required=True)
    p.add_argument("--brainmvp_inputs_dir", type=str, required=True)
    p.add_argument("--matched_labels_csv", type=str, required=True)
    p.add_argument("--out_dir",            type=str, required=True)
    p.add_argument("--batch_size",         type=int, default=8,
                   help="BrainMVP inputs are smaller (96x96x64 cropped); "
                        "batch=8 fits easily on A100.")
    p.add_argument("--num_workers",        type=int, default=4)
    p.add_argument("--max_subjects",       type=int, default=None)
    p.add_argument("--overwrite",          action="store_true")
    return p.parse_args()


def discover_scans(matched_csv: str, inputs_dir: Path) -> pd.DataFrame:
    matched_df = load_matched_labels(matched_csv)
    matched_df = matched_df[~matched_df["Patient_ID"].apply(
        is_excluded_subject)].copy()

    def _resolve_path(row):
        sub = patient_id_to_bids_sub(row["Patient_ID"])
        ses_label = f"ses-{row['bids_ses']}"
        return str(inputs_dir / sub / ses_label /
                   f"{sub}_{ses_label}_space-BrainMVP128x64_desc-preproc_T1w.nii.gz")

    matched_df["bids_sub"]   = matched_df["Patient_ID"].apply(patient_id_to_bids_sub)
    matched_df["image_path"] = matched_df.apply(_resolve_path, axis=1)
    have = matched_df["image_path"].apply(os.path.isfile)
    print(f"  Scans with preprocessed input on disk: {int(have.sum())}  "
          f"(missing {int((~have).sum())})")
    return matched_df[have].reset_index(drop=True)[
        ["bids_sub", "bids_ses", "image_path"]]


def build_loader(df, batch_size, num_workers):
    items = [{"image": row.image_path, "sub": row.bids_sub, "ses": row.bids_ses}
             for row in df.itertuples()]
    transform = mt.Compose([
        mt.LoadImaged(keys=["image"]),
        mt.EnsureChannelFirstd(keys=["image"]),
        # Paper Appendix I: classification ROI = 96x96x64 center crop
        mt.CenterSpatialCropd(keys=["image"], roi_size=BMVP_CROP),
        mt.ToTensord(keys=["image"]),
    ])
    ds = MonaiDataset(data=items, transform=transform)

    def _collate(batch):
        return {
            "image": torch.stack([b["image"] for b in batch], dim=0),
            "sub":   [b["sub"] for b in batch],
            "ses":   [b["ses"] for b in batch],
        }
    return DataLoader(ds, batch_size=batch_size, shuffle=False,
                      num_workers=num_workers, pin_memory=True,
                      collate_fn=_collate)


def main():
    args = parse_args()
    inputs_dir = Path(args.brainmvp_inputs_dir)
    out_dir    = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    parquet_path  = out_dir / "brainmvp_pooled.parquet"
    pt_path       = out_dir / "brainmvp_pooled.pt"
    manifest_path = out_dir / "manifest.csv"

    print("=" * 70)
    print(f"  BrainMVP frozen-encoder embedding extraction")
    print(f"  Inputs   : {inputs_dir}")
    print(f"  Out dir  : {out_dir}")
    print("=" * 70)

    df_all = discover_scans(args.matched_labels_csv, inputs_dir)
    if args.max_subjects:
        df_all = df_all.head(args.max_subjects)
    print(f"  Candidate scans      : {len(df_all)}")

    # Resume: skip scans already in parquet
    if parquet_path.exists() and not args.overwrite:
        existing = pd.read_parquet(parquet_path, columns=["bids_sub", "bids_ses"])
        existing["key"] = existing["bids_sub"] + "_" + existing["bids_ses"]
        df_all["key"] = df_all["bids_sub"] + "_" + df_all["bids_ses"]
        n_before = len(df_all)
        df_all = df_all[~df_all["key"].isin(existing["key"])]
        df_all = df_all.drop(columns="key").reset_index(drop=True)
        print(f"  Already in parquet   : {n_before - len(df_all)}  (skipping)")
    if len(df_all) == 0:
        print("  Nothing to extract. Done.")
        return
    print(f"  Scans to extract     : {len(df_all)}")

    # Build encoder (BrainMVPClassifier wraps uniformer_small + head;
    # we use only the encoder + pool, ignore the head)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = BrainMVPClassifier(n_classes=2, drop_rate=0.0)
    model = load_brainmvp_pretrained(model, args.pretrained_ckpt)
    model = model.to(device).eval()
    for p in model.parameters():
        p.requires_grad = False
    pool = nn.AdaptiveAvgPool3d(1).to(device)
    print(f"  Device   : {device}")

    loader = build_loader(df_all, args.batch_size, args.num_workers)
    embeddings_list, ids_list, manifest_rows = [], [], []
    n_done = 0

    with torch.no_grad():
        for batch in loader:
            x = batch["image"].to(device, non_blocking=True).float()
            features = model.encoder(x)
            x4 = features[-1]                                  # (B, 512, D', H', W')
            pooled = pool(x4).flatten(1).float().cpu().numpy() # (B, 512)
            for i in range(pooled.shape[0]):
                embeddings_list.append(pooled[i])
                ids_list.append((batch["sub"][i], batch["ses"][i]))
                manifest_rows.append({
                    "bids_sub":  batch["sub"][i],
                    "bids_ses":  batch["ses"][i],
                    "status":    "ok",
                    "source":    "uniformer.stage4_avg_pool",
                    "timestamp": datetime.datetime.now().isoformat(),
                })
            n_done += pooled.shape[0]
            if n_done % (10 * args.batch_size) == 0 or n_done >= len(df_all):
                print(f"  [{n_done:>4d}/{len(df_all)}] extracted")

    if not embeddings_list:
        print("  No embeddings produced. Aborting write.")
        return

    embeddings_arr = np.stack(embeddings_list, axis=0)
    print(f"  Embeddings shape: {embeddings_arr.shape}")

    new_df = pd.DataFrame({
        "bids_sub": [s for s, e in ids_list],
        "bids_ses": [e for s, e in ids_list],
        "status":   "ok",
    })
    for d in range(EMBED_DIM):
        new_df[f"dim_{d}"] = embeddings_arr[:, d]
    if parquet_path.exists() and not args.overwrite:
        old_df = pd.read_parquet(parquet_path)
        all_df = pd.concat([old_df, new_df], ignore_index=True)
    else:
        all_df = new_df
    all_df = all_df.sort_values(["bids_sub", "bids_ses"]).reset_index(drop=True)

    embed_cols = [f"dim_{d}" for d in range(EMBED_DIM)]
    embeddings_full = torch.from_numpy(all_df[embed_cols].to_numpy(dtype=np.float32))
    ids_full = list(zip(all_df["bids_sub"].tolist(), all_df["bids_ses"].tolist()))

    # Atomic writes
    tmp = parquet_path.with_suffix(".parquet.tmp")
    all_df.to_parquet(tmp, index=False); os.replace(tmp, parquet_path)
    tmp = pt_path.with_suffix(".pt.tmp")
    torch.save({
        "embeddings": embeddings_full,
        "ids":        ids_full,
        "embed_dim":  EMBED_DIM,
        "model_name": MODEL_NAME,
        "n_scans":    len(ids_full),
        "source_ckpt": str(args.pretrained_ckpt),
    }, tmp); os.replace(tmp, pt_path)

    new_manifest = pd.DataFrame(manifest_rows)
    if manifest_path.exists():
        all_manifest = pd.concat(
            [pd.read_csv(manifest_path), new_manifest], ignore_index=True)
    else:
        all_manifest = new_manifest
    all_manifest.to_csv(manifest_path, index=False)

    print("=" * 70)
    print(f"  Parquet : {parquet_path}  ({len(all_df)} rows total)")
    print(f"  PyTorch : {pt_path}        ({embeddings_full.shape[0]} embeddings, {EMBED_DIM}d)")
    print(f"  Manifest: {manifest_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()
