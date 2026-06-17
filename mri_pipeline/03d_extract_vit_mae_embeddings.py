"""
03d_extract_vit_mae_embeddings.py
==================================
One-time extraction of frozen ViT-MAE75 ViT-B/16 3D pooled-CLS
embeddings for every preprocessed (sub, ses) scan.

Output: paired Parquet + .pt keyed by (bids_sub, bids_ses) with the
row-aligned 768-d pooled CLS vector.

Forward (matches Vision_Transformer3D.forward without the head):
    x = patch_embed(volume)            # (B, n_patches, 768)
    x = concat([cls_token, x])
    x = x + pos_embed
    for block in blocks: x = block(x)
    x = norm(x)
    cls = x[:, 0]                      # (B, 768)

Implementation: swap model.head <- nn.Identity() so the existing
forward() returns the CLS feature directly. Avoids reimplementing the
forward chain.

Input: (128, 128, 128) z-scored NIfTI from vit_inputs/ (produced by
03_prepare_ViT.py). No augmentation -- deterministic eval transform.
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

THIS_DIR = Path(__file__).resolve().parent     # mri_pipeline/
REPO_ROOT = THIS_DIR.parent
sys.path.insert(0, str(THIS_DIR))
sys.path.insert(0, str(REPO_ROOT))

# Reuse ViT trainer's encoder builder + checkpoint loader
_vit_spec = importlib.util.spec_from_file_location(
    "_vit_pipeline", THIS_DIR / "04_supervised_finetuning_ViT.py")
_vit = importlib.util.module_from_spec(_vit_spec)
sys.modules["_vit_pipeline"] = _vit
_vit_spec.loader.exec_module(_vit)

load_pretrained_checkpoint = _vit.load_pretrained_checkpoint
load_matched_labels        = _vit.load_matched_labels
patient_id_to_bids_sub     = _vit.patient_id_to_bids_sub
from _vit_recipe.vit3d import Vision_Transformer3D
from bidsification.exclusions import is_excluded_subject

warnings.filterwarnings("ignore")

MODEL_NAME = "ViT-MAE75"
EMBED_DIM  = 768


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pretrained_ckpt",   type=str, required=True)
    p.add_argument("--vit_inputs_dir",    type=str, required=True)
    p.add_argument("--matched_labels_csv", type=str, required=True)
    p.add_argument("--out_dir",           type=str, required=True)
    p.add_argument("--batch_size",        type=int, default=4)
    p.add_argument("--num_workers",       type=int, default=4)
    p.add_argument("--max_subjects",      type=int, default=None)
    p.add_argument("--overwrite",         action="store_true")
    return p.parse_args()


def discover_scans(matched_csv: str, inputs_dir: Path) -> pd.DataFrame:
    matched_df = load_matched_labels(matched_csv)
    matched_df = matched_df[~matched_df["Patient_ID"].apply(
        is_excluded_subject)].copy()

    def _resolve_path(row):
        sub = patient_id_to_bids_sub(row["Patient_ID"])
        ses_label = f"ses-{row['bids_ses']}"
        return str(inputs_dir / sub / ses_label /
                   f"{sub}_{ses_label}_space-ViT128_desc-preproc_T1w.nii.gz")

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
    inputs_dir = Path(args.vit_inputs_dir)
    out_dir    = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    parquet_path  = out_dir / "vit_mae_pooled.parquet"
    pt_path       = out_dir / "vit_mae_pooled.pt"
    manifest_path = out_dir / "manifest.csv"

    print("=" * 70)
    print(f"  ViT-MAE75 frozen-encoder embedding extraction")
    print(f"  Inputs   : {inputs_dir}")
    print(f"  Out dir  : {out_dir}")
    print("=" * 70)

    df_all = discover_scans(args.matched_labels_csv, inputs_dir)
    if args.max_subjects:
        df_all = df_all.head(args.max_subjects)
    print(f"  Candidate scans      : {len(df_all)}")

    if parquet_path.exists() and not args.overwrite:
        existing = pd.read_parquet(parquet_path, columns=["bids_sub", "bids_ses"])
        existing["key"] = existing["bids_sub"] + "_" + existing["bids_ses"]
        df_all["key"] = df_all["bids_sub"] + "_" + df_all["bids_ses"]
        n_before = len(df_all)
        df_all = df_all[~df_all["key"].isin(existing["key"])]
        df_all = df_all.drop(columns="key").reset_index(drop=True)
        print(f"  Already in parquet   : {n_before - len(df_all)}")
    if len(df_all) == 0:
        print("  Nothing to extract. Done."); return
    print(f"  Scans to extract     : {len(df_all)}")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # Same architecture as in 04_supervised_finetuning_ViT.py line 985.
    # n_classes doesn't matter -- we replace the head with Identity anyway.
    model = Vision_Transformer3D(
        img_size=(128, 128, 128), patch_size=16, in_chans=1,
        n_classes=2, embed_dim=EMBED_DIM, depth=12, n_heads=12, mlp_ratio=4.0,
        qkv_bias=True, drop_path_rate=0.0, p=0.0, attn_p=0.0,
        global_avg_pool=False, pos_embed_type="learnable",
        patch_embed_fun="conv3d",
    )
    model = load_pretrained_checkpoint(model, args.pretrained_ckpt)
    # Swap head for Identity so model.forward returns the (B, 768) CLS feature.
    model.head = nn.Identity()
    model = model.to(device).eval()
    for p in model.parameters():
        p.requires_grad = False
    print(f"  Device   : {device}")

    loader = build_loader(df_all, args.batch_size, args.num_workers)
    embeddings_list, ids_list, manifest_rows = [], [], []
    n_done = 0

    with torch.no_grad():
        for batch in loader:
            x = batch["image"].to(device, non_blocking=True).float()
            cls = model(x)                                # (B, 768) since head=Identity
            cls_np = cls.float().cpu().numpy()
            for i in range(cls_np.shape[0]):
                embeddings_list.append(cls_np[i])
                ids_list.append((batch["sub"][i], batch["ses"][i]))
                manifest_rows.append({
                    "bids_sub":  batch["sub"][i],
                    "bids_ses":  batch["ses"][i],
                    "status":    "ok",
                    "source":    "vit3d.cls_token_final",
                    "timestamp": datetime.datetime.now().isoformat(),
                })
            n_done += cls_np.shape[0]
            if n_done % (10 * args.batch_size) == 0 or n_done >= len(df_all):
                print(f"  [{n_done:>4d}/{len(df_all)}] extracted "
                      f"(latest: {batch['sub'][-1]} ses-{batch['ses'][-1]})")

    if not embeddings_list:
        print("  No embeddings produced. Aborting write."); return

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

    # Atomic writes -- .pt FIRST (critical), parquet with CSV fallback.
    tmp = pt_path.with_suffix(".pt.tmp")
    torch.save({
        "embeddings": embeddings_full,
        "ids":        ids_full,
        "embed_dim":  EMBED_DIM,
        "model_name": MODEL_NAME,
        "n_scans":    len(ids_full),
        "source_ckpt": str(args.pretrained_ckpt),
    }, tmp); os.replace(tmp, pt_path)
    print(f"  Saved .pt: {pt_path}")
    try:
        tmp = parquet_path.with_suffix(".parquet.tmp")
        all_df.to_parquet(tmp, index=False); os.replace(tmp, parquet_path)
        print(f"  Saved parquet: {parquet_path}")
    except (ImportError, ValueError) as exc:
        csv_path = parquet_path.with_suffix(".csv.gz")
        print(f"  [WARN] parquet write failed ({exc.__class__.__name__}: "
              f"{exc}). Falling back to CSV: {csv_path}")
        all_df.to_csv(csv_path, index=False, compression="gzip")

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
