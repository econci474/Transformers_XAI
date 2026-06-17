"""
03_extract_braindino_embeddings.py
===================================
One-time extraction of frozen BrainDINO ViT-B/16 pooled-CLS embeddings
for every preprocessed (sub, ses) scan in the matched-labels CSV.

Output: paired Parquet (for inspection) + .pt (for fast trainer loads),
both keyed by (bids_sub, bids_ses), with the row-aligned 768-d pooled
CLS vector.

Why: the encoder forward (128 axial slices x 256^2 x ViT-B) is the
expensive part. When frozen + aug=none, the encoder output is fully
deterministic per scan. Cache once, then HP-sweep the head in seconds
per cell (see 04_head_finetune_from_embeddings.py + cached_head_sweep/).

Outputs (one shot, atomic):
  <out_dir>/braindino_pooled.parquet  - columns: bids_sub, bids_ses,
                                        status, dim_0..dim_767
  <out_dir>/braindino_pooled.pt       - dict: {embeddings (N, 768) f32,
                                        ids list[(sub, ses)],
                                        embed_dim 768, model_name 'BrainDINO',
                                        n_scans, source_ckpt}
  <out_dir>/manifest.csv              - per-scan status / source / timestamp

Resume: re-running skips (sub, ses) pairs already present in the parquet.
"""

import argparse
import datetime
import importlib.util
import json
import os
import sys
import warnings
from pathlib import Path

# torch first on Windows to dodge MKL DLL collision (see prep script).
import torch  # noqa: F401
import numpy as np
import pandas as pd
import monai.transforms as mt
from monai.data import Dataset as MonaiDataset
from torch.utils.data import DataLoader

THIS_DIR = Path(__file__).resolve().parent
MRI_DIR = THIS_DIR.parent
REPO_ROOT = MRI_DIR.parent
sys.path.insert(0, str(MRI_DIR))
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(THIS_DIR))

# Vendored DINOv3 + the BrainDINO classifier wrapper (we only need the
# encoder + checkpoint loader; the head is unused here).
from dinov3 import vit_base
_bdino_spec = importlib.util.spec_from_file_location(
    "_bdino_trainer", THIS_DIR / "02_supervised_finetuning_BrainDINO.py")
_bdino = importlib.util.module_from_spec(_bdino_spec)
sys.modules["_bdino_trainer"] = _bdino
_bdino_spec.loader.exec_module(_bdino)

BrainDINOClassifier      = _bdino.BrainDINOClassifier
load_braindino_pretrained = _bdino.load_braindino_pretrained
VIT_KWARGS               = _bdino.VIT_KWARGS
N_AXIAL                  = _bdino.N_AXIAL
EMBED_DIM                = _bdino.EMBED_DIM
SLICE_SIZE               = _bdino.SLICE_SIZE

# Shared ViT-pipeline utilities (load_matched_labels, patient_id_to_bids_sub,
# session_to_months, is_excluded_subject via _vit.is_excluded import chain).
_vit_spec = importlib.util.spec_from_file_location(
    "_vit_pipeline", MRI_DIR / "04_supervised_finetuning.py")
_vit = importlib.util.module_from_spec(_vit_spec)
sys.modules["_vit_pipeline"] = _vit
_vit_spec.loader.exec_module(_vit)

load_matched_labels      = _vit.load_matched_labels
patient_id_to_bids_sub   = _vit.patient_id_to_bids_sub
session_to_months        = _vit.session_to_months
from bidsification.exclusions import is_excluded_subject

warnings.filterwarnings("ignore")

MODEL_NAME = "BrainDINO"


# ── Args ──────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pretrained_ckpt", type=str, required=True,
                   help="Path to brain_dino_model.pth (DINO teacher dict).")
    p.add_argument("--braindino_inputs_dir", type=str, required=True,
                   help="Output dir of 01_prepare_braindino_inputs.py.")
    p.add_argument("--matched_labels_csv", type=str, required=True,
                   help="Per-(bids_sub, bids_ses) matched CSV (HPC "
                        "_extended_post_exclusion variant).")
    p.add_argument("--out_dir", type=str, required=True,
                   help="Output dir for braindino_pooled.{parquet,pt} + manifest.csv.")
    p.add_argument("--batch_size", type=int, default=4,
                   help="Volumes per forward pass (each volume = 128 axial "
                        "slices). batch=4 fits comfortably on A100-80GB at fp32.")
    p.add_argument("--num_workers", type=int, default=4)
    p.add_argument("--max_subjects", type=int, default=None,
                   help="Smoke test: cap to first N (sub, ses) pairs.")
    p.add_argument("--overwrite", action="store_true",
                   help="Re-extract scans already present in the parquet.")
    return p.parse_args()


# ── Discover scans (no train/val/test split -- embeddings are per-scan) ───────
def discover_scans(matched_csv: str, inputs_dir: Path) -> pd.DataFrame:
    """Returns df with bids_sub, bids_ses, image_path, for every scan in the
    matched CSV that (1) survives subject exclusion and (2) has a
    preprocessed NIfTI on disk."""
    matched_df = load_matched_labels(matched_csv)
    matched_df = matched_df[~matched_df["Patient_ID"].apply(
        is_excluded_subject)].copy()

    def _resolve_path(row):
        sub = patient_id_to_bids_sub(row["Patient_ID"])
        ses_label = f"ses-{row['bids_ses']}"
        return str(inputs_dir / sub / ses_label /
                   f"{sub}_{ses_label}_space-MNI128x256_desc-braindino_T1w.nii.gz")

    matched_df["bids_sub"] = matched_df["Patient_ID"].apply(patient_id_to_bids_sub)
    matched_df["image_path"] = matched_df.apply(_resolve_path, axis=1)
    have = matched_df["image_path"].apply(os.path.isfile)
    n_missing = int((~have).sum())
    print(f"  Scans with preprocessed input on disk: {int(have.sum())}  "
          f"(missing {n_missing})")
    return matched_df[have].reset_index(drop=True)[
        ["bids_sub", "bids_ses", "image_path"]]


# ── Dataset / dataloader (deterministic eval transform, no aug) ───────────────
def build_loader(df: pd.DataFrame, batch_size: int, num_workers: int) -> DataLoader:
    items = [{"image": row.image_path,
              "sub": row.bids_sub,
              "ses": row.bids_ses}
             for row in df.itertuples()]
    transform = mt.Compose([
        mt.LoadImaged(keys=["image"]),
        mt.EnsureChannelFirstd(keys=["image"]),
        mt.ToTensord(keys=["image"]),
    ])
    ds = MonaiDataset(data=items, transform=transform)
    # collate keys explicitly so the metadata strings pass through cleanly
    def _collate(batch):
        return {
            "image": torch.stack([b["image"] for b in batch], dim=0),
            "sub":   [b["sub"] for b in batch],
            "ses":   [b["ses"] for b in batch],
        }
    return DataLoader(ds, batch_size=batch_size, shuffle=False,
                      num_workers=num_workers, pin_memory=True,
                      collate_fn=_collate)


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    inputs_dir = Path(args.braindino_inputs_dir)
    out_dir    = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    parquet_path  = out_dir / "braindino_pooled.parquet"
    pt_path       = out_dir / "braindino_pooled.pt"
    manifest_path = out_dir / "manifest.csv"

    print("=" * 70)
    print(f"  BrainDINO embedding extraction")
    print(f"  Inputs   : {inputs_dir}")
    print(f"  Out dir  : {out_dir}")
    print(f"  Ckpt     : {args.pretrained_ckpt}")
    print("=" * 70)

    # ── Discover scans ────────────────────────────────────────────────────────
    df_all = discover_scans(args.matched_labels_csv, inputs_dir)
    if args.max_subjects:
        df_all = df_all.head(args.max_subjects)
    print(f"  Candidate scans      : {len(df_all)}")

    # ── Resume: drop scans already in the parquet ─────────────────────────────
    if parquet_path.exists() and not args.overwrite:
        existing = pd.read_parquet(parquet_path, columns=["bids_sub", "bids_ses"])
        existing["key"] = existing["bids_sub"] + "_" + existing["bids_ses"]
        df_all["key"] = df_all["bids_sub"] + "_" + df_all["bids_ses"]
        n_before = len(df_all)
        df_all = df_all[~df_all["key"].isin(existing["key"])]
        df_all = df_all.drop(columns="key").reset_index(drop=True)
        print(f"  Already in parquet   : {n_before - len(df_all)}  "
              f"(skipping; re-run with --overwrite to re-extract)")
    if len(df_all) == 0:
        print("  Nothing to extract. Done.")
        return

    print(f"  Scans to extract     : {len(df_all)}")

    # ── Build model: encoder only (head is unused for extraction) ─────────────
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = BrainDINOClassifier(n_classes=2, drop_rate=0.0)
    model = load_braindino_pretrained(model, args.pretrained_ckpt)
    model = model.to(device).eval()
    for p in model.parameters():
        p.requires_grad = False
    print(f"  Device   : {device}  "
          f"(encoder: {sum(p.numel() for p in model.encoder.parameters()):,} params)")

    # ── Extraction loop ───────────────────────────────────────────────────────
    loader = build_loader(df_all, batch_size=args.batch_size,
                          num_workers=args.num_workers)
    embeddings_list = []
    ids_list = []
    manifest_rows = []
    n_done = 0

    with torch.no_grad():
        for batch in loader:
            x = batch["image"].to(device, non_blocking=True).float()
            # Reuse BrainDINOClassifier.forward up to mean-pooled CLS (skip head).
            # x: (B, 1, 256, 256, 128) -> (B, S=128, 256, 256) -> (B*S, 1, 256, 256)
            B = x.shape[0]
            if x.shape[-1] != N_AXIAL:
                raise RuntimeError(
                    f"Expected (B, 1, 256, 256, {N_AXIAL}); got {tuple(x.shape)}")
            slices = x.squeeze(1).permute(0, 3, 2, 1).contiguous()
            slices = slices.reshape(B * N_AXIAL, 1, SLICE_SIZE, SLICE_SIZE)
            cls = model.encode_slices(slices)                  # (B*S, 768)
            cls = cls.view(B, N_AXIAL, EMBED_DIM)              # (B, S, 768)
            pooled = cls.mean(dim=1).float().cpu().numpy()     # (B, 768)

            for i in range(B):
                embeddings_list.append(pooled[i])
                ids_list.append((batch["sub"][i], batch["ses"][i]))
                manifest_rows.append({
                    "bids_sub":  batch["sub"][i],
                    "bids_ses":  batch["ses"][i],
                    "status":    "ok",
                    "source":    "encoder",
                    "timestamp": datetime.datetime.now().isoformat(),
                })

            n_done += B
            if n_done % (10 * args.batch_size) == 0 or n_done >= len(df_all):
                print(f"  [{n_done:>4d}/{len(df_all)}] extracted "
                      f"(latest: {batch['sub'][-1]} ses-{batch['ses'][-1]})")

    if not embeddings_list:
        print("  No embeddings produced. Aborting write.")
        return

    embeddings_arr = np.stack(embeddings_list, axis=0)        # (N, 768) float32
    print(f"  Embeddings shape: {embeddings_arr.shape}  dtype={embeddings_arr.dtype}")

    # ── Merge with any existing parquet for resumed runs ──────────────────────
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

    # Build the row-aligned .pt from the full merged set
    embed_cols = [f"dim_{d}" for d in range(EMBED_DIM)]
    embeddings_full = torch.from_numpy(
        all_df[embed_cols].to_numpy(dtype=np.float32))
    ids_full = list(zip(all_df["bids_sub"].tolist(), all_df["bids_ses"].tolist()))

    # ── Atomic writes -- .pt FIRST (the critical artifact for the head ──────
    # trainer). Parquet is nice-to-have for inspection; if pyarrow isn't
    # installed we fall back to CSV so the run never gets wasted again.
    tmp_pt = pt_path.with_suffix(".pt.tmp")
    torch.save({
        "embeddings": embeddings_full,                 # (N, D) float32
        "ids":        ids_full,                        # list[(sub, ses)]
        "embed_dim":  EMBED_DIM,
        "model_name": MODEL_NAME,
        "n_scans":    len(ids_full),
        "source_ckpt": str(args.pretrained_ckpt),
    }, tmp_pt)
    os.replace(tmp_pt, pt_path)
    print(f"  Saved .pt: {pt_path}")

    try:
        tmp_pq = parquet_path.with_suffix(".parquet.tmp")
        all_df.to_parquet(tmp_pq, index=False)
        os.replace(tmp_pq, parquet_path)
        print(f"  Saved parquet: {parquet_path}")
    except (ImportError, ValueError) as exc:
        csv_path = parquet_path.with_suffix(".csv.gz")
        print(f"  [WARN] parquet write failed ({exc.__class__.__name__}: "
              f"{exc}). Falling back to CSV: {csv_path}")
        all_df.to_csv(csv_path, index=False, compression="gzip")

    # ── Manifest (append) ─────────────────────────────────────────────────────
    new_manifest = pd.DataFrame(manifest_rows)
    if manifest_path.exists():
        old_manifest = pd.read_csv(manifest_path)
        all_manifest = pd.concat([old_manifest, new_manifest], ignore_index=True)
    else:
        all_manifest = new_manifest
    all_manifest.to_csv(manifest_path, index=False)

    print("=" * 70)
    print(f"  Parquet : {parquet_path}  ({len(all_df)} rows total)")
    print(f"  PyTorch : {pt_path}        ({embeddings_full.shape[0]} embeddings)")
    print(f"  Manifest: {manifest_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()
