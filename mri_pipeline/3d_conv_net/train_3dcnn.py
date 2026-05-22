"""
train_3dcnn.py
==============
Train a Spasov-style 3D CNN baseline (vanilla `MRI3DCNN` or separable
`MRI3DSeparableCNN`) on ADNI MRI for the diagnostic classification tasks.

Tasks (same definitions as the ViT pipeline, by-visit diagnosis, all sessions)
-----------------------------------------------------------------------------
  T1_binary     CN vs MCI+AD        (binary)
  T1b_binary    CN+MCI vs AD        (binary)
  T1c_binary    CN vs AD            (binary, MCI scans dropped)
  T2_multiclass CN / MCI / AD       (3-class)

Every scan is labelled by the diagnosis recorded at the time of that scan
(baseline scans use `Label_bl_multi`, follow-ups use `Label_visit_diag`).
A non-foundation, non-transformer baseline for the thesis comparison.

Inputs
------
  --cnn_inputs_dir     Per-(sub, ses) z-scored MNI volumes from
                       00_prepare_CNN_inputs.py (193x229x193, 1 mm, RAS).
  --matched_labels_csv MRI<->clinical matched labels CSV (per (bids_sub,
                       bids_ses) with Label_bl_multi / Label_visit_diag).
  --data_dir           Post-exclusion clinical splits root:
                       {data_dir}/seed_{seed}/{train,val,test}.csv. Subject
                       partitioning only (subject-level — all visits of a
                       subject stay in one split). Reused from the clinical
                       pipeline so test subjects stay consistent across
                       modalities.

Model
-----
  --model vanilla   -> MRI3DCNN          (3DCNN.py,  ~0.5 M params)
  --model separable -> MRI3DSeparableCNN (3DSCNN.py, ~0.03 M params)
Loaded via importlib (filenames start with a digit). The classifier head is
sized `n_outputs = num_labels` (2 or 3) and trained with CrossEntropyLoss
(class-weighted for the CN/MCI/AD imbalance) — uniform across binary and
3-class tasks, matching 04_supervised_finetuning_ViT.py.

Method (mirrors mri_pipeline/04_supervised_finetuning_ViT.py)
-------------------------------------------------------------
- Best checkpoint / early stopping on **val balanced accuracy** (tie-break on
  lower val_loss).
- Resumable per-epoch checkpoint (survives the 12 h SLURM/Colab cap).
- Per-subject evaluation (mean predicted probability per Patient_ID).
- metrics.json / train_log.csv schema matches the ViT pipeline.

Output layout
-------------
  {out_dir}/Spasov3DCNN_{model}/{task}/seed_{seed}/
    best_model.pt  last_checkpoint.pt  train_log.csv
    test_predictions.csv  dataset_manifest.csv  metrics.json

Smoke test
----------
  python train_3dcnn.py --model vanilla --task T1c_binary --seed 0 \
      --cnn_inputs_dir ... --matched_labels_csv ... --data_dir ... \
      --out_dir ... --max_subjects 30 --epochs 3
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import re
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.metrics import (
    accuracy_score, average_precision_score, balanced_accuracy_score,
    confusion_matrix, f1_score, precision_recall_fscore_support,
    precision_score, recall_score, roc_auc_score,
)
from sklearn.utils.class_weight import compute_class_weight
from torch.utils.data import DataLoader

import monai.transforms as mt
from monai.data import Dataset

warnings.filterwarnings("ignore")

THIS_DIR = Path(__file__).resolve().parent          # mri_pipeline/3d_conv_net
REPO_ROOT = THIS_DIR.parent.parent
sys.path.insert(0, str(REPO_ROOT))
from bidsification.exclusions import is_excluded_subject

# ── Task definitions (the 4 diagnostic tasks; same as the ViT pipeline) ───────
# By-visit diagnosis: bl rows use Label_bl_multi, follow-ups use Label_visit_diag.
# A diagnosis absent from label_map (e.g. MCI under T1c) -> the scan is dropped.
TASK_CONFIG = {
    "T1_binary": {
        "label_map":   {"CN": 0, "MCI": 1, "AD": 1},
        "num_labels":  2,
        "task_type":   "binary",
        "description": "Binary: CN vs MCI+AD",
    },
    "T1b_binary": {
        "label_map":   {"CN": 0, "MCI": 0, "AD": 1},
        "num_labels":  2,
        "task_type":   "binary",
        "description": "Binary: CN+MCI vs AD",
    },
    "T1c_binary": {
        "label_map":   {"CN": 0, "AD": 1},          # MCI dropped
        "num_labels":  2,
        "task_type":   "binary",
        "description": "Binary: CN vs AD (MCI excluded)",
    },
    "T2_multiclass": {
        "label_map":   {"CN": 0, "MCI": 1, "AD": 2},
        "num_labels":  3,
        "task_type":   "multiclass",
        "description": "Multiclass: CN / MCI / AD",
    },
}


def class_names_for(task_cfg: dict) -> list:
    """Per-index human-readable class names. label_map may be many-to-one
    (e.g. T1_binary CN->0, MCI->1, AD->1) -> 'MCI+AD'."""
    n, lm = task_cfg["num_labels"], task_cfg["label_map"]
    names = []
    for i in range(n):
        members = [k for k, v in lm.items() if v == i]
        names.append("+".join(members) if members else f"class_{i}")
    return names


# ── Model loading (filenames start with a digit -> importlib) ─────────────────
def load_model_class(model_kind: str):
    """Import MRI3DCNN / MRI3DSeparableCNN from 3DCNN.py / 3DSCNN.py by path."""
    fname, cls_name = {
        "vanilla":   ("3DCNN.py",  "MRI3DCNN"),
        "separable": ("3DSCNN.py", "MRI3DSeparableCNN"),
    }[model_kind]
    spec = importlib.util.spec_from_file_location(
        f"_spasov_{model_kind}", THIS_DIR / fname)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return getattr(mod, cls_name)


# ── Label / split / path helpers ──────────────────────────────────────────────
def bids_sub_to_ptid(bids_sub: str):
    """'sub-002S0413' or '002S0413' -> '002_S_0413'."""
    s = bids_sub[4:] if bids_sub.startswith("sub-") else bids_sub
    m = re.match(r"^(\d+)S(\d+)$", s)
    return f"{m.group(1)}_S_{m.group(2)}" if m else None


def cnn_image_path(patient_id: str, viscode: str, cnn_dir: Path) -> str:
    sub = "sub-" + patient_id.replace("_", "")
    ses = f"ses-{viscode}"
    return str(cnn_dir / sub / ses /
               f"{sub}_{ses}_space-CNN193_desc-preproc_T1w.nii.gz")


def load_matched_labels(csv_path: str) -> pd.DataFrame:
    """Load the matched-labels master; ensure a Patient_ID column exists."""
    df = pd.read_csv(csv_path, dtype=str)
    if "Patient_ID" not in df.columns:
        df["Patient_ID"] = df["bids_sub"].apply(bids_sub_to_ptid)
    return df.dropna(subset=["Patient_ID"])


def load_split(split_csv: Path, matched_df: pd.DataFrame, cnn_dir: Path,
               task_cfg: dict) -> tuple[pd.DataFrame, int]:
    """Build one split's labelled scan list: intersect the subject-level split
    CSV with the matched-labels master, resolve the by-visit label via the
    task's label_map (scans whose diagnosis is unmapped are dropped), and
    resolve the cnn_inputs NIfTI path. Returns (df, n_missing)."""
    subjects = pd.read_csv(split_csv, usecols=["Patient_ID"])["Patient_ID"].unique()
    df = matched_df[matched_df["Patient_ID"].isin(subjects)].copy()

    # Subject-level exclusions (defensive — post-exclusion splits already drop
    # these, but a non-post-exclusion split CSV would not).
    df = df[~df["Patient_ID"].apply(
        lambda p: is_excluded_subject(p, include_diagnostic_reversion=True))]

    # By-visit label: bl rows use Label_bl_multi, follow-ups use Label_visit_diag.
    df = df.dropna(subset=["Label_bl_multi", "Label_visit_diag"], how="any").copy()
    label_map = task_cfg["label_map"]

    def _label(row):
        dx = (row["Label_bl_multi"] if row["bids_ses"] == "bl"
              else row["Label_visit_diag"])
        return label_map.get(dx)        # unmapped diagnosis -> None -> dropped

    df["label"] = df.apply(_label, axis=1)
    df = df.dropna(subset=["label"]).copy()
    df["label"] = df["label"].astype(int)

    df["image_path"] = df.apply(
        lambda r: cnn_image_path(r["Patient_ID"], r["bids_ses"], cnn_dir), axis=1)
    have = df["image_path"].apply(os.path.isfile)
    n_missing = int((~have).sum())
    df = df[have].reset_index(drop=True)
    return df[["Patient_ID", "bids_ses", "label", "image_path"]], n_missing


# ── Datasets ──────────────────────────────────────────────────────────────────
def _train_transform() -> mt.Compose:
    """Load + light augmentation (random flips) for training."""
    return mt.Compose([
        mt.LoadImaged(keys=["image"]),
        mt.EnsureChannelFirstd(keys=["image"]),
        mt.RandFlipd(keys=["image"], prob=0.5, spatial_axis=0),
        mt.RandFlipd(keys=["image"], prob=0.5, spatial_axis=1),
        mt.RandFlipd(keys=["image"], prob=0.5, spatial_axis=2),
        mt.ToTensord(keys=["image"]),
    ])


def _eval_transform() -> mt.Compose:
    """Load only (no augmentation) for val / test."""
    return mt.Compose([
        mt.LoadImaged(keys=["image"]),
        mt.EnsureChannelFirstd(keys=["image"]),
        mt.ToTensord(keys=["image"]),
    ])


def build_dataset(df: pd.DataFrame, train: bool) -> Dataset:
    items = [{"image": r.image_path, "label": int(r.label)}
             for r in df.itertuples()]
    return Dataset(data=items,
                   transform=_train_transform() if train else _eval_transform())


# ── Metrics (mirror 04_supervised_finetuning_ViT.py) ──────────────────────────
def compute_test_metrics(y_true: np.ndarray, logits: np.ndarray, task_type: str):
    """Returns (scalar_metrics_dict, preds, probs). logits: [N, num_labels]."""
    probs = torch.softmax(torch.from_numpy(np.asarray(logits, dtype=np.float64)),
                          dim=-1).numpy()
    preds = probs.argmax(axis=-1)
    y_true = np.asarray(y_true).astype(int)
    multi = len(np.unique(y_true)) > 1

    out = {
        "accuracy":     float(accuracy_score(y_true, preds)),
        "balanced_acc": float(balanced_accuracy_score(y_true, preds)),
    }
    if task_type == "binary":
        positive = probs[:, 1]
        out.update({
            "precision":   float(precision_score(y_true, preds, zero_division=0)),
            "recall":      float(recall_score(y_true, preds, zero_division=0)),
            "sensitivity": float(recall_score(y_true, preds, pos_label=1, zero_division=0)),
            "specificity": float(recall_score(y_true, preds, pos_label=0, zero_division=0)),
            "f1":          float(f1_score(y_true, preds, zero_division=0)),
            "auc_roc":     float(roc_auc_score(y_true, positive)) if multi else float("nan"),
            "auc_pr":      float(average_precision_score(y_true, positive)) if multi else float("nan"),
        })
    else:
        out.update({
            "precision_macro": float(precision_score(y_true, preds, average="macro", zero_division=0)),
            "recall_macro":    float(recall_score(y_true, preds, average="macro", zero_division=0)),
            "macro_f1":        float(f1_score(y_true, preds, average="macro", zero_division=0)),
            "auc_roc_ovr":     float(roc_auc_score(y_true, probs, multi_class="ovr", average="macro")) if multi else float("nan"),
            "auc_pr_macro":    float(average_precision_score(
                np.eye(probs.shape[1])[y_true], probs, average="macro")) if multi else float("nan"),
        })
    metrics = {k: (round(v, 4) if isinstance(v, float) and not np.isnan(v) else v)
               for k, v in out.items()}
    return metrics, preds, probs


def compute_diagnostics(y_true: np.ndarray, preds: np.ndarray,
                        num_labels: int, label_names: list) -> dict:
    """Confusion matrix + per-class P/R/F1/support. String keys (wandb-safe)."""
    labels = list(range(num_labels))
    cm = confusion_matrix(y_true, preds, labels=labels)
    pr, rc, f1, sup = precision_recall_fscore_support(
        y_true, preds, labels=labels, zero_division=0)
    total = int(cm.sum())

    def _specificity(i):
        tp = int(cm[i, i])
        fp = int(cm[:, i].sum()) - tp
        fn = int(cm[i, :].sum()) - tp
        tn = total - tp - fp - fn
        return tn / (tn + fp) if (tn + fp) > 0 else 0.0

    per_class = {
        label_names[i]: {
            "precision":   round(float(pr[i]), 4),
            "recall":      round(float(rc[i]), 4),
            "sensitivity": round(float(rc[i]), 4),
            "specificity": round(_specificity(i), 4),
            "f1":          round(float(f1[i]), 4),
            "support":     int(sup[i]),
        }
        for i in labels
    }
    return {"labels": {str(i): label_names[i] for i in labels},
            "confusion_matrix": cm.tolist(),     # rows = true, cols = pred
            "per_class": per_class}


# ── Train / eval loop ─────────────────────────────────────────────────────────
def run_one_epoch(model, loader, criterion, optimizer, scaler, device,
                  train: bool):
    """One pass over `loader`. Returns (mean_loss, accuracy, logits, labels)."""
    model.train(mode=train)
    total_loss = total_correct = total_n = 0.0
    all_logits, all_labels = [], []

    ctx = torch.enable_grad() if train else torch.no_grad()
    with ctx:
        for batch in loader:
            x = batch["image"].to(device, non_blocking=True).float()
            y = batch["label"].to(device, non_blocking=True).long()   # [B]

            with torch.cuda.amp.autocast(enabled=(device.type == "cuda")):
                logits = model(x)                       # [B, num_labels]
                loss = criterion(logits, y)

            if train:
                optimizer.zero_grad(set_to_none=True)
                scaler.scale(loss).backward()
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                scaler.step(optimizer)
                scaler.update()

            preds = logits.detach().argmax(dim=-1)
            total_loss += loss.item() * x.size(0)
            total_correct += (preds == y).sum().item()
            total_n += x.size(0)
            all_logits.append(logits.detach().float().cpu().numpy())
            all_labels.append(y.detach().cpu().numpy())

    return (total_loss / max(total_n, 1),
            total_correct / max(total_n, 1),
            np.concatenate(all_logits, axis=0) if all_logits else np.empty((0,)),
            np.concatenate(all_labels) if all_labels else np.empty((0,)))


def lr_at(epoch_idx: int, args) -> float:
    """Linear warmup then cosine decay to a small floor."""
    eta_min = 1e-6
    if epoch_idx < args.warmup_epochs:
        return args.lr * (epoch_idx + 1) / max(args.warmup_epochs, 1)
    progress = ((epoch_idx - args.warmup_epochs)
                / max(args.epochs - args.warmup_epochs, 1))
    return eta_min + 0.5 * (args.lr - eta_min) * (1 + np.cos(np.pi * progress))


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--model", required=True, choices=["vanilla", "separable"])
    p.add_argument("--task", required=True, choices=list(TASK_CONFIG.keys()))
    p.add_argument("--seed", type=int, required=True, choices=[0, 1, 2])
    p.add_argument("--cnn_inputs_dir", required=True,
                   help="Per-(sub,ses) volumes from 00_prepare_CNN_inputs.py.")
    p.add_argument("--matched_labels_csv", required=True,
                   help="MRI<->clinical matched labels CSV.")
    p.add_argument("--data_dir", required=True,
                   help="Post-exclusion splits root ({data_dir}/seed_N/*.csv).")
    p.add_argument("--out_dir", required=True)
    p.add_argument("--epochs", type=int, default=60)
    p.add_argument("--batch_size", type=int, default=8)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--weight_decay", type=float, default=1e-4)
    p.add_argument("--warmup_epochs", type=int, default=5)
    p.add_argument("--patience", type=int, default=15,
                   help="Early-stopping patience on val balanced accuracy.")
    p.add_argument("--dropout", type=float, default=0.1)
    p.add_argument("--label_smoothing", type=float, default=0.0)
    p.add_argument("--num_workers", type=int, default=2)
    p.add_argument("--max_subjects", type=int, default=None,
                   help="Cap rows per split (smoke test).")
    p.add_argument("--no_resume", action="store_true",
                   help="Ignore any last_checkpoint.pt and start fresh.")
    p.add_argument("--wandb", action="store_true",
                   help="Enable Weights & Biases logging (no-op if absent).")
    p.add_argument("--wandb_project", type=str, default="cnn3d-mri")
    return p.parse_args()


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    task_cfg = TASK_CONFIG[args.task]
    num_labels = task_cfg["num_labels"]
    label_names = class_names_for(task_cfg)

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model_slug = f"Spasov3DCNN_{args.model}"
    out_dir = Path(args.out_dir) / model_slug / args.task / f"seed_{args.seed}"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print(f"  train_3dcnn — {args.model} | {args.task} | seed={args.seed}")
    print(f"  Device: {device}   Output: {out_dir}")
    print("=" * 70)

    # ── Load splits ───────────────────────────────────────────────────────────
    if not os.path.isfile(args.matched_labels_csv):
        raise FileNotFoundError(f"--matched_labels_csv not found: {args.matched_labels_csv}")
    matched_df = load_matched_labels(args.matched_labels_csv)
    print(f"  Matched labels CSV: {len(matched_df)} (subject, visit) rows")

    seed_dir = Path(args.data_dir) / f"seed_{args.seed}"
    cnn_dir = Path(args.cnn_inputs_dir)
    train_df, n_miss_tr = load_split(seed_dir / "train.csv", matched_df, cnn_dir, task_cfg)
    val_df,   n_miss_va = load_split(seed_dir / "val.csv",   matched_df, cnn_dir, task_cfg)
    test_df,  n_miss_te = load_split(seed_dir / "test.csv",  matched_df, cnn_dir, task_cfg)

    if args.max_subjects:
        train_df = train_df.head(args.max_subjects)
        val_df = val_df.head(min(args.max_subjects, len(val_df)))
        test_df = test_df.head(min(args.max_subjects, len(test_df)))

    print(f"  Splits — train: {len(train_df)} ({n_miss_tr} NIfTI missing)  "
          f"val: {len(val_df)} ({n_miss_va})  test: {len(test_df)} ({n_miss_te})")
    print(f"           subjects — train: {train_df['Patient_ID'].nunique()}, "
          f"val: {val_df['Patient_ID'].nunique()}, "
          f"test: {test_df['Patient_ID'].nunique()}")
    for name, df in [("train", train_df), ("val", val_df), ("test", test_df)]:
        print(f"           {name} class balance: {df['label'].value_counts().sort_index().to_dict()}")
    if min(len(train_df), len(val_df), len(test_df)) == 0:
        raise RuntimeError("Empty split — check cnn_inputs_dir / matched_labels_csv.")
    if train_df["label"].nunique() < 2 or val_df["label"].nunique() < 2:
        raise RuntimeError("train/val split has <2 classes — cannot train.")

    # Dataset manifest.
    manifest = pd.concat(
        [d.assign(split=s) for s, d in
         [("train", train_df), ("val", val_df), ("test", test_df)]],
        ignore_index=True)
    manifest[["split", "Patient_ID", "bids_ses", "label", "image_path"]].to_csv(
        out_dir / "dataset_manifest.csv", index=False)

    # ── Dataloaders ───────────────────────────────────────────────────────────
    train_loader = DataLoader(build_dataset(train_df, train=True),
                              batch_size=args.batch_size, shuffle=True,
                              num_workers=args.num_workers, pin_memory=True,
                              drop_last=True)          # BatchNorm needs B > 1
    val_loader = DataLoader(build_dataset(val_df, train=False), batch_size=1,
                            shuffle=False, num_workers=args.num_workers,
                            pin_memory=True)
    test_loader = DataLoader(build_dataset(test_df, train=False), batch_size=1,
                             shuffle=False, num_workers=args.num_workers,
                             pin_memory=True)

    # ── Model / loss / optimizer ──────────────────────────────────────────────
    ModelCls = load_model_class(args.model)
    model = ModelCls(in_channels=1, n_outputs=num_labels,
                     dropout=args.dropout).to(device)
    n_params = sum(p.numel() for p in model.parameters())
    print(f"  Model: {ModelCls.__name__}  n_outputs={num_labels}  "
          f"({n_params:,} parameters)")

    # Balanced class weights from the train split (padded to num_labels so the
    # weight tensor always matches the head — a small split may miss a class).
    classes_present = np.unique(train_df["label"].values).astype(int)
    cw_present = compute_class_weight("balanced", classes=classes_present,
                                      y=train_df["label"].values)
    cw = np.ones(num_labels, dtype=np.float32)
    for c, w in zip(classes_present, cw_present):
        cw[int(c)] = float(w)
    class_weights = torch.tensor(cw, dtype=torch.float, device=device)
    print(f"  Class weights: {dict(enumerate(cw.round(3).tolist()))}")
    criterion = nn.CrossEntropyLoss(weight=class_weights,
                                    label_smoothing=args.label_smoothing)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr,
                                  weight_decay=args.weight_decay)
    scaler = torch.cuda.amp.GradScaler(enabled=(device.type == "cuda"))

    # ── wandb (gated; no-op without --wandb) ──────────────────────────────────
    wb = None
    if args.wandb:
        try:
            import wandb
            wb = wandb.init(project=args.wandb_project,
                            name=f"{model_slug}-{args.task}-s{args.seed}",
                            config={"model": args.model, "seed": args.seed,
                                    "task": args.task, "num_labels": num_labels,
                                    "epochs": args.epochs, "lr": args.lr,
                                    "batch_size": args.batch_size,
                                    "n_params": n_params}, reinit=True)
        except ImportError:
            print("  [WARN] --wandb set but wandb not installed; local-only.")

    # ── Resumable checkpoint state ────────────────────────────────────────────
    ckpt_path = out_dir / "last_checkpoint.pt"
    best_metric = -1.0                  # best val balanced accuracy
    best_val_loss_at_best = float("inf")
    best_epoch = -1
    epochs_since_improve = 0
    log_rows: list = []
    start_epoch = 0

    if (not args.no_resume) and ckpt_path.exists():
        ck = torch.load(ckpt_path, map_location=device)
        model.load_state_dict(ck["model"])
        optimizer.load_state_dict(ck["optimizer"])
        scaler.load_state_dict(ck["scaler"])
        best_metric = ck["best_metric"]
        best_val_loss_at_best = ck["best_val_loss_at_best"]
        best_epoch = ck["best_epoch"]
        epochs_since_improve = ck["epochs_since_improve"]
        log_rows = ck["log_rows"]
        start_epoch = ck["epoch"] + 1
        print(f"  Resuming from epoch {start_epoch + 1} "
              f"(best val balanced acc={best_metric:.4f}).")

    def _save_checkpoint(ep: int):
        tmp = ckpt_path.with_suffix(".pt.tmp")
        torch.save({"model": model.state_dict(),
                    "optimizer": optimizer.state_dict(),
                    "scaler": scaler.state_dict(),
                    "epoch": ep, "best_metric": best_metric,
                    "best_val_loss_at_best": best_val_loss_at_best,
                    "best_epoch": best_epoch,
                    "epochs_since_improve": epochs_since_improve,
                    "log_rows": log_rows, "args": vars(args)}, tmp)
        os.replace(tmp, ckpt_path)

    train_log_path = out_dir / "train_log.csv"

    try:
        # ── Training loop ─────────────────────────────────────────────────────
        for epoch in range(start_epoch, args.epochs):
            cur_lr = lr_at(epoch, args)
            for g in optimizer.param_groups:
                g["lr"] = cur_lr

            tr_loss, tr_acc, _, _ = run_one_epoch(
                model, train_loader, criterion, optimizer, scaler, device, True)
            va_loss, va_acc, va_logits, va_labels = run_one_epoch(
                model, val_loader, criterion, optimizer, scaler, device, False)
            va_bacc = float(balanced_accuracy_score(
                va_labels.astype(int), va_logits.argmax(axis=-1)))

            log_rows.append({"epoch": epoch + 1, "lr": cur_lr,
                             "train_loss": tr_loss, "train_acc": tr_acc,
                             "val_loss": va_loss, "val_acc": va_acc,
                             "val_bacc": va_bacc})
            print(f"  [epoch {epoch+1:>3}/{args.epochs}] lr={cur_lr:.2e}  "
                  f"train_loss={tr_loss:.4f} train_acc={tr_acc:.4f}  "
                  f"val_loss={va_loss:.4f} val_acc={va_acc:.4f} "
                  f"val_bacc={va_bacc:.4f}")

            # Select on val balanced accuracy; tie-break on lower val_loss.
            improved = (va_bacc > best_metric + 1e-6) or (
                abs(va_bacc - best_metric) <= 1e-6
                and va_loss < best_val_loss_at_best)
            if improved:
                best_metric = va_bacc
                best_val_loss_at_best = va_loss
                best_epoch = epoch + 1
                epochs_since_improve = 0
                torch.save({"net": model.state_dict(), "epoch": best_epoch},
                           out_dir / "best_model.pt")
            else:
                epochs_since_improve += 1

            pd.DataFrame(log_rows).to_csv(train_log_path, index=False)
            _save_checkpoint(epoch)

            if wb is not None:
                import wandb
                wandb.log({"epoch": epoch + 1, "lr": cur_lr,
                           "train_loss": tr_loss, "train_acc": tr_acc,
                           "val_loss": va_loss, "val_acc": va_acc,
                           "val_bacc": va_bacc,
                           "best_val_balanced_acc": best_metric})

            if epochs_since_improve >= args.patience:
                print(f"  Early stopping at epoch {epoch+1} "
                      f"(no val balanced-acc improvement for {args.patience}).")
                break

        pd.DataFrame(log_rows).to_csv(train_log_path, index=False)

        # ── Final test on the best weights ────────────────────────────────────
        best_state = torch.load(out_dir / "best_model.pt", map_location=device)
        model.load_state_dict(best_state["net"])
        _, _, test_logits, test_labels = run_one_epoch(
            model, test_loader, criterion, optimizer, scaler, device, False)
        test_metrics, test_preds, test_probs = compute_test_metrics(
            test_labels, test_logits, task_cfg["task_type"])
        print(f"  Test metrics (image-level): {test_metrics}")
        test_diagnostics = compute_diagnostics(
            test_labels.astype(int), test_preds, num_labels, label_names)

        # Per-sample predictions (test_loader: shuffle=False, batch_size=1, so
        # row i aligns with test_df row i).
        n_te = len(test_labels)
        pred_df = pd.DataFrame({
            "Patient_ID": test_df["Patient_ID"].values[:n_te],
            "bids_ses":   test_df["bids_ses"].values[:n_te],
            "y_true": test_labels.astype(int),
            "y_pred": test_preds.astype(int),
        })
        for c in range(num_labels):
            pred_df[f"prob_{c}"] = test_probs[:, c]
        pred_df.to_csv(out_dir / "test_predictions.csv", index=False)

        # ── Per-subject evaluation (mean class probability per Patient_ID) ────
        prob_cols = [f"prob_{c}" for c in range(num_labels)]
        subj = pred_df.groupby("Patient_ID").agg(
            {**{pc: "mean" for pc in prob_cols}, "y_true": "first"})
        subj_probs = subj[prob_cols].to_numpy()
        # compute_test_metrics softmaxes its input; log(probs) round-trips back
        # to the (already mean-aggregated) probabilities.
        subj_logits = np.log(np.clip(subj_probs, 1e-8, 1.0))
        test_metrics_subject, _, _ = compute_test_metrics(
            subj["y_true"].to_numpy().astype(int), subj_logits,
            task_cfg["task_type"])
        print(f"  Test metrics (subject-level, n={len(subj)}): "
              f"{test_metrics_subject}")

        # ── metrics.json (schema matches the ViT pipeline) ────────────────────
        config = {
            "model_id":              model_slug,
            "model_kind":            args.model,
            "task":                  args.task,
            "task_description":      task_cfg["description"],
            "num_labels":            num_labels,
            "seed":                  args.seed,
            "n_params":              int(n_params),
            "epochs":                args.epochs,
            "best_epoch":            best_epoch,
            "best_val_balanced_acc": round(float(best_metric), 4),
            "lr":                    args.lr,
            "weight_decay":          args.weight_decay,
            "batch_size":            args.batch_size,
            "warmup_epochs":         args.warmup_epochs,
            "patience":              args.patience,
            "dropout":               args.dropout,
            "label_smoothing":       args.label_smoothing,
            "class_weights":         {i: round(float(w), 4) for i, w in enumerate(cw.tolist())},
            "n_train":               int(len(train_df)),
            "n_val":                 int(len(val_df)),
            "n_test":                int(len(test_df)),
            "n_test_subjects":       int(len(subj)),
            "timestamp":             datetime.now().isoformat(),
        }
        with open(out_dir / "metrics.json", "w") as f:
            json.dump({"config": config, "test_metrics": test_metrics,
                       "test_metrics_subject": test_metrics_subject,
                       "test_diagnostics": test_diagnostics}, f, indent=2)

        if wb is not None:
            import wandb
            wandb.log({f"test/{k}": v for k, v in test_metrics.items()
                       if isinstance(v, (int, float))})
            wandb.log({f"test_subject/{k}": v
                       for k, v in test_metrics_subject.items()
                       if isinstance(v, (int, float))})
            wandb.run.summary["best_epoch"] = best_epoch
            wandb.run.summary["best_val_balanced_acc"] = best_metric

        print(f"  Saved: {out_dir}/metrics.json")
        print("=" * 70)
    finally:
        if wb is not None:
            import wandb
            wandb.finish()


if __name__ == "__main__":
    main()
