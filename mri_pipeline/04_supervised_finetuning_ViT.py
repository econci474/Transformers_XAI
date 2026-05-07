"""
04_supervised_finetuning_ViT.py
================================
Supervised fine-tuning of the MAE-pretrained ViT-B/3D
(qasymjomart/ViT_recipe_for_AD) on ADNI baseline T1w MRIs preprocessed by
03_prepare_ViT.py (128 x 128 x 128 @ 1.75 mm RAS, z-scored on nonzero voxels).

Tasks (mirroring clinical_pipeline/03_encoder_finetune.py)
----------------------------------------------------------
  T1_binary     : CN vs MCI+AD          (binary)
  T2_multiclass : CN / MCI / AD         (3-class)
  T3a_conv3y    : Conversion to AD <=3y (binary, non-AD at baseline)
  T3b_conv5y    : Conversion to AD <=5y (binary, non-AD at baseline)

Strategies
----------
  full_ft : all parameters trained (lr=1e-4)
  frozen  : encoder frozen, classification head only (lr=1e-3)

Inputs
------
  --pretrained_ckpt    Path to the MAE-pretrained .pth (ckpt['net'] state_dict)
  --vit_inputs_dir     Per-(sub, ses) 128^3 NIfTIs from 03_prepare_ViT.py:
                       {vit_inputs_dir}/sub-XYZ/ses-{session}/sub-XYZ_ses-{session}_space-ViT128_desc-preproc_T1w.nii.gz
  --matched_labels_csv MRI<->clinical matched labels CSV from 03b. Sole source of
                       per-(subject, visit) labels and image-path resolution. Cohort
                       is constrained to scans matched via Strategy A (VISCODE2
                       direct join) or B (nearest clinical visit within +-14d).
  --session            Single session to fine-tune on. Default 'bl'.
  --long N             Longitudinal mode: subject-level splits expanded to all
                       sessions up to ses-m{12*N}. 'all' = no cap.
  --data_dir           Clinical splits root: {data_dir}/seed_{seed}/{train,val,test}.csv
                       Used for SUBJECT PARTITIONING ONLY (which Patient_IDs go to
                       train/val/test). Labels and image paths come from the
                       matched_labels_csv. Splits are reused from the clinical
                       pipeline so test subjects stay consistent across modalities
                       (per tool.md cross-modality consistency).

Output layout (matches clinical pipeline)
-----------------------------------------
  {out_dir}/ViT_B_mae75/{task}/seed_{seed}/{strategy}/
    metrics.json     {"config": {...}, "test_metrics": {...}}
    best_model.pt    weights with lowest val_loss
    train_log.csv    per-epoch train_loss, val_loss, val_acc, lr
    dataset_manifest.csv  which (Patient_ID, bids_ses, label, image_path) tuples
                          went into each split

Smoke-test usage (local Windows, mri conda env, 28 actual baseline scans)
--------------------------------------------------------------------------
  python mri_pipeline/04_supervised_finetuning_ViT.py \
      --task T1_binary --seed 0 --strategy full_ft \
      --pretrained_ckpt "D:/ViT_B_pretrained_.../p_noaug_..._077000.pth" \
      --epochs 2 --num_workers 0
"""

import argparse
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
    f1_score, precision_score, recall_score, roc_auc_score,
)
from sklearn.utils.class_weight import compute_class_weight
from torch.utils.data import DataLoader

import monai.transforms as mt
from monai.data import Dataset

# Vendored ViT_recipe_for_AD subset
THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))
from _vit_recipe.vit3d import Vision_Transformer3D
from _vit_recipe.checkpoint import load_pretrained_checkpoint

# Cross-pipeline subject-level exclusions (single source of truth in this repo)
sys.path.insert(0, str(THIS_DIR.parent))
from bidsification.exclusions import is_excluded_subject

warnings.filterwarnings("ignore")


# ── Task definitions (mirrors clinical_pipeline/03_encoder_finetune.py + session policy) ──
# session_policy controls how labels are resolved when fine-tuning on multiple
# (Patient_ID, VISCODE_long) rows under --long mode:
#   "current"           -> all sessions valid; label_col = Label_bl_multi for ses-bl rows,
#                          Label_visit_diag for non-bl rows. label_map applied.
#   "baseline_anchored" -> Label_Ny is constant per subject (broadcast on every visit
#                          row by 01b_build_clinical_csv.py via the prog dict). Sessions
#                          are capped to <= label_anchor_max_months to prevent leakage:
#                          at m12 within a 3y window, the 24 months of headroom keep
#                          the prediction non-trivial. With --long unset the session_filter
#                          collapses to ses-bl (matching the legacy 'baseline_only' case).
TASK_CONFIG = {
    "T1_binary": {
        "label_col":      "Label_bl_multi",
        "num_labels":     2,
        "task_type":      "binary",
        "label_map":      {"CN": 0, "MCI": 1, "AD": 1},
        "filter_non_ad":  False,
        "session_policy": "current",
        "description":    "Binary: CN vs MCI+AD",
    },
    "T2_multiclass": {
        "label_col":      "Label_bl_multi",
        "num_labels":     3,
        "task_type":      "multiclass",
        "label_map":      {"CN": 0, "MCI": 1, "AD": 2},
        "filter_non_ad":  False,
        "session_policy": "current",
        "description":    "Multiclass: CN / MCI / AD",
    },
    "T3a_conv3y": {
        "label_col":              "Label_3y",
        "num_labels":             2,
        "task_type":              "binary",
        "label_map":              None,
        "filter_non_ad":          True,
        "session_policy":         "baseline_anchored",
        "label_anchor_max_months": 12,
        "description":            "Prognosis: conversion to AD within 3 years",
    },
    "T3b_conv5y": {
        "label_col":              "Label_5y",
        "num_labels":             2,
        "task_type":              "binary",
        "label_map":              None,
        "filter_non_ad":          True,
        "session_policy":         "baseline_anchored",
        "label_anchor_max_months": 12,
        "description":            "Prognosis: conversion to AD within 5 years",
    },
}


def session_to_months(ses_label: str):
    """Convert VISCODE_long 'bl' -> 0, 'm24' -> 24. None for unsupported labels."""
    if ses_label == "bl":
        return 0
    m = re.match(r"^m(\d+)$", str(ses_label))
    if m:
        return int(m.group(1))
    return None

MODEL_SLUG = "ViT_B_mae75"


# ── CLI ────────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--task",            type=str, required=True,
                   choices=list(TASK_CONFIG.keys()))
    p.add_argument("--seed",            type=int, required=True, choices=[0, 1, 2])
    p.add_argument("--strategy",        type=str, default="full_ft",
                   choices=["full_ft", "frozen"])
    p.add_argument("--pretrained_ckpt", type=str, required=True)
    session_group = p.add_mutually_exclusive_group()
    session_group.add_argument("--session", type=str, default="bl",
                   help="Single session to pair with each Patient_ID (default 'bl'). "
                        "For non-bl sessions the master clinical CSV is consulted for "
                        "per-visit labels.")
    session_group.add_argument("--long", type=str, default=None,
                   help="Longitudinal mode (mutually exclusive with --session): "
                        "'all' = every available session per subject; "
                        "integer N = sessions from baseline up to ses-m{12*N} inclusive "
                        "(e.g. --long 1 for bl + m06 + m12). Subject-level splits are "
                        "preserved (all visits of a subject stay in the same split).")
    p.add_argument("--matched_labels_csv", type=str,
                   default=r"D:\ADNI_BIDS_project\derivatives\mri_clinical_matched\master_mri_clinical_matched_labels.csv",
                   help="Per-(bids_sub, bids_ses) MRI<->clinical matched labels, produced "
                        "by mri_pipeline/03b_match_mri_to_clinical.py. Used both for label "
                        "resolution and to constrain the cohort to scans with valid clinical "
                        "matches (match_status in viscode_exact or nearest_within_14d).")
    p.add_argument("--data_dir",        type=str,
                   default=r"D:\ADNI_BIDS_project\derivatives\clinical\verbose\baseline")
    p.add_argument("--vit_inputs_dir",  type=str,
                   default=r"D:\ADNI_BIDS_project\derivatives\vit_inputs")
    p.add_argument("--out_dir",         type=str,
                   default=str(THIS_DIR / "outputs"))
    p.add_argument("--epochs",          type=int, default=None,
                   help="Default: 50 (full_ft) or 70 (frozen)")
    p.add_argument("--batch_size",      type=int, default=4)
    p.add_argument("--lr",              type=float, default=None,
                   help="Default: 1e-4 (full_ft) or 1e-3 (frozen)")
    p.add_argument("--weight_decay",    type=float, default=1e-4)
    p.add_argument("--warmup_epochs",   type=int, default=5)
    p.add_argument("--patience",        type=int, default=10,
                   help="Early-stopping patience on val_loss")
    p.add_argument("--num_workers",     type=int, default=2)
    p.add_argument("--max_subjects",    type=int, default=None,
                   help="Cap subjects per split (smoke test)")
    args = p.parse_args()

    if args.epochs is None:
        args.epochs = 50 if args.strategy == "full_ft" else 70
    if args.lr is None:
        args.lr = 1e-4 if args.strategy == "full_ft" else 1e-3

    # Resolve session selection mode
    if args.long is None:
        args.long_mode = None
        args.max_months = None
    else:
        if args.long.lower() == "all" or args.long == "0":
            args.long_mode = "all"
            args.max_months = None
        else:
            try:
                n = int(args.long)
                if n <= 0:
                    raise ValueError
            except ValueError:
                p.error(f"--long must be 'all' or a positive integer (got {args.long!r})")
            args.long_mode = "cutoff"
            args.max_months = n * 12
        args.session = None  # cleared while in --long mode (cohort_label() handles display)
    return args


def cohort_label(args, task_cfg) -> str:
    """Human-readable cohort string for metrics.json + tables.

    Examples: 'bl', 'bl + m12', 'bl..m36', 'bl..mAll'.
    Legacy session_policy='baseline_only' always collapses to 'bl' regardless
    of any --long flag (non-bl scans were dropped by the loader).
    """
    sp = task_cfg.get("session_policy", "current")
    if sp == "baseline_only":
        return "bl"
    if args.long_mode is None:
        return args.session or "bl"
    if args.long_mode == "all":
        return "bl..mAll"
    eff = args.max_months
    if sp == "baseline_anchored":
        cap = task_cfg.get("label_anchor_max_months")
        if cap is not None:
            eff = min(eff, int(cap)) if eff is not None else int(cap)
    if eff is None:
        return "bl..mAll"
    if eff == 12:
        return "bl + m12"
    return f"bl..m{int(eff)}"


MATCHED_STATUSES = ("viscode_exact", "nearest_within_14d")


def load_matched_labels(path: str) -> pd.DataFrame:
    """Load the master MRI<->clinical matched labels CSV from 03b.
    Returns DataFrame with Patient_ID added (derived from bids_sub).
    Filtered to scans with a valid clinical match (Strategy A or B)."""
    df = pd.read_csv(path)
    if "match_status" in df.columns:
        df = df[df["match_status"].isin(MATCHED_STATUSES)].copy()

    def _bids_to_pid(s):
        s = str(s)
        # '002S0413' -> '002_S_0413'
        m = re.match(r"^(\d+)S(\d+)$", s)
        return f"{m.group(1)}_S_{m.group(2)}" if m else None

    df["Patient_ID"] = df["bids_sub"].apply(_bids_to_pid)
    df = df.dropna(subset=["Patient_ID"])
    return df


# ── Data loading ───────────────────────────────────────────────────────────────
def patient_id_to_bids_sub(pid: str) -> str:
    """Convert ADNI Patient_ID '002_S_0413' to BIDS subject 'sub-002S0413'."""
    return "sub-" + pid.replace("_", "")


def _resolve_image_path(patient_id: str, viscode: str, vit_inputs_dir: Path) -> str:
    sub = patient_id_to_bids_sub(patient_id)
    ses_label = f"ses-{viscode}"
    return str(vit_inputs_dir / sub / ses_label /
               f"{sub}_{ses_label}_space-ViT128_desc-preproc_T1w.nii.gz")


def _resolve_labels(df: pd.DataFrame, task_cfg: dict) -> pd.DataFrame:
    """Apply session_policy + label_map. Drops rows with NaN labels."""
    if task_cfg["session_policy"] == "current":
        df = df.dropna(subset=["Label_bl_multi", "Label_visit_diag"], how="any").copy()
        def _label(row):
            col = "Label_bl_multi" if row["bids_ses"] == "bl" else "Label_visit_diag"
            v = row[col]
            return task_cfg["label_map"].get(v) if task_cfg["label_map"] else v
        df["label"] = df.apply(_label, axis=1)
    else:
        df = df.dropna(subset=[task_cfg["label_col"]]).copy()
        df["label"] = df[task_cfg["label_col"]]
    df = df.dropna(subset=["label"])
    df["label"] = df["label"].astype(int)
    return df


def load_split_from_matched(baseline_csv: Path, matched_df: pd.DataFrame,
                            task_cfg: dict, vit_inputs_dir: Path,
                            session_filter: str = None, max_months=None
                            ) -> tuple[pd.DataFrame, int]:
    """
    Build a split's labelled scan list by intersecting:
      (a) baseline split CSV (Patient_IDs in this train/val/test partition)
      (b) the master matched_labels CSV (per-scan labels + bids_ses).

    session_filter: 'bl' / 'm12' / ... to keep one session only;
                    None to keep all sessions (subject-level expansion).
    max_months: cap sessions to <= N months from baseline (None = no cap).

    Returns (DataFrame[Patient_ID, bids_ses, label, image_path], n_missing_files).
    """
    subjects = pd.read_csv(baseline_csv, usecols=["Patient_ID"])["Patient_ID"].unique().tolist()
    df = matched_df[matched_df["Patient_ID"].isin(subjects)].copy()

    # Session selection
    if session_filter is not None:
        df = df[df["bids_ses"] == session_filter].copy()
    elif max_months is not None:
        df["_months"] = df["bids_ses"].map(session_to_months)
        df = df.dropna(subset=["_months"])
        df = df[df["_months"] <= max_months].drop(columns=["_months"])

    # Per-task session policy.
    #   baseline_anchored: Label_Ny is constant per subject. With --session bl
    #     (the default) session_filter has already collapsed to ses-bl above, so
    #     this block is a no-op. Under --long, the cap protects against leakage
    #     by clipping any session beyond label_anchor_max_months (e.g. m12 for
    #     T3a/T3b — far enough inside the 3y/5y window).
    if task_cfg["session_policy"] == "baseline_anchored":
        cap = int(task_cfg.get("label_anchor_max_months", 0))
        df["_months"] = df["bids_ses"].map(session_to_months)
        df = df.dropna(subset=["_months"])
        df = df[df["_months"] <= cap].drop(columns=["_months"]).copy()
    elif task_cfg["session_policy"] == "baseline_only":
        # Legacy alias preserved for any older configs still referencing this name.
        df = df[df["bids_ses"] == "bl"].copy()

    # filter_non_ad: drop baseline-AD subjects (T3a/T3b cohort definition)
    if task_cfg["filter_non_ad"]:
        df = df[df["Label_bl_multi"] != "AD"].copy()

    # Subject-level exclusions
    df = df[~df["Patient_ID"].apply(is_excluded_subject)].copy()

    # Resolve labels per session_policy
    df = _resolve_labels(df, task_cfg)

    # Resolve NIfTI paths and drop missing
    df["image_path"] = df.apply(
        lambda r: _resolve_image_path(r["Patient_ID"], r["bids_ses"], vit_inputs_dir),
        axis=1,
    )
    have = df["image_path"].apply(os.path.isfile)
    n_missing = int((~have).sum())
    df = df[have].copy()

    out_cols = ["Patient_ID", "bids_ses", "label", "image_path"]
    return df.reset_index(drop=True)[out_cols], n_missing




def build_dataset(df: pd.DataFrame, train: bool) -> Dataset:
    """
    MONAI Dataset over already-preprocessed 128^3 z-scored RAS NIfTIs.
    Train: light augmentation (flip/rot90/intensity).
    Val/Test: load + ToTensor only.
    """
    items = [{"image": row.image_path, "label": int(row.label)}
             for row in df.itertuples()]

    if train:
        tfm = mt.Compose([
            mt.LoadImaged(keys=["image"]),
            mt.EnsureChannelFirstd(keys=["image"]),
            mt.RandFlipd(keys=["image"], prob=0.2, spatial_axis=0),
            mt.RandFlipd(keys=["image"], prob=0.2, spatial_axis=1),
            mt.RandFlipd(keys=["image"], prob=0.2, spatial_axis=2),
            mt.RandRotate90d(keys=["image"], prob=0.2, max_k=3),
            mt.RandScaleIntensityd(keys=["image"], factors=0.1, prob=0.1),
            mt.RandShiftIntensityd(keys=["image"], offsets=0.1, prob=0.1),
            mt.ToTensord(keys=["image"]),
        ])
    else:
        tfm = mt.Compose([
            mt.LoadImaged(keys=["image"]),
            mt.EnsureChannelFirstd(keys=["image"]),
            mt.ToTensord(keys=["image"]),
        ])
    return Dataset(data=items, transform=tfm)


# ── Metrics (mirror clinical_pipeline/03_encoder_finetune.py) ─────────────────
def compute_test_metrics(y_true: np.ndarray, logits: np.ndarray, task_type: str) -> dict:
    probs = torch.softmax(torch.from_numpy(logits), dim=-1).numpy()
    preds = probs.argmax(axis=-1)

    out = {
        "accuracy":     float(accuracy_score(y_true, preds)),
        "balanced_acc": float(balanced_accuracy_score(y_true, preds)),
    }
    if task_type == "binary":
        positive = probs[:, 1]
        out.update({
            "precision": float(precision_score(y_true, preds, zero_division=0)),
            "recall":    float(recall_score(y_true, preds, zero_division=0)),
            "f1":        float(f1_score(y_true, preds, zero_division=0)),
            "auc_roc":   float(roc_auc_score(y_true, positive))
                         if len(np.unique(y_true)) > 1 else float("nan"),
            "auc_pr":    float(average_precision_score(y_true, positive)),
        })
    else:
        out.update({
            "precision_macro": float(precision_score(y_true, preds, average="macro", zero_division=0)),
            "recall_macro":    float(recall_score(y_true, preds, average="macro", zero_division=0)),
            "macro_f1":        float(f1_score(y_true, preds, average="macro", zero_division=0)),
            "auc_roc_ovr":     float(roc_auc_score(y_true, probs, multi_class="ovr", average="macro"))
                               if len(np.unique(y_true)) > 1 else float("nan"),
            "auc_pr_macro":    float(average_precision_score(
                np.eye(probs.shape[1])[y_true], probs, average="macro")),
        })
    return {k: round(v, 4) if isinstance(v, float) and not np.isnan(v) else v
            for k, v in out.items()}


# ── Train / eval loops ────────────────────────────────────────────────────────
def run_one_epoch(model, loader, criterion, optimizer, scaler, device, train: bool):
    model.train(mode=train)
    total_loss = 0.0
    total_correct = 0
    total_n = 0
    all_logits, all_labels = [], []

    ctx = torch.enable_grad() if train else torch.no_grad()
    with ctx:
        for batch in loader:
            x = batch["image"].to(device, non_blocking=True).float()
            y = batch["label"].to(device, non_blocking=True).long()

            with torch.cuda.amp.autocast(enabled=(device.type == "cuda")):
                logits = model(x)
                loss = criterion(logits, y)

            if train:
                optimizer.zero_grad(set_to_none=True)
                scaler.scale(loss).backward()
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                scaler.step(optimizer)
                scaler.update()

            total_loss += loss.item() * x.size(0)
            preds = logits.argmax(dim=-1)
            total_correct += (preds == y).sum().item()
            total_n += x.size(0)
            all_logits.append(logits.detach().float().cpu().numpy())
            all_labels.append(y.detach().cpu().numpy())

    return (total_loss / max(total_n, 1),
            total_correct / max(total_n, 1),
            np.concatenate(all_logits, axis=0) if all_logits else np.empty((0,)),
            np.concatenate(all_labels, axis=0) if all_labels else np.empty((0,)))


def lr_at(epoch_idx: int, args) -> float:
    """Linear warmup for warmup_epochs, then cosine decay to 5e-6."""
    eta_min = 5e-6
    if epoch_idx < args.warmup_epochs:
        return args.lr * (epoch_idx + 1) / max(args.warmup_epochs, 1)
    progress = (epoch_idx - args.warmup_epochs) / max(args.epochs - args.warmup_epochs, 1)
    return eta_min + 0.5 * (args.lr - eta_min) * (1 + np.cos(np.pi * progress))


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    task_cfg = TASK_CONFIG[args.task]

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    out_dir = Path(args.out_dir) / MODEL_SLUG / args.task / f"seed_{args.seed}" / args.strategy
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print(f"  04_supervised_finetuning_ViT — {args.task} | seed={args.seed} | {args.strategy}")
    print(f"  Device: {device}")
    print(f"  Output: {out_dir}")
    print("=" * 70)

    # ── Load splits ──────────────────────────────────────────────────────────
    seed_dir = Path(args.data_dir) / f"seed_{args.seed}"
    vit_dir = Path(args.vit_inputs_dir)

    if not os.path.isfile(args.matched_labels_csv):
        raise FileNotFoundError(
            f"--matched_labels_csv not found: {args.matched_labels_csv}\n"
            f"Run mri_pipeline/03b_match_mri_to_clinical.py first to produce it.")
    matched_df = load_matched_labels(args.matched_labels_csv)
    print(f"  Loaded matched labels CSV: {len(matched_df)} matched (subject, visit) rows")

    if args.long_mode is None:
        # Single-session mode (default: ses-bl).
        train_df, n_miss_tr = load_split_from_matched(
            seed_dir / "train.csv", matched_df, task_cfg, vit_dir, session_filter=args.session)
        val_df,   n_miss_va = load_split_from_matched(
            seed_dir / "val.csv",   matched_df, task_cfg, vit_dir, session_filter=args.session)
        test_df,  n_miss_te = load_split_from_matched(
            seed_dir / "test.csv",  matched_df, task_cfg, vit_dir, session_filter=args.session)
        cohort_note = f"single-session ses-{args.session}"
    else:
        # Longitudinal mode: subject-level splits expanded to multiple visits per subject.
        if args.long_mode == "cutoff":
            cohort_note = f"longitudinal up to ses-m{args.max_months} ({args.long}y)"
        else:
            cohort_note = "longitudinal: all sessions"
        if task_cfg["session_policy"] == "baseline_anchored":
            cap = task_cfg.get("label_anchor_max_months", 0)
            print(f"  [info] task {args.task} has session_policy='baseline_anchored' "
                  f"(Label_{task_cfg['label_col'].split('_')[-1]} is constant per subject; "
                  f"sessions capped to <= m{cap}).")
        elif task_cfg["session_policy"] == "baseline_only":
            print(f"  [info] task {args.task} has session_policy='baseline_only' "
                  "(non-bl scans dropped — legacy behavior).")
        train_df, n_miss_tr = load_split_from_matched(
            seed_dir / "train.csv", matched_df, task_cfg, vit_dir, max_months=args.max_months)
        val_df,   n_miss_va = load_split_from_matched(
            seed_dir / "val.csv",   matched_df, task_cfg, vit_dir, max_months=args.max_months)
        test_df,  n_miss_te = load_split_from_matched(
            seed_dir / "test.csv",  matched_df, task_cfg, vit_dir, max_months=args.max_months)

    if args.max_subjects:
        train_df = train_df.head(args.max_subjects)
        val_df   = val_df.head(min(args.max_subjects, len(val_df)))
        test_df  = test_df.head(min(args.max_subjects, len(test_df)))

    print(f"  Cohort: {cohort_note}")
    print(f"  Splits — train: {len(train_df)} rows (NIfTI missing dropped: {n_miss_tr})  "
          f"val: {len(val_df)} ({n_miss_va})  test: {len(test_df)} ({n_miss_te})")
    if args.long_mode is not None:
        print(f"           train: {train_df['Patient_ID'].nunique()} unique subjects, "
              f"val: {val_df['Patient_ID'].nunique()}, test: {test_df['Patient_ID'].nunique()}")
    if min(len(train_df), len(val_df), len(test_df)) == 0:
        raise RuntimeError("Empty split after filtering — check vit_inputs_dir / master_clinical_csv.")

    # Save dataset manifest for the run (which (sub, ses, label, image_path) tuples were used)
    manifest_rows = []
    for split_name, sdf in [("train", train_df), ("val", val_df), ("test", test_df)]:
        m = sdf.copy()
        m["split"] = split_name
        manifest_rows.append(m)
    manifest_df = pd.concat(manifest_rows, ignore_index=True)
    manifest_cols = ["split", "Patient_ID", "bids_ses", "label", "image_path"]
    manifest_cols = [c for c in manifest_cols if c in manifest_df.columns]
    manifest_df[manifest_cols].to_csv(out_dir / "dataset_manifest.csv", index=False)

    # ── Class weights from train ─────────────────────────────────────────────
    classes = np.unique(train_df["label"].values)
    cw = compute_class_weight("balanced", classes=classes, y=train_df["label"].values)
    class_weights = torch.tensor(cw, dtype=torch.float, device=device)
    print(f"  Class weights: {dict(zip(classes.tolist(), cw.round(3).tolist()))}")

    # ── Datasets / dataloaders ───────────────────────────────────────────────
    train_loader = DataLoader(build_dataset(train_df, train=True),
                              batch_size=args.batch_size, shuffle=True,
                              num_workers=args.num_workers, pin_memory=True,
                              drop_last=True)
    val_loader = DataLoader(build_dataset(val_df, train=False),
                            batch_size=1, shuffle=False,
                            num_workers=args.num_workers, pin_memory=True)
    test_loader = DataLoader(build_dataset(test_df, train=False),
                             batch_size=1, shuffle=False,
                             num_workers=args.num_workers, pin_memory=True)

    # ── Model + checkpoint ───────────────────────────────────────────────────
    model = Vision_Transformer3D(
        img_size=(128, 128, 128), patch_size=16, in_chans=1,
        n_classes=task_cfg["num_labels"],
        embed_dim=768, depth=12, n_heads=12, mlp_ratio=4.0,
        qkv_bias=True, drop_path_rate=0.1, p=0.0, attn_p=0.1,
        global_avg_pool=False, pos_embed_type="learnable",
        patch_embed_fun="conv3d",
    )
    model = load_pretrained_checkpoint(model, args.pretrained_ckpt)
    model = model.to(device)

    if args.strategy == "frozen":
        for n, p in model.named_parameters():
            if not n.startswith("head"):
                p.requires_grad = False
        n_train = sum(p.numel() for p in model.parameters() if p.requires_grad)
        print(f"  Frozen strategy: {n_train:,} trainable params (head only)")

    # ── Optimizer / scheduler / scaler / criterion ───────────────────────────
    trainable = [p for p in model.parameters() if p.requires_grad]
    optimizer = torch.optim.AdamW(trainable, lr=args.lr,
                                  weight_decay=args.weight_decay,
                                  betas=(0.9, 0.999))
    scaler = torch.cuda.amp.GradScaler(enabled=(device.type == "cuda"))
    criterion = nn.CrossEntropyLoss(weight=class_weights)

    # ── Training loop ────────────────────────────────────────────────────────
    best_val = float("inf")
    best_epoch = -1
    epochs_since_improve = 0
    log_rows = []

    for epoch in range(args.epochs):
        cur_lr = lr_at(epoch, args)
        for g in optimizer.param_groups:
            g["lr"] = cur_lr

        tr_loss, tr_acc, _, _ = run_one_epoch(
            model, train_loader, criterion, optimizer, scaler, device, train=True)
        va_loss, va_acc, _, _ = run_one_epoch(
            model, val_loader, criterion, optimizer, scaler, device, train=False)

        log_rows.append({"epoch": epoch + 1, "lr": cur_lr,
                         "train_loss": tr_loss, "train_acc": tr_acc,
                         "val_loss": va_loss, "val_acc": va_acc})
        print(f"  [epoch {epoch+1:>3}/{args.epochs}] "
              f"lr={cur_lr:.2e}  train_loss={tr_loss:.4f}  train_acc={tr_acc:.4f}  "
              f"val_loss={va_loss:.4f}  val_acc={va_acc:.4f}")

        if va_loss < best_val:
            best_val = va_loss
            best_epoch = epoch + 1
            epochs_since_improve = 0
            torch.save({"net": model.state_dict(), "epoch": best_epoch},
                       out_dir / "best_model.pt")
        else:
            epochs_since_improve += 1
            if epochs_since_improve >= args.patience:
                print(f"  Early stopping at epoch {epoch+1} "
                      f"(no val improvement for {args.patience}).")
                break

    pd.DataFrame(log_rows).to_csv(out_dir / "train_log.csv", index=False)

    # ── Final test on best weights ────────────────────────────────────────────
    best_state = torch.load(out_dir / "best_model.pt", map_location=device)
    model.load_state_dict(best_state["net"])
    _, _, test_logits, test_labels = run_one_epoch(
        model, test_loader, criterion, optimizer, scaler, device, train=False)
    test_metrics = compute_test_metrics(test_labels, test_logits, task_cfg["task_type"])
    print(f"  Test metrics: {test_metrics}")

    # ── Save metrics.json (matches clinical pipeline schema) ─────────────────
    config = {
        "model_id":         MODEL_SLUG,
        "pretrained_ckpt":  args.pretrained_ckpt,
        "task":             args.task,
        "task_description": task_cfg["description"],
        "session":          cohort_label(args, task_cfg),
        "long_mode":        args.long_mode,
        "max_months":       args.max_months,
        "session_policy":   task_cfg["session_policy"],
        "seed":             args.seed,
        "strategy":         args.strategy,
        "epochs":           args.epochs,
        "best_epoch":       best_epoch,
        "lr":               args.lr,
        "weight_decay":     args.weight_decay,
        "batch_size":       args.batch_size,
        "warmup_epochs":    args.warmup_epochs,
        "patience":         args.patience,
        "n_train":          int(len(train_df)),
        "n_val":            int(len(val_df)),
        "n_test":           int(len(test_df)),
        "class_weights":    {int(c): round(float(w), 4) for c, w in zip(classes, cw)},
        "timestamp":        datetime.now().isoformat(),
    }
    with open(out_dir / "metrics.json", "w") as f:
        json.dump({"config": config, "test_metrics": test_metrics}, f, indent=2)

    print(f"  Saved: {out_dir}/metrics.json")
    print("=" * 70)


if __name__ == "__main__":
    main()
