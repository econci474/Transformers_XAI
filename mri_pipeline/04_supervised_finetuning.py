"""
04_supervised_finetuning.py
================================
Supervised fine-tuning of the MAE-pretrained ViT-B/3D
(qasymjomart/ViT_recipe_for_AD) on ADNI baseline T1w MRIs preprocessed by
03_prepare_ViT.py (128 x 128 x 128 @ 1.75 mm RAS, z-scored on nonzero voxels).

Tasks (mirroring clinical_pipeline/03_encoder_finetune.py)
----------------------------------------------------------
  T1_binary     : CN vs MCI+AD          (binary)
  T1b_binary    : CN+MCI vs AD          (binary)
  T1c_binary    : CN vs AD              (binary, MCI excluded)
  T2_multiclass : CN / MCI / AD         (3-class)
  T3a_conv3y    : Conversion to AD <=3y (binary, non-AD at baseline)
  T3b_conv5y    : Conversion to AD <=5y (binary, non-AD at baseline)

Strategies
----------
  full_ft : MAE-pretrained encoder + head, all parameters trained (lr=1e-4)
  frozen  : MAE-pretrained encoder frozen, classification head only (lr=1e-3)
  scratch : random init, NO pretrained checkpoint, all parameters trained
            (lr=1e-3) — the ablation isolating the value of MAE pretraining

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
                       (cross-modality consistency; see mri_pipeline/README.md).

Output layout (matches clinical pipeline)
-----------------------------------------
  {out_dir}/ViT_B_mae75/{task}/seed_{seed}/{strategy}/
    metrics.json      {"config": {...}, "test_metrics": {...},
                       "test_diagnostics": {confusion_matrix, per_class
                       precision/recall/sensitivity/specificity/f1/support,
                       labels}}
    best_model.pt     weights with highest val balanced accuracy
    last_checkpoint.pt  atomic per-epoch resume state (model+optim+scaler+
                        counters); auto-resumed on rerun unless --no_resume.
                        Safety net for the 12 h SLURM/Colab wall-time cap.
    train_log.csv     per-epoch train_loss, val_loss, val_acc, lr (flushed
                      every epoch, not just at the end)
    test_predictions.csv  per-sample y_true, y_pred, prob_0..prob_{K-1}
    dataset_manifest.csv  which (Patient_ID, bids_ses, label, image_path) tuples
                          went into each split

  Augmentation (--augment): none | random (default; current behaviour) |
  plus_original (originals + --aug_copies always-on augmented copies).
  Optional wandb logging via --wandb (no-op otherwise; honours WANDB_MODE).

Smoke-test usage (local Windows, mri conda env, 28 actual baseline scans)
--------------------------------------------------------------------------
  python mri_pipeline/04_supervised_finetuning.py \
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

import torch  # import torch before numpy/pandas (Windows MKL DLL load order)
import torch.nn as nn
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score, average_precision_score, balanced_accuracy_score,
    confusion_matrix, f1_score, precision_recall_fscore_support,
    precision_score, recall_score, roc_auc_score,
)
from sklearn.utils.class_weight import compute_class_weight
from torch.utils.data import ConcatDataset, DataLoader

import monai.transforms as mt
import nibabel as nib
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
    "T1b_binary": {
        "label_col":      "Label_bl_multi",
        "num_labels":     2,
        "task_type":      "binary",
        "label_map":      {"CN": 0, "MCI": 0, "AD": 1},
        "filter_non_ad":  False,
        "session_policy": "current",
        "description":    "Binary: CN+MCI vs AD",
    },
    "T1c_binary": {
        # CN vs AD only — the cleanest binary contrast (no ambiguous MCI middle
        # group). MCI is deliberately absent from label_map, so _resolve_labels
        # maps MCI scans to None and drops them (see _resolve_labels, which ends
        # with dropna(subset=["label"])). Cohort ~216 AD + ~507 CN scans.
        "label_col":      "Label_bl_multi",
        "num_labels":     2,
        "task_type":      "binary",
        "label_map":      {"CN": 0, "AD": 1},
        "filter_non_ad":  False,
        "session_policy": "current",
        "description":    "Binary: CN vs AD (MCI excluded)",
    },
    "T1d_binary": {
        # Conversion task: among MCI subjects, predict whether they
        # progressed to AD (pMCI=1) vs stayed stable (sMCI=1). Labels come
        # from explicit pMCI/sMCI columns in the matched CSV (present on
        # the HPC `_extended_post_exclusion` variant; not in the local
        # master CSV, so T1d only runs on HPC). The earliest scan per
        # subject is kept -- conversion is a subject-level property, so
        # multiple sessions would over-weight the well-imaged subjects.
        "label_col":      None,
        "num_labels":     2,
        "task_type":      "binary",
        "label_map":      None,
        "filter_non_ad":  False,
        "session_policy": "conversion",
        "pos_col":        "pMCI",
        "neg_col":        "sMCI",
        "description":    "Binary: pMCI vs sMCI (MCI->AD conversion, earliest scan/subject)",
    },
    "T1e_pcn_vs_scn": {
        # Conversion task: among CN-baseline subjects, predict whether they
        # progressed (pCN -- to MCI OR to AD) vs stayed stable (sCN). Labels
        # come from pCN_to_AD + pCN_to_MCI + sCN columns in the matched CSV
        # (present on the HPC `_extended_post_exclusion` variant only). The
        # earliest scan per subject is kept (same rationale as T1d).
        # Cohort: 60 pCN + 152 sCN = 212 subjects, ~555 scans.
        "label_col":      None,
        "num_labels":     2,
        "task_type":      "binary",
        "label_map":      None,
        "filter_non_ad":  False,
        "session_policy": "conversion",
        "pos_cols":       ["pCN_to_AD", "pCN_to_MCI"],
        "neg_col":        "sCN",
        "description":    "Binary: pCN (CN->MCI or CN->AD) vs sCN, earliest scan/subject",
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
    "T3c_conv7y": {
        # Extends T3a/T3b to a 7-year prognostic horizon. Label_7y is
        # already in the master (01b_build_clinical_csv.py emits Label_1y
        # through Label_10y); test-cohort positive rate ~29%.
        # Cohort caveat: subjects need >= 7 y of follow-up to be labelled;
        # ADNI3-late subjects with shorter FU will fall out (Label_7y NaN
        # -> dropped by _resolve_labels).
        "label_col":              "Label_7y",
        "num_labels":             2,
        "task_type":              "binary",
        "label_map":              None,
        "filter_non_ad":          True,
        "session_policy":         "baseline_anchored",
        "label_anchor_max_months": 12,
        "description":            "Prognosis: conversion to AD within 7 years",
    },
    "T3d_conv10y": {
        # 10-year prognostic horizon. Label_10y is already in the master.
        # Cohort caveat: subjects need >= 10 y of follow-up to be labelled;
        # T3d is effectively an ADNI1/GO/2 analysis since most ADNI3
        # subjects only have 3-5 y of FU. Test-cohort positive rate ~33%.
        "label_col":              "Label_10y",
        "num_labels":             2,
        "task_type":              "binary",
        "label_map":              None,
        "filter_non_ad":          True,
        "session_policy":         "baseline_anchored",
        "label_anchor_max_months": 12,
        "description":            "Prognosis: conversion to AD within 10 years",
    },
    "T4_conv_horizon": {
        # Exploratory 3-class ordinal: predict the AD-conversion horizon
        # bucket among confirmed converters (pMCI + pCN_to_AD). Buckets:
        #   class 0: years_to_AD < 3
        #   class 1: 3 <= years_to_AD < 7
        #   class 2: years_to_AD >= 7
        # Cohort: 146 subjects (121 pMCI + 25 pCN_to_AD). Label_T4 column
        # is added to the MRI master by 01e_build_T4_labels_and_splits.py.
        # The T4 splits live under
        # derivatives/clinical/no_cdr_stratified_post_exclusion/tabular/
        # baseline_T4/seed_{0,1,2}/{train,val,test}.csv (point --data_dir
        # at this when submitting). Subjects outside the T4 cohort have
        # Label_T4=NaN and are dropped automatically by _resolve_labels.
        "label_col":              "Label_T4",
        "num_labels":             3,
        "task_type":              "multiclass",
        "label_map":              None,
        "filter_non_ad":          False,
        "session_policy":         "baseline_anchored",
        # bl scans are nearly absent in the matched T4 cohort (~11 total), which
        # left seed-0's test fold empty. Use 12 like T3a-d so T4 MRI is anchored
        # on the abundant m12 scans (subject-level Label_T4 is broadcast to every
        # visit row), matching the CL_bl + MRI_m12 late-fusion framing.
        "label_anchor_max_months": 12,
        "description":            "Multiclass: AD-conversion horizon "
                                  "(3-bucket: <3y / 3-7y / >=7y) on pMCI + pCN_to_AD",
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

# ── Architecture sizings (Vaswani / DeiT canonical) ───────────────────────────
# embed_dim, depth, n_heads, mlp_ratio
VIT_SIZES = {
    "tiny":  dict(embed_dim=192, depth=12, n_heads=3,  mlp_ratio=4.0),
    "small": dict(embed_dim=384, depth=12, n_heads=6,  mlp_ratio=4.0),
    "base":  dict(embed_dim=768, depth=12, n_heads=12, mlp_ratio=4.0),
}
VIT_SIZE_SLUG = {"tiny": "ViT_T", "small": "ViT_S", "base": "ViT_B"}


# ── CLI ────────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--task",            type=str, required=True,
                   choices=list(TASK_CONFIG.keys()))
    p.add_argument("--seed",            type=int, required=True, choices=[0, 1, 2])
    p.add_argument("--strategy",        type=str, default="full_ft",
                   choices=["full_ft", "frozen", "scratch"])
    p.add_argument("--vit_size",        type=str, default="base",
                   choices=["tiny", "small", "base"],
                   help="ViT architecture size. 'base' = paper-default ViT-B "
                        "(embed_dim=768, depth=12, heads=12, ~86M params). "
                        "'small' = ViT-S (384, 12, 6, ~22M). 'tiny' = ViT-Ti "
                        "(192, 12, 3, ~5.5M). The MAE-pretrained checkpoint "
                        "is ViT-B only; --vit_size != base requires "
                        "--strategy scratch.")
    p.add_argument("--pretrained_ckpt", type=str, default=None,
                   help="MAE-pretrained .pth (ckpt['net'] state_dict). Required "
                        "for full_ft / frozen; ignored for scratch.")
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
                   help="Default: 50 full_ft / 70 frozen / 100 scratch")
    p.add_argument("--batch_size",      type=int, default=4)
    p.add_argument("--lr",              type=float, default=None,
                   help="Default: 1e-4 full_ft / 1e-3 frozen / 1e-3 scratch")
    p.add_argument("--weight_decay",    type=float, default=1e-4)
    p.add_argument("--warmup_epochs",   type=int, default=None,
                   help="Default: 10 full_ft / 5 frozen / 10 scratch")
    p.add_argument("--patience",        type=int, default=10,
                   help="Early-stopping patience on val balanced accuracy")
    p.add_argument("--grad_accum_steps", type=int, default=8,
                   help="Micro-batches accumulated before an optimiser step. "
                        "Effective batch = batch_size * grad_accum_steps "
                        "(default 8 -> 32 at batch_size=4).")
    p.add_argument("--llrd_gamma",      type=float, default=None,
                   help="Layer-wise LR decay for full_ft: each ViT block deeper "
                        "from the head gets lr *= gamma. Default 0.70 for "
                        "full_ft/frozen, 1.0 for scratch (no decay — a "
                        "random-init net trains all layers at uniform LR).")
    p.add_argument("--num_workers",     type=int, default=2)
    p.add_argument("--max_subjects",    type=int, default=None,
                   help="Cap subjects per split (smoke test)")
    # ── Augmentation control ─────────────────────────────────────────────────
    p.add_argument("--augment",         type=str, default="random",
                   choices=["none", "random", "plus_original"],
                   help="none = no aug (eval transform on train); "
                        "random = on-the-fly stochastic aug, originals not retained "
                        "(default — current behaviour); "
                        "plus_original = originals (no aug) + K always-on augmented "
                        "copies (set size = (1+K)*N_train).")
    p.add_argument("--aug_copies",      type=int, default=1,
                   help="plus_original only: K always-on augmented copies per "
                        "original (default 1 -> 2x train set).")
    # ── Sweepable regularisation hyperparameters ─────────────────────────────
    p.add_argument("--drop_path_rate",  type=float, default=0.1,
                   help="ViT stochastic-depth rate (model drop_path_rate).")
    p.add_argument("--attn_dropout",    type=float, default=0.1,
                   help="ViT attention dropout (model attn_p).")
    p.add_argument("--label_smoothing", type=float, default=0.0,
                   help="CrossEntropyLoss label smoothing.")
    # ── Resume ───────────────────────────────────────────────────────────────
    p.add_argument("--no_resume",       action="store_true",
                   help="Ignore any last_checkpoint.pt and start fresh.")
    p.add_argument("--val_test",        action="store_true",
                   help="No training: load best_model.pt and recompute val + test, "
                        "patching a `val_metrics` block into metrics.json. In this "
                        "mode --out_dir IS the run dir. Pass the same "
                        "--strategy/--vit_size/--augment the run was trained with.")
    # ── wandb (optional; no-op unless --wandb) ───────────────────────────────
    p.add_argument("--wandb",           action="store_true",
                   help="Enable Weights & Biases logging (no-op if absent).")
    p.add_argument("--wandb_project",   type=str, default="vit-mri-finetune")
    p.add_argument("--wandb_entity",    type=str, default=None)
    p.add_argument("--wandb_run_name",  type=str, default=None)
    args = p.parse_args()

    # Strategy-dependent defaults. scratch trains a random-init ViT-B, which
    # needs longer and a higher LR than fine-tuning a pretrained net (tunable).
    if args.epochs is None:
        args.epochs = {"full_ft": 50, "frozen": 70, "scratch": 100}[args.strategy]
    if args.lr is None:
        args.lr = {"full_ft": 1e-4, "frozen": 1e-3, "scratch": 1e-3}[args.strategy]
    if args.warmup_epochs is None:
        args.warmup_epochs = {"full_ft": 10, "frozen": 5, "scratch": 10}[args.strategy]
    if args.llrd_gamma is None:
        # LLRD is a fine-tuning technique; a from-scratch net uses uniform LR.
        args.llrd_gamma = 1.0 if args.strategy == "scratch" else 0.70
    # --pretrained_ckpt is required for full_ft / frozen TRAINING, ignored for
    # scratch and for --val_test (which loads best_model.pt, not the pretrained).
    if (args.strategy in ("full_ft", "frozen") and not args.pretrained_ckpt
            and not args.val_test):
        p.error(f"--pretrained_ckpt is required for strategy '{args.strategy}'")
    # The MAE-pretrained .pth is ViT-B only; non-base sizes are scratch-only.
    if args.vit_size != "base" and args.strategy != "scratch":
        p.error(f"--vit_size={args.vit_size} only supports --strategy scratch "
                f"(MAE checkpoint is ViT-B only).")

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


def class_names_for(task_cfg: dict) -> list:
    """Per-index human-readable class names for confusion-matrix labelling.

    label_map may be many-to-one (e.g. T1_binary CN->0, MCI->1, AD->1), so each
    index collects every source name mapped to it ('MCI+AD'). label_map=None
    falls back to pos_col/neg_col (conversion tasks) or generic 'class_{i}'.
    """
    n = task_cfg["num_labels"]
    # Conversion tasks: pos_col(s)=class 1, neg_col(s)=class 0.
    if task_cfg.get("session_policy") == "conversion":
        pos = task_cfg.get("pos_cols") or [task_cfg.get("pos_col", "class_1")]
        neg = task_cfg.get("neg_cols") or [task_cfg.get("neg_col", "class_0")]
        return ["|".join(c for c in neg if c is not None) or "class_0",
                "|".join(c for c in pos if c is not None) or "class_1"]
    lm = task_cfg.get("label_map")
    if not lm:
        return [f"class_{i}" for i in range(n)]
    names = []
    for i in range(n):
        members = [k for k, v in lm.items() if v == i]
        names.append("+".join(members) if members else f"class_{i}")
    return names


def init_wandb(args, task_cfg: dict, extra=None):
    """Start a wandb run if --wandb is set; otherwise return None (full no-op).

    Honours WANDB_API_KEY / WANDB_MODE (offline|disabled) / WANDB_DIR from the
    environment automatically. Degrades to local-only if wandb is not installed.
    """
    if not args.wandb:
        return None
    try:
        import wandb
    except ImportError:
        print("  [WARN] --wandb set but wandb not installed; continuing local-only.")
        return None
    name = args.wandb_run_name or (
        f"{args.task}-s{args.seed}-{args.strategy}-{args.augment}")
    cfg = {
        "task":            args.task,
        "seed":            args.seed,
        "strategy":        args.strategy,
        "vit_size":        args.vit_size,
        "augment":         args.augment,
        "aug_copies":      args.aug_copies,
        "epochs":          args.epochs,
        "lr":              args.lr,
        "weight_decay":    args.weight_decay,
        "batch_size":      args.batch_size,
        "grad_accum_steps": args.grad_accum_steps,
        "llrd_gamma":      args.llrd_gamma,
        "warmup_epochs":   args.warmup_epochs,
        "patience":        args.patience,
        "drop_path_rate":  args.drop_path_rate,
        "attn_dropout":    args.attn_dropout,
        "label_smoothing": args.label_smoothing,
        "long_mode":       args.long_mode,
        "max_months":      args.max_months,
        "session_policy":  task_cfg["session_policy"],
        "num_labels":      task_cfg["num_labels"],
        "pretrained_ckpt": args.pretrained_ckpt,
    }
    if extra:
        cfg.update(extra)
    return wandb.init(project=args.wandb_project, entity=args.wandb_entity,
                      name=name, config=cfg, reinit=True)


MATCHED_STATUSES = ("viscode_exact", "viscode2_exact", "nearest_within_14d")


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
    """Apply session_policy + label_map. Drops rows with NaN labels.

    Session policies:
      current            label depends on bids_ses (Label_bl_multi for bl;
                         Label_visit_diag otherwise). T1, T1b, T1c, T2.
      baseline_anchored  label_col is constant per subject (Label_3y etc.).
                         T3a, T3b.
      conversion         label = 1 if pos_col==1, 0 if neg_col==1; rows
                         missing both are dropped; then collapse to one
                         (earliest) session per subject. T1d.
    """
    if task_cfg["session_policy"] == "current":
        df = df.dropna(subset=["Label_bl_multi", "Label_visit_diag"], how="any").copy()
        def _label(row):
            col = "Label_bl_multi" if row["bids_ses"] == "bl" else "Label_visit_diag"
            v = row[col]
            return task_cfg["label_map"].get(v) if task_cfg["label_map"] else v
        df["label"] = df.apply(_label, axis=1)
    elif task_cfg["session_policy"] == "conversion":
        # Tasks support EITHER single pos_col/neg_col OR pos_cols/neg_cols
        # (plural, list) for OR semantics across multiple positive/negative
        # indicator columns. T1e uses pos_cols=["pCN_to_AD", "pCN_to_MCI"]
        # because pCN is the union of two conversion subtypes.
        pos_cols = list(task_cfg.get("pos_cols") or [task_cfg.get("pos_col")])
        neg_cols = list(task_cfg.get("neg_cols") or [task_cfg.get("neg_col")])
        pos_cols = [c for c in pos_cols if c is not None]
        neg_cols = [c for c in neg_cols if c is not None]
        if not pos_cols or not neg_cols:
            raise KeyError(f"Conversion task {task_cfg.get('description', '?')!r} "
                           f"needs at least one pos_col/pos_cols and neg_col/neg_cols.")
        required = pos_cols + neg_cols
        for c in required:
            if c not in df.columns:
                raise KeyError(
                    f"matched_labels CSV is missing the {c!r} column "
                    f"required by conversion task '{task_cfg.get('description', '?')}'. "
                    f"This task only works on the HPC '_extended_post_exclusion' CSV.")
        # _is_one tolerates "1", "1.0", True, etc.
        def _is_one(v):
            return str(v).strip() in ("1", "1.0", "True", "true")
        df = df.dropna(subset=required, how="any").copy()
        df["label"] = df.apply(
            lambda r: (1 if any(_is_one(r[c]) for c in pos_cols)
                       else (0 if any(_is_one(r[c]) for c in neg_cols) else None)),
            axis=1)
        df = df.dropna(subset=["label"]).copy()
        df["label"] = df["label"].astype(int)
        # Earliest session per subject (conversion is subject-level).
        df["_months"] = df["bids_ses"].map(session_to_months)
        df = df.dropna(subset=["_months"])
        df = (df.sort_values(["Patient_ID", "_months"])
                .groupby("Patient_ID", as_index=False).first()
                .drop(columns="_months"))
        return df  # already int-cast, return early to skip the int recast below
    else:
        label_col = task_cfg["label_col"]
        if label_col not in df.columns:
            # An empty split collapses to a 0-column frame, which makes dropna
            # raise a misleading KeyError on the label column. Surface the real
            # cause (empty split vs genuinely-missing column) instead.
            if len(df) == 0:
                raise ValueError(
                    f"Split is empty (0 rows) before resolving label_col="
                    f"{label_col!r}. Likely no scans survived the session filter "
                    f"(e.g. baseline_anchored with too small label_anchor_max_months).")
            raise KeyError(
                f"label_col {label_col!r} not in matched-labels columns: "
                f"{list(df.columns)}")
        df = df.dropna(subset=[label_col]).copy()
        df["label"] = df[label_col]
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




def _eval_transform() -> mt.Compose:
    """Load + channel-first + tensor only (no augmentation)."""
    return mt.Compose([
        mt.LoadImaged(keys=["image"]),
        mt.EnsureChannelFirstd(keys=["image"]),
        mt.ToTensord(keys=["image"]),
    ])


def _train_transform(p_flip=0.2, p_rot=0.2, p_scale=0.1, p_shift=0.1) -> mt.Compose:
    """Stochastic train aug. Default probs == original behaviour; pass 1.0 for
    an always-on (deterministic-fire) augmented branch in plus_original mode."""
    return mt.Compose([
        mt.LoadImaged(keys=["image"]),
        mt.EnsureChannelFirstd(keys=["image"]),
        mt.RandFlipd(keys=["image"], prob=p_flip, spatial_axis=0),
        mt.RandFlipd(keys=["image"], prob=p_flip, spatial_axis=1),
        mt.RandFlipd(keys=["image"], prob=p_flip, spatial_axis=2),
        mt.RandRotate90d(keys=["image"], prob=p_rot, max_k=3),
        mt.RandScaleIntensityd(keys=["image"], factors=0.1, prob=p_scale),
        mt.RandShiftIntensityd(keys=["image"], offsets=0.1, prob=p_shift),
        mt.ToTensord(keys=["image"]),
    ])


def build_dataset(df: pd.DataFrame, train: bool,
                  augment: str = "random", aug_copies: int = 1):
    """
    MONAI Dataset over already-preprocessed 128^3 z-scored RAS NIfTIs.

    augment (train only; val/test always use the no-aug eval transform):
      "none"          -> eval transform on train (no augmentation at all)
      "random"        -> on-the-fly stochastic aug, originals NOT retained
                         (exact original behaviour; the default)
      "plus_original" -> ConcatDataset of every original (no aug) PLUS
                         `aug_copies` always-on (prob=1.0) augmented copies, so
                         originals are retained and the augmented copies are
                         genuinely distinct. Train set size = (1+aug_copies)*N.
    """
    items = [{"image": row.image_path, "label": int(row.label)}
             for row in df.itertuples()]

    if not train or augment == "none":
        return Dataset(data=items, transform=_eval_transform())
    if augment == "random":
        return Dataset(data=items, transform=_train_transform())
    if augment == "plus_original":
        k = max(1, int(aug_copies))
        orig = Dataset(data=items, transform=_eval_transform())
        aug = [Dataset(data=items, transform=_train_transform(1.0, 1.0, 1.0, 1.0))
               for _ in range(k)]
        return ConcatDataset([orig] + aug)
    raise ValueError(f"Unknown augment mode: {augment!r}")


# ── Metrics (mirror clinical_pipeline/03_encoder_finetune.py) ─────────────────
def compute_test_metrics(y_true: np.ndarray, logits: np.ndarray, task_type: str):
    """Returns (scalar_metrics_dict, preds, probs)."""
    probs = torch.softmax(torch.from_numpy(logits), dim=-1).numpy()
    preds = probs.argmax(axis=-1)

    out = {
        "accuracy":     float(accuracy_score(y_true, preds)),
        "balanced_acc": float(balanced_accuracy_score(y_true, preds)),
    }
    if task_type == "binary":
        positive = probs[:, 1]
        out.update({
            "precision":   float(precision_score(y_true, preds, zero_division=0)),
            "recall":      float(recall_score(y_true, preds, zero_division=0)),
            # sensitivity == recall/TPR for the positive class (class 1);
            # specificity == TNR == recall for the negative class (class 0).
            "sensitivity": float(recall_score(y_true, preds, pos_label=1,
                                               zero_division=0)),
            "specificity": float(recall_score(y_true, preds, pos_label=0,
                                               zero_division=0)),
            "f1":          float(f1_score(y_true, preds, zero_division=0)),
            "auc_roc":     float(roc_auc_score(y_true, positive))
                           if len(np.unique(y_true)) > 1 else float("nan"),
            "auc_pr":      float(average_precision_score(y_true, positive)),
        })
    else:
        # Multiclass OVR-AUC throws if a (small) fold is missing a class — the
        # len(unique)>1 guard is insufficient (it passes with 2 of 3 classes but
        # probs has all 3 columns). Guard both AUCs so a degenerate val/test fold
        # (common for tiny conversion cohorts, e.g. T4 seed folds) yields nan
        # rather than crashing the whole run after training.
        n_cls = probs.shape[1]
        try:
            auc_ovr = float(roc_auc_score(y_true, probs, multi_class="ovr",
                                          average="macro", labels=list(range(n_cls)))) \
                if len(np.unique(y_true)) > 1 else float("nan")
        except Exception:
            auc_ovr = float("nan")
        try:
            auc_pr = float(average_precision_score(
                np.eye(n_cls)[y_true], probs, average="macro"))
        except Exception:
            auc_pr = float("nan")
        out.update({
            "precision_macro": float(precision_score(y_true, preds, average="macro", zero_division=0)),
            "recall_macro":    float(recall_score(y_true, preds, average="macro", zero_division=0)),
            "macro_f1":        float(f1_score(y_true, preds, average="macro", zero_division=0)),
            "auc_roc_ovr":     auc_ovr,
            "auc_pr_macro":    auc_pr,
        })
    metrics = {k: round(v, 4) if isinstance(v, float) and not np.isnan(v) else v
               for k, v in out.items()}
    return metrics, preds, probs


def compute_diagnostics(y_true: np.ndarray, preds: np.ndarray,
                        num_labels: int, label_names: list) -> dict:
    """Confusion matrix + per-class P/R/F1/support (everything needed to rebuild
    a confusion matrix). Kept as a sibling of test_metrics so 05/06 aggregators
    (which flatten only test_metrics) stay scalar-clean."""
    labels = list(range(num_labels))
    cm = confusion_matrix(y_true, preds, labels=labels)
    pr, rc, f1, sup = precision_recall_fscore_support(
        y_true, preds, labels=labels, zero_division=0)
    total = int(cm.sum())

    def _specificity(i: int) -> float:
        # one-vs-rest TNR: TN / (TN + FP)
        tp = int(cm[i, i])
        fp = int(cm[:, i].sum()) - tp
        fn = int(cm[i, :].sum()) - tp
        tn = total - tp - fp - fn
        return tn / (tn + fp) if (tn + fp) > 0 else 0.0

    per_class = {
        int(i): {
            "name":        label_names[i],
            "precision":   round(float(pr[i]), 4),
            "recall":      round(float(rc[i]), 4),
            # recall == sensitivity (TPR) for this class one-vs-rest
            "sensitivity": round(float(rc[i]), 4),
            "specificity": round(_specificity(i), 4),
            "f1":          round(float(f1[i]), 4),
            "support":     int(sup[i]),
        }
        for i in range(num_labels)
    }
    return {
        "labels": {int(i): label_names[i] for i in range(num_labels)},
        "confusion_matrix": cm.tolist(),  # rows = true, cols = pred
        "per_class": per_class,
    }


# ── Train / eval loops ────────────────────────────────────────────────────────
def run_one_epoch(model, loader, criterion, optimizer, scaler, device,
                  train: bool, grad_accum_steps: int = 1):
    model.train(mode=train)
    total_loss = 0.0
    total_correct = 0
    total_n = 0
    all_logits, all_labels = [], []

    accum = max(1, int(grad_accum_steps)) if train else 1
    n_batches = len(loader)

    ctx = torch.enable_grad() if train else torch.no_grad()
    with ctx:
        if train:
            optimizer.zero_grad(set_to_none=True)
        for i, batch in enumerate(loader):
            x = batch["image"].to(device, non_blocking=True).float()
            y = batch["label"].to(device, non_blocking=True).long()

            with torch.cuda.amp.autocast(enabled=(device.type == "cuda")):
                logits = model(x)
                loss = criterion(logits, y)

            if train:
                # Divide by `accum` so accumulated grads average the
                # micro-batches; step once every `accum` batches (and on the
                # final batch, to flush any trailing partial group).
                scaler.scale(loss / accum).backward()
                if ((i + 1) % accum == 0) or ((i + 1) == n_batches):
                    scaler.unscale_(optimizer)
                    torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                    scaler.step(optimizer)
                    scaler.update()
                    optimizer.zero_grad(set_to_none=True)

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


def lr_multiplier_at(epoch_idx: int, args) -> float:
    """LR schedule as a multiplier in [0, 1] applied to every param group's base
    LR: linear warmup for warmup_epochs, then cosine decay to a small floor.

    Returning a multiplier (not an absolute LR) lets layer-wise LR decay work —
    each param group keeps its own base LR and is scaled by the same schedule.
    """
    eta_min_frac = 5e-6 / max(args.lr, 1e-12)  # floor as a fraction of peak LR
    if epoch_idx < args.warmup_epochs:
        return (epoch_idx + 1) / max(args.warmup_epochs, 1)
    progress = (epoch_idx - args.warmup_epochs) / max(args.epochs - args.warmup_epochs, 1)
    return eta_min_frac + 0.5 * (1.0 - eta_min_frac) * (1 + np.cos(np.pi * progress))


def build_param_groups(model, args) -> list:
    """AdamW parameter groups with layer-wise LR decay (LLRD).

    Each group carries a 'base_lr'; the scheduler scales every group's base_lr by
    the same [0,1] multiplier each epoch. The classifier head trains at the full
    --lr; ViT blocks get lr *= llrd_gamma per step away from the head; the patch
    / position embeddings get the smallest LR. With --llrd_gamma 1.0, or for the
    frozen strategy (only the head is trainable), this collapses to one group.
    """
    depth = len(model.blocks) if hasattr(model, "blocks") else 12
    max_idx = depth + 1  # 0=patch/pos embeds .. depth=last block .. depth+1=head

    def layer_index(name: str) -> int:
        if name.startswith(("patch_embed", "pos_embed", "cls_token")):
            return 0
        m = re.match(r"blocks\.(\d+)\.", name)
        if m:
            return int(m.group(1)) + 1
        return max_idx  # head, final norm, anything else late

    groups: dict = {}
    for name, p in model.named_parameters():
        if not p.requires_grad:
            continue
        base_lr = args.lr * (args.llrd_gamma ** (max_idx - layer_index(name)))
        key = round(base_lr, 12)
        groups.setdefault(
            key, {"params": [], "base_lr": base_lr, "lr": base_lr}
        )["params"].append(p)
    return list(groups.values())


def preflight_check_inputs(splits: dict, sample_k: int = 6) -> None:
    """Fail fast BEFORE training if the preprocessed inputs look wrong.

    For each split: assert it is non-empty and (train/val) has >=2 classes. Then
    load a random sample of volumes and assert each — shape (1,128,128,128),
    all-finite, non-empty brain (>0.8M nonzero voxels), and intensity within
    generous z-scored bands (nonzero mean in [-0.2,0.2], std in [0.75,1.15]).
    Bands are loose on purpose: catch gross failures (raw un-normalised
    intensities, all-zero volumes, NaN, wrong shape), not enforce exact 0/1
    (real volumes sit near 0.03 / 0.92). Raises RuntimeError on any failure.
    """
    import random
    rng = random.Random(0)
    expect = (128, 128, 128)
    for name, df in splits.items():
        if len(df) == 0:
            raise RuntimeError(f"[preflight] split '{name}' is empty.")
        if name in ("train", "val") and df["label"].nunique() < 2:
            raise RuntimeError(
                f"[preflight] split '{name}' has <2 classes "
                f"({sorted(df['label'].unique())}) — cannot train/validate.")
        paths = df["image_path"].tolist()
        for p in rng.sample(paths, min(sample_k, len(paths))):
            if not os.path.isfile(p):
                raise RuntimeError(f"[preflight] missing volume: {p}")
            arr = np.asarray(nib.load(p).get_fdata(dtype=np.float32))
            arr = np.squeeze(arr)
            if tuple(arr.shape) != expect:
                raise RuntimeError(
                    f"[preflight] {p}: shape {arr.shape}, expected {expect}.")
            if not np.isfinite(arr).all():
                raise RuntimeError(f"[preflight] {p}: contains NaN/Inf.")
            nz = arr[arr != 0]
            if nz.size < 800_000:
                raise RuntimeError(
                    f"[preflight] {p}: only {nz.size} nonzero voxels — "
                    f"brain mask looks empty/wrong.")
            m, s = float(nz.mean()), float(nz.std())
            if not (-0.2 <= m <= 0.2 and 0.75 <= s <= 1.15):
                raise RuntimeError(
                    f"[preflight] {p}: nonzero mean={m:.3f} std={s:.3f} outside "
                    f"expected z-scored bands — was NormalizeIntensity applied "
                    f"in 03_prepare_ViT.py?")
    print(f"  [preflight] inputs OK — sampled up to {sample_k}/split; "
          f"shape (1,128,128,128), finite, z-scored.")


# ── Main ───────────────────────────────────────────────────────────────────────
def patch_val_test_metrics(metrics_path, val_metrics, test_metrics,
                           verify=True, tol=0.02):
    """Recompute-mode writer: add a `val_metrics` block (and refresh
    `test_metrics`) on an existing metrics.json. If the recomputed test
    balanced_acc differs from the stored one by > tol, keep the stored
    test_metrics (likely an env/transform mismatch) but still write
    val_metrics, and warn — so a checkpoint recompute never silently changes
    the published test numbers."""
    metrics_path = Path(metrics_path)
    if metrics_path.exists():
        with open(metrics_path) as f:
            d = json.load(f)
    else:
        d = {"config": {}}
    old = (d.get("test_metrics") or {}).get("balanced_acc")
    new = test_metrics.get("balanced_acc")
    if verify and old is not None and new is not None and abs(old - new) > tol:
        print(f"  [WARN] recomputed test bACC {new:.4f} != stored {old:.4f} "
              f"(>|{tol}|); keeping stored test_metrics, writing val_metrics only.")
    else:
        d["test_metrics"] = test_metrics
    d["val_metrics"] = val_metrics
    with open(metrics_path, "w") as f:
        json.dump(d, f, indent=2)
    print(f"  [val_test] patched {metrics_path}")


def main():
    args = parse_args()
    task_cfg = TASK_CONFIG[args.task]

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # scratch runs land under a distinct model slug so 05/06 aggregators treat
    # the from-scratch ViT as a separate model from the MAE-pretrained one.
    # Sized scratch variants (Tiny / Small / Base) get their own slug so the
    # cross-model table can compare them as separate ablations.
    if args.strategy == "scratch":
        model_slug = f"{VIT_SIZE_SLUG[args.vit_size]}_scratch"
    else:
        model_slug = MODEL_SLUG  # full_ft / frozen: pretrained ViT-B
    out_dir = Path(args.out_dir) / model_slug / args.task / f"seed_{args.seed}" / args.strategy
    # --val_test targets an existing run directly: --out_dir IS the run dir.
    if args.val_test:
        out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print(f"  04_supervised_finetuning — {args.task} | seed={args.seed} | {args.strategy}")
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

    # ── Pre-flight: validate the preprocessed inputs before training ─────────
    preflight_check_inputs({"train": train_df, "val": val_df, "test": test_df})

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
    # Pad missing classes with weight=1.0 so the weight tensor always matches
    # the model's num_labels — small-N splits (e.g. T2 bl-only, N≈22 train)
    # may not contain every class.
    # INVARIANT: class weights are computed from train_df, NOT the Dataset. In
    # augment='plus_original' the ConcatDataset duplicates every original
    # uniformly across all classes (orig + K aug branches all iterate the same
    # items), so the class ratio is unchanged — do NOT recompute from the
    # (larger) dataset.
    num_labels = task_cfg["num_labels"]
    classes_present = np.unique(train_df["label"].values).astype(int)
    cw_present = compute_class_weight("balanced", classes=classes_present,
                                      y=train_df["label"].values)
    cw = np.ones(num_labels, dtype=np.float32)
    for c, w in zip(classes_present, cw_present):
        cw[int(c)] = float(w)
    class_weights = torch.tensor(cw, dtype=torch.float, device=device)
    missing = sorted(set(range(num_labels)) - set(classes_present.tolist()))
    if missing:
        print(f"  [WARN] Train split missing classes {missing}; padded with weight=1.0")
    print(f"  Class weights: {dict(enumerate(cw.round(3).tolist()))}")

    # ── Datasets / dataloaders ───────────────────────────────────────────────
    train_loader = DataLoader(build_dataset(train_df, train=True,
                                            augment=args.augment,
                                            aug_copies=args.aug_copies),
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
    size_cfg = VIT_SIZES[args.vit_size]
    model = Vision_Transformer3D(
        img_size=(128, 128, 128), patch_size=16, in_chans=1,
        n_classes=task_cfg["num_labels"],
        embed_dim=size_cfg["embed_dim"], depth=size_cfg["depth"],
        n_heads=size_cfg["n_heads"], mlp_ratio=size_cfg["mlp_ratio"],
        qkv_bias=True, drop_path_rate=args.drop_path_rate, p=0.0,
        attn_p=args.attn_dropout,
        global_avg_pool=False, pos_embed_type="learnable",
        patch_embed_fun="conv3d",
    )
    if args.strategy == "scratch":
        print("  Strategy 'scratch': random initialisation — no pretrained "
              "checkpoint loaded.")
    elif args.val_test:
        # --val_test loads best_model.pt next, which fully overwrites the weights,
        # so the MAE pretrained ckpt is unnecessary here (and need not be present —
        # lets the recompute run locally where only best_model.pt exists).
        print("  --val_test: skipping pretrained load (best_model.pt overwrites it).")
    else:
        model = load_pretrained_checkpoint(model, args.pretrained_ckpt)
    model = model.to(device)

    if args.strategy == "frozen":
        for n, p in model.named_parameters():
            if not n.startswith("head"):
                p.requires_grad = False
        n_train = sum(p.numel() for p in model.parameters() if p.requires_grad)
        print(f"  Frozen strategy: {n_train:,} trainable params (head only)")

    # ── Optimizer / scheduler / scaler / criterion ───────────────────────────
    param_groups = build_param_groups(model, args)
    optimizer = torch.optim.AdamW(param_groups, lr=args.lr,
                                  weight_decay=args.weight_decay,
                                  betas=(0.9, 0.999))
    _blr = [g["base_lr"] for g in param_groups]
    print(f"  Optimizer: AdamW, {len(param_groups)} LR group(s) "
          f"(LLRD gamma={args.llrd_gamma}); base LR "
          f"[{min(_blr):.2e}, {max(_blr):.2e}]")
    scaler = torch.cuda.amp.GradScaler(enabled=(device.type == "cuda"))
    criterion = nn.CrossEntropyLoss(weight=class_weights,
                                    label_smoothing=args.label_smoothing)

    # ── val/test recompute from a saved checkpoint (no training) ──────────────
    if args.val_test:
        bm = out_dir / "best_model.pt"
        if not bm.exists():
            print(f"  [ERROR] --val_test: no best_model.pt at {out_dir}")
            return
        model.load_state_dict(torch.load(bm, map_location=device)["net"])
        _, _, va_logits, va_labels = run_one_epoch(
            model, val_loader, criterion, optimizer, scaler, device, train=False)
        _, _, te_logits, te_labels = run_one_epoch(
            model, test_loader, criterion, optimizer, scaler, device, train=False)
        val_metrics,  _, _ = compute_test_metrics(
            va_labels, va_logits, task_cfg["task_type"])
        test_metrics, _, _ = compute_test_metrics(
            te_labels, te_logits, task_cfg["task_type"])
        print(f"  [val_test] val ={val_metrics}")
        print(f"  [val_test] test={test_metrics}")
        patch_val_test_metrics(out_dir / "metrics.json", val_metrics, test_metrics)
        return

    # ── wandb (gated; no-op without --wandb) ─────────────────────────────────
    wb = init_wandb(args, task_cfg, extra={
        "n_train": int(len(train_df)), "n_val": int(len(val_df)),
        "n_test": int(len(test_df))})

    # ── Resumable per-epoch checkpoint state ─────────────────────────────────
    ckpt_path = out_dir / "last_checkpoint.pt"
    train_log_path = out_dir / "train_log.csv"
    best_metric = -1.0                     # best val balanced accuracy seen
    best_val_loss_at_best = float("inf")   # val_loss at that epoch (tie-break)
    best_epoch = -1
    epochs_since_improve = 0
    log_rows = []
    start_epoch = 0

    if (not args.no_resume) and ckpt_path.exists():
        ck = torch.load(ckpt_path, map_location=device)
        model.load_state_dict(ck["model"])
        optimizer.load_state_dict(ck["optimizer"])
        scaler.load_state_dict(ck["scaler"])
        best_metric = ck.get("best_metric", -1.0)
        best_val_loss_at_best = ck.get("best_val_loss_at_best", float("inf"))
        best_epoch = ck["best_epoch"]
        epochs_since_improve = ck["epochs_since_improve"]
        log_rows = ck["log_rows"]
        start_epoch = ck["epoch"] + 1  # ck["epoch"] = last completed (0-based)
        print(f"  Resuming from epoch {start_epoch + 1} "
              f"(last completed epoch {ck['epoch'] + 1}; "
              f"best val balanced acc={best_metric:.4f}).")

    def _save_checkpoint(ep: int):
        """Atomic per-epoch checkpoint (survives a 12 h SLURM/Colab kill)."""
        state = {
            "model":                 model.state_dict(),
            "optimizer":             optimizer.state_dict(),
            "scaler":                scaler.state_dict(),
            "epoch":                 ep,  # 0-based, just completed
            "best_metric":           best_metric,
            "best_val_loss_at_best": best_val_loss_at_best,
            "best_epoch":            best_epoch,
            "epochs_since_improve":  epochs_since_improve,
            "log_rows":              log_rows,
            "args":                  vars(args),
        }
        tmp = ckpt_path.with_suffix(".pt.tmp")
        torch.save(state, tmp)
        os.replace(tmp, ckpt_path)

    try:
        # ── Training loop ────────────────────────────────────────────────────
        for epoch in range(start_epoch, args.epochs):
            lr_mult = lr_multiplier_at(epoch, args)
            for g in optimizer.param_groups:
                g["lr"] = g["base_lr"] * lr_mult
            cur_lr = args.lr * lr_mult  # head-group LR, for logging

            tr_loss, tr_acc, _, _ = run_one_epoch(
                model, train_loader, criterion, optimizer, scaler, device,
                train=True, grad_accum_steps=args.grad_accum_steps)
            va_loss, va_acc, va_logits, va_labels = run_one_epoch(
                model, val_loader, criterion, optimizer, scaler, device,
                train=False)

            # Balanced accuracy on val — the metric that selects the best
            # checkpoint. val_loss alone, under class imbalance + warmup,
            # selects an untrained epoch-1 model (see the debug plan).
            # Also pull the full val metric set via compute_test_metrics
            # so W&B has val_auc / val_f1 / val_prec / val_recall (binary:
            # + val_sens / val_spec + val_TP/FP/TN/FN) per epoch, mirroring
            # the test-side report.
            _val_full, va_preds, _ = compute_test_metrics(
                va_labels.astype(int), va_logits, task_cfg["task_type"])
            va_bacc = float(_val_full["balanced_acc"])
            val_extra = {
                "val_auc":    float(_val_full.get(
                    "auc_roc", _val_full.get("auc_roc_ovr", float("nan")))),
                "val_auc_pr": float(_val_full.get(
                    "auc_pr",  _val_full.get("auc_pr_macro", float("nan")))),
                "val_f1":     float(_val_full.get(
                    "f1",      _val_full.get("macro_f1", float("nan")))),
                "val_prec":   float(_val_full.get(
                    "precision", _val_full.get("precision_macro", float("nan")))),
                "val_recall": float(_val_full.get(
                    "recall",  _val_full.get("recall_macro", float("nan")))),
            }
            if task_cfg["task_type"] == "binary":
                val_extra["val_sens"] = float(_val_full["sensitivity"])
                val_extra["val_spec"] = float(_val_full["specificity"])
                _cm = confusion_matrix(
                    va_labels.astype(int), va_preds, labels=[0, 1])
                val_extra["val_TN"] = int(_cm[0, 0])
                val_extra["val_FP"] = int(_cm[0, 1])
                val_extra["val_FN"] = int(_cm[1, 0])
                val_extra["val_TP"] = int(_cm[1, 1])

            log_rows.append({"epoch": epoch + 1, "lr": cur_lr,
                             "train_loss": tr_loss, "train_acc": tr_acc,
                             "val_loss": va_loss, "val_acc": va_acc,
                             "val_bacc": va_bacc, **val_extra})
            print(f"  [epoch {epoch+1:>3}/{args.epochs}] "
                  f"lr={cur_lr:.2e}  train_loss={tr_loss:.4f}  train_acc={tr_acc:.4f}  "
                  f"val_loss={va_loss:.4f}  val_acc={va_acc:.4f}  val_bacc={va_bacc:.4f}  "
                  f"val_auc={val_extra['val_auc']:.3f}  val_f1={val_extra['val_f1']:.3f}")

            # Selection: maximise val balanced accuracy, tie-break on lower
            # val_loss (ties at 0.5 are common when a model collapses).
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

            # Incremental flush + atomic checkpoint every epoch (timeout safety)
            pd.DataFrame(log_rows).to_csv(train_log_path, index=False)
            _save_checkpoint(epoch)

            if wb is not None:
                import wandb
                wandb.log({"epoch": epoch + 1, "lr": cur_lr,
                           "train_loss": tr_loss, "train_acc": tr_acc,
                           "val_loss": va_loss, "val_acc": va_acc,
                           "val_bacc": va_bacc,
                           "best_val_balanced_acc": best_metric,
                           **val_extra})

            if epochs_since_improve >= args.patience:
                print(f"  Early stopping at epoch {epoch+1} "
                      f"(no val balanced-acc improvement for {args.patience}).")
                break

        pd.DataFrame(log_rows).to_csv(train_log_path, index=False)

        # ── Final test on best weights ────────────────────────────────────────
        best_state = torch.load(out_dir / "best_model.pt", map_location=device)
        model.load_state_dict(best_state["net"])
        _, _, test_logits, test_labels = run_one_epoch(
            model, test_loader, criterion, optimizer, scaler, device, train=False)
        test_metrics, test_preds, test_probs = compute_test_metrics(
            test_labels, test_logits, task_cfg["task_type"])
        print(f"  Test metrics (image-level): {test_metrics}")

        # Val-split metrics on the SAME best weights — a separate held-out cohort
        # from test (disjoint 80/10/10), so it must be evaluated on its own loader.
        # Saving it here populates the VALIDATION tables without a second pass.
        _, _, val_logits, val_labels = run_one_epoch(
            model, val_loader, criterion, optimizer, scaler, device, train=False)
        val_metrics, _, _ = compute_test_metrics(
            val_labels, val_logits, task_cfg["task_type"])
        print(f"  Val metrics (image-level): {val_metrics}")

        # Confusion matrix + per-class P/R/F1/sensitivity/specificity/support
        label_names = class_names_for(task_cfg)
        test_diagnostics = compute_diagnostics(
            test_labels, test_preds, num_labels, label_names)

        # Per-sample predictions. test_loader is shuffle=False / batch_size=1, so
        # row i of test_logits aligns with row i of test_df — attach the IDs.
        n_te = len(test_labels)
        pred_df = pd.DataFrame({
            "Patient_ID": test_df["Patient_ID"].values[:n_te],
            "bids_ses":   test_df["bids_ses"].values[:n_te],
            "y_true": test_labels.astype(int),
            "y_pred": test_preds.astype(int),
        })
        for c in range(test_probs.shape[1]):
            pred_df[f"prob_{c}"] = test_probs[:, c]
        pred_df.to_csv(out_dir / "test_predictions.csv", index=False)

        # ── Per-subject evaluation ────────────────────────────────────────────
        # In --long mode a subject contributes several scans; aggregate to one
        # prediction per Patient_ID by mean predicted probability. Steadier and
        # more honest than the image-level metric (test images are correlated).
        prob_cols = [f"prob_{c}" for c in range(test_probs.shape[1])]
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

        # ── Save metrics.json (matches clinical pipeline schema) ─────────────
        config = {
            "model_id":              model_slug,
            "pretrained_ckpt":       args.pretrained_ckpt,
            "task":                  args.task,
            "task_description":      task_cfg["description"],
            "session":               cohort_label(args, task_cfg),
            "long_mode":             args.long_mode,
            "max_months":            args.max_months,
            "session_policy":        task_cfg["session_policy"],
            "seed":                  args.seed,
            "strategy":              args.strategy,
            "vit_size":              args.vit_size,
            "augment":               args.augment,
            "aug_copies":            args.aug_copies,
            "epochs":                args.epochs,
            "best_epoch":            best_epoch,
            "best_val_balanced_acc": round(float(best_metric), 4),
            "lr":                    args.lr,
            "weight_decay":          args.weight_decay,
            "batch_size":            args.batch_size,
            "grad_accum_steps":      args.grad_accum_steps,
            "effective_batch_size":  args.batch_size * args.grad_accum_steps,
            "llrd_gamma":            args.llrd_gamma,
            "warmup_epochs":         args.warmup_epochs,
            "patience":              args.patience,
            "drop_path_rate":        args.drop_path_rate,
            "attn_dropout":          args.attn_dropout,
            "label_smoothing":       args.label_smoothing,
            "n_train":               int(len(train_df)),
            "n_val":                 int(len(val_df)),
            "n_test":                int(len(test_df)),
            "n_test_subjects":       int(len(subj)),
            "class_weights":         {i: round(float(w), 4) for i, w in enumerate(cw.tolist())},
            "timestamp":             datetime.now().isoformat(),
        }
        with open(out_dir / "metrics.json", "w") as f:
            json.dump({"config": config, "test_metrics": test_metrics,
                       "val_metrics": val_metrics,
                       "test_metrics_subject": test_metrics_subject,
                       "test_diagnostics": test_diagnostics}, f, indent=2)

        if wb is not None:
            import wandb
            wandb.log({f"test/{k}": v for k, v in test_metrics.items()
                       if isinstance(v, (int, float))})
            wandb.log({f"test_subject/{k}": v
                       for k, v in test_metrics_subject.items()
                       if isinstance(v, (int, float))})
            wandb.log({"test/confusion_matrix": wandb.plot.confusion_matrix(
                y_true=test_labels.astype(int).tolist(),
                preds=test_preds.astype(int).tolist(),
                class_names=label_names)})
            wandb.run.summary["best_val_balanced_acc"] = best_metric
            wandb.run.summary["best_epoch"] = best_epoch
            for k, v in test_metrics.items():
                if isinstance(v, (int, float)):
                    wandb.run.summary[f"test_{k}"] = v
            for k, v in test_metrics_subject.items():
                if isinstance(v, (int, float)):
                    wandb.run.summary[f"test_subject_{k}"] = v
            wandb.run.summary["confusion_matrix"] = test_diagnostics["confusion_matrix"]
            # per_class has integer class keys; wandb's summary encoder only
            # handles string keys in a nested dict — stringify them.
            wandb.run.summary["per_class"] = {
                str(k): v for k, v in test_diagnostics["per_class"].items()}

        print(f"  Saved: {out_dir}/metrics.json")
        print("=" * 70)
    finally:
        if wb is not None:
            import wandb
            wandb.finish()


if __name__ == "__main__":
    main()
