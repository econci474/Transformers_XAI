"""
04_supervised_finetuning_BrainMVP.py
=====================================
Supervised fine-tuning of BrainMVP UniFormer-Small (CVPR 2025 Highlight)
on ADNI T1w MRIs preprocessed by 03_prepare_BrainMVP.py (128×128×64 @ 1mm).

Input images are 128×128×64 NIfTI volumes. Following the BrainMVP paper
(Appendix I, Classification), a spatial crop to 96×96×64 is applied:
  - Training: RandSpatialCrop (data augmentation)
  - Inference: CenterSpatialCrop (deterministic)

Augmentation arms (--augment)
-----------------------------
  "none"           → CenterSpatialCrop 96×96×64 only (paper inference spec).
                     No augmentation. Baseline control.
  "stochastic"     → RandSpatialCrop 96×96×64 + stochastic flips/intensity
                     (p=0.5). Each epoch sees a different random crop + maybe
                     augmented version. This is the default.
  "plus_original"  → ConcatDataset of center-cropped originals (no aug) PLUS
                     aug_copies of (RandSpatialCrop + always-on aug with p=1.0).
                     Train set size = (1 + aug_copies) × N.

Architecture
------------
  uniformer_small(in_chans=1) → global avg pool of stage-4 features (512-d)
  → nn.Linear(512, n_classes) classification head.

Checkpoint
----------
  The BrainMVP pretrained checkpoint contains the full RecModel (encoder +
  decoder + template). We extract only encoder.uniformer.* weights.
"""

import argparse, json, os, re, sys, warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.utils.class_weight import compute_class_weight
from sklearn.metrics import balanced_accuracy_score
from torch.utils.data import DataLoader, ConcatDataset
import monai.transforms as mt
from monai.data import Dataset as MonaiDataset

# Reuse ALL shared infrastructure from the ViT pipeline
THIS_DIR = Path(__file__).resolve().parent
MRI_DIR = THIS_DIR.parent
sys.path.insert(0, str(MRI_DIR))

from _vit_recipe.checkpoint import load_pretrained_checkpoint  # noqa: unused but keeps import parity
# Import shared functions from ViT pipeline
import importlib.util
_vit_spec = importlib.util.spec_from_file_location(
    "_vit_pipeline", MRI_DIR / "04_supervised_finetuning.py")
_vit = importlib.util.module_from_spec(_vit_spec)
# Prevent main() from firing on import
sys.modules["_vit_pipeline"] = _vit
_vit_spec.loader.exec_module(_vit)

# Shared imports from ViT pipeline
TASK_CONFIG = _vit.TASK_CONFIG
cohort_label = _vit.cohort_label
class_names_for = _vit.class_names_for
load_matched_labels = _vit.load_matched_labels
session_to_months = _vit.session_to_months
patient_id_to_bids_sub = _vit.patient_id_to_bids_sub
compute_test_metrics = _vit.compute_test_metrics
compute_diagnostics = _vit.compute_diagnostics
run_one_epoch = _vit.run_one_epoch
lr_multiplier_at = _vit.lr_multiplier_at
# init_wandb is NOT imported from _vit: the ViT version references ViT-only
# argparse args (grad_accum_steps, llrd_gamma, drop_path_rate, attn_dropout).
# A local init_wandb is defined below, using only this script's own args.

# ── BrainMVP spatial crop dimensions ──────────────────────────────────────────
BMVP_CROP = (96, 96, 64)   # paper Appendix I: classification ROI

# BrainMVP UniFormer
from brain_mvp.uniformer_blocks import uniformer_small

warnings.filterwarnings("ignore")

MODEL_SLUG = "BrainMVP_uniformer"


def init_wandb(args, task_cfg, extra=None):
    """Start a wandb run if --wandb is set; otherwise return None (full no-op).

    Local to the BrainMVP pipeline — the wandb config dict uses only this
    script's own argparse args, so it cannot break when the ViT script's args
    change (see project_brainmvp_imports_vit_symbols). Honours WANDB_MODE /
    WANDB_DIR from the environment; degrades to local-only if wandb is absent.
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
        "augment":         args.augment,
        "aug_copies":      args.aug_copies,
        "epochs":          args.epochs,
        "lr":              args.lr,
        "weight_decay":    args.weight_decay,
        "batch_size":      args.batch_size,
        "warmup_epochs":   args.warmup_epochs,
        "patience":        args.patience,
        "drop_rate":       args.drop_rate,
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


# ── BrainMVP-specific transforms ──────────────────────────────────────────────
def _bmvp_eval_transform() -> mt.Compose:
    """Inference: load → center crop 96×96×64 (paper spec)."""
    return mt.Compose([
        mt.LoadImaged(keys=["image"]),
        mt.EnsureChannelFirstd(keys=["image"]),
        mt.CenterSpatialCropd(keys=["image"], roi_size=BMVP_CROP),
        mt.ToTensord(keys=["image"]),
    ])


def _bmvp_stochastic_transform(p_flip=0.5, p_scale=0.5, p_shift=0.5) -> mt.Compose:
    """Training: random crop 96×96×64 + stochastic augmentation."""
    return mt.Compose([
        mt.LoadImaged(keys=["image"]),
        mt.EnsureChannelFirstd(keys=["image"]),
        mt.RandSpatialCropd(keys=["image"], roi_size=BMVP_CROP, random_center=True),
        mt.RandFlipd(keys=["image"], prob=p_flip, spatial_axis=0),
        mt.RandFlipd(keys=["image"], prob=p_flip, spatial_axis=1),
        mt.RandFlipd(keys=["image"], prob=p_flip, spatial_axis=2),
        mt.RandScaleIntensityd(keys=["image"], factors=0.1, prob=p_scale),
        mt.RandShiftIntensityd(keys=["image"], offsets=0.1, prob=p_shift),
        mt.ToTensord(keys=["image"]),
    ])


def build_bmvp_dataset(df, train: bool,
                       augment: str = "stochastic", aug_copies: int = 1):
    """
    BrainMVP dataset with spatial crops per the paper.

    Input volumes: 128×128×64 (from 03_prepare_BrainMVP.py).
    Model input  : 96×96×64  (cropped here).

    augment (train only; val/test always use CenterSpatialCrop):
      "none"          → CenterSpatialCrop only, no augmentation (baseline)
      "stochastic"    → RandSpatialCrop + stochastic aug (p=0.5 flips, etc.)
      "plus_original" → ConcatDataset of center-cropped originals (no aug)
                        + aug_copies of (RandSpatialCrop + always-on aug).
                        Train set size = (1 + aug_copies) × N.
    """
    items = [{"image": row.image_path, "label": int(row.label)}
             for row in df.itertuples()]

    if not train or augment == "none":
        return MonaiDataset(data=items, transform=_bmvp_eval_transform())
    if augment == "stochastic":
        return MonaiDataset(data=items, transform=_bmvp_stochastic_transform())
    if augment == "plus_original":
        k = max(1, int(aug_copies))
        orig = MonaiDataset(data=items, transform=_bmvp_eval_transform())
        aug = [MonaiDataset(data=items,
                            transform=_bmvp_stochastic_transform(
                                p_flip=1.0, p_scale=1.0, p_shift=1.0))
               for _ in range(k)]
        return ConcatDataset([orig] + aug)
    raise ValueError(f"Unknown augment mode: {augment!r}")


def _resolve_brainmvp_path(patient_id: str, viscode: str, inputs_dir: Path) -> str:
    """Resolve image path for BrainMVP 128x64 preprocessed NIfTIs."""
    sub = patient_id_to_bids_sub(patient_id)
    ses_label = f"ses-{viscode}"
    return str(inputs_dir / sub / ses_label /
               f"{sub}_{ses_label}_space-BrainMVP128x64_desc-preproc_T1w.nii.gz")


def load_split_brainmvp(baseline_csv, matched_df, task_cfg, inputs_dir,
                        session_filter=None, max_months=None):
    """Same as ViT load_split_from_matched but uses BrainMVP96 file paths."""
    # Get the labelled rows from the ViT pipeline loader (which checks file existence)
    # We can't reuse it directly since it hardcodes ViT128 filenames,
    # so we replicate the logic with our path resolver.
    import re as _re
    from bidsification.exclusions import is_excluded_subject

    subjects = pd.read_csv(baseline_csv, usecols=["Patient_ID"])["Patient_ID"].unique().tolist()
    df = matched_df[matched_df["Patient_ID"].isin(subjects)].copy()

    if session_filter is not None:
        df = df[df["bids_ses"] == session_filter].copy()
    elif max_months is not None:
        df["_months"] = df["bids_ses"].map(session_to_months)
        df = df.dropna(subset=["_months"])
        df = df[df["_months"] <= max_months].drop(columns=["_months"])

    if task_cfg["session_policy"] == "baseline_anchored":
        cap = int(task_cfg.get("label_anchor_max_months", 0))
        df["_months"] = df["bids_ses"].map(session_to_months)
        df = df.dropna(subset=["_months"])
        df = df[df["_months"] <= cap].drop(columns=["_months"]).copy()
    elif task_cfg["session_policy"] == "baseline_only":
        df = df[df["bids_ses"] == "bl"].copy()

    if task_cfg["filter_non_ad"]:
        df = df[df["Label_bl_multi"] != "AD"].copy()

    df = df[~df["Patient_ID"].apply(is_excluded_subject)].copy()
    df = _vit._resolve_labels(df, task_cfg)

    df["image_path"] = df.apply(
        lambda r: _resolve_brainmvp_path(r["Patient_ID"], r["bids_ses"], inputs_dir),
        axis=1)
    have = df["image_path"].apply(os.path.isfile)
    n_missing = int((~have).sum())
    df = df[have].copy()

    out_cols = ["Patient_ID", "bids_ses", "label", "image_path"]
    return df.reset_index(drop=True)[out_cols], n_missing


# ── BrainMVP classifier wrapper ──────────────────────────────────────────────
class BrainMVPClassifier(nn.Module):
    """UniFormer-Small encoder + global avg pool + linear head."""

    def __init__(self, n_classes: int, drop_rate: float = 0.1):
        super().__init__()
        self.encoder = uniformer_small(in_chans=1)
        self.pool = nn.AdaptiveAvgPool3d(1)
        self.dropout = nn.Dropout(drop_rate)
        self.head = nn.Linear(512, n_classes)

    def forward(self, x):
        # encoder returns (x0, x1, x2, x3, x4); we use deepest (512-ch)
        features = self.encoder(x)
        x4 = features[-1]                      # (B, 512, D', H', W')
        pooled = self.pool(x4).flatten(1)       # (B, 512)
        return self.head(self.dropout(pooled))


def load_brainmvp_pretrained(model: BrainMVPClassifier, ckpt_path: str):
    """Load BrainMVP pretrained weights into the encoder only.
    
    The checkpoint contains the full RecModel state dict with keys like:
      encoder.uniformer.patch_embed1.proj.weight
      decoder.decoder5.transp_conv.conv.weight
      rep_template
    We extract only 'encoder.uniformer.*' keys, strip that prefix,
    and load into model.encoder (which IS the uniformer_small).
    """
    print(f"  Loading BrainMVP checkpoint: {ckpt_path}")
    ckpt = torch.load(ckpt_path, map_location="cpu")

    # Handle different checkpoint formats
    if isinstance(ckpt, dict):
        state = ckpt.get("state_dict", ckpt.get("model", ckpt.get("net", ckpt)))
    else:
        state = ckpt

    # Extract encoder.uniformer.* keys
    # Handle 'module.' prefix from DataParallel wrapping
    prefix_candidates = ["encoder.uniformer.", "module.encoder.uniformer."]
    encoder_state = {}
    matched_prefix = None
    for prefix in prefix_candidates:
        for k, v in state.items():
            if k.startswith(prefix):
                encoder_state[k[len(prefix):]] = v
                matched_prefix = prefix
        if encoder_state:
            break

    if not encoder_state:
        # Maybe keys don't have any prefix — try loading directly
        print(f"  [WARN] No encoder.uniformer keys found; attempting direct load")
        encoder_state = state
    else:
        print(f"  Matched prefix: '{matched_prefix}'")

    loaded = model.encoder.load_state_dict(encoder_state, strict=False)
    n_loaded = len(encoder_state) - len(loaded.unexpected_keys)
    print(f"  Loaded {n_loaded} encoder weight tensors")
    if loaded.missing_keys:
        print(f"  [WARN] Missing keys: {loaded.missing_keys[:5]}...")
    if loaded.unexpected_keys:
        print(f"  [WARN] Unexpected keys: {loaded.unexpected_keys[:5]}...")
    return model


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="BrainMVP UniFormer fine-tuning on ADNI MRI")
    p.add_argument("--task", type=str, required=True,
                   choices=list(TASK_CONFIG.keys()))
    p.add_argument("--seed", type=int, required=True, choices=[0, 1, 2])
    p.add_argument("--strategy", type=str, default="full_ft",
                   choices=["full_ft", "frozen"])
    p.add_argument("--pretrained_ckpt", type=str, required=True,
                   help="Path to BrainMVP_uniformer.pt")
    session_group = p.add_mutually_exclusive_group()
    session_group.add_argument("--session", type=str, default="bl")
    session_group.add_argument("--long", type=str, default=None)
    p.add_argument("--matched_labels_csv", type=str,
                   default=r"D:\ADNI_BIDS_project\derivatives\mri_clinical_matched"
                           r"\master_mri_clinical_matched_labels.csv")
    p.add_argument("--data_dir", type=str,
                   default=r"D:\ADNI_BIDS_project\derivatives\clinical\verbose\baseline")
    p.add_argument("--brainmvp_inputs_dir", type=str,
                   default=r"D:\ADNI_BIDS_project\derivatives\brainmvp_inputs")
    p.add_argument("--out_dir", type=str, required=True,
                   help="Output root (e.g. /path/to/brain_mvp_uniformer/aug_none)")
    p.add_argument("--epochs", type=int, default=None)
    p.add_argument("--batch_size", type=int, default=4)
    p.add_argument("--lr", type=float, default=None)
    p.add_argument("--weight_decay", type=float, default=1e-4)
    p.add_argument("--warmup_epochs", type=int, default=5)
    p.add_argument("--patience", type=int, default=50)
    p.add_argument("--num_workers", type=int, default=0,
                   help="DataLoader workers (0=main process; avoids bus errors)")
    p.add_argument("--max_subjects", type=int, default=None)
    p.add_argument("--augment", type=str, default="stochastic",
                   choices=["none", "stochastic", "plus_original"],
                   help="none = center crop only; "
                        "stochastic = random crop + p=0.5 aug; "
                        "plus_original = center crop originals + K always-on aug copies")
    p.add_argument("--aug_copies", type=int, default=1)
    p.add_argument("--drop_rate", type=float, default=0.1,
                   help="Dropout before classification head")
    p.add_argument("--label_smoothing", type=float, default=0.0)
    p.add_argument("--no_resume", action="store_true")
    p.add_argument("--val_test", action="store_true",
                   help="No training: load best_model.pt and recompute val + "
                        "test metrics, patching a `val_metrics` block into the "
                        "existing metrics.json (adds val AUC/F1 to runs that "
                        "only logged val_bacc). Requires best_model.pt present.")
    p.add_argument("--wandb", action="store_true")
    p.add_argument("--wandb_project", type=str, default="brainmvp-mri-finetune")
    p.add_argument("--wandb_entity", type=str, default=None)
    p.add_argument("--wandb_run_name", type=str, default=None)
    args = p.parse_args()

    # BrainMVP paper Appendix I: 200 epochs, lr=3e-4, cosine schedule
    if args.epochs is None:
        args.epochs = 200 if args.strategy == "full_ft" else 100
    if args.lr is None:
        # full_ft: 5e-5 -- a standard foundation-model full-fine-tune LR.
        # The paper's 3e-4 diverged to NaN (~epoch 27) on every full_ft run.
        args.lr = 5e-5 if args.strategy == "full_ft" else 1e-3

    # Resolve session selection mode (same logic as ViT)
    if args.long is None:
        args.long_mode = None
        args.max_months = None
    else:
        if args.long.lower() == "all" or args.long == "0":
            args.long_mode = "all"
            args.max_months = None
        else:
            n = int(args.long)
            args.long_mode = "cutoff"
            args.max_months = n * 12
        args.session = None

    # Compatibility attrs for init_wandb / cohort_label
    args.drop_path_rate = 0.0
    args.attn_dropout = 0.0
    return args


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


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    task_cfg = TASK_CONFIG[args.task]
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    out_dir = (Path(args.out_dir) / MODEL_SLUG / args.task
               / f"seed_{args.seed}" / args.strategy)
    # --val_test targets an existing run directly: --out_dir IS the run dir
    # (avoids slug/nesting reconstruction mismatches; the driver passes the
    # dir holding best_model.pt + metrics.json).
    if args.val_test:
        out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Skip if already done (but --val_test deliberately re-reads an existing run)
    if (out_dir / "metrics.json").exists() and not args.val_test:
        print(f"  metrics.json already exists at {out_dir} — skipping.")
        return

    print("=" * 70)
    print(f"  BrainMVP UniFormer — {args.task} | seed={args.seed} | {args.strategy}")
    print(f"  Device: {device}")
    print(f"  Output: {out_dir}")
    print("=" * 70)

    # ── Load splits (identical to ViT) ────────────────────────────────────────
    seed_dir = Path(args.data_dir) / f"seed_{args.seed}"
    mvp_dir = Path(args.brainmvp_inputs_dir)
    matched_df = load_matched_labels(args.matched_labels_csv)

    if args.long_mode is None:
        train_df, n_miss_tr = load_split_brainmvp(
            seed_dir / "train.csv", matched_df, task_cfg, mvp_dir,
            session_filter=args.session)
        val_df, n_miss_va = load_split_brainmvp(
            seed_dir / "val.csv", matched_df, task_cfg, mvp_dir,
            session_filter=args.session)
        test_df, n_miss_te = load_split_brainmvp(
            seed_dir / "test.csv", matched_df, task_cfg, mvp_dir,
            session_filter=args.session)
    else:
        train_df, n_miss_tr = load_split_brainmvp(
            seed_dir / "train.csv", matched_df, task_cfg, mvp_dir,
            max_months=args.max_months)
        val_df, n_miss_va = load_split_brainmvp(
            seed_dir / "val.csv", matched_df, task_cfg, mvp_dir,
            max_months=args.max_months)
        test_df, n_miss_te = load_split_brainmvp(
            seed_dir / "test.csv", matched_df, task_cfg, mvp_dir,
            max_months=args.max_months)

    if args.max_subjects:
        train_df = train_df.head(args.max_subjects)
        val_df = val_df.head(min(args.max_subjects, len(val_df)))
        test_df = test_df.head(min(args.max_subjects, len(test_df)))

    print(f"  Splits — train: {len(train_df)}  val: {len(val_df)}  "
          f"test: {len(test_df)}")

    # Save manifest
    manifest = pd.concat([
        train_df.assign(split="train"), val_df.assign(split="val"),
        test_df.assign(split="test")], ignore_index=True)
    manifest.to_csv(out_dir / "dataset_manifest.csv", index=False)

    # ── Class weights ─────────────────────────────────────────────────────────
    num_labels = task_cfg["num_labels"]
    classes_present = np.unique(train_df["label"].values).astype(int)
    cw_present = compute_class_weight("balanced", classes=classes_present,
                                      y=train_df["label"].values)
    cw = np.ones(num_labels, dtype=np.float32)
    for c, w in zip(classes_present, cw_present):
        cw[int(c)] = float(w)
    class_weights = torch.tensor(cw, dtype=torch.float, device=device)
    print(f"  Class weights: {dict(enumerate(cw.round(3).tolist()))}")

    # ── Dataloaders ───────────────────────────────────────────────────────────
    train_ds = build_bmvp_dataset(train_df, train=True, augment=args.augment,
                                  aug_copies=args.aug_copies)
    print(f"  Train dataset size: {len(train_ds)} "
          f"(augment={args.augment}, aug_copies={args.aug_copies})")
    train_loader = DataLoader(
        train_ds, batch_size=args.batch_size, shuffle=True,
        num_workers=args.num_workers, pin_memory=True, drop_last=True)
    val_loader = DataLoader(
        build_bmvp_dataset(val_df, train=False),
        batch_size=1, shuffle=False,
        num_workers=args.num_workers, pin_memory=True)
    test_loader = DataLoader(
        build_bmvp_dataset(test_df, train=False),
        batch_size=1, shuffle=False,
        num_workers=args.num_workers, pin_memory=True)

    # ── Model ─────────────────────────────────────────────────────────────────
    model = BrainMVPClassifier(n_classes=num_labels, drop_rate=args.drop_rate)
    model = load_brainmvp_pretrained(model, args.pretrained_ckpt)
    model = model.to(device)

    if args.strategy == "frozen":
        for n, p in model.named_parameters():
            if not n.startswith("head"):
                p.requires_grad = False
        n_train = sum(p.numel() for p in model.parameters() if p.requires_grad)
        print(f"  Frozen: {n_train:,} trainable params (head only)")

    # ── Optimizer / criterion ─────────────────────────────────────────────────
    trainable = [p for p in model.parameters() if p.requires_grad]
    optimizer = torch.optim.AdamW(trainable, lr=args.lr,
                                  weight_decay=args.weight_decay)
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
            model, val_loader, criterion, optimizer, scaler, device, False)
        _, _, te_logits, te_labels = run_one_epoch(
            model, test_loader, criterion, optimizer, scaler, device, False)
        val_metrics, _, _ = compute_test_metrics(
            va_labels, va_logits, task_cfg["task_type"])
        test_metrics, _, _ = compute_test_metrics(
            te_labels, te_logits, task_cfg["task_type"])
        print(f"  [val_test] val ={val_metrics}")
        print(f"  [val_test] test={test_metrics}")
        patch_val_test_metrics(out_dir / "metrics.json", val_metrics, test_metrics)
        return

    wb = init_wandb(args, task_cfg, extra={
        "n_train": len(train_df), "n_val": len(val_df),
        "n_test": len(test_df), "model": MODEL_SLUG})

    # ── Resume ────────────────────────────────────────────────────────────────
    ckpt_path = out_dir / "last_checkpoint.pt"
    best_metric = -1.0                    # best val balanced accuracy
    best_val_loss_at_best = float("inf")  # tie-breaker
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
        best_epoch = ck.get("best_epoch", -1)
        epochs_since_improve = ck.get("epochs_since_improve", 0)
        log_rows = ck.get("log_rows", [])
        start_epoch = ck["epoch"] + 1
        print(f"  Resuming from epoch {start_epoch + 1}")

    def _save_checkpoint(ep):
        state = {"model": model.state_dict(), "optimizer": optimizer.state_dict(),
                 "scaler": scaler.state_dict(), "epoch": ep,
                 "best_metric": best_metric,
                 "best_val_loss_at_best": best_val_loss_at_best,
                 "best_epoch": best_epoch,
                 "epochs_since_improve": epochs_since_improve,
                 "log_rows": log_rows}
        tmp = ckpt_path.with_suffix(".pt.tmp")
        torch.save(state, tmp)
        os.replace(tmp, ckpt_path)

    try:
        for epoch in range(start_epoch, args.epochs):
            # lr_multiplier_at returns a 0-1 multiplier; scale by the peak LR.
            cur_lr = args.lr * lr_multiplier_at(epoch, args)
            for g in optimizer.param_groups:
                g["lr"] = cur_lr

            tr_loss, tr_acc, _, _ = run_one_epoch(
                model, train_loader, criterion, optimizer, scaler, device, True)
            va_loss, va_acc, va_logits, va_labels = run_one_epoch(
                model, val_loader, criterion, optimizer, scaler, device, False)

            # Abort on a non-finite loss: once weights are NaN they stay NaN, so
            # there is nothing to recover -- stop and keep the last good
            # best_model.pt (full_ft at the old lr=3e-4 diverged here ~epoch 27).
            if not (np.isfinite(tr_loss) and np.isfinite(va_loss)):
                print(f"  [ABORT] non-finite loss at epoch {epoch+1} "
                      f"(tr_loss={tr_loss}, va_loss={va_loss}) -- stopping.")
                break

            va_bacc = float(balanced_accuracy_score(
                va_labels.astype(int), va_logits.argmax(axis=-1)))

            log_rows.append({"epoch": epoch + 1, "lr": cur_lr,
                             "train_loss": tr_loss, "train_acc": tr_acc,
                             "val_loss": va_loss, "val_acc": va_acc,
                             "val_bacc": va_bacc})
            print(f"  [epoch {epoch+1:>3}/{args.epochs}] "
                  f"lr={cur_lr:.2e}  tr_loss={tr_loss:.4f}  "
                  f"va_loss={va_loss:.4f}  va_acc={va_acc:.4f}  "
                  f"va_bacc={va_bacc:.4f}")

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

            pd.DataFrame(log_rows).to_csv(out_dir / "train_log.csv", index=False)
            _save_checkpoint(epoch)

            if wb is not None:
                import wandb
                wandb.log({"epoch": epoch + 1, "lr": cur_lr,
                           "train_loss": tr_loss, "val_loss": va_loss,
                           "val_acc": va_acc, "val_bacc": va_bacc,
                           "best_val_balanced_acc": best_metric})

            if epochs_since_improve >= args.patience:
                print(f"  Early stopping at epoch {epoch+1}")
                break

        pd.DataFrame(log_rows).to_csv(out_dir / "train_log.csv", index=False)

        # ── Test ──────────────────────────────────────────────────────────────
        if not (out_dir / "best_model.pt").exists():
            print("  [ERROR] no best_model.pt -- training produced no usable "
                  "epoch (NaN abort before the first checkpoint). "
                  "No metrics.json written.")
            return
        best_state = torch.load(out_dir / "best_model.pt", map_location=device)
        model.load_state_dict(best_state["net"])
        _, _, test_logits, test_labels = run_one_epoch(
            model, test_loader, criterion, optimizer, scaler, device, False)
        test_metrics, test_preds, test_probs = compute_test_metrics(
            test_labels, test_logits, task_cfg["task_type"])
        print(f"  Test metrics: {test_metrics}")

        # Val-split metrics on the SAME best weights — a separate held-out cohort
        # from test (disjoint 80/10/10), evaluated on its own loader. Saved here so
        # the VALIDATION tables fill without a second --val_test pass.
        _, _, val_logits, val_labels = run_one_epoch(
            model, val_loader, criterion, optimizer, scaler, device, False)
        val_metrics, _, _ = compute_test_metrics(
            val_labels, val_logits, task_cfg["task_type"])
        print(f"  Val metrics: {val_metrics}")

        label_names = class_names_for(task_cfg)
        test_diagnostics = compute_diagnostics(
            test_labels, test_preds, num_labels, label_names)

        pred_df = pd.DataFrame({"y_true": test_labels.astype(int),
                                 "y_pred": test_preds.astype(int)})
        for c in range(test_probs.shape[1]):
            pred_df[f"prob_{c}"] = test_probs[:, c]
        pred_df.to_csv(out_dir / "test_predictions.csv", index=False)

        config = {
            "model_id": MODEL_SLUG,
            "pretrained_ckpt": args.pretrained_ckpt,
            "task": args.task,
            "task_description": task_cfg["description"],
            "session": cohort_label(args, task_cfg),
            "long_mode": args.long_mode,
            "max_months": args.max_months,
            "session_policy": task_cfg["session_policy"],
            "seed": args.seed,
            "strategy": args.strategy,
            "augment": args.augment,
            "aug_copies": args.aug_copies,
            "epochs": args.epochs,
            "best_epoch": best_epoch,
            "best_val_balanced_acc": round(float(best_metric), 4),
            "lr": args.lr,
            "weight_decay": args.weight_decay,
            "batch_size": args.batch_size,
            "warmup_epochs": args.warmup_epochs,
            "patience": args.patience,
            "drop_rate": args.drop_rate,
            "label_smoothing": args.label_smoothing,
            "n_train": len(train_df),
            "n_val": len(val_df),
            "n_test": len(test_df),
            "class_weights": {i: round(float(w), 4) for i, w in enumerate(cw)},
            "timestamp": datetime.now().isoformat(),
        }
        with open(out_dir / "metrics.json", "w") as f:
            json.dump({"config": config, "test_metrics": test_metrics,
                       "val_metrics": val_metrics,
                       "test_diagnostics": test_diagnostics}, f, indent=2)
        print(f"  Saved: {out_dir}/metrics.json")
        print("=" * 70)
    finally:
        if wb is not None:
            import wandb
            wandb.finish()


if __name__ == "__main__":
    main()
