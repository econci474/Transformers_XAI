"""
04_supervised_finetuning_BrainMVP.py
=====================================
Supervised fine-tuning of BrainMVP UniFormer-Small (CVPR 2025 Highlight)
on ADNI T1w MRIs preprocessed by 03_prepare_ViT.py (128^3 @ 1.75mm RAS).

Mirrors the ViT pipeline (04_supervised_finetuning_ViT.py) exactly — same
tasks, splits, augmentation modes, metrics, output format — but swaps the
ViT-B encoder for BrainMVP's pretrained UniFormer-Small encoder.

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
from torch.utils.data import DataLoader

# Reuse ALL shared infrastructure from the ViT pipeline
THIS_DIR = Path(__file__).resolve().parent
MRI_DIR = THIS_DIR.parent
sys.path.insert(0, str(MRI_DIR))

from _vit_recipe.checkpoint import load_pretrained_checkpoint  # noqa: unused but keeps import parity
# Import shared functions from ViT pipeline
import importlib.util
_vit_spec = importlib.util.spec_from_file_location(
    "_vit_pipeline", MRI_DIR / "04_supervised_finetuning_ViT.py")
_vit = importlib.util.module_from_spec(_vit_spec)
# Prevent main() from firing on import
sys.modules["_vit_pipeline"] = _vit
_vit_spec.loader.exec_module(_vit)

# Shared imports from ViT pipeline
TASK_CONFIG = _vit.TASK_CONFIG
cohort_label = _vit.cohort_label
class_names_for = _vit.class_names_for
load_matched_labels = _vit.load_matched_labels
load_split_from_matched = _vit.load_split_from_matched
session_to_months = _vit.session_to_months
build_dataset = _vit.build_dataset
compute_test_metrics = _vit.compute_test_metrics
compute_diagnostics = _vit.compute_diagnostics
run_one_epoch = _vit.run_one_epoch
lr_at = _vit.lr_at
init_wandb = _vit.init_wandb

# BrainMVP UniFormer
from brain_mvp.uniformer_blocks import uniformer_small

warnings.filterwarnings("ignore")

MODEL_SLUG = "BrainMVP_uniformer"


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
    prefix = "encoder.uniformer."
    encoder_state = {}
    for k, v in state.items():
        if k.startswith(prefix):
            encoder_state[k[len(prefix):]] = v

    if not encoder_state:
        # Maybe keys don't have the prefix — try loading directly
        print(f"  [WARN] No '{prefix}' keys found; attempting direct load")
        encoder_state = state

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
    p.add_argument("--vit_inputs_dir", type=str,
                   default=r"D:\ADNI_BIDS_project\derivatives\vit_inputs")
    p.add_argument("--out_dir", type=str, required=True,
                   help="Output root (e.g. /path/to/brain_mvp_uniformer/aug_none)")
    p.add_argument("--epochs", type=int, default=None)
    p.add_argument("--batch_size", type=int, default=4)
    p.add_argument("--lr", type=float, default=None)
    p.add_argument("--weight_decay", type=float, default=1e-4)
    p.add_argument("--warmup_epochs", type=int, default=5)
    p.add_argument("--patience", type=int, default=10)
    p.add_argument("--num_workers", type=int, default=2)
    p.add_argument("--max_subjects", type=int, default=None)
    p.add_argument("--augment", type=str, default="random",
                   choices=["none", "random", "plus_original"])
    p.add_argument("--aug_copies", type=int, default=1)
    p.add_argument("--drop_rate", type=float, default=0.1,
                   help="Dropout before classification head")
    p.add_argument("--label_smoothing", type=float, default=0.0)
    p.add_argument("--no_resume", action="store_true")
    p.add_argument("--wandb", action="store_true")
    p.add_argument("--wandb_project", type=str, default="brainmvp-mri-finetune")
    p.add_argument("--wandb_entity", type=str, default=None)
    p.add_argument("--wandb_run_name", type=str, default=None)
    args = p.parse_args()

    if args.epochs is None:
        args.epochs = 50 if args.strategy == "full_ft" else 70
    if args.lr is None:
        args.lr = 1e-4 if args.strategy == "full_ft" else 1e-3

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


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    task_cfg = TASK_CONFIG[args.task]
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    out_dir = (Path(args.out_dir) / MODEL_SLUG / args.task
               / f"seed_{args.seed}" / args.strategy)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Skip if already done
    if (out_dir / "metrics.json").exists():
        print(f"  metrics.json already exists at {out_dir} — skipping.")
        return

    print("=" * 70)
    print(f"  BrainMVP UniFormer — {args.task} | seed={args.seed} | {args.strategy}")
    print(f"  Device: {device}")
    print(f"  Output: {out_dir}")
    print("=" * 70)

    # ── Load splits (identical to ViT) ────────────────────────────────────────
    seed_dir = Path(args.data_dir) / f"seed_{args.seed}"
    vit_dir = Path(args.vit_inputs_dir)
    matched_df = load_matched_labels(args.matched_labels_csv)

    if args.long_mode is None:
        train_df, n_miss_tr = load_split_from_matched(
            seed_dir / "train.csv", matched_df, task_cfg, vit_dir,
            session_filter=args.session)
        val_df, n_miss_va = load_split_from_matched(
            seed_dir / "val.csv", matched_df, task_cfg, vit_dir,
            session_filter=args.session)
        test_df, n_miss_te = load_split_from_matched(
            seed_dir / "test.csv", matched_df, task_cfg, vit_dir,
            session_filter=args.session)
    else:
        train_df, n_miss_tr = load_split_from_matched(
            seed_dir / "train.csv", matched_df, task_cfg, vit_dir,
            max_months=args.max_months)
        val_df, n_miss_va = load_split_from_matched(
            seed_dir / "val.csv", matched_df, task_cfg, vit_dir,
            max_months=args.max_months)
        test_df, n_miss_te = load_split_from_matched(
            seed_dir / "test.csv", matched_df, task_cfg, vit_dir,
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
    train_loader = DataLoader(
        build_dataset(train_df, train=True, augment=args.augment,
                      aug_copies=args.aug_copies),
        batch_size=args.batch_size, shuffle=True,
        num_workers=args.num_workers, pin_memory=True, drop_last=True)
    val_loader = DataLoader(
        build_dataset(val_df, train=False),
        batch_size=1, shuffle=False,
        num_workers=args.num_workers, pin_memory=True)
    test_loader = DataLoader(
        build_dataset(test_df, train=False),
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

    wb = init_wandb(args, task_cfg, extra={
        "n_train": len(train_df), "n_val": len(val_df),
        "n_test": len(test_df), "model": MODEL_SLUG})

    # ── Resume ────────────────────────────────────────────────────────────────
    ckpt_path = out_dir / "last_checkpoint.pt"
    best_val = float("inf")
    best_epoch = -1
    epochs_since_improve = 0
    log_rows = []
    start_epoch = 0

    if (not args.no_resume) and ckpt_path.exists():
        ck = torch.load(ckpt_path, map_location=device)
        model.load_state_dict(ck["model"])
        optimizer.load_state_dict(ck["optimizer"])
        scaler.load_state_dict(ck["scaler"])
        best_val = ck["best_val"]
        best_epoch = ck["best_epoch"]
        epochs_since_improve = ck["epochs_since_improve"]
        log_rows = ck["log_rows"]
        start_epoch = ck["epoch"] + 1
        print(f"  Resuming from epoch {start_epoch + 1}")

    def _save_checkpoint(ep):
        state = {"model": model.state_dict(), "optimizer": optimizer.state_dict(),
                 "scaler": scaler.state_dict(), "epoch": ep,
                 "best_val": best_val, "best_epoch": best_epoch,
                 "epochs_since_improve": epochs_since_improve,
                 "log_rows": log_rows}
        tmp = ckpt_path.with_suffix(".pt.tmp")
        torch.save(state, tmp)
        os.replace(tmp, ckpt_path)

    try:
        for epoch in range(start_epoch, args.epochs):
            cur_lr = lr_at(epoch, args)
            for g in optimizer.param_groups:
                g["lr"] = cur_lr

            tr_loss, tr_acc, _, _ = run_one_epoch(
                model, train_loader, criterion, optimizer, scaler, device, True)
            va_loss, va_acc, _, _ = run_one_epoch(
                model, val_loader, criterion, optimizer, scaler, device, False)

            log_rows.append({"epoch": epoch + 1, "lr": cur_lr,
                             "train_loss": tr_loss, "train_acc": tr_acc,
                             "val_loss": va_loss, "val_acc": va_acc})
            print(f"  [epoch {epoch+1:>3}/{args.epochs}] "
                  f"lr={cur_lr:.2e}  tr_loss={tr_loss:.4f}  "
                  f"va_loss={va_loss:.4f}  va_acc={va_acc:.4f}")

            if va_loss < best_val:
                best_val = va_loss
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
                           "val_acc": va_acc, "best_val_loss": best_val})

            if epochs_since_improve >= args.patience:
                print(f"  Early stopping at epoch {epoch+1}")
                break

        pd.DataFrame(log_rows).to_csv(out_dir / "train_log.csv", index=False)

        # ── Test ──────────────────────────────────────────────────────────────
        best_state = torch.load(out_dir / "best_model.pt", map_location=device)
        model.load_state_dict(best_state["net"])
        _, _, test_logits, test_labels = run_one_epoch(
            model, test_loader, criterion, optimizer, scaler, device, False)
        test_metrics, test_preds, test_probs = compute_test_metrics(
            test_labels, test_logits, task_cfg["task_type"])
        print(f"  Test metrics: {test_metrics}")

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
                       "test_diagnostics": test_diagnostics}, f, indent=2)
        print(f"  Saved: {out_dir}/metrics.json")
        print("=" * 70)
    finally:
        if wb is not None:
            import wandb
            wandb.finish()


if __name__ == "__main__":
    main()
