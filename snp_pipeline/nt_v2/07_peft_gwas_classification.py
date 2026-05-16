"""
07_peft_gwas_classification.py
================================
LoRA fine-tuning of NT v2 for GWAS binned z-score classification (7 classes).
Uses the same input CSVs as BMFM 10_ pipeline (from 08h).

Input: train.csv / dev.csv / test.csv with columns: sequence, z_bin, z_score, seq_id
       z_bin ∈ {-3, -2, -1, 0, +1, +2, +3}

Usage:
  python nt_v2/07_peft_gwas_classification.py \
      --data-dir "D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/bmfm_gwas_classification/forward_and_reverse" \
      --outdir results/ntv2_peft_gwas_classification
"""

import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                              f1_score, classification_report)
from torch.utils.data import DataLoader
from transformers import AutoTokenizer

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from peft_utils import (
    MODEL_ID, NTv2ForSequenceTask, SequenceLabelDataset,
    FocalLoss, train_epoch, eval_epoch, save_lora_checkpoint, save_metrics,
)

NUM_CLASSES = 7
ZBIN_TO_IDX = {-3: 0, -2: 1, -1: 2, 0: 3, 1: 4, 2: 5, 3: 6}
IDX_TO_ZBIN = {v: k for k, v in ZBIN_TO_IDX.items()}


def load_csv(csv_path, limit=None, max_null=None):
    df = pd.read_csv(csv_path)
    if max_null is not None:
        null_mask = df["z_bin"] == 0
        null_df = df[null_mask]
        sig_df = df[~null_mask]
        if len(null_df) > max_null:
            null_df = null_df.sample(n=max_null, random_state=42)
        df = pd.concat([sig_df, null_df]).sample(frac=1, random_state=42)
    if limit:
        df = df.head(limit)
    # Map z_bin to contiguous indices
    labels = df["z_bin"].map(ZBIN_TO_IDX).values.astype(int).tolist()
    return df["sequence"].tolist(), labels


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True,
                        help="Dir with train.csv, dev.csv, test.csv")
    parser.add_argument("--outdir", default="results/ntv2_peft_gwas_classification")
    parser.add_argument("--epochs", type=int, default=5)
    parser.add_argument("--lr", type=float, default=1e-5)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--grad-accum", type=int, default=16,
                        help="Effective batch = batch_size * grad_accum = 64")
    parser.add_argument("--max-length", type=int, default=2048)
    parser.add_argument("--max-null", type=int, default=None,
                        help="Cap null (z_bin=0) rows in training set")
    parser.add_argument("--focal-gamma", type=float, default=2.0,
                        help="Focal loss gamma (0 = standard CE)")
    parser.add_argument("--lora-r", type=int, default=8)
    parser.add_argument("--lora-alpha", type=int, default=16)
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    out_dir = pathlib.Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data_dir = pathlib.Path(args.data_dir)

    # ── Tokenizer ─────────────────────────────────────────────────────────────
    print(f"Loading tokenizer: {MODEL_ID}")
    tokenizer = AutoTokenizer.from_pretrained(MODEL_ID, trust_remote_code=True)

    # ── Datasets ──────────────────────────────────────────────────────────────
    print("Loading datasets...")
    train_seqs, train_labels = load_csv(data_dir / "train.csv", args.limit, args.max_null)
    val_seqs,   val_labels   = load_csv(data_dir / "dev.csv",   args.limit)
    test_seqs,  test_labels  = load_csv(data_dir / "test.csv",  args.limit)

    # Print class distribution
    train_counts = np.bincount(train_labels, minlength=NUM_CLASSES)
    print(f"  Train: {len(train_labels)} rows  class dist: {dict(zip(range(NUM_CLASSES), train_counts.tolist()))}")
    print(f"  Val:   {len(val_labels)},  Test: {len(test_labels)}")

    train_ds = SequenceLabelDataset(train_seqs, train_labels, tokenizer,
                                     args.max_length, label_dtype=torch.long)
    val_ds   = SequenceLabelDataset(val_seqs,   val_labels,   tokenizer,
                                     args.max_length, label_dtype=torch.long)
    test_ds  = SequenceLabelDataset(test_seqs,  test_labels,  tokenizer,
                                     args.max_length, label_dtype=torch.long)

    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True,
                              num_workers=0, pin_memory=True)
    val_loader   = DataLoader(val_ds,   batch_size=args.batch_size, num_workers=0)
    test_loader  = DataLoader(test_ds,  batch_size=args.batch_size, num_workers=0)

    # ── Model ─────────────────────────────────────────────────────────────────
    print(f"\nBuilding {NUM_CLASSES}-class model with LoRA (r={args.lora_r})")
    model = NTv2ForSequenceTask(
        num_labels=NUM_CLASSES, is_regression=False,
        lora_r=args.lora_r, lora_alpha=args.lora_alpha,
    ).to(args.device)

    # Focal loss with class weights
    class_weights = torch.tensor(
        [len(train_labels) / (NUM_CLASSES * max(c, 1)) for c in train_counts],
        dtype=torch.float32
    ).to(args.device)
    print(f"  Class weights: {[f'{w:.2f}' for w in class_weights.tolist()]}")
    criterion = FocalLoss(gamma=args.focal_gamma, weight=class_weights)

    optimizer = torch.optim.AdamW(
        [p for p in model.parameters() if p.requires_grad],
        lr=args.lr, weight_decay=0.01,
    )
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)

    # ── Train ─────────────────────────────────────────────────────────────────
    best_val_loss = float("inf")

    for epoch in range(1, args.epochs + 1):
        train_loss = train_epoch(model, train_loader, optimizer, criterion,
                                  args.device, args.grad_accum)
        val_loss, val_preds, val_labels_arr = eval_epoch(
            model, val_loader, criterion, args.device
        )
        scheduler.step()

        val_acc = accuracy_score(val_labels_arr, val_preds)
        val_bal = balanced_accuracy_score(val_labels_arr, val_preds)
        print(f"  Epoch {epoch}  train_loss={train_loss:.4f}  val_loss={val_loss:.4f}  "
              f"val_acc={val_acc:.3f}  val_bal_acc={val_bal:.3f}")

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            save_lora_checkpoint(model, out_dir, "best")
            torch.save(model.head.state_dict(), out_dir / "best_head.pt")

    # ── Evaluate ──────────────────────────────────────────────────────────────
    print("\nEvaluating on test set...")
    model.backbone.load_adapter(out_dir / "best_lora", adapter_name="default")
    model.head.load_state_dict(torch.load(out_dir / "best_head.pt", weights_only=True))

    test_loss, test_preds, test_labels_arr = eval_epoch(
        model, test_loader, criterion, args.device
    )

    # Map back to z-bin labels for readability
    zbin_preds = [IDX_TO_ZBIN[int(p)] for p in test_preds]
    zbin_labels = [IDX_TO_ZBIN[int(l)] for l in test_labels_arr]

    metrics = {
        "test_loss": float(test_loss),
        "accuracy": float(accuracy_score(test_labels_arr, test_preds)),
        "balanced_accuracy": float(balanced_accuracy_score(test_labels_arr, test_preds)),
        "f1_weighted": float(f1_score(test_labels_arr, test_preds, average="weighted")),
        "f1_macro": float(f1_score(test_labels_arr, test_preds, average="macro")),
    }
    print(f"\n  Test metrics: {json.dumps(metrics, indent=2)}")
    target_names = [f"z_bin={IDX_TO_ZBIN[i]}" for i in range(NUM_CLASSES)]
    print(classification_report(test_labels_arr, test_preds, target_names=target_names))
    save_metrics(metrics, out_dir)
    print(f"  Done → {out_dir}")


if __name__ == "__main__":
    main()
