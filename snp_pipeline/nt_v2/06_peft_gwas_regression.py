"""
06_peft_gwas_regression.py
===========================
LoRA fine-tuning of NT v2 for GWAS signed z-score regression.
Uses the same input CSVs as BMFM 09_ pipeline.

Input: train.csv / dev.csv / test.csv with columns: sequence, z_score, seq_id
       (from bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom_combos)

Usage:
  python nt_v2/06_peft_gwas_regression.py \
      --data-dir "D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom_combos/forward_and_reverse" \
      --outdir results/ntv2_peft_gwas_regression
"""

import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd
import torch
from scipy.stats import pearsonr, spearmanr
from torch.utils.data import DataLoader
from transformers import AutoTokenizer

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from peft_utils import (
    MODEL_ID, NTv2ForSequenceTask, SequenceLabelDataset,
    train_epoch, eval_epoch, save_lora_checkpoint, save_metrics,
)


def load_csv(csv_path, limit=None):
    df = pd.read_csv(csv_path)
    if limit:
        df = df.head(limit)
    return df["sequence"].tolist(), df["z_score"].values.astype(np.float32).tolist()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True,
                        help="Dir with train.csv, dev.csv, test.csv")
    parser.add_argument("--outdir", default="results/ntv2_peft_gwas_regression")
    parser.add_argument("--epochs", type=int, default=5)
    parser.add_argument("--lr", type=float, default=1e-5)
    parser.add_argument("--batch-size", type=int, default=2)
    parser.add_argument("--grad-accum", type=int, default=32,
                        help="Effective batch = batch_size * grad_accum = 64")
    parser.add_argument("--max-length", type=int, default=2048)
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
    train_seqs, train_labels = load_csv(data_dir / "train.csv", args.limit)
    val_seqs,   val_labels   = load_csv(data_dir / "dev.csv",   args.limit)
    test_seqs,  test_labels  = load_csv(data_dir / "test.csv",  args.limit)
    print(f"  Train: {len(train_labels)},  Val: {len(val_labels)},  Test: {len(test_labels)}")

    train_ds = SequenceLabelDataset(train_seqs, train_labels, tokenizer,
                                     args.max_length, label_dtype=torch.float32)
    val_ds   = SequenceLabelDataset(val_seqs,   val_labels,   tokenizer,
                                     args.max_length, label_dtype=torch.float32)
    test_ds  = SequenceLabelDataset(test_seqs,  test_labels,  tokenizer,
                                     args.max_length, label_dtype=torch.float32)

    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True,
                              num_workers=0, pin_memory=True)
    val_loader   = DataLoader(val_ds,   batch_size=args.batch_size, num_workers=0)
    test_loader  = DataLoader(test_ds,  batch_size=args.batch_size, num_workers=0)

    # ── Model ─────────────────────────────────────────────────────────────────
    print(f"\nBuilding regression model with LoRA (r={args.lora_r})")
    model = NTv2ForSequenceTask(
        num_labels=1, is_regression=True,
        lora_r=args.lora_r, lora_alpha=args.lora_alpha,
    ).to(args.device)

    criterion = torch.nn.MSELoss()
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

        r, _ = pearsonr(val_labels_arr, val_preds) if len(val_preds) > 1 else (0, 1)
        print(f"  Epoch {epoch}  train_MSE={train_loss:.4f}  "
              f"val_MSE={val_loss:.4f}  val_r={r:.4f}")

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
    r_pearson, _ = pearsonr(test_labels_arr, test_preds)
    r_spearman, _ = spearmanr(test_labels_arr, test_preds)
    rmse = float(np.sqrt(test_loss))

    metrics = {
        "test_MSE": float(test_loss),
        "test_RMSE": rmse,
        "pearson_r": float(r_pearson),
        "spearman_r": float(r_spearman),
    }
    print(f"\n  Test metrics: {json.dumps(metrics, indent=2)}")
    save_metrics(metrics, out_dir)
    print(f"  Done → {out_dir}")


if __name__ == "__main__":
    main()
