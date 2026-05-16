"""
05_peft_patient_classification.py
==================================
LoRA fine-tuning of NT v2 on patient DNA sequences for "ever converter"
classification (binary: 0/1).

Input:
  - Patient window sequences CSV (PTID, window_id, sequence)
  - Clinical labels CSV (Patient_ID, Label_ever_convert)
  - Uses no-CDR stratified train/val/test splits

Usage:
  python nt_v2/05_peft_patient_classification.py \
      --sequences  D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/patient_window_sequences/all_patients.csv \
      --labels-dir "D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified_ever_convert/tabular/baseline" \
      --seed 0 --outdir results/ntv2_peft_patient_clf
"""

import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                              f1_score, roc_auc_score, classification_report)
from torch.utils.data import DataLoader
from transformers import AutoTokenizer

# Add parent to path so we can import peft_utils
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from peft_utils import (
    MODEL_ID, NTv2ForSequenceTask, SequenceLabelDataset,
    train_epoch, eval_epoch, save_lora_checkpoint, save_metrics,
)


def build_patient_dataset(sequences_csv, labels_csv, tokenizer, max_length,
                           task_col="Label_ever_convert"):
    """Join patient sequences with labels. Mean-pool windows → one row per patient
    is NOT done here; instead we replicate the patient label for each window."""
    seq_df = pd.read_csv(sequences_csv)
    lab_df = pd.read_csv(labels_csv, usecols=["Patient_ID", task_col])

    # Rename for merge
    lab_df = lab_df.rename(columns={"Patient_ID": "PTID"})
    merged = seq_df.merge(lab_df, on="PTID", how="inner")
    print(f"    {len(merged)} rows ({merged['PTID'].nunique()} patients, "
          f"pos={merged[task_col].sum()})")

    sequences = merged["sequence"].tolist()
    labels = merged[task_col].values.astype(int).tolist()
    return SequenceLabelDataset(sequences, labels, tokenizer, max_length)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--labels-dir", required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--task", default="Label_ever_convert")
    parser.add_argument("--outdir", default="results/ntv2_peft_patient_clf")
    parser.add_argument("--epochs", type=int, default=10)
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--batch-size", type=int, default=2)
    parser.add_argument("--grad-accum", type=int, default=16,
                        help="Gradient accumulation steps (eff_batch = batch_size * grad_accum)")
    parser.add_argument("--max-length", type=int, default=2048)
    parser.add_argument("--patience", type=int, default=5)
    parser.add_argument("--lora-r", type=int, default=8)
    parser.add_argument("--lora-alpha", type=int, default=16)
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--limit", type=int, default=None, help="Limit rows for dry-run")
    args = parser.parse_args()

    out_dir = pathlib.Path(args.outdir) / f"seed_{args.seed}"
    out_dir.mkdir(parents=True, exist_ok=True)

    split_dir = pathlib.Path(args.labels_dir) / f"seed_{args.seed}"

    # ── Tokenizer ─────────────────────────────────────────────────────────────
    print(f"Loading tokenizer: {MODEL_ID}")
    tokenizer = AutoTokenizer.from_pretrained(MODEL_ID, trust_remote_code=True)

    # ── Datasets ──────────────────────────────────────────────────────────────
    print("Building datasets...")
    print("  Train:")
    train_ds = build_patient_dataset(args.sequences, split_dir / "train.csv",
                                      tokenizer, args.max_length, args.task)
    print("  Val:")
    val_ds = build_patient_dataset(args.sequences, split_dir / "val.csv",
                                    tokenizer, args.max_length, args.task)
    print("  Test:")
    test_ds = build_patient_dataset(args.sequences, split_dir / "test.csv",
                                     tokenizer, args.max_length, args.task)

    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True,
                              num_workers=0, pin_memory=True)
    val_loader = DataLoader(val_ds, batch_size=args.batch_size, num_workers=0)
    test_loader = DataLoader(test_ds, batch_size=args.batch_size, num_workers=0)

    # ── Model ─────────────────────────────────────────────────────────────────
    print(f"\nBuilding model with LoRA (r={args.lora_r}, alpha={args.lora_alpha})")
    model = NTv2ForSequenceTask(
        num_labels=2, is_regression=False,
        lora_r=args.lora_r, lora_alpha=args.lora_alpha,
    ).to(args.device)

    # ── Class weights ─────────────────────────────────────────────────────────
    all_labels = [train_ds.labels[i] for i in range(len(train_ds))]
    class_counts = np.bincount(all_labels)
    class_weights = torch.tensor(
        [len(all_labels) / (2 * c) for c in class_counts],
        dtype=torch.float32
    ).to(args.device)
    print(f"  Class weights: {class_weights.tolist()}")
    criterion = torch.nn.CrossEntropyLoss(weight=class_weights)

    optimizer = torch.optim.AdamW(
        [p for p in model.parameters() if p.requires_grad],
        lr=args.lr, weight_decay=0.01,
    )
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)

    # ── Train ─────────────────────────────────────────────────────────────────
    best_val_loss = float("inf")
    patience_counter = 0

    for epoch in range(1, args.epochs + 1):
        train_loss = train_epoch(model, train_loader, optimizer, criterion,
                                  args.device, args.grad_accum)
        val_loss, val_preds, val_labels = eval_epoch(model, val_loader, criterion, args.device)
        scheduler.step()

        val_acc = accuracy_score(val_labels, val_preds)
        print(f"  Epoch {epoch:2d}  train_loss={train_loss:.4f}  "
              f"val_loss={val_loss:.4f}  val_acc={val_acc:.3f}")

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
            save_lora_checkpoint(model, out_dir, "best")
            torch.save(model.head.state_dict(), out_dir / "best_head.pt")
        else:
            patience_counter += 1
            if patience_counter >= args.patience:
                print(f"  Early stopping at epoch {epoch}")
                break

    # ── Evaluate ──────────────────────────────────────────────────────────────
    print("\nEvaluating on test set...")
    model.backbone.load_adapter(out_dir / "best_lora", adapter_name="default")
    model.head.load_state_dict(torch.load(out_dir / "best_head.pt", weights_only=True))

    test_loss, test_preds, test_labels = eval_epoch(model, test_loader, criterion, args.device)

    # Get probabilities for AUROC
    model.eval()
    all_probs = []
    with torch.no_grad():
        for batch in test_loader:
            logits = model(batch["input_ids"].to(args.device),
                          batch["attention_mask"].to(args.device))
            probs = torch.softmax(logits, dim=1)[:, 1]
            all_probs.extend(probs.cpu().numpy())

    metrics = {
        "test_loss": float(test_loss),
        "accuracy": float(accuracy_score(test_labels, test_preds)),
        "balanced_accuracy": float(balanced_accuracy_score(test_labels, test_preds)),
        "f1_weighted": float(f1_score(test_labels, test_preds, average="weighted")),
        "auroc": float(roc_auc_score(test_labels, np.array(all_probs))),
    }
    print(f"\n  Test metrics: {json.dumps(metrics, indent=2)}")
    print(classification_report(test_labels, test_preds))
    save_metrics(metrics, out_dir)
    print(f"  Done → {out_dir}")


if __name__ == "__main__":
    main()
