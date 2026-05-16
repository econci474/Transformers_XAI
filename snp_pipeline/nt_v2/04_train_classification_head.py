"""
04_train_classification_head.py
================================
Train a simple MLP classification head on pre-extracted NT v2 embeddings.

Usage:
  python nt_v2/04_train_classification_head.py \
      --embeddings D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/patient_embeddings/ntv2/embeddings.npz \
      --labels-dir "D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified_ever_convert/tabular/baseline" \
      --seed 0 --task Label_ever_convert --outdir results/ntv2_clf_head
"""

import argparse
import json
import pathlib

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                              f1_score, roc_auc_score, classification_report)
from torch.utils.data import DataLoader, TensorDataset


class ClassificationHead(nn.Module):
    """MLP: Linear → ReLU → Dropout → Linear."""

    def __init__(self, input_dim: int, hidden_dim: int = 256,
                 num_classes: int = 2, dropout: float = 0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, num_classes),
        )

    def forward(self, x):
        return self.net(x)


def load_split(npz_path, csv_path, task_col):
    """Load embeddings for patients in the given clinical CSV split."""
    data = np.load(npz_path, allow_pickle=True)
    df = pd.read_csv(csv_path)
    ptids = df["Patient_ID"].tolist()
    labels = df[task_col].values

    embs = []
    valid_labels = []
    for ptid, label in zip(ptids, labels):
        if ptid in data:
            arr = data[ptid]  # (n_windows, hidden_dim)
            # Mean-pool across windows (ignoring NaN windows)
            mask = ~np.isnan(arr).any(axis=1)
            if mask.sum() > 0:
                embs.append(np.nanmean(arr[mask], axis=0))
                valid_labels.append(label)

    return np.stack(embs), np.array(valid_labels, dtype=np.int64)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--embeddings", required=True, help="Path to embeddings.npz")
    parser.add_argument("--labels-dir", required=True,
                        help="Dir containing seed_0/, seed_1/, seed_2/ with train/val/test.csv")
    parser.add_argument("--seed", type=int, default=0, help="Which seed split to use")
    parser.add_argument("--task", default="Label_ever_convert", help="Label column name")
    parser.add_argument("--outdir", default="results/ntv2_clf_head")
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--patience", type=int, default=15)
    parser.add_argument("--hidden-dim", type=int, default=256)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = parser.parse_args()

    out_dir = pathlib.Path(args.outdir) / f"seed_{args.seed}"
    out_dir.mkdir(parents=True, exist_ok=True)

    split_dir = pathlib.Path(args.labels_dir) / f"seed_{args.seed}"

    # ── Load data ─────────────────────────────────────────────────────────────
    print(f"Loading embeddings + labels (seed={args.seed}, task={args.task})")
    X_train, y_train = load_split(args.embeddings, split_dir / "train.csv", args.task)
    X_val,   y_val   = load_split(args.embeddings, split_dir / "val.csv",   args.task)
    X_test,  y_test  = load_split(args.embeddings, split_dir / "test.csv",  args.task)
    print(f"  Train: {len(y_train)} (pos={y_train.sum()})")
    print(f"  Val:   {len(y_val)}   (pos={y_val.sum()})")
    print(f"  Test:  {len(y_test)}  (pos={y_test.sum()})")

    input_dim = X_train.shape[1]
    num_classes = len(np.unique(y_train))

    # ── Class weights ─────────────────────────────────────────────────────────
    class_counts = np.bincount(y_train)
    class_weights = torch.tensor(
        [len(y_train) / (num_classes * c) for c in class_counts],
        dtype=torch.float32
    ).to(args.device)
    print(f"  Class weights: {class_weights.tolist()}")

    # ── Dataloaders ───────────────────────────────────────────────────────────
    def make_loader(X, y, shuffle=False):
        ds = TensorDataset(torch.tensor(X, dtype=torch.float32),
                           torch.tensor(y, dtype=torch.long))
        return DataLoader(ds, batch_size=args.batch_size, shuffle=shuffle)

    train_loader = make_loader(X_train, y_train, shuffle=True)
    val_loader = make_loader(X_val, y_val)
    test_loader = make_loader(X_test, y_test)

    # ── Model ─────────────────────────────────────────────────────────────────
    model = ClassificationHead(input_dim, args.hidden_dim, num_classes).to(args.device)
    criterion = nn.CrossEntropyLoss(weight=class_weights)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=0.01)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)

    # ── Training ──────────────────────────────────────────────────────────────
    best_val_loss = float("inf")
    patience_counter = 0

    for epoch in range(1, args.epochs + 1):
        model.train()
        train_loss = 0
        for X_batch, y_batch in train_loader:
            X_batch, y_batch = X_batch.to(args.device), y_batch.to(args.device)
            optimizer.zero_grad()
            logits = model(X_batch)
            loss = criterion(logits, y_batch)
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * len(y_batch)
        train_loss /= len(y_train)
        scheduler.step()

        # Validate
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for X_batch, y_batch in val_loader:
                X_batch, y_batch = X_batch.to(args.device), y_batch.to(args.device)
                logits = model(X_batch)
                val_loss += criterion(logits, y_batch).item() * len(y_batch)
        val_loss /= len(y_val)

        if epoch % 10 == 0 or epoch == 1:
            print(f"  Epoch {epoch:3d}  train_loss={train_loss:.4f}  val_loss={val_loss:.4f}")

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
            torch.save(model.state_dict(), out_dir / "best_model.pt")
        else:
            patience_counter += 1
            if patience_counter >= args.patience:
                print(f"  Early stopping at epoch {epoch}")
                break

    # ── Evaluate ──────────────────────────────────────────────────────────────
    model.load_state_dict(torch.load(out_dir / "best_model.pt", weights_only=True))
    model.eval()

    all_preds, all_probs, all_labels = [], [], []
    with torch.no_grad():
        for X_batch, y_batch in test_loader:
            X_batch = X_batch.to(args.device)
            logits = model(X_batch)
            probs = torch.softmax(logits, dim=1)
            preds = logits.argmax(dim=1)
            all_preds.extend(preds.cpu().numpy())
            all_probs.extend(probs[:, 1].cpu().numpy())
            all_labels.extend(y_batch.numpy())

    all_preds = np.array(all_preds)
    all_probs = np.array(all_probs)
    all_labels = np.array(all_labels)

    metrics = {
        "accuracy": float(accuracy_score(all_labels, all_preds)),
        "balanced_accuracy": float(balanced_accuracy_score(all_labels, all_preds)),
        "f1": float(f1_score(all_labels, all_preds, average="weighted")),
        "auroc": float(roc_auc_score(all_labels, all_probs)) if num_classes == 2 else None,
    }
    print(f"\n  Test metrics: {json.dumps(metrics, indent=2)}")
    print(classification_report(all_labels, all_preds))

    with open(out_dir / "metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"  Saved to {out_dir}")


if __name__ == "__main__":
    main()
