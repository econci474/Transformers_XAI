"""
15_train_classification_head.py
================================
Train a small aggregation + MLP head on per-patient BMFM embeddings for the
T1_binary task (CN vs MCI/AD) using the no_cdr_stratified splits.

Inputs
------
  patient_embeddings/<upstream_tag>/embeddings.npz  (from 14_extract_patient_embeddings.py)
  patient_window_sequences/windows.tsv               (window_id → chrom mapping)
  derivatives/clinical/no_cdr_stratified/tabular/baseline/seed_{0,1,2}/{train,val,test}.csv

Aggregation strategies (all output a single 768-dim or 22*768 vector per patient):
  - mean_all                : mean(all windows)
  - chrom_mean_then_mean    : mean of per-chrom means
  - attn_all                : softmax-weighted sum over windows
  - chrom_mean_then_attn    : per-chrom mean, then attention across chrom-means
  - chrom_attn_then_attn    : per-chrom attention, then attention across chrom pools
  - concat_chrom_means      : concat(chrom-means), padding missing chroms to zero

MLP head: hidden -> 256 -> 1, BCE-with-logits, AdamW, early stopping on val AUC.

Output
------
  classification_results/<upstream_tag>/<aggregation>/seed_<N>/
    metrics.json     {auc, balanced_accuracy, f1, sens, spec, ...}
    head.pt          state dict
    val_curve.csv    epoch, val_auc, val_loss

Usage
-----
  python 15_train_classification_head.py \\
      --embeddings patient_embeddings/base_ref/embeddings.npz \\
      --windows patient_window_sequences/windows.tsv \\
      --splits-root D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified/tabular/baseline \\
      --upstream-tag base_ref \\
      --aggregation mean_all \\
      --seed 0 \\
      --output-root classification_results/
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import (
    balanced_accuracy_score,
    f1_score,
    roc_auc_score,
    recall_score,
    confusion_matrix,
)
from torch.utils.data import DataLoader, TensorDataset

AGGREGATIONS = (
    "mean_all",
    "chrom_mean_then_mean",
    "attn_all",
    "chrom_mean_then_attn",
    "chrom_attn_then_attn",
    "concat_chrom_means",
)


# ── Aggregation modules ──────────────────────────────────────────────────────
class AdditiveAttention(nn.Module):
    """Single-headed additive attention, mean-pooling fallback when input has one row."""

    def __init__(self, dim: int, attn_hidden: int = 128):
        super().__init__()
        self.v = nn.Linear(dim, attn_hidden, bias=False)
        self.u = nn.Linear(attn_hidden, 1, bias=False)

    def forward(self, x: torch.Tensor, mask: torch.Tensor | None = None) -> torch.Tensor:
        # x: (B, N, D); mask: (B, N) bool (True=valid)
        scores = self.u(torch.tanh(self.v(x))).squeeze(-1)  # (B, N)
        if mask is not None:
            scores = scores.masked_fill(~mask, float("-inf"))
        weights = torch.softmax(scores, dim=-1)             # (B, N)
        return (weights.unsqueeze(-1) * x).sum(dim=1)        # (B, D)


class Aggregator(nn.Module):
    def __init__(self, mode: str, dim: int, n_chroms: int):
        super().__init__()
        if mode not in AGGREGATIONS:
            raise ValueError(f"unknown aggregation: {mode}")
        self.mode = mode
        self.dim = dim
        self.n_chroms = n_chroms
        if mode == "attn_all":
            self.attn = AdditiveAttention(dim)
        elif mode == "chrom_mean_then_attn":
            self.chrom_attn = AdditiveAttention(dim)
        elif mode == "chrom_attn_then_attn":
            self.intra = AdditiveAttention(dim)
            self.inter = AdditiveAttention(dim)
        # mean_all, chrom_mean_then_mean, concat_chrom_means : no params

    @property
    def output_dim(self) -> int:
        return self.n_chroms * self.dim if self.mode == "concat_chrom_means" else self.dim

    def forward(self, x: torch.Tensor, mask: torch.Tensor, chrom_idx: torch.Tensor) -> torch.Tensor:
        """
        x          : (B, N, D)  windows, padded
        mask       : (B, N) bool — True where window is valid for the patient
        chrom_idx  : (B, N) int  — chrom-of-window index in [0, n_chroms); -1 where padding
        """
        if self.mode == "mean_all":
            x_masked = x * mask.unsqueeze(-1)
            return x_masked.sum(1) / mask.sum(1, keepdim=True).clamp_min(1)

        if self.mode == "attn_all":
            return self.attn(x, mask)

        # All chrom-aware modes group windows by chrom_idx
        B, N, D = x.shape
        # one-hot: (B, N, C)
        one_hot = F.one_hot(chrom_idx.clamp_min(0), num_classes=self.n_chroms).float()
        one_hot = one_hot * mask.unsqueeze(-1)
        counts = one_hot.sum(1).clamp_min(1e-6)              # (B, C)
        chrom_present = (one_hot.sum(1) > 0)                  # (B, C)

        if self.mode == "concat_chrom_means":
            # Sum then divide → mean per chrom; zero where absent.
            chrom_means = torch.einsum("bnc,bnd->bcd", one_hot, x) / counts.unsqueeze(-1)
            return chrom_means.reshape(B, self.n_chroms * D)

        if self.mode == "chrom_mean_then_mean":
            chrom_means = torch.einsum("bnc,bnd->bcd", one_hot, x) / counts.unsqueeze(-1)
            cnts = chrom_present.float().sum(1, keepdim=True).clamp_min(1)
            return (chrom_means * chrom_present.unsqueeze(-1)).sum(1) / cnts

        if self.mode == "chrom_mean_then_attn":
            chrom_means = torch.einsum("bnc,bnd->bcd", one_hot, x) / counts.unsqueeze(-1)
            return self.chrom_attn(chrom_means, chrom_present)

        if self.mode == "chrom_attn_then_attn":
            # For each chrom c, run attention over the windows where chrom_idx==c.
            # We do this with a vectorised softmax: score every window, then mask by chrom.
            B, N, D = x.shape
            # window-level attention scores (B, N) computed once per chrom via the same head
            scores = self.intra.u(torch.tanh(self.intra.v(x))).squeeze(-1)  # (B, N)

            # For each (b, c), select windows with chrom_idx==c and softmax over them.
            chrom_pools = torch.zeros(B, self.n_chroms, D, device=x.device, dtype=x.dtype)
            for c in range(self.n_chroms):
                m = (chrom_idx == c) & mask  # (B, N)
                masked_scores = scores.masked_fill(~m, float("-inf"))
                # Skip chrom-rows with no windows (avoid NaN from -inf softmax).
                has_any = m.any(dim=1)
                if not has_any.any():
                    continue
                w = torch.softmax(masked_scores, dim=-1)                       # (B, N)
                pool = (w.unsqueeze(-1) * x).sum(1)                            # (B, D)
                chrom_pools[has_any, c] = pool[has_any]
            return self.inter(chrom_pools, chrom_present)

        raise RuntimeError(f"unhandled mode: {self.mode}")


# ── Classifier ───────────────────────────────────────────────────────────────
class MLPHead(nn.Module):
    def __init__(self, in_dim: int, hidden: int = 256, dropout: float = 0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x).squeeze(-1)


class FullModel(nn.Module):
    def __init__(self, aggregator: Aggregator, head: MLPHead):
        super().__init__()
        self.aggregator = aggregator
        self.head = head

    def forward(self, x: torch.Tensor, mask: torch.Tensor, chrom_idx: torch.Tensor) -> torch.Tensor:
        pooled = self.aggregator(x, mask, chrom_idx)
        return self.head(pooled)


# ── Data assembly ────────────────────────────────────────────────────────────
def build_chrom_map(windows_tsv: Path) -> tuple[dict[str, int], list[str]]:
    """Return (window_id → chrom_idx, ordered chrom list)."""
    df = pd.read_csv(windows_tsv, sep="\t")
    df["chrom"] = df["chrom"].astype(str)
    chroms = sorted(df["chrom"].unique(), key=lambda c: int(c) if c.isdigit() else 99)
    chrom_idx_map = {c: i for i, c in enumerate(chroms)}
    return {row["window_id"]: chrom_idx_map[row["chrom"]] for _, row in df.iterrows()}, chroms


def load_embeddings(npz_path: Path):
    npz = np.load(npz_path, allow_pickle=True)
    window_ids = list(npz["window_ids"])
    embeddings = {k: npz[k] for k in npz.files if k != "window_ids"}
    return embeddings, window_ids


def load_labels(splits_root: Path, seed: int) -> dict[str, pd.DataFrame]:
    out = {}
    for split in ("train", "val", "test"):
        df = pd.read_csv(splits_root / f"seed_{seed}" / f"{split}.csv")
        df["Patient_ID"] = df["Patient_ID"].astype(str)
        df["y"] = (df["Label_bl_multi"] != "CN").astype(int)
        out[split] = df[["Patient_ID", "y"]]
    return out


def stack_patient_arrays(
    embeddings: dict[str, np.ndarray],
    window_ids: list[str],
    chrom_of_window: dict[str, int],
    patient_ids: list[str],
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Return (x, mask, chrom_idx) tensors (B, N, D), (B, N), (B, N)."""
    N = len(window_ids)
    D = next(iter(embeddings.values())).shape[1]
    B = len(patient_ids)
    x = np.zeros((B, N, D), dtype=np.float32)
    mask = np.zeros((B, N), dtype=bool)
    chrom_idx = np.full((B, N), -1, dtype=np.int64)

    for j, w in enumerate(window_ids):
        chrom_idx[:, j] = chrom_of_window[w]

    for i, ptid in enumerate(patient_ids):
        arr = embeddings.get(ptid)
        if arr is None:
            continue
        x[i] = arr
        mask[i] = ~np.isnan(arr).any(axis=1)
    return torch.from_numpy(x), torch.from_numpy(mask), torch.from_numpy(chrom_idx)


# ── Training loop ────────────────────────────────────────────────────────────
def train_one(
    model: FullModel, train_data, val_data,
    epochs: int, batch_size: int, lr: float, weight_decay: float, device: torch.device,
    patience: int,
) -> tuple[FullModel, dict, list]:
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    loss_fn = nn.BCEWithLogitsLoss()

    x_tr, m_tr, c_tr, y_tr = (t.to(device) for t in train_data)
    x_va, m_va, c_va, y_va = (t.to(device) for t in val_data)

    loader = DataLoader(TensorDataset(x_tr, m_tr, c_tr, y_tr),
                        batch_size=batch_size, shuffle=True)

    best_state, best_auc, since = None, -1.0, 0
    curve = []
    for ep in range(epochs):
        model.train()
        for xb, mb, cb, yb in loader:
            opt.zero_grad()
            logits = model(xb, mb, cb)
            loss = loss_fn(logits, yb.float())
            loss.backward()
            opt.step()
        # eval
        model.eval()
        with torch.no_grad():
            val_logits = model(x_va, m_va, c_va).cpu().numpy()
        val_auc = roc_auc_score(y_va.cpu().numpy(), val_logits) if y_va.unique().numel() > 1 else float("nan")
        val_loss = float(loss_fn(torch.from_numpy(val_logits), y_va.cpu().float()).item())
        curve.append({"epoch": ep, "val_auc": val_auc, "val_loss": val_loss})
        if val_auc > best_auc:
            best_auc = val_auc
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            since = 0
        else:
            since += 1
        if since >= patience:
            break

    if best_state is not None:
        model.load_state_dict(best_state)
    return model, {"best_val_auc": best_auc, "n_epochs": len(curve)}, curve


def evaluate(model: FullModel, x, m, c, y, device: torch.device) -> dict:
    model.eval()
    with torch.no_grad():
        logits = model(x.to(device), m.to(device), c.to(device)).cpu().numpy()
    y_np = y.cpu().numpy() if hasattr(y, "cpu") else y
    pred = (logits > 0).astype(int)
    out = {
        "auc": float(roc_auc_score(y_np, logits)) if len(np.unique(y_np)) > 1 else float("nan"),
        "balanced_accuracy": float(balanced_accuracy_score(y_np, pred)),
        "f1": float(f1_score(y_np, pred, zero_division=0)),
        "sensitivity": float(recall_score(y_np, pred, pos_label=1, zero_division=0)),
        "specificity": float(recall_score(y_np, pred, pos_label=0, zero_division=0)),
    }
    cm = confusion_matrix(y_np, pred, labels=[0, 1])
    out["confusion_matrix"] = {"tn": int(cm[0, 0]), "fp": int(cm[0, 1]),
                                "fn": int(cm[1, 0]), "tp": int(cm[1, 1])}
    return out


# ── Main ─────────────────────────────────────────────────────────────────────
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--embeddings", type=Path, required=True)
    ap.add_argument("--windows", type=Path, required=True,
                    help="patient_window_sequences/windows.tsv")
    ap.add_argument("--splits-root", type=Path, required=True,
                    help="…/no_cdr_stratified/tabular/baseline")
    ap.add_argument("--upstream-tag", type=str, required=True)
    ap.add_argument("--aggregation", choices=AGGREGATIONS, required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--output-root", type=Path, required=True)
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--weight-decay", type=float, default=1e-4)
    ap.add_argument("--patience", type=int, default=20)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    device = torch.device(args.device)
    chrom_of_window, chroms = build_chrom_map(args.windows)
    n_chroms = len(chroms)
    print(f"Chromosomes present: {chroms} ({n_chroms} groups)")

    embeddings, window_ids = load_embeddings(args.embeddings)
    print(f"Embeddings: {len(embeddings)} patients, {len(window_ids)} windows")

    labels = load_labels(args.splits_root, args.seed)
    for split, df in labels.items():
        n_with_emb = df["Patient_ID"].isin(embeddings).sum()
        print(f"  {split}: {len(df)} patients, {n_with_emb} with embeddings  "
              f"(y mean={df['y'].mean():.3f})")

    # Filter to patients with embeddings.
    splits = {}
    for split, df in labels.items():
        df = df[df["Patient_ID"].isin(embeddings)].reset_index(drop=True)
        x, m, c = stack_patient_arrays(embeddings, window_ids, chrom_of_window,
                                       df["Patient_ID"].tolist())
        y = torch.from_numpy(df["y"].values.astype(np.int64))
        splits[split] = (x, m, c, y)

    # Build model.
    hidden_dim = next(iter(embeddings.values())).shape[1]
    aggregator = Aggregator(args.aggregation, hidden_dim, n_chroms)
    head = MLPHead(aggregator.output_dim, hidden=256, dropout=0.3)
    model = FullModel(aggregator, head).to(device)
    print(f"\nModel:  aggregator={args.aggregation}  in_dim={aggregator.output_dim}")

    model, fit_info, curve = train_one(
        model, splits["train"], splits["val"],
        epochs=args.epochs, batch_size=args.batch_size,
        lr=args.lr, weight_decay=args.weight_decay,
        device=device, patience=args.patience,
    )
    print(f"\nFit:  {fit_info}")

    test_metrics = evaluate(model, *splits["test"], device=device)
    val_metrics = evaluate(model, *splits["val"], device=device)
    print(f"\nVAL  : {val_metrics}")
    print(f"TEST : {test_metrics}")

    out_dir = args.output_root / args.upstream_tag / args.aggregation / f"seed_{args.seed}"
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "metrics.json", "w") as f:
        json.dump({
            "upstream_tag": args.upstream_tag,
            "aggregation": args.aggregation,
            "seed": args.seed,
            "fit_info": fit_info,
            "val": val_metrics,
            "test": test_metrics,
        }, f, indent=2)
    pd.DataFrame(curve).to_csv(out_dir / "val_curve.csv", index=False)
    torch.save(model.state_dict(), out_dir / "head.pt")
    print(f"\nWrote {out_dir / 'metrics.json'}")


if __name__ == "__main__":
    main()
