"""
18_patient_classification_finetune.py
======================================
End-to-end per-patient classifier: raw BMFM-DNA-REF + LoRA + aggregator + MLP.

Tests whether the GWAS-regression intermediate task is necessary by
fine-tuning BMFM **directly** on the patient classification task (no
intermediate pretext). Compare against script 15's frozen-cached-embedding
sweep that uses GWAS-regression-finetuned upstreams.

Pipeline per training step (one patient = one batch sample):

    one patient ─► 183 windows ─► tokenize ─► BMFM (raw + LoRA adapters)
                                                          │
                                              183 × 768 per-window CLS
                                                          │
                                          Aggregator (one of 6 strategies)
                                                          │
                                              1 × 768 patient vector
                                                          │
                                                       MLP head
                                                          │
                                                  focal-loss(y_patient)

Output directory mirrors script 15 with an extra ``mode`` segment:

    classification_results/<label_mode>/<upstream>/<mode>/<aggregation>/seed_<N>/metrics.json

where ``upstream`` is fixed to ``raw_bmfm_ref`` for this script and ``mode``
is ``lora`` (or ``frozen_live`` if --no-lora is passed).

Aggregator / MLPHead / metrics implementations are imported from script 15.

Usage on Colab A100
-------------------
After mounting Drive + installing BMFM + cloning the repo:

    python snp_pipeline/18_patient_classification_finetune.py \\
        --sequences /content/drive/MyDrive/ADNI_SNP/inputs/patient_sequences/all_patients.csv \\
        --windows   /content/drive/MyDrive/ADNI_SNP/inputs/windows/windows.tsv \\
        --splits-root /content/drive/MyDrive/ADNI_SNP/splits/no_cdr_stratified_ever_convert/baseline \\
        --label-mode ever_convert \\
        --aggregation chrom_mean_then_attn \\
        --seed 0 \\
        --output-root /content/drive/MyDrive/ADNI_SNP/classification_results \\
        --lora-r 16 --epochs 5 --batch-patients 1 --grad-checkpoint

The sweep driver (a wrapper bash / notebook cell) iterates this over
{cls, mean}-equivalent pools is **not** applicable here because we don't
pre-pool inside BMFM — the aggregator picks one pool implicitly via its
``cls`` vs ``mean`` choice. To match script 15's per-window pool dimension,
this script exposes ``--per-window-pool {cls, mean}`` controlling how each
window's hidden states get reduced to its 768-dim vector before aggregation.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader

# Reuse aggregator + MLPHead + helpers from script 15
# (script 15's filename starts with a digit, so we can't `import` it directly —
# load it via spec_from_file_location.)
import importlib.util as _il
_SCRIPT15 = Path(__file__).parent / "15_train_classification_head.py"
_spec = _il.spec_from_file_location("script15", _SCRIPT15)
_script15 = _il.module_from_spec(_spec)
_spec.loader.exec_module(_script15)
Aggregator      = _script15.Aggregator
MLPHead         = _script15.MLPHead
AGGREGATIONS    = _script15.AGGREGATIONS
build_chrom_map = _script15.build_chrom_map
load_labels     = _script15.load_labels
evaluate        = _script15.evaluate

# ── DNA helpers ───────────────────────────────────────────────────────────────
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")
def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


# ── Dataset: one item = one patient's 183 windows + label ────────────────────
class PatientWindowsDataset(Dataset):
    """Yields per-patient stacks of tokenised window sequences + label.

    Each __getitem__ returns:
        input_ids:      (n_windows, max_len) int64
        attention_mask: (n_windows, max_len) int64
        chrom_idx:      (n_windows,)         int64
        y:              scalar               int64
        ptid:           str

    Random fwd/rc sampling per window applied only when ``training=True``.
    """
    def __init__(self, *, sequences_df: pd.DataFrame, patient_ids: list[str],
                 labels: dict[str, int], window_order: list[str],
                 chrom_of_window: dict[str, int],
                 tokenizer, max_len: int, training: bool,
                 rc_prob: float = 0.5):
        self.patient_ids = list(patient_ids)
        self.labels = labels
        self.window_order = window_order
        self.chrom_idx_vec = torch.tensor(
            [chrom_of_window[w] for w in window_order], dtype=torch.long)
        self.tokenizer = tokenizer
        self.max_len = max_len
        self.training = training
        self.rc_prob = rc_prob

        # Group sequences by (PTID, window_id) for fast lookup
        df = sequences_df[sequences_df["PTID"].astype(str).isin(set(patient_ids))]
        self.seq_map: dict[tuple[str, str], str] = {}
        for row in df.itertuples(index=False):
            self.seq_map[(str(row.PTID), row.window_id)] = row.sequence

    def __len__(self) -> int:
        return len(self.patient_ids)

    def __getitem__(self, idx: int):
        ptid = self.patient_ids[idx]
        y = int(self.labels[ptid])
        seqs = []
        for w in self.window_order:
            s = self.seq_map.get((ptid, w))
            if s is None:
                # Pad missing window with N's (should not happen for our data)
                s = "N" * self.max_len
            if self.training and random.random() < self.rc_prob:
                s = reverse_complement(s)
            seqs.append(s)
        encs = self.tokenizer.encode_batch(seqs, add_special_tokens=True)
        ids = torch.tensor([e.ids for e in encs], dtype=torch.long)
        mask = torch.tensor([e.attention_mask for e in encs], dtype=torch.long)
        return {
            "input_ids": ids,
            "attention_mask": mask,
            "chrom_idx": self.chrom_idx_vec,
            "y": torch.tensor(y, dtype=torch.long),
            "ptid": ptid,
        }


def collate_patients(batch: list[dict]) -> dict:
    """Stack a list of patient items into a batch tensor.

    All patients are assumed to have the same n_windows (= 183). The output
    tensors have shape (B, n_windows, ...).
    """
    out = {
        "input_ids":      torch.stack([b["input_ids"]      for b in batch], dim=0),
        "attention_mask": torch.stack([b["attention_mask"] for b in batch], dim=0),
        "chrom_idx":      torch.stack([b["chrom_idx"]      for b in batch], dim=0),
        "y":              torch.stack([b["y"]              for b in batch], dim=0),
        "ptid":           [b["ptid"] for b in batch],
    }
    # All patient-windows are valid (no padding at the window level)
    B, N = out["input_ids"].shape[:2]
    out["window_mask"] = torch.ones(B, N, dtype=torch.bool)
    return out


# ── Model: BMFM(+LoRA) → per-window pool → Aggregator → MLP ─────────────────
class PatientClassifier(nn.Module):
    """End-to-end patient classifier.

    Forward: takes batched per-patient (input_ids, attention_mask, chrom_idx)
    and returns one logit per patient.
    """
    def __init__(self, bmfm: nn.Module, aggregator: Aggregator, head: MLPHead,
                 per_window_pool: str = "cls", micro_batch: int = 16):
        super().__init__()
        self.bmfm = bmfm
        self.aggregator = aggregator
        self.head = head
        if per_window_pool not in ("cls", "mean"):
            raise ValueError(f"per_window_pool must be cls|mean, got {per_window_pool!r}")
        self.per_window_pool = per_window_pool
        # Micro-batch size for chunking windows within a patient's forward.
        # 183 windows × 2048 tokens × 768 H × 22 L is too much to fit a single
        # forward+backward on A100 40 GB even with gradient checkpointing;
        # chunk into micro_batch-sized groups so only ~micro_batch/183 of the
        # per-window activations are live at any moment.
        self.micro_batch = micro_batch
        # Special tokens for BMFM-DNA-REF tokenizer
        self.CLS_ID, self.SEP_ID, self.PAD_ID = 3, 1, 2

    def _bmfm_forward(self, input_ids: torch.Tensor, attention_mask: torch.Tensor) -> torch.Tensor:
        """Forward through BMFM and return per-row 768-dim embedding.

        input_ids:      (R, max_len)    where R = B * n_windows
        attention_mask: (R, max_len)
        returns:        (R, hidden_dim)
        """
        # BMFM expects (R, num_fields, max_len) — num_fields=1 for dna_chunks
        ids3d = input_ids.unsqueeze(1)
        out = self.bmfm(input_ids=ids3d, attention_mask=attention_mask)
        hidden = out.last_hidden_state if hasattr(out, "last_hidden_state") else out[0]
        # hidden: (R, max_len, H)
        if self.per_window_pool == "cls":
            return hidden[:, 0, :]
        # mean over real tokens (exclude CLS/SEP/PAD)
        real = ((input_ids != self.CLS_ID) &
                (input_ids != self.SEP_ID) &
                (input_ids != self.PAD_ID)).to(hidden.dtype)
        denom = real.sum(dim=1, keepdim=True).clamp_min(1.0)
        return (hidden * real.unsqueeze(-1)).sum(dim=1) / denom

    def forward(self, input_ids: torch.Tensor, attention_mask: torch.Tensor,
                chrom_idx: torch.Tensor, window_mask: torch.Tensor) -> torch.Tensor:
        # input_ids: (B, N, S); attention_mask: (B, N, S); chrom_idx: (B, N); window_mask: (B, N)
        B, N, S = input_ids.shape
        ids_flat  = input_ids.reshape(B * N, S)
        mask_flat = attention_mask.reshape(B * N, S)
        # Chunk windows in micro-batches so the live activation set is
        # ~ micro_batch worth of windows, not the whole (B*N) batch. Gradients
        # still flow through every micro_batch via the cat below.
        chunks = []
        total = B * N
        mb = max(1, self.micro_batch)
        for start in range(0, total, mb):
            end = min(start + mb, total)
            chunks.append(self._bmfm_forward(ids_flat[start:end], mask_flat[start:end]))
        per_window = torch.cat(chunks, dim=0)                        # (B*N, H)
        per_window = per_window.view(B, N, -1)                       # (B, N, H)
        pooled = self.aggregator(per_window, window_mask, chrom_idx) # (B, output_dim)
        return self.head(pooled)                                     # (B,)


# ── Focal loss with class weighting ──────────────────────────────────────────
class FocalLoss(nn.Module):
    def __init__(self, gamma: float = 2.0, pos_weight: float = 1.0):
        super().__init__()
        self.gamma = gamma
        self.pos_weight = pos_weight

    def forward(self, logits: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
        # logits: (B,) , targets: (B,) in {0, 1}
        bce = F.binary_cross_entropy_with_logits(
            logits, targets.float(),
            pos_weight=torch.tensor(self.pos_weight, device=logits.device),
            reduction="none",
        )
        p = torch.sigmoid(logits)
        # p_t = p if y=1, 1-p if y=0
        p_t = p * targets + (1 - p) * (1 - targets)
        focal = (1 - p_t) ** self.gamma * bce
        return focal.mean()


# ── BMFM model loading + LoRA injection ─────────────────────────────────────
def load_bmfm_ref_with_lora(*, lora_r: int = 16, lora_alpha: int = 32,
                             lora_dropout: float = 0.1,
                             enable_lora: bool = True,
                             grad_checkpoint: bool = True):
    """Load raw BMFM-DNA-REF base encoder and inject LoRA adapters.

    Returns the wrapped model (PeftModel if enable_lora=True else the base
    SCModernBert encoder).
    """
    # Force SDPA (avoid FA2 multi-field bugs)
    import bmfm_targets.models.predictive.scmodernbert.modeling_scmodernbert as mb
    _Orig = mb.SCModernBertForMultiTaskModeling
    if not getattr(_Orig, "_sdpa_patched", False):
        _orig_init = _Orig.__init__
        def _sdpa_init(self, config, *a, **kw):
            config._attn_implementation = "sdpa"
            return _orig_init(self, config, *a, **kw)
        _Orig.__init__ = _sdpa_init
        _Orig._sdpa_patched = True

    from bmfm_targets.models.model_utils import (
        register_configs_and_models,
        download_ckpt_from_huggingface,
        migrate_checkpoint_if_needed,
        get_base_model_from_config,
    )
    register_configs_and_models()
    ckpt_path = download_ckpt_from_huggingface(
        "ibm-research/biomed.dna.ref.modernbert.113m.v1")
    state = migrate_checkpoint_if_needed(ckpt_path)
    model = get_base_model_from_config(state["hyper_parameters"]["model_config"])

    # Load encoder weights only (drop scmodernbert. prefix)
    sd_stripped = {k.removeprefix("scmodernbert."): v
                   for k, v in state["state_dict"].items()
                   if k.startswith("scmodernbert.")}
    r = model.load_state_dict(sd_stripped, strict=False)
    if r.missing_keys:
        print(f"[bmfm] warn missing keys: {len(r.missing_keys)} (first 3: {r.missing_keys[:3]})")

    # BMFM's get_input_embeddings() is hardcoded for the scRNA case
    # (returns self.embeddings.genes_embeddings, which doesn't exist on DNA
    # models — they have dna_chunks_embeddings instead). Patch it so peft's
    # prepare_model_for_gradient_checkpointing doesn't blow up.
    def _patched_get_input_embeddings(self):
        emb = self.embeddings
        for attr_name in ("dna_chunks_embeddings", "genes_embeddings"):
            attr = getattr(emb, attr_name, None)
            if attr is not None:
                return attr
        # Last-resort: find any nn.Embedding child on the embeddings layer
        for name, child in emb.named_children():
            if isinstance(child, nn.Embedding):
                return child
        raise AttributeError(
            f"No input embedding layer found on {type(emb).__name__}; "
            "expected dna_chunks_embeddings or genes_embeddings."
        )
    type(model).get_input_embeddings = _patched_get_input_embeddings

    # Gradient checkpointing — required to fit 183 windows × 2048 tokens × LoRA bp
    if grad_checkpoint and hasattr(model, "gradient_checkpointing_enable"):
        model.gradient_checkpointing_enable()
        print("[bmfm] gradient checkpointing: enabled")

    if not enable_lora:
        # Freeze BMFM entirely — used by --no-lora ("frozen_live") mode
        for p in model.parameters():
            p.requires_grad_(False)
        return model

    # LoRA via peft
    from peft import LoraConfig, get_peft_model
    # ModernBERT's attention block uses `Wqkv` and `Wo` (per BMFM source)
    target_modules = ["Wqkv", "Wo", "Wi"]  # Wi is the MLP input projection
    lora_cfg = LoraConfig(
        r=lora_r, lora_alpha=lora_alpha, lora_dropout=lora_dropout,
        bias="none", target_modules=target_modules,
    )
    model = get_peft_model(model, lora_cfg)
    n_train = sum(p.numel() for p in model.parameters() if p.requires_grad)
    n_total = sum(p.numel() for p in model.parameters())
    print(f"[lora] trainable params: {n_train/1e6:.2f} M / total {n_total/1e6:.1f} M  "
          f"({100*n_train/n_total:.2f}%)")
    return model


def load_tokenizer():
    """Load the BMFM-DNA-REF ref2vec tokenizer from the HF cache."""
    from tokenizers import Tokenizer
    import os
    hf_home = os.environ.get("HF_HOME") or os.environ.get("HF_HUB_CACHE")
    if hf_home:
        roots = [Path(hf_home), Path(hf_home)/"hub"]
    else:
        roots = [Path.home()/".cache/huggingface", Path.home()/".cache/huggingface/hub"]
    for root in roots:
        candidates = list(root.glob("**/models--ibm-research--biomed.dna.ref.modernbert.113m.v1/snapshots/*/tokenizers/dna_chunks/tokenizer.json"))
        if candidates:
            return Tokenizer.from_file(str(candidates[0]))
    raise FileNotFoundError("Could not find dna_chunks/tokenizer.json in HF cache")


# ── Patient-level metrics (reuses script 15's evaluate via wrapper) ─────────
def patient_metrics(logits: np.ndarray, y: np.ndarray) -> dict:
    """Compute patient-level metrics on (n_patients,) logits + labels."""
    from sklearn.metrics import (
        balanced_accuracy_score, f1_score, roc_auc_score,
        recall_score, precision_score, confusion_matrix,
    )
    pred = (logits > 0).astype(int)
    out = {
        "auc": float(roc_auc_score(y, logits)) if len(np.unique(y)) > 1 else float("nan"),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "f1": float(f1_score(y, pred, zero_division=0)),
        "precision": float(precision_score(y, pred, pos_label=1, zero_division=0)),
        "recall": float(recall_score(y, pred, pos_label=1, zero_division=0)),
        "sensitivity": float(recall_score(y, pred, pos_label=1, zero_division=0)),
        "specificity": float(recall_score(y, pred, pos_label=0, zero_division=0)),
    }
    cm = confusion_matrix(y, pred, labels=[0, 1])
    out["confusion_matrix"] = {"tn": int(cm[0, 0]), "fp": int(cm[0, 1]),
                                "fn": int(cm[1, 0]), "tp": int(cm[1, 1])}
    return out


# ── Train / eval loop ────────────────────────────────────────────────────────
def run_inference(model: PatientClassifier, loader: DataLoader, device: torch.device
                  ) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return (logits, y, ptids) on the given loader."""
    model.eval()
    all_logits, all_y, all_ptids = [], [], []
    with torch.no_grad():
        for batch in loader:
            ids  = batch["input_ids"].to(device)
            mask = batch["attention_mask"].to(device)
            cidx = batch["chrom_idx"].to(device)
            wmsk = batch["window_mask"].to(device)
            logits = model(ids, mask, cidx, wmsk)
            all_logits.append(logits.float().cpu().numpy())
            all_y.append(batch["y"].cpu().numpy())
            all_ptids.extend(batch["ptid"])
    return np.concatenate(all_logits), np.concatenate(all_y), all_ptids


def train_one(model: PatientClassifier, train_loader: DataLoader,
              val_loader: DataLoader, *, epochs: int, lr: float,
              weight_decay: float, focal_gamma: float, pos_weight: float,
              patience: int, device: torch.device,
              out_dir: Path | None = None):
    """Train end-to-end with per-epoch partial-progress checkpointing.

    If ``out_dir`` is given, after every epoch we write a fresh ``val_curve.csv``
    and ``metrics_partial.json`` capturing the best val AUC so far. Lets us
    survive a Colab disconnect mid-training; the next run sees the partial
    artefacts (though it currently still re-trains from scratch — true
    optimizer-state resume is overkill at 5 epochs).
    """
    opt = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad],
                            lr=lr, weight_decay=weight_decay)
    loss_fn = FocalLoss(gamma=focal_gamma, pos_weight=pos_weight)
    best_auc, since = -1.0, 0
    best_state = None
    curve = []
    for ep in range(epochs):
        t0 = time.time()
        model.train()
        running, n_batches = 0.0, 0
        for batch in train_loader:
            ids  = batch["input_ids"].to(device, non_blocking=True)
            mask = batch["attention_mask"].to(device, non_blocking=True)
            cidx = batch["chrom_idx"].to(device, non_blocking=True)
            wmsk = batch["window_mask"].to(device, non_blocking=True)
            y    = batch["y"].to(device, non_blocking=True)
            logits = model(ids, mask, cidx, wmsk)
            loss = loss_fn(logits, y)
            opt.zero_grad()
            loss.backward()
            opt.step()
            running += float(loss.item())
            n_batches += 1
        # eval
        val_logits, val_y, _ = run_inference(model, val_loader, device)
        val_m = patient_metrics(val_logits, val_y)
        ep_t = time.time() - t0
        curve.append({"epoch": ep, "train_loss": running/max(1,n_batches),
                       "val_auc": val_m["auc"], "val_bACC": val_m["balanced_accuracy"],
                       "elapsed_s": ep_t})
        print(f"  ep {ep:>2}  train_loss={running/max(1,n_batches):.4f}  "
              f"val_auc={val_m['auc']:.3f}  val_bACC={val_m['balanced_accuracy']:.3f}  "
              f"({ep_t:.1f}s)")
        if val_m["auc"] > best_auc:
            best_auc = val_m["auc"]
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            since = 0
        else:
            since += 1

        # Per-epoch partial-progress checkpoint (for disconnect recovery)
        if out_dir is not None:
            try:
                pd.DataFrame(curve).to_csv(out_dir / "val_curve.csv", index=False)
                with open(out_dir / "metrics_partial.json", "w") as f:
                    json.dump({
                        "epoch_completed": ep,
                        "best_val_auc_so_far": best_auc,
                        "last_val": val_m,
                    }, f, indent=2)
            except Exception as e:
                print(f"  [warn] partial-save failed: {e}")

        if since >= patience:
            print(f"  early stop at epoch {ep} (no val_auc improvement in {patience} epochs)")
            break
    if best_state is not None:
        model.load_state_dict(best_state)
    return model, {"best_val_auc": best_auc, "n_epochs": len(curve)}, curve


# ── Main ─────────────────────────────────────────────────────────────────────
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    # Data
    ap.add_argument("--sequences", type=Path, required=True,
                    help="all_patients.csv from 13_build_patient_window_sequences.py")
    ap.add_argument("--windows",   type=Path, required=True)
    ap.add_argument("--splits-root", type=Path, required=True)
    # Experiment selectors
    ap.add_argument("--label-mode", choices=["bl_multi", "ever_convert"], required=True)
    ap.add_argument("--aggregation", choices=AGGREGATIONS, required=True)
    ap.add_argument("--per-window-pool", choices=["cls", "mean"], default="cls")
    ap.add_argument("--seed", type=int, default=0)
    # Modes
    ap.add_argument("--no-lora", action="store_true",
                    help="Disable LoRA, keep BMFM fully frozen (mode=frozen_live). "
                         "By default LoRA is on (mode=lora).")
    ap.add_argument("--lora-r", type=int, default=16)
    ap.add_argument("--lora-alpha", type=int, default=32)
    ap.add_argument("--lora-dropout", type=float, default=0.1)
    # Training
    ap.add_argument("--epochs", type=int, default=5)
    ap.add_argument("--batch-patients", type=int, default=1,
                    help="Patients per gradient step (memory-limited; 1 on A100 40 GB).")
    ap.add_argument("--micro-batch", type=int, default=16,
                    help="Micro-batch size for chunking a patient's 183 windows "
                         "inside BMFM forward. Smaller = less memory but slower.")
    ap.add_argument("--lr", type=float, default=1e-4)
    ap.add_argument("--weight-decay", type=float, default=1e-3)
    ap.add_argument("--focal-gamma", type=float, default=2.0)
    ap.add_argument("--pos-weight", type=float, default=None,
                    help="Class weight for the positive class. If None, computed "
                         "from train data as n_neg/n_pos.")
    ap.add_argument("--rc-prob", type=float, default=0.5,
                    help="Probability of reverse-complementing each window per epoch (training only).")
    ap.add_argument("--patience", type=int, default=3)
    ap.add_argument("--max-len", type=int, default=2048,
                    help="Tokenizer max_length / padding (default: 2048 = pre-training limit).")
    ap.add_argument("--grad-checkpoint", action="store_true", default=True,
                    help="Enable gradient checkpointing in BMFM (memory saver).")
    ap.add_argument("--no-grad-checkpoint", action="store_false", dest="grad_checkpoint")
    # Output
    ap.add_argument("--output-root", type=Path, required=True)
    ap.add_argument("--force", action="store_true",
                    help="Re-run even if metrics.json already exists in output_dir.")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    random.seed(args.seed)
    device = torch.device(args.device)

    mode = "frozen_live" if args.no_lora else "lora"
    upstream = "raw_bmfm_ref"
    out_dir = (args.output_root / args.label_mode / upstream / mode /
               args.aggregation / f"seed_{args.seed}")
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n[run] label={args.label_mode}  upstream={upstream}  mode={mode}  "
          f"agg={args.aggregation}  pool={args.per_window_pool}  seed={args.seed}")
    print(f"[run] output_dir: {out_dir}")

    # ── Skip-if-exists: respect prior completed runs ────────────────────────
    final_metrics = out_dir / "metrics.json"
    if final_metrics.exists() and not args.force:
        print(f"[run] SKIP — final metrics already exist at {final_metrics}")
        print(f"      (pass --force to overwrite)")
        return

    # ── Chrom map + window order ────────────────────────────────────────────
    chrom_of_window, chroms = build_chrom_map(args.windows)
    n_chroms = len(chroms)
    # window_order: stable list following biological (chrom-sorted) order from windows.tsv
    wins_df = pd.read_csv(args.windows, sep="\t")
    window_order = wins_df["window_id"].tolist()

    # ── Labels per split ────────────────────────────────────────────────────
    labels = load_labels(args.splits_root, args.seed, label_mode=args.label_mode)
    print("[data] split sizes (and y mean):")
    for split, df in labels.items():
        print(f"  {split}: n={len(df)}  y_mean={df['y'].mean():.3f}")

    # ── Read patient-window sequences ───────────────────────────────────────
    seq_df = pd.read_csv(args.sequences, usecols=["PTID", "window_id", "sequence"])
    seq_df["PTID"] = seq_df["PTID"].astype(str)
    print(f"[data] sequences: {len(seq_df)} rows ({seq_df.PTID.nunique()} patients × "
          f"{seq_df.window_id.nunique()} windows)")

    # ── Tokenizer + Model ───────────────────────────────────────────────────
    print("[model] loading BMFM-DNA-REF + tokenizer …")
    tok = load_tokenizer()
    tok.enable_padding(length=args.max_len, pad_id=2, pad_token="[PAD]")
    tok.enable_truncation(max_length=args.max_len)

    bmfm = load_bmfm_ref_with_lora(
        lora_r=args.lora_r, lora_alpha=args.lora_alpha,
        lora_dropout=args.lora_dropout,
        enable_lora=not args.no_lora,
        grad_checkpoint=args.grad_checkpoint,
    )

    hidden_dim = 768
    aggregator = Aggregator(args.aggregation, hidden_dim, n_chroms)
    head = MLPHead(aggregator.output_dim, hidden=256, dropout=0.3)
    model = PatientClassifier(bmfm, aggregator, head,
                              per_window_pool=args.per_window_pool,
                              micro_batch=args.micro_batch).to(device)
    # All aggregator+head params are trainable; LoRA params already marked trainable.

    # ── DataLoaders ─────────────────────────────────────────────────────────
    def make_loader(split: str, training: bool) -> DataLoader:
        df = labels[split]
        df = df[df["Patient_ID"].astype(str).isin(seq_df["PTID"].astype(str).unique())]
        pids = df["Patient_ID"].astype(str).tolist()
        lbl  = dict(zip(pids, df["y"].astype(int).tolist()))
        ds = PatientWindowsDataset(
            sequences_df=seq_df, patient_ids=pids, labels=lbl,
            window_order=window_order, chrom_of_window=chrom_of_window,
            tokenizer=tok, max_len=args.max_len, training=training,
            rc_prob=args.rc_prob,
        )
        return DataLoader(ds, batch_size=args.batch_patients,
                          shuffle=training, num_workers=0,
                          collate_fn=collate_patients)
    train_loader = make_loader("train", training=True)
    val_loader   = make_loader("val",   training=False)
    test_loader  = make_loader("test",  training=False)

    # ── Positive-class weight (computed from train if not provided) ─────────
    if args.pos_weight is None:
        n_pos = labels["train"]["y"].sum()
        n_neg = len(labels["train"]) - n_pos
        args.pos_weight = float(n_neg / max(n_pos, 1))
    print(f"[loss] focal γ={args.focal_gamma}  pos_weight={args.pos_weight:.3f}")

    # ── Train ───────────────────────────────────────────────────────────────
    t_train = time.time()
    model, fit_info, curve = train_one(
        model, train_loader, val_loader,
        epochs=args.epochs, lr=args.lr, weight_decay=args.weight_decay,
        focal_gamma=args.focal_gamma, pos_weight=args.pos_weight,
        patience=args.patience, device=device, out_dir=out_dir,
    )
    train_min = (time.time() - t_train) / 60
    print(f"[train] done in {train_min:.1f} min  best_val_auc={fit_info['best_val_auc']:.3f}")

    # ── Final eval ──────────────────────────────────────────────────────────
    val_logits, val_y, _   = run_inference(model, val_loader,  device)
    test_logits, test_y, _ = run_inference(model, test_loader, device)
    val_m  = patient_metrics(val_logits,  val_y)
    test_m = patient_metrics(test_logits, test_y)
    print(f"\nVAL  : {val_m}")
    print(f"TEST : {test_m}")

    with open(out_dir / "metrics.json", "w") as f:
        json.dump({
            "label_mode": args.label_mode,
            "upstream_tag": upstream,
            "mode": mode,
            "pool": args.per_window_pool,   # for column-compat with script 15
            "aggregation": args.aggregation,
            "seed": args.seed,
            "fit_info": fit_info,
            "lora": {"r": args.lora_r, "alpha": args.lora_alpha, "dropout": args.lora_dropout},
            "train_minutes": train_min,
            "val": val_m,
            "test": test_m,
        }, f, indent=2)
    pd.DataFrame(curve).to_csv(out_dir / "val_curve.csv", index=False)
    # Save only the trainable head + aggregator + LoRA params, not the frozen base
    trainable_sd = {k: v.detach().cpu() for k, v in model.state_dict().items()
                    if model.get_parameter(k).requires_grad if k in dict(model.named_parameters())}
    # Simpler: save head + aggregator state separately
    torch.save({
        "head": model.head.state_dict(),
        "aggregator": model.aggregator.state_dict(),
    }, out_dir / "head_aggregator.pt")
    if not args.no_lora:
        try:
            model.bmfm.save_pretrained(out_dir / "lora_adapter")
        except Exception as e:
            print(f"[warn] could not save_pretrained LoRA: {e}")
    print(f"\nWrote {out_dir / 'metrics.json'}")


if __name__ == "__main__":
    main()
