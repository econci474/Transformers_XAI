"""
peft_utils.py
==============
Shared utilities for NT v2 LoRA fine-tuning scripts (05, 06, 07).
"""

import json
import pathlib
from typing import Optional

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset
from transformers import AutoTokenizer, AutoModelForMaskedLM
from peft import LoraConfig, get_peft_model, TaskType


MODEL_ID = "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species"


# ── Focal Loss ────────────────────────────────────────────────────────────────
class FocalLoss(nn.Module):
    """Focal loss for multi-class classification (Lin et al., 2017).
    
    Addresses class imbalance by down-weighting well-classified examples.
    FL(p_t) = -alpha_t * (1 - p_t)^gamma * log(p_t)
    """
    def __init__(self, gamma: float = 2.0, weight: Optional[torch.Tensor] = None):
        super().__init__()
        self.gamma = gamma
        self.weight = weight

    def forward(self, logits: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
        ce = nn.functional.cross_entropy(logits, targets, weight=self.weight, reduction="none")
        pt = torch.exp(-ce)
        focal = ((1 - pt) ** self.gamma) * ce
        return focal.mean()


# ── NT v2 wrapper with task head ──────────────────────────────────────────────
class NTv2ForSequenceTask(nn.Module):
    """NT v2 backbone + mean-pooling + task-specific head.
    
    The backbone is the encoder from AutoModelForMaskedLM (we discard the LM head).
    LoRA is applied to the backbone via peft.
    """
    def __init__(self, num_labels: int, is_regression: bool = False,
                 lora_r: int = 8, lora_alpha: int = 16, lora_dropout: float = 0.1):
        super().__init__()

        # Load backbone
        self.backbone = AutoModelForMaskedLM.from_pretrained(
            MODEL_ID, trust_remote_code=True
        )

        # Identify target modules for LoRA
        # NT v2 uses standard transformer attention with q_proj, v_proj
        lora_config = LoraConfig(
            r=lora_r,
            lora_alpha=lora_alpha,
            lora_dropout=lora_dropout,
            target_modules=["query", "value"],  # NT v2 naming
            bias="none",
            task_type=TaskType.FEATURE_EXTRACTION,
        )
        self.backbone = get_peft_model(self.backbone, lora_config)
        self.backbone.print_trainable_parameters()

        # Get hidden dim from config
        self.hidden_dim = self.backbone.config.hidden_size

        # Task head
        self.is_regression = is_regression
        if is_regression:
            self.head = nn.Linear(self.hidden_dim, 1)
        else:
            self.head = nn.Sequential(
                nn.Linear(self.hidden_dim, 256),
                nn.ReLU(),
                nn.Dropout(0.1),
                nn.Linear(256, num_labels),
            )

    def forward(self, input_ids, attention_mask):
        out = self.backbone(
            input_ids,
            attention_mask=attention_mask,
            encoder_attention_mask=attention_mask,
            output_hidden_states=True,
        )
        hidden = out.hidden_states[-1]  # (B, L, H)

        # Mean-pool
        mask = attention_mask.unsqueeze(-1).float()
        pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1)  # (B, H)

        return self.head(pooled)


# ── Datasets ──────────────────────────────────────────────────────────────────
class SequenceLabelDataset(Dataset):
    """Dataset for sequence → label (classification or regression)."""

    def __init__(self, sequences: list, labels: list, tokenizer,
                 max_length: int = 2048, label_dtype=torch.long):
        self.sequences = sequences
        self.labels = labels
        self.tokenizer = tokenizer
        self.max_length = max_length
        self.label_dtype = label_dtype

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        seq = self.sequences[idx]
        label = self.labels[idx]

        encoded = self.tokenizer(
            seq,
            max_length=self.max_length,
            padding="max_length",
            truncation=True,
            return_tensors="pt",
        )
        input_ids = encoded["input_ids"].squeeze(0)
        attention_mask = (input_ids != self.tokenizer.pad_token_id).long()

        return {
            "input_ids": input_ids,
            "attention_mask": attention_mask,
            "label": torch.tensor(label, dtype=self.label_dtype),
        }


# ── Training loop ─────────────────────────────────────────────────────────────
def train_epoch(model, loader, optimizer, criterion, device, grad_accum=1):
    model.train()
    total_loss = 0
    total_samples = 0
    optimizer.zero_grad()

    for step, batch in enumerate(loader):
        input_ids = batch["input_ids"].to(device)
        attention_mask = batch["attention_mask"].to(device)
        labels = batch["label"].to(device)

        logits = model(input_ids, attention_mask)
        if logits.dim() == 2 and logits.size(1) == 1:
            logits = logits.squeeze(-1)
        loss = criterion(logits, labels) / grad_accum

        loss.backward()

        if (step + 1) % grad_accum == 0:
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            optimizer.zero_grad()

        total_loss += loss.item() * grad_accum * len(labels)
        total_samples += len(labels)

    # Flush remaining grads
    if (step + 1) % grad_accum != 0:
        optimizer.step()
        optimizer.zero_grad()

    return total_loss / total_samples


@torch.no_grad()
def eval_epoch(model, loader, criterion, device):
    model.eval()
    total_loss = 0
    total_samples = 0
    all_preds = []
    all_labels = []

    for batch in loader:
        input_ids = batch["input_ids"].to(device)
        attention_mask = batch["attention_mask"].to(device)
        labels = batch["label"].to(device)

        logits = model(input_ids, attention_mask)
        if logits.dim() == 2 and logits.size(1) == 1:
            logits = logits.squeeze(-1)
        loss = criterion(logits, labels)

        total_loss += loss.item() * len(labels)
        total_samples += len(labels)

        if logits.dim() == 1:
            all_preds.extend(logits.cpu().numpy().tolist())
        else:
            all_preds.extend(logits.argmax(dim=1).cpu().numpy().tolist())
        all_labels.extend(labels.cpu().numpy().tolist())

    return total_loss / total_samples, np.array(all_preds), np.array(all_labels)


def save_lora_checkpoint(model, out_dir, name="best"):
    """Save only the LoRA adapter weights (not the full model)."""
    save_dir = pathlib.Path(out_dir) / f"{name}_lora"
    model.backbone.save_pretrained(save_dir)
    print(f"  Saved LoRA adapter to {save_dir}")


def save_metrics(metrics: dict, out_dir: pathlib.Path, filename="metrics.json"):
    with open(out_dir / filename, "w") as f:
        json.dump(metrics, f, indent=2)
