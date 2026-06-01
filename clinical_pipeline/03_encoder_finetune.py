"""
03_encoder_finetune.py
======================
Fine-tunes encoder-only transformer models (ModernBERT / BioClinical-ModernBERT)
for sequence classification on ADNI baseline clinical summaries.

Uses HuggingFace Trainer with:
  - Class-weighted cross-entropy loss (via Trainer subclass)
  - Early stopping (patience=5 by default)
  - fp16 mixed precision on GPU
  - Best checkpoint selection by val loss

Tasks
-----
  T1_binary     : CN vs MCI+AD  (binary)
  T2_multiclass : CN / MCI / AD (3-class)
  T3a_conv3y    : Conversion to AD within 3 years (binary, non-AD at baseline)
  T3b_conv5y    : Conversion to AD within 5 years (binary, non-AD at baseline)

Strategies
----------
  Full fine-tune  (default) : all parameters trained, lr=2e-5
  Frozen backbone            : --freeze_backbone flag, head only, lr=1e-3

Hyperparameters (per BioClinical-ModernBERT paper)
--------------------------------------------------
  Epochs: 10 (full FT) / 20 (frozen)
  Batch size: 16
  Weight decay: 1e-5
  Early stopping patience: 5

Usage
-----
  # Full fine-tune
  python 03_encoder_finetune.py --model_id answerdotai/ModernBERT-base --task T1_binary --seed 0

  # Frozen backbone
  python 03_encoder_finetune.py --model_id answerdotai/ModernBERT-base --task T1_binary --seed 0 --freeze_backbone
"""

import os, json, argparse, warnings
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

import torch
import torch.nn as nn
from torch.utils.data import Dataset

from transformers import (
    AutoTokenizer,
    AutoModelForSequenceClassification,
    Trainer,
    TrainingArguments,
    EarlyStoppingCallback,
)

from sklearn.utils.class_weight import compute_class_weight
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, roc_auc_score,
    f1_score, average_precision_score, precision_score, recall_score,
    confusion_matrix,
)

warnings.filterwarnings("ignore")


# ── Task definitions ──────────────────────────────────────────────────────────
TASK_CONFIG = {
    "T1_binary": {
        "label_col":     "Label_bl_multi",
        "num_labels":    2,
        "task_type":     "binary",
        "label_map":     {"CN": 0, "MCI": 1, "AD": 1},
        "filter_non_ad": False,
        "description":   "Binary: CN vs MCI+AD",
    },
    "T2_multiclass": {
        "label_col":     "Label_bl_multi",
        "num_labels":    3,
        "task_type":     "multiclass",
        "label_map":     {"CN": 0, "MCI": 1, "AD": 2},
        "filter_non_ad": False,
        "description":   "Multiclass: CN / MCI / AD",
    },
    "T3a_conv3y": {
        "label_col":     "Label_3y",
        "num_labels":    2,
        "task_type":     "binary",
        "label_map":     None,
        "filter_non_ad": True,
        "description":   "Prognosis: conversion to AD within 3 years",
    },
    "T3b_conv5y": {
        "label_col":     "Label_5y",
        "num_labels":    2,
        "task_type":     "binary",
        "label_map":     None,
        "filter_non_ad": True,
        "description":   "Prognosis: conversion to AD within 5 years",
    },
    "T3c_conv10y": {
        "label_col":     "Label_10y",
        "num_labels":    2,
        "task_type":     "binary",
        "label_map":     None,
        "filter_non_ad": True,
        "description":   "Prognosis: conversion to AD within 10 years",
    },
    "T3d_conv7y": {
        "label_col":     "Label_7y",
        "num_labels":    2,
        "task_type":     "binary",
        "label_map":     None,
        "filter_non_ad": True,
        "description":   "Prognosis: conversion to AD within 7 years",
    },
    # Stratified binary tasks: label derived from `conversion_group`. The label_map
    # only matches the relevant baseline cohort (other groups -> NaN -> dropped in
    # load_split), so no extra cohort filter is needed.
    "T1d_pmci_smci": {
        "label_col":     "conversion_group",
        "num_labels":    2,
        "task_type":     "binary",
        "label_map":     {"sMCI": 0, "pMCI": 1},
        "filter_non_ad": False,
        "description":   "Binary: pMCI vs sMCI (baseline MCI)",
    },
    "T1e_scn_pcn": {
        "label_col":     "conversion_group",
        "num_labels":    2,
        "task_type":     "binary",
        "label_map":     {"sCN": 0, "pCN_to_MCI": 1, "pCN_to_AD": 1},
        "filter_non_ad": False,
        "description":   "Binary: sCN vs pCN (baseline CN, progression to MCI or AD)",
    },
    "T1b_cnmci_ad": {
        "label_col":     "Label_bl_multi",
        "num_labels":    2,
        "task_type":     "binary",
        "label_map":     {"CN": 0, "MCI": 0, "AD": 1},
        "filter_non_ad": False,
        "description":   "Binary: CN+MCI vs AD",
    },
    "T1c_scn_prog": {
        "label_col":     "conversion_group",
        "num_labels":    2,
        "task_type":     "binary",
        "label_map":     {"sCN": 0, "pCN_to_MCI": 1, "pCN_to_AD": 1, "pMCI": 1, "AD_bl": 1},
        "filter_non_ad": False,
        "description":   "Binary: sCN vs progressors+AD (pCN/pMCI/AD, excl. sMCI)",
    },
}


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--model_id",        type=str,   required=True)
    p.add_argument("--task",            type=str,   required=True,
                   choices=list(TASK_CONFIG.keys()))
    p.add_argument("--seed",            type=int,   required=True, choices=[0, 1, 2])
    p.add_argument("--freeze_backbone", action="store_true",
                   help="Freeze encoder; train classifier head only")
    p.add_argument("--text_col",        type=str,   default="Generated_Text")
    p.add_argument("--max_length",      type=int,   default=1024)
    p.add_argument("--epochs",          type=int,   default=None,
                   help="Default: 10 (full FT) or 20 (frozen), per BioClinical-ModernBERT paper")
    p.add_argument("--patience",        type=int,   default=5,
                   help="Early stopping patience (epochs without val_loss improvement)")
    p.add_argument("--batch_size",      type=int,   default=16)
    p.add_argument("--lr",              type=float, default=None,
                   help="Default: 2e-5 (full FT) or 1e-3 (frozen)")
    p.add_argument("--warmup_ratio",    type=float, default=0.1)
    p.add_argument("--weight_decay",    type=float, default=1e-5,
                   help="Weight decay (default: 1e-5, per BioClinical-ModernBERT paper)")
    p.add_argument("--data_dir",        type=str,
                   default=r"D:\ADNI_BIDS_project\derivatives\clinical\verbose\baseline")
    p.add_argument("--out_dir",         type=str,
                   default=r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\clinical_pipeline\outputs")
    p.add_argument("--hf_cache",        type=str,   default=None,
                   help="HuggingFace cache dir (overrides HF_HOME)")
    p.add_argument("--wandb",           action="store_true",
                   help="Log this run to Weights & Biases in OFFLINE mode (sync later)")
    p.add_argument("--wandb_project",   type=str,   default="clinical-encoder-post-exclusion")
    return p.parse_args()


# ── Dataset ───────────────────────────────────────────────────────────────────
class ClinicalDataset(Dataset):
    def __init__(self, texts, labels, tokenizer, max_length):
        self.encodings = tokenizer(
            texts,
            truncation=True,
            max_length=max_length,
            padding="max_length",
            return_tensors="pt",
        )
        self.labels = torch.tensor(labels, dtype=torch.long)

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        return {
            "input_ids":      self.encodings["input_ids"][idx],
            "attention_mask": self.encodings["attention_mask"][idx],
            "labels":         self.labels[idx],
        }


# ── Data loading ──────────────────────────────────────────────────────────────
def load_split(csv_path, task_cfg, text_col):
    """Load a CSV split, apply label mapping and baseline-AD filter."""
    df = pd.read_csv(csv_path, low_memory=False)
    label_col = task_cfg["label_col"]

    if task_cfg["label_map"] is not None:
        df[label_col] = df[label_col].map(task_cfg["label_map"])
    else:
        df[label_col] = pd.to_numeric(df[label_col], errors="coerce")

    if task_cfg["filter_non_ad"] and "Label_bl_multi" in df.columns:
        df = df[df["Label_bl_multi"].isin(["CN", "MCI"])]

    df = df.dropna(subset=[label_col, text_col])
    df[label_col] = df[label_col].astype(int)
    return df[text_col].tolist(), df[label_col].tolist()


# ── Weighted-loss Trainer ─────────────────────────────────────────────────────
class WeightedTrainer(Trainer):
    """Trainer subclass that uses class-weighted cross-entropy loss."""

    def __init__(self, class_weights, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.class_weights = class_weights   # tensor on CPU; moved to device in compute_loss

    def compute_loss(self, model, inputs, return_outputs=False, **kwargs):
        labels = inputs.pop("labels")
        outputs = model(**inputs)
        logits = outputs.logits
        loss_fn = nn.CrossEntropyLoss(
            weight=self.class_weights.to(logits.device)
        )
        loss = loss_fn(logits, labels)
        return (loss, outputs) if return_outputs else loss


# ── Metrics computation ──────────────────────────────────────────────────────
def make_compute_metrics(task_type):
    """Returns a compute_metrics function compatible with HuggingFace Trainer."""
    def compute_metrics(eval_pred):
        logits, labels = eval_pred
        probs = torch.softmax(torch.tensor(logits), dim=-1).numpy()
        preds = np.argmax(logits, axis=-1)

        m = {
            "accuracy":     float(accuracy_score(labels, preds)),
            "balanced_acc": float(balanced_accuracy_score(labels, preds)),
        }
        if task_type == "binary":
            m["precision"] = float(precision_score(labels, preds, zero_division=0))
            m["recall"]    = float(recall_score(labels, preds, zero_division=0))
            m["f1"]        = float(f1_score(labels, preds, zero_division=0))
            # Sensitivity = recall of the positive class; specificity = recall of the
            # negative class. Measured directly via recall_score per class.
            m["sensitivity"] = float(recall_score(labels, preds, pos_label=1, zero_division=0))
            m["specificity"] = float(recall_score(labels, preds, pos_label=0, zero_division=0))
            try:
                m["auc_roc"] = float(roc_auc_score(labels, probs[:, 1]))
                m["auc_pr"]  = float(average_precision_score(labels, probs[:, 1]))
            except Exception:
                m["auc_roc"] = float("nan")
                m["auc_pr"]  = float("nan")
        else:
            m["precision_macro"] = float(precision_score(labels, preds, average="macro", zero_division=0))
            m["recall_macro"]    = float(recall_score(labels, preds, average="macro", zero_division=0))
            m["macro_f1"]        = float(f1_score(labels, preds, average="macro", zero_division=0))
            try:
                m["auc_roc_ovr"]  = float(roc_auc_score(labels, probs, multi_class="ovr", average="macro"))
                m["auc_pr_macro"] = float(average_precision_score(labels, probs, average="macro"))
            except Exception:
                m["auc_roc_ovr"]  = float("nan")
                m["auc_pr_macro"] = float("nan")
        return m
    return compute_metrics


def full_metrics(labels, logits, task_type):
    """
    Full metric set for a final (val/test) evaluation, including the confusion matrix.
    Used for reporting + wandb logging (not for per-epoch Trainer eval).
    Returns (metrics_dict, preds, labels) so callers can log a wandb confusion matrix.
    """
    labels = np.asarray(labels)
    probs = torch.softmax(torch.tensor(logits), dim=-1).numpy()
    preds = np.argmax(logits, axis=-1)

    m = {
        "accuracy":     float(accuracy_score(labels, preds)),
        "balanced_acc": float(balanced_accuracy_score(labels, preds)),
    }
    if task_type == "binary":
        m["precision"]   = float(precision_score(labels, preds, zero_division=0))
        m["recall"]      = float(recall_score(labels, preds, zero_division=0))
        m["f1"]          = float(f1_score(labels, preds, zero_division=0))
        # Directly measured per-class recalls.
        m["sensitivity"] = float(recall_score(labels, preds, pos_label=1, zero_division=0))
        m["specificity"] = float(recall_score(labels, preds, pos_label=0, zero_division=0))
        try:
            m["auc_roc"] = float(roc_auc_score(labels, probs[:, 1]))
            m["auc_pr"]  = float(average_precision_score(labels, probs[:, 1]))
        except Exception:
            m["auc_roc"] = float("nan")
            m["auc_pr"]  = float("nan")
        cm = confusion_matrix(labels, preds, labels=[0, 1])
        tn, fp, fn, tp = (int(x) for x in cm.ravel())
        m["confusion_matrix"] = {"tn": tn, "fp": fp, "fn": fn, "tp": tp}
    else:
        m["macro_precision"] = float(precision_score(labels, preds, average="macro", zero_division=0))
        m["macro_recall"]    = float(recall_score(labels, preds, average="macro", zero_division=0))
        m["macro_f1"]        = float(f1_score(labels, preds, average="macro", zero_division=0))
        m["per_class_recall"] = [float(x) for x in
                                 recall_score(labels, preds, average=None, zero_division=0)]
        try:
            m["auc_roc_ovr"]  = float(roc_auc_score(labels, probs, multi_class="ovr", average="macro"))
            m["auc_pr_macro"] = float(average_precision_score(labels, probs, average="macro"))
        except Exception:
            m["auc_roc_ovr"]  = float("nan")
            m["auc_pr_macro"] = float("nan")
        m["confusion_matrix"] = confusion_matrix(labels, preds).tolist()

    m = {k: (round(v, 4) if isinstance(v, float) else v) for k, v in m.items()}
    return m, preds, labels


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()

    # Strategy-dependent defaults
    if args.lr     is None: args.lr     = 1e-3 if args.freeze_backbone else 2e-5
    if args.epochs is None: args.epochs = 20   if args.freeze_backbone else 10

    strategy   = "frozen" if args.freeze_backbone else "full_ft"
    task_cfg   = TASK_CONFIG[args.task]
    model_slug = args.model_id.split("/")[-1]
    hf_cache   = args.hf_cache or os.environ.get("HF_HOME", None)

    # Output dir: outputs/<model>/<task>/seed_<N>/<strategy>/
    out_dir = Path(args.out_dir) / model_slug / args.task / f"seed_{args.seed}" / strategy
    out_dir.mkdir(parents=True, exist_ok=True)

    # Trainer writes checkpoints into a subdir — we'll clean up after
    trainer_dir = out_dir / "trainer_ckpts"

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    print(f"\n{'='*65}")
    print(f"  Model    : {args.model_id}")
    print(f"  Task     : {task_cfg['description']}")
    print(f"  Seed     : {args.seed}  |  Strategy: {strategy}")
    print(f"  Device   : {device}  |  Epochs: {args.epochs}  |  LR: {args.lr}")
    print(f"  Patience : {args.patience}")
    print(f"  Output   : {out_dir}")
    print(f"{'='*65}\n")

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    # ── Load data ─────────────────────────────────────────────────────────────
    data_dir = Path(args.data_dir) / f"seed_{args.seed}"
    train_texts, train_labels = load_split(data_dir / "train.csv", task_cfg, args.text_col)
    val_texts,   val_labels   = load_split(data_dir / "val.csv",   task_cfg, args.text_col)
    test_texts,  test_labels  = load_split(data_dir / "test.csv",  task_cfg, args.text_col)

    print(f"  Train: {len(train_labels)}  Val: {len(val_labels)}  Test: {len(test_labels)}")
    unique, counts = np.unique(train_labels, return_counts=True)
    print(f"  Train class dist: {dict(zip(unique.tolist(), counts.tolist()))}")

    # ── Tokenise ──────────────────────────────────────────────────────────────
    print(f"\n  Loading tokeniser ...")
    tokenizer = AutoTokenizer.from_pretrained(args.model_id, cache_dir=hf_cache)

    train_ds = ClinicalDataset(train_texts, train_labels, tokenizer, args.max_length)
    val_ds   = ClinicalDataset(val_texts,   val_labels,   tokenizer, args.max_length)
    test_ds  = ClinicalDataset(test_texts,  test_labels,  tokenizer, args.max_length)

    # ── Model ─────────────────────────────────────────────────────────────────
    print(f"  Loading model ...")
    model = AutoModelForSequenceClassification.from_pretrained(
        args.model_id,
        num_labels=task_cfg["num_labels"],
        cache_dir=hf_cache,
        ignore_mismatched_sizes=True,
    )

    if args.freeze_backbone:
        n_frozen = 0
        for name, param in model.named_parameters():
            if "classifier" not in name:
                param.requires_grad = False
                n_frozen += 1
        trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
        print(f"  Backbone FROZEN ({n_frozen} param tensors frozen). Trainable: {trainable:,}")
    else:
        trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
        print(f"  Full fine-tune. Trainable parameters: {trainable:,}")

    # ── Class weights ─────────────────────────────────────────────────────────
    classes = np.unique(train_labels)
    cw = compute_class_weight("balanced", classes=classes, y=np.array(train_labels))
    class_weights = torch.tensor(cw, dtype=torch.float)
    print(f"  Class weights: {dict(zip(classes.tolist(), cw.round(3).tolist()))}")

    # ── Weights & Biases (offline) ────────────────────────────────────────────
    use_wandb = args.wandb
    if use_wandb:
        try:
            import wandb
        except ImportError:
            print("  [WARN] wandb not installed — disabling wandb logging.")
            use_wandb = False
    if use_wandb:
        os.environ["WANDB_MODE"] = "offline"           # sync on HPC later
        run_name = f"{model_slug}_{strategy}_{args.task}_seed{args.seed}"
        wandb.init(
            project=args.wandb_project,
            name=run_name,
            group=args.task,
            job_type=strategy,
            tags=[model_slug, args.task, strategy, f"seed{args.seed}"],
            dir=os.environ.get("WANDB_DIR", str(out_dir)),
            reinit=True,
            config={
                "model_id": args.model_id, "model": model_slug, "task": args.task,
                "task_description": task_cfg["description"], "seed": args.seed,
                "strategy": strategy, "epochs": args.epochs, "lr": args.lr,
                "weight_decay": args.weight_decay, "batch_size": args.batch_size,
                "max_length": args.max_length, "warmup_ratio": args.warmup_ratio,
                "patience": args.patience, "num_labels": task_cfg["num_labels"],
                "n_train": len(train_labels), "n_val": len(val_labels),
                "n_test": len(test_labels),
            },
        )
        print(f"  wandb offline run: {run_name}")

    # ── Training arguments ────────────────────────────────────────────────────
    training_args = TrainingArguments(
        output_dir=str(trainer_dir),
        num_train_epochs=args.epochs,
        per_device_train_batch_size=args.batch_size,
        per_device_eval_batch_size=args.batch_size,
        learning_rate=args.lr,
        weight_decay=args.weight_decay,
        warmup_ratio=args.warmup_ratio,
        # Evaluation & checkpointing every epoch
        eval_strategy="epoch",
        save_strategy="epoch",
        load_best_model_at_end=True,
        # eval_dataset is a dict below, so HF prefixes metrics as eval_<key>_* .
        # Select the best checkpoint on the VALIDATION loss specifically.
        metric_for_best_model="eval_val_loss",
        greater_is_better=False,
        save_total_limit=2,           # keep only best + latest checkpoint
        # Mixed precision
        fp16=torch.cuda.is_available(),
        # Logging
        logging_steps=10,
        logging_dir=str(out_dir / "tb_logs"),
        report_to=(["wandb"] if use_wandb else "none"),   # offline wandb if enabled
        # Reproducibility
        seed=args.seed,
        data_seed=args.seed,
        # Misc
        remove_unused_columns=False,
        dataloader_num_workers=2,
        dataloader_pin_memory=True,
    )

    # ── Trainer ───────────────────────────────────────────────────────────────
    trainer = WeightedTrainer(
        class_weights=class_weights,
        model=model,
        args=training_args,
        train_dataset=train_ds,
        # Evaluate BOTH train and val every epoch → wandb logs true train-vs-val
        # learning curves (eval_train_* and eval_val_*) so over/under-fitting is
        # visible. VAL drives early stopping + best-checkpoint selection; TRAIN is
        # logged for the curve only. (Test is held out — evaluated once at the end.)
        eval_dataset={"val": val_ds, "train": train_ds},
        compute_metrics=make_compute_metrics(task_cfg["task_type"]),
        callbacks=[EarlyStoppingCallback(early_stopping_patience=args.patience)],
    )

    # ── Train ─────────────────────────────────────────────────────────────────
    print(f"\n  Starting training ...")
    train_result = trainer.train()

    # ── Final evaluation: full metrics on val + test (incl. confusion matrix) ──
    print(f"\n  Final evaluation (val + test) ...")
    val_out  = trainer.predict(val_ds,  metric_key_prefix="val")
    test_out = trainer.predict(test_ds, metric_key_prefix="test")

    val_metrics,  _,         _          = full_metrics(
        val_out.label_ids,  val_out.predictions,  task_cfg["task_type"])
    test_metrics, test_preds, test_true = full_metrics(
        test_out.label_ids, test_out.predictions, task_cfg["task_type"])
    val_metrics["loss"]  = round(float(val_out.metrics.get("val_loss", float("nan"))), 4)
    test_metrics["loss"] = round(float(test_out.metrics.get("test_loss", float("nan"))), 4)

    print(f"\n  Test metrics:")
    for k, v in test_metrics.items():
        print(f"    {k}: {v}")

    # ── wandb: log final val/test metrics + confusion matrices ─────────────────
    if use_wandb:
        wandb.summary.update({f"val/{k}":  v for k, v in val_metrics.items()
                              if isinstance(v, (int, float))})
        wandb.summary.update({f"test/{k}": v for k, v in test_metrics.items()
                              if isinstance(v, (int, float))})
        try:
            cls = [str(c) for c in sorted(set(test_true.tolist()))]
            wandb.log({"test/confusion_matrix": wandb.plot.confusion_matrix(
                y_true=test_true.tolist(), preds=test_preds.tolist(), class_names=cls)})
        except Exception as e:
            print(f"  [WARN] wandb confusion-matrix log failed: {e}")
        wandb.finish()

    # ── Save outputs ──────────────────────────────────────────────────────────
    # Training log from Trainer
    log_history = pd.DataFrame(trainer.state.log_history)
    log_history.to_csv(out_dir / "training_log.csv", index=False)

    # Run config + val/test metrics
    config = {
        "model_id":         args.model_id,
        "task":             args.task,
        "task_description": task_cfg["description"],
        "seed":             args.seed,
        "strategy":         strategy,
        "epochs":           args.epochs,
        "best_epoch":       int(trainer.state.best_model_checkpoint.split("-")[-1])
                            if trainer.state.best_model_checkpoint else args.epochs,
        "lr":               args.lr,
        "weight_decay":     args.weight_decay,
        "batch_size":       args.batch_size,
        "max_length":       args.max_length,
        "patience":         args.patience,
        "n_train":          len(train_labels),
        "n_val":            len(val_labels),
        "n_test":           len(test_labels),
        "class_weights":    dict(zip(classes.tolist(), cw.round(4).tolist())),
        "timestamp":        datetime.now().isoformat(),
    }

    results = {"config": config, "val_metrics": val_metrics, "test_metrics": test_metrics}
    with open(out_dir / "metrics.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n  Saved → {out_dir / 'metrics.json'}")
    print(f"  Saved → {out_dir / 'training_log.csv'}")

    # Save best model checkpoint
    best_ckpt_dir = out_dir / "best_checkpoint"
    model.save_pretrained(best_ckpt_dir)
    tokenizer.save_pretrained(best_ckpt_dir)
    print(f"  Best model saved → {best_ckpt_dir}")

    # Clean up intermediate Trainer checkpoints to save disk space
    import shutil
    if trainer_dir.exists():
        shutil.rmtree(trainer_dir, ignore_errors=True)
        print(f"  Cleaned up trainer checkpoints.")

    print(f"\n  Done.")


if __name__ == "__main__":
    main()
