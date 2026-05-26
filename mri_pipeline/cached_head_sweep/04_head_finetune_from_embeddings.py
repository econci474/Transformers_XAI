"""
04_head_finetune_from_embeddings.py
====================================
Model-agnostic head-only fine-tuner: trains a 2-layer MLP classification
head on cached frozen-encoder embeddings (BrainDINO / BrainMVP / ViT-MAE
all use the same head architecture, just different embed_dim).

Head architecture (mirrors BrainDINOClassifier.head):
    LayerNorm(D) -> Linear(D, D) -> GELU -> Dropout -> Linear(D, n_classes)

Training:
- Adam, weight_decay=1e-5 (paper default; not swept).
- ReduceLROnPlateau(mode=max, factor=0.5, patience=10) on val_bacc.
- EarlyStopping at --patience epochs without val_bacc improvement.
- Class-weighted CE with --label_smoothing.
- Max --epochs (default 50 -- head converges fast on cached features).

Metrics.json schema: identical to 02_supervised_finetuning_BrainDINO.py
output, so 05_aggregate_mri_results.py and 06_render_cross_model_table.py
ingest unchanged.

Usage example:
    python 04_head_finetune_from_embeddings.py \\
        --model_name BrainDINO \\
        --embeddings_pt <out>/braindino/braindino_pooled.pt \\
        --embed_dim 768 \\
        --task T1c_binary --seed 0 \\
        --lr 1e-4 --weight_decay 1e-5 \\
        --drop_rate 0.2 --label_smoothing 0.0 \\
        --data_dir <clinical_splits> \\
        --matched_labels_csv <...> \\
        --out_dir <braindino_outputs/aug_none/BrainDINO_vitb16_frozen_cached> \\
        --wandb --wandb_project braindino_frozen_cached
"""

import argparse
import importlib.util
import json
import os
import sys
import warnings
from datetime import datetime
from pathlib import Path

import torch  # noqa: F401  (Windows MKL DLL order)
import numpy as np
import pandas as pd
import torch.nn as nn
from sklearn.metrics import balanced_accuracy_score
from sklearn.utils.class_weight import compute_class_weight
from torch.utils.data import Dataset, DataLoader

THIS_DIR = Path(__file__).resolve().parent
MRI_DIR = THIS_DIR.parent
REPO_ROOT = MRI_DIR.parent
sys.path.insert(0, str(MRI_DIR))
sys.path.insert(0, str(REPO_ROOT))

_vit_spec = importlib.util.spec_from_file_location(
    "_vit_pipeline", MRI_DIR / "04_supervised_finetuning_ViT.py")
_vit = importlib.util.module_from_spec(_vit_spec)
sys.modules["_vit_pipeline"] = _vit
_vit_spec.loader.exec_module(_vit)

TASK_CONFIG            = _vit.TASK_CONFIG
cohort_label           = _vit.cohort_label
class_names_for        = _vit.class_names_for
load_matched_labels    = _vit.load_matched_labels
session_to_months      = _vit.session_to_months
patient_id_to_bids_sub = _vit.patient_id_to_bids_sub
compute_test_metrics   = _vit.compute_test_metrics
compute_diagnostics    = _vit.compute_diagnostics
from bidsification.exclusions import is_excluded_subject

warnings.filterwarnings("ignore")


# ── Head + dataset ────────────────────────────────────────────────────────────
class HeadOnly(nn.Module):
    """LayerNorm -> Linear -> GELU -> Dropout -> Linear. Matches the
    architecture used by BrainDINOClassifier / BrainMVPClassifier /
    ViT-MAE downstream heads."""
    def __init__(self, embed_dim: int, n_classes: int, drop_rate: float = 0.2):
        super().__init__()
        self.net = nn.Sequential(
            nn.LayerNorm(embed_dim),
            nn.Linear(embed_dim, embed_dim),
            nn.GELU(),
            nn.Dropout(drop_rate),
            nn.Linear(embed_dim, n_classes),
        )
    def forward(self, x):
        return self.net(x)


class EmbeddingDataset(Dataset):
    """In-memory dataset: takes a (sub, ses, label) DataFrame and the
    cached embedding tensor + ids index, returns (embedding, label)."""
    def __init__(self, df: pd.DataFrame, embeddings: torch.Tensor,
                 id_to_row: dict):
        self.embeddings = embeddings
        self.df = df.reset_index(drop=True)
        # Pre-compute row indices into the embedding tensor
        rows = []
        for r in df.itertuples():
            key = (r.bids_sub, r.bids_ses)
            if key not in id_to_row:
                raise KeyError(
                    f"Embedding missing for {key} -- re-run 03_extract_*_embeddings.py")
            rows.append(id_to_row[key])
        self.rows = torch.tensor(rows, dtype=torch.long)
        self.labels = torch.tensor(df["label"].astype(int).to_numpy(),
                                   dtype=torch.long)

    def __len__(self):
        return len(self.df)

    def __getitem__(self, i):
        return self.embeddings[self.rows[i]], self.labels[i]


# ── Load splits + join with embeddings ────────────────────────────────────────
def load_split_cached(baseline_csv: Path, matched_df: pd.DataFrame,
                      task_cfg: dict, id_to_row: dict,
                      session_filter=None, max_months=None):
    """Same load_split logic as the BrainMVP / BrainDINO trainers, but
    keeps only rows whose embedding is present in id_to_row."""
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

    # Derive bids_sub and keep only rows with cached embeddings
    df["bids_sub"] = df["Patient_ID"].apply(patient_id_to_bids_sub)
    key_series = df.apply(lambda r: (r["bids_sub"], r["bids_ses"]), axis=1)
    df["_has_embedding"] = key_series.map(lambda k: k in id_to_row)
    n_missing = int((~df["_has_embedding"]).sum())
    df = df[df["_has_embedding"]].drop(columns="_has_embedding").copy()

    return df.reset_index(drop=True)[
        ["Patient_ID", "bids_sub", "bids_ses", "label"]], n_missing


# ── wandb (local stub) ────────────────────────────────────────────────────────
def init_wandb(args, task_cfg, extra=None):
    if not args.wandb:
        return None
    try:
        import wandb
    except ImportError:
        print("  [WARN] --wandb set but wandb not installed; continuing local-only.")
        return None
    # Associate this run with a registered W&B Sweep, if --sweep_id was
    # passed (or WANDB_SWEEP_ID is already in env). Must be set BEFORE
    # wandb.init -- the run dir stamps the sweep_id and `wandb sync`
    # later uploads the run as a member of the sweep, so the sweep page
    # gets native parallel-coords / HP-vs-metric visualisation. The
    # sweep itself must be pre-registered once via `wandb sweep
    # sweep_<model>.yaml` on a login node (needs internet).
    if args.sweep_id:
        os.environ["WANDB_SWEEP_ID"] = args.sweep_id
        print(f"  W&B sweep: associating with {args.sweep_id}")
    # Group: clusters all cells of one model/sweep together in the W&B UI
    # (works even without a formal sweep registered).
    group = args.wandb_group or f"head_sweep_{args.model_name}"
    name = args.wandb_run_name or (
        f"{args.task}-s{args.seed}-cached"
        f"-lr{args.lr:.0e}_d{args.drop_rate}_ls{args.label_smoothing}")
    cfg = {
        "model_name":      args.model_name,
        "task":            args.task,
        "seed":            args.seed,
        "lr":              args.lr,
        "weight_decay":    args.weight_decay,
        "drop_rate":       args.drop_rate,
        "label_smoothing": args.label_smoothing,
        "epochs":          args.epochs,
        "patience":        args.patience,
        "embed_dim":       args.embed_dim,
        "strategy":        "frozen_cached",
        "augment":         "none",
        "embeddings_pt":   args.embeddings_pt,
    }
    if extra:
        cfg.update(extra)
    return wandb.init(project=args.wandb_project, entity=args.wandb_entity,
                      name=name, group=group, config=cfg, reinit=True)


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--model_name", type=str, required=True,
                   choices=["BrainDINO", "BrainMVP", "ViT-MAE75"],
                   help="Only used in metrics.json + wandb naming.")
    p.add_argument("--embeddings_pt", type=str, required=True,
                   help="Path to <model>_pooled.pt from 03_extract_*.py")
    p.add_argument("--embed_dim", type=int, required=True,
                   help="Embedding dim (768 for BrainDINO/ViT-MAE, 512 for BrainMVP).")
    p.add_argument("--task", type=str, required=True,
                   choices=list(TASK_CONFIG.keys()))
    p.add_argument("--seed", type=int, required=True, choices=[0, 1, 2])
    p.add_argument("--lr", type=float, default=1e-4)
    p.add_argument("--weight_decay", type=float, default=1e-5)
    p.add_argument("--drop_rate", type=float, default=0.2)
    p.add_argument("--label_smoothing", type=float, default=0.0)
    p.add_argument("--epochs", type=int, default=50)
    p.add_argument("--patience", type=int, default=15,
                   help="Early-stopping epochs without val_bacc improvement.")
    p.add_argument("--scheduler_patience", type=int, default=10)
    p.add_argument("--scheduler_factor",   type=float, default=0.5)
    p.add_argument("--batch_size", type=int, default=32,
                   help="Embeddings are tiny (768-d); large batches are fine.")
    session_group = p.add_mutually_exclusive_group()
    session_group.add_argument("--session", type=str, default="bl")
    session_group.add_argument("--long", type=str, default=None)
    p.add_argument("--matched_labels_csv", type=str, required=True)
    p.add_argument("--data_dir", type=str, required=True)
    p.add_argument("--out_dir", type=str, required=True,
                   help="Output dir for metrics.json / best_model.pt / etc.")
    p.add_argument("--no_resume", action="store_true")
    p.add_argument("--wandb", action="store_true")
    p.add_argument("--wandb_project", type=str, default=None)
    p.add_argument("--wandb_entity",  type=str, default=None)
    p.add_argument("--wandb_run_name", type=str, default=None)
    p.add_argument("--wandb_group", type=str, default=None,
                   help="W&B run group label (defaults to head_sweep_{model_name}). "
                        "All cells of one model cluster under this group in the UI.")
    p.add_argument("--sweep_id", type=str, default=None,
                   help="W&B Sweep ID returned by `wandb sweep <yaml>`. If set, "
                        "WANDB_SWEEP_ID env var is exported before wandb.init so the "
                        "offline run is stamped with the sweep_id and `wandb sync` "
                        "later uploads it as a member of the registered sweep.")
    args = p.parse_args()

    if args.long is None:
        args.long_mode = None; args.max_months = None
    else:
        if args.long.lower() == "all" or args.long == "0":
            args.long_mode = "all"; args.max_months = None
        else:
            n = int(args.long); args.long_mode = "cutoff"; args.max_months = n * 12
        args.session = None
    return args


# ── Eval helper ───────────────────────────────────────────────────────────────
def evaluate(model, loader, criterion, device):
    model.eval()
    total_loss, total_n = 0.0, 0
    all_logits, all_labels = [], []
    with torch.no_grad():
        for emb, lab in loader:
            emb = emb.to(device, non_blocking=True).float()
            lab = lab.to(device, non_blocking=True).long()
            logits = model(emb)
            loss = criterion(logits, lab)
            total_loss += loss.item() * emb.size(0)
            total_n += emb.size(0)
            all_logits.append(logits.float().cpu().numpy())
            all_labels.append(lab.cpu().numpy())
    return (total_loss / max(total_n, 1),
            np.concatenate(all_logits, axis=0) if all_logits else np.empty((0,)),
            np.concatenate(all_labels, axis=0) if all_labels else np.empty((0,)))


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    task_cfg = TASK_CONFIG[args.task]
    torch.manual_seed(args.seed); np.random.seed(args.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    out_dir = Path(args.out_dir) / args.task / f"seed_{args.seed}" / (
        f"lr{args.lr:.0e}_d{args.drop_rate}_ls{args.label_smoothing}")
    out_dir.mkdir(parents=True, exist_ok=True)

    if (out_dir / "metrics.json").exists():
        print(f"  metrics.json already exists at {out_dir} -- skipping.")
        return

    print("=" * 70)
    print(f"  {args.model_name} cached-head sweep -- {args.task} | seed={args.seed}")
    print(f"  HP: lr={args.lr}  wd={args.weight_decay}  drop={args.drop_rate}  "
          f"ls={args.label_smoothing}")
    print(f"  Output: {out_dir}")
    print("=" * 70)

    # ── Load embeddings ───────────────────────────────────────────────────────
    ck = torch.load(args.embeddings_pt, map_location="cpu", weights_only=False)
    embeddings = ck["embeddings"]                     # (N, D) float32
    ids        = ck["ids"]                            # list[(sub, ses)]
    if embeddings.shape[1] != args.embed_dim:
        raise RuntimeError(
            f"--embed_dim {args.embed_dim} mismatches loaded embeddings "
            f"{embeddings.shape}. Did you pass the right --embeddings_pt for "
            f"--model_name {args.model_name}?")
    id_to_row = {(s, e): i for i, (s, e) in enumerate(ids)}
    print(f"  Cached embeddings: {embeddings.shape}  ({len(ids)} scans)")

    # ── Splits ────────────────────────────────────────────────────────────────
    seed_dir = Path(args.data_dir) / f"seed_{args.seed}"
    matched_df = load_matched_labels(args.matched_labels_csv)

    common = dict(matched_df=matched_df, task_cfg=task_cfg, id_to_row=id_to_row)
    if args.long_mode is None:
        sf = args.session
        train_df, n_miss_tr = load_split_cached(seed_dir / "train.csv", session_filter=sf, **common)
        val_df,   n_miss_va = load_split_cached(seed_dir / "val.csv",   session_filter=sf, **common)
        test_df,  n_miss_te = load_split_cached(seed_dir / "test.csv",  session_filter=sf, **common)
    else:
        mm = args.max_months
        train_df, n_miss_tr = load_split_cached(seed_dir / "train.csv", max_months=mm, **common)
        val_df,   n_miss_va = load_split_cached(seed_dir / "val.csv",   max_months=mm, **common)
        test_df,  n_miss_te = load_split_cached(seed_dir / "test.csv",  max_months=mm, **common)

    print(f"  Splits -- train: {len(train_df)}  val: {len(val_df)}  test: {len(test_df)}")
    print(f"  Missing embedding -- train: {n_miss_tr}  val: {n_miss_va}  test: {n_miss_te}")
    manifest = pd.concat([train_df.assign(split="train"),
                          val_df.assign(split="val"),
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
    train_ds = EmbeddingDataset(train_df, embeddings, id_to_row)
    val_ds   = EmbeddingDataset(val_df,   embeddings, id_to_row)
    test_ds  = EmbeddingDataset(test_df,  embeddings, id_to_row)
    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True,
                              num_workers=0, drop_last=False)
    val_loader   = DataLoader(val_ds,   batch_size=args.batch_size, shuffle=False)
    test_loader  = DataLoader(test_ds,  batch_size=args.batch_size, shuffle=False)

    # ── Model / opt / scheduler / loss ────────────────────────────────────────
    model = HeadOnly(args.embed_dim, num_labels, drop_rate=args.drop_rate).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr,
                                 weight_decay=args.weight_decay)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="max",
        factor=args.scheduler_factor, patience=args.scheduler_patience)
    criterion = nn.CrossEntropyLoss(weight=class_weights,
                                    label_smoothing=args.label_smoothing)

    wb = init_wandb(args, task_cfg, extra={
        "n_train": len(train_df), "n_val": len(val_df), "n_test": len(test_df),
        "n_params": sum(p.numel() for p in model.parameters())})

    # ── Resume ────────────────────────────────────────────────────────────────
    ckpt_path = out_dir / "last_checkpoint.pt"
    best_metric = -1.0
    best_val_loss_at_best = float("inf")
    best_epoch = -1
    epochs_since_improve = 0
    log_rows = []
    start_epoch = 0
    if (not args.no_resume) and ckpt_path.exists():
        ck2 = torch.load(ckpt_path, map_location=device)
        model.load_state_dict(ck2["model"])
        optimizer.load_state_dict(ck2["optimizer"])
        if "scheduler" in ck2:
            scheduler.load_state_dict(ck2["scheduler"])
        best_metric = ck2.get("best_metric", -1.0)
        best_val_loss_at_best = ck2.get("best_val_loss_at_best", float("inf"))
        best_epoch = ck2.get("best_epoch", -1)
        epochs_since_improve = ck2.get("epochs_since_improve", 0)
        log_rows = ck2.get("log_rows", [])
        start_epoch = ck2["epoch"] + 1
        print(f"  Resuming from epoch {start_epoch + 1}")

    def _save_checkpoint(ep):
        tmp = ckpt_path.with_suffix(".pt.tmp")
        torch.save({"model": model.state_dict(),
                    "optimizer": optimizer.state_dict(),
                    "scheduler": scheduler.state_dict(),
                    "epoch": ep, "best_metric": best_metric,
                    "best_val_loss_at_best": best_val_loss_at_best,
                    "best_epoch": best_epoch,
                    "epochs_since_improve": epochs_since_improve,
                    "log_rows": log_rows}, tmp)
        os.replace(tmp, ckpt_path)

    # ── Train ─────────────────────────────────────────────────────────────────
    try:
        for epoch in range(start_epoch, args.epochs):
            cur_lr = optimizer.param_groups[0]["lr"]
            model.train()
            total_loss, total_correct, total_n = 0.0, 0, 0
            for emb, lab in train_loader:
                emb = emb.to(device, non_blocking=True).float()
                lab = lab.to(device, non_blocking=True).long()
                logits = model(emb)
                loss = criterion(logits, lab)
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                optimizer.step()
                total_loss += loss.item() * emb.size(0)
                total_correct += (logits.argmax(-1) == lab).sum().item()
                total_n += emb.size(0)
            tr_loss = total_loss / max(total_n, 1)
            tr_acc  = total_correct / max(total_n, 1)

            va_loss, va_logits, va_labels = evaluate(model, val_loader, criterion, device)
            if not (np.isfinite(tr_loss) and np.isfinite(va_loss)):
                print(f"  [ABORT] non-finite loss at epoch {epoch+1} -- stopping.")
                break
            va_bacc = float(balanced_accuracy_score(
                va_labels.astype(int), va_logits.argmax(axis=-1)))
            scheduler.step(va_bacc)

            log_rows.append({"epoch": epoch + 1, "lr": cur_lr,
                             "train_loss": tr_loss, "train_acc": tr_acc,
                             "val_loss": va_loss, "val_bacc": va_bacc})
            print(f"  [epoch {epoch+1:>3}/{args.epochs}] lr={cur_lr:.2e}  "
                  f"tr_loss={tr_loss:.4f}  va_loss={va_loss:.4f}  va_bacc={va_bacc:.4f}")

            improved = (va_bacc > best_metric + 1e-6) or (
                abs(va_bacc - best_metric) <= 1e-6 and va_loss < best_val_loss_at_best)
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
                           "val_bacc": va_bacc, "best_val_bacc": best_metric})

            if epochs_since_improve >= args.patience:
                print(f"  Early stopping at epoch {epoch+1}")
                break

        pd.DataFrame(log_rows).to_csv(out_dir / "train_log.csv", index=False)

        # ── Test ──────────────────────────────────────────────────────────────
        if not (out_dir / "best_model.pt").exists():
            print("  [ERROR] no best_model.pt -- no metrics.json written.")
            return
        best_state = torch.load(out_dir / "best_model.pt", map_location=device)
        model.load_state_dict(best_state["net"])
        _, test_logits, test_labels = evaluate(model, test_loader, criterion, device)
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
            "model_id":              f"{args.model_name}_frozen_cached",
            "model_name":            args.model_name,
            "embeddings_pt":         args.embeddings_pt,
            "embed_dim":             args.embed_dim,
            "task":                  args.task,
            "task_description":      task_cfg["description"],
            "session":               cohort_label(args, task_cfg),
            "long_mode":             args.long_mode,
            "max_months":            args.max_months,
            "session_policy":        task_cfg["session_policy"],
            "seed":                  args.seed,
            "strategy":              "frozen_cached",
            "augment":               "none",
            "epochs":                args.epochs,
            "best_epoch":            best_epoch,
            "best_val_balanced_acc": round(float(best_metric), 4),
            "lr":                    args.lr,
            "weight_decay":          args.weight_decay,
            "drop_rate":             args.drop_rate,
            "label_smoothing":       args.label_smoothing,
            "batch_size":            args.batch_size,
            "patience":              args.patience,
            "scheduler_patience":    args.scheduler_patience,
            "scheduler_factor":      args.scheduler_factor,
            "scheduler":             "ReduceLROnPlateau",
            "optimizer":             "Adam",
            "n_train":               len(train_df),
            "n_val":                 len(val_df),
            "n_test":                len(test_df),
            "class_weights":         {i: round(float(w), 4) for i, w in enumerate(cw)},
            "timestamp":             datetime.now().isoformat(),
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
