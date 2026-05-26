r"""
30_train_diff_attention_v2.py   (COLAB+WANDB or local; env: snp)
================================================================
Single-cell trainer for the diff-attention v2 sweep over the
128-SNP recover_all_pool. Driven by WandB sweep params; one process =
one (seq_length, model, functional_concat, aggregation, head, seed)
combination. Reuses BiasedAttention / GlobalAttnPool / ChromHierPool /
MLPHead / MLP3Head from fm_diff_lib.py.

Per-SNP input feature (per user 2026-05-26):
  x[p, i, :] = concat(diff_emb, func_feat (if on), β, dosage, β·dosage)

Attention bias (always): score + γ·|β| + δ·|β·dosage|  (existing BiasedAttention).

Heads:
  - mlp2 / mlp3 : end-to-end AdamW, weighted BCE, early-stop on val_loss
  - rf         : sklearn RandomForestClassifier, class_weight='balanced'
  - xgb        : xgboost.XGBClassifier, scale_pos_weight=neg/pos
  For tree heads the aggregator is a FIXED softmax over `γ·|β| + δ·|β·dosage|`
  with γ=δ=1 (no learnable params).

Sweep target metric: val_balanced_accuracy (logged + posted to WandB run summary).

Usage:
  python snp_pipeline/30_train_diff_attention_v2.py \\
      --base /content/drive/MyDrive/ADNI_SNP \\
      --splits-root /content/drive/MyDrive/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified_post_exclusion/tabular/baseline \\
      --seq-length 1001bp --model bmfm_ref \\
      --functional-concat on --aggregation global_attn \\
      --head mlp2 --seed 0 --wandb-project adni-snp-diff-attn-v2-mlp2
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import (balanced_accuracy_score, confusion_matrix,
                              f1_score, precision_score, recall_score,
                              roc_auc_score, log_loss)
from torch.utils.data import DataLoader, TensorDataset

# Import fm_diff_lib (numeric-leading filename safe)
_FM = Path(__file__).parent / "fm_diff_lib.py"
_spec = importlib.util.spec_from_file_location("fm_diff_lib", _FM)
fl = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(fl)

# Import 29b's loader
_FF = Path(__file__).parent / "29b_load_functional_features.py"
_spec_ff = importlib.util.spec_from_file_location("ff_loader", _FF)
ff = importlib.util.module_from_spec(_spec_ff); _spec_ff.loader.exec_module(ff)

POS_GROUPS = ("AD_bl", "pMCI", "CN_to_AD")
NEG_GROUPS = ("sCN",)

DEFAULT_BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
DEFAULT_SPLITS = Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                       "no_cdr_stratified_post_exclusion/tabular/baseline")


def _load_diff(base: Path, seq_length: str, model: str) -> tuple[np.ndarray, list[str]]:
    """Load (S, H) diff_emb + rsIDs for the given (seq_length, model)."""
    parent_name = ("fm_embeddings_short_seq_1kb_with_alphagenome" if seq_length == "1001bp"
                    else "fm_embeddings_short_seq_101bp")
    # Resolve model dir (caduceus_ph_d256 alias)
    for alias in [model, f"{model}_d256"]:
        d = base / parent_name / alias
        if d.exists():
            p = d / "recover_all_pool_snp_diff.npz"
            if p.exists():
                z = np.load(p, allow_pickle=True)
                return z["diff_emb"].astype(np.float32), [str(x) for x in z["rsids"]]
    raise FileNotFoundError(f"recover_all_pool_snp_diff.npz not found for "
                              f"{seq_length} × {model} under {base}")


def _load_dosage(base: Path) -> tuple[pd.DataFrame, list[str]]:
    """Returns (df with PTID column + per-rsID dosages, rsID col list)."""
    for rel in ("GWAS_comprehensive_v2/recover_all_pool/recover_all_pool_dosage.tsv",
                "tables/recover_all_pool_dosage.tsv",
                "recover_all_pool_dosage.tsv"):
        p = base / rel
        if p.exists():
            df = pd.read_csv(p, sep="\t")
            cols = [c for c in df.columns if c != "PTID"]
            return df, cols
    raise FileNotFoundError(f"recover_all_pool_dosage.tsv not found under {base}")


def _load_beta(base: Path) -> pd.DataFrame:
    for rel in ("GWAS_comprehensive_v2/recover_all_pool/recover_all_pool_beta_A1.tsv",
                "tables/recover_all_pool_beta_A1.tsv",
                "recover_all_pool_beta_A1.tsv"):
        p = base / rel
        if p.exists():
            return pd.read_csv(p, sep="\t")
    raise FileNotFoundError(f"recover_all_pool_beta_A1.tsv not found under {base}")


def _load_labels_for_split(splits_root: Path, seed: int, split: str,
                              base: Path | None = None) -> pd.DataFrame:
    """Returns DataFrame with PTID + y (POS=1 / NEG=0); rows of other groups dropped."""
    p = splits_root / f"seed_{seed}/{split}.csv"
    if not p.exists() and base is not None:
        fallback = base / "splits" / f"seed_{seed}/{split}.csv"
        if fallback.exists():
            p = fallback
    df = pd.read_csv(p, dtype=str)
    df["y_pos"] = ((df["AD_bl"].astype(int) == 1) |
                    (df["pMCI"].astype(int) == 1) |
                    (df["CN_to_AD"].astype(int) == 1)).astype(int)
    df["y_neg"] = (df["sCN"].astype(int) == 1).astype(int)
    df = df[(df["y_pos"] == 1) | (df["y_neg"] == 1)].copy()
    df["y"] = df["y_pos"].astype(int)
    df = df.rename(columns={"Patient_ID": "PTID"})[["PTID", "y"]]
    return df


def _build_x_tensor(diff: np.ndarray, func: np.ndarray | None,
                     beta: np.ndarray, dosage_p: np.ndarray) -> np.ndarray:
    """For one patient, build (S, F) feature matrix.
       Concat: diff (S,H) ‖ func (S,F_func) ‖ β (S,1) ‖ dosage (S,1) ‖ β·dosage (S,1)."""
    S, H = diff.shape
    blocks = [diff]
    if func is not None:
        blocks.append(func)
    beta_col = beta.reshape(S, 1)
    dosage_col = dosage_p.reshape(S, 1)
    blocks.append(beta_col)
    blocks.append(dosage_col)
    blocks.append(beta_col * dosage_col)
    return np.concatenate(blocks, axis=1).astype(np.float32)


def _per_patient_features(diff: np.ndarray, func: np.ndarray | None,
                           beta: np.ndarray, dosage: np.ndarray) -> np.ndarray:
    """Build (P, S, F) tensor."""
    P, S = dosage.shape
    return np.stack([_build_x_tensor(diff, func, beta, dosage[p]) for p in range(P)], axis=0)


def _fixed_aggregator(x: np.ndarray, beta: np.ndarray, dosage: np.ndarray) -> np.ndarray:
    """For tree heads: deterministic softmax over γ·|β| + δ·|β·dosage| (γ=δ=1)."""
    P, S, F = x.shape
    abs_beta = np.abs(beta)                          # (S,)
    abs_bd = np.abs(beta[None, :] * dosage)          # (P, S)
    score = abs_beta[None, :] + abs_bd               # (P, S)
    # Softmax over SNP axis per patient
    s = score - score.max(axis=1, keepdims=True)
    attn = np.exp(s)
    attn /= attn.sum(axis=1, keepdims=True)
    return (attn[:, :, None] * x).sum(axis=1)        # (P, F)


# ───────────────────── Tree-head MLP-emulating model ────────────────────

class TreePatientEmb(nn.Module):
    """Dummy for symmetry — never used for tree heads (numpy path)."""


# ───────────────────── MLP-end-to-end model ─────────────────────────────

class DiffAttnV2(nn.Module):
    """Per-SNP features (B, S, F) → biased attention pool → MLP head."""
    def __init__(self, in_dim: int, n_chroms: int, aggregation: str,
                 head: str, attn_hidden: int = 128):
        super().__init__()
        self.aggregation = aggregation
        if aggregation == "global_attn":
            self.pool = fl.GlobalAttnPool(in_dim, use_delta=True, attn_hidden=attn_hidden)
        else:
            self.pool = fl.ChromHierPool(in_dim, n_chroms=n_chroms, use_delta=True,
                                          attn_hidden=attn_hidden)
        self.norm = nn.LayerNorm(in_dim)
        if head == "mlp2":
            self.head = fl.MLPHead(in_dim)
        elif head == "mlp3":
            self.head = fl.MLP3Head(in_dim)
        else:
            raise ValueError(f"Unsupported head for end-to-end: {head}")

    def forward(self, x, beta, dosage_x_beta, chrom_idx=None):
        # x: (B, S, F); beta: (S,); dosage_x_beta: (B, S); chrom_idx: (S,)
        B, S, _ = x.shape
        mask = torch.ones(B, S, dtype=torch.bool, device=x.device)
        abs_beta = beta.abs().unsqueeze(0).expand(B, S)
        abs_dosbeta = dosage_x_beta.abs()
        if chrom_idx is None:
            chrom_idx = torch.zeros(S, dtype=torch.long, device=x.device)
        emb, _ = self.pool(x, mask, chrom_idx, abs_beta, abs_dosbeta)
        return self.head(self.norm(emb)).squeeze(-1)


# ───────────────────── Train / eval loops ───────────────────────────────

def _metrics(y_true, prob, pred, loss_val=None) -> dict:
    cm = confusion_matrix(y_true, pred, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel() if cm.size == 4 else (0, 0, 0, 0)
    sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    m = {
        "balanced_accuracy": balanced_accuracy_score(y_true, pred),
        "roc_auc": roc_auc_score(y_true, prob) if len(set(y_true)) > 1 else 0.5,
        "f1": f1_score(y_true, pred, zero_division=0),
        "precision": precision_score(y_true, pred, zero_division=0),
        "recall": sens,
        "sensitivity": sens,
        "specificity": spec,
        "tp": int(tp), "tn": int(tn), "fp": int(fp), "fn": int(fn),
    }
    if loss_val is not None:
        m["loss"] = float(loss_val)
    return m


def _train_mlp(model: nn.Module, train_data, val_data, *,
               epochs: int, batch_size: int, lr: float, weight_decay: float,
               patience: int, device: str, pos_weight: float,
               wandb_run=None) -> dict:
    """End-to-end SGD with early-stop on val_loss; per-epoch wandb log."""
    model = model.to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    Xtr, beta_tr, dxb_tr, chrom_tr, ytr = train_data
    Xva, beta_va, dxb_va, chrom_va, yva = val_data
    pos_w = torch.tensor([pos_weight], device=device)
    loss_fn = nn.BCEWithLogitsLoss(pos_weight=pos_w)

    n_train = Xtr.shape[0]
    best_val_loss = float("inf"); best_state = None; best_epoch = 0; bad = 0
    history = []
    for ep in range(epochs):
        model.train()
        idx = torch.randperm(n_train)
        train_loss_sum = 0.0; train_loss_n = 0
        # Also track train balAcc on the same epoch for a fairness curve
        tr_probs_chunks = []; tr_y_chunks = []
        for s in range(0, n_train, batch_size):
            b = idx[s:s+batch_size]
            opt.zero_grad()
            logits = model(Xtr[b].to(device), beta_tr.to(device),
                            dxb_tr[b].to(device),
                            chrom_tr.to(device) if chrom_tr is not None else None)
            loss = loss_fn(logits, ytr[b].float().to(device))
            loss.backward()
            opt.step()
            train_loss_sum += loss.item() * len(b)
            train_loss_n += len(b)
            tr_probs_chunks.append(torch.sigmoid(logits.detach()).cpu().numpy())
            tr_y_chunks.append(ytr[b].cpu().numpy())
        train_loss = train_loss_sum / max(1, train_loss_n)
        tr_probs = np.concatenate(tr_probs_chunks)
        tr_y = np.concatenate(tr_y_chunks)
        tr_pred = (tr_probs > 0.5).astype(int)
        tr_balacc = balanced_accuracy_score(tr_y, tr_pred)

        # Val pass
        model.eval()
        with torch.no_grad():
            v_logits = model(Xva.to(device), beta_va.to(device),
                              dxb_va.to(device),
                              chrom_va.to(device) if chrom_va is not None else None)
            v_loss = loss_fn(v_logits, yva.float().to(device)).item()
            v_prob = torch.sigmoid(v_logits).cpu().numpy()
            v_pred = (v_prob > 0.5).astype(int)
            v_m = _metrics(yva.cpu().numpy(), v_prob, v_pred, v_loss)
        history.append({"epoch": ep, "train_loss": train_loss,
                          "train_balanced_accuracy": tr_balacc,
                          **{f"val_{k}": v for k, v in v_m.items()}})

        # Per-epoch wandb log → learning curves visible in WandB UI
        if wandb_run is not None:
            wandb_run.log({
                "epoch": ep,
                "train_loss": train_loss,
                "train_balanced_accuracy": tr_balacc,
                "val_loss": v_loss,
                "val_balanced_accuracy": v_m["balanced_accuracy"],
                "val_roc_auc": v_m["roc_auc"],
            }, step=ep)

        if v_loss < best_val_loss - 1e-5:
            best_val_loss = v_loss; bad = 0
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            best_epoch = ep
        else:
            bad += 1
            if bad >= patience:
                break
    if best_state:
        model.load_state_dict(best_state)
    return {"epochs_run": ep + 1, "best_epoch": best_epoch, "history": history}


def _eval_mlp(model, data, device) -> tuple[dict, np.ndarray]:
    Xe, beta_e, dxb_e, chrom_e, y_e = data
    model.eval()
    with torch.no_grad():
        logits = model(Xe.to(device), beta_e.to(device), dxb_e.to(device),
                        chrom_e.to(device) if chrom_e is not None else None)
        prob = torch.sigmoid(logits).cpu().numpy()
        pred = (prob > 0.5).astype(int)
    y = y_e.cpu().numpy()
    pos_w = max(1.0, (y == 0).sum() / max(1, (y == 1).sum()))
    loss = log_loss(y, np.clip(prob, 1e-7, 1 - 1e-7),
                     sample_weight=np.where(y == 1, pos_w, 1.0))
    return _metrics(y, prob, pred, loss_val=loss), prob


def _fit_tree(head: str, train_emb, train_y, val_emb, val_y, seed: int) -> dict:
    pos = int((train_y == 1).sum()); neg = int((train_y == 0).sum())
    pos_w = neg / max(1, pos)
    if head == "rf":
        from sklearn.ensemble import RandomForestClassifier
        clf = RandomForestClassifier(n_estimators=500, n_jobs=-1, random_state=seed,
                                       class_weight="balanced")
        clf.fit(train_emb, train_y)
    elif head == "xgb":
        from xgboost import XGBClassifier
        clf = XGBClassifier(n_estimators=500, max_depth=4, learning_rate=0.05,
                             tree_method="hist", scale_pos_weight=pos_w,
                             eval_metric="logloss", random_state=seed)
        clf.fit(train_emb, train_y,
                eval_set=[(val_emb, val_y)], verbose=False)
    else:
        raise ValueError(head)
    prob = clf.predict_proba(val_emb)[:, 1]
    pred = (prob > 0.5).astype(int)
    loss = log_loss(val_y, np.clip(prob, 1e-7, 1-1e-7),
                     sample_weight=np.where(val_y == 1, pos_w, 1.0))
    return {"val_metrics": _metrics(val_y, prob, pred, loss),
            "val_prob": prob, "val_pred": pred, "classifier": clf}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", type=Path, default=DEFAULT_BASE)
    ap.add_argument("--splits-root", "--splits_root", type=Path, default=DEFAULT_SPLITS,
                    dest="splits_root")
    ap.add_argument("--seq-length", "--seq_length", choices=["1001bp", "101bp"],
                    required=True, dest="seq_length")
    ap.add_argument("--model", choices=["bmfm_ref", "bmfm_snp", "ntv2",
                                          "caduceus_ph", "caduceus_ps"], required=True)
    ap.add_argument("--functional-concat", "--functional_concat",
                    default="off", dest="functional_concat",
                    help="off/on (or True/False — WandB sweep may cast).")
    ap.add_argument("--aggregation", choices=["global_attn", "chrom_hier"], default="global_attn")
    ap.add_argument("--head", choices=["mlp2", "mlp3", "rf", "xgb"], required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--epochs", type=int, default=80)
    ap.add_argument("--batch-size", "--batch_size", type=int, default=32, dest="batch_size")
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--weight-decay", "--weight_decay", type=float, default=1e-4,
                    dest="weight_decay")
    ap.add_argument("--patience", type=int, default=15)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--wandb-project", "--wandb_project", default=None,
                    help="If set, init WandB with this project name",
                    dest="wandb_project")
    ap.add_argument("--wandb-entity", "--wandb_entity", default=None, dest="wandb_entity")
    ap.add_argument("--output-root", "--output_root", type=Path,
                    default=Path("outputs/diff_attn_v2"), dest="output_root")
    args = ap.parse_args()
    # Normalize functional_concat to {'off', 'on'} (WandB may send True/False/"True")
    fc_str = str(args.functional_concat).strip().lower()
    if fc_str in ("true", "on", "1", "yes"):
        args.functional_concat = "on"
    elif fc_str in ("false", "off", "0", "no", ""):
        args.functional_concat = "off"
    else:
        raise ValueError(f"Bad --functional-concat: {args.functional_concat!r}")
    np.random.seed(args.seed); torch.manual_seed(args.seed)

    # ── WandB init ─────────────────────────────────────────────────────
    wb = None
    if args.wandb_project:
        import wandb
        wb = wandb.init(project=args.wandb_project, entity=args.wandb_entity,
                         config=vars(args),
                         name=(f"{args.seq_length}_{args.model}_func{args.functional_concat}"
                                f"_{args.aggregation}_{args.head}_s{args.seed}"))

    # ── Load data ──────────────────────────────────────────────────────
    print(f"[load] diff: {args.seq_length} × {args.model}")
    diff, diff_rsids = _load_diff(args.base, args.seq_length, args.model)
    print(f"  diff shape: {diff.shape}")

    print(f"[load] dosage matrix")
    dos_df, dos_rsids = _load_dosage(args.base)
    print(f"  dosage shape: {dos_df.shape}")

    print(f"[load] beta_A1 table")
    beta_df = _load_beta(args.base)

    # Align rsIDs: diff_rsids should match dos_rsids order (same recover_all_pool order)
    target_rsids = diff_rsids
    # Reorder dosage cols
    dos_arr = dos_df[target_rsids].values.astype(np.float32)
    # Mean-impute NA per train fold later
    # Beta in target order
    beta_df_idx = beta_df.set_index("rsID")
    beta = np.array([float(beta_df_idx.loc[rs, "beta_A1"]) if rs in beta_df_idx.index else 0.0
                     for rs in target_rsids], dtype=np.float32)
    # CHR index
    chrom_arr = np.array([int(beta_df_idx.loc[rs, "CHR"]) if rs in beta_df_idx.index else 0
                          for rs in target_rsids], dtype=np.int64)
    chrom_uniq = sorted(set(chrom_arr.tolist()))
    chrom_idx_map = {c: i for i, c in enumerate(chrom_uniq)}
    chrom_idx = np.array([chrom_idx_map[c] for c in chrom_arr], dtype=np.int64)
    n_chroms = len(chrom_uniq)
    print(f"  n_chroms: {n_chroms}")

    # ── Optional functional features ───────────────────────────────────
    func = None
    if args.functional_concat == "on":
        parent_func = args.base / "fm_embeddings_short_seq_1kb_with_alphagenome"
        func, _ = ff.load_functional_features(parent_func, target_rsids)
        print(f"  functional features: {func.shape}")

    # ── Labels per split ───────────────────────────────────────────────
    print(f"[labels] splits_root={args.splits_root}, seed={args.seed}")
    labels = {sp: _load_labels_for_split(args.splits_root, args.seed, sp, base=args.base)
               for sp in ("train", "val", "test")}
    for sp, df in labels.items():
        print(f"  {sp}: n={len(df)}, pos={(df.y==1).sum()}, neg={(df.y==0).sum()}")

    # Build per-split (X, y) tensors. Align dosage_df by PTID.
    dos_by_ptid = dos_df.set_index("PTID")
    splits = {}
    for sp, df in labels.items():
        ptids = df["PTID"].tolist()
        # Get dosage rows
        present_ptids = [p for p in ptids if p in dos_by_ptid.index]
        missing = set(ptids) - set(present_ptids)
        if missing:
            print(f"  [warn] {sp}: {len(missing)} patients have no dosage row")
        df_ok = df[df["PTID"].isin(present_ptids)].reset_index(drop=True)
        dos_split = dos_by_ptid.loc[df_ok["PTID"], target_rsids].values.astype(np.float32)
        splits[sp] = {"ptid": df_ok["PTID"].tolist(),
                       "y": df_ok["y"].values.astype(np.int64),
                       "dosage": dos_split}

    # Mean-impute dosage NA using TRAIN means only
    train_means = np.nanmean(splits["train"]["dosage"], axis=0)
    for sp in splits:
        d = splits[sp]["dosage"]
        nan_mask = np.isnan(d)
        for i in range(d.shape[1]):
            d[nan_mask[:, i], i] = train_means[i]
        splits[sp]["dosage"] = d

    # Build per-(p, i, F) feature tensors
    def _build(sp):
        d = splits[sp]
        X = _per_patient_features(diff, func, beta, d["dosage"])  # (P, S, F)
        beta_t = beta.copy()                                      # (S,)
        dxb = beta[None, :] * d["dosage"]                         # (P, S) for use_delta
        return X, beta_t, dxb, d["y"]

    Xtr, beta_tr, dxb_tr, ytr = _build("train")
    Xva, beta_va, dxb_va, yva = _build("val")
    Xte, beta_te, dxb_te, yte = _build("test")
    F = Xtr.shape[2]
    print(f"[feat] in_dim per SNP F={F}; train n={Xtr.shape[0]}, val n={Xva.shape[0]}")

    # ── Train + eval ───────────────────────────────────────────────────
    out_dir = args.output_root / args.head / (f"{args.seq_length}_{args.model}_"
                                                f"func{args.functional_concat}_"
                                                f"{args.aggregation}_s{args.seed}")
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.head in ("mlp2", "mlp3"):
        # Build tensors for MLP path
        Xtr_t = torch.from_numpy(Xtr)
        Xva_t = torch.from_numpy(Xva)
        Xte_t = torch.from_numpy(Xte)
        beta_t = torch.from_numpy(beta_tr)
        dxb_tr_t = torch.from_numpy(dxb_tr.astype(np.float32))
        dxb_va_t = torch.from_numpy(dxb_va.astype(np.float32))
        dxb_te_t = torch.from_numpy(dxb_te.astype(np.float32))
        chrom_t = torch.from_numpy(chrom_idx) if args.aggregation == "chrom_hier" else None
        ytr_t = torch.from_numpy(ytr)
        yva_t = torch.from_numpy(yva)
        yte_t = torch.from_numpy(yte)

        pos = int((ytr == 1).sum()); neg = int((ytr == 0).sum())
        pos_weight = neg / max(1, pos)
        print(f"  pos_weight (neg/pos): {pos_weight:.3f}")

        model = DiffAttnV2(in_dim=F, n_chroms=n_chroms,
                            aggregation=args.aggregation, head=args.head)
        fit_info = _train_mlp(model,
                                (Xtr_t, beta_t, dxb_tr_t, chrom_t, ytr_t),
                                (Xva_t, beta_t, dxb_va_t, chrom_t, yva_t),
                                epochs=args.epochs, batch_size=args.batch_size,
                                lr=args.lr, weight_decay=args.weight_decay,
                                patience=args.patience, device=args.device,
                                pos_weight=pos_weight, wandb_run=wb)
        val_m, val_prob = _eval_mlp(model, (Xva_t, beta_t, dxb_va_t, chrom_t, yva_t), args.device)
        test_m, test_prob = _eval_mlp(model, (Xte_t, beta_t, dxb_te_t, chrom_t, yte_t), args.device)
        torch.save(model.state_dict(), out_dir / "model.pt")
    else:
        # Tree path: fixed aggregator → patient embedding → sklearn/xgboost
        emb_tr = _fixed_aggregator(Xtr, beta_tr, splits["train"]["dosage"])
        emb_va = _fixed_aggregator(Xva, beta_va, splits["val"]["dosage"])
        emb_te = _fixed_aggregator(Xte, beta_te, splits["test"]["dosage"])
        print(f"  tree emb dims: train {emb_tr.shape}, val {emb_va.shape}")
        result = _fit_tree(args.head, emb_tr, ytr, emb_va, yva, args.seed)
        val_m = result["val_metrics"]; val_prob = result["val_prob"]
        # Test eval
        clf = result["classifier"]
        test_prob = clf.predict_proba(emb_te)[:, 1]
        test_pred = (test_prob > 0.5).astype(int)
        pos_w = (yte == 0).sum() / max(1, (yte == 1).sum())
        test_loss = log_loss(yte, np.clip(test_prob, 1e-7, 1-1e-7),
                              sample_weight=np.where(yte == 1, pos_w, 1.0))
        test_m = _metrics(yte, test_prob, test_pred, test_loss)
        fit_info = {"head": args.head}

    # ── Log to WandB + write metrics ───────────────────────────────────
    metrics = {"val": val_m, "test": test_m, "config": vars(args)}
    (out_dir / "metrics.json").write_text(json.dumps(metrics, default=str, indent=2))

    log_payload = {}
    for k, v in val_m.items():
        log_payload[f"val_{k}"] = v
    for k, v in test_m.items():
        log_payload[f"test_{k}"] = v
    print(f"\n=== val balAcc={val_m['balanced_accuracy']:.4f}  ROC-AUC={val_m['roc_auc']:.4f}  loss={val_m.get('loss', float('nan')):.4f} ===")
    if wb:
        wb.log(log_payload)
        # Set summary metrics so the sweep ranks by these
        for k, v in log_payload.items():
            wb.summary[k] = v
        wb.finish()


if __name__ == "__main__":
    main()
