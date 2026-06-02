r"""
30v3_train_diff_attention_func.py   (COLAB+WANDB; env: snp/colab)
================================================================
v3 of the diff-attention sweep trainer — adds the
`--func-integration-mode` axis (7 values) to compare ways of injecting
AlphaGenome functional features into the SNP-attention pool / per-SNP
value vector. Everything else (cohort, splits, β harmonisation, dosage
imputation, weighted BCE, val_loss early-stop, per-epoch wandb logging)
is verbatim from v2 / `30_train_diff_attention_v2.py`.

Modes (passed via --func-integration-mode):
  - none                                : no func; value = [β, dosage, β·dosage]
  - per_snp_zscored_23                  : 23-D raw AG features (z-scored) concat
  - per_snp_summaries_8                 : 8 modality summaries (4 abs + 4 signed)
                                          z-scored concat
  - attn_bias_single                    : value = baseline 3-D; bias += ε · func_imp_abs
  - attn_bias_per_modality              : value = baseline 3-D; bias += Σ_m ε_m · m_abs
  - value_mult                          : value = 11-D multiplicative block
  - attn_bias_per_modality_plus_value_mult : combine the two

For tree heads (`rf`/`xgb`) all ε's are FIXED at 1; the `_fixed_aggregator`
is extended to honour the same bias formulation.

Z-scoring: all func columns z-scored using train-fold SNP-row stats only
(μ, σ saved per-run for reproducibility).

Splice signing: bias uses splice_abs (magnitude only — never reduces
attention); value-side multiplicative interactions use splice_signed.
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
                              f1_score, log_loss, precision_score, recall_score,
                              roc_auc_score)

# Import fm_diff_lib + functional feature loader by file path
_HERE = Path(__file__).parent
_FM = _HERE / "fm_diff_lib.py"
_spec = importlib.util.spec_from_file_location("fm_diff_lib", _FM)
fl = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(fl)

_FF = _HERE / "29b_load_functional_features.py"
_spec_ff = importlib.util.spec_from_file_location("ff_loader", _FF)
ff = importlib.util.module_from_spec(_spec_ff); _spec_ff.loader.exec_module(ff)

# ── Constants ─────────────────────────────────────────────────────────────

POS_GROUPS = ("AD_bl", "pMCI", "CN_to_AD")
NEG_GROUPS = ("sCN",)

DEFAULT_BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
DEFAULT_SPLITS = Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                       "no_cdr_stratified_post_exclusion/tabular/baseline")

ALL_MODES = [
    "none",
    "per_snp_zscored_23",
    "per_snp_summaries_8",
    "attn_bias_single",
    "attn_bias_per_modality",
    "value_mult",
    "attn_bias_per_modality_plus_value_mult",
]

MODALITIES = ("eqtl", "acc", "tss", "splice")  # canonical order


# ── Data loading (paths kept v2-compatible; modes resolved later) ─────────

def _load_diff(base: Path, seq_length: str, model: str) -> tuple[np.ndarray, list[str]]:
    parent_name = ("fm_embeddings_short_seq_1kb_with_alphagenome" if seq_length == "1001bp"
                    else "fm_embeddings_short_seq_101bp")
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
                              base: Path | None = None,
                              exclude_cn_to_ad: bool = False) -> pd.DataFrame:
    p = splits_root / f"seed_{seed}/{split}.csv"
    if not p.exists() and base is not None:
        fallback = base / "splits" / f"seed_{seed}/{split}.csv"
        if fallback.exists():
            p = fallback
    df = pd.read_csv(p, dtype=str)
    pos_terms = [df["AD_bl"].astype(int) == 1, df["pMCI"].astype(int) == 1]
    if not exclude_cn_to_ad:
        pos_terms.append(df["CN_to_AD"].astype(int) == 1)
    y_pos = pos_terms[0]
    for t in pos_terms[1:]:
        y_pos = y_pos | t
    df["y_pos"] = y_pos.astype(int)
    df["y_neg"] = (df["sCN"].astype(int) == 1).astype(int)
    df = df[(df["y_pos"] == 1) | (df["y_neg"] == 1)].copy()
    df["y"] = df["y_pos"].astype(int)
    df = df.rename(columns={"Patient_ID": "PTID"})[["PTID", "y"]]
    return df


# ── Per-task label loaders (added 2026-06-02 for task-conditioned re-extract) ─
#
# Each loader takes (splits_root, seed, split, base) and returns DataFrame
# with columns [PTID, y]. y is int (binary or multiclass class index).
#
# `splits_root` semantics match the original `--splits-root` flag:
#   - For sCN-vs-pCN and T3 horizons: pass `splits_root` pointing at
#     baseline/seed_{s}/ (or `/nonexistent/force-fallback` to fall back to
#     `<base>/splits/`).
#   - For T4 multiclass: pass `splits_root` pointing at baseline_T4/seed_{s}/
#     (different seeds; fall back to `<base>/splits_T4/` if the path doesn't
#     exist locally).

POS_GROUPS_SCN_VS_PCN = ("pCN_to_AD", "pCN_to_MCI")
T3_DROPPED_GROUPS = ("AD_bl", "Excluded")


def _read_split_csv(splits_root: Path, seed: int, split: str,
                      base: Path | None, fallback_subdir: str) -> pd.DataFrame:
    """Locate the split CSV. Tries `splits_root/seed_{s}/{split}.csv` first;
    falls back to `<base>/<fallback_subdir>/seed_{s}/{split}.csv`."""
    p = splits_root / f"seed_{seed}/{split}.csv"
    if not p.exists() and base is not None:
        fb = base / fallback_subdir / f"seed_{seed}/{split}.csv"
        if fb.exists():
            p = fb
    return pd.read_csv(p, dtype=str)


def _load_labels_sCN_vs_pCN(splits_root: Path, seed: int, split: str,
                                base: Path | None = None) -> pd.DataFrame:
    """sCN (neg) vs pCN_to_AD ∪ pCN_to_MCI (pos). Strict baseline-CN cohort."""
    df = _read_split_csv(splits_root, seed, split, base, "splits")
    for col in ("sCN", "pCN_to_AD", "pCN_to_MCI"):
        df[col] = pd.to_numeric(df.get(col, 0), errors="coerce").fillna(0).astype(int)
    y_pos = ((df["pCN_to_AD"] == 1) | (df["pCN_to_MCI"] == 1)).astype(int)
    y_neg = (df["sCN"] == 1).astype(int)
    df = df[(y_pos == 1) | (y_neg == 1)].copy()
    df["y"] = y_pos.loc[df.index].astype(int)
    return df.rename(columns={"Patient_ID": "PTID"})[["PTID", "y"]]


def _load_conv_meta(base: Path) -> pd.DataFrame:
    """Read conversion_labels.tsv. Looks at <base>/conversion_labels.tsv first
    (Colab convention), then <base>/conversion_labels/conversion_labels.tsv
    (local D: layout)."""
    for rel in ("conversion_labels.tsv",
                "conversion_labels/conversion_labels.tsv"):
        p = base / rel
        if p.exists():
            return pd.read_csv(p, sep="\t")
    raise FileNotFoundError(
        f"conversion_labels.tsv not found under {base}. Push the file to "
        f"{base / 'conversion_labels.tsv'} before running T3 tasks.")


def _load_labels_T3_horizon(splits_root: Path, seed: int, split: str,
                                base: Path, N: int) -> pd.DataFrame:
    """T3 binary horizon at threshold N years. Cohort drops AD_bl + Excluded.
    y=1 if `ever_conversion_AD AND years_to_AD <= N`; y=0 if (no AD AND
    FU_years >= N) or (years_to_AD > N); drop if (no AD AND FU_years < N)."""
    df = _read_split_csv(splits_root, seed, split, base, "splits")
    conv = _load_conv_meta(base)
    conv["Patient_ID"] = conv["Patient_ID"].astype(str)
    df["Patient_ID"] = df["Patient_ID"].astype(str)
    df = df.merge(conv[["Patient_ID", "conversion_group", "ever_conversion_AD",
                          "years_to_AD", "FU_years"]],
                    on="Patient_ID", how="left")
    df = df[~df["conversion_group"].isin(T3_DROPPED_GROUPS)].copy()
    ev = pd.to_numeric(df["ever_conversion_AD"], errors="coerce")
    yta = pd.to_numeric(df["years_to_AD"], errors="coerce")
    fu  = pd.to_numeric(df["FU_years"], errors="coerce")
    label = np.full(len(df), -1, dtype=np.int64)   # -1 = drop sentinel
    pos_mask = (ev == 1) & yta.notna() & (yta <= N)
    label[pos_mask.values] = 1
    neg_mask_a = (ev == 0) & fu.notna() & (fu >= N)
    neg_mask_b = (ev == 1) & yta.notna() & (yta > N)
    label[(neg_mask_a | neg_mask_b).values] = 0
    df["y"] = label
    df = df[df["y"] != -1].copy()
    return df.rename(columns={"Patient_ID": "PTID"})[["PTID", "y"]]


def _load_labels_T4_3class(splits_root: Path, seed: int, split: str,
                                base: Path | None = None) -> pd.DataFrame:
    """T4 3-class horizon (0=<3y, 1=3-7y, 2=>=7y). Cohort already restricted
    to converters (pMCI ∪ pCN_to_AD) by baseline_T4 split builder."""
    df = _read_split_csv(splits_root, seed, split, base, "splits_T4")
    df["y"] = pd.to_numeric(df["Label_T4"], errors="coerce").astype("Int64")
    df = df[df["y"].notna()].copy()
    df["y"] = df["y"].astype(int)
    return df.rename(columns={"Patient_ID": "PTID"})[["PTID", "y"]]


TASK_REGISTRY = {
    "ad_vs_cn":         {"num_classes": 1,
                          "loader": lambda sr, s, sp, base, exclude_cn_to_ad=False:
                              _load_labels_for_split(sr, s, sp, base, exclude_cn_to_ad)},
    "sCN_vs_pCN":       {"num_classes": 1,
                          "loader": lambda sr, s, sp, base, **_:
                              _load_labels_sCN_vs_pCN(sr, s, sp, base)},
    "T3_horizon_3y":    {"num_classes": 1,
                          "loader": lambda sr, s, sp, base, **_:
                              _load_labels_T3_horizon(sr, s, sp, base, 3)},
    "T3_horizon_5y":    {"num_classes": 1,
                          "loader": lambda sr, s, sp, base, **_:
                              _load_labels_T3_horizon(sr, s, sp, base, 5)},
    "T3_horizon_7y":    {"num_classes": 1,
                          "loader": lambda sr, s, sp, base, **_:
                              _load_labels_T3_horizon(sr, s, sp, base, 7)},
    "T3_horizon_10y":   {"num_classes": 1,
                          "loader": lambda sr, s, sp, base, **_:
                              _load_labels_T3_horizon(sr, s, sp, base, 10)},
    "T4_3class":        {"num_classes": 3,
                          "loader": lambda sr, s, sp, base, **_:
                              _load_labels_T4_3class(sr, s, sp, base)},
}


def _resolve_ag_parent(base: Path) -> Path:
    """Find the dir containing alphagenome_<modality>/ subdirs and tidy_long.tsv.gz."""
    for rel in ("fm_embeddings_short_seq_1kb_with_alphagenome",
                "fm_embeddings_short_seq_1kb"):
        p = base / rel
        if (p / "alphagenome_eqtl").exists():
            return p
    raise FileNotFoundError(f"alphagenome_eqtl/ not found under {base}")


# ── Z-scoring ─────────────────────────────────────────────────────────────

def _zscore_train(matrix: np.ndarray, train_snp_mask: np.ndarray) -> tuple[np.ndarray, dict]:
    """matrix: (N_SNPs, K)  train_snp_mask: (N_SNPs,) bool
    Returns (z-scored matrix, stats dict with mean/std arrays per column)."""
    if matrix.ndim == 1:
        matrix = matrix.reshape(-1, 1)
    train_rows = matrix[train_snp_mask]
    mu = train_rows.mean(axis=0)
    sigma = train_rows.std(axis=0) + 1e-8
    out = (matrix - mu) / sigma
    return out.astype(np.float32), {"mean": mu.tolist(), "std": sigma.tolist()}


# ── Mode-aware per-SNP feature builder ────────────────────────────────────

def _build_value_block(mode: str, beta: np.ndarray, dosage_p: np.ndarray,
                         func_block: np.ndarray | None,
                         summaries_z: dict | None) -> np.ndarray:
    """For one patient: (S, K) value-feature matrix per mode.
    beta: (S,)   dosage_p: (S,)
    func_block: (S, F_func) z-scored — used for per_snp_* modes
    summaries_z: dict of (S,) signed/abs z-scored — used for value_mult
    """
    S = beta.shape[0]
    beta_col = beta.reshape(S, 1).astype(np.float32)
    dosage_col = dosage_p.reshape(S, 1).astype(np.float32)
    bd_col = (beta_col * dosage_col).astype(np.float32)

    if mode in ("none", "attn_bias_single", "attn_bias_per_modality"):
        return np.concatenate([beta_col, dosage_col, bd_col], axis=1)

    if mode == "per_snp_zscored_23":
        return np.concatenate([beta_col, dosage_col, bd_col, func_block], axis=1)

    if mode == "per_snp_summaries_8":
        cols = [summaries_z[f"{m}_{kind}"].reshape(S, 1)
                 for m in MODALITIES for kind in ("abs", "signed")]
        return np.concatenate([beta_col, dosage_col, bd_col, *cols], axis=1)

    if mode in ("value_mult", "attn_bias_per_modality_plus_value_mult"):
        # 11-D multiplicative block (user spec 2026-05-26):
        # [β, dosage, β·dosage,
        #  dosage·{m}_signed for each m, β·dosage·{m}_signed for each m]
        blocks = [beta_col, dosage_col, bd_col]
        for m in MODALITIES:
            sgn = summaries_z[f"{m}_signed"].reshape(S, 1)
            blocks.append((dosage_col * sgn).astype(np.float32))
        for m in MODALITIES:
            sgn = summaries_z[f"{m}_signed"].reshape(S, 1)
            blocks.append((bd_col * sgn).astype(np.float32))
        return np.concatenate(blocks, axis=1)

    raise ValueError(f"Unknown func_integration_mode: {mode}")


def _per_patient_values(mode: str, beta: np.ndarray, dosage: np.ndarray,
                          func_block: np.ndarray | None,
                          summaries_z: dict | None,
                          diff: np.ndarray) -> np.ndarray:
    """Build (P, S, F_full) tensor = diff_emb || value_block (mode-specific)."""
    P, S = dosage.shape
    rows = []
    for p in range(P):
        v = _build_value_block(mode, beta, dosage[p], func_block, summaries_z)
        x_full = np.concatenate([diff, v], axis=1).astype(np.float32)
        rows.append(x_full)
    return np.stack(rows, axis=0)


# ── Extra-bias tensor builder ─────────────────────────────────────────────

def _build_extra_bias_static(mode: str, summaries_z: dict | None,
                                P: int, S: int) -> np.ndarray | None:
    """For tree-head fixed aggregator: deterministic (S,) bias (all ε=1).
    Returns (S,) array OR None for modes without bias terms."""
    if mode == "attn_bias_single":
        # func_importance = max across 4 abs modalities
        return np.maximum.reduce([summaries_z[f"{m}_abs"] for m in MODALITIES])
    if mode in ("attn_bias_per_modality",
                "attn_bias_per_modality_plus_value_mult"):
        return sum(summaries_z[f"{m}_abs"] for m in MODALITIES)
    return None


# ── DiffAttn v3 model (mlp heads) ─────────────────────────────────────────

class DiffAttnV3(nn.Module):
    """v3 model — uses fm_diff_lib pools but supplies per-modality learnable ε's
    via the new extra_bias kwarg."""

    def __init__(self, in_dim: int, n_chroms: int, aggregation: str,
                 head: str, mode: str, attn_hidden: int = 128,
                 mlp_width: int = 256, mlp_dropout: float = 0.3,
                 out_dim: int = 1):
        super().__init__()
        self.mode = mode
        self.aggregation = aggregation
        self.out_dim = out_dim
        if aggregation == "global_attn":
            self.pool = fl.GlobalAttnPool(in_dim, use_delta=True, attn_hidden=attn_hidden)
        else:
            self.pool = fl.ChromHierPool(in_dim, n_chroms=n_chroms, use_delta=True,
                                          attn_hidden=attn_hidden)
        self.norm = nn.LayerNorm(in_dim)
        # MLP-head width/dropout now configurable for the HP sweep. out_dim=K
        # for multiclass (T4); MLPHead emits (B, K) logits when out_dim>1.
        if head == "mlp2":
            self.head = fl.MLPHead(in_dim, hidden=mlp_width, dropout=mlp_dropout,
                                     out_dim=out_dim)
        elif head == "mlp3":
            self.head = fl.MLP3Head(in_dim, h1=mlp_width, h2=max(mlp_width // 2, 32),
                                     dropout=mlp_dropout, out_dim=out_dim)
        else:
            raise ValueError(f"Unsupported head for end-to-end: {head}")
        # Learnable ε's
        if mode == "attn_bias_single":
            self.eps = nn.Parameter(torch.zeros(1))
            self._eps_shape = "single"
        elif mode in ("attn_bias_per_modality",
                       "attn_bias_per_modality_plus_value_mult"):
            self.eps = nn.Parameter(torch.zeros(len(MODALITIES)))
            self._eps_shape = "per_modality"
        else:
            self.register_parameter("eps", None)
            self._eps_shape = "none"

    def _compute_extra_bias(self, modality_abs: torch.Tensor | None,
                              func_imp_abs: torch.Tensor | None) -> torch.Tensor | None:
        """modality_abs: (S, n_modalities) z-scored abs values per SNP
        func_imp_abs:  (S,) z-scored max-abs across modalities per SNP
        Returns (S,) tensor scaled by learnable ε(s), or None."""
        if self._eps_shape == "single" and func_imp_abs is not None:
            return self.eps * func_imp_abs                       # (S,)
        if self._eps_shape == "per_modality" and modality_abs is not None:
            return (modality_abs * self.eps.view(1, -1)).sum(dim=-1)  # (S,)
        return None

    def forward(self, x, beta, dosage_x_beta, chrom_idx=None,
                  modality_abs=None, func_imp_abs=None):
        B, S, _ = x.shape
        mask = torch.ones(B, S, dtype=torch.bool, device=x.device)
        abs_beta = beta.abs().unsqueeze(0).expand(B, S)
        abs_dosbeta = dosage_x_beta.abs()
        if chrom_idx is None:
            chrom_idx = torch.zeros(S, dtype=torch.long, device=x.device)
        # Extra bias is per-SNP (S,); broadcast to (B, S)
        extra_bias_snp = self._compute_extra_bias(modality_abs, func_imp_abs)
        extra_bias = extra_bias_snp.unsqueeze(0).expand(B, S) if extra_bias_snp is not None else None
        emb, _ = self.pool(x, mask, chrom_idx, abs_beta, abs_dosbeta,
                            extra_bias=extra_bias)
        return self.head(self.norm(emb)).squeeze(-1)


# ── Tree-head fixed aggregator (ε's = 1) ──────────────────────────────────

def _fixed_aggregator_v3(x: np.ndarray, beta: np.ndarray, dosage: np.ndarray,
                          extra_bias_snp: np.ndarray | None) -> np.ndarray:
    """Softmax over γ·|β| + δ·|β·dosage| + (extra_bias if not None), γ=δ=1.
    Extra bias is (S,) — broadcast across patients."""
    P, S, F = x.shape
    abs_beta = np.abs(beta)
    abs_bd = np.abs(beta[None, :] * dosage)
    score = abs_beta[None, :] + abs_bd                        # (P, S)
    if extra_bias_snp is not None:
        score = score + extra_bias_snp[None, :]
    s = score - score.max(axis=1, keepdims=True)
    attn = np.exp(s)
    attn /= attn.sum(axis=1, keepdims=True)
    return (attn[:, :, None] * x).sum(axis=1)                 # (P, F)


# ── Metrics + train + eval (copied from v2; minimal changes) ─────────────

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


def _metrics_multi(y_true, proba, pred, loss_val=None,
                     num_classes: int = 3) -> dict:
    """K-class metrics: macro BalAcc / AUC / F1 / precision / recall + per-class
    confusion-matrix counts (flatten KxK)."""
    labels = list(range(num_classes))
    bacc = balanced_accuracy_score(y_true, pred)
    try:
        auc = roc_auc_score(y_true, proba, multi_class="ovr", labels=labels)
    except Exception:
        auc = 0.5
    f1m   = f1_score(y_true, pred, average="macro", zero_division=0, labels=labels)
    precm = precision_score(y_true, pred, average="macro", zero_division=0,
                             labels=labels)
    recm  = recall_score(y_true, pred, average="macro", zero_division=0,
                          labels=labels)
    cm = confusion_matrix(y_true, pred, labels=labels)
    m = {
        "balanced_accuracy": bacc,    # MACRO BalAcc (sklearn handles this natively)
        "roc_auc":           auc,     # macro-OVR AUC
        "f1":                f1m,
        "precision":         precm,
        "recall":            recm,
    }
    for i in range(num_classes):
        for j in range(num_classes):
            m[f"cm_{i}{j}"] = int(cm[i, j]) if cm.size == num_classes*num_classes else 0
    if loss_val is not None:
        m["loss"] = float(loss_val)
    return m


def _train_mlp(model, train_data, val_data, *,
                epochs: int, batch_size: int, lr: float, weight_decay: float,
                patience: int, device: str, pos_weight: float,
                modality_abs_t, func_imp_abs_t,
                wandb_run=None, num_classes: int = 1,
                class_weights: np.ndarray | None = None) -> dict:
    model = model.to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    Xtr, beta_tr, dxb_tr, chrom_tr, ytr = train_data
    Xva, beta_va, dxb_va, chrom_va, yva = val_data
    # Binary path: BCEWithLogitsLoss with pos_weight (original behaviour).
    # Multiclass path: CrossEntropyLoss with per-class weights.
    is_multiclass = (num_classes > 1)
    if is_multiclass:
        cw = torch.tensor(class_weights, dtype=torch.float32, device=device) \
                if class_weights is not None else None
        loss_fn = nn.CrossEntropyLoss(weight=cw)
    else:
        pos_w = torch.tensor([pos_weight], device=device)
        loss_fn = nn.BCEWithLogitsLoss(pos_weight=pos_w)
    if modality_abs_t is not None: modality_abs_t = modality_abs_t.to(device)
    if func_imp_abs_t is not None: func_imp_abs_t = func_imp_abs_t.to(device)

    n_train = Xtr.shape[0]
    best_val_loss = float("inf"); best_state = None; best_epoch = 0; bad = 0
    history = []
    for ep in range(epochs):
        model.train()
        idx = torch.randperm(n_train)
        train_loss_sum = 0.0; train_loss_n = 0
        tr_probs_chunks, tr_y_chunks = [], []
        for s in range(0, n_train, batch_size):
            b = idx[s:s+batch_size]
            opt.zero_grad()
            logits = model(Xtr[b].to(device), beta_tr.to(device),
                            dxb_tr[b].to(device),
                            chrom_tr.to(device) if chrom_tr is not None else None,
                            modality_abs=modality_abs_t,
                            func_imp_abs=func_imp_abs_t)
            if is_multiclass:
                loss = loss_fn(logits, ytr[b].long().to(device))
            else:
                loss = loss_fn(logits, ytr[b].float().to(device))
            loss.backward()
            opt.step()
            train_loss_sum += loss.item() * len(b)
            train_loss_n += len(b)
            if is_multiclass:
                tr_probs_chunks.append(
                    torch.softmax(logits.detach(), dim=-1).cpu().numpy())
            else:
                tr_probs_chunks.append(torch.sigmoid(logits.detach()).cpu().numpy())
            tr_y_chunks.append(ytr[b].cpu().numpy())
        train_loss = train_loss_sum / max(1, train_loss_n)
        tr_probs = np.concatenate(tr_probs_chunks)
        tr_y = np.concatenate(tr_y_chunks)
        if is_multiclass:
            tr_pred = tr_probs.argmax(axis=1)
            tr_balacc = balanced_accuracy_score(tr_y, tr_pred)
        else:
            tr_balacc = balanced_accuracy_score(tr_y, (tr_probs > 0.5).astype(int))

        model.eval()
        with torch.no_grad():
            v_logits = model(Xva.to(device), beta_va.to(device),
                              dxb_va.to(device),
                              chrom_va.to(device) if chrom_va is not None else None,
                              modality_abs=modality_abs_t,
                              func_imp_abs=func_imp_abs_t)
            if is_multiclass:
                v_loss = loss_fn(v_logits, yva.long().to(device)).item()
                v_prob = torch.softmax(v_logits, dim=-1).cpu().numpy()
                v_pred = v_prob.argmax(axis=1)
                v_m = _metrics_multi(yva.cpu().numpy(), v_prob, v_pred,
                                       v_loss, num_classes)
            else:
                v_loss = loss_fn(v_logits, yva.float().to(device)).item()
                v_prob = torch.sigmoid(v_logits).cpu().numpy()
                v_pred = (v_prob > 0.5).astype(int)
                v_m = _metrics(yva.cpu().numpy(), v_prob, v_pred, v_loss)
        history.append({"epoch": ep, "train_loss": train_loss,
                          "train_balanced_accuracy": tr_balacc,
                          **{f"val_{k}": v for k, v in v_m.items()}})

        if wandb_run is not None:
            log_entry = {
                "epoch": ep,
                "train_loss": train_loss,
                "train_balanced_accuracy": tr_balacc,
                "val_loss": v_loss,
                "val_balanced_accuracy": v_m["balanced_accuracy"],
                "val_roc_auc": v_m["roc_auc"],
            }
            # Log learnable ε's per epoch for interpretation
            if model.eps is not None:
                eps_vals = model.eps.detach().cpu().numpy()
                if model._eps_shape == "single":
                    log_entry["eps_func_importance"] = float(eps_vals[0])
                else:
                    for i, mod_name in enumerate(MODALITIES):
                        log_entry[f"eps_{mod_name}"] = float(eps_vals[i])
            wandb_run.log(log_entry, step=ep)

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


def _eval_mlp(model, data, device, modality_abs_t, func_imp_abs_t,
                num_classes: int = 1) -> tuple[dict, np.ndarray]:
    Xe, beta_e, dxb_e, chrom_e, y_e = data
    if modality_abs_t is not None: modality_abs_t = modality_abs_t.to(device)
    if func_imp_abs_t is not None: func_imp_abs_t = func_imp_abs_t.to(device)
    model.eval()
    with torch.no_grad():
        logits = model(Xe.to(device), beta_e.to(device), dxb_e.to(device),
                        chrom_e.to(device) if chrom_e is not None else None,
                        modality_abs=modality_abs_t, func_imp_abs=func_imp_abs_t)
        if num_classes == 1:
            prob = torch.sigmoid(logits).cpu().numpy()
            pred = (prob > 0.5).astype(int)
        else:
            prob = torch.softmax(logits, dim=-1).cpu().numpy()
            pred = prob.argmax(axis=1)
    y = y_e.cpu().numpy()
    if num_classes == 1:
        pos_w = max(1.0, (y == 0).sum() / max(1, (y == 1).sum()))
        loss = log_loss(y, np.clip(prob, 1e-7, 1 - 1e-7),
                         sample_weight=np.where(y == 1, pos_w, 1.0))
        return _metrics(y, prob, pred, loss_val=loss), prob
    else:
        # Macro log-loss (no class weighting at eval time).
        loss = log_loss(y, np.clip(prob, 1e-7, 1 - 1e-7),
                         labels=list(range(num_classes)))
        return _metrics_multi(y, prob, pred, loss_val=loss,
                                 num_classes=num_classes), prob


def _fit_tree(head: str, train_emb, train_y, val_emb, val_y, seed: int,
                xgb_n_estimators: int = 500, xgb_max_depth: int = 4,
                xgb_lr: float = 0.05, wandb_run=None,
                out_dir: Path | None = None) -> dict:
    pos = int((train_y == 1).sum()); neg = int((train_y == 0).sum())
    pos_w = neg / max(1, pos)
    if head == "rf":
        from sklearn.ensemble import RandomForestClassifier
        clf = RandomForestClassifier(n_estimators=500, n_jobs=-1, random_state=seed,
                                       class_weight="balanced")
        clf.fit(train_emb, train_y)
    elif head == "xgb":
        from xgboost import XGBClassifier
        clf = XGBClassifier(n_estimators=xgb_n_estimators, max_depth=xgb_max_depth,
                             learning_rate=xgb_lr,
                             tree_method="hist", scale_pos_weight=pos_w,
                             eval_metric="logloss", random_state=seed)
        # Pass BOTH train and val as eval sets → train_logloss + val_logloss curves
        clf.fit(train_emb, train_y,
                eval_set=[(train_emb, train_y), (val_emb, val_y)], verbose=False)
        # Per-round logging to wandb → train/val curves visible in Charts tab
        if wandb_run is not None:
            evr = clf.evals_result_  # {'validation_0': {'logloss': [...]}, 'validation_1': {...}}
            train_curve = evr.get("validation_0", {}).get("logloss", [])
            val_curve = evr.get("validation_1", {}).get("logloss", [])
            for i in range(max(len(train_curve), len(val_curve))):
                payload = {"xgb_round": i}
                if i < len(train_curve): payload["train_logloss"] = float(train_curve[i])
                if i < len(val_curve):   payload["val_logloss"] = float(val_curve[i])
                wandb_run.log(payload, step=i)
        # Persist the XGB model for later interrogation / drop-splice ablation
        if out_dir is not None:
            clf.save_model(str(out_dir / "xgb_model.json"))
    else:
        raise ValueError(head)
    prob = clf.predict_proba(val_emb)[:, 1]
    pred = (prob > 0.5).astype(int)
    loss = log_loss(val_y, np.clip(prob, 1e-7, 1-1e-7),
                     sample_weight=np.where(val_y == 1, pos_w, 1.0))
    return {"val_metrics": _metrics(val_y, prob, pred, loss),
            "val_prob": prob, "val_pred": pred, "classifier": clf}


# ── Main ──────────────────────────────────────────────────────────────────

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
    ap.add_argument("--func-integration-mode", "--func_integration_mode",
                    choices=ALL_MODES, default="none",
                    dest="func_integration_mode")
    ap.add_argument("--aggregation", choices=["global_attn", "chrom_hier"],
                    default="global_attn")
    ap.add_argument("--head", choices=["mlp2", "mlp3", "rf", "xgb"], required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--epochs", type=int, default=80)
    ap.add_argument("--batch-size", "--batch_size", type=int, default=32, dest="batch_size")
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--weight-decay", "--weight_decay", type=float, default=1e-4,
                    dest="weight_decay")
    ap.add_argument("--patience", type=int, default=15)
    # MLP head HPs (Sweep A)
    ap.add_argument("--mlp-width", "--mlp_width", type=int, default=256,
                    dest="mlp_width", help="Hidden width of MLP head (mlp2: hidden; mlp3: h1).")
    ap.add_argument("--mlp-dropout", "--mlp_dropout", type=float, default=0.3,
                    dest="mlp_dropout", help="Dropout in MLP head.")
    # XGB head HPs (Sweep B)
    ap.add_argument("--xgb-n-estimators", "--xgb_n_estimators", type=int, default=500,
                    dest="xgb_n_estimators", help="XGB n_estimators.")
    ap.add_argument("--xgb-max-depth", "--xgb_max_depth", type=int, default=4,
                    dest="xgb_max_depth", help="XGB max_depth.")
    ap.add_argument("--xgb-lr", "--xgb_lr", type=float, default=0.05,
                    dest="xgb_lr", help="XGB learning_rate.")
    ap.add_argument("--drop-splice", "--drop_splice", action="store_true",
                    dest="drop_splice",
                    help="Zero the splice modality (46/128 SNPs missing). "
                          "Affects bias path + value-mult path.")
    ap.add_argument("--exclude-cn-to-ad", "--exclude_cn_to_ad",
                    action="store_true", dest="exclude_cn_to_ad",
                    help="Drop CN_to_AD subgroup from POS class "
                          "(seed-1 val has 0 CN_to_AD subjects — noisy comparison).")
    ap.add_argument("--save-predictions", "--save_predictions",
                    action="store_true", dest="save_predictions",
                    help="Write per-patient val + test predictions to "
                          "<out_dir>/predictions.tsv (cols: Patient_ID, split, y, prob).")
    # ── New-task re-extraction args (2026-06-02) ───────────────────────────
    ap.add_argument("--task", default="ad_vs_cn",
                    choices=list(TASK_REGISTRY.keys()),
                    help="Which task's labels to train on. Default 'ad_vs_cn' "
                          "matches the original v3 sweep. Other choices "
                          "(sCN_vs_pCN, T3_horizon_{3,5,7,10}y, T4_3class) "
                          "select task-specific label loaders + cohort filters.")
    ap.add_argument("--num-classes", "--num_classes", type=int, default=None,
                    dest="num_classes",
                    help="Override the task registry's num_classes (1=binary, "
                          ">1=multiclass softmax). Defaults to the registry value.")
    ap.add_argument("--save-embeddings", "--save_embeddings",
                    action="store_true", dest="save_embeddings",
                    help="After training, run forward passes with hooks and "
                          "save embeddings.npz (per-PTID embedding_771 + "
                          "embedding_256 + attn_weights) alongside model.pt.")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--wandb-project", "--wandb_project", default=None,
                    dest="wandb_project")
    ap.add_argument("--wandb-entity", "--wandb_entity", default=None, dest="wandb_entity")
    ap.add_argument("--output-root", "--output_root", type=Path,
                    default=Path("outputs/diff_attn_v3"), dest="output_root")
    args = ap.parse_args()
    np.random.seed(args.seed); torch.manual_seed(args.seed)

    # ── WandB init ───────────────────────────────────────────────────────
    wb = None
    if args.wandb_project:
        import wandb
        _name = (f"{args.seq_length}_{args.model}_{args.func_integration_mode}"
                  f"_{args.aggregation}_{args.head}_s{args.seed}")
        if args.drop_splice:
            _name += "_no_splice"
        if args.exclude_cn_to_ad:
            _name += "_no_cn_to_ad"
        wb = wandb.init(project=args.wandb_project, entity=args.wandb_entity,
                         config=vars(args), name=_name)

    # ── Load core tensors ────────────────────────────────────────────────
    print(f"[load] diff: {args.seq_length} × {args.model}")
    diff, diff_rsids = _load_diff(args.base, args.seq_length, args.model)
    print(f"  diff shape: {diff.shape}")

    print(f"[load] dosage matrix")
    dos_df, _ = _load_dosage(args.base)
    print(f"  dosage shape: {dos_df.shape}")

    print(f"[load] beta_A1 table")
    beta_df = _load_beta(args.base)

    target_rsids = diff_rsids
    beta_df_idx = beta_df.set_index("rsID")
    beta = np.array([float(beta_df_idx.loc[rs, "beta_A1"]) if rs in beta_df_idx.index else 0.0
                     for rs in target_rsids], dtype=np.float32)
    chrom_arr = np.array([int(beta_df_idx.loc[rs, "CHR"]) if rs in beta_df_idx.index else 0
                          for rs in target_rsids], dtype=np.int64)
    chrom_uniq = sorted(set(chrom_arr.tolist()))
    chrom_idx_map = {c: i for i, c in enumerate(chrom_uniq)}
    chrom_idx = np.array([chrom_idx_map[c] for c in chrom_arr], dtype=np.int64)
    n_chroms = len(chrom_uniq)
    print(f"  n_chroms: {n_chroms}")

    # ── Load func resources (only what this mode needs) ─────────────────
    mode = args.func_integration_mode
    needs_per_snp_23  = (mode == "per_snp_zscored_23")
    needs_summaries   = mode in {"per_snp_summaries_8", "attn_bias_single",
                                   "attn_bias_per_modality", "value_mult",
                                   "attn_bias_per_modality_plus_value_mult"}
    ag_parent = _resolve_ag_parent(args.base) if (needs_per_snp_23 or needs_summaries) else None

    func_raw23 = None
    summaries = None
    if needs_per_snp_23:
        print(f"[load] 23-D functional features ({ag_parent})")
        func_raw23, _ = ff.load_functional_features(ag_parent, target_rsids)
        print(f"  func shape: {func_raw23.shape}")
    if needs_summaries:
        print(f"[load] modality summaries ({ag_parent})")
        summaries = ff.load_modality_summaries(ag_parent, target_rsids)
        for k, v in summaries.items():
            print(f"  {k}: shape={v.shape}, nonzero={int((v!=0).sum())}")

    # ── Labels per split (per-task dispatcher; ad_vs_cn is the original) ──
    print(f"[labels] task={args.task} splits_root={args.splits_root}, seed={args.seed}")
    task_entry = TASK_REGISTRY[args.task]
    task_loader = task_entry["loader"]
    num_classes = (args.num_classes if args.num_classes is not None
                    else task_entry["num_classes"])
    labels = {sp: task_loader(args.splits_root, args.seed, sp,
                                args.base, exclude_cn_to_ad=args.exclude_cn_to_ad)
               for sp in ("train", "val", "test")}
    for sp, df in labels.items():
        if num_classes == 1:
            print(f"  {sp}: n={len(df)}, pos={(df.y==1).sum()}, neg={(df.y==0).sum()}")
        else:
            counts = df["y"].value_counts().sort_index().to_dict()
            print(f"  {sp}: n={len(df)}, by_class={counts}")

    dos_by_ptid = dos_df.set_index("PTID")
    splits = {}
    for sp, df in labels.items():
        ptids = df["PTID"].tolist()
        present = [p for p in ptids if p in dos_by_ptid.index]
        df_ok = df[df["PTID"].isin(present)].reset_index(drop=True)
        dos_split = dos_by_ptid.loc[df_ok["PTID"], target_rsids].values.astype(np.float32)
        splits[sp] = {"ptid": df_ok["PTID"].tolist(),
                       "y": df_ok["y"].values.astype(np.int64),
                       "dosage": dos_split}

    train_means = np.nanmean(splits["train"]["dosage"], axis=0)
    for sp in splits:
        d = splits[sp]["dosage"]
        nan_mask = np.isnan(d)
        for i in range(d.shape[1]):
            d[nan_mask[:, i], i] = train_means[i]
        splits[sp]["dosage"] = d

    # ── Z-score func features (train-fold rows = all SNPs present in train) ──
    # All target SNPs participate in train; so train_snp_mask is all-True.
    train_snp_mask = np.ones(len(target_rsids), dtype=bool)
    zscore_stats = {}
    func_raw23_z = None
    summaries_z = None
    if needs_per_snp_23:
        func_raw23_z, stats = _zscore_train(func_raw23, train_snp_mask)
        zscore_stats["func_raw23"] = stats
    if needs_summaries:
        summaries_z = {}
        for k, v in summaries.items():
            zv, st = _zscore_train(v, train_snp_mask)
            summaries_z[k] = zv.reshape(-1)
            zscore_stats[k] = st
        if args.drop_splice:
            for k in ("splice_abs", "splice_signed"):
                if k in summaries_z:
                    summaries_z[k] = np.zeros_like(summaries_z[k])
            print("  [drop-splice] zeroed splice_abs + splice_signed "
                  "(modality contribution → 0)")

    # ── Build per-(p, i, F_full) value tensors ──────────────────────────
    def _build(sp):
        d = splits[sp]
        X = _per_patient_values(mode, beta, d["dosage"],
                                  func_raw23_z, summaries_z, diff)
        beta_t = beta.copy()
        dxb = beta[None, :] * d["dosage"]
        return X, beta_t, dxb, d["y"]

    Xtr, beta_tr, dxb_tr, ytr = _build("train")
    Xva, beta_va, dxb_va, yva = _build("val")
    Xte, beta_te, dxb_te, yte = _build("test")
    F = Xtr.shape[2]
    print(f"[feat] in_dim per SNP F={F}; train n={Xtr.shape[0]}, val n={Xva.shape[0]}")

    # ── Bias tensors (MLP path: pass to model; tree path: pre-apply to aggregator) ──
    # Per-modality abs (S, 4) and func_importance abs (S,) z-scored.
    modality_abs_np = None
    func_imp_abs_np = None
    extra_bias_tree_np = None
    if needs_summaries:
        modality_abs_np = np.stack([summaries_z[f"{m}_abs"] for m in MODALITIES],
                                      axis=1).astype(np.float32)         # (S, 4)
        func_imp_abs_np = np.maximum.reduce([summaries_z[f"{m}_abs"]
                                               for m in MODALITIES]).astype(np.float32)
        extra_bias_tree_np = _build_extra_bias_static(mode, summaries_z,
                                                        P=splits["train"]["dosage"].shape[0],
                                                        S=len(target_rsids))

    out_dir = args.output_root / args.head / (f"{args.seq_length}_{args.model}_"
                                                f"{mode}_{args.aggregation}_s{args.seed}")
    out_dir.mkdir(parents=True, exist_ok=True)
    if zscore_stats:
        (out_dir / "func_zscore_stats.json").write_text(json.dumps(zscore_stats, indent=2))

    # ── Train + eval ────────────────────────────────────────────────────
    if args.head in ("mlp2", "mlp3"):
        Xtr_t = torch.from_numpy(Xtr)
        Xva_t = torch.from_numpy(Xva)
        Xte_t = torch.from_numpy(Xte)
        beta_t = torch.from_numpy(beta_tr)
        dxb_tr_t = torch.from_numpy(dxb_tr.astype(np.float32))
        dxb_va_t = torch.from_numpy(dxb_va.astype(np.float32))
        dxb_te_t = torch.from_numpy(dxb_te.astype(np.float32))
        chrom_t = torch.from_numpy(chrom_idx) if args.aggregation == "chrom_hier" else None
        ytr_t = torch.from_numpy(ytr); yva_t = torch.from_numpy(yva); yte_t = torch.from_numpy(yte)
        modality_abs_t = torch.from_numpy(modality_abs_np) if modality_abs_np is not None else None
        func_imp_abs_t = torch.from_numpy(func_imp_abs_np) if func_imp_abs_np is not None else None

        if num_classes == 1:
            pos = int((ytr == 1).sum()); neg = int((ytr == 0).sum())
            pos_weight = neg / max(1, pos)
            print(f"  pos_weight (neg/pos): {pos_weight:.3f}")
            class_weights = None
        else:
            # Inverse-frequency class weights for CrossEntropyLoss.
            pos_weight = 1.0  # unused
            ytr_arr = np.asarray(ytr)
            counts = np.bincount(ytr_arr, minlength=num_classes).astype(float)
            counts = np.where(counts == 0, 1.0, counts)
            class_weights = (counts.sum() / (num_classes * counts)).astype(np.float32)
            print(f"  class_weights (1/freq, normalised): {class_weights.tolist()}")

        model = DiffAttnV3(in_dim=F, n_chroms=n_chroms,
                            aggregation=args.aggregation, head=args.head, mode=mode,
                            mlp_width=args.mlp_width, mlp_dropout=args.mlp_dropout,
                            out_dim=num_classes)
        fit_info = _train_mlp(model,
                                (Xtr_t, beta_t, dxb_tr_t, chrom_t, ytr_t),
                                (Xva_t, beta_t, dxb_va_t, chrom_t, yva_t),
                                epochs=args.epochs, batch_size=args.batch_size,
                                lr=args.lr, weight_decay=args.weight_decay,
                                patience=args.patience, device=args.device,
                                pos_weight=pos_weight,
                                modality_abs_t=modality_abs_t,
                                func_imp_abs_t=func_imp_abs_t,
                                wandb_run=wb,
                                num_classes=num_classes,
                                class_weights=class_weights)
        val_m, val_prob = _eval_mlp(model, (Xva_t, beta_t, dxb_va_t, chrom_t, yva_t),
                              args.device, modality_abs_t, func_imp_abs_t,
                              num_classes=num_classes)
        test_m, test_prob = _eval_mlp(model, (Xte_t, beta_t, dxb_te_t, chrom_t, yte_t),
                              args.device, modality_abs_t, func_imp_abs_t,
                              num_classes=num_classes)
        torch.save(model.state_dict(), out_dir / "model.pt")
    else:
        # Tree path: build patient embeddings with fixed γ=δ=ε=1 aggregator
        emb_tr = _fixed_aggregator_v3(Xtr, beta_tr, splits["train"]["dosage"],
                                        extra_bias_tree_np)
        emb_va = _fixed_aggregator_v3(Xva, beta_va, splits["val"]["dosage"],
                                        extra_bias_tree_np)
        emb_te = _fixed_aggregator_v3(Xte, beta_te, splits["test"]["dosage"],
                                        extra_bias_tree_np)
        print(f"  tree emb dims: train {emb_tr.shape}, val {emb_va.shape}")
        result = _fit_tree(args.head, emb_tr, ytr, emb_va, yva, args.seed,
                              xgb_n_estimators=args.xgb_n_estimators,
                              xgb_max_depth=args.xgb_max_depth,
                              xgb_lr=args.xgb_lr,
                              wandb_run=wb, out_dir=out_dir)
        val_m = result["val_metrics"]
        val_prob = result.get("val_prob")
        clf = result["classifier"]
        test_prob = clf.predict_proba(emb_te)[:, 1]
        test_pred = (test_prob > 0.5).astype(int)
        pos_w = (yte == 0).sum() / max(1, (yte == 1).sum())
        test_loss = log_loss(yte, np.clip(test_prob, 1e-7, 1-1e-7),
                              sample_weight=np.where(yte == 1, pos_w, 1.0))
        test_m = _metrics(yte, test_prob, test_pred, test_loss)
        fit_info = {"head": args.head}

    metrics = {"val": val_m, "test": test_m, "config": vars(args)}
    (out_dir / "metrics.json").write_text(json.dumps(metrics, default=str, indent=2))

    # Per-patient predictions. Binary: single `prob` column (original).
    # Multiclass: `prob_class0..K-1` columns.
    if args.save_predictions and val_prob is not None and test_prob is not None:
        def _common_meta():
            return {"seed": args.seed, "model": args.model, "task": args.task,
                     "num_classes": num_classes,
                     "func_integration_mode": args.func_integration_mode,
                     "head": args.head, "seq_length": args.seq_length,
                     "aggregation": args.aggregation,
                     "mlp_width": args.mlp_width,
                     "mlp_dropout": args.mlp_dropout, "lr": args.lr}
        pred_rows = []
        for sp_name, sp_ptids, sp_y, sp_prob in [
                ("val",  splits["val"]["ptid"],  yva, np.asarray(val_prob)),
                ("test", splits["test"]["ptid"], yte, np.asarray(test_prob))]:
            for i, (ptid, y) in enumerate(zip(sp_ptids, sp_y)):
                row = {"Patient_ID": str(ptid), "split": sp_name,
                        "y": int(y), **_common_meta()}
                if num_classes == 1:
                    row["prob"] = float(np.asarray(sp_prob).ravel()[i])
                else:
                    for k in range(num_classes):
                        row[f"prob_class{k}"] = float(sp_prob[i, k])
                pred_rows.append(row)
        pd.DataFrame(pred_rows).to_csv(out_dir / "predictions.tsv",
                                          sep="\t", index=False)

    # ── Post-training embedding extraction (--save-embeddings) ──────────
    # Runs forward passes with hooks on the trained pool + MLP penultimate.
    # Saves per-PTID embedding_771 + embedding_256 (None for tree heads).
    if args.save_embeddings and args.head in ("mlp2", "mlp3"):
        print(f"[embeddings] extracting per-PTID embedding_771 + embedding_256...")
        captured = {sp: {"pool_out": None, "pool_attn": None,
                          "head_hidden": None, "prob": None}
                     for sp in ("train", "val", "test")}
        device_t = args.device
        cur_sp = {"value": None}    # closure-mutable

        def _h_pool(_m, _inp, output):
            sp = cur_sp["value"]
            if isinstance(output, tuple):
                captured[sp]["pool_out"]  = output[0].detach().cpu().numpy()
                attn = output[1]
                if isinstance(attn, dict):
                    if "a_snp" in attn and attn["a_snp"] is not None:
                        captured[sp]["pool_attn"] = attn["a_snp"].detach().cpu().numpy()
                elif attn is not None:
                    captured[sp]["pool_attn"] = attn.detach().cpu().numpy()
            else:
                captured[sp]["pool_out"] = output.detach().cpu().numpy()

        def _h_head(_m, _inp, output):
            sp = cur_sp["value"]
            captured[sp]["head_hidden"] = output.detach().cpu().numpy()

        h1 = model.pool.register_forward_hook(_h_pool)
        h2 = model.head.net[2].register_forward_hook(_h_head)
        try:
            model.eval()
            for sp_name, (Xe, be, dxbe, ye) in (
                    ("train", (Xtr_t, beta_t, dxb_tr_t, ytr_t)),
                    ("val",   (Xva_t, beta_t, dxb_va_t, yva_t)),
                    ("test",  (Xte_t, beta_t, dxb_te_t, yte_t))):
                cur_sp["value"] = sp_name
                with torch.no_grad():
                    logits = model(Xe.to(device_t), be.to(device_t),
                                     dxbe.to(device_t),
                                     chrom_t.to(device_t) if chrom_t is not None else None,
                                     modality_abs=modality_abs_t,
                                     func_imp_abs=func_imp_abs_t)
                    if num_classes == 1:
                        captured[sp_name]["prob"] = torch.sigmoid(logits).cpu().numpy()
                    else:
                        captured[sp_name]["prob"] = torch.softmax(
                            logits, dim=-1).cpu().numpy()
        finally:
            h1.remove(); h2.remove()
        # Concat all folds + write embeddings.npz
        ptids_all, fold_all, y_all, prob_all = [], [], [], []
        emb771_all, emb256_all, attn_all = [], [], []
        for sp_name, y_arr in (("train", ytr), ("val", yva), ("test", yte)):
            cap = captured[sp_name]
            n = cap["pool_out"].shape[0]
            ptids_all.extend(splits[sp_name]["ptid"])
            fold_all.extend([sp_name] * n)
            y_all.append(np.asarray(y_arr, dtype=np.int8))
            prob_all.append(np.asarray(cap["prob"], dtype=np.float32))
            emb771_all.append(cap["pool_out"].astype(np.float32))
            emb256_all.append(cap["head_hidden"].astype(np.float32))
            if cap["pool_attn"] is not None:
                attn_all.append(cap["pool_attn"].astype(np.float32))
        save_kwargs = dict(
            ptids=np.array(ptids_all),
            fold=np.array(fold_all),
            y_true=np.concatenate(y_all),
            proba=np.concatenate(prob_all, axis=0),
            embedding_771=np.concatenate(emb771_all, axis=0),
            embedding_256=np.concatenate(emb256_all, axis=0),
        )
        if attn_all:
            save_kwargs["attn_weights"] = np.concatenate(attn_all, axis=0)
        np.savez_compressed(out_dir / "embeddings.npz", **save_kwargs)
        print(f"  wrote {out_dir/'embeddings.npz'}  "
              f"e771={save_kwargs['embedding_771'].shape}  "
              f"e256={save_kwargs['embedding_256'].shape}")

    log_payload = {f"val_{k}": v for k, v in val_m.items()}
    log_payload.update({f"test_{k}": v for k, v in test_m.items()})
    print(f"\n=== val balAcc={val_m['balanced_accuracy']:.4f}  "
          f"ROC-AUC={val_m['roc_auc']:.4f}  loss={val_m.get('loss', float('nan')):.4f} ===")
    if wb:
        wb.log(log_payload)
        for k, v in log_payload.items():
            wb.summary[k] = v
        wb.finish()


if __name__ == "__main__":
    main()
