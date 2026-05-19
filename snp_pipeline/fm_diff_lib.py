r"""
fm_diff_lib.py
==============
Shared library for the frozen-DNA-FM **diff-embedding → β-biased attention →
AD/CN classifier** pipeline (scripts 26 → 27 → 28 → 28b). It is the SINGLE
source of truth for: canonical SNP ordering, allele/dosage QC, the
per-patient signed-effect tensor

    d[p, i, :] = (ALT_emb_i − REF_emb_i) · beta_A1_i · dosage_{p,i}

the β-biased attention pooling modules (global + chromosome-hierarchical),
and a minimal fork of script-15's train/eval loop (the only change vs 15 is
the 4-tensor model forward + a `train` metrics block + per-patient
predictions). Importable directly (no leading digit):

    from fm_diff_lib import (load_snp_table, canonical_snp_order, build_d,
                             DiffAttentionModel, train_one_diff, evaluate_diff)

Correctness invariants (enforced here, see asserts):
  * `effect_allele == ALT == A1` is what `beta_A1` refers to AND what the
    dosage counts (script 25 export convention). diff = ALT_emb − REF_emb,
    so `d = diff · beta_A1 · dosage` is sign-consistent — never re-sign.
  * SNP order is canonical (CHR-numeric → BP → rsID) everywhere; the diff
    .npz carries `rsids` so column order is self-describing and asserted.
  * Missing dosage is cohort-mean imputed using **train-fold** means only,
    propagated to val/test (no leakage) — matches scripts 20/25.

Reuses verbatim (cited): MLPHead — 15_train_classification_head.py L172-183;
MLP3Head — 19_compare_frozen_models.py L112-124; AdditiveAttention scoring
shape — 15 L74-88; chrom skip-empty trick — 15 L147-166; metric dict — 15
L350-368; early-stop loop — 15 L288-347.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.metrics import (
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from torch.utils.data import DataLoader, TensorDataset

# ── Canonical paths (Windows local; overridable via fn args) ─────────────────
BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
GWAS_FILTERED = BASE / "GWAS_filtered"
CONV_LABELS = BASE / "conversion_labels" / "conversion_labels.tsv"
SPLITS_ROOT = Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                    "no_cdr_stratified/tabular/baseline")
SNP_SETS = {"9W27B": 32, "9W27B5D": 37}          # nominal (pre-QC) sizes
MODEL_DIMS = {"bmfm_ref": 768, "bmfm_snp": 768, "ntv2": 1024,
              "caduceus_ph": 256}
VALID_BASES = set("ACGT")
_MISSING_AL = {"", ".", "NA", "NAN", "NONE", "<NA>"}

SNP_TSV_COLS = ["rsID", "gene", "lead_rsID", "source", "origin", "CHR",
                "BP_GRCh38", "effect_allele", "other_allele", "REF", "ALT",
                "beta_A1"]


# ════════════════════════════════════════════════════════════════════════════
# Data: SNP table, ordering, alleles, dosage, labels, splits
# ════════════════════════════════════════════════════════════════════════════
def _chr_key(c) -> int:
    s = str(c).replace("chr", "").strip().upper()
    if s.isdigit():
        return int(s)
    return {"X": 23, "Y": 24, "MT": 25, "M": 25}.get(s, 99)


def load_snp_table(snps_tsv: Path) -> pd.DataFrame:
    """Read a `<SET>_snps.tsv` (script-25 export). Alleles as str; beta_A1
    float (NaN on coerce). No filtering here."""
    df = pd.read_csv(snps_tsv, sep="\t", dtype=str)
    for c in SNP_TSV_COLS:
        if c not in df.columns:
            df[c] = np.nan
    df["beta_A1"] = pd.to_numeric(df["beta_A1"], errors="coerce")
    df["BP_GRCh38"] = pd.to_numeric(df["BP_GRCh38"], errors="coerce")
    for a in ("effect_allele", "other_allele", "REF", "ALT"):
        df[a] = df[a].astype(str).str.strip().str.upper()
    df["rsID"] = df["rsID"].astype(str).str.strip()
    df["CHR"] = df["CHR"].astype(str).str.strip()
    return df


def canonical_snp_order(df: pd.DataFrame) -> pd.DataFrame:
    """Deterministic order: CHR numeric → BP_GRCh38 → rsID. Adds int
    `snp_idx`. The one ordering used by 26/27/28 so column order is fixed."""
    d = df.copy()
    d["_ck"] = d["CHR"].map(_chr_key)
    d = (d.sort_values(["_ck", "BP_GRCh38", "rsID"], kind="mergesort")
           .drop(columns="_ck").reset_index(drop=True))
    d.insert(0, "snp_idx", np.arange(len(d), dtype=int))
    return d


def validate_alleles(ea: str, oa: str) -> tuple[bool, str]:
    ea = ("" if ea is None else str(ea)).strip().upper()
    oa = ("" if oa is None else str(oa)).strip().upper()
    if ea in _MISSING_AL or oa in _MISSING_AL:
        return False, "missing_allele"
    if ea not in VALID_BASES or oa not in VALID_BASES:
        return False, "non_ACGT_allele"
    if ea == oa:
        return False, "ea_eq_oa"
    return True, "ok"


def load_dosage(dosage_tsv: Path) -> pd.DataFrame:
    """`<SET>_patient_dosage.tsv` → index=PTID(str), cols=rsID, float;
    blank/'' → NaN. No imputation here."""
    df = pd.read_csv(dosage_tsv, sep="\t", dtype={"PTID": str})
    df["PTID"] = df["PTID"].astype(str).str.strip()
    df = df.set_index("PTID")
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def cohort_mean_impute(dosage: pd.DataFrame, snp_cols: list[str],
                       means: pd.Series | None = None
                       ) -> tuple[pd.DataFrame, pd.Series, dict]:
    """Per-SNP mean-impute missing dosage. If `means` is None they are
    computed from `dosage` (call on the TRAIN fold) and returned; pass the
    returned Series back in for val/test so no test info leaks (scripts
    20/25 convention). Returns (imputed_df, means, {rsID: n_imputed})."""
    sub = dosage.reindex(columns=snp_cols).astype(float)
    if means is None:
        means = sub.mean(axis=0, skipna=True)
    means = means.reindex(snp_cols)
    means = means.fillna(0.0)                     # SNP all-missing in train
    n_imp = {c: int(sub[c].isna().sum()) for c in snp_cols}
    return sub.fillna(means), means, n_imp


def load_cohort_labels(conv_labels: Path = CONV_LABELS) -> pd.DataFrame:
    """conversion_labels.tsv → df[Patient_ID, y] for the AD/CN cohort.
    Cohort = ad_case==1 (y=1) OR stable_CN==1 (y=0); asserts mutual
    exclusivity (script-22 definitions)."""
    df = pd.read_csv(conv_labels, sep="\t", dtype={"Patient_ID": str})
    df["Patient_ID"] = df["Patient_ID"].astype(str).str.strip()
    for c in ("ad_case", "stable_CN"):
        if c not in df.columns:
            raise ValueError(f"{conv_labels} lacks column {c!r}")
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    both = int(((df["ad_case"] == 1) & (df["stable_CN"] == 1)).sum())
    if both:
        raise ValueError(f"{both} patients have ad_case==1 AND stable_CN==1 "
                         "— labels not mutually exclusive")
    coh = df[(df["ad_case"] == 1) | (df["stable_CN"] == 1)].copy()
    coh["y"] = coh["ad_case"].astype(int)         # 1=AD case, 0=stable CN
    return coh[["Patient_ID", "y"]].reset_index(drop=True)


def split_membership(seed: int, splits_root: Path = SPLITS_ROOT
                     ) -> dict[str, set]:
    """seed_<seed>/{train,val,test}.csv → {split: set(Patient_ID)}."""
    out = {}
    for sp in ("train", "val", "test"):
        f = splits_root / f"seed_{seed}" / f"{sp}.csv"
        ids = pd.read_csv(f, usecols=["Patient_ID"])["Patient_ID"]
        out[sp] = set(ids.astype(str).str.strip())
    return out


def load_diff(diff_npz: Path, snp_meta: pd.DataFrame
              ) -> tuple[np.ndarray, np.ndarray]:
    """Load per-SNP diff (= ALT_emb − REF_emb) aligned to `snp_meta` row
    order. The .npz (script 27) has `rsids` (canonical order), `ref_emb`,
    `alt_emb`. Hard-asserts rsID alignment (the single load-bearing
    invariant). Returns (diff (S,H) float32, rsids (S,))."""
    z = np.load(diff_npz, allow_pickle=True)
    rsids = np.array([str(r) for r in z["rsids"]])
    ref = np.asarray(z["ref_emb"], dtype=np.float32)
    alt = np.asarray(z["alt_emb"], dtype=np.float32)
    want = snp_meta["rsID"].astype(str).tolist()
    have = list(rsids)
    if have != want:
        # allow npz to be a superset / different order → reindex by rsID
        pos = {r: i for i, r in enumerate(have)}
        missing = [r for r in want if r not in pos]
        if missing:
            raise ValueError(f"{diff_npz}: missing {len(missing)} rsIDs "
                             f"e.g. {missing[:5]}")
        order = [pos[r] for r in want]
        ref, alt, rsids = ref[order], alt[order], np.array(want)
    diff = (alt - ref).astype(np.float32)
    if not np.isfinite(diff).all():
        raise ValueError(f"{diff_npz}: non-finite values in ALT−REF diff")
    return diff, rsids


def build_d(diff: np.ndarray, beta: np.ndarray,
            dosage: np.ndarray) -> np.ndarray:
    """d[p,i,:] = diff[i,:] · beta[i] · dosage[p,i].  Shapes: diff (S,H),
    beta (S,), dosage (P,S) → (P,S,H) float32. beta sign already in
    beta_A1 (never re-sign). Zero dosage ⇒ zero d_i (definitional)."""
    diff = np.asarray(diff, np.float32)
    beta = np.asarray(beta, np.float32)
    dosage = np.asarray(dosage, np.float32)
    S, H = diff.shape
    assert beta.shape == (S,), f"beta {beta.shape} vs S={S}"
    assert dosage.shape[1] == S, f"dosage {dosage.shape} vs S={S}"
    for nm, arr in (("diff", diff), ("beta", beta), ("dosage", dosage)):
        if not np.isfinite(arr).all():
            raise ValueError(f"build_d: non-finite in {nm}")
    d = dosage[:, :, None] * beta[None, :, None] * diff[None, :, :]
    return d.astype(np.float32)


# ════════════════════════════════════════════════════════════════════════════
# Models: β-biased attention pooling (NEW) + reused MLP heads
# ════════════════════════════════════════════════════════════════════════════
class MLPHead(nn.Module):
    """2-layer head — verbatim 15_train_classification_head.py L172-183."""

    def __init__(self, in_dim: int, hidden: int = 256, dropout: float = 0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden), nn.GELU(), nn.Dropout(dropout),
            nn.Linear(hidden, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x).squeeze(-1)


class MLP3Head(nn.Module):
    """3-layer head — verbatim 19_compare_frozen_models.py L112-124."""

    def __init__(self, in_dim: int, h1: int = 256, h2: int = 128,
                 dropout: float = 0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, h1), nn.GELU(), nn.Dropout(dropout),
            nn.Linear(h1, h2), nn.GELU(), nn.Dropout(dropout),
            nn.Linear(h2, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x).squeeze(-1)


class SNPEncoder(nn.Module):
    """Optional per-SNP MLP d_i → h_i (identity if disabled)."""

    def __init__(self, in_dim: int, hidden: int = 128,
                 dropout: float = 0.1, enabled: bool = False):
        super().__init__()
        self.enabled = enabled
        self._in = in_dim
        if enabled:
            self.net = nn.Sequential(
                nn.Linear(in_dim, hidden), nn.GELU(), nn.Dropout(dropout),
                nn.Linear(hidden, in_dim),
            )

    @property
    def out_dim(self) -> int:
        return self._in

    def forward(self, d: torch.Tensor) -> torch.Tensor:
        return self.net(d) if self.enabled else d


class BiasedAttention(nn.Module):
    """Additive attention (scoring shape mirrors 15 AdditiveAttention
    L74-88) whose logits are biased by the genetic prior:

        score_i = u·tanh(v·x_i) + γ·|β_i| [+ δ·|dosage_i·β_i|]

    γ, δ are learnable scalars (nn.Parameter) initialised 0 so the model
    starts as plain learned attention and must *learn* to use the prior
    (debuggable: inspect γ/δ each epoch). Returns (pooled (B,D),
    weights (B,N)). All-masked rows → zero vector / zero weights (no
    -inf-softmax NaN)."""

    def __init__(self, dim: int, attn_hidden: int = 128,
                 use_delta: bool = False):
        super().__init__()
        self.v = nn.Linear(dim, attn_hidden, bias=False)
        self.u = nn.Linear(attn_hidden, 1, bias=False)
        self.gamma = nn.Parameter(torch.zeros(()))
        self.use_delta = use_delta
        self.delta = nn.Parameter(torch.zeros(())) if use_delta else None

    def forward(self, x, mask, abs_beta, abs_dosbeta=None):
        # x:(B,N,D) mask:(B,N) bool  abs_beta:(B,N)  abs_dosbeta:(B,N)|None
        score = self.u(torch.tanh(self.v(x))).squeeze(-1)        # (B,N)
        score = score + self.gamma * abs_beta
        if self.use_delta and abs_dosbeta is not None:
            score = score + self.delta * abs_dosbeta
        score = score.masked_fill(~mask, float("-inf"))
        a = torch.softmax(score, dim=-1)
        a = torch.nan_to_num(a, nan=0.0)                          # empty row
        pooled = (a.unsqueeze(-1) * x).sum(dim=1)                 # (B,D)
        return pooled, a


class GlobalAttnPool(nn.Module):
    """Mode A — one β-biased softmax over all of a patient's SNPs."""

    def __init__(self, dim: int, use_delta: bool, attn_hidden: int = 128):
        super().__init__()
        self.dim = dim
        self.attn = BiasedAttention(dim, attn_hidden, use_delta)

    @property
    def output_dim(self) -> int:
        return self.dim

    def forward(self, x, mask, chrom_idx, abs_beta, abs_dosbeta):
        emb, a = self.attn(x, mask, abs_beta, abs_dosbeta)
        return emb, {"a_snp": a}


class ChromHierPool(nn.Module):
    """Mode B — β-biased attention pools SNPs within each chromosome →
    one vector per chrom; a second (plain) attention pools chromosomes
    (numeric order). Skip-empty-chrom trick from 15 L147-166."""

    def __init__(self, dim: int, n_chroms: int, use_delta: bool,
                 attn_hidden: int = 128):
        super().__init__()
        self.dim = dim
        self.n_chroms = n_chroms
        self.intra = BiasedAttention(dim, attn_hidden, use_delta)
        self.inter = BiasedAttention(dim, attn_hidden, use_delta=False)

    @property
    def output_dim(self) -> int:
        return self.dim

    def forward(self, x, mask, chrom_idx, abs_beta, abs_dosbeta):
        B, S, D = x.shape
        chrom_emb = x.new_zeros(B, self.n_chroms, D)
        chrom_pres = torch.zeros(B, self.n_chroms, dtype=torch.bool,
                                 device=x.device)
        a_snp_full = x.new_zeros(B, S)
        for c in range(self.n_chroms):
            m = (chrom_idx == c) & mask                          # (B,S)
            has = m.any(dim=1)
            if not has.any():
                continue
            emb_c, a_c = self.intra(x, m, abs_beta, abs_dosbeta)
            chrom_emb[has, c] = emb_c[has]
            chrom_pres[has, c] = True
            a_snp_full[has] += a_c[has]            # disjoint chroms → add ok
        zeros = chrom_emb.new_zeros(B, self.n_chroms)
        patient, a_chrom = self.inter(chrom_emb, chrom_pres, zeros)
        return patient, {"a_snp": a_snp_full, "a_chrom": a_chrom}


class DiffAttentionModel(nn.Module):
    """SNP-encoder → β-biased pool → LayerNorm → MLP head → 1 logit.
    `forward` returns logits; `last_attn` / `last_emb` hold the most
    recent attention dict and pooled patient embedding (for saving)."""

    def __init__(self, *, model_dim: int, n_chroms: int, aggregation: str,
                 mlp_depth: str, use_delta: bool, snp_encoder: bool = False,
                 attn_hidden: int = 128):
        super().__init__()
        self.aggregation = aggregation
        self.encoder = SNPEncoder(model_dim, enabled=snp_encoder)
        d = self.encoder.out_dim
        if aggregation == "global_attn":
            self.pool = GlobalAttnPool(d, use_delta, attn_hidden)
        elif aggregation == "chrom_hier":
            self.pool = ChromHierPool(d, n_chroms, use_delta, attn_hidden)
        else:
            raise ValueError(f"unknown aggregation {aggregation!r}")
        ind = self.pool.output_dim
        self.norm = nn.LayerNorm(ind)             # scale guard (β·dosage)
        self.head = MLPHead(ind) if mlp_depth == "2L" else MLP3Head(ind)
        self.last_attn: dict | None = None
        self.last_emb: torch.Tensor | None = None

    def forward(self, d, mask, chrom_idx, abs_beta, abs_dosbeta):
        h = self.encoder(d)
        emb, attn = self.pool(h, mask, chrom_idx, abs_beta, abs_dosbeta)
        self.last_attn = attn
        self.last_emb = emb
        return self.head(self.norm(emb))


# ════════════════════════════════════════════════════════════════════════════
# Train / eval — minimal fork of 15 train_one/evaluate (4-tensor forward)
# ════════════════════════════════════════════════════════════════════════════
def _metrics(y_np: np.ndarray, logits: np.ndarray) -> dict:
    """Verbatim metric dict from 15 evaluate L356-368."""
    pred = (logits > 0).astype(int)
    out = {
        "auc": (float(roc_auc_score(y_np, logits))
                if len(np.unique(y_np)) > 1 else float("nan")),
        "balanced_accuracy": float(balanced_accuracy_score(y_np, pred)),
        "f1": float(f1_score(y_np, pred, zero_division=0)),
        "precision": float(precision_score(y_np, pred, pos_label=1,
                                           zero_division=0)),
        "recall": float(recall_score(y_np, pred, pos_label=1,
                                     zero_division=0)),
        "sensitivity": float(recall_score(y_np, pred, pos_label=1,
                                          zero_division=0)),
        "specificity": float(recall_score(y_np, pred, pos_label=0,
                                          zero_division=0)),
    }
    cm = confusion_matrix(y_np, pred, labels=[0, 1])
    out["confusion_matrix"] = {"tn": int(cm[0, 0]), "fp": int(cm[0, 1]),
                               "fn": int(cm[1, 0]), "tp": int(cm[1, 1])}
    return out


def train_one_diff(model: DiffAttentionModel, train_data, val_data,
                    *, epochs: int, batch_size: int, lr: float,
                    weight_decay: float, device: torch.device,
                    patience: int, loss_name: str = "weighted_bce"
                    ) -> tuple[DiffAttentionModel, dict, list]:
    """Fork of 15.train_one L288-347: data tuples are 5-tensor
    (d, mask, chrom_idx, abs_beta, abs_dosbeta, y); AdamW; (weighted-)BCE
    with pos_weight=n_neg/n_pos from train; early-stop on val AUC."""
    opt = torch.optim.AdamW(model.parameters(), lr=lr,
                            weight_decay=weight_decay)
    y_tr_all = train_data[-1]
    if loss_name == "weighted_bce":
        n_pos = float((y_tr_all == 1).sum().item())
        n_neg = float((y_tr_all == 0).sum().item())
        pw = max(n_neg / max(n_pos, 1.0), 1e-6)
        loss_fn = nn.BCEWithLogitsLoss(
            pos_weight=torch.tensor([pw], device=device))
        print(f"[loss] weighted_bce pos_weight={pw:.4f} "
              f"(n_pos={int(n_pos)}, n_neg={int(n_neg)})")
    elif loss_name == "bce":
        loss_fn = nn.BCEWithLogitsLoss()
        print("[loss] bce (unweighted)")
    else:
        raise ValueError(f"unknown loss {loss_name!r}")

    tr = [t.to(device) for t in train_data]
    va = [t.to(device) for t in val_data]
    loader = DataLoader(TensorDataset(*tr), batch_size=batch_size,
                        shuffle=True)
    best_state, best_auc, since, curve = None, -1.0, 0, []
    for ep in range(epochs):
        model.train()
        for db, mb, cb, bb, ab, yb in loader:
            opt.zero_grad()
            logits = model(db, mb, cb, bb, ab)
            loss_fn(logits, yb.float()).backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            vl_t = model(va[0], va[1], va[2], va[3], va[4])
            vloss = float(loss_fn(vl_t, va[5].float()).item())
            vl = vl_t.cpu().numpy()
        yv = va[5].cpu().numpy()
        vauc = (float(roc_auc_score(yv, vl)) if len(np.unique(yv)) > 1
                else float("nan"))
        g = float(model.pool.attn.gamma.detach().cpu()) \
            if hasattr(model.pool, "attn") \
            else float(model.pool.intra.gamma.detach().cpu())
        curve.append({"epoch": ep, "val_auc": vauc, "val_loss": vloss,
                      "gamma": g})
        if np.isfinite(vauc) and vauc > best_auc:
            best_auc, since = vauc, 0
            best_state = {k: v.detach().cpu().clone()
                          for k, v in model.state_dict().items()}
        else:
            since += 1
        if since >= patience:
            break
    if best_state is not None:
        model.load_state_dict(best_state)
    return model, {"best_val_auc": best_auc, "n_epochs": len(curve)}, curve


@torch.no_grad()
def evaluate_diff(model: DiffAttentionModel, data, device: torch.device
                  ) -> tuple[dict, dict]:
    """Returns (metrics, extras). metrics = 15's dict; extras has logits,
    prob, pred, y, a_snp (+a_chrom), patient_emb — for predictions.csv /
    attention.npz / embeddings.npz."""
    model.eval()
    d, m, c, b, ab, y = [t.to(device) for t in data]
    logits = model(d, m, c, b, ab)
    attn = {k: v.detach().cpu().numpy()
            for k, v in (model.last_attn or {}).items()}
    emb = model.last_emb.detach().cpu().numpy()
    lg = logits.cpu().numpy()
    yn = y.cpu().numpy()
    prob = 1.0 / (1.0 + np.exp(-lg))
    extras = {"logits": lg, "prob": prob, "pred": (lg > 0).astype(int),
              "y": yn, "patient_emb": emb, **attn}
    return _metrics(yn, lg), extras
