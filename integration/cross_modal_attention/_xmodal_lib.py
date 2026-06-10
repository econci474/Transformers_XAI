"""
_xmodal_lib.py
==============
Shared library for the cross-modal *feature*-fusion experiment: predict the 12-month (m12)
three-way diagnosis (CN/MCI/AD) by fusing a CLINICAL embedding (BioClinical-ModernBERT-large,
T2, full_ft) with an MRI embedding (m12 scan) using either a cross-modal attention block or a
plain MLP-concat head, optionally with SigLIP auxiliary losses.

Two clinical arms:
  - BL   : single baseline embedding  (encoder_outputs_no_cdr_post_exclusion, T2_multiclass)
  - LONG : {bl,m06,m12} per-visit embeddings from the single T2_long_multiclass model
           (trained on ALL visits, contemporaneous label; one shared embedding space)

Two MRI variants:
  - A : BrainDINO T2 (768-d)        mri_embeddings/T2_multiclass/braindino_frozen_none_cached
  - B : BrainMVP ft/stochastic T1b (512-d)  brainmvp_embeddings/T1b_binary/aug_stochastic

Target = the m12 concurrent label (taken from the CL_long m12 row; asserted == MRI m12 label).
Cohort = COMMON INTERSECTION: patients with bl+m06+m12 clinical AND an m12 MRI scan, so the BL
and LONG arms are compared on identical patients. Canonical per-seed splits (clinical split is
authoritative; MRI split asserted to match).

Embeddings are precomputed, so everything trains on CPU in seconds (env: snp, torch 2.12 cpu).

Used by 01_train_xmodal_fusion.py.
"""
import json
import os

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, confusion_matrix,
    f1_score, roc_auc_score,
)

# ── Paths / constants ──────────────────────────────────────────────────────────
# Local Windows default; override on HPC via env vars. The CSD3 clinical and MRI derivative trees
# live under DIFFERENT hpc-work roots, so they are configured separately (see the SL3-GPU submit
# script 02c_run_long_mamba_submit_csd3.sh). os.path.join works for both Windows and POSIX roots.
CL_ROOT = os.environ.get("XMODAL_CL_ROOT", r"D:\ADNI_BIDS_project\derivatives")
MRI_ROOT = os.environ.get("XMODAL_MRI_ROOT", r"D:\ADNI_BIDS_project\derivatives")
CL_BL_TMPL = os.path.join(
    CL_ROOT, "encoder_outputs_no_cdr_post_exclusion", "BioClinical-ModernBERT-large",
    "T2_multiclass", "seed_{s}", "full_ft", "embeddings", "embeddings.parquet")
CL_LONG_TMPL = os.path.join(
    CL_ROOT, "encoder_outputs_no_cdr_post_exclusion_longitudinal", "BioClinical-ModernBERT-large",
    "T2_long_multiclass", "seed_{s}", "full_ft", "embeddings", "embeddings.parquet")
CL_M12_TMPL = os.path.join(
    CL_ROOT, "encoder_outputs_no_cdr_post_exclusion_m12", "BioClinical-ModernBERT-large",
    "T2_m12_multiclass", "seed_{s}", "full_ft", "embeddings", "embeddings.parquet")
MRI_TMPL = {
    "A": os.path.join(MRI_ROOT, "mri_embeddings", "T2_multiclass",
                      "braindino_frozen_none_cached", "embeddings_seed_{s}.npz"),
    "B": os.path.join(MRI_ROOT, "brainmvp_embeddings", "T1b_binary", "aug_stochastic",
                      "seed_{s}", "embeddings_seed_{s}.npz"),
}
MRI_DIM = {"A": 768, "B": 512}
MRI_NAME = {"A": "BrainDINO-T2", "B": "BrainMVP-T1b"}
CLIN_DIM = 1024
N_CLASSES = 3
CLASS_NAMES = ["CN", "MCI", "AD"]
LONG_VISITS = ["bl", "m06", "m12"]
TYPE_IDX = {"bl": 0, "m06": 1, "m12": 2, "mri": 3}    # learned token-type embeddings
EMB_PREFIX = "emb_pool_"


# ── Data loading ────────────────────────────────────────────────────────────────
def _emb_cols(df):
    return sorted([c for c in df.columns if c.startswith(EMB_PREFIX)])


def _load_clin_bl(seed):
    df = pd.read_parquet(CL_BL_TMPL.format(s=seed))
    cols = _emb_cols(df)
    out = df[["Patient_ID", "split", "label"]].copy()
    out["emb"] = list(df[cols].to_numpy(np.float32))
    return out.set_index("Patient_ID")


def _load_clin_long(seed):
    """Return {visit: DataFrame indexed by Patient_ID with emb/split/label} for bl,m06,m12."""
    df = pd.read_parquet(CL_LONG_TMPL.format(s=seed))
    cols = _emb_cols(df)
    df = df[df["VISCODE_2"].isin(LONG_VISITS)].copy()
    by_visit = {}
    for v in LONG_VISITS:
        sub = df[df["VISCODE_2"] == v].drop_duplicates("Patient_ID")
        rec = sub[["Patient_ID", "split", "label"]].copy()
        rec["emb"] = list(sub[cols].to_numpy(np.float32))
        by_visit[v] = rec.set_index("Patient_ID")
    return by_visit


def _load_clin_m12(seed):
    """The dedicated 1-year model (T2_m12_multiclass): one m12 embedding per patient."""
    df = pd.read_parquet(CL_M12_TMPL.format(s=seed))
    cols = _emb_cols(df)
    out = df[["Patient_ID", "split", "label"]].copy()
    out["emb"] = list(df[cols].to_numpy(np.float32))
    return out.set_index("Patient_ID")


def _load_mri(variant, seed):
    z = np.load(MRI_TMPL[variant].format(s=seed), allow_pickle=True)
    df = pd.DataFrame({
        "Patient_ID": z["Patient_ID"].astype(str),
        "VISCODE_2": z["VISCODE_2"].astype(str),
        "split": z["split"].astype(str),
        "label": z["label"].astype(float),
    })
    E = z["embeddings"].astype(np.float32)
    df["emb"] = list(E)
    df = df[df["VISCODE_2"] == "m12"].drop_duplicates("Patient_ID")
    # Drop failed extractions (a single all-NaN m12 scan exists in both variants/all seeds).
    keep = ~df["emb"].apply(lambda v: bool(np.isnan(v).any()))
    df = df[keep]
    return df.set_index("Patient_ID")


def build_cohort(seed, clin_arm, variant, clin_only=False):
    """Assemble the common-intersection cohort for (seed, clin_arm, variant).

    Returns dict: pids, split (array), y (int array, m12 label), clin {visit: (N,1024)},
    mri (N, Dmri) or None when clin_only, and an integrity report.
    """
    cl_bl = _load_clin_bl(seed)
    cl_long = _load_clin_long(seed)
    mri = _load_mri(variant, seed)

    cl_m12 = _load_clin_m12(seed) if clin_arm == "M12" else None

    # Common intersection: need bl+m06+m12 clinical AND (always) the m12 MRI scan, so BL and
    # LONG arms share identical patients (n=478). m12 label (target) comes from CL_long m12 row.
    # The M12 ("1-year model") arm additionally requires the dedicated T2_m12_multiclass embedding
    # (covers 451/478 — 27 patients lack a dedicated-m12 embedding), so its cohort is a subset.
    pids = set(cl_bl.index) & set(mri.index)
    for v in LONG_VISITS:
        pids &= set(cl_long[v].index)
    if clin_arm == "M12":
        pids &= set(cl_m12.index)
    pids = sorted(pids)
    if not pids:
        raise RuntimeError(f"empty cohort for seed={seed} variant={variant}")

    tgt = cl_long["m12"].loc[pids]
    y = tgt["label"].to_numpy(int)
    split = cl_bl.loc[pids, "split"].to_numpy(str)

    # ── integrity checks ──
    rep = {}
    rep["split_mismatch_mri"] = int((mri.loc[pids, "split"].to_numpy(str) != split).sum())
    for v in LONG_VISITS:
        rep[f"split_mismatch_clin_{v}"] = int(
            (cl_long[v].loc[pids, "split"].to_numpy(str) != split).sum())
    # Only variant A's MRI carries the T2 3-way label; variant B's MRI label is the T1b binary
    # label, so this concurrency check is only meaningful (and asserted) for A. Target y always
    # comes from CL_long m12 regardless of variant.
    rep["label_mismatch_mri_vs_clin_m12"] = (
        int((mri.loc[pids, "label"].to_numpy(int) != y).sum()) if variant == "A" else "n/a (T1b)")
    if clin_arm == "M12":
        rep["split_mismatch_clin_m12model"] = int(
            (cl_m12.loc[pids, "split"].to_numpy(str) != split).sum())
    rep["n"] = len(pids)

    clin = {}
    if clin_arm == "LONG":
        visits = LONG_VISITS
        for v in LONG_VISITS:
            clin[v] = np.stack(cl_long[v].loc[pids, "emb"].to_numpy())
    elif clin_arm == "M12":
        visits = ["m12"]                       # single m12 token from the dedicated 1-year model
        clin["m12"] = np.stack(cl_m12.loc[pids, "emb"].to_numpy())
    else:
        visits = ["bl"]
        clin["bl"] = np.stack(cl_bl.loc[pids, "emb"].to_numpy())

    mri_arr = None if clin_only else np.stack(mri.loc[pids, "emb"].to_numpy())
    return {
        "pids": pids, "split": split, "y": y, "clin": clin, "visits": visits,
        "mri": mri_arr, "variant": variant, "clin_arm": clin_arm,
        "clin_only": clin_only, "report": rep,
    }


# ── Z-score (train-fold stats), tensor packing ──────────────────────────────────
def _zfit(x):
    mu = x.mean(0, keepdims=True)
    sd = x.std(0, keepdims=True)
    sd[sd < 1e-6] = 1.0
    return mu, sd


def pack_tensors(data):
    """Z-score each token source on the TRAIN fold and return per-split torch tensors."""
    tr = data["split"] == "train"
    stats = {}
    clin_t = {}
    for v, arr in data["clin"].items():
        mu, sd = _zfit(arr[tr]); stats[("clin", v)] = (mu, sd)
        clin_t[v] = torch.tensor((arr - mu) / sd, dtype=torch.float32)
    mri_t = None
    if not data["clin_only"]:
        mu, sd = _zfit(data["mri"][tr]); stats["mri"] = (mu, sd)
        mri_t = torch.tensor((data["mri"] - mu) / sd, dtype=torch.float32)
    y_t = torch.tensor(data["y"], dtype=torch.long)
    masks = {sp: torch.tensor(data["split"] == sp) for sp in ("train", "val", "test")}
    return clin_t, mri_t, y_t, masks, stats


# ── Model ───────────────────────────────────────────────────────────────────────
def _make_mamba_encoder():
    """Lazily build the frozen-Mamba pooler from the T4 mamba_long_clinical lib.

    Reduces the 3 visit CLS vectors (bl/m06/m12, 1024-d) to a single clinical vector via a frozen
    `state-spaces/mamba-130m-hf` backbone (requires_grad=False) with a trainable projection,
    last-present-step pooling -> 768-d. REQUIRES `transformers` (+ HF mamba weights), so this is
    only importable in the `clinical` env — never the `snp` env (Phase 2 only). Imported lazily so
    Phase-1 mean-pool / mlp_concat / cross_attn runs never touch transformers.
    """
    import sys
    t4 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "T4", "mamba_long_clinical")
    if t4 not in sys.path:
        sys.path.insert(0, t4)
    import _mamba_seq_lib as M  # noqa: E402  (transformers imported inside make_encoder)
    return M.make_encoder(backbone="mamba1", frozen=True, emb_dim=CLIN_DIM, d_model=768)


class FusionModel(nn.Module):
    def __init__(self, mri_dim, visits, arch, clin_only, d=256, dropout=0.1, clin_pool="mean"):
        super().__init__()
        self.visits, self.arch, self.clin_only, self.d = visits, arch, clin_only, d
        self.clin_pool = clin_pool
        self.proj_clin = nn.Sequential(nn.Linear(CLIN_DIM, d), nn.LayerNorm(d), nn.GELU())
        if not clin_only:
            self.proj_mri = nn.Sequential(nn.Linear(mri_dim, d), nn.LayerNorm(d), nn.GELU())
        self.type_emb = nn.Embedding(4, d)
        n_tok = len(visits) + (0 if clin_only else 1)
        if arch == "mlp_concat":
            self.head1 = nn.Sequential(nn.Linear(n_tok * d, d), nn.GELU(), nn.Dropout(dropout))
            self.head2 = nn.Linear(d, N_CLASSES)
        elif arch == "cross_attn":
            layer = nn.TransformerEncoderLayer(d_model=d, nhead=4, dim_feedforward=512,
                                               dropout=dropout, batch_first=True)
            self.encoder = nn.TransformerEncoder(layer, num_layers=1)
            self.head2 = nn.Linear(d, N_CLASSES)
        elif arch == "cross_attn_x":
            # True bidirectional cross-attention over a SINGLE CLS vector per modality
            # (BidirectionalCrossAttentionCLS): clinical<-MRI and MRI<-clinical, residual+LN, concat
            # the two updated CLS -> classifier. No mean-pool over tokens (there is one token each).
            if clin_only:
                raise ValueError("cross_attn_x requires the MRI modality (clin_only not supported)")
            self.mamba = None
            clin_in = CLIN_DIM
            if clin_pool == "mamba":                 # Phase-2: frozen-Mamba pooling of 3 visit CLS
                self.mamba = _make_mamba_encoder()
                clin_in = self.mamba.d_out           # 768
                self.register_buffer(
                    "_mamba_time", torch.tensor([0.0, 0.5, 1.0])[: len(visits)], persistent=False)
            self.clin_proj_x = nn.Sequential(nn.Linear(clin_in, d), nn.LayerNorm(d),
                                             nn.GELU(), nn.Dropout(dropout))
            self.mri_proj_x = nn.Sequential(nn.Linear(mri_dim, d), nn.LayerNorm(d),
                                            nn.GELU(), nn.Dropout(dropout))
            self.c2m = nn.MultiheadAttention(d, num_heads=4, dropout=dropout, batch_first=True)
            self.m2c = nn.MultiheadAttention(d, num_heads=4, dropout=dropout, batch_first=True)
            self.clin_norm = nn.LayerNorm(d); self.mri_norm = nn.LayerNorm(d)
            self.classifier = nn.Sequential(nn.Linear(2 * d, d), nn.GELU(), nn.Dropout(dropout),
                                            nn.Linear(d, N_CLASSES))
        else:
            raise ValueError(arch)
        # SigLIP learnable temperature (log space) + bias, per loss term
        self.log_t_pat = nn.Parameter(torch.tensor(0.0)); self.b_pat = nn.Parameter(torch.tensor(0.0))
        self.log_t_lab = nn.Parameter(torch.tensor(0.0)); self.b_lab = nn.Parameter(torch.tensor(0.0))

    def _forward_xattn(self, clin_t, mri_t):
        """Single-CLS bidirectional cross-attention. Returns logits, fused, c_cls, m_cls
        (c_cls/m_cls = pre-attention projected CLS, used by the patient-SigLIP loss)."""
        if self.clin_pool == "mamba":
            X = torch.stack([clin_t[v] for v in self.visits], 1)        # (N,Lc,1024) temporal order
            n, Lc = X.size(0), len(self.visits)
            # The sequential CPU Mamba scan (no fast kernels) segfaults above ~256 rows/call, so
            # chunk to keep each call (and its backward) small; torch.cat stays differentiable.
            CH = 128
            outs = []
            for i in range(0, n, CH):
                xb = X[i:i + CH]; nb = xb.size(0)
                Tt = self._mamba_time.unsqueeze(0).expand(nb, -1)      # (nb,Lc) months/12
                Msk = torch.ones(nb, Lc, device=xb.device)            # all 3 visits present in cohort
                Ln = torch.full((nb,), Lc, dtype=torch.long, device=xb.device)
                zb, _ = self.mamba(xb, Tt, Msk, Ln)                    # (nb,768) last-present step
                outs.append(zb)
            c_raw = torch.cat(outs, 0)                                  # (N,768)
        else:                                                          # mean of the visit CLS (1 for BL/M12)
            c_raw = torch.stack([clin_t[v] for v in self.visits], 0).mean(0)   # (N,1024)
        c = self.clin_proj_x(c_raw).unsqueeze(1)                       # (N,1,d)
        m = self.mri_proj_x(mri_t).unsqueeze(1)                        # (N,1,d)
        c_cross, _ = self.c2m(c, m, m, need_weights=False)            # clinical <- MRI
        m_cross, _ = self.m2c(m, c, c, need_weights=False)            # MRI <- clinical (original c)
        c_upd = self.clin_norm(c + c_cross)
        m_upd = self.mri_norm(m + m_cross)
        fused = torch.cat([c_upd.squeeze(1), m_upd.squeeze(1)], -1)    # (N,2d) — no token pooling
        return self.classifier(fused), fused, c.squeeze(1), m.squeeze(1)

    def forward(self, clin_t, mri_t):
        if self.arch == "cross_attn_x":
            return self._forward_xattn(clin_t, mri_t)
        toks, clin_proj = [], []
        for v in self.visits:
            p = self.proj_clin(clin_t[v]) + self.type_emb.weight[TYPE_IDX[v]]
            toks.append(p); clin_proj.append(p)
        c_pool = torch.stack(clin_proj, 0).mean(0)          # (B,d) clinical pooled (patient SigLIP)
        m_proj = None
        if not self.clin_only:
            m_proj = self.proj_mri(mri_t) + self.type_emb.weight[TYPE_IDX["mri"]]
            toks.append(m_proj)
        seq = torch.stack(toks, 1)                           # (B,K,d)
        if self.arch == "mlp_concat":
            f = self.head1(seq.reshape(seq.size(0), -1))      # (B,d) fused (label SigLIP)
        else:
            f = self.encoder(seq).mean(1)                     # (B,d)
        return self.head2(f), f, c_pool, m_proj


# ── Losses ──────────────────────────────────────────────────────────────────────
def siglip_cross(a, b, log_t, bias):
    """Cross-modal SigLIP: positives on the diagonal (same patient), negatives off-diagonal."""
    a = F.normalize(a, dim=-1); b = F.normalize(b, dim=-1)
    s = log_t.exp() * (a @ b.t()) + bias
    z = 2.0 * torch.eye(a.size(0), device=a.device) - 1.0
    return -F.logsigmoid(z * s).mean()


def siglip_label(f, y, log_t, bias):
    """Hierarchical SigLIP on the fused embedding: positives = same diagnosis (off-diagonal)."""
    f = F.normalize(f, dim=-1)
    n = f.size(0)
    s = log_t.exp() * (f @ f.t()) + bias
    same = (y[:, None] == y[None, :]).float()
    z = 2.0 * same - 1.0
    off = ~torch.eye(n, dtype=torch.bool, device=f.device)
    return -F.logsigmoid(z * s)[off].mean()


def class_weights(y):
    cnt = np.bincount(y, minlength=N_CLASSES).astype(float)
    cnt[cnt == 0] = 1.0
    w = len(y) / (N_CLASSES * cnt)
    return torch.tensor(w, dtype=torch.float32)


# ── Metrics ─────────────────────────────────────────────────────────────────────
def _safe_auc(y, probs):
    try:
        if len(np.unique(y)) < 2:
            return float("nan")
        return float(roc_auc_score(y, probs, multi_class="ovr", average="macro",
                                   labels=list(range(N_CLASSES))))
    except Exception:
        return float("nan")


def compute_metrics(y, probs):
    pred = probs.argmax(1)
    return {
        "accuracy": round(float(accuracy_score(y, pred)), 4),
        "balanced_acc": round(float(balanced_accuracy_score(y, pred)), 4),
        "auc_roc_ovr": round(_safe_auc(y, probs), 4),
        "macro_f1": round(float(f1_score(y, pred, average="macro",
                                         labels=list(range(N_CLASSES)), zero_division=0)), 4),
        "confusion_matrix": confusion_matrix(
            y, pred, labels=list(range(N_CLASSES))).tolist(),
    }


# ── Train / eval ────────────────────────────────────────────────────────────────
def seed_everything(seed):
    import random
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)


def _subset(clin_t, mri_t, mask):
    cs = {v: x[mask] for v, x in clin_t.items()}
    ms = None if mri_t is None else mri_t[mask]
    return cs, ms


def train_eval(data, arch, loss_cfg, seed, lam=0.01, gam=0.01,
               max_epochs=200, patience=30, lr=1e-3, wd=1e-4,
               warmup=0, min_epochs=0, clin_pool="mean"):
    """Full-batch train on TRAIN, early-stop on VAL weighted-CE, report TEST.

    loss_cfg ∈ {'ce','ce_patient','ce_patient_label'}. Patient/label SigLIP are added only when
    the config requests them and MRI is present (no patient term for clin_only).

    warmup    : first `warmup` epochs train on the AUX (SigLIP) loss ONLY — no CE — so the
                contrastive terms can organize the embedding space before CE dominates. Val-CE
                model-selection and early-stop counting only begin after warmup. (With loss_cfg
                'ce' there is no aux term, so warmup epochs are no-ops.)
    min_epochs: early-stopping (`wait >= patience`) cannot fire before epoch `min_epochs`, so
                later, contrastive-shaped epochs get a chance to win selection.
    """
    seed_everything(seed)
    clin_t, mri_t, y_t, masks, _ = pack_tensors(data)
    use_pat = (not data["clin_only"]) and loss_cfg in ("ce_patient", "ce_patient_label")
    use_lab = loss_cfg == "ce_patient_label"
    warmup = min(warmup, max_epochs - 1)

    # Auto-GPU: the frozen-Mamba pooler (clin_pool='mamba') MUST run on CUDA — its sequential CPU
    # scan is pathologically slow and segfaults on long runs. Harmless no-op on CPU for other arches.
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = FusionModel(MRI_DIM[data["variant"]], data["visits"], arch, data["clin_only"],
                        clin_pool=clin_pool).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=wd)
    cw = class_weights(y_t[masks["train"]].numpy()).to(dev)

    clin_t = {v: x.to(dev) for v, x in clin_t.items()}
    mri_t = None if mri_t is None else mri_t.to(dev)
    y_t = y_t.to(dev)
    masks = {k: m.to(dev) for k, m in masks.items()}

    ctr, mtr = _subset(clin_t, mri_t, masks["train"]); ytr = y_t[masks["train"]]
    cva, mva = _subset(clin_t, mri_t, masks["val"]);   yva = y_t[masks["val"]]

    best_val, best_state, best_ep, wait = float("inf"), None, -1, 0
    for ep in range(max_epochs):
        in_warmup = ep < warmup
        model.train(); opt.zero_grad()
        logits, f, c_pool, m_proj = model(ctr, mtr)
        aux = torch.zeros((), dtype=logits.dtype, device=logits.device)
        if use_pat:
            aux = aux + lam * siglip_cross(c_pool, m_proj, model.log_t_pat, model.b_pat)
        if use_lab:
            aux = aux + gam * siglip_label(f, ytr, model.log_t_lab, model.b_lab)
        loss = aux if in_warmup else F.cross_entropy(logits, ytr, weight=cw) + aux
        if loss.requires_grad:                 # skip no-op warmup steps when there is no aux term
            loss.backward(); opt.step()

        model.eval()
        with torch.no_grad():
            vl = F.cross_entropy(model(cva, mva)[0], yva, weight=cw).item()
        if in_warmup:                          # don't select / count during the aux-only warmup
            continue
        if vl < best_val - 1e-5:
            best_val, best_state, best_ep, wait = vl, {k: v.clone() for k, v in model.state_dict().items()}, ep, 0
        else:
            wait += 1
            if wait >= patience and ep >= min_epochs:
                break

    model.load_state_dict(best_state)
    model.eval()
    out = {}
    with torch.no_grad():
        for sp in ("val", "test"):
            cs, ms = _subset(clin_t, mri_t, masks[sp])
            probs = F.softmax(model(cs, ms)[0], dim=1).cpu().numpy()
            out[sp] = (y_t[masks[sp]].cpu().numpy(), probs)
    return out, {"best_epoch": best_ep, "best_val_ce": round(best_val, 4),
                 "warmup": warmup, "min_epochs": min_epochs, "clin_pool": clin_pool,
                 "class_weights": {str(i): round(float(w), 4) for i, w in enumerate(cw)}}


# ── metrics.json writer (matches clinical_pipeline/03 schema) ───────────────────
def write_metrics(out_dir, config, val_pack, test_pack, fit_info, pids_split, preds):
    os.makedirs(out_dir, exist_ok=True)
    blob = {
        "config": {**config, **fit_info, "objective": "m12_3way"},
        "val_metrics": compute_metrics(*val_pack),
        "test_metrics": compute_metrics(*test_pack),
    }
    with open(os.path.join(out_dir, "metrics.json"), "w") as fh:
        json.dump(blob, fh, indent=2)
    preds.to_parquet(os.path.join(out_dir, "predictions.parquet"), index=False)
    return blob
