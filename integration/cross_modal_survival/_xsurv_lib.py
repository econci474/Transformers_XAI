r"""
_xsurv_lib.py — cross-modal SURVIVAL: bidirectional CLS cross-attention (T2-long clinical ⊕ T1d MRI)
=====================================================================================================
Fuse the longitudinal T2 clinical representation (bl/m06/m12 CLS, mean-pooled) with a T1d MRI embedding
(m12 scan) via the bidirectional single-CLS cross-attention (clinical←MRI, MRI←clinical, residual+LN,
concat), then feed the fused vector into the T4 Weibull-piecewise survival head to parametrise
time-to-AD-conversion. Optionally align the two modalities with a patient-SigLIP auxiliary loss before
or during fusion.

Reuse (single source of truth):
  - T4 `_mamba_seq_lib`: survival labels/cohort/splits (`load_master_labels`, `assemble_sequences`,
    task="survival"), the `WeibullPiecewiseHead` (`make_head`), and `seed_everything`/`split_indices`.
  - Export schema MATCHES T4 `01_train_survival_mamba.py` (Patient_ID, split, event, time, Label_T4,
    S_000..S_060 on a 0..15y / 61-pt grid) so T4 `02`-style scoring (`02l.evaluate`) works unchanged.

Three alignment arms (fused only):
  none     : survival-NLL only.
  joint    : survival-NLL + λ·patient-SigLIP throughout.
  prealign : warmup epochs of patient-SigLIP ONLY (organise the cross-modal space), then NLL (+λ·SigLIP).

References on the SAME cohort: clinical-only (mean-pool → head) and MRI-only (T1d → head).

Env: clinical or snp (torch; mean-pooling needs NO transformers — the mamba backbone is never built).
Cohort is tiny (~25 events in val+test) → keep the fusion small (d=256, dropout) and read with caveats.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

HERE = Path(__file__).resolve().parent
T4DIR = HERE.parent / "T4" / "mamba_long_clinical"
sys.path.insert(0, str(T4DIR))
import _mamba_seq_lib as T4  # noqa: E402  (torch-free at import; mamba backbone never built here)

# ── Paths / constants ────────────────────────────────────────────────────────────
CLIN_ROOT = Path(r"D:\ADNI_BIDS_project\derivatives\encoder_outputs_no_cdr_post_exclusion_longitudinal"
                 r"\BioClinical-ModernBERT-large\T2_long_multiclass")
MASTER = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion"
              r"\verbose\longitudinal\master_clinical_verbose.csv")
MRI_ROOT = Path(r"D:\ADNI_BIDS_project\derivatives")

# T1d MRI variants (m12 scan). T1d is a binary sMCI/pMCI task — used ONLY as a frozen MRI embedding.
MRI_TMPL = {
    "mvp": MRI_ROOT / "brainmvp_embeddings" / "T1d_binary" / "aug_plus_original"
           / "seed_{s}" / "embeddings_seed_{s}.npz",
    "dino": MRI_ROOT / "mri_embeddings" / "T1d_binary" / "braindino_frozen_none_cached"
            / "embeddings_seed_{s}.npz",
}
MRI_DIM = {"mvp": 512, "dino": 768}
MRI_NAME = {"mvp": "BrainMVP-T1d", "dino": "BrainDINO-T1d"}

CLIN_DIM = 1024
GRID_TIMES = np.round(np.linspace(0.0, 15.0, 61), 4)     # matches T4 01 export grid


def clin_parquet(seed):
    return CLIN_ROOT / f"seed_{seed}" / "full_ft" / "embeddings" / "embeddings.parquet"


# ── Cohort: T4 survival labels/splits ∩ T1d m12 MRI ───────────────────────────────
def _load_mri(variant, seed):
    z = np.load(str(MRI_TMPL[variant]).format(s=seed), allow_pickle=True)
    df = pd.DataFrame({"Patient_ID": z["Patient_ID"].astype(str),
                       "VISCODE_2": z["VISCODE_2"].astype(str),
                       "split": z["split"].astype(str)})
    df["emb"] = list(z["embeddings"].astype(np.float32))
    df = df[df["VISCODE_2"] == "m12"].drop_duplicates("Patient_ID")
    df = df[~df["emb"].apply(lambda v: bool(np.isnan(v).any()))]    # drop failed extractions
    return df.set_index("Patient_ID")


def build_surv_cohort(seed, variant, master_labels=None):
    """T4 survival sequence (clinical bl/m06/m12) intersected with the T1d m12 MRI, mean-pooled
    clinical → one vector. Returns aligned numpy arrays + an integrity report."""
    if master_labels is None:
        master_labels = T4.load_master_labels(MASTER)
    seq = T4.assemble_sequences(clin_parquet(seed), master_labels, "survival")
    mri = _load_mri(variant, seed)

    keep = np.array([p in mri.index for p in seq["pid"]])
    idx = np.where(keep)[0]
    pid = seq["pid"][idx]
    # masked mean of the present visit CLS → one clinical vector
    X, M = seq["X"][idx], seq["M"][idx]                              # [n,3,1024], [n,3]
    c_pool = (X * M[..., None]).sum(1) / np.clip(M.sum(1, keepdims=True), 1.0, None)
    mri_emb = np.stack(mri.loc[pid, "emb"].to_numpy()).astype(np.float32)
    split = seq["split"][idx]
    rep = {"n": len(pid), "events": int(seq["event"][idx].sum()),
           "split_mismatch_mri": int((mri.loc[pid, "split"].to_numpy(str) != split).sum()),
           "n_clin_total": int(len(seq["pid"])), "n_mri_total": int(len(mri))}
    return {"pid": pid, "split": split, "variant": variant,
            "clin": c_pool.astype(np.float32), "mri": mri_emb, "mri_dim": MRI_DIM[variant],
            "event": seq["event"][idx].astype(np.float32), "time": seq["time"][idx].astype(np.float32),
            "y": seq["y"][idx].astype(np.int64), "report": rep}


def _zfit(x):
    mu = x.mean(0, keepdims=True); sd = x.std(0, keepdims=True); sd[sd < 1e-6] = 1.0
    return mu, sd


# ── Model: clinical / MRI projections + bidirectional cross-attention ─────────────
def siglip_cross(a, b, log_t, bias):
    """Patient-level cross-modal SigLIP: positives on the diagonal (same patient)."""
    a = F.normalize(a, dim=-1); b = F.normalize(b, dim=-1)
    s = log_t.exp() * (a @ b.t()) + bias
    z = 2.0 * torch.eye(a.size(0), device=a.device) - 1.0
    return -F.logsigmoid(z * s).mean()


class FusedSurv(nn.Module):
    """mode ∈ {fused, clinical, mri}. fused → bidirectional CLS cross-attention → concat (2d);
    clinical/mri → single projected CLS (d). Emits (emb, c_cls, m_cls) for the survival head + SigLIP."""
    def __init__(self, mode, mri_dim, d=256, dropout=0.1):
        super().__init__()
        self.mode = mode
        if mode in ("fused", "clinical"):
            self.clin_proj = nn.Sequential(nn.Linear(CLIN_DIM, d), nn.LayerNorm(d), nn.GELU(),
                                           nn.Dropout(dropout))
        if mode in ("fused", "mri"):
            self.mri_proj = nn.Sequential(nn.Linear(mri_dim, d), nn.LayerNorm(d), nn.GELU(),
                                          nn.Dropout(dropout))
        if mode == "fused":
            self.c2m = nn.MultiheadAttention(d, num_heads=4, dropout=dropout, batch_first=True)
            self.m2c = nn.MultiheadAttention(d, num_heads=4, dropout=dropout, batch_first=True)
            self.clin_norm = nn.LayerNorm(d); self.mri_norm = nn.LayerNorm(d)
            self.out_dim = 2 * d
        else:
            self.out_dim = d
        self.log_t_pat = nn.Parameter(torch.tensor(0.0)); self.b_pat = nn.Parameter(torch.tensor(0.0))

    def forward(self, clin_t, mri_t):
        if self.mode == "clinical":
            c = self.clin_proj(clin_t)
            return c, c, None
        if self.mode == "mri":
            m = self.mri_proj(mri_t)
            return m, None, m
        c = self.clin_proj(clin_t).unsqueeze(1)                     # (N,1,d)
        m = self.mri_proj(mri_t).unsqueeze(1)
        c_cross, _ = self.c2m(c, m, m, need_weights=False)          # clinical ← MRI
        m_cross, _ = self.m2c(m, c, c, need_weights=False)          # MRI ← clinical
        c_upd = self.clin_norm(c + c_cross); m_upd = self.mri_norm(m + m_cross)
        fused = torch.cat([c_upd.squeeze(1), m_upd.squeeze(1)], -1)  # (N,2d) — no token pooling
        return fused, c.squeeze(1), m.squeeze(1)


# ── Train one (mode, align) on one cohort/seed; export survival curves (T4 schema) ─
def train_export(cohort, mode, align, seed, out_dir, *, head_kind="weibull_piecewise",
                 d=256, lam=0.5, warmup=20, lr=1e-3, wd=1e-4, max_epochs=300, patience=30,
                 device="cpu"):
    """align ∈ {none, joint, prealign} (fused only; ignored for clinical/mri)."""
    T4.seed_everything(seed)
    dev = torch.device(device)
    split = cohort["split"]
    tr = split == "train"

    # z-score each modality on the TRAIN fold
    cmu, csd = _zfit(cohort["clin"][tr]); clin = (cohort["clin"] - cmu) / csd
    mmu, msd = _zfit(cohort["mri"][tr]);  mri = (cohort["mri"] - mmu) / msd
    clin_t = torch.tensor(clin, dtype=torch.float32, device=dev)
    mri_t = torch.tensor(mri, dtype=torch.float32, device=dev)
    t_all = torch.tensor(cohort["time"], dtype=torch.float32, device=dev)
    e_all = torch.tensor(cohort["event"], dtype=torch.float32, device=dev)
    masks = {sp: torch.tensor(split == sp, device=dev) for sp in ("train", "val", "test")}

    model = FusedSurv(mode, cohort["mri_dim"], d=d).to(dev)
    head = T4.make_head(head_kind, model.out_dim).to(dev)
    opt = torch.optim.AdamW([q for q in list(model.parameters()) + list(head.parameters())
                             if q.requires_grad], lr=lr, weight_decay=wd)

    use_aux = (mode == "fused") and (align in ("joint", "prealign"))
    warmup = warmup if (mode == "fused" and align == "prealign") else 0

    def fwd(mask):
        emb, c_cls, m_cls = model(clin_t[mask], mri_t[mask])
        return emb, c_cls, m_cls

    best_val, best_state, wait = float("inf"), None, 0
    for ep in range(max_epochs):
        in_warm = ep < warmup
        model.train(); head.train(); opt.zero_grad()
        emb, c_cls, m_cls = fwd(masks["train"])
        nll = head.nll(emb, t_all[masks["train"]], e_all[masks["train"]])
        aux = siglip_cross(c_cls, m_cls, model.log_t_pat, model.b_pat) if use_aux else emb.new_zeros(())
        loss = lam * aux if in_warm else (nll + (lam * aux if use_aux else 0.0))
        if loss.requires_grad:
            loss.backward(); opt.step()

        model.eval(); head.eval()
        with torch.no_grad():
            ev, _, _ = fwd(masks["val"])
            vloss = head.nll(ev, t_all[masks["val"]], e_all[masks["val"]]).item()
        if in_warm:
            continue
        if vloss < best_val - 1e-5:
            best_val, wait = vloss, 0
            best_state = {k: v.detach().cpu().clone() for k, v in
                          list(model.state_dict().items()) + [(f"H.{k}", v) for k, v in head.state_dict().items()]}
        else:
            wait += 1
            if wait >= patience:
                break
    if best_state is not None:
        model.load_state_dict({k: v for k, v in best_state.items() if not k.startswith("H.")})
        head.load_state_dict({k[2:]: v for k, v in best_state.items() if k.startswith("H.")})

    # export per-patient S on GRID_TIMES (all splits — train needed as IPCW reference downstream)
    model.eval(); head.eval()
    with torch.no_grad():
        emb_all, _, _ = model(clin_t, mri_t)
        S = head.survival(emb_all, torch.tensor(GRID_TIMES, dtype=torch.float32, device=dev)).cpu().numpy()
    df = pd.DataFrame({"Patient_ID": cohort["pid"], "split": cohort["split"],
                       "event": cohort["event"], "time": cohort["time"],
                       "Label_T4": np.where(cohort["y"] >= 0, cohort["y"], np.nan)})
    for j in range(len(GRID_TIMES)):
        df[f"S_{j:03d}"] = S[:, j]
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_dir / "predictions.parquet", index=False)
    json.dump({"variant": cohort["variant"], "mode": mode, "align": align, "seed": seed,
               "head": head_kind, "d": d, "lam": lam, "warmup": warmup, "lr": lr,
               "out_dim": model.out_dim, "val_nll": float(best_val),
               "grid_times": GRID_TIMES.tolist(), "report": cohort["report"]},
              open(out_dir / "meta.json", "w"), indent=2)
    return float(best_val), cohort["report"]
