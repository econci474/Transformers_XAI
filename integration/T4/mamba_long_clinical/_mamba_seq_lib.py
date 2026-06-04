r"""
_mamba_seq_lib.py — shared library for the T4 longitudinal-MAMBA sequence model (Task 2)
=========================================================================================
Consumes the single-space per-visit clinical embeddings (Task 1:
encoder_outputs_no_cdr_post_exclusion_longitudinal/.../embeddings.parquet) and builds a short
per-patient trajectory (bl/m06/m12) fed into a configurable sequence encoder with a frozen/finetuned
MAMBA backbone, plus swappable heads (classifier; survival: Weibull-AFT / Weibull-piecewise / Cox-PH /
SODEN) and an optional SigLIP-style patient-alignment auxiliary loss.

Design notes
------------
* Sequences are PACKED contiguous in temporal order (present visits only, length L∈{1,2,3}), padded at the
  END to 3 and masked; the per-visit time scalar t=months/12 encodes the actual visit time. We pool the
  LAST PRESENT step (index L-1) — for a causal SSM the trailing pad steps don't influence it.
* The projection is ALWAYS trained (random init; it is the 1024+1 -> d_model adapter). Only the MAMBA
  backbone is frozen-vs-finetuned (ablation).
* The Weibull-piecewise hazard mirrors 02l._h0/_H0 exactly (numpy-testable).
* Env: clinical (torch). SODEN needs torchdiffeq (guarded). No torch on the local box — math helpers are
  written so the Weibull/Cox pieces can be unit-tested in numpy (see _ref_* mirrors at bottom).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

# torch imported lazily so the pure-pandas assembly + numpy reference math import torch-free.
EPS = 1e-2
VISIT_ORDER = {"bl": 0.0, "m06": 6.0, "m12": 12.0}
CONVERTER_GROUPS = {"pCN_to_AD", "pMCI"}


# ─────────────────────────────────────────────────────────────────────────────
# 1. Sequence assembly (torch-free)
# ─────────────────────────────────────────────────────────────────────────────
def load_master_labels(master_csv):
    """Patient-level survival + T4 labels from the longitudinal master (one row/patient)."""
    m = pd.read_csv(master_csv, low_memory=False,
                    usecols=["Patient_ID", "conversion_group", "bl_dx", "years_to_AD", "FU_years"]) \
        .drop_duplicates("Patient_ID")
    m["Patient_ID"] = m["Patient_ID"].astype(str)
    ev = m["conversion_group"].isin(CONVERTER_GROUPS)
    t = np.where(ev, pd.to_numeric(m["years_to_AD"], errors="coerce"),
                 pd.to_numeric(m["FU_years"], errors="coerce"))
    m["surv_event"] = ev.astype(int).values
    m["surv_time"] = pd.Series(t).fillna(EPS).clip(lower=EPS).values
    m["surv_eligible"] = m["bl_dx"].isin(["CN", "MCI"]).values     # 02l cohort
    yA = pd.to_numeric(m["years_to_AD"], errors="coerce")
    lab = np.where(ev, np.select([yA < 3, yA < 7], [0, 1], default=2), np.nan)
    m["Label_T4"] = lab                                            # converters only, else NaN
    return m.set_index("Patient_ID")


def assemble_sequences(emb_parquet, master_labels, task):
    """Build per-patient padded sequences from one seed's embeddings parquet.

    task: 'survival' (cohort = bl_dx CN/MCI) or 'classifier' (cohort = converters w/ Label_T4).
    Returns a dict of numpy arrays:
      X [N,3,D] emb_pool, T [N,3] time(months/12), M [N,3] presence mask, L [N] lengths,
      pid [N], split [N], and label fields (survival: event,time ; classifier: y).
    """
    df = pd.read_parquet(emb_parquet)
    df["Patient_ID"] = df["Patient_ID"].astype(str)
    pool_cols = sorted(c for c in df.columns if c.startswith("emb_pool_"))
    D = len(pool_cols)
    rows, X, Tt, M, L, pid, split, ev, tt, y = [], [], [], [], [], [], [], [], [], []
    for p, g in df.groupby("Patient_ID"):
        if p not in master_labels.index:
            continue
        ml = master_labels.loc[p]
        if task == "survival" and not bool(ml["surv_eligible"]):
            continue
        if task == "classifier" and not np.isfinite(ml["Label_T4"]):
            continue
        g = g[g["VISCODE_2"].isin(VISIT_ORDER)].copy()
        if len(g) == 0:
            continue
        g["__ord"] = g["VISCODE_2"].map(VISIT_ORDER)
        g = g.sort_values("__ord")
        emb = g[pool_cols].to_numpy(np.float32)               # [L,D] contiguous, temporal order
        times = (g["__ord"].to_numpy(np.float32)) / 12.0       # months/12 -> {0,0.5,1.0}
        Lp = len(g)
        xp = np.zeros((3, D), np.float32); xp[:Lp] = emb
        tp = np.zeros((3,), np.float32); tp[:Lp] = times
        mp = np.zeros((3,), np.float32); mp[:Lp] = 1.0
        X.append(xp); Tt.append(tp); M.append(mp); L.append(Lp)
        pid.append(p); split.append(str(g["split"].iloc[0]))
        ev.append(int(ml["surv_event"])); tt.append(float(ml["surv_time"]))
        y.append(int(ml["Label_T4"]) if np.isfinite(ml["Label_T4"]) else -1)
    out = dict(X=np.stack(X), T=np.stack(Tt), M=np.stack(M), L=np.array(L, np.int64),
               pid=np.array(pid), split=np.array(split), emb_dim=D,
               event=np.array(ev, np.float32), time=np.array(tt, np.float32),
               y=np.array(y, np.int64))
    return out


def split_indices(seq, sp):
    return np.where(seq["split"] == sp)[0]


# ─────────────────────────────────────────────────────────────────────────────
# 2. Reference Weibull-piecewise math (numpy; mirrors 02l._h0/_H0) — unit-testable
# ─────────────────────────────────────────────────────────────────────────────
def ref_H0(p, g1, T, t, g2=None):
    g2 = (p * g1 * T ** (p - 1.0)) if g2 is None else g2
    t = np.asarray(t, float)
    return g1 * np.minimum(t, T) ** p + g2 * np.maximum(t - T, 0.0)


def ref_h0(p, g1, T, t, g2=None):
    g2 = (p * g1 * T ** (p - 1.0)) if g2 is None else g2
    t = np.asarray(t, float)
    return np.where(t < T, g1 * p * np.clip(t, 1e-8, None) ** (p - 1.0), g2)


def ref_nll_piecewise(p, g1, T, t, e, g2=None):
    """Per-sample-mean negative log-likelihood (scalar) for a parametric (non-PH) Weibull-CP."""
    h = np.clip(ref_h0(p, g1, T, t, g2), 1e-300, None)
    H = ref_H0(p, g1, T, t, g2)
    return -np.mean(np.asarray(e) * np.log(h) - H)


# ═════════════════════════════════════════════════════════════════════════════
# Everything below needs torch (imported lazily inside the factory).
# ═════════════════════════════════════════════════════════════════════════════
def _torch():
    import torch  # noqa
    return torch


def seed_everything(seed):
    import os, random
    torch = _torch()
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)


# ── Backbone factory ─────────────────────────────────────────────────────────
MAMBA1_ID = "state-spaces/mamba-130m-hf"
MAMBA2_ID = "AntonV/mamba2-130m-hf"


def build_backbone(kind, d_model, frozen, hf_cache=None):
    """Return (module, d_out, takes_inputs_embeds). kind ∈ {mamba1,mamba2,gru,meanpool}."""
    torch = _torch()
    import torch.nn as nn
    if kind in ("mamba1", "mamba2"):
        from transformers import AutoModel
        mid = MAMBA1_ID if kind == "mamba1" else MAMBA2_ID
        backbone = AutoModel.from_pretrained(mid, cache_dir=hf_cache)
        d_out = backbone.config.hidden_size
        if frozen:
            for q in backbone.parameters():
                q.requires_grad = False
        return backbone, d_out, True
    if kind == "gru":
        return nn.GRU(d_model, d_model, batch_first=True), d_model, False
    if kind == "meanpool":
        return nn.Identity(), d_model, False
    raise ValueError(kind)


# ── Sequence encoder ─────────────────────────────────────────────────────────
def make_encoder(*, backbone="mamba1", frozen=True, time_concat="append", emb_dim=1024,
                 d_model=768, proj_hidden=True, hf_cache=None):
    torch = _torch()
    import torch.nn as nn
    import torch.nn.functional as F

    class SeqEncoder(nn.Module):
        def __init__(self):
            super().__init__()
            self.time_concat = time_concat
            self.backbone_kind = backbone
            in_dim = emb_dim + 1
            layers = [nn.LayerNorm(in_dim), nn.Linear(in_dim, d_model)]
            if proj_hidden:
                layers += [nn.GELU(), nn.Linear(d_model, d_model)]
            self.proj = nn.Sequential(*layers)
            self.backbone, self.d_out, self.embeds = build_backbone(backbone, d_model, frozen, hf_cache)
            self.frozen = frozen
            # learnable SigLIP temperature/bias for the patient-align aux loss
            self.align_t = nn.Parameter(torch.tensor(10.0))
            self.align_b = nn.Parameter(torch.tensor(-10.0))

        def project(self, X, Tt):
            t = Tt.unsqueeze(-1)                                  # [N,3,1]
            x = torch.cat([X, t], -1) if self.time_concat == "append" else torch.cat([t, X], -1)
            return self.proj(x)                                   # [N,3,d_model]

        def forward(self, X, Tt, M, L):
            h = self.project(X, Tt)                               # [N,3,d_model]  (also used by align loss)
            if self.backbone_kind in ("mamba1", "mamba2"):
                # NB: a FROZEN backbone is frozen via requires_grad=False on its params (build_backbone),
                # NOT torch.no_grad() here — gradients MUST still flow THROUGH it to the trained
                # projection. Wrapping this in no_grad() would detach the graph and the proj never learns.
                out = self.backbone(inputs_embeds=h).last_hidden_state          # [N,3,d_out]
            elif self.backbone_kind == "gru":
                out, _ = self.backbone(h)
            else:                                                 # meanpool control (no SSM)
                out = h
            if self.backbone_kind == "meanpool":
                m = M.unsqueeze(-1)                                # masked MEAN over present visits
                z = (out * m).sum(1) / m.sum(1).clamp_min(1.0)
            else:                                                 # sequential: last-PRESENT step
                idx = (L - 1).clamp(min=0).view(-1, 1, 1).expand(-1, 1, out.size(-1))
                z = out.gather(1, idx).squeeze(1)                 # [N,d_out]
            return z, h

        def align_loss(self, h, M, pid_int):
            """SigLIP pairwise-sigmoid loss over all PRESENT projected visit embeddings in the batch."""
            N, S, Dm = h.shape
            flat = h.reshape(N * S, Dm)
            mask = M.reshape(-1) > 0.5
            pids = pid_int.view(-1, 1).expand(-1, S).reshape(-1)
            hv, pv = flat[mask], pids[mask]
            if hv.size(0) < 2:
                return hv.new_zeros(())
            hn = F.normalize(hv, dim=-1)
            sim = hn @ hn.t()                                     # cosine [n,n]
            s = self.align_t * sim + self.align_b
            z = torch.where(pv.view(-1, 1) == pv.view(1, -1), 1.0, -1.0)
            return -F.logsigmoid(z * s).mean()

    return SeqEncoder()


def _nullctx():
    import contextlib
    return contextlib.nullcontext()


# ── Heads (classifier + survival parametrisations) ───────────────────────────
def make_head(kind, d_out, *, n_classes=3, free_g2=False):
    torch = _torch()
    import torch.nn as nn
    import torch.nn.functional as F

    sp = F.softplus

    class ClassifierHead(nn.Module):
        kind = "classifier"
        def __init__(self):
            super().__init__(); self.fc = nn.Linear(d_out, n_classes)
        def logits(self, z): return self.fc(z)

    class WeibullAFTHead(nn.Module):
        kind = "weibull_aft"
        def __init__(self):
            super().__init__(); self.fc = nn.Linear(d_out, 2)
        def params(self, z):
            o = self.fc(z); k = sp(o[:, 0]) + 1e-3; lam = sp(o[:, 1]) + 1e-3
            return k, lam
        def nll(self, z, t, e):
            k, lam = self.params(z); t = t.clamp_min(1e-3)
            logh = torch.log(k / lam) + (k - 1) * torch.log(t / lam)
            H = (t / lam) ** k
            return -(e * logh - H).mean()
        def survival(self, z, times):
            k, lam = self.params(z)
            tt = times.view(1, -1)
            return torch.exp(-(tt / lam.view(-1, 1)) ** k.view(-1, 1))

    class WeibullPiecewiseHead(nn.Module):
        kind = "weibull_piecewise"
        def __init__(self):
            super().__init__(); self.fc = nn.Linear(d_out, 4 if free_g2 else 3)
        def params(self, z):
            o = self.fc(z)
            p = sp(o[:, 0]) + 1e-2; g1 = sp(o[:, 1]) + 1e-4; T = sp(o[:, 2]) + 1e-2
            g2 = (sp(o[:, 3]) + 1e-4) if free_g2 else p * g1 * T ** (p - 1.0)
            return p, g1, T, g2
        def _hH(self, z, t):
            p, g1, T, g2 = self.params(z); t = t.clamp_min(1e-6)
            h = torch.where(t < T, g1 * p * t ** (p - 1.0), g2)
            H = g1 * torch.minimum(t, T) ** p + g2 * torch.clamp(t - T, min=0.0)
            return h, H
        def nll(self, z, t, e):
            h, H = self._hH(z, t)
            return -(e * torch.log(h.clamp_min(1e-30)) - H).mean()
        def survival(self, z, times):
            p, g1, T, g2 = self.params(z)
            tt = times.view(1, -1)
            Tn, pn, g1n, g2n = T.view(-1, 1), p.view(-1, 1), g1.view(-1, 1), g2.view(-1, 1)
            H = g1n * torch.minimum(tt, Tn) ** pn + g2n * torch.clamp(tt - Tn, min=0.0)
            return torch.exp(-H.clamp(max=60.0))

    class CoxPHHead(nn.Module):
        """DeepSurv: linear risk on z, Cox partial-likelihood; Breslow baseline -> S at eval."""
        kind = "cox"
        def __init__(self):
            super().__init__(); self.fc = nn.Linear(d_out, 1, bias=False)
            self.register_buffer("bh_t", torch.zeros(0))         # Breslow grid (set at fit)
            self.register_buffer("bh_H0", torch.zeros(0))
        def risk(self, z): return self.fc(z).squeeze(-1)
        def nll(self, z, t, e):                                  # full-batch partial likelihood
            r = self.risk(z)
            order = torch.argsort(t, descending=True)
            r_s, e_s = r[order], e[order]
            log_cum = torch.logcumsumexp(r_s, dim=0)             # risk set R(t_i) (t sorted desc)
            return -((r_s - log_cum) * e_s).sum() / e_s.sum().clamp_min(1.0)
        @torch.no_grad()
        def fit_baseline(self, z, t, e):
            r = self.risk(z); ex = torch.exp(r)
            ut = torch.unique(t[e > 0.5])
            ut, _ = torch.sort(ut)
            H0, cum = [], 0.0
            for tk in ut:
                d = ((t == tk) & (e > 0.5)).sum().float()
                denom = ex[t >= tk].sum().clamp_min(1e-8)
                cum = cum + (d / denom).item(); H0.append(cum)
            self.bh_t = ut.to(z.device); self.bh_H0 = torch.tensor(H0, device=z.device)
        @torch.no_grad()
        def survival(self, z, times):
            r = self.risk(z); ex = torch.exp(r)
            if self.bh_t.numel() == 0:
                return torch.ones(z.size(0), times.numel(), device=z.device)
            H0_at = torch.stack([self.bh_H0[torch.searchsorted(self.bh_t, tt, right=True).clamp(max=self.bh_t.numel()) - 1]
                                 if tt >= self.bh_t[0] else torch.zeros((), device=z.device) for tt in times])
            return torch.exp(-H0_at.view(1, -1) * ex.view(-1, 1))

    heads = {"classifier": ClassifierHead, "weibull_aft": WeibullAFTHead,
             "weibull_piecewise": WeibullPiecewiseHead, "cox": CoxPHHead}
    if kind == "soden":
        return _make_soden_head(d_out)
    return heads[kind]()


def _make_soden_head(d_out):
    """SODEN neural-ODE cumulative-hazard head (guarded: needs torchdiffeq). Adapted from
    jiaqima/SODEN (MIT): H(t|z)=∫_0^t f(s,z) ds with f>0 (monotone H); S=exp(-H); NLL=-(e·log f(t)-H)."""
    torch = _torch()
    import torch.nn as nn
    import torch.nn.functional as F
    try:
        from torchdiffeq import odeint
    except Exception as exc:                                     # guarded
        raise ImportError(f"SODEN head needs torchdiffeq ({exc}); skip the soden config.")

    class SODENHead(nn.Module):
        kind = "soden"
        def __init__(self, hid=64):
            super().__init__()
            self.f = nn.Sequential(nn.Linear(d_out + 1, hid), nn.Tanh(),
                                   nn.Linear(hid, hid), nn.Tanh(), nn.Linear(hid, 1), nn.Softplus())
        def hazard(self, z, t):                                  # f(t,z) > 0  [N]
            return self.f(torch.cat([z, t.view(-1, 1)], -1)).squeeze(-1)
        def cum_hazard(self, z, t, n=16):                        # ∫_0^t f ds via grid quadrature
            ts = torch.linspace(0, 1, n, device=z.device).view(1, -1) * t.view(-1, 1)
            fs = torch.stack([self.hazard(z, ts[:, k]) for k in range(n)], -1)
            return torch.trapz(fs, ts, dim=-1)
        def nll(self, z, t, e):
            f = self.hazard(z, t.clamp_min(1e-4)); H = self.cum_hazard(z, t.clamp_min(1e-4))
            return -(e * torch.log(f.clamp_min(1e-30)) - H).mean()
        def survival(self, z, times):
            S = []
            for tt in times:
                H = self.cum_hazard(z, torch.full((z.size(0),), float(tt), device=z.device))
                S.append(torch.exp(-H.clamp(max=60.0)))
            return torch.stack(S, -1)
    return SODENHead()


# ── Training loop (full-batch; early stop on val task loss) ───────────────────
def train_survival(encoder, head, seq, *, device="cpu", lr=1e-3, max_epochs=300, patience=30,
                   weight_decay=1e-4, align_lambda=0.0, log_every=0):
    torch = _torch()
    enc, hd = encoder.to(device), head.to(device)
    tr, va = split_indices(seq, "train"), split_indices(seq, "val")

    def batch(idx):
        T = _torch()
        to = lambda a, d=torch.float32: T.tensor(a[idx], dtype=d, device=device)
        pid_int = T.tensor(pd.factorize(seq["pid"][idx])[0], dtype=torch.long, device=device)
        return (to(seq["X"]), to(seq["T"]), to(seq["M"]),
                T.tensor(seq["L"][idx], dtype=torch.long, device=device),
                to(seq["time"]), to(seq["event"]), pid_int)

    Xtr, Ttr, Mtr, Ltr, ttr, etr, ptr = batch(tr)
    Xva, Tva, Mva, Lva, tva, eva, pva = batch(va)
    params = [q for q in list(enc.parameters()) + list(hd.parameters()) if q.requires_grad]
    opt = torch.optim.AdamW(params, lr=lr, weight_decay=weight_decay)

    best, best_state, bad = float("inf"), None, 0
    for ep in range(max_epochs):
        enc.train(); hd.train(); opt.zero_grad()
        z, h = enc(Xtr, Ttr, Mtr, Ltr)
        if hd.kind == "cox" and hasattr(hd, "fit_baseline"):
            loss = hd.nll(z, ttr, etr)
        else:
            loss = hd.nll(z, ttr, etr)
        if align_lambda > 0:
            loss = loss + align_lambda * enc.align_loss(h, Mtr, ptr)
        loss.backward(); opt.step()

        enc.eval(); hd.eval()
        with torch.no_grad():
            zv, _ = enc(Xva, Tva, Mva, Lva)
            vloss = hd.nll(zv, tva, eva).item()
        if log_every and ep % log_every == 0:
            print(f"    ep{ep:3d} train {loss.item():.4f} val {vloss:.4f}")
        if vloss < best - 1e-5:
            best, bad = vloss, 0
            best_state = {k: v.detach().cpu().clone() for k, v in
                          list(enc.state_dict().items()) + [(f"head.{k}", v) for k, v in hd.state_dict().items()]}
        else:
            bad += 1
            if bad >= patience:
                break
    # restore best
    if best_state is not None:
        enc.load_state_dict({k: v for k, v in best_state.items() if not k.startswith("head.")})
        hd.load_state_dict({k[5:]: v for k, v in best_state.items() if k.startswith("head.")})
    if hd.kind == "cox" and hasattr(hd, "fit_baseline"):
        enc.eval(); hd.eval()
        with torch.no_grad():
            ztr, _ = enc(Xtr, Ttr, Mtr, Ltr)
        hd.fit_baseline(ztr, ttr, etr)
    return enc, hd, best


def predict_survival_curves(encoder, head, seq, idx, times, device="cpu"):
    """Return [len(idx), len(times)] survival matrix for the given patient indices."""
    torch = _torch()
    enc, hd = encoder.to(device).eval(), head.to(device).eval()
    to = lambda a, d=torch.float32: torch.tensor(a[idx], dtype=d, device=device)
    with torch.no_grad():
        z, _ = enc(to(seq["X"]), to(seq["T"]), to(seq["M"]),
                   torch.tensor(seq["L"][idx], dtype=torch.long, device=device))
        S = hd.survival(z, torch.tensor(times, dtype=torch.float32, device=device))
    return S.cpu().numpy()
