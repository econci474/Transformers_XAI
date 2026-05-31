"""
02j_weibull_changepoint_vs_cox.py
=================================
Piecewise Weibull-exponential CHANGE-POINT hazard model vs a matching Cox PH model,
on the post-exclusion pooled CN + MCI -> AD endpoint (baseline CN+MCI, event = pCN_to_AD
or pMCI), on both the time-from-baseline and chronological-age (left-truncated) axes.

Change-point hazard (continuity-constrained; covariates enter the change-point T):
    h(·) = p·γ₁··^(p-1)         for · < T        (Weibull)
    h(·) = γ₂ = p·γ₁·T^(p-1)    for · ≥ T        (constant — continuous at T)
    log T = b₀ + b·age_z(time only) + Σ_k c_k·X_k
Cox comparator (matching covariates):
    TIME: h(t|x) = h₀(t)·exp(a·age_z + b·age_z·log t + Σ c_k·X_k)   (time-varying age coef)
    AGE : h(a|x) = h₀(a)·exp(Σ c_k·X_k)                              (left-truncated)

The script loops over COVARIATE COMBINATIONS (add APOE2, then Sex) and writes every
combination to its own sub-directory:
    outputs/weibull/<combo>/   and   outputs/cox_ph/<combo>/
Every figure carries the fitted model equations underneath the number-at-risk rows, so
the model being plotted is always self-documented as more covariates are added.

Comparison metrics (per axis): Harrell C-index at a fixed horizon and a LEFT-TRUNCATION-
aware IPCW integrated Brier score (IBS). AIC is reported but flagged: parametric
full-likelihood AIC and Cox partial-likelihood AIC are NOT on the same scale.

Env: clinical (+ lifelines).  Run:  python clinical_pipeline/02j_weibull_changepoint_vs_cox.py
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import minimize

try:                                    # Windows console is cp1252; allow unicode prints
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

from lifelines import CoxPHFitter, CoxTimeVaryingFitter
from lifelines.utils import to_episodic_format, concordance_index

matplotlib.rcParams["font.family"] = "DejaVu Serif"

# ── Paths + covariate combinations ──────────────────────────────────────────────
REPO_ROOT  = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
MASTER_CSV = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion"
                  r"\tabular\longitudinal\master_clinical_tabular.csv")
WEIBULL_BASE = REPO_ROOT / "clinical_pipeline" / "outputs" / "weibull"
COX_BASE     = REPO_ROOT / "clinical_pipeline" / "outputs" / "cox_ph"

# (subdir name, covariate columns, display labels). APOE4 (covcols[0]) is the stratifier.
COMBOS = [
    ("apoe4",            ["apoe"],                 ["APOE4"]),
    ("apoe4_apoe2",      ["apoe", "apoe2"],        ["APOE4", "APOE2"]),
    ("apoe4_apoe2_sex",  ["apoe", "apoe2", "sex_M"], ["APOE4", "APOE2", "Sex(M)"]),
]

EPS = 0.01
BLUE, RED = "#2c7bb6", "#d7191c"


# ── Cohort (pooled CN+MCI -> AD; mirrors 02g/02i) ───────────────────────────────
def load_cohort():
    cols = ["Patient_ID", "VISCODE_long", "Date", "YOB", "APOE4_Dosage", "APOE_Alleles",
            "Sex", "bl_dx", "conversion_group", "FU_years", "years_to_AD"]
    df = pd.read_csv(MASTER_CSV, usecols=cols, low_memory=False)
    df["Date"] = pd.to_datetime(df["Date"], errors="coerce")
    bl = (df[df["VISCODE_long"] == "bl"].dropna(subset=["Date"])
          .sort_values("Date").drop_duplicates("Patient_ID").copy())
    birth = pd.to_datetime(pd.DataFrame(
        {"year": bl["YOB"].astype(int), "month": 7, "day": 1}).set_index(bl.index))
    bl["age_at_baseline"] = (bl["Date"] - birth) / np.timedelta64(1, "D") / 365.25

    c = bl[bl["bl_dx"].isin(["CN", "MCI"])].copy()
    c["event"] = c["conversion_group"].isin(["pCN_to_AD", "pMCI"]).astype(int)
    t = np.where(c["event"] == 1, c["years_to_AD"], c["FU_years"])
    c["t"] = pd.to_numeric(pd.Series(t, index=c.index), errors="coerce").fillna(EPS).clip(lower=EPS)
    c["entry_age"] = c["age_at_baseline"]
    c["exit_age"] = c["entry_age"] + c["t"]
    c["apoe"] = (pd.to_numeric(c["APOE4_Dosage"], errors="coerce") >= 1).astype(float)
    c["apoe2"] = c["APOE_Alleles"].fillna("").apply(
        lambda s: 1.0 if "2" in str(s).split("/") else 0.0)              # APOE2 carrier
    c["sex_M"] = (c["Sex"].astype(str) == "Male").astype(float)
    c["age_z"] = (c["age_at_baseline"] - c["age_at_baseline"].mean()) / c["age_at_baseline"].std(ddof=0)
    return c.dropna(subset=["t", "entry_age", "exit_age", "apoe", "apoe2", "sex_M", "age_z"]).reset_index(drop=True)


# ── Weibull-exponential change-point model (continuity-constrained): direct MLE ──
def _H_changepoint(theta, t, X):
    p = np.exp(theta[0]); g1 = np.exp(theta[1]); beta = theta[2:]
    Tc = np.exp(X @ beta)
    tmin = np.minimum(t, Tc)
    return g1 * tmin ** p + p * g1 * Tc ** (p - 1.0) * np.maximum(t - Tc, 0.0)


def _h_changepoint(theta, t, X):
    p = np.exp(theta[0]); g1 = np.exp(theta[1]); beta = theta[2:]
    Tc = np.exp(X @ beta)
    return np.where(t < Tc, g1 * p * np.clip(t, 1e-8, None) ** (p - 1.0), p * g1 * Tc ** (p - 1.0))


def fit_changepoint(dur, event, X, entry=None, t_bounds=(0.5, None), p_bounds=(0.3, 6.0),
                    n_starts=24):
    """MLE; covariates enter log T via design X (last column = intercept). p and the
    change-point are bounded to keep T inside the data range (else the model collapses
    toward exponential and p runs away)."""
    dur = np.asarray(dur, float); event = np.asarray(event, float); X = np.asarray(X, float)
    entry = None if entry is None else np.asarray(entry, float)
    tot_t, n_ev = dur.sum(), event.sum()
    t_lo, t_hi = t_bounds[0], (t_bounds[1] if t_bounds[1] else 0.9 * dur.max())
    lp_lo, lp_hi = np.log(p_bounds[0]), np.log(p_bounds[1])
    lt_lo, lt_hi = np.log(t_lo), np.log(t_hi)

    def nll(theta):
        if not np.all(np.isfinite(theta)) or not (lp_lo <= theta[0] <= lp_hi):
            return 1e12
        lp_pred = X @ theta[2:]
        if lp_pred.min() < lt_lo or lp_pred.max() > lt_hi:
            return 1e12
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            h = _h_changepoint(theta, dur, X)
            Ht = _H_changepoint(theta, dur, X)
            He = 0.0 if entry is None else _H_changepoint(theta, entry, X)
            ll = np.sum(event * np.log(np.clip(h, 1e-300, None)) - (Ht - He))
        return -ll if np.isfinite(ll) else 1e12

    med = np.median(dur[event == 1]) if n_ev > 0 else np.median(dur)
    base = np.zeros(2 + X.shape[1])
    base[1] = np.log(max(n_ev / max(tot_t, 1e-6), 1e-4))
    base[-1] = np.clip(np.log(max(med, t_lo)), lt_lo, lt_hi)
    rng = np.random.default_rng(20240601)
    best = None
    for s in range(n_starts):
        x0 = base.copy()
        if s:
            x0 += rng.normal(0, 0.4, base.shape); x0[0] = rng.uniform(lp_lo, lp_hi)
            x0[-1] = np.clip(x0[-1], lt_lo, lt_hi)
        try:
            res = minimize(nll, x0, method="Nelder-Mead",
                           options={"maxiter": 12000, "xatol": 1e-7, "fatol": 1e-7})
        except Exception:
            continue
        if np.isfinite(res.fun) and (best is None or res.fun < best.fun):
            best = res
    theta = best.x; k = len(theta); ll = -best.fun
    return {"theta": theta, "k": k, "loglik": ll, "aic": 2 * k - 2 * ll,
            "p": float(np.exp(theta[0])), "g1": float(np.exp(theta[1])), "beta": theta[2:],
            "predict_S": lambda Xs, tt, ent=0.0: np.exp(
                -(_H_changepoint(theta, tt, Xs) - _H_changepoint(theta, ent, Xs)))}


# ── Discontinuous change-point (free γ₂): place the covariate BLOCK on T, γ₁ and/or
#    γ₂ to ask WHERE the covariates act. ──────────────────────────────────────────
def _split_theta(theta, shapes):
    n1, n2, n3 = shapes
    return theta[0], theta[1:1+n1], theta[1+n1:1+n1+n2], theta[1+n1+n2:1+n1+n2+n3]


def _H_disc(theta, t, D, shapes):
    log_p, bg1, bg2, bT = _split_theta(theta, shapes)
    p = np.exp(log_p); g1 = np.exp(D["g1"] @ bg1); g2 = np.exp(D["g2"] @ bg2); Tc = np.exp(D["T"] @ bT)
    return g1 * np.minimum(t, Tc) ** p + g2 * np.maximum(t - Tc, 0.0)


def _h_disc(theta, t, D, shapes):
    log_p, bg1, bg2, bT = _split_theta(theta, shapes)
    p = np.exp(log_p); g1 = np.exp(D["g1"] @ bg1); g2 = np.exp(D["g2"] @ bg2); Tc = np.exp(D["T"] @ bT)
    return np.where(t < Tc, g1 * p * np.clip(t, 1e-8, None) ** (p - 1.0), g2)


def build_designs(c, axis, variant, covcols):
    n = len(c); ones = np.ones(n)
    Xcov = np.column_stack([c[cc].values for cc in covcols])    # (n, kc)
    Tcols = ([c["age_z"].values] if axis == "time" else [])
    if variant in ("T_only", "T_g1", "T_g2", "all"):
        Tcols.append(Xcov)
    Tcols.append(ones)
    XT = np.column_stack(Tcols)
    Xg1 = np.column_stack([Xcov, ones]) if variant in ("T_g1", "all") else ones.reshape(-1, 1)
    Xg2 = np.column_stack([Xcov, ones]) if variant in ("T_g2", "g2_only", "all") else ones.reshape(-1, 1)
    return {"g1": Xg1, "g2": Xg2, "T": XT}


def fit_cp_general(dur, event, D, entry=None, t_bounds=(0.5, None), p_bounds=(0.3, 6.0), n_starts=24):
    dur = np.asarray(dur, float); event = np.asarray(event, float)
    entry = None if entry is None else np.asarray(entry, float)
    n1, n2, n3 = D["g1"].shape[1], D["g2"].shape[1], D["T"].shape[1]
    shapes = (n1, n2, n3)
    tot, nev = dur.sum(), event.sum()
    t_lo, t_hi = t_bounds[0], (t_bounds[1] or 0.9 * dur.max())
    lp_lo, lp_hi = np.log(p_bounds[0]), np.log(p_bounds[1])
    lt_lo, lt_hi = np.log(t_lo), np.log(t_hi)

    def nll(theta):
        if not np.all(np.isfinite(theta)) or not (lp_lo <= theta[0] <= lp_hi):
            return 1e12
        _, _, _, bT = _split_theta(theta, shapes)
        lpT = D["T"] @ bT
        if lpT.min() < lt_lo or lpT.max() > lt_hi:
            return 1e12
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            h = _h_disc(theta, dur, D, shapes); Ht = _H_disc(theta, dur, D, shapes)
            He = 0.0 if entry is None else _H_disc(theta, entry, D, shapes)
            ll = np.sum(event * np.log(np.clip(h, 1e-300, None)) - (Ht - He))
        return -ll if np.isfinite(ll) else 1e12

    med = np.median(dur[event == 1]) if nev > 0 else np.median(dur)
    rate = np.log(max(nev / max(tot, 1e-6), 1e-4))
    base = np.zeros(1 + n1 + n2 + n3)
    base[1 + n1 - 1] = rate; base[1 + n1 + n2 - 1] = rate
    base[-1] = np.clip(np.log(max(med, t_lo)), lt_lo, lt_hi)
    rng = np.random.default_rng(20240601)
    best = None
    for s in range(n_starts):
        x0 = base.copy()
        if s:
            x0 += rng.normal(0, 0.4, base.shape); x0[0] = rng.uniform(lp_lo, lp_hi)
            x0[-1] = np.clip(x0[-1], lt_lo, lt_hi)
        try:
            res = minimize(nll, x0, method="Nelder-Mead",
                           options={"maxiter": 16000, "xatol": 1e-7, "fatol": 1e-7})
        except Exception:
            continue
        if np.isfinite(res.fun) and (best is None or res.fun < best.fun):
            best = res
    theta = best.x; k = len(theta); ll = -best.fun
    return {"theta": theta, "shapes": shapes, "k": k, "loglik": ll, "aic": 2 * k - 2 * ll,
            "predict_S": lambda Dd, tt, ent=0.0: np.exp(
                -(_H_disc(theta, tt, Dd, shapes) - _H_disc(theta, ent, Dd, shapes)))}


VARIANTS = [("T_only", "T only", "T"), ("T_g1", "T + γ₁", "T, γ₁"),
            ("T_g2", "T + γ₂", "T, γ₂"), ("all", "T + γ₁ + γ₂", "T, γ₁, γ₂"),
            ("g2_only", "γ₂ only", "γ₂")]


def fit_variants(c, axis, dur, event, entry, grid, hor, t_bounds, covcols):
    rows = []; models = []
    for key, lab, onstr in VARIANTS:
        D = build_designs(c, axis, key, covcols)
        res = fit_cp_general(dur, event, D, entry=entry, t_bounds=t_bounds)
        ent = 0.0 if entry is None else entry
        Sg = np.column_stack([res["predict_S"](D, g, ent) for g in grid])
        ibs = ipcw_ibs(Sg, dur, event, grid, entry=entry)
        cidx = cindex_at(res["predict_S"](D, hor, ent), dur, event)
        rows.append(dict(axis=axis, model=f"Weibull-cp: {lab}", placement=onstr, n_params=res["k"],
                         logLik=round(res["loglik"], 2), AIC=round(res["aic"], 2), AIC_kind="full",
                         c_index=round(cidx, 4), IBS=round(ibs, 4)))
        models.append((lab, res, D))
    return rows, models


# ── Soft (smoothly-differentiable) change-point: logistic blend of the two hazards ─
def _h_soft(a, p, g1, T, s):
    a = np.asarray(a, float)
    hw = p * g1 * np.clip(a, 1e-8, None) ** (p - 1.0)        # Weibull hazard
    g2 = p * g1 * T ** (p - 1.0)                              # plateau level (continuity value)
    w = 1.0 / (1.0 + np.exp(np.clip((a - T) / s, -50, 50)))   # ~1 below T, ~0 above T
    return g2 + (hw - g2) * w


def fit_soft(dur, event, entry, t_bounds, n_starts=20):
    """Soft change-point on the age axis (left-truncated); H computed numerically."""
    dur = np.asarray(dur, float); event = np.asarray(event, float); entry = np.asarray(entry, float)
    a_grid = np.linspace(50.0, float(dur.max()) + 1.0, 700)
    lt_lo, lt_hi = np.log(t_bounds[0]), np.log(t_bounds[1])
    ls_lo, ls_hi = np.log(0.25), np.log(10.0)
    lp_lo, lp_hi = np.log(0.3), np.log(6.0)

    def Hcum(p, g1, T, s):
        h = _h_soft(a_grid, p, g1, T, s)
        return np.concatenate([[0.0], np.cumsum(np.diff(a_grid) * (h[:-1] + h[1:]) / 2.0)])

    def nll(theta):
        if not np.all(np.isfinite(theta)):
            return 1e12
        if not (lp_lo <= theta[0] <= lp_hi) or not (lt_lo <= theta[2] <= lt_hi) or not (ls_lo <= theta[3] <= ls_hi):
            return 1e12
        p, g1, T, s = np.exp(theta[0]), np.exp(theta[1]), np.exp(theta[2]), np.exp(theta[3])
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            Hc = Hcum(p, g1, T, s)
            hev = _h_soft(dur, p, g1, T, s)
            ll = np.sum(event * np.log(np.clip(hev, 1e-300, None))
                        - (np.interp(dur, a_grid, Hc) - np.interp(entry, a_grid, Hc)))
        return -ll if np.isfinite(ll) else 1e12

    nev = event.sum(); med = np.median(dur[event == 1])
    base = np.array([0.0, np.log(max(nev / dur.sum(), 1e-4)), np.clip(np.log(med), lt_lo, lt_hi), np.log(2.0)])
    rng = np.random.default_rng(20240601); best = None
    for k in range(n_starts):
        x0 = base.copy()
        if k:
            x0 += rng.normal(0, 0.4, 4); x0[0] = rng.uniform(lp_lo, lp_hi)
            x0[2] = np.clip(x0[2], lt_lo, lt_hi); x0[3] = np.clip(x0[3], ls_lo, ls_hi)
        try:
            res = minimize(nll, x0, method="Nelder-Mead", options={"maxiter": 12000, "xatol": 1e-7, "fatol": 1e-7})
        except Exception:
            continue
        if np.isfinite(res.fun) and (best is None or res.fun < best.fun):
            best = res
    th = best.x; p, g1, T, s = (float(np.exp(th[0])), float(np.exp(th[1])), float(np.exp(th[2])), float(np.exp(th[3])))
    Hc = Hcum(p, g1, T, s); ll = -best.fun
    return {"p": p, "g1": g1, "T": T, "s": s, "loglik": ll, "aic": 2 * 4 - 2 * ll, "k": 4,
            "S": lambda a, ref: np.exp(-(np.interp(a, a_grid, Hc) - np.interp(ref, a_grid, Hc))),
            "h": lambda a: _h_soft(a, p, g1, T, s)}


def kink_experiment(c, outdir):
    """Compare 3 smoothness settings of the change-point (age axis, no covariates)."""
    exitg = c["exit_age"].values; entry = c["entry_age"].values; ev = c["event"].values
    tb = (float(np.quantile(exitg, 0.05)), float(np.quantile(exitg, 0.95)))
    n = len(c); onescol = np.ones((n, 1)); X1 = np.ones((1, 1))
    cont = fit_changepoint(exitg, ev, onescol, entry=entry, t_bounds=tb)
    hard = fit_cp_general(exitg, ev, {"g1": onescol, "g2": onescol, "T": onescol}, entry=entry, t_bounds=tb)
    soft = fit_soft(exitg, ev, entry, tb)
    Tc, pc, g1c = float(np.exp(cont["beta"][0])), cont["p"], cont["g1"]; g2c = pc * g1c * Tc ** (pc - 1)
    lp, bg1, bg2, bT = _split_theta(hard["theta"], hard["shapes"])
    ph, g1h, g2h, Th = float(np.exp(lp)), float(np.exp(bg1[0])), float(np.exp(bg2[0])), float(np.exp(bT[0]))

    ref = 65.0; a_hi = float(np.quantile(exitg, 0.95)); grid = np.linspace(ref, a_hi, 300)
    Xg = np.ones((len(grid), 1)); Dg = {"g1": Xg, "g2": Xg, "T": Xg}; D1 = {"g1": X1, "g2": X1, "T": X1}
    h_cont = _h_changepoint(cont["theta"], grid, Xg)
    h_hard = _h_disc(hard["theta"], grid, Dg, hard["shapes"])
    h_soft = soft["h"](grid)
    S_cont = np.exp(-(_H_changepoint(cont["theta"], grid, Xg) - _H_changepoint(cont["theta"], np.array([ref]), X1)[0]))
    S_hard = np.exp(-(_H_disc(hard["theta"], grid, Dg, hard["shapes"])
                      - _H_disc(hard["theta"], np.array([ref]), D1, hard["shapes"])[0]))
    S_soft = soft["S"](grid, ref)
    xs, ss = _km_lt_cond(entry, exitg, ev, ref)

    HARD_C, CONT_C, SOFT_C = "#d7191c", "#2c7bb6", "#1a9641"
    fig, (axh, axs) = plt.subplots(1, 2, figsize=(14, 6))
    axh.plot(grid, h_hard, color=HARD_C, lw=2, label=f"non-differentiable (free γ₂), T={Th:.1f}")
    axh.plot(grid, h_cont, color=CONT_C, lw=2, label=f"differentiable / current, T={Tc:.1f}")
    axh.plot(grid, h_soft, color=SOFT_C, lw=2, label=f"soft-differentiable, T={soft['T']:.1f}, s={soft['s']:.1f}")
    for T, col in ((Th, HARD_C), (Tc, CONT_C), (soft["T"], SOFT_C)):
        if ref <= T <= a_hi:
            axh.axvline(T, color=col, ls=":", lw=1, alpha=0.6)
    axh.set_title("Hazard h(a) at the change-point\n(jump vs kink vs smooth)", fontsize=11, fontweight="bold")
    axh.set_xlabel("Age (years)", fontweight="bold"); axh.set_ylabel("hazard", fontweight="bold")
    axs.step(xs, ss, where="post", color="black", lw=2.3, label="Kaplan–Meier (observed)")
    axs.plot(grid, S_hard, color=HARD_C, lw=2, ls="--", label=f"non-differentiable (AIC={hard['aic']:.0f})")
    axs.plot(grid, S_cont, color=CONT_C, lw=2, ls="--", label=f"differentiable / current (AIC={cont['aic']:.0f})")
    axs.plot(grid, S_soft, color=SOFT_C, lw=2, ls="--", label=f"soft-differentiable (AIC={soft['aic']:.0f})")
    axs.set_title("Survival S(a|65) vs observed KM", fontsize=11, fontweight="bold")
    axs.set_xlabel("Age (years)", fontweight="bold"); axs.set_ylabel("Survival (AD-free)", fontweight="bold")
    axs.set_ylim(0, 1.02)
    for ax in (axh, axs):
        ax.grid(alpha=0.3); ax.legend(fontsize=8.5)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    cap = (
        "Change-point transition — three smoothness settings (age axis, no covariates, left-truncated; lower AIC = better fit):\n"
        f"  non-differentiable (free γ₂):  h(a)=p·γ₁·a^(p−1) [a<T],  γ₂ [a≥T]  (hazard JUMPS).  "
        f"p={ph:.2f}, γ₁={g1h:.2g}, γ₂={g2h:.3g}, T={Th:.1f}.   logLik={hard['loglik']:.1f}, AIC={hard['aic']:.1f}\n"
        f"  differentiable / current (continuity):  γ₂=p·γ₁·T^(p−1)  (hazard CONTINUOUS, kinked).  "
        f"p={pc:.2f}, γ₁={g1c:.2g}, T={Tc:.1f}.   logLik={cont['loglik']:.1f}, AIC={cont['aic']:.1f}\n"
        f"  soft-differentiable:  h(a)=γ₂+(p·γ₁·a^(p−1)−γ₂)·σ((T−a)/s),  σ=logistic, width s  (hazard SMOOTH).  "
        f"p={soft['p']:.2f}, γ₁={soft['g1']:.2g}, T={soft['T']:.1f}, s={soft['s']:.2f}.   logLik={soft['loglik']:.1f}, AIC={soft['aic']:.1f}"
    )
    fig.subplots_adjust(bottom=0.30)
    fig.text(0.012, 0.012, cap, ha="left", va="bottom", fontsize=8.0, family="DejaVu Sans", linespacing=1.5)
    for ext in ("png", "pdf"):
        fig.savefig(outdir / f"kink_smoothness_experiment_age.{ext}", bbox_inches="tight", dpi=300)
    plt.close()
    print(f"\n[kink experiment]  AIC: non-diff={hard['aic']:.1f}  current={cont['aic']:.1f}  soft={soft['aic']:.1f}")


# ── Soft change-point WITH covariates on T (per-subject T), fixed s ─────────────
def _interp_rows(Hcum, q, a_grid):
    q = np.asarray(q, float)
    idx = np.clip(np.searchsorted(a_grid, q) - 1, 0, len(a_grid) - 2)
    rows = np.arange(Hcum.shape[0])
    a0 = a_grid[idx]; a1 = a_grid[idx + 1]
    H0 = Hcum[rows, idx]; H1 = Hcum[rows, idx + 1]
    return H0 + (q - a0) / (a1 - a0) * (H1 - H0)


def _interp_scalar(Hcum, a, a_grid):
    idx = int(np.clip(np.searchsorted(a_grid, a) - 1, 0, len(a_grid) - 2))
    a0, a1 = a_grid[idx], a_grid[idx + 1]
    return Hcum[:, idx] + (a - a0) / (a1 - a0) * (Hcum[:, idx + 1] - Hcum[:, idx])


def fit_soft_cov(dur, event, entry, X, s, t_bounds, n_starts=10):
    """Soft change-point with covariates on log T. s=None -> fit s (free); else s fixed."""
    dur = np.asarray(dur, float); event = np.asarray(event, float)
    entry = np.asarray(entry, float); X = np.asarray(X, float)
    n, kT = len(dur), X.shape[1]
    free_s = s is None
    a_grid = np.linspace(50.0, float(dur.max()) + 1.0, 900)   # fine enough to resolve sharp (small-s) transitions
    lt_lo, lt_hi = np.log(t_bounds[0]), np.log(t_bounds[1]); lp_lo, lp_hi = np.log(0.3), np.log(6.0)
    ls_lo, ls_hi = np.log(0.25), np.log(10.0)

    def Hcum_of(p, g1, bT, sv):
        T = np.exp(X @ bT); g2 = p * g1 * T ** (p - 1.0)
        hw = p * g1 * np.clip(a_grid, 1e-8, None) ** (p - 1.0)
        W = 1.0 / (1.0 + np.exp(np.clip((a_grid[None, :] - T[:, None]) / sv, -50, 50)))
        hmat = g2[:, None] + (hw[None, :] - g2[:, None]) * W
        dH = np.diff(a_grid)[None, :] * (hmat[:, :-1] + hmat[:, 1:]) / 2.0
        return T, g2, np.concatenate([np.zeros((n, 1)), np.cumsum(dH, axis=1)], axis=1)

    def nll(theta):
        if not np.all(np.isfinite(theta)) or not (lp_lo <= theta[0] <= lp_hi):
            return 1e12
        sv = np.exp(theta[-1]) if free_s else s
        if free_s and not (ls_lo <= theta[-1] <= ls_hi):
            return 1e12
        bT = theta[2:2 + kT]; lpT = X @ bT
        if lpT.min() < lt_lo or lpT.max() > lt_hi:
            return 1e12
        p, g1 = np.exp(theta[0]), np.exp(theta[1])
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            T, g2, Hcum = Hcum_of(p, g1, bT, sv)
            H_exit = _interp_rows(Hcum, dur, a_grid); H_entry = _interp_rows(Hcum, entry, a_grid)
            hw_e = p * g1 * np.clip(dur, 1e-8, None) ** (p - 1.0)
            W_e = 1.0 / (1.0 + np.exp(np.clip((dur - T) / sv, -50, 50)))
            hev = g2 + (hw_e - g2) * W_e
            ll = np.sum(event * np.log(np.clip(hev, 1e-300, None)) - (H_exit - H_entry))
        return -ll if np.isfinite(ll) else 1e12

    med = np.median(dur[event == 1])
    base = np.zeros(2 + kT + (1 if free_s else 0))
    base[1] = np.log(max(event.sum() / dur.sum(), 1e-4))
    base[2 + kT - 1] = np.clip(np.log(med), lt_lo, lt_hi)
    if free_s:
        base[-1] = np.log(2.0)
    rng = np.random.default_rng(20240601); best = None
    for k in range(n_starts):
        x0 = base.copy()
        if k:
            x0 += rng.normal(0, 0.4, base.shape); x0[0] = rng.uniform(lp_lo, lp_hi)
            x0[2 + kT - 1] = np.clip(x0[2 + kT - 1], lt_lo, lt_hi)
            if free_s:
                x0[-1] = np.clip(x0[-1], ls_lo, ls_hi)
        try:
            res = minimize(nll, x0, method="Nelder-Mead", options={"maxiter": 12000, "xatol": 1e-6, "fatol": 1e-6})
        except Exception:
            continue
        if np.isfinite(res.fun) and (best is None or res.fun < best.fun):
            best = res
    th = best.x; p, g1, bT = float(np.exp(th[0])), float(np.exp(th[1])), th[2:2 + kT]
    sv = float(np.exp(th[-1])) if free_s else s; ll = -best.fun; kk = len(th)
    _, _, Hcum = Hcum_of(p, g1, bT, sv)

    def predS(grid, ref):
        Href = _interp_scalar(Hcum, ref, a_grid)
        return np.column_stack([np.exp(-(_interp_scalar(Hcum, a, a_grid) - Href)) for a in grid])
    return {"p": p, "g1": g1, "bT": bT, "s": sv, "loglik": ll, "aic": 2 * kk - 2 * ll, "predict_S": predS}


# ── IPCW Brier (left-truncation-aware) + c-index ────────────────────────────────
def km_censoring(times, events):
    d = pd.DataFrame({"t": times, "c": 1 - np.asarray(events)}).sort_values("t")
    G, ts, gs = 1.0, [0.0], [1.0]
    for t in sorted(d["t"].unique()):
        n = int((d["t"] >= t).sum()); m = int(((d["t"] == t) & (d["c"] == 1)).sum())
        if n > 0 and m > 0:
            G *= (1 - m / n)
        ts.append(t); gs.append(G)
    ts = np.array(ts); gs = np.array(gs)
    return lambda q: gs[np.clip(np.searchsorted(ts, q, side="right") - 1, 0, len(gs) - 1)]


def ipcw_ibs(S_grid, dur, event, grid, entry=None):
    """Integrated Brier score (Graf IPCW). At each age only subjects AT RISK (entry ≤ age)
    are scored; survival is clipped to [0,1]."""
    G = km_censoring(dur, event)
    dur = np.asarray(dur, float); event = np.asarray(event, float)
    ent = None if entry is None else np.asarray(entry, float)
    bs = []
    for j, t in enumerate(grid):
        S = np.clip(S_grid[:, j], 0.0, 1.0)
        atrisk = np.ones(len(dur), bool) if ent is None else (ent <= t)
        had_event = (dur <= t) & (event == 1) & atrisk
        survived = (dur > t) & atrisk
        gt = max(G(t), 1e-6)
        w_ev = had_event / np.clip(G(dur), 1e-6, None)
        term = (S ** 2) * w_ev + ((1 - S) ** 2) * survived / gt
        bs.append(term.sum() / max(int(atrisk.sum()), 1))
    bs = np.array(bs)
    return float(np.sum(np.diff(grid) * (bs[:-1] + bs[1:]) / 2) / (grid[-1] - grid[0]))


def cindex_at(S_h, dur, event):
    return float(concordance_index(dur, np.clip(S_h, 0.0, 1.0), event))


# ── Cox: time axis (time-varying age coef), Breslow S(t|x) for any covariates ────
def cox_time_fit_predict(cohort, covcols, grid):
    df = cohort[["t", "event", "age_z"] + covcols].rename(columns={"event": "e"}).copy()
    df["id"] = np.arange(len(df))
    long = to_episodic_format(df, duration_col="t", event_col="e", id_col="id", time_gaps=0.5)
    long["age_x_logt"] = long["age_z"] * np.log(long["stop"].clip(lower=1e-3))
    feat = ["age_z", "age_x_logt"] + covcols
    ctv = CoxTimeVaryingFitter(penalizer=1e-3)
    ctv.fit(long[["id", "start", "stop", "e"] + feat], id_col="id", event_col="e",
            start_col="start", stop_col="stop", show_progress=False)
    ca = float(ctv.params_["age_z"]); ci = float(ctv.params_["age_x_logt"])
    ccov = ctv.params_[covcols].values.astype(float)

    start = long["start"].values; stop = long["stop"].values; e = long["e"].values
    az = long["age_z"].values; covconst_long = long[covcols].values.astype(float) @ ccov
    ev_times = np.array(sorted(stop[e == 1]))
    dH0 = []
    for s in ev_times:
        inrisk = (start < s) & (stop >= s)
        denom = np.exp(ca * az[inrisk] + ci * az[inrisk] * np.log(s) + covconst_long[inrisk]).sum()
        dH0.append(np.sum((stop == s) & (e == 1)) / denom if denom > 0 else 0.0)
    dH0 = np.array(dH0); logev = np.log(ev_times)
    saz = cohort["age_z"].values; covconst_s = cohort[covcols].values.astype(float) @ ccov
    S = np.ones((len(cohort), len(grid)))
    for j, t0 in enumerate(grid):
        m = ev_times <= t0
        if m.sum() == 0:
            continue
        LP = np.outer(saz, ca + ci * logev[m]) + covconst_s[:, None]
        S[:, j] = np.exp(-(np.exp(LP) * dH0[m][None, :]).sum(axis=1))
    return ctv, S


# ── Cox: age axis (left-truncated), per-subject survival curves ─────────────────
def cox_age_fit_predict(cohort, covcols, grid):
    df = cohort[["entry_age", "exit_age", "event"] + covcols].rename(columns={"event": "e"}).copy()
    cph = CoxPHFitter(penalizer=1e-3)
    cph.fit(df, duration_col="exit_age", event_col="e", entry_col="entry_age")
    fine = np.unique(np.concatenate([
        np.linspace(float(cohort["entry_age"].min()), float(cohort["exit_age"].max()), 400),
        np.asarray(grid, float)]))
    M = cph.predict_survival_function(cohort[covcols], times=list(fine)).values   # (len(fine), n)
    ent = cohort["entry_age"].values
    S = np.ones((len(cohort), len(grid)))
    for i in range(len(cohort)):
        col = M[:, i]
        s_ent = max(np.interp(ent[i], fine, col), 1e-8)
        S[i, :] = np.clip(np.interp(grid, fine, col) / s_ent, 0, 1)
    return cph, S, fine, M


# ── Model-equation strings (covariate-aware) ────────────────────────────────────
def eq_weibull(wt, axis, covlabels):
    b = wt["beta"]; i = 0; terms = []
    if axis == "time":
        terms.append(f"{b[i]:+.2f}·age_z"); i += 1
    for lab in covlabels:
        terms.append(f"{b[i]:+.2f}·{lab}"); i += 1
    sym = "t" if axis == "time" else "a"
    logT = f"log T = {b[i]:+.2f} " + " ".join(terms)
    return (f"Weibull→exp change-point:  h({sym}|x)=p·γ₁·{sym}^(p−1) [{sym}<T];  "
            f"γ₂=p·γ₁·T^(p−1) [{sym}≥T];  {logT}\n"
            f"   fitted:  p={wt['p']:.2f},  γ₁={wt['g1']:.3g}")


def eq_cox_time(ctv, covcols, covlabels):
    p = ctv.params_
    cterms = " ".join(f"{float(p[cc]):+.2f}·{lab}" for cc, lab in zip(covcols, covlabels))
    return (f"Cox PH (time-varying age):  h(t|x)=h₀(t)·exp({float(p['age_z']):+.2f}·age_z "
            f"{float(p['age_x_logt']):+.2f}·age_z·log t {cterms})")


def eq_cox_age(cph, covcols, covlabels):
    p = cph.params_
    cterms = " ".join(f"{float(p[cc]):+.2f}·{lab}" for cc, lab in zip(covcols, covlabels))
    return f"Cox PH:  h(a|x)=h₀(a)·exp({cterms})"


def _profile_T(axis, covcols, apoe_val):
    """T-design row with APOE4 = apoe_val, other covariates at 0, mean age (age_z=0)."""
    vals = [(apoe_val if cc == covcols[0] else 0.0) for cc in covcols]
    return np.array([([0.0] if axis == "time" else []) + vals + [1.0]], float)


# ── Number-at-risk row (returns lowest y) + equation block ──────────────────────
def _n_at_risk_row(ax, groups, kind):
    xlim = ax.get_xlim()
    ticks = [x for x in ax.get_xticks() if xlim[0] <= x <= xlim[1]]
    tr = ax.get_xaxis_transform()
    ax.text(0.0, -0.17, "Number at risk", transform=ax.transAxes,
            fontsize=7.5, fontweight="bold", va="top")
    y = -0.24
    for j, (lab, col, t2, ev, ent) in enumerate(groups):
        y = -0.24 - 0.058 * j
        for x in ticks:
            if kind == "time":
                nn = int(np.sum(np.asarray(t2, float) >= x))
            else:
                nn = int(np.sum((np.asarray(ent, float) <= x) & (np.asarray(t2, float) >= x)))
            ax.text(x, y, str(nn), transform=tr, ha="center", va="top", fontsize=6.8, color=col)
        ax.text(xlim[0], y, lab + "  ", transform=tr, ha="right", va="top",
                fontsize=6.8, color=col, fontweight="bold")
    return y


def _add_eqs(ax, y_low, eq_text):
    ax.text(0.0, y_low - 0.07, eq_text, transform=ax.transAxes, fontsize=6.6,
            family="DejaVu Sans", va="top", linespacing=1.5)


def _finish_calibration(ax, fig, groups, kind, eq_text, fname, dirs):
    fig.tight_layout()
    y_low = _n_at_risk_row(ax, groups, kind)
    _add_eqs(ax, y_low, eq_text)
    for d in dirs:
        for ext in ("png", "pdf"):
            fig.savefig(d / f"{fname}.{ext}", bbox_inches="tight", dpi=300)
    plt.close()


# ── Hazard + survival figures (APOE4 carrier vs non-carrier; equations underneath) ─
def _plot_weibull_haz_surv(wt, eq_text, axis, covcols, c, wdir):
    is_t = axis == "time"
    if is_t:
        x_hi = float(np.quantile(c["t"], 0.92)); xx = np.linspace(EPS, x_hi, 300)
        ref = 0.0; xlabel = "Time (years)"; sym = "t"
        dur, ev, ent = c["t"].values, c["event"].values, None
    else:
        x_lo, x_hi = 56.0, float(np.quantile(c["exit_age"], 0.95)); xx = np.linspace(x_lo, x_hi, 300)
        ref = x_lo; xlabel = "Age (years)"; sym = "a"
        dur, ev, ent = c["exit_age"].values, c["event"].values, c["entry_age"].values
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.8))
    for ap, lab, col in ((0, "Non-carrier", BLUE), (1, "APOE4 carrier", RED)):
        row = _profile_T(axis, covcols, ap); Xr = np.repeat(row, len(xx), axis=0)
        h = _h_changepoint(wt["theta"], xx, Xr)
        S = wt["predict_S"](Xr, xx, ref)
        T = float(np.exp(row[0] @ wt["beta"]))
        ax1.plot(xx, h, color=col, lw=2, label=f"{lab} (T={T:.1f}y)")
        ax2.plot(xx, S, color=col, lw=2, label=f"{lab} (T={T:.1f}y)")
        if xx[0] <= T <= xx[-1]:
            ax1.axvline(T, color=col, ls=":", lw=1.2, alpha=0.7)
            ax2.axvline(T, color=col, ls=":", lw=1.2, alpha=0.7)
    ax1.set_title(f"Hazard h({sym}|x) — Weibull→exponential change-point",
                  fontsize=11, fontweight="bold")
    ax1.set_ylabel("hazard", fontweight="bold")
    surv_ttl = "Survival S(t|x) (AD-free)" if is_t else f"Survival S(a|x), conditioned at age {ref:.0f}"
    ax2.set_title(surv_ttl, fontsize=11, fontweight="bold"); ax2.set_ylabel(f"S({sym})", fontweight="bold")
    for ax in (ax1, ax2):
        ax.set_xlabel(xlabel, fontweight="bold"); ax.legend(fontsize=9); ax.grid(alpha=0.3)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    _n_at_risk_row(ax2, [("Cohort", "black", dur, ev, ent)], axis)
    fig.subplots_adjust(bottom=0.40)
    fig.text(0.012, 0.012, eq_text, ha="left", va="bottom", fontsize=7.4,
             family="DejaVu Sans", linespacing=1.5)
    for ext in ("png", "pdf"):
        fig.savefig(wdir / f"weibull_{axis}_hazard_survival.{ext}", bbox_inches="tight", dpi=300)
    plt.close()


# ── Kaplan–Meier helpers ────────────────────────────────────────────────────────
def _km_marginal(times, events):
    times = np.asarray(times, float); events = np.asarray(events, int)
    S, xs, ss = 1.0, [], []
    for t in sorted(set(times[events == 1])):
        n = int((times >= t).sum()); d = int(((times == t) & (events == 1)).sum())
        if n > 0:
            S *= (1 - d / n)
        xs.append(t); ss.append(S)
    return np.array(xs), np.array(ss)


def _km_lt_cond(entry, exit_, event, ref):
    entry = np.asarray(entry, float); exit_ = np.asarray(exit_, float); event = np.asarray(event, int)
    S, xs, ss = 1.0, [ref], [1.0]
    for a in sorted(a for a in set(exit_[event == 1]) if a >= ref):
        n = int(np.sum((entry <= a) & (exit_ >= a))); d = int(np.sum((exit_ == a) & (event == 1)))
        if n > 0 and d > 0:
            S *= (1 - d / n)
        xs.append(a); ss.append(S)
    return np.array(xs), np.array(ss)


def _style_legend(ax):
    for sty, lab in (("-", "Kaplan–Meier (observed)"), ("--", "Weibull change-point"), (":", "Cox PH")):
        ax.plot([], [], color="grey", lw=2, ls=sty, label=lab)


# ── Calibration: time axis ──────────────────────────────────────────────────────
def _cox_age_cond(M, fine, idx, grid, ref):
    """Mean over subjects `idx` of Cox survival conditioned at age `ref`."""
    out = []
    for i in idx:
        col = M[:, i]
        out.append(np.interp(grid, fine, col) / max(np.interp(ref, fine, col), 1e-8))
    return np.mean(out, axis=0)


def _plot_cal_time(grid, S_w, S_c, dur, event, eq_text, wdir, cdir, strat=None):
    fig, ax = plt.subplots(figsize=(8.6, 6))
    if strat is not None:
        arr, lab0, lab1, short = strat
        for g, lab, col in ((0.0, lab0, BLUE), (1.0, lab1, RED)):
            m = arr == g
            if m.sum() == 0:
                continue
            xs, ss = _km_marginal(dur[m], event[m])
            ax.step(xs, ss, where="post", color=col, lw=2.3)
            ax.plot(grid, S_w[m].mean(axis=0), color=col, lw=1.7, ls="--")
            ax.plot(grid, S_c[m].mean(axis=0), color=col, lw=1.7, ls=":")
            ax.plot([], [], color=col, lw=3, label=f"{lab} (n={int(m.sum())}, ev={int(event[m].sum())})")
        _style_legend(ax)
        groups = [(lab0, BLUE, dur[arr == 0.0], event[arr == 0.0], None),
                  (lab1, RED, dur[arr == 1.0], event[arr == 1.0], None)]
        ttl = f"Calibration by {short.upper()} (time axis)"; fname = f"calibration_time_by_{short}"
    else:
        xs, ss = _km_marginal(dur, event)
        ax.step(xs, ss, where="post", color="black", lw=2, label="Kaplan–Meier (observed)")
        ax.plot(grid, S_w.mean(axis=0), color=RED, lw=2, ls="--", label="Weibull change-point (mean S)")
        ax.plot(grid, S_c.mean(axis=0), color=BLUE, lw=2, ls="-.", label="Cox PH (mean S)")
        groups = [("All cohort", "black", dur, event, None)]
        ttl = "Calibration (time axis)"; fname = "calibration_time"
    ax.set_title(f"{ttl} — predicted vs observed survival\nCN+MCI → AD", fontsize=11, fontweight="bold")
    ax.set_xlabel("Time (years)", fontweight="bold"); ax.set_ylabel("Survival (AD-free)", fontweight="bold")
    ax.set_ylim(0, 1.02); ax.grid(alpha=0.3); ax.legend(fontsize=8.3)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    _finish_calibration(ax, fig, groups, "time", eq_text, fname, (wdir, cdir))


def _subset_design(design, m):
    """Subset a Weibull design by a boolean mask (array design or per-param dict design)."""
    if isinstance(design, dict):
        return {k: v[m] for k, v in design.items()}
    return design[m]


def _plot_cal_age(c, wmodel, wdesign, fine, M, eq_text, wdir, cdir, ref, strat=None):
    entry = c["entry_age"].values; exitg = c["exit_age"].values; ev = c["event"].values
    a_hi = float(np.quantile(exitg, 0.92)); grid = np.linspace(ref, a_hi, 30)
    n_ref = int(np.sum((entry <= ref) & (exitg >= ref)))
    fig, ax = plt.subplots(figsize=(8.6, 6))
    if strat is not None:
        arr, lab0, lab1, short = strat
        for g, lab, col in ((0.0, lab0, BLUE), (1.0, lab1, RED)):
            m = arr == g
            if m.sum() == 0:
                continue
            xs, ss = _km_lt_cond(entry[m], exitg[m], ev[m], ref)
            ax.step(xs, ss, where="post", color=col, lw=2.3)
            Sw = np.column_stack([wmodel["predict_S"](_subset_design(wdesign, m), a, ref) for a in grid]).mean(axis=0)
            Sc = _cox_age_cond(M, fine, np.where(m)[0], grid, ref)
            ax.plot(grid, Sw, color=col, lw=1.7, ls="--"); ax.plot(grid, Sc, color=col, lw=1.7, ls=":")
            ax.plot([], [], color=col, lw=3, label=f"{lab} (n={int(m.sum())}, ev={int(ev[m].sum())})")
        _style_legend(ax)
        groups = [(lab0, BLUE, exitg[arr == 0.0], ev[arr == 0.0], entry[arr == 0.0]),
                  (lab1, RED, exitg[arr == 1.0], ev[arr == 1.0], entry[arr == 1.0])]
        ttl = f"Calibration by {short.upper()} (age axis)"; fname = f"calibration_age_by_{short}_cond{int(ref)}"
    else:
        xs, ss = _km_lt_cond(entry, exitg, ev, ref)
        ax.step(xs, ss, where="post", color="black", lw=2, label=f"Kaplan–Meier (cond. {ref:.0f}y)")
        Sw = np.column_stack([wmodel["predict_S"](wdesign, a, ref) for a in grid]).mean(axis=0)
        Sc = _cox_age_cond(M, fine, np.arange(len(c)), grid, ref)
        ax.plot(grid, Sw, color=RED, lw=2, ls="--", label="Weibull change-point (mean S)")
        ax.plot(grid, Sc, color=BLUE, lw=2, ls="-.", label="Cox PH (mean S)")
        groups = [("All cohort", "black", exitg, ev, entry)]
        ttl = "Calibration (age axis)"; fname = f"calibration_age_cond{int(ref)}"
    ax.set_title(f"{ttl} — CN+MCI → AD, conditioned at age {ref:.0f} (n at risk = {n_ref})",
                 fontsize=11, fontweight="bold")
    ax.set_xlabel("Age (years)", fontweight="bold"); ax.set_ylabel("Survival (AD-free)", fontweight="bold")
    ax.set_ylim(0, 1.02); ax.grid(alpha=0.3); ax.legend(fontsize=8.3)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    _finish_calibration(ax, fig, groups, "age", eq_text, fname, (wdir, cdir))


# ── Observed KM stratified by a covariate (model-free) ──────────────────────────
def _plot_km_strat(c, axis, strat, wdir, cdir, ref=None):
    arr, lab0, lab1, short = strat
    ev = c["event"].values
    fig, ax = plt.subplots(figsize=(8.4, 6))
    if axis == "time":
        dur = c["t"].values
        for g, lab, col in ((0.0, lab0, BLUE), (1.0, lab1, RED)):
            m = arr == g
            if m.sum() == 0:
                continue
            xs, ss = _km_marginal(dur[m], ev[m])
            ax.step(xs, ss, where="post", color=col, lw=2.3,
                    label=f"{lab} (n={int(m.sum())}, ev={int(ev[m].sum())})")
        groups = [(lab0, BLUE, dur[arr == 0.0], ev[arr == 0.0], None),
                  (lab1, RED, dur[arr == 1.0], ev[arr == 1.0], None)]
        xlabel = "Time (years)"; fname = f"km_time_by_{short}"
        ttl = f"Observed Kaplan–Meier by {short.upper()} (time axis)"
    else:
        entry = c["entry_age"].values; exitg = c["exit_age"].values
        for g, lab, col in ((0.0, lab0, BLUE), (1.0, lab1, RED)):
            m = arr == g
            if m.sum() == 0:
                continue
            xs, ss = _km_lt_cond(entry[m], exitg[m], ev[m], ref)
            ax.step(xs, ss, where="post", color=col, lw=2.3,
                    label=f"{lab} (n={int(m.sum())}, ev={int(ev[m].sum())})")
        groups = [(lab0, BLUE, exitg[arr == 0.0], ev[arr == 0.0], entry[arr == 0.0]),
                  (lab1, RED, exitg[arr == 1.0], ev[arr == 1.0], entry[arr == 1.0])]
        xlabel = "Age (years)"; fname = f"km_age_by_{short}_cond{int(ref)}"
        ttl = f"Observed Kaplan–Meier by {short.upper()} (age axis, conditioned at {ref:.0f})"
    ax.set_title(f"{ttl}\nCN+MCI → AD", fontsize=11, fontweight="bold")
    ax.set_xlabel(xlabel, fontweight="bold"); ax.set_ylabel("Survival (AD-free)", fontweight="bold")
    ax.set_ylim(0, 1.02); ax.grid(alpha=0.3); ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    _n_at_risk_row(ax, groups, axis)
    for d in (wdir, cdir):
        for ext in ("png", "pdf"):
            fig.savefig(d / f"{fname}.{ext}", bbox_inches="tight", dpi=300)
    plt.close()


# ── Calibration overlaying ALL Weibull variants (T / T+γ₁ / T+γ₂ / all / γ₂) ────
VAR_COLORS = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#666666"]


def _plot_variant_cal(axis, c, cont_res, cont_X, var_models, ref, wdir, cdir):
    event = c["event"].values
    if axis == "time":
        dur = c["t"].values; entry = None; ent_pred = 0.0
        x_hi = float(np.quantile(dur, 0.92)); grid = np.linspace(0.5, x_hi, 30)
        xs, ss = _km_marginal(dur, event); xlabel = "Time (years)"
        ttl = "Calibration — all Weibull variants (time axis)"; fname = "calibration_time_variants"
    else:
        dur = c["exit_age"].values; entry = c["entry_age"].values; ent_pred = ref
        a_hi = float(np.quantile(dur, 0.92)); grid = np.linspace(ref, a_hi, 30)
        xs, ss = _km_lt_cond(entry, dur, event, ref); xlabel = "Age (years)"
        ttl = f"Calibration — all Weibull variants (age axis, cond {ref:.0f})"
        fname = f"calibration_age_variants_cond{int(ref)}"
    fig, ax = plt.subplots(figsize=(9, 6.3))
    ax.step(xs, ss, where="post", color="black", lw=2.4, label="Kaplan–Meier (observed)")
    Sc = np.column_stack([cont_res["predict_S"](cont_X, g, ent_pred) for g in grid]).mean(axis=0)
    ax.plot(grid, Sc, color="#2c7bb6", lw=2.2, ls="--", label="T only (continuity, γ₂ tied)")
    for (lab, res, D), col in zip(var_models, VAR_COLORS):
        Sv = np.column_stack([res["predict_S"](D, g, ent_pred) for g in grid]).mean(axis=0)
        ax.plot(grid, Sv, color=col, lw=1.8, label=f"covariates on {lab}")
    ax.set_title(f"{ttl}\nCN+MCI → AD  (mean predicted S vs observed KM)", fontsize=11, fontweight="bold")
    ax.set_xlabel(xlabel, fontweight="bold"); ax.set_ylabel("Survival (AD-free)", fontweight="bold")
    ax.set_ylim(0, 1.02); ax.grid(alpha=0.3); ax.legend(fontsize=8.2)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    _n_at_risk_row(ax, [("All cohort", "black", dur, event, entry)], axis)
    for d in (wdir, cdir):
        for ext in ("png", "pdf"):
            fig.savefig(d / f"{fname}.{ext}", bbox_inches="tight", dpi=300)
    plt.close()


def _plot_variant_cal_strat(axis, c, cont_res, cont_X, var_models, strat, ref, wdir, cdir):
    """Variant-overlay calibration split into two panels by the stratifying covariate."""
    arr, lab0, lab1, short = strat
    event = c["event"].values
    if axis == "time":
        dur = c["t"].values; entry = None; ent_pred = 0.0
        x_hi = float(np.quantile(dur, 0.92)); grid = np.linspace(0.5, x_hi, 30); xlabel = "Time (years)"
        ttl = f"Calibration — all Weibull variants by {short.upper()} (time axis)"
        fname = f"calibration_time_variants_by_{short}"
        km = lambda m: _km_marginal(dur[m], event[m])
    else:
        dur = c["exit_age"].values; entry = c["entry_age"].values; ent_pred = ref
        a_hi = float(np.quantile(dur, 0.92)); grid = np.linspace(ref, a_hi, 30); xlabel = "Age (years)"
        ttl = f"Calibration — all Weibull variants by {short.upper()} (age axis, cond {ref:.0f})"
        fname = f"calibration_age_variants_by_{short}_cond{int(ref)}"
        km = lambda m: _km_lt_cond(entry[m], dur[m], event[m], ref)
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.4), sharey=True)
    for ax, (g, lab) in zip(axes, ((0.0, lab0), (1.0, lab1))):
        m = arr == g
        if m.sum() == 0:
            ax.set_visible(False); continue
        xs, ss = km(m)
        ax.step(xs, ss, where="post", color="black", lw=2.2, label="KM (observed)")
        Sc = np.column_stack([cont_res["predict_S"](cont_X[m], gg, ent_pred) for gg in grid]).mean(axis=0)
        ax.plot(grid, Sc, color="#2c7bb6", lw=2, ls="--", label="T only (continuity)")
        for (vlab, res, D), col in zip(var_models, VAR_COLORS):
            Dm = {k: v[m] for k, v in D.items()}
            Sv = np.column_stack([res["predict_S"](Dm, gg, ent_pred) for gg in grid]).mean(axis=0)
            ax.plot(grid, Sv, color=col, lw=1.7, label=f"on {vlab}")
        ax.set_title(f"{lab} (n={int(m.sum())}, ev={int(event[m].sum())})", fontsize=11, fontweight="bold")
        ax.set_xlabel(xlabel, fontweight="bold"); ax.set_ylim(0, 1.02); ax.grid(alpha=0.3)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        _n_at_risk_row(ax, [(lab, "black", dur[m], event[m], None if entry is None else entry[m])], axis)
    axes[0].set_ylabel("Survival (AD-free)", fontweight="bold")
    axes[1].legend(fontsize=8, loc="upper right")
    fig.suptitle(f"{ttl} — CN+MCI → AD", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    for d in (wdir, cdir):
        for ext in ("png", "pdf"):
            fig.savefig(d / f"{fname}.{ext}", bbox_inches="tight", dpi=300)
    plt.close()


def _plot_variant_grid_apoe_sex(c, wa, Xa, va_models, cox_eq, wdir, cdir):
    """2×2 (APOE4 × Sex): KM vs Weibull with the covariate block on T-only (continuity),
    T+γ₁, T+γ₂, and T+γ₁+γ₂ — to show which placement separates the strata."""
    entry = c["entry_age"].values; exitg = c["exit_age"].values; ev = c["event"].values
    apoe = c["apoe"].values; sexm = c["sex_M"].values
    ref = 65.0; a_hi = float(np.quantile(exitg, 0.92)); grid = np.linspace(ref, a_hi, 30)
    show = [("continuity (T only)", wa, Xa, "grey", "--"),
            ("T + γ₁", va_models[1][1], va_models[1][2], "#1b9e77", "-"),
            ("T + γ₂", va_models[2][1], va_models[2][2], "#d95f02", "-"),
            ("T + γ₁ + γ₂", va_models[3][1], va_models[3][2], "#7570b3", "-")]
    Smats = {lab: np.column_stack([model["predict_S"](design, a, ref) for a in grid])
             for lab, model, design, _c, _ls in show}
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), sharex=True, sharey=True)
    for r, (ga, alab) in enumerate([(0.0, "APOE4 non-carrier"), (1.0, "APOE4 carrier")]):
        for ci, (gs, slab) in enumerate([(0.0, "Female"), (1.0, "Male")]):
            ax = axes[r][ci]; mm = (apoe == ga) & (sexm == gs)
            if mm.sum() == 0:
                ax.set_visible(False); continue
            xs, ss = _km_lt_cond(entry[mm], exitg[mm], ev[mm], ref)
            ax.step(xs, ss, where="post", color="black", lw=2.2, label="KM (observed)")
            for lab, _m, _d, ccol, ls in show:
                ax.plot(grid, Smats[lab][mm].mean(0), color=ccol, lw=1.9, ls=ls, label=lab)
            ax.set_title(f"{alab} · {slab} (n={int(mm.sum())}, ev={int(ev[mm].sum())})",
                         fontsize=10, fontweight="bold")
            ax.set_ylim(0, 1.02); ax.grid(alpha=0.3)
            ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
            if r == 1:
                ax.set_xlabel("Age (years)", fontweight="bold")
            if ci == 0:
                ax.set_ylabel("Survival (AD-free)", fontweight="bold")
    axes[0][1].legend(fontsize=8.5, loc="upper right")
    fig.suptitle("Weibull change-point — covariate placement (T-only / γ₁ / γ₂ / all) on the hazard\n"
                 "age axis (cond 65); stratified by APOE4 × Sex", fontsize=12, fontweight="bold")
    cap = ("Variants place the covariate block on different parameters of  h(a|x)=p·γ₁·a^(p−1) [a<T], γ₂ [a≥T]:\n"
           "  'continuity (T only)' shifts only the change-point T → barely separates the strata.\n"
           "  'T + γ₁' adds covariates to the early Weibull level γ₁; 'T + γ₂' to the post-change level γ₂; "
           "'T + γ₁ + γ₂' to both → these scale the hazard like Cox and separate the strata.\n"
           "  " + cox_eq)
    fig.subplots_adjust(bottom=0.15, top=0.91)
    fig.text(0.012, 0.01, cap, ha="left", va="bottom", fontsize=8, family="DejaVu Sans", linespacing=1.5)
    for d in (wdir, cdir):
        for ext in ("png", "pdf"):
            fig.savefig(d / "variant_placement_grid_apoe4_sex.{}".format(ext), bbox_inches="tight", dpi=300)
    plt.close()


# ── Master comparison table across covariate combinations ───────────────────────
def render_master_table(comp, outdirs):
    sections = [("time", "Time-from-baseline axis:   h(t | x)"),
                ("age", "Chronological-age axis (left-truncated):   h(a | x)")]
    headers = ["Covariates", "Model", "k", "logLik", "AIC", "AIC type", "C-index", "IBS"]
    body = []
    for axis, title in sections:
        body.append(([title] + [""] * 7, "section"))
        for _, r in comp[comp.axis == axis].iterrows():
            body.append(([r.covariates, r.model, str(int(r.n_params)), f"{r.logLik:.1f}",
                          f"{r.AIC:.1f}", r.AIC_kind, f"{r.c_index:.3f}", f"{r.IBS:.3f}"], "data"))
    fig, ax = plt.subplots(figsize=(13.5, 0.40 * (len(body) + 1) + 1.5)); ax.axis("off")
    ax.set_title("Model comparison — Weibull change-point vs Cox PH across covariate sets"
                 "   (CN+MCI → AD, N=535, 146 events)\n"
                 "C-index (higher better) & left-truncation-aware IBS (lower better);"
                 "  AIC 'full' (parametric) vs 'partial' (Cox) NOT on the same scale",
                 fontsize=10, fontweight="bold", loc="left", pad=12)
    tbl = ax.table(cellText=[c for c, _ in body], colLabels=headers, loc="center", cellLoc="center",
                   colWidths=[0.20, 0.22, 0.045, 0.09, 0.09, 0.095, 0.085, 0.075])
    tbl.auto_set_font_size(False); tbl.set_fontsize(8.5); tbl.scale(1, 1.35)
    for j in range(len(headers)):
        cc = tbl[(0, j)]; cc.set_facecolor("#222"); cc.set_text_props(color="white", weight="bold")
    ridx = 1
    for axis, _t in sections:
        for j in range(len(headers)):
            cc = tbl[(ridx, j)]; cc.set_facecolor("#d9e3f0"); cc.set_text_props(weight="bold", style="italic")
        ridx += 1
        sub = comp[comp.axis == axis]
        best_c, best_ibs = sub.c_index.max(), sub.IBS.min()
        for _, r in sub.iterrows():
            if abs(r.c_index - best_c) < 1e-9:
                tbl[(ridx, 6)].set_text_props(weight="bold")
            if abs(r.IBS - best_ibs) < 1e-9:
                tbl[(ridx, 7)].set_text_props(weight="bold")
            ridx += 1
    for d in outdirs:
        for ext in ("png", "pdf"):
            fig.savefig(d / f"model_comparison_master.{ext}", bbox_inches="tight", dpi=200)
    plt.close()


# label0, label1, short-name for stratified calibration of each binary covariate
STRAT = {
    "apoe":  ("APOE4 non-carrier", "APOE4 carrier", "apoe4"),
    "apoe2": ("APOE2 non-carrier", "APOE2 carrier", "apoe2"),
    "sex_M": ("Female", "Male", "sex"),
}


def render_variant_table(comp, outdirs, covstr):
    sections = [("time", "Time axis:   h(t | x)"), ("age", "Age axis (left-truncated):   h(a | x)")]
    headers = ["Model", "Covariates on", "k", "logLik", "AIC", "AIC type", "C-index", "IBS"]
    body = []
    for axis, title in sections:
        body.append(([title] + [""] * 7, "section"))
        for _, r in comp[comp.axis == axis].iterrows():
            body.append(([r.model, r.placement, str(int(r.n_params)), f"{r.logLik:.1f}", f"{r.AIC:.1f}",
                          r.AIC_kind, f"{r.c_index:.3f}", f"{r.IBS:.3f}"], "data"))
    fig, ax = plt.subplots(figsize=(12.6, 0.40 * (len(body) + 1) + 1.5)); ax.axis("off")
    ax.set_title(f"Covariate-placement variants — where do the covariates ({covstr}) act?\n"
                 "discontinuous Weibull→exp; covariate block on T / γ₁ (early level) / γ₂ (post-change level)."
                 "   C-index higher better, IBS lower better; AIC full vs partial NOT comparable",
                 fontsize=10, fontweight="bold", loc="left", pad=12)
    tbl = ax.table(cellText=[c for c, _ in body], colLabels=headers, loc="center", cellLoc="center",
                   colWidths=[0.26, 0.14, 0.045, 0.09, 0.09, 0.095, 0.085, 0.075])
    tbl.auto_set_font_size(False); tbl.set_fontsize(8.5); tbl.scale(1, 1.35)
    for j in range(len(headers)):
        cc = tbl[(0, j)]; cc.set_facecolor("#222"); cc.set_text_props(color="white", weight="bold")
    ridx = 1
    for axis, _t in sections:
        for j in range(len(headers)):
            cc = tbl[(ridx, j)]; cc.set_facecolor("#d9e3f0"); cc.set_text_props(weight="bold", style="italic")
        ridx += 1
        sub = comp[comp.axis == axis]
        best_c, best_ibs = sub.c_index.max(), sub.IBS.min()
        for _, r in sub.iterrows():
            if abs(r.c_index - best_c) < 1e-9:
                tbl[(ridx, 6)].set_text_props(weight="bold")
            if abs(r.IBS - best_ibs) < 1e-9:
                tbl[(ridx, 7)].set_text_props(weight="bold")
            ridx += 1
    for d in outdirs:
        for ext in ("png", "pdf"):
            fig.savefig(d / f"variant_comparison_table.{ext}", bbox_inches="tight", dpi=200)
    plt.close()


def plot_soft_s_grid(c, covcols, covlabels, s_values, wdir, cdir):
    """2×2 grid (APOE4 × Sex): observed KM vs constrained (continuity), free-s soft,
    and 3 fixed-s soft change-point models. Covariates on T."""
    exitg = c["exit_age"].values; entry = c["entry_age"].values; ev = c["event"].values
    X = np.column_stack([c[cc].values for cc in covcols] + [np.ones(len(c))])
    tb = (float(np.quantile(exitg, 0.05)), float(np.quantile(exitg, 0.95)))
    ref = 65.0; a_hi = float(np.quantile(exitg, 0.92)); grid = np.linspace(ref, a_hi, 30)
    cont = fit_changepoint(exitg, ev, X, entry=entry, t_bounds=tb)
    S_cont = np.column_stack([cont["predict_S"](X, a, ref) for a in grid])
    free = fit_soft_cov(exitg, ev, entry, X, None, tb)
    S_free = free["predict_S"](grid, ref)
    fixed = []
    for s in s_values:
        m = fit_soft_cov(exitg, ev, entry, X, s, tb)
        fixed.append((s, m, m["predict_S"](grid, ref)))
    apoe = c["apoe"].values; sexm = c["sex_M"].values
    s_colors = ["#abd9e9", "#4575b4", "#08306b"]
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), sharex=True, sharey=True)
    for r, (ga, alab) in enumerate([(0.0, "APOE4 non-carrier"), (1.0, "APOE4 carrier")]):
        for col, (gs, slab) in enumerate([(0.0, "Female"), (1.0, "Male")]):
            ax = axes[r][col]; mm = (apoe == ga) & (sexm == gs)
            if mm.sum() == 0:
                ax.set_visible(False); continue
            xs, ss = _km_lt_cond(entry[mm], exitg[mm], ev[mm], ref)
            ax.step(xs, ss, where="post", color="black", lw=2.2, label="KM (observed)")
            ax.plot(grid, S_cont[mm].mean(0), color="#d7191c", lw=2, ls="--", label="constrained (continuity)")
            ax.plot(grid, S_free[mm].mean(0), color="#1a9641", lw=2, ls="-.", label=f"free s (s={free['s']:.1f})")
            for (s, m, S), scol in zip(fixed, s_colors):
                ax.plot(grid, S[mm].mean(0), color=scol, lw=1.7, label=f"soft s={s:g}")
            ax.set_title(f"{alab} · {slab} (n={int(mm.sum())}, ev={int(ev[mm].sum())})",
                         fontsize=10, fontweight="bold")
            ax.set_ylim(0, 1.02); ax.grid(alpha=0.3)
            ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
            if r == 1:
                ax.set_xlabel("Age (years)", fontweight="bold")
            if col == 0:
                ax.set_ylabel("Survival (AD-free)", fontweight="bold")
    axes[0][1].legend(fontsize=8.5, loc="upper right")
    fig.suptitle("Soft change-point — effect of smoothing width s, vs free-s and constrained (continuity)\n"
                 "age axis (cond 65); covariates APOE4+APOE2+Sex on T; stratified by APOE4 × Sex",
                 fontsize=12, fontweight="bold")
    cap = ("Soft: h(a|x)=γ₂+(p·γ₁·a^(p−1)−γ₂)·σ((T−a)/s),  γ₂=p·γ₁·T^(p−1),  log T=b₀+Σc·X.   "
           "Lower AIC = better fit (all parametric, comparable):\n"
           f"   constrained (continuity, s→0 sharp kink):  AIC={cont['aic']:.1f}\n"
           f"   free s (fitted s={free['s']:.2f}):  AIC={free['aic']:.1f}\n"
           + "".join(f"   soft s={s:g} (fixed):  AIC={m['aic']:.1f}\n" for s, m, _ in fixed))
    fig.subplots_adjust(bottom=0.16, top=0.90)
    fig.text(0.012, 0.01, cap.rstrip(), ha="left", va="bottom", fontsize=8, family="DejaVu Sans", linespacing=1.5)
    for d in (wdir, cdir):
        for ext in ("png", "pdf"):
            fig.savefig(d / "soft_s_effect_by_apoe4_sex.{}".format(ext), bbox_inches="tight", dpi=300)
    plt.close()
    print("[soft-s] AIC: continuity=%.1f  free(s=%.2f)=%.1f  " % (cont["aic"], free["s"], free["aic"])
          + " ".join(f"s={s:g}->{m['aic']:.1f}" for s, m, _ in fixed))


# ── Run one covariate combination ───────────────────────────────────────────────
def run_combo(c, name, covcols, covlabels, strat_cols):
    wdir = WEIBULL_BASE / name; cdir = COX_BASE / name
    wdir.mkdir(parents=True, exist_ok=True); cdir.mkdir(parents=True, exist_ok=True)
    covstr = "+".join(covlabels)
    dur_t, ev = c["t"].values, c["event"].values
    exitg, entry = c["exit_age"].values, c["entry_age"].values
    rows = []

    # ---- TIME axis ----
    Xt = np.column_stack([c["age_z"].values] + [c[cc].values for cc in covcols] + [np.ones(len(c))])
    wt = fit_changepoint(dur_t, ev, Xt, entry=None, t_bounds=(0.5, float(0.9 * dur_t.max())))
    t_hi = float(np.quantile(dur_t, 0.92)); grid_t = np.linspace(0.5, t_hi, 24); HOR_T = min(5.0, t_hi)
    St_w = np.column_stack([wt["predict_S"](Xt, t0, 0.0) for t0 in grid_t])
    ctv, St_c = cox_time_fit_predict(c, covcols, grid_t)
    rows.append(dict(axis="time", covariates=covstr, model="Weibull change-point", n_params=wt["k"],
                     logLik=round(wt["loglik"], 2), AIC=round(wt["aic"], 2), AIC_kind="full",
                     c_index=round(cindex_at(wt["predict_S"](Xt, HOR_T, 0.0), dur_t, ev), 4),
                     IBS=round(ipcw_ibs(St_w, dur_t, ev, grid_t), 4)))
    rows.append(dict(axis="time", covariates=covstr, model="Cox PH (age×log t)",
                     n_params=int(ctv.params_.shape[0]), logLik=round(float(ctv.log_likelihood_), 2),
                     AIC=round(float(ctv.AIC_partial_), 2), AIC_kind="partial",
                     c_index=round(cindex_at(St_c[:, np.argmin(np.abs(grid_t - HOR_T))], dur_t, ev), 4),
                     IBS=round(ipcw_ibs(St_c, dur_t, ev, grid_t), 4)))

    # ---- AGE axis ----
    Xa = np.column_stack([c[cc].values for cc in covcols] + [np.ones(len(c))])
    wa = fit_changepoint(exitg, ev, Xa, entry=entry,
                         t_bounds=(float(np.quantile(exitg, 0.05)), float(np.quantile(exitg, 0.95))))
    a_lo = float(np.quantile(entry, 0.10)); a_hi = float(np.quantile(exitg, 0.92))
    grid_a = np.linspace(a_lo, a_hi, 24); HOR_A = min(80.0, a_hi)
    Sa_w = np.column_stack([wa["predict_S"](Xa, a, entry) for a in grid_a])
    cph_a, Sa_c, fine_a, M_a = cox_age_fit_predict(c, covcols, grid_a)
    rows.append(dict(axis="age", covariates=covstr, model="Weibull change-point", n_params=wa["k"],
                     logLik=round(wa["loglik"], 2), AIC=round(wa["aic"], 2), AIC_kind="full",
                     c_index=round(cindex_at(wa["predict_S"](Xa, HOR_A, entry), exitg, ev), 4),
                     IBS=round(ipcw_ibs(Sa_w, exitg, ev, grid_a, entry=entry), 4)))
    rows.append(dict(axis="age", covariates=covstr, model="Cox PH",
                     n_params=int(cph_a.params_.shape[0]), logLik=round(float(cph_a.log_likelihood_), 2),
                     AIC=round(float(cph_a.AIC_partial_), 2), AIC_kind="partial",
                     c_index=round(cindex_at(Sa_c[:, np.argmin(np.abs(grid_a - HOR_A))], exitg, ev), 4),
                     IBS=round(ipcw_ibs(Sa_c, exitg, ev, grid_a, entry=entry), 4)))

    # ---- covariate-placement variants (T / T+γ₁ / T+γ₂ / all / γ₂) ----
    vt_rows, vt_models = fit_variants(c, "time", dur_t, ev, None, grid_t, HOR_T,
                                      (0.5, float(0.9 * dur_t.max())), covcols)
    va_rows, va_models = fit_variants(c, "age", exitg, ev, entry, grid_a, HOR_A,
                                      (float(np.quantile(exitg, 0.05)), float(np.quantile(exitg, 0.95))), covcols)
    var = ([dict(axis="time", model="Weibull-cp: T only (continuity)", placement="T (γ₂ tied)",
                 n_params=wt["k"], logLik=round(wt["loglik"], 2), AIC=round(wt["aic"], 2),
                 AIC_kind="full", c_index=rows[0]["c_index"], IBS=rows[0]["IBS"])]
           + vt_rows
           + [dict(axis="time", model="Cox PH (age×log t)", placement="linear pred.",
                   n_params=rows[1]["n_params"], logLik=rows[1]["logLik"], AIC=rows[1]["AIC"],
                   AIC_kind="partial", c_index=rows[1]["c_index"], IBS=rows[1]["IBS"])]
           + [dict(axis="age", model="Weibull-cp: T only (continuity)", placement="T (γ₂ tied)",
                   n_params=wa["k"], logLik=round(wa["loglik"], 2), AIC=round(wa["aic"], 2),
                   AIC_kind="full", c_index=rows[2]["c_index"], IBS=rows[2]["IBS"])]
           + va_rows
           + [dict(axis="age", model="Cox PH", placement="linear pred.",
                   n_params=rows[3]["n_params"], logLik=rows[3]["logLik"], AIC=rows[3]["AIC"],
                   AIC_kind="partial", c_index=rows[3]["c_index"], IBS=rows[3]["IBS"])])
    _plot_variant_cal("time", c, wt, Xt, vt_models, None, wdir, cdir)
    _plot_variant_cal("age", c, wa, Xa, va_models, 65.0, wdir, cdir)
    var_df = pd.DataFrame(var)
    var_df.to_csv(wdir / "variant_comparison.csv", index=False)
    var_df.to_csv(cdir / "variant_comparison.csv", index=False)
    render_variant_table(var_df, (wdir, cdir), covstr)

    # ---- equations + figures ----
    eq_t = eq_weibull(wt, "time", covlabels) + "\n" + eq_cox_time(ctv, covcols, covlabels)
    eq_a = eq_weibull(wa, "age", covlabels) + "\n" + eq_cox_age(cph_a, covcols, covlabels)
    _plot_weibull_haz_surv(wt, eq_t, "time", covcols, c, wdir)
    _plot_weibull_haz_surv(wa, eq_a, "age", covcols, c, wdir)
    # 'all' variant (covariates on T, γ₁ AND γ₂) — scales the hazard like Cox, so it
    # SEPARATES strata; use it for the stratified calibration (overall stays continuity).
    all_t, all_a = vt_models[3], va_models[3]            # ("T + γ₁ + γ₂", res, D)
    St_w_all = np.column_stack([all_t[1]["predict_S"](all_t[2], t0, 0.0) for t0 in grid_t])
    eq_t_var = ("Weibull = 'T+γ₁+γ₂' variant (covariates on the change-point AND both hazard "
                "levels → separates strata)\n" + eq_cox_time(ctv, covcols, covlabels))
    eq_a_var = ("Weibull = 'T+γ₁+γ₂' variant (covariates on the change-point AND both hazard "
                "levels → separates strata)\n" + eq_cox_age(cph_a, covcols, covlabels))

    _plot_cal_time(grid_t, St_w, St_c, dur_t, ev, eq_t, wdir, cdir, strat=None)
    for ref in (56.0, 65.0):
        _plot_cal_age(c, wa, Xa, fine_a, M_a, eq_a, wdir, cdir, ref, strat=None)
    for cc in strat_cols:                               # stratify ONLY by covariates new to this combo
        strat = (c[cc].values, *STRAT[cc])
        _plot_cal_time(grid_t, St_w_all, St_c, dur_t, ev, eq_t_var, wdir, cdir, strat=strat)
        _plot_km_strat(c, "time", strat, wdir, cdir)    # clean observed KM by covariate
        _plot_variant_cal_strat("time", c, wt, Xt, vt_models, strat, None, wdir, cdir)
        _plot_variant_cal_strat("age", c, wa, Xa, va_models, strat, 65.0, wdir, cdir)
        for ref in (56.0, 65.0):
            _plot_cal_age(c, all_a[1], all_a[2], fine_a, M_a, eq_a_var, wdir, cdir, ref, strat=strat)
            _plot_km_strat(c, "age", strat, wdir, cdir, ref=ref)

    if "sex_M" in covcols:                               # 2×2 APOE4×Sex variant grid (γ₁ / γ₂ / all)
        _plot_variant_grid_apoe_sex(c, wa, Xa, va_models, eq_cox_age(cph_a, covcols, covlabels), wdir, cdir)

    cmp = pd.DataFrame(rows)
    cmp.to_csv(wdir / "model_comparison.csv", index=False)
    cmp.to_csv(cdir / "model_comparison.csv", index=False)
    print(f"\n[{name}]  covariates = {covstr}")
    print(cmp[["axis", "model", "n_params", "logLik", "AIC", "AIC_kind", "c_index", "IBS"]].to_string(index=False))
    return rows


def main():
    c = load_cohort()
    print(f"Pooled CN+MCI → AD: N={len(c)}, events={int(c['event'].sum())}, "
          f"APOE4={int(c['apoe'].sum())}, APOE2={int(c['apoe2'].sum())}, male={int(c['sex_M'].sum())}")
    kink_experiment(c, WEIBULL_BASE)        # hard vs current vs soft change-point transition
    all_rows = []
    seen = set()
    for name, covcols, covlabels in COMBOS:
        new = [cc for cc in covcols if cc not in seen]   # only stratify by covariates new here
        seen.update(covcols)
        all_rows += run_combo(c, name, covcols, covlabels, new)
    comp = pd.DataFrame(all_rows)
    comp.to_csv(WEIBULL_BASE / "model_comparison_master.csv", index=False)
    comp.to_csv(COX_BASE / "model_comparison_master.csv", index=False)
    render_master_table(comp, (WEIBULL_BASE, COX_BASE))
    plot_soft_s_grid(c, ["apoe", "apoe2", "sex_M"], ["APOE4", "APOE2", "Sex(M)"], [0.5, 3.0, 8.0],
                     WEIBULL_BASE / "apoe4_apoe2_sex", COX_BASE / "apoe4_apoe2_sex")
    print(f"\nDone. Per-combo outputs under {WEIBULL_BASE}\\<combo> and {COX_BASE}\\<combo>.")


if __name__ == "__main__":
    main()
