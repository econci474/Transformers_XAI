"""
02l_survival_model_comparison.py
================================
Head-to-head survival models on the SAME compact tabular features, fit per-seed on TRAIN and
evaluated on val / test / val+test (same framework as 02k), for the pooled CN+MCI → AD endpoint
(event = pCN_to_AD or pMCI; time = years_to_AD if event else FU_years).

Models (all on the identical median-imputed + standardised 8 features):
  • RSF              — Random Survival Forest (scikit-survival; config B: 1500 trees, leaf 5, m=3)
  • Cox PH           — semiparametric proportional hazards (lifelines, ridge penalizer)
  • Weibull AFT      — Weibull accelerated-failure-time (lifelines)
  • Weibull piecewise — Weibull→exponential change-point PH baseline with a THRESHOLD T:
        H(t|x) = H₀(t)·exp(βx),  H₀(t) = γ₁·min(t,T)^p + γ₂·max(t−T,0),  γ₂ = p·γ₁·T^(p−1)
        (continuity-constrained baseline; covariates proportional; threshold T a fitted parameter)

Unified, model-agnostic evaluation: risk = 1 − S(t_ref) (so RSF/Cox/AFT/piecewise are scored
identically) → Harrell C-index, Uno IPCW C-index, integrated Brier score, time-dependent AUC@3y/5y.

Env: conda env `survml` (scikit-survival from conda-forge + lifelines via pip):
    conda create -n survml -c conda-forge python=3.11 scikit-survival pandas matplotlib
    conda run -n survml python -m pip install lifelines
Run:
    conda run -n survml python clinical_pipeline/02l_survival_model_comparison.py

Outputs → clinical_pipeline/outputs/model_comparison/:
    master_comparison.csv (per model × seed × split), master_comparison_summary.csv,
    master_comparison_table.{png,pdf}, calibration_valtest_models.{png,pdf}.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import minimize

from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import (concordance_index_censored, concordance_index_ipcw,
                            integrated_brier_score, cumulative_dynamic_auc)
from sksurv.util import Surv
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from lifelines import CoxPHFitter, WeibullAFTFitter

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
matplotlib.rcParams["font.family"] = "DejaVu Serif"

REPO_ROOT = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
SPLIT_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion"
                 r"\tabular\baseline")
META_DIR = Path(r"D:\ADNI_SNP_Omni2.5M_20140220\source_prs")  # per-seed/split meta-PRS-EN
OUT_DIR = REPO_ROOT / "clinical_pipeline" / "outputs" / "model_comparison"
OUT_EXT = REPO_ROOT / "clinical_pipeline" / "outputs" / "model_comparison_extended"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_EXT.mkdir(parents=True, exist_ok=True)

SEEDS = (0, 1, 2)
EPS = 0.01
# Compact set (original 8) and extended set (+ meta-PRS, CSF, plasma pTau181, memory/fluency)
COMPACT_FEATURES = ["MMSE_Total", "ADAS_Total13", "MoCA_Total", "FAQ_Total",
                    "age_at_baseline", "sex_M", "Education_Years", "APOE4_Dosage"]
EXTENDED_FEATURES = COMPACT_FEATURES + [
    "meta_prs_EN_combined", "CSF_Abeta_bl", "CSF_Tau_bl", "CSF_pTau_bl",
    "Plasma_pTau181", "CategoryFluency_Animals_Total", "AVLT_ListB", "AVLT_Delay30"]
MODELS = ["RSF", "Cox PH", "Weibull AFT", "Weibull piecewise"]


# ── Data prep (feature set median-imputed + standardised on each seed's TRAIN) ────
def _build(df, features):
    df = df[df["bl_dx"].isin(["CN", "MCI"])].copy()
    ev = df["conversion_group"].isin(["pCN_to_AD", "pMCI"]).astype(bool).values
    t = np.where(ev, pd.to_numeric(df["years_to_AD"], errors="coerce"),
                 pd.to_numeric(df["FU_years"], errors="coerce"))
    t = pd.Series(t).fillna(EPS).clip(lower=EPS).values
    df["Date"] = pd.to_datetime(df["Date"], errors="coerce")
    birth = pd.to_datetime(pd.DataFrame(
        {"year": df["YOB"].astype(int), "month": 7, "day": 1}).set_index(df.index))
    df["age_at_baseline"] = (df["Date"] - birth) / np.timedelta64(1, "D") / 365.25
    df["sex_M"] = (df["Sex"].astype(str) == "Male").astype(float)
    return df[features].apply(pd.to_numeric, errors="coerce"), ev, t


def _merge_meta_prs(df, seed, split):
    """Attach the per-seed/split meta-PRS-EN-combined (train-fitted EN, leakage-safe)."""
    m = pd.read_csv(META_DIR / f"meta_prs_en_combined_seed{seed}_{split}.tsv", sep="\t")
    m["Patient_ID"] = m["Patient_ID"].astype(str)
    df = df.copy(); df["Patient_ID"] = df["Patient_ID"].astype(str)
    return df.merge(m, on="Patient_ID", how="left")


def load_seed(seed, features):
    sd = SPLIT_DIR / f"seed_{seed}"
    parts = {}
    for sp in ("train", "val", "test"):
        df = pd.read_csv(sd / f"{sp}.csv", low_memory=False)
        if "meta_prs_EN_combined" in features:
            df = _merge_meta_prs(df, seed, sp)
        parts[sp] = df
    Xtr, _, _ = _build(parts["train"], features)
    imp = SimpleImputer(strategy="median").fit(Xtr)
    scl = StandardScaler().fit(imp.transform(Xtr))
    out = {}
    for sp in ("train", "val", "test"):
        X, e, t = _build(parts[sp], features)
        Xs = pd.DataFrame(scl.transform(imp.transform(X)), columns=features, index=X.index)
        out[sp] = (Xs, e, t)
    return out


# ── Weibull piecewise-PH (continuity-constrained change-point baseline) ──────────
def _H0(theta, tt):
    p, g1, T = np.exp(theta[0]), np.exp(theta[1]), np.exp(theta[2])
    g2 = p * g1 * T ** (p - 1.0)
    tt = np.asarray(tt, float)
    return g1 * np.minimum(tt, T) ** p + g2 * np.maximum(tt - T, 0.0)


def _h0(theta, tt):
    p, g1, T = np.exp(theta[0]), np.exp(theta[1]), np.exp(theta[2])
    g2 = p * g1 * T ** (p - 1.0)
    tt = np.asarray(tt, float)
    return np.where(tt < T, g1 * p * np.clip(tt, 1e-8, None) ** (p - 1.0), g2)


def fit_wpph(t, e, X, ridge=0.02, n_starts=12):
    t = np.asarray(t, float); e = np.asarray(e, float); X = np.asarray(X, float)
    n, k = X.shape
    ev = t[e == 1]
    t_lo, t_hi = max(0.5, float(np.percentile(ev, 5))), float(np.percentile(ev, 95))
    lt_lo, lt_hi = np.log(t_lo), np.log(t_hi); lp_lo, lp_hi = np.log(0.3), np.log(6.0)

    def nll(theta):
        if not np.all(np.isfinite(theta)) or not (lp_lo <= theta[0] <= lp_hi) or not (lt_lo <= theta[2] <= lt_hi):
            return 1e12
        beta = theta[3:]; lp = X @ beta
        if not np.all(np.isfinite(lp)) or np.abs(lp).max() > 30:
            return 1e12
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            ll = np.sum(e * (np.log(np.clip(_h0(theta, t), 1e-300, None)) + lp)
                        - _H0(theta, t) * np.exp(lp)) - ridge * np.sum(beta ** 2)
        return -ll if np.isfinite(ll) else 1e12

    base = np.zeros(3 + k)
    base[1] = np.log(max(e.sum() / t.sum(), 1e-4))
    base[2] = np.clip(np.log(np.median(ev)), lt_lo, lt_hi)
    rng = np.random.default_rng(0); best = None
    for s in range(n_starts):
        x0 = base.copy()
        if s:
            x0[:3] += rng.normal(0, 0.3, 3); x0[0] = rng.uniform(lp_lo, lp_hi)
            x0[2] = np.clip(x0[2], lt_lo, lt_hi); x0[3:] += rng.normal(0, 0.2, k)
        try:
            res = minimize(nll, x0, method="Nelder-Mead", options={"maxiter": 20000, "xatol": 1e-6, "fatol": 1e-6})
        except Exception:
            continue
        if np.isfinite(res.fun) and (best is None or res.fun < best.fun):
            best = res
    return best.x


# ── Piecewise change-point VARIANTS: covariates placed on T / γ₁ / γ₂ / all ──────
#   h(t|x) = p·γ₁(x)·t^(p−1)  for t < T(x);   h(t|x) = γ₂(x)  for t ≥ T(x)
#   γ₁(x)=exp(a₁+β₁x), γ₂(x)=exp(a₂+β₂x), T(x)=exp(aT+βT x); γ₂ FREE (discontinuous).
#   A variant only attaches the SAME 8 covariates to the flagged parameter block(s).
PW_VARIANTS = {
    "T only":        (False, False, True),   # (γ₁, γ₂, T)
    "γ1 + T":        (True,  False, True),
    "γ2 + T":        (False, True,  True),
    "all (γ1,γ2,T)": (True,  True,  True),
}


def _pw_unpack(params, k, flags):
    fg1, fg2, fT = flags
    i = 4
    log_p, log_g1, log_g2, log_T = params[0], params[1], params[2], params[3]
    bg1 = params[i:i + k] if fg1 else np.zeros(k); i += k if fg1 else 0
    bg2 = params[i:i + k] if fg2 else np.zeros(k); i += k if fg2 else 0
    bT = params[i:i + k] if fT else np.zeros(k); i += k if fT else 0
    return log_p, log_g1, log_g2, log_T, bg1, bg2, bT


def _pw_comp(params, X, k, flags):
    log_p, log_g1, log_g2, log_T, bg1, bg2, bT = _pw_unpack(params, k, flags)
    p = np.exp(log_p)
    g1 = np.exp(np.clip(log_g1 + X @ bg1, -30, 30))
    g2 = np.exp(np.clip(log_g2 + X @ bg2, -30, 30))
    T = np.exp(np.clip(log_T + X @ bT, -30, 30))
    return p, g1, g2, T


def _pw_surv(params, X, times, k, flags):
    p, g1, g2, T = _pw_comp(params, X, k, flags)
    tt = np.asarray(times, float)[None, :]
    Tn, g1n, g2n = T[:, None], g1[:, None], g2[:, None]
    H = g1n * np.minimum(tt, Tn) ** p + g2n * np.maximum(tt - Tn, 0.0)
    return np.exp(-np.clip(H, 0, 700))


def fit_pw_variant(t, e, X, flags, ridge=0.05, n_starts=6):
    t = np.asarray(t, float); e = np.asarray(e, float); X = np.asarray(X, float)
    n, k = X.shape
    ev = t[e == 1]
    lt_lo, lt_hi = np.log(max(0.5, np.percentile(ev, 5))), np.log(np.percentile(ev, 95))
    lp_lo, lp_hi = np.log(0.3), np.log(6.0)
    fg1, fg2, fT = flags
    nb = 4 + k * (int(fg1) + int(fg2) + int(fT))

    def nll(params):
        if not np.all(np.isfinite(params)):
            return 1e12
        p, g1, g2, T = _pw_comp(params, X, k, flags)
        if not (lp_lo <= params[0] <= lp_hi):
            return 1e12
        early = t < T
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            h = np.where(early, g1 * p * np.clip(t, 1e-8, None) ** (p - 1.0), g2)
            H = g1 * np.minimum(t, T) ** p + g2 * np.maximum(t - T, 0.0)
            ll = np.sum(e * np.log(np.clip(h, 1e-300, None)) - H) - ridge * np.sum(params[4:] ** 2)
        return -ll if np.isfinite(ll) else 1e12

    base = np.zeros(nb)
    base[1] = base[2] = np.log(max(e.sum() / t.sum(), 1e-4))
    base[3] = np.clip(np.log(np.median(ev)), lt_lo, lt_hi)
    bnds = [(lp_lo, lp_hi), (-12, 4), (-12, 4), (lt_lo, lt_hi)] + [(-2.5, 2.5)] * (nb - 4)
    rng = np.random.default_rng(0); best = None
    for s in range(n_starts):
        x0 = base.copy()
        if s:
            x0[0] = rng.uniform(lp_lo, lp_hi)
            x0[1:4] += rng.normal(0, 0.3, 3); x0[3] = np.clip(x0[3], lt_lo, lt_hi)
            x0[4:] += rng.normal(0, 0.2, nb - 4)
        try:
            res = minimize(nll, x0, method="L-BFGS-B", bounds=bnds, options={"maxiter": 5000})
        except Exception:
            continue
        if np.isfinite(res.fun) and (best is None or res.fun < best.fun):
            best = res
    return best.x


def fit_pw_variant_model(flags, Xtr, etr, ttr):
    params = fit_pw_variant(ttr, etr, Xtr.values, flags)
    k = Xtr.shape[1]
    return lambda X, times: np.clip(_pw_surv(params, X.values, times, k, flags), 1e-8, 1.0)


# ── Model fitters → return a surv(X, times) function ────────────────────────────
def fit_model(name, Xtr, etr, ttr):
    if name == "RSF":
        rsf = RandomSurvivalForest(n_estimators=1500, min_samples_leaf=5, max_features=3,
                                   n_jobs=-1, random_state=42).fit(Xtr.values, Surv.from_arrays(etr, ttr))

        def _rsf_surv(X, times):
            fns = rsf.predict_survival_function(X.values, return_array=False)
            # StepFunction is flat beyond its last event time; clamp queries into its domain
            lo, hi = fns[0].domain
            tq = np.clip(np.asarray(times, float), lo, hi)
            return np.array([[float(fn(tt)) for tt in tq] for fn in fns])

        return _rsf_surv
    if name == "Cox PH":
        df = Xtr.copy(); df["_t"] = ttr; df["_e"] = etr.astype(int)
        cph = CoxPHFitter(penalizer=0.01).fit(df, duration_col="_t", event_col="_e")
        return lambda X, times: cph.predict_survival_function(X, times=list(times)).values.T
    if name == "Weibull AFT":
        df = Xtr.copy(); df["_t"] = ttr; df["_e"] = etr.astype(int)
        aft = WeibullAFTFitter(penalizer=0.01).fit(df, duration_col="_t", event_col="_e")
        return lambda X, times: aft.predict_survival_function(X, times=list(times)).values.T
    if name == "Weibull piecewise":
        theta = fit_wpph(ttr, etr, Xtr.values)
        return lambda X, times: np.exp(-np.outer(np.exp(X.values @ theta[3:]),
                                                  _H0(theta, np.asarray(times, float))))
    raise ValueError(name)


# ── Unified evaluation (risk = 1 − S(t_ref); identical for every model) ──────────
def _grid(ttr, etr, tev, eev):
    hi = min(ttr[etr].max() if etr.any() else ttr.max(),
             tev[eev].max() if eev.any() else tev.max(), tev.max(), ttr.max()) - 1e-3
    lo = max(0.5, float(np.percentile(tev[eev], 5)) if eev.any() else 0.5)
    lo = min(lo, hi * 0.5)
    return np.linspace(lo, hi, 12)


def evaluate(surv_fn, ytr, ttr, etr, X, e, t):
    yev = Surv.from_arrays(event=e, time=t)
    times = _grid(ttr, etr, t, e)
    t_ref = float(np.clip(4.0, times[0] + 1e-6, times[-1] - 1e-6))
    S = np.clip(surv_fn(X, times), 1e-8, 1.0)
    risk = 1.0 - np.clip(surv_fn(X, np.array([t_ref]))[:, 0], 1e-8, 1.0)
    out = {"n": int(len(t)), "events": int(e.sum())}
    out["c_index"] = float(concordance_index_censored(e, t, risk)[0])
    try:
        out["c_index_ipcw"] = float(concordance_index_ipcw(ytr, yev, risk, tau=times[-1])[0])
    except Exception:
        out["c_index_ipcw"] = np.nan
    try:
        out["ibs"] = float(integrated_brier_score(ytr, yev, S, times))
    except Exception:
        out["ibs"] = np.nan
    hmax = min(float(ttr.max()), float(t.max())) - 1e-3
    for h in (3.0, 5.0, 7.0, 10.0):
        out[f"auc_{int(h)}y"] = np.nan
        if (t[e] < h).any() and h < hmax:
            try:
                out[f"auc_{int(h)}y"] = float(np.atleast_1d(cumulative_dynamic_auc(ytr, yev, risk, h)[0])[0])
            except Exception:
                pass
    return out


def _km(t, e):
    d = pd.DataFrame({"t": t, "e": np.asarray(e).astype(int)}).sort_values("t")
    S, xs, ss = 1.0, [], []
    for tt in sorted(d["t"].unique()):
        n = int((d["t"] >= tt).sum()); dd = int(((d["t"] == tt) & (d["e"] == 1)).sum())
        if n > 0:
            S *= (1 - dd / n)
        xs.append(tt); ss.append(S)
    return np.array(xs), np.array(ss)


# ── Comparison table (shared by headline-model and piecewise-variant tables) ─────
def covar_note(features):
    return (f"Covariates ({len(features)}, identical for every model): " + ", ".join(features) +
            "  —  median-imputed & standardised on each seed's train fold.")


def render_table(met, models, row_label, title, fname, out_dir, note):
    mc = ["c_index", "c_index_ipcw", "ibs", "auc_3y", "auc_5y", "auc_7y", "auc_10y"]
    g = met.groupby([row_label, "split"])[mc].agg(["mean", "std"])
    nev = met.groupby([row_label, "split"])[["n", "events"]].mean()
    order = {"val": 0, "test": 1, "valtest": 2}
    headers = [row_label.capitalize(), "Held-out", "N (events)", "C-index", "IPCW C",
               "IBS", "AUC@3y", "AUC@5y", "AUC@7y", "AUC@10y"]
    best = {sp: {"c_index": met[met.split == sp].groupby(row_label).c_index.mean().max(),
                 "ibs": met[met.split == sp].groupby(row_label).ibs.mean().min()} for sp in order}
    body, bold = [], []
    for model in models:
        for split in sorted(order, key=lambda s: order[s]):
            r = g.loc[(model, split)]; ne = nev.loc[(model, split)]
            cell = lambda m: f"{r[(m, 'mean')]:.3f}±{r[(m, 'std')]:.3f}"
            body.append([model, split, f"{ne['n']:.0f} ({ne['events']:.0f})",
                         cell("c_index"), cell("c_index_ipcw"), cell("ibs"),
                         cell("auc_3y"), cell("auc_5y"), cell("auc_7y"), cell("auc_10y")])
            bold.append((abs(r[("c_index", "mean")] - best[split]["c_index"]) < 1e-9,
                         abs(r[("ibs", "mean")] - best[split]["ibs"]) < 1e-9))
    fig, ax = plt.subplots(figsize=(16.5, 0.46 * len(body) + 2.0)); ax.axis("off")
    ax.set_title(title, fontsize=10.5, fontweight="bold", loc="left", pad=12)
    tbl = ax.table(cellText=body, colLabels=headers, loc="center", cellLoc="center",
                   colWidths=[0.15, 0.09, 0.10, 0.115, 0.115, 0.10, 0.085, 0.085, 0.085, 0.085])
    tbl.auto_set_font_size(False); tbl.set_fontsize(8.2); tbl.scale(1, 1.45)
    for j in range(len(headers)):
        cc = tbl[(0, j)]; cc.set_facecolor("#222"); cc.set_text_props(color="white", weight="bold")
    for i, (b_c, b_i) in enumerate(bold, start=1):
        if body[i - 1][1] == "valtest":
            for j in range(len(headers)):
                tbl[(i, j)].set_facecolor("#eef4fb")
        if b_c:
            tbl[(i, 3)].set_text_props(weight="bold")
        if b_i:
            tbl[(i, 5)].set_text_props(weight="bold")
    fig.text(0.012, 0.012, note, fontsize=7.0, style="italic", color="#333", wrap=True)
    for ext in ("png", "pdf"):
        fig.savefig(out_dir / f"{fname}.{ext}", bbox_inches="tight", dpi=200)
    plt.close()


def plot_calibration(cal, cal_grid, pooled_te, cols, title, fname, out_dir, note):
    tt, ee = pooled_te
    tmax = float(np.ceil(tt.max()))
    # numbers-at-risk timepoints: ~10 evenly spaced ticks over the full follow-up
    step = max(1.0, round(tmax / 10.0))
    risk_t = np.arange(0, tmax + step / 2, step)
    at_risk = [int((tt >= rt).sum()) for rt in risk_t]

    fig, (ax, axr) = plt.subplots(2, 1, figsize=(10.5, 6.9),
                                  gridspec_kw={"height_ratios": [5, 1.0], "hspace": 0.08})
    xs, ss = _km(tt, ee)
    ax.step(np.append(xs, tmax), np.append(ss, ss[-1]), where="post",
            color="black", lw=2.4, label=f"Kaplan–Meier (observed, pooled val+test, n={len(tt)})")
    for model, c in cols.items():
        ax.plot(cal_grid, np.vstack(cal[model]).mean(axis=0), color=c, lw=2, label=model)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_ylabel("Survival (AD-free)", fontweight="bold")
    ax.set_xlim(0, tmax); ax.set_ylim(0, 1.02); ax.grid(alpha=0.3); ax.legend(fontsize=8.5)
    # YEAR tick labels live on the main plot's x-axis (directly under the curve)
    ax.set_xticks(risk_t)
    ax.set_xticklabels([f"{rt:.0f}" for rt in risk_t])
    ax.tick_params(axis="x", labelbottom=True, length=4, labelsize=9)
    ax.set_xlabel("Time (years)", fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    # numbers-at-risk strip BELOW the year axis (counts only; no overlapping label)
    axr.set_xlim(0, tmax); axr.set_ylim(0, 1)
    axr.text(-0.055 * tmax, 0.5, "At risk", ha="right", va="center",
             fontsize=8.5, fontweight="bold", clip_on=False)
    for rt, nr in zip(risk_t, at_risk):
        axr.text(rt, 0.5, str(nr), ha="center", va="center", fontsize=8.5)
    axr.set_xticks([]); axr.set_yticks([])
    for s in ("top", "left", "right", "bottom"):
        axr.spines[s].set_visible(False)
    fig.text(0.012, 0.005, note, fontsize=7.0, style="italic", color="#333", wrap=True)

    fig.tight_layout(rect=[0, 0.03, 1, 1])
    for ext in ("png", "pdf"):
        fig.savefig(out_dir / f"{fname}.{ext}", bbox_inches="tight", dpi=300)
    plt.close()


def _pooled_max():
    hi = 0.0
    for seed in SEEDS:
        sd = SPLIT_DIR / f"seed_{seed}"
        for sp in ("val", "test"):
            _, _, t = _build(pd.read_csv(sd / f"{sp}.csv", low_memory=False), ["age_at_baseline"])
            hi = max(hi, float(t.max()))
    return hi


MC = ["c_index", "c_index_ipcw", "ibs", "auc_3y", "auc_5y", "auc_7y", "auc_10y"]


def make_model_specs(names):
    """name -> fit_fn(Xtr,etr,ttr)->surv_fn. Includes headline models + named PW variants."""
    spec = {}
    for n in names:
        if n in ("RSF", "Cox PH", "Weibull AFT", "Weibull piecewise"):
            spec[n] = (lambda Xtr, etr, ttr, _n=n: fit_model(_n, Xtr, etr, ttr))
        elif n in PW_VARIANTS:
            spec[n] = (lambda Xtr, etr, ttr, _f=PW_VARIANTS[n]: fit_pw_variant_model(_f, Xtr, etr, ttr))
        else:
            raise ValueError(n)
    return spec


def run_comparison(features, out_dir, model_names, cols, title, do_variants):
    cal_grid = np.linspace(0.5, float(np.ceil(_pooled_max())), 60)
    specs = make_model_specs(model_names)
    rows, vrows = [], []
    cal = {m: [] for m in model_names}
    vcal = {v: [] for v in PW_VARIANTS}
    pooled_t, pooled_e = [], []
    for seed in SEEDS:
        d = load_seed(seed, features)
        Xtr, etr, ttr = d["train"]
        ytr = Surv.from_arrays(event=etr, time=ttr)
        d["valtest"] = (pd.concat([d["val"][0], d["test"][0]]),
                        np.concatenate([d["val"][1], d["test"][1]]),
                        np.concatenate([d["val"][2], d["test"][2]]))
        pooled_t.append(d["valtest"][2]); pooled_e.append(d["valtest"][1])
        Xvt = d["valtest"][0]
        for name in model_names:
            surv_fn = specs[name](Xtr, etr, ttr)
            for sp in ("val", "test", "valtest"):
                X, e, t = d[sp]
                rows.append({"model": name, "seed": seed, "split": sp,
                             **evaluate(surv_fn, ytr, ttr, etr, X, e, t)})
            cal[name].append(np.clip(surv_fn(Xvt, cal_grid), 0, 1).mean(axis=0))
            print(f"[{out_dir.name}][seed {seed}] {name:20s} valtest "
                  f"C={rows[-1]['c_index']:.3f} IBS={rows[-1]['ibs']:.3f}")
        if do_variants:
            for vname, flags in PW_VARIANTS.items():
                surv_fn = fit_pw_variant_model(flags, Xtr, etr, ttr)
                for sp in ("val", "test", "valtest"):
                    X, e, t = d[sp]
                    vrows.append({"variant": vname, "seed": seed, "split": sp,
                                  **evaluate(surv_fn, ytr, ttr, etr, X, e, t)})
                vcal[vname].append(np.clip(surv_fn(Xvt, cal_grid), 0, 1).mean(axis=0))

    pooled = (np.concatenate(pooled_t), np.concatenate(pooled_e))
    note = covar_note(features)

    met = pd.DataFrame(rows)
    met.to_csv(out_dir / "master_comparison.csv", index=False)
    summ = met.groupby(["model", "split"])[MC].agg(["mean", "std"]).round(4)
    summ.to_csv(out_dir / "master_comparison_summary.csv")
    print(f"\n[{out_dir.name}] Master comparison (mean ± std over seeds):\n{summ.to_string()}")
    render_table(met, model_names, "model", title, "master_comparison_table", out_dir, note)
    plot_calibration(cal, cal_grid, pooled, cols,
                     "Survival models — calibration on val+test (mean predicted Ŝ vs observed KM)\n"
                     f"CN+MCI → AD  ({len(features)} features)",
                     "calibration_valtest_models", out_dir, note)

    if do_variants:
        vmet = pd.DataFrame(vrows)
        vmet.to_csv(out_dir / "piecewise_variant_comparison.csv", index=False)
        vmet.groupby(["variant", "split"])[MC].agg(["mean", "std"]).round(4).to_csv(
            out_dir / "piecewise_variant_comparison_summary.csv")
        render_table(
            vmet, list(PW_VARIANTS), "variant",
            "Piecewise Weibull→exponential change-point — covariate-placement variants (CN+MCI → AD)\n"
            "h(t|x)=p·γ₁(x)·t^(p−1) for t<T(x), =γ₂(x) for t≥T(x); each variant attaches the SAME covariates\n"
            "to the flagged block(s). Fit on train, mean ± std over 3 seeds; bold = best per split (max C-index / min IBS)",
            "piecewise_variant_comparison_table", out_dir, note)
        plot_calibration(
            vcal, cal_grid, pooled,
            {"T only": "#1f78b4", "γ1 + T": "#33a02c", "γ2 + T": "#ff7f00", "all (γ1,γ2,T)": "#e7298a"},
            "Piecewise change-point variants — calibration on val+test (covariates on T / γ₁ / γ₂ / all)\nCN+MCI → AD",
            "calibration_valtest_piecewise_variants", out_dir, note)


def main():
    # 1) Compact 8-feature comparison (RSF / Cox / Weibull AFT / Weibull piecewise) + variant table
    run_comparison(
        COMPACT_FEATURES, OUT_DIR, MODELS,
        {"RSF": "#1b9e77", "Cox PH": "#d95f02", "Weibull AFT": "#7570b3", "Weibull piecewise": "#e7298a"},
        "Survival model comparison — CN+MCI → AD  (8 compact features; fit on train, mean ± std over 3 seeds)\n"
        "risk = 1 − S(t_ref) for all models; bold = best per held-out split (max C-index / min IBS)",
        do_variants=True)

    # 2) Extended 16-feature comparison — best two PW variants (γ₁+T, γ₂+T) vs RSF / Cox / Weibull AFT
    ext_models = ["RSF", "Cox PH", "Weibull AFT", "γ1 + T", "γ2 + T"]
    run_comparison(
        EXTENDED_FEATURES, OUT_EXT, ext_models,
        {"RSF": "#1b9e77", "Cox PH": "#d95f02", "Weibull AFT": "#7570b3",
         "γ1 + T": "#33a02c", "γ2 + T": "#ff7f00"},
        "Survival model comparison — CN+MCI → AD  (16 extended features incl. meta-PRS-EN, CSF, plasma pTau181)\n"
        "Piecewise variants = best from variant comparison (γ₁+T, γ₂+T); fit on train, mean ± std over 3 seeds;\n"
        "risk = 1 − S(t_ref) for all models; bold = best per held-out split (max C-index / min IBS)",
        do_variants=False)

    print(f"\nDone. Outputs → {OUT_DIR}  and  {OUT_EXT}")


if __name__ == "__main__":
    main()
