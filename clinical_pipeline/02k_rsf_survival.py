"""
02k_rsf_survival.py
===================
Random Survival Forest on the tabular dataset — per-seed train → val/test — for the pooled
**CN+MCI → AD** endpoint (event = pCN_to_AD or pMCI; time = years_to_AD if event else FU_years,
floored at 0.01), time-from-baseline.

This is the first model in a per-seed train/val/test survival-ML framework; `load_seed()` and the
metric helpers are written to be REUSED by the parametric model (02j Weibull change-point / Cox)
for a head-to-head RSF-vs-parametric comparison.

Env: dedicated conda env `survml` — scikit-survival is not pip-installable on Windows without a
C compiler (ecos has no wheel), so use conda-forge prebuilt binaries:
    conda create -n survml -c conda-forge python=3.11 scikit-survival pandas matplotlib
Run:
    conda run -n survml python clinical_pipeline/02k_rsf_survival.py
    (or:  C:\\Users\\elena\\miniforge3\\envs\\survml\\python.exe clinical_pipeline\\02k_rsf_survival.py)

Data : D:\\ADNI_BIDS_project\\derivatives\\clinical\\no_cdr_stratified_post_exclusion
       \\tabular\\baseline\\seed_{0,1,2}\\{train,val,test}.csv
Cohort: baseline CN+MCI (bl_dx ∈ {CN,MCI}); N = 427/53/55 train/val/test per seed.
Features (compact, median-imputed on TRAIN, applied to val/test; RSF needs no scaling):
   MMSE_Total, ADAS_Total13, MoCA_Total, FAQ_Total, age_at_baseline (Date−YOB), sex_M,
   Education_Years, APOE4_Dosage.
Eval (val AND test; fit on train only): Harrell C-index, Uno IPCW C-index, integrated Brier
   score, time-dependent AUC at 3 y / 5 y. Held-out event counts are small (≈10–16) ⇒ noisy.

Outputs → clinical_pipeline/outputs/rsf/:
   rsf_metrics.csv (per seed × split), rsf_metrics_summary.csv (mean±std over seeds),
   rsf_feature_importance.{csv,png}, calibration_val.{png,pdf}.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import (concordance_index_censored, concordance_index_ipcw,
                            integrated_brier_score, cumulative_dynamic_auc)
from sksurv.util import Surv
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

matplotlib.rcParams["font.family"] = "DejaVu Serif"

REPO_ROOT = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
SPLIT_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion"
                 r"\tabular\baseline")
OUT_DIR = REPO_ROOT / "clinical_pipeline" / "outputs" / "rsf"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SEEDS = (0, 1, 2)
EPS = 0.01
FEATURES = ["MMSE_Total", "ADAS_Total13", "MoCA_Total", "FAQ_Total",
            "age_at_baseline", "sex_M", "Education_Years", "APOE4_Dosage"]
# RSF hyperparameter configs are defined in main(): A = previous, B = paper.
# (scikit-survival RSF uses the log-rank split rule by default; it has NO "number of
#  split points"/nsplit parameter — all candidate split points are evaluated.)


# ── Data prep (reusable by the parametric model) ────────────────────────────────
def _build(df):
    """Filter to CN+MCI, build (X, event, time) for the pooled CN+MCI→AD endpoint."""
    df = df[df["bl_dx"].isin(["CN", "MCI"])].copy()
    ev = df["conversion_group"].isin(["pCN_to_AD", "pMCI"]).astype(bool).values
    yad = pd.to_numeric(df["years_to_AD"], errors="coerce")
    fu = pd.to_numeric(df["FU_years"], errors="coerce")
    t = np.where(ev, yad, fu)
    t = pd.Series(t).fillna(EPS).clip(lower=EPS).values
    df["Date"] = pd.to_datetime(df["Date"], errors="coerce")
    birth = pd.to_datetime(pd.DataFrame(
        {"year": df["YOB"].astype(int), "month": 7, "day": 1}).set_index(df.index))
    df["age_at_baseline"] = (df["Date"] - birth) / np.timedelta64(1, "D") / 365.25
    df["sex_M"] = (df["Sex"].astype(str) == "Male").astype(float)
    X = df[FEATURES].apply(pd.to_numeric, errors="coerce")
    return X, ev, t


def load_seed(seed):
    """Returns {split: (X_imputed, y_struct, event, time)} for train/val/test of `seed`."""
    sd = SPLIT_DIR / f"seed_{seed}"
    parts = {sp: pd.read_csv(sd / f"{sp}.csv", low_memory=False) for sp in ("train", "val", "test")}
    Xtr, _, _ = _build(parts["train"])
    imp = SimpleImputer(strategy="median").fit(Xtr)
    out = {}
    for sp in ("train", "val", "test"):
        X, e, t = _build(parts[sp])
        Xi = pd.DataFrame(imp.transform(X), columns=FEATURES, index=X.index)
        out[sp] = (Xi, Surv.from_arrays(event=e, time=t), e, t)
    return out


# ── Metric helpers ──────────────────────────────────────────────────────────────
def _grid(ttr, etr, tev, eev):
    """Time grid for IBS/AUC, kept strictly inside the observed (train ∩ eval) event range."""
    hi = min(ttr[etr].max() if etr.any() else ttr.max(),
             tev[eev].max() if eev.any() else tev.max(), tev.max(), ttr.max()) - 1e-3
    lo = max(0.5, float(np.percentile(tev[eev], 5)) if eev.any() else 0.5)
    lo = min(lo, hi * 0.5)
    return np.linspace(lo, hi, 12)


def _surv_at(rsf, X, times):
    """Predicted survival S(t|x) at `times` via the per-subject step functions."""
    fns = rsf.predict_survival_function(X, return_array=False)     # array of StepFunction
    return np.array([[float(fn(t)) for t in times] for fn in fns])  # (n, len(times))


def eval_split(rsf, y_train, ttr, etr, X, y, e, t):
    risk = rsf.predict(X)                                          # higher = higher risk
    times = _grid(ttr, etr, t, e)
    out = {"n": int(len(t)), "events": int(e.sum())}
    out["c_index"] = float(concordance_index_censored(e, t, risk)[0])
    try:
        out["c_index_ipcw"] = float(concordance_index_ipcw(y_train, y, risk, tau=times[-1])[0])
    except Exception:
        out["c_index_ipcw"] = np.nan
    try:
        out["ibs"] = float(integrated_brier_score(y_train, y, _surv_at(rsf, X, times), times))
    except Exception:
        out["ibs"] = np.nan
    for h in (3.0, 5.0):
        out[f"auc_{int(h)}y"] = np.nan
        if times[0] < h < times[-1]:
            try:
                a, _ = cumulative_dynamic_auc(y_train, y, risk, h)
                out[f"auc_{int(h)}y"] = float(np.atleast_1d(a)[0])
            except Exception:
                pass
    return out, risk, times


def _km(t, e):
    d = pd.DataFrame({"t": t, "e": np.asarray(e).astype(int)}).sort_values("t")
    S, xs, ss = 1.0, [], []
    for tt in sorted(d["t"].unique()):
        n = int((d["t"] >= tt).sum()); dd = int(((d["t"] == tt) & (d["e"] == 1)).sum())
        if n > 0:
            S *= (1 - dd / n)
        xs.append(tt); ss.append(S)
    return np.array(xs), np.array(ss)


# ── Per-config run, feature importance, comparison table ────────────────────────
def _plot_importance(importances, tag, label):
    imp = pd.concat(importances, axis=1)
    s = pd.DataFrame({"mean": imp.mean(axis=1), "std": imp.std(axis=1)}).sort_values("mean")
    s.to_csv(OUT_DIR / f"rsf_feature_importance_{tag}.csv")
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(s.index, s["mean"], xerr=s["std"], color="#2c7bb6", ecolor="#888", capsize=3)
    ax.axvline(0, color="grey", lw=0.8)
    ax.set_xlabel("Permutation importance (drop in validation C-index)", fontweight="bold")
    ax.set_title(f"RSF feature importance — config {tag} ({label})\n"
                 "CN+MCI → AD (mean ± std over 3 seeds, validation)", fontsize=11, fontweight="bold")
    ax.grid(axis="x", alpha=0.3); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"rsf_feature_importance_{tag}.{ext}", bbox_inches="tight", dpi=300)
    plt.close()


def run_config(tag, label, kw):
    print(f"\n=== Config {tag}: {label} ===   {kw}")
    rows = []; importances = []; cal = {"val": {}, "test": {}, "valtest": {}}
    for seed in SEEDS:
        d = load_seed(seed)
        Xtr, ytr, etr, ttr = d["train"]
        rsf = RandomSurvivalForest(**kw).fit(Xtr, ytr)
        d["valtest"] = (pd.concat([d["val"][0], d["test"][0]]), None,
                        np.concatenate([d["val"][2], d["test"][2]]),
                        np.concatenate([d["val"][3], d["test"][3]]))
        d["valtest"] = (d["valtest"][0], Surv.from_arrays(event=d["valtest"][2], time=d["valtest"][3]),
                        d["valtest"][2], d["valtest"][3])
        for sp in ("val", "test", "valtest"):
            X, y, e, t = d[sp]
            m, risk, times = eval_split(rsf, ytr, ttr, etr, X, y, e, t)
            rows.append({"config": label, "tag": tag, "seed": seed, "split": sp, **m})
            print(f"   seed{seed} {sp:8s}: n={m['n']} ev={m['events']}  C={m['c_index']:.3f} "
                  f"C_ipcw={m['c_index_ipcw']:.3f} IBS={m['ibs']:.3f} AUC3y={m['auc_3y']:.3f} AUC5y={m['auc_5y']:.3f}")
            cal[sp][seed] = (*_km(t, e), times, _surv_at(rsf, X, times).mean(axis=0))
        Xv, yv, _, _ = d["val"]
        pi = permutation_importance(rsf, Xv, yv, n_repeats=15, random_state=42, n_jobs=1)
        importances.append(pd.Series(pi.importances_mean, index=FEATURES))
    for sp, name, f in (("val", "validation", f"calibration_val_{tag}"),
                        ("test", "test", f"calibration_test_{tag}"),
                        ("valtest", "val+test (combined)", f"calibration_valtest_{tag}")):
        _plot_calibration(cal[sp], f"{name} — config {tag}", f)
    _plot_importance(importances, tag, label)
    return rows


def render_hp_table(met):
    metric_cols = ["c_index", "c_index_ipcw", "ibs", "auc_3y", "auc_5y"]
    g = met.groupby(["tag", "config", "split"])[metric_cols].agg(["mean", "std"])
    order = {"val": 0, "test": 1, "valtest": 2}
    headers = ["Config", "Held-out", "C-index", "IPCW C", "IBS", "AUC@3y", "AUC@5y"]
    body = []
    for (tag, config) in sorted({(i[0], i[1]) for i in g.index}):
        for split in sorted({i[2] for i in g.index if i[0] == tag}, key=lambda s: order[s]):
            r = g.loc[(tag, config, split)]
            cell = lambda mc: f"{r[(mc, 'mean')]:.3f}±{r[(mc, 'std')]:.3f}"
            body.append([f"{tag}: {config.split(':', 1)[1].strip()}", split,
                         cell("c_index"), cell("c_index_ipcw"), cell("ibs"), cell("auc_3y"), cell("auc_5y")])
    fig, ax = plt.subplots(figsize=(13, 0.5 * len(body) + 1.6)); ax.axis("off")
    ax.set_title("Random Survival Forest — hyperparameter comparison  (CN+MCI → AD; mean ± std over 3 seeds)\n"
                 "A = previous (500 trees, min_leaf 15, m=√);  B = paper (1500 trees, min_leaf 5, m=3, log-rank split)\n"
                 "[scikit-survival has no 'number of split points' (nsplit) param → all splits evaluated]",
                 fontsize=10, fontweight="bold", loc="left", pad=10)
    tbl = ax.table(cellText=body, colLabels=headers, loc="center", cellLoc="center",
                   colWidths=[0.24, 0.12, 0.13, 0.13, 0.12, 0.12, 0.12])
    tbl.auto_set_font_size(False); tbl.set_fontsize(8.5); tbl.scale(1, 1.45)
    for j in range(len(headers)):
        cc = tbl[(0, j)]; cc.set_facecolor("#222"); cc.set_text_props(color="white", weight="bold")
    for i, b in enumerate(body, start=1):
        if b[1] == "valtest":                       # highlight the (least-noisy) combined row
            for j in range(len(headers)):
                tbl[(i, j)].set_facecolor("#eef4fb")
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"rsf_hp_comparison_table.{ext}", bbox_inches="tight", dpi=200)
    plt.close()


# ── Main ────────────────────────────────────────────────────────────────────────
def main():
    configs = [
        ("A", "prev: 500 trees, min_leaf 15, m=sqrt",
         dict(n_estimators=500, min_samples_leaf=15, max_features="sqrt", n_jobs=-1, random_state=42)),
        ("B", "paper: 1500 trees, min_leaf 5, m=3",
         dict(n_estimators=1500, min_samples_leaf=5, max_features=3, n_jobs=-1, random_state=42)),
    ]
    pd.DataFrame({tag: kw for tag, _l, kw in configs}).to_csv(OUT_DIR / "rsf_hyperparameters.csv")
    rows = []
    for tag, label, kw in configs:
        rows += run_config(tag, label, kw)
    met = pd.DataFrame(rows)
    met.to_csv(OUT_DIR / "rsf_metrics.csv", index=False)
    metric_cols = ["c_index", "c_index_ipcw", "ibs", "auc_3y", "auc_5y"]
    summ = met.groupby(["tag", "split"])[metric_cols].agg(["mean", "std"]).round(4)
    summ.to_csv(OUT_DIR / "rsf_metrics_summary.csv")
    print("\nSummary (mean ± std over seeds):"); print(summ.to_string())
    render_hp_table(met)
    print(f"\nDone. Outputs → {OUT_DIR}")


def _plot_calibration(cal_split, name, fname):
    fig, ax = plt.subplots(figsize=(8, 5.5))
    colors = ["#1a9641", "#2c7bb6", "#d7191c"]
    for seed, col in zip(SEEDS, colors):
        xs, ss, grid, S = cal_split[seed]
        ax.step(xs, ss, where="post", color=col, lw=1.8, alpha=0.85, label=f"seed {seed}: KM")
        ax.plot(grid, S, color=col, lw=1.8, ls="--", label=f"seed {seed}: RSF mean Ŝ")
    ax.set_title(f"RSF calibration on {name} — predicted vs observed survival\nCN+MCI → AD",
                 fontsize=12, fontweight="bold")
    ax.set_xlabel("Time (years)", fontweight="bold"); ax.set_ylabel("Survival (AD-free)", fontweight="bold")
    ax.set_ylim(0, 1.02); ax.grid(alpha=0.3); ax.legend(fontsize=8, ncol=2)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"{fname}.{ext}", bbox_inches="tight", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
