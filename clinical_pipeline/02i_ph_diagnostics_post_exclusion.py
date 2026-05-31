"""
02i_ph_diagnostics_post_exclusion.py
====================================
Proportional-hazards (PH) diagnostics for the post-exclusion no-CDR cohort, on the
pooled **CN + MCI -> AD** endpoint (the same one drawn by 02g as `km_cn_mci_to_ad`):
baseline CN and MCI pooled, event = conversion to AD (pCN_to_AD or pMCI),
CN->MCI-only converters censored at FU.

Runs on either timescale (mirrors 02g's time_to_event vs age_to_event):
  --timescale time  (default)  time-from-baseline; outputs -> outputs/ph_diagnostics
  --timescale age              chronological age, LEFT-TRUNCATED (delayed entry: each
                               subject enters the risk set at baseline age, exits at
                               event/censor age); outputs -> outputs/ph_diagnostics_age_to_event

Two complementary PH checks:

1. CATEGORICAL covariates -> complementary log-log plots + a FORMAL PH test
   (the log-log is only a visual; crossings in small/tail subgroups are noisy).
   sex, APOE4 carrier, APOE4 3-way genotype, APOE2 carrier. log(-log S) vs log(time):
   time axis S = ordinary KM; age axis S = left-truncated KM at a floor age. Parallel ⇒
   PH. The formal test (categorical_ph_test.csv) is scaled-Schoenfeld on the time axis
   and an indicator × log(age) interaction (left-truncated Cox) on the age axis; the
   latter targets a smooth monotone trend in log(age), so it has limited power against
   pure crossings — read it alongside subgroup sizes.

2. CONTINUOUS covariates -> PH test.
   TIME axis: scaled-Schoenfeld residual test (lifelines proportional_hazard_test) +
   residual-vs-time plots. AGE axis: lifelines does NOT implement Schoenfeld residuals
   under delayed entry, so we use the equivalent covariate × log(age) interaction in a
   left-truncated counting-process Cox (CoxTimeVaryingFitter) — a significant Wald p on
   the interaction means the effect is non-proportional (Therneau & Grambsch §6).
   a. OBSERVED covariates (full cohort): TIME axis = baseline age + MMSE/ADAS13/MoCA/FAQ;
      AGE axis = cognitive only (baseline age IS the timescale, fit with delayed entry).
      Median-impute + z-score (PH test is scale-invariant; z-scoring aids convergence /
      HR-per-SD reading and matches the tabular baseline in 02h).
   b. meta-PRS-EN-combined (PER SEED, train fold): the meta-PRS is a LEARNED score fit
      per-seed on each train fold (snp_pipeline/49 -> 46's EN Cox); a full-cohort fit
      would leak each subject's own outcome. We read each seed's train-fold scores,
      restrict to that fold, and run the PH test within the fold, per seed.

3. (TIME axis only) age × log(t) interaction — whether modelling a time-varying age
   coefficient explains age's PH violation. N/A on the age axis (age is the timescale).

Inputs
  D:\\ADNI_BIDS_project\\derivatives\\clinical\\no_cdr_stratified_post_exclusion
      \\tabular\\longitudinal\\master_clinical_tabular.csv          (covariates + survival)
  D:\\ADNI_SNP_Omni2.5M_20140220\\source_prs\\meta_prs_en_combined_seed{0,1,2}_train.tsv
      (per-seed train-fold meta-PRS; produced by snp_pipeline/49)

Env: clinical (+ lifelines).  Run:
  python clinical_pipeline/02i_ph_diagnostics_post_exclusion.py [--timescale time|age]
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path

from lifelines import CoxPHFitter, CoxTimeVaryingFitter
from lifelines.statistics import proportional_hazard_test
from lifelines.utils import to_episodic_format

matplotlib.rcParams["font.family"] = "DejaVu Serif"

# ── Paths ──────────────────────────────────────────────────────────────────────
REPO_ROOT  = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
MASTER_CSV = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion"
                  r"\tabular\longitudinal\master_clinical_tabular.csv")
PRS_DIR    = Path(r"D:\ADNI_SNP_Omni2.5M_20140220\source_prs")
OUT_BASE   = REPO_ROOT / "clinical_pipeline" / "outputs"
OUT_DIR    = None        # set in main() from --timescale

EPS = 0.01
AGE_LO = 65              # floor age for the left-truncated log-log (risk sets too small below)
SEEDS = (0, 1, 2)
COG_COVARS = ["MMSE_Total", "ADAS_Total13", "MoCA_Total", "FAQ_Total"]
DOSAGE_TO_3WAY = {0: "non-carrier", 1: "heterozygous", 2: "homozygous"}


# ── Load baseline (one row/subject) + build pooled CN+MCI -> AD survival ─────────
def load_cohort():
    cols = ["Patient_ID", "VISCODE_long", "Date", "YOB", "Sex", "APOE_Alleles",
            "APOE4_Dosage", "bl_dx", "conversion_group", "FU_years", "years_to_AD",
            *COG_COVARS]
    df = pd.read_csv(MASTER_CSV, usecols=cols, low_memory=False)
    df["Date"] = pd.to_datetime(df["Date"], errors="coerce")
    # Baseline = the 'bl' visit; keep the earliest-dated row if a subject has dupes.
    bl = (df[df["VISCODE_long"] == "bl"].dropna(subset=["Date"])
          .sort_values("Date").drop_duplicates("Patient_ID").copy())
    # Baseline age from YOB (birth = 1 July of YOB), as in 02g.
    birth = pd.to_datetime(pd.DataFrame(
        {"year": bl["YOB"].astype(int), "month": 7, "day": 1}).set_index(bl.index))
    bl["age_at_baseline"] = (bl["Date"] - birth) / np.timedelta64(1, "D") / 365.25

    # Pooled CN+MCI -> AD endpoint (mirrors 02g km_cn_mci_to_ad).
    cohort = bl[bl["bl_dx"].isin(["CN", "MCI"])].copy()
    cohort["event"] = cohort["conversion_group"].isin(["pCN_to_AD", "pMCI"]).astype(int)
    t = np.where(cohort["event"] == 1, cohort["years_to_AD"], cohort["FU_years"])
    cohort["t"] = pd.to_numeric(pd.Series(t, index=cohort.index),
                                errors="coerce").fillna(EPS).clip(lower=EPS)
    # Age timescale (left truncation): enter at baseline age, exit at event/censor age.
    cohort["entry_age"] = cohort["age_at_baseline"]
    cohort["exit_age"] = cohort["entry_age"] + cohort["t"]

    # Categorical strata.
    cohort["sex"] = cohort["Sex"]
    cohort["apoe4"] = np.where(cohort["APOE4_Dosage"] >= 1, "carrier", "non-carrier")
    cohort["apoe4_3way"] = cohort["APOE4_Dosage"].map(DOSAGE_TO_3WAY)
    cohort["apoe2"] = np.where(
        cohort["APOE_Alleles"].fillna("").apply(lambda s: "2" in str(s).split("/")),
        "carrier", "non-carrier")
    return cohort


# ── Hand-rolled KM for the log-log plots (matches 02g) ──────────────────────────
def kaplan_meier(times, events):
    """Ordinary KM S(t)."""
    d = pd.DataFrame({"time": times, "event": events}).sort_values("time")
    S, ts, ss = 1.0, [], []
    for t in sorted(d["time"].unique()):
        n = int((d["time"] >= t).sum())
        e = int(((d["time"] == t) & (d["event"] == 1)).sum())
        if n > 0 and e > 0:
            S *= (1 - e / n)
        ts.append(t); ss.append(S)
    return np.array(ts), np.array(ss)


def kaplan_meier_lt(entry, exit_, event, age_start):
    """Left-truncated KM on the age axis, conditioned at `age_start` (as in 02g)."""
    entry = np.asarray(entry, float); exit_ = np.asarray(exit_, float)
    event = np.asarray(event, int)
    event_ages = sorted(t for t in set(exit_[event == 1]) if t >= age_start)
    S, ages, ss = 1.0, [age_start], [1.0]
    for t in event_ages:
        n = int(np.sum((entry <= t) & (exit_ >= t)))
        d = int(np.sum((exit_ == t) & (event == 1)))
        if n > 0 and d > 0:
            S *= (1 - d / n)
        ages.append(t); ss.append(S)
    return np.array(ages), np.array(ss)


CATEGORICAL = {
    "sex":        [("Female", "Female", "#d7191c"), ("Male", "Male", "#2c7bb6")],
    "apoe4":      [("non-carrier", "Non-carrier", "#2c7bb6"), ("carrier", "Carrier", "#d7191c")],
    "apoe4_3way": [("non-carrier", "Non-carrier (0 e4)", "#2c7bb6"),
                   ("heterozygous", "Heterozygous (1 e4)", "#fdae61"),
                   ("homozygous", "Homozygous (2 e4)", "#d7191c")],
    "apoe2":      [("non-carrier", "Non-carrier", "#2c7bb6"), ("carrier", "Carrier (e2+)", "#1a9641")],
}
CAT_TITLE = {"sex": "Sex", "apoe4": "APOE4 carrier status",
             "apoe4_3way": "APOE4 genotype", "apoe2": "APOE2 carrier status"}


def loglog_plot(cohort, col, strata, age):
    """Complementary log-log: log(-log S) vs log(time). Parallel ⇒ PH holds."""
    fig, ax = plt.subplots(figsize=(8, 6))
    for grp, disp, color in strata:
        g = cohort[cohort[col] == grp]
        if len(g) < 5:
            continue
        if age:
            xs, ss = kaplan_meier_lt(g["entry_age"].values, g["exit_age"].values,
                                     g["event"].values, AGE_LO)
        else:
            xs, ss = kaplan_meier(g["t"].values, g["event"].values)
        m = (xs > 0) & (ss > 0) & (ss < 1)          # log domains
        if m.sum() < 2:
            continue
        x = np.log(xs[m]); y = np.log(-np.log(ss[m]))
        ax.step(x, y, where="post", color=color, linewidth=2.0,
                label=f"{disp}  (N={len(g)}, events={int(g['event'].sum())})")
    axis_lbl = "log(age)  [log years]" if age else "log(time)  [log years]"
    sub = ("CN+MCI → AD, age timescale (left-truncated)" if age
           else "CN+MCI → AD") + " (parallel ⇒ proportional hazards)"
    ax.set_xlabel(axis_lbl, fontsize=11, fontweight="bold")
    ax.set_ylabel(r"$\log(-\log \hat{S})$", fontsize=11, fontweight="bold")
    ax.set_title(f"Complementary log-log — {CAT_TITLE[col]}\n{sub}",
                 fontsize=12, fontweight="bold", pad=10)
    ax.grid(alpha=0.3, linewidth=0.5)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.legend(fontsize=9, frameon=True, framealpha=0.9, edgecolor="#cccccc")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"loglog_{col}.{ext}", bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved loglog_{col}")


# ── Schoenfeld residual test + plots ────────────────────────────────────────────
def zscore(frame, cols):
    out = frame.copy()
    for c in cols:
        x = pd.to_numeric(out[c], errors="coerce")
        x = x.fillna(x.median())
        sd = x.std(ddof=0) or 1.0
        out[c] = (x - x.mean()) / sd
    return out


def fit_cox(df_fit, age):
    """Fit Cox on the chosen timescale. Age axis uses delayed entry (entry_col)."""
    cph = CoxPHFitter(penalizer=1e-3)
    if age:
        cph.fit(df_fit, duration_col="exit_age", event_col="e", entry_col="entry_age")
    else:
        cph.fit(df_fit, duration_col="t", event_col="e")
    return cph


def schoenfeld_resid_plot(cph, df_fit, cov, fname, title, dur_col, xlabel):
    """Scatter scaled Schoenfeld residuals vs event time/age + rolling-mean trend."""
    res = cph.compute_residuals(df_fit, kind="scaled_schoenfeld")
    if cov not in res.columns or len(res) < 3:
        return
    times = df_fit.loc[res.index, dur_col].values
    order = np.argsort(times)
    tt, rr = times[order], res[cov].values[order]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.axhline(0.0, color="grey", linewidth=0.8, linestyle="--", alpha=0.7)
    ax.scatter(tt, rr, s=18, alpha=0.5, color="#2c7bb6", edgecolor="none")
    if len(tt) >= 8:                                   # rolling-mean trend (no statsmodels dep)
        w = max(5, len(tt) // 6)
        trend = pd.Series(rr).rolling(w, center=True, min_periods=max(2, w // 2)).mean()
        ax.plot(tt, trend.values, color="#d7191c", linewidth=2.0, label=f"rolling mean (w={w})")
        ax.legend(fontsize=9)
    ax.set_xlabel(xlabel, fontsize=11, fontweight="bold")
    ax.set_ylabel(f"Scaled Schoenfeld residual\n({cov})", fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=10)
    ax.grid(alpha=0.3, linewidth=0.5)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"{fname}.{ext}", bbox_inches="tight", dpi=300)
    plt.close()


# ── Age axis: PH test via time-varying coefficient (lifelines has no Schoenfeld
#    residuals under delayed entry; the covariate × log(age) interaction in a
#    counting-process Cox is an equivalent, left-truncation-aware PH test). ───────
def build_age_episodes(df, covars):
    """Split each subject into (start, stop] age intervals from entry_age to exit_age
    (delayed entry → left truncation is native to the counting-process format)."""
    recs = []
    for row in df.itertuples(index=False):
        e0, e1 = float(row.entry_age), float(row.exit_age)
        if e1 <= e0:
            e1 = e0 + EPS
        knots = sorted(set([e0] + [k for k in range(int(np.floor(e0)) + 1,
                                                     int(np.ceil(e1)) + 1) if e0 < k < e1] + [e1]))
        base = {c: getattr(row, c) for c in covars}
        for a, b in zip(knots[:-1], knots[1:]):
            if b <= a:
                continue
            rec = {"id": row.id, "start": a, "stop": b,
                   "e": int(row.e) if abs(b - e1) < 1e-9 else 0}
            rec.update(base)
            recs.append(rec)
    return pd.DataFrame(recs)


def tv_loghr_plot(b_const, b_int, cov2, fname, title, age_lo, age_hi):
    """Time-varying log-HR  β + β_int·log(age)  with a 95% CI band.

    cov2 = 2×2 covariance of (β, β_int); Var[g(age)] = V00 + (log age)²·V11
    + 2·log age·V01. The band makes the PH question honest: if a flat line fits
    inside it, the apparent slope is not distinguishable from no time trend."""
    aa = np.linspace(age_lo, age_hi, 200)
    L = np.log(aa)
    est = b_const + b_int * L
    var = cov2[0, 0] + L ** 2 * cov2[1, 1] + 2 * L * cov2[0, 1]
    se = np.sqrt(np.clip(var, 0.0, None))
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.axhline(0.0, color="grey", linewidth=0.9, linestyle="--", alpha=0.7, label="no effect")
    ax.fill_between(aa, est - 1.96 * se, est + 1.96 * se, color="#d7191c",
                    alpha=0.15, label="95% CI")
    ax.plot(aa, est, color="#d7191c", linewidth=2.2)
    ax.set_xlabel("Age (years)", fontsize=11, fontweight="bold")
    ax.set_ylabel("log-HR per +1 SD  (β + β_int·log age)", fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=10)
    ax.grid(alpha=0.3, linewidth=0.5)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.legend(fontsize=9)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"{fname}.{ext}", bbox_inches="tight", dpi=300)
    plt.close()


def age_ph_test(df_fit, covars, age_hi, plot=True):
    """Left-truncated counting-process Cox with covariate × log(age) interactions.
    Returns a per-covariate DataFrame (coef, coef_x_logage, PH Wald p, HR/SD)."""
    df = df_fit.copy()
    df["id"] = np.arange(len(df))
    long = build_age_episodes(df, covars)
    for c in covars:
        long[c + "__xlog"] = long[c] * np.log(long["stop"].clip(lower=1e-3))
    keep = ["id", "start", "stop", "e"] + covars + [c + "__xlog" for c in covars]
    ctv = CoxTimeVaryingFitter(penalizer=1e-3)
    ctv.fit(long[keep], id_col="id", event_col="e",
            start_col="start", stop_col="stop", show_progress=False)
    sm = ctv.summary
    try:
        V = ctv.variance_matrix_
    except Exception:
        V = None
    rows = []
    age_lo = max(AGE_LO, float(df["entry_age"].min()))
    for c in covars:
        b0 = float(sm.loc[c, "coef"]); bi = float(sm.loc[c + "__xlog", "coef"])
        pi = float(sm.loc[c + "__xlog", "p"])
        rows.append({"covariate": c, "coef": b0, "coef_x_logage": bi,
                     "p_ph_interaction": pi, "hr_per_sd_at_logage0": float(np.exp(b0))})
        if plot:
            try:
                cov2 = V.loc[[c, c + "__xlog"], [c, c + "__xlog"]].values
            except Exception:
                se = ctv.standard_errors_
                cov2 = np.diag([float(se[c]) ** 2, float(se[c + "__xlog"]) ** 2])
            tv_loghr_plot(b0, bi, cov2, f"phtv_{c}",
                          f"Time-varying {c} effect — CN+MCI → AD (age axis)\n"
                          f"β_int (·log age) = {bi:+.3f}, Wald p = {pi:.3f}  "
                          f"({'non-PH' if pi < 0.05 else 'PH ok'})",
                          age_lo, age_hi)
    return pd.DataFrame(rows)


def categorical_ph_test(cohort, age):
    """Formal PH test for the categorical covariates (the log-log plots are only a
    visual check). Time axis: scaled-Schoenfeld. Age axis: indicator × log(age)
    interaction in a left-truncated Cox. p < 0.05 ⇒ non-proportional / crossing."""
    enc = {
        "sex_Male":      (cohort["sex"] == "Male").astype(float),
        "apoe4_carrier": (cohort["apoe4"] == "carrier").astype(float),
        "apoe2_carrier": (cohort["apoe2"] == "carrier").astype(float),
        "apoe4_dosage":  pd.to_numeric(cohort["APOE4_Dosage"], errors="coerce"),
    }
    rows = []
    for name, x in enc.items():
        try:
            if not age:
                d = pd.DataFrame({"t": cohort["t"].values, "e": cohort["event"].values,
                                  name: x.values}).dropna()
                cph = CoxPHFitter(penalizer=1e-3); cph.fit(d, duration_col="t", event_col="e")
                s = proportional_hazard_test(cph, d, time_transform="rank").summary.loc[name]
                rows.append({"covariate": name, "test": "schoenfeld(rank)",
                             "statistic": float(s["test_statistic"]), "p": float(s["p"])})
            else:
                d = pd.DataFrame({"entry_age": cohort["entry_age"].values,
                                  "exit_age": cohort["exit_age"].values,
                                  "e": cohort["event"].values, name: x.values}).dropna()
                o = age_ph_test(d, [name], float(d["exit_age"].quantile(0.98)),
                                plot=False).iloc[0]
                rows.append({"covariate": name, "test": "X×log(age) Wald",
                             "statistic": float(o["coef_x_logage"]), "p": float(o["p_ph_interaction"])})
        except Exception as exc:
            rows.append({"covariate": name, "test": "failed", "statistic": float("nan"),
                         "p": float("nan")})
            print(f"  [{name}] PH test failed ({exc}); likely too few events.")
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "categorical_ph_test.csv", index=False)
    print("  Categorical PH test (p<0.05 ⇒ non-proportional / crossing log-log):")
    print(out.to_string(index=False))
    return out


def observed_continuous(cohort, age):
    """Continuous-covariate PH check, full pooled cohort.
    Time axis: scaled-Schoenfeld test (age + cognitive).
    Age axis: covariate × log(age) interaction test (cognitive only — age is the
    timescale; the Cox uses delayed entry at baseline age)."""
    if not age:
        covars = ["age_at_baseline"] + COG_COVARS
        df_fit = zscore(cohort[["t", "event"] + covars].rename(columns={"event": "e"}),
                        covars).dropna(subset=["t", "e"])
        cph = fit_cox(df_fit, age=False)
        summ = proportional_hazard_test(cph, df_fit, time_transform="rank").summary.copy()
        summ.index.name = "covariate"
        summ.to_csv(OUT_DIR / "schoenfeld_observed.csv")
        print(f"\nObserved-covariate Schoenfeld test (time-from-baseline; full cohort, "
              f"n={len(df_fit)}, events={int(df_fit['e'].sum())}):")
        print(summ[["test_statistic", "p"]].to_string())
        for cov in covars:
            schoenfeld_resid_plot(
                cph, df_fit, cov, f"schoenfeld_{cov}",
                f"Scaled Schoenfeld residuals — {cov}\n"
                f"CN+MCI → AD (flat ⇒ PH; p={summ.loc[cov, 'p']:.3f})",
                "t", "Event time (years)")
            print(f"  saved schoenfeld_{cov}")
        return

    covars = COG_COVARS
    df_fit = zscore(cohort[["entry_age", "exit_age", "event"] + covars]
                    .rename(columns={"event": "e"}), covars).dropna()
    age_hi = float(cohort["exit_age"].quantile(0.98))
    out = age_ph_test(df_fit, covars, age_hi)
    out.to_csv(OUT_DIR / "ph_timevarying_observed.csv", index=False)
    print(f"\nObserved-covariate PH test (age timescale, left-truncated; full cohort, "
          f"n={len(df_fit)}, events={int(df_fit['e'].sum())})")
    print("  covariate × log(age) interaction Wald test (p<0.05 ⇒ non-PH):")
    print(out[["covariate", "coef_x_logage", "p_ph_interaction"]].to_string(index=False))
    print("  saved phtv_<cov> plots + ph_timevarying_observed.csv")


def meta_prs_continuous(cohort, age):
    """Per-seed meta-PRS PH check within each train fold (Schoenfeld on time axis;
    meta-PRS × log(age) interaction on the age axis)."""
    surv_cols = ["entry_age", "exit_age"] if age else ["t"]
    surv = cohort.set_index("Patient_ID")[surv_cols + ["event"]]
    rows = []
    for seed in SEEDS:
        tsv = PRS_DIR / f"meta_prs_en_combined_seed{seed}_train.tsv"
        if not tsv.exists():
            print(f"  [seed {seed}] {tsv.name} missing — run snp_pipeline/49 first; skipping.")
            continue
        prs = pd.read_csv(tsv, sep="\t")
        prs["Patient_ID"] = prs["Patient_ID"].astype(str)
        merged = prs.merge(surv, left_on="Patient_ID", right_index=True, how="inner")
        df_fit = merged.rename(columns={"event": "e"})[surv_cols + ["e", "meta_prs_EN_combined"]]
        df_fit = zscore(df_fit, ["meta_prs_EN_combined"]).dropna()
        if len(df_fit) < 20 or df_fit["e"].sum() < 5:
            print(f"  [seed {seed}] too few (n={len(df_fit)}, ev={int(df_fit['e'].sum())}); skipping.")
            continue
        if not age:
            cph = fit_cox(df_fit, age=False)
            s = proportional_hazard_test(cph, df_fit, time_transform="rank").summary.loc["meta_prs_EN_combined"]
            p = float(s["p"])
            rows.append({"seed": seed, "n": len(df_fit), "events": int(df_fit["e"].sum()),
                         "test_statistic": float(s["test_statistic"]), "p": p,
                         "hr_per_sd": float(np.exp(cph.params_["meta_prs_EN_combined"]))})
            schoenfeld_resid_plot(
                cph, df_fit, "meta_prs_EN_combined", f"schoenfeld_meta_prs_seed{seed}",
                f"Scaled Schoenfeld residuals — meta-PRS-EN (seed {seed} train fold)\n"
                f"CN+MCI → AD (flat ⇒ PH; p={p:.3f}, n={len(df_fit)})",
                "t", "Event time (years)")
        else:
            age_hi = float(df_fit["exit_age"].quantile(0.98))
            o = age_ph_test(df_fit, ["meta_prs_EN_combined"], age_hi).iloc[0]
            p = float(o["p_ph_interaction"])
            rows.append({"seed": seed, "n": len(df_fit), "events": int(df_fit["e"].sum()),
                         "coef_x_logage": float(o["coef_x_logage"]), "p_ph_interaction": p,
                         "hr_per_sd_at_logage0": float(o["hr_per_sd_at_logage0"])})
            # age_ph_test already saved phtv_meta_prs_EN_combined; rename per seed.
            for ext in ("png", "pdf"):
                src = OUT_DIR / f"phtv_meta_prs_EN_combined.{ext}"
                if src.exists():
                    src.replace(OUT_DIR / f"phtv_meta_prs_seed{seed}.{ext}")
        print(f"  [seed {seed}] n={len(df_fit)} events={int(df_fit['e'].sum())} PH p={p:.3f}")
    if rows:
        fname = "schoenfeld_meta_prs_per_seed.csv" if not age else "ph_timevarying_meta_prs_per_seed.csv"
        pd.DataFrame(rows).to_csv(OUT_DIR / fname, index=False)
        print(f"  wrote {fname} ({len(rows)} seeds)")


def age_time_interaction(cohort):
    """TIME axis only: does an age × log(t) interaction explain age's PH violation?

    Two equivalent views of the SAME question:
      (A) Score test — the scaled-Schoenfeld test under time_transform g(t) IS the
          Grambsch-Therneau score test for beta2 = 0 in  h(t) = h0(t) exp(b1*age +
          b2*age*g(t)).  g(t)=log(t) corresponds to the age*log t interaction. We
          report age's Schoenfeld test under several transforms so g(t) is explicit.
      (B) Wald test — explicitly FIT the time-varying model with a constant age term
          plus an age*log(t) term (CoxTimeVaryingFitter on episode-split data) and
          report b1, b2 and the Wald p on b2. A significant positive b2 means age's
          log-HR rises with time (the rising Schoenfeld residual we saw).
    Cognitive covariates stay in as fixed adjusters, matching observed_schoenfeld.
    """
    covars = ["age_at_baseline"] + COG_COVARS
    df = cohort[["t", "event"] + covars].rename(columns={"event": "e"})
    df = zscore(df, covars).dropna(subset=["t", "e"]).reset_index(drop=True)

    # (A) Schoenfeld test for age under several time transforms.
    cph = CoxPHFitter(penalizer=1e-3)
    cph.fit(df, duration_col="t", event_col="e")
    rows = []
    for tf in ("rank", "identity", "log", "km"):
        s = proportional_hazard_test(cph, df, time_transform=tf).summary
        rows.append({"time_transform": tf,
                     "age_test_statistic": float(s.loc["age_at_baseline", "test_statistic"]),
                     "age_p": float(s.loc["age_at_baseline", "p"])})
    tdf = pd.DataFrame(rows)
    tdf.to_csv(OUT_DIR / "schoenfeld_age_transforms.csv", index=False)
    print("\n  (A) age Schoenfeld score test by time transform "
          "(time_transform='log' == age·log t interaction):")
    print(tdf.to_string(index=False))

    # (B) Explicit time-varying fit: b1*age + b2*age*log(t).
    df["id"] = np.arange(len(df))
    long = to_episodic_format(df, duration_col="t", event_col="e",
                              id_col="id", time_gaps=0.5)
    long["age_x_logt"] = long["age_at_baseline"] * np.log(long["stop"].clip(lower=1e-3))
    keep = ["id", "start", "stop", "e", "age_at_baseline", "age_x_logt"] + COG_COVARS
    ctv = CoxTimeVaryingFitter(penalizer=1e-3)
    ctv.fit(long[keep], id_col="id", event_col="e",
            start_col="start", stop_col="stop", show_progress=False)
    sm = ctv.summary
    b1 = float(sm.loc["age_at_baseline", "coef"])
    b2 = float(sm.loc["age_x_logt", "coef"])
    p2 = float(sm.loc["age_x_logt", "p"])
    sm[["coef", "exp(coef)", "p"]].to_csv(OUT_DIR / "age_timevarying_coef.csv")
    print(f"\n  (B) time-varying Cox  (β1·age + β2·age·log t, + cognitive adjusters):")
    print(f"      β1 (age, constant)   = {b1:+.4f}")
    print(f"      β2 (age × log t)     = {b2:+.4f}   Wald p = {p2:.2e}")
    print(f"      → age log-HR(t) = {b1:+.3f} {('+' if b2 >= 0 else '-')} "
          f"{abs(b2):.3f}·log t  "
          f"({'rises' if b2 > 0 else 'falls'} with time ⇒ confirms non-PH)")

    # Plot the implied time-varying age log-HR.
    tt = np.linspace(0.5, float(df['t'].quantile(0.98)), 200)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(tt, b1 + b2 * np.log(tt), color="#d7191c", linewidth=2.2)
    ax.axhline(b1, color="grey", linewidth=0.9, linestyle="--", alpha=0.7,
               label=f"constant-age model (β₁={b1:+.2f})")
    ax.set_xlabel("Time (years)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Age log-HR per +1 SD  (β₁ + β₂·log t)", fontsize=11, fontweight="bold")
    ax.set_title("Time-varying age effect — CN+MCI → AD\n"
                 f"β₂ (age × log t) = {b2:+.3f}, Wald p = {p2:.1e}",
                 fontsize=12, fontweight="bold", pad=10)
    ax.grid(alpha=0.3, linewidth=0.5)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.legend(fontsize=9)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"age_timevarying_loghr.{ext}", bbox_inches="tight", dpi=300)
    plt.close()
    print("      saved age_timevarying_loghr + schoenfeld_age_transforms.csv "
          "+ age_timevarying_coef.csv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--timescale", choices=["time", "age"], default="time",
                    help="time = time-from-baseline; age = chronological age (left-truncated)")
    args = ap.parse_args()
    age = args.timescale == "age"

    global OUT_DIR
    OUT_DIR = OUT_BASE / ("ph_diagnostics_age_to_event" if age else "ph_diagnostics")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    cohort = load_cohort()
    print(f"Pooled CN+MCI → AD cohort [{args.timescale} axis]: N={len(cohort)}, "
          f"events={int(cohort['event'].sum())}, median FU={cohort['t'].median():.1f} y")

    n_steps = 3 if age else 4
    print(f"\n[1/{n_steps}] Categorical log-log plots + formal PH test ...")
    for col, strata in CATEGORICAL.items():
        loglog_plot(cohort, col, strata, age)
    categorical_ph_test(cohort, age)

    cog_note = "cognitive only — age is the timescale" if age else "age + cognitive"
    test_note = "covariate × log(age) interaction" if age else "scaled-Schoenfeld"
    print(f"\n[2/{n_steps}] Observed-covariate {test_note} PH test ({cog_note}, full cohort) ...")
    observed_continuous(cohort, age)

    print(f"\n[3/{n_steps}] meta-PRS-EN PH test (per-seed train fold) ...")
    meta_prs_continuous(cohort, age)

    if not age:
        print(f"\n[4/{n_steps}] age × log(t) interaction (does it explain age's PH violation?) ...")
        age_time_interaction(cohort)

    print(f"\nDone. Outputs -> {OUT_DIR}")


if __name__ == "__main__":
    main()
