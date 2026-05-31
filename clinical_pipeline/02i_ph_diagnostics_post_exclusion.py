"""
02i_ph_diagnostics_post_exclusion.py
====================================
Proportional-hazards (PH) diagnostics for the post-exclusion no-CDR cohort, on the
pooled **CN + MCI -> AD** survival endpoint (the same one drawn by 02g as
`km_cn_mci_to_ad`): baseline CN and MCI pooled, event = conversion to AD
(pCN_to_AD or pMCI), CN->MCI-only converters censored at FU.

Two complementary PH checks:

1. CATEGORICAL covariates -> complementary log-log plots (full cohort).
   For sex, APOE4 carrier, APOE4 3-way genotype and APOE2 carrier we plot
   log(-log S(t)) against log(t) per stratum. Under proportional hazards the
   strata curves are roughly parallel (a constant vertical gap = constant log-HR).

2. CONTINUOUS covariates -> scaled Schoenfeld residual test (lifelines
   proportional_hazard_test) + residual-vs-time plots.
   a. OBSERVED covariates (full cohort): baseline age, MMSE, ADAS13, MoCA, FAQ.
      One multivariable Cox; median-impute + z-score (the test is invariant to
      linear scaling — z-scoring only aids convergence / HR-per-SD reading and
      matches the tabular-baseline convention in 02h).
   b. meta-PRS-EN-combined (PER SEED, train fold): the meta-PRS is a LEARNED score
      fit per-seed on each train fold (snp_pipeline/49 -> 46's EN Cox); a
      full-cohort fit would leak each subject's own outcome. So we read each seed's
      train-fold scores, restrict to that fold's subjects, and run the Schoenfeld
      test within the fold. Reported per seed so the meta-PRS PH behaviour is visible
      in each training seed independently.

Inputs
  D:\\ADNI_BIDS_project\\derivatives\\clinical\\no_cdr_stratified_post_exclusion
      \\tabular\\longitudinal\\master_clinical_tabular.csv          (covariates + survival)
  D:\\ADNI_SNP_Omni2.5M_20140220\\source_prs\\meta_prs_en_combined_seed{0,1,2}_train.tsv
      (per-seed train-fold meta-PRS; produced by snp_pipeline/49)

Outputs -> clinical_pipeline/outputs/ph_diagnostics/  (each plot as .png + .pdf)
  loglog_sex / loglog_apoe4 / loglog_apoe4_3way / loglog_apoe2          (categorical)
  schoenfeld_observed.csv  + schoenfeld_<cov>.{png,pdf}                 (age + cognitive)
  schoenfeld_meta_prs_per_seed.csv + schoenfeld_meta_prs_seed{N}.{png,pdf}

Env: clinical (+ lifelines).  Run:
  python clinical_pipeline/02i_ph_diagnostics_post_exclusion.py
"""
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path

from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

matplotlib.rcParams["font.family"] = "DejaVu Serif"

# ── Paths ──────────────────────────────────────────────────────────────────────
REPO_ROOT  = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
MASTER_CSV = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion"
                  r"\tabular\longitudinal\master_clinical_tabular.csv")
PRS_DIR    = Path(r"D:\ADNI_SNP_Omni2.5M_20140220\source_prs")
OUT_DIR    = REPO_ROOT / "clinical_pipeline" / "outputs" / "ph_diagnostics"
OUT_DIR.mkdir(parents=True, exist_ok=True)

EPS = 0.01
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

    # Categorical strata.
    cohort["sex"] = cohort["Sex"]
    cohort["apoe4"] = np.where(cohort["APOE4_Dosage"] >= 1, "carrier", "non-carrier")
    cohort["apoe4_3way"] = cohort["APOE4_Dosage"].map(DOSAGE_TO_3WAY)
    cohort["apoe2"] = np.where(
        cohort["APOE_Alleles"].fillna("").apply(lambda s: "2" in str(s).split("/")),
        "carrier", "non-carrier")
    return cohort


# ── Hand-rolled KM (S(t)) for the log-log plots (matches 02g) ───────────────────
def kaplan_meier(times, events):
    d = pd.DataFrame({"time": times, "event": events}).sort_values("time")
    S, ts, ss = 1.0, [], []
    for t in sorted(d["time"].unique()):
        n = int((d["time"] >= t).sum())
        e = int(((d["time"] == t) & (d["event"] == 1)).sum())
        if n > 0 and e > 0:
            S *= (1 - e / n)
        ts.append(t); ss.append(S)
    return np.array(ts), np.array(ss)


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


def loglog_plot(cohort, col, strata):
    """Complementary log-log: log(-log S(t)) vs log(t). Parallel => PH holds."""
    fig, ax = plt.subplots(figsize=(8, 6))
    for grp, disp, color in strata:
        g = cohort[cohort[col] == grp]
        if len(g) < 5:
            continue
        ts, ss = kaplan_meier(g["t"].values, g["event"].values)
        m = (ts > 0) & (ss > 0) & (ss < 1)          # log domains
        if m.sum() < 2:
            continue
        x = np.log(ts[m]); y = np.log(-np.log(ss[m]))
        ax.step(x, y, where="post", color=color, linewidth=2.0,
                label=f"{disp}  (N={len(g)}, events={int(g['event'].sum())})")
    ax.set_xlabel("log(time)  [log years]", fontsize=11, fontweight="bold")
    ax.set_ylabel(r"$\log(-\log \hat{S}(t))$", fontsize=11, fontweight="bold")
    ax.set_title(f"Complementary log-log — {CAT_TITLE[col]}\n"
                 f"CN+MCI → AD (parallel curves ⇒ proportional hazards)",
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


def schoenfeld_resid_plot(cph, df_fit, cov, fname, title):
    """Scatter scaled Schoenfeld residuals vs event time + rolling-mean trend."""
    res = cph.compute_residuals(df_fit, kind="scaled_schoenfeld")
    if cov not in res.columns or len(res) < 3:
        return
    times = df_fit.loc[res.index, "t"].values
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
    ax.set_xlabel("Event time (years)", fontsize=11, fontweight="bold")
    ax.set_ylabel(f"Scaled Schoenfeld residual\n({cov})", fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=10)
    ax.grid(alpha=0.3, linewidth=0.5)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"{fname}.{ext}", bbox_inches="tight", dpi=300)
    plt.close()


def observed_schoenfeld(cohort):
    """Multivariable Cox on age + cognitive covariates, full pooled cohort."""
    covars = ["age_at_baseline"] + COG_COVARS
    df_fit = cohort[["t", "event"] + covars].rename(columns={"event": "e"})
    df_fit = zscore(df_fit, covars)
    df_fit = df_fit.dropna(subset=["t", "e"])
    cph = CoxPHFitter(penalizer=1e-3)
    cph.fit(df_fit, duration_col="t", event_col="e")
    test = proportional_hazard_test(cph, df_fit, time_transform="rank")
    summ = test.summary.copy()
    summ.index.name = "covariate"
    summ.to_csv(OUT_DIR / "schoenfeld_observed.csv")
    print("\nObserved-covariate Schoenfeld test (full cohort, "
          f"n={len(df_fit)}, events={int(df_fit['e'].sum())}):")
    print(summ[["test_statistic", "p"]].to_string())
    for cov in covars:
        schoenfeld_resid_plot(
            cph, df_fit, cov, f"schoenfeld_{cov}",
            f"Scaled Schoenfeld residuals — {cov}\n"
            f"CN+MCI → AD (flat ⇒ PH; p={summ.loc[cov, 'p']:.3f})")
        print(f"  saved schoenfeld_{cov}")


def meta_prs_schoenfeld(cohort):
    """Per-seed meta-PRS Schoenfeld test on the pooled endpoint, within each train fold."""
    surv = cohort.set_index("Patient_ID")[["t", "event"]]
    rows = []
    for seed in SEEDS:
        tsv = PRS_DIR / f"meta_prs_en_combined_seed{seed}_train.tsv"
        if not tsv.exists():
            print(f"  [seed {seed}] {tsv.name} missing — run snp_pipeline/49 first; skipping.")
            continue
        prs = pd.read_csv(tsv, sep="\t")
        prs["Patient_ID"] = prs["Patient_ID"].astype(str)
        merged = prs.merge(surv, left_on="Patient_ID", right_index=True, how="inner")
        df_fit = merged.rename(columns={"event": "e"})[["t", "e", "meta_prs_EN_combined"]]
        df_fit = zscore(df_fit, ["meta_prs_EN_combined"]).dropna()
        if len(df_fit) < 20 or df_fit["e"].sum() < 5:
            print(f"  [seed {seed}] too few (n={len(df_fit)}, ev={int(df_fit['e'].sum())}); skipping.")
            continue
        cph = CoxPHFitter(penalizer=1e-3)
        cph.fit(df_fit, duration_col="t", event_col="e")
        test = proportional_hazard_test(cph, df_fit, time_transform="rank")
        s = test.summary.loc["meta_prs_EN_combined"]
        rows.append({"seed": seed, "n": len(df_fit), "events": int(df_fit["e"].sum()),
                     "test_statistic": float(s["test_statistic"]), "p": float(s["p"]),
                     "hr_per_sd": float(np.exp(cph.params_["meta_prs_EN_combined"]))})
        schoenfeld_resid_plot(
            cph, df_fit, "meta_prs_EN_combined", f"schoenfeld_meta_prs_seed{seed}",
            f"Scaled Schoenfeld residuals — meta-PRS-EN (seed {seed} train fold)\n"
            f"CN+MCI → AD (flat ⇒ PH; p={float(s['p']):.3f}, n={len(df_fit)})")
        print(f"  [seed {seed}] n={len(df_fit)} events={int(df_fit['e'].sum())} "
              f"PH p={float(s['p']):.3f} -> saved schoenfeld_meta_prs_seed{seed}")
    if rows:
        pd.DataFrame(rows).to_csv(OUT_DIR / "schoenfeld_meta_prs_per_seed.csv", index=False)
        print(f"  wrote schoenfeld_meta_prs_per_seed.csv ({len(rows)} seeds)")


def main():
    cohort = load_cohort()
    print(f"Pooled CN+MCI → AD cohort: N={len(cohort)}, "
          f"events={int(cohort['event'].sum())}, "
          f"median FU={cohort['t'].median():.1f} y")

    print("\n[1/3] Categorical log-log plots ...")
    for col, strata in CATEGORICAL.items():
        loglog_plot(cohort, col, strata)

    print("\n[2/3] Observed-covariate Schoenfeld (age + cognitive, full cohort) ...")
    observed_schoenfeld(cohort)

    print("\n[3/3] meta-PRS-EN Schoenfeld (per-seed train fold) ...")
    meta_prs_schoenfeld(cohort)

    print(f"\nDone. Outputs -> {OUT_DIR}")


if __name__ == "__main__":
    main()
