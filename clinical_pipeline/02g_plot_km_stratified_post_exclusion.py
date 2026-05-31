"""
02g_plot_km_stratified_post_exclusion.py
========================================
Kaplan-Meier survival curves for the *post-exclusion* no-CDR stratified cohort.

Five analyses (classic survival S(t), falling from 1.0):
  1. sCN vs pCN     — baseline CN, event = CN -> AD (pCN_to_AD),            AD-free survival
  2. sMCI vs pMCI   — baseline MCI, event = MCI -> AD (pMCI),              AD-free survival
  3. CN -> MCI      — baseline CN, event = CN -> MCI (pCN_to_MCI, stop at MCI),
                      MCI-free survival
  4. CN -> MCI/AD   — baseline CN, event = ANY progression (pCN_to_MCI OR pCN_to_AD),
                      progression-free survival; event time = first progression
                      (years_to_MCI, falling back to years_to_AD)
  5. CN+MCI -> AD   — POOLED baseline CN AND MCI, event = conversion to AD
                      (pCN_to_AD OR pMCI), AD-free survival. CN -> MCI-only converters
                      (pCN_to_MCI) are NOT events here — they remain at risk and are
                      censored at FU (i.e. "pcn_to_ad + pmci_to_ad, excluding pcn_to_mci").

CN -> MCI (stop at MCI) and CN -> AD (progress to AD) are mutually exclusive endpoints:
in the CN -> MCI analysis sCN and pCN_to_AD are right-censored; in CN -> AD, sCN and
pCN_to_MCI are censored. The CN -> MCI/AD analysis instead JOINS the two progressor
groups, so no CN progressor is censored as a non-event.

Each analysis is produced three ways:
  - pooled (one curve)
  - stratified by APOE4 carrier status      (carrier = >=1 e4 allele vs non-carrier)
  - stratified by APOE4 genotype, 3-way     (non-carrier / heterozygous / homozygous)
Stratified figures overlay the strata with per-stratum number-at-risk rows and a
k-sample log-rank test. APOE4 is fully observed in this cohort (no missing).

The "survival" curve is the proportion remaining event-free over time. Stable subjects
are right-censored at their total follow-up. AD-at-baseline (AD_bl) is prevalent AD and
excluded. Event membership follows `conversion_group` so the sustained-conversion logic
baked into that column is respected.

Data source:
  D:\\ADNI_BIDS_project\\derivatives\\clinical\\no_cdr_stratified_post_exclusion
      \\tabular\\longitudinal\\master_clinical_tabular.csv

Outputs (each as .png + .pdf):
  outputs/time_to_event/  — time-from-baseline x-axis:
    km_pcn_to_ad[ _by_apoe4 / _by_apoe4_3way ]        (CN -> AD;     events 25)
    km_pmci_to_ad[ _by_apoe4 / _by_apoe4_3way ]       (MCI -> AD;    events 121)
    km_pcn_to_mci[ _by_apoe4 / _by_apoe4_3way ]       (CN -> MCI;    events 35)
    km_pcn_to_mci_or_ad[ _by_apoe4 / _by_apoe4_3way ] (CN -> MCI/AD; events 60)
    km_cn_mci_to_ad[ _by_apoe4 / _by_apoe4_3way ]     (CN+MCI -> AD pooled; events 146)
  outputs/age_to_event/   — age x-axis, left-truncated (delayed entry):
    same five analyses x three views, each with an "_age" suffix
    (e.g. km_pcn_to_ad_age, km_pmci_to_ad_by_apoe4_3way_age, ...)

Run:
  python clinical_pipeline/02g_plot_km_stratified_post_exclusion.py
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import chi2

matplotlib.rcParams["font.family"] = "DejaVu Serif"

# ── Paths ──────────────────────────────────────────────────────────────────────
REPO_ROOT  = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
MASTER_CSV = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion"
                  r"\tabular\longitudinal\master_clinical_tabular.csv")
OUT_DIR     = REPO_ROOT / "clinical_pipeline" / "outputs" / "time_to_event"
AGE_OUT_DIR = REPO_ROOT / "clinical_pipeline" / "outputs" / "age_to_event"
OUT_DIR.mkdir(parents=True, exist_ok=True)
AGE_OUT_DIR.mkdir(parents=True, exist_ok=True)

EPS = 0.01  # floor for non-positive event/censor times (years)


# ── Load + dedupe to one row per subject ───────────────────────────────────────
print("Loading master CSV...")
clin = pd.read_csv(
    MASTER_CSV,
    usecols=["Patient_ID", "conversion_group", "bl_dx", "FU_years",
             "years_to_AD", "years_to_MCI", "APOE4_Dosage", "Date", "YOB"],
    low_memory=False,
)

# Per-subject baseline age (for the age-as-timescale analyses). Only YOB is recorded,
# so birth date is taken as 1 July of YOB; age at the earliest visit is the entry age.
clin["Date"] = pd.to_datetime(clin["Date"], errors="coerce")
_bl = clin.dropna(subset=["Date"]).sort_values("Date").groupby("Patient_ID").first()
_birth = pd.to_datetime(pd.DataFrame(
    {"year": _bl["YOB"].astype(int), "month": 7, "day": 1}).set_index(_bl.index))
BL_AGE = ((_bl["Date"] - _birth) / np.timedelta64(1, "D") / 365.25)  # Series, idx=Patient_ID

sub = clin.drop_duplicates("Patient_ID").reset_index(drop=True)
print(f"  {len(clin)} rows -> {len(sub)} unique subjects")
print("  conversion_group counts:")
print(sub["conversion_group"].value_counts(dropna=False).to_string())

# APOE4 stratum definitions: ordered list of (group_value, legend_label, color).
DOSAGE_TO_3WAY = {0: "non-carrier", 1: "heterozygous", 2: "homozygous"}
STRATA_BINARY = [
    ("non-carrier", "Non-carrier", "#2c7bb6"),
    ("carrier",     "Carrier",     "#d7191c"),
]
STRATA_3WAY = [
    ("non-carrier",  "Non-carrier (0 e4)",  "#2c7bb6"),
    ("heterozygous", "Heterozygous (1 e4)", "#fdae61"),
    ("homozygous",   "Homozygous (2 e4)",   "#d7191c"),
]


# ═══════════════════════════════════════════════════════════════════════════════
# BUILD SURVIVAL DATA
# ═══════════════════════════════════════════════════════════════════════════════
def build_survival(sub_df, bl_diag, event_groups, event_time_cols, label):
    """
    Build per-subject survival data for one cohort.

    bl_diag         : baseline diagnosis defining the at-risk cohort — a single value
                      ("CN" / "MCI") or an iterable to POOL cohorts (["CN", "MCI"]).
    event_groups    : conversion_group value(s) that count as the event (str or iterable)
    event_time_cols : event-time column(s) (str or iterable); when several are given they
                      are coalesced in order (first non-null wins) — e.g.
                      ["years_to_MCI", "years_to_AD"] = time of first progression.

    Event subjects get the (coalesced) event time; everyone else is censored at FU_years.
    Returns a per-subject DataFrame: Patient_ID, time_years, event (1/0), apoe4, apoe4_3way.
    """
    if isinstance(event_groups, str):
        event_groups = [event_groups]
    if isinstance(event_time_cols, str):
        event_time_cols = [event_time_cols]

    bl_list = [bl_diag] if isinstance(bl_diag, str) else list(bl_diag)
    cohort = sub_df[sub_df["bl_dx"].isin(bl_list)].copy()
    cohort["event"] = cohort["conversion_group"].isin(event_groups).astype(int)
    etime = cohort[event_time_cols[0]]
    for c in event_time_cols[1:]:
        etime = etime.fillna(cohort[c])
    cohort["time_years"] = np.where(cohort["event"] == 1, etime, cohort["FU_years"]).astype(float)
    # Guard against missing/non-positive times.
    cohort["time_years"] = cohort["time_years"].fillna(EPS).clip(lower=EPS)
    # APOE4 strata: binary carrier status and 3-way genotype.
    cohort["apoe4"] = np.where(cohort["APOE4_Dosage"] >= 1, "carrier", "non-carrier")
    cohort["apoe4_3way"] = cohort["APOE4_Dosage"].map(DOSAGE_TO_3WAY)
    # Age timescale: enter the risk set at baseline age, exit at event/censoring age.
    cohort["entry_age"] = cohort["Patient_ID"].map(BL_AGE)
    cohort["exit_age"] = cohort["entry_age"] + cohort["time_years"]

    out = cohort[["Patient_ID", "time_years", "event", "apoe4", "apoe4_3way",
                  "entry_age", "exit_age"]].reset_index(drop=True)
    n_events = int(out["event"].sum())
    n_cens = len(out) - n_events
    print(f"\n  [{label}] baseline {'+'.join(bl_list)}: N={len(out)}, "
          f"events({'|'.join(event_groups)})={n_events}, censored={n_cens}, "
          f"median follow-up={out['time_years'].median():.1f} y")
    return out


surv_cn        = build_survival(sub, "CN",  "pCN_to_AD",  "years_to_AD",  "sCN vs pCN (CN->AD)")
surv_mci       = build_survival(sub, "MCI", "pMCI",       "years_to_AD",  "sMCI vs pMCI (MCI->AD)")
surv_pcn_mci   = build_survival(sub, "CN",  "pCN_to_MCI", "years_to_MCI", "pCN -> MCI")
surv_pcn_mci_ad = build_survival(
    sub, "CN", ["pCN_to_MCI", "pCN_to_AD"], ["years_to_MCI", "years_to_AD"],
    "pCN -> MCI or AD (any progression)")
# Pooled baseline CN+MCI, event = conversion to AD (pCN_to_AD or pMCI). CN -> MCI-only
# converters (pCN_to_MCI) are NOT events here — they stay at risk and are censored at FU.
surv_cn_mci = build_survival(
    sub, ["CN", "MCI"], ["pCN_to_AD", "pMCI"], "years_to_AD",
    "CN+MCI -> AD (pooled, excl. CN->MCI-only)")


# ═══════════════════════════════════════════════════════════════════════════════
# KAPLAN-MEIER ESTIMATOR (manual implementation, reused from 02f)
# ═══════════════════════════════════════════════════════════════════════════════
def kaplan_meier(times, events):
    """
    Compute KM survival curve.
    Returns arrays: unique_times, survival_prob, n_at_risk, n_events
    """
    df = pd.DataFrame({"time": times, "event": events}).sort_values("time")
    unique_times = sorted(df["time"].unique())

    surv = 1.0
    km_times = [0.0]
    km_surv  = [1.0]
    km_risk  = [len(df)]
    km_events = [0]

    for t in unique_times:
        n_at_risk = int(df[df["time"] >= t].shape[0])
        d = int(df[(df["time"] == t) & (df["event"] == 1)].shape[0])
        if n_at_risk > 0 and d > 0:
            surv *= (1 - d / n_at_risk)
        km_times.append(t)
        km_surv.append(surv)
        km_risk.append(n_at_risk)
        km_events.append(d)

    return np.array(km_times), np.array(km_surv), np.array(km_risk), np.array(km_events)


def logrank_test(surv_df, group_col, groups):
    """
    k-sample log-rank test (Mantel-Cox) across the ordered `groups` in `group_col`.
    Returns (chi2_stat, p_value, dof). Uses the multivariate hypergeometric
    observed-minus-expected vector and its covariance at each distinct event time;
    the statistic uses the first k-1 components (the k-th is linearly dependent).
    """
    k = len(groups)
    subs = [surv_df[surv_df[group_col] == g] for g in groups]
    event_times = sorted(surv_df.loc[surv_df["event"] == 1, "time_years"].unique())

    O = np.zeros(k)
    E = np.zeros(k)
    V = np.zeros((k, k))
    for t in event_times:
        n_j = np.array([float((s["time_years"] >= t).sum()) for s in subs])
        d_j = np.array([float(((s["time_years"] == t) & (s["event"] == 1)).sum())
                        for s in subs])
        n, d = n_j.sum(), d_j.sum()
        if n < 2 or d == 0:
            continue
        p = n_j / n
        O += d_j
        E += d * p
        V += (d * (n - d) / (n - 1)) * (np.diag(p) - np.outer(p, p))

    OE = (O - E)[:k - 1]
    Vr = V[:k - 1, :k - 1]
    dof = k - 1
    try:
        stat = float(OE @ np.linalg.solve(Vr, OE))
    except np.linalg.LinAlgError:
        stat = float(OE @ np.linalg.pinv(Vr) @ OE)
    return stat, float(chi2.sf(stat, df=dof)), dof


# ═══════════════════════════════════════════════════════════════════════════════
# AGE-AS-TIMESCALE: left-truncated (delayed-entry) Kaplan-Meier
#   Each subject is at risk only between their baseline age (entry) and their
#   event/censoring age (exit). The risk set at age t is everyone enrolled who is
#   currently in the study at that age: {i : entry_i <= t <= exit_i}.
# ═══════════════════════════════════════════════════════════════════════════════
def kaplan_meier_lt(entry, exit_, event, age_start):
    """
    Left-truncated KM on the age axis, conditioned at `age_start`: S = 1 at age_start
    and steps only at event ages >= age_start (i.e. conditional on being event-free and
    still in study at age_start). This avoids the unstable tiny-risk-set region below it.
    """
    entry = np.asarray(entry, float)
    exit_ = np.asarray(exit_, float)
    event = np.asarray(event, int)
    event_ages = sorted(t for t in set(exit_[event == 1]) if t >= age_start)

    S = 1.0
    ages = [age_start]
    surv = [1.0]
    for t in event_ages:
        n_at_risk = int(np.sum((entry <= t) & (exit_ >= t)))
        d = int(np.sum((exit_ == t) & (event == 1)))
        if n_at_risk > 0 and d > 0:
            S *= (1 - d / n_at_risk)
        ages.append(t)
        surv.append(S)
    return np.array(ages), np.array(surv)


def n_at_risk_age(entry, exit_, a):
    """Number currently in the study at age `a` (left-truncated risk set)."""
    return int(np.sum((np.asarray(entry, float) <= a) & (np.asarray(exit_, float) >= a)))


def logrank_test_age(surv_df, group_col, groups, age_start):
    """k-sample log-rank on the age axis (left-truncated), conditioned at `age_start`."""
    k = len(groups)
    subs = [surv_df[surv_df[group_col] == g] for g in groups]
    event_ages = sorted(a for a in surv_df.loc[surv_df["event"] == 1, "exit_age"].unique()
                        if a >= age_start)

    O = np.zeros(k)
    E = np.zeros(k)
    V = np.zeros((k, k))
    for t in event_ages:
        n_j = np.array([float(((s["entry_age"] <= t) & (s["exit_age"] >= t)).sum())
                        for s in subs])
        d_j = np.array([float(((s["exit_age"] == t) & (s["event"] == 1)).sum())
                        for s in subs])
        n, d = n_j.sum(), d_j.sum()
        if n < 2 or d == 0:
            continue
        p = n_j / n
        O += d_j
        E += d * p
        V += (d * (n - d) / (n - 1)) * (np.diag(p) - np.outer(p, p))

    OE = (O - E)[:k - 1]
    Vr = V[:k - 1, :k - 1]
    dof = k - 1
    try:
        stat = float(OE @ np.linalg.solve(Vr, OE))
    except np.linalg.LinAlgError:
        stat = float(OE @ np.linalg.pinv(Vr) @ OE)
    return stat, float(chi2.sf(stat, df=dof)), dof


# ═══════════════════════════════════════════════════════════════════════════════
# PLOTTING — survival S(t), falling from 1.0
# ═══════════════════════════════════════════════════════════════════════════════
def plot_km_survival(surv_df, title, filename, out_dir, max_years, color="#2c7bb6",
                     ylabel="Survival probability (AD-free)"):
    """
    Kaplan-Meier SURVIVAL plot (event-free proportion) with number-at-risk table below.
    """
    times = surv_df["time_years"].values
    events = surv_df["event"].values

    km_t, km_s, km_r, km_e = kaplan_meier(times, events)

    fig, (ax_km, ax_risk) = plt.subplots(
        2, 1, figsize=(10, 6.5),
        gridspec_kw={"height_ratios": [4, 1], "hspace": 0.05})

    # Survival step curve
    ax_km.step(km_t, km_s, where="post", color=color, linewidth=2.0)
    ax_km.fill_between(km_t, km_s, 1.0, step="post", alpha=0.10, color=color)

    # Mark censoring ticks (on the curve, at each censored subject's time)
    cens = surv_df[surv_df["event"] == 0]["time_years"].values
    for ct in cens:
        if ct > max_years:
            continue
        idx = np.searchsorted(km_t, ct, side="right") - 1
        idx = min(max(idx, 0), len(km_s) - 1)
        ax_km.plot(ct, km_s[idx], "|", color="black", markersize=8, markeredgewidth=1.2)

    ax_km.set_ylabel(ylabel, fontsize=11, fontweight="bold")
    ax_km.set_ylim(-0.02, 1.05)
    ax_km.set_xlim(-0.3, max_years + 0.5)
    ax_km.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax_km.axhline(0.5, color="grey", linewidth=0.6, linestyle="--", alpha=0.5)
    ax_km.grid(axis="y", alpha=0.3, linewidth=0.5)
    ax_km.spines["top"].set_visible(False)
    ax_km.spines["right"].set_visible(False)

    # X-axis ticks: every year
    year_ticks = list(range(0, max_years + 1))
    ax_km.set_xticks(year_ticks)
    ax_km.set_xticklabels([str(y) for y in year_ticks], fontsize=9)
    ax_km.set_xlabel("")

    # Summary stats
    n_total = len(surv_df)
    n_events = int(surv_df["event"].sum())
    n_censored = n_total - n_events
    median_surv = None
    for i in range(len(km_s)):
        if km_s[i] <= 0.5:
            median_surv = km_t[i]
            break

    info_text = f"N = {n_total}   Events = {n_events}   Censored = {n_censored}"
    if median_surv is not None:
        info_text += f"\nMedian survival = {median_surv:.1f} years"
    else:
        info_text += "\nMedian survival: not reached"
    ax_km.text(0.02, 0.10, info_text, transform=ax_km.transAxes,
               ha="left", va="bottom", fontsize=9, fontstyle="italic",
               bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                         edgecolor="#cccccc", alpha=0.9))

    # Number-at-risk table (subjects still being followed at each year)
    risk_years = list(range(0, max_years + 1, 2))
    risk_counts = [int(surv_df[surv_df["time_years"] >= yr].shape[0]) for yr in risk_years]

    ax_risk.set_xlim(ax_km.get_xlim())
    ax_risk.set_ylim(0, 1)
    ax_risk.axis("off")
    ax_risk.text(-0.02, 0.5, "At risk", transform=ax_risk.transAxes,
                 ha="right", va="center", fontsize=9, fontweight="bold")
    for yr, n in zip(risk_years, risk_counts):
        ax_risk.text(yr, 0.5, str(n), ha="center", va="center", fontsize=9,
                     transform=ax_risk.get_xaxis_transform())

    fig.text(0.5, 0.01, "Time (years)", ha="center", fontsize=11, fontweight="bold")

    fig.subplots_adjust(bottom=0.08)
    fig.savefig(out_dir / f"{filename}.png", bbox_inches="tight", dpi=300)
    fig.savefig(out_dir / f"{filename}.pdf", bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved {filename} -> {out_dir / f'{filename}.png'}")


def plot_km_survival_stratified(surv_df, group_col, strata, title, filename, out_dir,
                                max_years, ylabel="Survival probability (AD-free)"):
    """
    Overlaid KM survival curves across `strata` (ordered list of
    (group_value, legend_label, color)) found in column `group_col`. Adds a
    per-stratum number-at-risk table and a k-sample log-rank p-value.
    """
    k = len(strata)
    fig, (ax_km, ax_risk) = plt.subplots(
        2, 1, figsize=(10, 6.6 + 0.25 * k),
        gridspec_kw={"height_ratios": [4, 0.45 * k + 0.45], "hspace": 0.18})

    risk_years = list(range(0, max_years + 1, 2))
    legend_handles = []

    for grp, disp, color in strata:
        g = surv_df[surv_df[group_col] == grp]
        if len(g) == 0:
            continue
        km_t, km_s, _, _ = kaplan_meier(g["time_years"].values, g["event"].values)

        line, = ax_km.step(km_t, km_s, where="post", color=color, linewidth=2.0)
        for ct in g[g["event"] == 0]["time_years"].values:
            if ct > max_years:
                continue
            idx = min(max(np.searchsorted(km_t, ct, side="right") - 1, 0), len(km_s) - 1)
            ax_km.plot(ct, km_s[idx], "|", color=color, markersize=7, markeredgewidth=1.1)

        legend_handles.append(
            (line, f"{disp}  (N={len(g)}, events={int(g['event'].sum())})"))

    ax_km.set_ylabel(ylabel, fontsize=11, fontweight="bold")
    ax_km.set_ylim(-0.02, 1.05)
    ax_km.set_xlim(-0.3, max_years + 0.5)
    ax_km.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax_km.axhline(0.5, color="grey", linewidth=0.6, linestyle="--", alpha=0.5)
    ax_km.grid(axis="y", alpha=0.3, linewidth=0.5)
    ax_km.spines["top"].set_visible(False)
    ax_km.spines["right"].set_visible(False)
    ax_km.set_xticks(list(range(0, max_years + 1)))
    ax_km.set_xticklabels([str(y) for y in range(0, max_years + 1)], fontsize=9)

    ax_km.legend([h for h, _ in legend_handles], [l for _, l in legend_handles],
                 loc="lower left", fontsize=9, frameon=True, framealpha=0.9,
                 edgecolor="#cccccc")

    # k-sample log-rank test across the strata
    stat, p, dof = logrank_test(surv_df, group_col, [s[0] for s in strata])
    p_str = "p < 0.001" if p < 1e-3 else f"p = {p:.3f}"
    ax_km.text(0.98, 0.95, f"Log-rank {p_str}\n($\\chi^2$ = {stat:.2f}, df = {dof})",
               transform=ax_km.transAxes, ha="right", va="top", fontsize=9,
               fontstyle="italic",
               bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                         edgecolor="#cccccc", alpha=0.9))

    # Number-at-risk table: one row per stratum
    ax_risk.set_xlim(ax_km.get_xlim())
    ax_risk.set_ylim(0, 1)
    ax_risk.axis("off")
    ax_risk.text(-0.02, 0.97, "Number at risk", transform=ax_risk.transAxes,
                 ha="right", va="top", fontsize=9, fontweight="bold")
    row_ys = np.linspace(0.72, 0.12, k)
    for (grp, disp, color), ry in zip(strata, row_ys):
        g = surv_df[surv_df[group_col] == grp]
        ax_risk.text(-0.02, ry, grp.capitalize(), transform=ax_risk.transAxes,
                     ha="right", va="center", fontsize=8.5, fontweight="bold", color=color)
        for yr in risk_years:
            n = int((g["time_years"] >= yr).sum())
            ax_risk.text(yr, ry, str(n), ha="center", va="center", fontsize=8.5,
                         color=color, transform=ax_risk.get_xaxis_transform())

    fig.text(0.5, 0.01, "Time (years)", ha="center", fontsize=11, fontweight="bold")
    fig.subplots_adjust(bottom=0.08)
    fig.savefig(out_dir / f"{filename}.png", bbox_inches="tight", dpi=300)
    fig.savefig(out_dir / f"{filename}.pdf", bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved {filename} -> {out_dir / f'{filename}.png'}  "
          f"(log-rank {p_str}, df={dof})")


AGE_NOTE = ("Age timescale with left truncation (delayed entry): at each age the risk set "
            "is the enrolled subjects still in study at that age. Curves are conditioned on "
            "being event-free and in study at the floor age (start of x-axis), since risk "
            "sets below it are too small to estimate reliably.")


def plot_km_age(surv_df, title, filename, out_dir, age_lo, age_hi, color="#2c7bb6",
                ylabel="Survival probability (AD-free)"):
    """Left-truncated KM survival on the AGE axis, with number-at-risk table below."""
    entry = surv_df["entry_age"].values
    exit_ = surv_df["exit_age"].values
    event = surv_df["event"].values

    km_a, km_s = kaplan_meier_lt(entry, exit_, event, age_lo)

    fig, (ax_km, ax_risk) = plt.subplots(
        2, 1, figsize=(10, 6.5),
        gridspec_kw={"height_ratios": [4, 1], "hspace": 0.05})

    ax_km.step(km_a, km_s, where="post", color=color, linewidth=2.0)
    ax_km.fill_between(km_a, km_s, 1.0, step="post", alpha=0.10, color=color)

    # Censoring ticks at each censored subject's exit age
    for ca in surv_df[surv_df["event"] == 0]["exit_age"].values:
        if ca < age_lo or ca > age_hi:
            continue
        idx = min(max(np.searchsorted(km_a, ca, side="right") - 1, 0), len(km_s) - 1)
        ax_km.plot(ca, km_s[idx], "|", color="black", markersize=8, markeredgewidth=1.2)

    ax_km.set_ylabel(ylabel, fontsize=11, fontweight="bold")
    ax_km.set_ylim(-0.02, 1.05)
    ax_km.set_xlim(age_lo - 0.5, age_hi + 0.5)
    ax_km.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax_km.axhline(0.5, color="grey", linewidth=0.6, linestyle="--", alpha=0.5)
    ax_km.grid(axis="y", alpha=0.3, linewidth=0.5)
    ax_km.spines["top"].set_visible(False)
    ax_km.spines["right"].set_visible(False)

    age_ticks = list(range(age_lo, age_hi + 1, 5))
    ax_km.set_xticks(age_ticks)
    ax_km.set_xticklabels([str(a) for a in age_ticks], fontsize=9)

    cond = surv_df[surv_df["exit_age"] >= age_lo]   # conditional cohort: in study at >= floor
    n_total = len(cond)
    n_events = int(cond["event"].sum())
    median_age = next((km_a[i] for i in range(len(km_s)) if km_s[i] <= 0.5), None)
    info_text = (f"N = {n_total} in study at {age_lo}+   Events = {n_events}   "
                 f"Censored = {n_total - n_events}")
    info_text += (f"\nMedian age at event = {median_age:.1f} y" if median_age is not None
                  else "\nMedian age at event: not reached")
    ax_km.text(0.02, 0.10, info_text, transform=ax_km.transAxes,
               ha="left", va="bottom", fontsize=9, fontstyle="italic",
               bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                         edgecolor="#cccccc", alpha=0.9))

    ax_risk.set_xlim(ax_km.get_xlim())
    ax_risk.set_ylim(0, 1)
    ax_risk.axis("off")
    ax_risk.text(-0.02, 0.5, "At risk", transform=ax_risk.transAxes,
                 ha="right", va="center", fontsize=9, fontweight="bold")
    for a in age_ticks:
        ax_risk.text(a, 0.5, str(n_at_risk_age(entry, exit_, a)),
                     ha="center", va="center", fontsize=9,
                     transform=ax_risk.get_xaxis_transform())

    fig.text(0.5, 0.015, "Age (years)", ha="center", fontsize=11, fontweight="bold")
    fig.text(0.5, -0.03, AGE_NOTE, ha="center", fontsize=7, fontstyle="italic",
             color="#555555", wrap=True)
    fig.subplots_adjust(bottom=0.08)
    fig.savefig(out_dir / f"{filename}.png", bbox_inches="tight", dpi=300)
    fig.savefig(out_dir / f"{filename}.pdf", bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved {filename} -> {out_dir / f'{filename}.png'}")


def plot_km_age_stratified(surv_df, group_col, strata, title, filename, out_dir,
                           age_lo, age_hi, ylabel="Survival probability (AD-free)"):
    """Left-truncated KM survival on the AGE axis, stratified, with k-sample log-rank."""
    k = len(strata)
    fig, (ax_km, ax_risk) = plt.subplots(
        2, 1, figsize=(10, 6.6 + 0.25 * k),
        gridspec_kw={"height_ratios": [4, 0.45 * k + 0.45], "hspace": 0.18})

    age_ticks = list(range(age_lo, age_hi + 1, 5))
    legend_handles = []

    for grp, disp, color in strata:
        g = surv_df[surv_df[group_col] == grp]
        if len(g) == 0:
            continue
        entry, exit_, event = g["entry_age"].values, g["exit_age"].values, g["event"].values
        km_a, km_s = kaplan_meier_lt(entry, exit_, event, age_lo)
        line, = ax_km.step(km_a, km_s, where="post", color=color, linewidth=2.0)
        for ca in g[g["event"] == 0]["exit_age"].values:
            if ca < age_lo or ca > age_hi:
                continue
            idx = min(max(np.searchsorted(km_a, ca, side="right") - 1, 0), len(km_s) - 1)
            ax_km.plot(ca, km_s[idx], "|", color=color, markersize=7, markeredgewidth=1.1)
        g_cond = g[g["exit_age"] >= age_lo]   # conditional cohort: in study at >= floor
        legend_handles.append(
            (line, f"{disp}  (N={len(g_cond)}, events={int(g_cond['event'].sum())})"))

    ax_km.set_ylabel(ylabel, fontsize=11, fontweight="bold")
    ax_km.set_ylim(-0.02, 1.05)
    ax_km.set_xlim(age_lo - 0.5, age_hi + 0.5)
    ax_km.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax_km.axhline(0.5, color="grey", linewidth=0.6, linestyle="--", alpha=0.5)
    ax_km.grid(axis="y", alpha=0.3, linewidth=0.5)
    ax_km.spines["top"].set_visible(False)
    ax_km.spines["right"].set_visible(False)
    ax_km.set_xticks(age_ticks)
    ax_km.set_xticklabels([str(a) for a in age_ticks], fontsize=9)

    ax_km.legend([h for h, _ in legend_handles], [l for _, l in legend_handles],
                 loc="lower left", fontsize=9, frameon=True, framealpha=0.9,
                 edgecolor="#cccccc")

    stat, p, dof = logrank_test_age(surv_df, group_col, [s[0] for s in strata], age_lo)
    p_str = "p < 0.001" if p < 1e-3 else f"p = {p:.3f}"
    ax_km.text(0.98, 0.95, f"Log-rank {p_str}\n($\\chi^2$ = {stat:.2f}, df = {dof})",
               transform=ax_km.transAxes, ha="right", va="top", fontsize=9,
               fontstyle="italic",
               bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                         edgecolor="#cccccc", alpha=0.9))

    ax_risk.set_xlim(ax_km.get_xlim())
    ax_risk.set_ylim(0, 1)
    ax_risk.axis("off")
    ax_risk.text(-0.02, 0.97, "Number at risk", transform=ax_risk.transAxes,
                 ha="right", va="top", fontsize=9, fontweight="bold")
    row_ys = np.linspace(0.72, 0.12, k)
    for (grp, disp, color), ry in zip(strata, row_ys):
        g = surv_df[surv_df[group_col] == grp]
        ax_risk.text(-0.02, ry, grp.capitalize(), transform=ax_risk.transAxes,
                     ha="right", va="center", fontsize=8.5, fontweight="bold", color=color)
        for a in age_ticks:
            ax_risk.text(a, ry, str(n_at_risk_age(g["entry_age"].values, g["exit_age"].values, a)),
                         ha="center", va="center", fontsize=8.5, color=color,
                         transform=ax_risk.get_xaxis_transform())

    fig.text(0.5, 0.015, "Age (years)", ha="center", fontsize=11, fontweight="bold")
    fig.text(0.5, -0.03, AGE_NOTE, ha="center", fontsize=7, fontstyle="italic",
             color="#555555", wrap=True)
    fig.subplots_adjust(bottom=0.08)
    fig.savefig(out_dir / f"{filename}.png", bbox_inches="tight", dpi=300)
    fig.savefig(out_dir / f"{filename}.pdf", bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved {filename} -> {out_dir / f'{filename}.png'}  "
          f"(log-rank {p_str}, df={dof})")


# ── Generate all KM curves ─────────────────────────────────────────────────────
all_times = pd.concat(
    [surv_cn, surv_mci, surv_pcn_mci, surv_pcn_mci_ad, surv_cn_mci])["time_years"]
MAX_YEARS = int(np.ceil(all_times.max()))
print(f"\nMax follow-up across cohorts: {MAX_YEARS} years")

# One spec per analysis: (surv_df, slug, pooled_color, ylabel, pooled_title, strat_title)
ANALYSES = [
    (surv_cn, "km_pcn_to_ad", "#1a9641", "Survival probability (AD-free)",
     "sCN vs pCN — AD-free survival (baseline CN)",
     "sCN vs pCN — AD-free survival (baseline CN)"),
    (surv_mci, "km_pmci_to_ad", "#2c7bb6", "Survival probability (AD-free)",
     "sMCI vs pMCI — AD-free survival (baseline MCI)",
     "sMCI vs pMCI — AD-free survival (baseline MCI)"),
    (surv_pcn_mci, "km_pcn_to_mci", "#8856a7", "Survival probability (MCI-free)",
     "CN → MCI conversion — MCI-free survival (baseline CN)",
     "CN → MCI conversion — MCI-free survival (baseline CN)"),
    (surv_pcn_mci_ad, "km_pcn_to_mci_or_ad", "#d95f0e",
     "Survival probability (progression-free)",
     "CN → MCI or AD — progression-free survival (baseline CN)",
     "CN → MCI or AD — progression-free survival (baseline CN)"),
    (surv_cn_mci, "km_cn_mci_to_ad", "#d7191c", "Survival probability (AD-free)",
     "CN + MCI → AD — AD-free survival (pooled baseline CN + MCI)",
     "CN + MCI → AD — AD-free survival (pooled baseline CN + MCI)"),
]

for surv_df, slug, color, ylabel, pooled_title, strat_title in ANALYSES:
    # Pooled
    plot_km_survival(surv_df, pooled_title, slug, OUT_DIR, MAX_YEARS,
                     color=color, ylabel=ylabel)
    # Stratified by APOE4 carrier status (binary)
    plot_km_survival_stratified(
        surv_df, "apoe4", STRATA_BINARY,
        f"{strat_title}\nby APOE4 status",
        f"{slug}_by_apoe4", OUT_DIR, MAX_YEARS, ylabel=ylabel)
    # Stratified by APOE4 genotype (3-way)
    plot_km_survival_stratified(
        surv_df, "apoe4_3way", STRATA_3WAY,
        f"{strat_title}\nby APOE4 genotype",
        f"{slug}_by_apoe4_3way", OUT_DIR, MAX_YEARS, ylabel=ylabel)


# ── Age-as-timescale versions (left-truncated / delayed entry) ─────────────────
# Condition at a floor age where risk sets are populated (ages below it have too few
# enrolled subjects to estimate reliably); cap the axis at the oldest observed exit age.
AGE_LO = 65
exit_max = max(s["exit_age"].max() for s, *_ in ANALYSES)
AGE_HI = int(np.ceil(exit_max / 5) * 5)
print(f"\nAge timescale: conditioned at {AGE_LO}, axis {AGE_LO}–{AGE_HI} y  -> {AGE_OUT_DIR}")

for surv_df, slug, color, ylabel, pooled_title, strat_title in ANALYSES:
    plot_km_age(surv_df, f"{pooled_title}\n(age timescale)",
                f"{slug}_age", AGE_OUT_DIR, AGE_LO, AGE_HI, color=color, ylabel=ylabel)
    plot_km_age_stratified(
        surv_df, "apoe4", STRATA_BINARY,
        f"{strat_title}\nby APOE4 status (age timescale)",
        f"{slug}_by_apoe4_age", AGE_OUT_DIR, AGE_LO, AGE_HI, ylabel=ylabel)
    plot_km_age_stratified(
        surv_df, "apoe4_3way", STRATA_3WAY,
        f"{strat_title}\nby APOE4 genotype (age timescale)",
        f"{slug}_by_apoe4_3way_age", AGE_OUT_DIR, AGE_LO, AGE_HI, ylabel=ylabel)

print("\nDone.")
