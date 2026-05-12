"""
02f_plot_time_to_event_no_cdr_stratified.py
============================================
Time-to-event analysis for the no-CDR stratified cohort.

Produces:
  1. Summary table: yearly cumulative conversions, censored, at risk
  2. Kaplan-Meier survival curves with number-at-risk table
  3. Nelson-Aalen cumulative hazard plots

Two analyses:
  A) MCI → AD conversion (baseline MCI subjects)
  B) CN → MCI conversion (baseline CN subjects)

Data source: no_cdr_stratified verbose longitudinal master
              (uses Label_visit_diag per visit to detect first conversion)

Outputs (to outputs/baseline_no_cdr/time_to_event/):
  time_to_event_summary.csv
  km_mci_to_ad.{png,pdf}
  km_cn_to_mci.{png,pdf}
  cumhaz_mci_to_ad.{png,pdf}
  cumhaz_cn_to_mci.{png,pdf}

Run:
  python clinical_pipeline/02f_plot_time_to_event_no_cdr_stratified.py
"""

import re
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path

matplotlib.rcParams["font.family"] = "DejaVu Serif"

# ── Paths ──────────────────────────────────────────────────────────────────────
REPO_ROOT   = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
MASTER_CSV  = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified"
                   r"\verbose\longitudinal\master_clinical_verbose.csv")
SPLIT_DIR   = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified\tabular\baseline")
OUT_DIR     = REPO_ROOT / "clinical_pipeline" / "outputs" / "baseline_no_cdr" / "time_to_event"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
print("Loading master CSV...")
clin = pd.read_csv(MASTER_CSV,
                    usecols=["Patient_ID", "VISCODE_long", "Date",
                             "Label_bl_multi", "Label_visit_diag"],
                    low_memory=False)
clin["Date"] = pd.to_datetime(clin["Date"], errors="coerce")
clin = clin.dropna(subset=["Date"]).copy()
print(f"  {len(clin)} rows, {clin['Patient_ID'].nunique()} subjects")


def viscode_to_months(v):
    v = str(v).strip().lower()
    if v == "sc": return -1
    if v == "bl": return 0
    m = re.match(r"^m(\d+)$", v)
    return int(m.group(1)) if m else None


clin["_months"] = clin["VISCODE_long"].apply(viscode_to_months)
clin = clin.dropna(subset=["_months"]).copy()
clin["_months"] = clin["_months"].astype(int)

# ── Load seed 0 split for test-set labelling ──────────────────────────────────
split_map = {}
for split_name in ["train", "val", "test"]:
    df = pd.read_csv(SPLIT_DIR / "seed_0" / f"{split_name}.csv",
                     usecols=["Patient_ID"], low_memory=False)
    for pid in df["Patient_ID"].unique():
        split_map[pid] = split_name

clin["_split"] = clin["Patient_ID"].map(split_map)


# ═══════════════════════════════════════════════════════════════════════════════
# BUILD SURVIVAL DATA
# ═══════════════════════════════════════════════════════════════════════════════
def build_survival(clin_df, bl_diag, target_diag, label):
    """
    Build survival data for subjects with baseline diagnosis `bl_diag`
    converting to `target_diag`.

    Returns DataFrame with columns:
        Patient_ID, split, time_years, event (1=converted, 0=censored)
    """
    # Subjects with the specified baseline diagnosis
    bl_subjects = clin_df[clin_df["Label_bl_multi"] == bl_diag]["Patient_ID"].unique()
    print(f"\n  [{label}] {len(bl_subjects)} subjects with baseline {bl_diag}")

    records = []
    for pid in bl_subjects:
        subj = clin_df[clin_df["Patient_ID"] == pid].sort_values("_months")
        split = split_map.get(pid, "unknown")

        # Find first visit where diagnosis matches target
        conversion_row = subj[subj["Label_visit_diag"] == target_diag]

        if not conversion_row.empty:
            first_conv = conversion_row.iloc[0]
            time_months = first_conv["_months"]
            if time_months <= 0:
                time_months = 1  # avoid zero/negative time
            records.append({
                "Patient_ID": pid, "split": split,
                "time_months": time_months,
                "time_years": time_months / 12.0,
                "event": 1
            })
        else:
            # Censored: last observed visit
            last_visit = subj[subj["_months"] >= 0].iloc[-1] if len(subj[subj["_months"] >= 0]) > 0 else subj.iloc[-1]
            time_months = max(last_visit["_months"], 1)
            records.append({
                "Patient_ID": pid, "split": split,
                "time_months": time_months,
                "time_years": time_months / 12.0,
                "event": 0
            })

    surv_df = pd.DataFrame(records)
    n_events = surv_df["event"].sum()
    n_censored = len(surv_df) - n_events
    print(f"    Events: {n_events}, Censored: {n_censored}, "
          f"Median follow-up: {surv_df['time_years'].median():.1f} years")
    return surv_df


surv_mci_ad = build_survival(clin, "MCI", "AD", "MCI -> AD")
surv_cn_mci = build_survival(clin, "CN", "MCI", "CN -> MCI")


# ═══════════════════════════════════════════════════════════════════════════════
# SUMMARY TABLE: yearly cumulative conversions, censored, at risk
# ═══════════════════════════════════════════════════════════════════════════════
def yearly_summary(surv_df, label, max_year=15):
    """Build yearly summary: N at risk, cumulative events, censored in interval."""
    rows = []
    years = range(1, max_year + 1)
    total = len(surv_df)

    for yr in years:
        # Cumulative events by year yr
        cum_events = int(surv_df[(surv_df["event"] == 1) & (surv_df["time_years"] <= yr)].shape[0])
        # Censored in interval (yr-1, yr]: left study without converting
        censored_in_interval = int(surv_df[
            (surv_df["event"] == 0) &
            (surv_df["time_years"] > (yr - 1)) &
            (surv_df["time_years"] <= yr)
        ].shape[0])
        # At risk at start of year: still being followed at year yr-1
        at_risk = int(surv_df[surv_df["time_years"] > (yr - 1)].shape[0])

        if at_risk == 0:
            break  # no one left to follow

        rows.append({
            "Analysis": label,
            "Year": yr,
            "N_total": total,
            "At_Risk": at_risk,
            "Events_in_Year": int(surv_df[
                (surv_df["event"] == 1) &
                (surv_df["time_years"] > (yr - 1)) &
                (surv_df["time_years"] <= yr)
            ].shape[0]),
            "Cum_Events": cum_events,
            "Censored_in_Year": censored_in_interval,
            "Cum_Censored": int(surv_df[
                (surv_df["event"] == 0) & (surv_df["time_years"] <= yr)
            ].shape[0]),
        })
    return pd.DataFrame(rows)


summary_mci = yearly_summary(surv_mci_ad, "MCI -> AD")
summary_cn  = yearly_summary(surv_cn_mci, "CN -> MCI")
summary_all = pd.concat([summary_mci, summary_cn], ignore_index=True)

csv_path = OUT_DIR / "time_to_event_summary.csv"
summary_all.to_csv(csv_path, index=False)
print(f"\nSaved summary CSV → {csv_path}")
print(summary_all.to_string(index=False))


# ═══════════════════════════════════════════════════════════════════════════════
# KAPLAN-MEIER ESTIMATOR (manual implementation)
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


def nelson_aalen(times, events):
    """
    Compute Nelson-Aalen cumulative hazard.
    Returns arrays: unique_times, cumulative_hazard
    """
    df = pd.DataFrame({"time": times, "event": events}).sort_values("time")
    unique_times = sorted(df["time"].unique())

    H = 0.0
    na_times = [0.0]
    na_H     = [0.0]

    for t in unique_times:
        n_at_risk = int(df[df["time"] >= t].shape[0])
        d = int(df[(df["time"] == t) & (df["event"] == 1)].shape[0])
        if n_at_risk > 0:
            H += d / n_at_risk
        na_times.append(t)
        na_H.append(H)

    return np.array(na_times), np.array(na_H)


# ═══════════════════════════════════════════════════════════════════════════════
# PLOTTING
# ═══════════════════════════════════════════════════════════════════════════════
def plot_km(surv_df, title, filename, out_dir, max_years, color="#2c7bb6"):
    """
    Kaplan-Meier plot showing CONVERSION probability (1 - survival)
    with number-at-risk table below.
    """
    times = surv_df["time_years"].values
    events = surv_df["event"].values

    km_t, km_s, km_r, km_e = kaplan_meier(times, events)

    # Conversion probability = 1 - survival
    km_conv = 1.0 - km_s

    fig, (ax_km, ax_risk) = plt.subplots(
        2, 1, figsize=(10, 6.5),
        gridspec_kw={"height_ratios": [4, 1], "hspace": 0.05})

    # Conversion probability step plot
    ax_km.step(km_t, km_conv, where="post", color=color, linewidth=2.0)
    ax_km.fill_between(km_t, km_conv, step="post", alpha=0.12, color=color)

    # Mark censoring ticks
    cens = surv_df[surv_df["event"] == 0]["time_years"].values
    for ct in cens:
        if ct > max_years:
            continue
        idx = np.searchsorted(km_t, ct, side="right") - 1
        idx = min(idx, len(km_conv) - 1)
        ax_km.plot(ct, km_conv[idx], "|", color="black", markersize=8, markeredgewidth=1.2)

    ax_km.set_ylabel("Conversion Probability", fontsize=11, fontweight="bold")
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
        info_text += f"\nMedian time to conversion = {median_surv:.1f} years"
    else:
        info_text += f"\nMedian time to conversion: not reached"
    ax_km.text(0.02, 0.95, info_text, transform=ax_km.transAxes,
               ha="left", va="top", fontsize=9, fontstyle="italic",
               bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                         edgecolor="#cccccc", alpha=0.9))

    # Number-at-risk table (subjects still being followed at each year)
    risk_years = list(range(0, max_years + 1, 2))
    risk_counts = []
    for yr in risk_years:
        n = int(surv_df[surv_df["time_years"] >= yr].shape[0])
        risk_counts.append(n)

    ax_risk.set_xlim(ax_km.get_xlim())
    ax_risk.set_ylim(0, 1)
    ax_risk.axis("off")

    # "At risk" label — use axes transform to avoid overlap
    ax_risk.text(-0.02, 0.5, "At risk", transform=ax_risk.transAxes,
                 ha="right", va="center", fontsize=9, fontweight="bold")

    for yr, n in zip(risk_years, risk_counts):
        ax_risk.text(yr, 0.5, str(n), ha="center", va="center", fontsize=9,
                     transform=ax_risk.get_xaxis_transform())

    # X-axis label below risk table
    fig.text(0.5, 0.01, "Time (years)", ha="center", fontsize=11, fontweight="bold")

    fig.subplots_adjust(bottom=0.08)
    fig.savefig(out_dir / f"{filename}.png", bbox_inches="tight", dpi=300)
    fig.savefig(out_dir / f"{filename}.pdf", bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved {filename} → {out_dir / f'{filename}.png'}")


def plot_cumhaz(surv_df, title, filename, out_dir, max_years, color="#d7191c"):
    """
    Nelson-Aalen cumulative hazard plot.
    """
    times = surv_df["time_years"].values
    events = surv_df["event"].values

    na_t, na_H = nelson_aalen(times, events)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.step(na_t, na_H, where="post", color=color, linewidth=2.0)
    ax.fill_between(na_t, na_H, step="post", alpha=0.12, color=color)

    ax.set_xlabel("Time (years)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Cumulative Hazard", fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax.set_xlim(-0.3, max_years + 0.5)
    ax.set_xticks(list(range(0, max_years + 1)))
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", alpha=0.3, linewidth=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Summary
    n_total = len(surv_df)
    n_events = int(surv_df["event"].sum())
    ax.text(0.02, 0.95,
            f"N = {n_total}   Events = {n_events}   "
            f"Final Λ(t) = {na_H[-1]:.3f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=9,
            fontstyle="italic",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="#cccccc", alpha=0.9))

    plt.tight_layout()
    fig.savefig(out_dir / f"{filename}.png", bbox_inches="tight", dpi=300)
    fig.savefig(out_dir / f"{filename}.pdf", bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved {filename} → {out_dir / f'{filename}.png'}")


# ── Generate all plots ────────────────────────────────────────────────────────
# Compute the actual max clinical follow-up in integer years
all_times = pd.concat([surv_mci_ad, surv_cn_mci])["time_years"]
clinical_max_years = int(np.ceil(all_times.max()))
MRI_WINDOW_YEARS = 11   # m132 = 11 years (last visit with MRI correspondence)

plot_configs = [
    ("imaging_window",  MRI_WINDOW_YEARS,   "within MRI window (≤m132)"),
    ("clinical_window", clinical_max_years,  f"full clinical follow-up (≤m{clinical_max_years*12})"),
]

for subdir, max_yr, desc in plot_configs:
    sub_out = OUT_DIR / subdir
    sub_out.mkdir(parents=True, exist_ok=True)
    print(f"\n── Generating plots: {desc}  (max {max_yr} years) ──")

    plot_km(surv_mci_ad,
            f"Time to AD Conversion (baseline MCI) — {desc}",
            "km_mci_to_ad", sub_out, max_yr, color="#2c7bb6")

    plot_km(surv_cn_mci,
            f"Time to MCI Conversion (baseline CN) — {desc}",
            "km_cn_to_mci", sub_out, max_yr, color="#1a9641")

    plot_cumhaz(surv_mci_ad,
                f"Cumulative Hazard: MCI -> AD — {desc}",
                "cumhaz_mci_to_ad", sub_out, max_yr, color="#d7191c")

    plot_cumhaz(surv_cn_mci,
                f"Cumulative Hazard: CN -> MCI — {desc}",
                "cumhaz_cn_to_mci", sub_out, max_yr, color="#e66101")


# ═══════════════════════════════════════════════════════════════════════════════
# SUMMARY TABLE (matplotlib rendering)
# ═══════════════════════════════════════════════════════════════════════════════
def render_summary_table(summary_df, surv_dfs, labels, colors):
    """Render a publication-quality summary table."""
    n_analyses = len(labels)

    # Column layout: Analysis | Year | At Risk | Events | Cum Events | Censored | Cum Censored
    col_headers = ["Year", "At Risk", "Events", "Cum.\nEvents", "Censored", "Cum.\nCensored"]
    n_cols = len(col_headers)

    COL_W = [0.60, 0.80, 0.80, 0.80, 0.80, 0.90]
    ROW_H = 0.28
    TITLE_H = 0.45
    SUBTITLE_H = 0.30
    GROUP_H = 0.32
    HDR_H = 0.32
    PAD = 0.15

    # Count data rows per analysis
    analyses_data = []
    for label in labels:
        adf = summary_df[summary_df["Analysis"] == label]
        analyses_data.append(adf)

    total_data_rows = sum(len(a) for a in analyses_data)
    fig_w = 0.25 + sum(COL_W) + 0.25
    fig_h = PAD + TITLE_H + n_analyses * (SUBTITLE_H + HDR_H) + total_data_rows * ROW_H + PAD

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    ax.set_xlim(0, fig_w)
    ax.set_ylim(0, fig_h)

    LEFT = 0.25
    col_left = [LEFT]
    for w in COL_W:
        col_left.append(col_left[-1] + w)
    RIGHT = col_left[-1]
    col_cx = [(col_left[i] + col_left[i+1]) / 2 for i in range(n_cols)]

    def hline(y, lw=1.0, ls="-"):
        ax.plot([LEFT, RIGHT], [y, y], color="black", linewidth=lw, linestyle=ls,
                solid_capstyle="butt", zorder=3)

    TOP = fig_h - PAD
    y = TOP

    # Title
    hline(y, lw=1.5)
    y -= TITLE_H
    ax.text((LEFT + RIGHT) / 2, y + TITLE_H / 2,
            "Time-to-Event Summary\n(no-CDR stratified cohort, all subjects)",
            ha="center", va="center", fontsize=10, fontweight="bold", color="black",
            linespacing=1.3)
    hline(y, lw=1.2)

    for a_idx, (label, adf, color) in enumerate(zip(labels, analyses_data, colors)):
        surv = surv_dfs[a_idx]
        n_total = len(surv)
        n_events = int(surv["event"].sum())

        # Group header
        y_top = y
        y -= SUBTITLE_H
        ax.text((LEFT + RIGHT) / 2, (y_top + y) / 2,
                f"{label}   (N={n_total}, events={n_events}, "
                f"censored={n_total - n_events})",
                ha="center", va="center", fontsize=9, fontweight="bold",
                color=color)
        hline(y, lw=0.8)

        # Column headers
        y_top = y
        y -= HDR_H
        for j, h in enumerate(col_headers):
            ax.text(col_cx[j], (y_top + y) / 2, h,
                    ha="center", va="center", fontsize=8,
                    fontstyle="italic", color="black")
        hline(y, lw=0.6, ls=(0, (3, 3)))

        # Data rows
        for _, row in adf.iterrows():
            y_top = y
            y -= ROW_H
            vals = [str(int(row["Year"])),
                    str(int(row["At_Risk"])),
                    str(int(row["Events_in_Year"])),
                    str(int(row["Cum_Events"])),
                    str(int(row["Censored_in_Year"])),
                    str(int(row["Cum_Censored"]))]
            for j, v in enumerate(vals):
                ax.text(col_cx[j], (y_top + y) / 2, v,
                        ha="center", va="center", fontsize=8, color="black")

        hline(y, lw=1.2)

    # Box
    rect = plt.Rectangle((LEFT, y), RIGHT - LEFT, TOP - y,
                          facecolor="none", edgecolor="black", linewidth=1.5, zorder=5)
    ax.add_patch(rect)

    plt.tight_layout(pad=0.1)
    fig.savefig(OUT_DIR / "time_to_event_table.png", bbox_inches="tight", dpi=300)
    fig.savefig(OUT_DIR / "time_to_event_table.pdf", bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved table → {OUT_DIR / 'time_to_event_table.png'}")


render_summary_table(
    summary_all,
    [surv_mci_ad, surv_cn_mci],
    ["MCI -> AD", "CN -> MCI"],
    ["#2c7bb6", "#1a9641"])

print("\nDone.")
