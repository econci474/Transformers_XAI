r"""
03_km_calibration_plots.py  (env: survml — lifelines + matplotlib)
==================================================================
Qualitative survival-curve / calibration figures for the BEST MAMBA survival heads, so we can SEE
what the parametrised curves look like and how well they track the data.

For each chosen head config (default: the two deep heads on the SAME mamba1-frozen backbone —
Weibull-piecewise [= A_default_mamba1_frozen] vs Cox-PH [= B_cox]) we overlay the model's MARGINAL
predicted survival S(t) (= mean over patients of the per-patient exported S, averaged across the 3
seed models) on the EMPIRICAL Kaplan-Meier curve (lifelines, 95% CI band):

  Fig 1 (average) : whole CN+MCI cohort — KM  vs  each head's marginal S(t).
  Fig 2 (by APOE4): faceted by APOE4 carrier status (0 copies vs >=1) — KM per stratum vs each head's
                    marginal S(t) within that stratum, so you can see whether the heads reproduce the
                    APOE4 risk separation the data show.

Population = whole cohort by default (535; smoothest, most representative for a QUALITATIVE shape view).
NB those predictions are largely IN-SAMPLE (mean over 3 seed models, each trained on ~427/535) — this is
a curve-shape / marginal-calibration figure, NOT a held-out generalisation claim (that's stage 02's
C-index/IPCW/IBS table). Pass --split test to redraw on the held-out test patients (noisier, n~55/seed).

Inputs : <pred_root>/<config>/seed_*/predictions.parquet (+ meta.json grid_times)   [stage 01 output]
         02l SPLIT_DIR/seed_0/{train,val,test}.csv  for APOE4_Dosage + Patient_ID   [canonical splits]
Outputs: <pred_root>/comparison/{km_calibration_average.png/pdf,
         km_calibration_by_apoe4.png/pdf, km_calibration_curves.csv}

Run:  conda run -n survml python integration/T4/mamba_long_clinical/03_km_calibration_plots.py \
        --pred_root D:/ADNI_BIDS_project/derivatives/mamba_survival/survival
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

HERE = Path(__file__).resolve().parent
PRED_ROOT = HERE / "outputs" / "survival"
SEEDS = (0, 1, 2)
TMAX = 15.0                                                        # model S grid is 0..15 yr
# canonical split CSVs (carry APOE4_Dosage + Patient_ID); seed_0 train/val/test partitions the cohort
SPLIT_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion"
                 r"\tabular\baseline")
CLIN = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\clinical_pipeline")


def _load_02l():
    spec = importlib.util.spec_from_file_location("l02", CLIN / "02l_survival_model_comparison.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


# the "best method" head comparison: same mamba1-frozen backbone, two parametrisations
DEFAULT_CONFIGS = ["A_default_mamba1_frozen", "B_cox"]
TAB_COX = "TAB_coxph"                                              # tabular Cox PH baseline (02l, 8 feat)
CONFIG_LABEL = {"A_default_mamba1_frozen": "MAMBA → Weibull-piecewise",
                "B_cox": "MAMBA → Cox-PH",
                "B_weibull_aft": "MAMBA → Weibull-AFT",
                "B_soden": "MAMBA → SODEN",
                "A_ctrl_meanpool": "mean-pool → Weibull-piecewise",
                TAB_COX: "Cox PH (MMSE, ADAS13, MoCA, FAQ, age, sex, educ, APOE4)"}
CONFIG_COLOR = {"A_default_mamba1_frozen": "#1b7837", "B_cox": "#762a83",
                "B_weibull_aft": "#b35806", "B_soden": "#2166ac",
                "A_ctrl_meanpool": "#d6604d", TAB_COX: "#e6550d"}


def load_apoe():
    """Patient_ID -> APOE4_Dosage from the canonical seed_0 splits (covers the whole cohort)."""
    frames = []
    for sp in ("train", "val", "test"):
        df = pd.read_csv(SPLIT_DIR / "seed_0" / f"{sp}.csv", low_memory=False,
                         usecols=lambda c: c in ("Patient_ID", "APOE4_Dosage"))
        frames.append(df)
    a = pd.concat(frames, ignore_index=True)
    a["Patient_ID"] = a["Patient_ID"].astype(str)
    a["APOE4_Dosage"] = pd.to_numeric(a["APOE4_Dosage"], errors="coerce")
    return a.dropna(subset=["APOE4_Dosage"]).drop_duplicates("Patient_ID").set_index("Patient_ID")["APOE4_Dosage"]


def load_config(name, split):
    """Per-patient mean predicted S (averaged across the 3 seeds) + event/time, for one head config."""
    grid = None
    per_seed = []
    meta_time = {}
    for seed in SEEDS:
        d = PRED_ROOT / name / f"seed_{seed}"
        if not (d / "predictions.parquet").exists():
            return None
        grid = np.asarray(json.load(open(d / "meta.json"))["grid_times"], float)
        df = pd.read_parquet(d / "predictions.parquet")
        df["Patient_ID"] = df["Patient_ID"].astype(str)
        if split != "all":
            keep = ["val", "test"] if split == "valtest" else [split]
            df = df[df["split"].isin(keep)]
        scol = [c for c in df.columns if c.startswith("S_")]
        per_seed.append(df.set_index("Patient_ID")[scol])
        for pid, tt, ee in zip(df["Patient_ID"], df["time"], df["event"]):
            meta_time[pid] = (float(tt), float(ee))                # patient-level, identical across seeds
    # average S across the seeds that contain each patient
    allS = pd.concat(per_seed).groupby(level=0).mean()
    pid = allS.index.to_numpy()
    S = allS.to_numpy(float)
    time = np.array([meta_time[p][0] for p in pid])
    event = np.array([meta_time[p][1] for p in pid]) > 0.5
    return {"pid": pid, "S": S, "grid": grid, "time": time, "event": event}


def km_curve(time, event):
    kmf = KaplanMeierFitter().fit(time, event_observed=event, label="km")
    t = kmf.survival_function_.index.to_numpy()
    s = kmf.survival_function_["km"].to_numpy()
    ci = kmf.confidence_interval_
    lo = ci.iloc[:, 0].to_numpy(); hi = ci.iloc[:, 1].to_numpy()
    m = t <= TMAX
    return t[m], s[m], lo[m], hi[m], kmf


def _draw(ax, time, event, curves, title):
    """curves = list of (label, color, grid, marginal_S, ls). Draws KM (black, CI) + each model curve."""
    t, s, lo, hi, _ = km_curve(time, event)
    ax.step(t, s, where="post", color="black", lw=2.4, label=f"Ground Truth KM (n={len(time)}, "
                                                             f"{int(event.sum())} ev)", zorder=3)
    ax.fill_between(t, lo, hi, step="post", color="black", alpha=0.12, lw=0, zorder=1,
                    label="95% CI of KM estimate (Greenwood)")
    for label, color, grid, mS, ls in curves:
        ax.plot(grid, mS, color=color, lw=2.2, ls=ls, label=label, zorder=4)
    ax.set_xlim(0, TMAX); ax.set_ylim(0, 1.02)
    ax.set_xlabel("Years from baseline"); ax.set_ylabel("Survival  S(t) = P(not yet AD)")
    ax.set_title(title, fontsize=11)
    ax.grid(alpha=0.25); ax.legend(loc="lower left", fontsize=8, framealpha=0.9)


def tabular_cox_meanS(pid_order, grid):
    """Vanilla tabular Cox PH (02l, 8 compact features) per-patient mean S over seeds, aligned to
    pid_order. Each patient is scored by every seed's Cox fit (train-fitted) and averaged."""
    L = _load_02l()
    spec = L.make_model_specs(["Cox PH"])["Cox PH"]
    acc = {}
    for seed in SEEDS:
        dd = L.load_seed(seed, L.COMPACT_FEATURES)
        Xtr, etr, ttr = dd["train"]
        sf = spec(Xtr, etr, ttr)                                  # surv_fn(X, times) -> [n, T]
        for sp in ("train", "val", "test"):
            X = dd[sp][0]
            S = np.asarray(sf(X, grid), float)
            csv = pd.read_csv(SPLIT_DIR / f"seed_{seed}" / f"{sp}.csv", low_memory=False)
            csv = csv[csv["bl_dx"].isin(["CN", "MCI"])].copy()
            csv["Patient_ID"] = csv["Patient_ID"].astype(str)
            pids = csv.loc[X.index, "Patient_ID"].to_numpy()       # X.index == filtered csv index
            for i, p in enumerate(pids):
                acc.setdefault(p, []).append(S[i])
    mean = {p: np.mean(v, axis=0) for p, v in acc.items()}
    return np.array([mean[p] if p in mean else np.full(len(grid), np.nan) for p in pid_order])


def main():
    global PRED_ROOT
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred_root", type=str, default=str(PRED_ROOT))
    ap.add_argument("--configs", type=str, nargs="*", default=DEFAULT_CONFIGS,
                    help="head configs to overlay (default: Weibull-piecewise + Cox, mamba1-frozen)")
    ap.add_argument("--split", type=str, default="all", choices=["all", "test", "valtest"],
                    help="patient population for BOTH KM and the marginal S (default whole cohort)")
    ap.add_argument("--no_tabular_cox", action="store_true",
                    help="omit the vanilla tabular Cox PH (8-feat 02l baseline) overlay")
    args = ap.parse_args()
    PRED_ROOT = Path(args.pred_root)
    out = PRED_ROOT / "comparison"; out.mkdir(parents=True, exist_ok=True)

    apoe = load_apoe()
    loaded = {}
    for name in args.configs:
        c = load_config(name, args.split)
        if c is None:
            print(f"  [skip] {name} (no predictions)"); continue
        loaded[name] = c
    if not loaded:
        print("No configs loaded — nothing to plot."); return
    ref = next(iter(loaded.values()))
    grid = ref["grid"]
    # KM time/event are patient-level: take from the first config (identical set across configs)
    pid, time, event = ref["pid"], ref["time"], ref["event"]
    csv_rows = []

    # vanilla tabular Cox PH baseline (02l, 8 feat), aligned to the SAME patient set
    if not args.no_tabular_cox:
        loaded[TAB_COX] = {"S": tabular_cox_meanS(pid, grid), "grid": grid}
        print(f"  + overlaid {CONFIG_LABEL[TAB_COX]}")

    def _ls(name):
        return ":" if name == TAB_COX else "--"

    # ── Fig 1: average / marginal ────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    curves = []
    for name, c in loaded.items():
        mS = np.nanmean(c["S"], 0)
        curves.append((CONFIG_LABEL.get(name, name), CONFIG_COLOR.get(name, None), c["grid"], mS, _ls(name)))
        for tg, sv in zip(c["grid"], mS):
            csv_rows.append({"figure": "average", "stratum": "all", "model": name,
                             "t": float(tg), "marginal_S": float(sv)})
    _draw(ax, time, event, curves,
          f"Marginal survival vs Ground Truth KM — CN+MCI → AD ({args.split} set)\n"
          "model curve = mean over patients of predicted S(t); shaded = 95% CI of KM estimate")
    fig.tight_layout(); fig.savefig(out / "km_calibration_average.png", dpi=150)
    fig.savefig(out / "km_calibration_average.pdf"); plt.close(fig)

    def stratified_fig(group_of, strata, fname, figkey, suptitle):
        """group_of: pid->code ; strata: list of (code, panel_title). One panel per stratum."""
        fig, axes = plt.subplots(1, len(strata), figsize=(6.5 * len(strata), 5.4), sharey=True)
        if len(strata) == 1:
            axes = [axes]
        for ax, (code, stitle) in zip(axes, strata):
            mask = np.array([group_of[p] == code for p in pid])
            if mask.sum() == 0:
                ax.set_title(f"{stitle}\n(no patients)"); ax.set_xlim(0, TMAX); ax.set_ylim(0, 1.02)
                continue
            curves = []
            for name, c in loaded.items():
                mS = np.nanmean(c["S"][mask], 0)
                curves.append((CONFIG_LABEL.get(name, name), CONFIG_COLOR.get(name, None), c["grid"], mS, _ls(name)))
                for tg, sv in zip(c["grid"], mS):
                    csv_rows.append({"figure": figkey, "stratum": stitle, "model": name,
                                     "t": float(tg), "marginal_S": float(sv)})
            _draw(ax, time[mask], event[mask], curves, stitle)
        fig.suptitle(suptitle, fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        fig.savefig(out / f"{fname}.png", dpi=150); fig.savefig(out / f"{fname}.pdf"); plt.close(fig)

    # ── Fig 2: stratified by APOE4 CARRIER status (2-way) ─────────────────────────
    carrier = {p: (1 if apoe.get(p, np.nan) >= 1 else (0 if apoe.get(p, np.nan) == 0 else -1)) for p in pid}
    stratified_fig(carrier, [(0, "APOE4 non-carrier (0 copies)"), (1, "APOE4 carrier (≥1 copy)")],
                   "km_calibration_by_apoe4", "by_apoe4",
                   f"Survival by APOE4 carrier status — model marginal S(t) vs Ground Truth KM ({args.split} set)")

    # ── Fig 3: stratified by APOE4 DOSAGE (3-way: 0 / 1 / 2 copies) ───────────────
    dose = {p: (int(apoe[p]) if (p in apoe.index and apoe.get(p) in (0, 1, 2)) else -1) for p in pid}
    stratified_fig(dose, [(0, "APOE4 = 0 copies"), (1, "APOE4 = 1 copy"), (2, "APOE4 = 2 copies")],
                   "km_calibration_by_apoe4_dosage", "by_apoe4_dosage",
                   f"Survival by APOE4 dosage — model marginal S(t) vs Ground Truth KM ({args.split} set)")

    pd.DataFrame(csv_rows).to_csv(out / "km_calibration_curves.csv", index=False)
    n0 = sum(v == 0 for v in carrier.values()); n1 = sum(v == 1 for v in carrier.values())
    d0 = sum(v == 0 for v in dose.values()); d1 = sum(v == 1 for v in dose.values()); d2 = sum(v == 2 for v in dose.values())
    print(f"APOE4 carrier (n): non-carrier={n0}  carrier={n1}")
    print(f"APOE4 dosage  (n): 0={d0}  1={d1}  2={d2}  unknown(dropped)={sum(v == -1 for v in dose.values())}")
    print(f"Done → {out}  (average / by_apoe4 / by_apoe4_dosage .png+.pdf , km_calibration_curves.csv)")


if __name__ == "__main__":
    main()
