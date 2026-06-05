r"""
02_survival_comparison.py  (env: survml — scikit-survival + lifelines)
======================================================================
Score every MAMBA survival config (exported by 01_train_survival_mamba.py) with the SAME
02l.evaluate() used for the tabular baselines, then table them together against 02l's Cox PH /
Weibull-AFT / RSF (8 tabular features) on the identical CN+MCI canonical split. Plus a derived
3-class T4-horizon evaluation on converters from each config's S(3)/S(7).

We reuse the MAMBA per-patient survival curves analytically: a config's surv_fn(X, times) linearly
interpolates that patient's exported dense S-grid to whatever times 02l.evaluate requests — so the
MAMBA heads are scored on identical C-index / IPCW C / IBS / time-AUC as the baselines.

Run:  conda run -n survml python integration/T4/mamba_long_clinical/02_survival_comparison.py
Out:  outputs/survival/comparison/{master_survival_comparison.csv, survival_comparison_table.{png,pdf},
      derived_horizon_from_survival.csv}
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score, f1_score
from sksurv.util import Surv

HERE = Path(__file__).resolve().parent
PRED_ROOT = HERE / "outputs" / "survival"
OUT = PRED_ROOT / "comparison"
OUT.mkdir(parents=True, exist_ok=True)
SEEDS = (0, 1, 2)
CLIN = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\clinical_pipeline")


def _load_02l():
    spec = importlib.util.spec_from_file_location("l02", CLIN / "02l_survival_model_comparison.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


L = _load_02l()
BASELINES = ["RSF", "Cox PH", "Weibull AFT"]


def interp_surv_fn(Sgrid, grid_times):
    """Return surv_fn(X, times) -> [n,len(times)] by per-row linear interp of the dense S grid."""
    def surv_fn(X, times):
        times = np.asarray(times, float)
        out = np.empty((Sgrid.shape[0], times.size))
        for i in range(Sgrid.shape[0]):
            out[i] = np.interp(times, grid_times, Sgrid[i])
        return np.clip(out, 1e-8, 1.0)
    return surv_fn


def eval_mamba_config(name):
    """All seeds × splits for one MAMBA config -> list of metric rows (schema = 02l.evaluate)."""
    rows = []
    for seed in SEEDS:
        d = PRED_ROOT / name / f"seed_{seed}"
        if not (d / "predictions.parquet").exists():
            return None                                          # config not run (e.g. soden skipped)
        import json
        gt = np.asarray(json.load(open(d / "meta.json"))["grid_times"], float)
        df = pd.read_parquet(d / "predictions.parquet")
        scol = [c for c in df.columns if c.startswith("S_")]
        parts = {sp: df[df.split == sp] for sp in ("train", "val", "test")}
        parts["valtest"] = df[df.split.isin(["val", "test"])]
        tr = parts["train"]
        etr = (tr["event"].to_numpy() > 0.5); ttr = tr["time"].to_numpy()
        ytr = Surv.from_arrays(event=etr, time=ttr)
        for sp in ("val", "test", "valtest"):
            g = parts[sp]
            e = (g["event"].to_numpy() > 0.5); t = g["time"].to_numpy()
            sf = interp_surv_fn(g[scol].to_numpy(float), gt)
            X = np.zeros((len(g), 1))
            rows.append({"model": name, "seed": seed, "split": sp,
                         **L.evaluate(sf, ytr, ttr, etr, X, e, t)})
    return rows


def eval_baselines():
    """02l tabular Cox/Weibull-AFT/RSF on the same canonical CN+MCI split."""
    specs = L.make_model_specs(BASELINES)
    rows = []
    for seed in SEEDS:
        dd = L.load_seed(seed, L.COMPACT_FEATURES)
        Xtr, etr, ttr = dd["train"]; ytr = Surv.from_arrays(event=etr, time=ttr)
        dd["valtest"] = (pd.concat([dd["val"][0], dd["test"][0]]),
                         np.concatenate([dd["val"][1], dd["test"][1]]),
                         np.concatenate([dd["val"][2], dd["test"][2]]))
        for name in BASELINES:
            sf = specs[name](Xtr, etr, ttr)
            for sp in ("val", "test", "valtest"):
                X, e, t = dd[sp]
                rows.append({"model": name, "seed": seed, "split": sp,
                             **L.evaluate(sf, ytr, ttr, etr, X, e, t)})
    return rows


def _median_time(S_row, gt):
    """Predicted conversion time = t where S(t)=0.5 (linear interp). inf if S never <0.5 on the grid.
    For the Weibull-piecewise head this equals its closed-form median (solve H0(t)=ln2)."""
    below = np.where(S_row <= 0.5)[0]
    if below.size == 0:
        return np.inf
    j = below[0]
    if j == 0:
        return float(gt[0])
    t0, t1, s0, s1 = gt[j - 1], gt[j], S_row[j - 1], S_row[j]
    return float(t1 if s1 == s0 else t0 + (0.5 - s0) * (t1 - t0) / (s1 - s0))


def _bin_time(tm):
    return 0 if tm < 3.0 else (1 if tm < 7.0 else 2)


def derived_horizon(name):
    """3-class T4 horizon on converters (test, pooled across seeds), by TWO mappings from the curve:
       • prob   : argmax of [1-S(3), S(3)-S(7), S(7)]  (uniform across all heads)
       • median : bin the predicted median conversion time (S(t)=0.5; closed form for Weibull).
    Returns one row per (model, method)."""
    yy, p_prob, p_med = [], [], []
    for seed in SEEDS:
        d = PRED_ROOT / name / f"seed_{seed}"
        if not (d / "predictions.parquet").exists():
            return None
        import json
        gt = np.asarray(json.load(open(d / "meta.json"))["grid_times"], float)
        df = pd.read_parquet(d / "predictions.parquet")
        g = df[(df.split == "test") & df["Label_T4"].notna()]
        if len(g) == 0:
            continue
        scol = [c for c in df.columns if c.startswith("S_")]
        S = g[scol].to_numpy(float)
        S3 = np.array([np.interp(3.0, gt, r) for r in S])
        S7 = np.array([np.interp(7.0, gt, r) for r in S])
        prob = np.clip(np.stack([1 - S3, S3 - S7, S7], 1), 0, None)
        prob = prob / prob.sum(1, keepdims=True)
        p_prob.append(prob.argmax(1))
        p_med.append(np.array([_bin_time(_median_time(r, gt)) for r in S]))
        yy.append(g["Label_T4"].to_numpy(int))
    if not yy:
        return None
    y = np.concatenate(yy)
    out = []
    for method, p in (("prob", np.concatenate(p_prob)), ("median_time", np.concatenate(p_med))):
        out.append({"model": name, "method": method, "n": int(len(y)),
                    "bacc": round(float(balanced_accuracy_score(y, p)), 4),
                    "macro_f1": round(float(f1_score(y, p, average="macro", zero_division=0)), 4)})
    return out


# ── APOE-stratified discrimination/calibration (honest per-subgroup reporting) ────
def apoe_dose_map():
    """Patient_ID -> APOE4 dosage (0/1/2) from the canonical seed_0 splits."""
    frames = []
    for sp in ("train", "val", "test"):
        d = pd.read_csv(L.SPLIT_DIR / "seed_0" / f"{sp}.csv", low_memory=False,
                        usecols=lambda c: c in ("Patient_ID", "APOE4_Dosage"))
        frames.append(d)
    a = pd.concat(frames, ignore_index=True)
    a["Patient_ID"] = a["Patient_ID"].astype(str)
    a["APOE4_Dosage"] = pd.to_numeric(a["APOE4_Dosage"], errors="coerce")
    a = a.dropna(subset=["APOE4_Dosage"]).drop_duplicates("Patient_ID")
    return {r.Patient_ID: int(r.APOE4_Dosage) for r in a.itertuples() if r.APOE4_Dosage in (0, 1, 2)}


def _pooled_mamba(name):
    """Pool the held-out TEST curves of a MAMBA config across the 3 seeds (+ pooled train reference)."""
    import json
    S, e, t, pid, trE, trT, grid = [], [], [], [], [], [], None
    for seed in SEEDS:
        d = PRED_ROOT / name / f"seed_{seed}"
        if not (d / "predictions.parquet").exists():
            return None
        grid = np.asarray(json.load(open(d / "meta.json"))["grid_times"], float)
        df = pd.read_parquet(d / "predictions.parquet"); df["Patient_ID"] = df["Patient_ID"].astype(str)
        sc = [c for c in df.columns if c.startswith("S_")]
        te, tr = df[df.split == "test"], df[df.split == "train"]
        S.append(te[sc].to_numpy(float)); e.append(te["event"].to_numpy() > 0.5); t.append(te["time"].to_numpy())
        pid.append(te["Patient_ID"].to_numpy()); trE.append(tr["event"].to_numpy() > 0.5); trT.append(tr["time"].to_numpy())
    return dict(S=np.vstack(S), e=np.concatenate(e), t=np.concatenate(t), pid=np.concatenate(pid),
                trE=np.concatenate(trE), trT=np.concatenate(trT), grid=grid)


def _pooled_baseline(name, grid):
    """Pool the held-out TEST curves of a tabular 02l baseline across the 3 seeds."""
    spec = L.make_model_specs([name])[name]
    S, e, t, pid, trE, trT = [], [], [], [], [], []
    for seed in SEEDS:
        dd = L.load_seed(seed, L.COMPACT_FEATURES)
        Xtr, etr, ttr = dd["train"]; sf = spec(Xtr, etr, ttr)
        X, ee, tt = dd["test"]; Sg = np.asarray(sf(X, grid), float)
        csv = pd.read_csv(L.SPLIT_DIR / f"seed_{seed}" / "test.csv", low_memory=False)
        csv = csv[csv["bl_dx"].isin(["CN", "MCI"])].copy(); csv["Patient_ID"] = csv["Patient_ID"].astype(str)
        pids = csv.loc[X.index, "Patient_ID"].to_numpy()
        S.append(Sg); e.append(ee > 0.5); t.append(tt); pid.append(pids)
        trE.append(etr > 0.5); trT.append(ttr)
    return dict(S=np.vstack(S), e=np.concatenate(e), t=np.concatenate(t), pid=np.concatenate(pid),
                trE=np.concatenate(trE), trT=np.concatenate(trT), grid=grid)


def stratified_rows(name, pooled, dose):
    """C-index (Harrell + IPCW) and IBS per APOE dosage stratum, on TEST pooled across seeds."""
    grid = pooled["grid"]; ytr = Surv.from_arrays(event=pooled["trE"], time=pooled["trT"])
    out = []
    for code, lab in [(-1, "all"), (0, "APOE_0"), (1, "APOE_1"), (2, "APOE_2")]:
        mask = (np.ones(len(pooled["pid"]), bool) if code == -1
                else np.array([dose.get(p, -9) == code for p in pooled["pid"]]))
        ev, tt = pooled["e"][mask], pooled["t"][mask]
        n, nev = int(mask.sum()), int(ev.sum())
        row = {"model": name, "stratum": lab, "n": n, "events": nev}
        if nev >= 3 and n >= 8 and int((~ev).sum()) >= 1:
            try:
                m = L.evaluate(interp_surv_fn(pooled["S"][mask], grid), ytr, pooled["trT"], pooled["trE"],
                               np.zeros((n, 1)), ev, tt)
                row.update({"c_index": round(m["c_index"], 4), "c_index_ipcw": round(m["c_index_ipcw"], 4),
                            "ibs": round(m["ibs"], 4)})
            except Exception as ex:
                row.update({"c_index": np.nan, "c_index_ipcw": np.nan, "ibs": np.nan, "note": str(ex)[:50]})
        else:
            row.update({"c_index": np.nan, "c_index_ipcw": np.nan, "ibs": np.nan, "note": "too few events"})
        out.append(row)
    return out


def main():
    global PRED_ROOT, OUT
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred_root", type=str, default=str(PRED_ROOT),
                    help="dir of 01_train predictions (<config>/seed_*/predictions.parquet); pass the "
                         "hpc-work path you pulled, so nothing is read from / written to the repo folder")
    args = ap.parse_args()
    PRED_ROOT = Path(args.pred_root); OUT = PRED_ROOT / "comparison"; OUT.mkdir(parents=True, exist_ok=True)

    configs = sorted(p.name for p in PRED_ROOT.iterdir()
                     if p.is_dir() and p.name != "comparison")
    print(f"MAMBA configs found: {configs}")
    rows = []
    for name in configs:
        r = eval_mamba_config(name)
        if r:
            rows += r
        else:
            print(f"  [skip] {name} (no predictions)")
    rows += eval_baselines()
    met = pd.DataFrame(rows)
    met.to_csv(OUT / "master_survival_comparison.csv", index=False)
    summ = met.groupby(["model", "split"])[L.MC].agg(["mean", "std"]).round(4)
    summ.to_csv(OUT / "master_survival_comparison_summary.csv")
    print(summ.to_string())

    models = [m for m in met["model"].unique() if m not in BASELINES] + BASELINES
    L.render_table(met, models, "model",
                   "T4 longitudinal-MAMBA survival heads vs tabular baselines (CN+MCI → AD; canonical split)\n"
                   "risk = 1 − S(t_ref); MAMBA = deep heads on bl/m06/m12 clinical embeddings; "
                   "baselines = 8 tabular features. mean ± std over 3 seeds; bold = best per split.",
                   "survival_comparison_table", OUT, "MAMBA S(t) interpolated from a 61-pt grid; "
                   "baselines via 02l.fit_model. Same CN+MCI canonical folds for all.")

    dh = [r for rows in (derived_horizon(n) for n in configs) if rows for r in rows]
    if dh:
        dhdf = pd.DataFrame(dh).sort_values(["method", "bacc"], ascending=[True, False])
        dhdf.to_csv(OUT / "derived_horizon_from_survival.csv", index=False)
        print("\nDerived 3-class T4 horizon on converters (test pooled) — from the SURVIVAL curve:")
        print("  prob = argmax[1-S(3),S(3)-S(7),S(7)] ; median_time = bin of predicted median (S=0.5)")
        print(dhdf.to_string(index=False))

    # ── APOE-dosage-stratified discrimination + calibration (honest subgroup table) ──
    dose = apoe_dose_map()
    GRID = np.round(np.linspace(0.0, 15.0, 61), 4)
    srows = []
    for name in configs:
        pm = _pooled_mamba(name)
        if pm is not None:
            srows += stratified_rows(name, pm, dose)
    for name in BASELINES:
        srows += stratified_rows(name, _pooled_baseline(name, GRID), dose)
    sdf = pd.DataFrame(srows)
    sdf.to_csv(OUT / "stratified_apoe_metrics.csv", index=False)
    print("\nAPOE-dosage-stratified C-index + IBS (held-out TEST pooled across 3 seeds):")
    print("  NB APOE_2 is small (~10-12 pooled) — read with the n/events columns.")
    print(sdf.to_string(index=False))
    print(f"\nDone → {OUT}")


if __name__ == "__main__":
    main()
