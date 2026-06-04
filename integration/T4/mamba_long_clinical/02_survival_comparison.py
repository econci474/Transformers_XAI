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


def derived_horizon(name):
    """3-class T4 horizon from S(3)/S(7) on converters (test, pooled across seeds)."""
    yy, pp = [], []
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
        pp.append(prob.argmax(1)); yy.append(g["Label_T4"].to_numpy(int))
    if not yy:
        return None
    y = np.concatenate(yy); p = np.concatenate(pp)
    return {"model": name, "n": int(len(y)),
            "bacc": round(float(balanced_accuracy_score(y, p)), 4),
            "macro_f1": round(float(f1_score(y, p, average="macro", zero_division=0)), 4)}


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

    dh = [r for r in (derived_horizon(n) for n in configs) if r]
    if dh:
        pd.DataFrame(dh).to_csv(OUT / "derived_horizon_from_survival.csv", index=False)
        print("\nDerived 3-class T4 horizon from S(3)/S(7) on converters (test pooled):")
        print(pd.DataFrame(dh).to_string(index=False))
    print(f"\nDone → {OUT}")


if __name__ == "__main__":
    main()
