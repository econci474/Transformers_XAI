r"""
02_score_xsurv.py  (env: survml — scikit-survival)
==================================================
Score every cross-modal survival run exported by 01_train_xsurv.py with the SAME 02l.evaluate() used for
the T4 baselines: risk = 1 − S(t_ref=4y); C-index, IPCW-C, IBS, AUC@{3,5,7,10}y; per split val/test/valtest.
A run's surv_fn linearly interpolates its exported dense S-grid, identical to T4 02_survival_comparison.py.

Run:  conda run -n survml python integration/cross_modal_survival/02_score_xsurv.py
Out:  outputs/xsurv/master_xsurv_comparison.csv (+ _summary.csv)
"""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sksurv.util import Surv

HERE = Path(__file__).resolve().parent
PRED_ROOT = HERE / "outputs" / "xsurv"
CLIN = HERE.parent.parent / "clinical_pipeline"


def _load_02l():
    spec = importlib.util.spec_from_file_location("l02", CLIN / "02l_survival_model_comparison.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


L = _load_02l()


def interp_surv_fn(Sgrid, grid_times):
    def surv_fn(X, times):
        times = np.asarray(times, float)
        out = np.empty((Sgrid.shape[0], times.size))
        for i in range(Sgrid.shape[0]):
            out[i] = np.interp(times, grid_times, Sgrid[i])
        return np.clip(out, 1e-8, 1.0)
    return surv_fn


def score_run(variant, arm):
    """All seeds × splits for one (variant, arm) → metric rows (schema = 02l.evaluate)."""
    rows = []
    for seed in (0, 1, 2):
        d = PRED_ROOT / variant / arm / f"seed_{seed}"
        if not (d / "predictions.parquet").exists():
            return None
        gt = np.asarray(json.load(open(d / "meta.json"))["grid_times"], float)
        df = pd.read_parquet(d / "predictions.parquet")
        scol = [c for c in df.columns if c.startswith("S_")]
        parts = {sp: df[df.split == sp] for sp in ("train", "val", "test")}
        parts["valtest"] = df[df.split.isin(["val", "test"])]
        tr = parts["train"]
        etr = tr["event"].to_numpy() > 0.5; ttr = tr["time"].to_numpy()
        ytr = Surv.from_arrays(event=etr, time=ttr)
        for sp in ("val", "test", "valtest"):
            g = parts[sp]
            e = g["event"].to_numpy() > 0.5; t = g["time"].to_numpy()
            sf = interp_surv_fn(g[scol].to_numpy(float), gt)
            Xz = np.zeros((len(g), 1))
            rows.append({"variant": variant, "arm": arm, "model": f"{variant}__{arm}",
                         "seed": seed, "split": sp, **L.evaluate(sf, ytr, ttr, etr, Xz, e, t)})
    return rows


def main():
    variants = sorted(p.name for p in PRED_ROOT.iterdir() if p.is_dir() and p.name != "comparison")
    rows = []
    for variant in variants:
        for armdir in sorted((PRED_ROOT / variant).iterdir()):
            if not armdir.is_dir():
                continue
            r = score_run(variant, armdir.name)
            if r:
                rows += r
            else:
                print(f"  [skip] {variant}/{armdir.name} (no predictions)")
    met = pd.DataFrame(rows)
    met.to_csv(PRED_ROOT / "master_xsurv_comparison.csv", index=False)
    summ = met.groupby(["model", "split"])[L.MC].agg(["mean", "std"]).round(4)
    summ.to_csv(PRED_ROOT / "master_xsurv_comparison_summary.csv")
    print(summ.to_string())
    print(f"\nWrote {PRED_ROOT / 'master_xsurv_comparison.csv'} ({len(met)} rows)")


if __name__ == "__main__":
    main()
