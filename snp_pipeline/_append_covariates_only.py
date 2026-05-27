"""
_append_covariates_only.py
==========================
Compute covariate-only baseline metrics (no PRS feature) and append them
into the existing per-task leaderboards. Lets us see by how much PRS lifts
performance vs the demographic + APOE-only baselines.

Two baseline rows are added per task:
  age+sex+apoe4
  age+sex+apoe4+apoe2

For raw-β leaderboard: append to all 4 LD-config × 3-task tables.
For PRS-CS leaderboard: append to 3-task tables (unpruned strict-QC pool).

Outputs (in-place updates):
  outputs/strict_qc_prs/ld_<cfg>/{task}/per_run_metrics.tsv      ← appended
  outputs/strict_qc_prs/ld_<cfg>/{task}/leaderboard.csv          ← re-aggregated
  outputs/strict_qc_prs_prscs/{task}/per_run_metrics.tsv          ← appended
  outputs/strict_qc_prs_prscs/{task}/leaderboard.csv              ← re-aggregated
"""
from __future__ import annotations
import math
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import (roc_auc_score, balanced_accuracy_score,
                              recall_score, precision_score, f1_score,
                              average_precision_score, confusion_matrix,
                              mean_absolute_error, mean_squared_error, r2_score)
from sklearn.pipeline import Pipeline
from scipy import stats
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _prs_strict_qc_lib import load_subject_covariates  # noqa: E402

SPLITS_BASELINE = Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                        "no_cdr_stratified_post_exclusion/tabular/baseline")
SPLITS_AAO = Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                   "no_cdr_stratified_post_exclusion_aao")
OUT_BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220/outputs/strict_qc_prs")
SEEDS = (0, 1, 2)

# Covariate sets to fit (no PRS feature in any of these)
COVAR_SETS = {
    "age+sex+apoe4":       ["age_at_baseline", "Sex_M", "APOE4_Dosage"],
    "age+sex+apoe4+apoe2": ["age_at_baseline", "Sex_M", "APOE4_Dosage", "APOE2_Dosage"],
}

LD_CONFIGS = ["ld_1000kb_r2_0.8", "ld_500kb_r2_0.2",
               "ld_250kb_r2_0.1", "ld_50_5_r2_0.5"]


# ── classification fit ────────────────────────────────────────────────────
def _cls_metrics(y_true, y_score):
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score).astype(float)
    pred = (y_score > 0.5).astype(int)
    try:    auc = roc_auc_score(y_true, y_score)
    except: auc = float("nan")
    try:    pr_auc = average_precision_score(y_true, y_score)
    except: pr_auc = float("nan")
    bacc = balanced_accuracy_score(y_true, pred)
    sens = recall_score(y_true, pred, zero_division=0)
    spec = recall_score(y_true, pred, pos_label=0, zero_division=0)
    prec = precision_score(y_true, pred, zero_division=0)
    f1   = f1_score(y_true, pred, zero_division=0)
    cm = confusion_matrix(y_true, pred, labels=[0, 1])
    tn, fp, fn, tp = (cm.ravel() if cm.size == 4 else (0, 0, 0, 0))
    return {"auc": auc, "bacc": bacc, "f1": f1, "sens": sens, "spec": spec,
             "prec": prec, "pr_auc": pr_auc,
             "tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp),
             "n": int(len(y_true)),
             "pos": int((y_true == 1).sum()), "neg": int((y_true == 0).sum())}


def _load_cls_labels(seed: int, split: str) -> pd.DataFrame:
    df = pd.read_csv(SPLITS_BASELINE / f"seed_{seed}/{split}.csv", dtype=str)
    df["Patient_ID"] = df["Patient_ID"].astype(str)
    df["y_pos"] = ((df["AD_bl"].astype(int) == 1)
                    | (df["pMCI"].astype(int) == 1)
                    | (df["CN_to_AD"].astype(int) == 1)).astype(int)
    df["y_neg"] = (df["sCN"].astype(int) == 1).astype(int)
    df = df[(df["y_pos"] == 1) | (df["y_neg"] == 1)].copy()
    df["y"] = df["y_pos"]
    return df[["Patient_ID", "y"]]


def classification_covars_only(cov: pd.DataFrame, covars: list, seed: int) -> dict | None:
    parts = {sp: _load_cls_labels(seed, sp) for sp in ("train", "val", "test")}
    fits = {sp: parts[sp].merge(cov, on="Patient_ID", how="left")
              for sp in ("train", "val", "test")}
    tr = fits["train"]
    if tr.empty or tr["y"].nunique() < 2: return None
    needed = ["y"] + covars
    tr_in = fits["train"][needed].dropna()
    val_in = fits["val"][needed].dropna()
    test_in = fits["test"][needed].dropna()
    if len(tr_in) < 20 or tr_in["y"].nunique() < 2: return None
    pipe = Pipeline([("imp", SimpleImputer(strategy="median")),
                      ("clf", LogisticRegression(max_iter=2000,
                                                   class_weight="balanced",
                                                   random_state=seed))])
    pipe.fit(tr_in[covars].astype(float), tr_in["y"].astype(int))
    val_p  = pipe.predict_proba(val_in[covars].astype(float))[:, 1]
    test_p = pipe.predict_proba(test_in[covars].astype(float))[:, 1]
    val_m  = _cls_metrics(val_in["y"], val_p)
    test_m = _cls_metrics(test_in["y"], test_p)
    out = {"source": "covariates_only",
            "covar_mode": f"covariates_only_{'+'.join(c.lower() for c in covars).replace('age_at_baseline','age').replace('sex_m','sex').replace('apoe4_dosage','apoe4').replace('apoe2_dosage','apoe2')}",
            "seed": seed, "or_per_sd": float("nan"), "n_snps_used": 0}
    for k, v in val_m.items():  out[f"val_{k}"] = v
    for k, v in test_m.items(): out[f"test_{k}"] = v
    return out


# ── cox fit ───────────────────────────────────────────────────────────────
EPS = 1e-3


def _load_cox_te(seed: int, split: str) -> pd.DataFrame:
    df = pd.read_csv(SPLITS_BASELINE / f"seed_{seed}/{split}.csv", dtype=str)
    df["Patient_ID"] = df["Patient_ID"].astype(str)
    df["AD_bl"] = pd.to_numeric(df["AD_bl"], errors="coerce")
    df["pMCI"]  = pd.to_numeric(df["pMCI"], errors="coerce")
    df["CN_to_AD"] = pd.to_numeric(df["CN_to_AD"], errors="coerce")
    df["AD_final"] = pd.to_numeric(df.get("AD_final", 0), errors="coerce").fillna(0)
    df["years_to_AD"] = pd.to_numeric(df["years_to_AD"], errors="coerce")
    df["FU_years"] = pd.to_numeric(df["FU_years"], errors="coerce")
    is_conv = ((df["AD_bl"] == 1) | (df["pMCI"] == 1) | (df["CN_to_AD"] == 1))
    t = np.where(df["AD_bl"] == 1, EPS,
                  np.where(is_conv, df["years_to_AD"], df["FU_years"]))
    e = (df["AD_final"] == 1).astype(int).values
    df["t"] = t.astype(float); df["e"] = e
    df = df[(df["t"] > 0) & df["t"].notna()].copy()
    return df[["Patient_ID","t","e"]]


def cox_covars_only(cov: pd.DataFrame, covars: list, seed: int) -> dict | None:
    parts = {sp: _load_cox_te(seed, sp) for sp in ("train","val","test")}
    fits = {sp: parts[sp].merge(cov, on="Patient_ID", how="left") for sp in ("train","val","test")}
    tr = fits["train"]
    if tr.empty or tr["e"].sum() < 5 or len(tr) < 30: return None
    needed = ["t","e"] + covars
    tr_in = fits["train"][needed].dropna()
    val_in = fits["val"][needed].dropna()
    test_in = fits["test"][needed].dropna()
    if tr_in.empty or tr_in["e"].sum() < 5: return None
    cph = CoxPHFitter(penalizer=1e-3)
    try:    cph.fit(tr_in, duration_col="t", event_col="e")
    except: return None

    def _cidx(p):
        if len(p) < 8 or p["e"].sum() < 2: return float("nan")
        pi = -cph.predict_partial_hazard(p[covars]).values
        try:    return float(concordance_index(p["t"], pi, p["e"]))
        except: return float("nan")
    return {"source": "covariates_only",
             "covar_mode": f"covariates_only_{'+'.join(c.lower() for c in covars).replace('age_at_baseline','age').replace('sex_m','sex').replace('apoe4_dosage','apoe4').replace('apoe2_dosage','apoe2')}",
             "seed": seed, "n_snps_used": 0,
             "n_train": len(tr_in), "events_train": int(tr_in["e"].sum()),
             "n_val": len(val_in), "events_val": int(val_in["e"].sum() if len(val_in) else 0),
             "n_test": len(test_in), "events_test": int(test_in["e"].sum() if len(test_in) else 0),
             "hr_per_sd": float("nan"), "hr_lo": float("nan"), "hr_hi": float("nan"), "hr_p": float("nan"),
             "val_cindex": _cidx(val_in), "test_cindex": _cidx(test_in)}


# ── AAO fit ───────────────────────────────────────────────────────────────
def _load_aao(seed: int, split: str) -> pd.DataFrame:
    df = pd.read_csv(SPLITS_AAO / f"seed_{seed}/{split}.csv", dtype=str)
    df["Patient_ID"] = df["Patient_ID"].astype(str)
    df["AAO"] = pd.to_numeric(df["AAO"], errors="coerce")
    return df[["Patient_ID","AAO"]].dropna()


def _aao_metrics(y_true, y_pred):
    if len(y_true) < 2:
        return {"mse":float("nan"),"mae":float("nan"),"rmse":float("nan"),
                 "r2":float("nan"),"pearson_r":float("nan"),"pearson_p":float("nan"),
                 "n":int(len(y_true))}
    mse = mean_squared_error(y_true, y_pred); mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred); rmse = math.sqrt(mse)
    if np.std(y_pred) > 1e-9: r, p = stats.pearsonr(y_true, y_pred)
    else: r, p = float("nan"), float("nan")
    return {"mse":float(mse),"mae":float(mae),"rmse":float(rmse),
             "r2":float(r2),"pearson_r":float(r),"pearson_p":float(p),
             "n":int(len(y_true))}


def aao_covars_only(cov: pd.DataFrame, covars: list, seed: int) -> dict | None:
    # AAO drops age_at_baseline (would be circular) — same convention as 47.
    covars_no_age = [c for c in covars if c != "age_at_baseline"]
    parts = {sp: _load_aao(seed, sp) for sp in ("train","val","test")}
    fits = {sp: parts[sp].merge(cov, on="Patient_ID", how="left") for sp in ("train","val","test")}
    tr = fits["train"]
    if len(tr) < 30 or not covars_no_age: return None
    needed = covars_no_age + ["AAO"]
    tr_in = fits["train"][needed].dropna()
    val_in = fits["val"][needed].dropna()
    test_in = fits["test"][needed].dropna()
    if len(tr_in) < 10 or len(val_in) < 3: return None
    pipe = Pipeline([("imp", SimpleImputer(strategy="median")), ("reg", LinearRegression())])
    pipe.fit(tr_in[covars_no_age].astype(float), tr_in["AAO"].astype(float))
    val_pred = pipe.predict(val_in[covars_no_age].astype(float))
    test_pred = pipe.predict(test_in[covars_no_age].astype(float))
    val_m = _aao_metrics(val_in["AAO"].values, val_pred)
    test_m = _aao_metrics(test_in["AAO"].values, test_pred)
    out = {"source": "covariates_only",
            "covar_mode": f"covariates_only_{'+'.join(c.lower() for c in covars).replace('age_at_baseline','age').replace('sex_m','sex').replace('apoe4_dosage','apoe4').replace('apoe2_dosage','apoe2')}",
            "seed": seed, "n_train": len(tr_in), "n_val": len(val_in), "n_test": len(test_in),
            "yrs_per_sd_prs": float("nan"), "n_snps_used": 0}
    for k, v in val_m.items():  out[f"val_{k}"] = v
    for k, v in test_m.items(): out[f"test_{k}"] = v
    return out


# ── aggregate (mirror what scripts 45/46/47 do) ───────────────────────────
def _agg_classification(df: pd.DataFrame) -> pd.DataFrame:
    metric_cols = ["val_auc","val_bacc","val_f1","val_sens","val_spec",
                    "val_prec","val_pr_auc","or_per_sd",
                    "test_auc","test_bacc","test_f1"]
    return df.groupby(["source","covar_mode"]).agg(
        n_snps=("n_snps_used","first"), n_seeds=("seed","count"),
        val_n=("val_n","mean"), val_pos=("val_pos","mean"), val_neg=("val_neg","mean"),
        val_tn=("val_tn","mean"), val_fp=("val_fp","mean"),
        val_fn=("val_fn","mean"), val_tp=("val_tp","mean"),
        **{f"{c}_mean": (c, "mean") for c in metric_cols},
        **{f"{c}_std":  (c, "std")  for c in metric_cols},
    ).reset_index().sort_values(["covar_mode","val_bacc_mean"], ascending=[True, False])


def _agg_cox(df: pd.DataFrame) -> pd.DataFrame:
    return df.groupby(["source","covar_mode"]).agg(
        n_snps=("n_snps_used","first"), n_seeds=("seed","count"),
        val_cindex_mean=("val_cindex","mean"), val_cindex_std=("val_cindex","std"),
        hr_per_sd_mean=("hr_per_sd","mean"), hr_per_sd_std=("hr_per_sd","std"),
        hr_lo_mean=("hr_lo","mean"), hr_hi_mean=("hr_hi","mean"),
        hr_p_min=("hr_p","min"),
        test_cindex_mean=("test_cindex","mean"),
    ).reset_index().sort_values(["covar_mode","val_cindex_mean"], ascending=[True, False])


def _agg_aao(df: pd.DataFrame) -> pd.DataFrame:
    return df.groupby(["source","covar_mode"]).agg(
        n_snps=("n_snps_used","first"), n_seeds=("seed","count"),
        val_r2_mean=("val_r2","mean"), val_r2_std=("val_r2","std"),
        val_rmse_mean=("val_rmse","mean"),
        val_pearson_r_mean=("val_pearson_r","mean"),
        val_pearson_p_min=("val_pearson_p","min"),
        yrs_per_sd_mean=("yrs_per_sd_prs","mean"), yrs_per_sd_std=("yrs_per_sd_prs","std"),
    ).reset_index().sort_values(["covar_mode","val_r2_mean"], ascending=[True, False])


# ── driver ────────────────────────────────────────────────────────────────
def update_leaderboard(task_dir: Path, fit_fn, agg_fn, cov: pd.DataFrame):
    per_run_p = task_dir / "per_run_metrics.tsv"
    lb_p      = task_dir / "leaderboard.csv"
    if not per_run_p.exists() or not lb_p.exists():
        print(f"  [skip] {task_dir} (missing per_run / leaderboard)")
        return
    existing = pd.read_csv(per_run_p, sep="\t")
    # Drop any existing covariates_only rows so we replace
    existing = existing[existing["source"] != "covariates_only"].copy()

    new_rows = []
    for covar_name, covars in COVAR_SETS.items():
        for seed in SEEDS:
            r = fit_fn(cov, covars, seed)
            if r is None: continue
            new_rows.append(r)
    if not new_rows:
        print(f"  [no rows] {task_dir}"); return
    combined = pd.concat([existing, pd.DataFrame(new_rows)], ignore_index=True)
    combined.to_csv(per_run_p, sep="\t", index=False)

    agg = agg_fn(combined)
    agg.to_csv(lb_p, index=False)
    print(f"  OK {task_dir.relative_to(OUT_BASE.parent)}: +{len(new_rows)} per-run rows, leaderboard agg -> {len(agg)}")


def main():
    cov = load_subject_covariates()

    # raw: 4 LD configs × 3 tasks
    for cfg in LD_CONFIGS:
        for task, fit_fn, agg_fn in [
            ("classification", classification_covars_only, _agg_classification),
            ("cox",            cox_covars_only,            _agg_cox),
            ("aao_regression", aao_covars_only,            _agg_aao),
        ]:
            update_leaderboard(OUT_BASE / cfg / task, fit_fn, agg_fn, cov)

    # prscs: 1 unpruned dir × 3 tasks
    prscs_base = OUT_BASE.parent / "strict_qc_prs_prscs"
    for task, fit_fn, agg_fn in [
        ("classification", classification_covars_only, _agg_classification),
        ("cox",            cox_covars_only,            _agg_cox),
        ("aao_regression", aao_covars_only,            _agg_aao),
    ]:
        update_leaderboard(prscs_base / task, fit_fn, agg_fn, cov)


if __name__ == "__main__":
    main()
