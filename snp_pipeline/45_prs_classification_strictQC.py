"""
45_prs_classification_strictQC.py
=================================
AAO Phase 2c.A — Per-source PRS classification on the LD-pruned strict-QC pool.

Per (source x covar_mode x seed):
  - PRS_p = Sigma beta_A1 * dosage_p (over LD-pruned strict-QC SNPs).
  - z-score on train fold, fit LogisticRegression, eval on val + test.
  - Covariates loaded via _prs_strict_qc_lib.load_subject_covariates().

Sources:
  16 individual sources (Bellenguez, Wightman, ... Kosteridis_MTAG_AD)
  + "prs_all_dedup" — for each rsID, pick beta from the source with the largest
    published GWAS N (per N_GWAS dict in the lib).

Covariate modes:
  prs_only                              (no covariates)
  prs+age+sex+apoe4                     (full clinical adjustment, current pipeline)
  prs+age+sex+apoe4+apoe2               (adds APOE2 protective allele)

Outputs (under outputs/strict_qc_prs/<ld_config>/classification/):
  per_run_metrics.tsv      — full metric set including confusion-matrix counts
  leaderboard.csv          — aggregate (mean +/- std) full metric set
  leaderboard.png          — concise per-cell view (n_snps, val_bacc, val_auc)
"""
from __future__ import annotations
import argparse
import math
from scipy.stats import norm
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (roc_auc_score, average_precision_score,
                              balanced_accuracy_score, confusion_matrix,
                              precision_score, recall_score, f1_score)
from sklearn.pipeline import Pipeline

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _prs_strict_qc_lib import (ALL_SOURCES, PRSCS_SOURCES, per_source_prs_table,
                                  load_subject_covariates,
                                  get_dedup_dosage_matrix,
                                  render_en_weights_png,
                                  META_FILTERED_SOURCES,
                                  DEFAULT_LD_CONFIG)  # noqa: E402
from _per_patient_io import (capture_predictions, save_per_patient_parquet,
                              concat_and_save_anchor)  # noqa: E402

# Per-seed EN coefficient capture for en_weights.png.
EN_COEFS = {"dosage": {}, "meta": {}, "meta_filtered": {}}  # flavor -> seed -> pd.Series

# Per-patient prediction capture for multimodal alignment.
# Only the 4 fusion-anchor sources get captured (matches plan 2026-06-02).
PER_PATIENT_SOURCES = {
    "meta_prs_EN_combined",
    "Kosteridis",
    "Kosteridis_MTAG_AD",
    "Kosteridis_shared_AD_CV",
}
PER_PATIENT_ROWS: list[pd.DataFrame] = []   # appended inside _run_one
# Single-element mutable containers so main() can publish run-context strings
# to _run_one without changing _run_one's call signature.
PER_PATIENT_BETA_SOURCE: list[str] = ["raw"]
PER_PATIENT_LD_CONFIG:   list[str] = ["na"]

SPLITS_ROOT = Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                    "no_cdr_stratified_post_exclusion/tabular/baseline")
OUT_BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220/outputs/strict_qc_prs")
SEEDS = (0, 1, 2)
COVAR_MODES = {
    "prs_only":                   [],
    "prs+age+sex+apoe4":          ["age_at_baseline", "Sex_M", "APOE4_Dosage"],
    "prs+age+sex+apoe4+apoe2":    ["age_at_baseline", "Sex_M", "APOE4_Dosage", "APOE2_Dosage"],
}


def _load_labels(seed: int, split: str) -> pd.DataFrame:
    df = pd.read_csv(SPLITS_ROOT / f"seed_{seed}/{split}.csv", dtype=str)
    df["Patient_ID"] = df["Patient_ID"].astype(str)
    df["y_pos"] = ((df["AD_bl"].astype(int) == 1)
                    | (df["pMCI"].astype(int) == 1)
                    | (df["CN_to_AD"].astype(int) == 1)).astype(int)
    df["y_neg"] = (df["sCN"].astype(int) == 1).astype(int)
    df = df[(df["y_pos"] == 1) | (df["y_neg"] == 1)].copy()
    df["y"] = df["y_pos"]
    return df[["Patient_ID", "y"]]


def _metrics(y_true, y_score) -> dict:
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score).astype(float)
    pred = (y_score > 0.5).astype(int)
    try:    auc = roc_auc_score(y_true, y_score)
    except: auc = float("nan")
    try:    pr_auc = average_precision_score(y_true, y_score)
    except: pr_auc = float("nan")
    bacc = balanced_accuracy_score(y_true, pred)
    sens = recall_score(y_true, pred, zero_division=0)        # TPR
    spec = recall_score(y_true, pred, pos_label=0, zero_division=0)  # TNR
    prec = precision_score(y_true, pred, zero_division=0)
    f1   = f1_score(y_true, pred, zero_division=0)
    cm = confusion_matrix(y_true, pred, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel() if cm.size == 4 else (0, 0, 0, 0)
    return {"auc": auc, "bacc": bacc, "f1": f1, "sens": sens, "spec": spec,
             "prec": prec, "pr_auc": pr_auc,
             "tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp),
             "n": int(len(y_true)),
             "pos": int((y_true == 1).sum()), "neg": int((y_true == 0).sum())}


def _fit_en_dosage_get_prs(dos_feat: pd.DataFrame, parts: dict,
                              cov: pd.DataFrame, seed: int):
    """Train ElasticNet LogisticRegression on TRAIN-fold dosage to predict
    AD. Return (model, log-odds-series). model=None when fit skipped/failed.

    Same EN model applied to val/test ⇒ no leakage."""
    train_pids = parts["train"]["Patient_ID"].tolist()
    train_labels = parts["train"].set_index("Patient_ID")["y"].astype(int)
    common = [p for p in train_pids if p in dos_feat.index]
    if len(common) < 30 or train_labels.loc[common].nunique() < 2:
        return None, pd.Series(np.nan, index=dos_feat.index)
    X_tr = dos_feat.loc[common].astype(float).fillna(dos_feat.mean()).values
    y_tr = train_labels.loc[common].values
    try:
        clf = LogisticRegression(penalty="elasticnet", solver="saga",
                                  l1_ratio=0.5, C=1.0, max_iter=5000,
                                  class_weight="balanced", random_state=seed)
        clf.fit(X_tr, y_tr)
        X_all = dos_feat.fillna(dos_feat.mean()).astype(float).values
        log_odds = clf.decision_function(X_all)
    except Exception:
        return None, pd.Series(np.nan, index=dos_feat.index)
    return clf, pd.Series(log_odds, index=dos_feat.index)


def _fit_en_meta_prs(source_prs: pd.DataFrame, parts: dict,
                       seed: int):
    """source_prs: PTID-indexed DataFrame with one column per source ("PRS_<src>").
    Train ElasticNet LR on TRAIN-fold source-PRS values to predict AD.
    Return (model, log-odds-series). model=None when fit skipped/failed."""
    train_pids = parts["train"]["Patient_ID"].tolist()
    train_labels = parts["train"].set_index("Patient_ID")["y"].astype(int)
    common = [p for p in train_pids if p in source_prs.index]
    if len(common) < 30 or train_labels.loc[common].nunique() < 2:
        return None, pd.Series(np.nan, index=source_prs.index)
    # Drop columns that are all-NaN
    X = source_prs.copy()
    keep_cols = [c for c in X.columns if X[c].notna().any()]
    X = X[keep_cols]
    X = X.fillna(X.mean())
    X_tr = X.loc[common].astype(float).values
    y_tr = train_labels.loc[common].values
    try:
        clf = LogisticRegression(penalty="elasticnet", solver="saga",
                                  l1_ratio=0.5, C=1.0, max_iter=5000,
                                  class_weight="balanced", random_state=seed)
        clf.fit(X_tr, y_tr)
        log_odds = clf.decision_function(X.astype(float).values)
    except Exception:
        return None, pd.Series(np.nan, index=source_prs.index)
    # Return both the model and the predictions; expose feature names so the
    # caller can build a labeled coef Series.
    clf._feature_names = keep_cols   # attached for downstream coef labeling
    return clf, pd.Series(log_odds, index=source_prs.index)


def _run_covariates_only(seed: int, mode: str, covars: list,
                           cov: pd.DataFrame) -> dict | None:
    """Clinical baseline: fit logistic(y ~ covariates) with NO PRS.

    All OR/PRS-related fields are NaN. Returns None when covars is empty
    (prs_only mode has no predictors here)."""
    if not covars: return None
    parts = {sp: _load_labels(seed, sp) for sp in ("train", "val", "test")}
    fits = {sp: parts[sp].merge(cov, on="Patient_ID", how="left")
             for sp in ("train", "val", "test")}
    tr_in   = fits["train"][["y"] + covars].dropna()
    val_in  = fits["val"][  ["y"] + covars].dropna()
    test_in = fits["test"][ ["y"] + covars].dropna()
    if len(tr_in) < 20 or tr_in["y"].nunique() < 2: return None
    pipe = Pipeline([("imp", SimpleImputer(strategy="median")),
                      ("clf", LogisticRegression(max_iter=2000,
                                                   class_weight="balanced",
                                                   random_state=seed))])
    pipe.fit(tr_in[covars].astype(float), tr_in["y"].astype(int))
    val_p  = (pipe.predict_proba(val_in[covars].astype(float))[:, 1]
                if len(val_in) else np.array([]))
    test_p = (pipe.predict_proba(test_in[covars].astype(float))[:, 1]
                if len(test_in) else np.array([]))
    val_m  = (_metrics(val_in["y"], val_p) if len(val_in)
                else _metrics(np.zeros(0), np.zeros(0)))
    test_m = (_metrics(test_in["y"], test_p) if len(test_in)
                else _metrics(np.zeros(0), np.zeros(0)))
    out = {"source": "covariates_only", "covar_mode": mode, "seed": seed,
            "or_per_sd":         float("nan"),
            "or_per_sd_lo":      float("nan"),
            "or_per_sd_hi":      float("nan"),
            "or_per_sd_p":       float("nan"),
            "or_per_sd_val":     float("nan"),
            "or_per_sd_val_lo":  float("nan"),
            "or_per_sd_val_hi":  float("nan"),
            "or_per_sd_val_p":   float("nan"),
            "or_per_sd_test":    float("nan"),
            "or_per_sd_test_lo": float("nan"),
            "or_per_sd_test_hi": float("nan"),
            "or_per_sd_test_p":  float("nan")}
    for k, v in val_m.items():  out[f"val_{k}"]  = v
    for k, v in test_m.items(): out[f"test_{k}"] = v
    return out


def _run_one(src: str, prs_full: pd.DataFrame, cov: pd.DataFrame,
              seed: int, mode: str, covars: list,
              dedup_dosage: pd.DataFrame | None = None) -> dict | None:
    if src == "covariates_only":
        return _run_covariates_only(seed, mode, covars, cov)
    parts = {sp: _load_labels(seed, sp) for sp in ("train", "val", "test")}
    # Build the "PRS" column per source. For non-EN sources, use the precomputed
    # PRS table. For EN sources, fit per-seed on train fold and produce log-odds
    # as the "PRS".
    if src == "prs_all_dedup_EN_dosage":
        if dedup_dosage is None: return None
        clf_en, en_prs = _fit_en_dosage_get_prs(dedup_dosage, parts, cov, seed)
        if en_prs.isna().all(): return None
        if clf_en is not None:
            EN_COEFS["dosage"][seed] = pd.Series(
                clf_en.coef_[0], index=list(dedup_dosage.columns))
        base = pd.DataFrame({"Patient_ID": en_prs.index.astype(str),
                              "PRS": en_prs.values})
    elif src == "meta_prs_EN_combined":
        # Drop the dedup columns (composed) and the Kosteridis SUBGROUPS
        # (MTAG_AD / shared_AD_CV / novel_AD): Kosteridis (full) = MTAG ∪ shared,
        # so including subgroups + full triple-counts the same Kosteridis SNPs.
        source_cols = [c for c in prs_full.columns if c.startswith("PRS_")
                        and c not in ("PRS_prs_all_dedup",
                                        "PRS_prs_all_dedup_ivw",
                                        "PRS_prs_all_dedup_filtered",
                                        "PRS_Kosteridis_MTAG_AD",
                                        "PRS_Kosteridis_shared_AD_CV",
                                        "PRS_Kosteridis_novel_AD")
                        and not c.endswith("EN_dosage")
                        and not c.endswith("EN_combined")]
        sp_df = prs_full[["PTID"] + source_cols].set_index("PTID")
        sp_df.index = sp_df.index.astype(str)
        clf_en, en_prs = _fit_en_meta_prs(sp_df, parts, seed)
        if en_prs.isna().all(): return None
        if clf_en is not None:
            feats = getattr(clf_en, "_feature_names", source_cols)
            EN_COEFS["meta"][seed] = pd.Series(
                clf_en.coef_[0],
                index=[c.replace("PRS_", "", 1) for c in feats])
        base = pd.DataFrame({"Patient_ID": en_prs.index,
                              "PRS": en_prs.values})
    elif src == "meta_prs_EN_filtered":
        # Same shape as meta_prs_EN_combined but restricted to the 7 curated
        # high-N European AD GWAS sources.
        source_cols = [f"PRS_{s}" for s in META_FILTERED_SOURCES
                        if f"PRS_{s}" in prs_full.columns
                        and prs_full[f"PRS_{s}"].notna().any()]
        if not source_cols: return None
        sp_df = prs_full[["PTID"] + source_cols].set_index("PTID")
        sp_df.index = sp_df.index.astype(str)
        clf_en, en_prs = _fit_en_meta_prs(sp_df, parts, seed)
        if en_prs.isna().all(): return None
        if clf_en is not None:
            feats = getattr(clf_en, "_feature_names", source_cols)
            EN_COEFS["meta_filtered"][seed] = pd.Series(
                clf_en.coef_[0],
                index=[c.replace("PRS_", "", 1) for c in feats])
        base = pd.DataFrame({"Patient_ID": en_prs.index,
                              "PRS": en_prs.values})
    else:
        prs_col = f"PRS_{src}"
        if prs_col not in prs_full.columns or prs_full[prs_col].isna().all():
            return None
        base = prs_full[["PTID", prs_col]].rename(columns={"PTID": "Patient_ID",
                                                              prs_col: "PRS"})
    base["Patient_ID"] = base["Patient_ID"].astype(str)
    fits = {}
    for sp in ("train", "val", "test"):
        f = parts[sp].merge(base, on="Patient_ID", how="inner")
        f = f.merge(cov, on="Patient_ID", how="left")
        fits[sp] = f
    tr = fits["train"]
    if tr.empty or len(tr["y"].unique()) < 2:
        return None
    mu = tr["PRS"].mean(); sd = tr["PRS"].std(ddof=0) or 1.0
    for sp in fits:
        fits[sp] = fits[sp].copy()
        fits[sp]["PRS_z"] = (fits[sp]["PRS"] - mu) / sd

    Xcols = ["PRS_z"] + covars
    needed = ["y"] + Xcols
    tr_in = fits["train"][needed].dropna()
    val_in = fits["val"][needed].dropna()
    test_in = fits["test"][needed].dropna()
    if len(tr_in) < 20 or len(tr_in["y"].unique()) < 2:
        return None
    pipe = Pipeline([("imp", SimpleImputer(strategy="median")),
                      ("clf", LogisticRegression(max_iter=2000,
                                                   class_weight="balanced",
                                                   random_state=seed))])
    pipe.fit(tr_in[Xcols].astype(float), tr_in["y"].astype(int))
    train_p = pipe.predict_proba(tr_in[Xcols].astype(float))[:, 1]
    val_p   = pipe.predict_proba(val_in[Xcols].astype(float))[:, 1]
    test_p  = pipe.predict_proba(test_in[Xcols].astype(float))[:, 1]
    val_m  = _metrics(val_in["y"], val_p)
    test_m = _metrics(test_in["y"], test_p)

    # Per-patient prediction capture (multimodal alignment).
    if src in PER_PATIENT_SOURCES:
        PER_PATIENT_ROWS.append(capture_predictions(
            fits=fits, tr_in=tr_in, val_in=val_in, test_in=test_in,
            train_p=train_p, val_p=val_p, test_p=test_p,
            source=src, seed=seed, covar_mode=mode,
            beta_source=PER_PATIENT_BETA_SOURCE[0],
            ld_config=PER_PATIENT_LD_CONFIG[0],
            task="classification"))

    # OR per 1-SD of PRS_z (which is already 1-SD scaled on train); 95% CI via
    # Wald from the train-fold Hessian (X^T W X, W = p(1-p)). PRS_z lives in
    # column 0 of Xcols, so its SE is the [1, 1] entry after prepending an
    # intercept column. sklearn's L2 regularization mildly biases the coef
    # toward 0 — CI is still reported as approximate Wald.
    coef = float(pipe.named_steps["clf"].coef_[0][0])
    or_per_sd = math.exp(coef)
    try:
        Xtr_imp = pipe.named_steps["imp"].transform(
            tr_in[Xcols].astype(float).values)
        ptr = pipe.predict_proba(tr_in[Xcols].astype(float))[:, 1]
        Wd  = ptr * (1.0 - ptr)
        Xa  = np.column_stack([np.ones(len(Xtr_imp)), Xtr_imp])
        H   = Xa.T @ (Wd[:, None] * Xa)
        cov_mat = np.linalg.inv(H)
        se = float(np.sqrt(np.diag(cov_mat))[1])  # PRS_z is feature 0
        or_lo = math.exp(coef - 1.96 * se)
        or_hi = math.exp(coef + 1.96 * se)
        or_p  = float(2 * (1 - norm.cdf(abs(coef / se)))) if se > 0 else float("nan")
    except (np.linalg.LinAlgError, ValueError, FloatingPointError):
        se, or_lo, or_hi, or_p = float("nan"), float("nan"), float("nan"), float("nan")

    # Val-fold OR per +1 SD + 95% CI: refit logistic on VAL with PRS_z + covars.
    # PRS_z is already z-scored using TRAIN mean/sd, so the val-fold slope is
    # comparable across folds. Removes train-fold data-snooping bias (esp. for
    # EN-derived PRS where the "PRS" is itself trained on TRAIN dosage).
    or_val, or_val_lo, or_val_hi, or_val_p = (float("nan"),) * 4
    try:
        if len(val_in) >= 8 and val_in["y"].nunique() >= 2:
            pipe_v = Pipeline([("imp", SimpleImputer(strategy="median")),
                                ("clf", LogisticRegression(max_iter=2000,
                                                             class_weight="balanced",
                                                             random_state=seed))])
            pipe_v.fit(val_in[Xcols].astype(float), val_in["y"].astype(int))
            coef_v = float(pipe_v.named_steps["clf"].coef_[0][0])
            or_val = math.exp(coef_v)
            Xv_imp = pipe_v.named_steps["imp"].transform(
                val_in[Xcols].astype(float).values)
            pv = pipe_v.predict_proba(val_in[Xcols].astype(float))[:, 1]
            Wv = pv * (1.0 - pv)
            Xa_v = np.column_stack([np.ones(len(Xv_imp)), Xv_imp])
            Hv = Xa_v.T @ (Wv[:, None] * Xa_v)
            covv = np.linalg.inv(Hv)
            se_v = float(np.sqrt(np.diag(covv))[1])
            or_val_lo = math.exp(coef_v - 1.96 * se_v)
            or_val_hi = math.exp(coef_v + 1.96 * se_v)
            or_val_p  = float(2 * (1 - norm.cdf(abs(coef_v / se_v)))) if se_v > 0 else float("nan")
    except (np.linalg.LinAlgError, ValueError, FloatingPointError):
        pass

    # Test-fold OR per +1 SD + 95% CI: same refit logic on TEST. Used by the
    # test-split master leaderboard for held-out effect-size + significance.
    or_test, or_test_lo, or_test_hi, or_test_p = (float("nan"),) * 4
    try:
        if len(test_in) >= 8 and test_in["y"].nunique() >= 2:
            pipe_t = Pipeline([("imp", SimpleImputer(strategy="median")),
                                ("clf", LogisticRegression(max_iter=2000,
                                                             class_weight="balanced",
                                                             random_state=seed))])
            pipe_t.fit(test_in[Xcols].astype(float), test_in["y"].astype(int))
            coef_t = float(pipe_t.named_steps["clf"].coef_[0][0])
            or_test = math.exp(coef_t)
            Xt_imp = pipe_t.named_steps["imp"].transform(
                test_in[Xcols].astype(float).values)
            pt = pipe_t.predict_proba(test_in[Xcols].astype(float))[:, 1]
            Wt = pt * (1.0 - pt)
            Xa_t = np.column_stack([np.ones(len(Xt_imp)), Xt_imp])
            Ht = Xa_t.T @ (Wt[:, None] * Xa_t)
            covt = np.linalg.inv(Ht)
            se_t = float(np.sqrt(np.diag(covt))[1])
            or_test_lo = math.exp(coef_t - 1.96 * se_t)
            or_test_hi = math.exp(coef_t + 1.96 * se_t)
            or_test_p  = float(2 * (1 - norm.cdf(abs(coef_t / se_t)))) if se_t > 0 else float("nan")
    except (np.linalg.LinAlgError, ValueError, FloatingPointError):
        pass

    out = {"source": src, "covar_mode": mode, "seed": seed,
            "or_per_sd": or_per_sd,
            "or_per_sd_lo": or_lo, "or_per_sd_hi": or_hi, "or_per_sd_p": or_p,
            "or_per_sd_val": or_val,
            "or_per_sd_val_lo": or_val_lo, "or_per_sd_val_hi": or_val_hi,
            "or_per_sd_val_p":  or_val_p,
            "or_per_sd_test": or_test,
            "or_per_sd_test_lo": or_test_lo, "or_per_sd_test_hi": or_test_hi,
            "or_per_sd_test_p": or_test_p}
    for k, v in val_m.items():
        out[f"val_{k}"] = v
    for k, v in test_m.items():
        out[f"test_{k}"] = v
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=None)
    ap.add_argument("--ld-config", default=DEFAULT_LD_CONFIG)
    ap.add_argument("--beta-source", default="raw", choices=["raw","prscs"],
                     help="raw = published lead-SNP β; prscs = PRS-CS posterior β.")
    args = ap.parse_args()
    if args.beta_source == "raw":
        out_dir = OUT_BASE / args.ld_config / "classification"
    else:
        # PRS-CS uses the unpruned 115-SNP pool; --ld-config is methodologically
        # irrelevant. Output goes to a single flat dir (no ld_<cfg>/ subdir).
        out_dir = OUT_BASE.parent / f"strict_qc_prs_{args.beta_source}" / "classification"
    out_dir.mkdir(parents=True, exist_ok=True)
    PER_PATIENT_BETA_SOURCE[0] = args.beta_source
    PER_PATIENT_LD_CONFIG[0]   = (args.ld_config if args.beta_source == "raw" else "na")

    print(f"Building per-source PRS table for {args.ld_config}  beta_source={args.beta_source}...")
    prs_full, snps_per_src = per_source_prs_table(ld_config=args.ld_config,
                                                     beta_source=args.beta_source)
    cov = load_subject_covariates()
    dedup_st, dedup_dosage = get_dedup_dosage_matrix(args.ld_config)
    n_snp_per_src = {s: len(st) for s, st in snps_per_src.items()}
    n_snp_per_src["prs_all_dedup_EN_dosage"] = dedup_dosage.shape[1]
    # meta_prs uses the union of source PRS values; "n_snps" is total
    # unique rsIDs in the LD-pruned pool (proxy — same as dedup size).
    n_snp_per_src["meta_prs_EN_combined"] = dedup_dosage.shape[1]
    # PRS-CS leaderboard: 6 individual PRS-CS sources only (no dedup, no EN).
    # Raw-β leaderboard: 16 individual sources + dedup + 2 EN ablation rows.
    if args.beta_source == "prscs":
        source_list = list(PRSCS_SOURCES) + ["covariates_only"]
    else:
        source_list = ALL_SOURCES + ["prs_all_dedup",
                                       "prs_all_dedup_ivw",
                                       "prs_all_dedup_filtered",
                                       "prs_all_dedup_EN_dosage",
                                       "meta_prs_EN_combined",
                                       "meta_prs_EN_filtered",
                                       "covariates_only"]
    # covariates_only carries no SNPs; mark explicitly.
    n_snp_per_src["covariates_only"] = 0
    # prs_all_dedup_ivw shares the LD-pruned union with the priority-pick dedup.
    if "prs_all_dedup_ivw" in snps_per_src:
        n_snp_per_src["prs_all_dedup_ivw"] = len(snps_per_src["prs_all_dedup_ivw"])
    # prs_all_dedup_filtered: subset of dedup SNPs by majority-positive EN coef.
    # Empty during the bootstrap pass (en_dosage_coefs.tsv files not yet written).
    if "prs_all_dedup_filtered" in snps_per_src:
        n_snp_per_src["prs_all_dedup_filtered"] = len(snps_per_src["prs_all_dedup_filtered"])
    # meta_prs_EN_filtered: same dosage matrix size as meta (proxy n_snps).
    n_snp_per_src["meta_prs_EN_filtered"] = dedup_dosage.shape[1]
    sources = [args.source] if args.source else [
        s for s in source_list
        if s == "covariates_only" or n_snp_per_src.get(s, 0) > 0]
    print(f"Sources with >=1 SNP after LD prune + orient: {len(sources)}")

    rows = []
    for src in sources:
        for mode, covs in COVAR_MODES.items():
            for seed in SEEDS:
                r = _run_one(src, prs_full, cov, seed, mode, covs,
                              dedup_dosage=dedup_dosage)
                if r is None: continue
                r["n_snps_used"] = n_snp_per_src.get(src, 0)
                rows.append(r)
    df = pd.DataFrame(rows)
    per_p = out_dir / "per_run_metrics.tsv"
    df.to_csv(per_p, sep="\t", index=False)
    print(f"wrote per-run metrics: {per_p}")

    # Aggregate full metric set
    metric_cols = ["val_auc","val_bacc","val_f1","val_sens","val_spec",
                    "val_prec","val_pr_auc","or_per_sd",
                    "or_per_sd_lo","or_per_sd_hi","or_per_sd_p",
                    "or_per_sd_val","or_per_sd_val_lo","or_per_sd_val_hi",
                    "or_per_sd_val_p",
                    "or_per_sd_test","or_per_sd_test_lo","or_per_sd_test_hi",
                    "or_per_sd_test_p",
                    "test_auc","test_bacc","test_f1"]
    g = df.groupby(["source","covar_mode"]).agg(
        n_snps=("n_snps_used","first"),
        n_seeds=("seed","count"),
        val_n=("val_n","mean"), val_pos=("val_pos","mean"), val_neg=("val_neg","mean"),
        val_tn=("val_tn","mean"), val_fp=("val_fp","mean"),
        val_fn=("val_fn","mean"), val_tp=("val_tp","mean"),
        **{f"{c}_mean": (c, "mean") for c in metric_cols},
        **{f"{c}_std":  (c, "std")  for c in metric_cols},
    ).reset_index().sort_values(["covar_mode","val_bacc_mean"], ascending=[True, False])
    lb_p = out_dir / "leaderboard.csv"
    g.to_csv(lb_p, index=False)
    print(f"wrote leaderboard: {lb_p}")

    # PNG — concise; primary metric (val_bacc) first
    show = pd.DataFrame({
        "source":     g["source"],
        "covar_mode": g["covar_mode"],
        "n_snps":     g["n_snps"].astype(int),
        "val_bacc":   [f"{m:.3f}+/-{s:.3f}" for m, s in zip(g["val_bacc_mean"], g["val_bacc_std"])],
        "val_auc":    [f"{m:.3f}+/-{s:.3f}" for m, s in zip(g["val_auc_mean"], g["val_auc_std"])],
        "val_f1":     [f"{m:.3f}+/-{s:.3f}" for m, s in zip(g["val_f1_mean"], g["val_f1_std"])],
        "OR/+1SD (val) [95% CI]": [
            (f"{m:.2f}" if (pd.isna(lo) or pd.isna(hi))
             else f"{m:.2f} [{lo:.2f}-{hi:.2f}]")
            for m, lo, hi in zip(g["or_per_sd_val_mean"],
                                   g["or_per_sd_val_lo_mean"],
                                   g["or_per_sd_val_hi_mean"])],
    })
    fig_h = max(2.5, 0.30 * len(show) + 1.3)
    fig, ax = plt.subplots(figsize=(14, fig_h)); ax.axis("off")
    ax.set_title(f"Classification — strict-QC + LD pruning [{args.ld_config}]   "
                 f"(VAL set; mean +/- std over 3 seeds)",
                 fontsize=11, fontweight="bold", pad=12, loc="left")
    col_widths = [0.26, 0.16, 0.06, 0.13, 0.13, 0.13, 0.18]
    tbl = ax.table(cellText=show.values.tolist(), colLabels=show.columns.tolist(),
                    loc="center", cellLoc="center", colWidths=col_widths)
    tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1, 1.2)
    for j in range(len(show.columns)):
        c = tbl[(0, j)]; c.set_facecolor("#222"); c.set_text_props(color="white", weight="bold")
    # Bold BalAcc > 0.5, and OR/+1SD (val) when its 95% CI excludes 1.0
    # (i.e. statistically significant departure from OR=1).
    g_lo = g["or_per_sd_val_lo_mean"].tolist()
    g_hi = g["or_per_sd_val_hi_mean"].tolist()
    for i, row in enumerate(show.itertuples(index=False), start=1):
        mean_bacc = float(row.val_bacc.split("+/-")[0])
        if mean_bacc > 0.5:
            tbl[(i, 3)].set_text_props(weight="bold")
        lo, hi = g_lo[i - 1], g_hi[i - 1]
        if not (pd.isna(lo) or pd.isna(hi)) and (lo > 1.0 or hi < 1.0):
            tbl[(i, 6)].set_text_props(weight="bold")
    png = out_dir / "leaderboard.png"
    fig.savefig(png, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"wrote PNG: {png}")

    # Per-patient prediction parquets (multimodal-alignment anchors).
    # Per-(source x covar_mode x seed) parquet in out_dir/per_patient/.
    # Per-(source x seed) anchor parquet in
    # outputs/strict_qc_prs/multimodal_anchors/prs/classification/
    # (all covar_modes long-format) — one canonical file per source x seed
    # that downstream fusion code can join on Patient_ID x seed x fold.
    if PER_PATIENT_ROWS:
        pp_dir = out_dir / "per_patient"
        anchor_dir = (OUT_BASE / "multimodal_anchors" / "prs" / "classification"
                       if args.beta_source == "raw"
                       else OUT_BASE.parent / f"strict_qc_prs_{args.beta_source}"
                            / "multimodal_anchors" / "prs" / "classification")
        all_rows = pd.concat(PER_PATIENT_ROWS, ignore_index=True)
        # Per-(source, covar_mode, seed) split for the in-task dir
        for (src, mode, seed), grp in all_rows.groupby(
                ["model", "covar_mode", "seed"], sort=False):
            save_per_patient_parquet(grp, pp_dir,
                                       source=src, seed=int(seed),
                                       covar_mode=mode)
        # Per-(source, seed) anchor for multimodal fusion
        for (src, seed), grp in all_rows.groupby(
                ["model", "seed"], sort=False):
            concat_and_save_anchor([grp], anchor_dir,
                                     source=src, seed=int(seed))
        n_pids = all_rows["Patient_ID"].nunique()
        print(f"wrote per-patient parquets: {pp_dir} "
              f"({all_rows.shape[0]:,} rows, {n_pids} unique PTIDs)")
        print(f"wrote multimodal anchors:   {anchor_dir}")

    # Save train-fold EN dosage coefs to disk for downstream consumers
    # (get_dedup_filtered_rsids() reads from all 4 LD configs to majority-vote
    # the prs_all_dedup_filtered SNP set).
    if EN_COEFS["dosage"]:
        dos_coefs = pd.DataFrame(EN_COEFS["dosage"]).rename(
            columns={s: f"seed_{s}" for s in EN_COEFS["dosage"]})
        dos_coefs.index.name = "rsID"
        coef_tsv = out_dir / "en_dosage_coefs.tsv"
        dos_coefs.to_csv(coef_tsv, sep="\t")
        print(f"wrote EN dosage coefs: {coef_tsv}")

    # EN coefficient bar chart (per-feature mean +/- std across seeds).
    dos_df  = pd.DataFrame(EN_COEFS["dosage"]) if EN_COEFS["dosage"] else None
    meta_df = pd.DataFrame(EN_COEFS["meta"])   if EN_COEFS["meta"]   else None
    if dos_df is not None or meta_df is not None:
        en_png = out_dir / "en_weights.png"
        render_en_weights_png(
            dos_df, meta_df, en_png,
            title=(f"Classification — EN coefficients (TRAIN-fold fit)   "
                    f"LD={args.ld_config}  beta_source={args.beta_source}\n"
                    f"top: per-SNP dosage-EN log-odds per +1 A1   ·   "
                    f"bottom: meta-EN log-odds per +1 source-PRS"))
        print(f"wrote EN weights PNG: {en_png}")

    # Filtered meta-EN coefs (7-source subset) as a separate PNG so the
    # comparison stays readable next to the full meta panel.
    meta_filt_df = (pd.DataFrame(EN_COEFS["meta_filtered"])
                     if EN_COEFS["meta_filtered"] else None)
    if meta_filt_df is not None:
        en_filt_png = out_dir / "en_weights_filtered.png"
        render_en_weights_png(
            None, meta_filt_df, en_filt_png,
            title=(f"Classification — meta_prs_EN_filtered coefficients   "
                    f"LD={args.ld_config}  beta_source={args.beta_source}\n"
                    f"7-source subset: Kosteridis, Leonenko, Lambert, "
                    f"Schwanzentruber, Kunkle, Bellenguez, Najar"))
        print(f"wrote filtered meta-EN PNG: {en_filt_png}")

    # Test-fold EN coefficient refit — mirrors the OR-on-test refit logic.
    # Fits a fresh EN on the TEST fold so the bar chart reflects what EN
    # selects when trained on held-out data. NB: small test sets (~30) mean
    # high coefficient variance — caveat reader.
    EN_COEFS_TEST = {"dosage": {}, "meta": {}}
    for seed in SEEDS:
        parts_seed = {sp: _load_labels(seed, sp) for sp in ("train","val","test")}
        # Swap: put test in the "train" slot so helpers fit on test data.
        parts_fit = {"train": parts_seed["test"],
                      "val":   parts_seed["val"],
                      "test":  parts_seed["train"]}
        clf_d, _ = _fit_en_dosage_get_prs(dedup_dosage, parts_fit, cov, seed)
        if clf_d is not None:
            EN_COEFS_TEST["dosage"][seed] = pd.Series(
                clf_d.coef_[0], index=list(dedup_dosage.columns))
        # Meta-EN: same exclusion list as the main loop's meta path.
        meta_source_cols = [
            c for c in prs_full.columns if c.startswith("PRS_")
            and c not in ("PRS_prs_all_dedup", "PRS_prs_all_dedup_ivw",
                            "PRS_Kosteridis_MTAG_AD",
                            "PRS_Kosteridis_shared_AD_CV",
                            "PRS_Kosteridis_novel_AD")
            and not c.endswith("EN_dosage") and not c.endswith("EN_combined")]
        sp_df = prs_full[["PTID"] + meta_source_cols].set_index("PTID")
        sp_df.index = sp_df.index.astype(str)
        clf_m, _ = _fit_en_meta_prs(sp_df, parts_fit, seed)
        if clf_m is not None:
            feats = getattr(clf_m, "_feature_names", meta_source_cols)
            EN_COEFS_TEST["meta"][seed] = pd.Series(
                clf_m.coef_[0],
                index=[c.replace("PRS_", "", 1) for c in feats])
    dos_df_t  = pd.DataFrame(EN_COEFS_TEST["dosage"]) if EN_COEFS_TEST["dosage"] else None
    meta_df_t = pd.DataFrame(EN_COEFS_TEST["meta"])   if EN_COEFS_TEST["meta"]   else None
    if dos_df_t is not None or meta_df_t is not None:
        en_png_t = out_dir / "en_weights_test.png"
        render_en_weights_png(
            dos_df_t, meta_df_t, en_png_t,
            title=(f"Classification — EN coefficients (TEST-fold refit)   "
                    f"LD={args.ld_config}  beta_source={args.beta_source}\n"
                    f"top: per-SNP dosage-EN log-odds per +1 A1   ·   "
                    f"bottom: meta-EN log-odds per +1 source-PRS\n"
                    f"caveat: small test sets (~30) -> high coef variance"))
        print(f"wrote EN weights PNG (TEST): {en_png_t}")


if __name__ == "__main__":
    main()
