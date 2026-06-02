"""
_covars_only_new_tasks.py
=========================
Covariates-only baseline arms (no PRS, no FM feature) for the new-task
leaderboards (49a sCN_vs_pCN, 49b T3 horizons, 49c T4 multiclass).

Lets the reader see by how much PRS / FM adds over age+sex+APOE alone --
without this baseline a "val_bacc=0.55" row is uncalibrated.

Each fit:
  - Pipeline = SimpleImputer(median) -> LogisticRegression(class_weight='balanced',
    multinomial for T4). Same head class as scripts 45/49a/49b's main PRS arm,
    just with the PRS column omitted.
  - Returns a single per-(seed, covar_set) row matching the host script's
    per_run_metrics.tsv schema, with `source='covariates_only'` and
    `covar_mode='covariates_only_<covar_set>'` so it appears as a distinct
    leaderboard row instead of overwriting a PRS row.

ALL FUNCTIONS HERE FIT ON THE TRAINING FOLD AND EVALUATE ON VAL + TEST,
mirroring the host script's evaluation convention.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (roc_auc_score, balanced_accuracy_score,
                              recall_score, precision_score, f1_score,
                              average_precision_score, confusion_matrix)
from sklearn.pipeline import Pipeline


COVAR_SETS_BINARY = {
    "covariates_only_age+sex+apoe4":
        ["age_at_baseline", "Sex_M", "APOE4_Dosage"],
    "covariates_only_age+sex+apoe4+apoe2":
        ["age_at_baseline", "Sex_M", "APOE4_Dosage", "APOE2_Dosage"],
}


def _cls_binary_metrics(y_true, y_score):
    y_true  = np.asarray(y_true).astype(int)
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
             "n": int(len(y_true)),
             "pos": int((y_true == 1).sum()), "neg": int((y_true == 0).sum()),
             "tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)}


def _cls_multi_metrics(y_true, y_proba, n_classes):
    y_true = np.asarray(y_true).astype(int)
    y_proba = np.asarray(y_proba).astype(float)
    pred = y_proba.argmax(axis=1)
    bacc = balanced_accuracy_score(y_true, pred)
    f1   = f1_score(y_true, pred, average="macro", zero_division=0)
    try:    auc = roc_auc_score(y_true, y_proba,
                                  multi_class="ovr", average="macro",
                                  labels=list(range(n_classes)))
    except: auc = float("nan")
    n_per_class = [int((y_true == k).sum()) for k in range(n_classes)]
    return {"bacc_macro": bacc, "auc_macro": auc, "f1_macro": f1,
             "n": int(len(y_true)),
             **{f"n_class{k}": n_per_class[k] for k in range(n_classes)}}


def fit_binary(cov_df: pd.DataFrame, covars: list[str], covar_mode: str,
                  seed: int,
                  fold_ptid: dict[str, list[str]],
                  fold_y: dict[str, np.ndarray]) -> dict | None:
    """Binary covariates-only fit on an arbitrary (fold_ptid, fold_y) split.
    Shared between 49a (sCN_vs_pCN) and 49b (per horizon).
    `fold_ptid` / `fold_y` are dicts keyed by 'train'/'val'/'test'."""
    cov_indexed = cov_df.set_index("Patient_ID")
    needed = covars

    def _slice(sp: str):
        pids = fold_ptid[sp]
        sub = cov_indexed.reindex(pids)[needed].copy()
        sub["y"] = fold_y[sp]
        return sub.dropna(subset=["y"])

    tr   = _slice("train")
    va   = _slice("val")
    te   = _slice("test")
    if tr["y"].nunique() < 2 or len(tr) < 20: return None
    pipe = Pipeline([("imp", SimpleImputer(strategy="median")),
                      ("clf", LogisticRegression(max_iter=2000,
                                                   class_weight="balanced",
                                                   random_state=seed))])
    pipe.fit(tr[needed].astype(float), tr["y"].astype(int))
    val_p  = pipe.predict_proba(va[needed].astype(float))[:, 1] if len(va) else np.array([])
    test_p = pipe.predict_proba(te[needed].astype(float))[:, 1] if len(te) else np.array([])
    val_m  = _cls_binary_metrics(va["y"], val_p) if len(va) else None
    test_m = _cls_binary_metrics(te["y"], test_p) if len(te) else None
    if val_m is None or test_m is None: return None
    row = {"source": "covariates_only", "covar_mode": covar_mode, "seed": seed}
    for k, v in val_m.items():  row[f"val_{k}"]  = v
    for k, v in test_m.items(): row[f"test_{k}"] = v
    return row


def fit_multiclass(cov_df: pd.DataFrame, covars: list[str], covar_mode: str,
                     seed: int,
                     fold_ptid: dict[str, list[str]],
                     fold_y: dict[str, np.ndarray],
                     n_classes: int) -> dict | None:
    """Multiclass covariates-only fit for 49c T4_3class."""
    cov_indexed = cov_df.set_index("Patient_ID")
    needed = covars

    def _slice(sp: str):
        pids = fold_ptid[sp]
        sub = cov_indexed.reindex(pids)[needed].copy()
        sub["y"] = fold_y[sp]
        return sub.dropna(subset=["y"])

    tr   = _slice("train")
    va   = _slice("val")
    te   = _slice("test")
    if tr["y"].nunique() < 2 or len(tr) < 20: return None
    pipe = Pipeline([("imp", SimpleImputer(strategy="median")),
                      ("clf", LogisticRegression(max_iter=2000,
                                                   class_weight="balanced",
                                                   random_state=seed))])
    pipe.fit(tr[needed].astype(float), tr["y"].astype(int))
    val_p  = pipe.predict_proba(va[needed].astype(float)) if len(va) else np.zeros((0, n_classes))
    test_p = pipe.predict_proba(te[needed].astype(float)) if len(te) else np.zeros((0, n_classes))
    val_m  = _cls_multi_metrics(va["y"], val_p,  n_classes) if len(va) else None
    test_m = _cls_multi_metrics(te["y"], test_p, n_classes) if len(te) else None
    if val_m is None or test_m is None: return None
    row = {"source": "covariates_only", "covar_mode": covar_mode, "seed": seed}
    for k, v in val_m.items():  row[f"val_{k}"]  = v
    for k, v in test_m.items(): row[f"test_{k}"] = v
    return row
