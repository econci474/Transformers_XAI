"""
01_fuse_T1e.py
==============
Late fusion for T1e_scn_pcn (binary: sCN=0, pCN=1) over per-patient class
probabilities from a CLINICAL encoder and a SNP-side model.

Per integration/T1e_classification/notes.txt: 2 clinical encoders x 4 SNP
models = 8 fusion combos.

  CL variants (no_cdr_post_exclusion / baseline parquet):
    - BioClinical-ModernBERT-base / frozen           ("BioBERT frozen")
    - ModernBERT-base             / full_ft          ("ModernBERT B ft")

  SNP variants (per-seed per-patient parquet, val+test only):
    - BMFM-SNP attn_bias_single chrom_hier MLP [attn=T1_AD_CN]   (multimodal_anchors/...)
    - BMFM-SNP none global_attn XGB                              (multimodal_anchors/...)
    - meta_prs_EN_combined                                       (per_patient/...)
    - Kosteridis_full_PRS                                        (per_patient/...)

Per the user pick, WEIGHTED AVERAGE ONLY (no stacked elastic-net here). For
each (CL, SNP, seed) we report 4 methods on val + test:
  clinical_only : argmax of CL probs                                 (reference)
  snp_only      : argmax of SNP probs                                (reference)
  equal_avg     : 0.5*CL + 0.5*SNP per-class prob
  best_alpha    : alpha*CL + (1-alpha)*SNP, alpha tuned on VAL bACC

Each base model already saw TRAIN; the meta-blend alpha is tuned on VAL only
(no train-fold leakage), and TEST is held out for reporting. SNP parquets
contain val+test rows only -- they never persist train-fold predictions.

Inputs (local D:):
  Clinical baseline parquet:
    encoder_outputs_no_cdr_post_exclusion/<model>/T1e_scn_pcn/seed_{s}/<strat>/
        embeddings/embeddings.parquet
    cols: Patient_ID, split, in_task_cohort, label, prob_0, prob_1
  SNP per-seed per-patient parquet (val+test only):
    cols: Patient_ID, seed, fold, outcome_proba, y_true (+ source metadata)

Outputs (outputs/baseline/):
  fusion_metrics_per_seed.csv  fusion_metrics.csv  coverage.csv
  fused_predictions.csv        failure_table.csv   fusion_table_T1e.png

Usage:
  python integration/T1e_classification/01_fuse_T1e.py
  python integration/T1e_classification/01_fuse_T1e.py \
      --clin-keys bioclin_mb_base_frozen \
      --snp-keys kosteridis_prs meta_prs
"""
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (balanced_accuracy_score, f1_score,
                              precision_recall_fscore_support,
                              recall_score, roc_auc_score)
from sklearn.model_selection import GridSearchCV, StratifiedKFold

warnings.filterwarnings("ignore")

CLASS_NAMES = ["sCN", "pCN"]
N_CLASSES = 2
SEEDS_DEFAULT = (0, 1, 2)

# --------------------------------------------------------------------------- #
# Paths and variant registries
# --------------------------------------------------------------------------- #
CLIN_ROOT = Path(r"D:/ADNI_BIDS_project/derivatives/encoder_outputs_no_cdr_post_exclusion")
SNP_ROOT = Path(r"D:/ADNI_SNP_Omni2.5M_20140220")
ANCHORS = SNP_ROOT / "outputs" / "multimodal_anchors"
NEW_TASKS_PERPT = (SNP_ROOT / "outputs" / "new_tasks" / "sCN_vs_pCN"
                   / "ld_1000kb_r2_0.8" / "per_patient")
OUT_DIR_DEFAULT = Path(__file__).resolve().parent / "outputs" / "baseline"

CLIN_VARIANTS = [
    {"key": "bioclin_mb_base_frozen",
     "display": "BioClinical-MB-base [frozen]",
     "short": "BioClin-B/frz",
     "model": "BioClinical-ModernBERT-base", "strat": "frozen"},
    {"key": "modernbert_base_full_ft",
     "display": "ModernBERT-base [full_ft]",
     "short": "ModBERT-B/ft",
     "model": "ModernBERT-base", "strat": "full_ft"},
]

SNP_VARIANTS = [
    {"key": "bmfm_snp_attn1_mlp",
     "display": "BMFM-SNP attn_bias_single chrom_hier MLP [attn=T1_AD_CN]",
     "short": "BMFM/attn1-MLP",
     "parquet_tmpl": str(ANCHORS / "bmfm_snp__attn_bias_single__chrom_hier__mlp2__101bp"
                         / "task_sCN_vs_pCN" / "seed_{seed}" / "prs+age+sex+apoe4"
                         / "predictions.parquet")},
    {"key": "bmfm_snp_none_xgb",
     "display": "BMFM-SNP none global_attn XGB",
     "short": "BMFM/none-XGB",
     "parquet_tmpl": str(ANCHORS / "bmfm_snp__none__global_attn__xgb__1001bp"
                         / "task_sCN_vs_pCN" / "seed_{seed}" / "prs+age+sex+apoe4"
                         / "predictions.parquet")},
    {"key": "meta_prs",
     "display": "meta_prs_EN_combined",
     "short": "meta_PRS",
     "parquet_tmpl": str(NEW_TASKS_PERPT
                         / "meta_prs_EN_combined__prs+age+sex+apoe4__seed{seed}.parquet")},
    {"key": "kosteridis_prs",
     "display": "Kosteridis_full_PRS",
     "short": "Kosteridis_PRS",
     "parquet_tmpl": str(NEW_TASKS_PERPT
                         / "Kosteridis__prs+age+sex+apoe4__seed{seed}.parquet")},
]


# --------------------------------------------------------------------------- #
# Loaders
# --------------------------------------------------------------------------- #
def load_clinical(seed: int, model_slug: str, strategy: str,
                  task: str = "T1e_scn_pcn",
                  clin_root: Path = CLIN_ROOT) -> pd.DataFrame:
    """Clinical probs -> Patient_ID, split, y_clin, cp0, cp1 (binary)."""
    p = (clin_root / model_slug / task / f"seed_{seed}" / strategy
         / "embeddings" / "embeddings.parquet")
    if not p.exists():
        raise FileNotFoundError(f"clinical parquet missing: {p}")
    cols = ["Patient_ID", "split", "in_task_cohort", "label", "prob_0", "prob_1"]
    d = pd.read_parquet(p, columns=cols)
    d = d[d["in_task_cohort"] & d["label"].notna()].copy()
    d["y_clin"] = d["label"].astype(int)
    d = d.rename(columns={"prob_0": "cp0", "prob_1": "cp1"})
    return d[["Patient_ID", "split", "y_clin", "cp0", "cp1"]]


def load_snp(seed: int, parquet_tmpl: str) -> pd.DataFrame:
    """SNP probs (val+test only) -> Patient_ID, split, y_snp, sp0, sp1."""
    p = Path(parquet_tmpl.format(seed=seed))
    if not p.exists():
        raise FileNotFoundError(f"SNP parquet missing: {p}")
    d = pd.read_parquet(p, columns=["Patient_ID", "fold", "outcome_proba", "y_true"])
    d = d.rename(columns={"fold": "split", "y_true": "y_snp"})
    d["sp1"] = d["outcome_proba"].astype(float)
    d["sp0"] = 1.0 - d["sp1"]
    d["y_snp"] = d["y_snp"].astype(int)
    return d[["Patient_ID", "split", "y_snp", "sp0", "sp1"]]


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #
def safe_auc(y_true: np.ndarray, prob_pos: np.ndarray) -> float:
    try:
        if len(np.unique(y_true)) < 2:
            return float("nan")
        return float(roc_auc_score(y_true, prob_pos))
    except Exception:
        return float("nan")


def metric_row(y_true: np.ndarray, y_pred: np.ndarray,
               probs: np.ndarray | None) -> dict:
    bacc = balanced_accuracy_score(y_true, y_pred)
    macro_f1 = f1_score(y_true, y_pred, average="macro",
                        labels=[0, 1], zero_division=0)
    _, rc, f1, _ = precision_recall_fscore_support(
        y_true, y_pred, labels=[0, 1], zero_division=0)
    out = {"n": int(len(y_true)), "bacc": float(bacc),
           "macro_f1": float(macro_f1)}
    out["auc"] = (safe_auc(y_true, probs[:, 1])
                  if probs is not None else float("nan"))
    for i, name in enumerate(CLASS_NAMES):
        out[f"f1_{name}"] = float(f1[i])
        out[f"recall_{name}"] = float(rc[i])
    return out


# --------------------------------------------------------------------------- #
# Fusion (weighted average, alpha tuned on VAL)
# --------------------------------------------------------------------------- #
def fit_weighted_avg(cl_v, sn_v, yv):
    """Pick alpha in [0,1] maximising VAL balanced accuracy at threshold 0.5.
    Returns alpha_star."""
    best_w, best_bacc = 0.5, -1.0
    for w in np.linspace(0.0, 1.0, 21):
        blend = w * cl_v + (1.0 - w) * sn_v
        b = balanced_accuracy_score(yv, blend.argmax(1))
        if b > best_bacc:
            best_bacc, best_w = b, float(w)
    return best_w


def fit_alpha_J(cp_v, sp_v, y_v, cp_t, sp_t):
    """Jointly grid-search (alpha, threshold) on VAL Youden's J = TPR - FPR.
    Per alpha, picks the threshold maximising J on VAL probabilities, then
    picks the alpha whose best-threshold J is highest. At test time, classifies
    with (alpha*, t*).

    Returns:
      (test_pred (int array), test_proba (N,2) blend (NO threshold applied),
       val_pred (int array), val_proba (N,2) blend,
       alpha_star, threshold_star, val_J_star)

    For the stored 2-column proba we use the blended distribution (cols sum
    to 1) so downstream metric_row/auc work; the threshold is reported but
    only used at the pred step. test_pred uses (alpha*, t*); val_pred too
    so the metrics block reflects the joint tuning."""
    alpha_grid = np.linspace(0.0, 1.0, 21)
    thr_grid = np.linspace(0.05, 0.95, 19)
    best = (-np.inf, 1.0, 0.5)
    for a in alpha_grid:
        blend_v = a * cp_v + (1.0 - a) * sp_v
        p_pos_v = blend_v[:, 1]
        for t in thr_grid:
            pred_v = (p_pos_v >= t).astype(int)
            tpr = recall_score(y_v, pred_v, pos_label=1, zero_division=0)
            tnr = recall_score(y_v, pred_v, pos_label=0, zero_division=0)
            j = tpr + tnr - 1.0
            if j > best[0]:
                best = (j, float(a), float(t))
    j_star, a_star, t_star = best

    blend_v = a_star * cp_v + (1.0 - a_star) * sp_v
    blend_t = a_star * cp_t + (1.0 - a_star) * sp_t
    val_pred = (blend_v[:, 1] >= t_star).astype(int)
    test_pred = (blend_t[:, 1] >= t_star).astype(int)
    return (test_pred, blend_t, val_pred, blend_v,
            a_star, t_star, j_star)


def fit_class_avg(cp_v, sp_v, y_v, cp_t, sp_t):
    """Per-class weighted average. Grid-search (w0, w1) on VAL bACC where
        p_unnorm[k] = w_k * P_CL[k] + (1 - w_k) * P_SNP[k]   for k in {0, 1}
        p_final     = p_unnorm / sum(p_unnorm)               (renormalise)
        pred        = argmax(p_final)

    Returns (test_pred, test_proba, val_pred, val_proba, w0_star, w1_star)."""
    grid = np.linspace(0.0, 1.0, 21)
    best = (-1.0, 1.0, 1.0)
    for w0 in grid:
        u0_v = w0 * cp_v[:, 0] + (1.0 - w0) * sp_v[:, 0]
        for w1 in grid:
            u1_v = w1 * cp_v[:, 1] + (1.0 - w1) * sp_v[:, 1]
            s = u0_v + u1_v
            s = np.where(s == 0, 1.0, s)
            p1 = u1_v / s
            p0 = u0_v / s
            pred = (p1 >= p0).astype(int)
            b = balanced_accuracy_score(y_v, pred)
            if b > best[0]:
                best = (b, float(w0), float(w1))
    _, w0s, w1s = best

    def _apply(cp, sp):
        u0 = w0s * cp[:, 0] + (1.0 - w0s) * sp[:, 0]
        u1 = w1s * cp[:, 1] + (1.0 - w1s) * sp[:, 1]
        s = u0 + u1
        s = np.where(s == 0, 1.0, s)
        proba = np.column_stack([u0 / s, u1 / s])
        return proba.argmax(1), proba

    val_pred,  val_proba  = _apply(cp_v, sp_v)
    test_pred, test_proba = _apply(cp_t, sp_t)
    return test_pred, test_proba, val_pred, val_proba, w0s, w1s


def fit_stacked_en(cp_v, sp_v, y_v, cp_t, sp_t, seed):
    """Stacked elastic-net meta-learner over [CL, SNP] class probs.
    Fits on VAL, predicts on val + test. CV grid on (C, l1_ratio) when VAL has
    >=2 samples per class; otherwise a fixed mild-EN fallback. Returns
    (val_pred, val_proba, test_pred, test_proba, info_dict) where info_dict
    carries best_C and best_l1_ratio (the model_info)."""
    Xv = np.column_stack([cp_v, sp_v])
    Xt = np.column_stack([cp_t, sp_t])
    if len(np.unique(y_v)) < 2:
        return None
    n_min = int(np.min(np.bincount(y_v, minlength=2)))
    if n_min < 2:
        clf = LogisticRegression(penalty="elasticnet", solver="saga", C=1.0,
                                 l1_ratio=0.5, max_iter=5000, random_state=seed)
        clf.fit(Xv, y_v)
        best = {"C": 1.0, "l1_ratio": 0.5}
    else:
        grid = {"C": [0.1, 1.0, 10.0], "l1_ratio": [0.1, 0.5, 0.9]}
        base = LogisticRegression(penalty="elasticnet", solver="saga",
                                  max_iter=5000, random_state=seed)
        cv = StratifiedKFold(n_splits=min(3, n_min), shuffle=True,
                             random_state=seed)
        clf = GridSearchCV(base, grid, scoring="balanced_accuracy",
                           cv=cv, n_jobs=-1)
        clf.fit(Xv, y_v)
        best = {"C": float(clf.best_params_["C"]),
                "l1_ratio": float(clf.best_params_["l1_ratio"])}

    def _full(proba, classes):
        full = np.zeros((proba.shape[0], 2))
        for j, c in enumerate(classes):
            full[:, int(c)] = proba[:, j]
        return full

    classes = (clf.best_estimator_.classes_ if hasattr(clf, "best_estimator_")
               else clf.classes_)
    pv_proba = _full(clf.predict_proba(Xv), classes)
    pt_proba = _full(clf.predict_proba(Xt), classes)
    return (pv_proba.argmax(1), pv_proba,
            pt_proba.argmax(1), pt_proba, best)


def build_frame(clin: pd.DataFrame, snp: pd.DataFrame,
                fill_missing_snp: bool = False) -> pd.DataFrame:
    """Outer-join CL + SNP on (Patient_ID, split). Tag presence flags. Canonical
    label = clinical's (NaN-filled with SNP's where missing).

    Default fill_missing_snp=False: only patients with BOTH modalities present
    are scored (complete-case). The integration script does NOT impute missing
    SNP probs -- if a patient is missing from the BMFM cohort, the proper fix
    is upstream re-extraction (see snp_pipeline/50b_bmfm_extract_missing_ptids.py),
    not a neutral 0.5 fill. The ``fill_missing_snp=True`` option remains for
    diagnostic / what-if analysis only."""
    fr = clin.merge(snp, on=["Patient_ID", "split"], how="outer", indicator=True)
    fr["clin_present"] = fr["_merge"].isin(["both", "left_only"])
    fr["snp_present"] = fr["_merge"].isin(["both", "right_only"])
    fr = fr.drop(columns=["_merge"])
    fr["y_true"] = fr["y_clin"]
    fill_mask = fr["y_true"].isna()
    fr.loc[fill_mask, "y_true"] = fr.loc[fill_mask, "y_snp"]
    fr = fr[fr["y_true"].notna()].copy()
    fr["y_true"] = fr["y_true"].astype(int)
    if fill_missing_snp:
        miss = fr["clin_present"] & ~fr["snp_present"]
        fr.loc[miss, "sp0"] = 0.5
        fr.loc[miss, "sp1"] = 0.5
    return fr


def _stash_preds(pred_rows, clin_key, snp_key, seed, method,
                 test_df, pred, proba):
    for i, (_, row) in enumerate(test_df.reset_index(drop=True).iterrows()):
        pred_rows.append({
            "clin_variant": clin_key, "snp_variant": snp_key,
            "seed": seed, "method": method,
            "Patient_ID": row["Patient_ID"], "y_true": int(row["y_true"]),
            "pred": int(pred[i]),
            "prob_sCN": float(proba[i, 0]),
            "prob_pCN": float(proba[i, 1]),
            "clin_present": bool(row["clin_present"]),
            "snp_present": bool(row["snp_present"]),
        })


def evaluate_combo(seed, clin_v, snp_v, clin, snp, pred_rows, fail_rows):
    """All methods x {val, test} for one (CL, SNP) combo at one seed. Returns
    list of metric dicts. Complete-case: evaluates only on patients where BOTH
    modalities are present. Missing SNP rows (e.g. pCN_to_MCI patients not in
    the Wave-1 BMFM training cohort) require upstream re-extraction; see
    snp_pipeline/50b_bmfm_extract_missing_ptids.py."""
    fr = build_frame(clin, snp, fill_missing_snp=False)
    val = fr[fr["split"] == "val"]
    test = fr[fr["split"] == "test"]
    val_b = val[val["clin_present"] & val["snp_present"]]
    test_b = test[test["clin_present"] & test["snp_present"]]
    rows = []
    if len(val_b) < 4 or len(test_b) == 0:
        return rows

    cp_v = val_b[["cp0", "cp1"]].to_numpy(float)
    sp_v = val_b[["sp0", "sp1"]].to_numpy(float)
    cp_t = test_b[["cp0", "cp1"]].to_numpy(float)
    sp_t = test_b[["sp0", "sp1"]].to_numpy(float)
    y_v = val_b["y_true"].to_numpy(int)
    y_t = test_b["y_true"].to_numpy(int)

    # alpha tuned on VAL bACC (threshold fixed at 0.5)
    alpha_star = fit_weighted_avg(cp_v, sp_v, y_v)
    blend_v_a1 = alpha_star * cp_v + (1 - alpha_star) * sp_v
    blend_t_a1 = alpha_star * cp_t + (1 - alpha_star) * sp_t

    # alpha + threshold jointly tuned on VAL Youden's J
    (a2_test_pred, a2_test_proba, a2_val_pred, a2_val_proba,
     a2_alpha, a2_thr, _) = fit_alpha_J(cp_v, sp_v, y_v, cp_t, sp_t)

    # per-class weighted average (w0, w1 jointly tuned on VAL bACC)
    (ca_test_pred, ca_test_proba, ca_val_pred, ca_val_proba,
     ca_w0, ca_w1) = fit_class_avg(cp_v, sp_v, y_v, cp_t, sp_t)

    # method_specs: method_name -> (val_pred, val_proba, test_pred, test_proba, info_dict)
    method_specs = {
        "clinical_only": (cp_v.argmax(1),     cp_v,
                          cp_t.argmax(1),     cp_t,    {"alpha": 1.0}),
        "snp_only":      (sp_v.argmax(1),     sp_v,
                          sp_t.argmax(1),     sp_t,    {"alpha": 0.0}),
        "equal_avg":     ((0.5*cp_v + 0.5*sp_v).argmax(1),  0.5*cp_v + 0.5*sp_v,
                          (0.5*cp_t + 0.5*sp_t).argmax(1),  0.5*cp_t + 0.5*sp_t,
                          {"alpha": 0.5}),
        "best_alpha":    (blend_v_a1.argmax(1), blend_v_a1,
                          blend_t_a1.argmax(1), blend_t_a1,
                          {"alpha": alpha_star}),
        "best_alpha_J":  (a2_val_pred, a2_val_proba,
                          a2_test_pred, a2_test_proba,
                          {"alpha": a2_alpha, "threshold": a2_thr}),
        "class_avg":     (ca_val_pred, ca_val_proba,
                          ca_test_pred, ca_test_proba,
                          {"w0": ca_w0, "w1": ca_w1}),
    }

    # stacked elastic-net (fit on VAL, predict val+test)
    en = fit_stacked_en(cp_v, sp_v, y_v, cp_t, sp_t, seed)
    if en is not None:
        pv_pred, pv_proba, pt_pred, pt_proba, en_info = en
        method_specs["stacked_en"] = (pv_pred, pv_proba, pt_pred, pt_proba,
                                       en_info)

    for method, (vp, vpr, tp, tpr, info) in method_specs.items():
        for split_name, pred, proba, df, y in [
            ("val",  vp, vpr, val_b,  y_v),
            ("test", tp, tpr, test_b, y_t),
        ]:
            r = metric_row(y, pred, proba)
            r.update({
                "clin_variant": clin_v["key"], "snp_variant": snp_v["key"],
                "seed": seed, "method": method, "split": split_name,
                "alpha":         info.get("alpha",     float("nan")),
                "threshold":     info.get("threshold", float("nan")),
                "w0":            info.get("w0",        float("nan")),
                "w1":            info.get("w1",        float("nan")),
                "best_C":        info.get("C",         float("nan")),
                "best_l1_ratio": info.get("l1_ratio",  float("nan")),
            })
            rows.append(r)
            if (method in ("equal_avg", "best_alpha", "best_alpha_J",
                            "class_avg", "stacked_en")
                    and split_name == "test"):
                _stash_preds(pred_rows, clin_v["key"], snp_v["key"],
                             seed, method, df, pred, proba)

    # failure table on test (cohort = full clin-present): one row per
    # misclassified-by-either patient; uses raw per-modality argmax
    cl_pred_t = cp_t.argmax(1)
    sn_pred_t = sp_t.argmax(1)
    for i, (_, row) in enumerate(test_b.reset_index(drop=True).iterrows()):
        y = int(row["y_true"])
        cp_ok = bool(cl_pred_t[i] == y)
        sp_ok = bool(sn_pred_t[i] == y)
        if cp_ok and sp_ok:
            continue
        fail_rows.append({
            "clin_variant": clin_v["key"], "snp_variant": snp_v["key"],
            "seed": seed,
            "Patient_ID": row["Patient_ID"],
            "y_true": CLASS_NAMES[y],
            "clin_pred": CLASS_NAMES[int(cl_pred_t[i])],
            "clin_correct": cp_ok,
            "snp_pred": CLASS_NAMES[int(sn_pred_t[i])],
            "snp_correct": sp_ok,
            "snp_filled": bool(not row["snp_present"]),
            "prob_pCN_clin": float(cp_t[i, 1]),
            "prob_pCN_snp":  float(sp_t[i, 1]),
        })
    return rows


# --------------------------------------------------------------------------- #
# Coverage
# --------------------------------------------------------------------------- #
def compute_coverage(seed, clin_v, snp_v, clin, snp):
    fr = build_frame(clin, snp, fill_missing_snp=False)
    te = fr[fr["split"] == "test"]
    miss = te[te["clin_present"] & ~te["snp_present"]]
    return {
        "clin_variant": clin_v["key"], "snp_variant": snp_v["key"], "seed": seed,
        "n_test": int(len(te)),
        "n_test_both": int((te["clin_present"] & te["snp_present"]).sum()),
        "n_test_clin_only": int(len(miss)),
        "n_test_snp_only":  int((~te["clin_present"] & te["snp_present"]).sum()),
        "n_test_eval": int(te["clin_present"].sum()),
        "missing_snp_ptids": ";".join(miss["Patient_ID"].astype(str).tolist()),
    }


# --------------------------------------------------------------------------- #
# PNG rendering -- mirror 01_fuse_T2.py:render_framing_png
# --------------------------------------------------------------------------- #
def _fmt(mean, std, n=None):
    if mean is None or (isinstance(mean, float) and np.isnan(mean)):
        return "-"
    s = f"{mean:.3f}"
    if std is not None and not (isinstance(std, float) and np.isnan(std)):
        s += f" ± {std:.3f}"
    if n is not None and not (isinstance(n, float) and np.isnan(n)):
        s += f" ({int(round(n))})"
    return s


METHODS = ["clinical_only", "snp_only", "equal_avg",
           "best_alpha", "best_alpha_J", "class_avg", "stacked_en"]
METHOD_SHORT = {"clinical_only": "CL", "snp_only": "SNP",
                "equal_avg": "avg",
                "best_alpha": "best-α₁", "best_alpha_J": "best-α₂",
                "class_avg": "class-avg", "stacked_en": "stack-EN"}
# Fusion methods displayed in section 3 (no single-modality refs there)
FUSION_METHODS = ["equal_avg", "best_alpha", "best_alpha_J",
                  "class_avg", "stacked_en"]
METRIC_KEYS = [("bacc", "bACC"), ("macro_f1", "macro-F1"), ("auc", "AUC")]


def _render_subtable(ax, body, col_labels, col_widths,
                      bold_top_cols=None, left_label_cols=(0,),
                      title=None):
    """Render one ax.table sub-section. bold_top_cols = list of column indices
    where the top-mean numeric value is bolded."""
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=10, weight="bold", loc="left", pad=4)
    tab = ax.table(cellText=body, colLabels=col_labels, loc="upper center",
                   cellLoc="center", colLoc="center", colWidths=col_widths)
    tab.auto_set_font_size(False)
    tab.set_fontsize(9)
    tab.scale(1.0, 1.4)
    n_rows = len(body)
    n_cols = len(col_labels)
    for j in range(n_cols):
        tab[0, j].set_facecolor("#ECECEC")
        tab[0, j].set_text_props(weight="bold")
    for i in range(n_rows + 1):
        for j in left_label_cols:
            tab[i, j].set_text_props(ha="left")
            tab[i, j].PAD = 0.03
    return tab


def _bold_top_in_table(tab, body, col_labels, numeric_offset: int):
    """Bold the top-numeric cell in each numeric column. numeric_offset = the
    column index where numeric data starts (after label columns)."""
    n_rows = len(body)
    n_cols = len(col_labels)
    nums = []
    for i in range(n_rows):
        row_nums = []
        for j in range(numeric_offset, n_cols):
            cell = str(body[i][j])
            try:
                val = float(cell.split("±")[0].strip())
            except (ValueError, AttributeError):
                val = float("nan")
            row_nums.append(val)
        nums.append(row_nums)
    nums = np.array(nums, dtype=float)
    if nums.size == 0:
        return
    for k in range(nums.shape[1]):
        col = nums[:, k]
        if np.all(np.isnan(col)):
            continue
        best = int(np.nanargmax(col))
        tab[best + 1, numeric_offset + k].set_text_props(weight="bold")


def render_combo_table_png(summary: pd.DataFrame, out_path: Path,
                           split: str = "test") -> None:
    """Three-section PNG:
      1. CL-only ablation (per CL variant)
      2. SNP-only ablation (per SNP variant)
      3. Fusion methods (per CL x SNP combo), no constant CL/SNP metric cols.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARN] matplotlib not installed -- skipping PNG.")
        return
    s = summary[summary["split"] == split]
    if s.empty:
        return

    # --- Build section 1: CL-only ablation -----------------------------------
    cl_body, cl_cols = [], ["CL", "bACC", "macro-F1", "AUC"]
    for cl in CLIN_VARIANTS:
        sub = s[(s["clin_variant"] == cl["key"])
                & (s["method"] == "clinical_only")]
        if sub.empty:
            continue
        m = sub.iloc[0]
        cells = [cl["display"]]
        for metric_key, _ in METRIC_KEYS:
            n = m.get("n_mean") if metric_key == "bacc" else None
            cells.append(_fmt(m.get(f"{metric_key}_mean"),
                              m.get(f"{metric_key}_std"), n))
        cl_body.append(cells)

    # --- Build section 2: SNP-only ablation ----------------------------------
    snp_body, snp_cols = [], ["SNP", "bACC", "macro-F1", "AUC"]
    seen = set()
    for sn in SNP_VARIANTS:
        sub = s[(s["snp_variant"] == sn["key"])
                & (s["method"] == "snp_only")]
        if sub.empty or sn["key"] in seen:
            continue
        seen.add(sn["key"])
        m = sub.iloc[0]
        cells = [sn["display"]]
        for metric_key, _ in METRIC_KEYS:
            n = m.get("n_mean") if metric_key == "bacc" else None
            cells.append(_fmt(m.get(f"{metric_key}_mean"),
                              m.get(f"{metric_key}_std"), n))
        snp_body.append(cells)

    # --- Build section 3: fusion --------------------------------------------
    fus_cols = ["CL", "SNP"]
    for method in FUSION_METHODS:
        for _, mlabel in METRIC_KEYS:
            fus_cols.append(f"{METHOD_SHORT[method]} {mlabel}")
    fus_body = []
    for cl in CLIN_VARIANTS:
        for sn in SNP_VARIANTS:
            sub = s[(s["clin_variant"] == cl["key"])
                    & (s["snp_variant"] == sn["key"])]
            if sub.empty:
                continue
            cells = [cl["short"], sn["short"]]
            for method in FUSION_METHODS:
                ms = sub[sub["method"] == method]
                for metric_key, _ in METRIC_KEYS:
                    if ms.empty:
                        cells.append("-")
                    else:
                        m = ms.iloc[0]
                        n = m.get("n_mean") if metric_key == "bacc" else None
                        cells.append(_fmt(m.get(f"{metric_key}_mean"),
                                          m.get(f"{metric_key}_std"), n))
            fus_body.append(cells)

    if not (cl_body and snp_body and fus_body):
        print("  [WARN] PNG sections missing -- aborting render.")
        return

    # --- Figure sizing -------------------------------------------------------
    CHAR_W = 0.10

    def col_widths_for(body, labels, total_chars_target=None):
        chars = []
        for j in range(len(labels)):
            texts = [labels[j]] + [str(b[j]) for b in body]
            chars.append(max(len(t) for t in texts) + 2)
        total = sum(chars)
        widths = [c / total for c in chars]
        return widths, total

    # Shared column widths for sections 1 + 2 so their column separators align.
    # Both sections have the same 4-column structure (label, bACC, F1, AUC).
    ab_chars = []
    for j in range(len(cl_cols)):
        candidates = ([cl_cols[j]] + [str(b[j]) for b in cl_body]
                       + [snp_cols[j]] + [str(b[j]) for b in snp_body])
        ab_chars.append(max(len(t) for t in candidates) + 2)
    ab_total = sum(ab_chars)
    ab_widths = [c / ab_total for c in ab_chars]
    fus_widths, fus_total = col_widths_for(fus_body, fus_cols)
    # Figure width driven by the widest section (Section 3 / fusion).
    fig_w = max(fus_total, ab_total) * CHAR_W
    # Tight per-section heights: just enough for header + rows.
    h_cl   = 0.45 * (len(cl_body)  + 1.5)
    h_snp  = 0.45 * (len(snp_body) + 1.5)
    h_fus  = 0.45 * (len(fus_body) + 1.5)
    h_foot = 1.4

    fig = plt.figure(figsize=(fig_w, h_cl + h_snp + h_fus + h_foot + 0.4))
    gs = fig.add_gridspec(nrows=4, ncols=1,
                          height_ratios=[h_cl, h_snp, h_fus, h_foot],
                          hspace=0.10)
    ax_cl  = fig.add_subplot(gs[0, 0])
    ax_snp = fig.add_subplot(gs[1, 0])
    ax_fus = fig.add_subplot(gs[2, 0])
    ax_ft  = fig.add_subplot(gs[3, 0])
    ax_ft.axis("off")

    tab_cl = _render_subtable(ax_cl, cl_body, cl_cols, ab_widths,
                              title="Section 1 -- CL ONLY ablation (TEST)")
    _bold_top_in_table(tab_cl, cl_body, cl_cols, numeric_offset=1)

    tab_snp = _render_subtable(ax_snp, snp_body, snp_cols, ab_widths,
                               title="Section 2 -- SNP ONLY ablation (TEST)")
    _bold_top_in_table(tab_snp, snp_body, snp_cols, numeric_offset=1)

    tab_fus = _render_subtable(ax_fus, fus_body, fus_cols, fus_widths,
                               left_label_cols=(0, 1),
                               title="Section 3 -- FUSION (TEST)")
    _bold_top_in_table(tab_fus, fus_body, fus_cols, numeric_offset=2)

    fig.suptitle("T1e (sCN vs pCN) late fusion -- TEST metrics\n"
                 "mean ± std across seeds 0,1,2",
                 fontsize=11, y=0.995)

    foot_lines = [
        "Classes: 0=sCN (stable CN), 1=pCN (progressive CN; converted to MCI or AD).",
        "Methods:",
        "  CL          = argmax of clinical-encoder class probs (Section 1).",
        "  SNP         = argmax of SNP-model class probs (Section 2).",
        "  avg         = 0.5*P_CL + 0.5*P_SNP, threshold 0.5.",
        "  best-α₁     = α*P_CL + (1-α)*P_SNP; α grid in [0,1] step 0.05;",
        "                tuned on VAL BALANCED ACCURACY at threshold 0.5.",
        "  best-α₂     = (α, threshold) jointly grid-searched;",
        "                tuned on VAL YOUDEN'S J = sensitivity + specificity - 1;",
        "                threshold grid [0.05, 0.95] step 0.05; both α, t* applied at test.",
        "  class-avg   = per-class weights (w₀, w₁) jointly grid-searched 21x21",
        "                on VAL bACC; p_fused[k] = w_k*P_CL[k] + (1-w_k)*P_SNP[k];",
        "                renormalised before argmax.",
        "  stack-EN    = LogisticRegression(elastic-net) over [P_CL, P_SNP] features;",
        "                (C, l1_ratio) GridSearchCV(cv=3) on VAL.",
        "(n) after bACC = mean #TEST patients across seeds. Bold = top per metric column.",
    ]
    ax_ft.text(0.0, 1.0, "\n".join(foot_lines), transform=ax_ft.transAxes,
               ha="left", va="top", fontsize=7.8, family="monospace",
               linespacing=1.4)
    fig.savefig(out_path, dpi=180, bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)
    print(f"  PNG  -> {out_path}")


def _tuned_params_str(method: str, row: pd.Series) -> str:
    """Compact representation of the tuned parameter(s) for one row."""
    if method == "best_alpha":
        a = row.get("alpha_mean", float("nan"))
        return f"α={a:.2f}" if not pd.isna(a) else "-"
    if method == "best_alpha_J":
        a = row.get("alpha_mean", float("nan"))
        t = row.get("threshold_mean", float("nan"))
        return f"α={a:.2f}, t={t:.2f}"
    if method == "class_avg":
        w0 = row.get("w0_mean", float("nan"))
        w1 = row.get("w1_mean", float("nan"))
        return f"w₀={w0:.2f}, w₁={w1:.2f}"
    if method == "stacked_en":
        C = row.get("best_C_mean", float("nan"))
        l1 = row.get("best_l1_ratio_mean", float("nan"))
        return f"C={C:.2f}, l1={l1:.2f}"
    if method == "equal_avg":
        return "α=0.50"
    if method == "clinical_only":
        return "α=1.00"
    if method == "snp_only":
        return "α=0.00"
    return "-"


def render_summary_table_png(summary: pd.DataFrame, met_per_seed: pd.DataFrame,
                              out_path: Path, top_n: int | None = None,
                              snp_filter: list[str] | None = None,
                              include_cl_refs: bool = False,
                              include_snp_refs: bool = False) -> None:
    """Leaderboard PNG: each (CL, SNP, method) as one row, sorted by TEST bACC.

    Args:
      snp_filter: if set, restrict the fusion-method ranking to these SNP
        variant keys. CL_only / SNP_only references are NOT in this pool.
      include_cl_refs: when True, append one CL-only reference row per CL
        variant (SNP column shown as "(no SNP)"). Included in the final sort.
      include_snp_refs: when True, append one SNP-only reference row per SNP
        variant in `snp_filter` (or all SNP variants if no filter). CL column
        shown as "(no CL)". Included in the final sort.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARN] matplotlib not installed -- skipping summary PNG.")
        return

    # Merge in tuned-param aggregates (alpha, threshold, w0, w1, C, l1_ratio)
    g = ["clin_variant", "snp_variant", "method", "split"]
    pmean = (met_per_seed.groupby(g)[
                ["alpha", "threshold", "w0", "w1", "best_C", "best_l1_ratio"]]
              .mean().round(3).reset_index())
    pmean.columns = (list(g)
                      + ["alpha_mean", "threshold_mean", "w0_mean", "w1_mean",
                          "best_C_mean", "best_l1_ratio_mean"])
    s_all = summary[summary["split"] == "test"].merge(pmean, on=g, how="left")

    # Main ranking pool: fusion methods only (drop the single-modality refs).
    pool = s_all[~s_all["method"].isin(["clinical_only", "snp_only"])]
    if snp_filter is not None:
        pool = pool[pool["snp_variant"].isin(snp_filter)]
    pool = pool.sort_values("bacc_mean", ascending=False).reset_index(drop=True)
    if top_n:
        pool = pool.head(top_n)

    # Optionally append CL-only reference rows (deduplicated per CL variant).
    if include_cl_refs:
        refs = (s_all[s_all["method"] == "clinical_only"]
                 .drop_duplicates(["clin_variant"])
                 .copy())
        refs["snp_variant"] = "(no SNP)"
        pool = pd.concat([pool, refs], ignore_index=True)

    # Optionally append SNP-only reference rows (deduplicated per SNP variant).
    if include_snp_refs:
        srefs = (s_all[s_all["method"] == "snp_only"]
                  .drop_duplicates(["snp_variant"])
                  .copy())
        if snp_filter is not None:
            srefs = srefs[srefs["snp_variant"].isin(snp_filter)]
        srefs["clin_variant"] = "(no CL)"
        pool = pd.concat([pool, srefs], ignore_index=True)

    s = pool.sort_values("bacc_mean", ascending=False).reset_index(drop=True)
    cl_short = {v["key"]: v["display"] for v in CLIN_VARIANTS}
    sn_short = {v["key"]: v["display"] for v in SNP_VARIANTS}

    col_labels = ["Rank", "CL", "SNP", "Method", "bACC",
                  "macro-F1", "AUC", "Tuned params"]
    body, numeric = [], []
    for rank, (_, row) in enumerate(s.iterrows(), start=1):
        method = row["method"]
        cells = [
            str(rank),
            cl_short.get(row["clin_variant"], row["clin_variant"]),
            sn_short.get(row["snp_variant"], row["snp_variant"]),
            METHOD_SHORT.get(method, method),
            _fmt(row.get("bacc_mean"),     row.get("bacc_std")),
            _fmt(row.get("macro_f1_mean"), row.get("macro_f1_std")),
            _fmt(row.get("auc_mean"),      row.get("auc_std")),
            _tuned_params_str(method, row),
        ]
        body.append(cells)
        numeric.append([row.get("bacc_mean", np.nan),
                         row.get("macro_f1_mean", np.nan),
                         row.get("auc_mean", np.nan)])
    if not body:
        return
    numeric = np.array(numeric, dtype=float)

    CHAR_W = 0.105
    col_chars = []
    for j in range(len(col_labels)):
        texts = [col_labels[j]] + [str(b[j]) for b in body]
        col_chars.append(max(len(t) for t in texts) + 2)
    total = sum(col_chars)
    col_widths = [c / total for c in col_chars]
    fig_w = total * CHAR_W
    h_tab = 0.45 * (len(body) + 1)            # ~ row height after scale=1.4
    h_foot = 1.1                                # tight band for footnote
    fig_h = h_tab + h_foot + 0.6                # 0.6 for title + margins

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(2, 1, height_ratios=[h_tab, h_foot], hspace=0.05)
    ax = fig.add_subplot(gs[0, 0])
    ax_ft = fig.add_subplot(gs[1, 0])
    ax.axis("off")
    ax_ft.axis("off")
    tab = ax.table(cellText=body, colLabels=col_labels,
                   loc="upper center", cellLoc="center", colLoc="center",
                   colWidths=col_widths)
    tab.auto_set_font_size(False)
    tab.set_fontsize(9)
    tab.scale(1.0, 1.4)

    # Header style (just bold text, no fill)
    for j in range(len(col_labels)):
        tab[0, j].set_text_props(weight="bold")
    # Left-align label columns (Rank, CL, SNP, Method, Tuned params)
    for i in range(len(body) + 1):
        for j in (0, 1, 2, 3, 7):
            tab[i, j].set_text_props(ha="left")
            tab[i, j].PAD = 0.03
    # Bold the top-1 per numeric column
    for j, label in enumerate(("bACC", "macro-F1", "AUC"), start=4):
        col = numeric[:, j - 4]
        if np.all(np.isnan(col)):
            continue
        best = int(np.nanargmax(col))
        tab[best + 1, j].set_text_props(weight="bold")

    title = "T1e (sCN vs pCN) late fusion -- LEADERBOARD"
    if snp_filter:
        title += " (SNP = " + ", ".join(snp_filter) + ")"
    if top_n:
        title += f" -- top {top_n}"
    refs_added = []
    if include_cl_refs:
        refs_added.append("CL-only")
    if include_snp_refs:
        refs_added.append("SNP-only")
    if refs_added:
        title += " + " + " & ".join(refs_added) + " references"
    ax.set_title(title + "\nsorted by TEST bACC (mean across seeds 0,1,2)",
                 pad=8, fontsize=11)
    foot_lines = [
        "One row per (CL, SNP, method); higher bACC = higher rank.",
        "Tuned params summarise the operating point chosen on VAL:",
        "  best-α₁ = α tuned on VAL bACC (threshold fixed = 0.5).",
        "  best-α₂ = (α, t) jointly tuned on VAL Youden's J = sens + spec − 1.",
        "  class-avg = (w₀, w₁) jointly tuned on VAL bACC; p_k = w_k·P_CL[k] + (1−w_k)·P_SNP[k], renormalised.",
        "  stack-EN = LogisticRegression(elastic-net); (C, l1_ratio) CV-tuned on VAL.",
        "  avg = unweighted 0.5/0.5 blend; CL / SNP = single-modality reference.",
    ]
    ax_ft.text(0.0, 1.0, "\n".join(foot_lines), transform=ax_ft.transAxes,
               ha="left", va="top", fontsize=7.8, family="monospace",
               linespacing=1.4)
    fig.savefig(out_path, dpi=180, bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)
    print(f"  PNG  -> {out_path}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clin-keys", nargs="+", default=None,
                    help=("subset CL variants. default=all: "
                          + ", ".join(v["key"] for v in CLIN_VARIANTS)))
    ap.add_argument("--snp-keys", nargs="+", default=None,
                    help=("subset SNP variants. default=all: "
                          + ", ".join(v["key"] for v in SNP_VARIANTS)))
    ap.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS_DEFAULT))
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR_DEFAULT)
    args = ap.parse_args()

    clin_keys = args.clin_keys or [v["key"] for v in CLIN_VARIANTS]
    snp_keys = args.snp_keys or [v["key"] for v in SNP_VARIANTS]
    clin_run = [v for v in CLIN_VARIANTS if v["key"] in clin_keys]
    snp_run = [v for v in SNP_VARIANTS if v["key"] in snp_keys]
    unknown_cl = set(clin_keys) - {v["key"] for v in CLIN_VARIANTS}
    unknown_sn = set(snp_keys) - {v["key"] for v in SNP_VARIANTS}
    if unknown_cl or unknown_sn:
        print(f"  [ERROR] unknown variant keys: CL={unknown_cl}, SNP={unknown_sn}")
        return 1
    if not clin_run or not snp_run:
        print("  [ERROR] empty CL or SNP set.")
        return 1

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    all_metrics, all_preds, all_fails, coverage = [], [], [], []
    for seed in args.seeds:
        # Cache CL loads per seed (re-used across SNP variants)
        clin_cache: dict[str, pd.DataFrame] = {}
        for clin_v in clin_run:
            try:
                clin_cache[clin_v["key"]] = load_clinical(
                    seed, clin_v["model"], clin_v["strat"])
            except FileNotFoundError as e:
                print(f"  [WARN] {e}; skipping CL '{clin_v['key']}' seed={seed}")
        for clin_v in clin_run:
            if clin_v["key"] not in clin_cache:
                continue
            clin = clin_cache[clin_v["key"]]
            for snp_v in snp_run:
                try:
                    snp = load_snp(seed, snp_v["parquet_tmpl"])
                except FileNotFoundError as e:
                    print(f"  [WARN] {e}; skipping SNP '{snp_v['key']}' seed={seed}")
                    continue
                rows = evaluate_combo(seed, clin_v, snp_v, clin, snp,
                                      all_preds, all_fails)
                all_metrics.extend(rows)
                coverage.append(compute_coverage(seed, clin_v, snp_v, clin, snp))

    if not all_metrics:
        print("  [ERROR] no metrics produced. Check input paths and split alignment.")
        return 1

    met = pd.DataFrame(all_metrics)
    cov = pd.DataFrame(coverage)
    preds = pd.DataFrame(all_preds)
    fails = pd.DataFrame(all_fails)

    # mean ± std over seeds per (clin_variant, snp_variant, method, split)
    agg_cols = ["bacc", "macro_f1", "auc",
                "f1_sCN", "f1_pCN", "recall_sCN", "recall_pCN"]
    g = ["clin_variant", "snp_variant", "method", "split"]
    summary = (met.groupby(g)[agg_cols].agg(["mean", "std"]).round(4))
    summary.columns = [f"{a}_{b}" for a, b in summary.columns]
    summary = summary.reset_index()
    n_mean = (met.groupby(g)["n"].mean().round(1)
              .rename("n_mean").reset_index())
    summary = summary.merge(n_mean, on=g)

    met.to_csv(out_dir / "fusion_metrics_per_seed.csv", index=False)
    summary.to_csv(out_dir / "fusion_metrics.csv", index=False)
    cov.to_csv(out_dir / "coverage.csv", index=False)
    preds.to_csv(out_dir / "fused_predictions.csv", index=False)
    fails.to_csv(out_dir / "failure_table.csv", index=False)
    render_combo_table_png(summary, out_dir / "fusion_table_T1e.png", split="test")
    render_summary_table_png(summary, met,
                              out_dir / "fusion_table_T1e_summary.png",
                              top_n=6,
                              snp_filter=["meta_prs"],
                              include_cl_refs=True,
                              include_snp_refs=True)

    # stdout summary
    print("=" * 78)
    print("  T1e late fusion -- CL x SNP combos (binary sCN vs pCN)")
    print("=" * 78)
    print("\n[coverage] test-set modality overlap per combo / seed:")
    print(cov.to_string(index=False))
    # best-alpha₁ (bACC-tuned) per combo
    a1 = met[met["method"] == "best_alpha"][
        ["clin_variant", "snp_variant", "seed", "split", "alpha"]].drop_duplicates()
    a1_mean = (a1[a1["split"] == "test"]
               .groupby(["clin_variant", "snp_variant"])["alpha"]
               .agg(["mean", "std"]).round(3).reset_index()
               .rename(columns={"mean": "alpha1_mean", "std": "alpha1_std"}))
    print("\n[best-alpha1 (VAL bACC, threshold=0.5) per combo, mean over seeds]:")
    print(a1_mean.to_string(index=False))

    # best-alpha₂ (Youden's J at optimal threshold) per combo
    a2 = met[met["method"] == "best_alpha_J"][
        ["clin_variant", "snp_variant", "seed", "split",
         "alpha", "threshold"]].drop_duplicates()
    a2_mean = (a2[a2["split"] == "test"]
               .groupby(["clin_variant", "snp_variant"])[["alpha", "threshold"]]
               .agg(["mean"]).round(3).reset_index())
    a2_mean.columns = ["clin_variant", "snp_variant", "alpha2_mean", "thr_mean"]
    print("\n[best-alpha2 (VAL Youden's J) per combo, mean over seeds]:")
    print(a2_mean.to_string(index=False))

    # class_avg weights per combo
    ca = met[met["method"] == "class_avg"][
        ["clin_variant", "snp_variant", "seed", "split", "w0", "w1"]].drop_duplicates()
    ca_mean = (ca[ca["split"] == "test"]
               .groupby(["clin_variant", "snp_variant"])[["w0", "w1"]]
               .agg(["mean"]).round(3).reset_index())
    ca_mean.columns = ["clin_variant", "snp_variant", "w0_mean", "w1_mean"]
    print("\n[class_avg (w0=sCN, w1=pCN) per combo, mean over seeds]:")
    print(ca_mean.to_string(index=False))
    show = ["clin_variant", "snp_variant", "method",
            "bacc_mean", "macro_f1_mean", "auc_mean", "n_mean"]
    test_sum = (summary[summary["split"] == "test"][show]
                .sort_values(["clin_variant", "snp_variant", "method"]))
    print("\n[summary] TEST metrics, mean over seeds:")
    print(test_sum.to_string(index=False))

    # snp_only parity check (each SNP variant should give the same snp_only
    # metric regardless of which CL is paired with it)
    print("\n[snp_only parity across CL pairings]:")
    snp_only_test = test_sum[test_sum["method"] == "snp_only"]
    for snp_key in snp_keys:
        sub = snp_only_test[snp_only_test["snp_variant"] == snp_key]
        if sub.empty:
            continue
        baccs = sub["bacc_mean"].round(4).unique()
        flag = " (CONSISTENT)" if len(baccs) == 1 else " (DIFFERS!)"
        print(f"  {snp_key}: snp_only test bACC across CL pairings = "
              f"{list(baccs)}{flag}")

    print(f"\n  wrote -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
