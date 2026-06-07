"""
01_fuse_T2.py
=============
Late-integration (late fusion) for T2_multiclass (CN/MCI/AD) over per-patient
class probabilities from the CLINICAL encoder and an MRI encoder.

Class index convention (enforced): 0=CN, 1=MCI, 2=AD.

Design (integration/T2_classification/notes.txt)
------------------------------------------------
Fuse output class probabilities with an elastic-net multinomial logistic
("stacking") meta-learner AND a tuned weighted-average baseline. Both are
**fit per-seed on the VAL split and reported on TEST** -- the base models were
trained on the TRAIN fold, so their train-fold probabilities are optimistic;
VAL is the first fold unseen by the base models, so the meta-learner is fit
there and TEST is held out for reporting.

Two TARGET FRAMINGS, reported side by side (baseline MRI is sparse: in the T2
test split bl=6 vs m12=48 scans, because ADNI MRI is overwhelmingly m12+):
  - "BL"  : predict the BASELINE diagnosis from CL_bl + MRI_bl probabilities.
            Target y comes from the clinical (baseline) label. MRI_bl is rare,
            so this exercises the missingness handling below.
  - "M12" : predict the 1-YEAR (m12) diagnosis from CL_bl + MRI_m12. Target y
            comes from the MRI m12 per-session label (Label_visit_diag at m12).
            CL_bl acts as a baseline prior. Viable now (~48 test patients).

Missingness (BL framing -- MRI_bl mostly absent), compared head to head:
  (A) complete-case + ZERO-FALLBACK: fit the EN on patients with BOTH modalities
      present; at inference a patient missing MRI gets its MRI probs set to 0, so
      only the clinical terms + intercept drive the prediction ("clinical only
      goes ahead"). Single fitted model, graceful fallback.
  (B) IMPUTE + INDICATOR: fill a missing modality's probs with the TRAIN-fold
      class prior and add a binary "modality present" feature; fit on all patients.

Single-modality references (clinical-only, MRI-only; raw argmax) are reported
for every framing so the fusion lift is explicit.

BLOCKED extension: the full pooled EN over [CL_bl, CL_m12, MRI_bl, MRI_m12]
needs a CLINICAL m12 run (does not exist yet -- raw text at verbose/by_visit/
visit_m12.csv, but no m12 splits / encoder run). Once that lands, add CL_m12 as
another feature block here; the harness is already timepoint-aware.

Inputs (local D:)
-----------------
Clinical per-seed parquet:
  <CLIN_ROOT>/<CLIN_MODEL>/T2_multiclass/seed_{s}/<CLIN_STRAT>/embeddings/embeddings.parquet
  cols: Patient_ID, split, in_task_cohort, label(0/1/2), pred, prob_0, prob_1, prob_2
MRI per-seed extraction CSV (05_extract_* schema), model-agnostic via template:
  cols: Patient_ID, adni_viscode, split, label(-1 if ineligible), pred,
        prob_class_0, prob_class_1, prob_class_2

Usage
-----
  # BrainDINO frozen-none (available locally now -- good for a smoke test):
  python integration/T2_classification/01_fuse_T2.py \
      --mri-name braindino_frozen_none \
      --mri-csv-template "D:/ADNI_BIDS_project/derivatives/mri_embeddings/T2_multiclass/braindino_frozen_none_cached/embeddings_seed_{seed}.csv"

  # BrainMVP ft stochastic (the notes' T2 MRI model; after the CSD3 extraction + scp):
  python integration/T2_classification/01_fuse_T2.py \
      --mri-name brainmvp_stochastic \
      --mri-csv-template "D:/ADNI_BIDS_project/derivatives/mri_embeddings_brainmvp/T2_multiclass/aug_stochastic/seed_{seed}/embeddings_seed_{seed}.csv"
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import (balanced_accuracy_score, f1_score,
                             roc_auc_score, confusion_matrix,
                             precision_recall_fscore_support)

warnings.filterwarnings("ignore")

CLASS_NAMES = ["CN", "MCI", "AD"]      # default 3-class (T2); index 0,1,2
N_CLASSES = 3
SEEDS_DEFAULT = (0, 1, 2)

# Class-count-derived column name lists (rebuilt by _set_classes for binary tasks
# like T1d). CP_/MP_ = renamed clinical/MRI prob cols; PROB_/MPROB_ = raw parquet/csv
# prob col names. METRIC_COLS = the rendered metric columns (3 headline + per-class F1).
def _set_classes(names):
    global CLASS_NAMES, N_CLASSES, CP_COLS, MP_COLS, PROB_COLS, MPROB_COLS, METRIC_COLS
    CLASS_NAMES = list(names)
    N_CLASSES = len(CLASS_NAMES)
    CP_COLS = [f"cp{i}" for i in range(N_CLASSES)]
    MP_COLS = [f"mp{i}" for i in range(N_CLASSES)]
    PROB_COLS = [f"prob_{i}" for i in range(N_CLASSES)]
    MPROB_COLS = [f"prob_class_{i}" for i in range(N_CLASSES)]
    METRIC_COLS = ([("bacc", "bACC", True), ("macro_f1", "macro-F1", False),
                    ("macro_auc", "macro-AUC", False)]
                   + [(f"f1_{n}", f"F1 {n}", False) for n in CLASS_NAMES])

_set_classes(CLASS_NAMES)

# --------------------------------------------------------------------------- #
# Default paths (LOCAL D:)
# --------------------------------------------------------------------------- #
CLIN_ROOT = Path(r"D:/ADNI_BIDS_project/derivatives/encoder_outputs_no_cdr_post_exclusion")
CLIN_MODEL = "BioClinical-ModernBERT-large"     # "BioBERT-L ft"
CLIN_STRAT = "full_ft"
OUT_DIR_DEFAULT = Path(__file__).resolve().parent / "outputs"


# --------------------------------------------------------------------------- #
# Loaders
# --------------------------------------------------------------------------- #
def load_clinical(seed: int, clin_root: Path, model: str, strat: str,
                  task: str = "T2_multiclass") -> pd.DataFrame:
    """Clinical probs -> Patient_ID, split, y_clin, cp0..cp2. `task` selects the
    parquet subdir: T2_multiclass (baseline) or T2_m12_multiclass (m12 concurrent)."""
    p = (clin_root / model / task / f"seed_{seed}" / strat /
         "embeddings" / "embeddings.parquet")
    if not p.exists():
        raise FileNotFoundError(f"clinical parquet missing: {p}")
    d = pd.read_parquet(p, columns=["Patient_ID", "split", "in_task_cohort",
                                    "label"] + PROB_COLS)
    d = d[d["in_task_cohort"] & d["label"].notna()].copy()
    d["y_clin"] = d["label"].astype(int)
    d = d.rename(columns=dict(zip(PROB_COLS, CP_COLS)))
    return d[["Patient_ID", "split", "y_clin"] + CP_COLS]


def load_mri(seed: int, template: str, timepoint: str) -> pd.DataFrame:
    """MRI probs at one session -> Patient_ID, split_mri, y_mri, mp0..mp2."""
    p = Path(template.format(seed=seed))
    if not p.exists():
        raise FileNotFoundError(f"MRI csv missing: {p}")
    d = pd.read_csv(p, low_memory=False)
    d = d[d["adni_viscode"] == timepoint].copy()
    # one row per patient at this session (master is unique by patient-scan, but guard)
    d = d.drop_duplicates(subset=["Patient_ID"], keep="first")
    rename = dict(zip(MPROB_COLS, MP_COLS))
    rename.update({"split": "split_mri", "label": "y_mri"})
    d = d.rename(columns=rename)
    d["y_mri"] = d["y_mri"].astype(int)
    return d[["Patient_ID", "split_mri", "y_mri"] + MP_COLS]


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #
def safe_macro_auc(y_true: np.ndarray, probs: np.ndarray) -> float:
    """AUC; NaN if a class is absent or it fails. Binary -> ROC AUC of the
    positive-class prob; multiclass -> macro one-vs-rest."""
    try:
        if len(np.unique(y_true)) < 2:
            return float("nan")
        if N_CLASSES == 2:
            return float(roc_auc_score(y_true, probs[:, 1]))
        return float(roc_auc_score(y_true, probs, multi_class="ovr",
                                   average="macro", labels=list(range(N_CLASSES))))
    except Exception:
        return float("nan")


def youden_j(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Macro Youden's J = mean_c(sensitivity_c + specificity_c - 1). For binary this
    is the usual sens+spec-1; for K classes it macro-averages the one-vs-rest J."""
    cm = confusion_matrix(y_true, y_pred, labels=list(range(N_CLASSES)))
    total = cm.sum()
    js = []
    for c in range(N_CLASSES):
        tp = cm[c, c]
        fn = cm[c, :].sum() - tp
        fp = cm[:, c].sum() - tp
        tn = total - tp - fn - fp
        sens = tp / (tp + fn) if (tp + fn) else 0.0
        spec = tn / (tn + fp) if (tn + fp) else 0.0
        js.append(sens + spec - 1.0)
    return float(np.mean(js))


def metric_row(y_true: np.ndarray, y_pred: np.ndarray,
               probs: np.ndarray | None) -> dict:
    bacc = balanced_accuracy_score(y_true, y_pred)
    macro_f1 = f1_score(y_true, y_pred, average="macro",
                        labels=list(range(N_CLASSES)), zero_division=0)
    pr, rc, f1, _ = precision_recall_fscore_support(
        y_true, y_pred, labels=list(range(N_CLASSES)), zero_division=0)
    out = {"n": int(len(y_true)), "bacc": float(bacc), "macro_f1": float(macro_f1)}
    out["macro_auc"] = safe_macro_auc(y_true, probs) if probs is not None else float("nan")
    for c, name in enumerate(CLASS_NAMES):
        out[f"f1_{name}"] = float(f1[c])
        out[f"recall_{name}"] = float(rc[c])
    return out


# --------------------------------------------------------------------------- #
# Meta-learners (fit on VAL, predict TEST)
# --------------------------------------------------------------------------- #
def fit_en(Xv, yv, Xt, seed):
    """Elastic-net multinomial logistic; small grid via CV on VAL.
    Returns (pred, full-proba) or None if VAL cannot support a fit (the
    both-present VAL subset is often tiny/single-class when baseline MRI is
    sparse -- skip rather than crash)."""
    if len(np.unique(yv)) < 2:
        return None
    # Enough rows per class for CV?
    n_splits = min(3, int(np.min(np.bincount(yv, minlength=N_CLASSES)[np.unique(yv)])))
    base = LogisticRegression(penalty="elasticnet", solver="saga",
                              max_iter=5000, random_state=seed)
    if n_splits < 2:
        # Degenerate VAL: fall back to a fixed mild EN, no CV.
        clf = LogisticRegression(penalty="elasticnet", solver="saga", C=1.0,
                                 l1_ratio=0.5, max_iter=5000, random_state=seed)
        clf.fit(Xv, yv)
    else:
        grid = {"C": [0.1, 1.0, 10.0], "l1_ratio": [0.1, 0.5, 0.9]}
        cv = StratifiedKFold(n_splits=int(n_splits), shuffle=True, random_state=seed)
        clf = GridSearchCV(base, grid, scoring="balanced_accuracy", cv=cv, n_jobs=-1)
        clf.fit(Xv, yv)
    proba = clf.predict_proba(Xt)
    # Align proba columns to full class set (LogReg uses clf.classes_).
    classes = clf.best_estimator_.classes_ if hasattr(clf, "best_estimator_") else clf.classes_
    full = np.zeros((proba.shape[0], N_CLASSES))
    for j, c in enumerate(classes):
        full[:, int(c)] = proba[:, j]
    return full.argmax(1), full


def fit_weighted_avg(cl_v, mr_v, yv, cl_t, mr_t, objective="bacc"):
    """Blend w*clinical + (1-w)*mri per-class prob; pick w on VAL by `objective`
    ('bacc' = balanced accuracy | 'youden' = macro Youden's J). Both present."""
    score = youden_j if objective == "youden" else balanced_accuracy_score
    best_w, best_s = 0.5, -2.0
    for w in np.linspace(0, 1, 21):
        s = score(yv, (w * cl_v + (1 - w) * mr_v).argmax(1))
        if s > best_s:
            best_s, best_w = s, w
    blend_t = best_w * cl_t + (1 - best_w) * mr_t
    return blend_t.argmax(1), blend_t, best_w


def fit_lr(Xv, yv, Xt, seed):
    """Plain logistic regression (default L2, no grid) over stacked [CL, MRI] probs.
    Returns (pred, full-proba) or None if VAL has <2 classes."""
    if len(np.unique(yv)) < 2:
        return None
    clf = LogisticRegression(max_iter=5000, random_state=seed)
    clf.fit(Xv, yv)
    proba = clf.predict_proba(Xt)
    full = np.zeros((proba.shape[0], N_CLASSES))
    for j, c in enumerate(clf.classes_):
        full[:, int(c)] = proba[:, j]
    return full.argmax(1), full


# --------------------------------------------------------------------------- #
# Detector-augmented fusion: sharpen CN with MRI T1 (CN-vs-rest) and AD with MRI
# T1b (rest-vs-AD); leave MCI to clinical. p(MCI)=CL_MCI, p(CN)/p(AD) blended,
# renormalized. The MRI 3-class (T2) is weak only on the MCI middle class, so we
# only borrow MRI's strong binary contrasts.
# --------------------------------------------------------------------------- #
def load_mri_detector(seed, template, timepoint, prob_col, out_name):
    p = Path(template.format(seed=seed))
    if not p.exists():
        raise FileNotFoundError(f"MRI detector csv missing: {p}")
    d = pd.read_csv(p, low_memory=False)
    d = d[d["adni_viscode"] == timepoint].drop_duplicates("Patient_ID")
    return d[["Patient_ID", prob_col]].rename(columns={prob_col: out_name})


def blend_detectors(cp, mri_cn, mri_ad, w_cn, w_ad):
    """cp:(N,3) clinical [CN,MCI,AD]; mri_cn/mri_ad:(N,) P(CN)/P(AD), may be NaN.
    Missing detector -> fall back to the clinical class prob (no augmentation).
    Returns renormalized (N,3) [CN,MCI,AD]."""
    pcn = np.where(np.isnan(mri_cn), cp[:, 0], w_cn * cp[:, 0] + (1 - w_cn) * mri_cn)
    pad = np.where(np.isnan(mri_ad), cp[:, 2], w_ad * cp[:, 2] + (1 - w_ad) * mri_ad)
    P = np.stack([pcn, cp[:, 1], pad], axis=1)
    P = np.clip(P, 1e-9, None)
    return P / P.sum(axis=1, keepdims=True)


def fit_detector_weights(cp, mri_cn, mri_ad, y):
    """Grid-search (w_cn, w_ad) in [0,1] maximizing VAL balanced accuracy."""
    grid = np.linspace(0, 1, 11)
    best, best_b = (0.5, 0.5), -1.0
    for wcn in grid:
        for wad in grid:
            b = balanced_accuracy_score(y, blend_detectors(cp, mri_cn, mri_ad, wcn, wad).argmax(1))
            if b > best_b:
                best_b, best = b, (wcn, wad)
    return best


def evaluate_detectors(seed, clin, t1, t1b, pred_rows):
    """CL_m12 3-class augmented by the AVAILABLE MRI binary detectors (T1->CN,
    T1b->AD). Either t1 or t1b may be None to use only one. Returns (rows, frame, fail_df)."""
    framing = "M12_DET"
    fr = clin.copy()
    fr["y_true"] = fr["y_clin"]
    if t1 is not None:
        fr = fr.merge(t1, on="Patient_ID", how="left")
    if "mri_cn" not in fr.columns:
        fr["mri_cn"] = np.nan
    if t1b is not None:
        fr = fr.merge(t1b, on="Patient_ID", how="left")
    if "mri_ad" not in fr.columns:
        fr["mri_ad"] = np.nan
    fr["clin_present"] = True
    det_cols = (["mri_cn"] if t1 is not None else []) + (["mri_ad"] if t1b is not None else [])
    fr["mri_present"] = fr[det_cols].notna().all(axis=1) if det_cols else False

    val, test = fr[fr["split"] == "val"], fr[fr["split"] == "test"]
    cpv = val[["cp0", "cp1", "cp2"]].to_numpy(float)
    cpt = test[["cp0", "cp1", "cp2"]].to_numpy(float)
    yv, yt = val["y_true"].to_numpy(int), test["y_true"].to_numpy(int)
    rows = []

    # clinical-only reference (3-class)
    r = metric_row(yt, cpt.argmax(1), cpt)
    r["val_bacc"] = _vbacc(yv, cpv.argmax(1))
    r.update(framing=framing, seed=seed, variant="clinical_only", missing="-")
    rows.append(r)

    # ---- MRI-only detector references: how each BINARY detector does on its OWN task,
    #      on this m12 test cohort. T1 = CN-vs-rest (target CN), T1b = rest-vs-AD (target AD).
    #      bACC/AUC are BINARY; the per-class F1 goes in the detector's target-class column. ----
    def _det_ref(prob_t, y3_t, prob_v, y3_v, target_idx):
        yb_t = (np.asarray(y3_t, int) == target_idx).astype(int)
        pb_t = (prob_t > 0.5).astype(int)
        rr = {"n": int(len(yb_t)),
              "bacc": float(balanced_accuracy_score(yb_t, pb_t)) if len(np.unique(yb_t)) > 1 else float("nan"),
              "macro_f1": float(f1_score(yb_t, pb_t, average="macro", zero_division=0)),
              "macro_auc": float(roc_auc_score(yb_t, prob_t)) if len(np.unique(yb_t)) > 1 else float("nan")}
        for nm in CLASS_NAMES:
            rr[f"f1_{nm}"] = float("nan"); rr[f"recall_{nm}"] = float("nan")
        rr[f"f1_{CLASS_NAMES[target_idx]}"] = float(f1_score(yb_t, pb_t, pos_label=1, zero_division=0))
        yb_v = (np.asarray(y3_v, int) == target_idx).astype(int)
        rr["val_bacc"] = _vbacc(yb_v, (prob_v > 0.5).astype(int))
        return rr

    if t1 is not None:
        te, ve = test[test["mri_cn"].notna()], val[val["mri_cn"].notna()]
        if len(te):
            r = _det_ref(te["mri_cn"].to_numpy(float), te["y_true"].to_numpy(int),
                         ve["mri_cn"].to_numpy(float), ve["y_true"].to_numpy(int), 0)
            r.update(framing=framing, seed=seed, variant="mri_only_T1", missing="-")
            rows.append(r)
    if t1b is not None:
        te, ve = test[test["mri_ad"].notna()], val[val["mri_ad"].notna()]
        if len(te):
            r = _det_ref(te["mri_ad"].to_numpy(float), te["y_true"].to_numpy(int),
                         ve["mri_ad"].to_numpy(float), ve["y_true"].to_numpy(int), 2)
            r.update(framing=framing, seed=seed, variant="mri_only_T1b", missing="-")
            rows.append(r)

    # detector-augmented structured blend (weights tuned on VAL, fallback on missing)
    wcn, wad = fit_detector_weights(cpv, val["mri_cn"].to_numpy(float),
                                    val["mri_ad"].to_numpy(float), yv)
    Pt = blend_detectors(cpt, test["mri_cn"].to_numpy(float),
                         test["mri_ad"].to_numpy(float), wcn, wad)
    r = metric_row(yt, Pt.argmax(1), Pt)
    r["val_bacc"] = _vbacc(yv, blend_detectors(cpv, val["mri_cn"].to_numpy(float),
                                               val["mri_ad"].to_numpy(float), wcn, wad).argmax(1))
    r.update(framing=framing, seed=seed, variant="detector_structured", missing="-",
             note=f"w_cn={wcn:.1f} w_ad={wad:.1f}")
    rows.append(r)
    _stash_preds(pred_rows, framing, seed, "detector_structured", "-", test, Pt.argmax(1), Pt)

    # data-driven EN over [CL3 + provided detector(s)] (complete-case)
    feats = ["cp0", "cp1", "cp2"] + det_cols
    vb, tb = val[val["mri_present"]], test[test["mri_present"]]
    if det_cols and len(vb) >= 4 and len(tb):
        en = fit_en(vb[feats].to_numpy(float), vb["y_true"].to_numpy(int),
                    tb[feats].to_numpy(float), seed)
        if en is not None:
            ep, epr = en
            r = metric_row(tb["y_true"].to_numpy(int), ep, epr)
            r.update(framing=framing, seed=seed, variant="elastic_net", missing="detectors")
            rows.append(r)
            _stash_preds(pred_rows, framing, seed, "elastic_net", "detectors", tb, ep, epr)

    # ---- per-patient failure table: expose each detector's BINARY call so the
    #      class-specific sensitivity (T1->CN, T1b->AD) is visible alongside the
    #      clinical 3-class and the fused (detector-augmented) prediction. ----
    fail = []
    test_r = test.reset_index(drop=True)
    for i in range(len(test_r)):
        y = int(test_r.loc[i, "y_true"])
        cl, fu = int(cpt[i].argmax()), int(Pt[i].argmax())
        rec = {"framing": framing, "seed": seed,
               "Patient_ID": test_r.loc[i, "Patient_ID"], "y_true": CLASS_NAMES[y],
               "clinical_pred": CLASS_NAMES[cl], "clinical_correct": cl == y,
               "fused_pred": CLASS_NAMES[fu], "fused_correct": fu == y}
        # T1 detector (CN-vs-rest): include only if provided. "correct" = CN/not-CN call matches CN membership.
        if t1 is not None:
            mc = test_r.loc[i, "mri_cn"]
            if pd.isna(mc):
                rec["t1_cn_pred"], rec["t1_cn_correct"] = np.nan, np.nan
            else:
                call = mc > 0.5
                rec["t1_cn_pred"] = "CN" if call else "not-CN"
                rec["t1_cn_correct"] = bool(call == (y == 0))
        # T1b detector (rest-vs-AD): include only if provided.
        if t1b is not None:
            ma = test_r.loc[i, "mri_ad"]
            if pd.isna(ma):
                rec["t1b_ad_pred"], rec["t1b_ad_correct"] = np.nan, np.nan
            else:
                call = ma > 0.5
                rec["t1b_ad_pred"] = "AD" if call else "not-AD"
                rec["t1b_ad_correct"] = bool(call == (y == 2))
        cc, fc = cl == y, fu == y
        rec["agreement"] = ("both_correct" if cc and fc else
                            "clinical_better" if cc and not fc else
                            "fused_better" if fc and not cc else "both_wrong")
        fail.append(rec)
    fail_df = pd.DataFrame(fail)
    return rows, fr, fail_df


# --------------------------------------------------------------------------- #
# Framing builder + evaluation
# --------------------------------------------------------------------------- #
def build_frame(framing: str, seed: int, clin: pd.DataFrame, mri: pd.DataFrame):
    """Join clinical(baseline) + MRI(at framing timepoint). Returns a per-patient frame
    with y_true, clinical present/probs, mri present/probs, and an authoritative split."""
    if framing in ("BL", "M12X"):
        # Clinical-as-base; target = clinical dx (baseline for BL, concurrent m12
        # for M12X when clinical is loaded from the T2_m12_multiclass parquet).
        # MRI is joined at this framing's timepoint (bl for BL, m12 for M12X).
        base = clin.copy()
        base["y_true"] = base["y_clin"]
        base["split"] = base["split"]
        m = base.merge(mri, on="Patient_ID", how="left")
    elif framing == "M12":
        # Target = m12 dx (MRI m12 label). Cohort = patients with an m12 MRI scan.
        base = mri[mri["y_mri"] >= 0].copy()
        base["y_true"] = base["y_mri"]
        m = base.merge(clin, on="Patient_ID", how="left")
        # split: prefer clinical's (consistent across modalities); else MRI's.
        m["split"] = m["split"].fillna(m["split_mri"])
    else:
        raise ValueError(framing)

    m["clin_present"] = m[CP_COLS].notna().all(axis=1)
    m["mri_present"] = m[MP_COLS].notna().all(axis=1)
    return m


def class_prior_from_train(clin: pd.DataFrame) -> np.ndarray:
    """Train-fold class prior (leakage-safe impute value), from clinical TRAIN labels."""
    tr = clin[clin["split"] == "train"]["y_clin"].astype(int)
    counts = np.bincount(tr, minlength=N_CLASSES).astype(float)
    return counts / counts.sum() if counts.sum() else np.ones(N_CLASSES) / N_CLASSES


def _vbacc(yv, pv):
    yv = np.asarray(yv, int)
    return (float(balanced_accuracy_score(yv, np.asarray(pv, int)))
            if len(yv) and len(np.unique(yv)) > 1 else float("nan"))


def evaluate_framing(framing, seed, clin, mri, prior, pred_rows):
    """Return list[metric dicts] for every variant in this framing/seed. Each row carries
    `val_bacc` (in-sample VAL balanced accuracy of the fitted strategy)."""
    fr = build_frame(framing, seed, clin, mri)
    rows = []

    def _probs(df, cols):
        return df[cols].to_numpy(dtype=float)

    def _fit2(fit_fn, Xv, yv_, Xt):
        """Fit on (Xv,yv); predict VAL+TEST in one fit -> (pred_v, pred_t, proba_t) or None."""
        res = fit_fn(Xv, yv_, np.vstack([Xv, Xt]), seed)
        if res is None:
            return None
        pa, proba = res
        nv = len(Xv)
        return pa[:nv], pa[nv:], proba[nv:]

    val = fr[fr["split"] == "val"]
    test = fr[fr["split"] == "test"]

    # ---- single-modality references (raw argmax on each modality's own cohort) ----
    for name, cols, presence in [("clinical_only", CP_COLS, "clin_present"),
                                  ("mri_only", MP_COLS, "mri_present")]:
        te = test[test[presence]]
        ve = val[val[presence]]
        if len(te):
            P = _probs(te, cols)
            r = metric_row(te["y_true"].to_numpy(int), P.argmax(1), P)
            r["val_bacc"] = _vbacc(ve["y_true"].to_numpy(int),
                                   _probs(ve, cols).argmax(1)) if len(ve) else float("nan")
            r.update(framing=framing, seed=seed, variant=name, missing="-")
            rows.append(r)
            _stash_preds(pred_rows, framing, seed, name, "-", te, P.argmax(1), P)

    # ---- complete-case (both present): weighted-avg(+J), LR, EN ----
    vb = val[val["clin_present"] & val["mri_present"]]
    tb = test[test["clin_present"] & test["mri_present"]]
    if len(vb) >= 4 and len(tb):
        Xv = _probs(vb, CP_COLS + MP_COLS)
        Xt = _probs(tb, CP_COLS + MP_COLS)
        yv, yt = vb["y_true"].to_numpy(int), tb["y_true"].to_numpy(int)
        clv, mrv = _probs(vb, CP_COLS), _probs(vb, MP_COLS)
        clt, mrt = _probs(tb, CP_COLS), _probs(tb, MP_COLS)
        for objective, vname in [("bacc", "weighted_avg"), ("youden", "weighted_avg_J")]:
            wpred, wproba, w = fit_weighted_avg(clv, mrv, yv, clt, mrt, objective=objective)
            r = metric_row(yt, wpred, wproba)
            r["val_bacc"] = _vbacc(yv, (w * clv + (1 - w) * mrv).argmax(1))
            r.update(framing=framing, seed=seed, variant=vname,
                     missing="complete_case", note=f"w_clin={w:.2f}")
            rows.append(r)
            _stash_preds(pred_rows, framing, seed, vname, "complete_case", tb, wpred, wproba)
        for fit, vname in [(fit_lr, "lr"), (fit_en, "elastic_net")]:
            out = _fit2(fit, Xv, yv, Xt)
            if out is None:
                continue
            pv, pt, prt = out
            r = metric_row(yt, pt, prt)
            r["val_bacc"] = _vbacc(yv, pv)
            r.update(framing=framing, seed=seed, variant=vname, missing="complete_case")
            rows.append(r)
            _stash_preds(pred_rows, framing, seed, vname, "complete_case", tb, pt, prt)
            if vname == "elastic_net":
                # ---- (A) zero-fallback: same fit, applied to clinical-present (missing MRI->0) ----
                ta = test[test["clin_present"]].copy()
                va = val[val["clin_present"]].copy()
                Xa = ta[CP_COLS + MP_COLS].copy(); Xa[MP_COLS] = Xa[MP_COLS].fillna(0.0)
                Xav = va[CP_COLS + MP_COLS].copy(); Xav[MP_COLS] = Xav[MP_COLS].fillna(0.0)
                res = fit_en(Xv, yv, np.vstack([Xav.to_numpy(float), Xa.to_numpy(float)]), seed)
                if res is not None:
                    pa_, pr_ = res; nv = len(Xav)
                    r = metric_row(ta["y_true"].to_numpy(int), pa_[nv:], pr_[nv:])
                    r["val_bacc"] = _vbacc(va["y_true"].to_numpy(int), pa_[:nv])
                    r.update(framing=framing, seed=seed, variant="elastic_net",
                             missing="zero_fallback")
                    rows.append(r)
                    _stash_preds(pred_rows, framing, seed, "elastic_net", "zero_fallback",
                                 ta, pa_[nv:], pr_[nv:])

    # ---- (B) impute + indicator: fill missing modality with train prior + flag ----
    def _impute_indicator(df):
        X = df[CP_COLS + MP_COLS].copy()
        X[CP_COLS] = X[CP_COLS].fillna(pd.Series(prior, index=CP_COLS))
        X[MP_COLS] = X[MP_COLS].fillna(pd.Series(prior, index=MP_COLS))
        X["clin_present"] = df["clin_present"].astype(float).values
        X["mri_present"] = df["mri_present"].astype(float).values
        return X.to_numpy(float)

    vi = val[val["clin_present"] | val["mri_present"]]
    ti = test[test["clin_present"] | test["mri_present"]]
    if len(vi) >= 4 and len(ti):
        out = _fit2(fit_en, _impute_indicator(vi), vi["y_true"].to_numpy(int),
                    _impute_indicator(ti))
        if out is not None:
            pv, pt, prt = out
            r = metric_row(ti["y_true"].to_numpy(int), pt, prt)
            r["val_bacc"] = _vbacc(vi["y_true"].to_numpy(int), pv)
            r.update(framing=framing, seed=seed, variant="elastic_net",
                     missing="impute_indicator")
            rows.append(r)
            _stash_preds(pred_rows, framing, seed, "elastic_net", "impute_indicator",
                         ti, pt, prt)

    return rows, fr


def _apply_en(Xv, yv, Xt, seed):
    """Thin wrapper: fit EN on (Xv,yv), predict Xt. Returns (pred, full proba)."""
    return fit_en(Xv, yv, Xt, seed)


def _stash_preds(pred_rows, framing, seed, variant, missing, test_df, pred, proba):
    for i, (_, row) in enumerate(test_df.reset_index(drop=True).iterrows()):
        rec = {
            "framing": framing, "seed": seed, "variant": variant, "missing": missing,
            "Patient_ID": row["Patient_ID"], "y_true": int(row["y_true"]),
            "pred": int(pred[i]),
            "clin_present": bool(row["clin_present"]),
            "mri_present": bool(row["mri_present"]),
        }
        for c, name in enumerate(CLASS_NAMES):
            rec[f"prob_{name}"] = float(proba[i, c])
        pred_rows.append(rec)


# --------------------------------------------------------------------------- #
# Failure-mode table (per-patient modality agreement on TEST)
# --------------------------------------------------------------------------- #
def failure_table(framing, seed, fr):
    """Per TEST patient: each modality's argmax correctness, to expose complementarity."""
    test = fr[fr["split"] == "test"].copy()
    out = []
    for _, r in test.iterrows():
        y = int(r["y_true"])
        rec = {"framing": framing, "seed": seed, "Patient_ID": r["Patient_ID"],
               "y_true": CLASS_NAMES[y]}
        for mod, cols, pres in [("clinical", CP_COLS, "clin_present"),
                                ("mri", MP_COLS, "mri_present")]:
            if r[pres]:
                p = int(np.argmax([r[c] for c in cols]))
                rec[f"{mod}_pred"] = CLASS_NAMES[p]
                rec[f"{mod}_correct"] = (p == y)
            else:
                rec[f"{mod}_pred"] = "absent"
                rec[f"{mod}_correct"] = None
        out.append(rec)
    df = pd.DataFrame(out)
    # keep patients where at least one modality is wrong (the interesting rows)
    if len(df):
        wrong = (df["clinical_correct"] == False) | (df["mri_correct"] == False)
        df = df[wrong]
    return df


def failure_stats(framing, seed, fr):
    """REAL per-modality performance on the FULL test cohort (NOT the wrong-enriched
    failure subset that `failure_table` returns). Lets the failure-table PNG state the
    true accuracy/bACC + how many both-correct rows were hidden, so the filtered green/
    red view is not mistaken for the real accuracy."""
    test = fr[fr["split"] == "test"]
    n_test = int(len(test))

    def _argmax_ok(row, cols):
        return int(np.argmax([row[c] for c in cols])) == int(row["y_true"])

    both_pres = test[test["clin_present"] & test["mri_present"]]
    n_both_present = int(len(both_pres))
    n_both_correct = int(sum(_argmax_ok(r, CP_COLS) and _argmax_ok(r, MP_COLS)
                             for _, r in both_pres.iterrows()))
    rows = []
    for mod, cols, pres in [("clinical", CP_COLS, "clin_present"),
                            ("mri", MP_COLS, "mri_present")]:
        sub = test[test[pres]]
        if not len(sub):
            continue
        y = sub["y_true"].to_numpy(int)
        pred = sub[cols].to_numpy(float).argmax(1)
        rows.append({
            "framing": framing, "seed": seed, "modality": mod,
            "n": int(len(sub)), "n_correct": int((y == pred).sum()),
            "raw_acc": float((y == pred).mean()),
            "balanced_acc": float(balanced_accuracy_score(y, pred))
                            if len(np.unique(y)) > 1 else float("nan"),
            "n_test": n_test, "n_both_present": n_both_present,
            "n_both_correct_hidden": n_both_correct,
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Rendered table (matplotlib) -- mirrors mri_pipeline/06_render_cross_model_table
# --------------------------------------------------------------------------- #
# Row display order + labels (only rows actually present are shown).
ROW_ORDER = [
    ("clinical_only",  "-",                "clinical-only (ref)"),
    ("mri_only",       "-",                "MRI-only (ref)"),
    ("weighted_avg",   "complete_case",    "weighted-avg (bACC-tuned)"),
    ("weighted_avg_J", "complete_case",    "weighted-avg (Youden's J)"),
    ("lr",             "complete_case",    "logistic regression"),
    ("elastic_net",    "complete_case",    "elastic-net (both present)"),
    ("elastic_net",    "zero_fallback",    "elastic-net (zero-fallback)"),
    ("elastic_net",    "impute_indicator", "elastic-net (impute+indicator)"),
    ("mri_only_T1",    "-",                "MRI-only: T1 detector (CN vs rest, binary)"),
    ("mri_only_T1b",   "-",                "MRI-only: T1b detector (rest vs AD, binary)"),
    ("detector_structured", "-",           "detector-augmented (CN/AD<-MRI binary; MCI=CL)"),
    ("elastic_net",    "detectors",        "elastic-net [CL3 + T1_CN + T1b_AD]"),
]
# METRIC_COLS is built by _set_classes (3 headline cols + one F1 per class) so it
# adapts to binary tasks (e.g. T1d sMCI/pMCI). Do not redefine statically here.
TP_LABEL = {"BL": "CL_bl + MRI_bl", "M12": "CL_bl + MRI_m12",
            "M12X": "CL_m12 + MRI_m12 (cross-sectional, concurrent label)",
            "M12_DET": "CL_m12 (3-class) augmented by MRI detectors T1->P(CN), T1b->P(AD); MCI=CL"}

FOOTNOTE_LINES = [
    "Meta-learners fit per-seed on VAL, reported on TEST.",
    "Variants:",
    "  clinical-only / MRI-only = raw argmax of that modality's class probs (reference).",
    "  weighted-avg  = w*CL + (1-w)*MRI per-class prob; w tuned on VAL bACC (both-present).",
    "  elastic-net   = multinomial logistic (L1/L2 mix) over stacked [CL, MRI] probs; CV grid on VAL.",
    "Missingness (baseline MRI is sparse):",
    "  both present  = complete-case (CL and MRI both available at this timepoint).",
    "  zero-fallback = EN fit on both-present; missing-MRI patient gets MRI probs=0 (clinical proceeds).",
    "  impute+indic. = missing modality filled with TRAIN-fold class prior + binary 'present' flag; fit on all.",
    "(n) after bACC = mean #TEST patients scored across seeds. '-' = variant not computable (too few/1-class VAL).",
]


def _fmt(mean, std, n=None):
    if mean is None or (isinstance(mean, float) and np.isnan(mean)):
        return "-"
    s = f"{mean:.3f}"
    if std is not None and not (isinstance(std, float) and np.isnan(std)):
        s += f" ± {std:.3f}"
    if n is not None and not np.isnan(n):
        s += f" ({int(round(n))})"
    return s


def render_confusions(preds, framing, picks, out_path, task_label="T2"):
    """One row of per-class confusion matrices (counts pooled across seeds) for the
    given (variant, missing, label) picks -- e.g. clinical_only, mri_only, best fusion.
    Each cell shows count and row-normalised %; the cohort n is printed per panel since
    references and complete-case fusion can be on different cohorts."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARN] matplotlib not installed -- skipping confusion PNG.")
        return
    sub_all = preds[preds["framing"] == framing]
    if sub_all.empty:
        return
    panels = []
    for variant, missing, label, *rest in picks:
        tp = rest[0] if rest else ""
        s = sub_all[(sub_all["variant"] == variant) & (sub_all["missing"] == missing)]
        if s.empty:
            continue
        y = s["y_true"].to_numpy(int)
        p = s["pred"].to_numpy(int)
        cm = confusion_matrix(y, p, labels=list(range(N_CLASSES)))
        bacc = balanced_accuracy_score(y, p) if len(np.unique(y)) > 1 else float("nan")
        n_per_seed = len(s) // s["seed"].nunique() if s["seed"].nunique() else len(s)
        panels.append((label, tp, cm, n_per_seed, len(s), bacc))
    if not panels:
        return
    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(3.6 * n, 4.3))
    if n == 1:
        axes = [axes]
    for ax, (label, tp, cm, n_per_seed, n_tot, bacc) in zip(axes, panels):
        row = cm.sum(1, keepdims=True)
        pct = np.divide(cm, np.where(row == 0, 1, row)) * 100
        ax.imshow(pct, cmap="Blues", vmin=0, vmax=100)
        ax.set_xticks(range(N_CLASSES)); ax.set_yticks(range(N_CLASSES))
        ax.set_xticklabels(CLASS_NAMES); ax.set_yticklabels(CLASS_NAMES)
        ax.set_xlabel("predicted"); ax.set_ylabel("true")
        for i in range(N_CLASSES):
            for j in range(N_CLASSES):
                ax.text(j, i, f"{cm[i, j]}\n{pct[i, j]:.0f}%", ha="center", va="center",
                        fontsize=8.5, color="white" if pct[i, j] > 55 else "black")
        tp_line = f"{tp}\n" if tp else ""
        ax.set_title(f"{label}\n{tp_line}{n_per_seed}/seed TEST (Σ{n_tot})  bACC {bacc:.3f}",
                     fontsize=9)
    fig.suptitle("Multimodal late integration: Confusion matrix", fontsize=13, fontweight="bold",
                 y=0.99)
    fig.text(0.5, 0.925, "counts pooled across seeds, row-normalised %",
             ha="center", va="center", fontsize=9.5, fontstyle="italic")
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    fig.savefig(out_path, dpi=170, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"  confusion: {out_path}")


def render_framing_png(summary, cov, framing, mri_name, clin_label, out_path,
                       task_label="T2", tp_label=None):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARN] matplotlib not installed -- skipping PNG.")
        return
    sub = summary[summary["framing"] == framing]
    if sub.empty:
        return
    idx = {(r["variant"], r["missing"]): r for _, r in sub.iterrows()}
    # Make the single-modality reference rows name the exact model used.
    label_override = {
        ("clinical_only", "-"): f"clinical-only [{clin_label}]",
        ("mri_only", "-"):      f"MRI-only [{mri_name}]",
    }
    rows, body, numeric = [], [], []
    for variant, missing, label in ROW_ORDER:
        if (variant, missing) not in idx:
            continue
        r = idx[(variant, missing)]
        disp = label_override.get((variant, missing), label)
        rows.append(disp)
        cells, nums = [], []
        for key, _, with_n in METRIC_COLS:
            n = r["n_test_mean"] if with_n else None
            cells.append(_fmt(r.get(f"{key}_mean"), r.get(f"{key}_std"), n))
            nums.append(r.get(f"{key}_mean", np.nan))
        body.append([disp] + cells)
        numeric.append(nums)
    if not body:
        return
    numeric = np.array(numeric, dtype=float)

    n_rows, n_cols = len(body), len(METRIC_COLS) + 1
    col_labels = ["Variant (fusion method / model)"] + [c[1] for c in METRIC_COLS]

    # Size EVERY column to its widest cell so nothing clips. Width is measured
    # in characters (header + all body cells), converted to inches at a fixed
    # per-char width; colWidths are fractions of the axes, fig width is set so
    # each column's fraction maps back to its char budget.
    CHAR_W = 0.105                       # inches per char at fontsize 9 (+margin)
    col_chars = []
    for j in range(n_cols):
        texts = [col_labels[j]] + [str(b[j]) for b in body]
        col_chars.append(max(len(t) for t in texts) + 2)   # +2 cell padding
    total_chars = sum(col_chars)
    col_widths = [c / total_chars for c in col_chars]
    fig_w = total_chars * CHAR_W

    fig, ax = plt.subplots(figsize=(fig_w, 0.55 + 0.32 * (n_rows + 2)))
    ax.axis("off")
    tab = ax.table(cellText=body, colLabels=col_labels, loc="upper center",
                   cellLoc="center", colLoc="center", colWidths=col_widths)
    tab.auto_set_font_size(False)
    tab.set_fontsize(9)
    tab.scale(1.0, 1.4)

    # bold best (max mean) per metric column
    for j in range(len(METRIC_COLS)):
        col = numeric[:, j]
        if np.all(np.isnan(col)):
            continue
        best = int(np.nanargmax(col))
        tab[best + 1, j + 1].set_text_props(weight="bold")
    for j in range(n_cols):
        tab[0, j].set_facecolor("#ECECEC")
        tab[0, j].set_text_props(weight="bold")
    for i in range(n_rows + 1):
        tab[i, 0].set_text_props(ha="left")
        tab[i, 0].PAD = 0.03

    cov_f = cov[cov["framing"] == framing]
    cov_line = ("  coverage (test, both modalities present per seed): "
                + ", ".join(f"s{int(c.seed)}={int(c.n_test_both)}/{int(c.n_test)}"
                            for c in cov_f.itertuples()))
    plt.title(f"{task_label} late fusion -- {framing} ({tp_label or TP_LABEL[framing]})\n"
              f"clinical={clin_label}   MRI={mri_name}   "
              f"(mean ± std across seeds, n)", pad=6, fontsize=11)
    # Detector-augmented: spell out the class-specific renormalised late-fusion formula.
    detector_lines = []
    if framing == "M12_DET":
        detector_lines = [
            "",
            "detector-augmented (class-specific late fusion, renormalised):",
            "  p(CN)_fused  = w_cn*CL_CN + (1-w_cn)*MRI_T1(CN)        [T1 = CN-vs-rest detector]",
            "  p(AD)_fused  = w_ad*CL_AD + (1-w_ad)*MRI_T1b(AD)       [T1b = rest-vs-AD detector]",
            "  p(MCI)_fused = CL_MCI                                  [MCI from clinical only]",
            "  then renormalise so p(CN)+p(MCI)+p(AD)=1. w_cn,w_ad tuned per-seed on VAL;",
            "  a missing detector -> that class falls back to CL (no augmentation).",
            "MRI-only T1/T1b rows = each BINARY detector's performance on its OWN task (same m12 test",
            "  cohort): T1 = CN-vs-rest, T1b = rest-vs-AD. bACC/macro-AUC are BINARY; F1 is shown in the",
            "  detector's target-class column (F1 CN for T1, F1 AD for T1b). MRI model named in the title.",
        ]
    foot = "\n".join(FOOTNOTE_LINES + detector_lines + ["", cov_line])
    ax.text(0.0, -0.04, foot, transform=ax.transAxes, ha="left", va="top",
            fontsize=7.5, family="monospace", linespacing=1.4)
    fig.savefig(out_path, dpi=180, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print(f"  PNG       : {out_path}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mri-name", required=True,
                    help="Label for the MRI model shown in the table "
                         "(e.g. BrainMVP_full_ft_stochastic).")
    ap.add_argument("--out-name", default=None,
                    help="Output subfolder name (default: --mri-name). Use to "
                         "name the folder by timeframe, e.g. baseline_only, "
                         "while --mri-name still labels the model in the table.")
    ap.add_argument("--mri-csv-template", required=False,
                    help="MRI extraction CSV path with a {seed} placeholder "
                         "(required unless --detectors).")
    ap.add_argument("--clin-root", type=Path, default=CLIN_ROOT)
    ap.add_argument("--clin-model", default=CLIN_MODEL)
    ap.add_argument("--clin-strat", default=CLIN_STRAT)
    ap.add_argument("--clin-task", default="T2_multiclass",
                    help="Clinical parquet task subdir. Use T2_m12_multiclass "
                         "with --framings M12X for the cross-sectional m12 fusion.")
    ap.add_argument("--clin-probs-root", type=Path, default=None,
                    help="Optional: load the clinical PROBS from a DIFFERENT parquet "
                         "tree than the cohort/label/split anchor (--clin-root/--clin-task). "
                         "Use to fuse CL_bl probs onto the m12 concurrent-label cohort "
                         "(present-only CL_bl + MRI_m12): anchor on T2_m12_multiclass, "
                         "probs from T2_multiclass. Defaults to --clin-root.")
    ap.add_argument("--clin-probs-task", default=None,
                    help="Task subdir for --clin-probs-root (defaults to --clin-task).")
    ap.add_argument("--mri-timepoint", default=None, choices=["bl", "m12"],
                    help="Override the MRI session joined for ALL framings (default: bl "
                         "for BL, m12 for M12/M12X). E.g. --framings BL --mri-timepoint m12 "
                         "= baseline-anchored present-only (full baseline clinical cohort) "
                         "fused with MRI at m12.")
    ap.add_argument("--tp-label", default=None,
                    help="Override the timepoint description in the table title (default "
                         "per-framing). Use to state the true clinical-probs source, e.g. "
                         "'CL_bl + MRI_m12 (baseline label, present-only)'.")
    ap.add_argument("--table-tag", default=None,
                    help="Name the rendered PNG fusion_table_<table_tag>.png instead of "
                         "fusion_table_<framing>.png (e.g. CLbl_Mm12), so the file is "
                         "named by clinical arm rather than the internal framing code.")
    ap.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS_DEFAULT))
    ap.add_argument("--n-classes", type=int, default=3,
                    help="Number of classes (3 for T2 CN/MCI/AD; 2 for T1d sMCI/pMCI).")
    ap.add_argument("--class-names", nargs="+", default=None,
                    help="Class names in index order (default CN MCI AD for 3, "
                         "class0..class{K-1} otherwise). E.g. --class-names sMCI pMCI.")
    ap.add_argument("--framings", nargs="+", default=["BL", "M12"],
                    choices=["BL", "M12", "M12X"])
    ap.add_argument("--detectors", action="store_true",
                    help="Detector-augmented mode: CL_m12 3-class sharpened by MRI "
                         "T1 (CN) + T1b (AD), MCI from clinical. Ignores --framings.")
    ap.add_argument("--mri-t1-template",
                    help="MRI T1 (CN-vs-rest) extraction CSV template with {seed}; "
                         "prob_class_0 is used as P(CN).")
    ap.add_argument("--mri-t1b-template",
                    help="MRI T1b (rest-vs-AD) extraction CSV template with {seed}; "
                         "prob_class_1 is used as P(AD).")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR_DEFAULT)
    ap.add_argument("--task-label", default="T2",
                    help="Task name shown in the table title (e.g. T1d).")
    args = ap.parse_args()

    # Configure class count / names (rebuilds CP_COLS/MP_COLS/PROB_COLS/METRIC_COLS).
    if args.class_names:
        names = args.class_names
    elif args.n_classes == 3:
        names = ["CN", "MCI", "AD"]
    else:
        names = [f"class{i}" for i in range(args.n_classes)]
    if len(names) != args.n_classes:
        ap.error(f"--class-names ({len(names)}) must match --n-classes ({args.n_classes})")
    _set_classes(names)

    if args.detectors:
        if args.n_classes != 3:
            ap.error("--detectors is the 3-class T2 CN/AD path; not valid for n_classes != 3")
        if not (args.mri_t1_template or args.mri_t1b_template):
            ap.error("--detectors requires at least one of --mri-t1-template / --mri-t1b-template")
    elif not args.mri_csv_template:
        ap.error("--mri-csv-template is required unless --detectors")

    out_dir = args.out_dir / (args.out_name or args.mri_name)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_metrics, all_preds, all_fail, all_stats = [], [], [], []
    coverage = []

    framings_run = ["M12_DET"] if args.detectors else args.framings
    for seed in args.seeds:
        clin = load_clinical(seed, args.clin_root, args.clin_model,
                             args.clin_strat, task=args.clin_task)
        # Optionally swap in clinical PROBS from a different parquet (e.g. CL_bl)
        # while keeping the anchor's cohort/label/split (e.g. m12 concurrent dx).
        # Splits are verified consistent across the bl/m12 trees, so this is
        # leakage-safe (a m12-test patient is never a bl-train patient).
        if args.clin_probs_root is not None or args.clin_probs_task is not None:
            cp_root = args.clin_probs_root or args.clin_root
            cp_task = args.clin_probs_task or args.clin_task
            clin_p = load_clinical(seed, cp_root, args.clin_model,
                                   args.clin_strat, task=cp_task)
            clin = (clin.drop(columns=CP_COLS)
                    .merge(clin_p[["Patient_ID"] + CP_COLS], on="Patient_ID", how="left"))
        prior = class_prior_from_train(clin)
        if args.detectors:
            t1 = (load_mri_detector(seed, args.mri_t1_template, "m12", "prob_class_0", "mri_cn")
                  if args.mri_t1_template else None)
            t1b = (load_mri_detector(seed, args.mri_t1b_template, "m12", "prob_class_1", "mri_ad")
                   if args.mri_t1b_template else None)
            rows, fr, fail_df = evaluate_detectors(seed, clin, t1, t1b, all_preds)
            all_metrics.extend(rows)
            all_fail.append(fail_df)
            te = fr[fr["split"] == "test"]
            coverage.append({
                "framing": "M12_DET", "seed": seed, "n_test": int(len(te)),
                "n_test_both": int(te["mri_present"].sum()),
                "n_test_clin_only": int((~te["mri_present"]).sum()),
                "n_test_mri_only": 0,
            })
            continue
        for framing in args.framings:
            tp = args.mri_timepoint or ("bl" if framing == "BL" else "m12")
            mri = load_mri(seed, args.mri_csv_template, tp)
            rows, fr = evaluate_framing(framing, seed, clin, mri, prior, all_preds)
            all_metrics.extend(rows)
            all_fail.append(failure_table(framing, seed, fr))
            all_stats.append(failure_stats(framing, seed, fr))
            te = fr[fr["split"] == "test"]
            coverage.append({
                "framing": framing, "seed": seed,
                "n_test": int(len(te)),
                "n_test_both": int((te["clin_present"] & te["mri_present"]).sum()),
                "n_test_clin_only": int((te["clin_present"] & ~te["mri_present"]).sum()),
                "n_test_mri_only": int((~te["clin_present"] & te["mri_present"]).sum()),
            })

    met = pd.DataFrame(all_metrics)
    cov = pd.DataFrame(coverage)
    preds = pd.DataFrame(all_preds)
    fails = pd.concat(all_fail, ignore_index=True) if all_fail else pd.DataFrame()

    # ---- aggregate across seeds: mean +/- std per (framing, variant, missing) ----
    agg_cols = ["val_bacc", "bacc", "macro_f1", "macro_auc"] + [f"f1_{n}" for n in CLASS_NAMES]
    summary = (met.groupby(["framing", "variant", "missing"])[agg_cols]
               .agg(["mean", "std"]).round(4))
    summary.columns = [f"{a}_{b}" for a, b in summary.columns]
    summary = summary.reset_index()
    n_per = met.groupby(["framing", "variant", "missing"])["n"].mean().round(1)
    summary = summary.merge(n_per.rename("n_test_mean").reset_index(),
                            on=["framing", "variant", "missing"])

    # ---- write ----
    met.to_csv(out_dir / "fusion_metrics_per_seed.csv", index=False)
    summary.to_csv(out_dir / "fusion_metrics.csv", index=False)
    cov.to_csv(out_dir / "coverage.csv", index=False)
    preds.to_csv(out_dir / "fused_predictions.csv", index=False)
    if len(fails):
        fails.to_csv(out_dir / "failure_table.csv", index=False)
    stats = pd.concat(all_stats, ignore_index=True) if all_stats else pd.DataFrame()
    if len(stats):
        stats.to_csv(out_dir / "failure_stats.csv", index=False)

    # ---- per-class confusion matrices: clinical_only, mri_only, best fusion ----
    if len(preds):
        # best fusion = highest mean test bACC among non-reference variants.
        fusion = summary[~summary["variant"].isin(["clinical_only", "mri_only"])]
        best = None
        if len(fusion):
            br = fusion.sort_values("bacc_mean", ascending=False).iloc[0]
            best = (br["variant"], br["missing"])
        _tp_word = lambda c: {"bl": "baseline", "m12": "month 12"}.get(c, c)
        clin_tp = _tp_word("m12" if "m12" in (args.clin_probs_task or args.clin_task) else "bl")
        for framing in framings_run:
            mri_tp = _tp_word(args.mri_timepoint or ("bl" if framing == "BL" else "m12"))
            picks = [("clinical_only", "-", f"clinical-only [{args.clin_model}]", clin_tp),
                     ("mri_only", "-", f"MRI-only [{args.mri_name}]", mri_tp)]
            if best is not None:
                picks.append((best[0], best[1], f"best fusion [{best[0]} / {best[1]}]",
                              f"CL {clin_tp} + MRI {mri_tp}"))
            render_confusions(preds, framing, picks,
                              out_dir / f"confusion_{args.table_tag or framing}.png",
                              task_label=args.task_label)

    # ---- rendered PNG table per framing (style of 06_render_cross_model_table) ----
    clin_label = f"{args.clin_model}/{args.clin_strat}"
    for framing in framings_run:
        tag = args.table_tag or framing
        render_framing_png(summary, cov, framing, args.mri_name, clin_label,
                           out_dir / f"fusion_table_{tag}.png",
                           task_label=args.task_label, tp_label=args.tp_label)

    print("=" * 72)
    print(f"  T2 late fusion -- clinical={args.clin_model}/{args.clin_strat}  "
          f"mri={args.mri_name}")
    print("=" * 72)
    print("\n[coverage] test-set modality overlap per framing/seed:")
    print(cov.to_string(index=False))
    print("\n[summary] mean over seeds (fit on VAL, reported on TEST):")
    show = ["framing", "variant", "missing", "n_test_mean",
            "bacc_mean", "bacc_std", "macro_f1_mean", "macro_auc_mean"]
    print(summary[show].to_string(index=False))
    print(f"\n  wrote -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
