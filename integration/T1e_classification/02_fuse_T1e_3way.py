"""
02_fuse_T1e_3way.py
===================
THREE-way late fusion for T1e_scn_pcn (binary: sCN=0, pCN=1) over per-patient
class probabilities from a CLINICAL encoder, a SNP-side model, and a frozen
MRI ViT encoder.

This extends the 2-way `01_fuse_T1e.py` (clinical + SNP) by adding an MRI arm.
It does NOT touch `01`'s published `outputs/baseline/` results -- everything
here is written to `outputs/3way/`.

Modalities & sources (all local D:):
  CL  clinical encoder parquet (no_cdr_post_exclusion / baseline):
        encoder_outputs_no_cdr_post_exclusion/<model>/T1e_scn_pcn/seed_{s}/<strat>/
            embeddings/embeddings.parquet
        cols: Patient_ID, split, in_task_cohort, label, prob_0, prob_1
  SNP per-seed per-patient parquet (val+test only):
        cols: Patient_ID, fold, outcome_proba, y_true
  MRI frozen cached-head per-scan probs (this repo's 05_extract_cached_head_probs.py,
        with the T1e branch):
        mri_embeddings/T1e_pcn_vs_scn/<arch>_frozen_none_cached/embeddings_seed_<s>.csv
        cols: Patient_ID, adni_viscode, split, label, prob_class_0, prob_class_1, ...
        -> the EARLIEST scan per subject is kept (T1e session_policy="conversion",
           mirroring 04_supervised_finetuning_ViT.py:_resolve_labels).

Each base model already saw TRAIN. Every fusion weight / meta-learner is tuned
on VAL only (no train-fold leakage); TEST is held out for reporting. SNP/MRI
predictions are per-seed; CL parquet carries its own split column.

Methods (per CL x SNP x MRI combo, reported on val + test):
  references : clinical_only, snp_only, mri_only        (argmax of one modality)
  pairwise   : cl+snp, cl+mri, snp+mri                  (alpha tuned on VAL bACC)
  three-way  : equal_avg3   (1/3 each)
               best_weights3 (w_cl,w_snp,w_mri on the simplex, step 0.1, VAL bACC)
               stacked_en3   (LogisticRegression elastic-net over the 6 prob cols)

Two cohorts per method:
  complete : subjects present in ALL three modalities (clinical AND SNP AND MRI).
  present  : clinical-anchored present-only -- keep every clinical subject; for a
             subject missing SNP and/or MRI, drop the absent term(s) and
             renormalise the weights (per-subject convex combination). stacked_en3
             needs a fixed-width feature vector, so it is complete-case only.

Outputs (outputs/3way/):
  fusion_metrics_per_seed.csv  fusion_metrics.csv  coverage.csv
  fused_predictions.csv        failure_contingency.csv
  fusion_table_T1e_3way.png    fusion_leaderboard_T1e_3way.png

Usage:
  python integration/T1e_classification/02_fuse_T1e_3way.py
  python integration/T1e_classification/02_fuse_T1e_3way.py \
      --clin-keys modernbert_base_full_ft --snp-keys meta_prs --mri-keys vit_mae
"""
from __future__ import annotations

import argparse
import re
import warnings
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (balanced_accuracy_score, f1_score,
                             precision_recall_fscore_support, recall_score,
                             roc_auc_score)
from sklearn.model_selection import GridSearchCV, StratifiedKFold

warnings.filterwarnings("ignore")

CLASS_NAMES = ["sCN", "pCN"]
MODS = ["clin", "snp", "mri"]          # canonical modality order for weight vectors
SEEDS_DEFAULT = (0, 1, 2)

# --------------------------------------------------------------------------- #
# Paths & variant registries
# --------------------------------------------------------------------------- #
CLIN_ROOT = Path(r"D:/ADNI_BIDS_project/derivatives/encoder_outputs_no_cdr_post_exclusion")
SNP_ROOT = Path(r"D:/ADNI_SNP_Omni2.5M_20140220")
MRI_ROOT = Path(r"D:/ADNI_BIDS_project/derivatives/mri_embeddings/T1e_pcn_vs_scn")
ANCHORS = SNP_ROOT / "outputs" / "multimodal_anchors"
NEW_TASKS_PERPT = (SNP_ROOT / "outputs" / "new_tasks" / "sCN_vs_pCN"
                   / "ld_1000kb_r2_0.8" / "per_patient")
OUT_DIR_DEFAULT = Path(__file__).resolve().parent / "outputs" / "3way"

CLIN_VARIANTS = [
    {"key": "bioclin_mb_base_frozen", "display": "BioClinical-MB-base [frozen]",
     "short": "BioClin-B/frz", "model": "BioClinical-ModernBERT-base", "strat": "frozen"},
    {"key": "modernbert_base_full_ft", "display": "ModernBERT-base [full_ft]",
     "short": "ModBERT-B/ft", "model": "ModernBERT-base", "strat": "full_ft"},
]

SNP_VARIANTS = [
    {"key": "bmfm_snp_attn1_mlp", "display": "BMFM-SNP attn_bias_single chrom_hier MLP",
     "short": "BMFM/attn1-MLP",
     "parquet_tmpl": str(ANCHORS / "bmfm_snp__attn_bias_single__chrom_hier__mlp2__101bp"
                         / "task_sCN_vs_pCN" / "seed_{seed}" / "prs+age+sex+apoe4"
                         / "predictions.parquet")},
    {"key": "bmfm_snp_none_xgb", "display": "BMFM-SNP none global_attn XGB",
     "short": "BMFM/none-XGB",
     "parquet_tmpl": str(ANCHORS / "bmfm_snp__none__global_attn__xgb__1001bp"
                         / "task_sCN_vs_pCN" / "seed_{seed}" / "prs+age+sex+apoe4"
                         / "predictions.parquet")},
    {"key": "meta_prs", "display": "meta_prs_EN_combined", "short": "meta_PRS",
     "parquet_tmpl": str(NEW_TASKS_PERPT
                         / "meta_prs_EN_combined__prs+age+sex+apoe4__seed{seed}.parquet")},
    {"key": "kosteridis_prs", "display": "Kosteridis_full_PRS", "short": "Kosteridis_PRS",
     "parquet_tmpl": str(NEW_TASKS_PERPT
                         / "Kosteridis__prs+age+sex+apoe4__seed{seed}.parquet")},
]

MRI_VARIANTS = [
    {"key": "vit_mae", "display": "ViT-MAE [frozen]", "short": "ViT-MAE/frz",
     "arch_dir": "vit_mae_frozen_none_cached"},
    {"key": "braindino", "display": "BrainDINO ViT-B/16 [frozen]", "short": "BrainDINO/frz",
     "arch_dir": "braindino_frozen_none_cached"},
]


# --------------------------------------------------------------------------- #
# Loaders
# --------------------------------------------------------------------------- #
def load_clinical(seed, model_slug, strategy, task="T1e_scn_pcn", clin_root=CLIN_ROOT):
    """Clinical probs -> Patient_ID, split, y_clin, cp0, cp1."""
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


def load_snp(seed, parquet_tmpl):
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


def _months(viscode):
    if viscode == "bl":
        return 0
    m = re.match(r"^m(\d+)$", str(viscode))
    return int(m.group(1)) if m else None


def load_mri_t1e(seed, arch_dir, mri_root=MRI_ROOT):
    """MRI frozen-cached probs -> Patient_ID, split, y_mri, mp0, mp1, at the
    EARLIEST scan per subject (T1e session_policy='conversion'; mirrors
    04_supervised_finetuning_ViT.py:_resolve_labels). The extractor writes one
    row per (subject, scan) for the whole master with a subject-level T1e label
    (label==-1 for out-of-cohort scans), so we filter to in-cohort then collapse
    to the earliest session per subject."""
    p = mri_root / arch_dir / f"embeddings_seed_{seed}.csv"
    if not p.exists():
        raise FileNotFoundError(f"MRI csv missing: {p}")
    d = pd.read_csv(p, low_memory=False)
    d = d[(d["label"] != -1) & d["split"].isin(["train", "val", "test"])].copy()
    d["_m"] = d["adni_viscode"].map(_months)
    d = d.dropna(subset=["_m"])
    d = (d.sort_values(["Patient_ID", "_m"])
           .groupby("Patient_ID", as_index=False).first())
    d = d.rename(columns={"prob_class_0": "mp0", "prob_class_1": "mp1",
                          "label": "y_mri"})
    d["y_mri"] = d["y_mri"].astype(int)
    return d[["Patient_ID", "split", "y_mri", "mp0", "mp1"]]


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #
def safe_auc(y_true, prob_pos):
    try:
        if len(np.unique(y_true)) < 2:
            return float("nan")
        return float(roc_auc_score(y_true, prob_pos))
    except Exception:
        return float("nan")


def metric_row(y_true, y_pred, probs):
    bacc = balanced_accuracy_score(y_true, y_pred)
    macro_f1 = f1_score(y_true, y_pred, average="macro", labels=[0, 1], zero_division=0)
    _, rc, f1, _ = precision_recall_fscore_support(
        y_true, y_pred, labels=[0, 1], zero_division=0)
    out = {"n": int(len(y_true)), "bacc": float(bacc), "macro_f1": float(macro_f1),
           "auc": safe_auc(y_true, probs[:, 1]) if probs is not None else float("nan")}
    for i, name in enumerate(CLASS_NAMES):
        out[f"f1_{name}"] = float(f1[i])
        out[f"recall_{name}"] = float(rc[i])
    return out


# --------------------------------------------------------------------------- #
# Frame builder (CL outer-join SNP outer-join MRI)
# --------------------------------------------------------------------------- #
def build_frame(clin, snp, mri):
    """Outer-join the three modalities on (Patient_ID, split); presence flags;
    canonical label from whichever modality is present (CL > SNP > MRI)."""
    fr = clin.merge(snp, on=["Patient_ID", "split"], how="outer")
    fr = fr.merge(mri, on=["Patient_ID", "split"], how="outer")
    fr["clin_present"] = fr["cp0"].notna()
    fr["snp_present"] = fr["sp0"].notna()
    fr["mri_present"] = fr["mp0"].notna()
    # canonical label
    fr["y_true"] = fr["y_clin"]
    for alt in ("y_snp", "y_mri"):
        fr["y_true"] = fr["y_true"].where(fr["y_true"].notna(), fr[alt])
    fr = fr[fr["y_true"].notna()].copy()
    fr["y_true"] = fr["y_true"].astype(int)
    # neutral 0.5/0.5 fill for absent modalities -- only ever used where the
    # corresponding presence flag is False, and `blend` zeroes those terms out.
    for c in ("cp0", "cp1", "sp0", "sp1", "mp0", "mp1"):
        fr[c] = fr[c].fillna(0.5)
    return fr


def _probs(df):
    """Return list of 3 (n,2) arrays in MODS order."""
    return [df[["cp0", "cp1"]].to_numpy(float),
            df[["sp0", "sp1"]].to_numpy(float),
            df[["mp0", "mp1"]].to_numpy(float)]


def _present(df):
    return [df["clin_present"].to_numpy(bool),
            df["snp_present"].to_numpy(bool),
            df["mri_present"].to_numpy(bool)]


# --------------------------------------------------------------------------- #
# Blend + tuning
# --------------------------------------------------------------------------- #
def blend(w, P, present):
    """Per-subject renormalised weighted average over PRESENT modalities.
        p_fused[i] = sum_m w_m * present_m[i] * P_m[i]  /  sum_m w_m * present_m[i]
    For the complete-case cohort all presence flags are True, so this reduces to
    the ordinary weighted average. Rows with zero usable weight (no modality with
    w>0 present) fall back to a uniform [0.5, 0.5]."""
    n = P[0].shape[0]
    num = np.zeros((n, 2))
    den = np.zeros((n, 1))
    for wi, Pi, pri in zip(w, P, present):
        mask = pri.astype(float)[:, None]
        num += wi * Pi * mask
        den += wi * mask
    out = np.divide(num, den, out=np.full((n, 2), 0.5), where=(den > 0))
    return out


def _simplex(step=0.1):
    """All (w_cl,w_snp,w_mri) >= 0 summing to 1 on a `step` grid."""
    k = int(round(1.0 / step))
    pts = []
    for a in range(k + 1):
        for b in range(k + 1 - a):
            c = k - a - b
            pts.append((a / k, b / k, c / k))
    return pts


def tune_weights(cands, P_v, present_v, y_v):
    """Pick the weight tuple in `cands` maximising VAL balanced accuracy."""
    best_w, best_b = cands[0], -1.0
    for w in cands:
        pred = blend(np.array(w), P_v, present_v).argmax(1)
        b = balanced_accuracy_score(y_v, pred)
        if b > best_b:
            best_b, best_w = b, w
    return np.array(best_w)


def fit_stacked_en(P_v, y_v, P_t, seed):
    """Elastic-net LR over the stacked 6-col prob feature [cp,sp,mp]. Complete-case
    only (fixed-width features). Returns (val_proba, test_proba, info)."""
    Xv = np.column_stack(P_v)
    Xt = np.column_stack(P_t)
    if len(np.unique(y_v)) < 2:
        return None
    n_min = int(np.min(np.bincount(y_v, minlength=2)))
    if n_min < 2:
        clf = LogisticRegression(penalty="elasticnet", solver="saga", C=1.0,
                                 l1_ratio=0.5, max_iter=5000, random_state=seed)
        clf.fit(Xv, y_v)
        info = {"C": 1.0, "l1_ratio": 0.5}
        est = clf
    else:
        grid = {"C": [0.1, 1.0, 10.0], "l1_ratio": [0.1, 0.5, 0.9]}
        base = LogisticRegression(penalty="elasticnet", solver="saga",
                                  max_iter=5000, random_state=seed)
        cv = StratifiedKFold(n_splits=min(3, n_min), shuffle=True, random_state=seed)
        clf = GridSearchCV(base, grid, scoring="balanced_accuracy", cv=cv, n_jobs=-1)
        clf.fit(Xv, y_v)
        info = {"C": float(clf.best_params_["C"]),
                "l1_ratio": float(clf.best_params_["l1_ratio"])}
        est = clf.best_estimator_

    def _full(X):
        proba = est.predict_proba(X)
        full = np.zeros((proba.shape[0], 2))
        for j, c in enumerate(est.classes_):
            full[:, int(c)] = proba[:, j]
        return full
    return _full(Xv), _full(Xt), info


# --------------------------------------------------------------------------- #
# Per-combo evaluation
# --------------------------------------------------------------------------- #
def _method_specs():
    """method_name -> (kind, candidate-weights). kind in {fixed, tune, en}."""
    e3 = (1 / 3, 1 / 3, 1 / 3)
    a = np.linspace(0.0, 1.0, 21)
    return {
        "clinical_only": ("fixed", [(1, 0, 0)]),
        "snp_only":      ("fixed", [(0, 1, 0)]),
        "mri_only":      ("fixed", [(0, 0, 1)]),
        "cl+snp":        ("tune", [(w, 1 - w, 0) for w in a]),
        "cl+mri":        ("tune", [(w, 0, 1 - w) for w in a]),
        "snp+mri":       ("tune", [(0, w, 1 - w) for w in a]),
        "equal_avg3":    ("fixed", [e3]),
        "best_weights3": ("tune", _simplex(0.1)),
        "stacked_en3":   ("en", None),
    }


# methods that are genuine fusion (>=2 modalities) -- stashed for prediction dumps
FUSION_METHODS = ["cl+snp", "cl+mri", "snp+mri",
                  "equal_avg3", "best_weights3", "stacked_en3"]


def evaluate_combo(seed, clin_v, snp_v, mri_v, clin, snp, mri,
                   pred_rows, contingency_rows):
    fr = build_frame(clin, snp, mri)
    rows = []
    specs = _method_specs()

    for cohort in ("complete", "present"):
        if cohort == "complete":
            sub = fr[fr["clin_present"] & fr["snp_present"] & fr["mri_present"]]
        else:                          # clinical-anchored present-only
            sub = fr[fr["clin_present"]]
        val = sub[sub["split"] == "val"]
        test = sub[sub["split"] == "test"]
        if len(val) < 4 or len(test) == 0:
            continue
        Pv, Pt = _probs(val), _probs(test)
        prv, prt = _present(val), _present(test)
        yv = val["y_true"].to_numpy(int)
        yt = test["y_true"].to_numpy(int)

        for method, (kind, cands) in specs.items():
            if kind == "en":
                if cohort != "complete":      # EN needs fixed-width features
                    continue
                res = fit_stacked_en(Pv, yv, Pt, seed)
                if res is None:
                    continue
                v_proba, t_proba, info = res
                v_pred, t_pred = v_proba.argmax(1), t_proba.argmax(1)
                wtuple = (float("nan"),) * 3
            else:
                if kind == "tune":
                    w = tune_weights(cands, Pv, prv, yv)
                else:                          # fixed
                    w = np.array(cands[0], dtype=float)
                v_proba, t_proba = blend(w, Pv, prv), blend(w, Pt, prt)
                v_pred, t_pred = v_proba.argmax(1), t_proba.argmax(1)
                wtuple = tuple(float(x) for x in w)
                info = {}

            for split_name, y, pred, proba, dframe in [
                ("val", yv, v_pred, v_proba, val),
                ("test", yt, t_pred, t_proba, test),
            ]:
                r = metric_row(y, pred, proba)
                r.update({
                    "clin_variant": clin_v["key"], "snp_variant": snp_v["key"],
                    "mri_variant": mri_v["key"], "method": method,
                    "cohort": cohort, "split": split_name,
                    "w_clin": wtuple[0], "w_snp": wtuple[1], "w_mri": wtuple[2],
                    "best_C": info.get("C", float("nan")),
                    "best_l1_ratio": info.get("l1_ratio", float("nan")),
                })
                rows.append(r)
                if method in FUSION_METHODS and split_name == "test":
                    for i, (_, row) in enumerate(dframe.reset_index(drop=True).iterrows()):
                        pred_rows.append({
                            "clin_variant": clin_v["key"], "snp_variant": snp_v["key"],
                            "mri_variant": mri_v["key"], "method": method,
                            "cohort": cohort, "seed": seed,
                            "Patient_ID": row["Patient_ID"], "y_true": int(row["y_true"]),
                            "pred": int(pred[i]),
                            "prob_sCN": float(proba[i, 0]), "prob_pCN": float(proba[i, 1]),
                            "clin_present": bool(row["clin_present"]),
                            "snp_present": bool(row["snp_present"]),
                            "mri_present": bool(row["mri_present"]),
                        })

    # failure contingency on the complete-case TEST cohort: per subject, the
    # TP/FP/FN/TN category of each single modality (do failure modes overlap?).
    comp_test = fr[fr["clin_present"] & fr["snp_present"] & fr["mri_present"]
                   & (fr["split"] == "test")]
    if len(comp_test):
        P = _probs(comp_test)
        for i, (_, row) in enumerate(comp_test.reset_index(drop=True).iterrows()):
            y = int(row["y_true"])
            rec = {"clin_variant": clin_v["key"], "snp_variant": snp_v["key"],
                   "mri_variant": mri_v["key"], "seed": seed,
                   "Patient_ID": row["Patient_ID"], "y_true": CLASS_NAMES[y]}
            for m, Pi in zip(MODS, P):
                pr = int(Pi[i].argmax())
                rec[f"{m}_pred"] = CLASS_NAMES[pr]
                rec[f"{m}_cat"] = ("TP" if (pr == 1 and y == 1) else
                                   "TN" if (pr == 0 and y == 0) else
                                   "FP" if (pr == 1 and y == 0) else "FN")
            contingency_rows.append(rec)
    return rows


# --------------------------------------------------------------------------- #
# Coverage
# --------------------------------------------------------------------------- #
def compute_coverage(seed, clin_v, snp_v, mri_v, clin, snp, mri):
    fr = build_frame(clin, snp, mri)
    te = fr[fr["split"] == "test"]
    return {
        "clin_variant": clin_v["key"], "snp_variant": snp_v["key"],
        "mri_variant": mri_v["key"], "seed": seed,
        "n_test_clin": int(te["clin_present"].sum()),
        "n_test_snp": int(te["snp_present"].sum()),
        "n_test_mri": int(te["mri_present"].sum()),
        "n_test_all3": int((te["clin_present"] & te["snp_present"]
                            & te["mri_present"]).sum()),
    }


# --------------------------------------------------------------------------- #
# Rendering
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


METHOD_SHORT = {"clinical_only": "CL", "snp_only": "SNP", "mri_only": "MRI",
                "cl+snp": "CL+SNP", "cl+mri": "CL+MRI", "snp+mri": "SNP+MRI",
                "equal_avg3": "avg3", "best_weights3": "best3", "stacked_en3": "stack-EN3"}


def render_leaderboard(summary, met_per_seed, out_path, cohort="complete", top_n=20):
    """Leaderboard PNG: one row per (CL,SNP,MRI,method) at `cohort`, sorted by
    TEST bACC. Single-modality refs included once each."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARN] matplotlib missing -- skipping leaderboard PNG.")
        return
    s = summary[(summary["split"] == "test") & (summary["cohort"] == cohort)].copy()
    if s.empty:
        print(f"  [WARN] no rows for cohort={cohort} -- skipping leaderboard.")
        return

    # weights aggregate for the tuned-param column
    g = ["clin_variant", "snp_variant", "mri_variant", "method", "cohort", "split"]
    wm = (met_per_seed.groupby(g)[["w_clin", "w_snp", "w_mri"]]
          .mean().round(2).reset_index())
    s = s.merge(wm, on=g, how="left")

    # One reference row per distinct single-modality model (so BOTH MRI encoders
    # and every SNP/CL variant show, not just the single best per method).
    refs = pd.concat([
        s[s["method"] == "clinical_only"].drop_duplicates("clin_variant"),
        s[s["method"] == "snp_only"].drop_duplicates("snp_variant"),
        s[s["method"] == "mri_only"].drop_duplicates("mri_variant"),
    ], ignore_index=True)
    fus = s[~s["method"].isin(["clinical_only", "snp_only", "mri_only"])]
    fus = fus.sort_values("bacc_mean", ascending=False)
    pool = pd.concat([fus.head(top_n), refs], ignore_index=True)
    pool = pool.sort_values("bacc_mean", ascending=False).reset_index(drop=True)

    cl_disp = {v["key"]: v["short"] for v in CLIN_VARIANTS}
    sn_disp = {v["key"]: v["short"] for v in SNP_VARIANTS}
    mr_disp = {v["key"]: v["short"] for v in MRI_VARIANTS}
    is_ref = pool["method"].isin(["clinical_only", "snp_only", "mri_only"])

    col_labels = ["Rank", "CL", "SNP", "MRI", "Method", "bACC", "macro-F1",
                  "AUC", "weights (cl/snp/mri)"]
    body, numeric = [], []
    for rank, (_, r) in enumerate(pool.iterrows(), start=1):
        ref = r["method"] in ("clinical_only", "snp_only", "mri_only")
        w = "-" if ref or any(pd.isna(r.get(k)) for k in ("w_clin", "w_snp", "w_mri")) \
            else f"{r['w_clin']:.2f}/{r['w_snp']:.2f}/{r['w_mri']:.2f}"
        if r["method"] == "stacked_en3":
            w = "EN"
        body.append([
            str(rank),
            "(no CL)" if r["method"] == "snp_only" or r["method"] == "mri_only" else cl_disp.get(r["clin_variant"], r["clin_variant"]),
            "(no SNP)" if r["method"] in ("clinical_only", "mri_only") else sn_disp.get(r["snp_variant"], r["snp_variant"]),
            "(no MRI)" if r["method"] in ("clinical_only", "snp_only") else mr_disp.get(r["mri_variant"], r["mri_variant"]),
            METHOD_SHORT.get(r["method"], r["method"]),
            _fmt(r.get("bacc_mean"), r.get("bacc_std"), r.get("n_mean")),
            _fmt(r.get("macro_f1_mean"), r.get("macro_f1_std")),
            _fmt(r.get("auc_mean"), r.get("auc_std")),
            w,
        ])
        numeric.append([r.get("bacc_mean", np.nan), r.get("macro_f1_mean", np.nan),
                        r.get("auc_mean", np.nan)])
    numeric = np.array(numeric, float)

    CHAR_W = 0.105
    col_chars = [max(len(col_labels[j]), *(len(str(b[j])) for b in body)) + 2
                 for j in range(len(col_labels))]
    widths = [c / sum(col_chars) for c in col_chars]
    fig_w = sum(col_chars) * CHAR_W
    h_tab = 0.45 * (len(body) + 1)
    fig = plt.figure(figsize=(fig_w, h_tab + 1.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[h_tab, 1.2], hspace=0.05)
    ax, ax_ft = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])
    ax.axis("off"); ax_ft.axis("off")
    tab = ax.table(cellText=body, colLabels=col_labels, loc="upper center",
                   cellLoc="center", colLoc="center", colWidths=widths)
    tab.auto_set_font_size(False); tab.set_fontsize(9); tab.scale(1.0, 1.4)
    for j in range(len(col_labels)):
        tab[0, j].set_text_props(weight="bold")
    for i in range(len(body) + 1):
        for j in (0, 1, 2, 3, 4, 8):
            tab[i, j].set_text_props(ha="left"); tab[i, j].PAD = 0.03
    for j in range(3):
        col = numeric[:, j]
        if not np.all(np.isnan(col)):
            tab[int(np.nanargmax(col)) + 1, j + 5].set_text_props(weight="bold")
    ax.set_title(f"T1e (sCN vs pCN) 3-way late fusion -- LEADERBOARD ({cohort}-case TEST)\n"
                 "clinical + SNP + MRI(ViT-frozen); sorted by TEST bACC, mean over seeds 0,1,2",
                 pad=8, fontsize=11)
    foot = [
        "Weights (cl/snp/mri) tuned on VAL balanced accuracy; refs/avg3 fixed; stack-EN3 = elastic-net LR over the 6 prob cols.",
        f"cohort={cohort}: " + ("subjects present in ALL three modalities."
                                if cohort == "complete"
                                else "clinical-anchored; absent SNP/MRI terms dropped & weights renormalised per subject."),
        "(n) after bACC = mean #TEST subjects across seeds. Bold = best per metric column.",
    ]
    ax_ft.text(0.0, 1.0, "\n".join(foot), transform=ax_ft.transAxes, ha="left",
               va="top", fontsize=7.8, family="monospace", linespacing=1.4)
    fig.savefig(out_path, dpi=180, bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)
    print(f"  PNG  -> {out_path}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clin-keys", nargs="+", default=None)
    ap.add_argument("--snp-keys", nargs="+", default=None)
    ap.add_argument("--mri-keys", nargs="+", default=None)
    ap.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS_DEFAULT))
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR_DEFAULT)
    args = ap.parse_args()

    clin_run = [v for v in CLIN_VARIANTS if not args.clin_keys or v["key"] in args.clin_keys]
    snp_run = [v for v in SNP_VARIANTS if not args.snp_keys or v["key"] in args.snp_keys]
    mri_run = [v for v in MRI_VARIANTS if not args.mri_keys or v["key"] in args.mri_keys]
    if not (clin_run and snp_run and mri_run):
        print("  [ERROR] empty CL/SNP/MRI set after filtering."); return 1
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    all_metrics, all_preds, contingency, coverage = [], [], [], []
    for seed in args.seeds:
        clin_cache, mri_cache = {}, {}
        for v in clin_run:
            try:
                clin_cache[v["key"]] = load_clinical(seed, v["model"], v["strat"])
            except FileNotFoundError as e:
                print(f"  [WARN] {e}; skip CL '{v['key']}' seed={seed}")
        for v in mri_run:
            try:
                mri_cache[v["key"]] = load_mri_t1e(seed, v["arch_dir"])
            except FileNotFoundError as e:
                print(f"  [WARN] {e}; skip MRI '{v['key']}' seed={seed}")
        for clin_v in clin_run:
            if clin_v["key"] not in clin_cache:
                continue
            for snp_v in snp_run:
                try:
                    snp = load_snp(seed, snp_v["parquet_tmpl"])
                except FileNotFoundError as e:
                    print(f"  [WARN] {e}; skip SNP '{snp_v['key']}' seed={seed}")
                    continue
                for mri_v in mri_run:
                    if mri_v["key"] not in mri_cache:
                        continue
                    rows = evaluate_combo(
                        seed, clin_v, snp_v, mri_v,
                        clin_cache[clin_v["key"]], snp, mri_cache[mri_v["key"]],
                        all_preds, contingency)
                    all_metrics.extend(rows)
                    coverage.append(compute_coverage(
                        seed, clin_v, snp_v, mri_v,
                        clin_cache[clin_v["key"]], snp, mri_cache[mri_v["key"]]))

    if not all_metrics:
        print("  [ERROR] no metrics produced -- check input paths / split alignment.")
        return 1

    met = pd.DataFrame(all_metrics)
    g = ["clin_variant", "snp_variant", "mri_variant", "method", "cohort", "split"]
    agg = ["bacc", "macro_f1", "auc", "f1_sCN", "f1_pCN", "recall_sCN", "recall_pCN"]
    summary = met.groupby(g)[agg].agg(["mean", "std"]).round(4)
    summary.columns = [f"{a}_{b}" for a, b in summary.columns]
    summary = summary.reset_index()
    n_mean = met.groupby(g)["n"].mean().round(1).rename("n_mean").reset_index()
    summary = summary.merge(n_mean, on=g)

    met.to_csv(out_dir / "fusion_metrics_per_seed.csv", index=False)
    summary.to_csv(out_dir / "fusion_metrics.csv", index=False)
    pd.DataFrame(coverage).to_csv(out_dir / "coverage.csv", index=False)
    pd.DataFrame(all_preds).to_csv(out_dir / "fused_predictions.csv", index=False)
    pd.DataFrame(contingency).to_csv(out_dir / "failure_contingency.csv", index=False)
    render_leaderboard(summary, met, out_dir / "fusion_leaderboard_T1e_3way.png",
                       cohort="complete")
    render_leaderboard(summary, met, out_dir / "fusion_leaderboard_T1e_3way_present.png",
                       cohort="present")

    # stdout digest: best per cohort vs single-modality refs
    print("=" * 78)
    print("  T1e 3-way late fusion (clinical + SNP + MRI ViT-frozen)")
    print("=" * 78)
    print("\n[coverage] mean test-set presence per combo (over seeds):")
    cov = pd.DataFrame(coverage).groupby(
        ["clin_variant", "snp_variant", "mri_variant"]).mean().round(1)
    print(cov.to_string())
    for cohort in ("complete", "present"):
        s = summary[(summary["split"] == "test") & (summary["cohort"] == cohort)]
        if s.empty:
            continue
        print(f"\n[{cohort}-case TEST] top rows by bACC:")
        show = ["clin_variant", "snp_variant", "mri_variant", "method",
                "bacc_mean", "macro_f1_mean", "auc_mean", "n_mean"]
        print(s.sort_values("bacc_mean", ascending=False)[show]
              .head(12).to_string(index=False))
        refs = s[s["method"].isin(["clinical_only", "snp_only", "mri_only"])]
        print(f"  [refs] best single-modality bACC: "
              + ", ".join(f"{m}={refs[refs['method']==m]['bacc_mean'].max():.3f}"
                          for m in ("clinical_only", "snp_only", "mri_only")
                          if not refs[refs['method'] == m].empty))
    print(f"\n  wrote -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
