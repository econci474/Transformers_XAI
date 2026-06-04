"""
01_fuse_temporal_T1d.py
=======================
TEMPORAL late fusion for T1d (pMCI vs sMCI; sMCI=0, pMCI=1). Fuse a patient's CLINICAL
trajectory up to 1 year (CL_bl + CL_m06 + CL_m12) with the MRI 1-year scan (MRI_m12),
three ways (per the user's notes in outputs/{flat,class_wise_fusion}/*.txt):

  - FLAT          : one elastic-net (and a plain logistic-regression) stacker over ALL
                    probability blocks at once -- EN([CL_bl, CL_m06, CL_m12, MRI_m12]).
  - HIERARCHICAL  : an EN within clinical across timepoints, then averaged with MRI --
                    w*EN([CL_bl, CL_m06, CL_m12]) + (1-w)*MRI_m12  (w on VAL bACC / Youden's J).
  - CLASS_WISE    : a CONVEX weighted average of the 4 sources' per-class probabilities,
                    weights w1..w4 >= 0 summing to 1 (same weights both classes -> stays on the
                    simplex), tuned on VAL (bACC / Youden's J).  argmax decision.

Binary task: clinical/MRI each output a 2-vector (sMCI,pMCI) summing to 1. Target = the
conversion label (concurrent m12 label where present, else MRI m12 label).

REUSE, don't duplicate: imports the MRI-arm harness `../T2_classification/01_fuse_T2.py`
(loaders, fit_en, fit_lr, youden_j, metric_row, class prior, METRIC_COLS, formatting) and sets
it to 2 classes via h._set_classes. Meta-learners fit per-seed on VAL, reported on TEST; mean ±
std over seeds. Missing blocks: impute+indicator (train-fold prior + presence flag) / zero-fallback.

Run (per MRI model):
  python integration/T1d_classification/01_fuse_temporal_T1d.py \
      --mri-name BrainMVP_plusorig \
      --mri-csv-template "D:/ADNI_BIDS_project/derivatives/brainmvp_embeddings/T1d_binary/aug_plus_original/seed_{seed}/embeddings_seed_{seed}.csv"
Outputs (per mode, per MRI model):
  integration/T1d_classification/outputs/{flat,hierarchical,class_wise_fusion}/<mri-name>/
    fusion_metrics.csv, fusion_metrics_per_seed.csv, coverage.csv,
    fused_predictions.csv, failure_table.csv, fusion_table_<mode>.png
"""
from __future__ import annotations

import argparse
import importlib.util
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score

# --------------------------------------------------------------------------- #
# Import the MRI-arm harness (01_fuse_T2.py) and switch it to binary sMCI/pMCI.
# --------------------------------------------------------------------------- #
_HARNESS = Path(__file__).resolve().parent.parent / "T2_classification" / "01_fuse_T2.py"
_spec = importlib.util.spec_from_file_location("fuse_harness", _HARNESS)
h = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(h)
h._set_classes(["sMCI", "pMCI"])                     # binary; rebuilds CP/MP/METRIC_COLS

K = h.N_CLASSES                                      # 2
CLASS_NAMES = h.CLASS_NAMES                          # [sMCI, pMCI]
SEEDS_DEFAULT = (0, 1, 2)

CL_BL_ROOT = h.CLIN_ROOT                             # encoder_outputs_no_cdr_post_exclusion
CL_M06_ROOT = Path(r"D:/ADNI_BIDS_project/derivatives/encoder_outputs_no_cdr_post_exclusion_m06")
CL_M12_ROOT = Path(r"D:/ADNI_BIDS_project/derivatives/encoder_outputs_no_cdr_post_exclusion_m12")
CLIN_MODEL = "ModernBERT-large"                      # T1d notes' clinical pick
CLIN_TASK = "T1d_pmci_smci"
OUT_ROOT = Path(__file__).resolve().parent / "outputs"

CLIN_BLOCKS = [("bl", CL_BL_ROOT), ("m06", CL_M06_ROOT), ("m12", CL_M12_ROOT)]
CLIN_PREFIXES = [b[0] for b in CLIN_BLOCKS]          # bl, m06, m12
FLAT_PREFIXES = CLIN_PREFIXES + ["mri"]              # + MRI m12


# --------------------------------------------------------------------------- #
# Per-seed frame
# --------------------------------------------------------------------------- #
def _cls_cols(pre):
    return [f"{pre}{i}" for i in range(K)]


def _load_clin_block(seed, root, prefix):
    d = h.load_clinical(seed, root, CLIN_MODEL, h.CLIN_STRAT, task=CLIN_TASK)
    ren = {f"cp{i}": f"{prefix}{i}" for i in range(K)}
    ren.update({"y_clin": f"y_{prefix}", "split": f"split_{prefix}"})
    return d.rename(columns=ren)


def build_seed_frame(seed, mri_template):
    frames = [_load_clin_block(seed, root, pre) for pre, root in CLIN_BLOCKS]
    f = frames[0]
    for nxt in frames[1:]:
        f = f.merge(nxt, on="Patient_ID", how="outer")
    mri = h.load_mri(seed, mri_template, "m12").rename(
        columns={f"mp{i}": f"mri{i}" for i in range(K)})
    f = f.merge(mri, on="Patient_ID", how="outer")

    # target = conversion label: clinical m12 where present, else MRI m12 (both = conversion_group).
    # MRI y_mri == -1 means the m12 scan was ineligible (not-MCI at m12) -> treat as missing label.
    ymri_valid = f["y_mri"].where(f["y_mri"] >= 0) if "y_mri" in f.columns else np.nan
    f["y_true"] = f["y_m12"].where(f["y_m12"].notna(), ymri_valid)
    both = f["y_m12"].notna() & (f["y_mri"] >= 0) if "y_mri" in f.columns else f["y_m12"].notna() & False
    n_mismatch = int((f.loc[both, "y_m12"] != f.loc[both, "y_mri"]).sum())
    f = f[f["y_true"].notna()].copy()
    f["y_true"] = f["y_true"].astype(int)

    f["split"] = (f["split_m12"].fillna(f["split_bl"]).fillna(f["split_m06"])
                  .fillna(f["split_mri"]))
    for pre in CLIN_PREFIXES:
        f[f"{pre}_present"] = f[_cls_cols(pre)].notna().all(axis=1)
    # MRI "present" only where the m12 scan was eligible (visit_diag==MCI, y_mri>=0); the T1d MRI
    # prediction is OOD on a converted (non-MCI) m12 scan, so treat those as MRI-absent.
    f["mri_present"] = (f[_cls_cols("mri")].notna().all(axis=1)
                        & (f["y_mri"].fillna(-1) >= 0))
    return f, n_mismatch


def _probs(df, pre):
    return df[_cls_cols(pre)].to_numpy(dtype=float)


def _X_raw(df, prefixes):
    return np.column_stack([_probs(df, p) for p in prefixes])


def _X_zero(df, prefixes):
    return np.column_stack([np.nan_to_num(_probs(df, p), nan=0.0) for p in prefixes])


def _X_impute(df, prefixes, prior):
    blocks = []
    for p in prefixes:
        X = _probs(df, p).copy()
        X[np.isnan(X).any(axis=1)] = prior
        blocks.append(X)
    ind = np.column_stack([df[f"{p}_present"].to_numpy(float) for p in prefixes])
    return np.column_stack(blocks + [ind])


def _simplex_grid(n_blocks, step=0.1):
    k = int(round(1 / step))
    for comp in itertools.product(range(k + 1), repeat=n_blocks):
        if sum(comp) == k:
            yield np.array(comp, dtype=float) / k


def _score(name):
    return h.youden_j if name == "youden" else balanced_accuracy_score


def _X_present(df, prefixes):
    """PRESENT-ONLY impute: fill a missing block with the patient's present-block mean
    (+ per-block presence indicators). Uses only that patient's available sources."""
    blocks = [_probs(df, p) for p in prefixes]            # list of (n,K), nan where absent
    stacked = np.stack(blocks, 0)                          # (P, n, K)
    pres = ~np.isnan(stacked).any(axis=2)                  # (P, n)
    bf = np.nan_to_num(stacked)
    cnt = pres.sum(0)                                      # (n,)
    pmean = (bf * pres[:, :, None]).sum(0) / np.clip(cnt[:, None], 1, None)  # (n, K)
    filled = [np.where(pres[i][:, None], bf[i], pmean) for i in range(len(prefixes))]
    ind = np.column_stack([pres[i].astype(float) for i in range(len(prefixes))])
    return np.column_stack(filled + [ind])


def _fuse_present(df, prefixes, w):
    """Per-patient convex average over PRESENT sources, weights renormalised over present
    (missing source dropped). sum_{i present} w_i*P_i / sum_{i present} w_i -- stays on simplex."""
    n = len(df)
    num = np.zeros((n, K)); den = np.zeros((n, 1))
    for pre, wi in zip(prefixes, w):
        m = df[f"{pre}_present"].to_numpy(bool).astype(float)[:, None]
        num += wi * m * np.nan_to_num(_probs(df, pre))
        den += wi * m
    return num / np.clip(den, 1e-12, None)


def _fit_vt(fit_fn, Xv, yv, Xt, seed):
    """Fit on (Xv,yv); predict VAL and TEST in one fit. Returns (pred_v, pred_t, proba_t) or None."""
    res = fit_fn(Xv, yv, np.vstack([Xv, Xt]), seed)
    if res is None:
        return None
    pa, proba = res
    nv = len(Xv)
    return pa[:nv], pa[nv:], proba[nv:]


# --------------------------------------------------------------------------- #
# Per-seed evaluation
# --------------------------------------------------------------------------- #
def _row(rows, mode, variant, missing, yt, pt, proba_t, yv=None, pv=None, **extra):
    r = h.metric_row(np.asarray(yt, int), np.asarray(pt, int), proba_t)
    r["val_bacc"] = (float(balanced_accuracy_score(np.asarray(yv, int), np.asarray(pv, int)))
                     if yv is not None and len(yv) and len(np.unique(np.asarray(yv, int))) > 1
                     else float("nan"))
    r.update(mode=mode, variant=variant, missing=missing, **extra)
    rows.append(r)


def _stash(preds, mode, variant, missing, seed, df, pred, proba):
    for i, (_, row) in enumerate(df.reset_index(drop=True).iterrows()):
        rec = {"mode": mode, "variant": variant, "missing": missing, "seed": seed,
               "Patient_ID": row["Patient_ID"], "y_true": int(row["y_true"]),
               "pred": int(pred[i]),
               **{f"{p}_present": bool(row[f"{p}_present"]) for p in FLAT_PREFIXES}}
        for c, name in enumerate(CLASS_NAMES):
            rec[f"prob_{name}"] = float(proba[i, c])
        preds.append(rec)


def evaluate_seed(seed, f, prior, preds, fails):
    val, test = f[f.split == "val"], f[f.split == "test"]
    yv_t = lambda d: d["y_true"].to_numpy(int)

    # ---- shared references (each carries VAL preds for val_bacc) ----
    refs = []  # (variant, missing, yv, pv, yt, pt, proba_t)
    for name, flag in [("clinical_only", "m12_present"), ("mri_only", "mri_present")]:
        src = "m12" if name == "clinical_only" else "mri"
        rv, rt = val[val[flag]], test[test[flag]]
        if len(rt):
            Pv, Pt = _probs(rv, src), _probs(rt, src)
            refs.append((name, "-", yv_t(rv), Pv.argmax(1), yv_t(rt), Pt.argmax(1), Pt))

    # clinical-temporal EN over [bl,m06,m12] (impute+indicator) -- fit VAL, predict VAL+TEST.
    clin_v = val[val[[f"{p}_present" for p in CLIN_PREFIXES]].any(axis=1)]
    clin_t = test[test[[f"{p}_present" for p in CLIN_PREFIXES]].any(axis=1)]
    P_clin_v = P_clin_t = None
    if len(clin_v) >= 4 and len(clin_t):
        res = h.fit_en(_X_impute(clin_v, CLIN_PREFIXES, prior), yv_t(clin_v),
                       np.vstack([_X_impute(clin_v, CLIN_PREFIXES, prior),
                                  _X_impute(clin_t, CLIN_PREFIXES, prior)]), seed)
        if res is not None:
            _, proba_all = res
            P_clin_v, P_clin_t = proba_all[:len(clin_v)], proba_all[len(clin_v):]
            refs.append(("clinical_temporal_EN", "impute_indicator", yv_t(clin_v),
                         P_clin_v.argmax(1), yv_t(clin_t), P_clin_t.argmax(1), P_clin_t))

    all4_v = val[val[[f"{p}_present" for p in FLAT_PREFIXES]].all(axis=1)]
    all4_t = test[test[[f"{p}_present" for p in FLAT_PREFIXES]].all(axis=1)]
    any_v = val[val[[f"{p}_present" for p in FLAT_PREFIXES]].any(axis=1)]
    any_t = test[test[[f"{p}_present" for p in FLAT_PREFIXES]].any(axis=1)]

    def _emit_refs(rows, mode):
        for v, miss, yv, pv, yt, pt, pr in refs:
            _row(rows, mode, v, miss, yt, pt, pr, yv=yv, pv=pv)

    # ========================= FLAT =========================
    flat_rows = []
    _emit_refs(flat_rows, "flat")
    if len(all4_v) >= 4 and len(all4_t):
        for fit, vname in [(h.fit_en, "elastic_net"), (h.fit_lr, "lr")]:
            out = _fit_vt(fit, _X_raw(all4_v, FLAT_PREFIXES), yv_t(all4_v),
                          _X_raw(all4_t, FLAT_PREFIXES), seed)
            if out is not None:
                pv, pt, prt = out
                _row(flat_rows, "flat", vname, "complete_case", yv_t(all4_t), pt, prt,
                     yv=yv_t(all4_v), pv=pv)
                _stash(preds, "flat", vname, "complete_case", seed, all4_t, pt, prt)
    if len(any_v) >= 4 and len(any_t):
        for fill, miss in [(_X_impute, "impute_indicator"), (_X_present, "present_only")]:
            Xv = fill(any_v, FLAT_PREFIXES, prior) if fill is _X_impute else fill(any_v, FLAT_PREFIXES)
            Xt = fill(any_t, FLAT_PREFIXES, prior) if fill is _X_impute else fill(any_t, FLAT_PREFIXES)
            out = _fit_vt(h.fit_en, Xv, yv_t(any_v), Xt, seed)
            if out is not None:
                pv, pt, prt = out
                _row(flat_rows, "flat", "elastic_net", miss, yv_t(any_t), pt, prt,
                     yv=yv_t(any_v), pv=pv)
                _stash(preds, "flat", "elastic_net", miss, seed, any_t, pt, prt)
                if miss == "impute_indicator":
                    _failrows(fails, "flat", seed, any_t, P_clin_t, clin_t, pt)

    # ===================== HIERARCHICAL =====================
    hier_rows = []
    _emit_refs(hier_rows, "hierarchical")

    def _blend(w, Pc, mri, mask):
        out = Pc.copy()
        out[mask] = w * Pc[mask] + (1 - w) * mri[mask]
        return out

    if P_clin_t is not None:
        mri_t = np.nan_to_num(_probs(clin_t, "mri")); tmask = clin_t["mri_present"].to_numpy(bool)
        mri_v = np.nan_to_num(_probs(clin_v, "mri")); vmask = clin_v["mri_present"].to_numpy(bool)
        yt, yv = yv_t(clin_t), yv_t(clin_v)
        Pt_eq, Pv_eq = _blend(0.5, P_clin_t, mri_t, tmask), _blend(0.5, P_clin_v, mri_v, vmask)
        _row(hier_rows, "hierarchical", "avg", "equal", yt, Pt_eq.argmax(1), Pt_eq,
             yv=yv, pv=Pv_eq.argmax(1))
        _stash(preds, "hierarchical", "avg", "equal", seed, clin_t, Pt_eq.argmax(1), Pt_eq)
        _failrows(fails, "hierarchical", seed, clin_t, P_clin_t, clin_t, Pt_eq.argmax(1))
        if vmask.sum() >= 4:
            for miss, sname in [("tuned_bacc", "bacc"), ("tuned_youden", "youden")]:
                sc = _score(sname)
                best_w, best_s = 0.5, -np.inf
                for w in np.linspace(0, 1, 21):
                    s = sc(yv[vmask], (w * P_clin_v[vmask] + (1 - w) * mri_v[vmask]).argmax(1))
                    if s > best_s:
                        best_s, best_w = s, float(w)
                Pt = _blend(best_w, P_clin_t, mri_t, tmask)
                Pv = _blend(best_w, P_clin_v, mri_v, vmask)
                _row(hier_rows, "hierarchical", "avg", miss, yt, Pt.argmax(1), Pt,
                     yv=yv, pv=Pv.argmax(1), note=f"w_clin={best_w:.2f}")
                _stash(preds, "hierarchical", "avg", miss, seed, clin_t, Pt.argmax(1), Pt)
        # present-only clinical-temporal EN (comparison) + present-only avg(equal)
        outp = _fit_vt(h.fit_en, _X_present(clin_v, CLIN_PREFIXES), yv,
                       _X_present(clin_t, CLIN_PREFIXES), seed)
        if outp is not None:
            pvp, ptp, prtp = outp
            prvp = h.fit_en(_X_present(clin_v, CLIN_PREFIXES), yv,
                            _X_present(clin_v, CLIN_PREFIXES), seed)
            Pcv = prvp[1] if prvp is not None else P_clin_v
            _row(hier_rows, "hierarchical", "clinical_temporal_EN", "present_only",
                 yt, ptp, prtp, yv=yv, pv=pvp)
            Pt_p = _blend(0.5, prtp, mri_t, tmask); Pv_p = _blend(0.5, Pcv, mri_v, vmask)
            _row(hier_rows, "hierarchical", "avg", "present_only_equal", yt, Pt_p.argmax(1), Pt_p,
                 yv=yv, pv=Pv_p.argmax(1))
            _stash(preds, "hierarchical", "avg", "present_only_equal", seed, clin_t,
                   Pt_p.argmax(1), Pt_p)

    # ====================== CLASS_WISE ======================
    cw_rows = []
    _emit_refs(cw_rows, "class_wise")
    grid = list(_simplex_grid(len(FLAT_PREFIXES), step=0.1))
    # ---- PRIMARY: present-only convex avg (any-present cohort; per-patient weight renorm) ----
    if len(any_v) >= 4 and len(any_t):
        yv_a, yt_a = yv_t(any_v), yv_t(any_t)
        w_eq = np.ones(len(FLAT_PREFIXES)) / len(FLAT_PREFIXES)
        Pt = _fuse_present(any_t, FLAT_PREFIXES, w_eq); Pv = _fuse_present(any_v, FLAT_PREFIXES, w_eq)
        _row(cw_rows, "class_wise", "convex_present", "equal", yt_a, Pt.argmax(1), Pt,
             yv=yv_a, pv=Pv.argmax(1))
        _stash(preds, "class_wise", "convex_present", "equal", seed, any_t, Pt.argmax(1), Pt)
        for miss, sname in [("tuned_bacc", "bacc"), ("tuned_youden", "youden")]:
            sc = _score(sname)
            best_w, best_s = None, -np.inf
            for w in grid:
                s = sc(yv_a, _fuse_present(any_v, FLAT_PREFIXES, w).argmax(1))
                if s > best_s:
                    best_s, best_w = s, w
            Pt = _fuse_present(any_t, FLAT_PREFIXES, best_w)
            Pv = _fuse_present(any_v, FLAT_PREFIXES, best_w)
            _row(cw_rows, "class_wise", "convex_present", miss, yt_a, Pt.argmax(1), Pt,
                 yv=yv_a, pv=Pv.argmax(1), note="w=[" + ",".join(f"{x:.1f}" for x in best_w) + "]")
            _stash(preds, "class_wise", "convex_present", miss, seed, any_t, Pt.argmax(1), Pt)
            _failrows(fails, "class_wise", seed, any_t, P_clin_t, clin_t, Pt.argmax(1))
    # ---- COMPARISON: complete-case convex avg (all-4 present) ----
    if len(all4_v) >= 4 and len(all4_t):
        Bv = [_probs(all4_v, p) for p in FLAT_PREFIXES]
        Bt = [_probs(all4_t, p) for p in FLAT_PREFIXES]
        yv4, yt4 = yv_t(all4_v), yv_t(all4_t)
        _fuse = lambda B, w: sum(wi * Bi for wi, Bi in zip(w, B))
        w_eq = np.ones(len(FLAT_PREFIXES)) / len(FLAT_PREFIXES)
        Pt, Pv = _fuse(Bt, w_eq), _fuse(Bv, w_eq)
        _row(cw_rows, "class_wise", "convex_complete", "equal", yt4, Pt.argmax(1), Pt,
             yv=yv4, pv=Pv.argmax(1))
        for miss, sname in [("tuned_bacc", "bacc")]:
            sc = _score(sname)
            best_w, best_s = None, -np.inf
            for w in grid:
                s = sc(yv4, _fuse(Bv, w).argmax(1))
                if s > best_s:
                    best_s, best_w = s, w
            Pt, Pv = _fuse(Bt, best_w), _fuse(Bv, best_w)
            _row(cw_rows, "class_wise", "convex_complete", miss, yt4, Pt.argmax(1), Pt,
                 yv=yv4, pv=Pv.argmax(1), note="w=[" + ",".join(f"{x:.1f}" for x in best_w) + "]")

    return flat_rows + hier_rows + cw_rows


def _failrows(fails, mode, seed, test_df, P_clin_t, clin_t, fused_pred):
    if P_clin_t is None:
        return
    clin_pred = {pid: int(p) for pid, p in
                 zip(clin_t["Patient_ID"].tolist(), P_clin_t.argmax(1).tolist())}
    td = test_df.reset_index(drop=True)
    for i, (_, r) in enumerate(td.iterrows()):
        y = int(r["y_true"]); pid = r["Patient_ID"]
        cp = clin_pred.get(pid)
        mp = int(np.argmax(_probs(td.iloc[[i]], "mri")[0])) if r["mri_present"] else None
        c_ok = (cp == y) if cp is not None else None
        m_ok = (mp == y) if mp is not None else None
        if c_ok and m_ok:
            agree = "both_correct"
        elif c_ok and m_ok is False:
            agree = "clinical_only"
        elif m_ok and c_ok is False:
            agree = "mri_only"
        elif c_ok is False and m_ok is False:
            agree = "both_wrong"
        else:
            agree = "mri_absent" if mp is None else "na"
        fails.append({
            "mode": mode, "seed": seed, "Patient_ID": pid, "y_true": CLASS_NAMES[y],
            "clin_temporal_pred": CLASS_NAMES[cp] if cp is not None else "absent",
            "mri_pred": CLASS_NAMES[mp] if mp is not None else "absent",
            "fused_pred": CLASS_NAMES[int(fused_pred[i])],
            "clin_correct": c_ok, "mri_correct": m_ok,
            "fused_correct": bool(int(fused_pred[i]) == y), "agreement": agree,
        })


# --------------------------------------------------------------------------- #
# Aggregate + render
# --------------------------------------------------------------------------- #
AGG_COLS = ["val_bacc", "bacc", "macro_f1", "macro_auc"] + [f"f1_{n}" for n in CLASS_NAMES]

_REFS = [
    ("clinical_only", "-", "clinical-only [CL_m12] (ref)"),
    ("mri_only", "-", "MRI-only [m12] (ref)"),
    ("clinical_temporal_EN", "impute_indicator", "clinical-temporal EN [bl+m06+m12] (ref)"),
]
ROW_ORDER = {
    "flat": _REFS + [
        ("elastic_net", "complete_case", "flat EN (all 4 present)"),
        ("lr", "complete_case", "flat logistic regression (all 4)"),
        ("elastic_net", "impute_indicator", "flat EN (impute+indicator)"),
        ("elastic_net", "present_only", "flat EN (present-only)"),
    ],
    "hierarchical": _REFS + [
        ("clinical_temporal_EN", "present_only", "clinical-temporal EN [present-only] (ref)"),
        ("avg", "equal", "hierarchical avg (0.5/0.5)"),
        ("avg", "tuned_bacc", "hierarchical avg (VAL-tuned w, bACC)"),
        ("avg", "tuned_youden", "hierarchical avg (VAL-tuned w, Youden's J)"),
        ("avg", "present_only_equal", "hierarchical avg present-only (0.5/0.5)"),
    ],
    "class_wise": _REFS + [
        ("convex_present", "equal", "class-wise present-only (equal 1/4)"),
        ("convex_present", "tuned_bacc", "class-wise present-only (VAL-tuned w, bACC)"),
        ("convex_present", "tuned_youden", "class-wise present-only (VAL-tuned w, Youden's J)"),
        ("convex_complete", "equal", "class-wise complete-case (equal 1/4)"),
        ("convex_complete", "tuned_bacc", "class-wise complete-case (VAL-tuned w, bACC)"),
    ],
}
TITLE = {
    "flat": "FLAT: EN/LR([CL_bl, CL_m06, CL_m12, MRI_m12])",
    "hierarchical": "HIERARCHICAL: w*EN([CL_bl,CL_m06,CL_m12]) + (1-w)*MRI_m12",
    "class_wise": "CLASS-WISE: convex avg  sum_i w_i * P_source_i  (w>=0, sum w=1)",
}
FOOT = {
    "flat": [
        "FLAT = a SINGLE binary meta-learner over the stacked 4-block prob vector:",
        "  [CL_bl(sMCI,pMCI), CL_m06, CL_m12, MRI_m12]  (impute+indicator adds per-block presence flags).",
        "  EN = elastic-net logistic; LR = plain L2 logistic. Missing block -> TRAIN-fold prior + flag.",
    ],
    "hierarchical": [
        "HIERARCHICAL = stage-1 EN over clinical [CL_bl,CL_m06,CL_m12] (impute+indicator) = P_clin;",
        "  then P_fused = w*P_clin + (1-w)*MRI_m12 (MRI absent -> P_clin). w on VAL: 0.5, or grid",
        "  {0,0.05,..,1} maximising balanced accuracy, or Youden's J (sens+spec-1).",
    ],
    "class_wise": [
        "CLASS-WISE = convex weighted average of the 4 sources' per-class probs (all 4 present):",
        "  p(sMCI)_fused = sum_i w_i * p(sMCI)_source_i ;  p(pMCI)_fused = sum_i w_i * p(pMCI)_source_i ;",
        "  sources = [CL_bl, CL_m06, CL_m12, MRI_m12];  w_i>=0, sum w_i=1 (same weights both classes,",
        "  so the result is already a valid distribution -- no renormalisation). w tuned on VAL over a",
        "  simplex grid (step 0.1) by balanced accuracy or Youden's J; argmax decision.",
    ],
}


def render_png(summary, cov, mode, mri_name, clin_label, out_path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARN] matplotlib not installed -- skipping PNG.")
        return
    idx = {(r["variant"], r["missing"]): r for _, r in summary.iterrows()}
    body, numeric = [], []
    for variant, missing, label in ROW_ORDER[mode]:
        if (variant, missing) not in idx:
            continue
        r = idx[(variant, missing)]
        cells, nums = [], []
        for key, _, with_n in h.METRIC_COLS:
            n = r["n_test_mean"] if with_n else None
            cells.append(h._fmt(r.get(f"{key}_mean"), r.get(f"{key}_std"), n))
            nums.append(r.get(f"{key}_mean", np.nan))
        body.append([label] + cells)
        numeric.append(nums)
    if not body:
        print(f"  [skip] no rows for mode={mode}")
        return
    numeric = np.array(numeric, dtype=float)
    n_rows, n_cols = len(body), len(h.METRIC_COLS) + 1
    col_labels = ["Variant (fusion method)"] + [c[1] for c in h.METRIC_COLS]
    CHAR_W = 0.105
    col_chars = [max(len(str(x)) for x in [col_labels[j]] + [b[j] for b in body]) + 2
                 for j in range(n_cols)]
    total = sum(col_chars)
    fig, ax = plt.subplots(figsize=(total * CHAR_W, 0.55 + 0.32 * (n_rows + 2)))
    ax.axis("off")
    tab = ax.table(cellText=body, colLabels=col_labels, loc="upper center",
                   cellLoc="center", colLoc="center",
                   colWidths=[c / total for c in col_chars])
    tab.auto_set_font_size(False); tab.set_fontsize(9); tab.scale(1.0, 1.4)
    for j in range(len(h.METRIC_COLS)):
        col = numeric[:, j]
        if not np.all(np.isnan(col)):
            tab[int(np.nanargmax(col)) + 1, j + 1].set_text_props(weight="bold")
    for j in range(n_cols):
        tab[0, j].set_facecolor("#ECECEC"); tab[0, j].set_text_props(weight="bold")
    for i in range(n_rows + 1):
        tab[i, 0].set_text_props(ha="left")
    cov_line = ("  coverage (test, per seed): "
                + ", ".join(f"s{int(c.seed)} all4={int(c.n_test_all4)}/{int(c.n_test)}"
                            for c in cov.itertuples()))
    plt.title(f"T1d temporal late fusion -- {TITLE[mode]}\n"
              f"clinical={clin_label}   MRI={mri_name}   target=conversion (pMCI vs sMCI)   "
              f"(mean ± std across seeds, n)", pad=6, fontsize=10.5)
    foot = (["Meta-learners fit per-seed on VAL, reported on TEST. Classes: 0=sMCI, 1=pMCI.",
             "References: clinical-only=argmax CL_m12; MRI-only=argmax MRI_m12;",
             "  clinical-temporal EN = EN over [CL_bl,CL_m06,CL_m12] (impute+indicator) -- no MRI."]
            + FOOT[mode]
            + ["(n) after bACC = mean #TEST patients across seeds. '-' = not computable (too few/1-class VAL).",
               "", cov_line])
    ax.text(0.0, -0.04, "\n".join(foot), transform=ax.transAxes, ha="left", va="top",
            fontsize=7.5, family="monospace", linespacing=1.4)
    fig.savefig(out_path, dpi=180, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print(f"  PNG       : {out_path}")


def write_mode(mode, folder, met, cov, preds, fails, mri_name, clin_label, out_root):
    out_dir = out_root / folder / mri_name
    out_dir.mkdir(parents=True, exist_ok=True)
    m = met[met["mode"] == mode]
    summary = (m.groupby(["variant", "missing"])[AGG_COLS].agg(["mean", "std"]).round(4))
    summary.columns = [f"{a}_{b}" for a, b in summary.columns]
    summary = summary.reset_index()
    n_per = m.groupby(["variant", "missing"])["n"].mean().round(1)
    summary = summary.merge(n_per.rename("n_test_mean").reset_index(), on=["variant", "missing"])
    m.to_csv(out_dir / "fusion_metrics_per_seed.csv", index=False)
    summary.to_csv(out_dir / "fusion_metrics.csv", index=False)
    cov.to_csv(out_dir / "coverage.csv", index=False)
    preds[preds["mode"] == mode].to_csv(out_dir / "fused_predictions.csv", index=False)
    fails[fails["mode"] == mode].to_csv(out_dir / "failure_table.csv", index=False)
    render_png(summary, cov, mode, mri_name, clin_label, out_dir / f"fusion_table_{mode}.png")
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mri-name", required=True, help="Label/subfolder for the MRI model.")
    ap.add_argument("--mri-csv-template", required=True,
                    help="MRI T1d extraction CSV with {seed} (adni_viscode==m12 rows used).")
    ap.add_argument("--modes", nargs="+", default=["flat", "hierarchical", "class_wise"],
                    choices=["flat", "hierarchical", "class_wise"])
    ap.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS_DEFAULT))
    ap.add_argument("--out-root", type=Path, default=OUT_ROOT)
    args = ap.parse_args()

    all_metrics, preds, fails, cov = [], [], [], []
    for seed in args.seeds:
        f, n_mis = build_seed_frame(seed, args.mri_csv_template)
        prior = h.class_prior_from_train(f.rename(columns={"y_true": "y_clin"}))
        rows = evaluate_seed(seed, f, prior, preds, fails)
        all_metrics.extend([{**r, "seed": seed} for r in rows])
        te = f[f.split == "test"]
        cov.append({"seed": seed, "n_test": int(len(te)),
                    **{f"n_test_{p}": int(te[f"{p}_present"].sum()) for p in FLAT_PREFIXES},
                    "n_test_all4": int(te[[f"{p}_present" for p in FLAT_PREFIXES]].all(axis=1).sum())})
        if n_mis:
            print(f"  [seed {seed}] {n_mis} CL_m12 vs MRI_m12 label mismatch (using clinical).")

    met = pd.DataFrame(all_metrics)
    covdf = pd.DataFrame(cov)
    preds_df, fails_df = pd.DataFrame(preds), pd.DataFrame(fails)
    clin_label = f"{CLIN_MODEL}/{h.CLIN_STRAT}"
    out_root = args.out_root

    # class_wise mode -> user's folder name 'class_wise_fusion'
    folder = {"flat": "flat", "hierarchical": "hierarchical", "class_wise": "class_wise_fusion"}
    print("=" * 74)
    print(f"  T1d TEMPORAL late fusion  clinical={clin_label}  MRI={args.mri_name}")
    print("=" * 74)
    print("\n[coverage] test-set block presence per seed:")
    print(covdf.to_string(index=False))
    for mode in args.modes:
        summary = write_mode(mode, folder[mode], met, covdf, preds_df, fails_df,
                             args.mri_name, clin_label, out_root)
        print(f"\n[{mode}] mean over seeds (fit VAL, report TEST):")
        show = ["variant", "missing", "n_test_mean", "bacc_mean", "bacc_std", "macro_auc_mean"]
        print(summary[show].to_string(index=False))
        print(f"  wrote -> {out_root / folder[mode] / args.mri_name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
