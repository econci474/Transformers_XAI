"""
02_fuse_temporal_T2.py
======================
TEMPORAL late fusion for T2_multiclass (CN/MCI/AD): fuse a patient's CLINICAL trajectory up to
1 year (CL_bl + CL_m06 + CL_m12) with the MRI 1-year scan (MRI_m12), two ways:

  - FLAT         : one elastic-net stacker over ALL probability blocks at once
                   EN([CL_bl, CL_m06, CL_m12, MRI_m12]).
  - HIERARCHICAL : an EN within clinical across timepoints, then averaged with MRI
                   average( EN([CL_bl, CL_m06, CL_m12]) , MRI_m12 ).

Motivation: clinical assessments are far denser than MRI (clinical at bl/m06/m12; MRI essentially
m12 only), so the clinical-temporal EN exploits more evidence than a single timepoint. Target = the
CONCURRENT m12 diagnosis (Label_visit_diag at m12). Class index 0=CN/1=MCI/2=AD.

REUSE, don't duplicate: this imports the MRI arm's harness `01_fuse_T2.py` (loaders, EN fitter,
metric_row, class prior, formatting) and never edits it. Meta-learners fit per-seed on VAL, reported
on TEST; mean ± std over seeds. Missing blocks handled by impute+indicator (train-fold class prior +
per-block presence flag) and a zero-fallback comparison.

Inputs (local D:):
  CL_bl  : encoder_outputs_no_cdr_post_exclusion/<model>/T2_multiclass/seed_*/full_ft/embeddings/embeddings.parquet
  CL_m06 : encoder_outputs_no_cdr_post_exclusion_m06/<model>/T2_m06_multiclass/seed_*/full_ft/...
  CL_m12 : encoder_outputs_no_cdr_post_exclusion_m12/<model>/T2_m12_multiclass/seed_*/full_ft/...
  MRI_m12: brainmvp_embeddings/T2_multiclass/aug_stochastic/seed_{seed}/embeddings_seed_{seed}.csv  (adni_viscode==m12)

Run:
  python integration/T2_classification/02_fuse_temporal_T2.py --mode both \
      --mri-name BrainMVP_full_ft_stochastic
Outputs (per mode):
  integration/T2_classification/outputs/up_to_m12/{flat,hierarchical}/
    fusion_metrics.csv, fusion_metrics_per_seed.csv, coverage.csv,
    fused_predictions.csv, failure_table.csv, fusion_table_<mode>.png
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score

# --------------------------------------------------------------------------- #
# Import the MRI-arm harness (01_fuse_T2.py) as a module -- reuse, don't edit.
# --------------------------------------------------------------------------- #
_HARNESS = Path(__file__).resolve().parent / "01_fuse_T2.py"
_spec = importlib.util.spec_from_file_location("fuse_harness", _HARNESS)
h = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(h)

N_CLASSES = h.N_CLASSES
CLASS_NAMES = h.CLASS_NAMES
SEEDS_DEFAULT = (0, 1, 2)

CL_BL_ROOT = h.CLIN_ROOT
CL_M06_ROOT = Path(r"D:/ADNI_BIDS_project/derivatives/encoder_outputs_no_cdr_post_exclusion_m06")
CL_M12_ROOT = Path(r"D:/ADNI_BIDS_project/derivatives/encoder_outputs_no_cdr_post_exclusion_m12")
MRI_TEMPLATE_DEFAULT = ("D:/ADNI_BIDS_project/derivatives/brainmvp_embeddings/T2_multiclass/"
                        "aug_stochastic/seed_{seed}/embeddings_seed_{seed}.csv")
OUT_ROOT_DEFAULT = Path(__file__).resolve().parent / "outputs" / "up_to_m12"

# clinical timepoint blocks (prefix, encoder root, parquet task)
CLIN_BLOCKS = [
    ("bl",  CL_BL_ROOT,  "T2_multiclass"),
    ("m06", CL_M06_ROOT, "T2_m06_multiclass"),
    ("m12", CL_M12_ROOT, "T2_m12_multiclass"),
]
CLIN_PREFIXES = [b[0] for b in CLIN_BLOCKS]            # bl, m06, m12
FLAT_PREFIXES = CLIN_PREFIXES + ["mri"]               # + MRI m12


# --------------------------------------------------------------------------- #
# Per-seed frame: join the 4 blocks, anchor the m12 target, tag presence/split
# --------------------------------------------------------------------------- #
def _load_clin_block(seed, root, task, prefix):
    d = h.load_clinical(seed, root, h.CLIN_MODEL, h.CLIN_STRAT, task=task)
    return d.rename(columns={"cp0": f"{prefix}0", "cp1": f"{prefix}1", "cp2": f"{prefix}2",
                             "y_clin": f"y_{prefix}", "split": f"split_{prefix}"})


def build_seed_frame(seed, mri_template):
    frames = [_load_clin_block(seed, root, task, pre) for pre, root, task in CLIN_BLOCKS]
    f = frames[0]
    for nxt in frames[1:]:
        f = f.merge(nxt, on="Patient_ID", how="outer")
    mri = h.load_mri(seed, mri_template, "m12").rename(
        columns={"mp0": "mri0", "mp1": "mri1", "mp2": "mri2"})
    f = f.merge(mri, on="Patient_ID", how="outer")

    # target = concurrent m12 dx: clinical m12 label where present, else MRI m12 label
    f["y_true"] = f["y_m12"].where(f["y_m12"].notna(), f.get("y_mri"))
    both = f["y_m12"].notna() & f["y_mri"].notna()
    n_mismatch = int((f.loc[both, "y_m12"] != f.loc[both, "y_mri"]).sum())
    f = f[f["y_true"].notna()].copy()
    f["y_true"] = f["y_true"].astype(int)

    # canonical split (clinical blocks are byte-identical; coalesce, prefer clinical)
    f["split"] = (f["split_m12"].fillna(f["split_bl"]).fillna(f["split_m06"])
                  .fillna(f["split_mri"]))
    for pre in FLAT_PREFIXES:
        f[f"{pre}_present"] = f[[f"{pre}0", f"{pre}1", f"{pre}2"]].notna().all(axis=1)
    return f, n_mismatch


def _probs(df, pre):
    return df[[f"{pre}0", f"{pre}1", f"{pre}2"]].to_numpy(dtype=float)


def _X_raw(df, prefixes):
    return np.column_stack([_probs(df, p) for p in prefixes])


def _X_zero(df, prefixes):
    return np.column_stack([np.nan_to_num(_probs(df, p), nan=0.0) for p in prefixes])


def _X_impute(df, prefixes, prior):
    """Missing block probs -> train-fold class prior; append per-block presence indicators."""
    blocks = []
    for p in prefixes:
        X = _probs(df, p)
        miss = np.isnan(X).any(axis=1)
        X = X.copy()
        X[miss] = prior
        blocks.append(X)
    ind = np.column_stack([df[f"{p}_present"].to_numpy(float) for p in prefixes])
    return np.column_stack(blocks + [ind])


def _row(rows, mode, variant, missing, y, pred, proba, **extra):
    r = h.metric_row(np.asarray(y, int), np.asarray(pred, int), proba)
    r.update(mode=mode, variant=variant, missing=missing, **extra)
    rows.append(r)


def _stash(preds, mode, variant, missing, seed, df, pred, proba):
    for i, (_, row) in enumerate(df.reset_index(drop=True).iterrows()):
        preds.append({
            "mode": mode, "variant": variant, "missing": missing, "seed": seed,
            "Patient_ID": row["Patient_ID"], "y_true": int(row["y_true"]),
            "pred": int(pred[i]), "prob_CN": float(proba[i, 0]),
            "prob_MCI": float(proba[i, 1]), "prob_AD": float(proba[i, 2]),
            **{f"{p}_present": bool(row[f"{p}_present"]) for p in FLAT_PREFIXES},
        })


# --------------------------------------------------------------------------- #
# Per-seed evaluation
# --------------------------------------------------------------------------- #
def macro_youden_j(y_true, y_pred, n_classes=3):
    """Mean over classes of one-vs-rest Youden's J = sensitivity_c + specificity_c - 1."""
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)
    js = []
    for c in range(n_classes):
        tp = int(np.sum((y_true == c) & (y_pred == c)))
        fn = int(np.sum((y_true == c) & (y_pred != c)))
        tn = int(np.sum((y_true != c) & (y_pred != c)))
        fp = int(np.sum((y_true != c) & (y_pred == c)))
        sens = tp / (tp + fn) if (tp + fn) else 0.0
        spec = tn / (tn + fp) if (tn + fp) else 0.0
        js.append(sens + spec - 1.0)
    return float(np.mean(js))


def evaluate_seed(seed, f, prior, preds, fails):
    """Return list[metric dict] tagged with mode in {flat, hierarchical} + shared references."""
    rows = []
    val, test = f[f.split == "val"], f[f.split == "test"]
    yv_all, yt_all = val["y_true"].to_numpy(int), test["y_true"].to_numpy(int)

    # ---------- shared single-modality / temporal references (added to BOTH modes) ----------
    refs = []  # (variant, missing, y, pred, proba)
    # clinical_only = argmax CL_m12 (best single clinical timepoint, concurrent)
    tcm = test[test.m12_present]
    if len(tcm):
        P = _probs(tcm, "m12")
        refs.append(("clinical_only", "-", tcm["y_true"].to_numpy(int), P.argmax(1), P))
    # mri_only = argmax MRI_m12
    tmr = test[test.mri_present]
    if len(tmr):
        P = _probs(tmr, "mri")
        refs.append(("mri_only", "-", tmr["y_true"].to_numpy(int), P.argmax(1), P))

    # clinical_temporal_EN = stage-1 EN over [bl,m06,m12] (impute+indicator).
    # Fit on VAL (any clinical present), predict VAL+TEST in one fit (need val proba for tuning).
    clin_any_v = val[val[[f"{p}_present" for p in CLIN_PREFIXES]].any(axis=1)]
    clin_any_t = test[test[[f"{p}_present" for p in CLIN_PREFIXES]].any(axis=1)]
    P_clin_v = P_clin_t = None
    if len(clin_any_v) >= 4 and len(clin_any_t):
        Xv = _X_impute(clin_any_v, CLIN_PREFIXES, prior)
        Xt = _X_impute(clin_any_t, CLIN_PREFIXES, prior)
        res = h.fit_en(Xv, clin_any_v["y_true"].to_numpy(int),
                       np.vstack([Xv, Xt]), seed)
        if res is not None:
            _, proba_all = res
            P_clin_v, P_clin_t = proba_all[:len(Xv)], proba_all[len(Xv):]
            refs.append(("clinical_temporal_EN", "impute_indicator",
                         clin_any_t["y_true"].to_numpy(int), P_clin_t.argmax(1), P_clin_t))

    # ---------- FLAT: EN over [CL_bl, CL_m06, CL_m12, MRI_m12] ----------
    all4_v = val[val[[f"{p}_present" for p in FLAT_PREFIXES]].all(axis=1)]
    all4_t = test[test[[f"{p}_present" for p in FLAT_PREFIXES]].all(axis=1)]
    any_t = test[test[[f"{p}_present" for p in FLAT_PREFIXES]].any(axis=1)]
    flat_rows = []
    for v, miss, y, pred, proba in refs:
        _row(flat_rows, "flat", v, miss, y, pred, proba)
    # complete-case
    if len(all4_v) >= 4 and len(all4_t):
        res = h.fit_en(_X_raw(all4_v, FLAT_PREFIXES), all4_v["y_true"].to_numpy(int),
                       _X_raw(all4_t, FLAT_PREFIXES), seed)
        if res is not None:
            pred, proba = res
            _row(flat_rows, "flat", "elastic_net", "complete_case",
                 all4_t["y_true"].to_numpy(int), pred, proba)
            _stash(preds, "flat", "elastic_net", "complete_case", seed, all4_t, pred, proba)
        # zero-fallback: same complete-case fit, applied to any-present test (missing blocks -> 0)
        if len(any_t):
            res = h.fit_en(_X_raw(all4_v, FLAT_PREFIXES), all4_v["y_true"].to_numpy(int),
                           _X_zero(any_t, FLAT_PREFIXES), seed)
            if res is not None:
                pred, proba = res
                _row(flat_rows, "flat", "elastic_net", "zero_fallback",
                     any_t["y_true"].to_numpy(int), pred, proba)
                _stash(preds, "flat", "elastic_net", "zero_fallback", seed, any_t, pred, proba)
    # impute + indicator (fit on any-present val)
    any_v = val[val[[f"{p}_present" for p in FLAT_PREFIXES]].any(axis=1)]
    if len(any_v) >= 4 and len(any_t):
        res = h.fit_en(_X_impute(any_v, FLAT_PREFIXES, prior), any_v["y_true"].to_numpy(int),
                       _X_impute(any_t, FLAT_PREFIXES, prior), seed)
        if res is not None:
            pred, proba = res
            _row(flat_rows, "flat", "elastic_net", "impute_indicator",
                 any_t["y_true"].to_numpy(int), pred, proba)
            _stash(preds, "flat", "elastic_net", "impute_indicator", seed, any_t, pred, proba)
            _failrows(fails, "flat", seed, any_t, P_clin_t, clin_any_t, pred)

    # ---------- HIERARCHICAL: average(P_clin, MRI_m12) ----------
    hier_rows = []
    for v, miss, y, pred, proba in refs:
        _row(hier_rows, "hierarchical", v, miss, y, pred, proba)
    if P_clin_t is not None:
        # align MRI probs onto the clinical-any test rows
        mri_t = _probs(clin_any_t, "mri")
        tmask = clin_any_t["mri_present"].to_numpy(bool)
        yt = clin_any_t["y_true"].to_numpy(int)

        def _blend(w, P_clin, mri, mask):
            out = P_clin.copy()
            out[mask] = w * P_clin[mask] + (1 - w) * mri[mask]
            return out

        # equal-weight average (headline)
        Pt_eq = _blend(0.5, P_clin_t, np.nan_to_num(mri_t), tmask)
        _row(hier_rows, "hierarchical", "avg", "equal", yt, Pt_eq.argmax(1), Pt_eq)
        _stash(preds, "hierarchical", "avg", "equal", seed, clin_any_t, Pt_eq.argmax(1), Pt_eq)
        _failrows(fails, "hierarchical", seed, clin_any_t, P_clin_t, clin_any_t, Pt_eq.argmax(1))

        # weight tuned on VAL (rows where both clinical-EN proba and MRI present), two criteria
        if P_clin_v is not None:
            vmask = clin_any_v["mri_present"].to_numpy(bool)
            if vmask.sum() >= 4:
                mri_v = np.nan_to_num(_probs(clin_any_v, "mri"))
                yv = clin_any_v["y_true"].to_numpy(int)

                def _tune(score_fn):
                    best_w, best_s = 0.5, -np.inf
                    for w in np.linspace(0, 1, 21):
                        pr = (w * P_clin_v[vmask] + (1 - w) * mri_v[vmask]).argmax(1)
                        s = score_fn(yv[vmask], pr)
                        if s > best_s:
                            best_s, best_w = s, float(w)
                    return best_w

                for miss, score_fn in [
                        ("tuned_bacc", balanced_accuracy_score),
                        ("tuned_youden", lambda a, b: macro_youden_j(a, b))]:
                    w = _tune(score_fn)
                    Pt = _blend(w, P_clin_t, np.nan_to_num(mri_t), tmask)
                    _row(hier_rows, "hierarchical", "avg", miss, yt, Pt.argmax(1), Pt,
                         note=f"w_clin={w:.2f}")
                    _stash(preds, "hierarchical", "avg", miss, seed, clin_any_t, Pt.argmax(1), Pt)

    rows.extend(flat_rows)
    rows.extend(hier_rows)
    return rows


def _failrows(fails, mode, seed, test_df, P_clin_t, clin_any_t, fused_pred):
    """Per TEST patient: clinical-temporal-EN vs MRI vs fused correctness + agreement."""
    if P_clin_t is None:
        return
    # map Patient_ID -> clinical-temporal-EN pred (computed on clin_any_t order)
    clin_pred = {pid: int(p) for pid, p in
                 zip(clin_any_t["Patient_ID"].tolist(), P_clin_t.argmax(1).tolist())}
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
AGG_COLS = ["bacc", "macro_f1", "macro_auc", "f1_CN", "f1_MCI", "f1_AD"]

ROW_ORDER = {
    "flat": [
        ("clinical_only", "-", "clinical-only [CL_m12] (ref)"),
        ("mri_only", "-", "MRI-only [m12] (ref)"),
        ("clinical_temporal_EN", "impute_indicator", "clinical-temporal EN [bl+m06+m12] (ref)"),
        ("elastic_net", "complete_case", "flat EN (all 4 present)"),
        ("elastic_net", "zero_fallback", "flat EN (zero-fallback)"),
        ("elastic_net", "impute_indicator", "flat EN (impute+indicator)"),
    ],
    "hierarchical": [
        ("clinical_only", "-", "clinical-only [CL_m12] (ref)"),
        ("mri_only", "-", "MRI-only [m12] (ref)"),
        ("clinical_temporal_EN", "impute_indicator", "clinical-temporal EN [bl+m06+m12] (ref)"),
        ("avg", "equal", "hierarchical avg (0.5/0.5)"),
        ("avg", "tuned_bacc", "hierarchical avg (VAL-tuned w, bACC)"),
        ("avg", "tuned_youden", "hierarchical avg (VAL-tuned w, Youden's J)"),
    ],
}
TITLE = {
    "flat": "FLAT: EN([CL_bl, CL_m06, CL_m12, MRI_m12])",
    "hierarchical": "HIERARCHICAL: average(EN([CL_bl, CL_m06, CL_m12]), MRI_m12)",
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
    cov_line = ("  coverage (test patients with all 4 blocks / total target, per seed): "
                + ", ".join(f"s{int(c.seed)}={int(c.n_test_all4)}/{int(c.n_test)}"
                            for c in cov.itertuples()))
    plt.title(f"T2 temporal late fusion -- {TITLE[mode]}\n"
              f"clinical={clin_label}   MRI={mri_name}   target=concurrent m12 dx   "
              f"(mean ± std across seeds, n)", pad=6, fontsize=10.5)
    foot_common = [
        "Meta-learners fit per-seed on VAL, reported on TEST. Class index 0=CN, 1=MCI, 2=AD.",
        "References: clinical-only=argmax CL_m12; MRI-only=argmax MRI_m12;",
        "  clinical-temporal EN = elastic net over [CL_bl,CL_m06,CL_m12] (impute+indicator) -- no MRI.",
    ]
    if mode == "flat":
        foot_mode = [
            "FLAT = a SINGLE multinomial elastic-net over the 12-dim T2-probability vector:",
            "  [CL_bl(CN,MCI,AD), CL_m06(CN,MCI,AD), CL_m12(CN,MCI,AD), MRI_m12(CN,MCI,AD)]",
            "  (16-dim with per-block presence flags in the impute+indicator variant); output = 3-class",
            "  softmax, joint over all classes (NOT class-wise).  All blocks are T2 3-class probabilities.",
            "Missing block -> TRAIN-fold class prior + presence flag (impute+indicator) / 0 (zero-fallback).",
        ]
    else:  # hierarchical
        foot_mode = [
            "HIERARCHICAL = average( clinical-temporal EN , MRI_m12 ) per class, then argmax;",
            "  P_fused = w*P_clin + (1-w)*MRI_m12 where MRI present, else P_clin (MRI absent -> clinical).",
            "  0.5/0.5 = equal weight;  VAL-tuned w = w on the 21-pt grid {0,0.05,..,1} maximising, on VAL,",
            "  balanced accuracy (bACC) or macro Youden's J = mean_c[sensitivity_c + specificity_c - 1] (CN/MCI/AD).",
        ]
    foot_lines = foot_common + foot_mode + ["(n) after bACC = mean #TEST patients across seeds."]
    ax.text(0.0, -0.04, "\n".join(foot_lines + ["", cov_line]), transform=ax.transAxes,
            ha="left", va="top", fontsize=7.5, family="monospace", linespacing=1.4)
    fig.savefig(out_path, dpi=180, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print(f"  PNG       : {out_path}")


def write_mode(mode, met, cov, preds, fails, mri_name, clin_label, out_root):
    out_dir = out_root / mode
    out_dir.mkdir(parents=True, exist_ok=True)
    m = met[met["mode"] == mode]
    summary = (m.groupby(["variant", "missing"])[AGG_COLS].agg(["mean", "std"]).round(4))
    summary.columns = [f"{a}_{b}" for a, b in summary.columns]
    summary = summary.reset_index()
    n_per = m.groupby(["variant", "missing"])["n"].mean().round(1)
    summary = summary.merge(n_per.rename("n_test_mean").reset_index(),
                            on=["variant", "missing"])
    m.to_csv(out_dir / "fusion_metrics_per_seed.csv", index=False)
    summary.to_csv(out_dir / "fusion_metrics.csv", index=False)
    cov.to_csv(out_dir / "coverage.csv", index=False)
    pm = preds[preds["mode"] == mode]
    pm.to_csv(out_dir / "fused_predictions.csv", index=False)
    fm = fails[fails["mode"] == mode]
    fm.to_csv(out_dir / "failure_table.csv", index=False)
    render_png(summary, cov, mode, mri_name, clin_label, out_dir / f"fusion_table_{mode}.png")
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["flat", "hierarchical", "both"], default="both")
    ap.add_argument("--mri-name", default="BrainMVP_full_ft_stochastic")
    ap.add_argument("--mri-csv-template", default=MRI_TEMPLATE_DEFAULT)
    ap.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS_DEFAULT))
    ap.add_argument("--out-root", type=Path, default=OUT_ROOT_DEFAULT)
    args = ap.parse_args()

    all_metrics, preds, fails, cov = [], [], [], []
    for seed in args.seeds:
        f, n_mis = build_seed_frame(seed, args.mri_csv_template)
        prior = h.class_prior_from_train(
            f.rename(columns={"y_true": "y_clin"}))  # train-fold class prior from this cohort
        rows = evaluate_seed(seed, f, prior, preds, fails)
        all_metrics.extend([{**r, "seed": seed} for r in rows])
        te = f[f.split == "test"]
        cov.append({
            "seed": seed, "n_test": int(len(te)),
            **{f"n_test_{p}": int(te[f"{p}_present"].sum()) for p in FLAT_PREFIXES},
            "n_test_all4": int(te[[f"{p}_present" for p in FLAT_PREFIXES]].all(axis=1).sum()),
        })
        if n_mis:
            print(f"  [seed {seed}] {n_mis} patients with CL_m12 vs MRI_m12 label mismatch "
                  f"(using clinical label).")

    met = pd.DataFrame(all_metrics)
    covdf = pd.DataFrame(cov)
    preds_df = pd.DataFrame(preds)
    fails_df = pd.DataFrame(fails)
    clin_label = f"{h.CLIN_MODEL}/{h.CLIN_STRAT}"

    modes = ["flat", "hierarchical"] if args.mode == "both" else [args.mode]
    print("=" * 74)
    print(f"  T2 TEMPORAL late fusion (up to m12)  clinical={clin_label}  MRI={args.mri_name}")
    print("=" * 74)
    print("\n[coverage] test-set block presence per seed:")
    print(covdf.to_string(index=False))
    for mode in modes:
        summary = write_mode(mode, met, covdf, preds_df, fails_df,
                             args.mri_name, clin_label, args.out_root)
        print(f"\n[{mode}] mean over seeds (fit VAL, report TEST):")
        show = ["variant", "missing", "n_test_mean", "bacc_mean", "bacc_std",
                "macro_f1_mean", "macro_auc_mean"]
        print(summary[show].to_string(index=False))
        print(f"  wrote -> {args.out_root / mode}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
