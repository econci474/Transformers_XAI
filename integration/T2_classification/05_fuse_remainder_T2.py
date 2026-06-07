"""
05_fuse_remainder_T2.py
=======================
T2 (CN/MCI/AD) via the REMAINDER-PRODUCT detector decomposition (the recreated "#11" comparison),
now WITH an MRI-only reference. Distinct from the class-wise fusion (03): here the 3-class is built
PURELY from two longitudinal binary detectors, with MCI as the residual (no clinical T2 in the fused
prediction; T2 appears only as reference rows).

Detectors (longitudinal binary elastic-net logistic, class-vs-rest, fit on VAL):
  CN-det:  CL      = EN([ P(CN) clinT1 @bl , @m06 , @m12 ])
           CL+MRI  = EN([ P(CN) clinT1 @bl , @m06 , @m12 , P(CN) BrainMVP-T1 @m12 ])
  AD-det:  CL      = EN([ P(AD) clinT1b @bl , @m06 , @m12 ])
           CL+MRI  = EN([ P(AD) clinT1b @bl , @m06 , @m12 , P(AD) BrainMVP-T1b @m12 ])
Compose to 3-class [CN, MCI, AD], argmax:
  structured : P = [ CN , (1-CN)(1-AD) , AD ]  then RENORMALISED to sum 1   (MCI = remainder).
  EN         : multinomial elastic-net over [ CN , AD , (1-CN)(1-AD) ], fit on VAL.
4 fusion variants = {structured, EN} x {CL, CL+MRI}.

References: clinical_only = argmax CL_m12 T2 3-class; temporal_CL_T2_EN = multinomial EN over CL T2
[bl,m06,m12]; MRI-only = the two MRI detectors composed alone, [t1,(1-t1)(1-t1b),t1b] renormalised
(BrainMVP T1+T1b, stochastic aug, @m12; MRI-present patients, n~44).

MRI model throughout = BrainMVP (uniformer, full_ft, stochastic aug). Clinical T1/T1b/T2 probs taken
DIRECTLY from the dedicated encoders. All metrics on TEST (fit on VAL); target = concurrent m12 dx;
cohort = TEST with an m12 visit (n=53). mean +/- std over seeds 0/1/2.

Output: integration/T2_classification/outputs/up_to_m12/detectors/remainder_product_fusion/
  fusion_table_detectors.png, fusion_metrics{,_per_seed}.csv, coverage.csv, fused_predictions.csv,
  failure_table.csv, failure_summary.csv, <method>_<src>/ {failure_table.csv, failure_seed{n}.png, failure_summary.png}
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV

_HARNESS = Path(__file__).resolve().parent / "01_fuse_T2.py"
_spec = importlib.util.spec_from_file_location("fuse_harness", _HARNESS)
h = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(h)

N_CLASSES, CLASS_NAMES = h.N_CLASSES, h.CLASS_NAMES
SEEDS_DEFAULT = (0, 1, 2)

DERIV = Path(r"D:/ADNI_BIDS_project/derivatives")
CL_ROOT = {"bl": h.CLIN_ROOT,
           "m06": DERIV / "encoder_outputs_no_cdr_post_exclusion_m06",
           "m12": DERIV / "encoder_outputs_no_cdr_post_exclusion_m12"}
T2 = {"model": "BioClinical-ModernBERT-large",
      "task": {"bl": "T2_multiclass", "m06": "T2_m06_multiclass", "m12": "T2_m12_multiclass"}}
T1 = {"model": "BioClinical-ModernBERT-base",
      "task": {"bl": "T1_binary", "m06": "T1_m06_binary", "m12": "T1_m12_binary"}}
T1B = {"model": "ModernBERT-base",
       "task": {"bl": "T1b_cnmci_ad", "m06": "T1b_m06_cnmci_ad", "m12": "T1b_m12_cnmci_ad"}}
MRI_T1_TMPL = ("D:/ADNI_BIDS_project/derivatives/brainmvp_embeddings/T1_binary/"
               "aug_stochastic/seed_{seed}/embeddings_seed_{seed}.csv")
MRI_T1B_TMPL = ("D:/ADNI_BIDS_project/derivatives/brainmvp_embeddings/T1b_binary/"
                "aug_stochastic/seed_{seed}/embeddings_seed_{seed}.csv")
MRI_NAME_DEFAULT = "BrainMVP T1+T1b stoch @m12"   # --mri-name; namespaces the MRI-only ref label
TPS = ["bl", "m06", "m12"]
OUT_DIR_DEFAULT = (Path(__file__).resolve().parent / "outputs" / "up_to_m12"
                   / "detectors" / "remainder_product_fusion")
SRC_SLUG = {"CL": "clinical_longitudinal", "CL+MRI": "clinical_longitudinal_plus_mri"}


def cn_cols(src):
    return [f"cn_{tp}" for tp in TPS] + (["mri_cn"] if src == "CL+MRI" else [])


def ad_cols(src):
    return [f"ad_{tp}" for tp in TPS] + (["mri_ad"] if src == "CL+MRI" else [])


# --------------------------------------------------------------------------- #
# Loaders (same as 03)
# --------------------------------------------------------------------------- #
def load_bin_prob(seed, tp, spec, pos_col, out):
    p = (CL_ROOT[tp] / spec["model"] / spec["task"][tp] / f"seed_{seed}" / "full_ft"
         / "embeddings" / "embeddings.parquet")
    d = pd.read_parquet(p, columns=["Patient_ID", "in_task_cohort", pos_col])
    return d[d["in_task_cohort"]][["Patient_ID", pos_col]].rename(columns={pos_col: out})


def load_t2_block(seed, tp):
    d = h.load_clinical(seed, CL_ROOT[tp], T2["model"], "full_ft", task=T2["task"][tp])
    return d.rename(columns={"cp0": f"t2_{tp}_0", "cp1": f"t2_{tp}_1", "cp2": f"t2_{tp}_2",
                             "y_clin": "y_m12" if tp == "m12" else f"_y_{tp}",
                             "split": "split" if tp == "m12" else f"_split_{tp}"})


def build_frame(seed, t1_tmpl=MRI_T1_TMPL, t1b_tmpl=MRI_T1B_TMPL):
    f = load_t2_block(seed, "m12")[["Patient_ID", "split", "y_m12",
                                    "t2_m12_0", "t2_m12_1", "t2_m12_2"]].copy()
    f["y_true"] = f["y_m12"].astype(int)
    for tp in ("bl", "m06"):
        f = f.merge(load_t2_block(seed, tp)[["Patient_ID", f"t2_{tp}_0", f"t2_{tp}_1", f"t2_{tp}_2"]],
                    on="Patient_ID", how="left")
    for tp in TPS:
        f = f.merge(load_bin_prob(seed, tp, T1, "prob_0", f"cn_{tp}"), on="Patient_ID", how="left")
        f = f.merge(load_bin_prob(seed, tp, T1B, "prob_1", f"ad_{tp}"), on="Patient_ID", how="left")
    f = f.merge(h.load_mri_detector(seed, t1_tmpl, "m12", "prob_class_0", "mri_cn"),
                on="Patient_ID", how="left")
    f = f.merge(h.load_mri_detector(seed, t1b_tmpl, "m12", "prob_class_1", "mri_ad"),
                on="Patient_ID", how="left")
    return f


def _bin_X(df, cols, impute=0.5):
    X = df[cols].copy()
    pres = X.notna().astype(float)
    return np.column_stack([X.fillna(impute).to_numpy(float), pres.to_numpy(float)])


def fit_binary_en(Xv, yv, Xall, seed):
    if len(np.unique(yv)) < 2:
        return None
    counts = np.bincount(yv, minlength=2)
    n_splits = min(3, int(counts[counts > 0].min()))
    if n_splits < 2:
        clf = LogisticRegression(penalty="elasticnet", solver="saga", C=1.0, l1_ratio=0.5,
                                 max_iter=5000, random_state=seed).fit(Xv, yv)
    else:
        clf = GridSearchCV(
            LogisticRegression(penalty="elasticnet", solver="saga", max_iter=5000, random_state=seed),
            {"C": [0.1, 1.0, 10.0], "l1_ratio": [0.1, 0.5, 0.9]}, scoring="balanced_accuracy",
            cv=StratifiedKFold(n_splits, shuffle=True, random_state=seed), n_jobs=-1).fit(Xv, yv)
    classes = clf.best_estimator_.classes_ if hasattr(clf, "best_estimator_") else clf.classes_
    return clf.predict_proba(Xall)[:, list(classes).index(1)]


def detector(val, test, cols, pos_class, seed):
    yv = (val["y_true"].to_numpy(int) == pos_class).astype(int)
    Xv, Xt = _bin_X(val, cols), _bin_X(test, cols)
    p = fit_binary_en(yv=yv, Xv=Xv, Xall=np.vstack([Xv, Xt]), seed=seed)
    return (None, None) if p is None else (p[:len(Xv)], p[len(Xv):])


def structured(cn, ad):
    P = np.clip(np.stack([cn, (1 - cn) * (1 - ad), ad], axis=1), 1e-9, None)
    return P / P.sum(axis=1, keepdims=True)


def _t2_impute(df, prior):
    blocks = []
    for tp in TPS:
        X = df[[f"t2_{tp}_0", f"t2_{tp}_1", f"t2_{tp}_2"]].to_numpy(float)
        miss = np.isnan(X).any(axis=1)
        X = X.copy(); X[miss] = prior
        blocks.append(X)
    ind = np.column_stack([df[f"t2_{tp}_0"].notna().to_numpy(float) for tp in TPS])
    return np.column_stack(blocks + [ind])


def variant_rows(mri_name=MRI_NAME_DEFAULT):
    out = [("clinical_only", "CL_m12", "REF  clinical T2 3-class @ m12 (argmax)"),
           ("temporal_CL_T2_EN", "CL bl+m06+m12", "REF  clinical T2 3-class EN @ bl+m06+m12"),
           ("MRI_only", mri_name,
            f"REF  MRI-only [{mri_name}]  (struct: [t1, (1-t1)(1-t1b), t1b])")]
    for method in ("structured", "EN"):
        for src in ("CL", "CL+MRI"):
            out.append((f"detector_{method}", src,
                        f"detector-{method} [CN<-clinT1, AD<-clinT1b @bl+m06+m12{' + BrainMVP T1/T1b' if src=='CL+MRI' else ''}]"))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS_DEFAULT))
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR_DEFAULT)
    ap.add_argument("--mri-name", default=MRI_NAME_DEFAULT,
                    help="MRI-model label for the MRI-only ref source (e.g. 'BrainDINO T1+T1b frozen @m12').")
    ap.add_argument("--mri-t1-template", default=MRI_T1_TMPL,
                    help="CSV template (with {seed}) for the T1/CN detector probs (prob_class_0).")
    ap.add_argument("--mri-t1b-template", default=MRI_T1B_TMPL,
                    help="CSV template (with {seed}) for the T1b/AD detector probs (prob_class_1).")
    args = ap.parse_args()

    rows, preds, cov, fails = [], [], [], []
    for seed in args.seeds:
        f = build_frame(seed, args.mri_t1_template, args.mri_t1b_template)
        prior = h.class_prior_from_train(f.rename(columns={"y_true": "y_clin"}))
        val, test = f[f.split == "val"], f[f.split == "test"]
        yv = val["y_true"].to_numpy(int)
        cpt = test[["t2_m12_0", "t2_m12_1", "t2_m12_2"]].to_numpy(float)

        def add(variant, source, pred, proba, df=None):
            d = test if df is None else df
            r = h.metric_row(d["y_true"].to_numpy(int), pred, proba)
            r.update(variant=variant, source=source, seed=seed)
            rows.append(r)
            for i, (_, rr) in enumerate(d.reset_index(drop=True).iterrows()):
                preds.append({"variant": variant, "source": source, "seed": seed,
                              "Patient_ID": rr["Patient_ID"], "y_true": int(rr["y_true"]),
                              "pred": int(pred[i]), "prob_CN": float(proba[i, 0]),
                              "prob_MCI": float(proba[i, 1]), "prob_AD": float(proba[i, 2])})

        def add_fail(method, src, fused_pred):
            tr = test.reset_index(drop=True)
            for i in range(len(tr)):
                y = int(tr.loc[i, "y_true"]); fu = int(fused_pred[i])
                rec = {"method": method, "src": src, "seed": seed,
                       "Patient_ID": tr.loc[i, "Patient_ID"], "y_true": CLASS_NAMES[y],
                       "fused_pred": CLASS_NAMES[fu], "fused_correct": fu == y}
                for tp in TPS:
                    v = tr.loc[i, [f"t2_{tp}_0", f"t2_{tp}_1", f"t2_{tp}_2"]].to_numpy(float)
                    if np.isnan(v).any():
                        rec[f"clinT2_{tp}_pred"], rec[f"clinT2_{tp}_correct"] = "absent", np.nan
                    else:
                        p = int(v.argmax())
                        rec[f"clinT2_{tp}_pred"], rec[f"clinT2_{tp}_correct"] = CLASS_NAMES[p], bool(p == y)
                base_ok = rec["clinT2_m12_correct"] is True
                rec["agreement"] = ("both_correct" if base_ok and fu == y else
                                    "base_better" if base_ok else
                                    "fused_better" if fu == y else "both_wrong")
                fails.append(rec)

        # references
        add("clinical_only", "CL_m12", cpt.argmax(1), cpt)
        base = h.fit_en(_t2_impute(val, prior), yv, _t2_impute(test, prior), seed)
        if base is not None:
            add("temporal_CL_T2_EN", "CL bl+m06+m12", base[0], base[1])
        tm = test[test["mri_cn"].notna() & test["mri_ad"].notna()]
        if len(tm):
            mcn, mad = tm["mri_cn"].to_numpy(float), tm["mri_ad"].to_numpy(float)
            add("MRI_only", args.mri_name, structured(mcn, mad).argmax(1),
                structured(mcn, mad), df=tm)

        # detector variants (remainder-product)
        for src in ("CL", "CL+MRI"):
            cn_v, cn_t = detector(val, test, cn_cols(src), 0, seed)
            ad_v, ad_t = detector(val, test, ad_cols(src), 2, seed)
            if cn_t is None or ad_t is None:
                continue
            Ps = structured(cn_t, ad_t)
            add("detector_structured", src, Ps.argmax(1), Ps)
            add_fail("structured", src, Ps.argmax(1))
            Xc_v = np.column_stack([cn_v, ad_v, (1 - cn_v) * (1 - ad_v)])
            Xc_t = np.column_stack([cn_t, ad_t, (1 - cn_t) * (1 - ad_t)])
            en = h.fit_en(Xc_v, yv, Xc_t, seed)
            if en is not None:
                add("detector_EN", src, en[0], en[1])
                add_fail("EN", src, en[0])

        cov.append({"seed": seed, "n_test": int(len(test)),
                    "n_test_cn_temporal": int(test[[f"cn_{tp}" for tp in TPS]].notna().any(axis=1).sum()),
                    "n_test_mri_cn": int(test["mri_cn"].notna().sum())})

    met = pd.DataFrame(rows)
    covdf = pd.DataFrame(cov)
    AGG = ["bacc", "macro_f1", "macro_auc", "f1_CN", "f1_MCI", "f1_AD"]
    summary = met.groupby(["variant", "source"], sort=False)[AGG].agg(["mean", "std"]).round(4)
    summary.columns = [f"{a}_{b}" for a, b in summary.columns]
    summary = summary.reset_index().merge(
        met.groupby(["variant", "source"], sort=False)["n"].mean().round(1)
        .rename("n_test_mean").reset_index(), on=["variant", "source"])

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    met.to_csv(out / "fusion_metrics_per_seed.csv", index=False)
    summary.to_csv(out / "fusion_metrics.csv", index=False)
    covdf.to_csv(out / "coverage.csv", index=False)
    pd.DataFrame(preds).to_csv(out / "fused_predictions.csv", index=False)
    fail_df = pd.DataFrame(fails)
    fail_df.to_csv(out / "failure_table.csv", index=False)
    if len(fail_df):
        for (method, src), sub in fail_df.groupby(["method", "src"]):
            vdir = out / f"{method}_{SRC_SLUG[src]}"
            vdir.mkdir(parents=True, exist_ok=True)
            sub.to_csv(vdir / "failure_table.csv", index=False)
            render_failure(sub, method, src, vdir)

    render_metrics_png(summary, covdf, out / "fusion_table_detectors.png", args.mri_name)
    print("=" * 92)
    print("  T2 REMAINDER-PRODUCT detector decomposition (recreated #11) + MRI-only ref")
    print("=" * 92)
    print(summary[["variant", "source", "n_test_mean", "bacc_mean", "bacc_std",
                   "macro_f1_mean", "macro_auc_mean"]].to_string(index=False))
    print(f"\n  wrote -> {out}")
    return 0


# --------------------------------------------------------------------------- #
FOOT = [
    "REMAINDER-PRODUCT detector decomposition: the 3-class is built PURELY from two longitudinal binary",
    "detectors (no clinical T2 in the fused prediction; T2 only as reference). Detectors use the dedicated",
    "binary clinical encoders directly: P(CN) from clinT1 (prob_0), P(AD) from clinT1b (prob_1).",
    "MRI = BrainMVP (uniformer, full_ft, stochastic aug): T1 -> P(CN), T1b -> P(AD) @ m12.",
    "CN-det / AD-det = binary elastic net over the task's bl+m06+m12 probs (fit on VAL; +present flags,",
    "  missing -> 0.5);  CL+MRI additionally feeds the MRI prob into the EN.",
    "COMPOSE to 3-class [CN, MCI, AD], argmax:",
    "  structured : P = [ CN , (1-CN)(1-AD) , AD ]  then RENORMALISED to sum 1   (MCI = remainder).",
    "  EN         : multinomial elastic net over [ CN , AD , (1-CN)(1-AD) ], fit on VAL.",
    "REFs: clinical_only / temporal_CL_T2_EN use the clinical T2 3-class head; MRI-only = the two MRI",
    "  detectors composed alone, [t1, (1-t1)(1-t1b), t1b] renormalised (MRI-present, n~44).",
    "All metrics on TEST (fit on VAL). Target = concurrent m12 dx; cohort = TEST with an m12 visit (n=53).",
]


def render_metrics_png(summary, cov, out_path, mri_name=MRI_NAME_DEFAULT):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARN] matplotlib missing."); return
    idx = {(r["variant"], r["source"]): r for _, r in summary.iterrows()}
    body, numeric = [], []
    for variant, source, label in variant_rows(mri_name):
        if (variant, source) not in idx:
            continue
        r = idx[(variant, source)]
        cells, nums = [], []
        for key, _, with_n in h.METRIC_COLS:
            n = r["n_test_mean"] if with_n else None
            cells.append(h._fmt(r.get(f"{key}_mean"), r.get(f"{key}_std"), n))
            nums.append(r.get(f"{key}_mean", np.nan))
        body.append([label] + cells); numeric.append(nums)
    if not body:
        return
    numeric = np.array(numeric, float)
    n_cols = len(h.METRIC_COLS) + 1
    col_labels = ["Variant  (CN<-clinT1, AD<-clinT1b, MCI=remainder)"] + [c[1] for c in h.METRIC_COLS]
    col_chars = [max(len(str(x)) for x in [col_labels[j]] + [b[j] for b in body]) + 2 for j in range(n_cols)]
    total = sum(col_chars)
    fig, ax = plt.subplots(figsize=(total * 0.102, 0.6 + 0.34 * (len(body) + 2)))
    ax.axis("off")
    tab = ax.table(cellText=body, colLabels=col_labels, loc="upper center", cellLoc="center",
                   colLoc="center", colWidths=[c / total for c in col_chars])
    tab.auto_set_font_size(False); tab.set_fontsize(9); tab.scale(1.0, 1.4)
    for j in range(len(h.METRIC_COLS)):
        col = numeric[:, j]
        if not np.all(np.isnan(col)):
            tab[int(np.nanargmax(col)) + 1, j + 1].set_text_props(weight="bold")
    for j in range(n_cols):
        tab[0, j].set_facecolor("#ECECEC"); tab[0, j].set_text_props(weight="bold")
    for i in range(len(body) + 1):
        tab[i, 0].set_text_props(ha="left")
    cov_line = ("  coverage (test, per seed): "
                + ", ".join(f"s{int(c.seed)} clinT1/T1b={int(c.n_test_cn_temporal)} "
                            f"MRI={int(c.n_test_mri_cn)}/{int(c.n_test)}" for c in cov.itertuples()))
    plt.title("T2 remainder-product detector decomposition (recreated #11)\n"
              "CN<-clinT1, AD<-clinT1b @ bl+m06+m12 (+ BrainMVP MRI); MCI = (1-CN)(1-AD).  "
              "(Test mean +/- std across seeds, n)", pad=6, fontsize=10)
    ax.text(0.0, -0.04, "\n".join(FOOT + ["", cov_line]), transform=ax.transAxes, ha="left",
            va="top", fontsize=7.3, family="monospace", linespacing=1.4)
    fig.savefig(out_path, dpi=180, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print(f"  PNG: {out_path}")


GREEN, RED, GREY, HEAD = "#cdebc5", "#f4c7c3", "#e6e6e6", "#d9d9d9"
FAIL_COLS = [("clin T2@bl", "clinT2_bl"), ("clin T2@m06", "clinT2_m06"),
             ("clin T2@m12", "clinT2_m12")]


def _color(c):
    return GREY if (c is None or (isinstance(c, float) and np.isnan(c))) else (GREEN if bool(c) else RED)


def render_failure(sub, method, src, vdir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    mtag = " + BrainMVP MRI" if src == "CL+MRI" else ""
    caption = [
        f"detector-{method}: CN-det = EN(clinT1 @bl,m06,m12{mtag}); AD-det = EN(clinT1b @bl,m06,m12{mtag}).",
        ("fused = [ CN-det , (1-CN-det)(1-AD-det) , AD-det ] renormalised (MCI = remainder)."
         if method == "structured" else
         "fused = multinomial EN( CN-det , AD-det , (1-CN-det)(1-AD-det) ), fit on VAL."),
        "Columns clin T2@bl / @m06 / @m12 = clinical T2 3-class argmax at each visit (diagnosis trajectory);",
        "these are CONTEXT only (not used by this fusion). green=correct, red=wrong, grey=visit absent.",
        "agreement compares clin T2 @m12 (reference) vs fused.",
    ]
    seeds = sorted(sub["seed"].unique())
    pats = ["both_correct", "base_better", "fused_better", "both_wrong"]
    srows = [[f"seed {s}"] + [int(sub[sub.seed == s]["agreement"].value_counts().get(p, 0)) for p in pats]
             for s in seeds]
    srows.append(["OVERALL"] + [int(sub["agreement"].value_counts().get(p, 0)) for p in pats])
    fig, ax = plt.subplots(figsize=(2 + 1.9 * (len(pats) + 1), 0.9 + 0.32 * (len(srows) + 4)))
    ax.axis("off")
    t = ax.table(cellText=srows, colLabels=["cohort"] + pats, loc="upper center", cellLoc="center", colLoc="center")
    t.auto_set_font_size(False); t.set_fontsize(9); t.scale(1.0, 1.4)
    for j in range(len(pats) + 1):
        t[0, j].set_facecolor(HEAD); t[0, j].set_text_props(weight="bold")
    plt.title(f"Failure-mode agreement (clin T2 @m12 vs fused)\ndetector-{method} | {src}", pad=6, fontsize=11)
    ax.text(0.0, -0.05, "\n".join(caption), transform=ax.transAxes, ha="left", va="top",
            fontsize=8, family="monospace", linespacing=1.4)
    fig.savefig(vdir / "failure_summary.png", dpi=170, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    cols = ["Patient_ID", "y_true"] + [c[0] for c in FAIL_COLS] + ["fused", "agreement"]
    for s in seeds:
        d = sub[sub.seed == s].reset_index(drop=True)
        body, colors = [], []
        for _, r in d.iterrows():
            row = [str(r["Patient_ID"]), str(r["y_true"])]
            crow = ["white", "white"]
            for _, base in FAIL_COLS:
                row.append(str(r[f"{base}_pred"])); crow.append(_color(r[f"{base}_correct"]))
            row.append(str(r["fused_pred"])); crow.append(GREEN if r["fused_correct"] else RED)
            row.append(str(r["agreement"])); crow.append("white")
            body.append(row); colors.append(crow)
        fig, ax = plt.subplots(figsize=(1.6 + 1.6 * len(cols), 0.8 + 0.27 * (len(body) + 4)))
        ax.axis("off")
        tab = ax.table(cellText=body, colLabels=cols, cellColours=colors, loc="upper center",
                       cellLoc="center", colLoc="center")
        tab.auto_set_font_size(False); tab.set_fontsize(8); tab.scale(1.0, 1.25)
        for j in range(len(cols)):
            tab[0, j].set_facecolor(HEAD); tab[0, j].set_text_props(weight="bold")
        plt.title(f"Failure table — detector-{method} | {src}  (seed {s})", pad=6, fontsize=11)
        ax.text(0.0, -0.03, "\n".join(caption), transform=ax.transAxes, ha="left", va="top",
                fontsize=7.5, family="monospace", linespacing=1.4)
        fig.savefig(vdir / f"failure_seed{s}.png", dpi=170, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)
    print(f"  PNG: {vdir} (summary + {len(seeds)} seed tables)")


if __name__ == "__main__":
    raise SystemExit(main())
