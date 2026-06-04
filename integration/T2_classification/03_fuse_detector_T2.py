"""
03_fuse_detector_T2.py
======================
T2 (CN/MCI/AD) via CLASS-WISE learned late fusion across clinical timepoints + MRI.

For each class c a small learned combiner g_c (LR or EN, class-c-vs-rest, fit on VAL) takes that
class's probability from several sources and outputs a score; the three scores are RENORMALISED to
sum to 1; argmax. Two clinical SOURCE formulations are compared:

  Formulation A  (clinical = T2 3-class marginal):
     P(CN)_fused  = g_CN([ P(CN)_T2·CL_bl , P(CN)_T2·CL_m06 , P(CN)_T2·CL_m12 (, P(CN)_MRI·T1_m12) ])
     P(AD)_fused  = g_AD([ P(AD)_T2·CL_bl , P(AD)_T2·CL_m06 , P(AD)_T2·CL_m12 (, P(AD)_MRI·T1b_m12) ])
  Formulation B  (clinical = dedicated binary detector):
     P(CN)_fused  = g_CN([ P(CN)_T1·CL_bl , P(CN)_T1·CL_m06 , P(CN)_T1·CL_m12 (, P(CN)_MRI·T1_m12) ])
     P(AD)_fused  = g_AD([ P(AD)_T1b·CL_bl, P(AD)_T1b·CL_m06, P(AD)_T1b·CL_m12(, P(AD)_MRI·T1b_m12)])
  MCI (both, no MRI; there is no MCI detector):
     P(MCI)_fused = g_MCI([ P(MCI)_T2·CL_bl , P(MCI)_T2·CL_m06 , P(MCI)_T2·CL_m12 ])
  then  P_fused = [CN, MCI, AD] / sum,  argmax.   "(, MRI)" terms present only in the +MRI variants.

Clinical probs come from the dedicated encoders DIRECTLY: T2 3-class marginals (prob_0/1/2), T1 (prob_0
= P(CN)), T1b (prob_1 = P(AD)); MRI = BrainMVP T1 (prob_class_0=P(CN)) / T1b (prob_class_1=P(AD)) @m12.
Two REFERENCE rows use the T2 3-class head only: clinical_only = argmax CL_m12; temporal_CL_T2_EN =
multinomial EN over CL T2 [bl,m06,m12].

Grid = 2 formulations x {CL-only, CL+MRI} x {LR, EN} = 8 fusion rows + 2 refs. Target = concurrent m12
dx; cohort = TEST with an m12 visit (n=53; MRI imputed 0.5 + present-flag where absent, so cohort kept).
Meta-learners fit per-seed on VAL, reported on TEST, mean +/- std over seeds. Reuses the MRI arm's
harness (import only; no edit).

Output: integration/T2_classification/outputs/up_to_m12/detectors/class_wise_detector_fusion/
  fusion_table_detectors.png, fusion_metrics.csv, fusion_metrics_per_seed.csv, coverage.csv,
  fused_predictions.csv, failure_table.csv, failure_summary.csv,
  <formulation>_<mri>_<combiner>/ {failure_table.csv, failure_seed{n}.png, failure_summary.png}
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
TPS = ["bl", "m06", "m12"]
OUT_DIR_DEFAULT = (Path(__file__).resolve().parent / "outputs" / "up_to_m12"
                   / "detectors" / "class_wise_detector_fusion")

FORMS = ["A", "B"]                      # A=clinical T2-marginal, B=clinical binary detector
MRIS = ["CL", "CL+MRI"]
KINDS = ["LR", "EN"]
MCI_COLS = [f"t2_{tp}_1" for tp in TPS]  # clinical T2 P(MCI) @ bl,m06,m12 (both formulations, no MRI)


def cn_cols(form, mri):
    base = [f"t2_{tp}_0" for tp in TPS] if form == "A" else [f"cn_{tp}" for tp in TPS]
    return base + (["mri_cn"] if mri == "CL+MRI" else [])


def ad_cols(form, mri):
    base = [f"t2_{tp}_2" for tp in TPS] if form == "A" else [f"ad_{tp}" for tp in TPS]
    return base + (["mri_ad"] if mri == "CL+MRI" else [])


# --------------------------------------------------------------------------- #
# Loaders
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


def build_frame(seed):
    f = load_t2_block(seed, "m12")[["Patient_ID", "split", "y_m12",
                                    "t2_m12_0", "t2_m12_1", "t2_m12_2"]].copy()
    f["y_true"] = f["y_m12"].astype(int)
    for tp in ("bl", "m06"):
        f = f.merge(load_t2_block(seed, tp)[["Patient_ID", f"t2_{tp}_0", f"t2_{tp}_1", f"t2_{tp}_2"]],
                    on="Patient_ID", how="left")
    for tp in TPS:
        f = f.merge(load_bin_prob(seed, tp, T1, "prob_0", f"cn_{tp}"), on="Patient_ID", how="left")
        f = f.merge(load_bin_prob(seed, tp, T1B, "prob_1", f"ad_{tp}"), on="Patient_ID", how="left")
    f = f.merge(h.load_mri_detector(seed, MRI_T1_TMPL, "m12", "prob_class_0", "mri_cn"),
                on="Patient_ID", how="left")
    f = f.merge(h.load_mri_detector(seed, MRI_T1B_TMPL, "m12", "prob_class_1", "mri_ad"),
                on="Patient_ID", how="left")
    return f


# --------------------------------------------------------------------------- #
# Per-class learned combiner g_c (LR or EN), class-c-vs-rest, fit on VAL -> P(c)
# --------------------------------------------------------------------------- #
def _bin_X(df, cols, impute=0.5):
    X = df[cols].copy()
    pres = X.notna().astype(float)
    return np.column_stack([X.fillna(impute).to_numpy(float), pres.to_numpy(float)])


def fit_binary(kind, Xv, yv, Xall, seed):
    if len(np.unique(yv)) < 2:
        return None
    if kind == "LR":
        clf = LogisticRegression(max_iter=5000, random_state=seed).fit(Xv, yv)   # plain L2 logistic
    else:  # EN
        counts = np.bincount(yv, minlength=2)
        n_splits = min(3, int(counts[counts > 0].min()))
        if n_splits < 2:
            clf = LogisticRegression(penalty="elasticnet", solver="saga", C=1.0, l1_ratio=0.5,
                                     max_iter=5000, random_state=seed).fit(Xv, yv)
        else:
            clf = GridSearchCV(
                LogisticRegression(penalty="elasticnet", solver="saga", max_iter=5000,
                                   random_state=seed),
                {"C": [0.1, 1.0, 10.0], "l1_ratio": [0.1, 0.5, 0.9]}, scoring="balanced_accuracy",
                cv=StratifiedKFold(n_splits, shuffle=True, random_state=seed), n_jobs=-1).fit(Xv, yv)
    classes = clf.best_estimator_.classes_ if hasattr(clf, "best_estimator_") else clf.classes_
    return clf.predict_proba(Xall)[:, list(classes).index(1)]


def class_score(val, test, cols, target_class, kind, seed):
    """g_c over `cols` predicting (y==target_class); returns (score_val, score_test) or (None,None)."""
    yv = (val["y_true"].to_numpy(int) == target_class).astype(int)
    Xv, Xt = _bin_X(val, cols), _bin_X(test, cols)
    p = fit_binary(kind, Xv, yv, np.vstack([Xv, Xt]), seed)
    return (None, None) if p is None else (p[:len(Xv)], p[len(Xv):])


def present_mean(df, cols):
    """Present-only: row-wise MEAN of the NON-missing probabilities in `cols` (no imputation)."""
    return df[cols].mean(axis=1, skipna=True).to_numpy(float)


def _t2_impute(df, prior):
    blocks = []
    for tp in TPS:
        X = df[[f"t2_{tp}_0", f"t2_{tp}_1", f"t2_{tp}_2"]].to_numpy(float)
        miss = np.isnan(X).any(axis=1)
        X = X.copy(); X[miss] = prior
        blocks.append(X)
    ind = np.column_stack([df[f"t2_{tp}_0"].notna().to_numpy(float) for tp in TPS])
    return np.column_stack(blocks + [ind])


def variant_rows():
    """(variant, source, label) for the metrics table, in display order."""
    out = [("clinical_only", "CL_m12", "REF  clinical T2 3-class @ m12 (argmax)"),
           ("temporal_CL_T2_EN", "CL bl+m06+m12", "REF  clinical T2 3-class EN @ bl+m06+m12"),
           ("MRI_only", "BrainMVP T1+T1b stoch @m12",
            "REF  MRI-only [BrainMVP T1+T1b stoch @m12]  (struct: [t1, (1-t1)(1-t1b), t1b])")]
    kd = {"LR": "LR", "EN": "EN", "meanP": "mean(present)"}
    for form in FORMS:
        for mri in MRIS:
            for kind in KINDS + ["meanP"]:
                clin = "T2-marg" if form == "A" else "T1/T1b-bin"
                out.append((f"{form}_{kind}", mri, f"{form}: clin={clin} | {mri} | {kd[kind]}"))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS_DEFAULT))
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR_DEFAULT)
    args = ap.parse_args()

    rows, preds, cov, fails = [], [], [], []
    for seed in args.seeds:
        f = build_frame(seed)
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

        def add_fail(form, mri, kind, fused_pred, proba):
            tr = test.reset_index(drop=True)
            for i in range(len(tr)):
                y = int(tr.loc[i, "y_true"]); fu = int(fused_pred[i])
                rec = {"form": form, "mri": mri, "kind": kind, "seed": seed,
                       "Patient_ID": tr.loc[i, "Patient_ID"], "y_true": CLASS_NAMES[y],
                       "fused_pred": CLASS_NAMES[fu], "fused_correct": fu == y,
                       "pCN": round(float(proba[i, 0]), 3), "pMCI": round(float(proba[i, 1]), 3),
                       "pAD": round(float(proba[i, 2]), 3)}
                # clinical T2 3-class argmax at EACH timepoint (the clinical diagnosis trajectory)
                for tp in TPS:
                    v = tr.loc[i, [f"t2_{tp}_0", f"t2_{tp}_1", f"t2_{tp}_2"]].to_numpy(float)
                    if np.isnan(v).any():
                        rec[f"clinT2_{tp}_pred"], rec[f"clinT2_{tp}_correct"] = "absent", np.nan
                    else:
                        p = int(v.argmax())
                        rec[f"clinT2_{tp}_pred"], rec[f"clinT2_{tp}_correct"] = CLASS_NAMES[p], bool(p == y)
                base_ok = rec["clinT2_m12_correct"] is True   # base reference = clin T2 @ m12
                rec["agreement"] = ("both_correct" if base_ok and fu == y else
                                    "base_better" if base_ok else
                                    "fused_better" if fu == y else "both_wrong")
                fails.append(rec)

        # ---- references ----
        add("clinical_only", "CL_m12", cpt.argmax(1), cpt)
        base = h.fit_en(_t2_impute(val, prior), yv, _t2_impute(test, prior), seed)
        if base is not None:
            add("temporal_CL_T2_EN", "CL bl+m06+m12", base[0], base[1])
        # MRI-only reference: 3-class from the two MRI detectors alone (BrainMVP, stochastic aug,
        # T1 -> P(CN) / T1b -> P(AD) @ m12), structured composition; on MRI-present patients (n~44).
        tm = test[test["mri_cn"].notna() & test["mri_ad"].notna()]
        if len(tm):
            mcn, mad = tm["mri_cn"].to_numpy(float), tm["mri_ad"].to_numpy(float)
            Pm = np.clip(np.stack([mcn, (1 - mcn) * (1 - mad), mad], axis=1), 1e-9, None)
            Pm = Pm / Pm.sum(axis=1, keepdims=True)
            add("MRI_only", "BrainMVP T1+T1b stoch @m12", Pm.argmax(1), Pm, df=tm)

        # ---- class-wise fusion variants (learned LR/EN, + parameter-free present-only mean) ----
        for form in FORMS:
            for mri in MRIS:
                for kind in KINDS + ["meanP"]:
                    if kind == "meanP":   # present-only: mean of the NON-missing probs, no imputation
                        cn_t = present_mean(test, cn_cols(form, mri))
                        ad_t = present_mean(test, ad_cols(form, mri))
                        mci_t = present_mean(test, MCI_COLS)
                    else:                 # learned per-class combiner (g_c fit on VAL)
                        _, cn_t = class_score(val, test, cn_cols(form, mri), 0, kind, seed)
                        _, ad_t = class_score(val, test, ad_cols(form, mri), 2, kind, seed)
                        _, mci_t = class_score(val, test, MCI_COLS, 1, kind, seed)
                    if any(x is None for x in (cn_t, ad_t, mci_t)):
                        continue
                    P = np.clip(np.stack([cn_t, mci_t, ad_t], axis=1), 1e-9, None)
                    P = P / P.sum(axis=1, keepdims=True)
                    add(f"{form}_{kind}", mri, P.argmax(1), P)
                    add_fail(form, mri, kind, P.argmax(1), P)

        cov.append({"seed": seed, "n_test": int(len(test)),
                    "n_test_cn_temporal": int(test[[f"cn_{tp}" for tp in TPS]].notna().any(axis=1).sum()),
                    "n_test_mri_cn": int(test["mri_cn"].notna().sum()),
                    "n_test_mri_ad": int(test["mri_ad"].notna().sum())})

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
        (fail_df.groupby(["form", "mri", "kind", "agreement"]).size().unstack(fill_value=0)
         .reindex(columns=["both_correct", "base_better", "fused_better", "both_wrong"],
                  fill_value=0).reset_index()).to_csv(out / "failure_summary.csv", index=False)
        for (form, mri, kind), sub in fail_df.groupby(["form", "mri", "kind"]):
            slug = f"{form}_{mri.replace('+', 'plus')}_{kind}"
            vdir = out / slug
            vdir.mkdir(parents=True, exist_ok=True)
            sub.to_csv(vdir / "failure_table.csv", index=False)
            render_failure(sub, form, mri, kind, vdir)

    render_metrics_png(summary, covdf, out / "fusion_table_detectors.png")
    print("=" * 94)
    print("  T2 class-wise learned late fusion (per-class LR/EN over clinical timepoints + MRI)")
    print("=" * 94)
    print("\n[coverage] test detector-input presence per seed:")
    print(covdf.to_string(index=False))
    print("\n[summary] mean over seeds (fit VAL, report TEST):")
    print(summary[["variant", "source", "n_test_mean", "bacc_mean", "bacc_std",
                   "macro_f1_mean", "macro_auc_mean"]].to_string(index=False))
    print(f"\n  wrote -> {out}")
    return 0


# --------------------------------------------------------------------------- #
# Metrics table PNG
# --------------------------------------------------------------------------- #
FOOT = [
    "CLASS-WISE learned late fusion: per class c, a combiner g_c (class-c-vs-rest, fit on VAL) maps that",
    "class's probabilities to a score; then P_fused = [P(CN),P(MCI),P(AD)] / (P(CN)+P(MCI)+P(AD)), argmax.",
    "g = LR (plain L2 logistic) or EN (elastic-net L1/L2, CV grid); both fit on VAL. 'mean(present)' is",
    "PRESENT-ONLY: g_c = the MEAN of only the non-missing terms (no imputation, no learned weights, n=53).",
    "All probabilities taken DIRECTLY from the dedicated encoders. Full per-class forms (written out in full):",
    "",
    "Formulation A  (clinical = T2 3-class marginal P(class)):",
    "  CL:     P(CN)  = g_CN([ P(CN) T2.CL @bl , P(CN) T2.CL @m06 , P(CN) T2.CL @m12 ])",
    "          P(AD)  = g_AD([ P(AD) T2.CL @bl , P(AD) T2.CL @m06 , P(AD) T2.CL @m12 ])",
    "  CL+MRI: P(CN)  = g_CN([ P(CN) T2.CL @bl , P(CN) T2.CL @m06 , P(CN) T2.CL @m12 , P(CN) MRI.T1 @m12 ])",
    "          P(AD)  = g_AD([ P(AD) T2.CL @bl , P(AD) T2.CL @m06 , P(AD) T2.CL @m12 , P(AD) MRI.T1b @m12 ])",
    "Formulation B  (clinical = dedicated binary detector: T1 -> P(CN), T1b -> P(AD)):",
    "  CL:     P(CN)  = g_CN([ P(CN) T1.CL @bl , P(CN) T1.CL @m06 , P(CN) T1.CL @m12 ])",
    "          P(AD)  = g_AD([ P(AD) T1b.CL @bl , P(AD) T1b.CL @m06 , P(AD) T1b.CL @m12 ])",
    "  CL+MRI: P(CN)  = g_CN([ P(CN) T1.CL @bl , P(CN) T1.CL @m06 , P(CN) T1.CL @m12 , P(CN) MRI.T1 @m12 ])",
    "          P(AD)  = g_AD([ P(AD) T1b.CL @bl , P(AD) T1b.CL @m06 , P(AD) T1b.CL @m12 , P(AD) MRI.T1b @m12 ])",
    "MCI (both formulations, every variant; no MRI):",
    "          P(MCI) = g_MCI([ P(MCI) T2.CL @bl , P(MCI) T2.CL @m06 , P(MCI) T2.CL @m12 ])",
    "MRI = BrainMVP (uniformer, full_ft, stochastic aug): T1 -> P(CN), T1b -> P(AD) @ m12 (the only MRI model used).",
    "MRI prob absent for a patient -> imputed 0.5 + a 'present' flag feature (cohort kept at n=53).",
    "REFs: clinical_only / temporal_CL_T2_EN use the clinical T2 3-class head; MRI-only = 3-class from the two",
    "  MRI detectors alone, [t1, (1-t1)(1-t1b), t1b] renormalised, on MRI-present patients (n~44).",
    "All metrics on TEST (meta-learners fit on VAL). Target = concurrent m12 dx; cohort = TEST with an m12 visit.",
]


def render_metrics_png(summary, cov, out_path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARN] matplotlib missing -- skip metrics PNG."); return
    idx = {(r["variant"], r["source"]): r for _, r in summary.iterrows()}
    body, numeric = [], []
    for variant, source, label in variant_rows():
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
    n_rows, n_cols = len(body), len(h.METRIC_COLS) + 1
    col_labels = ["Variant  (formulation: clin source | CL or CL+MRI | combiner)"] + \
                 [c[1] for c in h.METRIC_COLS]
    col_chars = [max(len(str(x)) for x in [col_labels[j]] + [b[j] for b in body]) + 2
                 for j in range(n_cols)]
    total = sum(col_chars)
    fig, ax = plt.subplots(figsize=(total * 0.102, 0.6 + 0.34 * (n_rows + 2)))
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
    for i in range(n_rows + 1):
        tab[i, 0].set_text_props(ha="left")
    cov_line = ("  coverage (test, per seed): "
                + ", ".join(f"s{int(c.seed)} clinT1/T1b={int(c.n_test_cn_temporal)} "
                            f"MRI={int(c.n_test_mri_cn)}/{int(c.n_test)}" for c in cov.itertuples()))
    plt.title("T2 class-wise learned late fusion — per-class LR/EN over clinical timepoints "
              "(+ MRI), renormalised\n(A: clinical=T2 marginal | B: clinical=T1/T1b binary detector).  "
              "(mean +/- std across seeds, n)", pad=6, fontsize=10)
    ax.text(0.0, -0.04, "\n".join(FOOT + ["", cov_line]), transform=ax.transAxes, ha="left",
            va="top", fontsize=7.2, family="monospace", linespacing=1.4)
    fig.savefig(out_path, dpi=180, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print(f"  PNG: {out_path}")


# --------------------------------------------------------------------------- #
# Failure-table PNGs (explicit per-class formula in the caption)
# --------------------------------------------------------------------------- #
GREEN, RED, GREY, HEAD = "#cdebc5", "#f4c7c3", "#e6e6e6", "#d9d9d9"


def _color(c):
    return GREY if (c is None or (isinstance(c, float) and np.isnan(c))) else (GREEN if bool(c) else RED)


def _caption(form, mri, kind):
    if form == "A":
        cn = ["P(CN) T2.CL @bl", "P(CN) T2.CL @m06", "P(CN) T2.CL @m12"]
        ad = ["P(AD) T2.CL @bl", "P(AD) T2.CL @m06", "P(AD) T2.CL @m12"]
        clin = "clinical = T2 3-class marginal"
    else:
        cn = ["P(CN) T1.CL @bl", "P(CN) T1.CL @m06", "P(CN) T1.CL @m12"]
        ad = ["P(AD) T1b.CL @bl", "P(AD) T1b.CL @m06", "P(AD) T1b.CL @m12"]
        clin = "clinical = T1/T1b binary detector"
    if mri == "CL+MRI":
        cn = cn + ["P(CN) MRI.T1 @m12"]
        ad = ad + ["P(AD) MRI.T1b @m12"]
    gdesc = ("mean of the PRESENT terms (present-only; no imputation, no learned weights)"
             if kind == "meanP" else f"{kind} (class-vs-rest, fit on VAL)")
    return [
        f"Formulation {form} ({clin}); combiner g = {gdesc}.",
        f"  P(CN)  = g_CN([ {' , '.join(cn)} ])",
        f"  P(AD)  = g_AD([ {' , '.join(ad)} ])",
        "  P(MCI) = g_MCI([ P(MCI) T2.CL @bl , P(MCI) T2.CL @m06 , P(MCI) T2.CL @m12 ])",
        "  P_fused = [P(CN), P(MCI), P(AD)] / (P(CN)+P(MCI)+P(AD)) ; argmax.",
        "Columns clin T2@bl / @m06 / @m12 = clinical T2 3-class argmax at each visit (diagnosis trajectory);",
        "P(CN)/P(MCI)/P(AD) = the fused class probabilities. green=correct, red=wrong, grey=visit absent.",
        "agreement compares clin T2 @m12 (base reference) vs fused.",
    ]


def render_failure(sub, form, mri, kind, vdir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    caption = _caption(form, mri, kind)
    tag = f"{form} | clin={'T2-marg' if form=='A' else 'T1/T1b-bin'} | {mri} | {kind}"
    seeds = sorted(sub["seed"].unique())
    pats = ["both_correct", "base_better", "fused_better", "both_wrong"]
    srows = []
    for s in seeds:
        vc = sub[sub.seed == s]["agreement"].value_counts()
        srows.append([f"seed {s}"] + [int(vc.get(p, 0)) for p in pats])
    vc = sub["agreement"].value_counts()
    srows.append(["OVERALL"] + [int(vc.get(p, 0)) for p in pats])
    fig, ax = plt.subplots(figsize=(2 + 1.9 * (len(pats) + 1), 0.9 + 0.32 * (len(srows) + 4)))
    ax.axis("off")
    t = ax.table(cellText=srows, colLabels=["cohort"] + pats, loc="upper center",
                 cellLoc="center", colLoc="center")
    t.auto_set_font_size(False); t.set_fontsize(9); t.scale(1.0, 1.4)
    for j in range(len(pats) + 1):
        t[0, j].set_facecolor(HEAD); t[0, j].set_text_props(weight="bold")
    plt.title(f"Failure-mode agreement (base clin T2 @m12 vs fused)\n{tag}", pad=6, fontsize=11)
    ax.text(0.0, -0.05, "\n".join(caption), transform=ax.transAxes, ha="left", va="top",
            fontsize=8, family="monospace", linespacing=1.4)
    fig.savefig(vdir / "failure_summary.png", dpi=170, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    cols = ["Patient_ID", "y_true", "clin T2@bl", "clin T2@m06", "clin T2@m12", "fused",
            "P(CN)", "P(MCI)", "P(AD)", "agreement"]
    for s in seeds:
        d = sub[sub.seed == s].reset_index(drop=True)
        body, colors = [], []
        for _, r in d.iterrows():
            body.append([str(r["Patient_ID"]), str(r["y_true"]),
                         str(r["clinT2_bl_pred"]), str(r["clinT2_m06_pred"]), str(r["clinT2_m12_pred"]),
                         str(r["fused_pred"]), f"{r['pCN']:.2f}", f"{r['pMCI']:.2f}",
                         f"{r['pAD']:.2f}", str(r["agreement"])])
            colors.append(["white", "white",
                           _color(r["clinT2_bl_correct"]), _color(r["clinT2_m06_correct"]),
                           _color(r["clinT2_m12_correct"]),
                           GREEN if r["fused_correct"] else RED,
                           "white", "white", "white", "white"])
        fig, ax = plt.subplots(figsize=(1.6 + 1.5 * len(cols), 0.8 + 0.27 * (len(body) + 4)))
        ax.axis("off")
        tab = ax.table(cellText=body, colLabels=cols, cellColours=colors, loc="upper center",
                       cellLoc="center", colLoc="center")
        tab.auto_set_font_size(False); tab.set_fontsize(8); tab.scale(1.0, 1.25)
        for j in range(len(cols)):
            tab[0, j].set_facecolor(HEAD); tab[0, j].set_text_props(weight="bold")
        plt.title(f"Failure table — {tag}  (seed {s})", pad=6, fontsize=11)
        ax.text(0.0, -0.03, "\n".join(caption), transform=ax.transAxes, ha="left", va="top",
                fontsize=7.5, family="monospace", linespacing=1.4)
        fig.savefig(vdir / f"failure_seed{s}.png", dpi=170, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)
    print(f"  PNG: {vdir} (summary + {len(seeds)} seed tables)")


if __name__ == "__main__":
    raise SystemExit(main())
