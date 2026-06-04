"""
01_fuse_T4_cumulative.py — T4 horizon via cumulative binary-T3 probabilities
============================================================================
The direct 3-class T4 models are too degenerate to fuse (see the cross-model
table: bACC≈chance, severe overfit). Per `integration/T4/notes.txt`, we instead
DERIVE the 3 horizon buckets (<3y / 3-7y / ≥7y) from the cumulative binary
conversion probabilities, separately for clinical (baseline) and MRI (m12),
then fuse the two modalities.

For modality m, using p_m(<3y) [3yr model] and p_m(<7y) [7yr model]:
    bucket_m = [ p_m(<3y),  p_m(<7y) − p_m(<3y),  1 − p_m(<7y) ]
Monotonicity guard (independently-trained binaries can give p7<p3):
    p7 ← max(p3, p7);  clip buckets to ≥0;  renormalise to sum 1.

Sources (leakage-safe: the baseline-split TEST patients are out-of-fold for
every T3 model, since clinical T3a/T3d and the MRI cached heads all train on
the same per-seed baseline splits):
  CL  p(<3y) = BioClinical-ModernBERT-large / T3a_conv3y  (baseline parquet, prob_1)
  CL  p(<7y) = ModernBERT-large            / T3d_conv7y  (baseline parquet, prob_1)
  MRI p(<3y) = ViT-MAE cached / T3a_conv3y  (m12, prob_class_1)
  MRI p(<7y) = ViT-MAE cached / T3c_conv7y  (m12, prob_class_1)
Eval cohort = T4 converters (Label_T4 not-NaN); val/test = the baseline split
(out-of-fold); y_true = Label_T4 (0=<3y,1=3-7y,2=≥7y).

Meta-learners reuse the engine `01_fuse_T2.py` (fit on VAL, report TEST):
clinical_only, mri_only, equal (0.5), weighted_avg (val bACC), weighted_avg_J
(Youden), elastic_net, lr. Small n (T4 test ≈ 12-15) → indicative only.

Run:  conda run -n snp python integration/T4/01_fuse_T4_cumulative.py
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# --- import the fusion engine as h (filename starts with a digit) ----------
HERE = Path(__file__).resolve().parent
ENGINE = HERE.parent / "T2_classification" / "01_fuse_T2.py"
_spec = importlib.util.spec_from_file_location("_fuse_engine", ENGINE)
h = importlib.util.module_from_spec(_spec)
sys.modules["_fuse_engine"] = h
_spec.loader.exec_module(h)
h._set_classes(["h_lt3", "h_3_7", "h_ge7"])     # 3-class horizon buckets

DERIV = Path("D:/ADNI_BIDS_project/derivatives")
CLIN_ROOT = DERIV / "encoder_outputs_no_cdr_post_exclusion"
MRI_ROOT = DERIV / "mri_embeddings"
MASTER = (DERIV / "mri_clinical_matched" / "VISCODE_2_aligned_extended_post_exclusion"
          / "master_mri_clinical_matched_viscode2_extended_post_exclusion.csv")
SEEDS = [0, 1, 2]

# (modality source) -> how to read the conversion probability p(convert ≤ horizon)
CL3 = ("BioClinical-ModernBERT-large", "T3a_conv3y")   # p_CL(<3y)
CL7 = ("ModernBERT-large", "T3d_conv7y")               # p_CL(<7y)  (clinical names 7y as T3d)
MRI3_TASK, MRI7_TASK = "T3a_conv3y", "T3c_conv7y"       # MRI names 7y as T3c
MRI_TP = "m12"
BMVP_FT_ROOT = DERIV / "brainmvp_embeddings"

# MRI source key -> (display label, template builder (task, seed) -> Path)
MRI_SOURCES = {
    "vit_mae_cached": ("ViT-MAE frozen/none",
                       lambda t, s: MRI_ROOT / t / "vit_mae_frozen_none_cached" / f"embeddings_seed_{s}.csv"),
    "brainmvp_ft":    ("BrainMVP full_ft plus_original",
                       lambda t, s: BMVP_FT_ROOT / t / "aug_plus_original" / f"seed_{s}" / f"embeddings_seed_{s}.csv"),
}


def _clin_prob(model, task, seed):
    """Return Patient_ID, split, p(=prob_1 conversion prob) from a clinical parquet."""
    p = CLIN_ROOT / model / task / f"seed_{seed}" / "full_ft" / "embeddings" / "embeddings.parquet"
    d = pd.read_parquet(p, columns=["Patient_ID", "split", "in_task_cohort", "prob_1"])
    d = d[d["in_task_cohort"].fillna(False)]
    return d[["Patient_ID", "split", "prob_1"]]


def _mri_prob(task, seed, tmpl):
    """Return Patient_ID, split_mri, p(=prob_class_1) from the m12 MRI CSV (one row/patient)."""
    d = pd.read_csv(tmpl(task, seed), low_memory=False)
    d = d[d["adni_viscode"] == MRI_TP].drop_duplicates("Patient_ID", keep="first")
    return d[["Patient_ID", "split", "prob_class_1"]].rename(
        columns={"split": "split_mri", "prob_class_1": "p"})


def _buckets(p3, p7):
    """Cumulative -> [<3y, 3-7y, ≥7y] with monotonicity guard + renormalise."""
    p3 = np.asarray(p3, float); p7 = np.maximum(np.asarray(p7, float), p3)
    b = np.clip(np.stack([p3, p7 - p3, 1.0 - p7], axis=1), 0.0, None)
    return b / b.sum(1, keepdims=True)


def _label_map():
    m = pd.read_csv(MASTER, low_memory=False)[["Patient_ID", "Label_T4"]].dropna(subset=["Label_T4"])
    return m.drop_duplicates("Patient_ID").set_index("Patient_ID")["Label_T4"].astype(int)


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mri-source", default="brainmvp_ft", choices=list(MRI_SOURCES),
                    help="MRI prob source for the cumulative buckets.")
    args = ap.parse_args()
    mri_label, mri_tmpl = MRI_SOURCES[args.mri_source]
    (HERE / "mri_source.txt").write_text(mri_label)     # sidecar for 02_summary title

    y_by_pid = _label_map()
    print(f"[T4-cumulative] T4 cohort (Label_T4 present): {len(y_by_pid)} subjects  "
          f"| MRI source: {mri_label}")
    rows, n_mismatch = [], 0
    pred_rows = []                 # per test patient/seed (clinical/mri/lr — for failure plots)
    all_preds = []                 # per test patient/seed/variant (for POOLED bACC)
    Xv_all, yv_all = [], []        # pooled VAL for the illustrative LR-weights fit
    FEATS = ["CL_<3y", "CL_3-7y", "CL_>=7y", "MRI_<3y", "MRI_3-7y", "MRI_>=7y"]

    for seed in SEEDS:
        c3 = _clin_prob(*CL3, seed).rename(columns={"prob_1": "p_cl3"})
        c7 = _clin_prob(*CL7, seed).rename(columns={"prob_1": "p_cl7"})[["Patient_ID", "p_cl7"]]
        m3 = _mri_prob(MRI3_TASK, seed, mri_tmpl).rename(columns={"p": "p_mri3"})
        m7 = _mri_prob(MRI7_TASK, seed, mri_tmpl)[["Patient_ID", "p"]].rename(columns={"p": "p_mri7"})

        df = (c3.merge(c7, on="Patient_ID").merge(m3, on="Patient_ID")
                .merge(m7, on="Patient_ID"))
        # harmonisation: clinical baseline split vs MRI split for the same patients
        n_mismatch += int((df["split"] != df["split_mri"]).sum())
        df["y"] = df["Patient_ID"].map(y_by_pid)
        df = df[df["y"].notna()].copy()
        df["y"] = df["y"].astype(int)

        bcl = _buckets(df["p_cl3"].values, df["p_cl7"].values)
        bmr = _buckets(df["p_mri3"].values, df["p_mri7"].values)
        vm, tm = (df["split"] == "val").values, (df["split"] == "test").values
        yv, yt = df["y"].values[vm], df["y"].values[tm]
        clv, clt = bcl[vm], bcl[tm]
        mrv, mrt = bmr[vm], bmr[tm]
        if len(yt) == 0 or len(yv) == 0:
            print(f"  seed {seed}: empty val/test (val={len(yv)} test={len(yt)}) — skip"); continue

        pid_t = df["Patient_ID"].values[tm]

        def add(variant, pred, proba, val_bacc):
            r = h.metric_row(yt, pred, proba)
            r.update(seed=seed, variant=variant, val_bacc=float(val_bacc))
            rows.append(r)
            for i in range(len(yt)):           # capture per-patient pred for POOLED bACC
                all_preds.append(dict(seed=seed, Patient_ID=pid_t[i], variant=variant,
                                      y_true=int(yt[i]), pred=int(pred[i])))

        from sklearn.metrics import balanced_accuracy_score as bacc
        cl_pred, mri_pred = clt.argmax(1), mrt.argmax(1)
        add("clinical_only", cl_pred, clt, bacc(yv, clv.argmax(1)))
        add("mri_only",      mri_pred, mrt, bacc(yv, mrv.argmax(1)))
        eq = 0.5 * clt + 0.5 * mrt
        add("equal_0.5", eq.argmax(1), eq, bacc(yv, (0.5 * clv + 0.5 * mrv).argmax(1)))
        wpr, wbl, ww = h.fit_weighted_avg(clv, mrv, yv, clt, mrt, "bacc")
        add("weighted_avg", wpr, wbl, bacc(yv, (ww * clv + (1 - ww) * mrv).argmax(1)))
        wavg_pred = wpr                                  # best fusion -> capture for plots
        pr, bl, w = h.fit_weighted_avg(clv, mrv, yv, clt, mrt, "youden")
        add("weighted_avg_J", pr, bl, bacc(yv, (w * clv + (1 - w) * mrv).argmax(1)))
        Xv, Xt = np.hstack([clv, mrv]), np.hstack([clt, mrt])
        Xv_all.append(Xv); yv_all.append(yv)
        en_res = h.fit_en(Xv, yv, Xt, seed)
        if en_res is not None:
            add("elastic_net", en_res[0], en_res[1], float("nan"))
        lr_res = h.fit_lr(Xv, yv, Xt, seed)
        if lr_res is not None:
            add("lr", lr_res[0], lr_res[1], float("nan"))
        lr_pred = lr_res[0] if lr_res is not None else np.full(len(yt), -1)
        for i in range(len(yt)):
            pred_rows.append(dict(
                seed=seed, Patient_ID=pid_t[i], y_true=int(yt[i]),
                clinical_pred=int(cl_pred[i]), clinical_correct=int(cl_pred[i] == yt[i]),
                mri_pred=int(mri_pred[i]),     mri_correct=int(mri_pred[i] == yt[i]),
                wavg_pred=int(wavg_pred[i]),   wavg_correct=int(wavg_pred[i] == yt[i]),
                lr_pred=int(lr_pred[i]),       lr_correct=int(lr_pred[i] == yt[i])))
        print(f"  seed {seed}: n_val={len(yv)} n_test={len(yt)} "
              f"clin_only_bACC={bacc(yt, cl_pred):.3f}  lr_bACC="
              f"{bacc(yt, lr_pred) if lr_res is not None else float('nan'):.3f}")

    if not rows:
        sys.exit("No T4 rows produced — check inputs.")
    per_seed = pd.DataFrame(rows)
    per_seed.to_csv(HERE / "fusion_metrics_per_seed.csv", index=False)

    # per-patient predictions (for failure plots) + all-variant (for pooled bACC)
    pd.DataFrame(pred_rows).to_csv(HERE / "fusion_predictions_T4.csv", index=False)
    pd.DataFrame(all_preds).to_csv(HERE / "fusion_predictions_all_T4.csv", index=False)

    # LR weights — fit one multinomial LR on POOLED VAL (illustrative; the
    # per-seed eval LRs vary, but pooled gives a clean 3×6 coefficient matrix).
    from sklearn.linear_model import LogisticRegression
    Xpool, ypool = np.vstack(Xv_all), np.concatenate(yv_all)
    coef_rows = []
    if len(np.unique(ypool)) >= 2:
        lr = LogisticRegression(max_iter=5000, random_state=0).fit(Xpool, ypool)
        coef = lr.coef_                      # (n_classes_seen, 6); binary -> (1,6)
        cls = list(lr.classes_)
        for ri in range(coef.shape[0]):
            cname = h.CLASS_NAMES[cls[ri]] if coef.shape[0] == len(cls) else f"cls{cls[-1]}"
            for fi, fname in enumerate(FEATS):
                coef_rows.append(dict(class_name=cname, feature=fname, coef=float(coef[ri, fi])))
        pd.DataFrame(coef_rows).to_csv(HERE / "lr_coefficients_T4.csv", index=False)
        print(f"  pooled-VAL LR fit on n={len(ypool)} (classes {cls}) -> lr_coefficients_T4.csv")

    metric_cols = ["bacc", "macro_auc", "macro_f1"] + [f"f1_{c}" for c in h.CLASS_NAMES] + ["val_bacc", "n"]
    agg = (per_seed.groupby("variant")[metric_cols]
           .agg(["mean", "std"]).reset_index())
    # flatten + write
    summary = pd.DataFrame({"variant": agg["variant"]})
    for c in metric_cols:
        summary[f"{c}_mean"] = agg[(c, "mean")]
        summary[f"{c}_std"] = agg[(c, "std")]
    order = ["clinical_only", "mri_only", "equal_0.5", "weighted_avg",
             "weighted_avg_J", "elastic_net", "lr"]
    summary["__o"] = summary["variant"].map({v: i for i, v in enumerate(order)}).fillna(99)
    summary = summary.sort_values("__o").drop(columns="__o").reset_index(drop=True)
    summary.to_csv(HERE / "fusion_metrics.csv", index=False)

    print(f"\n  harmonisation: clinical-vs-MRI split mismatches across seeds = {n_mismatch} "
          f"(expect 0 if baseline splits align)")
    print("\n[T4 cumulative-horizon fusion] mean over seeds (fit VAL / report TEST):")
    show = summary[["variant", "bacc_mean", "bacc_std", "macro_auc_mean",
                    "macro_f1_mean", "val_bacc_mean", "n_mean"]].copy()
    with pd.option_context("display.float_format", lambda v: f"{v:.3f}"):
        print(show.to_string(index=False))
    print(f"\n  wrote {HERE/'fusion_metrics.csv'} + fusion_metrics_per_seed.csv")
    print("  NOTE: T4 test n≈12-15 (converters only) — high variance, indicative only.")


if __name__ == "__main__":
    main()
