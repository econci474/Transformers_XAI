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
from _prs_strict_qc_lib import (ALL_SOURCES, per_source_prs_table,
                                  load_subject_covariates,
                                  get_dedup_dosage_matrix,
                                  DEFAULT_LD_CONFIG)  # noqa: E402

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
                              cov: pd.DataFrame, seed: int) -> pd.Series:
    """Train ElasticNet LogisticRegression on TRAIN-fold dosage to predict
    AD. Return decision function (log-odds) for all subjects in dos_feat,
    indexed by PTID. Same EN model applied to val/test ⇒ no leakage."""
    train_pids = parts["train"]["Patient_ID"].tolist()
    train_labels = parts["train"].set_index("Patient_ID")["y"].astype(int)
    common = [p for p in train_pids if p in dos_feat.index]
    if len(common) < 30 or train_labels.loc[common].nunique() < 2:
        return pd.Series(np.nan, index=dos_feat.index)
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
        return pd.Series(np.nan, index=dos_feat.index)
    return pd.Series(log_odds, index=dos_feat.index)


def _fit_en_meta_prs(source_prs: pd.DataFrame, parts: dict,
                       seed: int) -> pd.Series:
    """source_prs: PTID-indexed DataFrame with one column per source ("PRS_<src>").
    Train ElasticNet LR on TRAIN-fold source-PRS values to predict AD.
    Return decision function for all subjects."""
    train_pids = parts["train"]["Patient_ID"].tolist()
    train_labels = parts["train"].set_index("Patient_ID")["y"].astype(int)
    common = [p for p in train_pids if p in source_prs.index]
    if len(common) < 30 or train_labels.loc[common].nunique() < 2:
        return pd.Series(np.nan, index=source_prs.index)
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
        return pd.Series(np.nan, index=source_prs.index)
    return pd.Series(log_odds, index=source_prs.index)


def _run_one(src: str, prs_full: pd.DataFrame, cov: pd.DataFrame,
              seed: int, mode: str, covars: list,
              dedup_dosage: pd.DataFrame | None = None) -> dict | None:
    parts = {sp: _load_labels(seed, sp) for sp in ("train", "val", "test")}
    # Build the "PRS" column per source. For non-EN sources, use the precomputed
    # PRS table. For EN sources, fit per-seed on train fold and produce log-odds
    # as the "PRS".
    if src == "prs_all_dedup_EN_dosage":
        if dedup_dosage is None: return None
        en_prs = _fit_en_dosage_get_prs(dedup_dosage, parts, cov, seed)
        if en_prs.isna().all(): return None
        base = pd.DataFrame({"Patient_ID": en_prs.index.astype(str),
                              "PRS": en_prs.values})
    elif src == "meta_prs_EN_combined":
        source_cols = [c for c in prs_full.columns if c.startswith("PRS_")
                        and c not in ("PRS_prs_all_dedup",)
                        and not c.endswith("EN_dosage")
                        and not c.endswith("EN_combined")]
        sp_df = prs_full[["PTID"] + source_cols].set_index("PTID")
        sp_df.index = sp_df.index.astype(str)
        en_prs = _fit_en_meta_prs(sp_df, parts, seed)
        if en_prs.isna().all(): return None
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
    val_p  = pipe.predict_proba(val_in[Xcols].astype(float))[:, 1]
    test_p = pipe.predict_proba(test_in[Xcols].astype(float))[:, 1]
    val_m  = _metrics(val_in["y"], val_p)
    test_m = _metrics(test_in["y"], test_p)

    # OR per 1-SD of PRS_z (which is already 1-SD scaled on train); CI via Wald.
    coef = float(pipe.named_steps["clf"].coef_[0][0])
    or_per_sd = math.exp(coef)
    # Approximate Wald CI requires SE; for small-n we leave CI fields blank.
    out = {"source": src, "covar_mode": mode, "seed": seed,
            "or_per_sd": or_per_sd}
    for k, v in val_m.items():
        out[f"val_{k}"] = v
    for k, v in test_m.items():
        out[f"test_{k}"] = v
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=None)
    ap.add_argument("--ld-config", default=DEFAULT_LD_CONFIG)
    args = ap.parse_args()
    out_dir = OUT_BASE / args.ld_config / "classification"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Building per-source PRS table for {args.ld_config}...")
    prs_full, snps_per_src = per_source_prs_table(ld_config=args.ld_config)
    cov = load_subject_covariates()
    dedup_st, dedup_dosage = get_dedup_dosage_matrix(args.ld_config)
    n_snp_per_src = {s: len(st) for s, st in snps_per_src.items()}
    n_snp_per_src["prs_all_dedup_EN_dosage"] = dedup_dosage.shape[1]
    # meta_prs uses the union of source PRS values; "n_snps" is total
    # unique rsIDs in the LD-pruned pool (proxy — same as dedup size).
    n_snp_per_src["meta_prs_EN_combined"] = dedup_dosage.shape[1]
    all_with_dedup = (ALL_SOURCES + ["prs_all_dedup",
                                       "prs_all_dedup_EN_dosage",
                                       "meta_prs_EN_combined"])
    sources = [args.source] if args.source else [s for s in all_with_dedup
                                                    if n_snp_per_src.get(s, 0) > 0]
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
        "OR/+1SD":    [f"{x:.2f}" for x in g["or_per_sd_mean"]],
    })
    fig_h = max(2.5, 0.30 * len(show) + 1.3)
    fig, ax = plt.subplots(figsize=(14, fig_h)); ax.axis("off")
    ax.set_title(f"Classification — strict-QC + LD pruning [{args.ld_config}]   "
                 f"(VAL set; mean +/- std over 3 seeds)",
                 fontsize=11, fontweight="bold", pad=12, loc="left")
    col_widths = [0.30, 0.18, 0.07, 0.16, 0.16, 0.16, 0.10]
    tbl = ax.table(cellText=show.values.tolist(), colLabels=show.columns.tolist(),
                    loc="center", cellLoc="center", colWidths=col_widths)
    tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1, 1.2)
    for j in range(len(show.columns)):
        c = tbl[(0, j)]; c.set_facecolor("#222"); c.set_text_props(color="white", weight="bold")
    # Bold BalAcc > 0.5
    for i, row in enumerate(show.itertuples(index=False), start=1):
        mean_bacc = float(row.val_bacc.split("+/-")[0])
        if mean_bacc > 0.5:
            tbl[(i, 3)].set_text_props(weight="bold")
    png = out_dir / "leaderboard.png"
    fig.savefig(png, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"wrote PNG: {png}")


if __name__ == "__main__":
    main()
