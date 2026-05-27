"""
47_prs_aao_regression_strictQC.py
=================================
AAO Phase 2c.C — Per-source PRS -> canonical AAO linear regression on the
LD-pruned strict-QC pool.

Cohort: CONVERTERS ONLY — Phase 1 AAO splits.
Target: canonical AAO.
COVAR_MODES same as classification (script 45).
SOURCES: 16 individual + "prs_all_dedup".
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
from scipy import stats
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.pipeline import Pipeline

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from sklearn.linear_model import ElasticNet
from _prs_strict_qc_lib import (ALL_SOURCES, per_source_prs_table,
                                  load_subject_covariates,
                                  get_dedup_dosage_matrix,
                                  DEFAULT_LD_CONFIG)  # noqa: E402

AAO_SPLITS_ROOT = Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                        "no_cdr_stratified_post_exclusion_aao")
OUT_BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220/outputs/strict_qc_prs")
SEEDS = (0, 1, 2)
# Covariate modes are NAMED identically to scripts 45/46 (for join-able master
# leaderboard rows) but the AAO regression INTERNALLY drops `age_at_baseline`
# from each set — because AAO = age_at_baseline + years_to_AD, so using age
# as a feature is circular. The dropped age is noted on the table footer.
COVAR_MODES = {
    "prs_only":                          [],
    "prs+age+sex+apoe4":                 ["Sex_M", "APOE4_Dosage"],          # age dropped
    "prs+age+sex+apoe4+apoe2":           ["Sex_M", "APOE4_Dosage", "APOE2_Dosage"],  # age dropped
}


def _load_aao_split(seed: int, split: str) -> pd.DataFrame:
    df = pd.read_csv(AAO_SPLITS_ROOT / f"seed_{seed}/{split}.csv", dtype=str)
    df["Patient_ID"] = df["Patient_ID"].astype(str)
    df["AAO"] = pd.to_numeric(df["AAO"], errors="coerce")
    return df[["Patient_ID","AAO"]].dropna(subset=["AAO"])


def _metrics(y_true, y_pred) -> dict:
    if len(y_true) < 2:
        return {"mse":float("nan"),"mae":float("nan"),"rmse":float("nan"),
                 "r2":float("nan"),"pearson_r":float("nan"),"pearson_p":float("nan"),
                 "n":int(len(y_true))}
    mse = mean_squared_error(y_true, y_pred)
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    rmse = math.sqrt(mse)
    if np.std(y_pred) > 1e-9:
        r, p = stats.pearsonr(y_true, y_pred)
    else:
        r, p = float("nan"), float("nan")
    return {"mse":float(mse),"mae":float(mae),"rmse":float(rmse),
             "r2":float(r2),"pearson_r":float(r),"pearson_p":float(p),
             "n":int(len(y_true))}


def _fit_en_aao_dosage(dos_feat: pd.DataFrame, parts: dict,
                          seed: int) -> pd.Series:
    train_y = parts["train"].set_index("Patient_ID")["AAO"].astype(float)
    common = [p for p in train_y.index if p in dos_feat.index]
    if len(common) < 30: return pd.Series(np.nan, index=dos_feat.index)
    X = dos_feat.fillna(dos_feat.mean()).astype(float)
    try:
        en = ElasticNet(alpha=0.1, l1_ratio=0.5, max_iter=5000, random_state=seed)
        en.fit(X.loc[common].values, train_y.loc[common].values)
        scores = en.predict(X.values)
        return pd.Series(scores, index=X.index)
    except Exception:
        return pd.Series(np.nan, index=X.index)


def _fit_en_aao_meta_prs(source_prs: pd.DataFrame, parts: dict,
                            seed: int) -> pd.Series:
    train_y = parts["train"].set_index("Patient_ID")["AAO"].astype(float)
    common = [p for p in train_y.index if p in source_prs.index]
    if len(common) < 30: return pd.Series(np.nan, index=source_prs.index)
    X = source_prs.copy()
    keep_cols = [c for c in X.columns if X[c].notna().any()]
    X = X[keep_cols].fillna(X.mean()).astype(float)
    try:
        en = ElasticNet(alpha=0.1, l1_ratio=0.5, max_iter=5000, random_state=seed)
        en.fit(X.loc[common].values, train_y.loc[common].values)
        scores = en.predict(X.values)
        return pd.Series(scores, index=X.index)
    except Exception:
        return pd.Series(np.nan, index=X.index)


def _run_one(src: str, prs_full: pd.DataFrame, cov: pd.DataFrame,
              seed: int, mode: str, covars: list,
              dedup_dosage: pd.DataFrame | None = None) -> dict | None:
    parts = {sp: _load_aao_split(seed, sp) for sp in ("train","val","test")}
    if src == "prs_all_dedup_EN_dosage":
        if dedup_dosage is None: return None
        en_prs = _fit_en_aao_dosage(dedup_dosage, parts, seed)
        if en_prs.isna().all(): return None
        base = pd.DataFrame({"Patient_ID": en_prs.index.astype(str),
                              "PRS": en_prs.values})
    elif src == "meta_prs_EN_combined":
        source_cols = [c for c in prs_full.columns if c.startswith("PRS_")
                        and c not in ("PRS_prs_all_dedup",)]
        sp_df = prs_full[["PTID"] + source_cols].set_index("PTID")
        sp_df.index = sp_df.index.astype(str)
        en_prs = _fit_en_aao_meta_prs(sp_df, parts, seed)
        if en_prs.isna().all(): return None
        base = pd.DataFrame({"Patient_ID": en_prs.index, "PRS": en_prs.values})
    else:
        prs_col = f"PRS_{src}"
        if prs_col not in prs_full.columns or prs_full[prs_col].isna().all():
            return None
        base = prs_full[["PTID", prs_col]].rename(columns={"PTID":"Patient_ID",
                                                              prs_col:"PRS"})
    base["Patient_ID"] = base["Patient_ID"].astype(str)
    fits = {}
    for sp in ("train","val","test"):
        f = parts[sp].merge(base, on="Patient_ID", how="inner")
        f = f.merge(cov, on="Patient_ID", how="left")
        fits[sp] = f
    tr = fits["train"]
    if len(tr) < 30: return None
    mu = tr["PRS"].mean(); sd = tr["PRS"].std(ddof=0) or 1.0
    for sp in fits:
        fits[sp] = fits[sp].copy()
        fits[sp]["PRS_z"] = (fits[sp]["PRS"] - mu) / sd
    Xcols = ["PRS_z"] + covars
    needed = Xcols + ["AAO"]
    tr_in = fits["train"][needed].dropna()
    val_in = fits["val"][needed].dropna()
    test_in = fits["test"][needed].dropna()
    if len(tr_in) < 10 or len(val_in) < 3: return None
    pipe = Pipeline([("imp", SimpleImputer(strategy="median")),
                      ("reg", LinearRegression())])
    pipe.fit(tr_in[Xcols].astype(float), tr_in["AAO"].astype(float))
    val_pred = pipe.predict(val_in[Xcols].astype(float))
    test_pred = pipe.predict(test_in[Xcols].astype(float))
    val_m = _metrics(val_in["AAO"].values, val_pred)
    test_m = _metrics(test_in["AAO"].values, test_pred)
    coef = float(pipe.named_steps["reg"].coef_[0])  # PRS_z coef → years per SD
    out = {"source": src, "covar_mode": mode, "seed": seed,
            "n_train": len(tr_in), "n_val": len(val_in), "n_test": len(test_in),
            "yrs_per_sd_prs": coef}
    for k, v in val_m.items():
        out[f"val_{k}"] = v
    for k, v in test_m.items():
        out[f"test_{k}"] = v
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=None)
    ap.add_argument("--ld-config", default=DEFAULT_LD_CONFIG)
    ap.add_argument("--beta-source", default="raw", choices=["raw","prscs"],
                     help="raw = published lead-SNP β; prscs = PRS-CS posterior β.")
    args = ap.parse_args()
    if args.beta_source == "raw":
        out_dir = OUT_BASE / args.ld_config / "aao_regression"
    else:
        out_dir = OUT_BASE.parent / f"strict_qc_prs_{args.beta_source}" / args.ld_config / "aao_regression"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Building per-source PRS table for {args.ld_config}  beta_source={args.beta_source}...")
    prs_full, snps_per_src = per_source_prs_table(ld_config=args.ld_config,
                                                     beta_source=args.beta_source)
    cov = load_subject_covariates()
    dedup_st, dedup_dosage = get_dedup_dosage_matrix(args.ld_config)
    n_snp_per_src = {s: len(st) for s, st in snps_per_src.items()}
    n_snp_per_src["prs_all_dedup_EN_dosage"] = dedup_dosage.shape[1]
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

    g = df.groupby(["source","covar_mode"]).agg(
        n_snps=("n_snps_used","first"),
        n_seeds=("seed","count"),
        val_r2_mean=("val_r2","mean"), val_r2_std=("val_r2","std"),
        val_rmse_mean=("val_rmse","mean"),
        val_pearson_r_mean=("val_pearson_r","mean"),
        val_pearson_p_min=("val_pearson_p","min"),
        yrs_per_sd_mean=("yrs_per_sd_prs","mean"), yrs_per_sd_std=("yrs_per_sd_prs","std"),
    ).reset_index().sort_values(["covar_mode","val_r2_mean"], ascending=[True, False])
    lb_p = out_dir / "leaderboard.csv"
    g.to_csv(lb_p, index=False)
    print(f"wrote leaderboard: {lb_p}")

    show = pd.DataFrame({
        "source":     g["source"],
        "covar_mode": g["covar_mode"],
        "n_snps":     g["n_snps"].astype(int),
        "val_r2":     [f"{m:.3f}+/-{s:.3f}" for m, s in zip(g["val_r2_mean"], g["val_r2_std"])],
        "val_rmse":   [f"{x:.2f}" for x in g["val_rmse_mean"]],
        "yrs/+1SD":   [f"{m:+.2f}+/-{s:.2f}" for m, s in zip(g["yrs_per_sd_mean"], g["yrs_per_sd_std"])],
        "min_p":      [f"{p:.2f}" for p in g["val_pearson_p_min"]],
    })
    fig_h = max(2.5, 0.30 * len(show) + 1.3)
    fig, ax = plt.subplots(figsize=(14, fig_h)); ax.axis("off")
    ax.set_title(f"AAO regression — strict-QC + LD pruning [{args.ld_config}]   "
                 f"(VAL set; mean +/- std over 3 seeds)",
                 fontsize=11, fontweight="bold", pad=12, loc="left")
    col_widths = [0.28, 0.18, 0.07, 0.16, 0.10, 0.16, 0.08]
    tbl = ax.table(cellText=show.values.tolist(), colLabels=show.columns.tolist(),
                    loc="center", cellLoc="center", colWidths=col_widths)
    tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1, 1.2)
    for j in range(len(show.columns)):
        c = tbl[(0, j)]; c.set_facecolor("#222"); c.set_text_props(color="white", weight="bold")
    png = out_dir / "leaderboard.png"
    fig.savefig(png, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"wrote PNG: {png}")


if __name__ == "__main__":
    main()
