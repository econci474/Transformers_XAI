"""
46_prs_cox_strictQC.py
======================
AAO Phase 2c.B — Cox proportional-hazards survival on LD-pruned strict-QC PRS.

Cohort: ALL subjects with FU_years > 0 in the standard
no_cdr_stratified_post_exclusion splits.

Per (source x covar_mode x seed):
  PRS_p = Sigma beta_A1 * dosage_p  (over LD-pruned strict-QC SNPs)
  z-score on train fold
  lifelines CoxPHFitter on {PRS_z, *covariates}
  Eval on val: concordance, HR per 1-SD with 95% CI + p

COVAR_MODES same as classification (script 45).
SOURCES: 16 individual + "prs_all_dedup".
"""
from __future__ import annotations
import argparse
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

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
EPS = 1e-3  # epsilon for AD_bl baseline-converters


def _load_time_event(seed: int, split: str) -> pd.DataFrame:
    df = pd.read_csv(SPLITS_ROOT / f"seed_{seed}/{split}.csv", dtype=str)
    df["Patient_ID"] = df["Patient_ID"].astype(str)
    df["AD_bl"] = pd.to_numeric(df["AD_bl"], errors="coerce")
    df["pMCI"]  = pd.to_numeric(df["pMCI"], errors="coerce")
    df["CN_to_AD"] = pd.to_numeric(df["CN_to_AD"], errors="coerce")
    df["AD_final"] = pd.to_numeric(df.get("AD_final", 0), errors="coerce").fillna(0)
    df["years_to_AD"] = pd.to_numeric(df["years_to_AD"], errors="coerce")
    df["FU_years"] = pd.to_numeric(df["FU_years"], errors="coerce")
    is_conv = ((df["AD_bl"] == 1) | (df["pMCI"] == 1) | (df["CN_to_AD"] == 1))
    t = np.where(df["AD_bl"] == 1, EPS,
                  np.where(is_conv, df["years_to_AD"], df["FU_years"]))
    e = (df["AD_final"] == 1).astype(int).values
    df["t"] = t.astype(float); df["e"] = e
    df = df[(df["t"] > 0) & df["t"].notna()].copy()
    return df[["Patient_ID","t","e"]]


def _fit_en_cox_dosage_get_prs(dos_feat: pd.DataFrame, parts: dict,
                                  seed: int) -> pd.Series:
    """Fit ElasticNet Cox via lifelines (penalizer + l1_ratio) on TRAIN
    fold dosage with (t, e). Return log-partial-hazard for all subjects."""
    train_te = parts["train"].set_index("Patient_ID")
    common = [p for p in train_te.index if p in dos_feat.index]
    if len(common) < 30 or train_te.loc[common, "e"].sum() < 5:
        return pd.Series(np.nan, index=dos_feat.index)
    X = dos_feat.fillna(dos_feat.mean()).astype(float)
    df_tr = X.loc[common].copy()
    df_tr["t"] = train_te.loc[common, "t"].values
    df_tr["e"] = train_te.loc[common, "e"].astype(int).values
    try:
        cph = CoxPHFitter(penalizer=0.1, l1_ratio=0.5)
        cph.fit(df_tr, duration_col="t", event_col="e", show_progress=False)
        # log partial hazard = X @ beta — lifelines exposes predict_log_partial_hazard
        scores = cph.predict_log_partial_hazard(X)
        return pd.Series(scores.values, index=X.index)
    except Exception:
        return pd.Series(np.nan, index=X.index)


def _fit_en_cox_meta_prs(source_prs: pd.DataFrame, parts: dict,
                            seed: int) -> pd.Series:
    train_te = parts["train"].set_index("Patient_ID")
    common = [p for p in train_te.index if p in source_prs.index]
    if len(common) < 30 or train_te.loc[common, "e"].sum() < 5:
        return pd.Series(np.nan, index=source_prs.index)
    X = source_prs.copy()
    keep_cols = [c for c in X.columns if X[c].notna().any()]
    X = X[keep_cols].fillna(X.mean())
    df_tr = X.loc[common].copy()
    df_tr["t"] = train_te.loc[common, "t"].values
    df_tr["e"] = train_te.loc[common, "e"].astype(int).values
    try:
        cph = CoxPHFitter(penalizer=0.1, l1_ratio=0.5)
        cph.fit(df_tr, duration_col="t", event_col="e", show_progress=False)
        scores = cph.predict_log_partial_hazard(X)
        return pd.Series(scores.values, index=X.index)
    except Exception:
        return pd.Series(np.nan, index=X.index)


def _run_one(src: str, prs_full: pd.DataFrame, cov: pd.DataFrame,
              seed: int, mode: str, covars: list,
              dedup_dosage: pd.DataFrame | None = None) -> dict | None:
    parts = {sp: _load_time_event(seed, sp) for sp in ("train","val","test")}
    if src == "prs_all_dedup_EN_dosage":
        if dedup_dosage is None: return None
        en_prs = _fit_en_cox_dosage_get_prs(dedup_dosage, parts, seed)
        if en_prs.isna().all(): return None
        base = pd.DataFrame({"Patient_ID": en_prs.index.astype(str),
                              "PRS": en_prs.values})
    elif src == "meta_prs_EN_combined":
        source_cols = [c for c in prs_full.columns if c.startswith("PRS_")
                        and c not in ("PRS_prs_all_dedup",)]
        sp_df = prs_full[["PTID"] + source_cols].set_index("PTID")
        sp_df.index = sp_df.index.astype(str)
        en_prs = _fit_en_cox_meta_prs(sp_df, parts, seed)
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
    if tr.empty or tr["e"].sum() < 5 or len(tr) < 30:
        return None
    mu = tr["PRS"].mean(); sd = tr["PRS"].std(ddof=0) or 1.0
    for sp in fits:
        fits[sp] = fits[sp].copy()
        fits[sp]["PRS_z"] = (fits[sp]["PRS"] - mu) / sd
    Xcols = ["PRS_z"] + covars
    needed = ["t","e"] + Xcols
    tr_in = fits["train"][needed].dropna()
    val_in = fits["val"][needed].dropna()
    test_in = fits["test"][needed].dropna()
    if tr_in.empty or tr_in["e"].sum() < 5: return None
    cph = CoxPHFitter(penalizer=1e-3)
    try:
        cph.fit(tr_in, duration_col="t", event_col="e")
    except Exception:
        return None

    out = {"source": src, "covar_mode": mode, "seed": seed,
            "n_train": len(tr_in), "events_train": int(tr_in["e"].sum()),
            "n_val": len(val_in), "events_val": int(val_in["e"].sum() if len(val_in) else 0),
            "n_test": len(test_in), "events_test": int(test_in["e"].sum() if len(test_in) else 0)}
    try:
        sm = cph.summary
        out["hr_per_sd"] = float(sm.loc["PRS_z","exp(coef)"])
        out["hr_lo"]    = float(sm.loc["PRS_z","exp(coef) lower 95%"])
        out["hr_hi"]    = float(sm.loc["PRS_z","exp(coef) upper 95%"])
        out["hr_p"]     = float(sm.loc["PRS_z","p"])
    except Exception:
        for k in ("hr_per_sd","hr_lo","hr_hi","hr_p"): out[k] = float("nan")

    def _cidx(p):
        if len(p) < 8 or p["e"].sum() < 2: return float("nan")
        pi = -cph.predict_partial_hazard(p[Xcols]).values
        try:    return float(concordance_index(p["t"], pi, p["e"]))
        except: return float("nan")
    out["val_cindex"] = _cidx(val_in)
    out["test_cindex"] = _cidx(test_in)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=None)
    ap.add_argument("--ld-config", default=DEFAULT_LD_CONFIG)
    ap.add_argument("--beta-source", default="raw", choices=["raw","prscs"],
                     help="raw = published lead-SNP β; prscs = PRS-CS posterior β.")
    args = ap.parse_args()
    if args.beta_source == "raw":
        out_dir = OUT_BASE / args.ld_config / "cox"
    else:
        out_dir = OUT_BASE.parent / f"strict_qc_prs_{args.beta_source}" / args.ld_config / "cox"
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
        val_cindex_mean=("val_cindex","mean"), val_cindex_std=("val_cindex","std"),
        hr_per_sd_mean=("hr_per_sd","mean"), hr_per_sd_std=("hr_per_sd","std"),
        hr_lo_mean=("hr_lo","mean"), hr_hi_mean=("hr_hi","mean"),
        hr_p_min=("hr_p","min"),
        test_cindex_mean=("test_cindex","mean"),
    ).reset_index().sort_values(["covar_mode","val_cindex_mean"], ascending=[True, False])
    lb_p = out_dir / "leaderboard.csv"
    g.to_csv(lb_p, index=False)
    print(f"wrote leaderboard: {lb_p}")

    # PNG — concise
    show = pd.DataFrame({
        "source":     g["source"],
        "covar_mode": g["covar_mode"],
        "n_snps":     g["n_snps"].astype(int),
        "val_cindex": [f"{m:.3f}+/-{s:.3f}" for m, s in zip(g["val_cindex_mean"], g["val_cindex_std"])],
        "HR/+1SD":    [f"{m:.2f} [{lo:.2f}-{hi:.2f}]"
                       for m, lo, hi in zip(g["hr_per_sd_mean"], g["hr_lo_mean"], g["hr_hi_mean"])],
        "min_p":      [f"{p:.1e}" for p in g["hr_p_min"]],
    })
    fig_h = max(2.5, 0.30 * len(show) + 1.3)
    fig, ax = plt.subplots(figsize=(14, fig_h)); ax.axis("off")
    ax.set_title(f"Cox survival — strict-QC + LD pruning [{args.ld_config}]   "
                 f"(VAL set; mean +/- std over 3 seeds)",
                 fontsize=11, fontweight="bold", pad=12, loc="left")
    col_widths = [0.30, 0.18, 0.07, 0.16, 0.20, 0.12]
    tbl = ax.table(cellText=show.values.tolist(), colLabels=show.columns.tolist(),
                    loc="center", cellLoc="center", colWidths=col_widths)
    tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1, 1.2)
    for j in range(len(show.columns)):
        c = tbl[(0, j)]; c.set_facecolor("#222"); c.set_text_props(color="white", weight="bold")
    # Bold C-index > 0.5
    for i, row in enumerate(show.itertuples(index=False), start=1):
        mean_c = float(row.val_cindex.split("+/-")[0])
        if mean_c > 0.5:
            tbl[(i, 3)].set_text_props(weight="bold")
    png = out_dir / "leaderboard.png"
    fig.savefig(png, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"wrote PNG: {png}")


if __name__ == "__main__":
    main()
