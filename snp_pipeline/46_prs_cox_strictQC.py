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
from sklearn.metrics import roc_auc_score

TIME_HORIZONS = (3, 5, 10)  # years

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _prs_strict_qc_lib import (ALL_SOURCES, PRSCS_SOURCES, per_source_prs_table,
                                  load_subject_covariates,
                                  get_dedup_dosage_matrix,
                                  render_en_weights_png,
                                  META_FILTERED_SOURCES,
                                  DEFAULT_LD_CONFIG)  # noqa: E402
from _per_patient_io import (capture_predictions, save_per_patient_parquet,
                              concat_and_save_anchor)  # noqa: E402

EN_COEFS = {"dosage": {}, "meta": {}, "meta_filtered": {}}

# Per-patient capture targets (matches script 45). outcome_proba in Cox is
# the partial hazard, not a probability; consumer code uses it as a relative
# risk score for ranking + fusion.
PER_PATIENT_SOURCES = {
    "meta_prs_EN_combined",
    "Kosteridis",
    "Kosteridis_MTAG_AD",
    "Kosteridis_shared_AD_CV",
}
PER_PATIENT_ROWS: list[pd.DataFrame] = []
PER_PATIENT_BETA_SOURCE: list[str] = ["raw"]
PER_PATIENT_LD_CONFIG:   list[str] = ["na"]

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
                                  seed: int):
    """Fit EN Cox on TRAIN dosage. Return (model, log-partial-hazard-series).
    model=None when fit skipped/failed."""
    train_te = parts["train"].set_index("Patient_ID")
    common = [p for p in train_te.index if p in dos_feat.index]
    if len(common) < 30 or train_te.loc[common, "e"].sum() < 5:
        return None, pd.Series(np.nan, index=dos_feat.index)
    X = dos_feat.fillna(dos_feat.mean()).astype(float)
    df_tr = X.loc[common].copy()
    df_tr["t"] = train_te.loc[common, "t"].values
    df_tr["e"] = train_te.loc[common, "e"].astype(int).values
    try:
        cph = CoxPHFitter(penalizer=0.1, l1_ratio=0.5)
        cph.fit(df_tr, duration_col="t", event_col="e", show_progress=False)
        scores = cph.predict_log_partial_hazard(X)
        return cph, pd.Series(scores.values, index=X.index)
    except Exception:
        return None, pd.Series(np.nan, index=X.index)


def _fit_en_cox_meta_prs(source_prs: pd.DataFrame, parts: dict,
                            seed: int):
    train_te = parts["train"].set_index("Patient_ID")
    common = [p for p in train_te.index if p in source_prs.index]
    if len(common) < 30 or train_te.loc[common, "e"].sum() < 5:
        return None, pd.Series(np.nan, index=source_prs.index)
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
        cph._feature_names = keep_cols
        return cph, pd.Series(scores.values, index=X.index)
    except Exception:
        return None, pd.Series(np.nan, index=X.index)


def _run_covariates_only(seed: int, mode: str, covars: list,
                           cov: pd.DataFrame) -> dict | None:
    """Clinical baseline: fit Cox(t,e ~ covariates) with NO PRS.

    HR-per-SD-PRS fields are NaN; c-index from the covariates-only model is
    reported on val + test."""
    if not covars: return None
    parts = {sp: _load_time_event(seed, sp) for sp in ("train","val","test")}
    fits = {sp: parts[sp].merge(cov, on="Patient_ID", how="left")
             for sp in ("train","val","test")}
    needed = ["t","e"] + covars
    tr_in   = fits["train"][needed].dropna()
    val_in  = fits["val"][  needed].dropna()
    test_in = fits["test"][ needed].dropna()
    if tr_in.empty or tr_in["e"].sum() < 5: return None
    cph = CoxPHFitter(penalizer=1e-3)
    try:
        cph.fit(tr_in, duration_col="t", event_col="e")
    except Exception:
        return None
    out = {"source": "covariates_only", "covar_mode": mode, "seed": seed,
            "n_train": len(tr_in), "events_train": int(tr_in["e"].sum()),
            "n_val": len(val_in), "events_val": int(val_in["e"].sum() if len(val_in) else 0),
            "n_test": len(test_in), "events_test": int(test_in["e"].sum() if len(test_in) else 0)}
    for k in ("hr_per_sd","hr_lo","hr_hi","hr_p",
                "hr_per_sd_val","hr_lo_val","hr_hi_val","hr_p_val",
                "hr_per_sd_test","hr_lo_test","hr_hi_test","hr_p_test"):
        out[k] = float("nan")

    def _cidx(p):
        if len(p) < 8 or p["e"].sum() < 2: return float("nan")
        pi = -cph.predict_partial_hazard(p[covars]).values
        try:    return float(concordance_index(p["t"], pi, p["e"]))
        except: return float("nan")
    out["val_cindex"]  = _cidx(val_in)
    out["test_cindex"] = _cidx(test_in)
    # td-AUC + O/E (use covars-only predictions)
    def _td_auc_oe(df_in: pd.DataFrame, t_hor: float) -> tuple:
        if len(df_in) < 8: return float("nan"), float("nan")
        try:
            surv = cph.predict_survival_function(df_in[covars], times=[t_hor])
            p_event = 1.0 - surv.iloc[0].values
        except Exception:
            return float("nan"), float("nan")
        cases    = ((df_in["t"] <= t_hor) & (df_in["e"] == 1)).values
        controls = (df_in["t"]  > t_hor).values
        if cases.sum() < 1 or controls.sum() < 1:
            td_auc = float("nan")
        else:
            keep = cases | controls
            y    = np.where(cases, 1, 0)
            try:    td_auc = float(roc_auc_score(y[keep], p_event[keep]))
            except: td_auc = float("nan")
        n_obs = float(cases.sum()); n_exp = float(p_event.sum())
        oe    = (n_obs / n_exp) if n_exp > 1e-9 else float("nan")
        return td_auc, oe
    for split_name, df_in in [("val", val_in), ("test", test_in)]:
        for thor in TIME_HORIZONS:
            ta, oe = _td_auc_oe(df_in, float(thor))
            out[f"{split_name}_td_auc_{thor}y"] = ta
            out[f"{split_name}_oe_{thor}y"]    = oe
    return out


def _run_one(src: str, prs_full: pd.DataFrame, cov: pd.DataFrame,
              seed: int, mode: str, covars: list,
              dedup_dosage: pd.DataFrame | None = None) -> dict | None:
    if src == "covariates_only":
        return _run_covariates_only(seed, mode, covars, cov)
    parts = {sp: _load_time_event(seed, sp) for sp in ("train","val","test")}
    if src == "prs_all_dedup_EN_dosage":
        if dedup_dosage is None: return None
        cph_en, en_prs = _fit_en_cox_dosage_get_prs(dedup_dosage, parts, seed)
        if en_prs.isna().all(): return None
        if cph_en is not None:
            # cph.params_ is a Series indexed by feature name (rsID), already
            # excluding t/e duration/event columns.
            ps = cph_en.params_
            EN_COEFS["dosage"][seed] = ps[[c for c in ps.index
                                            if c not in ("t","e")]]
        base = pd.DataFrame({"Patient_ID": en_prs.index.astype(str),
                              "PRS": en_prs.values})
    elif src == "meta_prs_EN_combined":
        # Drop dedup columns and Kosteridis SUBGROUPS (Kosteridis full =
        # MTAG_AD ∪ shared_AD_CV) to avoid triple-counting the same Kosteridis
        # SNPs in the meta-EN.
        source_cols = [c for c in prs_full.columns if c.startswith("PRS_")
                        and c not in ("PRS_prs_all_dedup",
                                        "PRS_prs_all_dedup_ivw",
                                        "PRS_prs_all_dedup_filtered",
                                        "PRS_Kosteridis_MTAG_AD",
                                        "PRS_Kosteridis_shared_AD_CV",
                                        "PRS_Kosteridis_novel_AD")]
        sp_df = prs_full[["PTID"] + source_cols].set_index("PTID")
        sp_df.index = sp_df.index.astype(str)
        cph_en, en_prs = _fit_en_cox_meta_prs(sp_df, parts, seed)
        if en_prs.isna().all(): return None
        if cph_en is not None:
            ps = cph_en.params_
            ps = ps[[c for c in ps.index if c not in ("t","e")]]
            ps.index = [c.replace("PRS_", "", 1) for c in ps.index]
            EN_COEFS["meta"][seed] = ps
        base = pd.DataFrame({"Patient_ID": en_prs.index, "PRS": en_prs.values})
    elif src == "meta_prs_EN_filtered":
        source_cols = [f"PRS_{s}" for s in META_FILTERED_SOURCES
                        if f"PRS_{s}" in prs_full.columns
                        and prs_full[f"PRS_{s}"].notna().any()]
        if not source_cols: return None
        sp_df = prs_full[["PTID"] + source_cols].set_index("PTID")
        sp_df.index = sp_df.index.astype(str)
        cph_en, en_prs = _fit_en_cox_meta_prs(sp_df, parts, seed)
        if en_prs.isna().all(): return None
        if cph_en is not None:
            ps = cph_en.params_
            ps = ps[[c for c in ps.index if c not in ("t","e")]]
            ps.index = [c.replace("PRS_", "", 1) for c in ps.index]
            EN_COEFS["meta_filtered"][seed] = ps
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

    # Val-fold HR per +1 SD + 95% CI: refit Cox on VAL with PRS_z + covars.
    # Independent of TRAIN; removes train-fold data-snooping bias.
    for k in ("hr_per_sd_val","hr_lo_val","hr_hi_val","hr_p_val"):
        out[k] = float("nan")
    try:
        if len(val_in) >= 8 and val_in["e"].sum() >= 3:
            cph_v = CoxPHFitter(penalizer=1e-3)
            cph_v.fit(val_in, duration_col="t", event_col="e")
            smv = cph_v.summary
            out["hr_per_sd_val"] = float(smv.loc["PRS_z","exp(coef)"])
            out["hr_lo_val"]    = float(smv.loc["PRS_z","exp(coef) lower 95%"])
            out["hr_hi_val"]    = float(smv.loc["PRS_z","exp(coef) upper 95%"])
            out["hr_p_val"]     = float(smv.loc["PRS_z","p"])
    except Exception:
        pass

    # Test-fold HR per +1 SD + 95% CI: refit Cox on TEST. Used by the
    # test-split master leaderboard.
    for k in ("hr_per_sd_test","hr_lo_test","hr_hi_test","hr_p_test"):
        out[k] = float("nan")
    try:
        if len(test_in) >= 8 and test_in["e"].sum() >= 3:
            cph_t = CoxPHFitter(penalizer=1e-3)
            cph_t.fit(test_in, duration_col="t", event_col="e")
            smt = cph_t.summary
            out["hr_per_sd_test"] = float(smt.loc["PRS_z","exp(coef)"])
            out["hr_lo_test"]    = float(smt.loc["PRS_z","exp(coef) lower 95%"])
            out["hr_hi_test"]    = float(smt.loc["PRS_z","exp(coef) upper 95%"])
            out["hr_p_test"]     = float(smt.loc["PRS_z","p"])
    except Exception:
        pass

    def _cidx(p):
        if len(p) < 8 or p["e"].sum() < 2: return float("nan")
        pi = -cph.predict_partial_hazard(p[Xcols]).values
        try:    return float(concordance_index(p["t"], pi, p["e"]))
        except: return float("nan")
    out["val_cindex"] = _cidx(val_in)
    out["test_cindex"] = _cidx(test_in)

    # Per-patient prediction capture (multimodal alignment).
    # outcome_proba = predict_partial_hazard (relative risk score), NOT a
    # probability. Carries t, e in extra_cols for downstream survival fusion.
    if src in PER_PATIENT_SOURCES:
        try:
            train_ph = cph.predict_partial_hazard(tr_in[Xcols]).values.astype(float)
            val_ph   = (cph.predict_partial_hazard(val_in[Xcols]).values.astype(float)
                        if len(val_in) else np.array([], dtype=float))
            test_ph  = (cph.predict_partial_hazard(test_in[Xcols]).values.astype(float)
                        if len(test_in) else np.array([], dtype=float))
            cap_rows = capture_predictions(
                fits=fits, tr_in=tr_in, val_in=val_in, test_in=test_in,
                train_p=train_ph, val_p=val_ph, test_p=test_ph,
                source=src, seed=seed, covar_mode=mode,
                beta_source=PER_PATIENT_BETA_SOURCE[0],
                ld_config=PER_PATIENT_LD_CONFIG[0],
                task="cox", extra_cols={"t": "t", "e": "e"})
            # In Cox, y_true is the event flag; replace y_true with e to be unambiguous.
            cap_rows["y_true"] = cap_rows["e"].astype(int)
            PER_PATIENT_ROWS.append(cap_rows)
        except Exception:
            pass

    # Time-dependent AUC + cohort-level O/E calibration at fixed horizons.
    # td-AUC: discriminates cases (event by t_hor) vs controls (still at risk
    #         after t_hor) using predicted event probability 1 - S(t_hor).
    # O/E:    observed events by t_hor / expected events (sum of 1-S(t_hor))
    #         over the cohort. 1.0 = perfectly calibrated; >1 = under-predicts.
    def _td_auc_oe(df_in: pd.DataFrame, t_hor: float) -> tuple:
        if len(df_in) < 8:
            return float("nan"), float("nan")
        try:
            surv = cph.predict_survival_function(df_in[Xcols], times=[t_hor])
            p_event = 1.0 - surv.iloc[0].values
        except Exception:
            return float("nan"), float("nan")
        cases    = ((df_in["t"] <= t_hor) & (df_in["e"] == 1)).values
        controls = (df_in["t"]  > t_hor).values
        # td-AUC
        if cases.sum() < 1 or controls.sum() < 1:
            td_auc = float("nan")
        else:
            keep = cases | controls
            y    = np.where(cases, 1, 0)
            try:    td_auc = float(roc_auc_score(y[keep], p_event[keep]))
            except: td_auc = float("nan")
        # Cohort O/E
        n_obs = float(cases.sum())
        n_exp = float(p_event.sum())
        oe    = (n_obs / n_exp) if n_exp > 1e-9 else float("nan")
        return td_auc, oe
    for split_name, df_in in [("val", val_in), ("test", test_in)]:
        for thor in TIME_HORIZONS:
            ta, oe = _td_auc_oe(df_in, float(thor))
            out[f"{split_name}_td_auc_{thor}y"] = ta
            out[f"{split_name}_oe_{thor}y"]    = oe
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
        # PRS-CS already does LD-aware shrinkage; use unpruned 115-SNP pool only.
        out_dir = OUT_BASE.parent / f"strict_qc_prs_{args.beta_source}" / "cox"
    out_dir.mkdir(parents=True, exist_ok=True)
    PER_PATIENT_BETA_SOURCE[0] = args.beta_source
    PER_PATIENT_LD_CONFIG[0]   = (args.ld_config if args.beta_source == "raw" else "na")

    print(f"Building per-source PRS table for {args.ld_config}  beta_source={args.beta_source}...")
    prs_full, snps_per_src = per_source_prs_table(ld_config=args.ld_config,
                                                     beta_source=args.beta_source)
    cov = load_subject_covariates()
    dedup_st, dedup_dosage = get_dedup_dosage_matrix(args.ld_config)
    n_snp_per_src = {s: len(st) for s, st in snps_per_src.items()}
    n_snp_per_src["prs_all_dedup_EN_dosage"] = dedup_dosage.shape[1]
    n_snp_per_src["meta_prs_EN_combined"] = dedup_dosage.shape[1]
    if args.beta_source == "prscs":
        source_list = list(PRSCS_SOURCES) + ["covariates_only"]
    else:
        source_list = ALL_SOURCES + ["prs_all_dedup",
                                       "prs_all_dedup_ivw",
                                       "prs_all_dedup_filtered",
                                       "prs_all_dedup_EN_dosage",
                                       "meta_prs_EN_combined",
                                       "meta_prs_EN_filtered",
                                       "covariates_only"]
    n_snp_per_src["covariates_only"] = 0
    if "prs_all_dedup_ivw" in snps_per_src:
        n_snp_per_src["prs_all_dedup_ivw"] = len(snps_per_src["prs_all_dedup_ivw"])
    if "prs_all_dedup_filtered" in snps_per_src:
        n_snp_per_src["prs_all_dedup_filtered"] = len(snps_per_src["prs_all_dedup_filtered"])
    n_snp_per_src["meta_prs_EN_filtered"] = dedup_dosage.shape[1]
    sources = [args.source] if args.source else [
        s for s in source_list
        if s == "covariates_only" or n_snp_per_src.get(s, 0) > 0]
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

    td_oe_cols = {f"val_td_auc_{t}y": ("val_td_auc_" + str(t) + "y", "mean")
                    for t in TIME_HORIZONS}
    td_oe_cols.update({f"val_oe_{t}y": ("val_oe_" + str(t) + "y", "mean")
                        for t in TIME_HORIZONS})
    td_oe_cols.update({f"test_td_auc_{t}y": ("test_td_auc_" + str(t) + "y", "mean")
                        for t in TIME_HORIZONS})
    td_oe_cols.update({f"test_oe_{t}y": ("test_oe_" + str(t) + "y", "mean")
                        for t in TIME_HORIZONS})
    td_oe_cols_named = {f"{k}_mean": v for k, v in td_oe_cols.items()}
    g = df.groupby(["source","covar_mode"]).agg(
        n_snps=("n_snps_used","first"),
        n_seeds=("seed","count"),
        val_cindex_mean=("val_cindex","mean"), val_cindex_std=("val_cindex","std"),
        hr_per_sd_mean=("hr_per_sd","mean"), hr_per_sd_std=("hr_per_sd","std"),
        hr_lo_mean=("hr_lo","mean"), hr_hi_mean=("hr_hi","mean"),
        hr_p_min=("hr_p","min"),
        hr_per_sd_val_mean=("hr_per_sd_val","mean"),
        hr_per_sd_val_std =("hr_per_sd_val","std"),
        hr_lo_val_mean=("hr_lo_val","mean"), hr_hi_val_mean=("hr_hi_val","mean"),
        hr_p_val_min=("hr_p_val","min"),
        test_cindex_mean=("test_cindex","mean"),
        test_cindex_std=("test_cindex","std"),
        hr_per_sd_test_mean=("hr_per_sd_test","mean"),
        hr_per_sd_test_std =("hr_per_sd_test","std"),
        hr_lo_test_mean=("hr_lo_test","mean"), hr_hi_test_mean=("hr_hi_test","mean"),
        hr_p_test_min=("hr_p_test","min"),
        **td_oe_cols_named,
    ).reset_index().sort_values(["covar_mode","val_cindex_mean"], ascending=[True, False])
    lb_p = out_dir / "leaderboard.csv"
    g.to_csv(lb_p, index=False)
    print(f"wrote leaderboard: {lb_p}")

    # PNG — concise; time-dep AUC + O/E at 3/5/10y to show calibration & discrim.
    def _fmt(x, fmt="{:.3f}"):
        return "" if pd.isna(x) else fmt.format(x)
    show = pd.DataFrame({
        "source":     g["source"],
        "covar_mode": g["covar_mode"],
        "n_snps":     g["n_snps"].astype(int),
        "val_cindex": [f"{m:.3f}+/-{s:.3f}" for m, s in zip(g["val_cindex_mean"], g["val_cindex_std"])],
        "HR/+1SD (val)":
                      [(f"{m:.2f}" if (pd.isna(lo) or pd.isna(hi))
                        else f"{m:.2f} [{lo:.2f}-{hi:.2f}]")
                       for m, lo, hi in zip(g["hr_per_sd_val_mean"],
                                              g["hr_lo_val_mean"],
                                              g["hr_hi_val_mean"])],
        "td_AUC 3y": [_fmt(x) for x in g["val_td_auc_3y_mean"]],
        "td_AUC 5y": [_fmt(x) for x in g["val_td_auc_5y_mean"]],
        "td_AUC 10y":[_fmt(x) for x in g["val_td_auc_10y_mean"]],
        "O/E 3y":    [_fmt(x, "{:.2f}") for x in g["val_oe_3y_mean"]],
        "O/E 5y":    [_fmt(x, "{:.2f}") for x in g["val_oe_5y_mean"]],
        "O/E 10y":   [_fmt(x, "{:.2f}") for x in g["val_oe_10y_mean"]],
        "min_p":     [f"{p:.1e}" for p in g["hr_p_min"]],
    })
    fig_h = max(2.5, 0.30 * len(show) + 1.3)
    fig, ax = plt.subplots(figsize=(18, fig_h)); ax.axis("off")
    ax.set_title(f"Cox survival — strict-QC + LD pruning [{args.ld_config}]   "
                 f"(VAL set; mean +/- std over 3 seeds; td-AUC + O/E at 3/5/10y)",
                 fontsize=11, fontweight="bold", pad=12, loc="left")
    col_widths = [0.18, 0.13, 0.05, 0.12, 0.14,
                   0.06, 0.06, 0.06, 0.05, 0.05, 0.05, 0.07]
    tbl = ax.table(cellText=show.values.tolist(), colLabels=show.columns.tolist(),
                    loc="center", cellLoc="center", colWidths=col_widths)
    tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1, 1.2)
    for j in range(len(show.columns)):
        c = tbl[(0, j)]; c.set_facecolor("#222"); c.set_text_props(color="white", weight="bold")
    # Bold C-index > 0.5, and HR/+1SD (val) when its 95% CI excludes 1.0.
    g_lo = g["hr_lo_val_mean"].tolist()
    g_hi = g["hr_hi_val_mean"].tolist()
    for i, row in enumerate(show.itertuples(index=False), start=1):
        mean_c = float(row.val_cindex.split("+/-")[0])
        if mean_c > 0.5:
            tbl[(i, 3)].set_text_props(weight="bold")
        lo, hi = g_lo[i - 1], g_hi[i - 1]
        if not (pd.isna(lo) or pd.isna(hi)) and (lo > 1.0 or hi < 1.0):
            tbl[(i, 4)].set_text_props(weight="bold")
    png = out_dir / "leaderboard.png"
    fig.savefig(png, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"wrote PNG: {png}")

    # Per-patient prediction parquets (multimodal-alignment anchors).
    # outcome_proba is the partial hazard (not a probability); t, e carried
    # in extra columns so downstream survival fusion has the duration + event.
    if PER_PATIENT_ROWS:
        pp_dir = out_dir / "per_patient"
        anchor_dir = (OUT_BASE / "multimodal_anchors" / "prs" / "cox"
                       if args.beta_source == "raw"
                       else OUT_BASE.parent / f"strict_qc_prs_{args.beta_source}"
                            / "multimodal_anchors" / "prs" / "cox")
        all_rows = pd.concat(PER_PATIENT_ROWS, ignore_index=True)
        for (src, mode, seed), grp in all_rows.groupby(
                ["model", "covar_mode", "seed"], sort=False):
            save_per_patient_parquet(grp, pp_dir,
                                       source=src, seed=int(seed),
                                       covar_mode=mode)
        for (src, seed), grp in all_rows.groupby(
                ["model", "seed"], sort=False):
            concat_and_save_anchor([grp], anchor_dir,
                                     source=src, seed=int(seed))
        n_pids = all_rows["Patient_ID"].nunique()
        print(f"wrote per-patient parquets: {pp_dir} "
              f"({all_rows.shape[0]:,} rows, {n_pids} unique PTIDs)")
        print(f"wrote multimodal anchors:   {anchor_dir}")

    # EN coefficient bar chart for the Cox EN models.
    dos_df  = pd.DataFrame(EN_COEFS["dosage"]) if EN_COEFS["dosage"] else None
    meta_df = pd.DataFrame(EN_COEFS["meta"])   if EN_COEFS["meta"]   else None
    if dos_df is not None or meta_df is not None:
        en_png = out_dir / "en_weights.png"
        render_en_weights_png(
            dos_df, meta_df, en_png,
            title=(f"Cox survival — EN coefficients (TRAIN-fold fit)   "
                    f"LD={args.ld_config}  beta_source={args.beta_source}\n"
                    f"top: per-SNP dosage-EN log-hazard per +1 A1   ·   "
                    f"bottom: meta-EN log-hazard per +1 source-PRS"))
        print(f"wrote EN weights PNG: {en_png}")

    # Test-fold EN coefficient refit — mirror of train-fold for the test
    # leaderboard. Fits EN on test data (size ~30 with ~10-15 events).
    EN_COEFS_TEST = {"dosage": {}, "meta": {}}
    for seed in SEEDS:
        parts_seed = {sp: _load_time_event(seed, sp) for sp in ("train","val","test")}
        parts_fit = {"train": parts_seed["test"],
                      "val":   parts_seed["val"],
                      "test":  parts_seed["train"]}
        cph_d, _ = _fit_en_cox_dosage_get_prs(dedup_dosage, parts_fit, seed)
        if cph_d is not None:
            ps = cph_d.params_
            EN_COEFS_TEST["dosage"][seed] = ps[[c for c in ps.index
                                                  if c not in ("t","e")]]
        meta_source_cols = [
            c for c in prs_full.columns if c.startswith("PRS_")
            and c not in ("PRS_prs_all_dedup", "PRS_prs_all_dedup_ivw",
                            "PRS_Kosteridis_MTAG_AD",
                            "PRS_Kosteridis_shared_AD_CV",
                            "PRS_Kosteridis_novel_AD")]
        sp_df = prs_full[["PTID"] + meta_source_cols].set_index("PTID")
        sp_df.index = sp_df.index.astype(str)
        cph_m, _ = _fit_en_cox_meta_prs(sp_df, parts_fit, seed)
        if cph_m is not None:
            ps = cph_m.params_
            ps = ps[[c for c in ps.index if c not in ("t","e")]]
            ps.index = [c.replace("PRS_", "", 1) for c in ps.index]
            EN_COEFS_TEST["meta"][seed] = ps
    dos_df_t  = pd.DataFrame(EN_COEFS_TEST["dosage"]) if EN_COEFS_TEST["dosage"] else None
    meta_df_t = pd.DataFrame(EN_COEFS_TEST["meta"])   if EN_COEFS_TEST["meta"]   else None
    if dos_df_t is not None or meta_df_t is not None:
        en_png_t = out_dir / "en_weights_test.png"
        render_en_weights_png(
            dos_df_t, meta_df_t, en_png_t,
            title=(f"Cox survival — EN coefficients (TEST-fold refit)   "
                    f"LD={args.ld_config}  beta_source={args.beta_source}\n"
                    f"top: per-SNP dosage-EN log-hazard per +1 A1   ·   "
                    f"bottom: meta-EN log-hazard per +1 source-PRS\n"
                    f"caveat: small test sets (~30) -> high coef variance"))
        print(f"wrote EN weights PNG (TEST): {en_png_t}")


if __name__ == "__main__":
    main()
