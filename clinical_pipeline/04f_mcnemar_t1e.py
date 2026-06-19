r"""
04f_mcnemar_t1e.py   (env: snp — pandas/sklearn/xgboost/scipy/pyarrow)
=====================================================================
McNemar test (per-patient paired) of the best clinical LLM vs the best tabular baseline on the T1e
(sCN vs pCN) held-out TEST set, pooled across seeds 0/1/2 — replacing the per-seed 2-sided t-test that
currently annotates table_val_t1e.tex.

Why McNemar: the LLM and the tabular baseline are evaluated on the SAME held-out patients, so the right
test of "is model A better than model B" is the paired per-patient discordance (McNemar exact binomial),
not a t-test on 3 seed-level means (df=2).

The mcnemar() function is the SAME one used in
integration/cross_modal_attention/04_significance_unimodal.py (reused verbatim, not reinvented). Only the
T1e data loading (LLM preds from the encoder parquet; tabular XGBoost refit) is new here.

Best LLM (val bACC, full fine-tune): ModernBERT-base (0.87); ModernBERT-large (0.82) also reported.
Best baseline: XGBoost (the T1e tabular leaderboard winner).

Run:  conda run -n snp python clinical_pipeline/04f_mcnemar_t1e.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.pipeline import Pipeline
from sklearn.metrics import balanced_accuracy_score
from xgboost import XGBClassifier

ENC_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\encoder_outputs_no_cdr_post_exclusion")
TAB_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion"
               r"\tabular\baseline")
SEEDS = (0, 1, 2)
TASK_ENC = "T1e_scn_pcn"
LABEL_COL = "conversion_group"
LABEL_MAP = {"sCN": 0, "pCN_to_MCI": 1, "pCN_to_AD": 1}


# ── reused verbatim from integration/cross_modal_attention/04_significance_unimodal.py ──
def mcnemar(correct_a, correct_b):
    """Exact McNemar on paired per-sample correctness. Returns (b, c, p)."""
    b = int(np.sum(correct_a & ~correct_b))   # A right, B wrong
    c = int(np.sum(~correct_a & correct_b))   # A wrong, B right
    p = stats.binomtest(min(b, c), b + c, 0.5).pvalue if (b + c) > 0 else 1.0
    return b, c, p


# ── tabular feature prep (replicated from 02h / 03_baseline_shap_pdp, no shap dep) ──
DROP_ALWAYS = [
    "Patient_ID", "VISCODE_long", "VISCODE_2", "Date", "Generated_Text",
    "Label_bl_multi", "Label_visit_diag", "Label_visit_granular",
    "Label_1y", "Label_2y", "Label_3y", "Label_4y", "Label_5y",
    "Label_6y", "Label_7y", "Label_8y", "Label_9y", "Label_10y",
    "LogicalMemory_I_Total", "LogicalMemory_II_Total", "LogicalMemory_II_Cued",
    "Plasma_Abeta42", "Plasma_Abeta40", "Plasma_Abeta42_40",
    "Plasma_pTau217", "Plasma_pTau181", "Plasma_NfL", "Plasma_GFAP",
    "CDR_Global", "CDR_SumBoxes",
    "bl_dx", "last_dx", "conversion_group", "AD_bl", "AD_final",
    "pMCI", "sMCI", "pCN_to_AD", "pCN_to_MCI", "sCN", "CN_to_AD", "CN_to_MCI",
    "Excluded", "FU_years", "years_to_MCI", "years_to_AD",
]
CAT = ["Sex", "Handedness", "Marital_Status", "Retired", "Language", "Ethnicity",
       "APOE4_Status", "APOE_Alleles"]


def prepare_features(df, fit_encoders=None, feature_cols=None):
    df = df.copy().drop(columns=[c for c in DROP_ALWAYS if c in df.columns])
    enc = fit_encoders or {}
    for c in CAT:
        if c not in df.columns:
            continue
        df[c] = df[c].astype(str).str.strip().replace({"NaN": np.nan, "nan": np.nan})
        if fit_encoders is None:
            le = LabelEncoder(); le.fit(df[c].dropna().astype(str)); enc[c] = le
        le = enc[c]
        df[c] = df[c].apply(lambda v: le.transform([str(v)])[0]
                            if pd.notna(v) and str(v) in le.classes_ else np.nan)
    df = df.apply(pd.to_numeric, errors="coerce")
    if feature_cols is None:
        feature_cols = list(df.columns)
    for c in feature_cols:
        if c not in df.columns:
            df[c] = np.nan
    return df[feature_cols], enc, feature_cols


def make_xgb():
    return Pipeline([("imp", SimpleImputer(strategy="median")),
                     ("scl", StandardScaler()),
                     ("clf", XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.05,
                                           use_label_encoder=False, eval_metric="logloss",
                                           random_state=42, verbosity=0, scale_pos_weight=1))])


def baseline_preds(seed, split):
    """XGBoost refit on T1e train; predictions on `split` -> [Patient_ID, y, pred]."""
    def prep(sp):
        d = pd.read_csv(TAB_DIR / f"seed_{seed}" / f"{sp}.csv", low_memory=False)
        d[LABEL_COL] = d[LABEL_COL].map(LABEL_MAP)
        return d.dropna(subset=[LABEL_COL])
    tr, ev = prep("train"), prep(split)
    Xtr, enc, feat = prepare_features(tr)
    Xev, _, _ = prepare_features(ev, fit_encoders=enc, feature_cols=feat)
    pipe = make_xgb().fit(Xtr, tr[LABEL_COL].astype(int).to_numpy())
    return pd.DataFrame({"Patient_ID": ev["Patient_ID"].astype(str).to_numpy(),
                         "y": ev[LABEL_COL].astype(int).to_numpy(),
                         "pred": pipe.predict(Xev).astype(int)})


def llm_preds(model_slug, seed, split):
    """LLM predictions on the T1e cohort for `split` -> [Patient_ID, y, pred]."""
    p = ENC_DIR / model_slug / TASK_ENC / f"seed_{seed}" / "full_ft" / "embeddings" / "embeddings.parquet"
    d = pd.read_parquet(p, columns=["Patient_ID", "split", "in_task_cohort", "label", "pred"])
    d = d[(d["split"] == split) & (d["in_task_cohort"]) & d["label"].notna()]
    return pd.DataFrame({"Patient_ID": d["Patient_ID"].astype(str).to_numpy(),
                         "y": d["label"].astype(int).to_numpy(),
                         "pred": d["pred"].astype(int).to_numpy()})


def run(model_slug, split):
    ca_all, cb_all, baccs_a, baccs_b = [], [], [], []
    for s in SEEDS:
        L, B = llm_preds(model_slug, s, split), baseline_preds(s, split)
        m = L.merge(B, on="Patient_ID", suffixes=("_llm", "_bl"))
        assert (m["y_llm"] == m["y_bl"]).all(), f"label mismatch seed {s}"
        ca = (m["pred_llm"] == m["y_llm"]).to_numpy(); cb = (m["pred_bl"] == m["y_bl"]).to_numpy()
        ca_all.append(ca); cb_all.append(cb)
        baccs_a.append(balanced_accuracy_score(m["y_llm"], m["pred_llm"]))
        baccs_b.append(balanced_accuracy_score(m["y_bl"], m["pred_bl"]))
        print(f"    seed {s}: n={len(m)}  {model_slug} acc={ca.mean():.3f} bACC={baccs_a[-1]:.3f} | "
              f"XGB acc={cb.mean():.3f} bACC={baccs_b[-1]:.3f}")
    ca, cb = np.concatenate(ca_all), np.concatenate(cb_all)
    b, c, p = mcnemar(ca, cb)                              # b = LLM-right/XGB-wrong, c = LLM-wrong/XGB-right
    star = "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "ns"
    print(f"  [{split}] {model_slug} vs XGBoost  (pooled n={len(ca)}):  "
          f"{model_slug} acc {ca.mean():.3f}/bACC {np.mean(baccs_a):.3f}  vs  "
          f"XGB acc {cb.mean():.3f}/bACC {np.mean(baccs_b):.3f}")
    print(f"    discordant b(LLM>XGB)={b}, c(XGB>LLM)={c}  ->  McNemar exact p={p:.4g}  [{star}]\n")
    return dict(split=split, model=model_slug, n=int(len(ca)), b_llm_right=b, c_bl_right=c,
                mcnemar_p=float(p), star=star, llm_bacc=float(np.mean(baccs_a)),
                xgb_bacc=float(np.mean(baccs_b)))


def main():
    print("McNemar (exact, paired per-patient, pooled seeds 0/1/2) — T1e best LLM vs XGBoost:\n")
    res = pd.DataFrame([run(m, sp) for sp in ("val", "test")
                        for m in ("ModernBERT-base", "ModernBERT-large")])
    out = Path(__file__).resolve().parent / "outputs" / "encoder_only_post_exclusions" / "mcnemar_t1e.csv"
    res.to_csv(out, index=False)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
