"""
02h_baseline_post_exclusion.py
==============================
Tabular baseline models (LogReg / SVM / XGBoost) on the NO-CDR stratified
*POST-EXCLUSION* data, for the expanded results table (8 tasks).

Same training core as 02d (Part 1) but:
  - reads no_cdr_stratified_post_exclusion/tabular/baseline
  - DROP_ALWAYS extended with the post-exclusion OUTCOME columns
    (conversion_group, pMCI, sMCI, sCN, pCN_*, AD_final, FU_years, years_to_* ...)
    — these are present in the post-exclusion splits and would otherwise LEAK the
    labels into the features (catastrophic for T1d/T1e/conversion tasks).
  - 8 tasks: adds T1d (pMCI vs sMCI), T1e (sCN vs pCN -> MCI/AD), T3d (conv 7y).
  - writes baseline_model_comparison.csv + results_per_seed.csv into
    outputs/encoder_only_post_exclusions/ (same folder the table reads from).

Task labels for T1d/T1e are derived from `conversion_group` via label_map; rows in
other groups map to NaN and are dropped (so no separate cohort filter is needed).

Run (local Windows, env: clinical or snp — needs sklearn + xgboost):
  python clinical_pipeline/02h_baseline_post_exclusion.py
"""

import warnings
import numpy as np
import pandas as pd
from pathlib import Path

from sklearn.pipeline      import Pipeline
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.impute        import SimpleImputer
from sklearn.linear_model  import LogisticRegression
from sklearn.svm           import SVC
from sklearn.metrics       import (
    accuracy_score, balanced_accuracy_score, roc_auc_score,
    f1_score, precision_score, recall_score,
)
from xgboost import XGBClassifier

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────────
SPLIT_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical"
                 r"\no_cdr_stratified_post_exclusion\tabular\baseline")
REPO_ROOT = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
OUT_DIR   = REPO_ROOT / "clinical_pipeline" / "outputs" / "encoder_only_post_exclusions"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SEEDS  = [0, 1, 2]
MODELS = ["LogReg", "SVM", "XGBoost"]

# ── Feature columns ───────────────────────────────────────────────────────────
DROP_ALWAYS = [
    "Patient_ID", "VISCODE_long", "VISCODE_2", "Date", "Generated_Text",
    "Label_bl_multi", "Label_visit_diag", "Label_visit_granular",
    "Label_1y", "Label_2y", "Label_3y", "Label_4y", "Label_5y",
    "Label_6y", "Label_7y", "Label_8y", "Label_9y", "Label_10y",
    "LogicalMemory_I_Total", "LogicalMemory_II_Total", "LogicalMemory_II_Cued",
    "Plasma_Abeta42", "Plasma_Abeta40", "Plasma_Abeta42_40",
    "Plasma_pTau217", "Plasma_pTau181", "Plasma_NfL", "Plasma_GFAP",
    "CDR_Global", "CDR_SumBoxes",
    # ── POST-EXCLUSION OUTCOME COLUMNS — drop to prevent target leakage ──
    "bl_dx", "last_dx", "conversion_group", "AD_bl", "AD_final",
    "pMCI", "sMCI", "pCN_to_AD", "pCN_to_MCI", "sCN", "CN_to_AD", "CN_to_MCI",
    "Excluded", "FU_years", "years_to_MCI", "years_to_AD",
]

CATEGORICAL_ENCODE = [
    "Sex", "Handedness", "Marital_Status", "Retired", "Language", "Ethnicity",
    "APOE4_Status", "APOE_Alleles",
]


def prepare_features(df, fit_encoders=None, feature_cols=None):
    df = df.copy()
    df = df.drop(columns=[c for c in DROP_ALWAYS if c in df.columns])
    encoders = fit_encoders or {}
    for col in CATEGORICAL_ENCODE:
        if col not in df.columns:
            continue
        df[col] = df[col].astype(str).str.strip().replace({"NaN": np.nan, "nan": np.nan})
        if fit_encoders is None:
            le = LabelEncoder()
            le.fit(df[col].dropna().astype(str))
            encoders[col] = le
        le = encoders[col]
        df[col] = df[col].apply(
            lambda v: le.transform([str(v)])[0]
            if pd.notna(v) and str(v) in le.classes_ else np.nan)
    df = df.apply(pd.to_numeric, errors="coerce")
    if feature_cols is None:
        feature_cols = list(df.columns)
    for col in feature_cols:
        if col not in df.columns:
            df[col] = np.nan
    df = df[feature_cols]
    return df, encoders, feature_cols


def make_pipeline(model_name):
    imp = SimpleImputer(strategy="median")
    scl = StandardScaler()
    if model_name == "LogReg":
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", random_state=42)
    elif model_name == "SVM":
        clf = SVC(probability=True, class_weight="balanced", random_state=42)
    elif model_name == "XGBoost":
        clf = XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.05,
                            use_label_encoder=False, eval_metric="logloss",
                            random_state=42, verbosity=0, scale_pos_weight=1)
    else:
        raise ValueError(model_name)
    return Pipeline([("imp", imp), ("scl", scl), ("clf", clf)])


def evaluate(pipe, X_test, y_test, task_type="binary"):
    y_pred = pipe.predict(X_test)
    metrics = {
        "Accuracy":    round(accuracy_score(y_test, y_pred), 4),
        "BalancedAcc": round(balanced_accuracy_score(y_test, y_pred), 4),
    }
    if task_type == "binary":
        y_prob = pipe.predict_proba(X_test)[:, 1]
        metrics["AUC-ROC"]     = round(roc_auc_score(y_test, y_prob), 4)
        metrics["Precision"]   = round(precision_score(y_test, y_pred, zero_division=0), 4)
        metrics["Recall"]      = round(recall_score(y_test, y_pred, zero_division=0), 4)
        # Sensitivity = recall of positive class; specificity = recall of negative class.
        metrics["Sensitivity"] = round(recall_score(y_test, y_pred, pos_label=1, zero_division=0), 4)
        metrics["Specificity"] = round(recall_score(y_test, y_pred, pos_label=0, zero_division=0), 4)
        metrics["F1"]          = round(f1_score(y_test, y_pred, zero_division=0), 4)
    else:
        y_prob = pipe.predict_proba(X_test)
        try:
            metrics["AUC-ROC (OvR)"] = round(
                roc_auc_score(y_test, y_prob, multi_class="ovr", average="macro"), 4)
        except Exception:
            metrics["AUC-ROC (OvR)"] = np.nan
        metrics["MacroPrecision"] = round(precision_score(y_test, y_pred, average="macro", zero_division=0), 4)
        metrics["MacroRecall"]    = round(recall_score(y_test, y_pred, average="macro", zero_division=0), 4)
        metrics["MacroF1"]        = round(f1_score(y_test, y_pred, average="macro", zero_division=0), 4)
    return metrics


# ── Task definitions (8) ───────────────────────────────────────────────────────
# NOTE: descriptions must not share substrings used as the table's `bl_task` matcher
# (e.g. only T1 contains "Binary"; T1d/T1e drop that prefix).
TASKS = {
    "T1_CN_vs_MCIAD": {
        "label_col": "Label_bl_multi", "task_type": "binary",
        "label_map": {"CN": 0, "MCI": 1, "AD": 1}, "filter": None,
        "description": "Binary: CN vs MCI+AD",
    },
    "T1b_CNMCI_AD": {
        "label_col": "Label_bl_multi", "task_type": "binary",
        "label_map": {"CN": 0, "MCI": 0, "AD": 1}, "filter": None,
        "description": "CN+MCI vs AD",
    },
    "T2_Multiclass": {
        "label_col": "Label_bl_multi", "task_type": "multiclass",
        "label_map": {"CN": 0, "MCI": 1, "AD": 2}, "filter": None,
        "description": "Multi-class: CN / MCI / AD",
    },
    "T1c_sCN_prog": {
        "label_col": "conversion_group", "task_type": "binary",
        "label_map": {"sCN": 0, "pCN_to_MCI": 1, "pCN_to_AD": 1, "pMCI": 1, "AD_bl": 1},
        "filter": None,
        "description": "sCN vs progressors+AD (excl. sMCI)",
    },
    "T1d_pMCI_sMCI": {
        "label_col": "conversion_group", "task_type": "binary",
        "label_map": {"sMCI": 0, "pMCI": 1}, "filter": None,
        "description": "pMCI vs sMCI (baseline MCI)",
    },
    "T1e_sCN_pCN": {
        "label_col": "conversion_group", "task_type": "binary",
        "label_map": {"sCN": 0, "pCN_to_MCI": 1, "pCN_to_AD": 1}, "filter": None,
        "description": "sCN vs pCN (baseline CN, to MCI/AD)",
    },
    "T3a_Conv3y": {
        "label_col": "Label_3y", "task_type": "binary", "label_map": None,
        "filter": "non_ad", "description": "Conversion to AD within 3 years",
    },
    "T3b_Conv5y": {
        "label_col": "Label_5y", "task_type": "binary", "label_map": None,
        "filter": "non_ad", "description": "Conversion to AD within 5 years",
    },
    "T3d_Conv7y": {
        "label_col": "Label_7y", "task_type": "binary", "label_map": None,
        "filter": "non_ad", "description": "Conversion to AD within 7 years",
    },
    "T3c_Conv10y": {
        "label_col": "Label_10y", "task_type": "binary", "label_map": None,
        "filter": "non_ad", "description": "Conversion to AD within 10 years",
    },
}


# ── Sanity check: confirm leakage columns will be dropped ──────────────────────
_t0 = pd.read_csv(SPLIT_DIR / "seed_0" / "test.csv", low_memory=False)
_leak = [c for c in ("conversion_group", "pMCI", "sMCI", "years_to_AD", "AD_final")
         if c in _t0.columns]
print(f"Outcome columns present in splits (will be dropped from features): {_leak}")
print(f"Seed_0 test rows: {len(_t0)}")


# ═══════════════════════════════════════════════════════════════════════════════
# MODEL TRAINING
# ═══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*65}\n  BASELINE MODELS — NO CDR, POST-EXCLUSION\n{'='*65}")

all_seed_results = []

for seed in SEEDS:
    seed_dir = SPLIT_DIR / f"seed_{seed}"
    try:
        train_s = pd.read_csv(seed_dir / "train.csv", low_memory=False)
        val_s   = pd.read_csv(seed_dir / "val.csv",   low_memory=False)
        test_s  = pd.read_csv(seed_dir / "test.csv",  low_memory=False)
    except FileNotFoundError as e:
        print(f"  [seed {seed}] Missing: {e}"); continue

    trainval_s = pd.concat([train_s, val_s], ignore_index=True)
    print(f"\n{'#'*65}\n  SEED {seed}  —  Train+Val: {len(trainval_s)}  Test: {len(test_s)}\n{'#'*65}")

    for task_id, task_cfg in TASKS.items():
        print(f"\n  [{task_cfg['description']}]")
        lcol, lmap, filter_ = task_cfg["label_col"], task_cfg["label_map"], task_cfg["filter"]

        def prep_split(df, _lcol=lcol, _lmap=lmap, _filter=filter_):
            df = df.copy()
            if _lmap:
                df[_lcol] = df[_lcol].map(_lmap)
            else:
                df[_lcol] = pd.to_numeric(df[_lcol], errors="coerce")
            if _filter == "non_ad":
                df = df[df["Label_bl_multi"].map({"CN": True, "MCI": True, "AD": False})]
            return df.dropna(subset=[_lcol])

        tv_raw   = prep_split(trainval_s)
        test_raw = prep_split(test_s)
        if tv_raw.empty or test_raw.empty:
            print("    Insufficient data — skipping."); continue

        y_tv   = tv_raw[lcol].astype(int)
        y_test = test_raw[lcol].astype(int)
        print(f"    TV: {len(tv_raw)} | Test: {len(test_raw)} | "
              f"Classes: {dict(y_tv.value_counts().sort_index())}")

        X_tv,   enc, feat_cols = prepare_features(tv_raw)
        X_test, _,   _         = prepare_features(test_raw, fit_encoders=enc,
                                                  feature_cols=feat_cols)
        for model_name in MODELS:
            pipe = make_pipeline(model_name)
            try:
                pipe.fit(X_tv, y_tv)
                metrics = evaluate(pipe, X_test, y_test, task_type=task_cfg["task_type"])
                metrics.update({"Model": model_name, "Task": task_cfg["description"],
                                "Task_ID": task_id, "Seed": seed})
                all_seed_results.append(metrics)
                m_str = "  ".join(f"{k}: {v}" for k, v in metrics.items()
                                  if k not in ("Model", "Task", "Task_ID", "Seed"))
                print(f"    [{model_name:8s}] {m_str}")
            except Exception as e:
                print(f"    [{model_name:8s}] ERROR: {e}")

# ── Aggregate (mean ± std across seeds) ───────────────────────────────────────
if all_seed_results:
    raw_df    = pd.DataFrame(all_seed_results)
    meta_cols = ["Task", "Model", "Task_ID", "Seed"]
    metric_cols = [c for c in raw_df.columns if c not in meta_cols]

    raw_df[["Task", "Model", "Seed"] + metric_cols].to_csv(
        OUT_DIR / "results_per_seed.csv", index=False)

    agg = raw_df.groupby(["Task_ID", "Task", "Model"])[metric_cols].agg(["mean", "std"])
    agg.columns = [f"{m}_{stat}" for m, stat in agg.columns]
    agg = agg.reset_index()

    summary_rows = []
    for _, row in agg.iterrows():
        entry = {"Task": row["Task"], "Model": row["Model"]}
        for m in metric_cols:
            mu = row.get(f"{m}_mean", np.nan)
            sd = row.get(f"{m}_std",  np.nan)
            if pd.notna(mu):
                entry[m] = f"{mu:.4f} ± {sd:.4f}" if pd.notna(sd) else f"{mu:.4f}"
            else:
                entry[m] = ""
        summary_rows.append(entry)
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT_DIR / "baseline_model_comparison.csv", index=False)

    print(f"\n{'='*75}\n  SUMMARY  (mean ± std, no CDR, post-exclusion)\n{'='*75}")
    print(summary_df.to_string(index=False))
    print(f"\nSaved → {OUT_DIR / 'baseline_model_comparison.csv'}")
    print(f"Saved → {OUT_DIR / 'results_per_seed.csv'}")
else:
    print("\n[ERROR] No results produced.")
