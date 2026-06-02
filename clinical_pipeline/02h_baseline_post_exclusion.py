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

import argparse
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
OUTPUTS   = REPO_ROOT / "clinical_pipeline" / "outputs"

SEEDS  = [0, 1, 2]
MODELS = ["LogReg", "SVM", "XGBoost"]

# ── Feature-set selection (full battery vs reduced "totals/filtered") ──────────
_ap = argparse.ArgumentParser(description="Tabular baselines (no-CDR, post-exclusion).")
_ap.add_argument("--feature_set", choices=["full", "totals"], default="full",
                 help="'full' = whole battery (current, feeds the combined table); "
                      "'totals' = one summary score per neurocognitive test + "
                      "demographics/APOE/CSF (26 features). Outputs get a '_filtered' suffix.")
_args, _ = _ap.parse_known_args()
FEATURE_SET = _args.feature_set
SUFFIX = "" if FEATURE_SET == "full" else "_filtered"

# Reduced "totals/filtered" set (26 features): demographics + APOE + CSF + one summary score
# per neurocognitive test. AVLT/Trails have no single native total → user-selected primary
# components (AVLT T6 / List B / 30-min delayed recall; Trails A/B completion time). EXCLUDES
# CDR (no-cdr), Logical Memory (not at baseline), Plasma (sparse), GDS/Hachinski, and ALL
# per-trial / per-item / error-count breakdown columns.
TOTALS_INCLUDE = [
    "Sex", "YOB", "Handedness", "Marital_Status", "Education_Years",
    "Retired", "Residence", "Language", "Ethnicity",
    "APOE4_Status", "APOE_Alleles", "APOE4_Dosage",
    "CSF_Abeta_bl", "CSF_Tau_bl", "CSF_pTau_bl",
    "ADAS_Total13", "MMSE_Total", "MoCA_Total", "CDT_Total",
    "CategoryFluency_Animals_Total", "FAQ_Total",
    "AVLT_T6", "AVLT_ListB", "AVLT_Delay30",
    "TrailsA_Time", "TrailsB_Time",
]
FIT_FEATURE_COLS = TOTALS_INCLUDE if FEATURE_SET == "totals" else None

# full-feature run feeds the combined encoder table; reduced run gets its own table dir.
OUT_DIR = (OUTPUTS / "encoder_only_post_exclusions" if FEATURE_SET == "full"
           else OUTPUTS / "baseline_no_cdr_post_exclusion")
OUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"Feature set: {FEATURE_SET}"
      + (f"  ({len(TOTALS_INCLUDE)} features) → {OUT_DIR.name}/  (suffix '{SUFFIX}')"
         if FEATURE_SET == "totals" else f" → {OUT_DIR.name}/"))
RESOLVED_FEATURES = None   # captured from the first task's fit (same set across tasks)

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
if FEATURE_SET == "totals":
    _missing = [c for c in TOTALS_INCLUDE if c not in _t0.columns]
    if _missing:
        print(f"  [WARN] requested 'totals' features missing from the splits "
              f"(will be all-NaN → median-imputed): {_missing}")


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

    # Fit on TRAIN only; evaluate VAL and TEST as genuine held-out sets so a real
    # validation metric can be reported. Imputer/scaler live inside the pipeline and are
    # therefore fit on train only — no val/test leakage. (Previously train+val were merged
    # for fitting and only test was scored, leaving no held-out validation metric.)
    print(f"\n{'#'*65}\n  SEED {seed}  —  Train: {len(train_s)}  Val: {len(val_s)}  "
          f"Test: {len(test_s)}\n{'#'*65}")

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

        train_raw = prep_split(train_s)
        val_raw   = prep_split(val_s)
        test_raw  = prep_split(test_s)
        if train_raw.empty or test_raw.empty:
            print("    Insufficient data — skipping."); continue

        y_train = train_raw[lcol].astype(int)
        print(f"    Train: {len(train_raw)} | Val: {len(val_raw)} | Test: {len(test_raw)} | "
              f"Train classes: {dict(y_train.value_counts().sort_index())}")

        X_train, enc, feat_cols = prepare_features(train_raw, feature_cols=FIT_FEATURE_COLS)
        if RESOLVED_FEATURES is None:
            RESOLVED_FEATURES = list(feat_cols)
            print(f"    Feature columns ({len(feat_cols)}): {feat_cols}")

        # Build held-out design matrices once (independent of model).
        eval_splits = []
        if not val_raw.empty:
            X_val, _, _ = prepare_features(val_raw, fit_encoders=enc, feature_cols=feat_cols)
            eval_splits.append(("val", X_val, val_raw[lcol].astype(int)))
        X_test, _, _ = prepare_features(test_raw, fit_encoders=enc, feature_cols=feat_cols)
        eval_splits.append(("test", X_test, test_raw[lcol].astype(int)))

        for model_name in MODELS:
            pipe = make_pipeline(model_name)
            try:
                pipe.fit(X_train, y_train)
            except Exception as e:
                print(f"    [{model_name:8s}] FIT ERROR: {e}"); continue
            for split_name, X_split, y_split in eval_splits:
                try:
                    metrics = evaluate(pipe, X_split, y_split, task_type=task_cfg["task_type"])
                    metrics.update({"Model": model_name, "Task": task_cfg["description"],
                                    "Task_ID": task_id, "Seed": seed, "Split": split_name})
                    all_seed_results.append(metrics)
                    m_str = "  ".join(f"{k}: {v}" for k, v in metrics.items()
                                      if k not in ("Model", "Task", "Task_ID", "Seed", "Split"))
                    print(f"    [{model_name:8s}|{split_name:4s}] {m_str}")
                except Exception as e:
                    print(f"    [{model_name:8s}|{split_name:4s}] EVAL ERROR: {e}")

# ── Aggregate (mean ± std across seeds, per split) ────────────────────────────
if all_seed_results:
    raw_df    = pd.DataFrame(all_seed_results)
    meta_cols = ["Task", "Model", "Task_ID", "Seed", "Split"]
    metric_cols = [c for c in raw_df.columns if c not in meta_cols]

    raw_df[["Task", "Model", "Seed", "Split"] + metric_cols].to_csv(
        OUT_DIR / f"results_per_seed{SUFFIX}.csv", index=False)

    agg = raw_df.groupby(["Task_ID", "Task", "Model", "Split"])[metric_cols].agg(["mean", "std"])
    agg.columns = [f"{m}_{stat}" for m, stat in agg.columns]
    agg = agg.reset_index()

    summary_rows = []
    for _, row in agg.iterrows():
        entry = {"Task": row["Task"], "Model": row["Model"], "Split": row["Split"]}
        for m in metric_cols:
            mu = row.get(f"{m}_mean", np.nan)
            sd = row.get(f"{m}_std",  np.nan)
            if pd.notna(mu):
                entry[m] = f"{mu:.4f} ± {sd:.4f}" if pd.notna(sd) else f"{mu:.4f}"
            else:
                entry[m] = ""
        summary_rows.append(entry)
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT_DIR / f"baseline_model_comparison{SUFFIX}.csv", index=False)

    # Sidecar: the exact feature set used (for the table footnote + provenance).
    import json
    feat_meta = {"feature_set": FEATURE_SET,
                 "n_features": len(RESOLVED_FEATURES or []),
                 "features": RESOLVED_FEATURES or []}
    with open(OUT_DIR / f"feature_set{SUFFIX}.json", "w") as _fh:
        json.dump(feat_meta, _fh, indent=2)

    print(f"\n{'='*75}\n  SUMMARY  (mean ± std, no CDR, post-exclusion, {FEATURE_SET})\n{'='*75}")
    print(summary_df.to_string(index=False))
    print(f"\nSaved → {OUT_DIR / f'baseline_model_comparison{SUFFIX}.csv'}")
    print(f"Saved → {OUT_DIR / f'results_per_seed{SUFFIX}.csv'}")
    print(f"Saved → {OUT_DIR / f'feature_set{SUFFIX}.json'}  ({len(RESOLVED_FEATURES or [])} features)")
else:
    print("\n[ERROR] No results produced.")
