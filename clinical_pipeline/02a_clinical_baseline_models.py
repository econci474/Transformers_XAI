"""
clinical_baseline_models.py
============================
Trains and evaluates Logistic Regression, SVM, and XGBoost on tabular baseline
clinical data across three tasks:

  Task 1: Binary classification      — CN vs MCI+AD
  Task 2: Multi-class classification — CN / MCI / AD
  Task 3a: Conversion prediction     — conversion to AD within 3 years
  Task 3b: Conversion prediction     — conversion to AD within 5 years

Input:  D:\\ADNI_BIDS_project\\derivatives\\clinical\\tabular\\baseline\\
Output: D:\\ADNI_BIDS_project\\derivatives\\clinical\\tabular\\baseline_model_results\\
"""

import os, sys, warnings, argparse
import numpy as np
import pandas as pd
from pathlib import Path

from sklearn.pipeline         import Pipeline
from sklearn.preprocessing    import StandardScaler, LabelEncoder
from sklearn.impute            import SimpleImputer
from sklearn.linear_model     import LogisticRegression
from sklearn.svm              import SVC
from sklearn.metrics          import (
    accuracy_score, balanced_accuracy_score, roc_auc_score,
    f1_score, classification_report
)
from xgboost import XGBClassifier

warnings.filterwarnings("ignore")

# ── CLI flags ──────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--imputation", action="store_true",
                    help="Print imputation audit only, then exit (no model training)")
args = parser.parse_args()

# ── Paths ──────────────────────────────────────────────────────────────────────
SPLIT_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\tabular\baseline")
OUT_DIR   = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\clinical_pipeline\outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Load split (seed_0 used for imputation audit) ─────────────────────────────
_s0 = SPLIT_DIR / "seed_0"
train = pd.read_csv(_s0 / "train.csv", low_memory=False)
val   = pd.read_csv(_s0 / "val.csv",   low_memory=False)
test  = pd.read_csv(_s0 / "test.csv",  low_memory=False)
trainval = pd.concat([train, val], ignore_index=True)
print(f"Seed_0 — Train: {len(train)}  Val: {len(val)}  Test: {len(test)}")


# ── Feature columns ────────────────────────────────────────────────────────────
DROP_ALWAYS = [
    "Patient_ID", "VISCODE_long", "VISCODE_2", "Date", "Generated_Text",
    "Label_bl_multi", "Label_visit_diag", "Label_visit_granular",
    "Label_1y", "Label_2y", "Label_3y", "Label_4y", "Label_5y",
    "Label_6y", "Label_7y", "Label_8y", "Label_9y", "Label_10y",
    # Logical Memory not administered at the baseline visit:
    "LogicalMemory_I_Total", "LogicalMemory_II_Total", "LogicalMemory_II_Cued",
    # Plasma biomarkers — too sparse at baseline; not imputed:
    "Plasma_Abeta42", "Plasma_Abeta40", "Plasma_Abeta42_40",
    "Plasma_pTau217", "Plasma_pTau181", "Plasma_NfL", "Plasma_GFAP",
]

# Categorical columns to label-encode (string → integer)
CATEGORICAL_ENCODE = [
    "Sex", "Handedness", "Marital_Status", "Retired", "Language", "Ethnicity",
    "APOE4_Status", "APOE_Alleles",
]

def prepare_features(df, fit_encoders=None, feature_cols=None):
    """
    Return (X_numeric, encoder_dict, feature_cols).
    - fit_encoders=None  → fit new encoders (training call)
    - feature_cols       → re-index to this exact column list (test call)
    """
    df = df.copy()
    df = df.drop(columns=[c for c in DROP_ALWAYS if c in df.columns])

    # Encode categoricals
    encoders = fit_encoders or {}
    for col in CATEGORICAL_ENCODE:
        if col not in df.columns:
            continue
        df[col] = df[col].astype(str).str.strip().replace({"NaN": np.nan, "nan": np.nan})
        if fit_encoders is None:
            le = LabelEncoder()
            non_null = df[col].dropna()
            le.fit(non_null.astype(str))
            encoders[col] = le
        le = encoders[col]
        df[col] = df[col].apply(
            lambda v: le.transform([str(v)])[0]
            if pd.notna(v) and str(v) in le.classes_ else np.nan
        )

    # Convert everything to numeric (handles leftover string "NaN" values
    # and ensures all-null columns stay as float64, not object)
    df = df.apply(pd.to_numeric, errors="coerce")

    # On training call: record which columns survived
    if feature_cols is None:
        feature_cols = list(df.columns)

    # On test call: align to training columns exactly
    for col in feature_cols:
        if col not in df.columns:
            df[col] = np.nan
    df = df[feature_cols]

    return df, encoders, feature_cols


def make_pipeline(model_name):
    """Build a Pipeline with imputer + scaler + classifier."""
    imp = SimpleImputer(strategy="median")
    scl = StandardScaler()

    if model_name == "LogReg":
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", random_state=42)
    elif model_name == "SVM":
        clf = SVC(probability=True, class_weight="balanced", random_state=42)
    elif model_name == "XGBoost":
        clf = XGBClassifier(
            n_estimators=300, max_depth=4, learning_rate=0.05,
            use_label_encoder=False, eval_metric="logloss",
            random_state=42, verbosity=0,
            scale_pos_weight=1,   # overridden per task if needed
        )
    else:
        raise ValueError(model_name)

    return Pipeline([("imp", imp), ("scl", scl), ("clf", clf)])


def evaluate(pipe, X_test, y_test, task_type="binary", label_names=None):
    """Return metrics dict."""
    y_pred = pipe.predict(X_test)
    metrics = {
        "Accuracy":          round(accuracy_score(y_test, y_pred), 4),
        "BalancedAcc":       round(balanced_accuracy_score(y_test, y_pred), 4),
    }
    if task_type == "binary":
        y_prob = pipe.predict_proba(X_test)[:, 1]
        metrics["AUC-ROC"] = round(roc_auc_score(y_test, y_prob), 4)
        metrics["F1"]      = round(f1_score(y_test, y_pred, zero_division=0), 4)
    else:
        y_prob = pipe.predict_proba(X_test)
        try:
            metrics["AUC-ROC (OvR)"] = round(
                roc_auc_score(y_test, y_prob, multi_class="ovr", average="macro"), 4)
        except Exception:
            metrics["AUC-ROC (OvR)"] = np.nan
        metrics["MacroF1"] = round(f1_score(y_test, y_pred, average="macro", zero_division=0), 4)
    return metrics


# ── Task definitions ───────────────────────────────────────────────────────────
TASKS = {
    "T1_CN_vs_MCIAD": {
        "label_col":   "Label_bl_multi",
        "task_type":   "binary",
        "label_map":   {"CN": 0, "MCI": 1, "AD": 1},
        "description": "Binary: CN vs MCI+AD",
    },
    "T2_Multiclass": {
        "label_col":   "Label_bl_multi",
        "task_type":   "multiclass",
        "label_map":   {"CN": 0, "MCI": 1, "AD": 2},
        "description": "Multi-class: CN / MCI / AD",
    },
    "T3a_Conv3y": {
        "label_col":   "Label_3y",
        "task_type":   "binary",
        "label_map":   None,   # already 0/1
        "description": "Conversion to AD within 3 years",
        "filter":      "non_ad",   # exclude subjects already AD at baseline
    },
    "T3b_Conv5y": {
        "label_col":   "Label_5y",
        "task_type":   "binary",
        "label_map":   None,
        "description": "Conversion to AD within 5 years",
        "filter":      "non_ad",
    },
    "T3c_Conv10y": {
        "label_col":   "Label_10y",
        "task_type":   "binary",
        "label_map":   None,
        "description": "Conversion to AD within 10 years",
        "filter":      "non_ad",
    },
}

MODELS = ["LogReg", "SVM", "XGBoost"]

# ── Imputation audit (run once on full train+val split) ────────────────────────
print(f"\n{'='*65}")
print("  IMPUTATION AUDIT  (Train+Val baseline split, n={})".format(len(trainval)))
print(f"{'='*65}")

X_audit, _, _ = prepare_features(trainval)
n_rows      = len(X_audit)
nan_counts  = X_audit.isna().sum()
nan_audit   = (
    nan_counts[nan_counts > 0]
    .rename("N_Missing")
    .to_frame()
    .assign(Pct_Missing=lambda df: (df["N_Missing"] / n_rows * 100).round(1))
    .sort_values("N_Missing", ascending=False)
)
if nan_audit.empty:
    print("  No missing values \u2014 no imputation required.")
else:
    print(f"  Features with missing values ({len(nan_audit)} of {X_audit.shape[1]} total features):\n")
    print(f"  {'Feature':<40} {'N_Missing':>10} {'% Missing':>10}")
    print("  " + "-" * 62)
    for feat, row in nan_audit.iterrows():
        print(f"  {feat:<40} {int(row['N_Missing']):>10} {row['Pct_Missing']:>9.1f}%")
    nan_audit.to_csv(OUT_DIR / "imputation_audit.csv")
    print(f"\n  Saved audit \u2192 {OUT_DIR / 'imputation_audit.csv'}")

print(f"  Features with NO missing values: "
      f"{(nan_counts == 0).sum()} / {X_audit.shape[1]}")
print(f"  Imputation strategy: median (per-feature, fitted on train+val only)")

if args.imputation:
    print("\n--imputation flag set \u2014 exiting without running models.")
    sys.exit(0)
print()

# \u2500\u2500 Run all tasks across 3 seeds \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
SEEDS = [0, 1, 2]
all_seed_results = []

for seed in SEEDS:
    seed_dir = SPLIT_DIR / f"seed_{seed}"
    try:
        train_s = pd.read_csv(seed_dir / "train.csv", low_memory=False)
        val_s   = pd.read_csv(seed_dir / "val.csv",   low_memory=False)
        test_s  = pd.read_csv(seed_dir / "test.csv",  low_memory=False)
    except FileNotFoundError as e:
        print(f"  [seed {seed}] Missing file: {e} — skipping seed.")
        continue

    trainval_s = pd.concat([train_s, val_s], ignore_index=True)
    print(f"\n{'#'*65}")
    print(f"  SEED {seed}  —  Train+Val: {len(trainval_s)}  Test: {len(test_s)}")
    print(f"{'#'*65}")

    for task_id, task_cfg in TASKS.items():
        print(f"\n  [{task_cfg['description']}]")
        lcol    = task_cfg["label_col"]
        lmap    = task_cfg["label_map"]
        filter_ = task_cfg.get("filter")

        def prep_split(df):
            df = df.copy()
            if lmap:
                df[lcol] = df[lcol].map(lmap)
            else:
                df[lcol] = pd.to_numeric(df[lcol], errors="coerce")
            if filter_ == "non_ad":
                df = df[df["Label_bl_multi"].map({"CN": True, "MCI": True, "AD": False})]
            return df.dropna(subset=[lcol])

        tv_raw   = prep_split(trainval_s)
        test_raw = prep_split(test_s)

        if tv_raw.empty or test_raw.empty:
            print(f"    Insufficient data — skipping.")
            continue

        y_tv   = tv_raw[lcol].astype(int)
        y_test = test_raw[lcol].astype(int)
        print(f"    Train+Val: {len(tv_raw)} | Test: {len(test_raw)} | "
              f"Classes (TV): {dict(y_tv.value_counts().sort_index())}")

        X_tv,   enc, feat_cols = prepare_features(tv_raw)
        X_test, _,   _         = prepare_features(test_raw, fit_encoders=enc, feature_cols=feat_cols)

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

# ── Aggregate across seeds: mean ± std ────────────────────────────────────────
if all_seed_results:
    raw_df    = pd.DataFrame(all_seed_results)
    meta_cols = ["Task", "Model", "Task_ID", "Seed"]
    metric_cols = [c for c in raw_df.columns if c not in meta_cols]

    # Save per-seed results
    raw_df[["Task", "Model", "Seed"] + metric_cols].to_csv(
        OUT_DIR / "results_per_seed.csv", index=False)

    # Aggregate
    agg = raw_df.groupby(["Task_ID", "Task", "Model"])[metric_cols].agg(["mean", "std"])
    agg.columns = [f"{m}_{stat}" for m, stat in agg.columns]
    agg = agg.reset_index()

    # Build readable summary: "mean ± std" per metric
    summary_rows = []
    for _, row in agg.iterrows():
        entry = {"Task": row["Task"], "Model": row["Model"]}
        for m in metric_cols:
            mu  = row.get(f"{m}_mean", np.nan)
            sd  = row.get(f"{m}_std",  np.nan)
            if pd.notna(mu):
                entry[m] = f"{mu:.4f} ± {sd:.4f}" if pd.notna(sd) else f"{mu:.4f}"
            else:
                entry[m] = ""
        summary_rows.append(entry)
    summary_df = pd.DataFrame(summary_rows)

    out_path = OUT_DIR / "baseline_model_comparison.csv"
    summary_df.to_csv(out_path, index=False)

    print(f"\n{'='*75}")
    print("  SUMMARY  (mean ± std across seeds 0, 1, 2)")
    print(f"{'='*75}")
    print(summary_df.to_string(index=False))
    print(f"\nSaved per-seed results → {OUT_DIR / 'results_per_seed.csv'}")
    print(f"Saved summary          → {out_path}")
