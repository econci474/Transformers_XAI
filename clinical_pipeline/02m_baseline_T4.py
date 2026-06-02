"""
02m_baseline_T4.py
==================
Tabular baseline models (LogReg / SVM / XGBoost) for the **T4 AD-conversion horizon**
task (3-class: <3y / 3–7y / ≥7y) over the converter cohort (pMCI + pCN_to_AD), so the
T4 column of the combined table gets a fair baseline-vs-encoder comparison.

Why a separate script: `tabular/baseline_T4` (built by 01e from the MRI master) carries
only labels + conversion columns — NOT the cognitive feature battery. So we JOIN the
full no-CDR feature battery from `tabular/baseline` (per-subject baseline, seed-independent)
onto the baseline_T4 cohort by Patient_ID, then train on the SAME folds as the T4 encoder.

Protocol matches 02h (train-fit; val + test held-out; median-impute + standardise fit on
train; class-weighted). 3-class multiclass metrics.

Outputs (→ outputs/baseline_no_cdr_post_exclusion/):
  baseline_model_comparison_T4.csv   (Task, Model, Split, multiclass metrics; read by 04d)
  results_per_seed_T4.csv

Run (local Windows, env: snp/clinical):
  python clinical_pipeline/02m_baseline_T4.py
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
BASE_DERIV   = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion")
FEATURE_SRC  = BASE_DERIV / "tabular" / "baseline"      # 60-feature battery (per-subject)
T4_SPLIT_DIR = BASE_DERIV / "tabular" / "baseline_T4"   # converter folds + Label_T4
REPO_ROOT    = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
OUT_DIR      = REPO_ROOT / "clinical_pipeline" / "outputs" / "baseline_no_cdr_post_exclusion"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SEEDS  = [0, 1, 2]
MODELS = ["LogReg", "SVM", "XGBoost"]
TASK_LABEL = "AD horizon (3-class): pMCI + pCN_to_AD"   # MUST match 04d's T4 bl_task

# ── Feature columns (same no-CDR battery as 02h; + drop Label_T4, the label) ───
DROP_ALWAYS = [
    "Patient_ID", "VISCODE_long", "VISCODE_2", "Date", "Generated_Text",
    "bids_sub", "bids_ses", "adni_viscode", "mri_date", "clinical_date",
    "days_diff", "match_status",
    "Label_bl_multi", "Label_visit_diag", "Label_visit_granular",
    "Label_1y", "Label_2y", "Label_3y", "Label_4y", "Label_5y",
    "Label_6y", "Label_7y", "Label_8y", "Label_9y", "Label_10y",
    "Label_T4",                       # ← the T4 label; must NOT be a feature
    "LogicalMemory_I_Total", "LogicalMemory_II_Total", "LogicalMemory_II_Cued",
    "Plasma_Abeta42", "Plasma_Abeta40", "Plasma_Abeta42_40",
    "Plasma_pTau217", "Plasma_pTau181", "Plasma_NfL", "Plasma_GFAP",
    "CDR_Global", "CDR_SumBoxes",
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
            le = LabelEncoder(); le.fit(df[col].dropna().astype(str)); encoders[col] = le
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
                            use_label_encoder=False, eval_metric="mlogloss",
                            random_state=42, verbosity=0)
    else:
        raise ValueError(model_name)
    return Pipeline([("imp", imp), ("scl", scl), ("clf", clf)])


def evaluate(pipe, X, y):
    """3-class multiclass metrics."""
    y_pred = pipe.predict(X)
    metrics = {
        "Accuracy":    round(accuracy_score(y, y_pred), 4),
        "BalancedAcc": round(balanced_accuracy_score(y, y_pred), 4),
    }
    y_prob = pipe.predict_proba(X)
    try:
        metrics["AUC-ROC (OvR)"] = round(
            roc_auc_score(y, y_prob, multi_class="ovr", average="macro"), 4)
    except Exception:
        metrics["AUC-ROC (OvR)"] = np.nan
    metrics["MacroPrecision"] = round(precision_score(y, y_pred, average="macro", zero_division=0), 4)
    metrics["MacroRecall"]    = round(recall_score(y, y_pred, average="macro", zero_division=0), 4)
    metrics["MacroF1"]        = round(f1_score(y, y_pred, average="macro", zero_division=0), 4)
    return metrics


# ── Build Patient_ID → feature-row map (from tabular/baseline; baseline, seed-indep) ──
_feat_frames = [pd.read_csv(FEATURE_SRC / "seed_0" / f"{sp}.csv", low_memory=False)
                for sp in ("train", "val", "test")]
feat_map = pd.concat(_feat_frames, ignore_index=True)
feat_map["Patient_ID"] = feat_map["Patient_ID"].astype(str)
feat_map = feat_map.drop_duplicates("Patient_ID", keep="first").set_index("Patient_ID")
feat_map["_has_features"] = 1
print(f"Feature battery: {feat_map.shape[1]-1} cols for {len(feat_map)} subjects "
      f"(from {FEATURE_SRC.name})")


def load_t4_fold(seed, split):
    """Return (X-ready joined df, y) for one baseline_T4 fold, features joined in."""
    t4 = pd.read_csv(T4_SPLIT_DIR / f"seed_{seed}" / f"{split}.csv", low_memory=False)
    t4["Patient_ID"] = t4["Patient_ID"].astype(str)
    keep = t4[["Patient_ID", "Label_T4"]].dropna(subset=["Label_T4"]).copy()
    joined = keep.join(feat_map, on="Patient_ID")
    missing = joined["_has_features"].isna()
    if missing.any():
        print(f"    [seed {seed}/{split}] DROPPED {int(missing.sum())} converters with no "
              f"feature row: {sorted(joined.loc[missing, 'Patient_ID'])}")
    joined = joined[~missing].drop(columns=["_has_features"])
    y = joined["Label_T4"].astype(int)
    return joined, y


# ═══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*65}\n  BASELINE MODELS — T4 (3-class horizon), NO CDR, POST-EXCLUSION\n{'='*65}")
all_seed_results = []
RESOLVED_FEATURES = None

for seed in SEEDS:
    try:
        train_raw, y_train = load_t4_fold(seed, "train")
        val_raw,   y_val   = load_t4_fold(seed, "val")
        test_raw,  y_test  = load_t4_fold(seed, "test")
    except FileNotFoundError as e:
        print(f"  [seed {seed}] Missing: {e}"); continue

    print(f"\n{'#'*65}\n  SEED {seed}  —  Train: {len(train_raw)}  Val: {len(val_raw)}  "
          f"Test: {len(test_raw)}  |  Train classes: "
          f"{dict(y_train.value_counts().sort_index())}\n{'#'*65}")

    X_train, enc, feat_cols = prepare_features(train_raw)
    if RESOLVED_FEATURES is None:
        RESOLVED_FEATURES = list(feat_cols)
        print(f"  Feature columns ({len(feat_cols)}): {feat_cols}")
    X_val,  _, _ = prepare_features(val_raw,  fit_encoders=enc, feature_cols=feat_cols)
    X_test, _, _ = prepare_features(test_raw, fit_encoders=enc, feature_cols=feat_cols)

    for model_name in MODELS:
        pipe = make_pipeline(model_name)
        try:
            pipe.fit(X_train, y_train)
        except Exception as e:
            print(f"    [{model_name:8s}] FIT ERROR: {e}"); continue
        for split_name, X_split, y_split in (("val", X_val, y_val), ("test", X_test, y_test)):
            try:
                metrics = evaluate(pipe, X_split, y_split)
                metrics.update({"Model": model_name, "Task": TASK_LABEL,
                                "Task_ID": "T4_conv_horizon", "Seed": seed, "Split": split_name})
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
        OUT_DIR / "results_per_seed_T4.csv", index=False)

    agg = raw_df.groupby(["Task_ID", "Task", "Model", "Split"])[metric_cols].agg(["mean", "std"])
    agg.columns = [f"{m}_{stat}" for m, stat in agg.columns]
    agg = agg.reset_index()

    summary_rows = []
    for _, row in agg.iterrows():
        entry = {"Task": row["Task"], "Model": row["Model"], "Split": row["Split"]}
        for m in metric_cols:
            mu = row.get(f"{m}_mean", np.nan); sd = row.get(f"{m}_std", np.nan)
            entry[m] = (f"{mu:.4f} ± {sd:.4f}" if pd.notna(sd) else f"{mu:.4f}") if pd.notna(mu) else ""
        summary_rows.append(entry)
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT_DIR / "baseline_model_comparison_T4.csv", index=False)

    print(f"\n{'='*70}\n  T4 SUMMARY (mean ± std, 3-class, converters)\n{'='*70}")
    print(summary_df.to_string(index=False))
    print(f"\nSaved → {OUT_DIR / 'baseline_model_comparison_T4.csv'}")
    print(f"Saved → {OUT_DIR / 'results_per_seed_T4.csv'}")
else:
    print("\n[ERROR] No results produced.")
