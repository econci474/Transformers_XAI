"""
02d_baseline_models_and_explainability_no_cdr.py
=================================================
Baseline models + explainability on the NO-CDR stratified data.
Combines logic from 02a (training) and 02c (explainability) but reads
from no_cdr_stratified/ and outputs to outputs/baseline_no_cdr/.

Run:
  python clinical_pipeline/02d_baseline_models_and_explainability_no_cdr.py
  python clinical_pipeline/02d_baseline_models_and_explainability_no_cdr.py --imputation
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
    f1_score, precision_score, recall_score
)
from sklearn.inspection       import permutation_importance
from xgboost import XGBClassifier
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

matplotlib.rcParams["font.family"] = "DejaVu Serif"
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument("--imputation", action="store_true",
                    help="Print imputation audit only, then exit")
args = parser.parse_args()

# ── Paths ──────────────────────────────────────────────────────────────────────
SPLIT_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified\tabular\baseline")
REPO_ROOT = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
OUT_DIR   = REPO_ROOT / "clinical_pipeline" / "outputs" / "baseline_no_cdr"
XAI_DIR   = OUT_DIR / "explainability"
OUT_DIR.mkdir(parents=True, exist_ok=True)
XAI_DIR.mkdir(parents=True, exist_ok=True)

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
    # CDR should already be removed from the data, but drop just in case
    "CDR_Global", "CDR_SumBoxes",
]

CATEGORICAL_ENCODE = [
    "Sex", "Handedness", "Marital_Status", "Retired", "Language", "Ethnicity",
    "APOE4_Status", "APOE_Alleles",
]

DIAGNOSTIC_CRITERIA_FEATURES = {
    "MMSE_Total":            "Mini-Mental State Exam — directly used in ADNI dx algorithm",
    "LogicalMemory_II_Total": "Logical Memory II — directly used in ADNI dx algorithm",
}


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
        metrics["AUC-ROC"]   = round(roc_auc_score(y_test, y_prob), 4)
        metrics["Precision"] = round(precision_score(y_test, y_pred, zero_division=0), 4)
        metrics["Recall"]    = round(recall_score(y_test, y_pred, zero_division=0), 4)
        metrics["F1"]        = round(f1_score(y_test, y_pred, zero_division=0), 4)
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


def extract_importance(pipe, model_name, feature_names, X_test, y_test):
    n_feat = len(feature_names)
    clf = pipe.named_steps["clf"]
    if model_name == "LogReg":
        coefs = np.abs(clf.coef_)
        imp = coefs.mean(axis=0) if coefs.ndim == 2 and coefs.shape[0] > 1 else coefs.flatten()
    elif model_name == "XGBoost":
        imp = clf.feature_importances_
    elif model_name == "SVM":
        result = permutation_importance(pipe, X_test, y_test, n_repeats=10,
                                        random_state=42, n_jobs=1)
        imp = np.clip(result.importances_mean, 0, None)
    else:
        imp = np.zeros(n_feat)
    if len(imp) != n_feat:
        padded = np.zeros(n_feat)
        padded[:min(len(imp), n_feat)] = imp[:n_feat]
        return padded
    return imp


# ── Task definitions ───────────────────────────────────────────────────────────
TASKS = {
    "T1_CN_vs_MCIAD": {
        "label_col": "Label_bl_multi", "task_type": "binary",
        "label_map": {"CN": 0, "MCI": 1, "AD": 1},
        "description": "Binary: CN vs MCI+AD",
        "class_names": {0: "CN", 1: "MCI+AD"},
    },
    "T2_Multiclass": {
        "label_col": "Label_bl_multi", "task_type": "multiclass",
        "label_map": {"CN": 0, "MCI": 1, "AD": 2},
        "description": "Multi-class: CN / MCI / AD",
        "class_names": {0: "CN", 1: "MCI", 2: "AD"},
    },
    "T3a_Conv3y": {
        "label_col": "Label_3y", "task_type": "binary", "label_map": None,
        "description": "Conversion to AD within 3 years", "filter": "non_ad",
        "class_names": {0: "non-conv", 1: "AD_3yr"},
    },
    "T3b_Conv5y": {
        "label_col": "Label_5y", "task_type": "binary", "label_map": None,
        "description": "Conversion to AD within 5 years", "filter": "non_ad",
        "class_names": {0: "non-conv", 1: "AD_5yr"},
    },
    "T3c_Conv10y": {
        "label_col": "Label_10y", "task_type": "binary", "label_map": None,
        "description": "Conversion to AD within 10 years", "filter": "non_ad",
        "class_names": {0: "non-conv", 1: "AD_10yr"},
    },
}


# ── Verify CDR is absent ──────────────────────────────────────────────────────
_s0 = SPLIT_DIR / "seed_0"
train = pd.read_csv(_s0 / "train.csv", low_memory=False)
val   = pd.read_csv(_s0 / "val.csv",   low_memory=False)
test  = pd.read_csv(_s0 / "test.csv",  low_memory=False)
trainval = pd.concat([train, val], ignore_index=True)

cdr_cols = [c for c in train.columns if "CDR" in c.upper()]
print(f"CDR columns found in data: {cdr_cols if cdr_cols else 'NONE ✓'}")
print(f"Seed_0 — Train: {len(train)}  Val: {len(val)}  Test: {len(test)}")

# ── Imputation audit ──────────────────────────────────────────────────────────
print(f"\n{'='*65}")
print("  IMPUTATION AUDIT  (Train+Val baseline, no CDR)")
print(f"{'='*65}")

X_audit, _, _ = prepare_features(trainval)
n_rows     = len(X_audit)
nan_counts = X_audit.isna().sum()
nan_audit  = (nan_counts[nan_counts > 0].rename("N_Missing").to_frame()
              .assign(Pct_Missing=lambda df: (df["N_Missing"] / n_rows * 100).round(1))
              .sort_values("N_Missing", ascending=False))
if nan_audit.empty:
    print("  No missing values.")
else:
    print(f"  Features with missing values ({len(nan_audit)} of {X_audit.shape[1]}):\n")
    for feat, row in nan_audit.iterrows():
        print(f"  {feat:<40} {int(row['N_Missing']):>10} {row['Pct_Missing']:>9.1f}%")
    nan_audit.to_csv(OUT_DIR / "imputation_audit.csv")

print(f"  Total features: {X_audit.shape[1]}")

if args.imputation:
    print("\n--imputation flag set — exiting.")
    sys.exit(0)


# ═══════════════════════════════════════════════════════════════════════════════
# PART 1: MODEL TRAINING (same as 02a but on no-CDR data)
# ═══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*65}")
print("  BASELINE MODELS — NO CDR (stratified splits)")
print(f"{'='*65}")

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
    print(f"\n{'#'*65}")
    print(f"  SEED {seed}  —  Train+Val: {len(trainval_s)}  Test: {len(test_s)}")
    print(f"{'#'*65}")

    for task_id, task_cfg in TASKS.items():
        print(f"\n  [{task_cfg['description']}]")
        lcol, lmap = task_cfg["label_col"], task_cfg["label_map"]
        filter_ = task_cfg.get("filter")

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

# ── Aggregate ─────────────────────────────────────────────────────────────────
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

    print(f"\n{'='*75}")
    print("  SUMMARY  (mean ± std, no CDR)")
    print(f"{'='*75}")
    print(summary_df.to_string(index=False))
    print(f"\nSaved → {OUT_DIR / 'baseline_model_comparison.csv'}")


# ═══════════════════════════════════════════════════════════════════════════════
# PART 2: EXPLAINABILITY (same as 02c but on no-CDR data)
# ═══════════════════════════════════════════════════════════════════════════════
print(f"\n\n{'='*70}")
print("  EXPLAINABILITY ANALYSIS — NO CDR")
print(f"{'='*70}")

leakage_rows = []

for task_id, task_cfg in TASKS.items():
    print(f"\n{'─'*60}")
    print(f"  Task: {task_cfg['description']}")
    print(f"{'─'*60}")

    lcol, lmap = task_cfg["label_col"], task_cfg["label_map"]
    filter_ = task_cfg.get("filter")

    def prep_split(df, _lcol=lcol, _lmap=lmap, _filter=filter_):
        df = df.copy()
        if _lmap:
            df[_lcol] = df[_lcol].map(_lmap)
        else:
            df[_lcol] = pd.to_numeric(df[_lcol], errors="coerce")
        if _filter == "non_ad":
            df = df[df["Label_bl_multi"].map({"CN": True, "MCI": True, "AD": False})]
        return df.dropna(subset=[_lcol])

    seed_importances = {m: [] for m in MODELS}
    feature_names = None
    canonical_enc = None
    test_counts = {}

    # Establish canonical features from seed 0
    sd0 = SPLIT_DIR / "seed_0"
    t0 = pd.read_csv(sd0 / "train.csv", low_memory=False)
    v0 = pd.read_csv(sd0 / "val.csv",   low_memory=False)
    tv0 = prep_split(pd.concat([t0, v0], ignore_index=True))
    if not tv0.empty:
        _, canonical_enc, feature_names = prepare_features(tv0)
        print(f"    Canonical features: {len(feature_names)} columns")

    if feature_names is None:
        print("    No data — skipping."); continue

    for seed in SEEDS:
        sd = SPLIT_DIR / f"seed_{seed}"
        ts = pd.read_csv(sd / "train.csv", low_memory=False)
        vs = pd.read_csv(sd / "val.csv",   low_memory=False)
        es = pd.read_csv(sd / "test.csv",  low_memory=False)
        tvs = pd.concat([ts, vs], ignore_index=True)

        tv_raw   = prep_split(tvs)
        test_raw = prep_split(es)
        if tv_raw.empty or test_raw.empty:
            continue

        y_tv   = tv_raw[lcol].astype(int)
        y_test = test_raw[lcol].astype(int)
        test_counts[seed] = dict(y_test.value_counts().sort_index())

        X_tv,   _, _ = prepare_features(tv_raw,   fit_encoders=canonical_enc,
                                         feature_cols=feature_names)
        X_test, _, _ = prepare_features(test_raw, fit_encoders=canonical_enc,
                                         feature_cols=feature_names)

        for model_name in MODELS:
            pipe = make_pipeline(model_name)
            pipe.fit(X_tv, y_tv)
            imp = extract_importance(pipe, model_name, feature_names, X_test, y_test)
            seed_importances[model_name].append(imp)
            acc = accuracy_score(y_test, pipe.predict(X_test))
            print(f"    [{model_name:8s}] seed {seed}  acc={acc:.4f}")

    # Test set distribution
    print(f"\n    Test set class distribution:")
    for seed, counts in test_counts.items():
        labels = ", ".join(f"{task_cfg['class_names'].get(k, k)}={v}"
                           for k, v in sorted(counts.items()))
        print(f"      seed {seed}: N={sum(counts.values())}  {labels}")

    # Aggregate importances
    importance_dfs = []
    for model_name in MODELS:
        arrs = np.array(seed_importances[model_name], dtype=float)
        imp_df = pd.DataFrame({
            "Feature": list(feature_names), "Model": model_name,
            "Importance_mean": arrs.mean(axis=0).tolist(),
            "Importance_std":  arrs.std(axis=0).tolist(),
        })
        importance_dfs.append(imp_df)

    all_imp = pd.concat(importance_dfs, ignore_index=True)
    all_imp.to_csv(XAI_DIR / f"feature_importance_{task_id}.csv", index=False)
    print(f"\n    Saved CSV → {XAI_DIR / f'feature_importance_{task_id}.csv'}")

    # Leakage audit
    for model_name in MODELS:
        model_imp = all_imp[all_imp["Model"] == model_name].sort_values(
            "Importance_mean", ascending=False).reset_index(drop=True)
        total_imp = model_imp["Importance_mean"].sum()
        for feat, reason in DIAGNOSTIC_CRITERIA_FEATURES.items():
            row = model_imp[model_imp["Feature"] == feat]
            if not row.empty:
                imp_val = row.iloc[0]["Importance_mean"]
                rank = model_imp.index.get_loc(row.index[0]) + 1
                pct = imp_val / total_imp * 100 if total_imp > 0 else 0
                leakage_rows.append({
                    "Task": task_cfg["description"], "Model": model_name,
                    "Feature": feat, "Rank": rank,
                    "Importance": round(imp_val, 6), "% of total": round(pct, 2),
                    "Reason": reason,
                })

    # Plot top-20 features
    fig, axes = plt.subplots(1, len(MODELS), figsize=(6 * len(MODELS), 8), sharey=False)
    if len(MODELS) == 1:
        axes = [axes]

    for ax, model_name in zip(axes, MODELS):
        model_imp = all_imp[all_imp["Model"] == model_name].sort_values(
            "Importance_mean", ascending=True)
        top20 = model_imp.tail(20)
        colors = ["#d62728" if f in DIAGNOSTIC_CRITERIA_FEATURES else "#1f77b4"
                  for f in top20["Feature"]]
        ax.barh(range(len(top20)), top20["Importance_mean"],
                xerr=top20["Importance_std"], color=colors,
                edgecolor="black", linewidth=0.3, capsize=2, height=0.7)
        ax.set_yticks(range(len(top20)))
        ax.set_yticklabels(top20["Feature"], fontsize=8)
        ax.set_xlabel("Importance (mean ± std)", fontsize=9)
        imp_type = {"LogReg": "|coef|", "XGBoost": "gain", "SVM": "permutation"}
        ax.set_title(f"{model_name}\n({imp_type.get(model_name, '')})",
                     fontsize=10, fontweight="bold")

    fig.legend(handles=[
        Patch(facecolor="#d62728", edgecolor="black", label="Used in ADNI dx algorithm"),
        Patch(facecolor="#1f77b4", edgecolor="black", label="Other feature"),
    ], loc="lower center", ncol=2, fontsize=9, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle(f"Feature Importance (no CDR): {task_cfg['description']}",
                 fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()
    fig.savefig(XAI_DIR / f"feature_importance_{task_id}.png", bbox_inches="tight", dpi=200)
    fig.savefig(XAI_DIR / f"feature_importance_{task_id}.pdf", bbox_inches="tight", dpi=200)
    plt.close()
    print(f"    Saved plot → {XAI_DIR / f'feature_importance_{task_id}.png'}")

# Leakage summary
if leakage_rows:
    leak_df = pd.DataFrame(leakage_rows)
    leak_df.to_csv(XAI_DIR / "leakage_audit.csv", index=False)
    print(f"\n{'='*70}")
    print("  LEAKAGE AUDIT (remaining diagnostic-criterion features)")
    print(f"{'='*70}")
    print(leak_df.to_string(index=False))
    print(f"\nSaved → {XAI_DIR / 'leakage_audit.csv'}")

print("\nDone.")
