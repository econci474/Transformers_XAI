"""
02c_baseline_models_explainability.py
======================================
Investigates *why* baseline clinical models achieve near-perfect accuracy by
extracting feature importances and flagging potentially leaky features.

Outputs (to outputs/baseline/explainability/):
  - feature_importance_{task}.csv   : per-feature importance (mean ± std across seeds)
  - feature_importance_{task}.png   : horizontal bar chart of top-20 features
  - leakage_audit.csv              : summary of suspected leaky features per task

Feature-importance sources:
  - LogReg  : |coef_|  (absolute logistic regression coefficients)
  - XGBoost : feature_importances_ (gain-based)
  - SVM     : permutation importance (RBF kernel has no coef_)

Run:
  python clinical_pipeline/02c_baseline_models_explainability.py
"""

import os, sys, warnings
import numpy as np
import pandas as pd
from pathlib import Path

from sklearn.pipeline         import Pipeline
from sklearn.preprocessing    import StandardScaler, LabelEncoder
from sklearn.impute            import SimpleImputer
from sklearn.linear_model     import LogisticRegression
from sklearn.svm              import SVC
from sklearn.inspection       import permutation_importance
from sklearn.metrics          import accuracy_score
from xgboost import XGBClassifier

import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"
warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────────
SPLIT_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\tabular\baseline")
REPO_ROOT = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
OUT_DIR   = REPO_ROOT / "clinical_pipeline" / "outputs" / "baseline" / "explainability"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SEEDS  = [0, 1, 2]
MODELS_CFG = ["LogReg", "SVM", "XGBoost"]

# ── Feature columns (replicated from 02a) ─────────────────────────────────────
DROP_ALWAYS = [
    "Patient_ID", "VISCODE_long", "VISCODE_2", "Date", "Generated_Text",
    "Label_bl_multi", "Label_visit_diag", "Label_visit_granular",
    "Label_1y", "Label_2y", "Label_3y", "Label_4y", "Label_5y",
    "Label_6y", "Label_7y", "Label_8y", "Label_9y", "Label_10y",
    "LogicalMemory_I_Total", "LogicalMemory_II_Total", "LogicalMemory_II_Cued",
    "Plasma_Abeta42", "Plasma_Abeta40", "Plasma_Abeta42_40",
    "Plasma_pTau217", "Plasma_pTau181", "Plasma_NfL", "Plasma_GFAP",
]

CATEGORICAL_ENCODE = [
    "Sex", "Handedness", "Marital_Status", "Retired", "Language", "Ethnicity",
    "APOE4_Status", "APOE_Alleles",
]

# Features that are part of the ADNI diagnostic criteria — likely leaky
DIAGNOSTIC_CRITERIA_FEATURES = {
    "CDR_Global":   "Clinical Dementia Rating — directly used in ADNI dx algorithm",
    "CDR_SumBoxes": "CDR Sum of Boxes — granular CDR composite",
    "MMSE_Total":   "Mini-Mental State Exam — cognitive screening used in dx",
    "ADAS_Total13": "ADAS-Cog 13 — AD severity scale, correlated with dx",
    "GDS_Total":    "Geriatric Depression Scale — used in MCI/CN differentiation",
    "FAQ_Total":    "Functional Assessment Questionnaire — functional impairment",
}


def prepare_features(df, fit_encoders=None, feature_cols=None):
    """Return (X_numeric, encoder_dict, feature_cols)."""
    df = df.copy()
    df = df.drop(columns=[c for c in DROP_ALWAYS if c in df.columns])
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
    df = df.apply(pd.to_numeric, errors="coerce")
    if feature_cols is None:
        feature_cols = list(df.columns)
    for col in feature_cols:
        if col not in df.columns:
            df[col] = np.nan
    df = df[feature_cols]
    return df, encoders, feature_cols


def make_pipeline(model_name):
    """Build Pipeline with imputer + scaler + classifier."""
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
            random_state=42, verbosity=0, scale_pos_weight=1,
        )
    else:
        raise ValueError(model_name)
    return Pipeline([("imp", imp), ("scl", scl), ("clf", clf)])


# ── Task definitions ───────────────────────────────────────────────────────────
TASKS = {
    "T1_CN_vs_MCIAD": {
        "label_col":   "Label_bl_multi",
        "task_type":   "binary",
        "label_map":   {"CN": 0, "MCI": 1, "AD": 1},
        "description": "Binary: CN vs MCI+AD",
        "class_names": {0: "CN", 1: "MCI+AD"},
    },
    "T2_Multiclass": {
        "label_col":   "Label_bl_multi",
        "task_type":   "multiclass",
        "label_map":   {"CN": 0, "MCI": 1, "AD": 2},
        "description": "Multi-class: CN / MCI / AD",
        "class_names": {0: "CN", 1: "MCI", 2: "AD"},
    },
    "T3a_Conv3y": {
        "label_col":   "Label_3y",
        "task_type":   "binary",
        "label_map":   None,
        "description": "Conversion to AD within 3 years",
        "filter":      "non_ad",
        "class_names": {0: "non-conv", 1: "AD_3yr"},
    },
    "T3b_Conv5y": {
        "label_col":   "Label_5y",
        "task_type":   "binary",
        "label_map":   None,
        "description": "Conversion to AD within 5 years",
        "filter":      "non_ad",
        "class_names": {0: "non-conv", 1: "AD_5yr"},
    },
    "T3c_Conv10y": {
        "label_col":   "Label_10y",
        "task_type":   "binary",
        "label_map":   None,
        "description": "Conversion to AD within 10 years",
        "filter":      "non_ad",
        "class_names": {0: "non-conv", 1: "AD_10yr"},
    },
}


def extract_importance(pipe, model_name, feature_names, X_test, y_test):
    """Extract feature importance array aligned to feature_names."""
    n_feat = len(feature_names)
    clf = pipe.named_steps["clf"]

    if model_name == "LogReg":
        coefs = np.abs(clf.coef_)
        if coefs.ndim == 2 and coefs.shape[0] > 1:
            # Multi-class: mean absolute coefficient across classes
            imp = coefs.mean(axis=0)
        else:
            imp = coefs.flatten()

    elif model_name == "XGBoost":
        imp = clf.feature_importances_

    elif model_name == "SVM":
        # RBF kernel: use permutation importance (n_jobs=1 for Windows stability)
        result = permutation_importance(
            pipe, X_test, y_test, n_repeats=10, random_state=42, n_jobs=1)
        imp = np.clip(result.importances_mean, 0, None)

    else:
        imp = np.zeros(n_feat)

    # Safety: verify length matches expected features
    if len(imp) != n_feat:
        print(f"      WARNING: {model_name} returned {len(imp)} importances, "
              f"expected {n_feat}. Padding/truncating.")
        padded = np.zeros(n_feat)
        padded[:min(len(imp), n_feat)] = imp[:n_feat]
        return padded
    return imp


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN LOOP
# ═══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("  BASELINE MODEL EXPLAINABILITY ANALYSIS")
print("=" * 70)

leakage_rows = []

for task_id, task_cfg in TASKS.items():
    print(f"\n{'─'*60}")
    print(f"  Task: {task_cfg['description']}")
    print(f"{'─'*60}")

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

    # Accumulate importances across seeds: {model: [array_seed0, ...]}
    seed_importances = {m: [] for m in MODELS_CFG}
    feature_names = None
    canonical_enc = None
    test_counts = {}  # {seed: {class_label: count}}

    # First pass: establish canonical feature columns from seed 0
    seed_dir_0 = SPLIT_DIR / "seed_0"
    train_0 = pd.read_csv(seed_dir_0 / "train.csv", low_memory=False)
    val_0   = pd.read_csv(seed_dir_0 / "val.csv",   low_memory=False)
    tv_raw_0 = prep_split(pd.concat([train_0, val_0], ignore_index=True))
    if not tv_raw_0.empty:
        _, canonical_enc, feature_names = prepare_features(tv_raw_0)
        print(f"    Canonical features: {len(feature_names)} columns (from seed 0)")

    if feature_names is None:
        print(f"    No data for this task — skipping.")
        continue

    for seed in SEEDS:
        seed_dir = SPLIT_DIR / f"seed_{seed}"
        train_s = pd.read_csv(seed_dir / "train.csv", low_memory=False)
        val_s   = pd.read_csv(seed_dir / "val.csv",   low_memory=False)
        test_s  = pd.read_csv(seed_dir / "test.csv",  low_memory=False)
        trainval_s = pd.concat([train_s, val_s], ignore_index=True)

        tv_raw   = prep_split(trainval_s)
        test_raw = prep_split(test_s)
        if tv_raw.empty or test_raw.empty:
            continue

        y_tv   = tv_raw[lcol].astype(int)
        y_test = test_raw[lcol].astype(int)

        # Track test set class distribution
        test_counts[seed] = dict(y_test.value_counts().sort_index())

        # Use canonical encoders and feature columns for all seeds
        X_tv,   _, _ = prepare_features(tv_raw,   fit_encoders=canonical_enc,
                                         feature_cols=feature_names)
        X_test, _, _ = prepare_features(test_raw, fit_encoders=canonical_enc,
                                         feature_cols=feature_names)

        for model_name in MODELS_CFG:
            pipe = make_pipeline(model_name)
            pipe.fit(X_tv, y_tv)
            imp = extract_importance(pipe, model_name, feature_names, X_test, y_test)
            seed_importances[model_name].append(imp)
            acc = accuracy_score(y_test, pipe.predict(X_test))
            print(f"    [{model_name:8s}] seed {seed}  acc={acc:.4f}")

    if feature_names is None:
        continue

    # ── Print test set class distribution ─────────────────────────────────
    print(f"\n    Test set class distribution:")
    for seed, counts in test_counts.items():
        labels = ", ".join(f"{task_cfg['class_names'].get(k, k)}={v}"
                           for k, v in sorted(counts.items()))
        print(f"      seed {seed}: N={sum(counts.values())}  {labels}")

    # ── Aggregate importances: mean ± std across seeds ────────────────────
    importance_dfs = []
    for model_name in MODELS_CFG:
        imp_list = seed_importances[model_name]
        # Debug: print lengths
        lens = [len(a) for a in imp_list]
        print(f"    {model_name}: importance array lengths = {lens} "
              f"(expected {len(feature_names)})")
        arrs = np.array(imp_list, dtype=float)  # shape (n_seeds, n_features)
        imp_mean = arrs.mean(axis=0)
        imp_std  = arrs.std(axis=0)
        imp_df = pd.DataFrame({
            "Feature":         list(feature_names),
            "Model":           [model_name] * len(feature_names),
            "Importance_mean": imp_mean.tolist(),
            "Importance_std":  imp_std.tolist(),
        })
        importance_dfs.append(imp_df)

    all_imp = pd.concat(importance_dfs, ignore_index=True)
    csv_path = OUT_DIR / f"feature_importance_{task_id}.csv"
    all_imp.to_csv(csv_path, index=False)
    print(f"\n    Saved CSV → {csv_path}")

    # ── Leakage audit: flag diagnostic-criteria features ──────────────────
    for model_name in MODELS_CFG:
        model_imp = all_imp[all_imp["Model"] == model_name].copy()
        model_imp = model_imp.sort_values("Importance_mean", ascending=False)
        total_imp = model_imp["Importance_mean"].sum()

        for feat, reason in DIAGNOSTIC_CRITERIA_FEATURES.items():
            row = model_imp[model_imp["Feature"] == feat]
            if not row.empty:
                imp_val  = row.iloc[0]["Importance_mean"]
                rank     = int(model_imp.index.get_loc(row.index[0])) + 1
                pct_total = imp_val / total_imp * 100 if total_imp > 0 else 0
                leakage_rows.append({
                    "Task":       task_cfg["description"],
                    "Model":      model_name,
                    "Feature":    feat,
                    "Rank":       rank,
                    "Importance": round(imp_val, 6),
                    "% of total": round(pct_total, 2),
                    "Reason":     reason,
                })

    # ── Plot: top-20 features per model (side-by-side subplot) ────────────
    fig, axes = plt.subplots(1, len(MODELS_CFG), figsize=(6 * len(MODELS_CFG), 8),
                              sharey=False)
    if len(MODELS_CFG) == 1:
        axes = [axes]

    for ax, model_name in zip(axes, MODELS_CFG):
        model_imp = all_imp[all_imp["Model"] == model_name].copy()
        model_imp = model_imp.sort_values("Importance_mean", ascending=True)
        top20 = model_imp.tail(20)

        colors = []
        for feat in top20["Feature"]:
            if feat in DIAGNOSTIC_CRITERIA_FEATURES:
                colors.append("#d62728")  # red for leaky
            else:
                colors.append("#1f77b4")  # blue for normal
        
        ax.barh(range(len(top20)), top20["Importance_mean"],
                xerr=top20["Importance_std"], color=colors,
                edgecolor="black", linewidth=0.3, capsize=2, height=0.7)
        ax.set_yticks(range(len(top20)))
        ax.set_yticklabels(top20["Feature"], fontsize=8)
        ax.set_xlabel("Importance (mean ± std, 3 seeds)", fontsize=9)

        imp_type = {"LogReg": "|coef|", "XGBoost": "gain",
                     "SVM": "permutation"}
        ax.set_title(f"{model_name}\n({imp_type.get(model_name, '')})",
                      fontsize=10, fontweight="bold")
        ax.axvline(0, color="black", linewidth=0.5)

    # Legend
    from matplotlib.patches import Patch
    legend_patches = [
        Patch(facecolor="#d62728", edgecolor="black", label="Diagnostic criterion"),
        Patch(facecolor="#1f77b4", edgecolor="black", label="Other feature"),
    ]
    fig.legend(handles=legend_patches, loc="lower center", ncol=2,
               fontsize=9, frameon=True, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle(f"Feature Importance: {task_cfg['description']}",
                  fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()
    fig_path = OUT_DIR / f"feature_importance_{task_id}.png"
    fig.savefig(fig_path, bbox_inches="tight", dpi=200)
    fig.savefig(fig_path.with_suffix(".pdf"), bbox_inches="tight", dpi=200)
    plt.close()
    print(f"    Saved plot → {fig_path}")

# ── Leakage audit summary ────────────────────────────────────────────────────
if leakage_rows:
    leak_df = pd.DataFrame(leakage_rows)
    leak_path = OUT_DIR / "leakage_audit.csv"
    leak_df.to_csv(leak_path, index=False)
    print(f"\n{'='*70}")
    print("  LEAKAGE AUDIT: Diagnostic-criterion features in top importances")
    print(f"{'='*70}")
    print(leak_df.to_string(index=False))
    print(f"\nSaved → {leak_path}")

print("\nDone.")
