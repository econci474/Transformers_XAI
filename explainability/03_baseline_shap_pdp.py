r"""
03_baseline_shap_pdp.py   (env: clinical + `pip install shap`; CPU)
===================================================================
SHAP values + partial-dependence (PDP) plots for the BEST tabular baseline per explainability task
(T1d / T1e / T2). Best baseline = XGBoost for all three (post-exclusion leaderboard, val balanced
accuracy: T1d 0.755, T1e 0.442[near chance], T2 0.817).

Refits the baseline per seed exactly as in 02h_baseline_post_exclusion.py (same DROP_ALWAYS / categorical
encoding / median-impute + standardise pipeline / 60-feature set) so the explained model IS the leaderboard
model. Helpers replicated (not imported) because 02h trains at import time.

Outputs (per task) -> explainability/<T1d|T1e|T2>/: shap_beeswarm, shap_bar, pdp_grid (png/pdf),
shap_values_seed0.csv, shap_global_importance.csv. CAVEAT: T1e baseline near chance.

Run:  conda run -n clinical python explainability/03_baseline_shap_pdp.py
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.pipeline import Pipeline
from sklearn.inspection import PartialDependenceDisplay
from xgboost import XGBClassifier
import shap

matplotlib.rcParams["font.family"] = "DejaVu Serif"

HERE = Path(__file__).resolve().parent
# Figures/CSVs are OUTPUTS — they live on D: (off the code repo), like every other derivative.
OUT_BASE = Path(r"D:\ADNI_BIDS_project\derivatives\explainability")
SPLIT_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion"
                 r"\tabular\baseline")
SEEDS = (0, 1, 2)

TASKS = {
    "T1d": dict(label_col="conversion_group", label_map={"sMCI": 0, "pMCI": 1},
                task_type="binary", classes=["sMCI", "pMCI"],
                title="T2a - pMCI vs sMCI (baseline MCI)"),
    "T1e": dict(label_col="conversion_group", label_map={"sCN": 0, "pCN_to_MCI": 1, "pCN_to_AD": 1},
                task_type="binary", classes=["sCN", "pCN"],
                title="T2b - sCN vs pCN (baseline CN)"),
    "T2":  dict(label_col="Label_bl_multi", label_map={"CN": 0, "MCI": 1, "AD": 2},
                task_type="multiclass", classes=["CN", "MCI", "AD"],
                title="T1c - CN / MCI / AD"),
}
CHANCE_NOTE = {"T1e": "  (val bACC 0.44)"}

# Harmonised class colours (shared verbatim with 02_attention_rollout.py) for the multiclass T1c figure.
CLASS_HEX = {"CN": "#CCFFCC", "MCI": "#FFCC99", "AD": "#F4CCCC"}

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
CATEGORICAL_ENCODE = ["Sex", "Handedness", "Marital_Status", "Retired", "Language", "Ethnicity",
                      "APOE4_Status", "APOE_Alleles"]


def prepare_features(df, fit_encoders=None, feature_cols=None):
    df = df.copy().drop(columns=[c for c in DROP_ALWAYS if c in df.columns])
    encoders = fit_encoders or {}
    for col in CATEGORICAL_ENCODE:
        if col not in df.columns:
            continue
        df[col] = df[col].astype(str).str.strip().replace({"NaN": np.nan, "nan": np.nan})
        if fit_encoders is None:
            le = LabelEncoder(); le.fit(df[col].dropna().astype(str)); encoders[col] = le
        le = encoders[col]
        df[col] = df[col].apply(lambda v: le.transform([str(v)])[0]
                                if pd.notna(v) and str(v) in le.classes_ else np.nan)
    df = df.apply(pd.to_numeric, errors="coerce")
    if feature_cols is None:
        feature_cols = list(df.columns)
    for col in feature_cols:
        if col not in df.columns:
            df[col] = np.nan
    return df[feature_cols], encoders, feature_cols


def make_xgb():
    return Pipeline([("imp", SimpleImputer(strategy="median")),
                     ("scl", StandardScaler()),
                     ("clf", XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.05,
                                           use_label_encoder=False, eval_metric="logloss",
                                           random_state=42, verbosity=0, scale_pos_weight=1))])


def load_task_seed(seed, cfg):
    def prep(split):
        df = pd.read_csv(SPLIT_DIR / f"seed_{seed}" / f"{split}.csv", low_memory=False)
        df[cfg["label_col"]] = df[cfg["label_col"]].map(cfg["label_map"])
        return df.dropna(subset=[cfg["label_col"]])
    tr, te = prep("train"), prep("test")
    Xtr, enc, feat = prepare_features(tr)
    Xte, _, _ = prepare_features(te, fit_encoders=enc, feature_cols=feat)
    return Xtr, tr[cfg["label_col"]].astype(int).to_numpy(), Xte, te[cfg["label_col"]].astype(int).to_numpy(), feat


def shap_on(pipe, X):
    """SHAP on imputed+scaled X. Uses the pipeline's SURVIVING feature names (SimpleImputer drops
    all-NaN columns, so the model has <= len(feat) inputs)."""
    names = list(pipe[:-1].get_feature_names_out())
    Xt = pd.DataFrame(pipe[:-1].transform(X), columns=names, index=X.index)
    return shap.TreeExplainer(pipe.named_steps["clf"])(Xt), Xt


def _abs_importance(values):
    a = np.abs(values)
    return a.mean(axis=0).sum(axis=1) if a.ndim == 3 else a.mean(axis=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed_for_plots", type=int, default=0)
    ap.add_argument("--topk", type=int, default=6)
    ap.add_argument("--out_dir", type=str, default=str(OUT_BASE),
                    help="output root (default: D: derivatives/explainability; off the repo)")
    args = ap.parse_args()
    out_base = Path(args.out_dir)

    for sub, cfg in TASKS.items():
        out = out_base / sub; out.mkdir(parents=True, exist_ok=True)
        print(f"\n{'='*70}\n  {cfg['title']}{CHANCE_NOTE.get(sub, '')}\n{'='*70}")

        imp_per_seed, feat, plot_pack = [], None, None
        for seed in SEEDS:
            Xtr, ytr, Xte, yte, feat = load_task_seed(seed, cfg)
            pipe = make_xgb().fit(Xtr, ytr)
            sv, Xt = shap_on(pipe, Xte)
            imp_per_seed.append(pd.Series(_abs_importance(sv.values), index=Xt.columns))
            if seed == args.seed_for_plots:
                plot_pack = (pipe, sv, Xt, Xte)
            print(f"  seed {seed}: n_test={len(yte)} classes={np.bincount(yte)}")
        # align per-seed importances by feature name (all-NaN drops can differ across seeds)
        combo = pd.concat(imp_per_seed, axis=1).reindex(feat).fillna(0.0)
        gi = pd.DataFrame({"feature": feat, "mean_abs_shap": combo.mean(1).values,
                           "std_abs_shap": combo.std(1, ddof=1).values}) \
            .sort_values("mean_abs_shap", ascending=False).reset_index(drop=True)
        gi.to_csv(out / "shap_global_importance.csv", index=False)

        pipe, sv, Xt, Xraw = plot_pack
        feat_used = list(Xt.columns)
        topk = [f for f in gi["feature"].tolist() if f in feat_used][:args.topk]
        cap = f"Best baseline = XGBoost (60 features){CHANCE_NOTE.get(sub, '')}"

        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        top = gi.head(15).iloc[::-1]
        ax.barh(top["feature"], top["mean_abs_shap"], xerr=top["std_abs_shap"],
                color="#4576b5", error_kw=dict(ecolor="#888", lw=0.8))
        ax.set_xlabel("mean |SHAP| (mean +/- std over seeds 0/1/2)")
        ax.set_title(f"{cfg['title']} - XGBoost global feature importance\n{cap}", fontsize=9.5)
        fig.tight_layout(); fig.savefig(out / "shap_bar.png", dpi=200); fig.savefig(out / "shap_bar.pdf"); plt.close(fig)

        if sv.values.ndim == 3:
            for ci, cname in enumerate(cfg["classes"]):
                plt.figure()
                shap.summary_plot(sv.values[:, :, ci], Xt, feature_names=feat_used, show=False, max_display=15)
                plt.title(f"{cfg['title']} - XGBoost SHAP beeswarm (class {cname}, seed {args.seed_for_plots})",
                          fontsize=9)
                plt.tight_layout(); plt.savefig(out / f"shap_beeswarm_{cname}.png", dpi=200); plt.close()
            # all-classes stacked bar — colour each class by the harmonised palette (map by NAME), and keep
            # the legend in CN/MCI/AD order (class_inds="original") to match the attention figure.
            color_ok = all(c in CLASS_HEX for c in cfg["classes"])
            shap.summary_plot(sv.values, Xt, feature_names=feat_used, show=False, max_display=15,
                              class_names=cfg["classes"], plot_type="bar",
                              **(dict(class_inds="original", color=lambda i: CLASS_HEX[cfg["classes"][i]])
                                 if color_ok else {}))
            if color_ok:                                  # thin edge so the pale pastels read on white
                for p in plt.gca().patches:
                    p.set_edgecolor("#444444"); p.set_linewidth(0.4)
            plt.title(f"{cfg['title']} - XGBoost SHAP (all classes)", fontsize=9)
            plt.tight_layout(); plt.savefig(out / "shap_beeswarm.png", dpi=200); plt.close()
        else:
            plt.figure()
            shap.summary_plot(sv.values, Xt, feature_names=feat_used, show=False, max_display=15)
            plt.title(f"{cfg['title']} - XGBoost SHAP beeswarm (seed {args.seed_for_plots})\n{cap}", fontsize=9)
            plt.tight_layout(); plt.savefig(out / "shap_beeswarm.png", dpi=200); plt.close()

        if sv.values.ndim == 2:
            pd.DataFrame(sv.values, columns=feat_used).to_csv(out / "shap_values_seed0.csv", index=False)
        else:
            np.save(out / "shap_values_seed0.npy", sv.values)

        targets = list(range(len(cfg["classes"]))) if cfg["task_type"] == "multiclass" else [None]
        for tgt in targets:
            feats_1d = topk[:args.topk]
            pair = [(topk[0], topk[1])] if len(topk) >= 2 else []
            n = len(feats_1d) + len(pair)
            ncol = 3; nrow = int(np.ceil(n / ncol))
            fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.2 * nrow))
            axes = np.atleast_1d(axes).ravel()
            kw = dict(target=tgt) if tgt is not None else {}
            try:
                PartialDependenceDisplay.from_estimator(pipe, Xraw, feats_1d + pair, ax=axes[:n], **kw)
            except Exception as ex:
                print(f"    [PDP warn] {sub} tgt={tgt}: {str(ex)[:90]}")
            for a in axes[n:]:
                a.axis("off")
            ttl = cfg["title"] + (f" - PDP (class {cfg['classes'][tgt]})" if tgt is not None else " - PDP")
            fig.suptitle(f"{ttl}\ntop-{args.topk} features by mean|SHAP| - {cap}", fontsize=9.5)
            fig.tight_layout(rect=[0, 0, 1, 0.96])
            suff = f"_{cfg['classes'][tgt]}" if tgt is not None else ""
            fig.savefig(out / f"pdp_grid{suff}.png", dpi=200); fig.savefig(out / f"pdp_grid{suff}.pdf"); plt.close(fig)

        print(f"  top-{args.topk} features: {topk}\n  -> {out}")

    print("\nDone. SHAP + PDP written to explainability/<task>/.")


if __name__ == "__main__":
    main()
