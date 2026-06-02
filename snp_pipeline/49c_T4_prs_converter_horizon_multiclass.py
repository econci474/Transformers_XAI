"""
49c_T4_prs_converter_horizon_multiclass.py
==========================================
**T4-style 3-class horizon classification on the converter-only cohort**
(pMCI ∪ pCN_to_AD, n = 146). Labels per `years_to_AD`:
  0 → <3 years    1 → 3-7 years    2 → >=7 years

- Splits root: `D:/ADNI_BIDS_project/derivatives/clinical/
  no_cdr_stratified_post_exclusion/tabular/baseline_T4/seed_{0,1,2}/{train,val,test}.csv`
  **(DIFFERENT seeds from baseline/ — stratified on Label_T4.)**
- Labels are precomputed in the T4 split CSVs as `Label_T4 ∈ {0, 1, 2}`
  (built by `clinical_pipeline/01e_build_T4_labels_and_splits.py`).
- Sources (4): `meta_prs_EN_combined`, `Kosteridis`, `Kosteridis_MTAG_AD`,
  `Kosteridis_shared_AD_CV`.
- 3-class softmax LR (multi_class='multinomial', solver='lbfgs').

Outputs:
  outputs/strict_qc_prs/new_tasks/horizon_T4_multiclass/<ld_config>/
    leaderboard.{csv,png}, per_patient/{source}__{covar_mode}__seed{s}.parquet
  outputs/strict_qc_prs/multimodal_anchors/prs/T4_3class/
    {source}__seed{s}.parquet                          # all covar_modes
"""
from __future__ import annotations
import argparse
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (roc_auc_score, balanced_accuracy_score,
                              confusion_matrix, precision_score,
                              recall_score, f1_score)
from sklearn.pipeline import Pipeline

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _prs_strict_qc_lib import (per_source_prs_table,
                                  load_subject_covariates,
                                  DEFAULT_LD_CONFIG)  # noqa: E402
import importlib
mod_45 = importlib.import_module("45_prs_classification_strictQC")
_fit_en_meta_prs = mod_45._fit_en_meta_prs
from _per_patient_io import (capture_predictions, save_per_patient_parquet,
                              concat_and_save_anchor)  # noqa: E402
from _fm_arm import (FM_VARIANTS, FM_ATTENTION_DISCLAIMER,
                       run_fm_head_multiclass, save_per_task_anchor,
                       task_conditioned_variants_for)  # noqa: E402
from _covars_only_new_tasks import (COVAR_SETS_BINARY,
                                       fit_multiclass as _covars_only_fit_multi)  # noqa: E402
from _wandb_io import init_wandb, safe_log, safe_finish  # noqa: E402

SPLITS_ROOT_T4 = Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                       "no_cdr_stratified_post_exclusion/tabular/baseline_T4")
# Per-task outputs live at outputs/new_tasks/<task>/<ld>/ ; multimodal
# anchors live at outputs/multimodal_anchors/.  Both moved out from under
# the legacy outputs/strict_qc_prs/ tree per user request 2026-06-02.
NEW_TASKS_ROOT = Path("D:/ADNI_SNP_Omni2.5M_20140220/outputs/new_tasks")
ANCHOR_ROOT    = Path("D:/ADNI_SNP_Omni2.5M_20140220/outputs/multimodal_anchors")
SEEDS = (0, 1, 2)
COVAR_MODES = {
    "prs_only":                   [],
    "prs+age+sex+apoe4":          ["age_at_baseline", "Sex_M", "APOE4_Dosage"],
    "prs+age+sex+apoe4+apoe2":    ["age_at_baseline", "Sex_M", "APOE4_Dosage", "APOE2_Dosage"],
}
SOURCES = ("meta_prs_EN_combined", "Kosteridis",
            "Kosteridis_MTAG_AD", "Kosteridis_shared_AD_CV")
META_SOURCE = "meta_prs_EN_combined"
TASK_NAME = "T4_3class"
N_CLASSES = 3


def _load_labels(seed: int, split: str) -> pd.DataFrame:
    df = pd.read_csv(SPLITS_ROOT_T4 / f"seed_{seed}/{split}.csv", dtype=str)
    df["Patient_ID"] = df["Patient_ID"].astype(str)
    df["y"] = pd.to_numeric(df["Label_T4"], errors="coerce").astype("Int64")
    df = df[df["y"].notna()].copy()
    df["y"] = df["y"].astype(int)
    return df[["Patient_ID", "y"]]


def _metrics_multi(y_true, y_proba) -> dict:
    """3-class metrics; y_proba is (N, 3)."""
    y_true = np.asarray(y_true).astype(int)
    y_proba = np.asarray(y_proba)
    if y_proba.ndim == 1: y_proba = y_proba.reshape(-1, 1)
    pred = y_proba.argmax(axis=1)
    bacc = balanced_accuracy_score(y_true, pred)
    try:    auc_macro = roc_auc_score(y_true, y_proba, multi_class="ovr",
                                       labels=list(range(N_CLASSES)))
    except Exception:
        auc_macro = float("nan")
    f1_macro   = f1_score(y_true, pred, average="macro", zero_division=0,
                           labels=list(range(N_CLASSES)))
    prec_macro = precision_score(y_true, pred, average="macro", zero_division=0,
                                   labels=list(range(N_CLASSES)))
    rec_macro  = recall_score(y_true, pred, average="macro", zero_division=0,
                                labels=list(range(N_CLASSES)))
    cm = confusion_matrix(y_true, pred, labels=list(range(N_CLASSES)))
    out = {"bacc_macro": bacc, "auc_macro": auc_macro,
            "f1_macro": f1_macro, "prec_macro": prec_macro,
            "rec_macro": rec_macro,
            "n": int(len(y_true)),
            "n_class0": int((y_true == 0).sum()),
            "n_class1": int((y_true == 1).sum()),
            "n_class2": int((y_true == 2).sum())}
    # per-class confusion-matrix counts (flatten 3x3)
    for i in range(N_CLASSES):
        for j in range(N_CLASSES):
            out[f"cm_{i}{j}"] = int(cm[i, j])
    return out


def _build_prs_base(src: str, prs_full: pd.DataFrame,
                      parts: dict, seed: int) -> pd.DataFrame | None:
    if src == META_SOURCE:
        # meta_prs_EN_combined expects parts with "y" — _fit_en_meta_prs fits
        # logistic regression on the train fold's y. For multiclass, this still
        # produces a useful score (collapses to binary AD-yes/no via class>0?
        # Actually parts["train"]["y"] here is multiclass — _fit_en_meta_prs
        # will still fit and produce a continuous score; it generalizes to
        # multiclass scoring at decision_function level).
        # Safer: collapse to binary 'ever-AD-converter' (always 1 here) — but
        # then nunique=1 and helper bails. Instead pre-binarize: class 0 (early
        # converter, <3y) vs class 2 (late, >=7y) trained as 1 vs 0 — gives a
        # meaningful EN-meta direction. Class 1 sits between, as a tiebreaker.
        bin_parts = {}
        for sp, df_sp in parts.items():
            d = df_sp.copy()
            d = d[d["y"].isin([0, 2])].copy()
            d["y"] = (d["y"] == 0).astype(int)   # 1 = early converter
            bin_parts[sp] = d
        source_cols = [c for c in prs_full.columns if c.startswith("PRS_")
                        and c not in ("PRS_prs_all_dedup",
                                        "PRS_prs_all_dedup_ivw",
                                        "PRS_prs_all_dedup_filtered",
                                        "PRS_Kosteridis_MTAG_AD",
                                        "PRS_Kosteridis_shared_AD_CV",
                                        "PRS_Kosteridis_novel_AD")
                        and not c.endswith("EN_dosage")
                        and not c.endswith("EN_combined")]
        sp_df = prs_full[["PTID"] + source_cols].set_index("PTID")
        sp_df.index = sp_df.index.astype(str)
        _clf, en_prs = _fit_en_meta_prs(sp_df, bin_parts, seed)
        if en_prs.isna().all(): return None
        return pd.DataFrame({"Patient_ID": en_prs.index.astype(str),
                              "PRS": en_prs.values})
    prs_col = f"PRS_{src}"
    if prs_col not in prs_full.columns or prs_full[prs_col].isna().all():
        return None
    return prs_full[["PTID", prs_col]].rename(
        columns={"PTID": "Patient_ID", prs_col: "PRS"})


def _run_one(src: str, prs_full: pd.DataFrame, cov: pd.DataFrame,
              seed: int, mode: str, covars: list,
              beta_source: str, ld_config: str,
              capture_rows: list) -> dict | None:
    parts = {sp: _load_labels(seed, sp) for sp in ("train", "val", "test")}
    base = _build_prs_base(src, prs_full, parts, seed)
    if base is None: return None
    base["Patient_ID"] = base["Patient_ID"].astype(str)
    fits = {}
    for sp in ("train", "val", "test"):
        f = parts[sp].merge(base, on="Patient_ID", how="inner")
        f = f.merge(cov, on="Patient_ID", how="left")
        fits[sp] = f
    tr = fits["train"]
    if tr.empty or tr["y"].nunique() < 2: return None
    mu = tr["PRS"].mean(); sd = tr["PRS"].std(ddof=0) or 1.0
    for sp in fits:
        fits[sp] = fits[sp].copy()
        fits[sp]["PRS_z"] = (fits[sp]["PRS"] - mu) / sd
    Xcols = ["PRS_z"] + covars
    needed = ["y"] + Xcols
    tr_in   = fits["train"][needed].dropna()
    val_in  = fits["val"][needed].dropna()
    test_in = fits["test"][needed].dropna()
    # T4 cohort is tiny (~146 split 80/10/10 ≈ 117/14/15). Tolerate small folds.
    if len(tr_in) < 20 or tr_in["y"].nunique() < 2: return None
    # sklearn 1.5+ removed the `multi_class` kwarg; lbfgs handles multinomial
    # by default when n_classes > 2.
    pipe = Pipeline([("imp", SimpleImputer(strategy="median")),
                      ("clf", LogisticRegression(max_iter=4000,
                                                   solver="lbfgs",
                                                   class_weight="balanced",
                                                   random_state=seed))])
    pipe.fit(tr_in[Xcols].astype(float), tr_in["y"].astype(int))
    classes_ = pipe.named_steps["clf"].classes_  # sorted unique classes
    n_seen   = len(classes_)
    def _full_proba(p_arr):
        """Pad missing classes with 0.0 so output is always (N, 3)."""
        if n_seen == N_CLASSES: return p_arr
        full = np.zeros((p_arr.shape[0], N_CLASSES), dtype=float)
        for col_idx, c in enumerate(classes_):
            full[:, int(c)] = p_arr[:, col_idx]
        return full
    train_p_full = _full_proba(pipe.predict_proba(tr_in[Xcols].astype(float)))
    val_p_full   = (_full_proba(pipe.predict_proba(val_in[Xcols].astype(float)))
                     if len(val_in) else np.zeros((0, N_CLASSES)))
    test_p_full  = (_full_proba(pipe.predict_proba(test_in[Xcols].astype(float)))
                     if len(test_in) else np.zeros((0, N_CLASSES)))
    val_m  = (_metrics_multi(val_in["y"], val_p_full) if len(val_in)
                else {"bacc_macro": float("nan"), "auc_macro": float("nan"),
                       "f1_macro": float("nan"), "prec_macro": float("nan"),
                       "rec_macro": float("nan"), "n": 0,
                       "n_class0": 0, "n_class1": 0, "n_class2": 0})
    test_m = (_metrics_multi(test_in["y"], test_p_full) if len(test_in)
                else {"bacc_macro": float("nan"), "auc_macro": float("nan"),
                       "f1_macro": float("nan"), "prec_macro": float("nan"),
                       "rec_macro": float("nan"), "n": 0,
                       "n_class0": 0, "n_class1": 0, "n_class2": 0})

    # For per-patient capture: store proba_class0..2 + a scalar outcome_proba
    # (the predicted class label as float for convenience). Override with class
    # probabilities via extra columns added after capture.
    pred_train = train_p_full.argmax(axis=1).astype(float)
    pred_val   = (val_p_full.argmax(axis=1).astype(float)
                   if len(val_p_full) else np.array([]))
    pred_test  = (test_p_full.argmax(axis=1).astype(float)
                   if len(test_p_full) else np.array([]))
    cap = capture_predictions(
        fits=fits, tr_in=tr_in, val_in=val_in, test_in=test_in,
        train_p=pred_train, val_p=pred_val, test_p=pred_test,
        source=src, seed=seed, covar_mode=mode,
        beta_source=beta_source, ld_config=ld_config,
        task=TASK_NAME)
    # Attach class-probability columns. Order: train, val, test rows.
    proba_stack = np.concatenate(
        [train_p_full] +
        ([val_p_full]  if len(val_p_full)  else []) +
        ([test_p_full] if len(test_p_full) else []),
        axis=0)
    cap["proba_class0"] = proba_stack[:, 0]
    cap["proba_class1"] = proba_stack[:, 1]
    cap["proba_class2"] = proba_stack[:, 2]
    capture_rows.append(cap)

    out = {"source": src, "covar_mode": mode, "seed": seed}
    for k, v in val_m.items():  out[f"val_{k}"]  = v
    for k, v in test_m.items(): out[f"test_{k}"] = v
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ld-config", default=DEFAULT_LD_CONFIG)
    ap.add_argument("--beta-source", default="raw", choices=["raw","prscs"])
    ap.add_argument("--wandb-project", default=None, dest="wandb_project")
    ap.add_argument("--wandb-entity", default=None, dest="wandb_entity")
    args = ap.parse_args()
    out_dir = (NEW_TASKS_ROOT / "horizon_T4_multiclass" / args.ld_config
                if args.beta_source == "raw"
                else NEW_TASKS_ROOT / f"horizon_T4_multiclass__{args.beta_source}")
    out_dir.mkdir(parents=True, exist_ok=True)
    ld_tag = args.ld_config if args.beta_source == "raw" else "na"

    print(f"[49c T4 multiclass] LD={args.ld_config} beta_source={args.beta_source}")
    prs_full, _snps_per_src = per_source_prs_table(ld_config=args.ld_config,
                                                       beta_source=args.beta_source)
    cov = load_subject_covariates()

    print("Cohort sizes (T4 splits — converter-only):")
    for seed in SEEDS:
        for sp in ("train","val","test"):
            d = _load_labels(seed, sp)
            counts = d["y"].value_counts().sort_index().to_dict()
            print(f"  seed={seed} {sp:5s}: n={len(d):3d}  by class={counts}")

    wandb_run = init_wandb(args.wandb_project, args.wandb_entity,
                             name=f"T4_3class__{args.ld_config}__{args.beta_source}",
                             config=vars(args))

    rows, capture_rows = [], []
    for src in SOURCES:
        for mode, covars in COVAR_MODES.items():
            for seed in SEEDS:
                r = _run_one(src, prs_full, cov, seed, mode, covars,
                              args.beta_source, ld_tag, capture_rows)
                if r is not None:
                    rows.append(r)
                    safe_log(wandb_run, {"arm": "PRS", **{k: v for k, v in r.items()
                                                              if not isinstance(v, (dict, list))}})

    # ── Covariates-only baseline arms (multiclass) ───────────────────────
    print(f"\n[49c] Covariates-only baseline arms ({len(COVAR_SETS_BINARY)} sets)")
    for covar_name, covars_list in COVAR_SETS_BINARY.items():
        for seed in SEEDS:
            parts = {sp: _load_labels(seed, sp) for sp in ("train","val","test")}
            fold_ptid = {sp: parts[sp]["Patient_ID"].astype(str).tolist()
                          for sp in ("train","val","test")}
            fold_y = {sp: parts[sp]["y"].astype(int).values
                       for sp in ("train","val","test")}
            r = _covars_only_fit_multi(cov, covars_list, covar_name, seed,
                                          fold_ptid=fold_ptid, fold_y=fold_y,
                                          n_classes=N_CLASSES)
            if r is not None:
                rows.append(r)
                safe_log(wandb_run, {"arm": "covars_only",
                                      **{k: v for k, v in r.items()
                                          if not isinstance(v, (dict, list))}})

    # ── FM arms (multiclass; Task-1-conditioned attention; SEE DISCLAIMER) ─
    print(f"\n[49c] FM arms -- DISCLAIMER:\n  {FM_ATTENTION_DISCLAIMER}\n")
    for variant in FM_VARIANTS:
        for mode, covars in COVAR_MODES.items():
            for seed in SEEDS:
                parts = {sp: _load_labels(seed, sp) for sp in ("train","val","test")}
                fold_ptid = {sp: parts[sp]["Patient_ID"].astype(str).tolist()
                              for sp in ("train","val","test")}
                fold_y = {sp: parts[sp]["y"].astype(int).values
                           for sp in ("train","val","test")}
                fm = run_fm_head_multiclass(
                    variant=variant, seed=seed,
                    fold_ptid=fold_ptid, fold_y=fold_y,
                    cov_df=cov, covars=covars,
                    n_classes=N_CLASSES)
                if fm is None: continue
                if "val"  not in fm["fold"] or "test" not in fm["fold"]: continue
                val_m  = _metrics_multi(fm["fold"]["val"]["y"],  fm["fold"]["val"]["proba"])
                test_m = _metrics_multi(fm["fold"]["test"]["y"], fm["fold"]["test"]["proba"])
                row = {"source": variant["short_name"],
                        "covar_mode": mode, "seed": seed}
                for k, v in val_m.items():  row[f"val_{k}"]  = v
                for k, v in test_m.items(): row[f"test_{k}"] = v
                rows.append(row)
                safe_log(wandb_run, {"arm": "FM",
                                      "attention_source": variant.get("attention_source"),
                                      **{k: v for k, v in row.items()
                                          if not isinstance(v, (dict, list))}})
                save_per_task_anchor(variant=variant, task_name=TASK_NAME,
                                       seed=seed, covar_mode=mode,
                                       fold_capture=fm["fold"])

    # ── Wave 2 FM arms (multiclass; task-conditioned attention) ─────────
    wave2 = task_conditioned_variants_for(TASK_NAME)
    if wave2:
        print(f"\n[49c] Wave 2 FM arms ({len(wave2)} task-conditioned variants)")
        for variant in wave2:
            for mode, covars in COVAR_MODES.items():
                for seed in SEEDS:
                    parts = {sp: _load_labels(seed, sp) for sp in ("train","val","test")}
                    fold_ptid = {sp: parts[sp]["Patient_ID"].astype(str).tolist()
                                  for sp in ("train","val","test")}
                    fold_y = {sp: parts[sp]["y"].astype(int).values
                               for sp in ("train","val","test")}
                    fm = run_fm_head_multiclass(
                        variant=variant, seed=seed,
                        fold_ptid=fold_ptid, fold_y=fold_y,
                        cov_df=cov, covars=covars,
                        n_classes=N_CLASSES)
                    if fm is None: continue
                    if "val"  not in fm["fold"] or "test" not in fm["fold"]: continue
                    val_m  = _metrics_multi(fm["fold"]["val"]["y"],  fm["fold"]["val"]["proba"])
                    test_m = _metrics_multi(fm["fold"]["test"]["y"], fm["fold"]["test"]["proba"])
                    row = {"source": variant["short_name"],
                            "covar_mode": mode, "seed": seed}
                    for k, v in val_m.items():  row[f"val_{k}"]  = v
                    for k, v in test_m.items(): row[f"test_{k}"] = v
                    rows.append(row)
                    safe_log(wandb_run, {"arm": "FM_wave2",
                                          "attention_source": variant.get("attention_source"),
                                          **{k: v for k, v in row.items()
                                              if not isinstance(v, (dict, list))}})
                    save_per_task_anchor(variant=variant, task_name=TASK_NAME,
                                           seed=seed, covar_mode=mode,
                                           fold_capture=fm["fold"])

    df = pd.DataFrame(rows)
    per_p = out_dir / "per_run_metrics.tsv"
    df.to_csv(per_p, sep="\t", index=False)
    print(f"wrote per-run metrics: {per_p}")

    metric_cols = ["val_bacc_macro","val_auc_macro","val_f1_macro",
                    "test_bacc_macro","test_auc_macro","test_f1_macro"]
    g = df.groupby(["source","covar_mode"]).agg(
        n_seeds=("seed","count"),
        val_n=("val_n","mean"),  test_n=("test_n","mean"),
        val_n_class0=("val_n_class0","mean"),
        val_n_class1=("val_n_class1","mean"),
        val_n_class2=("val_n_class2","mean"),
        test_n_class0=("test_n_class0","mean"),
        test_n_class1=("test_n_class1","mean"),
        test_n_class2=("test_n_class2","mean"),
        **{f"{c}_mean": (c, "mean") for c in metric_cols},
        **{f"{c}_std":  (c, "std")  for c in metric_cols},
    ).reset_index().sort_values(["covar_mode","val_bacc_macro_mean"],
                                  ascending=[True, False])
    lb_p = out_dir / "leaderboard.csv"
    g.to_csv(lb_p, index=False)
    print(f"wrote leaderboard: {lb_p}")

    # val_n / test_n columns dropped per user request; the per-seed footer below
    # carries the cohort breakdown.
    show = pd.DataFrame({
        "source":     g["source"],
        "covar_mode": g["covar_mode"],
        "val_bacc":   [f"{m:.3f}+/-{s:.3f}" for m,s in zip(g["val_bacc_macro_mean"], g["val_bacc_macro_std"])],
        "val_auc":    [f"{m:.3f}+/-{s:.3f}" for m,s in zip(g["val_auc_macro_mean"],  g["val_auc_macro_std"])],
        "test_bacc":  [f"{m:.3f}+/-{s:.3f}" for m,s in zip(g["test_bacc_macro_mean"], g["test_bacc_macro_std"])],
        "test_auc":   [f"{m:.3f}+/-{s:.3f}" for m,s in zip(g["test_auc_macro_mean"],  g["test_auc_macro_std"])],
    })
    fig_h = max(3.5, 0.30 * len(show) + 2.5)
    fig, ax = plt.subplots(figsize=(18, fig_h)); ax.axis("off")
    ax.set_title(f"T4 horizon multiclass -- converter cohort   "
                  f"LD={args.ld_config}  beta={args.beta_source}  "
                  f"(macro mean +/- std over 3 seeds)\n"
                  f"FM rows: attention provenance annotated in [..] "
                  f"-- see footer for the labels each attention pool was trained on",
                  fontsize=10, fontweight="bold", pad=10, loc="left")
    col_widths = [0.42, 0.18, 0.10, 0.10, 0.10, 0.10]
    tbl = ax.table(cellText=show.values.tolist(), colLabels=show.columns.tolist(),
                    loc="center", cellLoc="center", colWidths=col_widths)
    tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1, 1.2)
    for j in range(len(show.columns)):
        c = tbl[(0, j)]; c.set_facecolor("#222"); c.set_text_props(color="white", weight="bold")
    # Left-align source column so the variant name reads cleanly.
    src_col = show.columns.get_loc("source")
    for row_idx in range(len(show)):
        cell = tbl[(row_idx + 1, src_col)]
        cell.set_text_props(ha="left")
        cell.PAD = 0.02
    # Bold top-3 means for val/test bacc + auc (tie-aware at 3 dp).
    BOLD_COLS = {
        "val_bacc":  ("val_bacc_macro_mean",  show.columns.get_loc("val_bacc")),
        "test_bacc": ("test_bacc_macro_mean", show.columns.get_loc("test_bacc")),
        "val_auc":   ("val_auc_macro_mean",   show.columns.get_loc("val_auc")),
        "test_auc":  ("test_auc_macro_mean",  show.columns.get_loc("test_auc")),
    }
    g_indexed = g.reset_index(drop=True)
    for _key, (mean_col, col_pos) in BOLD_COLS.items():
        vals_3dp = g_indexed[mean_col].round(3)
        top3 = sorted(vals_3dp.dropna().unique(), reverse=True)[:3]
        if not top3: continue
        cutoff = top3[-1]
        for row_idx in g_indexed.index:
            v = vals_3dp.iloc[row_idx]
            if pd.notna(v) and v >= cutoff:
                tbl[(row_idx + 1, col_pos)].set_text_props(weight="bold")
    # Per-seed cohort breakdown (use PRS source as reference; FM rows may
    # differ in cohort if some PTIDs are missing from the FM anchor).
    src0 = next((s for s in df["source"].unique()
                  if "BMFM" not in s), df["source"].iloc[0])
    mode0 = df["covar_mode"].iloc[0]
    per_seed = df[(df["source"] == src0) & (df["covar_mode"] == mode0)]\
                 .sort_values("seed")
    footer_lines = [
        ("Per-seed val and test cohort sizes (T4 splits are stratified on "
         "Label_T4, so class proportions are stable across seeds):"),
        "  Class counts in brackets: (n_<3y / n_3-7y / n_>=7y)",
    ]
    for _, r in per_seed.iterrows():
        footer_lines.append(
            f"  seed {int(r['seed'])}: "
            f"val n={int(r['val_n'])} "
            f"({int(r['val_n_class0'])}/{int(r['val_n_class1'])}/{int(r['val_n_class2'])})    "
            f"test n={int(r['test_n'])} "
            f"({int(r['test_n_class0'])}/{int(r['test_n_class1'])}/{int(r['test_n_class2'])})")
    footer_lines.append("")
    footer_lines.append("FM attention provenance (the labels each ChromHierPool"
                         " was trained on -- attn=... in row source):")
    footer_lines.append("  [attn=T1_AD_CN]          ChromHierPool trained on the "
                         "ORIGINAL Task 1 labels: AD_bl + pMCI + CN_to_AD (positive) "
                         "vs sCN (negative).  Softly biased toward AD-onset signal.")
    footer_lines.append("  [attn=fixed_aggregator]   Deterministic _fixed_aggregator_v3 "
                         "(gamma=delta=1). No training -- task-independent. The XGB "
                         "head still sees raw 771-d features.")
    footer_lines.append("  [attn=T4_3class]          ChromHierPool RE-TRAINED on this "
                         "task's labels: converters partitioned by years-to-AD bins  "
                         "(class 0: <3 y;  class 1: 3-7 y;  class 2: >=7 y).")
    footer_lines.append("")
    footer_lines.append("Bold cells: top-3 within each of val_bacc / val_auc / "
                         "test_bacc / test_auc (tie-aware at 3 dp).")
    fig.text(0.01, 0.01, "\n".join(footer_lines), fontsize=8,
              family="monospace", va="bottom")
    png = out_dir / "leaderboard.png"
    fig.savefig(png, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"wrote PNG: {png}")
    if wandb_run is not None:
        try:
            import wandb as _wb
            wandb_run.log({"leaderboard": _wb.Table(dataframe=g),
                            "leaderboard_png": _wb.Image(str(png))})
        except Exception as e:
            print(f"[wandb] table/image log skipped ({e})")

    if capture_rows:
        all_rows = pd.concat(capture_rows, ignore_index=True)
        pp_dir = out_dir / "per_patient"
        anchor_dir = (ANCHOR_ROOT / "prs" / TASK_NAME
                       if args.beta_source == "raw"
                       else ANCHOR_ROOT / f"prs_{args.beta_source}" / TASK_NAME)
        for (src, mode, seed), grp in all_rows.groupby(
                ["model","covar_mode","seed"], sort=False):
            save_per_patient_parquet(grp, pp_dir,
                                       source=src, seed=int(seed),
                                       covar_mode=mode)
        for (src, seed), grp in all_rows.groupby(["model","seed"], sort=False):
            concat_and_save_anchor([grp], anchor_dir,
                                     source=src, seed=int(seed))
        n_pids = all_rows["Patient_ID"].nunique()
        print(f"wrote per-patient parquets: {pp_dir} "
              f"({all_rows.shape[0]:,} rows, {n_pids} unique PTIDs)")
        print(f"wrote multimodal anchors:   {anchor_dir}")

    safe_finish(wandb_run)


if __name__ == "__main__":
    main()
