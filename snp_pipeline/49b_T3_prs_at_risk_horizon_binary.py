"""
49b_T3_prs_at_risk_horizon_binary.py
====================================
**T3-style binary horizon classification on the whole CN+MCI at-risk
population**, for thresholds N ∈ {3, 5, 7, 10} years.

Cohort + label rule (per-N right-censoring; matches the survival-derived
formulation locked 2026-06-02):
  - Drop `conversion_group ∈ {AD_bl, Excluded}` → base cohort n ≈ 536.
  - For each N:
      y = 1   if `ever_conversion_AD == 1 AND years_to_AD ≤ N`
      y = 0   if `(ever_conversion_AD == 0 AND FU_years ≥ N)`
              OR `(ever_conversion_AD == 1 AND years_to_AD > N)`
      drop    if `ever_conversion_AD == 0 AND FU_years < N`  (right-censored)

- Splits root: `D:/ADNI_BIDS_project/derivatives/clinical/
  no_cdr_stratified_post_exclusion/tabular/baseline/seed_{0,1,2}/{train,val,test}.csv`
  (SAME seeds as scripts 45/46/49a — distinct from baseline_T4/.)
- Sources (4): `meta_prs_EN_combined`, `Kosteridis`, `Kosteridis_MTAG_AD`,
  `Kosteridis_shared_AD_CV`.

Outputs:
  outputs/strict_qc_prs/new_tasks/horizon_binary/<ld_config>/leaderboard.{csv,png}
  outputs/strict_qc_prs/new_tasks/horizon_binary/<ld_config>/per_patient/
    {source}__{covar_mode}__{task}__seed{s}.parquet
  outputs/strict_qc_prs/multimodal_anchors/prs/horizon_{3,5,7,10}y/
    {source}__seed{s}.parquet
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
from sklearn.metrics import (roc_auc_score, average_precision_score,
                              balanced_accuracy_score, confusion_matrix,
                              precision_score, recall_score, f1_score)
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
                       run_fm_head_binary, save_per_task_anchor,
                       task_conditioned_variants_for)  # noqa: E402
from _wandb_io import init_wandb, safe_log, safe_finish  # noqa: E402

SPLITS_ROOT = Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                    "no_cdr_stratified_post_exclusion/tabular/baseline")
# Per-task outputs live at outputs/new_tasks/<task>/<ld>/ ; multimodal
# anchors live at outputs/multimodal_anchors/.  Both moved out from under
# the legacy outputs/strict_qc_prs/ tree per user request 2026-06-02.
NEW_TASKS_ROOT = Path("D:/ADNI_SNP_Omni2.5M_20140220/outputs/new_tasks")
ANCHOR_ROOT    = Path("D:/ADNI_SNP_Omni2.5M_20140220/outputs/multimodal_anchors")
CONV_TSV    = Path("D:/ADNI_SNP_Omni2.5M_20140220/conversion_labels/"
                    "conversion_labels.tsv")
SEEDS = (0, 1, 2)
HORIZONS = (3, 5, 7, 10)
COVAR_MODES = {
    "prs_only":                   [],
    "prs+age+sex+apoe4":          ["age_at_baseline", "Sex_M", "APOE4_Dosage"],
    "prs+age+sex+apoe4+apoe2":    ["age_at_baseline", "Sex_M", "APOE4_Dosage", "APOE2_Dosage"],
}
SOURCES = ("meta_prs_EN_combined", "Kosteridis",
            "Kosteridis_MTAG_AD", "Kosteridis_shared_AD_CV")
META_SOURCE = "meta_prs_EN_combined"

CONV_DROP = {"AD_bl", "Excluded"}


def _load_conversion_meta() -> pd.DataFrame:
    df = pd.read_csv(CONV_TSV, sep="\t")
    df["Patient_ID"] = df["Patient_ID"].astype(str)
    keep = ["Patient_ID", "conversion_group",
             "ever_conversion_AD", "years_to_AD", "FU_years"]
    return df[keep].copy()


def _load_pids(seed: int, split: str) -> pd.DataFrame:
    """PIDs from baseline split csv (just the Patient_ID column)."""
    df = pd.read_csv(SPLITS_ROOT / f"seed_{seed}/{split}.csv", dtype=str)
    df["Patient_ID"] = df["Patient_ID"].astype(str)
    return df[["Patient_ID"]]


def _label_for_N(meta_row: pd.Series, N: int):
    """Return label (0/1/None) per the T3 survival-derived rule."""
    ev = pd.to_numeric(meta_row.get("ever_conversion_AD"), errors="coerce")
    yta = pd.to_numeric(meta_row.get("years_to_AD"), errors="coerce")
    fu  = pd.to_numeric(meta_row.get("FU_years"), errors="coerce")
    if pd.notna(ev) and int(ev) == 1:
        if pd.isna(yta): return None
        return 1 if yta <= N else 0
    # ever_conversion_AD == 0
    if pd.isna(fu): return None
    return 0 if fu >= N else None   # else: right-censored, drop


def _build_label_table(conv: pd.DataFrame, N: int) -> pd.DataFrame:
    """Apply cohort filter + per-N label rule. Returns
    DataFrame[Patient_ID, y] with y ∈ {0, 1}."""
    df = conv[~conv["conversion_group"].isin(CONV_DROP)].copy()
    df["y"] = df.apply(lambda r: _label_for_N(r, N), axis=1)
    df = df[df["y"].notna()].copy()
    df["y"] = df["y"].astype(int)
    return df[["Patient_ID", "y"]]


def _metrics(y_true, y_score) -> dict:
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score).astype(float)
    pred = (y_score > 0.5).astype(int)
    try:    auc = roc_auc_score(y_true, y_score)
    except: auc = float("nan")
    try:    pr_auc = average_precision_score(y_true, y_score)
    except: pr_auc = float("nan")
    bacc = balanced_accuracy_score(y_true, pred)
    sens = recall_score(y_true, pred, zero_division=0)
    spec = recall_score(y_true, pred, pos_label=0, zero_division=0)
    prec = precision_score(y_true, pred, zero_division=0)
    f1   = f1_score(y_true, pred, zero_division=0)
    cm = confusion_matrix(y_true, pred, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel() if cm.size == 4 else (0, 0, 0, 0)
    return {"auc": auc, "bacc": bacc, "f1": f1, "sens": sens, "spec": spec,
             "prec": prec, "pr_auc": pr_auc,
             "tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp),
             "n": int(len(y_true)),
             "pos": int((y_true == 1).sum()), "neg": int((y_true == 0).sum())}


def _build_prs_base(src: str, prs_full: pd.DataFrame,
                      parts: dict, seed: int) -> pd.DataFrame | None:
    if src == META_SOURCE:
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
        _clf, en_prs = _fit_en_meta_prs(sp_df, parts, seed)
        if en_prs.isna().all(): return None
        return pd.DataFrame({"Patient_ID": en_prs.index.astype(str),
                              "PRS": en_prs.values})
    prs_col = f"PRS_{src}"
    if prs_col not in prs_full.columns or prs_full[prs_col].isna().all():
        return None
    return prs_full[["PTID", prs_col]].rename(
        columns={"PTID": "Patient_ID", prs_col: "PRS"})


def _run_one(src: str, prs_full: pd.DataFrame, cov: pd.DataFrame,
              conv: pd.DataFrame, seed: int, mode: str, covars: list,
              N: int, beta_source: str, ld_config: str,
              capture_rows: list) -> dict | None:
    # Build (split -> [Patient_ID, y]) by joining baseline PIDs with the
    # per-N label table.
    label_df = _build_label_table(conv, N)
    parts = {}
    for sp in ("train", "val", "test"):
        pids = _load_pids(seed, sp)
        parts[sp] = pids.merge(label_df, on="Patient_ID", how="inner")
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
    if len(tr_in) < 20 or tr_in["y"].nunique() < 2: return None
    pipe = Pipeline([("imp", SimpleImputer(strategy="median")),
                      ("clf", LogisticRegression(max_iter=2000,
                                                   class_weight="balanced",
                                                   random_state=seed))])
    pipe.fit(tr_in[Xcols].astype(float), tr_in["y"].astype(int))
    train_p = pipe.predict_proba(tr_in[Xcols].astype(float))[:, 1]
    val_p   = (pipe.predict_proba(val_in[Xcols].astype(float))[:, 1]
                if len(val_in) else np.array([]))
    test_p  = (pipe.predict_proba(test_in[Xcols].astype(float))[:, 1]
                if len(test_in) else np.array([]))
    val_m  = (_metrics(val_in["y"], val_p) if len(val_in)
                else _metrics(np.zeros(0), np.zeros(0)))
    test_m = (_metrics(test_in["y"], test_p) if len(test_in)
                else _metrics(np.zeros(0), np.zeros(0)))

    task = f"horizon_{N}y"
    capture_rows.append(capture_predictions(
        fits=fits, tr_in=tr_in, val_in=val_in, test_in=test_in,
        train_p=train_p, val_p=val_p, test_p=test_p,
        source=src, seed=seed, covar_mode=mode,
        beta_source=beta_source, ld_config=ld_config,
        task=task))

    out = {"task": task, "source": src, "covar_mode": mode, "seed": seed}
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
    out_dir = (NEW_TASKS_ROOT / "horizon_binary" / args.ld_config
                if args.beta_source == "raw"
                else NEW_TASKS_ROOT / f"horizon_binary__{args.beta_source}")
    out_dir.mkdir(parents=True, exist_ok=True)
    ld_tag = args.ld_config if args.beta_source == "raw" else "na"

    print(f"[49b T3 binary horizon] LD={args.ld_config} beta_source={args.beta_source}")
    prs_full, _snps_per_src = per_source_prs_table(ld_config=args.ld_config,
                                                       beta_source=args.beta_source)
    cov  = load_subject_covariates()
    conv = _load_conversion_meta()

    # Per-N cohort-size diagnostics
    print("Cohort sizes per N (per-seed × per-fold):")
    for N in HORIZONS:
        label_df = _build_label_table(conv, N)
        print(f"  N={N}y: n={len(label_df):4d}  pos={int(label_df['y'].sum()):3d}  "
               f"neg={int((label_df['y']==0).sum()):3d}")
        for seed in SEEDS:
            for sp in ("train","val","test"):
                pids = _load_pids(seed, sp)
                merged = pids.merge(label_df, on="Patient_ID", how="inner")
                n_pos = int(merged["y"].sum())
                print(f"        seed={seed} {sp:5s}: n={len(merged):3d}  pos={n_pos}  neg={len(merged)-n_pos}")

    wandb_run = init_wandb(args.wandb_project, args.wandb_entity,
                             name=f"T3_horizon__{args.ld_config}__{args.beta_source}",
                             config=vars(args))

    rows, capture_rows = [], []
    for N in HORIZONS:
        for src in SOURCES:
            for mode, covars in COVAR_MODES.items():
                for seed in SEEDS:
                    r = _run_one(src, prs_full, cov, conv, seed, mode, covars,
                                  N, args.beta_source, ld_tag, capture_rows)
                    if r is not None:
                        rows.append(r)
                        safe_log(wandb_run, {"arm": "PRS", **{k: v for k, v in r.items()
                                                                  if not isinstance(v, (dict, list))}})

    # ── FM arms per horizon (Task-1-conditioned attention; SEE DISCLAIMER) ─
    print(f"\n[49b] FM arms -- DISCLAIMER:\n  {FM_ATTENTION_DISCLAIMER}\n")
    for N in HORIZONS:
        task_for_fm = f"horizon_{N}y"
        label_df = _build_label_table(conv, N)
        # Filter task labels to the per-N at-risk cohort once; the FM arms
        # share the same cohort across covar_modes per seed.
        for variant in FM_VARIANTS:
            for mode, covars in COVAR_MODES.items():
                for seed in SEEDS:
                    fold_ptid = {}; fold_y = {}
                    for sp in ("train","val","test"):
                        pids = _load_pids(seed, sp).merge(label_df,
                                                            on="Patient_ID",
                                                            how="inner")
                        fold_ptid[sp] = pids["Patient_ID"].astype(str).tolist()
                        fold_y[sp]    = pids["y"].astype(int).values
                    fm = run_fm_head_binary(
                        variant=variant, seed=seed,
                        fold_ptid=fold_ptid, fold_y=fold_y,
                        cov_df=cov, covars=covars)
                    if fm is None: continue
                    if "val"  not in fm["fold"] or "test" not in fm["fold"]: continue
                    val_m  = _metrics(fm["fold"]["val"]["y"],  fm["fold"]["val"]["proba"])
                    test_m = _metrics(fm["fold"]["test"]["y"], fm["fold"]["test"]["proba"])
                    row = {"task": task_for_fm,
                            "source": variant["short_name"],
                            "covar_mode": mode, "seed": seed}
                    for k, v in val_m.items():  row[f"val_{k}"]  = v
                    for k, v in test_m.items(): row[f"test_{k}"] = v
                    rows.append(row)
                    safe_log(wandb_run, {"arm": "FM",
                                          "attention_source": variant.get("attention_source"),
                                          **{k: v for k, v in row.items()
                                              if not isinstance(v, (dict, list))}})
                    save_per_task_anchor(variant=variant, task_name=task_for_fm,
                                           seed=seed, covar_mode=mode,
                                           fold_capture=fm["fold"])

    # ── Wave 2 FM arms per horizon: task-conditioned attention ──────────
    # One variant set per horizon (sweep 2026-06-02 trained a separate
    # ChromHierPool on each horizon's labels). Loaded from the local
    # task-conditioned sync via load_embedding_for_variant() inside
    # run_fm_head_binary().
    for N in HORIZONS:
        task_for_fm = f"horizon_{N}y"
        wave2 = task_conditioned_variants_for(task_for_fm)
        if not wave2: continue
        label_df = _build_label_table(conv, N)
        for variant in wave2:
            for mode, covars in COVAR_MODES.items():
                for seed in SEEDS:
                    fold_ptid = {}; fold_y = {}
                    for sp in ("train","val","test"):
                        pids = _load_pids(seed, sp).merge(label_df,
                                                            on="Patient_ID",
                                                            how="inner")
                        fold_ptid[sp] = pids["Patient_ID"].astype(str).tolist()
                        fold_y[sp]    = pids["y"].astype(int).values
                    fm = run_fm_head_binary(
                        variant=variant, seed=seed,
                        fold_ptid=fold_ptid, fold_y=fold_y,
                        cov_df=cov, covars=covars)
                    if fm is None: continue
                    if "val"  not in fm["fold"] or "test" not in fm["fold"]: continue
                    val_m  = _metrics(fm["fold"]["val"]["y"],  fm["fold"]["val"]["proba"])
                    test_m = _metrics(fm["fold"]["test"]["y"], fm["fold"]["test"]["proba"])
                    row = {"task": task_for_fm,
                            "source": variant["short_name"],
                            "covar_mode": mode, "seed": seed}
                    for k, v in val_m.items():  row[f"val_{k}"]  = v
                    for k, v in test_m.items(): row[f"test_{k}"] = v
                    rows.append(row)
                    safe_log(wandb_run, {"arm": "FM_wave2",
                                          "attention_source": variant.get("attention_source"),
                                          **{k: v for k, v in row.items()
                                              if not isinstance(v, (dict, list))}})
                    save_per_task_anchor(variant=variant, task_name=task_for_fm,
                                           seed=seed, covar_mode=mode,
                                           fold_capture=fm["fold"])

    df = pd.DataFrame(rows)
    per_p = out_dir / "per_run_metrics.tsv"
    df.to_csv(per_p, sep="\t", index=False)
    print(f"wrote per-run metrics: {per_p}")

    # Pin task ordering for both the leaderboard and the per-seed footer.
    task_order = ["horizon_3y", "horizon_5y", "horizon_7y", "horizon_10y"]
    df["task"] = pd.Categorical(df["task"], categories=task_order, ordered=True)
    metric_cols = ["val_auc","val_bacc","val_f1","val_sens","val_spec",
                    "val_prec","val_pr_auc",
                    "test_auc","test_bacc","test_f1"]
    g = df.groupby(["task","source","covar_mode"]).agg(
        n_seeds=("seed","count"),
        val_n=("val_n","mean"),  val_pos=("val_pos","mean"), val_neg=("val_neg","mean"),
        test_n=("test_n","mean"), test_pos=("test_pos","mean"), test_neg=("test_neg","mean"),
        **{f"{c}_mean": (c, "mean") for c in metric_cols},
        **{f"{c}_std":  (c, "std")  for c in metric_cols},
    ).reset_index().sort_values(["task","covar_mode","val_bacc_mean"],
                                  ascending=[True, True, False])
    lb_p = out_dir / "leaderboard.csv"
    g.to_csv(lb_p, index=False)
    print(f"wrote leaderboard: {lb_p}")

    # Compact PNG: one row per (task, source, covar_mode)
    def _np(n, p):
        return f"{n:.1f} ({p:.1f})"
    show = pd.DataFrame({
        "task":       g["task"],
        "source":     g["source"],
        "covar_mode": g["covar_mode"],
        "val_n (pos)":  [_np(n, p) for n, p in zip(g["val_n"],  g["val_pos"])],
        "test_n (pos)": [_np(n, p) for n, p in zip(g["test_n"], g["test_pos"])],
        "val_bacc":   [f"{m:.3f}+/-{s:.3f}" for m,s in zip(g["val_bacc_mean"], g["val_bacc_std"])],
        "val_auc":    [f"{m:.3f}+/-{s:.3f}" for m,s in zip(g["val_auc_mean"],  g["val_auc_std"])],
        "test_bacc":  [f"{m:.3f}+/-{s:.3f}" for m,s in zip(g["test_bacc_mean"], g["test_bacc_std"])],
        "test_auc":   [f"{m:.3f}+/-{s:.3f}" for m,s in zip(g["test_auc_mean"],  g["test_auc_std"])],
    })
    fig_h = max(3.5, 0.22 * len(show) + 2.5)
    fig, ax = plt.subplots(figsize=(17, fig_h)); ax.axis("off")
    ax.set_title(f"T3 horizon -- at-risk CN+MCI cohort   "
                  f"LD={args.ld_config}  beta={args.beta_source}  "
                  f"(mean +/- std over 3 seeds)",
                  fontsize=11, fontweight="bold", pad=10, loc="left")
    col_widths = [0.07, 0.20, 0.16, 0.08, 0.08, 0.10, 0.10, 0.10, 0.10]
    tbl = ax.table(cellText=show.values.tolist(), colLabels=show.columns.tolist(),
                    loc="center", cellLoc="center", colWidths=col_widths)
    tbl.auto_set_font_size(False); tbl.set_fontsize(7); tbl.scale(1, 1.15)
    for j in range(len(show.columns)):
        c = tbl[(0, j)]; c.set_facecolor("#222"); c.set_text_props(color="white", weight="bold")
    # Bold the top cell per (horizon, metric column) — best performance per
    # horizon. Tie-aware at 3 dp (display precision); bold every row that ties
    # the max within its horizon block.
    METRIC_COLS = {
        "val_bacc":  ("val_bacc_mean",  show.columns.get_loc("val_bacc")),
        "val_auc":   ("val_auc_mean",   show.columns.get_loc("val_auc")),
        "test_bacc": ("test_bacc_mean", show.columns.get_loc("test_bacc")),
        "test_auc":  ("test_auc_mean",  show.columns.get_loc("test_auc")),
    }
    g_indexed = g.reset_index(drop=True)
    for task_name in task_order:
        block_rows = g_indexed[g_indexed["task"] == task_name].index.tolist()
        if not block_rows: continue
        for _key, (mean_col, col_pos) in METRIC_COLS.items():
            vals = g_indexed.loc[block_rows, mean_col]
            top = round(float(vals.max()), 3)
            for row_idx in block_rows:
                if round(float(g_indexed.loc[row_idx, mean_col]), 3) == top:
                    # +1 because table row 0 is the header row.
                    tbl[(row_idx + 1, col_pos)].set_text_props(weight="bold")
    # Per-seed cohort breakdown by horizon (counts identical across source x
    # covar_mode for the same task+seed, so use the first PRS source — FM
    # rows may have a slightly different cohort if some PTIDs are missing
    # from the FM anchor, but the PRS source is the reference here).
    src0 = next((s for s in df["source"].unique()
                  if "BMFM" not in s), df["source"].iloc[0])
    mode0 = df["covar_mode"].iloc[0]
    per_seed = df[(df["source"] == src0) & (df["covar_mode"] == mode0)]\
                 .sort_values(["task", "seed"])
    footer_lines = [
        ("val_n / test_n columns above are means over 3 seeds; per-N cohort "
         "varies by seed because the 80/10/10 split assigns subjects "
         "differently and the right-censored drop is N-specific."),
        "Per-seed cohort sizes (val / test, with `pos` = AD converters within N):",
    ]
    for task_name in task_order:
        rows = per_seed[per_seed["task"] == task_name]
        if rows.empty: continue
        parts = []
        for _, r in rows.iterrows():
            parts.append(
                f"s{int(r['seed'])}: val={int(r['val_n'])}({int(r['val_pos'])}) "
                f"test={int(r['test_n'])}({int(r['test_pos'])})")
        footer_lines.append(f"  {task_name:<12s}  " + "    ".join(parts))
    footer_lines.append("")
    footer_lines.append("FM attention provenance:")
    import textwrap as _tw
    for ln in _tw.wrap(FM_ATTENTION_DISCLAIMER, width=170):
        footer_lines.append(f"  {ln}")
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
        for (src, mode, task, seed), grp in all_rows.groupby(
                ["model","covar_mode","task","seed"], sort=False):
            save_per_patient_parquet(grp, pp_dir,
                                       source=src, seed=int(seed),
                                       covar_mode=mode, task=task)
        # One anchor parquet per (source, seed, task): downstream fusion code
        # selects by horizon (task) independently.
        for (src, seed, task), grp in all_rows.groupby(
                ["model","seed","task"], sort=False):
            anchor_dir = (ANCHOR_ROOT / "prs" / task
                           if args.beta_source == "raw"
                           else ANCHOR_ROOT / f"prs_{args.beta_source}" / task)
            concat_and_save_anchor([grp], anchor_dir,
                                     source=src, seed=int(seed))
        n_pids = all_rows["Patient_ID"].nunique()
        print(f"wrote per-patient parquets: {pp_dir} "
              f"({all_rows.shape[0]:,} rows, {n_pids} unique PTIDs)")

    safe_finish(wandb_run)


if __name__ == "__main__":
    main()
