"""
49a_prs_sCN_vs_pCN_classification.py
====================================
T3-style PRS classification: **sCN (neg) vs pCN (pos)**,
where pCN = pCN_to_AD ∪ pCN_to_MCI (strict baseline-CN cohort).

- Splits root: `D:/ADNI_BIDS_project/derivatives/clinical/
  no_cdr_stratified_post_exclusion/tabular/baseline/seed_{0,1,2}/{train,val,test}.csv`
  (SAME seeds as scripts 45/46/49b — distinct from baseline_T4/.)
- Cohort filter at runtime: keep rows where
  `sCN==1` (negative)  OR  `pCN_to_AD==1 OR pCN_to_MCI==1` (positive).
- Sources (4): `meta_prs_EN_combined`, `Kosteridis`, `Kosteridis_MTAG_AD`,
  `Kosteridis_shared_AD_CV` (per multimodal-anchor plan).
- 3 covar modes (prs_only / +cov / +cov_apoe2) × 3 seeds.

Outputs:
  outputs/strict_qc_prs/new_tasks/sCN_vs_pCN/<ld_config>/leaderboard.{csv,png}
  outputs/strict_qc_prs/new_tasks/sCN_vs_pCN/<ld_config>/per_patient/
    {source}__{covar_mode}__seed{s}.parquet
  outputs/strict_qc_prs/multimodal_anchors/prs/sCN_vs_pCN/
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
# Reuse the meta-EN fit from script 45 (zero drift — same training-fold rule).
import importlib
mod_45 = importlib.import_module("45_prs_classification_strictQC")
_fit_en_meta_prs = mod_45._fit_en_meta_prs
from _per_patient_io import (capture_predictions, save_per_patient_parquet,
                              concat_and_save_anchor)  # noqa: E402
from _fm_arm import (FM_VARIANTS, FM_ATTENTION_DISCLAIMER,
                       run_fm_head_binary, save_per_task_anchor,
                       task_conditioned_variants_for)  # noqa: E402
from _covars_only_new_tasks import (COVAR_SETS_BINARY,
                                       fit_binary as _covars_only_fit_binary)  # noqa: E402
from _wandb_io import init_wandb, safe_log, safe_finish  # noqa: E402

SPLITS_ROOT = Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                    "no_cdr_stratified_post_exclusion/tabular/baseline")
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

TASK_NAME = "sCN_vs_pCN"


def _load_labels(seed: int, split: str) -> pd.DataFrame:
    """sCN vs (pCN_to_AD ∪ pCN_to_MCI). Cohort restricted to strict baseline-CN."""
    df = pd.read_csv(SPLITS_ROOT / f"seed_{seed}/{split}.csv", dtype=str)
    df["Patient_ID"] = df["Patient_ID"].astype(str)
    df["sCN"]        = pd.to_numeric(df["sCN"], errors="coerce").fillna(0).astype(int)
    df["pCN_to_AD"]  = pd.to_numeric(df["pCN_to_AD"], errors="coerce").fillna(0).astype(int)
    df["pCN_to_MCI"] = pd.to_numeric(df["pCN_to_MCI"], errors="coerce").fillna(0).astype(int)
    y_pos = ((df["pCN_to_AD"] == 1) | (df["pCN_to_MCI"] == 1)).astype(int)
    y_neg = (df["sCN"] == 1).astype(int)
    df = df[(y_pos == 1) | (y_neg == 1)].copy()
    df["y"] = y_pos.loc[df.index]
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
    """Return DataFrame with columns ['Patient_ID', 'PRS'] for `src`."""
    if src == META_SOURCE:
        # Same exclusion rule as script 45's meta_prs_EN_combined.
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
    tr_in  = fits["train"][needed].dropna()
    val_in = fits["val"][needed].dropna()
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

    capture_rows.append(capture_predictions(
        fits=fits, tr_in=tr_in, val_in=val_in, test_in=test_in,
        train_p=train_p, val_p=val_p, test_p=test_p,
        source=src, seed=seed, covar_mode=mode,
        beta_source=beta_source, ld_config=ld_config,
        task=TASK_NAME))

    out = {"source": src, "covar_mode": mode, "seed": seed}
    for k, v in val_m.items():  out[f"val_{k}"]  = v
    for k, v in test_m.items(): out[f"test_{k}"] = v
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ld-config", default=DEFAULT_LD_CONFIG)
    ap.add_argument("--beta-source", default="raw", choices=["raw","prscs"])
    ap.add_argument("--wandb-project", default=None, dest="wandb_project",
                    help="WandB project name (skip wandb if not set).")
    ap.add_argument("--wandb-entity", default=None, dest="wandb_entity",
                    help="WandB entity (default = ec474-university-of-cambridge).")
    args = ap.parse_args()
    out_dir = (NEW_TASKS_ROOT / TASK_NAME / args.ld_config
                if args.beta_source == "raw"
                else NEW_TASKS_ROOT / f"{TASK_NAME}__{args.beta_source}")
    out_dir.mkdir(parents=True, exist_ok=True)
    ld_tag = args.ld_config if args.beta_source == "raw" else "na"

    print(f"[49a sCN-vs-pCN] LD={args.ld_config} beta_source={args.beta_source}")
    wandb_run = init_wandb(args.wandb_project, args.wandb_entity,
                             name=f"{TASK_NAME}__{args.ld_config}__{args.beta_source}",
                             config=vars(args))

    prs_full, _snps_per_src = per_source_prs_table(ld_config=args.ld_config,
                                                       beta_source=args.beta_source)
    cov = load_subject_covariates()

    # Cohort-size diagnostics (per-seed × per-fold).
    print("Cohort sizes (after sCN-vs-pCN filter):")
    for seed in SEEDS:
        for sp in ("train","val","test"):
            d = _load_labels(seed, sp)
            n_pos = int(d["y"].sum())
            print(f"  seed={seed} {sp:5s}: n={len(d):3d}  pos={n_pos}  neg={len(d)-n_pos}")

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

    # ── Covariates-only baseline arms (no PRS, no FM feature) ───────────
    # Floor reference: what does age+sex+APOE alone get on this cohort?
    print(f"\n[49a] Covariates-only baseline arms ({len(COVAR_SETS_BINARY)} sets)")
    for covar_name, covars_list in COVAR_SETS_BINARY.items():
        for seed in SEEDS:
            parts = {sp: _load_labels(seed, sp) for sp in ("train","val","test")}
            fold_ptid = {sp: parts[sp]["Patient_ID"].astype(str).tolist()
                          for sp in ("train","val","test")}
            fold_y = {sp: parts[sp]["y"].astype(int).values
                       for sp in ("train","val","test")}
            r = _covars_only_fit_binary(cov, covars_list, covar_name, seed,
                                          fold_ptid=fold_ptid, fold_y=fold_y)
            if r is not None:
                rows.append(r)
                safe_log(wandb_run, {"arm": "covars_only",
                                      **{k: v for k, v in r.items()
                                          if not isinstance(v, (dict, list))}})

    # ── FM arms (Task-1-conditioned attention; SEE DISCLAIMER) ───────────
    print(f"\n[49a] FM arms — DISCLAIMER:\n  {FM_ATTENTION_DISCLAIMER}\n")
    for variant in FM_VARIANTS:
        for mode, covars in COVAR_MODES.items():
            for seed in SEEDS:
                parts = {sp: _load_labels(seed, sp) for sp in ("train","val","test")}
                fold_ptid = {sp: parts[sp]["Patient_ID"].astype(str).tolist()
                              for sp in ("train","val","test")}
                fold_y = {sp: parts[sp]["y"].astype(int).values
                           for sp in ("train","val","test")}
                fm = run_fm_head_binary(
                    variant=variant, seed=seed,
                    fold_ptid=fold_ptid, fold_y=fold_y,
                    cov_df=cov, covars=covars)
                if fm is None: continue
                val_m  = (_metrics(fm["fold"]["val"]["y"],   fm["fold"]["val"]["proba"])
                            if "val"  in fm["fold"] else None)
                test_m = (_metrics(fm["fold"]["test"]["y"],  fm["fold"]["test"]["proba"])
                            if "test" in fm["fold"] else None)
                if val_m is None or test_m is None: continue
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

    # ── Wave 2 FM arms: task-conditioned attention (Colab sweep 2026-06-02)
    # Same shape as the Wave 1 loop above, but `variant` carries
    # `variant_subdir` so load_embedding_for_variant() reads from the
    # task-conditioned Drive sync at TASK_CONDITIONED_ROOT rather than the
    # script-50 anchor. Skipped silently per-variant if the sync hasn't
    # landed yet (the file just isn't on disk).
    wave2 = task_conditioned_variants_for(TASK_NAME)
    if wave2:
        print(f"\n[49a] Wave 2 FM arms ({len(wave2)} task-conditioned variants)")
        for variant in wave2:
            for mode, covars in COVAR_MODES.items():
                for seed in SEEDS:
                    parts = {sp: _load_labels(seed, sp) for sp in ("train","val","test")}
                    fold_ptid = {sp: parts[sp]["Patient_ID"].astype(str).tolist()
                                  for sp in ("train","val","test")}
                    fold_y = {sp: parts[sp]["y"].astype(int).values
                               for sp in ("train","val","test")}
                    fm = run_fm_head_binary(
                        variant=variant, seed=seed,
                        fold_ptid=fold_ptid, fold_y=fold_y,
                        cov_df=cov, covars=covars)
                    if fm is None: continue
                    val_m  = (_metrics(fm["fold"]["val"]["y"],   fm["fold"]["val"]["proba"])
                                if "val"  in fm["fold"] else None)
                    test_m = (_metrics(fm["fold"]["test"]["y"],  fm["fold"]["test"]["proba"])
                                if "test" in fm["fold"] else None)
                    if val_m is None or test_m is None: continue
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

    metric_cols = ["val_auc","val_bacc","val_f1","val_sens","val_spec",
                    "val_prec","val_pr_auc",
                    "test_auc","test_bacc","test_f1"]
    g = df.groupby(["source","covar_mode"]).agg(
        n_seeds=("seed","count"),
        val_n=("val_n","mean"),  val_pos=("val_pos","mean"), val_neg=("val_neg","mean"),
        test_n=("test_n","mean"), test_pos=("test_pos","mean"), test_neg=("test_neg","mean"),
        **{f"{c}_mean": (c, "mean") for c in metric_cols},
        **{f"{c}_std":  (c, "std")  for c in metric_cols},
    ).reset_index().sort_values(["covar_mode","val_bacc_mean"], ascending=[True, False])
    lb_p = out_dir / "leaderboard.csv"
    g.to_csv(lb_p, index=False)
    print(f"wrote leaderboard: {lb_p}")

    # Compact PNG (val + test BalAcc/AUC).
    # `val_n N (P)` / `test_n N (P)` format — P = positive class (pCN converters).
    # N and P are the means across 3 seeds (one decimal). Per-seed counts can
    # differ slightly because the 80/10/10 split assigns subjects differently
    # by seed.
    def _np(n, p):
        return f"{n:.1f} ({p:.1f})"
    show = pd.DataFrame({
        "source":     g["source"],
        "covar_mode": g["covar_mode"],
        "val_n (pos)":  [_np(n, p) for n, p in zip(g["val_n"],  g["val_pos"])],
        "test_n (pos)": [_np(n, p) for n, p in zip(g["test_n"], g["test_pos"])],
        "val_bacc":   [f"{m:.3f}+/-{s:.3f}" for m,s in zip(g["val_bacc_mean"], g["val_bacc_std"])],
        "val_auc":    [f"{m:.3f}+/-{s:.3f}" for m,s in zip(g["val_auc_mean"],  g["val_auc_std"])],
        "test_bacc":  [f"{m:.3f}+/-{s:.3f}" for m,s in zip(g["test_bacc_mean"], g["test_bacc_std"])],
        "test_auc":   [f"{m:.3f}+/-{s:.3f}" for m,s in zip(g["test_auc_mean"],  g["test_auc_std"])],
    })
    fig_h = max(4.5, 0.30 * len(show) + 3.5)
    fig, ax = plt.subplots(figsize=(18, fig_h)); ax.axis("off")
    ax.set_title(f"sCN vs pCN -- strict baseline-CN cohort   "
                  f"LD={args.ld_config}  beta={args.beta_source}  "
                  f"(mean +/- std over 3 seeds)\n"
                  f"FM rows: attention provenance annotated in [..] "
                  f"-- see footer for the labels each attention pool was trained on",
                  fontsize=10, fontweight="bold", pad=10, loc="left")
    # First column widened from 0.24 -> 0.42 so the longest FM source name fits.
    col_widths = [0.42, 0.16, 0.08, 0.08, 0.085, 0.085, 0.085, 0.085]
    tbl = ax.table(cellText=show.values.tolist(), colLabels=show.columns.tolist(),
                    loc="center", cellLoc="center", colWidths=col_widths)
    tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1, 1.2)
    for j in range(len(show.columns)):
        c = tbl[(0, j)]; c.set_facecolor("#222"); c.set_text_props(color="white", weight="bold")
    # Left-align the long source column so the variant name reads cleanly.
    src_col = show.columns.get_loc("source")
    for row_idx in range(len(show)):
        cell = tbl[(row_idx + 1, src_col)]
        cell.set_text_props(ha="left")
        cell.PAD = 0.02
    # Bold the TOP-3 means per metric column (val/test BalAcc + AUC). Tie-aware
    # at display precision (3 dp): if a tied value sits at the boundary, all
    # tied rows are bolded together so we don't arbitrarily pick one.
    BOLD_COLS = {
        "val_bacc":  ("val_bacc_mean",  show.columns.get_loc("val_bacc")),
        "test_bacc": ("test_bacc_mean", show.columns.get_loc("test_bacc")),
        "val_auc":   ("val_auc_mean",   show.columns.get_loc("val_auc")),
        "test_auc":  ("test_auc_mean",  show.columns.get_loc("test_auc")),
    }
    g_indexed = g.reset_index(drop=True)
    for _key, (mean_col, col_pos) in BOLD_COLS.items():
        vals_3dp = g_indexed[mean_col].round(3)
        top3_threshold = sorted(vals_3dp.dropna().unique(), reverse=True)[:3]
        if not top3_threshold: continue
        cutoff = top3_threshold[-1]
        for row_idx in g_indexed.index:
            v = vals_3dp.iloc[row_idx]
            if pd.notna(v) and v >= cutoff:
                tbl[(row_idx + 1, col_pos)].set_text_props(weight="bold")
    # Per-seed cohort breakdown footer — same per source so use the first source.
    per_seed = (df[df["source"] == df["source"].iloc[0]]
                 [df["covar_mode"] == df["covar_mode"].iloc[0]]
                 .sort_values("seed"))
    footer_lines = [
        ("val_n / test_n columns above are the mean over 3 seeds; the "
         "positive class is pCN = pCN_to_AD union pCN_to_MCI."),
        "Per-seed cohort sizes (val / test, with pCN positives in brackets):",
    ]
    for _, r in per_seed.iterrows():
        footer_lines.append(
            f"  seed {int(r['seed'])}: "
            f"val n={int(r['val_n'])} (pCN={int(r['val_pos'])}, sCN={int(r['val_neg'])})    "
            f"test n={int(r['test_n'])} (pCN={int(r['test_pos'])}, sCN={int(r['test_neg'])})")
    footer_lines.append("")
    footer_lines.append("FM attention provenance (the labels each ChromHierPool"
                         " was trained on -- attn=... in row source):")
    footer_lines.append("  [attn=T1_AD_CN]          ChromHierPool trained on the "
                         "ORIGINAL Task 1 labels: AD_bl + pMCI + CN_to_AD (positive) "
                         "vs sCN (negative).")
    footer_lines.append("                            Softly biased toward AD-onset "
                         "signal; included as a reference / leakage-style baseline.")
    footer_lines.append("  [attn=fixed_aggregator]   Deterministic _fixed_aggregator_v3 "
                         "(gamma=delta=1). No training -- task-independent. The XGB "
                         "head still sees raw 771-d features.")
    footer_lines.append("  [attn=sCN_vs_pCN]         ChromHierPool RE-TRAINED on this "
                         "task's labels: sCN (negative) vs pCN_to_AD + pCN_to_MCI "
                         "(positive). The fair task-conditioned comparison.")
    footer_lines.append("")
    footer_lines.append("Bold cells: top-3 within each of val_bacc / val_auc / "
                         "test_bacc / test_auc (tie-aware at 3 dp).")
    fig.text(0.01, 0.01, "\n".join(footer_lines), fontsize=8,
              family="monospace", va="bottom")
    png = out_dir / "leaderboard.png"
    fig.savefig(png, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"wrote PNG: {png}")
    # WandB: log final leaderboard as a Table + the PNG as an Image.
    if wandb_run is not None:
        try:
            import wandb as _wb
            wandb_run.log({"leaderboard": _wb.Table(dataframe=g),
                            "leaderboard_png": _wb.Image(str(png))})
        except Exception as e:
            print(f"[wandb] table/image log skipped ({e})")

    # Per-patient parquets + multimodal anchors.
    if capture_rows:
        all_rows = pd.concat(capture_rows, ignore_index=True)
        pp_dir = out_dir / "per_patient"
        anchor_dir = (ANCHOR_ROOT / "prs" / TASK_NAME
                       if args.beta_source == "raw"
                       else ANCHOR_ROOT / f"prs_{args.beta_source}" / TASK_NAME)
        for (src, mode, seed), grp in all_rows.groupby(
                ["model","covar_mode","seed"], sort=False):
            save_per_patient_parquet(grp, pp_dir,
                                       source=src, seed=int(seed), covar_mode=mode)
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
