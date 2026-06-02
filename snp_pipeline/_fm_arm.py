"""
_fm_arm.py
==========
Shared helpers for adding BMFM-derived FM arms to the new-task scripts 49a/b/c.

Workflow per (variant, seed, covar_mode, task):
  1. Load the variant's per-PTID 771-d input embedding from script 50's anchor
     (`outputs/.../multimodal_anchors/<variant>/seed_{s}/embeddings.npz`).
  2. Subset the embedding to the task's cohort PTIDs.
  3. Build feature matrix (771-d input ± covariates) for train/val/test folds.
  4. Train a fresh head on the new task labels (variant 1 -> MLPClassifier;
     variant 2 -> XGBClassifier). Sample weights = "balanced".
  5. Return per-fold predictions + (for MLP variant) the 256-d penultimate of
     the new task-trained head.

Per-task anchors are written to
  outputs/strict_qc_prs/multimodal_anchors/<variant>/task_<task>/seed_{s}/
    - predictions.parquet  (canonical schema)
    - embeddings_taskhead.npz   (variant 1 only: embedding_256_taskhead, plus
                                  the 771-d input the head consumed for self-
                                  contained re-fittability)
    - embeddings_input.npz       (the 771-d input, both variants)

Attention bias caveat: variant 1's 771-d is produced by a TRAINED
ChromHierPool whose attention weights were fitted on the original AD/CN task.
The representation is therefore SOFTLY BIASED toward the original task signal.
Saving both the input 771-d and the new-task 256-d penultimate makes this
auditable; downstream analyses can re-fit anything from disk.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.utils.class_weight import compute_sample_weight

# ── BMFM variant registry (matches script 50's VARIANTS list) ──────────────
# Anchors live at `outputs/multimodal_anchors/<variant_desc_name>/seed_{s}/`.
ANCHOR_ROOT = Path("D:/ADNI_SNP_Omni2.5M_20140220/outputs/multimodal_anchors")

FM_VARIANTS = [
    dict(
        desc_name="bmfm_snp__attn_bias_single__chrom_hier__mlp2__101bp",
        # CRITICAL: attention is task-1-conditioned (trained on AD/CN).
        # short_name surfaces this in every leaderboard row.
        short_name="BMFM-SNP attn_bias_single chrom_hier MLP [attn=T1_AD_CN]",
        head_kind="mlp",
        attention_source="T1_AD_CN_trained",
    ),
    dict(
        desc_name="bmfm_snp__none__global_attn__xgb__1001bp",
        # Attention is the deterministic fixed-aggregator (γ=δ=1); no task bias.
        short_name="BMFM-SNP none global_attn XGB [attn=fixed_aggregator]",
        head_kind="xgb",
        attention_source="fixed_aggregator",
    ),
]


# ── Wave 2: task-conditioned BMFM variants (Colab sweep 2026-06-02) ─────────
# Each (task, variant) pair has its own attention pool fitted on that task's
# labels via the new --task arg in 30v3_train_diff_attention_func.py. After
# the Colab sweep completes, outputs land on Drive at
#   diff_attn_v2_upload/fm_new_tasks_top_predictions/<task>/mlp2/
#       101bp_bmfm_snp_<mode>_chrom_hier_s{0,1,2}/
#         {model.pt, predictions.tsv, embeddings.npz, metrics.json}
# To use locally, sync this tree to D:/ADNI_SNP_Omni2.5M_20140220/
#   diff_attn_v2_upload/fm_new_tasks_top_predictions/  (Drive Desktop, rclone,
# or manual copy). The variant registry below points at the local D: paths.
TASK_CONDITIONED_ROOT = Path(
    "D:/ADNI_SNP_Omni2.5M_20140220/diff_attn_v2_upload/"
    "fm_new_tasks_top_predictions")

_TASK_CONDITIONED_MODES = ("none", "attn_bias_single", "attn_bias_per_modality")
_TASK_CONDITIONED_MODELS = ("bmfm_snp", "ntv2")

# 49a/b/c TASK_NAME convention -> trainer task name used by 30v3
# (the Colab sweep writes its output dirs under the TRAINER name, so the
# variant_subdir below must use that form).
_TRAINER_TASK_FOR_49 = {
    "sCN_vs_pCN":  "sCN_vs_pCN",
    "horizon_3y":  "T3_horizon_3y",
    "horizon_5y":  "T3_horizon_5y",
    "horizon_7y":  "T3_horizon_7y",
    "horizon_10y": "T3_horizon_10y",
    "T4_3class":   "T4_3class",
}
_MODEL_LABEL = {"bmfm_snp": "BMFM-SNP", "ntv2": "NTv2"}


def _task_conditioned_variant(task_49: str, model: str, mode: str) -> dict:
    """One task-conditioned variant entry covering a (task, model, mode)
    triple. `task_49` is the 49a/b/c TASK_NAME (e.g. 'horizon_5y'); the
    Drive sub-path below uses the trainer's task name ('T3_horizon_5y')
    because that is what the Colab sweep writes."""
    trainer_task = _TRAINER_TASK_FOR_49.get(task_49, task_49)
    model_lbl = _MODEL_LABEL[model]
    short = f"{model_lbl} {mode} chrom_hier MLP [attn={task_49}]"
    return dict(
        desc_name=f"{model}__{mode}__chrom_hier__mlp2__101bp__attn_{task_49}",
        short_name=short,
        head_kind="mlp",
        attention_source=f"{task_49}_trained",
        # Drive sub-path where Colab saved this variant's outputs.
        variant_subdir=f"{trainer_task}/mlp2/101bp_{model}_{mode}_chrom_hier",
    )


FM_VARIANTS_TASK_CONDITIONED: dict[str, list[dict]] = {
    task: [_task_conditioned_variant(task, model, mode)
            for model in _TASK_CONDITIONED_MODELS
            for mode in _TASK_CONDITIONED_MODES]
    for task in ("sCN_vs_pCN",
                  "horizon_3y", "horizon_5y", "horizon_7y", "horizon_10y",
                  "T4_3class")
}


def task_conditioned_variants_for(task_name: str) -> list[dict]:
    """Return the 3 task-conditioned BMFM variants for `task_name`.
    `task_name` matches the 49a/b/c TASK_NAME convention (e.g. 'sCN_vs_pCN',
    'horizon_5y', 'T4_3class'). Returns empty list if the task isn't in
    the registry (e.g. 'classification' or 'cox' from the upstream PRS scripts).
    """
    return FM_VARIANTS_TASK_CONDITIONED.get(task_name, [])


def load_task_conditioned_embedding(variant: dict, seed: int) -> dict | None:
    """Load embeddings.npz from a Wave-2 task-conditioned anchor (Colab sweep
    output, synced to local D:). Same return shape as load_variant_embedding().
    Returns None if the file isn't on disk yet (sweep still running or not
    yet synced from Drive to D:)."""
    sub = variant.get("variant_subdir")
    if not sub:
        return None
    path = TASK_CONDITIONED_ROOT / sub / f"s{seed}" / "embeddings.npz"
    # The sweep saves under .../<mode>_chrom_hier_s{seed}/embeddings.npz; the
    # _s{seed} suffix is part of the directory name (no separate seed dir).
    direct = TASK_CONDITIONED_ROOT / (sub + f"_s{seed}") / "embeddings.npz"
    p = direct if direct.exists() else path
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=False)
    out = {
        "ptids":         np.array([str(x) for x in z["ptids"]]),
        "fold":          np.array([str(x) for x in z["fold"]]),
        "embedding_771": z["embedding_771"].astype(np.float32),
    }
    if "embedding_256" in z.files:
        out["embedding_256"] = z["embedding_256"].astype(np.float32)
    return out

# Disclaimer text reproduced in every per-task PNG and saved README. ASCII
# only so Windows cp1252 stdout/streams render the text correctly.
FM_ATTENTION_DISCLAIMER = (
    "FM rows fall into two waves. Wave 1: Task-1-conditioned attention - "
    "BMFM-SNP 771-d embeddings from script 50, whose ChromHierPool "
    "(gamma/delta/epsilon + Q/K/V) was trained on the original AD vs CN labels "
    "and is therefore SOFTLY BIASED toward that signal; plus the XGB variant "
    "whose deterministic _fixed_aggregator_v3 (gamma=delta=1) is task-"
    "independent. Wave 2: task-conditioned attention - BMFM-SNP and NTv2 "
    "ChromHierPools re-trained on the current task's labels (Colab sweep "
    "2026-06-02). Each row's source column annotates its attention provenance."
)


def load_embedding_for_variant(variant: dict, seed: int) -> dict | None:
    """Unified loader. Picks the task-conditioned loader if `variant` has the
    `variant_subdir` field (Wave-2 entries from FM_VARIANTS_TASK_CONDITIONED),
    or the original Task-1 loader otherwise."""
    if "variant_subdir" in variant:
        return load_task_conditioned_embedding(variant, seed)
    return load_variant_embedding(variant["desc_name"], seed)


def load_variant_embedding(desc_name: str, seed: int) -> dict | None:
    """Load embeddings.npz from a script-50 anchor.

    Returns dict with keys: ptids (str ndarray), fold (str ndarray),
    embedding_771 (np.float32 (N, 771)), optionally embedding_256 (N, 256).
    Returns None if the file is missing.
    """
    path = ANCHOR_ROOT / desc_name / f"seed_{seed}" / "embeddings.npz"
    if not path.exists():
        return None
    z = np.load(path, allow_pickle=False)
    out = {
        "ptids":         np.array([str(p) for p in z["ptids"]]),
        "fold":          np.array([str(f) for f in z["fold"]]),
        "embedding_771": z["embedding_771"].astype(np.float32),
    }
    if "embedding_256" in z.files:
        out["embedding_256"] = z["embedding_256"].astype(np.float32)
    return out


def _stack_with_covars(emb: np.ndarray, cov_df: pd.DataFrame | None,
                         ptids: np.ndarray, covars: list[str]) -> np.ndarray:
    """Stack a 771-d embedding with optional covariates (looked up by PTID).
    Covariate NaNs are mean-imputed using the input slice itself."""
    if not covars or cov_df is None:
        return emb
    sub = (cov_df.set_index("Patient_ID").loc[ptids, covars]
                 .reset_index(drop=True))
    arr = sub.astype(float).values
    col_means = np.nanmean(arr, axis=0)
    inds = np.where(np.isnan(arr))
    arr[inds] = np.take(col_means, inds[1])
    return np.concatenate([emb, arr.astype(np.float32)], axis=1)


def _mlp_penultimate(mlp: MLPClassifier, X: np.ndarray) -> np.ndarray:
    """Compute MLPClassifier's hidden-layer activation. Assumes
    `hidden_layer_sizes=(256,)` (single hidden layer) and the default
    `activation="relu"`. Returns (n, 256)."""
    Z = X @ mlp.coefs_[0] + mlp.intercepts_[0]
    act = mlp.activation
    if act == "relu":     return np.maximum(Z, 0).astype(np.float32)
    if act == "logistic": return (1.0 / (1.0 + np.exp(-Z))).astype(np.float32)
    if act == "tanh":     return np.tanh(Z).astype(np.float32)
    if act == "identity": return Z.astype(np.float32)
    raise ValueError(f"Unsupported MLPClassifier activation: {act}")


def run_fm_head_binary(
    *,
    variant: dict,
    seed: int,
    fold_ptid: dict,           # {"train": [ptids], "val": [ptids], "test": [ptids]}
    fold_y:    dict,           # {"train": np.ndarray int, ...}
    cov_df: pd.DataFrame | None,
    covars: list[str],
) -> dict | None:
    """Train a fresh head on a binary task using variant's 771-d input.

    Returns dict with per-fold predicted probability arrays AND (variant 1
    only) per-fold 256-d task-trained penultimate, plus the 771-d input
    matrices the head consumed:
      {
        "fold": {"train": {"ptids", "y", "X_input_771", "proba",
                              "embedding_256_taskhead"|None},
                 "val":   ...,
                 "test":  ...},
        "head_kind": "mlp"|"xgb",
      }
    Returns None if the variant's anchor is missing or insufficient data.
    """
    emb = load_embedding_for_variant(variant, seed)
    if emb is None:
        print(f"  [fm:{variant['head_kind']}] seed={seed}: missing anchor for "
               f"{variant['desc_name']} -- skipping")
        return None
    ptid_to_row = {p: i for i, p in enumerate(emb["ptids"].tolist())}

    def _slice(fold: str):
        wanted = [p for p in fold_ptid[fold] if p in ptid_to_row]
        if not wanted: return None
        idx = np.array([ptid_to_row[p] for p in wanted], dtype=np.int64)
        X771 = emb["embedding_771"][idx]
        y    = np.array([int(fold_y[fold][i]) for i, p in enumerate(fold_ptid[fold])
                          if p in ptid_to_row], dtype=np.int64)
        return {"ptids": wanted, "X_input_771": X771, "y": y, "idx": idx}

    slices = {sp: _slice(sp) for sp in ("train", "val", "test")}
    if slices["train"] is None or len(slices["train"]["y"]) < 20:
        print(f"  [fm:{variant['head_kind']}] seed={seed}: train slice empty/"
               f"too small -- skipping")
        return None
    if len(set(slices["train"]["y"])) < 2:
        return None

    # Build feature matrices (with optional covariates).
    Xs, ys = {}, {}
    for sp, sl in slices.items():
        if sl is None: continue
        Xfull = _stack_with_covars(sl["X_input_771"], cov_df,
                                     np.array(sl["ptids"]), covars)
        Xs[sp] = Xfull
        ys[sp] = sl["y"]

    head_kind = variant["head_kind"]
    sample_weight = compute_sample_weight("balanced", ys["train"])

    out_fold = {}
    embedding_256_per_fold = {}

    if head_kind == "mlp":
        scaler = StandardScaler().fit(Xs["train"])
        Xtr_s = scaler.transform(Xs["train"]).astype(np.float32)
        # Single hidden layer width=256 mirrors the original MLPHead.
        mlp = MLPClassifier(
            hidden_layer_sizes=(256,), activation="relu",
            alpha=1e-3, max_iter=2000, random_state=seed,
            early_stopping=False, solver="adam")
        try:
            mlp.fit(Xtr_s, ys["train"])
        except Exception as e:
            print(f"  [fm:mlp] seed={seed}: fit failed ({e})")
            return None
        for sp, X in Xs.items():
            Xs_t = scaler.transform(X).astype(np.float32)
            proba = mlp.predict_proba(Xs_t)
            pos_proba = proba[:, list(mlp.classes_).index(1)]
            penult = _mlp_penultimate(mlp, Xs_t)
            embedding_256_per_fold[sp] = penult
            out_fold[sp] = {
                "ptids":  slices[sp]["ptids"],
                "y":      ys[sp],
                "X_input_771": slices[sp]["X_input_771"],
                "proba":  pos_proba.astype(np.float32),
                "embedding_256_taskhead": penult,
            }
    elif head_kind == "xgb":
        import xgboost as xgb
        n_pos = int((ys["train"] == 1).sum())
        n_neg = int((ys["train"] == 0).sum())
        spw = (n_neg / max(1, n_pos))
        clf = xgb.XGBClassifier(
            n_estimators=200, max_depth=4, learning_rate=0.05,
            random_state=seed, scale_pos_weight=spw,
            objective="binary:logistic", eval_metric="logloss",
            tree_method="hist", verbosity=0)
        try:
            clf.fit(Xs["train"], ys["train"])
        except Exception as e:
            print(f"  [fm:xgb] seed={seed}: fit failed ({e})")
            return None
        for sp, X in Xs.items():
            proba = clf.predict_proba(X)
            pos_proba = proba[:, list(clf.classes_).index(1)]
            out_fold[sp] = {
                "ptids":  slices[sp]["ptids"],
                "y":      ys[sp],
                "X_input_771": slices[sp]["X_input_771"],
                "proba":  pos_proba.astype(np.float32),
                "embedding_256_taskhead": None,
            }
    else:
        raise ValueError(f"Unknown head_kind: {head_kind}")

    return {"fold": out_fold, "head_kind": head_kind}


def run_fm_head_multiclass(
    *,
    variant: dict,
    seed: int,
    fold_ptid: dict,           # {"train": [ptids], ...}
    fold_y:    dict,           # {"train": int ndarray with K classes}
    cov_df: pd.DataFrame | None,
    covars: list[str],
    n_classes: int = 3,
) -> dict | None:
    """K-class version of run_fm_head_binary. Variant 1 uses MLPClassifier
    (handles multiclass natively); variant 2 uses XGBClassifier with
    objective='multi:softprob'. Returns the same dict shape but the per-fold
    'proba' is now (N, K) and 'pred_class' is also included."""
    emb = load_embedding_for_variant(variant, seed)
    if emb is None:
        print(f"  [fm:{variant['head_kind']}/multi] seed={seed}: missing anchor "
               f"for {variant['desc_name']} -- skipping")
        return None
    ptid_to_row = {p: i for i, p in enumerate(emb["ptids"].tolist())}

    def _slice(fold: str):
        wanted = [p for p in fold_ptid[fold] if p in ptid_to_row]
        if not wanted: return None
        idx = np.array([ptid_to_row[p] for p in wanted], dtype=np.int64)
        X771 = emb["embedding_771"][idx]
        keep_mask = np.array([p in ptid_to_row for p in fold_ptid[fold]])
        y = np.array(fold_y[fold])[keep_mask].astype(np.int64)
        return {"ptids": wanted, "X_input_771": X771, "y": y, "idx": idx}

    slices = {sp: _slice(sp) for sp in ("train", "val", "test")}
    if slices["train"] is None or len(slices["train"]["y"]) < 20:
        return None
    if len(set(slices["train"]["y"])) < 2:
        return None

    Xs, ys = {}, {}
    for sp, sl in slices.items():
        if sl is None: continue
        Xfull = _stack_with_covars(sl["X_input_771"], cov_df,
                                     np.array(sl["ptids"]), covars)
        Xs[sp] = Xfull
        ys[sp] = sl["y"]

    head_kind = variant["head_kind"]
    sample_weight = compute_sample_weight("balanced", ys["train"])

    out_fold = {}
    if head_kind == "mlp":
        scaler = StandardScaler().fit(Xs["train"])
        Xtr_s = scaler.transform(Xs["train"]).astype(np.float32)
        mlp = MLPClassifier(
            hidden_layer_sizes=(256,), activation="relu",
            alpha=1e-3, max_iter=2000, random_state=seed,
            early_stopping=False, solver="adam")
        try:
            mlp.fit(Xtr_s, ys["train"])
        except Exception as e:
            print(f"  [fm:mlp/multi] seed={seed}: fit failed ({e})")
            return None
        seen = list(mlp.classes_)
        def _full_proba(arr):
            full = np.zeros((arr.shape[0], n_classes), dtype=float)
            for col_idx, c in enumerate(seen):
                full[:, int(c)] = arr[:, col_idx]
            return full
        for sp, X in Xs.items():
            Xs_t = scaler.transform(X).astype(np.float32)
            proba = _full_proba(mlp.predict_proba(Xs_t))
            penult = _mlp_penultimate(mlp, Xs_t)
            out_fold[sp] = {
                "ptids":  slices[sp]["ptids"],
                "y":      ys[sp],
                "X_input_771": slices[sp]["X_input_771"],
                "proba":  proba.astype(np.float32),       # (N, K)
                "embedding_256_taskhead": penult,
            }
    elif head_kind == "xgb":
        import xgboost as xgb
        clf = xgb.XGBClassifier(
            n_estimators=200, max_depth=4, learning_rate=0.05,
            random_state=seed, sample_weight=None,
            objective="multi:softprob", num_class=n_classes,
            eval_metric="mlogloss", tree_method="hist", verbosity=0)
        try:
            clf.fit(Xs["train"], ys["train"], sample_weight=sample_weight)
        except Exception as e:
            print(f"  [fm:xgb/multi] seed={seed}: fit failed ({e})")
            return None
        seen = list(clf.classes_)
        def _full_proba(arr):
            full = np.zeros((arr.shape[0], n_classes), dtype=float)
            for col_idx, c in enumerate(seen):
                full[:, int(c)] = arr[:, col_idx]
            return full
        for sp, X in Xs.items():
            proba = _full_proba(clf.predict_proba(X))
            out_fold[sp] = {
                "ptids":  slices[sp]["ptids"],
                "y":      ys[sp],
                "X_input_771": slices[sp]["X_input_771"],
                "proba":  proba.astype(np.float32),
                "embedding_256_taskhead": None,
            }
    else:
        raise ValueError(f"Unknown head_kind: {head_kind}")
    return {"fold": out_fold, "head_kind": head_kind, "n_classes": n_classes}


def save_per_task_anchor(
    *,
    variant: dict,
    task_name: str,
    seed: int,
    covar_mode: str,
    fold_capture: dict,
    anchor_root: Path = ANCHOR_ROOT,
) -> Path:
    """Write the per-task FM anchor (predictions.parquet + embeddings npz).

    Layout:
      <ANCHOR_ROOT>/<variant>/task_<task>/seed_{s}/<covar_mode>/
        - predictions.parquet
        - embeddings_input.npz       (variant 1 + 2: 771-d input)
        - embeddings_taskhead.npz    (variant 1 only: new MLP 256-d penultimate)
    """
    out_dir = (anchor_root / variant["desc_name"]
                / f"task_{task_name}"
                / f"seed_{seed}"
                / covar_mode)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Predictions parquet (canonical multimodal schema + attention_source
    # column so the bias provenance is visible at every downstream read).
    # For multiclass (`proba` shape (N, K)) we save outcome_proba = argmax-class's
    # probability and add proba_class0..K-1 columns. Binary stays scalar.
    pred_rows = []
    attn_src = variant.get("attention_source", "unknown")
    for fold in ("train", "val", "test"):
        if fold not in fold_capture: continue
        rec = fold_capture[fold]
        n = len(rec["ptids"])
        proba = np.asarray(rec["proba"], dtype=float)
        is_multiclass = (proba.ndim == 2 and proba.shape[1] > 1)
        if is_multiclass:
            pred_class = proba.argmax(axis=1)
            outcome_proba = proba[np.arange(n), pred_class]
        else:
            outcome_proba = proba.reshape(-1)
        cols = {
            "Patient_ID":       rec["ptids"],
            "seed":             np.full(n, seed, dtype=np.int8),
            "fold":             np.full(n, fold, dtype=object),
            "task":             np.full(n, task_name, dtype=object),
            "model":            np.full(n, variant["short_name"], dtype=object),
            "attention_source": np.full(n, attn_src, dtype=object),
            "covar_mode":       np.full(n, covar_mode, dtype=object),
            "beta_source":      np.full(n, "na", dtype=object),
            "ld_config":        np.full(n, "na", dtype=object),
            "prs_raw":          np.full(n, np.nan, dtype=float),
            "prs_z":            np.full(n, np.nan, dtype=float),
            "outcome_proba":    outcome_proba,
            "y_true":           np.asarray(rec["y"], dtype=int),
        }
        if is_multiclass:
            for k in range(proba.shape[1]):
                cols[f"proba_class{k}"] = proba[:, k]
        pred_rows.append(pd.DataFrame(cols))
    if pred_rows:
        pd.concat(pred_rows, ignore_index=True).to_parquet(
            out_dir / "predictions.parquet", index=False)
    # Also drop a README describing the attention provenance.
    (out_dir / "README.txt").write_text(
        f"variant: {variant['desc_name']}\n"
        f"short_name: {variant['short_name']}\n"
        f"attention_source: {attn_src}\n"
        f"task: {task_name}    seed: {seed}    covar_mode: {covar_mode}\n"
        f"\n{FM_ATTENTION_DISCLAIMER}\n",
        encoding="utf-8")

    # 771-d input npz (saved for both variants; supports auditability of
    # the soft task bias for variant 1 + makes the per-task anchor self-
    # contained and re-fittable).
    ptids_all, fold_all, y_all, prob_all = [], [], [], []
    in_all = []
    for fold in ("train", "val", "test"):
        if fold not in fold_capture: continue
        rec = fold_capture[fold]
        ptids_all.extend(rec["ptids"])
        fold_all.extend([fold] * len(rec["ptids"]))
        y_all.append(np.asarray(rec["y"], dtype=np.int8))
        # proba can be (N,) for binary or (N, K) for multiclass; keep the
        # native shape so downstream code can read either.
        prob_all.append(np.asarray(rec["proba"], dtype=np.float32))
        in_all.append(np.asarray(rec["X_input_771"], dtype=np.float32))
    if in_all:
        np.savez_compressed(
            out_dir / "embeddings_input.npz",
            ptids=np.array(ptids_all),
            fold=np.array(fold_all),
            y_true=np.concatenate(y_all),
            proba=np.concatenate(prob_all, axis=0),
            embedding_771=np.concatenate(in_all, axis=0),
        )

    # 256-d task-trained MLP penultimate (variant 1 only).
    th_all = []
    for fold in ("train", "val", "test"):
        if fold not in fold_capture: continue
        rec = fold_capture[fold]
        if rec.get("embedding_256_taskhead") is None:
            th_all = []
            break
        th_all.append(np.asarray(rec["embedding_256_taskhead"],
                                    dtype=np.float32))
    if th_all:
        np.savez_compressed(
            out_dir / "embeddings_taskhead.npz",
            ptids=np.array(ptids_all),
            fold=np.array(fold_all),
            y_true=np.concatenate(y_all),
            proba=np.concatenate(prob_all),
            embedding_256_taskhead=np.concatenate(th_all, axis=0),
        )
    return out_dir
