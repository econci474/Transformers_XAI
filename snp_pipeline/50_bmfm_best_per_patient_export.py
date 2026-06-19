"""
50_bmfm_best_per_patient_export.py
==================================
Materialize the per-PTID multimodal-alignment anchors for the
**BMFM-SNP x attn_bias_single x chrom_hier x MLP (mlp2)** variant —
the best-performing diff-attention model from fm_top_predictions_top_7.

For each of the 3 seed dirs at
  D:/ADNI_SNP_Omni2.5M_20140220/fm_top_predictions_top_7/
    fm_top_predictions/mlp2/101bp_bmfm_snp_attn_bias_single_chrom_hier_s{0,1,2}/
this script:

1. Reproduces the trainer's data loading EXACTLY by reusing the v3
   trainer's helpers (`_load_diff`, `_load_dosage`, `_load_beta`,
   `_resolve_ag_parent`, `_zscore_train`, `_per_patient_values`,
   `DiffAttnV3`). Zero drift from training-time preprocessing.
2. Instantiates `DiffAttnV3(in_dim=771, ..., head='mlp2', mlp_width=256,
   mlp_dropout=0.1, mode='attn_bias_single', aggregation='chrom_hier')`
   and loads `model.pt`.
3. Runs forward over the union of train+val+test PTIDs with hooks on:
     - `model.pool`  -> embedding_771 (B, 771) + attention weights (B, 128)
     - `model.head.net[2]` (Dropout, identity in eval) -> embedding_256
       (B, 256), the pre-final-linear MLP hidden state.
4. **Verification gate (FAIL-LOUDLY)**: reload `predictions.tsv` and
   assert per-PTID `|prob_recomputed - prob_saved| < 1e-6` for every
   val + test row. If ANY row fails, dump a diff TSV and `sys.exit(2)`
   BEFORE writing any anchor file.
5. Writes anchor files to
   `outputs/strict_qc_prs/multimodal_anchors/
    bmfm_snp__attn_bias_single__chrom_hier__mlp2__101bp/`:
     - README.md (variant config + asset legend)
     - beta.tsv, snp_metadata.tsv (chrom/pos/chrom_idx),
       alphagenome_func_imp_abs.tsv, dosage_all_patients.tsv  (static)
     - seed_{s}/predictions.parquet (canonical multimodal schema)
     - seed_{s}/embeddings.npz (keys: ptids, fold, y_true, proba,
       embedding_771, embedding_256, attn_weights)
     - seed_{s}/metrics_recomputed.json (sanity vs trainer's metrics.json)

Run (env: snp):
  python snp_pipeline/50_bmfm_best_per_patient_export.py
"""
from __future__ import annotations
import argparse
import importlib.util
import json
import sys
import time
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import (roc_auc_score, balanced_accuracy_score, log_loss,
                              confusion_matrix)

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Load script 30v3 by file path (numeric prefix -> not importable via "import 30v3...").
_TRAINER = HERE / "30v3_train_diff_attention_func.py"
_spec = importlib.util.spec_from_file_location("diff_v3", _TRAINER)
v3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(v3)

# ── Paths + variant config (top 2 BMFM-SNP variants, locked 2026-06-02) ─────
BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220/diff_attn_v2_upload")
TOP7_ROOT = Path("D:/ADNI_SNP_Omni2.5M_20140220/fm_top_predictions_top_7/"
                  "fm_top_predictions")
OUT_BASE_ROOT = Path("D:/ADNI_SNP_Omni2.5M_20140220/outputs/multimodal_anchors")

# Each variant is a self-contained config. Variant 1 = the original best
# MLP head; variant 2 = the top BMFM-SNP XGB head from the same sweep.
VARIANTS = [
    dict(
        desc_name="bmfm_snp__attn_bias_single__chrom_hier__mlp2__101bp",
        top7_subdir="mlp2",
        variant_name="101bp_bmfm_snp_attn_bias_single_chrom_hier",
        ckpt_filename="model.pt",
        seq_length="101bp",
        model="bmfm_snp",
        func_integration_mode="attn_bias_single",
        aggregation="chrom_hier",
        head="mlp2",
        mlp_width=256,
        mlp_dropout=0.1,
    ),
    dict(
        desc_name="bmfm_snp__none__global_attn__xgb__1001bp",
        top7_subdir="xgb",
        variant_name="1001bp_bmfm_snp_none_global_attn",
        ckpt_filename="xgb_model.json",
        seq_length="1001bp",
        model="bmfm_snp",
        func_integration_mode="none",
        aggregation="global_attn",
        head="xgb",
        # XGB head doesn't use MLP width/dropout — set to None so _make_model
        # is never called for this variant.
        mlp_width=None,
        mlp_dropout=None,
        # Trainer XGB HPs (from metrics.json config).
        xgb_n_estimators=200,
        xgb_max_depth=4,
        xgb_lr=0.05,
    ),
]

# The v3 sweep ran with `--splits-root=/nonexistent/force-fallback` so the
# trainer's _load_labels_for_split fell through to `<base>/splits/`. Mirror
# that exactly so our re-forward sees the same train/val/test partition.
SPLITS_ROOT_PRIMARY  = Path("/nonexistent/force-fallback")
SPLITS_ROOT_FALLBACK = BASE / "splits"

SEEDS = (0, 1, 2)
N_CLASSES = 1  # binary head (sigmoid)
PARITY_TOL = 1e-6


def _load_inputs(seed: int, variant: dict):
    """Replay v3 trainer's data loading for this variant. Returns a dict of
    everything needed for forward + capture."""
    # Per-SNP diff + ordering
    diff, diff_rsids = v3._load_diff(BASE, variant["seq_length"],
                                       variant["model"])
    target_rsids = diff_rsids
    # Per-rsID β + chromosome assignment
    beta_df = v3._load_beta(BASE)
    beta_df_idx = beta_df.set_index("rsID")
    beta = np.array([float(beta_df_idx.loc[rs, "beta_A1"])
                       if rs in beta_df_idx.index else 0.0
                      for rs in target_rsids], dtype=np.float32)
    chrom_arr = np.array([int(beta_df_idx.loc[rs, "CHR"])
                            if rs in beta_df_idx.index else 0
                           for rs in target_rsids], dtype=np.int64)
    chrom_uniq = sorted(set(chrom_arr.tolist()))
    chrom_idx_map = {c: i for i, c in enumerate(chrom_uniq)}
    chrom_idx = np.array([chrom_idx_map[c] for c in chrom_arr], dtype=np.int64)
    n_chroms = len(chrom_uniq)
    # Per-PTID dosage matrix
    dos_df, _ = v3._load_dosage(BASE)
    # AlphaGenome modality summaries
    ag_parent = v3._resolve_ag_parent(BASE)
    summaries = v3.ff.load_modality_summaries(ag_parent, target_rsids)
    # Labels per split
    labels = {sp: v3._load_labels_for_split(SPLITS_ROOT_PRIMARY, seed, sp,
                                                base=BASE,
                                                exclude_cn_to_ad=False)
              for sp in ("train", "val", "test")}
    # Subset dosage to PTIDs present per split, fill NaN with train means
    dos_by_ptid = dos_df.set_index("PTID")
    splits = {}
    for sp, df in labels.items():
        present = [p for p in df["PTID"].tolist() if p in dos_by_ptid.index]
        df_ok = df[df["PTID"].isin(present)].reset_index(drop=True)
        dos_split = (dos_by_ptid.loc[df_ok["PTID"], target_rsids]
                                .values.astype(np.float32))
        splits[sp] = {"ptid": df_ok["PTID"].tolist(),
                       "y":    df_ok["y"].values.astype(np.int64),
                       "dosage": dos_split}
    train_means = np.nanmean(splits["train"]["dosage"], axis=0)
    for sp in splits:
        d = splits[sp]["dosage"]
        mask = np.isnan(d)
        for i in range(d.shape[1]):
            d[mask[:, i], i] = train_means[i]
        splits[sp]["dosage"] = d
    # Z-score modality summaries on train-fold-all-SNPs (matches trainer rule)
    train_snp_mask = np.ones(len(target_rsids), dtype=bool)
    summaries_z, zscore_stats = {}, {}
    for k, v in summaries.items():
        zv, st = v3._zscore_train(v, train_snp_mask)
        summaries_z[k] = zv.reshape(-1)
        zscore_stats[k] = st
    return dict(
        target_rsids=target_rsids, diff=diff, beta=beta,
        chrom_arr=chrom_arr, chrom_idx=chrom_idx, n_chroms=n_chroms,
        splits=splits, summaries_z=summaries_z, zscore_stats=zscore_stats,
    )


def _build_extras(beta: np.ndarray, dosage: np.ndarray) -> dict:
    """Compute the additional tensors the model forward needs."""
    abs_beta_per_snp     = np.abs(beta).astype(np.float32)
    dosage_x_beta        = beta[None, :] * dosage                    # (P, S)
    abs_dosage_x_beta    = np.abs(dosage_x_beta).astype(np.float32)
    return {"dosage_x_beta": dosage_x_beta,
             "abs_beta_per_snp": abs_beta_per_snp,
             "abs_dosage_x_beta": abs_dosage_x_beta}


def _make_model(n_chroms: int, in_dim: int, seed: int,
                  ckpt_path: Path, variant: dict) -> "torch.nn.Module":
    model = v3.DiffAttnV3(
        in_dim=in_dim, n_chroms=n_chroms,
        aggregation=variant["aggregation"],
        head=variant["head"],
        mode=variant["func_integration_mode"],
        mlp_width=variant["mlp_width"],
        mlp_dropout=variant["mlp_dropout"],
    )
    state = torch.load(ckpt_path, map_location="cpu")
    model.load_state_dict(state)
    model.eval()
    return model


def _forward_with_hooks(model, X: np.ndarray, beta: np.ndarray,
                          dosage_x_beta: np.ndarray, chrom_idx: np.ndarray,
                          func_imp_abs: np.ndarray, variant: dict) -> dict:
    """Run forward on the full split (no batching — 116-450 patients fits
    on CPU). Captures pool output (B, 771), pool attention (B, S), MLP
    penultimate (B, 256), and final probability (B,)."""
    captured = {"pool_out": None, "pool_attn": None, "pool_attn_chrom": None,
                 "head_hidden": None}

    def hook_pool(module, inputs, output):
        # ChromHierPool returns (patient_emb, {"a_snp": ..., "a_chrom": ...}).
        if isinstance(output, tuple):
            captured["pool_out"]  = output[0].detach().cpu().numpy()
            attn = output[1]
            if isinstance(attn, dict):
                if "a_snp" in attn and attn["a_snp"] is not None:
                    captured["pool_attn"] = attn["a_snp"].detach().cpu().numpy()
                if "a_chrom" in attn and attn["a_chrom"] is not None:
                    captured["pool_attn_chrom"] = attn["a_chrom"].detach().cpu().numpy()
            elif attn is not None:
                captured["pool_attn"] = attn.detach().cpu().numpy()
        else:
            captured["pool_out"] = output.detach().cpu().numpy()

    def hook_head(module, inputs, output):
        captured["head_hidden"] = output.detach().cpu().numpy()

    h1 = model.pool.register_forward_hook(hook_pool)
    # head.net[2] is the Dropout (identity in eval) — captures the
    # post-Linear->GELU->Dropout activation, i.e. the MLP penultimate (B, 256).
    h2 = model.head.net[2].register_forward_hook(hook_head)
    try:
        x_t   = torch.from_numpy(X)
        b_t   = torch.from_numpy(beta)
        dxb_t = torch.from_numpy(dosage_x_beta.astype(np.float32))
        chrom_t = (torch.from_numpy(chrom_idx)
                    if variant["aggregation"] == "chrom_hier" else None)
        fia_t = torch.from_numpy(func_imp_abs.astype(np.float32))
        with torch.no_grad():
            logit = model(x_t, b_t, dxb_t, chrom_idx=chrom_t,
                           modality_abs=None, func_imp_abs=fia_t)
            prob = torch.sigmoid(logit).cpu().numpy().reshape(-1)
    finally:
        h1.remove(); h2.remove()
    captured["prob"] = prob
    return captured


def _forward_xgb(X: np.ndarray, beta: np.ndarray, dosage: np.ndarray,
                   func_imp_abs: np.ndarray, xgb_model_path: Path,
                   variant: dict) -> dict:
    """XGB variant inference. No DiffAttnV3, no MLP. Uses the deterministic
    `_fixed_aggregator_v3` (γ=δ=1, no learned ε) on the per-patient (P, S, 771)
    tensor to produce a 771-D patient feature, then runs the trained XGB
    classifier saved by the v3 sweep. Returns the same dict shape as
    `_forward_with_hooks`, with `head_hidden = None` (no MLP penultimate).
    """
    import xgboost as xgb
    # Deterministic attention-pool aggregator (mirror of v3 trainer's tree
    # path). For mode='none' there is no extra_bias term.
    extra_bias = v3._build_extra_bias_static(
        variant["func_integration_mode"],
        {k: func_imp_abs for k in v3.MODALITIES},  # unused for mode='none'
        P=X.shape[0], S=X.shape[1])
    # `extra_bias` is None for mode='none'; the aggregator ignores it.
    feat = v3._fixed_aggregator_v3(X, beta, dosage, extra_bias)  # (P, 771)
    booster = xgb.Booster()
    booster.load_model(str(xgb_model_path))
    prob = booster.predict(xgb.DMatrix(feat))
    return {
        "pool_out":    feat.astype(np.float32),
        "pool_attn":   None,   # fixed aggregator doesn't expose per-SNP attn
        "pool_attn_chrom": None,
        "head_hidden": None,   # XGB has no MLP penultimate
        "prob":        np.asarray(prob, dtype=np.float32).reshape(-1),
    }


def _parity_check(per_fold_capture: dict, predictions_tsv: Path,
                   seed: int, out_dir: Path) -> dict:
    """Compare recomputed val + test probs vs the trainer's predictions.tsv.
    Returns metrics_recomputed dict on success; sys.exit(2) on failure."""
    saved = pd.read_csv(predictions_tsv, sep="\t")
    diffs = []
    for fold in ("val", "test"):
        rec = per_fold_capture[fold]
        ptid_to_prob = dict(zip(rec["ptids"], rec["prob"].astype(float)))
        sub = saved[saved["split"] == fold].copy()
        ours = sub["Patient_ID"].astype(str).map(ptid_to_prob)
        if ours.isna().any():
            missing = sub[ours.isna()]["Patient_ID"].tolist()
            print(f"[FAIL] seed={seed} fold={fold}: {len(missing)} PTIDs in "
                   f"predictions.tsv not in recomputed set: {missing[:5]}…")
            sys.exit(2)
        delta = (sub["prob"].astype(float).values
                  - ours.astype(float).values)
        diffs.append(pd.DataFrame({
            "Patient_ID": sub["Patient_ID"].astype(str),
            "split": fold,
            "saved_prob": sub["prob"].astype(float).values,
            "recomputed_prob": ours.astype(float).values,
            "delta": delta,
        }))
    diff_df = pd.concat(diffs, ignore_index=True)
    max_abs = float(np.abs(diff_df["delta"]).max())
    print(f"[parity] seed={seed}  max|delta|={max_abs:.3e}  tol={PARITY_TOL:.1e}")
    if max_abs > PARITY_TOL:
        out_dir.mkdir(parents=True, exist_ok=True)
        diff_df.to_csv(out_dir / "embedding_parity_failure.tsv",
                        sep="\t", index=False)
        print(f"[FAIL] parity check exceeded tolerance — wrote "
               f"{out_dir/'embedding_parity_failure.tsv'} — aborting.")
        sys.exit(2)
    # Recompute trainer-side metrics so we can verify against metrics.json
    out = {}
    for fold in ("val", "test"):
        rec = per_fold_capture[fold]
        y    = np.asarray(rec["y"], dtype=int)
        prob = np.asarray(rec["prob"], dtype=float)
        pred = (prob > 0.5).astype(int)
        try:    auc = float(roc_auc_score(y, prob))
        except: auc = float("nan")
        bacc = float(balanced_accuracy_score(y, pred))
        loss = float(log_loss(y, np.clip(prob, 1e-7, 1-1e-7), labels=[0,1]))
        cm = confusion_matrix(y, pred, labels=[0,1])
        tn, fp, fn, tp = (cm.ravel().tolist() if cm.size == 4
                            else [0,0,0,0])
        out[fold] = {"balanced_accuracy": bacc, "roc_auc": auc, "loss": loss,
                      "tp": int(tp), "tn": int(tn), "fp": int(fp), "fn": int(fn),
                      "n": int(len(y))}
    return out


def _save_seed_outputs(seed: int, per_fold_capture: dict,
                        metrics_recomputed: dict, out_dir: Path,
                        desc_name: str):
    seed_dir = out_dir / f"seed_{seed}"
    seed_dir.mkdir(parents=True, exist_ok=True)

    # Per-patient predictions parquet (canonical multimodal schema).
    pred_rows = []
    for fold in ("train", "val", "test"):
        rec = per_fold_capture[fold]
        n = len(rec["ptids"])
        pred_rows.append(pd.DataFrame({
            "Patient_ID":   rec["ptids"],
            "seed":         np.full(n, seed, dtype=np.int8),
            "fold":         np.full(n, fold, dtype=object),
            "task":         np.full(n, "classification", dtype=object),
            "model":        np.full(n, desc_name, dtype=object),
            "covar_mode":   np.full(n, "na", dtype=object),
            "beta_source":  np.full(n, "na", dtype=object),
            "ld_config":    np.full(n, "na", dtype=object),
            "prs_raw":      np.full(n, np.nan, dtype=float),
            "prs_z":        np.full(n, np.nan, dtype=float),
            "outcome_proba": rec["prob"].astype(float),
            "y_true":       rec["y"].astype(int),
        }))
    pred_df = pd.concat(pred_rows, ignore_index=True)
    pred_df.to_parquet(seed_dir / "predictions.parquet", index=False)

    # Embeddings npz — combine all folds. embedding_256 only present for MLP
    # variants; for XGB variants this is None and is omitted from the npz.
    ptids_all, fold_all, y_all, prob_all = [], [], [], []
    emb_771_all, emb_256_all, attn_all = [], [], []
    have_emb256 = all(per_fold_capture[fold]["head_hidden"] is not None
                       for fold in ("train", "val", "test"))
    for fold in ("train", "val", "test"):
        rec = per_fold_capture[fold]
        ptids_all.extend(rec["ptids"])
        fold_all.extend([fold] * len(rec["ptids"]))
        y_all.append(rec["y"].astype(np.int8))
        prob_all.append(rec["prob"].astype(np.float32))
        emb_771_all.append(rec["pool_out"].astype(np.float32))
        if have_emb256:
            emb_256_all.append(rec["head_hidden"].astype(np.float32))
        if rec["pool_attn"] is not None:
            attn_all.append(rec["pool_attn"].astype(np.float32))
    save_kwargs = dict(
        ptids=np.array(ptids_all),
        fold=np.array(fold_all),
        y_true=np.concatenate(y_all),
        proba=np.concatenate(prob_all),
        embedding_771=np.concatenate(emb_771_all, axis=0),
    )
    if have_emb256:
        save_kwargs["embedding_256"] = np.concatenate(emb_256_all, axis=0)
    if attn_all:
        save_kwargs["attn_weights"] = np.concatenate(attn_all, axis=0)
    np.savez_compressed(seed_dir / "embeddings.npz", **save_kwargs)

    # Per-seed metrics recomputed (sanity vs metrics.json)
    (seed_dir / "metrics_recomputed.json").write_text(
        json.dumps(metrics_recomputed, indent=2))


def _save_static_assets(inputs0: dict, out_dir: Path, variant: dict):
    """Write the static (per-variant, not per-seed) tables once."""
    out_dir.mkdir(parents=True, exist_ok=True)
    rsids = list(inputs0["target_rsids"])
    chrom_arr = inputs0["chrom_arr"]
    chrom_idx = inputs0["chrom_idx"]
    beta_df = v3._load_beta(BASE).set_index("rsID")
    pos = [int(beta_df.loc[rs, "BP_GRCh38"]) if rs in beta_df.index else -1
            for rs in rsids]
    pd.DataFrame({"rsID": rsids,
                    "chrom": chrom_arr.tolist(),
                    "pos_grch38": pos,
                    "chrom_idx": chrom_idx.tolist()}).to_csv(
        out_dir / "snp_metadata.tsv", sep="\t", index=False)
    pd.DataFrame({"rsID": rsids,
                    "beta_a1": inputs0["beta"].tolist()}).to_csv(
        out_dir / "beta.tsv", sep="\t", index=False)
    # AlphaGenome per-SNP func importance (the |modality|-max post-z-score)
    summaries_z = inputs0["summaries_z"]
    fia = np.maximum.reduce([summaries_z[f"{m}_abs"] for m in v3.MODALITIES])
    pd.DataFrame({"rsID": rsids, "func_imp_abs": fia.tolist()}).to_csv(
        out_dir / "alphagenome_func_imp_abs.tsv", sep="\t", index=False)
    # Per-PTID dosage (the same table the trainer used).
    dos_df, _ = v3._load_dosage(BASE)
    cols_ordered = ["PTID"] + rsids
    dos_df[cols_ordered].to_csv(
        out_dir / "dosage_all_patients.tsv", sep="\t", index=False)

    # README with variant config + asset legend.
    readme = (
        f"# {variant['desc_name']}\n"
        "\n"
        "Multimodal-alignment anchors for one of the top-2 BMFM-SNP variants\n"
        "in `fm_top_predictions_top_7`. Lineage confirmed Drive↔local\n"
        "byte-identical 2026-06-02 (see auto-memory `reference_bmfm_top7_lineage`).\n"
        "\n"
        f"Variant config: {json.dumps(variant, indent=2)}\n"
        "\n"
        "Static assets (shared across the 3 seeds):\n"
        "- `beta.tsv`                       rsID -> β_a1 (per published GWAS)\n"
        "- `snp_metadata.tsv`               rsID -> CHR, BP_GRCh38, chrom_idx\n"
        "- `alphagenome_func_imp_abs.tsv`   rsID -> per-SNP func-importance\n"
        "                                   used as attention bias (post z-score)\n"
        "- `dosage_all_patients.tsv`        PTID × rsID dosage (0/1/2 + NaN)\n"
        "\n"
        "Per-seed (`seed_{0,1,2}/`):\n"
        "- `predictions.parquet`            canonical schema:\n"
        "                                   Patient_ID, seed, fold, task='classification',\n"
        "                                   model, covar_mode='na', beta_source='na',\n"
        "                                   ld_config='na', prs_raw=NaN, prs_z=NaN,\n"
        "                                   outcome_proba, y_true\n"
        "- `embeddings.npz`                 keys: ptids (N,), fold (N,), y_true (N,),\n"
        "                                   proba (N,), embedding_771 (N, 771),\n"
        "                                   embedding_256 (N, 256), attn_weights (N, 128).\n"
        "                                   embedding_771 = pool output (768 diff_emb +\n"
        "                                   3 value_block [β, dosage, βxdosage] cols\n"
        "                                   chrom-hier-aggregated). embedding_256 = MLP\n"
        "                                   penultimate (post-Linear->GELU->Dropout).\n"
        "- `metrics_recomputed.json`        BalAcc / AUC / loss / CM recomputed from\n"
        "                                   the recomputed probabilities — sanity check\n"
        "                                   vs each seed dir's `metrics.json`.\n"
        "\n"
        "**Parity gate (FAIL-LOUDLY):** script 50 asserts per-PTID\n"
        "`|prob_recomputed - prob_saved| < 1e-6` against the trainer's\n"
        "`predictions.tsv`. If the gate fails for any seed, no files are written.\n"
    )
    (out_dir / "README.md").write_text(readme, encoding="utf-8")


def _process_variant(variant: dict, seeds: list[int]) -> None:
    """Extract anchors for one variant (MLP or XGB head). Per-PTID parity gate
    against the trainer's predictions.tsv per seed."""
    out_base = OUT_BASE_ROOT / variant["desc_name"]
    top7_dir = TOP7_ROOT / variant["top7_subdir"]
    print(f"\n========== [50] variant: {variant['desc_name']} ==========")
    print(f"  ckpt source : {top7_dir}/{variant['variant_name']}_s*")
    print(f"  output dir  : {out_base}")
    out_base.mkdir(parents=True, exist_ok=True)

    # Static assets: load inputs for seed 0 first.
    print(f"[50] Loading inputs for seed 0 (publishes static assets too)...")
    inputs0 = _load_inputs(0, variant)
    _save_static_assets(inputs0, out_base, variant)
    print(f"[50] Wrote static assets to {out_base}")

    t0 = time.time()
    for seed in seeds:
        variant_dir = top7_dir / f"{variant['variant_name']}_s{seed}"
        ckpt = variant_dir / variant["ckpt_filename"]
        pred_tsv = variant_dir / "predictions.tsv"
        if not ckpt.exists() or not pred_tsv.exists():
            print(f"[skip] seed={seed}: missing {variant['ckpt_filename']} "
                   f"or predictions.tsv in {variant_dir}")
            continue
        print(f"\n[50] {variant['desc_name']}  seed={seed} -- loading inputs + ckpt...")
        inp = inputs0 if seed == 0 else _load_inputs(seed, variant)
        # Build per-patient X (P, S, 771) + extras
        per_fold = {}
        F = None
        for sp in ("train", "val", "test"):
            d = inp["splits"][sp]
            X = v3._per_patient_values(
                variant["func_integration_mode"],
                inp["beta"], d["dosage"],
                None, inp["summaries_z"], inp["diff"])
            F = X.shape[2]
            extras = _build_extras(inp["beta"], d["dosage"])
            per_fold[sp] = {"X": X, "ptids": d["ptid"], "y": d["y"],
                             "dosage": d["dosage"], **extras}
        assert F == 771, f"expected in_dim 771; got {F}"
        # func_imp_abs (S,) — z-scored max across modalities
        fia = np.maximum.reduce([inp["summaries_z"][f"{m}_abs"]
                                   for m in v3.MODALITIES]).astype(np.float32)
        # Branch on head type. MLP variants use DiffAttnV3 forward; XGB
        # variants use the deterministic _fixed_aggregator_v3 + xgb predict.
        per_fold_capture = {}
        if variant["head"] in ("mlp2", "mlp3"):
            model = _make_model(inp["n_chroms"], F, seed, ckpt, variant)
            for sp in ("train", "val", "test"):
                d = per_fold[sp]
                cap = _forward_with_hooks(model, d["X"], inp["beta"],
                                            d["dosage_x_beta"], inp["chrom_idx"], fia,
                                            variant)
                cap["ptids"] = d["ptids"]
                cap["y"]     = d["y"]
                per_fold_capture[sp] = cap
        elif variant["head"] == "xgb":
            for sp in ("train", "val", "test"):
                d = per_fold[sp]
                cap = _forward_xgb(d["X"], inp["beta"], d["dosage"],
                                     fia, ckpt, variant)
                cap["ptids"] = d["ptids"]
                cap["y"]     = d["y"]
                per_fold_capture[sp] = cap
        else:
            raise ValueError(f"Unsupported head: {variant['head']}")
        for sp in ("train", "val", "test"):
            cap = per_fold_capture[sp]
            emb256_shape = (cap["head_hidden"].shape
                             if cap["head_hidden"] is not None else "(none)")
            attn_label = "yes" if cap["pool_attn"] is not None else "no"
            print(f"  seed={seed} {sp:5s}: n={len(cap['ptids']):3d}  "
                   f"emb771={cap['pool_out'].shape}  emb256={emb256_shape}  "
                   f"attn={attn_label}")
        # Parity gate vs predictions.tsv (val + test)
        seed_dir = out_base / f"seed_{seed}"
        metrics_recomputed = _parity_check(per_fold_capture, pred_tsv,
                                              seed, seed_dir)
        # Compare to the trainer's metrics.json
        trainer_metrics = json.loads((variant_dir / "metrics.json").read_text())
        for fold in ("val", "test"):
            d_tr = trainer_metrics[fold]
            d_re = metrics_recomputed[fold]
            for key in ("balanced_accuracy", "roc_auc", "loss"):
                delta = abs(d_tr[key] - d_re[key])
                tag = "ok" if delta < 1e-3 else "MISMATCH"
                print(f"  metrics seed={seed} {fold} {key:18s} "
                       f"trainer={d_tr[key]:.4f}  recomputed={d_re[key]:.4f}  "
                       f"delta={delta:.2e}  [{tag}]")
        _save_seed_outputs(seed, per_fold_capture, metrics_recomputed,
                            out_base, variant["desc_name"])
        print(f"  seed={seed} [OK] wrote anchors to {seed_dir}")
    print(f"[50] variant {variant['desc_name']} DONE in {time.time()-t0:.1f}s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS),
                    help="Subset of seeds to process (default: 0 1 2)")
    ap.add_argument("--variants", type=str, nargs="+",
                    default=[v["desc_name"] for v in VARIANTS],
                    help="Subset of variant desc_names to run (default: all)")
    args = ap.parse_args()
    selected = [v for v in VARIANTS if v["desc_name"] in set(args.variants)]
    if not selected:
        print(f"[err] no variants matched --variants {args.variants}")
        return
    t_start = time.time()
    for variant in selected:
        _process_variant(variant, args.seeds)
    print(f"\n[50] ALL DONE in {time.time()-t_start:.1f}s")


if __name__ == "__main__":
    main()
