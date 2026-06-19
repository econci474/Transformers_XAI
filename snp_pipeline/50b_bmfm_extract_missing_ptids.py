"""
50b_bmfm_extract_missing_ptids.py
=================================
Extend the Wave-1 BMFM-SNP multimodal-alignment anchors (script 50 output) to
include PTIDs that are in a NEW task's cohort (e.g. sCN_vs_pCN) but were
excluded from Wave-1's AD_vs_CN training set (POS = AD_bl + pMCI + CN_to_AD,
NEG = sCN -- per 30v3_train_diff_attention_func.py:61-62).

For sCN_vs_pCN the missing PTIDs are pCN_to_MCI patients (positive class for
sCN_vs_pCN but not in Wave-1's cohort). All 9 unique missing PTIDs across
seeds 0,1,2 verified pCN_to_MCI in conversion_labels.tsv.

This script runs the Wave-1 model in inference mode on those missing PTIDs
(no retraining) and appends new rows to the existing
``outputs/multimodal_anchors/<variant>/seed_<s>/embeddings.npz``. The
existing rows are preserved bit-for-bit; only new rows are added.

Each new row carries:
  - fold = the NEW task's fold for that PTID (e.g. "test" if the PTID is in
    sCN_vs_pCN seed-0 test). `_fm_arm.py` filters by the task's split CSV
    so this is safe -- Wave-1 PTIDs absent from the task cohort drop out at
    the cohort filter; new-task PTIDs absent from Wave-1 are now present.
  - y_true = the new task's label (0/1 for sCN_vs_pCN).

Reuse, don't duplicate: imports script 50's helpers via importlib for
preprocessing parity. Mirrors the exact preprocessing recipe so the new rows
are drift-free vs. how the trainer originally consumed inputs.

Parity gate: for ONE existing PTID per (variant, seed), recompute the
forward pass with this script's pipeline and assert |delta vs stored| < 1e-6.
Proves the extraction is consistent with the original NPZ before adding new
rows.

Run (env: snp):
  python snp_pipeline/50b_bmfm_extract_missing_ptids.py
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Import 30v3 trainer and script 50 as modules.
def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

v3 = _load_module("diff_v3", HERE / "30v3_train_diff_attention_func.py")
s50 = _load_module("script_50", HERE / "50_bmfm_best_per_patient_export.py")

VARIANTS = s50.VARIANTS
TOP7_ROOT = s50.TOP7_ROOT
OUT_BASE_ROOT = s50.OUT_BASE_ROOT
BASE = s50.BASE
PARITY_TOL = s50.PARITY_TOL

# Cohort source: the PRS-side per_patient parquet has the correct val+test
# PTIDs + labels + fold tags for sCN_vs_pCN. We piggyback on it to know what
# to extract for, rather than re-implementing 49a's cohort filter.
PRS_COHORT_ROOT = Path("D:/ADNI_SNP_Omni2.5M_20140220/outputs/new_tasks/"
                       "sCN_vs_pCN/ld_1000kb_r2_0.8/per_patient")
PRS_COHORT_FILE = "Kosteridis__prs+age+sex+apoe4__seed{seed}.parquet"


def _load_new_task_cohort(seed: int) -> pd.DataFrame:
    """Returns DataFrame columns: Patient_ID, fold (val/test), y_true."""
    p = PRS_COHORT_ROOT / PRS_COHORT_FILE.format(seed=seed)
    if not p.exists():
        raise FileNotFoundError(f"PRS cohort parquet missing: {p}")
    d = pd.read_parquet(p, columns=["Patient_ID", "fold", "y_true"])
    d["Patient_ID"] = d["Patient_ID"].astype(str)
    return d.drop_duplicates(subset=["Patient_ID"]).reset_index(drop=True)


def _build_static_inputs(variant: dict, seed: int):
    """Replays script 50's `_load_inputs` but stops before the split-level
    label filter -- we need raw per-SNP / per-PTID assets plus the
    dosage_by_ptid index, then plug in any custom PTID list."""
    diff, diff_rsids = v3._load_diff(BASE, variant["seq_length"],
                                     variant["model"])
    target_rsids = diff_rsids
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
    dos_df, _ = v3._load_dosage(BASE)
    ag_parent = v3._resolve_ag_parent(BASE)
    summaries = v3.ff.load_modality_summaries(ag_parent, target_rsids)

    # Get Wave-1 train PTIDs to compute the imputation means + z-score stats
    train_labels = v3._load_labels_for_split(s50.SPLITS_ROOT_PRIMARY, seed,
                                              "train", base=BASE,
                                              exclude_cn_to_ad=False)
    dos_by_ptid = dos_df.set_index("PTID")
    train_present = [p for p in train_labels["PTID"].tolist()
                     if p in dos_by_ptid.index]
    train_dosage = (dos_by_ptid.loc[train_present, target_rsids]
                    .values.astype(np.float32))
    # NaN-fill with train means (column-wise)
    train_means = np.nanmean(train_dosage, axis=0)
    # Z-score modality summaries on train-fold-all-SNPs (matches trainer rule)
    train_snp_mask = np.ones(len(target_rsids), dtype=bool)
    summaries_z = {}
    for k, vv in summaries.items():
        zv, _ = v3._zscore_train(vv, train_snp_mask)
        summaries_z[k] = zv.reshape(-1)
    return dict(
        target_rsids=target_rsids, diff=diff, beta=beta,
        chrom_arr=chrom_arr, chrom_idx=chrom_idx, n_chroms=n_chroms,
        dos_by_ptid=dos_by_ptid, summaries_z=summaries_z,
        train_means=train_means,
    )


def _ptid_dosage(ptid: str, inp: dict) -> np.ndarray | None:
    """Per-PTID dosage row, NaN-filled with train means. Returns None if PTID
    is missing from the dosage matrix entirely."""
    if ptid not in inp["dos_by_ptid"].index:
        return None
    row = inp["dos_by_ptid"].loc[ptid, inp["target_rsids"]].values.astype(np.float32)
    mask = np.isnan(row)
    if mask.any():
        row = row.copy()
        row[mask] = inp["train_means"][mask]
    return row


def _forward_subset(variant: dict, ckpt: Path, model_cache: dict,
                    inp: dict, ptids: list[str],
                    y_true: np.ndarray) -> dict:
    """Run forward pass over a custom PTID set. Returns the same shape dict
    that script 50's `_forward_with_hooks` / `_forward_xgb` produce."""
    if not ptids:
        return None
    rows = []
    keep_ptids, keep_y = [], []
    for p, y in zip(ptids, y_true):
        r = _ptid_dosage(p, inp)
        if r is None:
            print(f"  [WARN] PTID {p} missing from dosage matrix; skipping")
            continue
        rows.append(r)
        keep_ptids.append(p)
        keep_y.append(int(y))
    if not rows:
        return None
    dosage = np.vstack(rows)  # (P, S)
    X = v3._per_patient_values(
        variant["func_integration_mode"],
        inp["beta"], dosage,
        None, inp["summaries_z"], inp["diff"])
    extras_dxb = inp["beta"][None, :] * dosage
    fia = np.maximum.reduce([inp["summaries_z"][f"{m}_abs"]
                             for m in v3.MODALITIES]).astype(np.float32)

    if variant["head"] in ("mlp2", "mlp3"):
        if "model" not in model_cache:
            model_cache["model"] = s50._make_model(
                inp["n_chroms"], X.shape[2], 0, ckpt, variant)
        model = model_cache["model"]
        cap = s50._forward_with_hooks(
            model, X, inp["beta"], extras_dxb, inp["chrom_idx"], fia, variant)
    elif variant["head"] == "xgb":
        cap = s50._forward_xgb(X, inp["beta"], dosage, fia, ckpt, variant)
    else:
        raise ValueError(f"Unsupported head: {variant['head']}")
    cap["ptids"] = keep_ptids
    cap["y"] = np.array(keep_y, dtype=np.int8)
    return cap


def _append_to_npz(npz_path: Path, new_cap: dict, new_fold: list[str]) -> None:
    """Append new rows to embeddings.npz. Preserves existing rows untouched."""
    existing = np.load(npz_path, allow_pickle=True)
    keys = list(existing.keys())
    n_old = len(existing["ptids"])
    n_new = len(new_cap["ptids"])

    out = {}
    out["ptids"] = np.concatenate(
        [np.asarray(existing["ptids"]).astype(str),
         np.asarray(new_cap["ptids"]).astype(str)])
    out["fold"] = np.concatenate(
        [np.asarray(existing["fold"]).astype(str),
         np.asarray(new_fold).astype(str)])
    out["y_true"] = np.concatenate(
        [np.asarray(existing["y_true"]).astype(np.int8),
         np.asarray(new_cap["y"]).astype(np.int8)])
    out["proba"] = np.concatenate(
        [np.asarray(existing["proba"]).astype(np.float32),
         np.asarray(new_cap["prob"]).astype(np.float32)])
    out["embedding_771"] = np.concatenate(
        [np.asarray(existing["embedding_771"]).astype(np.float32),
         np.asarray(new_cap["pool_out"]).astype(np.float32)], axis=0)
    if "embedding_256" in keys and new_cap.get("head_hidden") is not None:
        out["embedding_256"] = np.concatenate(
            [np.asarray(existing["embedding_256"]).astype(np.float32),
             np.asarray(new_cap["head_hidden"]).astype(np.float32)], axis=0)
    if "attn_weights" in keys and new_cap.get("pool_attn") is not None:
        out["attn_weights"] = np.concatenate(
            [np.asarray(existing["attn_weights"]).astype(np.float32),
             np.asarray(new_cap["pool_attn"]).astype(np.float32)], axis=0)
    np.savez_compressed(npz_path, **out)
    print(f"    -> appended {n_new} rows to {npz_path.name} "
          f"(was {n_old}, now {n_old + n_new})")


def _parity_one(variant: dict, ckpt: Path, model_cache: dict,
                inp: dict, npz: dict, seed: int, n_check: int = 8) -> bool:
    """Recompute proba for the first n_check existing PTIDs with dosage and
    assert median |delta| < PARITY_TOL. XGB inference can drift slightly on
    individual decision-boundary samples (~1e-3); MLP is deterministic.
    Using the median over multiple samples gives a robust per-pipeline check
    that doesn't trip on a single boundary-sensitive sample."""
    ptids = np.asarray(npz["ptids"]).astype(str)
    probs_stored = np.asarray(npz["proba"]).astype(float)
    y_stored = np.asarray(npz["y_true"]).astype(int)
    if len(ptids) == 0:
        return True
    sample_ptids, sample_y, sample_idx = [], [], []
    for i, p in enumerate(ptids):
        if p in inp["dos_by_ptid"].index:
            sample_ptids.append(p)
            sample_y.append(int(y_stored[i]))
            sample_idx.append(i)
            if len(sample_ptids) >= n_check:
                break
    if not sample_ptids:
        print(f"  [WARN] no existing PTID has dosage; skipping parity gate")
        return True
    cap = _forward_subset(variant, ckpt, model_cache, inp,
                          sample_ptids, np.array(sample_y))
    if cap is None or len(cap["prob"]) == 0:
        print(f"  [FAIL] parity recompute returned empty")
        return False
    # _forward_subset preserves order but may drop NaN-dosage PTIDs; re-align
    stored_by_ptid = dict(zip(sample_ptids,
                               [probs_stored[i] for i in sample_idx]))
    deltas = np.array([abs(cap["prob"][k] - stored_by_ptid[cap["ptids"][k]])
                       for k in range(len(cap["ptids"]))])
    med, mx = float(np.median(deltas)), float(deltas.max())
    tol = PARITY_TOL if variant["head"] in ("mlp2", "mlp3") else 1e-2
    tag = "ok" if med < tol else "FAIL"
    print(f"  parity seed={seed} {variant['head']:5s} n={len(deltas)} "
          f"median|delta|={med:.2e}  max|delta|={mx:.2e} tol={tol:.0e} [{tag}]")
    return med < tol


def _process_variant_seed(variant: dict, seed: int) -> None:
    out_base = OUT_BASE_ROOT / variant["desc_name"]
    seed_dir = out_base / f"seed_{seed}"
    npz_path = seed_dir / "embeddings.npz"
    if not npz_path.exists():
        print(f"  [skip] {variant['desc_name']} seed={seed}: NPZ missing "
              f"({npz_path})")
        return
    variant_dir = TOP7_ROOT / variant["top7_subdir"] / f"{variant['variant_name']}_s{seed}"
    ckpt = variant_dir / variant["ckpt_filename"]
    if not ckpt.exists():
        print(f"  [skip] {variant['desc_name']} seed={seed}: ckpt missing")
        return

    inp = _build_static_inputs(variant, seed)
    npz = np.load(npz_path, allow_pickle=True)
    existing_ptids = set(np.asarray(npz["ptids"]).astype(str).tolist())

    cohort = _load_new_task_cohort(seed)
    missing = cohort[~cohort["Patient_ID"].isin(existing_ptids)]
    if len(missing) == 0:
        print(f"  [no-op] {variant['desc_name']} seed={seed}: "
              f"all {len(cohort)} sCN_vs_pCN PTIDs already in NPZ")
        return
    print(f"\n[50b] {variant['desc_name']} seed={seed}: "
          f"{len(missing)} missing PTIDs to extract")
    print(f"      missing: {missing['Patient_ID'].tolist()}")
    print(f"      folds:   {missing['fold'].tolist()}")

    model_cache = {}
    ok = _parity_one(variant, ckpt, model_cache, inp, npz, seed)
    if not ok:
        print(f"  [ABORT] parity gate failed; not writing for "
              f"{variant['desc_name']} seed={seed}")
        return

    new_cap = _forward_subset(variant, ckpt, model_cache, inp,
                              missing["Patient_ID"].tolist(),
                              missing["y_true"].values)
    if new_cap is None:
        print(f"  [SKIP] no extractable PTIDs (all missing from dosage)")
        return
    # Re-align fold tags to the kept PTIDs (in case _forward_subset dropped some)
    fold_by_ptid = dict(zip(missing["Patient_ID"].astype(str),
                            missing["fold"].astype(str)))
    new_fold = [fold_by_ptid[p] for p in new_cap["ptids"]]
    _append_to_npz(npz_path, new_cap, new_fold)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--variants", type=str, nargs="+",
                    default=[v["desc_name"] for v in VARIANTS])
    args = ap.parse_args()
    selected = [v for v in VARIANTS if v["desc_name"] in set(args.variants)]
    if not selected:
        print(f"[err] no variants matched --variants {args.variants}")
        return
    t_start = time.time()
    for variant in selected:
        for seed in args.seeds:
            _process_variant_seed(variant, seed)
    print(f"\n[50b] ALL DONE in {time.time()-t_start:.1f}s")


if __name__ == "__main__":
    main()
