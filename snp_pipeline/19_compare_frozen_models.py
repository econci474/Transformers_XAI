"""
19_compare_frozen_models.py
============================
Frozen-embedding patient-classification bake-off across foundation models.

Compares per-patient frozen embeddings from several DNA foundation models on
the binary ``ever_convert`` task (stable-CN vs ever-MCI/AD), using a SINGLE
common per-window pool (``mean``) and aggregation, and FOUR classifier heads
to see which best extracts signal:

    2-layer MLP, 3-layer MLP, RandomForest, XGBoost

Aggregation design (see plan / README):
  * ``chrom_mean_then_mean`` — parameter-free per-chrom mean then mean. The
    SAME fixed per-patient vector feeds all 4 heads (apples-to-apples).
  * ``chrom_mean_then_attn`` — trainable per-chrom attention; only valid for
    the differentiable MLP heads (RF/XGBoost cannot backprop through it), so
    it is run for ``mlp2``/``mlp3`` only.

Per-window pooling is fixed to ``mean`` for ALL models (the only pool common
to BMFM/NT/Caduceus); it is not a swept axis.

All heavy lifting is reused from ``15_train_classification_head.py`` (imported
by file path because its name starts with a digit) — that script is NOT
modified. New here: ``MLP3Head``, sklearn RF/XGBoost heads, a model registry,
and a hard preflight that guards the BMFM window-ordering caveat.

Output tree (script-16 compatible — also writes ``upstream_tag``):

    <output-root>/<label_mode>/<model>/mean/<aggregation>/<head>/seed_<N>/metrics.json

Usage
-----
    python snp_pipeline/19_compare_frozen_models.py \\
        --output-root frozen_model_comparison \\
        --models all --heads all --seeds 0 1 2 --device cpu

    # fast smoke (1 model x xgb x seed0 x chrom_mean_then_mean)
    python snp_pipeline/19_compare_frozen_models.py --smoke \\
        --models bmfm_base_ref --heads xgb --seeds 0 --device cpu
"""
from __future__ import annotations

import argparse
import importlib.util as _il
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    balanced_accuracy_score, f1_score, roc_auc_score,
    recall_score, precision_score, confusion_matrix,
)

# ── Reuse script 15 (digit-leading filename → import by path) ───────────────
_SCRIPT15 = Path(__file__).parent / "15_train_classification_head.py"
_spec = _il.spec_from_file_location("script15", _SCRIPT15)
_s15 = _il.module_from_spec(_spec)
_spec.loader.exec_module(_s15)
build_chrom_map      = _s15.build_chrom_map
load_embeddings      = _s15.load_embeddings
load_labels          = _s15.load_labels
stack_patient_arrays = _s15.stack_patient_arrays
Aggregator           = _s15.Aggregator
MLPHead              = _s15.MLPHead
FullModel            = _s15.FullModel
train_one            = _s15.train_one
evaluate             = _s15.evaluate          # used for MLP heads (model,x,m,c,y,device)

# ── Paths / registry ────────────────────────────────────────────────────────
DATA      = Path("D:/ADNI_SNP_Omni2.5M_20140220")
STAGED    = DATA / "embeddings"                          # user-downloaded HPC npz
WINDOWS_8K   = DATA / "bmfm_inputs/patient_window_sequences/windows.tsv"        # post-fix (numeric)
WINDOWS_131K = DATA / "caduceus_inputs/patient_window_sequences/windows.tsv"    # ~57 windows
SPLITS = {
    "ever_convert": Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                          "no_cdr_stratified_ever_convert/tabular/baseline"),
    "bl_multi":     Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                          "no_cdr_stratified/tabular/baseline"),
}

# model -> (npz path, expected hidden dim, windows.tsv, is_8k_group)
MODEL_REGISTRY: dict[str, dict] = {
    "bmfm_base_ref":    {"npz": DATA / "bmfm_inputs/patient_embeddings/base_ref/embeddings.npz",
                          "hidden": 768,  "windows": WINDOWS_8K,   "grp8k": True},
    "ntv2":             {"npz": STAGED / "ntv2_8k/embeddings.npz",
                          "hidden": 1024, "windows": WINDOWS_8K,   "grp8k": True},
    "caduceus_ps_8k":   {"npz": STAGED / "caduceus_ps_8k/embeddings.npz",
                          "hidden": 512,  "windows": WINDOWS_8K,   "grp8k": True},
    "caduceus_ph_8k":   {"npz": STAGED / "caduceus_ph_8k/embeddings.npz",
                          "hidden": 256,  "windows": WINDOWS_8K,   "grp8k": True},
    "caduceus_ps_131k": {"npz": STAGED / "caduceus_ps_131k/embeddings.npz",
                          "hidden": 512,  "windows": WINDOWS_131K, "grp8k": False},
    "caduceus_ph_131k": {"npz": STAGED / "caduceus_ph_131k/embeddings.npz",
                          "hidden": 256,  "windows": WINDOWS_131K, "grp8k": False},
}

ALL_HEADS = ["mlp2", "mlp3", "rf", "xgb"]
# which aggregations each head may use
HEAD_AGGS = {
    "chrom_mean_then_mean": ["mlp2", "mlp3", "rf", "xgb"],   # param-free → all heads
    "chrom_mean_then_attn": ["mlp2", "mlp3"],                # trainable → MLP only
}


# ── 3-layer MLP head (only genuinely new model class) ───────────────────────
class MLP3Head(torch.nn.Module):
    """Linear→GELU→Dropout→Linear→GELU→Dropout→Linear (mirrors script 15
    MLPHead's (B,in)→(B,) contract so it drops into FullModel unchanged)."""
    def __init__(self, in_dim: int, h1: int = 256, h2: int = 128, dropout: float = 0.3):
        super().__init__()
        self.net = torch.nn.Sequential(
            torch.nn.Linear(in_dim, h1), torch.nn.GELU(), torch.nn.Dropout(dropout),
            torch.nn.Linear(h1, h2),     torch.nn.GELU(), torch.nn.Dropout(dropout),
            torch.nn.Linear(h2, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x).squeeze(-1)


# ── Metric helper for the sklearn (RF/XGB) heads — mirrors script 15
#    evaluate() exactly, but takes real-valued scores + a decision threshold
#    (MLP uses logits>0; RF/XGB use proba>0.5). Schema is identical. ──────────
def eval_metrics(y_true: np.ndarray, score: np.ndarray, thr: float) -> dict:
    y = np.asarray(y_true).astype(int)
    pred = (score > thr).astype(int)
    out = {
        "auc": float(roc_auc_score(y, score)) if len(np.unique(y)) > 1 else float("nan"),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "f1": float(f1_score(y, pred, zero_division=0)),
        "precision": float(precision_score(y, pred, pos_label=1, zero_division=0)),
        "recall": float(recall_score(y, pred, pos_label=1, zero_division=0)),
        "sensitivity": float(recall_score(y, pred, pos_label=1, zero_division=0)),
        "specificity": float(recall_score(y, pred, pos_label=0, zero_division=0)),
    }
    cm = confusion_matrix(y, pred, labels=[0, 1])
    out["confusion_matrix"] = {"tn": int(cm[0, 0]), "fp": int(cm[0, 1]),
                               "fn": int(cm[1, 0]), "tp": int(cm[1, 1])}
    return out


# ── Preflight: validate every requested model; guard the BMFM ordering caveat ─
def preflight(model_keys: list[str]) -> dict[str, dict]:
    """Returns {model: info} for models that pass. Missing npz → skipped with a
    warning (lets you run partial before the scp completes). The 8k group
    (base_ref/ntv2/caduceus_*_8k) must share an IDENTICAL window_id set with
    WINDOWS_8K — proving the alphanumeric→numeric windows.tsv fix was
    order-only, not a re-windowing (else BMFM vs NT/Caduceus is confounded)."""
    print("=" * 70)
    print("PREFLIGHT")
    print("=" * 70)
    w8k_ids = None
    if WINDOWS_8K.exists():
        w8k_ids = set(pd.read_csv(WINDOWS_8K, sep="\t")["window_id"].astype(str))

    ok: dict[str, dict] = {}
    grp8k_idsets: dict[str, set] = {}
    for mk in model_keys:
        reg = MODEL_REGISTRY[mk]
        npz = Path(reg["npz"])
        if not npz.exists():
            print(f"  [skip] {mk}: npz not found ({npz}) — stage it then re-run.")
            continue
        if not Path(reg["windows"]).exists():
            print(f"  [skip] {mk}: windows.tsv missing ({reg['windows']}).")
            continue
        try:
            emb, win_ids = load_embeddings(npz, pool="mean")
        except Exception as e:
            print(f"  [skip] {mk}: failed to load npz ({type(e).__name__}: {e}).")
            continue
        if len(emb) == 0:
            print(f"  [skip] {mk}: no 'mean__' embeddings in npz "
                  f"(needs mean-pooled per-patient arrays).")
            continue
        hid = next(iter(emb.values())).shape[1]
        if hid != reg["hidden"]:
            print(f"  [skip] {mk}: hidden dim {hid} != expected {reg['hidden']}.")
            continue
        chrom_of_window, chroms = build_chrom_map(Path(reg["windows"]))
        win_set = set(map(str, win_ids))
        missing = win_set - set(chrom_of_window)
        if missing:
            print(f"  [skip] {mk}: {len(missing)} npz window_ids absent from "
                  f"windows.tsv (e.g. {sorted(missing)[:3]}).")
            continue
        print(f"  [ok]   {mk}: {len(emb)} patients x {len(win_ids)} windows x "
              f"{hid}d  ({len(chroms)} chroms)")
        ok[mk] = {**reg, "n_windows": len(win_ids), "n_chroms": len(chroms),
                  "chrom_of_window": chrom_of_window, "win_ids": win_ids}
        if reg["grp8k"]:
            grp8k_idsets[mk] = win_set

    # 8k-group window-definition equality gate (the headline risk).
    if grp8k_idsets:
        ref = w8k_ids if w8k_ids is not None else next(iter(grp8k_idsets.values()))
        bad = [mk for mk, s in grp8k_idsets.items() if s != ref]
        if bad:
            print("\n  [HALT] 8k window-id sets disagree — the BMFM windows "
                  "ordering fix changed window DEFINITIONS, not just order:")
            for mk in bad:
                print(f"         {mk}: |Δ vs WINDOWS_8K| = "
                      f"{len(grp8k_idsets[mk] ^ ref)} window_ids")
            print("         Re-extract the offending model(s) with the "
                  "corrected 13_build_patient_window_sequences.py.\n"
                  "         (chrom_mean_then_* is order-invariant, so a pure "
                  "reorder would pass this check.)")
            for mk in bad:
                ok.pop(mk, None)
        else:
            print(f"\n  [ok] 8k group window_id sets identical "
                  f"({len(ref)} windows) — ordering fix was order-only; "
                  f"BMFM vs NT/Caduceus is a fair comparison.")
    if not ok:
        raise SystemExit("Preflight: no usable models. Stage the npz files "
                          "(see plan) and re-run.")
    print("=" * 70)
    return ok


# ── Assemble per-split tensors for one model (mirrors script 15 main) ────────
def assemble(info: dict, label_mode: str, seed: int):
    npz = Path(info["npz"])
    emb, win_ids = load_embeddings(npz, pool="mean")
    labels = load_labels(SPLITS[label_mode], seed, label_mode=label_mode)
    splits = {}
    for split, df in labels.items():
        df = df[df["Patient_ID"].isin(emb)].reset_index(drop=True)
        x, m, c = stack_patient_arrays(emb, win_ids, info["chrom_of_window"],
                                       df["Patient_ID"].tolist())
        y = torch.from_numpy(df["y"].values.astype(np.int64))
        splits[split] = (x, m, c, y)
    return splits


def fixed_features(agg_name: str, hidden: int, n_chroms: int, split):
    """Parameter-free Aggregator forward (no grad) → (n, D) numpy for RF/XGB."""
    x, m, c, y = split
    agg = Aggregator(agg_name, hidden, n_chroms).eval()
    with torch.no_grad():
        feats = agg(x, m, c).cpu().numpy()
    return feats, y.cpu().numpy().astype(int)


# ── Run one (model, label_mode, aggregation, head, seed) cell ───────────────
def run_cell(mk: str, info: dict, label_mode: str, agg: str, head: str,
             seed: int, args) -> None:
    out_dir = (Path(args.output_root) / label_mode / mk / "mean" / agg /
               head / f"seed_{seed}")
    out_dir.mkdir(parents=True, exist_ok=True)
    final = out_dir / "metrics.json"
    if final.exists() and not args.force:
        print(f"  SKIP (exists): {final}")
        return

    torch.manual_seed(seed)
    np.random.seed(seed)
    splits = assemble(info, label_mode, seed)
    hidden, n_chroms = info["hidden"], info["n_chroms"]
    device = torch.device(args.device)
    t0 = time.time()

    if head in ("mlp2", "mlp3"):
        aggregator = Aggregator(agg, hidden, n_chroms)
        head_mod = (MLPHead(aggregator.output_dim, hidden=256, dropout=0.3)
                    if head == "mlp2"
                    else MLP3Head(aggregator.output_dim))
        model = FullModel(aggregator, head_mod).to(device)
        model, fit_info, _ = train_one(
            model, splits["train"], splits["val"],
            epochs=args.epochs, batch_size=args.batch_size,
            lr=args.lr, weight_decay=args.weight_decay,
            device=device, patience=args.patience, loss_name="weighted_bce")
        val_m  = evaluate(model, *splits["val"],  device=device)
        test_m = evaluate(model, *splits["test"], device=device)
        clf_info = {"best_val_auc": fit_info["best_val_auc"],
                    "n_epochs": fit_info["n_epochs"]}
    else:
        Xtr, ytr = fixed_features(agg, hidden, n_chroms, splits["train"])
        Xva, yva = fixed_features(agg, hidden, n_chroms, splits["val"])
        Xte, yte = fixed_features(agg, hidden, n_chroms, splits["test"])
        n_pos = int((ytr == 1).sum()); n_neg = int((ytr == 0).sum())
        spw = max(n_neg / max(n_pos, 1), 1e-6)
        if head == "rf":
            best, best_auc, best_n = None, -1.0, None
            for n_est in (200, 500):
                clf = RandomForestClassifier(
                    n_estimators=n_est, class_weight="balanced",
                    random_state=seed, n_jobs=-1)
                clf.fit(Xtr, ytr)
                p = clf.predict_proba(Xva)[:, 1]
                a = (roc_auc_score(yva, p) if len(np.unique(yva)) > 1
                     else float("nan"))
                if np.nan_to_num(a, nan=-1.0) > best_auc:
                    best, best_auc, best_n = clf, a, n_est
            clf, clf_info = best, {"best_val_auc": float(best_auc),
                                   "n_estimators": best_n}
        else:  # xgb (lazy import — MLP/RF still work if xgboost absent)
            try:
                from xgboost import XGBClassifier
            except ImportError:
                print("  [skip] xgboost not installed; `pip install xgboost`.")
                return
            clf = XGBClassifier(
                n_estimators=400, max_depth=4, learning_rate=0.05,
                subsample=0.9, colsample_bytree=0.9,
                eval_metric="logloss", scale_pos_weight=spw,
                random_state=seed, verbosity=0, n_jobs=-1)
            try:
                clf.fit(Xtr, ytr, eval_set=[(Xva, yva)],
                        early_stopping_rounds=30, verbose=False)
                best_n = getattr(clf, "best_iteration", None)
            except TypeError:                       # older/newer xgboost API
                clf.fit(Xtr, ytr)
                best_n = 400
            pva = clf.predict_proba(Xva)[:, 1]
            clf_info = {"best_val_auc": float(
                roc_auc_score(yva, pva) if len(np.unique(yva)) > 1
                else float("nan")),
                "n_estimators": (int(best_n) if best_n is not None else 400)}
        val_m  = eval_metrics(yva, clf.predict_proba(Xva)[:, 1], thr=0.5)
        test_m = eval_metrics(yte, clf.predict_proba(Xte)[:, 1], thr=0.5)

    with open(final, "w") as f:
        json.dump({
            "label_mode": label_mode,
            "model": mk,
            "upstream_tag": mk,          # alias → script 16 stays compatible
            "pool": "mean",
            "aggregation": agg,
            "head": head,
            "clf": head,
            "seed": seed,
            "fit_info": clf_info,
            "val": val_m,
            "test": test_m,
        }, f, indent=2)
    print(f"  wrote {final}  "
          f"(val_auc={val_m['auc']:.3f} test_auc={test_m['auc']:.3f} "
          f"{time.time()-t0:.1f}s)")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output-root", type=Path, default=Path("frozen_model_comparison"))
    ap.add_argument("--models", nargs="+", default=["all"],
                    help="registry keys or 'all'")
    ap.add_argument("--heads", nargs="+", default=["all"],
                    help="subset of mlp2 mlp3 rf xgb, or 'all'")
    ap.add_argument("--aggregations", nargs="+",
                    default=["chrom_mean_then_mean", "chrom_mean_then_attn"])
    ap.add_argument("--label-modes", nargs="+", default=["ever_convert"],
                    choices=["ever_convert", "bl_multi"])
    ap.add_argument("--seeds", nargs="+", type=int, default=[0, 1, 2])
    ap.add_argument("--device", default="cpu")
    # MLP hyperparams = script 15 defaults
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--patience", type=int, default=20)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--weight-decay", type=float, default=1e-4)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--smoke", action="store_true",
                    help="tiny MLP budget (epochs=3, patience=2) for a fast end-to-end check")
    args = ap.parse_args()

    if args.smoke:
        args.epochs, args.patience = 3, 2

    model_keys = (list(MODEL_REGISTRY) if args.models == ["all"]
                  else [m for m in args.models if m in MODEL_REGISTRY])
    if not model_keys:
        raise SystemExit(f"No valid --models in {list(MODEL_REGISTRY)}")
    heads = ALL_HEADS if args.heads == ["all"] else args.heads

    info = preflight(model_keys)

    n_done = 0
    for mk in model_keys:
        if mk not in info:
            continue
        for label_mode in args.label_modes:
            if not SPLITS[label_mode].exists():
                print(f"[warn] splits for {label_mode} missing "
                      f"({SPLITS[label_mode]}) — skipping label_mode.")
                continue
            for agg in args.aggregations:
                allowed = HEAD_AGGS.get(agg)
                if allowed is None:
                    print(f"[warn] unknown aggregation {agg!r} — skipped.")
                    continue
                for head in heads:
                    if head not in allowed:
                        continue
                    for seed in args.seeds:
                        print(f"[{mk} | {label_mode} | {agg} | {head} | "
                              f"seed {seed}]")
                        run_cell(mk, info[mk], label_mode, agg, head,
                                 seed, args)
                        n_done += 1
    print(f"\nDone. {n_done} cells processed → {args.output_root}")


if __name__ == "__main__":
    main()
