"""
05_extract_embeddings_braindino_frozen.py
==========================================
Extract per-scan embeddings + per-class probabilities from the BrainDINO
frozen-cached classifier (HP-tuned winner per task) for use in downstream
multimodal models.

Why this script exists:
  The user's multimodal model needs per-scan (BrainDINO 768-d embedding,
  task-specific softmax probabilities) for tasks {T1, T1b, T1c, T1d, T2}.
  Because BrainDINO/frozen means the encoder is frozen at pretraining
  weights, the embedding is identical regardless of task/seed. The cached
  embeddings already on D: are reused; we only need to apply the trained
  task-specific HEAD to compute probabilities.

Inputs:
  D:/ADNI_BIDS_project/derivatives/cached_embeddings/braindino/
      braindino_pooled.pt   -- frozen-encoder embeddings (1447 scans x 768)
      manifest.csv          -- bids_sub + bids_ses per row

  D:/ADNI_BIDS_project/derivatives/braindino_outputs/
      aug_none_hp_tuned/BrainDINO_vitb16_frozen_cached/<task>/seed_<n>/
      <hp_winner>/best_model.pt   -- HeadOnly state_dict (LN -> L -> GELU -> Drop -> L)

  D:/ADNI_BIDS_project/derivatives/mri_clinical_matched/
      VISCODE_2_aligned_extended_post_exclusion/master_..._post_exclusion.csv
      -- the canonical MRI cohort (1456 scans) with conversion_group, labels.

  D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified_post_exclusion/
      tabular/baseline/seed_{0,1,2}/{train,val,test}.csv
      -- per-seed Patient_ID splits (for fold tagging).

HP-winner per task (from braindino_outputs/.../_winners/winners.csv):
  T1_binary    : lr1e-03_d0.2_ls0.1   (val_bacc 0.6677)
  T1b_binary   : lr1e-03_d0.2_ls0.1   (val_bacc 0.8123)
  T1c_binary   : lr1e-03_d0.2_ls0.0   (val_bacc 0.8442)
  T1d_binary   : lr1e-03_d0.2_ls0.1   (val_bacc 0.7827)
  T2_multiclass: lr1e-03_d0.2_ls0.1   (val_bacc 0.6019)

Output:
  D:/ADNI_BIDS_project/derivatives/mri_embeddings/<task>/braindino_frozen_none_cached/
      embeddings_seed_<n>.npz   -- {embeddings, logits, probs, ids, label,
                                    split_seed_<n>, model_meta}
      embeddings_seed_<n>.csv   -- flat row-per-scan format for inspection

Schema (NPZ):
  embeddings   : (N, 768) float32
  logits       : (N, K)   float32  -- pre-softmax
  probs        : (N, K)   float32  -- softmax (per-class probability)
  pred         : (N,)     int8     -- argmax
  bids_sub     : (N,)     str
  Patient_ID   : (N,)     str
  adni_viscode : (N,)     str
  label        : (N,)     int8     -- task label (-1 if not labelled in master)
  split        : (N,)     str      -- 'train'/'val'/'test'/'' (per the seed used)
  conversion_group : (N,) str
  model_meta   : dict     -- task, hp_winner, ckpt_path, n_classes, embed_dim

Leakage-safety note:
  The frozen encoder NEVER saw our labels (BrainDINO is pretrained
  externally), so embeddings are leakage-free.

  The HEAD however IS task-specific and trained on per-seed train folds.
  So per-scan probabilities are LEAKED for the seed where the scan was in
  the train fold. For downstream out-of-fold use the caller should pick,
  per scan, the seed where its `split == 'test'` (or 'val').

  This script writes one NPZ per seed with that seed's `split` column;
  downstream code can pick the right row per scan.

Usage:
  python mri_pipeline/05_extract_embeddings_braindino_frozen.py \
      [--task T1_binary,T1b_binary,T1c_binary,T1d_binary,T2_multiclass] \
      [--seeds 0,1,2] [--device cpu|cuda] [--out-dir <path>]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

# Work around MKL/libomp DLL clash before importing torch.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

# On Windows, torch must be imported BEFORE numpy/pandas for DLL ordering.
import torch  # noqa: E402
import torch.nn as nn  # noqa: E402
import torch.nn.functional as F  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# --------------------------------------------------------------------------- #
# Paths + constants
# --------------------------------------------------------------------------- #
CACHED_DIR = Path(r"D:/ADNI_BIDS_project/derivatives/cached_embeddings/braindino")
HEADS_ROOT = Path(r"D:/ADNI_BIDS_project/derivatives/braindino_outputs/"
                  r"aug_none_hp_tuned/BrainDINO_vitb16_frozen_cached")
MRI_MASTER = Path(r"D:/ADNI_BIDS_project/derivatives/mri_clinical_matched/"
                  r"VISCODE_2_aligned_extended_post_exclusion/"
                  r"master_mri_clinical_matched_viscode2_extended_post_exclusion.csv")
CLINICAL_SPLITS = Path(r"D:/ADNI_BIDS_project/derivatives/clinical/"
                       r"no_cdr_stratified_post_exclusion/tabular/baseline")
DEFAULT_OUT = Path(r"D:/ADNI_BIDS_project/derivatives/mri_embeddings")

# Per-task HP winners (from winners.csv).
HP_WINNERS = {
    "T1_binary":     "lr1e-03_d0.2_ls0.1",
    "T1b_binary":    "lr1e-03_d0.2_ls0.1",
    "T1c_binary":    "lr1e-03_d0.2_ls0.0",
    "T1d_binary":    "lr1e-03_d0.2_ls0.1",
    "T2_multiclass": "lr1e-03_d0.2_ls0.1",
}

# Per-task label-builder logic from the project's MRI master columns.
# Aligns with TASK_CONFIG in 04_supervised_finetuning.py.
def task_label(task: str, row: pd.Series) -> int | None:
    """Compute the task-specific class label for one MRI master row.
    Returns None if the row is not eligible for this task."""
    bl = row.get("bl_dx")
    last = row.get("last_dx")
    cg = row.get("conversion_group")
    visit_diag = row.get("Label_visit_diag")

    if task == "T1_binary":      # CN vs (MCI + AD) at visit
        if visit_diag in ("CN",):
            return 0
        if visit_diag in ("MCI", "AD"):
            return 1
        return None
    if task == "T1b_binary":     # CN vs AD at visit
        if visit_diag == "CN":
            return 0
        if visit_diag == "AD":
            return 1
        return None
    if task == "T1c_binary":     # MCI vs AD at visit
        if visit_diag == "MCI":
            return 0
        if visit_diag == "AD":
            return 1
        return None
    if task == "T1d_binary":     # pMCI vs sMCI (visit-time MCI subjects, by conversion_group)
        if visit_diag != "MCI":
            return None
        if cg == "sMCI":
            return 0
        if cg == "pMCI":
            return 1
        return None
    if task == "T2_multiclass":  # CN/MCI/AD at visit
        return {"CN": 0, "MCI": 1, "AD": 2}.get(visit_diag, None)
    return None


# --------------------------------------------------------------------------- #
# HeadOnly (mirrors cached_head_sweep/04_head_finetune_from_embeddings.py)
# --------------------------------------------------------------------------- #
class HeadOnly(nn.Module):
    """LayerNorm -> Linear -> GELU -> Dropout -> Linear."""
    def __init__(self, embed_dim: int, n_classes: int, drop_rate: float = 0.2):
        super().__init__()
        self.net = nn.Sequential(
            nn.LayerNorm(embed_dim),
            nn.Linear(embed_dim, embed_dim),
            nn.GELU(),
            nn.Dropout(drop_rate),
            nn.Linear(embed_dim, n_classes),
        )

    def forward(self, x):
        return self.net(x)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def n_classes_for_task(task: str) -> int:
    if task == "T2_multiclass":
        return 3
    return 2


def ptid_to_bids_sub(ptid: str) -> str:
    return "sub-" + ptid.replace("_", "")


def ptid_to_rid(ptid: str) -> int | None:
    """'041_S_4629' -> 4629. None if unparseable."""
    if not isinstance(ptid, str):
        return None
    parts = ptid.split("_")
    try:
        return int(parts[-1])
    except (ValueError, IndexError):
        return None


def load_cached_embeddings() -> tuple[torch.Tensor, list[tuple[str, str]]]:
    """Returns (embeddings (N, D), ids list of (bids_sub, bids_ses))."""
    obj = torch.load(CACHED_DIR / "braindino_pooled.pt",
                     map_location="cpu", weights_only=False)
    emb = obj["embeddings"]               # torch float32 (N, 768)
    ids = obj["ids"]                       # list of (bids_sub, bids_ses)
    return emb, ids


def load_master_cohort() -> pd.DataFrame:
    df = pd.read_csv(MRI_MASTER, low_memory=False)
    return df


def load_split_for_seed(seed: int, task: str) -> dict[str, str]:
    """Map Patient_ID -> 'train'/'val'/'test' for the given seed.
    Falls back to '' if Patient_ID not in any split (e.g. T1d eligibility-filtered).
    """
    seed_dir = CLINICAL_SPLITS / f"seed_{seed}"
    out: dict[str, str] = {}
    for fold in ("train", "val", "test"):
        f = seed_dir / f"{fold}.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f, low_memory=False)
        if "Patient_ID" not in df.columns:
            continue
        for p in df["Patient_ID"].dropna().unique():
            out[p] = fold
    return out


def load_head_state(task: str, seed: int) -> dict:
    hp = HP_WINNERS[task]
    p = HEADS_ROOT / task / f"seed_{seed}" / hp / "best_model.pt"
    if not p.exists():
        raise FileNotFoundError(f"Head checkpoint missing: {p}")
    sd = torch.load(p, map_location="cpu", weights_only=False)
    state = sd["net"] if isinstance(sd, dict) and "net" in sd else sd
    return state, str(p)


# --------------------------------------------------------------------------- #
# Per-task extraction
# --------------------------------------------------------------------------- #
def extract_one_task_seed(task: str, seed: int, emb: torch.Tensor,
                          ids: list[tuple[str, str]], master: pd.DataFrame,
                          device: str) -> dict:
    """Builds the per-scan output dict for one (task, seed)."""
    K = n_classes_for_task(task)
    head = HeadOnly(embed_dim=emb.shape[1], n_classes=K, drop_rate=0.2)
    state, ckpt_path = load_head_state(task, seed)
    missing, unexpected = head.load_state_dict(state, strict=False)
    if missing or unexpected:
        print(f"  [WARN] task={task} seed={seed}: missing={missing}, "
              f"unexpected={unexpected}", file=sys.stderr)
    head.to(device).eval()

    # Build index from (bids_sub, bids_ses) → row in cache.
    cache_idx = {(s, v): i for i, (s, v) in enumerate(ids)}

    # Iterate master scans (1456 rows). Map each scan to its cached embedding.
    embeds, hits, miss = [], [], []
    row_meta = []
    for _, r in master.iterrows():
        s = r["bids_sub"] if r["bids_sub"].startswith("sub-") else "sub-" + r["bids_sub"]
        v = r["adni_viscode"]
        key = (s, v)
        if key in cache_idx:
            embeds.append(emb[cache_idx[key]])
            hits.append(True)
        else:
            embeds.append(torch.full((emb.shape[1],), float("nan")))
            miss.append(key)
            hits.append(False)
        row_meta.append({
            "bids_sub":         r["bids_sub"],
            "Patient_ID":       r["Patient_ID"],
            "RID":              ptid_to_rid(r["Patient_ID"]),
            "adni_viscode":     r["adni_viscode"],
            "VISCODE_2":        r.get("VISCODE_2"),
            "mri_date":         r.get("mri_date"),
            "conversion_group": r.get("conversion_group"),
            "label_raw":        task_label(task, r),
        })

    X = torch.stack(embeds).to(device)                  # (N, D); rows for miss are NaN
    with torch.no_grad():
        # Mask NaN rows -- inference would produce NaN logits.
        nan_mask = torch.isnan(X).any(dim=1)
        X_safe = torch.where(nan_mask.unsqueeze(1), torch.zeros_like(X), X)
        logits = head(X_safe)                            # (N, K)
        probs = F.softmax(logits, dim=1)                 # (N, K)
        preds = probs.argmax(dim=1)                      # (N,)
        # Re-NaN out the bad rows.
        logits[nan_mask] = float("nan")
        probs[nan_mask] = float("nan")
        preds[nan_mask] = -1

    # Per-seed split mapping.
    ptid_to_split = load_split_for_seed(seed, task)
    splits = [ptid_to_split.get(rm["Patient_ID"], "") for rm in row_meta]
    labels = [rm["label_raw"] if rm["label_raw"] is not None else -1 for rm in row_meta]

    print(f"  task={task} seed={seed}: cache_hit={sum(hits)}/{len(hits)}, "
          f"miss={len(miss)}; split_counts={pd.Series(splits).value_counts().to_dict()}")

    return {
        "embeddings":     X.cpu().numpy().astype("float32"),
        "logits":         logits.cpu().numpy().astype("float32"),
        "probs":          probs.cpu().numpy().astype("float32"),
        "pred":           preds.cpu().numpy().astype("int8"),
        "bids_sub":       np.array([rm["bids_sub"] for rm in row_meta], dtype=object),
        "Patient_ID":     np.array([rm["Patient_ID"] for rm in row_meta], dtype=object),
        "RID":            np.array([rm["RID"] if rm["RID"] is not None else -1 for rm in row_meta], dtype="int32"),
        "adni_viscode":   np.array([rm["adni_viscode"] for rm in row_meta], dtype=object),
        "VISCODE_2":      np.array([rm["VISCODE_2"] for rm in row_meta], dtype=object),
        "mri_date":       np.array([str(rm["mri_date"]) for rm in row_meta], dtype=object),
        "label":          np.array(labels, dtype="int8"),
        "split":          np.array(splits, dtype=object),
        "conversion_group": np.array([rm["conversion_group"] for rm in row_meta], dtype=object),
        "model_meta": {
            "task":           task,
            "seed":           seed,
            "model":          "BrainDINO_vitb16_frozen_cached",
            "encoder_source": "DINO-pretrained ViT-B/16 (frozen)",
            "embed_dim":      int(X.shape[1]),
            "n_classes":      K,
            "hp_winner":      HP_WINNERS[task],
            "head_ckpt":      ckpt_path,
            "cached_embeddings_pt": str(CACHED_DIR / "braindino_pooled.pt"),
            "cache_hits":     int(sum(hits)),
            "cache_misses":   int(len(miss)),
            "miss_examples":  [list(k) for k in miss[:10]],
        },
    }


# --------------------------------------------------------------------------- #
# Output writers + verification
# --------------------------------------------------------------------------- #
def write_outputs(out_dir: Path, task: str, seed: int, payload: dict) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    npz_path = out_dir / f"embeddings_seed_{seed}.npz"
    np.savez(
        npz_path,
        embeddings=payload["embeddings"],
        logits=payload["logits"],
        probs=payload["probs"],
        pred=payload["pred"],
        bids_sub=payload["bids_sub"],
        Patient_ID=payload["Patient_ID"],
        RID=payload["RID"],
        adni_viscode=payload["adni_viscode"],
        VISCODE_2=payload["VISCODE_2"],
        mri_date=payload["mri_date"],
        label=payload["label"],
        split=payload["split"],
        conversion_group=payload["conversion_group"],
        model_meta=np.array(json.dumps(payload["model_meta"])),
    )
    # Flat CSV (sans embedding cols -- keeps file small + readable).
    K = payload["probs"].shape[1]
    flat = pd.DataFrame({
        "bids_sub":        payload["bids_sub"],
        "Patient_ID":      payload["Patient_ID"],
        "RID":             payload["RID"],
        "adni_viscode":    payload["adni_viscode"],
        "VISCODE_2":       payload["VISCODE_2"],
        "mri_date":        payload["mri_date"],
        "label":           payload["label"],
        "split":           payload["split"],
        "conversion_group": payload["conversion_group"],
        "pred":            payload["pred"],
    })
    for k in range(K):
        flat[f"prob_class_{k}"] = payload["probs"][:, k]
        flat[f"logit_class_{k}"] = payload["logits"][:, k]
    flat.to_csv(out_dir / f"embeddings_seed_{seed}.csv", index=False)
    print(f"    wrote {npz_path} ({npz_path.stat().st_size/1e6:.2f} MB) + .csv")


def verify_one(payload: dict) -> None:
    """Run a basic verification battery on the output."""
    emb = payload["embeddings"]
    probs = payload["probs"]
    n = emb.shape[0]
    n_nan = int(np.isnan(emb).any(axis=1).sum())
    n_ok = n - n_nan
    print(f"    verify: n={n}, ok_rows={n_ok}, nan_rows={n_nan}; "
          f"emb_dim={emb.shape[1]}; probs_sum_mean={np.nanmean(probs.sum(axis=1)):.6f} "
          f"(should be ~1.0)")
    # Class balance check.
    valid = payload["label"] >= 0
    if valid.sum() > 0:
        labs = payload["label"][valid]
        preds = payload["pred"][valid]
        agreement = (labs == preds).mean()
        print(f"    label_distribution_valid: {pd.Series(labs).value_counts().to_dict()}")
        print(f"    pred==label fraction (sanity, may be biased by overfit head): {agreement:.3f}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tasks", default=",".join(HP_WINNERS.keys()),
                   help="Comma-separated task list.")
    p.add_argument("--seeds", default="0,1,2")
    p.add_argument("--device", default="cpu", choices=["cpu", "cuda"])
    p.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    tasks = [t.strip() for t in args.tasks.split(",") if t.strip()]
    seeds = [int(s) for s in args.seeds.split(",")]

    print("[05] BrainDINO frozen-cached embedding extraction (local)")
    print(f"  tasks: {tasks}")
    print(f"  seeds: {seeds}")
    print(f"  device: {args.device}")
    print(f"  out:    {args.out_dir}")
    print()

    print("Loading cached BrainDINO frozen embeddings ...")
    emb, ids = load_cached_embeddings()
    print(f"  shape: {tuple(emb.shape)}; n_ids: {len(ids)}")

    print("Loading MRI master cohort ...")
    master = load_master_cohort()
    print(f"  shape: {master.shape}; unique Patient_ID: {master.Patient_ID.nunique()}")

    print()
    for task in tasks:
        out_task = args.out_dir / task / "braindino_frozen_none_cached"
        for seed in seeds:
            print(f"-- task={task} seed={seed} --")
            payload = extract_one_task_seed(task, seed, emb, ids, master, args.device)
            write_outputs(out_task, task, seed, payload)
            verify_one(payload)
            print()

    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
