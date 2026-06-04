"""
05_extract_cached_head_probs.py
===============================
Generalised per-scan probability extractor for the frozen cached-head models
(BrainDINO / BrainMVP / ViT-MAE). Companion to
`05_extract_embeddings_braindino_frozen.py` (BrainDINO-only); this one takes
`--arch` so the same code produces fusion-ready per-patient probability CSVs
for any of the three cached encoders, on any task — built to supply the
T3/T4 late-fusion arm with BrainMVP/ViT-MAE cached probs (no extractor for
those existed).

For each (arch, task, seed) it:
  1. loads the frozen cached embeddings `cached_embeddings/<arch>/<arch>_pooled.pt`
     ({embeddings (N,D), ids list[(bids_sub,bids_ses)]});
  2. selects the HP-winner per (arch,task) = highest mean
     `config.best_val_balanced_acc` across seeds (scanning the cached sweep
     metrics.json — same rule the cross-model table uses), loads that
     `best_model.pt["net"]` into HeadOnly;
  3. forwards every master scan, writes softmax probs for ALL splits/sessions.

Output (the exact schema the fusion engine `integration/T2_classification/
01_fuse_T2.py::load_mri` expects):
  D:/.../mri_embeddings/<task>/<arch>_frozen_none_cached/embeddings_seed_<n>.csv
  cols: Patient_ID, adni_viscode, split, label, prob_class_0, prob_class_1[, prob_class_2], ...

Usage:
  python mri_pipeline/05_extract_cached_head_probs.py --arch vit_mae \
      --tasks T3a_conv3y,T3b_conv5y,T3c_conv7y
  python mri_pipeline/05_extract_cached_head_probs.py --arch brainmvp \
      --tasks T3a_conv3y,T3b_conv5y,T3c_conv7y
  # T4 (3-class) once its cached sweeps are local:
  #   --arch vit_mae --tasks T4_conv_horizon --splits-dir baseline_T4
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import torch  # noqa: E402
import torch.nn as nn  # noqa: E402
import torch.nn.functional as F  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# --------------------------------------------------------------------------- #
DERIV = Path(r"D:/ADNI_BIDS_project/derivatives")
MRI_MASTER = (DERIV / "mri_clinical_matched" / "VISCODE_2_aligned_extended_post_exclusion"
              / "master_mri_clinical_matched_viscode2_extended_post_exclusion.csv")
CLINICAL_SPLITS_ROOT = DERIV / "clinical" / "no_cdr_stratified_post_exclusion" / "tabular"
DEFAULT_OUT = DERIV / "mri_embeddings"

# Per-arch cached-embedding file, cached-sweep head root, and embed dim.
ARCHS = {
    "braindino": dict(
        pooled=DERIV / "cached_embeddings/braindino/braindino_pooled.pt",
        heads=DERIV / "braindino_outputs/aug_none_hp_tuned/BrainDINO_vitb16_frozen_cached",
        embed_dim=768, model="BrainDINO_vitb16_frozen_cached"),
    "brainmvp": dict(
        pooled=DERIV / "cached_embeddings/brainmvp/brainmvp_pooled.pt",
        heads=DERIV / "brainmvp_debug/aug_none/BrainMVP_uniformer_frozen_cached",
        embed_dim=512, model="BrainMVP_uniformer_frozen_cached"),
    "vit_mae": dict(
        pooled=DERIV / "cached_embeddings/vit_mae/vit_mae_pooled.pt",
        heads=DERIV / "vit_outputs_debug/aug_none/ViT_B_mae75_frozen_cached",
        embed_dim=768, model="ViT_B_mae75_frozen_cached"),
}

# task -> (n_classes, label resolver). Conversion tasks read the subject-level
# Label_Ny / Label_T4 columns broadcast onto every visit row in the master.
THREE_CLASS = {"T2_multiclass", "T4_conv_horizon"}
T3_LABEL_COL = {"T3a_conv3y": "Label_3y", "T3b_conv5y": "Label_5y", "T3c_conv7y": "Label_7y"}


def n_classes_for_task(task: str) -> int:
    return 3 if task in THREE_CLASS else 2


def task_label(task: str, row: pd.Series):
    """Task-specific class label for one master row, or None if ineligible.
    Mirrors TASK_CONFIG in 04_supervised_finetuning_ViT.py."""
    visit_diag = row.get("Label_visit_diag")
    cg = row.get("conversion_group")
    # --- diagnostic tasks (kept for parity with the braindino extractor) ---
    if task == "T1_binary":
        return 0 if visit_diag == "CN" else (1 if visit_diag in ("MCI", "AD") else None)
    if task == "T1b_binary":
        return 0 if visit_diag == "CN" else (1 if visit_diag == "AD" else None)
    if task == "T1c_binary":
        return 0 if visit_diag == "MCI" else (1 if visit_diag == "AD" else None)
    if task == "T1d_binary":
        if visit_diag != "MCI":
            return None
        return 0 if cg == "sMCI" else (1 if cg == "pMCI" else None)
    if task == "T2_multiclass":
        return {"CN": 0, "MCI": 1, "AD": 2}.get(visit_diag, None)
    # --- conversion tasks (T3 binary, baseline-anchored, filter_non_ad) ---
    if task in T3_LABEL_COL:
        if row.get("Label_bl_multi") == "AD":          # filter_non_ad=True
            return None
        v = row.get(T3_LABEL_COL[task])
        return int(v) if pd.notna(v) else None
    # --- T4 3-class horizon (converters only; no AD filter) ---
    if task == "T4_conv_horizon":
        v = row.get("Label_T4")
        return int(v) if pd.notna(v) else None
    return None


class StandardizeLayer(nn.Module):
    """Frozen z-score baked as buffers (mirrors the trainer's StandardizeLayer)."""
    def __init__(self, dim: int):
        super().__init__()
        self.register_buffer("mean", torch.zeros(dim))
        self.register_buffer("std", torch.ones(dim))

    def forward(self, x):
        return (x - self.mean) / self.std


class HeadOnly(nn.Module):
    """Must match cached_head_sweep/04_head_finetune_from_embeddings.py::HeadOnly
    exactly (same layer order/keys) so the winner checkpoint loads. Defaults
    (head='mlp', standardize=False) reproduce the original architecture."""
    def __init__(self, embed_dim: int, n_classes: int, drop_rate: float = 0.2,
                 head: str = "mlp", standardize: bool = False):
        super().__init__()
        layers = []
        if standardize:
            layers.append(StandardizeLayer(embed_dim))
        if head == "linear":
            layers += [nn.LayerNorm(embed_dim), nn.Dropout(drop_rate),
                       nn.Linear(embed_dim, n_classes)]
        else:  # mlp
            layers += [nn.LayerNorm(embed_dim), nn.Linear(embed_dim, embed_dim),
                       nn.GELU(), nn.Dropout(drop_rate), nn.Linear(embed_dim, n_classes)]
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


def ptid_to_rid(ptid):
    if not isinstance(ptid, str):
        return None
    try:
        return int(ptid.split("_")[-1])
    except (ValueError, IndexError):
        return None


def select_hp_winner(heads_root: Path, task: str) -> str:
    """HP-winner dir name = highest mean config.best_val_balanced_acc across
    seeds (same rule as the cross-model table's _filter_cached_to_hp_winners)."""
    rows = []
    for mp in glob.glob(str(heads_root / task / "seed_*" / "*" / "metrics.json")):
        try:
            cfg = json.load(open(mp)).get("config", {})
            rows.append((os.path.basename(os.path.dirname(mp)),
                         cfg.get("best_val_balanced_acc")))
        except Exception:
            continue
    if not rows:
        raise FileNotFoundError(f"No metrics.json under {heads_root/task}")
    df = pd.DataFrame(rows, columns=["hp", "val"]).dropna()
    mean = df.groupby("hp")["val"].mean()
    winner = mean.idxmax()
    print(f"    HP-winner[{task}] = {winner}  (mean val_bacc {mean.max():.4f}, "
          f"{df['hp'].nunique()} HP combos)")
    return winner


def load_split_for_seed(splits_dir: Path, seed: int) -> dict:
    out = {}
    for fold in ("train", "val", "test"):
        f = splits_dir / f"seed_{seed}" / f"{fold}.csv"
        if not f.exists():
            continue
        d = pd.read_csv(f, low_memory=False)
        if "Patient_ID" in d.columns:
            for p in d["Patient_ID"].dropna().unique():
                out[p] = fold
    return out


def extract(arch, task, seed, emb, ids, master, splits_dir, device, heads_root=None):
    acfg = ARCHS[arch]
    K = n_classes_for_task(task)
    hroot = Path(heads_root) if heads_root else acfg["heads"]
    hp = select_hp_winner(hroot, task)
    hp_dir = hroot / task / f"seed_{seed}" / hp
    ckpt = hp_dir / "best_model.pt"
    if not ckpt.exists():
        raise FileNotFoundError(
            f"Head checkpoint missing: {ckpt}\n"
            f"        (a --metrics_only sweep tree has no best_model.pt — point "
            f"--heads-root at the finalised winner tree that re-trained the winner.)")
    # Rebuild the EXACT head the winner used (head type / standardize / drop) from
    # its metrics.json config; defaults reproduce the original mlp head for older runs.
    cfg = json.load(open(hp_dir / "metrics.json")).get("config", {}) if (hp_dir / "metrics.json").exists() else {}
    head_type = cfg.get("head", "mlp")
    standardize = bool(cfg.get("standardize", False))
    drop_rate = float(cfg.get("drop_rate", 0.2))
    sd = torch.load(ckpt, map_location="cpu", weights_only=False)
    state = sd["net"] if isinstance(sd, dict) and "net" in sd else sd
    head = HeadOnly(acfg["embed_dim"], K, drop_rate=drop_rate,
                    head=head_type, standardize=standardize)
    missing, unexpected = head.load_state_dict(state, strict=False)
    if missing or unexpected:
        print(f"    [WARN] load_state_dict missing={missing} unexpected={unexpected}",
              file=sys.stderr)
    print(f"    head={head_type}  standardize={standardize}  drop={drop_rate}")
    head.to(device).eval()

    cache_idx = {(s, v): i for i, (s, v) in enumerate(ids)}
    embeds, hits, meta = [], 0, []
    for _, r in master.iterrows():
        s = r["bids_sub"] if str(r["bids_sub"]).startswith("sub-") else "sub-" + str(r["bids_sub"])
        key = (s, r["adni_viscode"])
        if key in cache_idx:
            embeds.append(emb[cache_idx[key]]); hits += 1
        else:
            embeds.append(torch.full((acfg["embed_dim"],), float("nan")))
        meta.append(r)

    X = torch.stack(embeds).to(device)
    with torch.no_grad():
        nan = torch.isnan(X).any(dim=1)
        logits = head(torch.where(nan.unsqueeze(1), torch.zeros_like(X), X))
        probs = F.softmax(logits, dim=1)
        preds = probs.argmax(dim=1)
        logits[nan] = float("nan"); probs[nan] = float("nan"); preds[nan] = -1

    p2s = load_split_for_seed(splits_dir, seed)
    flat = pd.DataFrame({
        "bids_sub":        [r["bids_sub"] for r in meta],
        "Patient_ID":      [r["Patient_ID"] for r in meta],
        "RID":             [ptid_to_rid(r["Patient_ID"]) for r in meta],
        "adni_viscode":    [r["adni_viscode"] for r in meta],
        "VISCODE_2":       [r.get("VISCODE_2") for r in meta],
        "mri_date":        [str(r.get("mri_date")) for r in meta],
        "conversion_group": [r.get("conversion_group") for r in meta],
        "label":           [(lambda v: v if v is not None else -1)(task_label(task, r)) for r in meta],
        "split":           [p2s.get(r["Patient_ID"], "") for r in meta],
        "pred":            preds.cpu().numpy().astype("int8"),
    })
    P = probs.cpu().numpy()
    for k in range(K):
        flat[f"prob_class_{k}"] = P[:, k]
    print(f"    cache_hit={hits}/{len(meta)}  split={pd.Series(flat['split']).value_counts().to_dict()}  "
          f"probs_sum~{np.nanmean(P.sum(1)):.4f}")
    return flat


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arch", required=True, choices=list(ARCHS))
    ap.add_argument("--tasks", required=True, help="Comma-separated task list.")
    ap.add_argument("--seeds", default="0,1,2")
    ap.add_argument("--splits-dir", default="baseline",
                    help="Sub-dir under .../no_cdr_stratified_post_exclusion/tabular/ "
                         "(baseline for T3; baseline_T4 for T4).")
    ap.add_argument("--heads-root", default=None,
                    help="Override the cached-head tree (e.g. the finalised T2 wide-sweep "
                         "winner tree). Defaults to ARCHS[arch]['heads'].")
    ap.add_argument("--device", default="cpu", choices=["cpu", "cuda"])
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    tasks = [t.strip() for t in args.tasks.split(",") if t.strip()]
    seeds = [int(s) for s in args.seeds.split(",")]
    splits_dir = CLINICAL_SPLITS_ROOT / args.splits_dir

    print(f"[05] cached-head prob extraction  arch={args.arch}  tasks={tasks}  seeds={seeds}")
    obj = torch.load(ARCHS[args.arch]["pooled"], map_location="cpu", weights_only=False)
    emb, ids = obj["embeddings"], obj["ids"]
    print(f"  cached embeddings: {tuple(emb.shape)}  ({len(ids)} scans)")
    master = pd.read_csv(MRI_MASTER, low_memory=False)
    print(f"  master: {master.shape}")

    for task in tasks:
        out_task = args.out_dir / task / f"{args.arch}_frozen_none_cached"
        out_task.mkdir(parents=True, exist_ok=True)
        for seed in seeds:
            print(f"-- {args.arch} {task} seed={seed} --")
            flat = extract(args.arch, task, seed, emb, ids, master, splits_dir,
                           args.device, heads_root=args.heads_root)
            flat.to_csv(out_task / f"embeddings_seed_{seed}.csv", index=False)
            print(f"    wrote {out_task / f'embeddings_seed_{seed}.csv'}")
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
