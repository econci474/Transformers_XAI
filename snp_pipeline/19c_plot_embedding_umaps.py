"""
19c_plot_embedding_umaps.py
============================
UMAP of each foundation model's frozen per-patient embedding, coloured by
3-class baseline diagnosis (CN / MCI / AD) and by binary ``ever_convert``.

Per-patient vector = mean per-window pool → ``chrom_mean_then_mean``
(parameter-free; the SAME feature family the classifier bake-off in
``19_compare_frozen_models.py`` uses, so the UMAP visualises exactly what the
heads see). Model registry, paths and the script-15 helpers are reused from
script 19 (imported by path — digit-leading filename). UMAP params mirror
``snp_pipeline/02_pca.py``; styling matches the 16b/19b house style
(DejaVu Serif, white, no bold).

Outputs (under --out-root, default frozen_model_comparison/umaps):
  <model>_dx.png            scatter coloured CN/MCI/AD
  <model>_everconvert.png   scatter coloured stable/convert
  all_models_dx.png         combined grid (one panel per model)

Usage
-----
  python snp_pipeline/19c_plot_embedding_umaps.py [--models all] \\
      [--out-root frozen_model_comparison/umaps]
"""
from __future__ import annotations

import argparse
import importlib.util as _il
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

plt.rcParams.update({
    "font.family": "DejaVu Serif",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": "#cccccc",
})

# ── Reuse script 19 (registry + script-15 helpers); digit-leading filename ──
_S19 = Path(__file__).parent / "19_compare_frozen_models.py"
_spec = _il.spec_from_file_location("script19", _S19)
_m19 = _il.module_from_spec(_spec)
_spec.loader.exec_module(_m19)
MODEL_REGISTRY       = _m19.MODEL_REGISTRY
MODEL_ORDER          = list(MODEL_REGISTRY)
SPLITS               = _m19.SPLITS
build_chrom_map      = _m19.build_chrom_map
load_embeddings      = _m19.load_embeddings
stack_patient_arrays = _m19.stack_patient_arrays
Aggregator           = _m19.Aggregator

MODEL_PRETTY = {
    "bmfm_base_ref": "BMFM-DNA-REF (frozen)",
    "ntv2": "Nucleotide Transformer v2",
    "caduceus_ps_8k": "Caduceus-PS (8 kb)",
    "caduceus_ph_8k": "Caduceus-PH (8 kb)",
    "caduceus_ps_131k": "Caduceus-PS (131 kb)",
    "caduceus_ph_131k": "Caduceus-PH (131 kb)",
}
DX_COLORS  = {"CN": "#4fc3f7", "MCI": "#ffa726", "AD": "#ef5350", "Other": "#999999"}
DX_ORDER   = ["CN", "MCI", "AD", "Other"]
EVC_COLORS = {"stable": "#4fc3f7", "convert": "#ef5350"}
EVC_ORDER  = ["stable", "convert"]


def patient_label_table() -> pd.DataFrame:
    """Patient_ID → (Label_bl_multi, ever_convert) from the ever_convert
    seed_0 split CSVs (union of train/val/test)."""
    root = SPLITS["ever_convert"] / "seed_0"
    frames = []
    for split in ("train", "val", "test"):
        df = pd.read_csv(root / f"{split}.csv",
                         usecols=["Patient_ID", "Label_bl_multi",
                                  "Label_ever_convert"])
        frames.append(df)
    lab = pd.concat(frames, ignore_index=True).drop_duplicates("Patient_ID")
    lab["Patient_ID"] = lab["Patient_ID"].astype(str)
    lab["dx"] = lab["Label_bl_multi"].where(
        lab["Label_bl_multi"].isin(["CN", "MCI", "AD"]), "Other")
    lab["evc"] = lab["Label_ever_convert"].map({0: "stable", 1: "convert"})
    return lab.set_index("Patient_ID")


def model_umap(mk: str, lab: pd.DataFrame):
    """Return (df with UMAP1/2/dx/evc) for one model, or None if unavailable."""
    import umap  # local import: clean error if umap-learn missing
    reg = MODEL_REGISTRY[mk]
    npz = Path(reg["npz"])
    if not npz.exists():
        print(f"  [skip] {mk}: npz not found ({npz}).")
        return None
    if not Path(reg["windows"]).exists():
        print(f"  [skip] {mk}: windows.tsv missing ({reg['windows']}).")
        return None
    emb, win_ids = load_embeddings(npz, pool="mean")
    if not emb:
        print(f"  [skip] {mk}: no 'mean__' embeddings.")
        return None
    chrom_of_window, chroms = build_chrom_map(Path(reg["windows"]))
    pids = [p for p in map(str, emb.keys()) if p in lab.index]
    if len(pids) < 10:
        print(f"  [skip] {mk}: only {len(pids)} labelled patients.")
        return None
    x, m, c = stack_patient_arrays(emb, win_ids, chrom_of_window, pids)
    hidden = next(iter(emb.values())).shape[1]
    with torch.no_grad():
        vecs = Aggregator("chrom_mean_then_mean", hidden,
                          len(chroms)).eval()(x, m, c).cpu().numpy()
    reducer = umap.UMAP(n_components=2, n_neighbors=15, min_dist=0.3,
                        metric="euclidean", random_state=42)
    xy = reducer.fit_transform(vecs)
    sub = lab.loc[pids]
    return pd.DataFrame({"UMAP1": xy[:, 0], "UMAP2": xy[:, 1],
                         "dx": sub["dx"].to_numpy(),
                         "evc": sub["evc"].to_numpy()})


def _scatter(ax, df, key, colors, order, title):
    for lvl in order:
        mask = df[key] == lvl
        if not mask.any():
            continue
        ax.scatter(df.loc[mask, "UMAP1"], df.loc[mask, "UMAP2"],
                   c=colors[lvl], s=16, alpha=0.8, linewidths=0,
                   label=f"{lvl} (n={int(mask.sum())})")
    ax.set_xlabel("UMAP 1", fontsize=10)
    ax.set_ylabel("UMAP 2", fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(fontsize=8, framealpha=0.85, markerscale=1.3)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", nargs="+", default=["all"])
    ap.add_argument("--out-root", type=Path,
                    default=Path("frozen_model_comparison/umaps"))
    args = ap.parse_args()
    args.out_root.mkdir(parents=True, exist_ok=True)

    keys = (MODEL_ORDER if args.models == ["all"]
            else [k for k in args.models if k in MODEL_REGISTRY])
    lab = patient_label_table()
    print(f"Labels: {len(lab)} patients "
          f"({lab['dx'].value_counts().to_dict()})")

    results: dict[str, pd.DataFrame] = {}
    for mk in keys:
        print(f"[umap] {mk}")
        df = model_umap(mk, lab)
        if df is None:
            continue
        results[mk] = df
        for key, colors, order, tag, ttl in (
            ("dx",  DX_COLORS,  DX_ORDER,  "dx",
             f"{MODEL_PRETTY.get(mk, mk)} — baseline diagnosis"),
            ("evc", EVC_COLORS, EVC_ORDER, "everconvert",
             f"{MODEL_PRETTY.get(mk, mk)} — ever_convert"),
        ):
            fig, ax = plt.subplots(figsize=(7, 6), dpi=170)
            _scatter(ax, df, key, colors, order, ttl)
            fig.tight_layout()
            p = args.out_root / f"{mk}_{tag}.png"
            fig.savefig(p, dpi=170, bbox_inches="tight", facecolor="white")
            plt.close(fig)
            print(f"  saved {p}")

    if results:
        n = len(results)
        ncol = min(3, n)
        nrow = math.ceil(n / ncol)
        fig, axes = plt.subplots(nrow, ncol, figsize=(6.2 * ncol, 5.4 * nrow),
                                 dpi=160, squeeze=False)
        for ax in axes.flat:
            ax.set_axis_off()
        for i, (mk, df) in enumerate(results.items()):
            ax = axes[i // ncol][i % ncol]
            ax.set_axis_on()
            _scatter(ax, df, "dx", DX_COLORS, DX_ORDER,
                     MODEL_PRETTY.get(mk, mk))
        fig.suptitle("Frozen embedding UMAPs — baseline diagnosis (CN/MCI/AD)",
                     fontsize=14, y=1.005)
        fig.tight_layout()
        p = args.out_root / "all_models_dx.png"
        fig.savefig(p, dpi=160, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"  saved {p}")
    else:
        print("No models produced UMAPs (stage the npz files and re-run).")


if __name__ == "__main__":
    main()
