r"""
28_train_diff_attention_classifier.py   (LOCAL, env: snp, CPU)
=============================================================
One grid cell: train a β-biased-attention + MLP classifier for
**ad_case vs stable_CN** on the per-SNP signed-effect embeddings
d = (ALT_emb − REF_emb) · beta_A1 · dosage, for one
(snp_set, model, aggregation, mlp_depth, delta, seed).

Everything model/data is in fm_diff_lib (the shared contract with 26/27).
This script only orchestrates: assemble splits → train → evaluate → write
the full per-run output tree.

Inputs:
  GWAS_filtered/<SET>/<SET>_snps.tsv            (beta_A1, CHR, alleles)
  GWAS_filtered/<SET>/fm_sequences/<SET>_snp_manifest.tsv  (QC drop authority)
  GWAS_filtered/<SET>/<SET>_patient_dosage.tsv  (ALT dosage; cohort-mean impute)
  <diffs-root>/<model>/<SET>_snp_embeddings.npz (script 27: ref_emb, alt_emb)
  conversion_labels.tsv (ad_case=1 / stable_CN=0)
  no_cdr_stratified/.../seed_<seed>/{train,val,test}.csv  (membership)

Output tree:
  <out>/<SET>/<model>/<agg>/<depth>/delta_<off|on>/seed_<n>/
    config.json metrics.json predictions.csv model.pt val_curve.csv
    attention.npz embeddings.npz umaps/{patient_emb,dbg_diff,dbg_altref}.png

Usage (smoke):
  conda run -n snp python snp_pipeline/28_train_diff_attention_classifier.py \
    --snp-set 9W27B --model bmfm_ref --aggregation global_attn \
    --mlp-depth 2L --delta off --seed 0 --output-root \
    snp_pipeline/outputs/diff_attn --smoke
"""
from __future__ import annotations

import argparse
import importlib.util as _il
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

_LIB = Path(__file__).parent / "fm_diff_lib.py"
_spec = _il.spec_from_file_location("fm_diff_lib", _LIB)
fl = _il.module_from_spec(_spec)
_spec.loader.exec_module(fl)

WANDB_KEY_FILE = Path(__file__).parent.parent / "mri_pipeline" / "wandb_api.txt"


def init_wandb(args, config: dict):
    """ImportError-safe (pattern: 18 L433-450). Reads the gitignored
    mri_pipeline/wandb_api.txt → WANDB_API_KEY (never echoed)."""
    if not args.wandb:
        return None
    try:
        import wandb
    except ImportError:
        print("  [warn] --wandb set but wandb missing; local-only.")
        return None
    if not os.environ.get("WANDB_API_KEY") and WANDB_KEY_FILE.exists():
        os.environ["WANDB_API_KEY"] = WANDB_KEY_FILE.read_text().strip()
    name = (f"{args.snp_set}-{args.model}-{args.aggregation}-"
            f"{args.mlp_depth}-d{args.delta}-s{args.seed}")
    return wandb.init(project=args.wandb_project, entity=args.wandb_entity,
                      name=name, reinit=True, config=config)


def _umap(ax, X, lab, title, classes):
    import umap
    if len(X) < 4:
        ax.set_title(f"{title} (n<4)")
        return
    nn = min(15, len(X) - 1)
    z = umap.UMAP(n_components=2, n_neighbors=nn, min_dist=0.3,
                  metric="euclidean", random_state=42).fit_transform(X)
    for c in classes:
        m = np.asarray(lab) == c
        if m.any():
            ax.scatter(z[m, 0], z[m, 1], s=14, label=str(c), alpha=0.75)
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    ax.set_xticks([]); ax.set_yticks([])


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snp-set", required=True)
    ap.add_argument("--model", required=True, choices=list(fl.MODEL_DIMS))
    ap.add_argument("--aggregation", required=True,
                    choices=["global_attn", "chrom_hier"])
    ap.add_argument("--mlp-depth", required=True, choices=["2L", "3L"])
    ap.add_argument("--delta", required=True, choices=["off", "on"])
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--snp-encoder", action="store_true")
    ap.add_argument("--gwas-filtered", type=Path, default=fl.GWAS_FILTERED)
    ap.add_argument("--diffs-root", type=Path,
                    default=fl.BASE / "fm_embeddings_short_seq_1kb",
                    help="flat layout: <diffs-root>/<model>/<set>_snp_embeddings.npz "
                         "(matches the Drive structure script 27 writes; "
                         "downloaded locally to "
                         "D:/ADNI_SNP_Omni2.5M_20140220/fm_embeddings_short_seq_1kb/)")
    ap.add_argument("--conv-labels", type=Path, default=fl.CONV_LABELS)
    ap.add_argument("--splits-root", type=Path, default=fl.SPLITS_ROOT)
    ap.add_argument("--output-root", type=Path, required=True)
    ap.add_argument("--epochs", type=int, default=300)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--weight-decay", type=float, default=1e-4)
    ap.add_argument("--patience", type=int, default=30)
    ap.add_argument("--loss", choices=["bce", "weighted_bce"],
                    default="weighted_bce")
    ap.add_argument("--attn-hidden", type=int, default=128)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--wandb", action="store_true")
    ap.add_argument("--wandb-project", default="adni-snp-diff-attn")
    ap.add_argument("--wandb-entity", default=None)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()

    if args.smoke:
        args.epochs, args.patience = 3, 2

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    import random
    random.seed(args.seed)
    device = torch.device(args.device)

    sset = args.snp_set
    sdir = args.gwas_filtered / sset
    out_dir = (args.output_root / sset / args.model / args.aggregation /
               args.mlp_depth / f"delta_{args.delta}" / f"seed_{args.seed}")
    if (out_dir / "metrics.json").exists() and not args.force:
        print(f"[skip] exists: {out_dir}")
        return
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── SNP table → canonical order → drop QC-failed (manifest authority) ──
    snp = fl.canonical_snp_order(fl.load_snp_table(
        sdir / f"{sset}_snps.tsv"))
    man = pd.read_csv(
        sdir / "fm_sequences" / f"{sset}_snp_manifest.tsv", sep="\t")
    keep_rs = set(man.loc[~man["qc_verdict"].astype(str)
                          .str.startswith("drop_"), "rsID"].astype(str))
    snp = snp[snp["rsID"].astype(str).isin(keep_rs)].reset_index(drop=True)

    # Drive-folder alias: caduceus npz dirs carry a `_d256` suffix on Drive
    # (e.g. `caduceus_ph_d256/`) while the model key in MODEL_DIMS is the
    # bare `caduceus_ph`. Map driver-side names → on-disk dir names. Default
    # is identity so bmfm_ref / bmfm_snp / ntv2 are unaffected.
    _MODEL_DIR_ALIAS = {"caduceus_ph": "caduceus_ph_d256",
                         "caduceus_ps": "caduceus_ps_d256"}
    model_dir = _MODEL_DIR_ALIAS.get(args.model, args.model)
    diff_npz = (args.diffs_root / model_dir /
                f"{sset}_snp_embeddings.npz")
    diff, rs_emb = fl.load_diff(diff_npz, snp)             # asserts alignment
    H = diff.shape[1]
    assert H == fl.MODEL_DIMS[args.model], \
        f"H {H} != {fl.MODEL_DIMS[args.model]} for {args.model}"

    dosage = fl.load_dosage(sdir / f"{sset}_patient_dosage.tsv")
    snp_rs = snp["rsID"].astype(str).tolist()
    have = [r for r in snp_rs if r in dosage.columns]
    if len(have) != len(snp_rs):
        miss = sorted(set(snp_rs) - set(have))
        print(f"  [warn] {len(miss)} SNPs not in dosage matrix → dropped: "
              f"{miss[:5]}")
        keep = snp["rsID"].astype(str).isin(set(have))
        snp = snp[keep].reset_index(drop=True)
        diff = diff[keep.to_numpy()]
        snp_rs = have
    beta = snp["beta_A1"].to_numpy(np.float32)
    S = len(snp_rs)

    # chromosome grouping (numeric order; 15 build_chrom_map convention)
    chroms = sorted(snp["CHR"].astype(str).unique(), key=fl._chr_key)
    cmap = {c: i for i, c in enumerate(chroms)}
    chrom_of_snp = snp["CHR"].astype(str).map(cmap).to_numpy(np.int64)
    n_chroms = len(chroms)

    cohort = fl.load_cohort_labels(args.conv_labels)
    members = fl.split_membership(args.seed, args.splits_root)
    yb_map = dict(zip(cohort["Patient_ID"], cohort["y"]))

    def assemble(split: str, impute_means):
        pids = [p for p in members[split]
                if p in yb_map and p in dosage.index]
        pids = sorted(pids)
        dsub = dosage.loc[pids, snp_rs]
        dz, means, n_imp = fl.cohort_mean_impute(
            dsub, snp_rs, means=impute_means)
        dose = dz.to_numpy(np.float32)                     # (P,S)
        d = fl.build_d(diff, beta, dose)                   # (P,S,H)
        y = np.array([yb_map[p] for p in pids], np.int64)
        mask = np.ones((len(pids), S), dtype=bool)
        cidx = np.tile(chrom_of_snp, (len(pids), 1))
        ab = np.tile(np.abs(beta), (len(pids), 1)).astype(np.float32)
        adb = np.abs(dose * beta[None, :]).astype(np.float32)
        tens = (torch.from_numpy(d), torch.from_numpy(mask),
                torch.from_numpy(cidx), torch.from_numpy(ab),
                torch.from_numpy(adb), torch.from_numpy(y))
        return pids, tens, means, n_imp

    tr_pids, tr, means, n_imp = assemble("train", None)    # train-fold means
    va_pids, va, _, _ = assemble("val", means)             # no leakage
    te_pids, te, _, _ = assemble("test", means)
    cls = {sp: {"pos": int((t[-1] == 1).sum()),
                "neg": int((t[-1] == 0).sum())}
           for sp, t in (("train", tr), ("val", va), ("test", te))}
    print(f"[data] {sset}/{args.model} S={S} chroms={n_chroms} "
          f"train={len(tr_pids)}{cls['train']} val={len(va_pids)}"
          f"{cls['val']} test={len(te_pids)}{cls['test']}")
    for sp in ("train", "val", "test"):
        n = cls[sp]["pos"] + cls[sp]["neg"]
        if n < 30 or min(cls[sp].values()) < 5:
            print(f"  [warn] small {sp} fold: n={n} {cls[sp]}")

    model = fl.DiffAttentionModel(
        model_dim=H, n_chroms=n_chroms, aggregation=args.aggregation,
        mlp_depth=args.mlp_depth, use_delta=(args.delta == "on"),
        snp_encoder=args.snp_encoder, attn_hidden=args.attn_hidden).to(device)

    config = {"snp_set": sset, "model": args.model,
              "aggregation": args.aggregation, "mlp_depth": args.mlp_depth,
              "delta": args.delta, "snp_encoder": args.snp_encoder,
              "seed": args.seed, "model_dim": H, "n_snps": S,
              "n_chroms": n_chroms, "chroms": chroms,
              "n_patients": {k: len(v) for k, v in
                             (("train", tr_pids), ("val", va_pids),
                              ("test", te_pids))},
              "class_balance": cls, "loss": args.loss,
              "impute": "train_fold_per_snp_mean",
              "hparams": {"epochs": args.epochs, "lr": args.lr,
                          "batch_size": args.batch_size,
                          "weight_decay": args.weight_decay,
                          "patience": args.patience,
                          "attn_hidden": args.attn_hidden},
              "diff_npz": str(diff_npz),
              "seeds_set": {"torch": args.seed, "numpy": args.seed,
                            "python": args.seed}}
    run = init_wandb(args, config)

    model, fit, curve = fl.train_one_diff(
        model, tr, va, epochs=args.epochs, batch_size=args.batch_size,
        lr=args.lr, weight_decay=args.weight_decay, device=device,
        patience=args.patience, loss_name=args.loss)

    m_tr, x_tr = fl.evaluate_diff(model, tr, device)
    m_va, x_va = fl.evaluate_diff(model, va, device)
    m_te, x_te = fl.evaluate_diff(model, te, device)
    gpool = model.pool
    gamma = float((gpool.attn.gamma if hasattr(gpool, "attn")
                   else gpool.intra.gamma).detach().cpu())
    dpar = (gpool.attn.delta if hasattr(gpool, "attn") else gpool.intra.delta)
    delta_v = None if dpar is None else float(dpar.detach().cpu())
    print(f"[fit] {fit}  γ={gamma:.4f} δ={delta_v}  "
          f"VAL bAcc={m_va['balanced_accuracy']:.3f} "
          f"TEST bAcc={m_te['balanced_accuracy']:.3f} "
          f"AUC={m_te['auc']:.3f}")

    # ── write outputs ────────────────────────────────────────────────────
    (out_dir / "config.json").write_text(json.dumps(config, indent=2),
                                         encoding="utf-8")
    json.dump({**{k: config[k] for k in
                  ("snp_set", "model", "aggregation", "mlp_depth", "delta",
                   "seed", "loss")},
               "fit_info": fit, "learned": {"gamma": gamma,
                                            "delta": delta_v},
               "train": m_tr, "val": m_va, "test": m_te},
              open(out_dir / "metrics.json", "w"), indent=2)
    pd.DataFrame(curve).to_csv(out_dir / "val_curve.csv", index=False)
    torch.save(model.state_dict(), out_dir / "model.pt")

    rows = []
    for sp, pids, ex in (("train", tr_pids, x_tr), ("val", va_pids, x_va),
                         ("test", te_pids, x_te)):
        for i, p in enumerate(pids):
            rows.append({"patient_id": p, "split": sp,
                         "true_label": int(ex["y"][i]),
                         "predicted_probability": float(ex["prob"][i]),
                         "predicted_label": int(ex["pred"][i])})
    pd.DataFrame(rows).to_csv(out_dir / "predictions.csv", index=False)

    np.savez_compressed(
        out_dir / "attention.npz", rsids=np.array(snp_rs),
        chrom_of_snp=chrom_of_snp,
        **{f"a_snp_{sp}": ex.get("a_snp")
           for sp, ex in (("train", x_tr), ("val", x_va), ("test", x_te))
           if ex.get("a_snp") is not None},
        **{f"a_chrom_{sp}": ex["a_chrom"]
           for sp, ex in (("train", x_tr), ("val", x_va), ("test", x_te))
           if "a_chrom" in ex})
    np.savez_compressed(
        out_dir / "embeddings.npz",
        **{f"{sp}_pids": np.array(p) for sp, p in
           (("train", tr_pids), ("val", va_pids), ("test", te_pids))},
        **{f"{sp}_emb": ex["patient_emb"] for sp, ex in
           (("train", x_tr), ("val", x_va), ("test", x_te))},
        **{f"{sp}_y": ex["y"] for sp, ex in
           (("train", x_tr), ("val", x_va), ("test", x_te))})

    # ── debug UMAPs ──────────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        umd = out_dir / "umaps"
        umd.mkdir(exist_ok=True)
        # (3) patient embeddings AD/CN (all splits)
        Pe = np.concatenate([x_tr["patient_emb"], x_va["patient_emb"],
                             x_te["patient_emb"]])
        Py = np.concatenate([x_tr["y"], x_va["y"], x_te["y"]])
        fig, ax = plt.subplots(figsize=(5, 4.5))
        _umap(ax, Pe, ["AD" if v == 1 else "CN" for v in Py],
              "Patient embeddings (AD vs CN)", ["CN", "AD"])
        fig.savefig(umd / "patient_emb.png", dpi=140,
                    bbox_inches="tight"); plt.close(fig)
        # (2) per-SNP diff coloured by CHR
        fig, ax = plt.subplots(figsize=(5, 4.5))
        _umap(ax, diff, [str(c) for c in snp["CHR"]],
              "Per-SNP diff (ALT−REF) by CHR", sorted(set(snp["CHR"]),
              key=fl._chr_key))
        fig.savefig(umd / "dbg_diff.png", dpi=140,
                    bbox_inches="tight"); plt.close(fig)
        # (1) raw ALT vs REF (needs ref/alt in npz)
        z = np.load(diff_npz, allow_pickle=True)
        if "ref_emb" in z and "alt_emb" in z:
            X = np.concatenate([np.asarray(z["ref_emb"]),
                                np.asarray(z["alt_emb"])])
            lab = ["REF"] * S + ["ALT"] * S
            fig, ax = plt.subplots(figsize=(5, 4.5))
            _umap(ax, X, lab, "Raw REF vs ALT embeddings", ["REF", "ALT"])
            fig.savefig(umd / "dbg_altref.png", dpi=140,
                        bbox_inches="tight"); plt.close(fig)
    except Exception as e:
        print(f"  [warn] UMAP step skipped: {e!r}")

    if run is not None:
        import wandb
        wandb.log({f"test/{k}": v for k, v in m_te.items()
                   if isinstance(v, (int, float))})
        wandb.finish()
    print(f"[out] {out_dir}")


if __name__ == "__main__":
    main()
