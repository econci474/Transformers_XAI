"""
01_train_xmodal_fusion.py
=========================
Train ONE cross-modal feature-fusion run for the m12 three-way diagnosis (CN/MCI/AD).

Fuses a clinical embedding (BioClin-ModernBERT-L, T2, full_ft) with an m12 MRI embedding using
either a cross-modal attention block or an MLP-concat head, optionally with SigLIP auxiliary
losses (patient alignment + label/hierarchical). See _xmodal_lib.py for the data/model/loss core
and the approved plan for the full design.

Cohort = common intersection (bl+m06+m12 clinical AND m12 MRI), so BL and LONG clinical arms are
compared on identical patients. Target = m12 concurrent label.

Run (env: snp, local CPU):
  python integration/cross_modal_attention/01_train_xmodal_fusion.py \
      --clin_arm BL --variant A --arch mlp_concat --loss ce --seed 0
  # add --clin_only to drop the MRI token (clinical-only reference, same cohort)

Out: outputs/<clin_arm[/clin_only]>/<variant>/<arch>/<loss>/seed_<s>/metrics.json + predictions.parquet
"""
import argparse
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _xmodal_lib as L  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="Cross-modal feature fusion (m12 T2).")
    ap.add_argument("--clin_arm", choices=["BL", "M12", "LONG"], required=True,
                    help="BL=baseline model, M12=dedicated 1-year model, LONG=all-visits model")
    ap.add_argument("--variant", choices=["A", "B"], required=True,
                    help="A=BrainDINO-T2 (768d), B=BrainMVP-T1b (512d)")
    ap.add_argument("--arch", choices=["mlp_concat", "cross_attn", "cross_attn_x"], required=True,
                    help="cross_attn_x = true bidirectional cross-attention over single CLS vectors")
    ap.add_argument("--loss", choices=["ce", "ce_patient", "ce_patient_label"], required=True)
    ap.add_argument("--clin_pool", choices=["mean", "mamba"], default="mean",
                    help="LONG-arm pooling of the 3 visit CLS for cross_attn_x: mean (snp/CPU) or "
                         "frozen-Mamba (clinical env only). Ignored for BL/M12 and other arches.")
    ap.add_argument("--seed", type=int, choices=[0, 1, 2], required=True)
    ap.add_argument("--clin_only", action="store_true",
                    help="drop the MRI token (clinical-only reference on the same cohort)")
    ap.add_argument("--lam", type=float, default=0.01, help="patient-SigLIP weight")
    ap.add_argument("--gam", type=float, default=0.01, help="label-SigLIP weight")
    ap.add_argument("--warmup", type=int, default=0,
                    help="train aux(SigLIP)-only for this many epochs before adding CE")
    ap.add_argument("--min_epochs", type=int, default=0,
                    help="block early-stopping until this epoch (let contrastive-shaped epochs win)")
    ap.add_argument("--run_tag", default="",
                    help="optional subdir under out_root (e.g. 'sweep/lam0.1_wu10') to isolate a sweep")
    ap.add_argument("--out_root", default=os.path.join(HERE, "outputs"))
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    # clin_only forces loss=ce (no MRI -> no patient/label cross-modal terms)
    loss = "ce" if args.clin_only else args.loss
    # clin_pool only varies the LONG arm under cross_attn_x; pin to mean elsewhere so paths/configs
    # stay canonical (BL/M12 have a single visit; other arches ignore pooling).
    clin_pool = args.clin_pool if (args.arch == "cross_attn_x" and args.clin_arm == "LONG") else "mean"
    arm_dir = f"{args.clin_arm}_clinonly" if args.clin_only else args.clin_arm
    root = os.path.join(args.out_root, *args.run_tag.split("/")) if args.run_tag else args.out_root
    seg = [root, arm_dir, args.variant, args.arch]
    if args.arch == "cross_attn_x" and args.clin_arm == "LONG":
        seg.append(clin_pool)                       # …/LONG/<variant>/cross_attn_x/<pool>/<loss>/seed_<s>
    seg += [loss, f"seed_{args.seed}"]
    out_dir = os.path.join(*seg)
    mfile = os.path.join(out_dir, "metrics.json")
    if os.path.exists(mfile) and not args.overwrite:
        print(f"[skip] exists: {mfile}")
        return

    data = L.build_cohort(args.seed, args.clin_arm, args.variant, clin_only=args.clin_only)
    rep = data["report"]
    print(f"[cohort] {args.clin_arm}{'/clin_only' if args.clin_only else ''} "
          f"variant={args.variant} seed={args.seed} n={rep['n']} | integrity={rep}")
    # Hard fail on split/label leakage (must be 0) — see plan verification step 2.
    bad = sum(v for k, v in rep.items() if k.startswith("split_mismatch_"))
    if bad:
        raise RuntimeError(f"split mismatch across modalities ({bad}) — cohort not aligned")
    lm = rep["label_mismatch_mri_vs_clin_m12"]
    if isinstance(lm, int) and lm:
        print(f"  [warn] m12 MRI label != CL_long m12 label for {lm} patients "
              f"(target uses CL_long m12).")

    out, fit = L.train_eval(data, args.arch, loss, args.seed, lam=args.lam, gam=args.gam,
                            warmup=args.warmup, min_epochs=args.min_epochs, clin_pool=clin_pool)

    # predictions.parquet (test fold)
    yte, pte = out["test"]
    test_mask = data["split"] == "test"
    preds = pd.DataFrame({
        "Patient_ID": [p for p, m in zip(data["pids"], test_mask) if m],
        "split": "test", "y_true": yte, "pred": pte.argmax(1),
        "prob_0": pte[:, 0], "prob_1": pte[:, 1], "prob_2": pte[:, 2],
    })
    config = {
        "clin_arm": args.clin_arm, "clin_only": args.clin_only,
        "variant": args.variant, "mri_model": L.MRI_NAME[args.variant],
        "arch": args.arch, "loss": loss, "seed": args.seed, "clin_pool": clin_pool,
        "lambda": args.lam, "gam": args.gam,
        "warmup": args.warmup, "min_epochs": args.min_epochs, "run_tag": args.run_tag,
        "n_train": int((data["split"] == "train").sum()),
        "n_val": int((data["split"] == "val").sum()),
        "n_test": int((data["split"] == "test").sum()),
        "integrity": rep,
    }
    blob = L.write_metrics(out_dir, config, out["val"], out["test"], fit,
                           data["split"], preds)
    vm, tm = blob["val_metrics"], blob["test_metrics"]
    print(f"[done] {out_dir}")
    print(f"  VAL  bACC {vm['balanced_acc']}  AUC {vm['auc_roc_ovr']}  macroF1 {vm['macro_f1']}")
    print(f"  TEST bACC {tm['balanced_acc']}  AUC {tm['auc_roc_ovr']}  macroF1 {tm['macro_f1']}")


if __name__ == "__main__":
    main()
