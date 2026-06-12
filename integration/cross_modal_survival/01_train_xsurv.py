r"""
01_train_xsurv.py  (env: snp or clinical — torch, no transformers needed for mean pooling)
==========================================================================================
Train the cross-modal survival models and export per-patient survival curves (T4 schema) for scoring by
02_score_xsurv.py (env survml). Per T1d MRI variant × arm × seed:

  arms:  clinical  — clinical-only mean-pool → Weibull-piecewise (reference, same cohort)
         mri       — T1d MRI-only → Weibull-piecewise (reference)
         none      — fused (cross-attention), survival-NLL only
         joint     — fused + λ·patient-SigLIP (joint alignment)
         prealign  — fused, warmup patient-SigLIP only → then NLL (+λ·SigLIP)

Matrix: variants {mvp (BrainMVP-T1d 512), dino (BrainDINO-T1d 768)} × 5 arms × seeds {0,1,2}.
Output: outputs/xsurv/<variant>/<arm>/seed_{s}/{predictions.parquet, meta.json}.  Resumable.

Run:  conda run -n snp python integration/cross_modal_survival/01_train_xsurv.py [--seeds 0 1 2]
      [--variants mvp dino] [--arms clinical mri none joint prealign] [--overwrite]
"""
import argparse
import sys
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _xsurv_lib as X

OUT = HERE / "outputs" / "xsurv"
ARMS = ["clinical", "mri", "none", "joint", "prealign"]
ARM_MODE = {"clinical": ("clinical", "none"), "mri": ("mri", "none"),
            "none": ("fused", "none"), "joint": ("fused", "joint"), "prealign": ("fused", "prealign")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--variants", type=str, nargs="+", default=["mvp", "dino"])
    ap.add_argument("--arms", type=str, nargs="+", default=ARMS)
    ap.add_argument("--device", type=str, default="cpu")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    master = X.T4.load_master_labels(X.MASTER)
    print(f"variants={args.variants} arms={args.arms} seeds={args.seeds}")
    fails = 0
    for variant in args.variants:
        for seed in args.seeds:
            cohort = X.build_surv_cohort(seed, variant, master_labels=master)
            r = cohort["report"]
            print(f"[{variant} seed{seed}] cohort n={r['n']} events={r['events']} "
                  f"(clin_total={r['n_clin_total']}, mri_total={r['n_mri_total']}, "
                  f"split_mismatch_mri={r['split_mismatch_mri']})")
            for arm in args.arms:
                mode, align = ARM_MODE[arm]
                out_dir = OUT / variant / arm / f"seed_{seed}"
                if (out_dir / "predictions.parquet").exists() and not args.overwrite:
                    print(f"  [skip] {variant}/{arm}/seed{seed}")
                    continue
                try:
                    vnll, _ = X.train_export(cohort, mode, align, seed, out_dir, device=args.device)
                    print(f"  [done] {variant}/{arm}/seed{seed}  val_nll={vnll:.4f}")
                except Exception:
                    fails += 1
                    print(f"  [FAIL] {variant}/{arm}/seed{seed}\n{traceback.format_exc()}")
    print(f"\nDone. failures: {fails}")


if __name__ == "__main__":
    main()
