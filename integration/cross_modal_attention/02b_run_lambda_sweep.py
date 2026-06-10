"""
02b_run_lambda_sweep.py
=======================
Contrastive-weight sweep on the longitudinal (LONG) and 1-year (M12) clinical arms. Gives the
SigLIP auxiliary losses a fair chance via:
  - warmup=10     : train aux(SigLIP)-only for 10 epochs first (organize the embedding space)
  - min_epochs=30 : block early-stopping so contrastive-shaped epochs can win selection
  - lambda=gamma in {0.1, 0.5, 1.0}  (vs the lambda=0.01 main-grid baseline)

Grid (72): LONG x variant {A,B} x arch {mlp_concat,cross_attn} x loss {ce_patient,ce_patient_label}
           x lambda {0.1,0.5,1.0} x seed {0,1,2}.  (loss 'ce' has no aux term, so it is not swept.)
Outputs isolated under outputs/sweep/lam<λ>_wu10/ so the main grid is untouched.

Run:  python integration/cross_modal_attention/02b_run_lambda_sweep.py [--overwrite]
"""
import argparse
import itertools
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TRAIN = os.path.join(HERE, "01_train_xmodal_fusion.py")

ARMS = ["LONG", "M12"]
VARIANTS = ["A", "B"]
ARCHS = ["mlp_concat", "cross_attn", "cross_attn_x"]
LOSSES = ["ce_patient", "ce_patient_label"]
LAMBDAS = [0.1, 0.5, 1.0]
SEEDS = [0, 1, 2]
WARMUP, MIN_EPOCHS = 10, 30


def run(extra, overwrite):
    cmd = [sys.executable, TRAIN] + extra + (["--overwrite"] if overwrite else [])
    r = subprocess.run(cmd, capture_output=True, text=True)
    tail = [l for l in r.stdout.splitlines() if l.startswith(("[done]", "[skip]", "  TEST"))]
    print("  " + " | ".join(tail[-2:]) if tail else f"  (rc={r.returncode})")
    if r.returncode:
        print("  STDERR:", (r.stderr.strip().splitlines() or [""])[-1])
    return r.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()
    grid = list(itertools.product(ARMS, VARIANTS, ARCHS, LOSSES, LAMBDAS, SEEDS))
    print(f"Sweep runs: {len(grid)}  (arms {ARMS}, warmup={WARMUP}, min_epochs={MIN_EPOCHS})")
    fails = 0
    for i, (arm, var, arch, loss, lam, seed) in enumerate(grid, 1):
        tag = f"sweep/lam{lam}_wu{WARMUP}"
        print(f"[{i}/{len(grid)}] {arm} {var} {arch} {loss} lam{lam} seed{seed}")
        fails += bool(run(["--clin_arm", arm, "--variant", var, "--arch", arch,
                           "--loss", loss, "--seed", str(seed),
                           "--lam", str(lam), "--gam", str(lam),
                           "--warmup", str(WARMUP), "--min_epochs", str(MIN_EPOCHS),
                           "--run_tag", tag], args.overwrite))
    print(f"\nDone. failures: {fails}")


if __name__ == "__main__":
    main()
