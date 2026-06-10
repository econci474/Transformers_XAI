"""
02c_run_long_mamba.py
=====================
PHASE 2 (env: clinical, NOT snp) — the LONG-arm frozen-Mamba pooling variant of the true
bidirectional CLS cross-attention (`cross_attn_x`), to compare against the Phase-1 mean pooling.

Why a separate driver: the frozen-Mamba pooler reuses the T4 `SeqEncoder` (`integration/T4/
mamba_long_clinical/_mamba_seq_lib.py`), which loads `state-spaces/mamba-130m-hf` via
`transformers`. The `snp` env has no `transformers`, so this cannot run there. Run it in the
`clinical` env (transformers + HF Mamba); CPU is fine via the pure-PyTorch `slow_forward`
(no `mamba_ssm`/CUDA kernels needed). First run downloads the ~130M Mamba weights to HF_HOME.

Mirrors the LONG·cross_attn_x cells of 02 (main grid, lambda=0.01) and 02b (sweep, warmup=10,
min_epochs=30, lambda in {0.1,0.5,1.0}) but with --clin_pool mamba, so LONG·mean vs LONG·mamba are
directly comparable in 03/03b/03c (each writes a separate `…/LONG/<variant>/cross_attn_x/<pool>/…`
tree). Resumable: skips when metrics.json exists.

Run (env: clinical):
  conda run -n clinical python integration/cross_modal_attention/02c_run_long_mamba.py [--overwrite]
"""
import argparse
import itertools
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TRAIN = os.path.join(HERE, "01_train_xmodal_fusion.py")

SEEDS = [0, 1, 2]
MAIN_LOSSES = ["ce", "ce_patient", "ce_patient_label"]      # lambda=0.01, no warmup (matches 02)
SWEEP_LOSSES = ["ce_patient", "ce_patient_label"]           # lambda swept (matches 02b)
LAMBDAS = [0.1, 0.5, 1.0]
WARMUP, MIN_EPOCHS = 10, 30
COMMON = ["--clin_arm", "LONG", "--arch", "cross_attn_x", "--clin_pool", "mamba"]


def run(extra, overwrite):
    cmd = [sys.executable, TRAIN] + COMMON + extra + (["--overwrite"] if overwrite else [])
    r = subprocess.run(cmd, capture_output=True, text=True)
    tail = [l for l in r.stdout.splitlines() if l.startswith(("[done]", "[skip]", "  TEST"))]
    print("  " + " | ".join(tail[-2:]) if tail else f"  (rc={r.returncode})")
    if r.returncode:
        print("  STDERR:", (r.stderr.strip().splitlines() or [""])[-1])
    return r.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--variants", nargs="+", default=["A", "B"], choices=["A", "B"],
                    help="MRI variants to run (A=BrainDINO-T2, B=BrainMVP-T1b). Default both.")
    args = ap.parse_args()
    variants = args.variants

    main_grid = list(itertools.product(variants, MAIN_LOSSES, SEEDS))
    sweep = list(itertools.product(variants, SWEEP_LOSSES, LAMBDAS, SEEDS))
    print(f"LONG/mamba runs: main {len(main_grid)} (lambda=0.01) + sweep {len(sweep)} "
          f"(warmup={WARMUP}, min_epochs={MIN_EPOCHS}) = {len(main_grid) + len(sweep)}")

    fails = 0
    for i, (var, loss, seed) in enumerate(main_grid, 1):
        print(f"[main {i}/{len(main_grid)}] LONG mamba {var} {loss} seed{seed}")
        fails += bool(run(["--variant", var, "--loss", loss, "--seed", str(seed)], args.overwrite))
    for i, (var, loss, lam, seed) in enumerate(sweep, 1):
        tag = f"sweep/lam{lam}_wu{WARMUP}"
        print(f"[sweep {i}/{len(sweep)}] LONG mamba {var} {loss} lam{lam} seed{seed}")
        fails += bool(run(["--variant", var, "--loss", loss, "--seed", str(seed),
                           "--lam", str(lam), "--gam", str(lam),
                           "--warmup", str(WARMUP), "--min_epochs", str(MIN_EPOCHS),
                           "--run_tag", tag], args.overwrite))
    print(f"\nDone. failures: {fails}")


if __name__ == "__main__":
    main()
