"""
02_run_grid.py
==============
Sequential driver for the full cross-modal fusion grid (local CPU, env: snp).

Fusion grid: clin_arm {BL,M12,LONG} x variant {A,B} x arch {mlp_concat,cross_attn,cross_attn_x}
                  x loss {ce,ce_patient,ce_patient_label} x seed {0,1,2}.
cross_attn_x = true bidirectional cross-attention over single CLS vectors (LONG uses --clin_pool
mean here; the frozen-Mamba pool is a separate clinical-env Phase-2 driver).
Clinical-only references (same cohort): clin_arm {BL,M12,LONG} x variant {A,B}
                  x arch {mlp_concat,cross_attn} x seed {0,1,2}  (--clin_only, loss forced ce;
                  cross_attn_x excluded — it requires the MRI modality).
The clin_only refs are run per-variant because each variant's cohort drops its own NaN scans
(here identical), keeping refs comparable to their fusion runs.

Each sub-run is invoked as a subprocess of 01_train_xmodal_fusion.py and skips when its
metrics.json already exists (resumable). Tiny models -> the whole grid is minutes on CPU.

Run:  python integration/cross_modal_attention/02_run_grid.py [--overwrite]
"""
import argparse
import itertools
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TRAIN = os.path.join(HERE, "01_train_xmodal_fusion.py")

ARMS = ["BL", "M12", "LONG"]
VARIANTS = ["A", "B"]
ARCHS = ["mlp_concat", "cross_attn", "cross_attn_x"]
ARCHS_REF = ["mlp_concat", "cross_attn"]     # clin_only refs only — cross_attn_x requires MRI
LOSSES = ["ce", "ce_patient", "ce_patient_label"]
SEEDS = [0, 1, 2]


def run(extra, overwrite):
    cmd = [sys.executable, TRAIN] + extra + (["--overwrite"] if overwrite else [])
    r = subprocess.run(cmd, capture_output=True, text=True)
    tail = [l for l in r.stdout.splitlines() if l.startswith(("[done]", "[skip]", "  TEST"))]
    print("  " + " | ".join(tail[-2:]) if tail else f"  (rc={r.returncode})")
    if r.returncode != 0:
        print("  STDERR:", r.stderr.strip().splitlines()[-1] if r.stderr.strip() else "")
    return r.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    fusion = list(itertools.product(ARMS, VARIANTS, ARCHS, LOSSES, SEEDS))
    refs = list(itertools.product(ARMS, VARIANTS, ARCHS_REF, SEEDS))
    print(f"Fusion runs: {len(fusion)} | clinical-only refs: {len(refs)} | "
          f"total: {len(fusion) + len(refs)}")

    fails = 0
    for i, (arm, var, arch, loss, seed) in enumerate(fusion, 1):
        print(f"[{i}/{len(fusion)}] fusion {arm} {var} {arch} {loss} seed{seed}")
        fails += bool(run(["--clin_arm", arm, "--variant", var, "--arch", arch,
                           "--loss", loss, "--seed", str(seed)], args.overwrite))
    for i, (arm, var, arch, seed) in enumerate(refs, 1):
        print(f"[ref {i}/{len(refs)}] clin_only {arm} {var} {arch} seed{seed}")
        fails += bool(run(["--clin_arm", arm, "--variant", var, "--arch", arch,
                           "--loss", "ce", "--clin_only", "--seed", str(seed)], args.overwrite))
    print(f"\nDone. failures: {fails}")


if __name__ == "__main__":
    main()
