"""
audit_sweeps.py
===============
Completion audit for the cached-head HP sweeps + the BrainMVP T4 full_ft job. For each submit script
it enumerates the exact HP grid, checks whether each cell's `metrics.json` exists, and reports
found/expected/% complete plus the SLURM **array indices** of any missing cells (so you can resubmit
just those -- though every submit also has a `metrics.json` skip-guard, so a plain re-`sbatch` only
re-runs the missing cells anyway).

Run on CSD3 (where the outputs live):
    python mri_pipeline/cached_head_sweep/audit_sweeps.py
Or after scp'ing the trees to local D:, point --base at the local root:
    python mri_pipeline/cached_head_sweep/audit_sweeps.py --base D:/ADNI_BIDS_project/derivatives/ADNI_MRI

Array decode (matches the submit scripts):
  T3abcd : id = TASK + 4*(SEED + 3*(DROP + 3*(LS + 2*LR)))   [N_TASKS=4]
  T4/T1e : id = SEED + 3*(DROP + 3*(LS + 2*LR))              [single task]
  04i    : id = SEED                                          [per AUGMENT submission]
HP dir name: lr{lr:.0e}_d{drop}_ls{ls}  e.g. lr1e-03_d0.1_ls0.0
"""
from __future__ import annotations

import argparse
import os

# CSD3 base that all OUT_DIRs share; override with --base after scp to local.
DEFAULT_BASE = "/home/ec474/rds/hpc-work/ADNI_MRI"

SEEDS = [0, 1, 2]
DROPS = [0.1, 0.2, 0.3]
LSS = [0.0, 0.1]
LRS = [1e-3, 1e-4, 1e-5]
T3_TASKS = ["T3a_conv3y", "T3b_conv5y", "T3c_conv7y", "T3d_conv10y"]

# Cached-head OUT_DIRs (relative to --base)
CACHED = {
    "brainmvp": "brainmvp_debug/aug_none/BrainMVP_uniformer_frozen_cached",
    "braindino": "braindino_outputs/aug_none_hp_tuned/BrainDINO_vitb16_frozen_cached",
    "vit_mae": "vit_outputs_debug/aug_none/ViT_B_mae75_frozen_cached",
}


def _hp_dir(lr, drop, ls):
    return f"lr{lr:.0e}_d{drop}_ls{ls}"


def _cached_cells(tasks, n_tasks_for_decode):
    """Yield (array_id, rel_runpath) for a cached-head sweep grid."""
    for ti, task in enumerate(tasks):
        for si, seed in enumerate(SEEDS):
            for di, drop in enumerate(DROPS):
                for li, ls in enumerate(LSS):
                    for ri, lr in enumerate(LRS):
                        if n_tasks_for_decode == 1:
                            aid = si + 3 * (di + 3 * (li + 2 * ri))
                        else:
                            aid = ti + n_tasks_for_decode * (si + 3 * (di + 3 * (li + 2 * ri)))
                        rel = f"{task}/seed_{seed}/{_hp_dir(lr, drop, ls)}/metrics.json"
                        yield aid, rel


def _audit(label, out_dir, cells):
    found, missing = 0, []
    for aid, rel in cells:
        if os.path.isfile(os.path.join(out_dir, rel)):
            found += 1
        else:
            missing.append(aid)
    total = found + len(missing)
    pct = 100.0 * found / total if total else 0.0
    status = "COMPLETE" if not missing else f"{len(missing)} MISSING"
    print(f"\n[{label}]  {found}/{total}  ({pct:5.1f}%)  {status}")
    print(f"    dir: {out_dir}")
    if missing:
        ids = sorted(set(missing))
        print(f"    missing array ids ({len(ids)}): {','.join(map(str, ids))}")
        print(f"    resubmit: sbatch --array={','.join(map(str, ids))}%4 <that submit script>")
    return found, total


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=DEFAULT_BASE,
                    help="Root containing brainmvp_debug/, braindino_outputs/, vit_outputs_debug/ "
                         "(CSD3 default; override after scp).")
    ap.add_argument("--remote", default="ec474@login-cpu.hpc.cam.ac.uk",
                    help="user@host for the transfer commands printed after the audit.")
    ap.add_argument("--local-dest", default="/d/ADNI_BIDS_project/derivatives/ADNI_MRI",
                    help="Local destination root for the printed transfer commands.")
    args = ap.parse_args()
    b = args.base

    print("=" * 78)
    print(f"  Cached-head + full_ft sweep audit   base={b}")
    print("=" * 78)

    specs = []
    # T3abcd (216 each), 3 archs
    for arch, rel in CACHED.items():
        specs.append((f"T3abcd / {arch}", os.path.join(b, rel),
                      list(_cached_cells(T3_TASKS, 4))))
    # T4 (54 each), brainmvp + braindino
    for arch in ("brainmvp", "braindino"):
        specs.append((f"T4 / {arch}", os.path.join(b, CACHED[arch]),
                      list(_cached_cells(["T4_conv_horizon"], 1))))
    # T1e (54 each), brainmvp + braindino
    for arch in ("brainmvp", "braindino"):
        specs.append((f"T1e / {arch}", os.path.join(b, CACHED[arch]),
                      list(_cached_cells(["T1e_pcn_vs_scn"], 1))))
    # 04i BrainMVP T4 full_ft: 3 seeds per AUGMENT
    for aug in ("none", "stochastic", "plus_original"):
        cells = [(s, f"BrainMVP_uniformer/T4_conv_horizon/seed_{s}/full_ft/metrics.json")
                 for s in SEEDS]
        specs.append((f"T4 full_ft / brainmvp aug={aug}",
                      os.path.join(b, f"brainmvp_debug/aug_{aug}"), cells))

    grand_f = grand_t = 0
    for label, out_dir, cells in specs:
        f, t = _audit(label, out_dir, cells)
        grand_f += f; grand_t += t

    print("\n" + "=" * 78)
    print(f"  TOTAL: {grand_f}/{grand_t} cells complete "
          f"({100.0 * grand_f / grand_t if grand_t else 0:.1f}%)")
    print("=" * 78)

    # ---- transfer commands: pull each task subtree CSD3 -> local (run LOCALLY) ----
    host, dest = args.remote, args.local_dest
    print("\n" + "=" * 78)
    print("  TRANSFER COMMANDS  (run from your LOCAL machine; rsync = resumable, recommended)")
    print(f"  remote = {host}   local dest root = {dest}")
    print("=" * 78)

    def _xfer(src_rel, note=""):
        src = f"{b}/{src_rel}"
        dparent = f"{dest}/{os.path.dirname(src_rel)}"
        print(f"\n  # {src_rel}{('   ' + note) if note else ''}")
        print(f"  mkdir -p {dparent}")
        print(f"  rsync -av {host}:{src}/ {dest}/{src_rel}/")
        print(f"  # scp alt: scp -r {host}:{src} {dparent}/")

    # cached-head sweeps: one subtree per (arch, task) we need for fusion
    for arch, rel in CACHED.items():
        tasks = T3_TASKS + (["T4_conv_horizon", "T1e_pcn_vs_scn"] if arch != "vit_mae" else [])
        for task in tasks:
            _xfer(f"{rel}/{task}")
    # BrainMVP T4 full_ft, the three augmentation folders
    for aug in ("none", "stochastic", "plus_original"):
        _xfer(f"brainmvp_debug/aug_{aug}/BrainMVP_uniformer/T4_conv_horizon",
              note=f"(T4 full_ft, aug={aug})")
    print()


if __name__ == "__main__":
    main()
