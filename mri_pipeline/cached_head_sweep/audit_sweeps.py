"""
audit_sweeps.py
===============
Per-TASK completion audit for the cached-head HP sweeps + the BrainMVP T4 full_ft job. For each task
it checks whether every cell's `metrics.json` exists and prints, inline:
  - if INCOMPLETE -> the exact `sbatch --array=<missing ids>%4 <submit script>` to resubmit just those
    (the `metrics.json` skip-guard also makes a plain re-`sbatch` idempotent), OR
  - if COMPLETE   -> the `rsync`/`scp` command to pull that task's outputs CSD3 -> local D:.
So for each task you either RESUBMIT or TRANSFER.

Run on CSD3 (where the outputs live):
    python mri_pipeline/cached_head_sweep/audit_sweeps.py
Override paths if needed: --base (CSD3 root), --remote (user@host), --local-dest (local derivatives/).

Array decode (matches the submit scripts):
  T3abcd : id = TASK + 4*(SEED + 3*(DROP + 3*(LS + 2*LR)))   [N_TASKS=4]
  T4/T1e : id = SEED + 3*(DROP + 3*(LS + 2*LR))              [single task]
  04i    : id = SEED                                          [per AUGMENT submission]
HP dir name: lr{lr:.0e}_d{drop}_ls{ls}  e.g. lr1e-03_d0.1_ls0.0
"""
from __future__ import annotations

import argparse
import os

DEFAULT_BASE = "/home/ec474/rds/hpc-work/ADNI_MRI"     # CSD3 root shared by all OUT_DIRs
SEEDS = [0, 1, 2]
DROPS = [0.1, 0.2, 0.3]
LSS = [0.0, 0.1]
LRS = [1e-3, 1e-4, 1e-5]
T3_TASKS = ["T3a_conv3y", "T3b_conv5y", "T3c_conv7y", "T3d_conv10y"]

CACHED = {  # cached-head OUT_DIR (relative to --base) per arch
    "brainmvp": "brainmvp_debug/aug_none/BrainMVP_uniformer_frozen_cached",
    "braindino": "braindino_outputs/aug_none_hp_tuned/BrainDINO_vitb16_frozen_cached",
    "vit_mae": "vit_outputs_debug/aug_none/ViT_B_mae75_frozen_cached",
}
SUB = "mri_pipeline/cached_head_sweep"                  # cached-head submit-script dir


def _hp_dir(lr, drop, ls):
    return f"lr{lr:.0e}_d{drop}_ls{ls}"


def _cells_by_task(tasks, n_tasks):
    """Return {task: [(array_id, rel_metrics_path)]} for a cached-head sweep grid."""
    out = {t: [] for t in tasks}
    for ti, task in enumerate(tasks):
        for si, seed in enumerate(SEEDS):
            for di, drop in enumerate(DROPS):
                for li, ls in enumerate(LSS):
                    for ri, lr in enumerate(LRS):
                        inner = si + 3 * (di + 3 * (li + 2 * ri))
                        aid = inner if n_tasks == 1 else ti + n_tasks * inner
                        out[task].append((aid, f"{task}/seed_{seed}/{_hp_dir(lr, drop, ls)}/metrics.json"))
    return out


def _xfer(b, host, dest, src_rel, dest_rel=None):
    """rsync (+ scp alt) lines to pull a subtree CSD3 -> local. dest_rel overrides for odd layouts."""
    dest_rel = dest_rel or src_rel
    dparent = f"{dest}/{os.path.dirname(dest_rel)}"
    return [f"mkdir -p {dparent}",
            f"rsync -av {host}:{b}/{src_rel}/ {dest}/{dest_rel}/",
            f"# scp alt: scp -r {host}:{b}/{src_rel} {dparent}/"]


def _audit(label, out_dir, cells, submit, xfer_lines):
    found = sum(os.path.isfile(os.path.join(out_dir, rel)) for _, rel in cells)
    missing = [aid for aid, rel in cells if not os.path.isfile(os.path.join(out_dir, rel))]
    total = len(cells)
    pct = 100.0 * found / total if total else 0.0
    status = "COMPLETE" if not missing else f"{len(missing)} MISSING"
    print(f"\n[{label}]  {found}/{total}  ({pct:5.1f}%)  {status}")
    print(f"    dir: {out_dir}")
    if missing:
        ids = ",".join(map(str, sorted(set(missing))))
        print(f"    -> RESUBMIT: sbatch --array={ids}%4 {submit}")
    else:
        print(f"    -> TRANSFER:")
        for ln in xfer_lines:
            print(f"        {ln}")
    return found, total


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=DEFAULT_BASE)
    ap.add_argument("--remote", default="ec474@login-cpu.hpc.cam.ac.uk")
    ap.add_argument("--local-dest", default="/d/ADNI_BIDS_project/derivatives",
                    help="Local derivatives/ root. CSD3 ADNI_MRI/<tree> -> derivatives/<tree>.")
    args = ap.parse_args()
    b, host, dest = args.base, args.remote, args.local_dest

    print("=" * 78)
    print(f"  Per-task sweep audit (resubmit if incomplete / transfer if complete)   base={b}")
    print("=" * 78)

    specs = []  # (label, out_dir, cells, submit_cmd, xfer_lines)

    # T3a-d cached, 3 archs (one block per horizon)
    for arch, rel in CACHED.items():
        buckets = _cells_by_task(T3_TASKS, 4)
        submit = f"{SUB}/04_{arch}_head_sweep_T3abcd_submit_csd3.sh"
        for task in T3_TASKS:
            specs.append((f"{task} / {arch}", os.path.join(b, rel), buckets[task], submit,
                          _xfer(b, host, dest, f"{rel}/{task}")))
    # T4 cached + T1e cached, brainmvp + braindino
    for task, tag in [("T4_conv_horizon", "T4"), ("T1e_pcn_vs_scn", "T1e")]:
        for arch in ("brainmvp", "braindino"):
            rel = CACHED[arch]
            submit = f"{SUB}/04_{arch}_head_sweep_{tag}_submit_csd3.sh"
            cells = _cells_by_task([task], 1)[task]
            specs.append((f"{task} / {arch}", os.path.join(b, rel), cells, submit,
                          _xfer(b, host, dest, f"{rel}/{task}")))
    # BrainMVP T4 full_ft per AUGMENT (CSD3 single-nest; LOCAL double-nest brainmvp_debug/brainmvp_debug/)
    for aug in ("none", "stochastic", "plus_original"):
        cells = [(s, f"BrainMVP_uniformer/T4_conv_horizon/seed_{s}/full_ft/metrics.json") for s in SEEDS]
        src_rel = f"brainmvp_debug/aug_{aug}/BrainMVP_uniformer/T4_conv_horizon"
        dst_rel = f"brainmvp_debug/brainmvp_debug/aug_{aug}/BrainMVP_uniformer/T4_conv_horizon"
        submit = f"--export=ALL,AUGMENT={aug} mri_pipeline/brain_mvp/04i_finetune_BrainMVP_T4_full_ft_submit_csd3.sh"
        specs.append((f"T4_conv_horizon full_ft / aug={aug}",
                      os.path.join(b, f"brainmvp_debug/aug_{aug}"), cells, submit,
                      _xfer(b, host, dest, src_rel, dest_rel=dst_rel)))

    gf = gt = 0
    for label, out_dir, cells, submit, xfer in specs:
        f, t = _audit(label, out_dir, cells, submit, xfer)
        gf += f; gt += t

    print("\n" + "=" * 78)
    print(f"  TOTAL: {gf}/{gt} cells complete ({100.0 * gf / gt if gt else 0:.1f}%)")
    print("  Per task above: RESUBMIT line if incomplete, TRANSFER (rsync/scp) line if complete.")
    print("=" * 78)


if __name__ == "__main__":
    main()
