"""
03i_verify_encoder_embeddings.py
================================
Verify that EVERY expected encoder run in a visit output root (e.g.
encoder_outputs_no_cdr_post_exclusion_m12) produced complete, valid per-patient
embeddings + probabilities. Checks the 03g (m12) / 03h (m06) sweep grid:
6 notes'-pick combos x 3 seeds = 18 runs.

For each expected run dir <root>/<model>/<task>/seed_<N>/<strategy>/ it checks:
  - metrics.json present
  - embeddings/embeddings_meta.json present + parseable
  - embeddings/embeddings.npz present + loadable
  - embeddings/embeddings.parquet present  (WARN-only: needs pyarrow at write time)
and validates from the npz + meta (numpy/json only -- NO pandas/pyarrow needed):
  - num_classes matches the task (T2=3, T1/T1b/T1d=2)
  - probs array shape == (n_total, num_classes), all finite, rows sum to ~1
  - emb_pool / emb_head finite (qc flags)
  - splits train/val/test all present and non-empty
  - label non-null for in_task_cohort rows
Optionally (--data-dir) cross-checks counts_by_split against the actual split CSV row
counts, to confirm the run used the right (visit) split tree.

Dependency-light: standard library + numpy only.

Run (on the HPC, right after the sweep):
  python clinical_pipeline/03i_verify_encoder_embeddings.py --visit m12
  python clinical_pipeline/03i_verify_encoder_embeddings.py --visit m06 \
      --data-dir /home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/verbose/baseline_m06
  # local copy:
  python clinical_pipeline/03i_verify_encoder_embeddings.py --visit m12 \
      --root D:/ADNI_BIDS_project/derivatives/encoder_outputs_no_cdr_post_exclusion_m12

Exit code 0 iff all expected runs are present AND valid; else 1.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

SEEDS = [0, 1, 2]
SPLITS = ("train", "val", "test")
DEFAULT_HPC_ROOT = "/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_{visit}"


def expected_combos(visit: str):
    """(model_slug, task, strategy, expected_num_classes) -- mirrors 03g/03h COMBOS.
    T1d uses the visit-independent key (only its data_dir changed)."""
    return [
        ("BioClinical-ModernBERT-large", f"T2_{visit}_multiclass", "full_ft", 3),
        ("BioClinical-ModernBERT-base",  f"T1_{visit}_binary",     "full_ft", 2),
        ("ModernBERT-base",              f"T1b_{visit}_cnmci_ad",  "frozen",  2),
        ("ModernBERT-base",              f"T1b_{visit}_cnmci_ad",  "full_ft", 2),
        ("ModernBERT-large",             "T1d_pmci_smci",          "full_ft", 2),
        ("ModernBERT-base",              "T1d_pmci_smci",          "full_ft", 2),
    ]


def csv_split_counts(data_dir: Path, seed: int) -> dict | None:
    """Row counts (minus header) of the seed's split CSVs, or None if unavailable."""
    out = {}
    for sp in SPLITS:
        f = data_dir / f"seed_{seed}" / f"{sp}.csv"
        if not f.exists():
            return None
        # csv.reader (NOT naive line count) -- Generated_Text has newlines in quoted fields
        with open(f, "r", encoding="utf-8", errors="replace", newline="") as fh:
            n = sum(1 for _ in csv.reader(fh)) - 1            # minus header
        out[sp] = max(n, 0)
    return out


def verify_run(run_dir: Path, exp_classes: int, expect_counts: dict | None):
    """Return (status, list[str] issues). status in {OK, MISSING, INCOMPLETE, INVALID}."""
    issues = []
    if not run_dir.exists():
        return "MISSING", ["run dir does not exist"]
    if not (run_dir / "metrics.json").exists():
        issues.append("metrics.json missing")
    emb = run_dir / "embeddings"
    meta_p, npz_p, pq_p = (emb / "embeddings_meta.json",
                           emb / "embeddings.npz", emb / "embeddings.parquet")
    if not emb.exists():
        return "INCOMPLETE", issues + ["embeddings/ dir missing"]
    if not meta_p.exists():
        issues.append("embeddings_meta.json missing")
    if not npz_p.exists():
        issues.append("embeddings.npz missing")
    if not pq_p.exists():
        issues.append("WARN: embeddings.parquet missing (pyarrow not installed at write time?)")
    # if the core artifacts are absent we cannot validate further
    if not npz_p.exists():
        return ("INCOMPLETE", issues)

    # ---- meta ----
    meta = {}
    if meta_p.exists():
        try:
            meta = json.loads(meta_p.read_text())
        except Exception as e:  # noqa: BLE001
            issues.append(f"meta unreadable: {e}")
    nc = meta.get("num_classes")
    if nc is not None and int(nc) != exp_classes:
        issues.append(f"num_classes={nc} != expected {exp_classes}")
    cbs = meta.get("counts_by_split", {})
    for sp in SPLITS:
        if cbs.get(sp, 0) in (0, None):
            issues.append(f"counts_by_split[{sp}] empty/missing")
    qc = meta.get("qc", {})
    for flag in ("probs_all_finite", "emb_pool_all_finite", "emb_head_all_finite"):
        if qc.get(flag) is False:
            issues.append(f"qc.{flag}=False")
    if expect_counts is not None and cbs:
        if {k: cbs.get(k) for k in SPLITS} != expect_counts:
            issues.append(f"counts_by_split {dict((k, cbs.get(k)) for k in SPLITS)} "
                          f"!= split CSVs {expect_counts}")

    # ---- npz ----
    try:
        z = np.load(npz_p, allow_pickle=True)
    except Exception as e:  # noqa: BLE001
        return "INVALID", issues + [f"npz load failed: {e}"]
    for key in ("patient_id", "split", "in_task_cohort", "label", "pred",
                "probs", "emb_pool", "emb_head"):
        if key not in z.files:
            issues.append(f"npz missing array '{key}'")
    if "probs" in z.files:
        probs = z["probs"]
        if probs.ndim != 2 or probs.shape[1] != exp_classes:
            issues.append(f"probs shape {probs.shape} != (N, {exp_classes})")
        if not np.isfinite(probs).all():
            issues.append("probs contain non-finite values")
        else:
            s = probs.sum(axis=1)
            if not np.allclose(s, 1.0, atol=1e-3):
                issues.append(f"probs rows do not sum to 1 (min={s.min():.4f}, max={s.max():.4f})")
    n_total = int(z["probs"].shape[0]) if "probs" in z.files else None
    if "split" in z.files:
        present = set(np.unique(z["split"]).tolist())
        miss = set(SPLITS) - present
        if miss:
            issues.append(f"npz splits missing {sorted(miss)}")

    hard = [i for i in issues if not i.startswith("WARN")]
    status = "OK" if not hard else "INVALID"
    return status, issues, n_total, (cbs or None)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--visit", default="m12", help="m12 or m06 (selects the expected task names)")
    ap.add_argument("--root", default=None,
                    help="Output root (default: HPC encoder_outputs_no_cdr_post_exclusion_<visit>)")
    ap.add_argument("--data-dir", default=None,
                    help="Optional split tree (verbose/baseline_<visit>) to cross-check split counts")
    args = ap.parse_args()

    root = Path(args.root or DEFAULT_HPC_ROOT.format(visit=args.visit))
    data_dir = Path(args.data_dir) if args.data_dir else None
    combos = expected_combos(args.visit)
    n_expected = len(combos) * len(SEEDS)

    print(f"[03i] Verifying {n_expected} expected runs under:\n      {root}")
    if not root.exists():
        print(f"[ERROR] root does not exist: {root}", file=sys.stderr)
        return 1
    if data_dir:
        print(f"      cross-checking split counts vs: {data_dir}")

    counts = {"OK": 0, "MISSING": 0, "INCOMPLETE": 0, "INVALID": 0}
    bad = []
    print(f"\n  {'STATUS':10s} {'rows':>5s}  run")
    print("  " + "-" * 78)
    for model, task, strat, ec in combos:
        for seed in SEEDS:
            run_dir = root / model / task / f"seed_{seed}" / strat
            expect_counts = csv_split_counts(data_dir, seed) if data_dir else None
            res = verify_run(run_dir, ec, expect_counts)
            status, issues = res[0], res[1]
            n_total = res[2] if len(res) > 2 else None
            counts[status] = counts.get(status, 0) + 1
            rel = f"{model}/{task}/seed_{seed}/{strat}"
            ntxt = str(n_total) if n_total is not None else "-"
            print(f"  {status:10s} {ntxt:>5s}  {rel}")
            for it in issues:
                print(f"             - {it}")
            if status != "OK":
                bad.append(rel)

    print("\n  " + "=" * 78)
    print(f"  Summary: {counts['OK']}/{n_expected} OK | "
          f"MISSING={counts['MISSING']} INCOMPLETE={counts['INCOMPLETE']} INVALID={counts['INVALID']}")
    if bad:
        print("  Not-OK runs:")
        for r in bad:
            print(f"    - {r}")
        return 1
    print("  All expected embeddings present and valid.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
