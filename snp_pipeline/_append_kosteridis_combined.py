"""
_append_kosteridis_combined.py
==============================
Compute the new unified Kosteridis source's classification + Cox + AAO
metrics across the 4 LD configs and append them to the existing leaderboards
without re-running other sources.

Imports the _run_one functions from 45/46/47 directly (via importlib because
of the digit-leading filenames) so we re-use the exact same fit logic and
metric schema as the existing tables.

Updates in place:
  outputs/strict_qc_prs/ld_<cfg>/{classification,cox,aao_regression}/
      per_run_metrics.tsv   ← append/replace Kosteridis rows
      leaderboard.csv       ← re-aggregate

Then the master leaderboard at outputs/strict_qc_prs/master_leaderboard.{tsv,png}
should be re-rendered separately via script 48.
"""
from __future__ import annotations
import importlib.util as _il
import sys
import warnings
from pathlib import Path
import pandas as pd

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _prs_strict_qc_lib import (per_source_prs_table, load_subject_covariates,
                                  get_dedup_dosage_matrix)  # noqa: E402

OUT_BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220/outputs/strict_qc_prs")
LD_CONFIGS = ["ld_1000kb_r2_0.8", "ld_500kb_r2_0.2",
               "ld_250kb_r2_0.1", "ld_50_5_r2_0.5"]
SEEDS = (0, 1, 2)
SRC = "Kosteridis"


def _load_task_module(filename: str):
    """importlib loader for the digit-leading task scripts."""
    p = HERE / filename
    spec = _il.spec_from_file_location(filename.replace(".py","").replace("-","_"), p)
    mod  = _il.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def update_for_task(task_subdir: str, task_mod, agg_fn):
    print(f"\n=== task: {task_subdir} ===")
    for cfg in LD_CONFIGS:
        out_dir = OUT_BASE / cfg / task_subdir
        per_run_p = out_dir / "per_run_metrics.tsv"
        lb_p = out_dir / "leaderboard.csv"
        if not per_run_p.exists() or not lb_p.exists():
            print(f"  [skip] {out_dir} (missing per_run/leaderboard)"); continue
        prs_full, _ = per_source_prs_table(ld_config=cfg, beta_source="raw")
        cov = load_subject_covariates()
        dedup_st, dedup_dosage = get_dedup_dosage_matrix(cfg)
        # Determine n_snps for the Kosteridis column
        prs_col = f"PRS_{SRC}"
        n_snps_kost = int(prs_full[prs_col].notna().any())  # rough flag, recomputed below
        # Run all (mode × seed) for Kosteridis
        new_rows = []
        for mode, covs in task_mod.COVAR_MODES.items():
            for seed in SEEDS:
                r = task_mod._run_one(SRC, prs_full, cov, seed, mode, covs,
                                        dedup_dosage=dedup_dosage)
                if r is None: continue
                new_rows.append(r)
        if not new_rows:
            print(f"  [no rows] {cfg}"); continue
        # Set n_snps_used per row using the merged SNP table size
        # (need to re-query, since some scripts compute this differently).
        # We do a quick post-hoc: count the SNPs in the Kosteridis source table.
        from _prs_strict_qc_lib import load_source_beta, load_ldpruned_dosage
        bim, _ = load_ldpruned_dosage(cfg)
        st = load_source_beta(SRC, bim, beta_source="raw")
        n_snps_kost = len(st)
        for r in new_rows:
            r["n_snps_used"] = n_snps_kost

        # Merge into existing per_run + drop any prior Kosteridis rows.
        existing = pd.read_csv(per_run_p, sep="\t")
        existing = existing[existing["source"] != SRC].copy()
        combined = pd.concat([existing, pd.DataFrame(new_rows)], ignore_index=True)
        combined.to_csv(per_run_p, sep="\t", index=False)

        agg = agg_fn(combined)
        agg.to_csv(lb_p, index=False)
        print(f"  OK {cfg}: +{len(new_rows)} Kosteridis rows  "
              f"(n_snps_used={n_snps_kost}); leaderboard -> {len(agg)} rows")


# Aggregators (mirror what scripts 45/46/47 emit) ----------------------------
def _agg_classification(df: pd.DataFrame) -> pd.DataFrame:
    metric_cols = ["val_auc","val_bacc","val_f1","val_sens","val_spec",
                    "val_prec","val_pr_auc","or_per_sd",
                    "test_auc","test_bacc","test_f1"]
    return df.groupby(["source","covar_mode"]).agg(
        n_snps=("n_snps_used","first"), n_seeds=("seed","count"),
        val_n=("val_n","mean"), val_pos=("val_pos","mean"), val_neg=("val_neg","mean"),
        val_tn=("val_tn","mean"), val_fp=("val_fp","mean"),
        val_fn=("val_fn","mean"), val_tp=("val_tp","mean"),
        **{f"{c}_mean": (c, "mean") for c in metric_cols},
        **{f"{c}_std":  (c, "std")  for c in metric_cols},
    ).reset_index().sort_values(["covar_mode","val_bacc_mean"], ascending=[True, False])


def _agg_cox(df: pd.DataFrame) -> pd.DataFrame:
    return df.groupby(["source","covar_mode"]).agg(
        n_snps=("n_snps_used","first"), n_seeds=("seed","count"),
        val_cindex_mean=("val_cindex","mean"), val_cindex_std=("val_cindex","std"),
        hr_per_sd_mean=("hr_per_sd","mean"), hr_per_sd_std=("hr_per_sd","std"),
        hr_lo_mean=("hr_lo","mean"), hr_hi_mean=("hr_hi","mean"),
        hr_p_min=("hr_p","min"),
        test_cindex_mean=("test_cindex","mean"),
    ).reset_index().sort_values(["covar_mode","val_cindex_mean"], ascending=[True, False])


def _agg_aao(df: pd.DataFrame) -> pd.DataFrame:
    return df.groupby(["source","covar_mode"]).agg(
        n_snps=("n_snps_used","first"), n_seeds=("seed","count"),
        val_r2_mean=("val_r2","mean"), val_r2_std=("val_r2","std"),
        val_rmse_mean=("val_rmse","mean"),
        val_pearson_r_mean=("val_pearson_r","mean"),
        val_pearson_p_min=("val_pearson_p","min"),
        yrs_per_sd_mean=("yrs_per_sd_prs","mean"), yrs_per_sd_std=("yrs_per_sd_prs","std"),
    ).reset_index().sort_values(["covar_mode","val_r2_mean"], ascending=[True, False])


def main():
    mod_cls = _load_task_module("45_prs_classification_strictQC.py")
    mod_cox = _load_task_module("46_prs_cox_strictQC.py")
    mod_aao = _load_task_module("47_prs_aao_regression_strictQC.py")

    update_for_task("classification", mod_cls, _agg_classification)
    update_for_task("cox",            mod_cox, _agg_cox)
    update_for_task("aao_regression", mod_aao, _agg_aao)


if __name__ == "__main__":
    main()
