"""
49_export_meta_prs_en_combined.py
=================================
Materialise the per-subject **meta-PRS-EN-combined** score so downstream analyses
(e.g. the clinical PH-diagnostics in clinical_pipeline/02i) can read it as a column.

LEAKAGE-SAFE BY DESIGN — fit PER SEED on the TRAIN fold only, never on the whole
cohort. A full-cohort fit would let each subject's own AD-conversion outcome inform
their score, which is test-set leakage the moment those scores enter any evaluation.
This mirrors exactly the meta-EN Cox fit already used in 46_prs_cox_strictQC.py
(CoxPHFitter(penalizer=0.1, l1_ratio=0.5) over the per-source PRS table). The exported
score is `predict_log_partial_hazard` (a relative log-hazard / risk score) for the
TRAIN subjects of that seed — i.e. in-sample for the fold the EN was fit on, which is
the canonical setting for a PH-assumption (Schoenfeld) check of that fitted model.

Reuses 46's helpers via importlib (digit-prefixed module name) so the meta-EN logic
stays identical and cannot drift:
  - _load_time_event(seed, split)        -> t (years), e (AD_final==1)
  - _fit_en_cox_meta_prs(sp_df, parts, seed) -> (cph, log-partial-hazard Series by PTID)
  - the meta source-column selection (drop dedup variants + Kosteridis subgroups)
and per_source_prs_table() from _prs_strict_qc_lib for the PRS source matrix.

Inputs
  D:\\ADNI_BIDS_project\\derivatives\\clinical\\no_cdr_stratified_post_exclusion
      \\tabular\\baseline\\seed_{0,1,2}\\{train,val,test}.csv   (t,e; via 46)
  LD-pruned strict-QC per-source PRS table (built in memory by per_source_prs_table,
  ld_config = ld_1000kb_r2_0.8).

Outputs (TSV: Patient_ID <tab> meta_prs_EN_combined)
  D:\\ADNI_SNP_Omni2.5M_20140220\\source_prs\\meta_prs_en_combined_seed{0,1,2}_train.tsv

Env: the same env that runs 45/46/47 (needs lifelines + the strict-QC PRS stack).
Run:
  python snp_pipeline/49_export_meta_prs_en_combined.py
  python snp_pipeline/49_export_meta_prs_en_combined.py --ld-config ld_1000kb_r2_0.8
"""
from __future__ import annotations
import argparse
import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from _prs_strict_qc_lib import per_source_prs_table, DEFAULT_LD_CONFIG  # noqa: E402

# ── Import 46's helpers (module name starts with a digit -> importlib) ───────────
_spec = importlib.util.spec_from_file_location(
    "prs_cox_strictQC", HERE / "46_prs_cox_strictQC.py")
_cox = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_cox)
_load_time_event   = _cox._load_time_event
_fit_en_cox_meta_prs = _cox._fit_en_cox_meta_prs
SEEDS = _cox.SEEDS

OUT_DIR = Path("D:/ADNI_SNP_Omni2.5M_20140220/source_prs")

# Meta-EN source columns — identical rule to 46._run_one (meta_prs_EN_combined):
# every PRS_* source except the dedup variants and the Kosteridis SUBGROUPS
# (the full Kosteridis = MTAG_AD ∪ shared_AD_CV, so its subgroups would triple-count).
_META_DROP = ("PRS_prs_all_dedup", "PRS_prs_all_dedup_ivw",
              "PRS_prs_all_dedup_filtered", "PRS_Kosteridis_MTAG_AD",
              "PRS_Kosteridis_shared_AD_CV", "PRS_Kosteridis_novel_AD")


def meta_source_frame(prs_full: pd.DataFrame) -> pd.DataFrame:
    """PTID-indexed frame of the meta-EN source PRS columns (mirrors 46)."""
    source_cols = [c for c in prs_full.columns
                   if c.startswith("PRS_") and c not in _META_DROP]
    sp_df = prs_full[["PTID"] + source_cols].set_index("PTID")
    sp_df.index = sp_df.index.astype(str)
    return sp_df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ld-config", default=DEFAULT_LD_CONFIG)
    ap.add_argument("--beta-source", default="raw", choices=["raw", "prscs"])
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Building per-source PRS table (ld={args.ld_config}, beta={args.beta_source})...")
    prs_full, _ = per_source_prs_table(ld_config=args.ld_config,
                                       beta_source=args.beta_source)
    sp_df = meta_source_frame(prs_full)
    print(f"  meta-EN sources: {sp_df.shape[1]} columns, {sp_df.shape[0]} subjects with PRS")

    for seed in SEEDS:
        parts = {sp: _load_time_event(seed, sp) for sp in ("train", "val", "test")}
        cph, en_prs = _fit_en_cox_meta_prs(sp_df, parts, seed)
        if cph is None or en_prs.isna().all():
            print(f"  [seed {seed}] meta-EN fit skipped/failed — no output.")
            continue
        # Export TRAIN (in-sample for the EN fit) AND val/test (OUT-OF-SAMPLE
        # predictions from the SAME train-fitted EN — leakage-safe: the model never
        # saw val/test outcomes, so applying it to those subjects is a legitimate
        # downstream feature, exactly as in 45/46's held-out evaluation).
        for split in ("train", "val", "test"):
            split_ids = parts[split]["Patient_ID"].astype(str).tolist()
            scored = en_prs.reindex([p for p in split_ids if p in en_prs.index]).dropna()
            out = pd.DataFrame({"Patient_ID": scored.index,
                                "meta_prs_EN_combined": scored.values})
            out_path = OUT_DIR / f"meta_prs_en_combined_seed{seed}_{split}.tsv"
            out.to_csv(out_path, sep="\t", index=False)
            n_ev = int(parts[split].set_index("Patient_ID")
                       .reindex(scored.index)["e"].sum())
            print(f"  [seed {seed}] {split:5s}: wrote {len(out):3d} subjects "
                  f"(events={n_ev}, score range [{scored.min():.3f}, {scored.max():.3f}]) "
                  f"-> {out_path.name}")

    print("\nDone.")


if __name__ == "__main__":
    main()
