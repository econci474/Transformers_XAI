"""
42_strict_qc_recover_all_pool.py
================================
AAO Phase 2a — strict-QC re-do.

Apply tighter PLINK filters to the GRCh38 / 616-subject / pre-MAF cohort
and re-extract the 128 recover_all_pool target rsIDs from that strictly-
filtered superset.

Why not just filter the existing recover_all_pool.bed?  Because the existing
pool already passed loose --geno 0.02 (with several bypass-exception rescues
for high-missingness published PRS leads). To apply a stricter --geno 0.01
we have to back up to a BED that still contains the variants we'd be cutting.

Upstream pre-strict-QC BED (chosen 2026-05-27):
  D:/ADNI_SNP_Omni2.5M_20140220/SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr
  (2,248,376 variants, 616 subjects — has had --autosome --mind 0.02
   --geno 0.02 --hwe 1e-7, kgp->rsID, dedup palindromic, deprecated rsID
   patch, GRCh37->GRCh38 liftover, REF-correction; NO --maf yet.)

Strict thresholds (user-confirmed 2026-05-27):
  --geno 0.01     (variant missingness ≤ 1%; stricter than the existing 0.02)
  --mind 0.02     (sample missingness ≤ 2%; same as existing)
  --maf  0.01     (MAF ≥ 1%; same threshold but now applied with strict geno)
  --hwe  1e-6     (HWE p > 1e-6; stricter than existing 1e-7)

Output folder (sibling to recover_all_pool/, per `feedback-aao-qc-strict-folder`):
  D:/ADNI_SNP_Omni2.5M_20140220/GWAS_comprehensive_v2/QC_strict/

Outputs:
  full_strict_qc.bed/bim/fam              — full cohort post-strict-QC
  full_strict_qc.log                       — PLINK log
  recover_all_pool_strictQC.bed/bim/fam   — 128 targets ∩ post-strict-QC
  recover_all_pool_strictQC_dosage.tsv     — per-patient dosage matrix
  pre_strict_diag.frq                      — pre-filter MAF
  pre_strict_diag.lmiss                    — pre-filter variant missingness
  pre_strict_diag.hwe                      — pre-filter HWE
  strict_qc_report.tsv                     — per-target-SNP retention + drop reason
  strict_qc_summary.txt                    — terminal-style summary

Usage (env: snp):
  python snp_pipeline/42_strict_qc_recover_all_pool.py
  python snp_pipeline/42_strict_qc_recover_all_pool.py --dry-run
"""
from __future__ import annotations
import argparse
import shutil
import subprocess
import time
from pathlib import Path
import pandas as pd
import numpy as np

# ── Paths ─────────────────────────────────────────────────────────────────
ROOT = Path("D:/ADNI_SNP_Omni2.5M_20140220")
PLINK = ROOT / "plink.exe"
UPSTREAM = ROOT / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr"
RECOVER  = ROOT / "GWAS_comprehensive_v2/recover_all_pool/SNP_recover_all_GRCh38_PRS_only"
OUT_DIR  = ROOT / "GWAS_comprehensive_v2/QC_strict"

# ── Strict thresholds (do NOT change without user explicit approval) ──────
THR = {
    "geno": 0.01,
    "mind": 0.02,
    "maf":  0.01,
    "hwe":  1e-6,
}


def _run(cmd: list[str], log: str = ""):
    print(f"\n>>> {' '.join(cmd)}")
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    if r.returncode != 0:
        print("STDOUT:\n", r.stdout[-1500:])
        print("STDERR:\n", r.stderr[-1500:])
        raise SystemExit(f"PLINK failed (exit {r.returncode}) for: {log or cmd[0]}")
    print(f"    done in {dt:.1f}s")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    # 0. Sanity
    for ext in (".bed", ".bim", ".fam"):
        p = UPSTREAM.with_suffix(ext)
        if not p.exists():
            raise SystemExit(f"Missing upstream file: {p}")
    if not RECOVER.with_suffix(".bim").exists():
        raise SystemExit(f"Missing recover_all_pool BIM at {RECOVER}.bim")
    if not PLINK.exists():
        raise SystemExit(f"plink.exe not found at {PLINK}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Read target rsIDs from recover_all_pool.bim
    pool_bim = pd.read_csv(RECOVER.with_suffix(".bim"), sep=r"\s+",
                            header=None, names=["chrom","rsID","cM","bp","A1","A2"])
    target_rsids = pool_bim["rsID"].tolist()
    target_set = set(target_rsids)
    print(f"recover_all_pool target rsIDs: {len(target_rsids)}")

    extract_file = OUT_DIR / "target_rsids.txt"
    extract_file.write_text("\n".join(target_rsids) + "\n")

    if args.dry_run:
        print("[DRY] Would now run PLINK strict-QC + diagnostic passes.")
        return

    # 2. Diagnostic pass on the PRE-strict-QC BED — gives us per-target-SNP
    #    raw MAF / missingness / HWE so we can attribute drop reasons.
    diag = OUT_DIR / "pre_strict_diag"
    _run([str(PLINK),
           "--bfile", str(UPSTREAM),
           "--extract", str(extract_file),
           "--freq", "--missing", "--hardy",
           "--out", str(diag)],
          log="diagnostic stats for targets")

    # 3. Strict-QC pass on the FULL upstream cohort
    full = OUT_DIR / "full_strict_qc"
    _run([str(PLINK),
           "--bfile", str(UPSTREAM),
           "--geno", str(THR["geno"]),
           "--mind", str(THR["mind"]),
           "--maf",  str(THR["maf"]),
           "--hwe",  str(THR["hwe"]),
           "--make-bed",
           "--out", str(full)],
          log="full strict QC pass")

    # 4. Extract the 128 targets from the strict-QC cohort
    target_strict = OUT_DIR / "recover_all_pool_strictQC"
    _run([str(PLINK),
           "--bfile", str(full),
           "--extract", str(extract_file),
           "--make-bed",
           "--out", str(target_strict)],
          log="extract 128 targets from strict-QC cohort")

    # 5. Per-patient dosage matrix
    _run([str(PLINK),
           "--bfile", str(target_strict),
           "--recode", "A",
           "--out", str(target_strict.with_name(target_strict.name + "_dosage"))],
          log="dosage matrix")

    # 5b. Reshape .raw → tidy tsv (PTID × rsID)
    raw_p = target_strict.with_name(target_strict.name + "_dosage.raw")
    raw_df = pd.read_csv(raw_p, sep=r"\s+")
    snp_cols = [c for c in raw_df.columns if c not in
                  ("FID","IID","PAT","MAT","SEX","PHENOTYPE")]
    # PLINK suffixes _A1 to rsID columns — strip it for clean rsID headers
    rename = {c: c.split("_")[0] for c in snp_cols}
    tidy = raw_df[["IID"] + snp_cols].rename(columns={"IID": "PTID", **rename})
    tidy_p = target_strict.with_name(target_strict.name + "_dosage.tsv")
    tidy.to_csv(tidy_p, sep="\t", index=False)
    print(f"\nwrote dosage tsv: {tidy_p}  ({tidy.shape})")

    # 6. QC report — per target rsID: kept? if dropped, which filter?
    strict_bim = pd.read_csv(target_strict.with_suffix(".bim"), sep=r"\s+",
                              header=None, names=["chrom","rsID","cM","bp","A1","A2"])
    kept_in_strict = set(strict_bim["rsID"])

    frq = pd.read_csv(f"{diag}.frq", sep=r"\s+")            # CHR SNP A1 A2 MAF NCHROBS
    lmiss = pd.read_csv(f"{diag}.lmiss", sep=r"\s+")        # CHR SNP N_MISS N_GENO F_MISS
    hwe = pd.read_csv(f"{diag}.hwe", sep=r"\s+")            # CHR SNP TEST A1 A2 GENO O(HET) E(HET) P
    hwe = hwe[hwe["TEST"] == "ALL(NP)"] if "ALL(NP)" in hwe["TEST"].unique() else \
          hwe[hwe["TEST"].str.startswith("ALL")]

    frq_d   = frq.set_index("SNP")["MAF"].to_dict()
    lmiss_d = lmiss.set_index("SNP")["F_MISS"].to_dict()
    hwe_d   = hwe.set_index("SNP")["P"].to_dict()
    bim_full = pd.read_csv(UPSTREAM.with_suffix(".bim"), sep=r"\s+",
                            header=None, names=["chrom","rsID","cM","bp","A1","A2"])
    in_upstream = set(bim_full["rsID"])

    rows = []
    for rs in target_rsids:
        in_up = rs in in_upstream
        maf = frq_d.get(rs, np.nan)
        fmiss = lmiss_d.get(rs, np.nan)
        hwep  = hwe_d.get(rs, np.nan)
        kept = rs in kept_in_strict
        reason = ""
        if not in_up:
            reason = "not_in_upstream"
        elif kept:
            reason = "kept"
        else:
            # Diagnose primary drop reason from the diagnostic stats
            fail = []
            if not np.isnan(fmiss) and fmiss > THR["geno"]:
                fail.append(f"F_MISS={fmiss:.4f}")
            if not np.isnan(maf) and maf < THR["maf"]:
                fail.append(f"MAF={maf:.4f}")
            if not np.isnan(hwep) and hwep < THR["hwe"]:
                fail.append(f"HWE_p={hwep:.2e}")
            reason = " AND ".join(fail) if fail else "dropped_in_pipeline_(unknown_cause)"
        rows.append({"rsID": rs,
                       "in_upstream_pre_strict": in_up,
                       "kept_in_strict_QC": kept,
                       "pre_strict_MAF": maf,
                       "pre_strict_F_MISS": fmiss,
                       "pre_strict_HWE_p": hwep,
                       "drop_reason": reason})
    report = pd.DataFrame(rows)
    rep_p = OUT_DIR / "strict_qc_report.tsv"
    report.to_csv(rep_p, sep="\t", index=False)
    print(f"\nwrote report: {rep_p}")

    # 7. Summary text
    n_targets = len(target_rsids)
    n_kept    = int(report["kept_in_strict_QC"].sum())
    n_dropped = n_targets - n_kept
    breakdown = (report.loc[~report["kept_in_strict_QC"], "drop_reason"]
                       .value_counts())
    sum_p = OUT_DIR / "strict_qc_summary.txt"
    full_bim = pd.read_csv(full.with_suffix(".bim"), sep=r"\s+",
                            header=None, names=["chrom","rsID","cM","bp","A1","A2"])
    full_fam = pd.read_csv(full.with_suffix(".fam"), sep=r"\s+",
                            header=None,
                            names=["FID","IID","PID","MID","SEX","PHENO"])
    text = []
    text.append(f"Strict-QC summary — built {pd.Timestamp.now():%Y-%m-%d %H:%M}")
    text.append(f"Upstream BED: {UPSTREAM}")
    text.append(f"  variants: {len(bim_full):,}  subjects: 616 (MRI keep-cohort)")
    text.append(f"Strict thresholds: --geno {THR['geno']} --mind {THR['mind']} "
                f"--maf {THR['maf']} --hwe {THR['hwe']:.0e}")
    text.append(f"Post-strict full cohort: {len(full_bim):,} variants, {len(full_fam)} subjects")
    text.append("")
    text.append(f"recover_all_pool targets:  {n_targets}")
    text.append(f"  kept in strict QC:       {n_kept}")
    text.append(f"  dropped in strict QC:    {n_dropped}")
    text.append("")
    text.append("Drop reason breakdown:")
    for r, c in breakdown.items():
        text.append(f"  {c:>3d}  {r}")
    text.append("")
    text.append(f"Output BED:    {target_strict}.bed/.bim/.fam")
    text.append(f"Dosage TSV:    {tidy_p}")
    text.append(f"Report:        {rep_p}")
    sum_p.write_text("\n".join(text) + "\n")
    print("\n" + "\n".join(text))


if __name__ == "__main__":
    main()
