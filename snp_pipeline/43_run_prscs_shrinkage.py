"""
43_run_prscs_shrinkage.py
=========================
AAO Phase 2b — PRS-CS Bayesian β shrinkage with 1000G EUR LD reference.

For each PRS source with full sumstats available:
  1. Harmonise sumstats to PRS-CS format: SNP A1 A2 BETA SE (or BETA P).
  2. Run PRScs.py with the target BIM (strict-QC pool) and EUR LD ref.
  3. Concatenate per-chromosome posterior β outputs into a single TSV.

Sources WITHOUT full sumstats fall back to raw β only (PRS-CS skipped).

Output:  D:/ADNI_SNP_Omni2.5M_20140220/source_prs/prscs_posterior/
  <SOURCE>_prscs_posterior_beta.tsv          (rsID, A1, A2, posterior_beta)
  <SOURCE>_prscs_harmonised_sumstats.tsv     (sumstats input we fed PRS-CS)
  <SOURCE>_prscs_run.log                     (PRS-CS console log)
  prscs_summary.tsv                          (per-source status + N_matched)

PRS-CS settings (defaults):
  phi: auto (Bayesian; ok with ≥100k GWAS samples)
  iterations: 1000  (500 burnin, thin 5) — default
  ref_dir: D:/ADNI_SNP_Omni2.5M_20140220/PRScs_ref/ldblk_1kg_eur
  bim_prefix: D:/ADNI_SNP_Omni2.5M_20140220/GWAS_comprehensive_v2/QC_strict/recover_all_pool_strictQC

Usage (env: snp):
  python snp_pipeline/43_run_prscs_shrinkage.py                 # all eligible sources
  python snp_pipeline/43_run_prscs_shrinkage.py --source Bellenguez
  python snp_pipeline/43_run_prscs_shrinkage.py --dry-run        # harmonise + smoke-test, skip PRS-CS
  python snp_pipeline/43_run_prscs_shrinkage.py --chrom 22       # one chrom only (faster smoke test)
"""
from __future__ import annotations
import argparse
import gzip
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
import pandas as pd
import numpy as np

# ── Paths ─────────────────────────────────────────────────────────────────
ROOT = Path("D:/ADNI_SNP_Omni2.5M_20140220")
PRSCS_PY = ROOT / "PRScs_tool" / "PRScs.py"
LD_REF   = ROOT / "PRScs_ref" / "ldblk_1kg_eur"
TARGET_BIM_PREFIX = (ROOT / "GWAS_comprehensive_v2" / "QC_strict"
                      / "recover_all_pool_strictQC")
OUT_DIR  = ROOT / "source_prs" / "prscs_posterior"
TMP_DIR  = OUT_DIR / "_tmp"

# ── Per-source sumstats config ────────────────────────────────────────────
# Each entry says: where the file is, how to map its columns to
# {SNP, A1, A2, BETA, SE, P}, and the GWAS sample size to pass.
# n_gwas chosen as the total effective sample size from each paper's abstract.
SOURCES = {
    "Bellenguez": {
        "path": ROOT / "GWAS/Bellenguez_2022/GCST90027158_buildGRCh38.tsv.gz",
        "sep": "\t",
        "snp":  "variant_id",
        "a1":   "effect_allele",
        "a2":   "other_allele",
        "beta": "beta",
        "se":   "standard_error",
        "p":    "p_value",
        "n_gwas": 487511,                     # 39106 cases + 401577 controls + 46828 proxy_cases (from abstract)
        "build": "GRCh38",
    },
    "Wightman": {
        # Use the user's already-lifted Wightman without-UKB file (3,121 rows,
        # has rsID + BETA + SE + p resolved against ADNI hg19 BIM).
        "path": ROOT / "GWAS/Wightman_2021/wightman_without_ukb_GRCh38.tsv",
        "sep": "\t",
        "snp":  "rsID",
        "a1":   "effect_allele",
        "a2":   "other_allele",
        "beta": "beta",
        "se":   "standard_error",
        "p":    "p_value",
        "n_gwas": 472868,
        "build": "GRCh37 (lifted to GRCh38 in this file)",
    },
    "Kunkle": {
        "path": ROOT / "GWAS/Kunkle_2019/NG00075_Kunkle_etal_Stage1_P-val_only_results.txt",
        "sep":  r"\s+",
        "snp":  "MarkerName",
        "a1":   "Effect_allele",
        "a2":   "Non_Effect_allele",
        "beta": "Beta",
        "se":   "SE",
        "p":    "Pvalue",
        "n_gwas": 63926,                      # 21,982 cases + 41,944 controls (Stage 1)
        "build": "GRCh37",
    },
    "Kosteridis": {
        "path": ROOT / "GWAS/Kosteridis_2024/GCST90449060.tsv.gz",
        "sep": "\t",
        "snp":  "rs_id",
        "a1":   "effect_allele",
        "a2":   "other_allele",
        "beta": "beta",
        "se":   "standard_error",
        "p":    "p_value",
        "n_gwas": 436000,                     # MTAG paper: 436k effective
        "build": "GRCh37",
    },
    "Schwanzentruber": {
        # Schwanzentruber has both MTAG `beta` and marginal `GWAS_BETA`. For
        # PRS-CS we want the marginal (un-MTAG-boosted) beta.
        "path": ROOT / "GWAS/Schwanzentruber_2021/GCST90012877_buildGRCh37.tsv.gz",
        "sep": "\t",
        "snp":  "variant_id",
        "a1":   "effect_allele",
        "a2":   "other_allele",
        "beta": "GWAS_BETA",
        "se":   "GWAS_SE",
        "p":    "GWAS_P",
        "n_gwas": 446092,                     # paper Table 1: 446,092 effective
        "build": "GRCh37",
    },
    "DeRojas": {
        # SPIGAPUK2 meta (Spanish + IGAP + UKB Phase 2). Extracted from
        # Sumstats_SPIGAPUK2_20190625.zip on 2026-05-27.
        "path": ROOT / "GWAS/DeRojas_2021/Sumstats_SPIGAPUK2_20190625.txt",
        "sep":  "\t",
        "snp":  "RS",
        "a1":   "A1",
        "a2":   "A2",
        "beta": "Beta",
        "se":   "SE",
        "p":    "P",
        "n_gwas": 88000,                       # 21,982 SP + 25,580 IGAP + ~40k UKB Phase 2
        "build": "GRCh37",
    },
}

# Sources without full sumstats — skipped from PRS-CS; raw β retained downstream.
# Lambert: OpenGWAS ieu-a-298 only exposes top-hits at p<0.01, not full sumstats.
SKIPPED_SOURCES = [
    "Desikan", "Ebenau", "Felsky_IT", "Felsky_MF", "Huang",
    "Lambert", "Leonenko", "Najar", "ONeil_GHR", "ONeil_NPY", "Vesilievick",
    "Zhang", "Kosteridis_novel_AD", "Kosteridis_shared_AD_CV",
]


# ── Sumstats harmonisation ────────────────────────────────────────────────
def _read_sumstats(cfg: dict) -> pd.DataFrame:
    """Read raw sumstats and return a DataFrame with at least columns:
      SNP A1 A2 BETA SE P (some may be NaN for sources where we derive them)."""
    p = cfg["path"]
    print(f"  reading {p} ({p.stat().st_size // (1024*1024)} MB)...", flush=True)
    if str(p).endswith(".gz"):
        df = pd.read_csv(p, sep=cfg["sep"], low_memory=False, compression="gzip")
    else:
        df = pd.read_csv(p, sep=cfg["sep"], low_memory=False)
    cols = {}
    if cfg.get("snp"):
        cols["SNP"] = df[cfg["snp"]]
    cols["A1"] = df[cfg["a1"]]
    cols["A2"] = df[cfg["a2"]]
    if "beta" in cfg:
        cols["BETA"] = pd.to_numeric(df[cfg["beta"]], errors="coerce")
    if "se" in cfg:
        cols["SE"]   = pd.to_numeric(df[cfg["se"]],   errors="coerce")
    if "p" in cfg:
        cols["P"]    = pd.to_numeric(df[cfg["p"]],    errors="coerce")
    # Wightman z-only path
    if "z" in cfg:
        cols["Z"] = pd.to_numeric(df[cfg["z"]], errors="coerce")
        if "eaf" in cfg:
            cols["EAF"] = pd.to_numeric(df[cfg["eaf"]], errors="coerce")
        if "n" in cfg:
            cols["N"] = pd.to_numeric(df[cfg["n"]], errors="coerce")
        # chr:pos for rsID lookup
        if "chromosome" in df.columns: cols["chrom"] = df["chromosome"]
        if "base_pair_location" in df.columns: cols["bp"] = df["base_pair_location"]
    out = pd.DataFrame(cols)
    print(f"    raw rows: {len(out):,}", flush=True)
    return out


def _derive_beta_se_from_z(df: pd.DataFrame) -> pd.DataFrame:
    """For Wightman: BETA = Z / sqrt(2*N*EAF*(1-EAF) + Z²), SE = 1 / sqrt(...).
    Standard binary-trait log-OR scale derivation."""
    z, eaf, n = df["Z"], df["EAF"], df["N"]
    denom = np.sqrt(2.0 * n * eaf * (1.0 - eaf) + z**2)
    df["SE"]   = 1.0 / denom
    df["BETA"] = z / denom
    return df


def _attach_rsid_via_chrpos(df: pd.DataFrame, target_bim: pd.DataFrame) -> pd.DataFrame:
    """Wightman has no rsID column. Match rsID via chrom:bp using the target BIM
    (GRCh38 positions). Wightman is GRCh37 → matches require liftover.
    For 115 target SNPs, simplest is to look them up in dbSNP, but
    we cheat: only the rsIDs in target_bim matter, so we attempt
    chr:pos join in BOTH builds via the target BIM's CHR:BP."""
    # Target BIM has GRCh38 positions. Wightman has GRCh37. We can join via
    # GRCh37 positions only if we have a GRCh37 version of the target BIM.
    # For now, just skip Wightman if no rsID column — flag and warn.
    print("    [WARN] Wightman has no rsID; skipping chr:pos rsID join (GRCh37/38 mismatch)")
    df["SNP"] = "rs_unknown"
    return df


def _format_for_prscs(df: pd.DataFrame, use_se: bool) -> pd.DataFrame:
    """Return df with columns SNP A1 A2 BETA SE (or BETA P), dropping rows
    with missing required values."""
    if use_se:
        out = df[["SNP", "A1", "A2", "BETA", "SE"]].dropna()
    else:
        out = df[["SNP", "A1", "A2", "BETA", "P"]].dropna()
    # uppercase alleles, drop non-ATGC
    for c in ("A1", "A2"):
        out = out[out[c].astype(str).str.upper().isin({"A","T","G","C"})]
        out[c] = out[c].astype(str).str.upper()
    # drop duplicated rsIDs (PRS-CS requires unique)
    out = out.drop_duplicates(subset=["SNP"], keep="first")
    return out.reset_index(drop=True)


def _run_prscs(src_label: str, sst_path: Path, n_gwas: int, out_prefix: Path,
                chrom: int | None = None) -> dict:
    """Invoke PRS-CS subprocess. Returns dict with per-chrom output paths."""
    if not LD_REF.exists():
        raise SystemExit(f"LD ref not found at {LD_REF}. Did download/extract finish?")
    cmd = [
        sys.executable, str(PRSCS_PY),
        f"--ref_dir={LD_REF}",
        f"--bim_prefix={TARGET_BIM_PREFIX}",
        f"--sst_file={sst_path}",
        f"--n_gwas={n_gwas}",
        f"--out_dir={out_prefix}",
        "--seed=42",
    ]
    if chrom is not None:
        cmd.append(f"--chrom={chrom}")
    log_path = out_prefix.with_name(out_prefix.name + "_run.log")
    print(f"  >>> {' '.join(cmd)}", flush=True)
    with open(log_path, "w") as lf:
        t0 = time.time()
        r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, text=True)
        dt = time.time() - t0
    if r.returncode != 0:
        print(f"    PRS-CS exit {r.returncode} (see {log_path})")
        return {"ok": False, "log": log_path}
    print(f"    PRS-CS done in {dt:.1f}s", flush=True)
    return {"ok": True, "log": log_path}


def _aggregate_per_chrom(out_prefix: Path) -> pd.DataFrame:
    """Concatenate _pst_eff_*.txt files across chromosomes into one DataFrame."""
    parent = out_prefix.parent
    pattern = out_prefix.name + "_pst_eff_a1_b0.5_phiauto_chr"
    rows = []
    for f in sorted(parent.glob(pattern + "*.txt")):
        df = pd.read_csv(f, sep=r"\s+", header=None,
                          names=["chrom","rsID","bp","A1","A2","posterior_beta"])
        rows.append(df)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", help="Run only this source.")
    ap.add_argument("--dry-run", action="store_true",
                    help="Harmonise sumstats + write to disk; skip PRS-CS subprocess.")
    ap.add_argument("--chrom", type=int, help="Run only this chromosome (smoke test).")
    args = ap.parse_args()

    if not PRSCS_PY.exists():
        raise SystemExit(f"PRScs.py not found at {PRSCS_PY}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    # Target BIM rsIDs (used to filter sumstats for size + accept-rate diagnostics)
    bim = pd.read_csv(TARGET_BIM_PREFIX.with_suffix(".bim"), sep=r"\s+", header=None,
                       names=["chrom","rsID","cM","bp","A1","A2"])
    target_rsids = set(bim["rsID"])
    print(f"Target BIM: {TARGET_BIM_PREFIX}.bim  ({len(target_rsids)} SNPs)")

    sources_to_run = [args.source] if args.source else list(SOURCES.keys())

    summary = []
    for src in sources_to_run:
        if src not in SOURCES:
            print(f"\n*** {src} not in SOURCES registry — skipping ***")
            continue
        print(f"\n========== {src} ==========")
        cfg = SOURCES[src]
        try:
            df = _read_sumstats(cfg)
        except Exception as e:
            print(f"  failed to read sumstats: {e}")
            summary.append({"source": src, "status": "read_failed", "error": str(e)})
            continue

        # Wightman: derive BETA/SE from Z then attempt rsID join
        if "Z" in df.columns:
            df = _derive_beta_se_from_z(df)
            if cfg.get("needs_rsid_join"):
                df = _attach_rsid_via_chrpos(df, bim)

        # PRS-CS prefers BETA + SE; fall back to BETA + P
        use_se = "SE" in df.columns and df["SE"].notna().any()
        sst = _format_for_prscs(df, use_se=use_se)
        n_intersect = sst["SNP"].isin(target_rsids).sum()
        print(f"  harmonised sumstats: {len(sst):,} SNPs  "
              f"(intersect target: {n_intersect}/{len(target_rsids)})")

        if len(sst) == 0:
            summary.append({"source": src, "status": "empty_sumstats", "n_intersect": 0})
            continue

        sst_path = TMP_DIR / f"{src}_prscs_sst.txt"
        sst.to_csv(sst_path, sep="\t", index=False)
        print(f"  wrote {sst_path}")

        if args.dry_run:
            summary.append({"source": src, "status": "dry_run",
                             "n_sst": len(sst), "n_intersect": int(n_intersect)})
            continue

        out_prefix = OUT_DIR / src
        result = _run_prscs(src, sst_path, cfg["n_gwas"], out_prefix, chrom=args.chrom)
        if not result["ok"]:
            summary.append({"source": src, "status": "prscs_failed",
                             "log": str(result["log"])})
            continue
        posterior = _aggregate_per_chrom(out_prefix)
        post_path = OUT_DIR / f"{src}_prscs_posterior_beta.tsv"
        posterior.to_csv(post_path, sep="\t", index=False)
        summary.append({"source": src, "status": "ok",
                         "n_posteriors": len(posterior),
                         "n_intersect": int(n_intersect),
                         "out": str(post_path)})
        print(f"  posterior β: {len(posterior)} SNPs → {post_path}")

    sum_df = pd.DataFrame(summary)
    sum_p = OUT_DIR / "prscs_summary.tsv"
    sum_df.to_csv(sum_p, sep="\t", index=False)
    print(f"\nwrote summary: {sum_p}")
    print(sum_df.to_string(index=False))

    # Skipped sources note
    if SKIPPED_SOURCES:
        print(f"\nSkipped (no full sumstats; raw β retained downstream):")
        for s in SKIPPED_SOURCES:
            print(f"  - {s}")


if __name__ == "__main__":
    main()
