"""
parse_gwas_results.py
=====================
Post-GWAS significance report for the ADNI CN vs AD/LMCI logistic regression.

Reads the PLINK .assoc.logistic output and applies four thresholds:

  1. Genome-wide significance   p  < 5e-8     (conventional, no correction)
  2. Suggestive significance    p  < 1e-5     (conventional, no correction)
  3. Bonferroni correction      p  < 0.05 / N_tested
  4. Benjamini-Hochberg FDR     q  < 0.05    (standard GWAS FDR; controls expected
                                              fraction of false positives among hits)

Outputs
-------
  results/gwas/gwas_significance_report.txt   human-readable summary
  results/gwas/gwas_significant_snps.tsv      all thresholds per SNP
  results/gwas/gwas_bmfm_labels.tsv           for BMFM fine-tuning
                                               (label=1: GW-significant;
                                                label=0: p > 0.5, null SNPs)

Usage
-----
  python parse_gwas_results.py
  python parse_gwas_results.py --assoc results/gwas/gwas_CN_vs_AD.assoc.logistic
  python parse_gwas_results.py --fdr-threshold 0.05   # relax FDR
"""

import argparse
import io
import pathlib
import sys
from datetime import datetime

import numpy as np
import pandas as pd

# Force UTF-8 on Windows
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ── CLI ───────────────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser(
    description="Post-GWAS significance analysis with multiple correction methods.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
_p.add_argument("--assoc", default="results/gwas/gwas_CN_vs_AD.assoc.logistic",
                help="PLINK .assoc.logistic file")
_p.add_argument("--out-dir", default="results/gwas", dest="out_dir",
                help="Output directory (default: results/gwas)")
_p.add_argument("--alpha", type=float, default=0.05,
                help="Family-wise error rate for Bonferroni (default: 0.05)")
_p.add_argument("--fdr-threshold", type=float, default=0.05, dest="fdr_threshold",
                help="BH-FDR q-value threshold (default: 0.05 — standard GWAS FDR)")
_p.add_argument("--null-p-min", type=float, default=0.5, dest="null_p_min",
                help="Minimum p-value for a SNP to be labelled 'null' (label=0) "
                     "for BMFM fine-tuning (default: 0.5)")
args = _p.parse_args()

ASSOC_FILE    = pathlib.Path(args.assoc)
OUT_DIR       = pathlib.Path(args.out_dir)
REPORT_FILE   = OUT_DIR / "gwas_significance_report.txt"
SIGSNPS_FILE  = OUT_DIR / "gwas_significant_snps.tsv"
BMFM_FILE     = OUT_DIR / "gwas_bmfm_labels.tsv"

OUT_DIR.mkdir(parents=True, exist_ok=True)

if not ASSOC_FILE.exists():
    print(f"[ERROR] ASSOC file not found: {ASSOC_FILE}")
    sys.exit(1)

# ── Load PLINK assoc.logistic ─────────────────────────────────────────────────
# Columns: CHR  SNP  BP  A1  TEST  NMISS  OR  SE  L95  U95  STAT  P
print(f"Loading: {ASSOC_FILE}")
raw = pd.read_csv(ASSOC_FILE, sep=r"\s+")
raw.columns = raw.columns.str.strip()

# Keep only the ADD (allelic association) rows; drop covariate rows
df = raw[raw["TEST"] == "ADD"].copy()
df["P"] = pd.to_numeric(df["P"], errors="coerce")
df = df.dropna(subset=["P"])

N = len(df)
print(f"  {N:,} SNPs with valid p-values (ADD test)")

# ── Thresholds ────────────────────────────────────────────────────────────────
GW_THRESH        = 5e-8
SUGGESTIVE_THRESH = 1e-5
BONF_THRESH      = args.alpha / N
FDR_THRESH       = args.fdr_threshold

# Additional BH-FDR thresholds to always report
FDR_EXTRA = [0.20, 0.10, 0.05, 0.01]

print(f"\n  Genome-wide threshold : {GW_THRESH:.0e}")
print(f"  Suggestive threshold  : {SUGGESTIVE_THRESH:.0e}")
print(f"  Bonferroni threshold  : {BONF_THRESH:.3e}  (0.05 / {N:,})")
print(f"  BH-FDR threshold      : {FDR_THRESH:.0e}")

# ── Benjamini-Hochberg FDR ────────────────────────────────────────────────────
try:
    from statsmodels.stats.multitest import multipletests
    _reject, q_values, _, _ = multipletests(df["P"].values, alpha=0.05, method="fdr_bh")
    df["q_BH"] = q_values
    print("  BH-FDR computed via statsmodels.")
except ImportError:
    # Manual BH implementation
    print("  statsmodels not found — computing BH-FDR manually.")
    p_sorted_idx = np.argsort(df["P"].values)
    p_sorted     = df["P"].values[p_sorted_idx]
    m            = len(p_sorted)
    q_manual     = np.minimum.accumulate((p_sorted * m / np.arange(1, m+1))[::-1])[::-1]
    q_back       = np.empty(m)
    q_back[p_sorted_idx] = q_manual
    df["q_BH"] = np.minimum(q_back, 1.0)

# ── Flag each SNP under each threshold ────────────────────────────────────────
df["sig_GW"]        = df["P"] < GW_THRESH
df["sig_suggestive"] = df["P"] < SUGGESTIVE_THRESH
df["sig_bonferroni"] = df["P"] < BONF_THRESH
df["sig_FDR"]       = df["q_BH"] < FDR_THRESH

def get_subset(flag_col):
    cols = ["SNP", "CHR", "BP", "A1", "OR", "L95", "U95", "STAT", "P", "q_BH"]
    return (df[df[flag_col]][cols]
            .sort_values("P")
            .reset_index(drop=True))

gw_hits        = get_subset("sig_GW")
sugg_hits      = get_subset("sig_suggestive")
bonf_hits      = get_subset("sig_bonferroni")
fdr_hits       = get_subset("sig_FDR")

# Extra FDR levels
fdr_extra_counts = {q: int((df["q_BH"] < q).sum()) for q in FDR_EXTRA}

# ── Write significance report ─────────────────────────────────────────────────
SEP  = "═" * 72
SEP2 = "─" * 72

def fmt_snp_table(subset: pd.DataFrame) -> list[str]:
    if subset.empty:
        return ["  (none)"]
    lines = []
    hdr = f"  {'rsID':<20} {'CHR':>4} {'BP':>12} {'A1':>3} {'OR':>8} {'95%CI':>20} {'P-value':>12} {'q_BH':>12}"
    lines.append(hdr)
    lines.append("  " + "-" * 70)
    for _, r in subset.iterrows():
        ci = f"[{r['L95']:.3f}, {r['U95']:.3f}]"
        lines.append(
            f"  {str(r['SNP']):<20} {int(r['CHR']):>4} {int(r['BP']):>12,} "
            f"{str(r['A1']):>3} {r['OR']:>8.3f} {ci:>20} {r['P']:>12.3e} {r['q_BH']:>12.3e}"
        )
    return lines

report_lines = [
    SEP,
    "  GWAS Significance Report — ADNI CN vs AD/LMCI",
    f"  Generated : {datetime.now().strftime('%Y-%m-%d %H:%M')}",
    f"  Input     : {ASSOC_FILE}",
    f"  SNPs tested: {N:,}",
    SEP,
    "",
    "  CORRECTION THRESHOLDS APPLIED",
    SEP2,
    f"  1. Genome-wide significance  : p < {GW_THRESH:.0e}  (conventional, no correction)",
    f"  2. Suggestive significance   : p < {SUGGESTIVE_THRESH:.0e}  (conventional, no correction)",
    f"  3. Bonferroni correction     : p < {BONF_THRESH:.3e}  (alpha=0.05 / {N:,} tests)",
    f"  4. Benjamini-Hochberg FDR    : q < {FDR_THRESH:.0e}",
    "",
]

# Section 1: GW
report_lines += [
    SEP2,
    f"  [1] GENOME-WIDE SIGNIFICANT  (p < {GW_THRESH:.0e})   —  {len(gw_hits)} SNP(s)",
    SEP2,
] + fmt_snp_table(gw_hits) + [""]

# Section 2: Suggestive
report_lines += [
    SEP2,
    f"  [2] SUGGESTIVE               (p < {SUGGESTIVE_THRESH:.0e})   —  {len(sugg_hits)} SNP(s)",
    SEP2,
] + fmt_snp_table(sugg_hits) + [""]

# Section 3: Bonferroni
report_lines += [
    SEP2,
    f"  [3] BONFERRONI               (p < {BONF_THRESH:.3e})  —  {len(bonf_hits)} SNP(s)",
    f"      (alpha = {args.alpha} / {N:,} tests)",
    SEP2,
] + fmt_snp_table(bonf_hits) + [""]

# Section 4: BH-FDR
report_lines += [
    SEP2,
    f"  [4] BENJAMINI-HOCHBERG FDR   (q < {FDR_THRESH:.0e})  —  {len(fdr_hits)} SNP(s)",
    SEP2,
] + fmt_snp_table(fdr_hits) + [""]

# BH-FDR at extra thresholds (informational)
report_lines += [
    "  BH-FDR significant counts at other q thresholds (informational):",
]
for q, cnt in sorted(fdr_extra_counts.items(), reverse=True):
    report_lines.append(f"    q < {q:.0e}  :  {cnt} SNP(s)")
report_lines += [""]

# Summary table
report_lines += [
    SEP,
    "  SUMMARY",
    SEP2,
    f"  {'Correction method':<40} {'Threshold':>12}  {'Hits':>6}",
    f"  {'-'*40} {'----------':>12}  {'----':>6}",
    f"  {'Genome-wide significance':<40} {GW_THRESH:>12.0e}  {len(gw_hits):>6}",
    f"  {'Suggestive significance':<40} {SUGGESTIVE_THRESH:>12.0e}  {len(sugg_hits):>6}",
    f"  {'Bonferroni (FWER < 0.05)':<40} {BONF_THRESH:>12.3e}  {len(bonf_hits):>6}",
    f"  {'Benjamini-Hochberg FDR':<40} {FDR_THRESH:>12.0e}  {len(fdr_hits):>6}",
    SEP,
    "",
    "  BMFM LABELLING NOTE",
    SEP2,
    f"  label=1 : genome-wide significant SNPs (p < {GW_THRESH:.0e})  →  {len(gw_hits)} SNPs",
    f"  label=0 : null SNPs (p >= {args.null_p_min})  →  {int((df['P'] >= args.null_p_min).sum()):,} SNPs",
    f"  (All intermediate-p SNPs are excluded from fine-tuning labels)",
    f"  Output  : {BMFM_FILE}",
    SEP,
]

report_text = "\n".join(report_lines)
print("\n" + report_text)

REPORT_FILE.write_text(report_text, encoding="utf-8")
print(f"\nReport written: {REPORT_FILE}")

# ── Write significant SNPs TSV ────────────────────────────────────────────────
# All SNPs significant under at least one threshold
any_sig = df["sig_GW"] | df["sig_suggestive"] | df["sig_bonferroni"] | df["sig_FDR"]
sig_all = df[any_sig].copy()
sig_all["thresholds_passed"] = (
    sig_all["sig_GW"].map({True: "GW;", False: ""}) +
    sig_all["sig_suggestive"].map({True: "suggestive;", False: ""}) +
    sig_all["sig_bonferroni"].map({True: "bonferroni;", False: ""}) +
    sig_all["sig_FDR"].map({True: f"FDR_{FDR_THRESH:.0e};", False: ""})
).str.rstrip(";")

out_cols = ["SNP", "CHR", "BP", "A1", "OR", "L95", "U95", "STAT", "P", "q_BH",
            "sig_GW", "sig_suggestive", "sig_bonferroni", "sig_FDR", "thresholds_passed"]
sig_all[out_cols].sort_values("P").to_csv(SIGSNPS_FILE, sep="\t", index=False)
print(f"Significant SNPs TSV: {SIGSNPS_FILE}  ({len(sig_all)} rows)")

# ── Write BMFM labels TSV ─────────────────────────────────────────────────────
# label=1 : genome-wide significant
# label=0 : p >= null_p_min (clearly non-associated)
# Excluded: 1e-5 < p < null_p_min (ambiguous)
gw_label  = df[df["sig_GW"]].copy()
gw_label["label"] = 1
gw_label["label_basis"] = f"GW_sig_p<{GW_THRESH:.0e}"

null_label = df[df["P"] >= args.null_p_min].copy()
null_label["label"] = 0
null_label["label_basis"] = f"null_p>={args.null_p_min}"

bmfm = pd.concat([gw_label, null_label], ignore_index=True)
bmfm_cols = ["SNP", "CHR", "BP", "A1", "P", "q_BH", "label", "label_basis"]
bmfm[bmfm_cols].to_csv(BMFM_FILE, sep="\t", index=False)

print(f"BMFM labels TSV     : {BMFM_FILE}")
print(f"  label=1 (case)    : {int((bmfm['label']==1).sum()):,} SNPs (GW-significant)")
print(f"  label=0 (null)    : {int((bmfm['label']==0).sum()):,} SNPs (p >= {args.null_p_min})")
print(f"  excluded          : {N - len(bmfm):,} SNPs (ambiguous p-value range)")
