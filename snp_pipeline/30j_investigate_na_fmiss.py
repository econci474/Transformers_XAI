r"""
30j_investigate_na_fmiss.py   (LOCAL, env: snp)
==============================================
For the recoverable SNPs whose `filter_reason` in 30c reports F_MISS=NA
(reporting-limitation of the QC-stage probe-name keying), run PLINK
--missing on the unfiltered 812-subject BED to get their TRUE per-SNP
F_MISS. Emit a verdict on whether they're recoverable at --geno 0.05.

Outputs → `GWAS_comprehensive_v2/LD_report/`:
  fmiss_na_investigation.tsv  one row per recoverable rsID with
                              F_MISS_true + verdict
  fmiss_na_investigation.txt  human summary

Usage:
  conda run -n snp python snp_pipeline/30j_investigate_na_fmiss.py
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
UNFILT_BED = BASE / "WGS_Omni25_BIN_wo_ConsentsIssues"
PLINK = BASE / "plink.exe"
DELTA_RECOV = BASE / "GWAS_comprehensive_v2" / "delta_recoverable.tsv"
RECOV_ALL = (BASE / "source_prs" / "unfiltered_SNP_reconciliation"
              / "recoverable_snps.tsv")
OUT_DIR = BASE / "GWAS_comprehensive_v2" / "LD_report"


def _verdict(f_miss: float) -> str:
    if f_miss < 0.02:
        return "should_have_passed_geno_0.02_unclear_why_filtered"
    if f_miss < 0.05:
        return "recoverable_at_geno_0.05"
    if f_miss < 0.10:
        return "recoverable_at_geno_0.10_borderline"
    return "not_recoverable_too_sparse"


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    # Read both delta_recoverable.tsv (v2+v3 NEW) AND the full recoverable
    # list, then dedupe by rsID — investigate ALL F_MISS=NA cases.
    recov_all = pd.read_csv(RECOV_ALL, sep="\t", dtype=str,
                              keep_default_na=False)
    na_rows = recov_all[
        recov_all["filter_reason"].str.contains("F_MISS=NA",
                                                  na=False, regex=False)].copy()
    print(f"[30j] F_MISS=NA rows in recoverable_snps.tsv: {len(na_rows)}")
    target_rsids = sorted(set(na_rows["rsID_pub"]))
    print(f"[30j] unique rsIDs to investigate: {len(target_rsids)}")
    for rs in target_rsids:
        srcs = na_rows.loc[na_rows["rsID_pub"] == rs, "source"].tolist()
        print(f"  {rs}  ←  sources: {','.join(sorted(set(srcs)))}")

    if not target_rsids:
        print("[30j] No F_MISS=NA rsIDs found.")
        return

    # Run PLINK --missing on the unfiltered BED for these rsIDs
    tmp = Path(tempfile.mkdtemp(prefix="fmiss_inv_"))
    extract_path = tmp / "rsids.txt"
    extract_path.write_text("\n".join(target_rsids) + "\n", "utf-8")
    out_prefix = tmp / "missing"
    cmd = [str(PLINK),
           "--bfile", str(UNFILT_BED),
           "--extract", str(extract_path),
           "--missing",
           "--out", str(out_prefix)]
    print(f"[30j] {' '.join(cmd)}")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("[30j] PLINK STDERR:", r.stderr[-2000:])
        sys.exit(r.returncode)

    lmiss = out_prefix.with_suffix(".lmiss")
    if not lmiss.exists():
        sys.exit(f"[30j] expected {lmiss}, got nothing")

    miss = pd.read_csv(lmiss, sep=r"\s+", dtype=str)
    # PLINK .lmiss columns: CHR SNP N_MISS N_GENO F_MISS
    miss["F_MISS"] = pd.to_numeric(miss["F_MISS"], errors="coerce")
    print(f"[30j] PLINK returned {len(miss)} rsIDs in .lmiss")

    # Per-rsID resolution row (merge with source membership)
    src_map = (na_rows.groupby("rsID_pub")["source"]
                 .apply(lambda s: ",".join(sorted(set(s)))).to_dict())
    fr_map = (na_rows.groupby("rsID_pub")["filter_reason"]
                .first().to_dict())

    rows = []
    for _, r in miss.iterrows():
        rs = r["SNP"]
        fm = r["F_MISS"]
        rows.append({
            "rsID": rs,
            "sources": src_map.get(rs, ""),
            "old_filter_reason": fr_map.get(rs, ""),
            "F_MISS_true_812subj": round(float(fm), 5) if pd.notna(fm) else "NA",
            "N_miss": r.get("N_MISS", ""),
            "N_geno": r.get("N_GENO", ""),
            "verdict": _verdict(float(fm)) if pd.notna(fm) else "no_data",
        })

    out_df = pd.DataFrame(rows).sort_values("F_MISS_true_812subj")
    out_path = args.out / "fmiss_na_investigation.tsv"
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"  [out] {out_path}")
    print()
    print(out_df.to_string(index=False))

    # Text summary
    lines = [
        "F_MISS=NA recoverable SNPs — true missingness investigation",
        f"Source data: WGS_Omni25_BIN_wo_ConsentsIssues "
        f"(812 subjects, pre-QC). PLINK --missing.",
        f"Reason for NA in 30c: QC reports keyed by Illumina kgp probe "
        f"names; these rsIDs sit directly in the unfiltered BIM under "
        f"their rs* identifier so the rsID-keyed lookup misses them.",
        "",
        out_df.to_string(index=False),
        "",
        "Verdict thresholds:",
        "  < 0.02   should_have_passed_geno_0.02_unclear_why_filtered",
        "  < 0.05   recoverable_at_geno_0.05",
        "  < 0.10   recoverable_at_geno_0.10_borderline",
        "  >= 0.10  not_recoverable_too_sparse",
    ]
    (args.out / "fmiss_na_investigation.txt").write_text(
        "\n".join(lines) + "\n", "utf-8")
    print(f"  [out] {args.out/'fmiss_na_investigation.txt'}")


if __name__ == "__main__":
    main()
