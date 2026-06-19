r"""
30k_recover_all_pool.py   (LOCAL, env: snp)
==========================================
Build the **recover-all** PRS BED by extending the on-MAF pool (116 SNPs
across 22 PRS sources) with all 12 recoverable rsIDs (non-monomorphic,
on-chip-but-filtered) from `recoverable_snps.tsv` + the F_MISS=NA cases
resolved by 30j_investigate_na_fmiss.

User-directed (2026-05-25 night): recover ALL non-monomorphic SNPs from
the 22-source PRS pool. Result is a single PRS-ready BED:

  Recovery classes:
    standard (--geno 0.05 passes):
      rs117618017  F_MISS=0.022 (Bellenguez+Wightman+others)
      rs34674752   F_MISS=0.023 (DeRojas+Ebenau+Vesilievick)
      rs2405442    F_MISS=0.027 (Yang_Mic+Yang_Oli, PTK2B)
      rs7401792    F_MISS=0.034 (Bellenguez)
      rs7617515    F_MISS=0.042 (Najar, ZNF852)
    relaxed (--geno 0.10 needed):
      rs17057043   F_MISS=0.060 (Zhang, PTK2B)
      rs1859788    F_MISS=0.074 (DeRojas+Ebenau+Schwanzentruber+Vesilievick, PILRA)
      rs744373     F_MISS=0.083 (Zhang, BIN1)
      rs12325539   F_MISS=0.089 (Yang_Ex+Yang_In)
    bypass exceptions (no QC; published-PRS lead retention):
      rs4266886    F_MISS=0.120 (Desikan+Leonenko, CR1)
      rs114105899  MAF=0.0080  (Felsky_IT)
      rs115186657  MAF=0.0037  (Wightman)
  Monomorphic — NOT recoverable (no within-cohort variation):
      rs429358     APOE ε4   (Schwanzentruber+Wightman+Zhang)
      rs75932628   TREM2 R47H (Bellenguez+DeRojas+Ebenau+Vesilievick)

Two-pass PLINK from the 812-subject unfiltered BED:
  Pass A: 9 standard rescues with `--geno 0.10 --hwe 1e-7 --maf 0.001`
  Pass B: 3 bypass exceptions extracted without --geno/--maf/--hwe
Merge with the existing 22-source MAF BIM (already GRCh38), then output.

Outputs → `D:\…\GWAS_comprehensive_v2\recover_all_pool\`:
  recovery_report.tsv                     per-rsID recovery class + source + true F_MISS/MAF
  recovery_summary.txt                    human summary
  SNP_recover_all_GRCh38_PRS_only.{bed,bim,fam}   final BED (116 + 12 = 128 SNPs × 616 patients)

Usage:
  conda run -n snp python snp_pipeline/30k_recover_all_pool.py
"""
from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
UNFILT = BASE / "WGS_Omni25_BIN_wo_ConsentsIssues"
KEEP_FAM = BASE / "SNP_filtered_with_mri.fam"
MAF_BIM_STEM = BASE / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr_maf001"
CHAIN = BASE / "liftover" / "hg19ToHg38.over.chain.gz"
PLINK = BASE / "plink.exe"
RECON_DIR = BASE / "source_prs" / "unfiltered_SNP_reconciliation"
FMISS_INV = (BASE / "GWAS_comprehensive_v2" / "LD_report"
              / "fmiss_na_investigation.tsv")
OUT_DEFAULT = BASE / "GWAS_comprehensive_v2" / "recover_all_pool"

# Recovery classification — locked from 30j + recoverable_snps.tsv investigation
STANDARD_RESCUE = {
    "rs117618017": ("F_MISS=0.0222", "Bellenguez+Wightman+others", ""),
    "rs34674752":  ("F_MISS=0.0234", "DeRojas+Ebenau+Vesilievick", "SHARPIN"),
    "rs2405442":   ("F_MISS=0.0271", "Yang_Mic+Yang_Oli", "PTK2B region"),
    "rs7401792":   ("F_MISS=0.0345", "Bellenguez", ""),
    "rs7617515":   ("F_MISS=0.0419", "Najar", "ZNF852"),
}
RELAXED_RESCUE = {
    "rs17057043":  ("F_MISS=0.0603", "Zhang", "PTK2B"),
    "rs1859788":   ("F_MISS=0.0739", "DeRojas+Ebenau+Schwanzentruber+Vesilievick", "PILRA"),
    "rs744373":    ("F_MISS=0.0825", "Zhang", "BIN1"),
    "rs12325539":  ("F_MISS=0.0887", "Yang_Ex+Yang_In", ""),
}
BYPASS_EXCEPTIONS = {
    "rs4266886":   ("F_MISS=0.1195", "Desikan+Leonenko", "CR1 (published PRS lead)"),
    "rs114105899": ("MAF=0.0080",    "Felsky_IT",        "(published PRS lead)"),
    "rs115186657": ("MAF=0.0037",    "Wightman",         "(published PRS lead)"),
}
MONOMORPHIC_LOST = {
    "rs429358":    ("MAF=0.0000",    "Schwanzentruber+Wightman+Zhang", "APOE ε4 — no variation"),
    "rs75932628":  ("MAF=0.0000",    "Bellenguez+DeRojas+Ebenau+Vesilievick", "TREM2 R47H — no variation"),
}


def _run(cmd: list[str], desc: str) -> str:
    print(f"  [plink] {desc}")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"    STDERR (tail): {r.stderr[-1500:]}")
        raise RuntimeError(f"plink failed for: {desc}")
    return r.stdout


def _build_unfilt_name_map() -> dict[str, str]:
    """Return {rsID_pub: unfilt_name} for all 12 recoverable rsIDs by
    looking up `unfilt_name` from recoverable_snps.tsv. Some are direct
    rsIDs (rs2405442, rs1859788, etc.); others are Illumina kgp probe
    IDs (rs4266886→kgp5825939) that need rename after extract."""
    rec_path = (BASE / "source_prs" / "unfiltered_SNP_reconciliation"
                  / "recoverable_snps.tsv")
    rec = pd.read_csv(rec_path, sep="\t", dtype=str, keep_default_na=False)
    target = set(STANDARD_RESCUE) | set(RELAXED_RESCUE) | set(BYPASS_EXCEPTIONS)
    out = {}
    for rs in target:
        rows = rec[rec["rsID_pub"] == rs]
        if rows.empty:
            print(f"  [warn] {rs} not in recoverable_snps.tsv")
            out[rs] = rs   # fallback to direct rsID
            continue
        out[rs] = str(rows.iloc[0]["unfilt_name"])
    return out


def _load_on_maf_pool() -> set[str]:
    """116 unique on-MAF rsIDs across 22 PRS sources."""
    PRS = ["Bellenguez","Wightman","Kunkle","Desikan","Lambert","DeRojas",
            "Schwanzentruber","Najar","Ebenau","Leonenko","Vesilievick",
            "Zhang","Felsky_MF","Felsky_IT","ONeil_NPY","ONeil_GHR",
            "Kosteridis_novel_AD","Kosteridis_shared_AD_CV",
            "Yang_Ex","Yang_In","Yang_Ast","Yang_Mic","Yang_Oli","Yang_Opc"]
    pool: set[str] = set()
    for src in PRS:
        p = RECON_DIR / f"{src}_unfiltered_reconciliation.tsv"
        if not p.exists():
            continue
        df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
        df["in_maf_bim"] = df["in_maf_bim"].str.lower().isin({"true","1"})
        pool |= set(df.loc[df["in_maf_bim"], "rsID_pub"])
    return pool


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--keep-tmp", action="store_true")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    print(f"[pool] loading 22-source on-MAF rsID pool …")
    on_maf = _load_on_maf_pool()
    print(f"[pool] 22-source on-MAF: {len(on_maf)} rsIDs")

    recover_rsids = (set(STANDARD_RESCUE) | set(RELAXED_RESCUE)
                       | set(BYPASS_EXCEPTIONS))
    print(f"[recover] recoverable rsIDs: {len(recover_rsids)}")
    print(f"[recover] monomorphic lost: {len(MONOMORPHIC_LOST)}")

    tmp = Path(tempfile.mkdtemp(prefix="recover_all_"))
    print(f"[tmp] {tmp}")

    # Build unfilt-name map: rsID_pub → actual name in unfiltered BIM
    # (rsID for some, kgp probe ID for others)
    rsid2unfilt = _build_unfilt_name_map()
    print(f"[map] rsID → unfilt name resolution:")
    for rs, un in sorted(rsid2unfilt.items()):
        print(f"    {rs:14s} → {un}")

    # Write rename file for kgp → rsID after extract
    rename_path = tmp / "rename_kgp_to_rsid.txt"
    with rename_path.open("w", encoding="utf-8") as f:
        for rs, un in rsid2unfilt.items():
            if un != rs:                       # only kgp rows need rename
                f.write(f"{un}\t{rs}\n")
    print(f"[map] rename file: {rename_path} "
          f"({sum(1 for rs,un in rsid2unfilt.items() if un != rs)} entries)")

    # ── Pass A: standard rescue (--geno 0.10) ───────────────────────────
    passA_rsids = sorted(set(STANDARD_RESCUE) | set(RELAXED_RESCUE))
    passA_unfilt_names = [rsid2unfilt[rs] for rs in passA_rsids]
    (tmp / "passA_extract.txt").write_text(
        "\n".join(passA_unfilt_names) + "\n", "utf-8")
    _run([str(PLINK),
           "--bfile", str(UNFILT),
           "--extract", str(tmp / "passA_extract.txt"),
           "--keep", str(KEEP_FAM),
           "--geno", "0.10",
           "--hwe", "1e-7",
           "--maf", "0.001",
           "--make-bed",
           "--out", str(tmp / "passA")],
          "passA — standard rescue (--geno 0.10)")

    # ── Pass B: bypass exceptions (no QC filters) ───────────────────────
    passB_rsids = sorted(BYPASS_EXCEPTIONS)
    passB_unfilt_names = [rsid2unfilt[rs] for rs in passB_rsids]
    (tmp / "passB_extract.txt").write_text(
        "\n".join(passB_unfilt_names) + "\n", "utf-8")
    _run([str(PLINK),
           "--bfile", str(UNFILT),
           "--extract", str(tmp / "passB_extract.txt"),
           "--keep", str(KEEP_FAM),
           "--make-bed",
           "--out", str(tmp / "passB")],
          "passB — bypass exceptions (no --geno/--maf/--hwe)")

    # ── Merge passA + passB → passAB (recovered pool only, GRCh37) ──────
    _run([str(PLINK),
           "--bfile", str(tmp / "passA"),
           "--bmerge", str(tmp / "passB"),
           "--make-bed",
           "--out", str(tmp / "passAB_kgp")],
          "merge passA + passB")

    # ── Rename kgp probe IDs → rsIDs ─────────────────────────────────────
    _run([str(PLINK),
           "--bfile", str(tmp / "passAB_kgp"),
           "--update-name", str(rename_path),
           "--make-bed",
           "--out", str(tmp / "passAB")],
          "rename kgp → rsID")

    # ── Liftover passAB GRCh37 → GRCh38 via pyliftover ──────────────────
    print(f"[liftover] applying hg19 → hg38 chain to recovered SNPs")
    try:
        from pyliftover import LiftOver
    except ImportError:
        sys.exit("[ERROR] pyliftover not installed (need it for GRCh37→38)")
    lo = LiftOver(str(CHAIN))
    bim37 = pd.read_csv(tmp / "passAB.bim", sep="\t", header=None,
                          names=["CHR", "rsID", "cM", "BP", "A1", "A2"],
                          dtype=str)
    new_pos = []
    bim37["BP_int"] = bim37["BP"].astype(int)
    for _, r in bim37.iterrows():
        hits = lo.convert_coordinate(f"chr{r['CHR']}", r["BP_int"] - 1)
        if hits:
            new_pos.append((r["rsID"], r["CHR"], int(hits[0][1]) + 1))
        else:
            print(f"  [warn] no liftover for {r['rsID']} chr{r['CHR']}:{r['BP']}")
            new_pos.append((r["rsID"], r["CHR"], int(r["BP"])))
    update_map = pd.DataFrame(new_pos, columns=["rsID", "CHR", "BP_GRCh38"])
    update_map[["rsID", "BP_GRCh38"]].to_csv(
        tmp / "update_map.txt", sep="\t", index=False, header=False)
    _run([str(PLINK),
           "--bfile", str(tmp / "passAB"),
           "--update-map", str(tmp / "update_map.txt"),
           "--make-bed",
           "--out", str(tmp / "passAB_GRCh38")],
          "apply lifted positions (GRCh38)")

    # ── Merge with the 22-source MAF BIM ────────────────────────────────
    # Restrict MAF BIM to the 116 on-MAF pool first (saves disk + speed)
    (tmp / "maf_pool.txt").write_text("\n".join(sorted(on_maf)) + "\n", "utf-8")
    _run([str(PLINK),
           "--bfile", str(MAF_BIM_STEM),
           "--extract", str(tmp / "maf_pool.txt"),
           "--keep", str(KEEP_FAM),
           "--make-bed",
           "--out", str(tmp / "maf_pool")],
          f"restrict MAF BIM to 22-source pool ({len(on_maf)} rsIDs)")

    final_stem = args.out / "SNP_recover_all_GRCh38_PRS_only"
    _run([str(PLINK),
           "--bfile", str(tmp / "maf_pool"),
           "--bmerge", str(tmp / "passAB_GRCh38"),
           "--make-bed",
           "--out", str(final_stem)],
          f"merge on-MAF pool ({len(on_maf)}) + recovered "
          f"({len(recover_rsids)}) → final BED")

    # Report on the final
    fin_bim = pd.read_csv(final_stem.with_suffix(".bim"), sep="\t", header=None,
                            names=["CHR","rsID","cM","BP","A1","A2"], dtype=str)
    fin_fam = pd.read_csv(final_stem.with_suffix(".fam"), sep=r"\s+",
                            header=None, dtype=str)
    print(f"\n[final] {final_stem}.bed")
    print(f"  SNPs: {len(fin_bim)}")
    print(f"  Samples: {len(fin_fam)}")
    rec_in = sum(1 for rs in recover_rsids if rs in set(fin_bim["rsID"]))
    pool_in = sum(1 for rs in on_maf if rs in set(fin_bim["rsID"]))
    print(f"  Recovered (of 12): {rec_in}")
    print(f"  On-MAF pool (of {len(on_maf)}): {pool_in}")

    # ── recovery_report.tsv ─────────────────────────────────────────────
    final_rsids = set(fin_bim["rsID"])
    rows = []
    for cat, mapping in [("standard_geno_0.05",   STANDARD_RESCUE),
                          ("relaxed_geno_0.10",    RELAXED_RESCUE),
                          ("bypass_exception",     BYPASS_EXCEPTIONS),
                          ("monomorphic_lost",     MONOMORPHIC_LOST)]:
        for rs, (stat, srcs, note) in mapping.items():
            rows.append({
                "rsID": rs,
                "recovery_class": cat,
                "statistic": stat,
                "sources": srcs,
                "gene_or_note": note,
                "in_final_BED": rs in final_rsids,
            })
    rep = pd.DataFrame(rows).sort_values(["recovery_class", "rsID"])
    rep.to_csv(args.out / "recovery_report.tsv", sep="\t", index=False)
    print(f"\n  [out] {args.out/'recovery_report.tsv'}")

    # ── Summary text ─────────────────────────────────────────────────────
    lines = [
        "Recover-all PRS pool — 22-source extension (2026-05-25 night)",
        f"On-MAF pool (22 sources, already in MAF BIM): {len(on_maf)} rsIDs",
        f"Recoverables added: {len(recover_rsids)} "
        f"({len(STANDARD_RESCUE)} standard + {len(RELAXED_RESCUE)} relaxed "
        f"+ {len(BYPASS_EXCEPTIONS)} bypass)",
        f"Monomorphic — TRULY lost: {len(MONOMORPHIC_LOST)} "
        f"({', '.join(MONOMORPHIC_LOST)})",
        f"Final BED SNP count: {len(fin_bim)}",
        f"Final BED sample count: {len(fin_fam)} (MRI-keep cohort)",
        "",
        "QC thresholds:",
        "  Pass A (standard rescue): --geno 0.10 --hwe 1e-7 --maf 0.001",
        "  Pass B (bypass exceptions): no --geno / --maf / --hwe filters",
        "  Final MAF retained: --maf 0.001 (lets rare PRS-published leads through)",
        "",
        "Recovery breakdown:",
        rep.to_string(index=False),
    ]
    (args.out / "recovery_summary.txt").write_text(
        "\n".join(lines) + "\n", "utf-8")
    print(f"  [out] {args.out/'recovery_summary.txt'}")

    if not args.keep_tmp:
        import shutil
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
