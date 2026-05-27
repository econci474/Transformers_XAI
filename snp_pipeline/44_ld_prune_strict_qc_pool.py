"""
44_ld_prune_strict_qc_pool.py
=============================
LD-prune the strict-QC 115-SNP recover_all_pool at r² > 0.8 within 1 MB
physical-distance windows. Produces an LD-independent subset for downstream
classification + Cox + AAO regression (per the user's existing methodology
preference — see feedback_strict_qc_ld_pruning_first).

PLINK call: --indep-pairwise 1000kb 1 0.8
  - 1000kb : 1 MB sliding window (physical distance, not SNP count)
  - 1      : step = 1 variant (slide by one SNP at a time)
  - 0.8    : keep pairwise r² < 0.8 (drop variants tagged by an earlier kept SNP)

Outputs (under QC_strict/):
  recover_all_pool_strictQC_LDpruned.bed/bim/fam   — survivors only
  recover_all_pool_strictQC_LDpruned.prune.in       — kept rsIDs
  recover_all_pool_strictQC_LDpruned.prune.out      — dropped rsIDs
  recover_all_pool_strictQC_LDpruned.log            — PLINK log
  recover_all_pool_strictQC_LDpruned_dosage.tsv     — PTID × kept-rsID matrix
  ld_pruning_report.tsv                              — rsID-level kept/dropped + LD-tag (if dropped)
  ld_pruning_summary.txt                             — terminal-style summary

Usage (env: snp):
  python snp_pipeline/44_ld_prune_strict_qc_pool.py
  python snp_pipeline/44_ld_prune_strict_qc_pool.py --dry-run
"""
from __future__ import annotations
import argparse
import subprocess
import time
from pathlib import Path
import pandas as pd

ROOT     = Path("D:/ADNI_SNP_Omni2.5M_20140220")
PLINK    = ROOT / "plink.exe"
QC_DIR   = ROOT / "GWAS_comprehensive_v2" / "QC_strict"
SRC_BFILE = QC_DIR / "recover_all_pool_strictQC"

# Default config (user's stated preference) — see feedback_strict_qc_ld_pruning_first.
DEFAULT_LD_CONFIG = "ld_1000kb_r2_0.8"
LD_CONFIGS = {
    "ld_1000kb_r2_0.8": ("1000kb", "1", "0.8"),
    "ld_500kb_r2_0.2":  ("500kb",  "1", "0.2"),
    "ld_250kb_r2_0.1":  ("250kb",  "1", "0.1"),
    "ld_50_5_r2_0.5":   ("50",     "5", "0.5"),   # PLINK tutorial default
}


def _run(cmd: list[str], label: str = ""):
    print(f"\n>>> {' '.join(cmd)}")
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    if r.returncode != 0:
        print("STDOUT:\n", r.stdout[-1500:])
        print("STDERR:\n", r.stderr[-1500:])
        raise SystemExit(f"PLINK failed (exit {r.returncode}) at: {label or cmd[0]}")
    print(f"    done in {dt:.1f}s")


def _do_one_config(ld_config: str, dry_run: bool = False):
    LD_WINDOW, LD_STEP, LD_R2 = LD_CONFIGS[ld_config]
    out_dir = QC_DIR / ld_config
    out_dir.mkdir(parents=True, exist_ok=True)
    OUT_BFILE = out_dir / "recover_all_pool_strictQC_LDpruned"
    print(f"\n========== {ld_config}  (--indep-pairwise {LD_WINDOW} {LD_STEP} {LD_R2}) ==========")

    if dry_run:
        print(f"[DRY] would write to {out_dir}")
        return

    # 1. Generate prune.in / prune.out via --indep-pairwise
    _run([str(PLINK),
           "--bfile", str(SRC_BFILE),
           "--indep-pairwise", LD_WINDOW, LD_STEP, LD_R2,
           "--out", str(OUT_BFILE)],
          label=f"LD pruning ({ld_config})")
    prune_in_p  = OUT_BFILE.with_suffix(".prune.in")
    prune_out_p = OUT_BFILE.with_suffix(".prune.out")

    kept    = pd.read_csv(prune_in_p,  header=None, names=["rsID"])
    dropped = pd.read_csv(prune_out_p, header=None, names=["rsID"])
    print(f"\nKept: {len(kept)}  Dropped: {len(dropped)}  "
          f"(of {len(kept)+len(dropped)} strict-QC SNPs)")

    # 2. Materialise the kept-SNP BED via --extract
    _run([str(PLINK),
           "--bfile", str(SRC_BFILE),
           "--extract", str(prune_in_p),
           "--make-bed",
           "--out", str(OUT_BFILE)],
          label="extract kept SNPs into LD-pruned BED")

    # 3. Dosage matrix on the LD-pruned set
    raw_prefix = OUT_BFILE.with_name(OUT_BFILE.name + "_dosage")
    _run([str(PLINK),
           "--bfile", str(OUT_BFILE),
           "--recode", "A",
           "--out", str(raw_prefix)],
          label="dosage matrix --recode A")
    raw_p = raw_prefix.with_suffix(".raw")
    raw = pd.read_csv(raw_p, sep=r"\s+")
    snp_cols = [c for c in raw.columns if c not in
                  ("FID","IID","PAT","MAT","SEX","PHENOTYPE")]
    rename = {c: c.split("_")[0] for c in snp_cols}
    dosage = raw[["IID"] + snp_cols].rename(columns={"IID": "PTID", **rename})
    dosage_p = raw_prefix.with_suffix(".tsv")
    dosage.to_csv(dosage_p, sep="\t", index=False)
    print(f"wrote dosage TSV: {dosage_p}  shape={dosage.shape}")

    # 4. Pruning report (per-rsID kept/dropped + which kept SNP "tagged" each dropped)
    strict_bim = pd.read_csv(SRC_BFILE.with_suffix(".bim"), sep=r"\s+",
                              header=None, names=["chrom","rsID","cM","bp","A1","A2"])
    kept_set = set(kept["rsID"])
    rows = []
    for _, r in strict_bim.iterrows():
        rs = r["rsID"]
        if rs in kept_set:
            rows.append({"rsID": rs, "chrom": r["chrom"], "bp": r["bp"],
                          "status": "kept", "tagged_by": ""})
        else:
            same_chrom_kept = strict_bim[(strict_bim["chrom"] == r["chrom"]) &
                                          (strict_bim["rsID"].isin(kept_set))]
            if len(same_chrom_kept) == 0:
                tag = ""
            else:
                nearest = same_chrom_kept.iloc[
                    (same_chrom_kept["bp"] - r["bp"]).abs().argmin()]
                tag = f"{nearest['rsID']} (dist={abs(nearest['bp']-r['bp'])} bp)"
            rows.append({"rsID": rs, "chrom": r["chrom"], "bp": r["bp"],
                          "status": "dropped", "tagged_by": tag})
    report = pd.DataFrame(rows)
    rep_p = out_dir / "ld_pruning_report.tsv"
    report.to_csv(rep_p, sep="\t", index=False)
    print(f"wrote pruning report: {rep_p}")

    # 5. Terminal summary
    sum_p = out_dir / "ld_pruning_summary.txt"
    text = [
        f"LD pruning of strict-QC pool [{ld_config}] — {pd.Timestamp.now():%Y-%m-%d %H:%M}",
        f"Input BED: {SRC_BFILE}.bed  ({len(strict_bim)} variants, 616 subj)",
        f"PLINK args: --indep-pairwise {LD_WINDOW} {LD_STEP} {LD_R2}",
        "",
        f"  Kept    : {len(kept)}",
        f"  Dropped : {len(dropped)}",
        "",
        f"Per-chromosome breakdown:",
    ]
    chrom_summary = report.groupby(["chrom","status"]).size().unstack(fill_value=0)
    for c in sorted(chrom_summary.index):
        k = int(chrom_summary.loc[c].get("kept", 0))
        d = int(chrom_summary.loc[c].get("dropped", 0))
        text.append(f"  chr {c:>2d}: kept={k:>3d}  dropped={d:>3d}")
    text.append("")
    text.append(f"Output BED:    {OUT_BFILE}.bed/.bim/.fam")
    text.append(f"Dosage TSV:    {dosage_p}")
    text.append(f"Report TSV:    {rep_p}")
    sum_p.write_text("\n".join(text) + "\n")
    print("\n" + "\n".join(text))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--ld-config", default=None,
                     help=f"Run only this config. Choices: {list(LD_CONFIGS)}. "
                          f"Default: all 4 configs.")
    args = ap.parse_args()

    for ext in (".bed", ".bim", ".fam"):
        if not SRC_BFILE.with_suffix(ext).exists():
            raise SystemExit(f"Missing strict-QC BED: {SRC_BFILE}{ext}")
    if not PLINK.exists():
        raise SystemExit(f"plink.exe not found at {PLINK}")
    QC_DIR.mkdir(parents=True, exist_ok=True)

    if args.ld_config:
        if args.ld_config not in LD_CONFIGS:
            raise SystemExit(f"Unknown ld_config {args.ld_config!r}")
        _do_one_config(args.ld_config, dry_run=args.dry_run)
    else:
        for cfg in LD_CONFIGS:
            _do_one_config(cfg, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
