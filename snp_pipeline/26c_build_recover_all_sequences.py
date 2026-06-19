r"""
26c_build_recover_all_sequences.py   (LOCAL, env: snp)
=====================================================
Build the SNP-centred 1001 bp REF/ALT (+ reverse-complement) sequences
for the full 128-SNP recover-all pool:
  - 116 on-MAF SNPs (union across 22 PRS sources)
  - 12 recovered SNPs (5 standard --geno 0.05 + 4 relaxed --geno 0.10
    + 3 bypass-exception published PRS leads)
  - excludes audit-only sources (Huang, ONeil_SST_candidates)
  - excludes the 2 monomorphic SNPs (rs429358 APOE ε4, rs75932628 TREM2
    R47H) — no within-cohort variation so cannot have ALT dosage anyway

The recover-all BED already has all 128 SNPs on GRCh38 with correct
A1/A2 orientation (--ref-correct earlier in the pipeline). We pull
CHR/BP/A1/A2 from its BIM, gene + β from the per-source resolution
TSVs (priority: B > D > W > K > others), then invoke script 26 to
build the sequences.

Pipeline:
  1. Load recover_all BIM for CHR/BP_GRCh38/A1/A2 of 128 rsIDs.
  2. For each rsID, pull a representative β + gene from the
     per-source resolution TSVs (any source with non-empty β works —
     the sequence doesn't depend on β; β is needed only for script 26's
     QC validation).
  3. Write `recover_all_pool_snps.tsv` to
     `diff_attn_drive_upload/recover_all_pool/`.
  4. Subprocess-invoke `26_build_fm_snp_sequences.py
     --set recover_all_pool --gwas-root diff_attn_drive_upload
     --emit-rc` to produce 4 × 128 = 512 sequence rows.

Usage:
  conda run -n snp python snp_pipeline/26c_build_recover_all_sequences.py
"""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path

import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
RECOVER_BED = (BASE / "GWAS_comprehensive_v2" / "recover_all_pool"
                / "SNP_recover_all_GRCh38_PRS_only")
RECON_DIR = BASE / "source_prs" / "unfiltered_SNP_reconciliation"
RECOVERY_REPORT = (BASE / "GWAS_comprehensive_v2" / "recover_all_pool"
                     / "recovery_report.tsv")
UPLOAD_ROOT = BASE / "diff_attn_drive_upload"
SCRIPT_26 = Path(__file__).with_name("26_build_fm_snp_sequences.py")

# Priority order for picking β/gene per rsID. Bellenguez/Wightman/Kunkle/
# Desikan have the longest-validated β values; the new sources fall after.
SOURCE_PRIORITY = [
    "Bellenguez", "Wightman", "Kunkle", "Desikan",
    "Lambert", "DeRojas", "Schwanzentruber", "Najar",
    "Ebenau", "Leonenko", "Vesilievick", "Zhang",
    "Kosteridis_novel_AD", "Kosteridis_shared_AD_CV",
    "Yang_Mic", "Yang_Ast", "Yang_Oli", "Yang_Opc", "Yang_Ex", "Yang_In",
    "Felsky_MF", "Felsky_IT",
    "ONeil_NPY", "ONeil_GHR",
]


def _load_per_source_lookups() -> dict[str, dict]:
    """Return {rsID: {source, gene, beta, ea_orig, oa_orig, recovery_class}}
    picking the highest-priority source per rsID. β orientation from this
    source is reported in `beta_A1` (matches the effect_allele column from
    that source — not necessarily the recover_all BIM's A1, which is why
    the script 26 validates and strand-flips as needed)."""
    pool: dict[str, dict] = {}
    for src in SOURCE_PRIORITY:
        p = RECON_DIR / f"{src}_unfiltered_reconciliation.tsv"
        if not p.exists():
            continue
        df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
        for _, r in df.iterrows():
            rs = str(r["rsID_pub"])
            if rs in pool:
                continue
            beta = str(r.get("OR_or_beta_pub", "")).strip()
            if not beta:
                continue
            pool[rs] = {
                "source": src,
                "gene": str(r.get("locus_name", "")).strip(),
                "ea_orig": str(r.get("effect_allele_pub", "")).strip(),
                "oa_orig": str(r.get("other_allele_pub", "")).strip(),
                "beta": beta,
            }
    return pool


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set-name", default="recover_all_pool")
    ap.add_argument("--upload-root", type=Path, default=UPLOAD_ROOT)
    args = ap.parse_args()

    out_dir = args.upload_root / args.set_name
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Stage 1: read recover-all BIM → 128 rsIDs with CHR/BP/A1/A2 ───
    print(f"[bim] loading {RECOVER_BED.name}.bim …")
    bim = pd.read_csv(RECOVER_BED.with_suffix(".bim"), sep=r"\s+",
                       header=None, dtype=str,
                       names=["CHR", "rsID", "cM", "BP_GRCh38", "A1", "A2"])
    bim["A1"] = bim["A1"].str.upper()
    bim["A2"] = bim["A2"].str.upper()
    print(f"  {len(bim)} SNPs in recover-all BED")

    # ── Stage 2: pick β/gene/source per rsID ───────────────────────────
    print(f"[source] loading per-source resolutions for β + gene …")
    lookup = _load_per_source_lookups()
    print(f"  {len(lookup)} rsIDs have a usable β across the 22 sources")

    # ── Stage 3: recovery_class for `origin` column ─────────────────────
    print(f"[recov] loading recovery_report.tsv …")
    rec = pd.read_csv(RECOVERY_REPORT, sep="\t", dtype=str, keep_default_na=False)
    rec_class = dict(zip(rec["rsID"], rec["recovery_class"]))

    # ── Stage 4: assemble snps.tsv rows in script-26 schema ────────────
    rows = []
    missing_beta = []
    for snp_idx, r in enumerate(bim.itertuples(index=False)):
        rs = str(r.rsID)
        info = lookup.get(rs, {})
        beta = info.get("beta", "")
        if not beta:
            # Recoverable SNPs (esp. the F_MISS=NA ones from sumstats
            # fallback) may have empty β in some sources but still need
            # a sequence built. Fall back to 1.0 — the sequence itself
            # is β-invariant; script 26 needs only a numeric value.
            beta = "1.0"
            missing_beta.append(rs)
        try:
            beta_val = float(beta)
        except ValueError:
            beta_val = 1.0
            missing_beta.append(rs)
        # A1 = ALT = effect_allele (PLINK convention after ref-correct);
        # A2 = REF = other_allele
        rows.append({
            "rsID": rs,
            "gene": info.get("gene", "") or "",
            "lead_rsID": rs,
            "source": info.get("source", "recovery_only"),
            "origin": rec_class.get(rs, "on_MAF_pool"),
            "CHR": r.CHR,
            "BP_GRCh38": int(r.BP_GRCh38),
            "effect_allele": r.A1,
            "other_allele": r.A2,
            "REF": r.A2,
            "ALT": r.A1,
            "beta_A1": beta_val,
        })

    snps_df = pd.DataFrame(rows)
    snps_path = out_dir / f"{args.set_name}_snps.tsv"
    snps_df.to_csv(snps_path, sep="\t", index=False)
    print(f"\n[snps] {snps_path}  ({len(snps_df)} SNPs)")
    if missing_beta:
        print(f"  [warn] {len(missing_beta)} rsIDs had no β in any source — "
              f"using β=1.0 placeholder (sequence is β-invariant): "
              f"{missing_beta[:5]}{'…' if len(missing_beta) > 5 else ''}")
    print(snps_df[["rsID", "gene", "source", "origin", "CHR",
                    "BP_GRCh38", "effect_allele", "other_allele",
                    "beta_A1"]].head(10).to_string(index=False))

    # Also write a no-op patient_dosage.tsv with just a header row so
    # script 26's `in_dosage` check finds it (but doesn't dictate drops).
    # (script 26 only WARNS on not_in_dosage, doesn't drop.)
    dos_path = out_dir / f"{args.set_name}_patient_dosage.tsv"
    if not dos_path.exists():
        cols = ["PTID"] + list(snps_df["rsID"])
        pd.DataFrame(columns=cols).to_csv(dos_path, sep="\t", index=False)
        print(f"[dos] {dos_path}  (header-only placeholder, 0 patients)")

    # ── Stage 5: invoke 26 to build sequences ────────────────────────────
    print(f"\n[26] building 1001-bp REF/ALT fwd+rc sequences …")
    cmd = [sys.executable, str(SCRIPT_26),
           "--set", args.set_name,
           "--gwas-root", str(args.upload_root),
           "--emit-rc"]
    print("  " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=False, text=True)
    if r.returncode != 0:
        raise SystemExit(f"script 26 failed (rc={r.returncode})")

    # ── Stage 6: report bundle contents ──────────────────────────────────
    seq_dir = out_dir / "fm_sequences"
    print(f"\n[bundle] {out_dir}")
    for p in sorted(out_dir.rglob("*")):
        if p.is_file():
            print(f"  {p.relative_to(out_dir)}  ({p.stat().st_size:,} B)")

    # Quick verify
    seq_tsv = seq_dir / f"{args.set_name}_snp_sequences.tsv"
    if seq_tsv.exists():
        sdf = pd.read_csv(seq_tsv, sep="\t", dtype=str)
        n_snps = sdf["rsID"].nunique()
        n_rows = len(sdf)
        seqlens = sdf["sequence"].str.len().value_counts() if "sequence" in sdf.columns else {}
        print(f"\n[verify] sequence TSV: {n_rows} rows / {n_snps} unique SNPs")
        if len(seqlens):
            print(f"  sequence length distribution: {dict(seqlens.head())}")


if __name__ == "__main__":
    main()
