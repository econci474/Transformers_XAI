r"""
25e_export_missingness_resolved_sets.py   (LOCAL, env: snp)
===========================================================
Build the **missingness-resolved** PRS bundles by re-running PLINK on the
unfiltered Omni2.5M BED with `--geno 0.05` (relaxed from the canonical
`--geno 0.02`), plus two published-PRS exceptions whose probes bypass the
QC filters:

  Standard rescues (pass --geno 0.05):
    rs7401792    Bellenguez chr14    F_MISS=0.0345
    rs117618017  Bellenguez+Wightman F_MISS=0.0222
  Filter-bypass exceptions (published-PRS leads kept regardless):
    rs4266886    Desikan PGS CR1     F_MISS=0.1195
    rs115186657  Wightman            MAF=0.0037

The 2 monomorphic probes (rs75932628 TREM2, rs429358 APOE ε4) remain
excluded — no variation in the consent-cleared cohort means they cannot
contribute to a PRS.

Mirrors `25d_export_maf_resolved_sets.py` for output schema; differs only
in the dosage source (a freshly built BED at relaxed missingness) and
the explicit handling of the bypass list.

Output → `D:\ADNI_SNP_Omni2.5M_20140220\GWAS_filtered_missingness_resolved_geno005\`:
  sourceB/, sourceW/, sourceD/, sourceK/  (per-source PRS bundles)
  5W26B/, 5W26B14D/, 5W26B13D/, 5W26B14D2K/  (union combos)
  Each combo: <combo>_snps.tsv, <combo>_patient_dosage.tsv,
              <combo>_patient_dosage.bim, README.txt
  combos_summary.tsv   tally per combo (n SNPs by source)
  rescue_log.tsv       which SNPs rescued / still dropped / why

Intermediate BED:
  D:\…\SNP_filtered_with_mri_geno005_PRS_only_GRCh38_maf001.{bed,bim,fam}

Usage:
  conda run -n snp python snp_pipeline/25e_export_missingness_resolved_sets.py
"""
from __future__ import annotations

import argparse
import importlib.util
import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
UNFILT = BASE / "WGS_Omni25_BIN_wo_ConsentsIssues"
KEEP_FAM = BASE / "SNP_filtered_with_mri.fam"
CHAIN = BASE / "liftover" / "hg19ToHg38.over.chain.gz"
SRC_PRS = BASE / "source_prs"
GBF = BASE / "genetic_baselines_filtered"
PLINK = BASE / "plink.exe"

# GWAS-Catalogue Wightman 2021 sumstats (source-of-truth for β)
WIGHT_GCST_DIR = BASE / "GWAS" / "Wightman_2021" / "GWAS_Catalogue_summary_stats"
WIGHT_GCST196 = WIGHT_GCST_DIR / "GCST013196.tsv.gz"   # no-UKB, β direct
WIGHT_GCST197 = WIGHT_GCST_DIR / "GCST013197.tsv.gz"   # with-UKB, Z + N

OUT_DEFAULT = BASE / "GWAS_filtered_missingness_resolved_geno005"
NEW_BED_STEM = BASE / "SNP_filtered_with_mri_geno005_PRS_only_GRCh38_maf001"

# Published-PRS exceptions — extracted without --geno / --hwe / --maf
BYPASS_FILTERS = {"rs4266886", "rs115186657"}

# Reuse helpers from 25d (β harmonisation, write_combo) and 30c (Illumina map)
SCRIPT_25D = Path(__file__).with_name("25d_export_maf_resolved_sets.py")
_s25d_spec = importlib.util.spec_from_file_location("_s25d", SCRIPT_25D)
_s25d = importlib.util.module_from_spec(_s25d_spec)
_s25d_spec.loader.exec_module(_s25d)

SCRIPT_30C = Path(__file__).with_name("30c_unfiltered_prs_reconciliation.py")
_s30c_spec = importlib.util.spec_from_file_location("_s30c", SCRIPT_30C)
_s30c = importlib.util.module_from_spec(_s30c_spec)
_s30c_spec.loader.exec_module(_s30c)

harmonise = _s25d.harmonise
dedupe_union = _s25d.dedupe_union
write_combo = _s25d.write_combo
SRC_IS_OR = _s25d.SRC_IS_OR
PRIORITY = _s25d.PRIORITY
WIGHT_ALLOWLIST = _s25d.WIGHT_ALLOWLIST
WFILT_ENR = _s25d.WFILT_ENR
_load_illumina_rsid_to_kgp = _s30c._load_illumina_rsid_to_kgp

SOURCES = ("Bellenguez", "Wightman", "Desikan", "Kunkle")


def _run(cmd: list[str], desc: str) -> None:
    print(f"  [plink] {desc}")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"    STDERR (tail): {r.stderr[-1500:]}")
        raise RuntimeError(f"plink failed for: {desc}")


def _build_kgp_extract_and_update(pub_rsids: set[str],
                                  rsid2kgp: dict[str, list[tuple[str, str]]]
                                  ) -> tuple[list[str], dict[str, str]]:
    """For each published rsID, find a kgp predecessor (if any). Return
    (extract_list = published_rsids + winning_kgp_per_rsID,
     update_name = {kgp: rsID} for plink --update-name).
    Pick the first kgp predecessor (typically only one per rsID)."""
    extract = set(pub_rsids)
    update: dict[str, str] = {}
    for rs in pub_rsids:
        cands = rsid2kgp.get(rs, [])
        if not cands:
            continue
        kgp = cands[0][0]   # winning kgp = first listed
        extract.add(kgp)
        update[kgp] = rs
    return sorted(extract), update


def _bim_index(bim_path: Path) -> pd.DataFrame:
    return pd.read_csv(bim_path, sep="\t", header=None,
                       names=["CHR", "rsID", "cM", "BP", "A1", "A2"],
                       dtype=str).set_index("rsID")


def _liftover_bim(bim_path: Path, chain_path: Path) -> None:
    """In-place liftover of a small BIM file GRCh37 → GRCh38 using
    pyliftover. Mirrors `06_liftover_hg19_to_hg38.py` minus the
    industrial-scale handling — for <500 SNPs this is fast."""
    from pyliftover import LiftOver
    print(f"  [liftover] {bim_path.name}  GRCh37 → GRCh38")
    lo = LiftOver(str(chain_path))
    bim = pd.read_csv(bim_path, sep="\t", header=None,
                      names=["CHR", "rsID", "cM", "BP", "A1", "A2"],
                      dtype=str)
    new_chr, new_bp, dropped = [], [], 0
    for chrom, bp in zip(bim["CHR"], bim["BP"]):
        try:
            bp0 = int(bp) - 1
        except ValueError:
            new_chr.append(chrom); new_bp.append(bp); dropped += 1
            continue
        hits = lo.convert_coordinate(f"chr{chrom}", bp0)
        if not hits:
            new_chr.append(chrom); new_bp.append(bp); dropped += 1
            continue
        h_chr = hits[0][0].lstrip("chr")
        h_bp = hits[0][1] + 1
        new_chr.append(h_chr); new_bp.append(str(h_bp))
    bim["CHR"] = new_chr
    bim["BP"] = new_bp
    bim.to_csv(bim_path, sep="\t", header=False, index=False)
    print(f"  [liftover] {len(bim)} SNPs, {dropped} unmapped (positions left unchanged)")


def _build_new_bed(out_root: Path, rsid2kgp: dict) -> Path:
    """Build the missingness-resolved BED via the two-pass + merge +
    liftover + MAF + bypass-restore pipeline. Returns the final BED stem."""
    tmp = out_root / "_build_tmp"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)

    # ── union of all published rsIDs (collect from 4 resolution TSVs) ────
    pub_rs: set[str] = set()
    for src in SOURCES:
        df = pd.read_csv(SRC_PRS / f"{src}_full_snps_resolution.tsv",
                         sep="\t", dtype=str, keep_default_na=False)
        for rs in pd.concat([df["rsID_pub"], df["rsID_canonical"]]):
            if isinstance(rs, str) and rs.startswith("rs"):
                pub_rs.add(rs)
    print(f"  union of published rsIDs (incl. canonical): {len(pub_rs)}")

    extract_all, update_name_map = _build_kgp_extract_and_update(pub_rs, rsid2kgp)
    print(f"  extract list size: {len(extract_all)} "
          f"(rsIDs + kgp probes); update_name pairs: {len(update_name_map)}")

    # ── build helper files ──────────────────────────────────────────────
    bypass_kgp = {kgp for kgp, rs in update_name_map.items() if rs in BYPASS_FILTERS}
    bypass_extract = sorted({rs for rs in BYPASS_FILTERS} | bypass_kgp)
    standard_extract = sorted(set(extract_all) - set(bypass_extract))

    (tmp / "extract_standard.txt").write_text(
        "\n".join(standard_extract) + "\n", "utf-8")
    (tmp / "extract_bypass.txt").write_text(
        "\n".join(bypass_extract) + "\n", "utf-8")
    update_lines = "\n".join(f"{old}\t{new}" for old, new
                              in sorted(update_name_map.items()))
    (tmp / "update_name.txt").write_text(update_lines + "\n", "utf-8")

    # ── (4a) Standard pass — full QC filters ────────────────────────────
    # NOTE: --mind 0.02 is INTENTIONALLY OMITTED. The canonical cohort
    # (SNP_filtered_with_mri.fam) was already produced post --mind 0.02
    # at step 4 of 01_qc.sh. Re-applying --mind on a tiny ~50-SNP subset
    # would spuriously drop patients (1 missing genotype = 2% on a 50-SNP
    # set; 2 missing = 4%) — confirmed in the v1 run which dropped
    # 072_S_4383 and 941_S_4420.
    _run([str(PLINK), "--bfile", str(UNFILT),
          "--extract", str(tmp / "extract_standard.txt"),
          "--autosome", "--geno", "0.05",
          "--hwe", "1e-7",
          "--keep", str(KEEP_FAM),
          "--update-name", str(tmp / "update_name.txt"),
          "--keep-allele-order",
          "--make-bed", "--out", str(tmp / "step1_standard")],
         "4a) standard pass --geno 0.05 --hwe 1e-7")

    # ── (4b) Bypass pass — same cohort, NO --geno/--hwe checks ──────────
    _run([str(PLINK), "--bfile", str(UNFILT),
          "--extract", str(tmp / "extract_bypass.txt"),
          "--autosome",
          "--keep", str(KEEP_FAM),
          "--update-name", str(tmp / "update_name.txt"),
          "--keep-allele-order",
          "--make-bed", "--out", str(tmp / "step1_bypass")],
         "4b) bypass pass (no filters)")

    # ── (4c) Merge ──────────────────────────────────────────────────────
    _run([str(PLINK), "--bfile", str(tmp / "step1_standard"),
          "--bmerge", str(tmp / "step1_bypass"),
          "--keep-allele-order",
          "--make-bed", "--out", str(tmp / "step1_merged")],
         "4c) merge standard + bypass")

    # ── (5) Liftover GRCh37 → GRCh38 (in-place on BIM only) ─────────────
    _liftover_bim(tmp / "step1_merged.bim", CHAIN)

    # ── (6a) Apply --maf 0.01 ──────────────────────────────────────────
    _run([str(PLINK), "--bfile", str(tmp / "step1_merged"),
          "--maf", "0.01",
          "--keep-allele-order",
          "--make-bed", "--out", str(tmp / "step6_maf")],
         "6a) --maf 0.01")

    # ── (6b) Bypass restore if any exception got MAF-dropped ───────────
    maf_bim_rs = set(_bim_index(tmp / "step6_maf.bim").index)
    missing_bypass = sorted(BYPASS_FILTERS - maf_bim_rs)
    if missing_bypass:
        print(f"  MAF dropped exception(s): {missing_bypass} — restoring")
        (tmp / "bypass_restore.txt").write_text(
            "\n".join(missing_bypass) + "\n", "utf-8")
        _run([str(PLINK), "--bfile", str(tmp / "step1_merged"),
              "--extract", str(tmp / "bypass_restore.txt"),
              "--keep-allele-order",
              "--make-bed", "--out", str(tmp / "step6_bypass_restore")],
             "6b-extract) bypass restore")
        _run([str(PLINK), "--bfile", str(tmp / "step6_maf"),
              "--bmerge", str(tmp / "step6_bypass_restore"),
              "--keep-allele-order",
              "--make-bed", "--out", str(NEW_BED_STEM)],
             "6b-merge) final BED")
    else:
        for ext in (".bed", ".bim", ".fam"):
            shutil.copy(tmp / f"step6_maf{ext}",
                        Path(str(NEW_BED_STEM) + ext))
        print("  no bypass restore needed — copied step6_maf as final BED")

    final_n = len(_bim_index(Path(str(NEW_BED_STEM) + ".bim")))
    print(f"  [final BED] {NEW_BED_STEM}  ({final_n} SNPs)")
    return NEW_BED_STEM


def plink_extract_dosage(rsids: list[str], bfile: Path, outdir: Path
                          ) -> pd.DataFrame:
    """plink --bfile <new_BED> --extract <rsids> --recode A → wide dosage."""
    outdir.mkdir(parents=True, exist_ok=True)
    snplist = outdir / "rsid_extract.txt"
    snplist.write_text("\n".join(rsids) + "\n", "utf-8")
    stem = outdir / "patient_dosage"
    _run([str(PLINK), "--bfile", str(bfile),
          "--extract", str(snplist),
          "--keep-allele-order",
          "--recode", "A",
          "--out", str(stem)],
         f"--recode A (n={len(rsids)})")
    raw = pd.read_csv(stem.with_suffix(".raw"), sep=r"\s+", engine="python")
    pid = raw["IID"].astype(str)
    geno = raw.drop(columns=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"])
    geno.columns = [c.rsplit("_", 1)[0] for c in geno.columns]
    geno.index = pid
    geno.index.name = "PTID"
    return geno


def _load_gcst_wightman_beta(targets_chrpos: set[tuple[str, str]]
                              ) -> tuple[dict, dict]:
    """Stream GCST013196 (no-UKB, β direct) and GCST013197 (with-UKB,
    Z + N, SE=NA) once, caching by (chr, pos_GRCh37) for the target
    Wightman PRS positions.

    Returns (gcst196, gcst197) dicts:
      gcst196: {(chr, pos): {beta, SE, EA, OA, EAF, N, source='GCST013196'}}
      gcst197: {(chr, pos): {Z, SE_computed, beta_derived, EA, OA, EAF, N,
                              source='GCST013197_Z*SE'}}"""
    import gzip
    gcst196: dict = {}
    gcst197: dict = {}

    if WIGHT_GCST196.exists() and targets_chrpos:
        print(f"  [gcst] streaming {WIGHT_GCST196.name} ({len(targets_chrpos)} targets)")
        with gzip.open(WIGHT_GCST196, "rt") as fh:
            next(fh, None)
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 8:
                    continue
                k = (p[0], p[1])
                if k in targets_chrpos:
                    try:
                        gcst196[k] = {
                            "EA": p[2].upper(), "OA": p[3].upper(),
                            "beta": float(p[4]), "SE": float(p[5]),
                            "EAF": float(p[6]), "p": p[7],
                            "N": int(float(p[8])) if len(p) > 8 else None,
                            "source": "GCST013196_no-UKB_beta",
                        }
                    except (ValueError, IndexError):
                        pass
                if len(gcst196) == len(targets_chrpos):
                    break
        print(f"  [gcst] GCST013196 hits: {len(gcst196)}/{len(targets_chrpos)}")

    # Tier 2: streamed only for SNPs still missing
    remaining = targets_chrpos - set(gcst196.keys())
    if WIGHT_GCST197.exists() and remaining:
        print(f"  [gcst] streaming {WIGHT_GCST197.name} ({len(remaining)} remaining)")
        with gzip.open(WIGHT_GCST197, "rt") as fh:
            next(fh, None)
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 8:
                    continue
                k = (p[0], p[1])
                if k in remaining:
                    try:
                        z = float(p[4])
                        eaf = float(p[6])
                        n = int(float(p[8])) if len(p) > 8 else None
                        maf = min(eaf, 1 - eaf)
                        se = (1.0 / (2 * maf * (1 - maf) * n)) ** 0.5
                        gcst197[k] = {
                            "EA": p[2].upper(), "OA": p[3].upper(),
                            "Z": z, "SE_computed": se,
                            "beta": z * se, "EAF": eaf, "N": n,
                            "p": p[7],
                            "source": "GCST013197_with-UKB_Z*SE_computed",
                        }
                    except (ValueError, IndexError, ZeroDivisionError):
                        pass
                if len(gcst197) == len(remaining):
                    break
        print(f"  [gcst] GCST013197 hits: {len(gcst197)}/{len(remaining)}")
    return gcst196, gcst197


def _wightman_targets_chrpos() -> dict[str, tuple[str, str]]:
    """rsID → (CHR_GRCh37, BP_GRCh37) for every published Wightman SNP
    in the resolution TSV (regardless of drop_reason)."""
    p = SRC_PRS / "Wightman_full_snps_resolution.tsv"
    df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    out: dict[str, tuple[str, str]] = {}
    for _, r in df.iterrows():
        rs = str(r["rsID_pub"])
        chrom = str(r["CHR_pub"]).strip()
        bp = str(r["BP_pub"]).strip()
        if rs and chrom and bp:
            out[rs] = (chrom, bp)
    return out


def load_resolved_with_rescue(source: str, new_bim: pd.DataFrame
                                ) -> pd.DataFrame:
    """Like 25d.load_resolved but uses the NEW BED's CHR/BP/A1/A2 (after
    --geno 0.05 + liftover + maf-with-bypass). For every published rsID
    that survives in `new_bim`, build the harmonised row.
    Marks `filter_exception=True` for the bypass list."""
    p = SRC_PRS / f"{source}_full_snps_resolution.tsv"
    df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)

    # Wightman: GCST primary chain (per user 2026-05-23):
    #   Tier 1: GCST013196 (no-UKB, β direct)
    #   Tier 2: GCST013197 (with-UKB, Z + N, SE computed) → β = Z × SE
    #   Cross-check legacy wmap (allowlist + wfilt_enriched); WARN on
    #   |Δβ| > 0.02 (absolute, allele-flip-normalised)
    if source == "Wightman":
        targets = _wightman_targets_chrpos()
        target_chrpos = set(targets.values())
        gcst196, gcst197 = _load_gcst_wightman_beta(target_chrpos)

        # legacy wmap (for cross-check only — NOT a source any more)
        wmap_legacy = {}
        allow = pd.read_csv(WIGHT_ALLOWLIST, sep="\t")
        wmap_legacy.update(dict(zip(allow["pipeline_rsID"].astype(str),
                                     allow["beta"].astype(float))))
        wfilt = pd.read_csv(WFILT_ENR, sep="\t")
        wfilt = wfilt[(wfilt["novel"] == True) &                 # noqa: E712
                      wfilt["source"].astype(str).str.contains("Wightman")]
        wmap_legacy.update(dict(zip(wfilt["rsID"].astype(str),
                                     wfilt["beta_A1"].astype(float))))

        rows = []
        for _, r in df.iterrows():
            rs_pub = str(r["rsID_pub"])
            rs_can = str(r["rsID_canonical"])
            rs = next((x for x in (rs_pub, rs_can) if x in new_bim.index), None)
            if rs is None:
                continue
            chrpos = targets.get(rs_pub) or targets.get(rs_can)
            if chrpos is None:
                continue
            # Tier 1 → Tier 2 lookup
            hit, tier_tag = None, ""
            if chrpos in gcst196:
                hit = gcst196[chrpos]; tier_tag = "GCST013196_no-UKB_beta"
            elif chrpos in gcst197:
                hit = gcst197[chrpos]; tier_tag = "GCST013197_with-UKB_Z*SE"
            if hit is None:
                continue
            # harmonise β to the new BED's A1
            br = new_bim.loc[rs]
            if isinstance(br, pd.DataFrame):
                br = br.iloc[0]
            a1, a2 = str(br["A1"]).upper(), str(br["A2"]).upper()
            beta_gcst, flip = harmonise(hit["beta"], hit["EA"], a1, a2,
                                         is_or=False)
            if not math.isfinite(beta_gcst):
                print(f"  [{source}] {rs}: harmonise failed: {flip} "
                      f"(EA={hit['EA']} vs A1={a1}/A2={a2})")
                continue

            # cross-check vs legacy wmap (|Δβ| > 0.02 → WARNING)
            delta_log = ""
            legacy_beta = wmap_legacy.get(rs_pub) or wmap_legacy.get(rs_can)
            if legacy_beta is not None and math.isfinite(float(legacy_beta)):
                d_abs = abs(abs(beta_gcst) - abs(float(legacy_beta)))
                if d_abs > 0.02:
                    delta_log = (f"|Δβ|={d_abs:.4f} (GCST={beta_gcst:+.4f}, "
                                  f"wmap={float(legacy_beta):+.4f})")
                    print(f"  [{source}] WARN {rs}: β mismatch — {delta_log}")

            rows.append({
                "rsID": rs,
                "CHR": str(br["CHR"]),
                "BP_GRCh38": int(br["BP"]),
                "A1": a1, "A2": a2,
                "effect_allele_pub": str(r["effect_allele_pub"]).upper(),
                "beta_A1": float(beta_gcst),
                "source": "Wightman2021",
                "harmonisation": flip,
                "beta_source": tier_tag,
                "wmap_delta_beta": delta_log,
                "lead_rsID": rs_pub,
                "gene": "",
                "filter_exception": rs in BYPASS_FILTERS,
            })
        return pd.DataFrame(rows)

    # Bellenguez / Kunkle / Desikan: harmonise() based on OR_or_beta_pub
    is_or = SRC_IS_OR[source]
    src_label = {"Bellenguez": "Bellenguez2022",
                 "Desikan": "Desikan2019_PGS",
                 "Kunkle": "Kunkle2019"}[source]
    rows, drops = [], []
    for _, r in df.iterrows():
        rs_pub = str(r["rsID_pub"])
        rs_can = str(r["rsID_canonical"])
        rs = next((x for x in (rs_pub, rs_can) if x in new_bim.index), None)
        if rs is None:
            continue
        # APOE haplotype rows (Desikan) — skip (handled in APOE control)
        if str(r.get("is_haplotype", "")).upper() == "TRUE":
            continue
        try:
            val = float(r["OR_or_beta_pub"]) if r["OR_or_beta_pub"] else float("nan")
        except ValueError:
            val = float("nan")
        br = new_bim.loc[rs]
        if isinstance(br, pd.DataFrame):
            br = br.iloc[0]
        a1, a2 = str(br["A1"]).upper(), str(br["A2"]).upper()
        beta, flip = harmonise(val, str(r["effect_allele_pub"]), a1, a2, is_or)
        if not math.isfinite(beta):
            drops.append((rs, flip))
            continue
        rows.append({
            "rsID": rs, "CHR": str(br["CHR"]),
            "BP_GRCh38": int(br["BP"]),
            "A1": a1, "A2": a2,
            "effect_allele_pub": str(r["effect_allele_pub"]).upper(),
            "beta_A1": beta, "source": src_label,
            "harmonisation": flip,
            "lead_rsID": rs_pub,
            "gene": str(r.get("locus_name", "")),
            "filter_exception": rs in BYPASS_FILTERS,
        })
    if drops:
        print(f"  [{source}] {len(drops)} dropped during harmonisation:")
        for rs, why in drops[:5]:
            print(f"      {rs}: {why}")
    return pd.DataFrame(rows)


def _write_rescue_log(per_source: dict[str, pd.DataFrame],
                       new_bim_rs: set[str], out_root: Path) -> None:
    """For every PRS rsID across all sources, log whether it landed in the
    new BED or not, and (for the bypass list) flag it as an exception."""
    rows = []
    for src in SOURCES:
        df_res = pd.read_csv(SRC_PRS / f"{src}_full_snps_resolution.tsv",
                              sep="\t", dtype=str, keep_default_na=False)
        for _, r in df_res.iterrows():
            rs_pub = str(r["rsID_pub"])
            rs_can = str(r["rsID_canonical"])
            in_new = (rs_pub in new_bim_rs) or (rs_can in new_bim_rs)
            in_old = (str(r["drop_reason"]) == "ok")
            new_status = ("kept" if in_new else "still_dropped")
            if not in_old and in_new:
                outcome = "rescued"
            elif in_old and in_new:
                outcome = "carried_over"
            elif in_old and not in_new:
                outcome = "REGRESSED"           # shouldn't happen
            else:
                outcome = "still_dropped"
            rows.append({
                "source": src,
                "rsID_pub": rs_pub, "rsID_canonical": rs_can,
                "is_filter_exception": rs_pub in BYPASS_FILTERS or rs_can in BYPASS_FILTERS,
                "in_old_maf_bim": in_old,
                "in_new_geno005_bim": in_new,
                "outcome": outcome,
                "drop_reason_pre": str(r.get("drop_reason", "")),
            })
    out = pd.DataFrame(rows)
    out.to_csv(out_root / "rescue_log.tsv", sep="\t", index=False)
    n_rescued = (out["outcome"] == "rescued").sum()
    n_carried = (out["outcome"] == "carried_over").sum()
    print(f"  rescue_log: {n_carried} carried over, {n_rescued} rescued")
    if (out["outcome"] == "REGRESSED").any():
        regr = out[out["outcome"] == "REGRESSED"]
        print(f"  [WARN] {len(regr)} SNPs REGRESSED (in old MAF, not in new):")
        print(regr.to_string())


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--skip-build-bed", action="store_true",
                    help="Skip the plink rebuild — useful when iterating "
                         "on the export step only.")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    # ── Stage A: build the missingness-resolved BED ──────────────────────
    if not args.skip_build_bed or not (Path(str(NEW_BED_STEM) + ".bed")).exists():
        print("\n=== Stage A — build new BED (--geno 0.05 + bypass list) ===")
        print(f"[illumina] loading kgp→rsID manifests …")
        rsid2kgp = _load_illumina_rsid_to_kgp()
        _build_new_bed(args.out, rsid2kgp)
    else:
        print(f"\n[skip-build-bed] reusing existing {NEW_BED_STEM}.bed")

    # ── Load the new BIM for CHR/BP/A1/A2 lookups ────────────────────────
    new_bim = _bim_index(Path(str(NEW_BED_STEM) + ".bim"))
    new_bim_rs = set(new_bim.index)
    print(f"\nNew BED contains {len(new_bim_rs)} SNPs after all filters.")

    # ── Stage B: per-source harmonised tables (using new BIM coords) ─────
    print("\n=== Stage B — per-source β harmonisation against the new BED ===")
    bell = load_resolved_with_rescue("Bellenguez", new_bim)
    wight = load_resolved_with_rescue("Wightman", new_bim)
    desi = load_resolved_with_rescue("Desikan", new_bim)
    kunk = load_resolved_with_rescue("Kunkle", new_bim)
    print(f"  Bellenguez: {len(bell)} SNPs  (was 26 in MAF-resolved)")
    print(f"  Wightman  : {len(wight)} SNPs  (was 8)")
    print(f"  Desikan   : {len(desi)} SNPs  (was 15)")
    print(f"  Kunkle    : {len(kunk)} SNPs  (was 3)")

    # ── Stage C: dosage extraction for the union ─────────────────────────
    union_rs = sorted(set(bell["rsID"]) | set(wight["rsID"]) |
                       set(desi["rsID"]) | set(kunk["rsID"]))
    print(f"\n=== Stage C — extract dosage for union ({len(union_rs)} rsIDs) ===")
    dos = plink_extract_dosage(union_rs, NEW_BED_STEM, args.out / "_dosage_union")
    print(f"  dosage shape: {dos.shape}")
    if dos.shape[1] != len(union_rs):
        print(f"  [WARN] dosage cols={dos.shape[1]} != union={len(union_rs)}")

    # ── Stage D: per-source combos ───────────────────────────────────────
    print("\n=== Stage D — write per-source combos (4) ===")
    summary = []
    summary.append(write_combo(
        "sourceB", bell, dos, args.out,
        f"Full Bellenguez 2022 PRS — {len(bell)} missingness-resolved "
        f"(--geno 0.05) leads from the 89-SNP catalog."))
    summary.append(write_combo(
        "sourceW", wight, dos, args.out,
        f"Full Wightman 2021 PRS — {len(wight)} missingness-resolved "
        f"leads from the 38 no-UKBB filtered loci. Includes rs115186657 "
        f"as a published-PRS MAF-bypass exception."))
    summary.append(write_combo(
        "sourceD", desi, dos, args.out,
        f"Full Desikan 2017 PGS000026 — {len(desi)} missingness-resolved "
        f"leads. Includes rs4266886 as a published-PRS missingness-bypass "
        f"exception (F_MISS=0.12 → ~74 of 616 patients require cohort-mean "
        f"imputation during scoring)."))
    summary.append(write_combo(
        "sourceK", kunk, dos, args.out,
        f"Full Kunkle 2019 PRS — {len(kunk)} missingness-resolved leads "
        f"from the 24-SNP GWAS-Catalog GCST007511."))

    # ── Stage E: union combos ────────────────────────────────────────────
    print("\n=== Stage E — write union combos (4) ===")
    u_wb = dedupe_union([("Bellenguez", bell), ("Wightman", wight)])
    summary.append(write_combo(
        "5W26B", u_wb, dos, args.out,
        f"W ∪ B = {len(u_wb)} SNPs (missingness-resolved). B wins B∩W "
        f"overlaps."))

    u_wbd = dedupe_union([("Bellenguez", bell), ("Desikan", desi),
                           ("Wightman", wight)])
    summary.append(write_combo(
        "5W26B14D", u_wbd, dos, args.out,
        f"W ∪ B ∪ D = {len(u_wbd)} SNPs (missingness-resolved). "
        f"B wins B∩D rs6733839; no LD prune."))

    desi_pruned = desi[desi["rsID"] != "rs2597283"].copy()
    u_wbdp = dedupe_union([("Bellenguez", bell), ("Desikan", desi_pruned),
                            ("Wightman", wight)])
    summary.append(write_combo(
        "5W26B13D", u_wbdp, dos, args.out,
        f"W ∪ B ∪ D with intra-Desikan LD prune = {len(u_wbdp)} SNPs "
        f"(missingness-resolved)."))

    u_wbdk = dedupe_union([("Bellenguez", bell), ("Desikan", desi),
                            ("Wightman", wight), ("Kunkle", kunk)])
    summary.append(write_combo(
        "5W26B14D2K", u_wbdk, dos, args.out,
        f"W ∪ B ∪ D ∪ K = {len(u_wbdk)} SNPs (missingness-resolved). "
        f"D wins D∩K rs7920721."))

    # ── Stage F: top-level summary + rescue log ──────────────────────────
    pd.DataFrame(summary).to_csv(args.out / "combos_summary.tsv",
                                  sep="\t", index=False)
    _write_rescue_log({"Bellenguez": bell, "Wightman": wight,
                        "Desikan": desi, "Kunkle": kunk},
                       new_bim_rs, args.out)

    (args.out / "README.txt").write_text(
        "GWAS_filtered_missingness_resolved_geno005 — missingness-resolved "
        "(--geno 0.05) 4-source PRS combinations.\n\n"
        "Differs from GWAS_filtered_maf_resolved/ (which used --geno 0.02) "
        "in that the per-SNP missingness threshold is relaxed to 5%. Two "
        "published-PRS leads are additionally kept as filter-bypass "
        "exceptions:\n"
        "  rs4266886    Desikan PGS CR1 — bypasses --geno 0.05 "
        "(F_MISS=0.12)\n"
        "  rs115186657  Wightman      — bypasses --maf 0.01 (MAF=0.0037)\n\n"
        "The 2 monomorphic probes (rs75932628 TREM2 R47H, rs429358 APOE ε4) "
        "remain excluded — no variation in cohort.\n\n"
        f"Per-source counts:\n"
        f"  sourceB     {len(bell)} SNPs  Bellenguez 2022 (+2 vs MAF-resolved)\n"
        f"  sourceW     {len(wight)} SNPs  Wightman 2021 (+2 vs MAF-resolved)\n"
        f"  sourceD     {len(desi)} SNPs  Desikan 2017 PGS (+1 exception)\n"
        f"  sourceK     {len(kunk)} SNPs  Kunkle 2019\n\n"
        f"Union counts (winner's β kept on overlaps):\n"
        f"  5W26B        {len(u_wb)} SNPs  W ∪ B\n"
        f"  5W26B14D     {len(u_wbd)} SNPs  + D (no LD prune)\n"
        f"  5W26B13D     {len(u_wbdp)} SNPs  + intra-D LD prune\n"
        f"  5W26B14D2K   {len(u_wbdk)} SNPs  + K\n\n"
        "Each <COMBO>/ has <COMBO>_snps.tsv, <COMBO>_patient_dosage.tsv, "
        "<COMBO>_patient_dosage.bim, README.txt. combos_summary.tsv tallies "
        "all 8. rescue_log.tsv documents which SNPs were rescued / still "
        "dropped vs the MAF-resolved bundle.\n\n"
        "snps.tsv contains the same 12 columns as the MAF-resolved set; "
        "filter-exception rows are flagged via downstream loading code by "
        "rsID membership in {rs4266886, rs115186657}.\n",
        encoding="utf-8")

    print("\n=== Summary ===")
    for s in summary:
        print(f"  {s['combo']:14s}  n_snps={s['n_snps']:>3}  "
              f"B={s['n_bellenguez']}, W={s['n_wightman']}, "
              f"D={s['n_desikan']}, K={s['n_kunkle']}")
    print(f"\n[out] {args.out}")


if __name__ == "__main__":
    main()
