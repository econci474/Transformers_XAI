r"""
25d_export_maf_resolved_sets.py   (LOCAL, env: snp)
===================================================
Export the MAF>0.01-resolved PRS SNP combinations as self-contained,
training-ready bundles mirroring the script-25 layout
(`D:\…\GWAS_filtered\<COMBO>\`) but rebuilt against the **MAF-resolved
4-source per-SNP resolution tables** (script 30b output).

For each combo →
  D:\ADNI_SNP_Omni2.5M_20140220\GWAS_filtered_maf_resolved\<COMBO>\
    <COMBO>_snps.tsv            rsID, gene, lead_rsID, source, origin,
                                CHR, BP_GRCh38, effect_allele,
                                other_allele, REF, ALT, beta_A1
    <COMBO>_patient_dosage.tsv  PTID + one column per combo rsID
                                (additively coded 0/1/2 = A1=ALT dosage)
    <COMBO>_patient_dosage.bim  plink-style rsID CHR BP A1 A2
    README.txt

Combos (8 total) — naming follows the user's dedup convention
(component count = "after rsID-dedup against previously-named sources";
LD redundancies at distinct rsIDs are NOT deducted unless intra-source).
Hierarchy: B>D>W>K (larger source wins shared rsIDs; matches user math):

  Per-source PRS — each source uses all its own genotyped SNPs / its own β:
    sourceB        26 SNPs  Bellenguez 2022 (MAF>0.01 resolved)
    sourceW         8 SNPs  Wightman 2021
    sourceD        15 SNPs  Desikan 2017 PGS000026
    sourceK         3 SNPs  Kunkle 2019 catalog

  Union PRS — winner's β kept on overlaps (B wins B∩W,B∩D; D wins D∩K):
    5W26B          31 SNPs  W ∪ B
    5W26B14D       45 SNPs  + D (no LD prune; both PICALM rsIDs kept)
    5W26B13D       44 SNPs  + intra-Desikan LD prune (drop rs2597283)
    5W26B14D2K     47 SNPs  + K

β harmonisation:
  - Bellenguez/Kunkle: published OR  → β = ln(OR), then flip to A1.
  - Desikan: published β/effect_weight → flip to A1.
  - Wightman: β loaded from `wightman_resolved_allowlist.tsv` (5 in-430
    SNPs, already A1-aligned via script 24) + `wfilt_enriched_weights.tsv`
    (3 re-extracted, `beta_A1` already harmonised).

Dosage: re-extracted via PLINK from the canonical MAF BED
(`SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr_maf001`) into
a single 47-SNP × 616-patient matrix, then subset per combo.

Usage:  conda run -n snp python snp_pipeline/25d_export_maf_resolved_sets.py
"""
from __future__ import annotations

import argparse
import math
import subprocess
from pathlib import Path

import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
MAF_BED = BASE / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr_maf001"
SRC_PRS = BASE / "source_prs"
GBF = BASE / "genetic_baselines_filtered"
PLINK = BASE / "plink.exe"

# Wightman β sources (allowlist for in-430 + wfilt_enriched_weights for re-extracted)
WIGHT_ALLOWLIST = GBF / "wightman_resolved_allowlist.tsv"
WFILT_ENR = GBF / "wfilt_enriched_weights.tsv"

OUT_DEFAULT = BASE / "GWAS_filtered_maf_resolved"

# Sources where OR_or_beta_pub is an OR (need ln() conversion)
SRC_IS_OR = {"Bellenguez": True, "Kunkle": True,
              "Desikan": False, "Wightman": False}

SNP_COLS = ["rsID", "gene", "lead_rsID", "source", "origin", "CHR",
            "BP_GRCh38", "effect_allele", "other_allele", "REF", "ALT",
            "beta_A1"]

# Winner of shared rsIDs (B > D > W > K)
PRIORITY = {"Bellenguez": 0, "Desikan": 1, "Wightman": 2, "Kunkle": 3}

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def harmonise(or_or_beta: float, effect_allele: str, a1: str, a2: str,
              is_or: bool) -> tuple[float, str]:
    """Return (beta_A1, flip_label) or (nan, reason)."""
    if not (or_or_beta is not None and math.isfinite(or_or_beta)):
        return float("nan"), "missing_or_or_beta"
    if is_or:
        if or_or_beta <= 0:
            return float("nan"), "non_positive_OR"
        b = math.log(or_or_beta)
    else:
        b = float(or_or_beta)
    ea = effect_allele.upper().strip()
    a1, a2 = a1.upper(), a2.upper()
    if ea == a1:
        return b, "direct"
    if ea == a2:
        return -b, "flip"
    if COMP.get(ea) == a1:
        return b, "strand"
    if COMP.get(ea) == a2:
        return -b, "strand_flip"
    return float("nan"), f"allele_mismatch_EA={ea}_A1={a1}_A2={a2}"


def load_resolved(source: str) -> pd.DataFrame:
    """Load a source's resolution TSV, filter to drop_reason='ok',
    harmonise β. For Wightman, β comes from the allowlist/wfilt files."""
    p = SRC_PRS / f"{source}_full_snps_resolution.tsv"
    df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    df = df[df["drop_reason"] == "ok"].copy()
    df["BP_GRCh38"] = pd.to_numeric(df["BP_GRCh38"], errors="coerce").astype("Int64")
    df["A1"] = df["A1_in_bim"].str.upper()
    df["A2"] = df["A2_in_bim"].str.upper()

    if source == "Wightman":
        allow = pd.read_csv(WIGHT_ALLOWLIST, sep="\t")
        wmap = dict(zip(allow["pipeline_rsID"].astype(str),
                         allow["beta"].astype(float)))
        wfilt = pd.read_csv(WFILT_ENR, sep="\t")
        wfilt = wfilt[(wfilt["novel"] == True) &                 # noqa: E712
                      wfilt["source"].astype(str).str.contains("Wightman")]
        wmap.update(dict(zip(wfilt["rsID"].astype(str),
                              wfilt["beta_A1"].astype(float))))
        rows = []
        for _, r in df.iterrows():
            rs = str(r["rsID_in_bim"])
            beta = wmap.get(rs)
            if beta is None or not math.isfinite(beta):
                continue
            rows.append({
                "rsID": rs, "CHR": str(r["CHR_GRCh38"]),
                "BP_GRCh38": int(r["BP_GRCh38"]),
                "A1": r["A1"], "A2": r["A2"],
                "effect_allele_pub": r["effect_allele_pub"].upper(),
                "beta_A1": beta, "source": "Wightman2021",
                "harmonisation": "wightman_pre_aligned",
                "lead_rsID": str(r["rsID_pub"]),
                "gene": ""})
        return pd.DataFrame(rows)

    is_or = SRC_IS_OR[source]
    src_label = {"Bellenguez": "Bellenguez2022", "Desikan": "Desikan2019_PGS",
                  "Kunkle": "Kunkle2019"}[source]
    rows = []
    drops = []
    for _, r in df.iterrows():
        try:
            val = float(r["OR_or_beta_pub"]) if r["OR_or_beta_pub"] else float("nan")
        except ValueError:
            val = float("nan")
        beta, flip = harmonise(val, r["effect_allele_pub"], r["A1"], r["A2"],
                                is_or)
        if not math.isfinite(beta):
            drops.append((r["rsID_in_bim"], flip))
            continue
        rows.append({
            "rsID": str(r["rsID_in_bim"]), "CHR": str(r["CHR_GRCh38"]),
            "BP_GRCh38": int(r["BP_GRCh38"]),
            "A1": r["A1"], "A2": r["A2"],
            "effect_allele_pub": r["effect_allele_pub"].upper(),
            "beta_A1": beta, "source": src_label,
            "harmonisation": flip,
            "lead_rsID": str(r["rsID_pub"]),
            "gene": str(r.get("locus_name", ""))})
    if drops:
        print(f"  [{source}] dropped {len(drops)} during harmonisation:")
        for rs, why in drops[:5]:
            print(f"      {rs}: {why}")
    return pd.DataFrame(rows)


def plink_extract_dosage(rsids: list[str], outdir: Path) -> pd.DataFrame:
    """plink --bfile MAF_BED --extract <rsids> --recode A → return wide
    dosage DataFrame indexed by PTID with 0/1/2 A1 counts."""
    outdir.mkdir(parents=True, exist_ok=True)
    snplist = outdir / "rsid_extract.txt"
    snplist.write_text("\n".join(rsids) + "\n", encoding="utf-8")
    stem = outdir / "patient_dosage"
    cmd = [str(PLINK), "--bfile", str(MAF_BED),
           "--extract", str(snplist),
           "--keep-allele-order",
           "--recode", "A",
           "--out", str(stem)]
    print(f"  [plink] {' '.join(cmd[1:])}")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"plink failed:\n{r.stderr}")
    raw = pd.read_csv(stem.with_suffix(".raw"), sep=r"\s+", engine="python")
    pid = raw["IID"].astype(str)
    geno = raw.drop(columns=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"])
    # Column names in PLINK .raw are `rsID_<A1>`; strip to rsID
    geno.columns = [c.rsplit("_", 1)[0] for c in geno.columns]
    geno.index = pid
    geno.index.name = "PTID"
    return geno


def dedupe_union(parts: list[tuple[str, pd.DataFrame]]) -> pd.DataFrame:
    """Concat per-source frames; winner = lowest PRIORITY value (B<D<W<K)."""
    frames = []
    for src, df in parts:
        df = df.copy()
        df["__src"] = src
        df["__prio"] = PRIORITY[src]
        frames.append(df)
    all_df = pd.concat(frames, ignore_index=True)
    all_df = all_df.sort_values(["rsID", "__prio"], kind="stable")
    deduped = all_df.drop_duplicates("rsID", keep="first")
    return deduped.drop(columns=["__src", "__prio"]).reset_index(drop=True)


def write_combo(combo: str, snp_df: pd.DataFrame, dos: pd.DataFrame,
                out_root: Path, description: str) -> dict:
    d = out_root / combo
    d.mkdir(parents=True, exist_ok=True)
    cols = [c for c in snp_df["rsID"] if c in dos.columns]
    miss = sorted(set(snp_df["rsID"]) - set(cols))
    if miss:
        print(f"  [WARN] {combo}: {len(miss)} SNPs not in dosage matrix: {miss}")
    # snps.tsv (sorted CHR/BP, with REF/ALT and full provenance)
    rows = []
    for s in cols:
        r = snp_df[snp_df["rsID"] == s].iloc[0]
        rows.append({"rsID": s, "gene": r.get("gene", ""),
                      "lead_rsID": r.get("lead_rsID", s),
                      "source": r["source"],
                      "origin": f"{r['source']}_resolved",
                      "CHR": str(r["CHR"]),
                      "BP_GRCh38": int(r["BP_GRCh38"]),
                      "effect_allele": r["A1"],
                      "other_allele": r["A2"],
                      "REF": r["A2"], "ALT": r["A1"],
                      "beta_A1": float(r["beta_A1"])})
    snp_out = (pd.DataFrame(rows, columns=SNP_COLS)
                 .sort_values(["CHR", "BP_GRCh38"],
                               key=lambda c: c.map(
                                   lambda v: (int(v) if str(v).isdigit() else 99)
                                   if c.name == "CHR" else v))
                 .reset_index(drop=True))
    snp_out.to_csv(d / f"{combo}_snps.tsv", sep="\t", index=False)

    dsub = dos[cols].copy().reset_index()
    dsub.to_csv(d / f"{combo}_patient_dosage.tsv", sep="\t", index=False)

    # .bim — A1 = effect = ALT
    bim = snp_out[["rsID", "CHR", "BP_GRCh38", "effect_allele",
                   "other_allele"]].rename(
                       columns={"BP_GRCh38": "BP",
                                 "effect_allele": "A1",
                                 "other_allele": "A2"})
    bim.to_csv(d / f"{combo}_patient_dosage.bim", sep="\t", index=False)

    nB = int((snp_out["source"] == "Bellenguez2022").sum())
    nW = int((snp_out["source"] == "Wightman2021").sum())
    nD = int((snp_out["source"] == "Desikan2019_PGS").sum())
    nK = int((snp_out["source"] == "Kunkle2019").sum())
    (d / "README.txt").write_text(
        f"{combo} — MAF>0.01-resolved PRS set.\n\n"
        f"{description}\n\n"
        f"SNP count: {len(snp_out)}  "
        f"(Bellenguez={nB}, Wightman={nW}, Desikan={nD}, Kunkle={nK})\n"
        f"Patients : {len(dsub)}\n\n"
        f"{combo}_snps.tsv columns:\n"
        "  rsID          dbSNP id (GRCh38; MAF-resolved via 30b)\n"
        "  gene          mapped gene (catalog/PGS locus name)\n"
        "  lead_rsID     original published rsID (= rsID before "
        "RsMergeArch canonicalisation)\n"
        "  source        Bellenguez2022 | Wightman2021 | "
        "Desikan2019_PGS | Kunkle2019\n"
        "  origin        <source>_resolved\n"
        "  CHR,BP_GRCh38 chromosome and GRCh38 base-pair position\n"
        "  effect_allele = A1 = ALT (the allele the dosage counts)\n"
        "  other_allele  = A2 = REF\n"
        "  REF/ALT       as in refcorr GRCh38 BIM\n"
        "  beta_A1       harmonised effect size on A1 "
        "(Bellenguez/Kunkle: ln(OR); Desikan/Wightman: published β)\n\n"
        f"{combo}_patient_dosage.tsv: PTID + one column per rsID, "
        "0/1/2 = count of effect_allele (A1=ALT). Missing → blank "
        "(downstream models impute with cohort mean of the training "
        "fold).\n\n"
        f"{combo}_patient_dosage.bim: plink-style (rsID CHR BP A1 A2); "
        "A1 = counted/effect/ALT, A2 = other/REF.\n\n"
        "Weighted PRS:  PRS = Σ_i beta_A1_i × dosage_i  over the "
        "rsIDs in *_snps.tsv (cohort-mean impute missing dosage, then "
        "z-score on the training fold).\n",
        encoding="utf-8")
    print(f"  [out] {d}  ({len(cols)} SNPs: B={nB} W={nW} D={nD} K={nK}; "
          f"{len(dsub)} patients)")
    return {"combo": combo, "n_snps": len(cols), "n_bellenguez": nB,
            "n_wightman": nW, "n_desikan": nD, "n_kunkle": nK,
            "n_patients": len(dsub)}


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    print("=== Load resolved per-source β tables ===")
    bell = load_resolved("Bellenguez")
    wight = load_resolved("Wightman")
    desi = load_resolved("Desikan")
    kunk = load_resolved("Kunkle")
    print(f"  Bellenguez: {len(bell)} SNPs")
    print(f"  Wightman  : {len(wight)} SNPs")
    print(f"  Desikan   : {len(desi)} SNPs")
    print(f"  Kunkle    : {len(kunk)} SNPs")

    # Sanity: counts should match 26 / 8 / 15 / 3
    assert len(bell) == 26, f"Bellenguez = {len(bell)} (expected 26)"
    assert len(wight) == 8,  f"Wightman = {len(wight)} (expected 8)"
    assert len(desi) == 15, f"Desikan = {len(desi)} (expected 15)"
    assert len(kunk) == 3,  f"Kunkle = {len(kunk)} (expected 3)"

    # Union rsIDs across all 4 sources for the plink extract (47)
    union_rs = sorted(set(bell["rsID"]) | set(wight["rsID"]) |
                       set(desi["rsID"]) | set(kunk["rsID"]))
    print(f"\n=== Union {len(union_rs)} rsIDs → plink extract from MAF BED ===")
    dos = plink_extract_dosage(union_rs, args.out / "_dosage_union")
    print(f"  dosage shape: {dos.shape}")
    assert dos.shape[1] == len(union_rs), \
        f"plink returned {dos.shape[1]} cols vs expected {len(union_rs)}"

    print("\n=== Write per-source combos (4) ===")
    summary = []
    summary.append(write_combo(
        "sourceB", bell, dos, args.out,
        "Full Bellenguez 2022 PRS — 26 MAF>0.01-resolved leads from the "
        "89-SNP catalog (Belloy et al. 2022 supplementary)."))
    summary.append(write_combo(
        "sourceW", wight, dos, args.out,
        "Full Wightman 2021 PRS — 8 MAF>0.01-resolved leads from the 38 "
        "no-UKBB filtered loci (5 in-430 + 3 re-extracted via plink)."))
    summary.append(write_combo(
        "sourceD", desi, dos, args.out,
        "Full Desikan 2017 PGS000026 — 15 MAF>0.01-resolved leads "
        "(2 APOE haplotype rows ε4/ε2 are excluded — handled in the "
        "APOE control arm)."))
    summary.append(write_combo(
        "sourceK", kunk, dos, args.out,
        "Full Kunkle 2019 PRS — 3 MAF>0.01-resolved leads from the "
        "24-SNP GWAS-Catalog GCST007511 (the other 21 are below MAF "
        "0.01 or absent from Omni2.5M)."))

    print("\n=== Write union combos (4) ===")
    # 5W26B = W ∪ B  (B wins overlaps)
    u_wb = dedupe_union([("Bellenguez", bell), ("Wightman", wight)])
    summary.append(write_combo(
        "5W26B", u_wb, dos, args.out,
        "W ∪ B (Wightman + Bellenguez). 31 unique SNPs (26 Bellenguez + "
        "5 Wightman after removing 3 B∩W rsID overlaps; B retains "
        "rs871269, rs1582763, rs7912495).  ≈ legacy 9W27B set, 1 fewer "
        "SNP due to MAF>0.01 filter dropping 1 Bellenguez variant."))

    # 5W26B14D = W ∪ B ∪ D
    u_wbd = dedupe_union([("Bellenguez", bell), ("Desikan", desi),
                           ("Wightman", wight)])
    summary.append(write_combo(
        "5W26B14D", u_wbd, dos, args.out,
        "W ∪ B ∪ D (Wightman + Bellenguez + Desikan). 45 unique SNPs "
        "(B wins B∩D rs6733839; both PICALM rsIDs kept — Wightman "
        "rs561655 and Desikan rs543293, r²=0.84 LD-redundant but "
        "distinct rsIDs).  No LD-prune."))

    # 5W26B13D = remove intra-Desikan LD pair (drop rs2597283, keep rs17265593)
    desi_pruned = desi[desi["rsID"] != "rs2597283"].copy()
    u_wbdp = dedupe_union([("Bellenguez", bell), ("Desikan", desi_pruned),
                            ("Wightman", wight)])
    summary.append(write_combo(
        "5W26B13D", u_wbdp, dos, args.out,
        "W ∪ B ∪ D with intra-Desikan LD prune (drop rs2597283; keep "
        "rs17265593 as the labelled BC043356 lead). 44 unique SNPs."))

    # 5W26B14D2K = W ∪ B ∪ D ∪ K  (D wins D∩K)
    u_wbdk = dedupe_union([("Bellenguez", bell), ("Desikan", desi),
                            ("Wightman", wight), ("Kunkle", kunk)])
    summary.append(write_combo(
        "5W26B14D2K", u_wbdk, dos, args.out,
        "W ∪ B ∪ D ∪ K (all 4 sources). 47 unique SNPs.  D wins D∩K "
        "rs7920721; B wins B∩W (3) and B∩D (1). No LD-prune."))

    # Top-level summary
    pd.DataFrame(summary).to_csv(args.out / "combos_summary.tsv",
                                  sep="\t", index=False)
    (args.out / "README.txt").write_text(
        "GWAS_filtered_maf_resolved — MAF>0.01-resolved 4-source PRS "
        "combinations.\n\n"
        "8 combos (each its own sub-directory), reconstructed from the "
        "per-source resolution tables at D:\\…\\source_prs\\ (script 30b "
        "output, MAF>0.01 BIM resolution) and dosage extracted from the "
        "canonical MAF BED via plink.\n\n"
        "Per-source PRS (use all native genotyped SNPs per source):\n"
        "  sourceB        26 SNPs  Bellenguez 2022\n"
        "  sourceW         8 SNPs  Wightman 2021\n"
        "  sourceD        15 SNPs  Desikan 2017 PGS000026\n"
        "  sourceK         3 SNPs  Kunkle 2019\n\n"
        "Union PRS (winner's β kept on shared rsIDs; user's dedup naming):\n"
        "  5W26B          31 SNPs  W ∪ B\n"
        "  5W26B14D       45 SNPs  + D (no LD-prune)\n"
        "  5W26B13D       44 SNPs  + intra-Desikan LD-prune\n"
        "  5W26B14D2K     47 SNPs  + K\n\n"
        "Each <COMBO>/ has <COMBO>_snps.tsv, <COMBO>_patient_dosage.tsv, "
        "<COMBO>_patient_dosage.bim, and README.txt. combos_summary.tsv "
        "tallies all 8.\n",
        encoding="utf-8")
    print(f"\n[out] {args.out / 'combos_summary.tsv'}")
    print("\n=== combos summary ===")
    for s in summary:
        print(f"  {s['combo']:<14}  n_snps={s['n_snps']:>3}  "
              f"B={s['n_bellenguez']} W={s['n_wightman']} "
              f"D={s['n_desikan']} K={s['n_kunkle']}  "
              f"n_pat={s['n_patients']}")


if __name__ == "__main__":
    main()
