r"""
25c_export_source_specific_prs.py   (LOCAL, env: snp)
====================================================
Build per-source PRS for each of the 4 GWAS studies using **all native
genotyped SNPs of that source** (no inter-source de-duplication — SNPs
overlapping multiple studies contribute to each PRS with their own
study-specific β). Source SNP sets:

  Bellenguez    27 leads  (from external_gwas_labels.tsv source=='Bellenguez2022')
  Wightman      9 leads   (script-24 resolved: 6 in-430 + 3 re-extracted from 38 published)
  Kunkle        ~3-10     (24 GWAS-Catalog leads → genotyped subset on MAF>0.01 array;
                          tried rsID + positional + RsMergeArch matching)
  Desikan       ~10       (desikan_pgs_resolved.tsv, in_430 subset)

For each patient, computes one PRS per source:
    PRS_<src>(p) = Σ_i  beta_A1_i · dosage_{p,i}   over source's genotyped SNPs.
Missing dosage imputed with cohort mean (per-SNP), matching script 25 / 20.

Output → D:\…\source_prs\:
  source_snps.tsv         per-SNP table (rsID, source, CHR, BP, A1, A2,
                          beta_A1, gene, in_each_source_flag)
  patient_source_prs.tsv  PTID + PRS_Bellenguez + PRS_Wightman + PRS_Kunkle
                          + PRS_Desikan (raw, NOT z-scored)
  source_summary.txt      per-source n_SNPs, overlap counts, Kunkle drop log

Usage:
  conda run -n snp python snp_pipeline/25c_export_source_specific_prs.py
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util as _il
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

# ── reuse script 20 (build_beta_table needs DOSAGE_TSV/_BIM/etc.) ────────
_S20 = Path(__file__).parent / "20_genetic_baselines.py"
_spec = _il.spec_from_file_location("s20", _S20)
_s20 = _il.module_from_spec(_spec)
_spec.loader.exec_module(_s20)
build_beta_table = _s20.build_beta_table
_load_wightman_allowlist = _s20._load_wightman_allowlist

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
MAF_BED = BASE / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr_maf001"
GBF = BASE / "genetic_baselines_filtered"
GB = BASE / "genetic_baselines"
LABELS = BASE / "bmfm_inputs" / "external_gwas_labels.tsv"

# Kunkle GWAS-Catalog (24 leads — accession GCST007511, GRCh38)
KUNKLE_CAT = (BASE / "GWAS" / "Kunkle_2019" /
              "Kunkle_2019_gwas-association-downloaded_2026-05-18-"
              "accessionId_GCST007511.tsv")
KUNKLE_STAGE1 = (BASE / "GWAS" / "Kunkle_2019" /
                 "NG00075_Kunkle_etal_Stage1_P-val_only_results.txt")
RSMERGE = BASE / "liftover" / "NCBI" / "RsMergeArch.bcp.gz"
CHAIN_HG19_TO_HG38 = BASE / "liftover" / "hg19ToHg38.over.chain.gz"

OUT = BASE / "source_prs"


def _rs_canonical(rsids: set[str]) -> dict[str, str]:
    """RsMergeArch low→high merge map for the given rsIDs (only)."""
    out = {}
    needles = {r.lstrip("rs").lstrip("0") for r in rsids if r.startswith("rs")}
    if not RSMERGE.exists():
        return out
    with gzip.open(RSMERGE, "rt") as f:
        for line in f:
            t = line.rstrip("\n").split("\t")
            if len(t) < 3:
                continue
            lo, hi = t[0], t[1]
            if lo in needles:
                out[f"rs{lo}"] = f"rs{hi}"
    return out


def parse_kunkle_catalog(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    # column-name lookup robust to spaces / hyphens in GWAS-Catalog headers
    cols = {c.upper().strip(): c for c in df.columns}

    def col(*aliases):
        for a in aliases:
            k = a.upper().strip()
            if k in cols:
                return cols[k]
        return None
    c_rs = col("SNPS")
    c_chr = col("CHR_ID", "CHR")
    c_bp = col("CHR_POS", "BP")
    c_risk = col("STRONGEST SNP-RISK ALLELE", "STRONGEST_SNP_RISK_ALLELE")
    c_or = col("OR or BETA", "OR_or_BETA", "OR")
    c_p = col("P-VALUE", "P_VALUE", "PVAL", "P")
    if c_rs is None:
        raise ValueError(f"No SNPS column in {path.name}")
    keep = []
    for _, r in df.iterrows():
        rs = (str(r[c_rs]) if c_rs else "").strip()
        if not rs.startswith("rs"):
            continue
        chrom = (str(r[c_chr]) if c_chr else "").strip()
        bp = (str(r[c_bp]) if c_bp else "").strip()
        risk = (str(r[c_risk]) if c_risk else "").strip()
        ea = risk.split("-")[-1].upper() if "-" in risk else ""
        if ea in {"", "?"}:
            ea = ""
        or_str = (str(r[c_or]) if c_or else "").strip()
        pval = (str(r[c_p]) if c_p else "").strip()
        try:
            or_v = float(or_str)
        except ValueError:
            or_v = float("nan")
        try:
            p_v = float(pval)
        except ValueError:
            p_v = float("nan")
        keep.append({"rsID_pub": rs, "chrom_pub": chrom, "bp_pub": bp,
                     "effect_allele_pub": ea, "OR": or_v, "P": p_v})
    return pd.DataFrame(keep)


def resolve_kunkle(maf_bim: Path) -> pd.DataFrame:
    """Resolve Kunkle catalog leads onto the MAF-filtered BIM via:
      1. direct rsID match
      2. RsMergeArch canonicalisation (low→high)
      3. positional match (CHR:BP exact, GRCh38)
      4. **hg19→hg38 liftover** then positional (handles GRCh37 coords —
         Kunkle 2019 paper is in GRCh37 and the GWAS-Catalog row's
         CHR_POS may carry the original paper position).
    Returns one row per resolved + harmonised lead with beta_A1=ln(OR)
    harmonised to the BIM A1."""
    cat = parse_kunkle_catalog(KUNKLE_CAT)
    print(f"[kunkle] catalog leads: {len(cat)}")

    # need EA but some rows have '?' — drop those (can't harmonise)
    n_no_ea = int((cat["effect_allele_pub"] == "").sum())
    print(f"[kunkle] missing effect_allele in catalog: {n_no_ea}")
    rsmap = _rs_canonical(set(cat["rsID_pub"]))
    cat["rsID_canonical"] = cat["rsID_pub"].map(
        lambda r: rsmap.get(r, r))
    if rsmap:
        print(f"[kunkle] RsMergeArch hits: {len(rsmap)} "
              f"(e.g. {list(rsmap.items())[:3]})")

    # hg19→hg38 liftover (Kunkle paper is GRCh37)
    lo = None
    try:
        from pyliftover import LiftOver
        if CHAIN_HG19_TO_HG38.exists():
            lo = LiftOver(str(CHAIN_HG19_TO_HG38))
            print(f"[kunkle] hg19→hg38 chain loaded from "
                  f"{CHAIN_HG19_TO_HG38.name}")
    except ImportError:
        print("[kunkle] pyliftover not available — skipping hg19→hg38")

    def _lift(chrom: str, bp_str: str) -> int | None:
        if lo is None or not chrom or not bp_str:
            return None
        try:
            bp = int(bp_str) - 1                     # liftover wants 0-based
        except ValueError:
            return None
        hits = lo.convert_coordinate(f"chr{chrom}", bp)
        if not hits:
            return None
        return int(hits[0][1]) + 1                   # 1-based back

    bim = pd.read_csv(maf_bim, sep="\t", header=None,
                       names=["CHR", "rsID", "cM", "BP", "A1", "A2"],
                       dtype=str)
    bim["BP"] = pd.to_numeric(bim["BP"], errors="coerce").astype("Int64")
    bim_by_rs = bim.set_index("rsID")
    bim_by_pos = bim.set_index(bim["CHR"].astype(str) + ":" + bim["BP"].astype(str))

    rows, drops = [], []
    for r in cat.itertuples(index=False):
        rs_pub = r.rsID_pub
        rs_can = r.rsID_canonical
        chrom = r.chrom_pub
        bp = r.bp_pub
        ea = r.effect_allele_pub.upper().strip()
        method = None
        bim_row = None
        if rs_pub in bim_by_rs.index:
            bim_row = bim_by_rs.loc[rs_pub]
            method = "rsID"
        elif rs_can in bim_by_rs.index:
            bim_row = bim_by_rs.loc[rs_can]
            method = "rsID_merged"
        elif chrom and bp:
            # try (a) catalog BP as GRCh38 first
            key = f"{chrom}:{int(bp)}" if bp.lstrip("-").isdigit() else None
            if key and key in bim_by_pos.index:
                bim_row = bim_by_pos.loc[key]
                if isinstance(bim_row, pd.DataFrame):
                    bim_row = bim_row.iloc[0]
                method = "position_grch38"
            else:
                # (b) treat catalog BP as GRCh37, liftover to GRCh38
                bp38 = _lift(chrom, bp)
                if bp38 is not None:
                    key2 = f"{chrom}:{int(bp38)}"
                    if key2 in bim_by_pos.index:
                        bim_row = bim_by_pos.loc[key2]
                        if isinstance(bim_row, pd.DataFrame):
                            bim_row = bim_row.iloc[0]
                        method = "position_lifted_grch37_to_38"
        if bim_row is None:
            drops.append({"rsID_pub": rs_pub, "chrom_pub": chrom,
                          "bp_pub": bp, "reason": "not_in_MAF_BIM"})
            continue
        if not ea:
            drops.append({"rsID_pub": rs_pub, "chrom_pub": chrom,
                          "bp_pub": bp, "reason": "missing_effect_allele"})
            continue
        if not (r.OR > 0 and math.isfinite(r.OR)):
            drops.append({"rsID_pub": rs_pub, "chrom_pub": chrom,
                          "bp_pub": bp, "reason": "missing_OR"})
            continue
        a1, a2 = str(bim_row["A1"]).upper(), str(bim_row["A2"]).upper()
        beta_raw = math.log(r.OR)
        comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
        if ea == a1:
            beta_a1 = beta_raw
            flip = False
        elif ea == a2:
            beta_a1 = -beta_raw
            flip = True
        elif comp.get(ea) == a1:
            beta_a1 = beta_raw
            flip = "strand"
        elif comp.get(ea) == a2:
            beta_a1 = -beta_raw
            flip = "strand"
        else:
            drops.append({"rsID_pub": rs_pub, "chrom_pub": chrom,
                          "bp_pub": bp,
                          "reason": f"allele_mismatch_EA={ea}_A1={a1}_A2={a2}"})
            continue
        rows.append({"rsID": str(bim_row.name),
                     "rsID_pub": rs_pub,
                     "CHR": str(bim_row["CHR"]),
                     "BP_GRCh38": int(bim_row["BP"]),
                     "A1": a1, "A2": a2,
                     "effect_allele_pub": ea,
                     "beta_A1": beta_a1, "OR": r.OR, "P": r.P,
                     "match_method": method, "allele_flip": flip,
                     "source": "Kunkle2019"})
    res = pd.DataFrame(rows)
    print(f"[kunkle] resolved+harmonised: {len(res)} / {len(cat)} "
          f"(methods: {res['match_method'].value_counts().to_dict() if len(res) else '{}'})")
    print(f"[kunkle] dropped: {len(drops)}")
    return res, pd.DataFrame(drops)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    # ── 1. Bellenguez (27) — from build_beta_table (no Wightman filter) ──
    print("\n[bellenguez] building β-table (full unfiltered) …")
    bt = build_beta_table()
    bell = bt[bt["source"].astype(str).str.contains("Bellenguez")].copy()
    bell = bell.reset_index().rename(columns={"SNP": "rsID"})
    bell["source"] = "Bellenguez2022"
    bell["match_method"] = "external_gwas_labels"
    bell["BP_GRCh38"] = bell["BP"].astype(int)
    bell = bell[["rsID", "CHR", "BP_GRCh38", "A1", "A2",
                 "beta_A1", "source", "match_method"]]
    print(f"[bellenguez] {len(bell)} SNPs")

    # ── 2. Wightman (9 of 38) — allowlist β directly + 3 re-extracted ────
    # `wightman_resolved_allowlist.tsv` is the script-24 output: 6 in-430
    # Wightman leads with β (the allowlist's own β column — not all are
    # picked up by build_beta_table's trusted-source matching, hence we
    # load directly). + 3 re-extracted from wfilt_enriched_weights (β_A1
    # harmonised). MAF-BIM gives A1/A2/CHR/BP_GRCh38 for harmonisation.
    print("\n[wightman] resolving (allowlist β direct + re-extracted) …")
    allow = pd.read_csv(GBF / "wightman_resolved_allowlist.tsv", sep="\t")
    bim_full = pd.read_csv(MAF_BED.with_suffix(".bim"), sep="\t",
                           header=None, names=["CHR", "rsID", "cM", "BP",
                                                "A1", "A2"], dtype=str)
    bim_full["BP"] = pd.to_numeric(bim_full["BP"], errors="coerce").astype("Int64")
    bim_by_rs = bim_full.set_index("rsID")
    w_rows = []
    for _, r in allow.iterrows():
        rs = str(r["pipeline_rsID"]).strip()
        if rs not in bim_by_rs.index:
            print(f"  [skip] {rs}: not in MAF BIM")
            continue
        br = bim_by_rs.loc[rs]
        # allowlist β is from Wightman raw; we trust it as A1-oriented for
        # the array's A1 (the allowlist was built strand-aware in script 24)
        beta = float(r["beta"]) if pd.notna(r["beta"]) else float("nan")
        if not math.isfinite(beta):
            print(f"  [skip] {rs}: no beta in allowlist")
            continue
        w_rows.append({"rsID": rs, "CHR": str(br["CHR"]),
                       "BP_GRCh38": int(br["BP"]),
                       "A1": str(br["A1"]).upper(),
                       "A2": str(br["A2"]).upper(),
                       "beta_A1": beta, "source": "Wightman2021",
                       "match_method": "wightman_allowlist_in_430"})
    wfilt = pd.DataFrame(w_rows)
    re_ex = pd.read_csv(GBF / "wfilt_enriched_weights.tsv", sep="\t")
    re_ex = re_ex[(re_ex["novel"] == True)
                   & re_ex["source"].astype(str).str.contains("Wightman")]
    re_ex = re_ex.rename(columns={"BP": "BP_GRCh38"})
    re_ex["match_method"] = "re_extracted_plink"
    re_ex = re_ex[["rsID", "CHR", "BP_GRCh38", "A1", "A2",
                   "beta_A1", "source", "match_method"]]
    wight = pd.concat([wfilt, re_ex], ignore_index=True).drop_duplicates("rsID")
    wight["source"] = "Wightman2021"
    print(f"[wightman] {len(wight)} SNPs "
          f"(allowlist in-430 {len(wfilt)} + re-extracted {len(re_ex)})")

    # ── 3. Kunkle (catalog 24 → genotyped subset) ────────────────────────
    kun, kun_drops = resolve_kunkle(MAF_BED.with_suffix(".bim"))
    if len(kun):
        kun = kun.rename(columns={"rsID_pub": "_pub"})
        kun_min = kun[["rsID", "CHR", "BP_GRCh38", "A1", "A2",
                       "beta_A1", "source", "match_method"]].copy()
    else:
        kun_min = pd.DataFrame(columns=["rsID", "CHR", "BP_GRCh38",
                                         "A1", "A2", "beta_A1", "source",
                                         "match_method"])

    # ── 4. Desikan — from desikan_pgs_resolved.tsv (script-21 output) ────
    des_path = GB / "desikan_pgs_resolved.tsv"
    if des_path.exists():
        des = pd.read_csv(des_path, sep="\t")
        des = des[des["beta_A1"].notna() & (des["in_430"] == True)].copy()
        des = des.rename(columns={"BP": "BP_GRCh38"})
        des["source"] = "Desikan2019_PGS"
        des["match_method"] = "desikan_pgs_resolved"
        des = des[["rsID", "CHR", "BP_GRCh38", "A1", "A2", "beta_A1",
                    "source", "match_method"]]
        print(f"[desikan] {len(des)} SNPs (from {des_path.name})")
    else:
        print(f"[desikan] {des_path} MISSING")
        des = pd.DataFrame(columns=bell.columns)

    # ── Consolidated per-SNP source table ────────────────────────────────
    all_snps = pd.concat([bell, wight, kun_min, des], ignore_index=True)
    # overlap flag: a SNP can appear in multiple sources with different β
    counts = all_snps["rsID"].value_counts()
    overlap = set(counts[counts > 1].index)
    all_snps["overlap_multi_source"] = all_snps["rsID"].isin(overlap)
    all_snps.to_csv(args.out / "source_snps.tsv", sep="\t", index=False)
    print(f"\n[out] {args.out/'source_snps.tsv'}  ({len(all_snps)} rows; "
          f"{len(overlap)} rsIDs in ≥2 sources)")

    # ── Per-patient PRS for each source ──────────────────────────────────
    # Use the FULL ADNI dosage (need the 1.49M-SNP plink-extracted dosage).
    # For now load the existing dosage_TSV used by script 20 (430-SNP set)
    # — that covers Bellenguez 27, Wightman 9, Desikan ≤10; for Kunkle's 3
    # genotyped SNPs we also need them. Check which are in DOSAGE_TSV.
    dosage_tsv = _s20.DOSAGE_TSV
    dos = pd.read_csv(dosage_tsv, sep="\t")
    dos["PTID"] = dos["PTID"].astype(str)
    dos = dos.set_index("PTID")
    dos = dos.apply(lambda c: c.fillna(c.mean()))    # cohort-mean impute
    # wfilt dosage covers Wightman re-extracted SNPs
    wfilt_dos_path = GBF / "wfilt_dosage" / "patient_dosage.tsv"
    if wfilt_dos_path.exists():
        wdos = pd.read_csv(wfilt_dos_path, sep="\t")
        wdos["PTID"] = wdos["PTID"].astype(str)
        wdos = wdos.set_index("PTID")
        wdos = wdos.apply(lambda c: c.fillna(c.mean()))
        new_cols = [c for c in wdos.columns if c not in dos.columns]
        if new_cols:
            dos = dos.join(wdos[new_cols], how="left")
            print(f"[dosage] +{len(new_cols)} wfilt SNPs joined "
                  f"(re-extracted)")
    print(f"[dosage] {dos.shape[0]} patients × {dos.shape[1]} SNP cols")

    def _prs(src_df: pd.DataFrame, label: str) -> pd.Series:
        cc = [r for r in src_df["rsID"] if r in dos.columns]
        miss = sorted(set(src_df["rsID"]) - set(cc))
        if miss:
            print(f"  [{label}] {len(miss)} of {len(src_df)} SNPs missing "
                  f"in dosage matrix: {miss[:8]}…")
        bw = src_df.set_index("rsID").loc[cc, "beta_A1"].astype(float)
        sub = dos[cc].astype(float)
        # impute any residual NaN (rare — wfilt cols are mean-imputed above)
        sub = sub.fillna(sub.mean())
        return sub.mul(bw.values, axis=1).sum(axis=1), len(cc)

    prs_df = pd.DataFrame(index=dos.index)
    counts_used = {}
    for label, sub in [("Bellenguez", bell), ("Wightman", wight),
                       ("Kunkle", kun_min), ("Desikan", des)]:
        col, n = _prs(sub, label)
        prs_df[f"PRS_{label}"] = col
        counts_used[label] = n
    prs_df.index.name = "PTID"
    prs_df.reset_index().to_csv(args.out / "patient_source_prs.tsv",
                                  sep="\t", index=False)
    print(f"\n[out] {args.out/'patient_source_prs.tsv'}  "
          f"({len(prs_df)} patients × 4 PRS)")

    # ── Summary ──────────────────────────────────────────────────────────
    n_pub_kunkle = len(parse_kunkle_catalog(KUNKLE_CAT))
    lines = [
        "Source-specific PRS — per-SNP source attribution (overlap allowed)",
        "",
        f"Bellenguez 2022     : 27 published →  {counts_used['Bellenguez']} "
        "SNPs used (all in 430-SNP dosage matrix)",
        f"Wightman   2021     : 38 published →  {counts_used['Wightman']} "
        f"SNPs used  ({len(wfilt)} resolved in-430 + {len(re_ex)} re-extracted "
        "via plink; remaining 29 not genotyped on Omni2.5M)",
        f"Kunkle     2019     : {n_pub_kunkle} catalog leads → "
        f"{counts_used['Kunkle']} SNPs used "
        f"(genotyped+harmonised: {len(kun)}; the rest are below MAF 0.01 "
        "or absent from the Omni2.5M array — rare suggestive hits in the "
        "Kunkle paper at P=1e-6 to 1e-7)",
        f"Desikan PGS         : (PGS000026) →  {counts_used['Desikan']} "
        "SNPs used (from desikan_pgs_resolved.tsv, in_430 subset)",
        "",
        f"Overlap (same rsID in ≥2 sources): {len(overlap)} rsIDs",
    ]
    if len(overlap):
        sample = sorted(overlap)[:10]
        lines.append(f"  example overlapping rsIDs (first 10): {sample}")
    lines += ["", "Kunkle drop log:"]
    if len(kun_drops):
        for _, row in kun_drops.iterrows():
            lines.append(f"  - {row['rsID_pub']:<14} chr{row['chrom_pub']}:"
                         f"{row['bp_pub']}   {row['reason']}")
    (args.out / "source_summary.txt").write_text("\n".join(lines), "utf-8")
    print(f"[out] {args.out/'source_summary.txt'}")
    print("\n=== Counts ===")
    for s, n in counts_used.items():
        print(f"  PRS_{s:<10} = {n} SNPs")


if __name__ == "__main__":
    main()
