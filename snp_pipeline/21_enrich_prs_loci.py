r"""
21_enrich_prs_loci.py
======================
Check whether the Kunkle-2019 and Desikan-2019 (PGS000026) AD GWAS contribute
genome-wide-significant **lead** loci that are genotyped in ADNI but NOT in
the current 430-SNP PRS, and build the artefacts to enrich the β-PRS +
run the published Desikan PHS as an external comparator (script 20).

Sources
-------
- Kunkle 2019: the 24 GWAS-Catalog **lead** loci (GCST007511.tsv); β/SE/p
  fetched from the raw Stage-1 sumstats (keep p<5e-8). Independent by
  construction (lead SNPs) → no LD-clumping.
- Desikan 2019: PGS000026.txt — the published polygenic-hazard score
  (31 SNPs + APOE_e2/APOE_e4 haplotype terms; effect_weight = log-OR-scale).

Matching is by **rsID** (assembly-independent; the GWAS are GRCh37, ADNI is
GRCh38) against the FULL QC'd ADNI BIM; RsMergeArch is a graceful fallback
for deprecated rsIDs. Liftover is NOT needed for SNP identity.

Outputs (under D:\...\genetic_baselines\)
  enriched_snp_weights.tsv  — every candidate lead/PHS SNP: source, β_A1,
                              A1/A2, p, genotyped, in_430, novel
  enriched_rsids.txt        — union(430, novel-genotyped) for plink --extract
  desikan_pgs_resolved.tsv  — Desikan SNP weights (A1-oriented, genotyped)
                              + APOE_e2/APOE_e4 weights
Plus, with --extract-dosages, a wide dosage matrix
  D:\...\bmfm_inputs\patient_genotypes_enriched\patient_dosage.{tsv,bim}
(via plink --recode A --keep-allele-order; no imputation).

Usage
-----
  conda run -n snp python snp_pipeline/21_enrich_prs_loci.py
  conda run -n snp python snp_pipeline/21_enrich_prs_loci.py --extract-dosages
"""
from __future__ import annotations

import argparse
import gzip
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("D:/ADNI_SNP_Omni2.5M_20140220")
KUNKLE_CAT = DATA / "GWAS/Kunkle_2019/Kunkle_2019_gwas-association-downloaded_2026-05-18-accessionId_GCST007511.tsv"
KUNKLE_S1  = DATA / "GWAS/Kunkle_2019/NG00075_Kunkle_etal_Stage1_P-val_only_results.txt"
DESIKAN    = DATA / "GWAS/Desikan_2019/PGS000026.txt"
FULL_BIM   = DATA / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr.bim"
PLINK_STEM = DATA / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr"
LABELS     = DATA / "bmfm_inputs/external_gwas_labels.tsv"
DOSAGE430  = DATA / "bmfm_inputs/patient_genotypes/patient_dosage.tsv"
RSMERGE    = DATA / "liftover/NCBI/RsMergeArch.bcp.gz"
OUT        = DATA / "genetic_baselines"
ENR_DIR    = DATA / "bmfm_inputs/patient_genotypes_enriched"
GW_SIG     = 5e-8


def parse_kunkle_catalog() -> pd.DataFrame:
    """24 lead loci → rsID, risk_allele, p, or  (GWAS-Catalog schema)."""
    b = pd.read_csv(KUNKLE_CAT, sep="\t", low_memory=False)
    risk = (b["STRONGEST SNP-RISK ALLELE"].astype(str)
            .str.rsplit("-", n=1).str[-1].str.upper().str.strip())
    p = pd.to_numeric(b["P-VALUE"], errors="coerce")
    orb = pd.to_numeric(b["OR or BETA"], errors="coerce")
    df = pd.DataFrame({"rsID": b["SNPS"].astype(str).str.strip(),
                       "risk_allele": risk, "p": p, "cat_or": orb})
    df = df[df["rsID"].str.startswith("rs") & (df["p"] < GW_SIG)]
    return df.drop_duplicates("rsID").reset_index(drop=True)


def fetch_kunkle_stage1(rsids: set[str]) -> dict:
    """Stream the 544 MB Stage-1 file; for the target leads return
    rsID → (effect_allele, beta_logOR, se, p).  Keep p<5e-8."""
    cols = ["MarkerName", "Effect_allele", "Non_Effect_allele",
            "Beta", "SE", "Pvalue"]
    out: dict[str, tuple] = {}
    for ch in pd.read_csv(KUNKLE_S1, sep=r"\s+", usecols=cols,
                          chunksize=1_000_000, low_memory=False):
        hit = ch[ch["MarkerName"].astype(str).isin(rsids)]
        for r in hit.itertuples(index=False):
            try:
                p = float(r.Pvalue)
                if p >= GW_SIG:
                    continue
                out.setdefault(str(r.MarkerName), (
                    str(r.Effect_allele).upper().strip(),
                    float(r.Beta), float(r.SE), p))
            except (ValueError, TypeError):
                continue
        if len(out) >= len(rsids):
            break
    return out


def parse_desikan() -> tuple[pd.DataFrame, dict]:
    """SNP weights (is_haplotype FALSE) + APOE term weights (TRUE)."""
    d = pd.read_csv(DESIKAN, sep="\t", comment="#")
    d.columns = [c.strip() for c in d.columns]
    is_hap = d["is_haplotype"].astype(str).str.upper() == "TRUE"
    apoe = {}
    for r in d[is_hap].itertuples(index=False):
        apoe[str(r.effect_allele).strip()] = float(r.effect_weight)  # APOE_e2/e4
    snp = d[~is_hap & d["rsID"].astype(str).str.startswith("rs")].copy()
    snp = pd.DataFrame({
        "rsID": snp["rsID"].astype(str).str.strip(),
        "effect_allele": snp["effect_allele"].astype(str).str.upper().str.strip(),
        "weight": pd.to_numeric(snp["effect_weight"], errors="coerce"),
        "gene": snp["locus_name"].astype(str)})
    return snp.dropna(subset=["weight"]).reset_index(drop=True), apoe


def load_full_bim() -> pd.DataFrame:
    b = pd.read_csv(FULL_BIM, sep=r"\s+", header=None,
                    names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
                    usecols=[0, 1, 3, 4, 5], dtype={"SNP": str})
    b["SNP"] = b["SNP"].astype(str)
    return b.drop_duplicates("SNP").set_index("SNP")


def rsmerge_map(query: set[str]) -> dict:
    """Best-effort deprecated→current rsID map for the few unmatched leads.
    RsMergeArch.bcp: col0=rsHigh(merged-away) col1=rsLow(kept)."""
    want = {q[2:] for q in query if q.startswith("rs") and q[2:].isdigit()}
    if not want:
        return {}
    m: dict[str, str] = {}
    try:
        with gzip.open(RSMERGE, "rt") as fh:
            for line in fh:
                p = line.split("\t")
                if len(p) >= 2 and p[0] in want:
                    m["rs" + p[0]] = "rs" + p[1]
    except Exception as e:
        print(f"[rsmerge] skipped ({type(e).__name__}: {e})")
    return m


def harmonise(rsid, eff_allele, beta, bim):
    """β oriented to BIM A1 (dosage-counted). None if SNP absent / ambiguous."""
    if rsid not in bim.index:
        return None
    a1 = str(bim.loc[rsid, "A1"]).upper()
    a2 = str(bim.loc[rsid, "A2"]).upper()
    ea = str(eff_allele).upper()
    if ea == a1:
        flip = False
    elif ea == a2:
        flip = True
    else:
        return None                       # strand/allele mismatch → drop
    return dict(CHR=str(bim.loc[rsid, "CHR"]), BP=int(bim.loc[rsid, "BP"]),
                A1=a1, A2=a2, genotyped=(a1 not in (".", "0")),
                beta_A1=(-beta if flip else beta), allele_flip=flip)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--extract-dosages", action="store_true",
                    help="Run plink --recode A for union(430, novel) → "
                         "patient_genotypes_enriched/patient_dosage.{tsv,bim}")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    # current 430 (genotyped GW-sig already in the PRS)
    lab = pd.read_csv(LABELS, sep="\t", low_memory=False)
    cur430 = set(pd.read_csv(DOSAGE430, sep="\t", nrows=0).columns) - {"PTID"}
    print(f"[base] current PRS SNPs (patient_dosage): {len(cur430)}")

    bim = load_full_bim()
    print(f"[bim] full ADNI panel rsIDs: {len(bim):,}")

    # ── Kunkle: 24 catalog leads, β/SE/p from Stage-1 ───────────────────────
    kcat = parse_kunkle_catalog()
    print(f"[kunkle] catalog GW-sig lead loci: {len(kcat)}")
    s1 = fetch_kunkle_stage1(set(kcat["rsID"]))
    print(f"[kunkle] of those, Stage-1 p<5e-8 with β: {len(s1)}")

    # ── Desikan PGS000026 ───────────────────────────────────────────────────
    dsnp, apoe_w = parse_desikan()
    print(f"[desikan] PGS000026 SNPs: {len(dsnp)}  APOE terms: {apoe_w}")

    # candidate (source, rsID, effect_allele, beta) before harmonisation
    cand = []
    for r in kcat.itertuples(index=False):
        if r.rsID in s1:
            ea, be = s1[r.rsID][0], s1[r.rsID][1]; pv = s1[r.rsID][3]
        else:                                  # fall back to catalog OR
            ea, be, pv = r.risk_allele, math.log(float(r.cat_or)) if r.cat_or > 0 else np.nan, r.p
        if np.isfinite(be):
            cand.append(("Kunkle2019", r.rsID, ea, be, pv))
    for r in dsnp.itertuples(index=False):
        cand.append(("Desikan2019", r.rsID, r.effect_allele, float(r.weight), np.nan))

    # rsID canonicalisation fallback for any not in the BIM
    unmatched = {c[1] for c in cand if c[1] not in bim.index}
    rmap = rsmerge_map(unmatched) if unmatched else {}
    if rmap:
        print(f"[rsmerge] resolved {sum(v in bim.index for v in rmap.values())}"
              f"/{len(unmatched)} deprecated rsIDs")

    rows = []
    for src, rsid, ea, be, pv in cand:
        rid = rsid if rsid in bim.index else rmap.get(rsid, rsid)
        h = harmonise(rid, ea, be, bim)
        rec = {"rsID": rid, "query_rsID": rsid, "source": src,
               "effect_allele_raw": ea, "beta_raw": be, "p": pv,
               "in_bim": rid in bim.index}
        if h:
            rec.update(h)
        rec["in_430"] = rid in cur430
        rec["novel"] = bool(h) and h["genotyped"] and rid not in cur430
        rows.append(rec)
    enr = pd.DataFrame(rows).drop_duplicates("rsID")
    for _c in ("CHR", "BP", "A1", "A2", "beta_A1", "allele_flip", "genotyped"):
        if _c not in enr.columns:
            enr[_c] = np.nan
    enr["genotyped"] = enr["genotyped"].fillna(False).astype(bool)
    enr.to_csv(args.out / "enriched_snp_weights.tsv", sep="\t", index=False)

    # Desikan resolved (genotyped, A1-oriented) + APOE terms for the PHS
    dres = enr[(enr["source"] == "Desikan2019") & enr["genotyped"]]
    dres[["rsID", "CHR", "BP", "A1", "A2", "beta_A1", "in_430"]].to_csv(
        args.out / "desikan_pgs_resolved.tsv", sep="\t", index=False)
    pd.Series(apoe_w).to_csv(args.out / "desikan_apoe_weights.tsv",
                             sep="\t", header=["weight"])

    # ── Report (this is the answer to "did we miss loci?") ──────────────────
    def _rep(src):
        s = enr[enr["source"] == src]
        g = s[s["genotyped"]]
        nov = s[s["novel"]]
        print(f"\n[{src}] candidates={len(s)}  genotyped_in_ADNI={len(g)}  "
              f"already_in_430={int(s['in_430'].sum())}  NOVEL={len(nov)}")
        if len(nov):
            cols = [c for c in ["rsID", "CHR", "BP", "beta_A1", "p"] if c in nov]
            print(nov[cols].to_string(index=False))
    for src in ("Kunkle2019", "Desikan2019"):
        _rep(src)
    novel = enr[enr["novel"]]["rsID"].tolist()
    print(f"\n[ENRICHMENT] {len(novel)} novel genotyped GW-sig lead SNP(s) "
          f"beyond the current {len(cur430)}: {novel}")

    union = sorted(cur430 | set(novel))
    (args.out / "enriched_rsids.txt").write_text("\n".join(union) + "\n")
    print(f"[out] enriched_rsids.txt: {len(union)} rsIDs "
          f"({len(cur430)} + {len(novel)} novel)")
    print(f"[out] {args.out/'enriched_snp_weights.tsv'}")
    print(f"[out] {args.out/'desikan_pgs_resolved.tsv'}  "
          f"({len(dres)} genotyped PGS000026 SNPs)")

    # ── Optional: re-extract the wider dosage matrix via plink ──────────────
    if args.extract_dosages:
        ENR_DIR.mkdir(parents=True, exist_ok=True)
        raw_stem = ENR_DIR / "enriched"
        cmd = ["plink", "--bfile", str(PLINK_STEM),
               "--extract", str(args.out / "enriched_rsids.txt"),
               "--recode", "A", "--keep-allele-order",
               "--out", str(raw_stem)]
        print(f"\n[plink] {' '.join(cmd)}", flush=True)
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout[-1500:]); print(r.stderr[-1500:])
            sys.exit("plink --recode A failed")
        raw = pd.read_csv(str(raw_stem) + ".raw", sep=r"\s+")
        keep = [c for c in raw.columns if c not in
                ("FID", "PAT", "MAT", "SEX", "PHENOTYPE")]
        raw = raw[keep].rename(columns={"IID": "PTID"})
        raw.columns = ["PTID"] + [c.rsplit("_", 1)[0] for c in raw.columns[1:]]
        raw.to_csv(ENR_DIR / "patient_dosage.tsv", sep="\t", index=False)
        snpcols = [c for c in raw.columns if c != "PTID"]
        bim.loc[[s for s in snpcols if s in bim.index],
                ["CHR", "BP", "A1", "A2"]].reset_index().rename(
            columns={"index": "rsID", "SNP": "rsID"}).to_csv(
            ENR_DIR / "patient_dosage.bim", sep="\t", index=False)
        print(f"[plink] wrote {ENR_DIR/'patient_dosage.tsv'} "
              f"({len(raw)} patients x {len(snpcols)} SNPs)")


if __name__ == "__main__":
    main()
