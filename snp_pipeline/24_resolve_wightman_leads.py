r"""
24_resolve_wightman_leads.py
============================
Robustly resolve the 38 PUBLISHED Wightman-2021 genome-wide-significant lead
variants against our GRCh38 ADNI pipeline, so the PRS can be restricted to
true published leads instead of the ~407 unfiltered sub-threshold Wightman
variants currently in `external_gwas_labels.tsv` (label==1).

Naïve rsID matching loses leads to rsID-merge / build drift, so for each
lead we:
  1. RsMergeArch-canonicalise the lead rsID (deprecated → current),
  2. pyliftover GRCh37 → GRCh38 the lead position (cached UCSC chain),
  3. find β / effect-allele in the TRUSTED ADNI-hg19 Wightman file first
     (`…_trusted_adni_hg19.tsv`, has BP_37+BP_38+beta), else the full
     GRCh38 Wightman file,
  4. determine whether the lead is in `external_gwas_labels[label==1]`
     (the 430), the dosage matrix (genotyped now), and the full GRCh38
     array — matching on {rsID, rsID_canonical} OR GRCh38 CHR:BP.

Reuses the liftover + RsMergeArch conventions from
`08a_liftover_hg19_to_hg38_wightman_2021.py`.

Outputs (under --out, default genetic_baselines_filtered/):
  * wightman_lead_resolution_report.tsv  (38 rows, full diagnostics)
  * wightman_resolved_allowlist.tsv      (resolved_genotyped subset →
                                          the rsIDs to keep as the filtered
                                          Wightman set for scripts 20/23)
  + a console summary (counts by status; explicit no-β / not-on-array list).

This is a CHECKPOINT: review the report before the filtered re-run.

Usage:  conda run -n snp python snp_pipeline/24_resolve_wightman_leads.py
"""
from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
GWAS = BASE / "GWAS" / "Wightman_2021"
FILTERED = GWAS / "wightman_2021_filtered.tsv"
TRUSTED = GWAS / "wightman_without_ukb_GRCh38_trusted_adni_hg19.tsv"
FULL = GWAS / "wightman_without_ukb_GRCh38.tsv"
# with-UKB = full PGC-ALZ2 meta (incl. UKB/23andMe). ADNI is NOT in
# UKB/23andMe ⇒ no sample overlap with our target; used to RECOVER leads
# absent from the (more independent) without-UKB source.
TRUSTED_W = GWAS / "wightman_with_ukb_GRCh38_trusted_adni_hg19.tsv"
FULL_W = GWAS / "wightman_with_ukb_GRCh38.tsv"
# RAW genome-wide PGC-ALZ2 (without UKB&23andMe) — 12.2M variants, has
# β+EAF, GRCh37, chr:pos-keyed (no rsID). The lifted *_GRCh38.tsv above are
# only ~3k-row derived subsets; the published leads live here. Recover by
# GRCh37 chr:pos (same provenance as the 5 already resolved → consistent β).
RAW_NOUKB = (GWAS /
             "PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis (1).txt" /
             "PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis (1).txt")
CHAIN = BASE / "liftover" / "hg19ToHg38.over.chain.gz"
MERGE_ARCH = BASE / "liftover" / "NCBI" / "RsMergeArch.bcp.gz"
EXT_LABELS = BASE / "bmfm_inputs" / "external_gwas_labels.tsv"
DOSAGE_BIM = BASE / "bmfm_inputs" / "patient_genotypes" / "patient_dosage.bim"
FULL_BIM = BASE / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr.bim"
PLINK_STEM = BASE / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr"
DOSAGE_TSV = BASE / "bmfm_inputs" / "patient_genotypes" / "patient_dosage.tsv"
OUT_DEFAULT = BASE / "genetic_baselines_filtered"
CANON = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}


def _rsnum(x) -> int | None:
    s = str(x).strip().lower().replace("rs", "")
    return int(s) if s.isdigit() else None


def merge_canonical(lead_rsids: list[str]) -> dict:
    """{lead_rsid → canonical rsID} via RsMergeArch (low→high)."""
    ints = {n for n in (_rsnum(r) for r in lead_rsids) if n is not None}
    mm: dict[int, int] = {}
    if MERGE_ARCH.exists() and ints:
        with gzip.open(MERGE_ARCH, "rt", encoding="utf-8",
                       errors="replace") as fh:
            for line in fh:
                p = line.split("\t")
                if len(p) >= 2:
                    try:
                        high, low = int(p[0]), int(p[1])
                    except ValueError:
                        continue
                    if low in ints:
                        mm[low] = high
    out = {}
    for r in lead_rsids:
        n = _rsnum(r)
        out[r] = f"rs{mm[n]}" if (n is not None and n in mm) else str(r)
    return out


def _read_wightman(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    for c in ("rsID", "rsID_canonical"):
        if c in df.columns:
            df[c + "_l"] = df[c].astype(str).str.strip().str.lower()
    return df


def _lookup_beta(df: pd.DataFrame, ids: set, chrom: str,
                 pos37: int, pos38, require_lifted: bool) -> tuple | None:
    """Return (effect_allele, beta, eaf, by) or None."""
    m = pd.Series(False, index=df.index)
    if "rsID_l" in df:
        m |= df["rsID_l"].isin(ids)
    if "rsID_canonical_l" in df:
        m |= df["rsID_canonical_l"].isin(ids)
    by = "rsID/canonical"
    if not m.any() and {"CHR_37", "BP_37"} <= set(df.columns):
        pm = ((df["CHR_37"].astype(str) == str(chrom))
              & (pd.to_numeric(df["BP_37"], errors="coerce") == pos37))
        if pm.any():
            m, by = pm, "pos_GRCh37"
    if not m.any() and pos38 is not None and {"CHR_38", "BP_38"} <= set(df.columns):
        pm = ((df["CHR_38"].astype(str) == str(chrom))
              & (pd.to_numeric(df["BP_38"], errors="coerce") == pos38))
        if require_lifted and "lifted" in df.columns:
            pm &= df["lifted"].astype(str).str.lower().isin(["true", "1"])
        if pm.any():
            m, by = pm, "pos_GRCh38(lifted)"
    if not m.any():
        return None
    r = df[m].iloc[0]
    try:
        oa = str(r.get("other_allele", "")).upper().strip()
        return (str(r["effect_allele"]).upper().strip(), float(r["beta"]),
                float(r["effect_allele_frequency"]), oa, by)
    except (ValueError, TypeError, KeyError):
        return None


def scan_raw_noukb(targets: set) -> dict:
    """Stream the 12.2M-row RAW without-UKB PGC-ALZ2 sumstats; return
    {(chr,pos37_str) → (effect_allele, other_allele, beta, eaf)} for the
    lead targets. Cols: chromosome base_pair_location effect_allele
    other_allele beta standard_error effect_allele_frequency p_value SNPID
    N Neffective Build. GRCh37, chr:pos-keyed (no rsID)."""
    out: dict = {}
    if not RAW_NOUKB.exists() or not targets:
        return out
    with open(RAW_NOUKB, "r", encoding="utf-8", errors="replace") as fh:
        next(fh, None)                                   # header
        for line in fh:
            p = line.split()
            if len(p) < 7:
                continue
            k = (p[0], p[1])
            if k in targets:
                try:
                    out[k] = (p[2].upper().strip(), p[3].upper().strip(),
                              float(p[4]), float(p[6]))
                except (ValueError, IndexError):
                    pass
                if len(out) == len(targets):
                    break
    return out


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--extract-dosages", action="store_true",
                    help="plink-extract the filtered SNP set (Bellenguez "
                         "430∩dosage + 6 in-430 Wightman leads + the 3 "
                         "on_array_not_in_430 leads) from the full GRCh38 "
                         "set; write the filtered dosage matrix + the "
                         "3-lead enriched-weights for script-20 injection.")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    leads = pd.read_csv(FILTERED, sep="\t")
    leads["Lead_variant"] = leads["Lead_variant"].astype(str).str.strip()
    print(f"[in] {len(leads)} published Wightman-2021 leads", flush=True)

    canon = merge_canonical(leads["Lead_variant"].tolist())
    n_merged = sum(1 for k, v in canon.items() if v.lower() != k.lower())
    print(f"[RsMergeArch] {n_merged} lead rsIDs canonicalised (deprecated→current)",
          flush=True)

    from pyliftover import LiftOver
    print("[liftover] loading chain (first load builds index ~30s)…",
          flush=True)
    lo = LiftOver(str(CHAIN))

    trusted = _read_wightman(TRUSTED)
    full = _read_wightman(FULL)
    trustedW = _read_wightman(TRUSTED_W)
    fullW = _read_wightman(FULL_W)
    print(f"[β src] without-UKB trusted={len(trusted)}/full={len(full)}; "
          f"with-UKB trusted={len(trustedW)}/full={len(fullW)}", flush=True)

    _extcols = pd.read_csv(EXT_LABELS, sep="\t", nrows=0).columns.tolist()
    _src = [c for c in ("z_source", "source") if c in _extcols]
    ext = pd.read_csv(EXT_LABELS, sep="\t",
                      usecols=["SNP", "CHR", "BP", "label"] + _src,
                      low_memory=False)
    ext = ext[ext["label"] == 1].copy()
    ext["SNP_l"] = ext["SNP"].astype(str).str.strip().str.lower()
    ext["_srcjoin"] = ext[_src].astype(str).agg(" ".join, axis=1) if _src \
        else ""
    ext_pos = {(str(c), int(b)): str(s) for s, c, b in
               zip(ext["SNP"], ext["CHR"], pd.to_numeric(ext["BP"],
                   errors="coerce").fillna(-1).astype(int))}
    ext_ids = set(ext["SNP_l"])
    print(f"[ext] external_gwas_labels[label==1] = {len(ext)} SNPs", flush=True)

    dbim = pd.read_csv(DOSAGE_BIM, sep="\t")
    dbim["rsID_l"] = dbim["rsID"].astype(str).str.strip().str.lower()
    dos_ids = set(dbim["rsID_l"])
    dos_pos = {(str(c), int(b)) for c, b in zip(dbim["CHR"],
               pd.to_numeric(dbim["BP"], errors="coerce").fillna(-1).astype(int))}

    fbim = pd.read_csv(FULL_BIM, sep=r"\s+", header=None,
                       usecols=[0, 1, 3], names=["CHR", "rsID", "BP"],
                       low_memory=False)
    fbim["rsID_l"] = fbim["rsID"].astype(str).str.strip().str.lower()
    farr_ids = set(fbim["rsID_l"])
    farr_idmap = dict(zip(fbim["rsID_l"], fbim["rsID"].astype(str)))
    farr_pos = {(str(c), int(b)): str(r) for c, b, r in zip(fbim["CHR"],
                pd.to_numeric(fbim["BP"], errors="coerce").fillna(-1)
                .astype(int), fbim["rsID"])}
    print(f"[array] full GRCh38 array BIM = {len(fbim)} SNPs", flush=True)

    raw_targets = {(str(L["Chromosome"]).strip(),
                    str(int(L["Position_GRCh37"])))
                   for _, L in leads.iterrows()}
    print(f"[raw] scanning RAW without-UKB sumstats (12.2M rows) for "
          f"{len(raw_targets)} lead positions…", flush=True)
    raw_map = scan_raw_noukb(raw_targets)
    print(f"[raw] matched {len(raw_map)}/{len(raw_targets)} leads by "
          f"GRCh37 chr:pos in the raw genome-wide without-UKB meta",
          flush=True)

    rows = []
    for _, L in leads.iterrows():
        rsid = str(L["Lead_variant"]).strip()
        cano = canon.get(rsid, rsid)
        chrom = str(L["Chromosome"]).strip()
        pos37 = int(L["Position_GRCh37"])
        ids = {rsid.lower(), cano.lower()}

        res = lo.convert_coordinate(f"chr{chrom}", pos37 - 1)
        pos38 = None
        if res:
            nc, npn, *_ = res[0]
            if nc.replace("chr", "") == chrom:
                pos38 = int(npn) + 1

        # β precedence: independent without-UKB lifted subsets → with-UKB
        # lifted subsets → RAW genome-wide without-UKB sumstats (same
        # provenance as the resolved 5, consistent β, no ADNI/UKB overlap).
        b, bsrc = None, ""
        for name, df_, reqL in (("trusted_noUKB", trusted, False),
                                 ("full_noUKB", full, True),
                                 ("trusted_withUKB", trustedW, False),
                                 ("full_withUKB", fullW, True)):
            hit = _lookup_beta(df_, ids, chrom, pos37, pos38, reqL)
            if hit:
                b, bsrc = hit, name
                break
        if b is None and (chrom, str(pos37)) in raw_map:
            ea_, oa_, be_, ef_ = raw_map[(chrom, str(pos37))]
            b = (ea_, be_, ef_, oa_, "raw_GRCh37_posmatch")
            bsrc = "raw_noUKB_fullmeta"
        beta = b[1] if b else np.nan
        ea = b[0] if b else ""
        oa = b[3] if b else ""

        pos38_key = (chrom, pos38) if pos38 is not None else None
        in_ext = bool(ids & ext_ids) or (pos38_key in ext_pos)
        pipe_rsid = ""
        if ids & ext_ids:
            pipe_rsid = next(s for s in (rsid, cano)
                             if s.lower() in ext_ids)
        elif pos38_key in ext_pos:
            pipe_rsid = ext_pos[pos38_key]
        genotyped = bool(ids & dos_ids) or (pos38_key in dos_pos)
        on_array = bool(ids & farr_ids) or (pos38_key in farr_pos)
        array_rsID = ""
        _hit_id = ids & farr_ids
        if _hit_id:
            array_rsID = farr_idmap[next(iter(_hit_id))]
        elif pos38_key in farr_pos:
            array_rsID = farr_pos[pos38_key]

        if not b:
            status = "no_beta"
        elif in_ext and genotyped:
            status = "resolved_genotyped"
        elif on_array:
            status = "on_array_not_in_430"
        else:
            status = "beta_but_not_on_array"

        rows.append({
            "Genomic_locus": L.get("Genomic_locus"), "Gene": L.get("Gene"),
            "lead_rsID": rsid, "lead_rsID_canonical": cano,
            "chr": chrom, "pos_GRCh37": pos37,
            "pos_GRCh38": pos38 if pos38 is not None else "",
            "lifted": pos38 is not None, "P_value": L.get("P_value"),
            "beta": beta, "effect_allele": ea, "other_allele": oa,
            "beta_source": bsrc,
            "in_ext_labels_430": in_ext, "pipeline_rsID": pipe_rsid,
            "genotyped_dosage": genotyped, "on_full_array": on_array,
            "array_rsID": array_rsID, "status": status,
        })

    rep = pd.DataFrame(rows)
    rfile = args.out / "wightman_lead_resolution_report.tsv"
    rep.to_csv(rfile, sep="\t", index=False)
    print(f"\n[out] {rfile}  ({len(rep)} rows)")

    allow = rep[rep["status"] == "resolved_genotyped"][
        ["pipeline_rsID", "lead_rsID", "lead_rsID_canonical",
         "beta", "beta_source"]].copy()
    afile = args.out / "wightman_resolved_allowlist.tsv"
    allow.to_csv(afile, sep="\t", index=False)
    print(f"[out] {afile}  ({len(allow)} resolved-genotyped leads)")

    print("\n[QC] status counts:")
    for s, n in rep["status"].value_counts().items():
        print(f"  {s:24s} {n}")
    rg = rep[rep["status"] == "resolved_genotyped"]
    print(f"\n[QC] resolved & genotyped (→ filtered Wightman set): {len(rg)}")
    print("[QC] resolved-genotyped β provenance:")
    for s, n in rg["beta_source"].value_counts().items():
        print(f"  {s:18s} {n}")
    wuk = rep[rep["beta_source"].astype(str).str.contains("withUKB")]
    print(f"[QC] β RECOVERED from with-UKB meta: {len(wuk)} leads "
          f"({int((wuk['status']=='resolved_genotyped').sum())} of them "
          f"genotyped → enter the PRS)")
    nb = rep[rep["status"] == "no_beta"]["lead_rsID"].tolist()
    print(f"\n[QC] leads with NO β found ({len(nb)}): {nb}")
    na = rep[~rep["on_full_array"]]["lead_rsID"].tolist()
    print(f"\n[QC] leads NOT on the GRCh38 array ({len(na)}): {na}")
    oa = rep[rep["status"] == "on_array_not_in_430"][
        ["lead_rsID", "chr", "pos_GRCh38"]]
    print(f"\n[QC] on the array but DROPPED from the 430 "
          f"({len(oa)}) — relevant to scope decision:")
    if len(oa):
        print(oa.to_string(index=False))

    if not args.extract_dosages:
        return

    # ── gene annotation map (Bellenguez GWAS-Catalog: covers Bellenguez +
    #    the AD Desikan loci; one reproducible source) ───────────────────────
    BELL_CAT = (BASE / "GWAS" / "Bellenguez_2022" /
                "gwas-association-downloaded_2026-04-28-pubmedId_35379992.tsv")
    gene_map: dict = {}
    if BELL_CAT.exists():
        gc = pd.read_csv(BELL_CAT, sep="\t", low_memory=False)
        for _, x in gc.iterrows():
            g = str(x.get("MAPPED_GENE", "")).strip()
            if not g or g.lower() == "nan":
                g = str(x.get("REPORTED GENE(S)", "")).strip()
            for tok in str(x.get("SNPS", "")).replace(";", " ").replace(
                    ",", " ").split():
                t = tok.strip().lower()
                if t.startswith("rs") and t not in gene_map and g \
                        and g.lower() != "nan":
                    gene_map[t] = g

    # Desikan locus names from the PGS Catalog scoring file (authoritative
    # for the Desikan loci; Bellenguez catalog is the fallback).
    DESIKAN_PGS = BASE / "GWAS" / "Desikan_2019" / "PGS000026.txt"
    dlocus: dict = {}
    if DESIKAN_PGS.exists():
        dp = pd.read_csv(DESIKAN_PGS, sep="\t", comment="#", low_memory=False)
        for _, x in dp.iterrows():
            rid = str(x.get("rsID", "")).strip().lower()
            ln = str(x.get("locus_name", "")).strip()
            if rid.startswith("rs") and ln and ln.lower() != "nan":
                dlocus.setdefault(rid, ln)

    def _gene(rsid, fallback=""):
        k = str(rsid).strip().lower()
        return dlocus.get(k) or gene_map.get(k, fallback)

    # ── filtered-set members (positions only; pre-plink) ───────────────────
    bell_rows = ext[ext["_srcjoin"].str.contains("Bellenguez", case=False,
                                                 na=False)
                    & ext["SNP_l"].isin(dos_ids)]
    bell_geno = sorted({s for s in bell_rows["SNP"].astype(str)})
    w6 = sorted({s for s in rg["pipeline_rsID"].astype(str) if s})
    reext = rep[rep["status"] == "on_array_not_in_430"].copy()
    w3 = sorted({s for s in reext["array_rsID"].astype(str) if s})
    fset = ([(str(w["chr"]), int(float(w["pos_GRCh38"])))
             for _, w in rep[rep["status"].isin(["resolved_genotyped",
                 "on_array_not_in_430"])].iterrows()
             if str(w["pos_GRCh38"]).strip() not in ("", "nan")]
            + [(str(c), int(b)) for c, b in zip(bell_rows["CHR"],
               pd.to_numeric(bell_rows["BP"], errors="coerce")
               .fillna(-1).astype(int))])

    # ── Desikan genotyped loci + non-overlap classification ────────────────
    dres_f = BASE / "genetic_baselines" / "desikan_pgs_resolved.tsv"
    dd, desikan5 = pd.DataFrame(), []
    if dres_f.exists():
        dres = pd.read_csv(dres_f, sep="\t")
        WIN = 250_000                                  # ±250 kb LD proxy
        drows = []
        for _, d in dres.iterrows():
            dc, db = str(d["CHR"]), int(float(d["BP"]))
            near = [abs(db - b) for c, b in fset if c == dc]
            mind = min(near) if near else None
            ov = mind is not None and mind <= WIN
            drows.append({"rsID": d["rsID"], "gene": _gene(d["rsID"]),
                          "CHR": dc, "BP_GRCh38": db, "A1": d["A1"],
                          "A2": d["A2"], "beta_A1": d["beta_A1"],
                          "in_orig_430": d["in_430"],
                          "nearest_filtered_bp": (mind if mind is not None
                                                  else ""),
                          "overlaps_filtered_pm250kb": ov,
                          "classification": ("overlapping" if ov
                                             else "additional_nonoverlap")})
        dd = pd.DataFrame(drows)
        desikan5 = dd.loc[~dd["overlaps_filtered_pm250kb"],
                          "rsID"].astype(str).tolist()

    # ── plink-extract: Bellenguez ∪ 6 in-430 W ∪ 3 re-extract ∪ 5 Desikan ──
    keep = sorted(set(bell_geno) | set(w6) | set(w3) | set(desikan5))
    rfile = args.out / "wfilt_rsids.txt"
    rfile.write_text("\n".join(keep) + "\n", encoding="utf-8")
    print(f"\n[wfilt] extract set = {len(bell_geno)} Bellenguez + {len(w6)} "
          f"in-430 W + {len(w3)} re-extract W + {len(desikan5)} novel "
          f"Desikan = {len(keep)} SNPs → {rfile}")

    wdir = args.out / "wfilt_dosage"
    wdir.mkdir(parents=True, exist_ok=True)
    stem = wdir / "wfilt"
    cmd = ["plink", "--bfile", str(PLINK_STEM), "--extract", str(rfile),
           "--recode", "A", "--keep-allele-order", "--out", str(stem)]
    print(f"[plink] {' '.join(cmd)}", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-1500:]); print(r.stderr[-1500:])
        sys.exit("plink --recode A failed")
    raw = pd.read_csv(str(stem) + ".raw", sep=r"\s+")
    keepc = [c for c in raw.columns
             if c not in ("FID", "PAT", "MAT", "SEX", "PHENOTYPE")]
    raw = raw[keepc].rename(columns={"IID": "PTID"})
    raw.columns = ["PTID"] + [c.rsplit("_", 1)[0] for c in raw.columns[1:]]
    raw.to_csv(wdir / "patient_dosage.tsv", sep="\t", index=False)
    snpcols = [c for c in raw.columns if c != "PTID"]
    fb = pd.read_csv(str(PLINK_STEM) + ".bim", sep=r"\s+", header=None,
                     names=["CHR", "rsID", "cM", "BP", "A1", "A2"])
    fb["rsID"] = fb["rsID"].astype(str)
    fb = fb.drop_duplicates("rsID").set_index("rsID")
    fb.loc[[s for s in snpcols if s in fb.index],
           ["CHR", "BP", "A1", "A2"]].reset_index().rename(
        columns={"index": "rsID"}).to_csv(wdir / "patient_dosage.bim",
                                           sep="\t", index=False)
    print(f"[plink] wrote {wdir/'patient_dosage.tsv'} "
          f"({len(raw)} patients x {len(snpcols)} SNPs)")

    # Strand-aware biallelic harmonisation to the EXTRACTED array A1
    # (reverse-complement handled; palindromic A/T & C/G NOT strand-flipped).
    COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

    def _orient(ea, oa, a1, a2):
        S = {a1, a2}
        if a1 in (".", "0") or a2 in (".", "0") or len(S) < 2:
            return None
        pal = S == {"A", "T"} or S == {"C", "G"}
        if oa and {ea, oa} == S:
            return 1.0 if ea == a1 else -1.0
        if (oa and not pal and ea in COMP and oa in COMP
                and {COMP[ea], COMP[oa]} == S):
            return 1.0 if COMP[ea] == a1 else -1.0
        if not oa:
            if ea == a1:
                return 1.0
            if ea == a2:
                return -1.0
            if not pal and ea in COMP:
                if COMP[ea] == a1:
                    return 1.0
                if COMP[ea] == a2:
                    return -1.0
        return None

    def _harm(rsid, ea, oa, beta, src, gene, lead=""):
        if rsid not in fb.index:
            return None, (rsid, "absent_in_plink_bim")
        a1 = str(fb.loc[rsid, "A1"]).upper()
        a2 = str(fb.loc[rsid, "A2"]).upper()
        s = _orient(ea.upper(), str(oa).upper().strip(), a1, a2)
        if s is None:
            return None, (rsid, f"unresolvable ea={ea} oa={oa} "
                                 f"a1={a1} a2={a2}")
        return ({"rsID": rsid, "lead_rsID": lead, "Gene": gene,
                 "source": src, "CHR": str(fb.loc[rsid, "CHR"]),
                 "BP": int(fb.loc[rsid, "BP"]), "A1": a1, "A2": a2,
                 "beta_A1": s * float(beta), "allele_flip": s < 0,
                 "genotyped": a1 not in (".", "0"), "novel": True}), None

    ew_w, ew_d, dropped = [], [], []
    for _, L in reext.iterrows():               # 3 re-extract Wightman leads
        row, err = _harm(str(L["array_rsID"]), str(L["effect_allele"]),
                         L.get("other_allele", ""), L["beta"],
                         "Wightman2021_lead", L["Gene"], L["lead_rsID"])
        (ew_w.append(row) if row else dropped.append(err))
    for _, D in dd[~dd["overlaps_filtered_pm250kb"]].iterrows() \
            if len(dd) else []:                 # 5 novel Desikan loci
        row, err = _harm(str(D["rsID"]), str(D["A1"]), str(D["A2"]),
                         D["beta_A1"], "Desikan2019_PGS",
                         _gene(D["rsID"]), str(D["rsID"]))
        (ew_d.append(row) if row else dropped.append(err))

    ew = pd.DataFrame(ew_w)
    ewf = args.out / "wfilt_enriched_weights.tsv"
    ew.to_csv(ewf, sep="\t", index=False)
    ewd = pd.DataFrame(ew_w + ew_d)
    ewdf = args.out / "wfilt_desikan_enriched_weights.tsv"
    ewd.to_csv(ewdf, sep="\t", index=False)
    print(f"[wfilt] {ewf} ({len(ew_w)} W re-extract); "
          f"{ewdf} ({len(ew_w)} W + {len(ew_d)} novel Desikan); "
          f"{len(dropped)} dropped: {dropped}")

    # ── consolidated FINAL filtered SNP list (9 W + Bellenguez), w/ genes ───
    hw = {r["rsID"]: r for r in ew_w}
    fin = []
    for _, w in rep[rep["status"].isin(["resolved_genotyped",
                                        "on_array_not_in_430"])].iterrows():
        rid = w["pipeline_rsID"] or w["array_rsID"]
        bA1 = (hw[rid]["beta_A1"] if rid in hw else w["beta"])
        fin.append({"rsID": rid, "source": "Wightman2021_lead",
                    "gene": w["Gene"], "CHR": w["chr"],
                    "BP_GRCh38": w["pos_GRCh38"], "beta_A1": bA1,
                    "beta_source": w["beta_source"],
                    "origin": ("in_430" if w["status"] == "resolved_genotyped"
                               else "re_extracted")})
    for _, b in bell_rows.iterrows():
        fin.append({"rsID": b["SNP"], "source": "Bellenguez2022",
                    "gene": _gene(b["SNP"]), "CHR": str(b["CHR"]),
                    "BP_GRCh38": int(b["BP"]), "beta_A1": "",
                    "beta_source": "build_beta_table(parse_bellenguez)",
                    "origin": "bellenguez_in_430"})
    finf = args.out / "filtered_prs_snps.tsv"
    findf = pd.DataFrame(fin)
    _cn = pd.to_numeric(findf["CHR"].astype(str).str.replace("X", "23")
                        .str.replace("Y", "24").str.replace("MT", "25"),
                        errors="coerce")
    _bp = pd.to_numeric(findf["BP_GRCh38"], errors="coerce")
    findf = (findf.assign(_c=_cn, _b=_bp)
             .sort_values(["_c", "_b"]).drop(columns=["_c", "_b"]))
    findf.to_csv(finf, sep="\t", index=False)
    nw = sum(1 for r in fin if r["source"] != "Bellenguez2022")
    print(f"[out] {finf}  (FINAL: {nw} Wightman + {len(fin)-nw} "
          f"Bellenguez = {len(fin)} SNPs)")

    if len(dd):
        ddf = args.out / "desikan_nonoverlap_report.tsv"
        dd.to_csv(ddf, sep="\t", index=False)
        print(f"[out] {ddf}  ({len(dd)} genotyped Desikan; "
              f"{len(desikan5)} additional_nonoverlap: {desikan5})")

    print(f"\n[NEXT] filtered:    --wightman-filtered {afile} "
          f"--enriched-weights {ewf} --enriched-dosage "
          f"{wdir/'patient_dosage.tsv'}")
    print(f"[NEXT] withdesikan: --wightman-filtered {afile} "
          f"--enriched-weights {ewdf} --enriched-dosage "
          f"{wdir/'patient_dosage.tsv'}  (+{len(ew_d)} novel Desikan)")


if __name__ == "__main__":
    main()
