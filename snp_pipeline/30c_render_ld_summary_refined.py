r"""
30c_render_ld_summary_refined.py   (LOCAL, env: snp)
====================================================
S1b refinement (user-directed): re-render `ld_summary.txt` from the
existing `pairwise_r2.tsv` (script 30) with:
  • partners flagged against the union of ALL 4 PRS sources (Bellenguez
    27 + Wightman 9 + Desikan 9 + Kunkle 3) + the 430 GW-sig set;
  • top-3 partners by r² per query with rsIDs + source-membership tags;
  • 4×4 cross-source overlap matrix (rsID + positional within 1 Mb).

No new PLINK runs — purely a re-render. Updates the same output file
`genetic_baselines_ld/ld_summary.txt` and adds:
  `ld_summary_topK.tsv`        machine-readable top-K partners with tags
  `cross_source_overlap.tsv`   pairwise rsID & positional overlap counts

Usage:
  conda run -n snp python snp_pipeline/30c_render_ld_summary_refined.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
LD_DIR = BASE / "genetic_baselines_ld"
LABELS = BASE / "bmfm_inputs" / "external_gwas_labels.tsv"
GBF = BASE / "genetic_baselines_filtered"
GB = BASE / "genetic_baselines"

APOE_CHR, APOE_LO, APOE_HI = "19", 44_400_000, 46_500_000
TOP_K = 3


def _gather_sources() -> dict[str, dict[str, dict]]:
    """rsID → source meta for each of the 4 PRS sources + 430-set.
    Returns {source: {rsID: {CHR, BP, gene}}}."""
    out: dict[str, dict[str, dict]] = {
        "Bellenguez": {}, "Wightman": {}, "Desikan": {}, "Kunkle": {},
        "in_430": {}}
    lab = pd.read_csv(LABELS, sep="\t", low_memory=False)
    lab = lab[lab["label"] == 1]
    # 430-set
    for _, r in lab.iterrows():
        out["in_430"][str(r["SNP"])] = {"CHR": str(r["CHR"]),
                                          "BP": int(r["BP"])
                                          if pd.notna(r["BP"]) else 0,
                                          "gene": ""}
    # Bellenguez (subset of in_430)
    b = lab[lab["source"].astype(str).str.contains("Bellenguez")]
    for _, r in b.iterrows():
        out["Bellenguez"][str(r["SNP"])] = out["in_430"][str(r["SNP"])]
    # Wightman 9
    allow = pd.read_csv(GBF / "wightman_resolved_allowlist.tsv", sep="\t")
    for _, r in allow.iterrows():
        rs = str(r["pipeline_rsID"])
        meta = out["in_430"].get(rs, {"CHR": "", "BP": 0, "gene": ""})
        out["Wightman"][rs] = meta
    re_ex = pd.read_csv(GBF / "wfilt_enriched_weights.tsv", sep="\t")
    re_ex = re_ex[(re_ex["novel"] == True)
                   & re_ex["source"].astype(str).str.contains("Wightman")]
    for _, r in re_ex.iterrows():
        out["Wightman"][str(r["rsID"])] = {"CHR": str(r["CHR"]),
                                            "BP": int(float(r["BP"])),
                                            "gene": str(r.get("Gene", ""))}
    # Desikan 9
    des = pd.read_csv(GB / "desikan_pgs_resolved.tsv", sep="\t")
    des = des[des["beta_A1"].notna() & (des["in_430"] == True)]
    for _, r in des.iterrows():
        out["Desikan"][str(r["rsID"])] = {"CHR": str(r["CHR"]),
                                           "BP": int(r["BP"]) if pd.notna(r["BP"]) else 0,
                                           "gene": ""}
    de_novel = pd.read_csv(GBF / "wfilt_desikan_enriched_weights.tsv",
                            sep="\t")
    de_novel = de_novel[(de_novel["novel"] == True)
                          & de_novel["source"].astype(str)
                          .str.contains("Desikan")]
    for _, r in de_novel.iterrows():
        out["Desikan"][str(r["rsID"])] = {"CHR": str(r["CHR"]),
                                           "BP": int(float(r["BP"])),
                                           "gene": str(r.get("Gene", ""))}
    # Kunkle 3 (genotyped)
    out["Kunkle"]["rs7920721"] = {"CHR": "10", "BP": 11678309, "gene": "ECHDC3"}
    out["Kunkle"]["rs593742"] = {"CHR": "15", "BP": 58753575, "gene": ""}
    out["Kunkle"]["rs2830500"] = {"CHR": "21", "BP": 26784537, "gene": "ADAMTS1"}
    return out


def _tags_for(rs: str, src: dict[str, dict[str, dict]]) -> str:
    t = []
    if rs in src["Bellenguez"]:
        t.append("B")
    if rs in src["Wightman"]:
        t.append("W")
    if rs in src["Desikan"]:
        t.append("D")
    if rs in src["Kunkle"]:
        t.append("K")
    if rs in src["in_430"] and not (t):                    # only-430
        t.append("in_430")
    return "/".join(t) if t else "-"


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ld-dir", type=Path, default=LD_DIR)
    args = ap.parse_args()

    pwr = pd.read_csv(args.ld_dir / "pairwise_r2.tsv", sep="\t")
    print(f"[in ] {args.ld_dir/'pairwise_r2.tsv'}  ({len(pwr)} pairs)")
    src = _gather_sources()
    print(f"[src] Bellenguez={len(src['Bellenguez'])}  "
          f"Wightman={len(src['Wightman'])}  Desikan={len(src['Desikan'])}  "
          f"Kunkle={len(src['Kunkle'])}  in_430={len(src['in_430'])}")

    pwr["partner_tags"] = pwr["partner_rsID"].map(lambda r: _tags_for(r, src))
    pwr["partner_in_PRS_source"] = pwr["partner_tags"].apply(
        lambda t: any(x in t for x in ("B", "W", "D", "K")))

    # Refined summary
    queries = pwr["query_rsID"].drop_duplicates().tolist()
    lines = [
        "REFINED LD r² summary — pairwise r² in 1 Mb (in-cohort MAF>0.01 BED,",
        "  616 ADNI ~European samples). Re-render of `pairwise_r2.tsv` with",
        "  multi-source flagging.",
        f"APOE cluster: chr{APOE_CHR}:{APOE_LO:,}–{APOE_HI:,}.",
        f"Partner tags: B=Bellenguez, W=Wightman, D=Desikan, K=Kunkle, "
        f"in_430=label==1 GW-sig (broader set), '-'=anonymous LD-block",
        f"  variant (not in any PRS source).",
        "",
    ]
    topk_rows = []
    for q in queries:
        sub = pwr[pwr["query_rsID"] == q].sort_values("r2", ascending=False)
        n_part = len(sub)
        max_r2 = float(sub["r2"].max()) if n_part else 0.0
        # tagged partners (in any PRS source)
        tagged = sub[sub["partner_in_PRS_source"]]
        in430 = sub[sub["partner_tags"].str.contains("in_430")]
        gt8 = sub[sub["flag_r2_gt_080"]]
        gt9 = sub[sub["flag_r2_gt_090"]]
        # top-K with tags
        top = sub.head(TOP_K)
        q_src = sub["query_source"].iloc[0] if n_part else "?"
        lines += [
            f"## {q}  ({q_src})",
            f"  n_partners (1 Mb)              : {n_part}",
            f"  max r² (any partner)           : {max_r2:.4f}",
            f"  partners in ANY PRS source     : {len(tagged)}  "
            f"max r²={tagged['r2'].max() if len(tagged) else 0:.4f}",
            f"  partners in 430-set (broader)  : {len(in430)}  "
            f"max r²={in430['r2'].max() if len(in430) else 0:.4f}",
            f"  partners at r²>0.8 (any)       : {len(gt8)}",
            f"    of those, in any PRS source  : "
            f"{int(gt8['partner_in_PRS_source'].sum())}",
            f"    of those, in 430-set         : "
            f"{int(gt8['partner_tags'].str.contains('in_430').sum())}",
            f"  partners at r²>0.9             : {len(gt9)}",
            f"  top-{TOP_K} partners by r²:",
        ]
        for _, r in top.iterrows():
            lines.append(
                f"    {r['partner_rsID']:<14} tags={r['partner_tags']:<12}"
                f" dist={int(r['dist_bp']):>9} bp   r²={r['r2']:.4f}")
            topk_rows.append({
                "query_rsID": q, "query_source": q_src,
                "rank": int(top["partner_rsID"].tolist()
                              .index(r["partner_rsID"])) + 1,
                "partner_rsID": r["partner_rsID"],
                "partner_tags": r["partner_tags"],
                "partner_CHR": r["partner_CHR"],
                "partner_BP": int(r["partner_BP"]),
                "dist_bp": int(r["dist_bp"]),
                "r2": float(r["r2"]),
                "flag_r2_gt_080": bool(r["flag_r2_gt_080"]),
                "flag_r2_gt_090": bool(r["flag_r2_gt_090"])})
        # r²>0.8 partners that hit a PRS source / 430-set (the interesting ones)
        hot = gt8[gt8["partner_in_PRS_source"]
                   | gt8["partner_tags"].str.contains("in_430")]
        if len(hot):
            lines.append(f"  ⚠ r²>0.8 partners hitting PRS-source or 430-set "
                         f"({len(hot)}):")
            for _, r in hot.iterrows():
                lines.append(
                    f"      {r['partner_rsID']:<14} "
                    f"tags={r['partner_tags']:<12} dist={int(r['dist_bp']):>9}"
                    f" r²={r['r2']:.4f}")
        lines.append("")

    # 4×4 cross-source overlap matrix (rsID + positional within 1 Mb)
    sources = ["Bellenguez", "Wightman", "Desikan", "Kunkle"]
    rsid_mat = pd.DataFrame(0, index=sources, columns=sources)
    pos_mat = pd.DataFrame(0, index=sources, columns=sources)
    for sa in sources:
        for sb in sources:
            sa_rs = set(src[sa])
            sb_rs = set(src[sb])
            rsid_mat.loc[sa, sb] = len(sa_rs & sb_rs)
            # positional within 1 Mb (any rsID in sa within 1Mb of any rsID in sb)
            n_pos = 0
            for ra, ma in src[sa].items():
                if ra in sb_rs:
                    continue
                for rb, mb in src[sb].items():
                    if (str(ma["CHR"]) == str(mb["CHR"])
                            and ma["CHR"] != "" and ma["BP"] > 0
                            and mb["BP"] > 0
                            and abs(int(ma["BP"]) - int(mb["BP"])) <= 1_000_000):
                        n_pos += 1
                        break
            pos_mat.loc[sa, sb] = n_pos

    lines += ["---", "## 4×4 cross-source overlap matrices", "",
              "### rsID exact overlap (diagonal = own count)",
              rsid_mat.to_string(),
              "",
              "### positional overlap (same chr, |dist|≤1 Mb; off-diagonal "
              "only; excludes rsID-identical pairs to count distinct loci)",
              pos_mat.to_string(),
              ""]

    (args.ld_dir / "ld_summary.txt").write_text("\n".join(lines), "utf-8")
    pd.DataFrame(topk_rows).to_csv(
        args.ld_dir / "ld_summary_topK.tsv", sep="\t", index=False)
    # cross_source_overlap
    ov = []
    for sa in sources:
        for sb in sources:
            ov.append({"sourceA": sa, "sourceB": sb,
                       "rsID_overlap": int(rsid_mat.loc[sa, sb]),
                       "positional_within_1Mb": int(pos_mat.loc[sa, sb])})
    pd.DataFrame(ov).to_csv(args.ld_dir / "cross_source_overlap.tsv",
                              sep="\t", index=False)
    print(f"[out] {args.ld_dir/'ld_summary.txt'}  (refined)")
    print(f"[out] {args.ld_dir/'ld_summary_topK.tsv'}  ({len(topk_rows)} rows)")
    print(f"[out] {args.ld_dir/'cross_source_overlap.tsv'}")

    print("\n=== quick verdict ===")
    n_hot_q = pwr.groupby("query_rsID").apply(
        lambda g: int(((g["flag_r2_gt_080"]) &
                        (g["partner_in_PRS_source"]
                         | g["partner_tags"].str.contains("in_430"))).any()),
        include_groups=False).sum()
    print(f"  query SNPs with ≥1 r²>0.8 partner in PRS-source OR 430-set: "
          f"{n_hot_q} / {len(queries)}")


if __name__ == "__main__":
    main()
