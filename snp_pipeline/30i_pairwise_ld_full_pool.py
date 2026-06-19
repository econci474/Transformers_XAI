r"""
30i_pairwise_ld_full_pool.py   (LOCAL, env: snp)
================================================
Pairwise LD r² (1 Mb window) on the full **22-PRS-source, 116-on-MAF**
SNP pool produced by 30c (post v2+v3 extension). Mirrors 30e but:
- 22 sources instead of 4 (uses full source names, not BWDK letters)
- Standard MAF BIM (not the missingness-resolved geno005 BED)
- Adds per-phenoclass LD aggregation
- Generates a thesis-grade PNG heatmap

Outputs → `D:\…\GWAS_comprehensive_v2\LD_report\`:
  pairwise_r2.tsv                 one row per intra-chr pair
  ld_summary.txt                  human-readable summary
  ld_clusters.tsv                 connected components at r²>0.8
  per_phenoclass_ld_matrix.tsv    phenoclass×phenoclass pair counts
  ld_thesis_heatmap.png           22×22 source heatmap @ r²>0.8

Usage:
  conda run -n snp python snp_pipeline/30i_pairwise_ld_full_pool.py
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
MAF_BIM = BASE / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr_maf001"
RECON = BASE / "source_prs" / "unfiltered_SNP_reconciliation"
PLINK = BASE / "plink.exe"
OUT_DEFAULT = BASE / "GWAS_comprehensive_v2" / "LD_report"

PRS_SOURCES = [
    # 8 baseline (existing + v1)
    "Bellenguez", "Wightman", "Kunkle", "Desikan",
    "Lambert", "DeRojas", "Schwanzentruber", "Najar",
    # 8 v2 (2026-05-25 evening)
    "Ebenau", "Leonenko", "Vesilievick", "Zhang",
    "Felsky_MF", "Felsky_IT", "ONeil_NPY", "ONeil_GHR",
    "Kosteridis_novel_AD", "Kosteridis_shared_AD_CV",
    # 6 v3 — Yang (2026-05-25 night)
    "Yang_Ex", "Yang_In", "Yang_Ast",
    "Yang_Mic", "Yang_Oli", "Yang_Opc",
]

PHENO_TYPE = {
    "clinical_AD":                 "clinical",
    "clinical_AD_MTAG":            "clinical",
    "clinical_AAO":                "clinical",
    "pathology_microglia":         "pathology",
    "pathology_amyloid":           "pathology",
    "pathology_plaque":            "pathology",
    "pathology_mixed":             "pathology",
    "pathology_astrocyte":         "pathology",
    "pathology_neuron_excitatory": "pathology",
    "pathology_neuron_inhibitory": "pathology",
    "pathology_oligodendrocyte":   "pathology",
    "pathology_opc":               "pathology",
}

# Priority order for multi-class rsID main_class selection (matches 30f)
_PHENO_PRIORITY = ["clinical_AD", "clinical_AD_MTAG",
                   "pathology_microglia",
                   "pathology_astrocyte",
                   "pathology_neuron_excitatory",
                   "pathology_neuron_inhibitory",
                   "pathology_oligodendrocyte",
                   "pathology_opc",
                   "pathology_plaque", "pathology_amyloid",
                   "clinical_AAO", "pathology_mixed"]

APOE_CHR, APOE_LO, APOE_HI = "19", 44_400_000, 46_500_000


def _in_apoe(chrom: str, bp: int) -> bool:
    return str(chrom) == APOE_CHR and APOE_LO <= int(bp) <= APOE_HI


def _pick_main_class(classes: set[str]) -> str:
    for p in _PHENO_PRIORITY:
        if p in classes:
            return p
    return next(iter(classes), "")


def _load_on_maf_pool() -> dict[str, dict]:
    """{rsID: {sources: set, pheno_classes: set, chr, bp}} for the
    116-SNP union of on-MAF rsIDs across the 22 PRS sources."""
    pool: dict[str, dict] = defaultdict(
        lambda: {"sources": set(), "pheno_classes": set(),
                  "chr": "", "bp": ""})
    for src in PRS_SOURCES:
        p = RECON / f"{src}_unfiltered_reconciliation.tsv"
        if not p.exists():
            print(f"  [warn] missing {p}")
            continue
        df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
        df["in_maf_bim"] = df["in_maf_bim"].str.lower().isin({"true", "1"})
        df = df[df["in_maf_bim"]]
        pheno = df["pheno_class"].iloc[0] if len(df) and "pheno_class" in df.columns else ""
        for _, r in df.iterrows():
            rs = str(r["rsID_pub"])
            pool[rs]["sources"].add(src)
            if pheno:
                pool[rs]["pheno_classes"].add(pheno)
            # use unfilt CHR/BP if present, else CHR_pub/BP_pub
            if not pool[rs]["chr"]:
                pool[rs]["chr"] = (str(r.get("unfilt_CHR_grch37", ""))
                                    or str(r.get("CHR_pub", "")))
            if not pool[rs]["bp"]:
                pool[rs]["bp"] = (str(r.get("unfilt_BP_grch37", ""))
                                   or str(r.get("BP_pub", "")))
    return dict(pool)


def _connected_components(edges: list[tuple[str, str]],
                            nodes: set[str]) -> list[set[str]]:
    """Union-find: build connected components from a list of edges."""
    parent = {n: n for n in nodes}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    for a, b in edges:
        if a in parent and b in parent: union(a, b)
    comps: dict[str, set[str]] = defaultdict(set)
    for n in nodes:
        comps[find(n)].add(n)
    return [c for c in comps.values() if len(c) >= 2]


def _render_heatmap(pair_matrix: pd.DataFrame, png_path: Path,
                     pool: dict[str, dict]) -> None:
    """22×22 heatmap of per-(srcA × srcB) pair counts at r²>0.8."""
    fig, ax = plt.subplots(figsize=(11, 10), dpi=300)
    data = pair_matrix.fillna(0).astype(int).values
    im = ax.imshow(data, cmap="OrRd", aspect="auto")
    ax.set_xticks(range(len(PRS_SOURCES)))
    ax.set_yticks(range(len(PRS_SOURCES)))
    ax.set_xticklabels(PRS_SOURCES, rotation=70, ha="right", fontsize=8)
    ax.set_yticklabels(PRS_SOURCES, fontsize=8)
    # Show cell counts (skip zeros)
    for i in range(len(PRS_SOURCES)):
        for j in range(len(PRS_SOURCES)):
            v = data[i, j]
            if v > 0:
                col = "white" if v >= data.max() * 0.6 else "black"
                ax.text(j, i, str(v), ha="center", va="center",
                         fontsize=7, color=col)
    cb = plt.colorbar(im, ax=ax, label="n pairs at r² > 0.8")
    ax.set_title("Cross-source LD on ADNI Omni2.5M MAF chip\n"
                  "(22 PRS sources, 116 unique on-MAF SNPs, r² > 0.8, "
                  "1 Mb window)", fontsize=12, pad=12)
    fig.text(0.04, 0.02,
              "Diagonal = intra-source LD-redundant pairs; "
              "off-diagonal = cross-source pairs at same locus. "
              "Counts are pairs of rsIDs (bidirectional); a rsID "
              "shared by 2 sources is NOT counted (it's one node, not "
              "a pair). APOE cluster (chr19:44.4–46.5 Mb) dominates "
              "most off-diagonal cells.",
              fontsize=7, ha="left", va="bottom", wrap=True)
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bfile", type=Path, default=MAF_BIM)
    ap.add_argument("--plink", type=Path, default=PLINK)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--window-kb", type=int, default=1000)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    pool = _load_on_maf_pool()
    union = sorted(pool)
    print(f"[pool] {len(union)} unique on-MAF rsIDs across "
          f"{len(PRS_SOURCES)} PRS sources")

    # Per-source counts
    per_src_n = {s: sum(1 for r in pool.values() if s in r["sources"])
                  for s in PRS_SOURCES}
    print("[pool] per-source on-MAF counts:")
    for s in PRS_SOURCES:
        print(f"  {s:24s}: {per_src_n[s]:3d}")

    tmp = Path(tempfile.mkdtemp(prefix="ld_full_pool_"))
    extract_path = tmp / "rsids.txt"
    extract_path.write_text("\n".join(union) + "\n", "utf-8")

    print(f"[plink] pairwise --r2 within {args.window_kb} kb")
    out_prefix = tmp / "pairwise"
    cmd = [str(args.plink),
           "--bfile", str(args.bfile),
           "--extract", str(extract_path),
           "--r2",
           "--ld-window", "99999",
           "--ld-window-kb", str(args.window_kb),
           "--ld-window-r2", "0",
           "--out", str(out_prefix)]
    print("  " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("[plink] STDERR:", r.stderr[-2000:])
        sys.exit(r.returncode)
    ld_file = out_prefix.with_suffix(".ld")
    if not ld_file.exists():
        sys.exit(f"[ERROR] no .ld file at {ld_file}")

    df = pd.read_csv(ld_file, sep=r"\s+")
    df["snp_A"] = df["SNP_A"].astype(str)
    df["snp_B"] = df["SNP_B"].astype(str)
    df = df[df["snp_A"] != df["snp_B"]].copy()
    swap = df["snp_A"] > df["snp_B"]
    df.loc[swap, ["snp_A", "snp_B", "CHR_A", "BP_A", "CHR_B", "BP_B"]] = \
        df.loc[swap, ["snp_B", "snp_A", "CHR_B", "BP_B", "CHR_A", "BP_A"]].values
    df = df.drop_duplicates(subset=["snp_A", "snp_B"])
    df["chr_A"] = df["CHR_A"].astype(str)
    df["bp_A"] = df["BP_A"].astype(int)
    df["chr_B"] = df["CHR_B"].astype(str)
    df["bp_B"] = df["BP_B"].astype(int)
    df["dist_bp"] = (df["bp_B"] - df["bp_A"]).abs().astype(int)
    df["r2"] = df["R2"].astype(float)
    df["sources_A"] = df["snp_A"].map(
        lambda r: "/".join(sorted(pool.get(r, {"sources": set()})["sources"])))
    df["sources_B"] = df["snp_B"].map(
        lambda r: "/".join(sorted(pool.get(r, {"sources": set()})["sources"])))
    df["pheno_classes_A"] = df["snp_A"].map(
        lambda r: "|".join(sorted(pool.get(r, {"pheno_classes": set()})["pheno_classes"])))
    df["pheno_classes_B"] = df["snp_B"].map(
        lambda r: "|".join(sorted(pool.get(r, {"pheno_classes": set()})["pheno_classes"])))
    df["A_in_apoe"] = df.apply(lambda r: _in_apoe(r["chr_A"], r["bp_A"]), axis=1)
    df["B_in_apoe"] = df.apply(lambda r: _in_apoe(r["chr_B"], r["bp_B"]), axis=1)
    df["flag_r2_gt_080"] = df["r2"] > 0.8
    df["flag_r2_gt_060"] = (df["r2"] > 0.6) & ~df["flag_r2_gt_080"]
    df["same_source"] = (df["sources_A"] == df["sources_B"])
    df["cross_source"] = ~df["same_source"]
    # phenoclass aggregation
    df["main_class_A"] = df["snp_A"].map(
        lambda r: _pick_main_class(pool.get(r, {"pheno_classes": set()})["pheno_classes"]))
    df["main_class_B"] = df["snp_B"].map(
        lambda r: _pick_main_class(pool.get(r, {"pheno_classes": set()})["pheno_classes"]))
    df["same_phenoclass"] = (df["main_class_A"] == df["main_class_B"])
    df["cross_phenoclass"] = ~df["same_phenoclass"]

    out_cols = ["snp_A", "sources_A", "pheno_classes_A", "main_class_A",
                "chr_A", "bp_A", "A_in_apoe",
                "snp_B", "sources_B", "pheno_classes_B", "main_class_B",
                "chr_B", "bp_B", "B_in_apoe",
                "dist_bp", "r2", "same_source", "cross_source",
                "same_phenoclass", "cross_phenoclass",
                "flag_r2_gt_080", "flag_r2_gt_060"]
    df[out_cols].to_csv(args.out / "pairwise_r2.tsv", sep="\t", index=False)
    print(f"[out] {args.out/'pairwise_r2.tsv'}  ({len(df)} pairs)")

    # ── Per-(source × source) pair count + max r² ─────────────────────
    src_pair_n = pd.DataFrame(0, index=PRS_SOURCES, columns=PRS_SOURCES)
    src_pair_maxr2 = pd.DataFrame(0.0, index=PRS_SOURCES, columns=PRS_SOURCES)
    hot = df[df["flag_r2_gt_080"]]
    for _, rr in hot.iterrows():
        sa_set = pool[rr["snp_A"]]["sources"]
        sb_set = pool[rr["snp_B"]]["sources"]
        r2v = float(rr["r2"])
        for sa in sa_set:
            for sb in sb_set:
                src_pair_n.loc[sa, sb] += 1
                src_pair_n.loc[sb, sa] += 1
                if r2v > src_pair_maxr2.loc[sa, sb]:
                    src_pair_maxr2.loc[sa, sb] = r2v
                    src_pair_maxr2.loc[sb, sa] = r2v

    # ── Per-phenoclass aggregation ───────────────────────────────────────
    classes = sorted({c for v in pool.values() for c in v["pheno_classes"]})
    pc_pair_n = pd.DataFrame(0, index=classes, columns=classes)
    for _, rr in hot.iterrows():
        ca, cb = rr["main_class_A"], rr["main_class_B"]
        if ca in classes and cb in classes:
            pc_pair_n.loc[ca, cb] += 1
            if ca != cb:
                pc_pair_n.loc[cb, ca] += 1
    pc_pair_n.to_csv(args.out / "per_phenoclass_ld_matrix.tsv", sep="\t")
    print(f"[out] {args.out/'per_phenoclass_ld_matrix.tsv'}  "
          f"({len(classes)} phenoclasses)")

    # ── Connected-component clusters at r²>0.8 ───────────────────────────
    edges = list(zip(hot["snp_A"], hot["snp_B"]))
    clusters = _connected_components(edges, set(union))
    cluster_rows = []
    for cid, comp in enumerate(sorted(clusters, key=lambda c: -len(c)), start=1):
        comp_sorted = sorted(comp)
        srcs = set()
        phenos = set()
        chroms, bps = set(), []
        for rs in comp:
            srcs |= pool[rs]["sources"]
            phenos |= pool[rs]["pheno_classes"]
            chroms.add(pool[rs]["chr"])
            try:
                bps.append(int(pool[rs]["bp"]))
            except (ValueError, TypeError):
                pass
        hot_in_comp = hot[hot["snp_A"].isin(comp) & hot["snp_B"].isin(comp)]
        max_r2 = float(hot_in_comp["r2"].max()) if len(hot_in_comp) else 0.0
        cluster_rows.append({
            "cluster_id": cid,
            "n_snps": len(comp),
            "snps": "/".join(comp_sorted),
            "sources": "/".join(sorted(srcs)),
            "pheno_classes": "|".join(sorted(phenos)),
            "chroms": ",".join(sorted(chroms)),
            "bp_min": min(bps) if bps else "",
            "bp_max": max(bps) if bps else "",
            "max_r2": round(max_r2, 4),
        })
    pd.DataFrame(cluster_rows).to_csv(
        args.out / "ld_clusters.tsv", sep="\t", index=False)
    print(f"[out] {args.out/'ld_clusters.tsv'}  ({len(cluster_rows)} clusters at r²>0.8)")

    # ── Human-readable summary ─────────────────────────────────────────
    lines = [
        "Pairwise LD r² — full PRS pool (22 sources, 116 on-MAF SNPs)",
        f"Cohort: 616 ADNI ~European samples on the MAF BIM ({MAF_BIM.name}).",
        f"Window: {args.window_kb} kb (--ld-window-r2 0 → every pair).",
        f"Intra-chromosomal pairs reported: {len(df)}.",
        f"Pairs at r²>0.8: {int(df['flag_r2_gt_080'].sum())}.",
        f"Pairs at 0.6<r²≤0.8: {int(df['flag_r2_gt_060'].sum())}.",
        f"Cross-source pairs at r²>0.8: "
        f"{int((df['flag_r2_gt_080'] & df['cross_source']).sum())}.",
        f"Connected-component clusters at r²>0.8: {len(cluster_rows)}.",
        f"APOE cluster: chr{APOE_CHR}:{APOE_LO:,}–{APOE_HI:,}.",
        "",
        "=" * 70,
        f"## Per-source on-MAF counts",
        "",
    ]
    for s in PRS_SOURCES:
        lines.append(f"  {s:24s}: {per_src_n[s]:>3} on-MAF")
    lines.append("")
    lines.append("=" * 70)
    lines.append("## Top pairs at r² > 0.8")
    lines.append("")
    if len(hot) == 0:
        lines.append("  (none)")
    else:
        for _, rr in hot.sort_values("r2", ascending=False).iterrows():
            kind = "cross-source" if rr["cross_source"] else "intra-source"
            apoe = " [APOE]" if rr["A_in_apoe"] and rr["B_in_apoe"] else ""
            lines.append(
                f"  {rr['snp_A']:<13} [{rr['sources_A']:<40}]  chr{rr['chr_A']:<2}:"
                f"{rr['bp_A']:>10}\n"
                f"   ↔  {rr['snp_B']:<13} [{rr['sources_B']:<40}]  chr{rr['chr_B']:<2}:"
                f"{rr['bp_B']:>10}\n"
                f"      dist={int(rr['dist_bp']):>9} bp   r²={rr['r2']:.4f}   "
                f"({kind}){apoe}")
    lines.append("")
    lines.append("## Pairs at 0.6 < r² ≤ 0.8 (moderate)")
    lines.append("")
    mid = df[df["flag_r2_gt_060"]].sort_values("r2", ascending=False)
    if len(mid) == 0:
        lines.append("  (none)")
    else:
        for _, rr in mid.iterrows():
            kind = "cross-source" if rr["cross_source"] else "intra-source"
            apoe = " [APOE]" if rr["A_in_apoe"] and rr["B_in_apoe"] else ""
            lines.append(
                f"  {rr['snp_A']:<13} ↔ {rr['snp_B']:<13}  "
                f"dist={int(rr['dist_bp']):>9} bp  r²={rr['r2']:.4f}  "
                f"({kind}){apoe}")
    lines.append("")
    lines.append("=" * 70)
    lines.append("## LD clusters at r²>0.8 (connected components)")
    lines.append("")
    for c in cluster_rows:
        lines.append(
            f"  cluster {c['cluster_id']:>2} ({c['n_snps']} SNPs, "
            f"chr{c['chroms']} {c['bp_min']:,}-{c['bp_max']:,}, "
            f"max r²={c['max_r2']:.3f}):")
        lines.append(f"     snps: {c['snps']}")
        lines.append(f"     sources: {c['sources']}")
        lines.append(f"     pheno: {c['pheno_classes']}")
    lines.append("")
    lines.append("=" * 70)
    lines.append("## Per-(phenoclass × phenoclass) pair counts at r²>0.8")
    lines.append("")
    lines.append(pc_pair_n.to_string())
    lines.append("")
    lines.append("=" * 70)
    lines.append("## Per-(source × source) pair counts at r²>0.8")
    lines.append("(see ld_thesis_heatmap.png for a visual)")
    lines.append("")
    lines.append(src_pair_n.to_string())
    lines.append("")

    (args.out / "ld_summary.txt").write_text("\n".join(lines), "utf-8")
    print(f"[out] {args.out/'ld_summary.txt'}")

    # ── PNG heatmap ─────────────────────────────────────────────────────
    png_path = args.out / "ld_thesis_heatmap.png"
    _render_heatmap(src_pair_n, png_path, pool)
    print(f"[out] {png_path}")

    print(f"\n=== verdict ===")
    print(f"  Total pairs within 1 Mb: {len(df)}")
    print(f"  Pairs at r²>0.8: {int(df['flag_r2_gt_080'].sum())}")
    print(f"  Cross-source pairs at r²>0.8: "
          f"{int((df['flag_r2_gt_080'] & df['cross_source']).sum())}")
    print(f"  LD clusters at r²>0.8: {len(cluster_rows)}")


if __name__ == "__main__":
    main()
