r"""
30d_pairwise_ld_resolved_only.py   (LOCAL, env: snp)
====================================================
S1e (user-directed): pairwise LD r² **restricted to the resolved-set
partner pool** — the union of the 4 PRS sources' rsIDs that were
recovered against the MAF>0.01 BIM (script 30b output:
`Bellenguez_full_snps_resolution.tsv` etc., `drop_reason=='ok'`).

Discards (i) the wider 430 GW-sig set and (ii) the ~1,400 anonymous
LD-block neighbours that PLINK reports by default for each query —
only PRS-source × PRS-source pairs.

PLINK: `--extract <resolved_union.txt> --r2 --ld-window 99999
--ld-window-kb 999999 --ld-window-r2 0` so every intra-chromosomal pair
within the resolved set is reported (cross-chromosome r² is undefined and
naturally excluded). Each pair annotated with source-membership tags for
both rsIDs.

Output → `D:\ADNI_SNP_Omni2.5M_20140220\genetic_baselines_ld\`:
  pairwise_r2_resolved.tsv   one row per intra-chr pair within resolved set
  ld_summary_resolved.txt    per-source intra-LD + cross-source LD-pair
                              listings; max r² per source-source cell

Usage:
  conda run -n snp python snp_pipeline/30d_pairwise_ld_resolved_only.py
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
MAF_BED = BASE / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr_maf001"
PLINK = BASE / "plink.exe"
SOURCE_DIR = BASE / "source_prs"
OUT_DIR = BASE / "genetic_baselines_ld"

SOURCES = ("Bellenguez", "Wightman", "Desikan", "Kunkle")
APOE_CHR, APOE_LO, APOE_HI = "19", 44_400_000, 46_500_000


def _load_resolved() -> tuple[dict[str, set[str]], dict[str, dict]]:
    """Per-source resolved rsID set + per-rsID metadata (gene, BP, etc.)."""
    src_rs: dict[str, set[str]] = {}
    meta: dict[str, dict] = {}
    for src in SOURCES:
        p = SOURCE_DIR / f"{src}_full_snps_resolution.tsv"
        if not p.exists():
            print(f"[warn] {p} not found")
            src_rs[src] = set()
            continue
        df = pd.read_csv(p, sep="\t")
        ok = df[df["drop_reason"] == "ok"].copy()
        src_rs[src] = set(ok["rsID_in_bim"].astype(str))
        for _, r in ok.iterrows():
            rs = str(r["rsID_in_bim"])
            if rs in meta:
                meta[rs]["sources"].add(src)
            else:
                meta[rs] = {
                    "sources": {src},
                    "CHR": str(r["CHR_GRCh38"]),
                    "BP": int(r["BP_GRCh38"]) if pd.notna(r["BP_GRCh38"])
                          else 0,
                    "gene": str(r.get("locus_name", "")) or "",
                }
    return src_rs, meta


def _tags(rs: str, meta: dict[str, dict]) -> str:
    if rs not in meta:
        return "-"
    return "/".join(sorted("BWDK"[SOURCES.index(s)] for s in meta[rs]["sources"]))


def _in_apoe(chrom: str, bp: int) -> bool:
    return str(chrom) == APOE_CHR and APOE_LO <= int(bp) <= APOE_HI


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bfile", type=Path, default=MAF_BED)
    ap.add_argument("--plink", type=Path, default=PLINK)
    ap.add_argument("--out", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    src_rs, meta = _load_resolved()
    union = sorted({rs for st in src_rs.values() for rs in st})
    print(f"[resolved] per-source counts:  "
          + "  ".join(f"{s}={len(src_rs[s])}" for s in SOURCES))
    print(f"[resolved] union (dedup): {len(union)} rsIDs")

    # write extract list
    tmp = Path(tempfile.mkdtemp(prefix="ld_resolved_"))
    extract_path = tmp / "resolved_rsids.txt"
    extract_path.write_text("\n".join(union) + "\n", "utf-8")

    print("[plink] pairwise --r2 within resolved set …")
    out_prefix = tmp / "resolved"
    cmd = [str(args.plink),
           "--bfile", str(args.bfile),
           "--extract", str(extract_path),
           "--r2",
           "--ld-window", "99999",
           "--ld-window-kb", "999999",
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
    # PLINK .ld columns: CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
    df["snp_A"] = df["SNP_A"].astype(str)
    df["snp_B"] = df["SNP_B"].astype(str)
    df = df[df["snp_A"] != df["snp_B"]].copy()       # drop self-pairs
    # canonical order so pair (A,B) == (B,A) appears once
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
    df["tags_A"] = df["snp_A"].map(lambda r: _tags(r, meta))
    df["tags_B"] = df["snp_B"].map(lambda r: _tags(r, meta))
    df["A_in_apoe"] = df.apply(lambda r: _in_apoe(r["chr_A"], r["bp_A"]),
                                axis=1)
    df["B_in_apoe"] = df.apply(lambda r: _in_apoe(r["chr_B"], r["bp_B"]),
                                axis=1)
    df["flag_r2_gt_080"] = df["r2"] > 0.8
    df["flag_r2_gt_090"] = df["r2"] > 0.9
    df["same_source"] = (df["tags_A"] == df["tags_B"])
    df["cross_source"] = ~df["same_source"]

    out_cols = ["snp_A", "tags_A", "chr_A", "bp_A", "A_in_apoe",
                "snp_B", "tags_B", "chr_B", "bp_B", "B_in_apoe",
                "dist_bp", "r2", "same_source", "cross_source",
                "flag_r2_gt_080", "flag_r2_gt_090"]
    df[out_cols].to_csv(args.out / "pairwise_r2_resolved.tsv",
                         sep="\t", index=False)
    print(f"[out] {args.out/'pairwise_r2_resolved.tsv'}  ({len(df)} pairs)")

    # ── Per-source intra-LD + cross-source LD pairs ──────────────────────
    lines = [
        "Pairwise LD r² — RESOLVED-set only (partner pool = the union of",
        "the 4 PRS sources' rsIDs after MAF>0.01 + liftover resolution).",
        f"Cohort: 616 ADNI ~European samples on MAF>0.01 BIM.",
        f"Resolved-set sizes:  "
        + "  ".join(f"{s}={len(src_rs[s])}" for s in SOURCES),
        f"Union (dedup): {len(union)} rsIDs.",
        f"Intra-chromosomal pairs reported by PLINK: {len(df)}.",
        f"APOE cluster: chr{APOE_CHR}:{APOE_LO:,}–{APOE_HI:,}.",
        "Tags: B=Bellenguez W=Wightman D=Desikan K=Kunkle (concatenated with '/'",
        "  if an rsID is in multiple sources).",
        "",
        "=" * 70,
        "## Top pairs at r² > 0.8",
    ]
    hot = df[df["flag_r2_gt_080"]].sort_values("r2", ascending=False)
    if len(hot) == 0:
        lines.append("  (none)")
    else:
        for _, r in hot.iterrows():
            kind = "cross-source" if r["cross_source"] else "intra-source"
            lines.append(
                f"  {r['snp_A']:<14} [{r['tags_A']:<4}] chr{r['chr_A']:<2}:"
                f"{r['bp_A']:>10}   ↔   {r['snp_B']:<14} [{r['tags_B']:<4}]"
                f" chr{r['chr_B']:<2}:{r['bp_B']:>10}   dist="
                f"{int(r['dist_bp']):>9} bp   r²={r['r2']:.4f}   ({kind})")
    lines.append("")
    lines.append("## Pairs at 0.6 < r² ≤ 0.8 (moderate, for context)")
    mid = df[(df["r2"] > 0.6) & ~df["flag_r2_gt_080"]
              ].sort_values("r2", ascending=False)
    if len(mid) == 0:
        lines.append("  (none)")
    else:
        for _, r in mid.iterrows():
            kind = "cross-source" if r["cross_source"] else "intra-source"
            lines.append(
                f"  {r['snp_A']:<14} [{r['tags_A']:<4}]   ↔   "
                f"{r['snp_B']:<14} [{r['tags_B']:<4}]   "
                f"dist={int(r['dist_bp']):>9} bp   r²={r['r2']:.4f}   ({kind})")
    lines.append("")
    lines.append("=" * 70)
    lines.append("## Per-(sourceA × sourceB) cell — n_pairs / max r²")
    for sa in SOURCES:
        for sb in SOURCES:
            sub = df[df["tags_A"].str.contains(SOURCES.index(sa) + 0
                                                 and 'X' or 'X')]
            # simpler: enumerate cells based on tag containment of single letters
    lines.append("")
    lines.append("  (intra-source diagonal = pairs where BOTH rsIDs map to "
                  "the same source; off-diagonal = pairs spanning two sources)")
    # build the cells with a simple containment test
    letters = {"Bellenguez": "B", "Wightman": "W",
                "Desikan": "D", "Kunkle": "K"}

    def _has(tag: str, letter: str) -> bool:
        return letter in tag.split("/")
    cell_n = pd.DataFrame(0, index=SOURCES, columns=SOURCES)
    cell_maxr2 = pd.DataFrame(0.0, index=SOURCES, columns=SOURCES)
    for _, row in df.iterrows():
        ta, tb = row["tags_A"], row["tags_B"]
        r2 = float(row["r2"])
        for sa in SOURCES:
            for sb in SOURCES:
                if _has(ta, letters[sa]) and _has(tb, letters[sb]):
                    cell_n.loc[sa, sb] += 1
                    if r2 > cell_maxr2.loc[sa, sb]:
                        cell_maxr2.loc[sa, sb] = r2
                # symmetrise off-diagonal
                if sa != sb and _has(ta, letters[sb]) and _has(tb, letters[sa]):
                    cell_n.loc[sa, sb] += 1
                    if r2 > cell_maxr2.loc[sa, sb]:
                        cell_maxr2.loc[sa, sb] = r2
    lines.append("")
    lines.append("### Pair counts (n; bidirectional counted)")
    lines.append(cell_n.to_string())
    lines.append("")
    lines.append("### Max r² per cell")
    lines.append(cell_maxr2.round(4).to_string())
    lines.append("")
    lines.append("=" * 70)
    lines.append("## Per-source intra-LD summary")
    for src in SOURCES:
        intra = df[df["tags_A"].str.contains(letters[src])
                    & df["tags_B"].str.contains(letters[src])]
        # filter to pairs where BOTH have this letter (intra-source)
        n_pairs = len(intra)
        n_gt8 = int(intra["flag_r2_gt_080"].sum())
        max_r2 = float(intra["r2"].max()) if n_pairs else 0.0
        lines.append(f"  {src:<12} n_intra_pairs={n_pairs:<4} "
                     f"r²>0.8={n_gt8}  max r²={max_r2:.4f}")
    lines.append("")
    (args.out / "ld_summary_resolved.txt").write_text("\n".join(lines),
                                                       "utf-8")
    print(f"[out] {args.out/'ld_summary_resolved.txt'}")
    print(f"\n=== verdict ===")
    print(f"  Pairs at r²>0.8 (any tag combination): {int(df['flag_r2_gt_080'].sum())}")
    print(f"  Pairs at r²>0.9: {int(df['flag_r2_gt_090'].sum())}")
    print(f"  Cross-source pairs at r²>0.8: "
          f"{int((df['flag_r2_gt_080'] & df['cross_source']).sum())}")


if __name__ == "__main__":
    main()
