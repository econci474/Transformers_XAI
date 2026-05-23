r"""
30e_pairwise_ld_missingness_resolved.py   (LOCAL, env: snp)
===========================================================
Pairwise LD r² (1 Mb window, threshold 0.8) on the **missingness-resolved**
BED produced by `25e_export_missingness_resolved_sets.py`. Mirrors the
`30d` script in structure but reads the per-source rsID sets from
`GWAS_filtered_missingness_resolved_geno005/<sourceX>/<sourceX>_snps.tsv`
and uses the new --geno 0.05 BED.

Reports:
  - Every intra-chromosomal pair within 1 Mb (no filter at PLINK time;
    threshold applied at the summary step).
  - r²>0.8 hot list (intra- and inter-source).
  - r²>0.6 moderate list.
  - Per-(source × source) cell counts and max r².
  - Per-source intra-LD totals.

Output → `D:\ADNI_SNP_Omni2.5M_20140220\genetic_baselines_ld_missingness_resolved\`:
  pairwise_r2_resolved.tsv    one row per intra-chr pair
  ld_summary_resolved.txt     human-readable summary

Usage:
  conda run -n snp python snp_pipeline/30e_pairwise_ld_missingness_resolved.py
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
BED = BASE / "SNP_filtered_with_mri_geno005_PRS_only_GRCh38_maf001"
PLINK = BASE / "plink.exe"
SOURCE_ROOT = BASE / "GWAS_filtered_missingness_resolved_geno005"
OUT_DEFAULT = BASE / "genetic_baselines_ld_missingness_resolved"

SOURCES = ("Bellenguez", "Wightman", "Desikan", "Kunkle")
LETTERS = {"Bellenguez": "B", "Wightman": "W",
           "Desikan": "D", "Kunkle": "K"}
SOURCE_DIRS = {"Bellenguez": "sourceB", "Wightman": "sourceW",
               "Desikan": "sourceD", "Kunkle": "sourceK"}

APOE_CHR, APOE_LO, APOE_HI = "19", 44_400_000, 46_500_000


def _in_apoe(chrom: str, bp: int) -> bool:
    return str(chrom) == APOE_CHR and APOE_LO <= int(bp) <= APOE_HI


def _load_source_rsids() -> tuple[dict[str, set[str]], dict[str, dict]]:
    src_rs: dict[str, set[str]] = {}
    meta: dict[str, dict] = {}
    for src in SOURCES:
        p = SOURCE_ROOT / SOURCE_DIRS[src] / f"{SOURCE_DIRS[src]}_snps.tsv"
        if not p.exists():
            print(f"[warn] {p} missing")
            src_rs[src] = set()
            continue
        df = pd.read_csv(p, sep="\t", dtype=str)
        src_rs[src] = set(df["rsID"].astype(str))
        for _, r in df.iterrows():
            rs = str(r["rsID"])
            if rs in meta:
                meta[rs]["sources"].add(src)
            else:
                meta[rs] = {
                    "sources": {src},
                    "CHR": str(r["CHR"]),
                    "BP": int(r["BP_GRCh38"]),
                    "gene": str(r.get("gene", "")) or "",
                }
    return src_rs, meta


def _tags(rs: str, meta: dict[str, dict]) -> str:
    if rs not in meta:
        return "-"
    return "/".join(sorted("BWDK"[SOURCES.index(s)] for s in meta[rs]["sources"]))


def _has(tag: str, letter: str) -> bool:
    return letter in tag.split("/")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bfile", type=Path, default=BED)
    ap.add_argument("--plink", type=Path, default=PLINK)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--window-kb", type=int, default=1000,
                    help="--ld-window-kb (default 1000 = 1 Mb)")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    src_rs, meta = _load_source_rsids()
    union = sorted({rs for st in src_rs.values() for rs in st})
    print(f"[resolved] per-source counts: "
          + "  ".join(f"{s}={len(src_rs[s])}" for s in SOURCES))
    print(f"[resolved] union (dedup): {len(union)} rsIDs")

    tmp = Path(tempfile.mkdtemp(prefix="ld_geno005_"))
    extract_path = tmp / "resolved_rsids.txt"
    extract_path.write_text("\n".join(union) + "\n", "utf-8")

    print(f"[plink] pairwise --r2 within 1 Mb on {args.bfile.name}")
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
    # PLINK .ld columns: CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
    df["snp_A"] = df["SNP_A"].astype(str)
    df["snp_B"] = df["SNP_B"].astype(str)
    df = df[df["snp_A"] != df["snp_B"]].copy()
    # canonical order so (A,B) == (B,A) appears once
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
    df["A_in_apoe"] = df.apply(lambda r: _in_apoe(r["chr_A"], r["bp_A"]), axis=1)
    df["B_in_apoe"] = df.apply(lambda r: _in_apoe(r["chr_B"], r["bp_B"]), axis=1)
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

    # ── Per-(source × source) cell counts + max r² ───────────────────────
    cell_n = pd.DataFrame(0, index=SOURCES, columns=SOURCES)
    cell_maxr2 = pd.DataFrame(0.0, index=SOURCES, columns=SOURCES)
    for _, row in df.iterrows():
        ta, tb = row["tags_A"], row["tags_B"]
        r2 = float(row["r2"])
        for sa in SOURCES:
            for sb in SOURCES:
                if _has(ta, LETTERS[sa]) and _has(tb, LETTERS[sb]):
                    cell_n.loc[sa, sb] += 1
                    if r2 > cell_maxr2.loc[sa, sb]:
                        cell_maxr2.loc[sa, sb] = r2
                if sa != sb and _has(ta, LETTERS[sb]) and _has(tb, LETTERS[sa]):
                    cell_n.loc[sa, sb] += 1
                    if r2 > cell_maxr2.loc[sa, sb]:
                        cell_maxr2.loc[sa, sb] = r2

    # ── Human-readable summary ──────────────────────────────────────────
    lines = [
        "Pairwise LD r² — missingness-resolved (--geno 0.05) PRS bundle",
        f"Cohort: 616 ADNI ~European samples on geno005 GRCh38 BIM",
        f"Per-source counts:  "
        + "  ".join(f"{s}={len(src_rs[s])}" for s in SOURCES),
        f"Union (dedup): {len(union)} rsIDs.",
        f"Window: 1 Mb (user-set); --ld-window-r2 0 (every pair).",
        f"Intra-chromosomal pairs reported by PLINK: {len(df)}.",
        f"APOE cluster: chr{APOE_CHR}:{APOE_LO:,}–{APOE_HI:,}.",
        "Tags: B=Bellenguez W=Wightman D=Desikan K=Kunkle "
        "(concatenated with '/' if rsID is in multiple sources).",
        "",
        "=" * 70,
        "## Top pairs at r² > 0.8",
    ]
    hot = df[df["flag_r2_gt_080"]].sort_values("r2", ascending=False)
    if len(hot) == 0:
        lines.append("  (none)")
    else:
        for _, rr in hot.iterrows():
            kind = "cross-source" if rr["cross_source"] else "intra-source"
            lines.append(
                f"  {rr['snp_A']:<14} [{rr['tags_A']:<4}] chr{rr['chr_A']:<2}:"
                f"{rr['bp_A']:>10}   ↔   {rr['snp_B']:<14} [{rr['tags_B']:<4}]"
                f" chr{rr['chr_B']:<2}:{rr['bp_B']:>10}   dist="
                f"{int(rr['dist_bp']):>9} bp   r²={rr['r2']:.4f}   ({kind})")
    lines.append("")
    lines.append("## Pairs at 0.6 < r² ≤ 0.8 (moderate, for context)")
    mid = df[(df["r2"] > 0.6) & ~df["flag_r2_gt_080"]].sort_values(
        "r2", ascending=False)
    if len(mid) == 0:
        lines.append("  (none)")
    else:
        for _, rr in mid.iterrows():
            kind = "cross-source" if rr["cross_source"] else "intra-source"
            lines.append(
                f"  {rr['snp_A']:<14} [{rr['tags_A']:<4}]   ↔   "
                f"{rr['snp_B']:<14} [{rr['tags_B']:<4}]   "
                f"dist={int(rr['dist_bp']):>9} bp   r²={rr['r2']:.4f}   ({kind})")
    lines.append("")
    lines.append("=" * 70)
    lines.append("## Per-(sourceA × sourceB) cell counts  (bidirectional)")
    lines.append(cell_n.to_string())
    lines.append("")
    lines.append("## Max r² per cell")
    lines.append(cell_maxr2.round(4).to_string())
    lines.append("")
    lines.append("=" * 70)
    lines.append("## Per-source intra-LD summary  (pairs where BOTH rsIDs map "
                  "to the source)")
    for src in SOURCES:
        intra = df[df["tags_A"].str.contains(LETTERS[src])
                    & df["tags_B"].str.contains(LETTERS[src])]
        n_pairs = len(intra)
        n_gt8 = int(intra["flag_r2_gt_080"].sum())
        max_r2 = float(intra["r2"].max()) if n_pairs else 0.0
        lines.append(f"  {src:<12} n_intra_pairs={n_pairs:<5} "
                     f"r²>0.8={n_gt8}  max r²={max_r2:.4f}")
    lines.append("")
    lines.append("=" * 70)
    lines.append("Note: source membership is post-dedup-union (B > D > W > K "
                  "priority); a single rsID surfaces once per source it belongs "
                  "to, joined with '/' in the tag column.")

    (args.out / "ld_summary_resolved.txt").write_text("\n".join(lines), "utf-8")
    print(f"[out] {args.out/'ld_summary_resolved.txt'}")

    print(f"\n=== verdict ===")
    print(f"  Pairs at r²>0.8 (any tag combination): "
          f"{int(df['flag_r2_gt_080'].sum())}")
    print(f"  Pairs at r²>0.9: {int(df['flag_r2_gt_090'].sum())}")
    print(f"  Cross-source pairs at r²>0.8: "
          f"{int((df['flag_r2_gt_080'] & df['cross_source']).sum())}")


if __name__ == "__main__":
    main()
