r"""
30_ld_wightman_desikan_kunkle_vs_bellenguez.py   (LOCAL, env: snp)
==================================================================
Pairwise LD r² within 1 Mb between the **3 re-extracted Wightman + 5 novel
Desikan + 2 novel Kunkle** query SNPs (the ones added on top of the
27 Bellenguez leads) and every neighbouring SNP in the cohort, restricted
to **MAF>0.01** (per user: file `…_GRCh38_refcorr_maf001.{bed,bim,fam}`).
Then filter to Bellenguez partners and flag pairs at r²>0.8 / r²>0.9 —
detects whether any of our "added" SNPs are actually colinear with a
Bellenguez lead and would double-count in a combined PRS.

User-set: **fixed 1 Mb window**, **`--ld-window-r2 0`** (integer; print
every pair, no threshold). Flag `query_in_apoe_cluster` and
`partner_in_apoe_cluster` for chr19:44.4–46.5 Mb (the script-20
`APOE_CLUSTER_DEFAULT`; sanity guardrail — queries are noAPOE by
construction so the flag should be False for all 10 query SNPs).

In-cohort LD (616 samples; ADNI is predominantly European — documented
caveat).

Output → `D:\ADNI_SNP_Omni2.5M_20140220\genetic_baselines_ld\`:
  pairwise_r2.tsv   one row per (query, partner) pair within 1 Mb
  ld_summary.txt    per-query: n_partners, n_Bellenguez_partners, max_r²,
                    top-3 partners, count r²>0.8, APOE-cluster overlap.

Usage:
  conda run -n snp python snp_pipeline/30_ld_wightman_desikan_kunkle_vs_bellenguez.py
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
BFILE = BASE / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr_maf001"
PLINK = BASE / "plink.exe"
LABELS = BASE / "bmfm_inputs" / "external_gwas_labels.tsv"
OUT_DEFAULT = BASE / "genetic_baselines_ld"

# 10 query SNPs: 3 re-extracted Wightman + 5 novel Desikan + 2 novel Kunkle
QUERIES = {
    "rs4504245":  "Wightman2021_lead_re_extracted",
    "rs708382":   "Wightman2021_lead_re_extracted",
    "rs1354106":  "Wightman2021_lead_re_extracted",
    "rs1109581":  "Desikan2019_PGS_novel",
    "rs17265593": "Desikan2019_PGS_novel",
    "rs1476679":  "Desikan2019_PGS_novel",
    "rs12590273": "Desikan2019_PGS_novel",
    "rs74615166": "Desikan2019_PGS_novel",
    "rs593742":   "Kunkle2019_novel",
    "rs2830500":  "Kunkle2019_novel",
}
# APOE cluster window — same as script-20 APOE_CLUSTER_DEFAULT
APOE_CHR, APOE_LO, APOE_HI = "19", 44_400_000, 46_500_000


def _in_apoe(chrom: str, bp: int) -> bool:
    return str(chrom) == APOE_CHR and APOE_LO <= int(bp) <= APOE_HI


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bfile", type=Path, default=BFILE)
    ap.add_argument("--plink", type=Path, default=PLINK)
    ap.add_argument("--labels", type=Path, default=LABELS)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--window-kb", type=int, default=1000,
                    help="--ld-window-kb (default 1000 = 1 Mb)")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    # ── Bellenguez partner set (rsID, CHR, BP, gene) ─────────────────────
    print(f"[labels] {args.labels}")
    lab = pd.read_csv(args.labels, sep="\t", low_memory=False)
    bell = lab[(lab["label"] == 1) & lab["source"].astype(str)
               .str.contains("Bellenguez", na=False)].copy()
    bell["CHR"] = bell["CHR"].astype(str)
    bell["BP"] = pd.to_numeric(bell["BP"], errors="coerce").astype("Int64")
    bell_rs = set(bell["SNP"].astype(str))
    bell_meta = bell.set_index(bell["SNP"].astype(str))[["CHR", "BP"]]
    # gene column may or may not exist; tolerate missing
    bell_meta["gene"] = (lab.set_index(lab["SNP"].astype(str))
                          .get("gene", pd.Series("", index=bell_meta.index))
                          .reindex(bell_meta.index))
    print(f"  Bellenguez partners: {len(bell_rs)}  (chr distribution: "
          + ", ".join(f"chr{c}:{n}" for c, n in
                       bell["CHR"].value_counts().sort_index().items())
          + ")")

    # ── Per-query PLINK r² ───────────────────────────────────────────────
    tmp = Path(tempfile.mkdtemp(prefix="ld_"))
    rows = []
    per_query = {}
    for q, src in QUERIES.items():
        print(f"[plink] --ld-snp {q}  (source: {src})", flush=True)
        out_prefix = tmp / q
        cmd = [str(args.plink),
               "--bfile", str(args.bfile),
               "--ld-snp", q,
               "--r2",
               "--ld-window", "99999",
               "--ld-window-kb", str(args.window_kb),
               "--ld-window-r2", "0",        # integer 0 — print every pair
               "--out", str(out_prefix)]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  [FAIL] rc={r.returncode}: {r.stderr[-500:]}")
            per_query[q] = {"n_partners": 0, "error": r.stderr[-200:]}
            continue
        ld_file = out_prefix.with_suffix(".ld")
        if not ld_file.exists():
            print(f"  [skip] no .ld file for {q}")
            per_query[q] = {"n_partners": 0, "error": "no .ld"}
            continue
        df = pd.read_csv(ld_file, sep=r"\s+")
        # PLINK .ld columns: CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
        df["query_rsID"] = q
        df["query_source"] = src
        df["partner_rsID"] = df["SNP_B"].astype(str)
        # drop self (SNP_A == SNP_B)
        df = df[df["partner_rsID"] != q].copy()
        df["query_CHR"] = df["CHR_A"].astype(str)
        df["query_BP"] = df["BP_A"].astype(int)
        df["partner_CHR"] = df["CHR_B"].astype(str)
        df["partner_BP"] = df["BP_B"].astype(int)
        df["dist_bp"] = (df["partner_BP"] - df["query_BP"]).abs().astype(int)
        df["r2"] = df["R2"].astype(float)
        # is partner a Bellenguez lead?
        df["partner_source"] = df["partner_rsID"].map(
            lambda r: "Bellenguez2022" if r in bell_rs else "")
        df["partner_gene"] = df["partner_rsID"].map(
            lambda r: bell_meta.loc[r, "gene"] if r in bell_meta.index else "")
        df["query_in_apoe_cluster"] = df.apply(
            lambda r: _in_apoe(r["query_CHR"], r["query_BP"]), axis=1)
        df["partner_in_apoe_cluster"] = df.apply(
            lambda r: _in_apoe(r["partner_CHR"], r["partner_BP"]), axis=1)
        df["flag_r2_gt_080"] = df["r2"] > 0.8
        df["flag_r2_gt_090"] = df["r2"] > 0.9
        rows.append(df[["query_rsID", "query_source", "query_CHR",
                         "query_BP", "query_in_apoe_cluster",
                         "partner_rsID", "partner_source", "partner_gene",
                         "partner_CHR", "partner_BP",
                         "partner_in_apoe_cluster", "dist_bp", "r2",
                         "flag_r2_gt_080", "flag_r2_gt_090"]])

        # per-query summary
        bell_partners = df[df["partner_source"] == "Bellenguez2022"]
        top3 = (df.nlargest(3, "r2")
                  [["partner_rsID", "partner_source", "dist_bp", "r2"]]
                  .to_dict("records"))
        per_query[q] = {
            "n_partners": int(len(df)),
            "n_bellenguez_partners": int(len(bell_partners)),
            "max_r2_overall": float(df["r2"].max()) if len(df) else 0.0,
            "max_r2_bellenguez": (float(bell_partners["r2"].max())
                                   if len(bell_partners) else 0.0),
            "n_r2_gt_080": int(df["flag_r2_gt_080"].sum()),
            "n_r2_gt_080_bellenguez":
                int(bell_partners["flag_r2_gt_080"].sum()),
            "n_r2_gt_090_bellenguez":
                int(bell_partners["flag_r2_gt_090"].sum()),
            "query_in_apoe_cluster": bool(df["query_in_apoe_cluster"].any()),
            "n_partners_in_apoe_cluster":
                int(df["partner_in_apoe_cluster"].sum()),
            "top3_partners": top3,
        }

    if not rows:
        sys.exit("[ERROR] no LD output produced")

    all_pairs = pd.concat(rows, ignore_index=True)
    all_pairs.to_csv(args.out / "pairwise_r2.tsv", sep="\t", index=False)
    print(f"\n[out] {args.out/'pairwise_r2.tsv'}  ({len(all_pairs)} pairs)")

    # ── Human-readable summary ───────────────────────────────────────────
    lines = [
        "Pairwise LD r² — 3 re-extracted Wightman + 5 novel Desikan + 2 "
        "novel Kunkle query SNPs vs neighbouring SNPs in MAF>0.01 cohort "
        "(in-cohort LD, 616 ADNI samples ~European)",
        "Window: 1 Mb (user-set); --ld-window-r2 0 (every pair printed).",
        f"APOE cluster: chr{APOE_CHR}:{APOE_LO:,}–{APOE_HI:,}",
        f"Bellenguez partner set: {len(bell_rs)} SNPs",
        "",
    ]
    for q, src in QUERIES.items():
        s = per_query.get(q, {})
        if "error" in s:
            lines.append(f"## {q}  ({src})  [PLINK FAILED: {s['error']}]")
            continue
        lines.append(f"## {q}  ({src})")
        lines.append(f"  n_partners (within 1 Mb)       : {s['n_partners']}")
        lines.append(f"  n_Bellenguez partners          : "
                     f"{s['n_bellenguez_partners']}")
        lines.append(f"  max r² (any partner)           : "
                     f"{s['max_r2_overall']:.4f}")
        lines.append(f"  max r² to Bellenguez           : "
                     f"{s['max_r2_bellenguez']:.4f}")
        lines.append(f"  r² > 0.8 (any partner)         : {s['n_r2_gt_080']}")
        lines.append(f"  r² > 0.8 to Bellenguez (★ LD)  : "
                     f"{s['n_r2_gt_080_bellenguez']}")
        lines.append(f"  r² > 0.9 to Bellenguez (★★)   : "
                     f"{s['n_r2_gt_090_bellenguez']}")
        lines.append(f"  query_in_apoe_cluster          : "
                     f"{s['query_in_apoe_cluster']}")
        lines.append(f"  n partners in APOE cluster     : "
                     f"{s['n_partners_in_apoe_cluster']}")
        lines.append("  top 3 partners by r²:")
        for t in s["top3_partners"]:
            src_tag = (t["partner_source"] or "-")
            lines.append(f"    {t['partner_rsID']:<14} {src_tag:<18} "
                         f"dist={t['dist_bp']:>9} bp   r²={t['r2']:.4f}")
        lines.append("")

    n_high = int(all_pairs[(all_pairs["partner_source"] == "Bellenguez2022")
                            & all_pairs["flag_r2_gt_080"]]
                  ["query_rsID"].nunique())
    lines.append(f"---\nQuery SNPs with ≥1 Bellenguez partner at r²>0.8: "
                 f"{n_high} / {len(QUERIES)}")
    n_apoe_overlap = int(all_pairs["partner_in_apoe_cluster"].sum())
    lines.append(f"Partners in APOE cluster (any query): {n_apoe_overlap}")
    (args.out / "ld_summary.txt").write_text("\n".join(lines), "utf-8")
    print(f"[out] {args.out/'ld_summary.txt'}")
    shutil.rmtree(tmp, ignore_errors=True)
    print(f"\n=== verdict ===")
    print(f"{n_high} / {len(QUERIES)} query SNPs are LD-redundant (r²>0.8) "
          f"with ≥1 Bellenguez lead.")


if __name__ == "__main__":
    main()
