r"""
25f_extend_combos.py   (LOCAL, env: snp)
========================================
Extend the missingness-resolved combos at
`D:\…\GWAS_filtered_missingness_resolved_geno005\` to the full 17-combo
table: 4 source-only + 6 two-way unions + 4 three-way + 1 four-way + 2
LD-pruned variants. Also rename the 4 existing union folders to the
correct missingness-resolved counts (`5W26B*` → `6W28B*` per the
B > D > W > K priority dedup), rewrite `combos_summary.tsv`, and emit a
`duplicates_summary.txt` documenting the 6 multi-source rsIDs + the 2
LD-redundant pairs from 30e.

Reads only the per-source `<sourceX>/<sourceX>_snps.tsv` files (already
β-harmonised by 25e) and the 51-SNP union dosage; does NOT re-extract
from PLINK so the dosage of any rsID is bit-identical across combos.

Combo naming (W → B → D → K order, smallest priority first within each
multi-letter combo, matching the user's `5W26B14D2K` convention):

  source-only:   sourceB (28), sourceW (10), sourceD (16), sourceK (3)
  2-way (6):     6W28B, 28B15D, 28B3K, 10W16D, 10W3K, 16D2K
  3-way (4):     6W28B15D, 6W28B3K, 28B15D2K, 10W16D2K
  4-way (1):     6W28B15D2K
  D-pruned (2):  6W28B14D, 6W28B14D2K  (drop rs2597283; intra-D LD r²=0.91)

Usage:
  conda run -n snp python snp_pipeline/25f_extend_combos.py
"""
from __future__ import annotations

import argparse
import importlib.util
import shutil
from pathlib import Path

import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
OUT_ROOT = BASE / "GWAS_filtered_missingness_resolved_geno005"
PLINK = BASE / "plink.exe"

# Reuse 25d's write_combo + dedupe_union + helpers via importlib
SCRIPT_25D = Path(__file__).with_name("25d_export_maf_resolved_sets.py")
_s25d_spec = importlib.util.spec_from_file_location("_s25d", SCRIPT_25D)
_s25d = importlib.util.module_from_spec(_s25d_spec)
_s25d_spec.loader.exec_module(_s25d)

SCRIPT_25E = Path(__file__).with_name("25e_export_missingness_resolved_sets.py")
_s25e_spec = importlib.util.spec_from_file_location("_s25e", SCRIPT_25E)
_s25e = importlib.util.module_from_spec(_s25e_spec)
_s25e_spec.loader.exec_module(_s25e)

write_combo = _s25d.write_combo
dedupe_union = _s25d.dedupe_union
PRIORITY = _s25d.PRIORITY
plink_extract_dosage = _s25e.plink_extract_dosage
NEW_BED = _s25e.NEW_BED_STEM

# Source label → directory name (existing convention)
SOURCE_DIR = {"Bellenguez": "sourceB", "Wightman": "sourceW",
              "Desikan": "sourceD", "Kunkle": "sourceK"}

# Old → new union folder renames
RENAMES = {
    "5W26B": "6W28B",
    "5W26B14D": "6W28B15D",
    "5W26B13D": "6W28B14D",
    "5W26B14D2K": "6W28B15D2K",
}

INTRA_D_LD_PAIR = ("rs17265593", "rs2597283")   # r²=0.91; keep rs17265593


def _load_source_df(src_dir: str) -> pd.DataFrame:
    """Load a per-source snps.tsv into the format write_combo expects
    (renames effect_allele → A1, other_allele → A2)."""
    p = OUT_ROOT / src_dir / f"{src_dir}_snps.tsv"
    df = pd.read_csv(p, sep="\t", dtype=str)
    df["BP_GRCh38"] = pd.to_numeric(df["BP_GRCh38"], errors="coerce").astype("Int64")
    df["beta_A1"] = pd.to_numeric(df["beta_A1"], errors="coerce")
    return df.rename(columns={"effect_allele": "A1", "other_allele": "A2"})


def _rename_union(old_name: str, new_name: str) -> None:
    """Rename a union folder and the files inside that embed the combo name."""
    old_dir = OUT_ROOT / old_name
    new_dir = OUT_ROOT / new_name
    if not old_dir.exists():
        print(f"  [skip rename] {old_name} not found")
        return
    if new_dir.exists():
        print(f"  [skip rename] {new_name} already exists — removing old")
        shutil.rmtree(old_dir)
        return
    print(f"  [rename] {old_name}/ → {new_name}/")
    old_dir.rename(new_dir)
    # Rename embedded files
    for ext in ("_snps.tsv", "_patient_dosage.tsv", "_patient_dosage.bim"):
        old_file = new_dir / f"{old_name}{ext}"
        new_file = new_dir / f"{new_name}{ext}"
        if old_file.exists():
            old_file.rename(new_file)
    # Patch README to mention the new name (best-effort string replace)
    readme = new_dir / "README.txt"
    if readme.exists():
        txt = readme.read_text("utf-8")
        readme.write_text(txt.replace(old_name, new_name), "utf-8")


def _write_duplicates_summary(per_source: dict[str, pd.DataFrame]) -> None:
    """Compute multi-source rsIDs + emit the duplicates_summary.txt."""
    membership: dict[str, set[str]] = {}
    for src, df in per_source.items():
        for rs in df["rsID"]:
            membership.setdefault(rs, set()).add(src)

    multi = {rs: srcs for rs, srcs in membership.items() if len(srcs) > 1}
    print(f"  multi-source rsIDs: {len(multi)}")

    # Group by source-set
    pair_groups: dict[frozenset, list[str]] = {}
    for rs, srcs in multi.items():
        pair_groups.setdefault(frozenset(srcs), []).append(rs)

    def _winner(srcs: set[str]) -> str:
        return min(srcs, key=lambda s: PRIORITY[s])

    lines = [
        "Missingness-resolved (geno 0.05) PRS multi-source duplicates summary",
        "",
        "Priority hierarchy (winner kept when an rsID is in ≥2 sources):",
        "  B > D > W > K",
        "  i.e. Bellenguez (largest post-2022 GWAS) > Desikan PGS >",
        "       Wightman (no-UKB, independent replication) > Kunkle 2019.",
        "  The winner's β/OR is the one used in the dedup_union; the loser's",
        "  β is logged but not applied. Use the per-source folders",
        "  (sourceB/W/D/K) if you specifically want a source's own β even",
        "  for shared rsIDs.",
        "",
        f"Multi-source rsIDs ({len(multi)} total):",
        "",
    ]
    for srcs_set, rsids in sorted(pair_groups.items(),
                                    key=lambda kv: tuple(sorted(kv[0]))):
        srcs = set(srcs_set)
        winner = _winner(srcs)
        loser = sorted(s for s in srcs if s != winner)
        lines.append(f"  ## {' ∩ '.join(sorted(srcs))}  →  winner: {winner}")
        for rs in sorted(rsids):
            wb = per_source[winner].loc[per_source[winner]["rsID"] == rs,
                                          "beta_A1"]
            wb_str = (f"{float(wb.iloc[0]):+.4f}" if len(wb) else "—")
            loser_str_parts = []
            for ls in loser:
                lb = per_source[ls].loc[per_source[ls]["rsID"] == rs,
                                          "beta_A1"]
                lb_str = (f"{float(lb.iloc[0]):+.4f}" if len(lb) else "—")
                loser_str_parts.append(f"{ls}={lb_str}")
            lines.append(f"    {rs:14s}  "
                          f"{winner}_β={wb_str}   "
                          f"({', '.join(loser_str_parts)})")
        lines.append("")

    # LD-redundant pairs (from 30e LD analysis)
    lines.extend([
        "=" * 70,
        "LD-redundant pairs (distinct rsIDs in high LD; from 30e):",
        "",
        "  Intra-Desikan (chr7 BC043356 locus):",
        f"    rs17265593  ↔  rs2597283   r²=0.9059   dist=70,585 bp",
        "    → keep rs17265593 (labelled BC043356 lead in PGS000026);",
        "      drop rs2597283 in the LD-pruned variant `6W28B14D` (and",
        "      `6W28B14D2K`).",
        "",
        "  Cross-source W↔D (chr11 PICALM locus):",
        f"    rs561655 (W)  ↔  rs543293 (D)   r²=0.8414   dist=19,798 bp",
        "    → BOTH kept in default combos (distinct rsIDs at distinct",
        "      positions; downstream LD-aware PRS aggregation can decide).",
        "      User has not requested a `*_pruneWD` variant.",
        "",
        "=" * 70,
        "Combo naming convention (W → B → D → K order, matching the existing",
        "`5W26B14D2K` precedent). Each digit is the count of rsIDs CONTRIBUTED",
        "BY that source AFTER the priority dedup:",
        "",
        "  28B = 28 Bellenguez   (all of Bellenguez when B is in the combo)",
        "  6W  = 6 Wightman      (10 − 4 B∩W when both B and W are in the combo)",
        "  15D = 15 Desikan      (16 − 1 B∩D when B is in the combo)",
        "  16D = 16 Desikan      (full Desikan when no B in the combo)",
        "  14D = 14 Desikan      (15D − 1 intra-D LD prune rs2597283)",
        "  2K  = 2 Kunkle        (3 − 1 D∩K when D is in the combo)",
        "  3K  = 3 Kunkle        (full Kunkle when no D in the combo)",
        "  10W = 10 Wightman     (full Wightman when no B in the combo)",
    ])

    (OUT_ROOT / "duplicates_summary.txt").write_text("\n".join(lines), "utf-8")
    print(f"  [out] {OUT_ROOT/'duplicates_summary.txt'}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    args = ap.parse_args()

    # ── Load the 4 per-source DataFrames ──────────────────────────────────
    print("=== Stage 1 — load per-source SNP tables ===")
    bell = _load_source_df("sourceB")
    wight = _load_source_df("sourceW")
    desi = _load_source_df("sourceD")
    kunk = _load_source_df("sourceK")
    print(f"  Bellenguez={len(bell)}  Wightman={len(wight)}  "
          f"Desikan={len(desi)}  Kunkle={len(kunk)}")

    # ── Extract union dosage from the 51-SNP BED ─────────────────────────
    print("\n=== Stage 2 — extract union dosage from missingness-resolved BED ===")
    union_rs = sorted(set(bell["rsID"]) | set(wight["rsID"])
                       | set(desi["rsID"]) | set(kunk["rsID"]))
    dos = plink_extract_dosage(union_rs, NEW_BED, OUT_ROOT / "_dosage_union")
    print(f"  dosage shape: {dos.shape}")

    # ── Stage 3 — rename existing 4 union folders ────────────────────────
    print("\n=== Stage 3 — rename existing 4 union folders ===")
    for old, new in RENAMES.items():
        _rename_union(old, new)

    # ── Stage 4 — write all 17 combos (rewriting renamed + adding new) ───
    print("\n=== Stage 4 — write all 17 combos ===")
    desi_pruned = desi[desi["rsID"] != INTRA_D_LD_PAIR[1]].copy()

    combos_def = [
        # (name, [(src, df), ...], description)
        ("sourceB", [("Bellenguez", bell)],
         "Bellenguez 2022 alone (28 SNPs, missingness-resolved)."),
        ("sourceW", [("Wightman", wight)],
         "Wightman 2021 alone (10 SNPs; 2 published-PRS rescues: "
         "rs117618017 from GCST013196 direct β, rs115186657 from "
         "GCST013197 Z×SE_computed)."),
        ("sourceD", [("Desikan", desi)],
         "Desikan 2017 PGS000026 alone (16 SNPs; +rs4266886 as a "
         "published-PRS missingness-bypass exception)."),
        ("sourceK", [("Kunkle", kunk)],
         "Kunkle 2019 alone (3 SNPs)."),

        # 2-way unions
        ("6W28B", [("Bellenguez", bell), ("Wightman", wight)],
         "W ∪ B (34 SNPs). B wins B∩W overlaps "
         "(rs871269, rs1582763, rs7912495, rs117618017)."),
        ("28B15D", [("Bellenguez", bell), ("Desikan", desi)],
         "B ∪ D (43 SNPs). B wins B∩D rs6733839; no W or K."),
        ("28B3K", [("Bellenguez", bell), ("Kunkle", kunk)],
         "B ∪ K (31 SNPs). B∩K = ∅; no W or D."),
        ("10W16D", [("Desikan", desi), ("Wightman", wight)],
         "W ∪ D (26 SNPs). W∩D = ∅ (rs561655/rs543293 at PICALM are "
         "distinct rsIDs in r²=0.84 LD; both kept). No B or K."),
        ("10W3K", [("Wightman", wight), ("Kunkle", kunk)],
         "W ∪ K (13 SNPs). W∩K = ∅; no B or D."),
        ("16D2K", [("Desikan", desi), ("Kunkle", kunk)],
         "D ∪ K (18 SNPs). D wins D∩K rs7920721; no B or W."),

        # 3-way unions
        ("6W28B15D", [("Bellenguez", bell), ("Desikan", desi),
                       ("Wightman", wight)],
         "W ∪ B ∪ D (49 SNPs). No LD prune; no K."),
        ("6W28B3K", [("Bellenguez", bell), ("Wightman", wight),
                      ("Kunkle", kunk)],
         "W ∪ B ∪ K (37 SNPs). No D."),
        ("28B15D2K", [("Bellenguez", bell), ("Desikan", desi),
                        ("Kunkle", kunk)],
         "B ∪ D ∪ K (45 SNPs). D wins D∩K rs7920721; no W."),
        ("10W16D2K", [("Desikan", desi), ("Wightman", wight),
                        ("Kunkle", kunk)],
         "W ∪ D ∪ K (28 SNPs). D wins D∩K rs7920721; no B."),

        # 4-way + pruned
        ("6W28B15D2K", [("Bellenguez", bell), ("Desikan", desi),
                          ("Wightman", wight), ("Kunkle", kunk)],
         "W ∪ B ∪ D ∪ K (51 SNPs). All 4 sources; no LD prune."),
        ("6W28B14D", [("Bellenguez", bell), ("Desikan", desi_pruned),
                       ("Wightman", wight)],
         "W ∪ B ∪ D with intra-D LD prune (48 SNPs; drop rs2597283, "
         "keep rs17265593 as the BC043356 lead; intra-D r²=0.91)."),
        ("6W28B14D2K", [("Bellenguez", bell), ("Desikan", desi_pruned),
                          ("Wightman", wight), ("Kunkle", kunk)],
         "W ∪ B ∪ D ∪ K with intra-D LD prune (50 SNPs; drop rs2597283)."),
    ]

    summary = []
    for name, parts, desc in combos_def:
        if len(parts) == 1:
            # source-only: write_combo already handles single-source df
            combo_df = parts[0][1]
        else:
            combo_df = dedupe_union(parts)
        summary.append(write_combo(name, combo_df, dos, OUT_ROOT, desc))

    # ── Stage 5 — rewrite combos_summary.tsv ─────────────────────────────
    sdf = pd.DataFrame(summary)
    sdf.to_csv(OUT_ROOT / "combos_summary.tsv", sep="\t", index=False)
    print(f"\n[out] {OUT_ROOT/'combos_summary.tsv'}  ({len(sdf)} rows)")

    # ── Stage 6 — duplicates_summary.txt ─────────────────────────────────
    print("\n=== Stage 6 — duplicates_summary.txt ===")
    _write_duplicates_summary({"Bellenguez": bell, "Wightman": wight,
                                 "Desikan": desi, "Kunkle": kunk})

    print("\n=== Summary ===")
    for s in summary:
        print(f"  {s['combo']:14s}  n={s['n_snps']:>3}  "
              f"B={s['n_bellenguez']}, W={s['n_wightman']}, "
              f"D={s['n_desikan']}, K={s['n_kunkle']}")


if __name__ == "__main__":
    main()
