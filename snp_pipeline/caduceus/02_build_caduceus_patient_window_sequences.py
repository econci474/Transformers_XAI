"""
02_build_caduceus_patient_window_sequences.py
==============================================
Build per-patient × per-window DNA sequences sized for the Caduceus 131k
context window (131,072 bp), rather than the BMFM 8 kb window.

This is a direct adaptation of snp_pipeline/13_build_patient_window_sequences.py
with MAX_LEN raised from 8,192 bp to 131,072 bp.  With single-nucleotide
tokenisation (1 bp = 1 token), the full 131k context is available.

Effect:
  - Fewer, much wider windows per chromosome (many SNPs grouped together).
  - Far more flanking genomic context around each SNP cluster.
  - Typically yields ~5–10 windows total (vs ~25 at 8 kb).

Input
-----
  Same as 13_build_patient_window_sequences.py:
    patient_genotypes/patient_dosage.tsv   (from 12_extract_patient_genotypes.py)
    patient_genotypes/patient_dosage.bim   (allele table; A1 = patient alt)
    external_gwas_labels.tsv               (BP positions for the SNPs)
    Homo_sapiens.GRCh38.dna.primary_assembly.fa

Output
------
  D:/ADNI_SNP_Omni2.5M_20140220/caduceus_inputs/patient_window_sequences/
    windows.tsv
    sequences/all_patients.csv

Usage
-----
  conda activate snp
  python snp_pipeline/caduceus/02_build_caduceus_patient_window_sequences.py
  python snp_pipeline/caduceus/02_build_caduceus_patient_window_sequences.py --check
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR = Path(r"D:/ADNI_SNP_Omni2.5M_20140220")
DOSAGE_TSV = BASE_DIR / "bmfm_inputs" / "patient_genotypes" / "patient_dosage.tsv"
DOSAGE_BIM = BASE_DIR / "bmfm_inputs" / "patient_genotypes" / "patient_dosage.bim"
LABELS_TSV = BASE_DIR / "bmfm_inputs" / "external_gwas_labels.tsv"
FASTA = BASE_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Caduceus-specific output directory
OUT_DIR = BASE_DIR / "caduceus_inputs" / "patient_window_sequences"

# Caduceus 131k context: 131,072 bp per window (1 bp = 1 token)
MAX_LEN = 131_072
MIN_FLANK = 500    # more generous flanking for the wider context


def pack_windows(snp_rows: list, max_inner_span: int) -> list[list]:
    """Greedy 1-D packing of SNP rows (sorted by BP) into ≤ max_inner_span groups.

    Identical to the algorithm in 13_build_patient_window_sequences.py.
    """
    windows = []
    cur = [snp_rows[0]]
    for snp in snp_rows[1:]:
        if snp.BP - cur[0].BP >= max_inner_span:
            windows.append(cur)
            cur = [snp]
        else:
            cur.append(snp)
    windows.append(cur)
    return windows


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dosage", type=Path, default=DOSAGE_TSV)
    p.add_argument("--bim", type=Path, default=DOSAGE_BIM)
    p.add_argument("--labels", type=Path, default=LABELS_TSV)
    p.add_argument("--fasta", type=Path, default=FASTA)
    p.add_argument("--outdir", type=Path, default=OUT_DIR)
    p.add_argument("--max-len", type=int, default=MAX_LEN, dest="max_len",
                    help="Max window length in bp (default: 131,072 for Caduceus)")
    p.add_argument("--min-flank", type=int, default=MIN_FLANK, dest="min_flank")
    p.add_argument("--check", action="store_true",
                   help="Dry-run: print window stats only.")
    args = p.parse_args()

    for f, name in [(args.dosage, "dosage TSV"), (args.bim, "BIM"),
                    (args.labels, "labels TSV"), (args.fasta, "FASTA")]:
        if not f.exists():
            sys.exit(f"[ERROR] {name} not found: {f}")

    if 2 * args.min_flank >= args.max_len:
        sys.exit(f"[ERROR] min_flank ({args.min_flank}) * 2 >= max_len ({args.max_len})")

    inner_span_limit = args.max_len - 2 * args.min_flank

    # ── Load dosage matrix ────────────────────────────────────────────────────
    print(f"Loading dosage matrix: {args.dosage}")
    dosage = pd.read_csv(args.dosage, sep="\t")
    dosage["PTID"] = dosage["PTID"].astype(str)
    snp_cols = [c for c in dosage.columns if c != "PTID"]
    print(f"  {len(dosage)} patients × {len(snp_cols)} SNPs")

    # ── Load BIM (A1 = patient alt; A2 = REF) ─────────────────────────────────
    bim = pd.read_csv(args.bim, sep=r"\s+", engine="python", dtype=str)
    bim["CHR"] = bim["CHR"].astype(str)
    bim["BP"] = bim["BP"].astype(int)
    bim_lookup = {row.rsID: row for row in bim.itertuples(index=False)}
    print(f"  BIM: {len(bim)} rows")

    # ── Load labels TSV ───────────────────────────────────────────────────────
    labels = pd.read_csv(args.labels, sep="\t", low_memory=False)
    sig = labels.loc[labels["label"] == 1, ["SNP", "CHR", "BP"]].copy()
    sig["CHR"] = sig["CHR"].astype(str)
    sig["BP"] = sig["BP"].astype(int)
    sig = sig.sort_values(["CHR", "BP"]).reset_index(drop=True)
    sig = sig[sig["SNP"].isin(bim_lookup)]
    print(f"  Label=1 SNPs with BIM match: {len(sig)}")

    # ── Open FASTA ────────────────────────────────────────────────────────────
    try:
        from pyfaidx import Fasta
    except ImportError:
        sys.exit("[ERROR] pyfaidx not installed:  pip install pyfaidx")
    fa = Fasta(str(args.fasta), rebuild=False)
    chrom_names = set(fa.keys())

    def resolve_chrom(c: str) -> str | None:
        for candidate in (c, f"chr{c}"):
            if candidate in chrom_names:
                return candidate
        return None

    # ── Pack windows ──────────────────────────────────────────────────────────
    print(f"\nPacking windows (max_len={args.max_len:,} bp, "
          f"min_flank={args.min_flank:,} bp, "
          f"inner_span_limit={inner_span_limit:,} bp)")
    windows_meta = []
    win_refseqs = {}
    win_snp_info = {}
    win_id = 0
    skipped_chroms = []

    for chrom, grp in sig.groupby("CHR"):
        chrom_fa = resolve_chrom(chrom)
        if chrom_fa is None:
            skipped_chroms.append(chrom)
            continue
        snp_rows = list(grp.itertuples(index=False))
        windows = pack_windows(snp_rows, inner_span_limit)

        chrom_len = len(fa[chrom_fa])
        for wi, win in enumerate(windows):
            first_bp = win[0].BP
            last_bp = win[-1].BP
            inner = last_bp - first_bp
            extra = args.max_len - inner - 2 * args.min_flank
            left_flank = args.min_flank + extra // 2
            right_flank = args.min_flank + (extra - extra // 2)

            win_start0 = max(0, first_bp - 1 - left_flank)
            win_end0 = min(chrom_len, last_bp + right_flank)
            ref_seq = str(fa[chrom_fa][win_start0:win_end0]).upper()

            snp_info = []
            for row in win:
                bim_row = bim_lookup[row.SNP]
                offset = row.BP - 1 - win_start0
                snp_info.append({
                    "rsID": row.SNP,
                    "BP": row.BP,
                    "offset": offset,
                    "A1": bim_row.A1.upper(),
                    "A2": bim_row.A2.upper(),
                })
            window_id = f"chr{chrom}_w{wi:03d}"
            win_id += 1

            windows_meta.append({
                "window_id": window_id,
                "chrom": chrom,
                "start_1based": win_start0 + 1,
                "end_1based": win_end0,
                "length_bp": len(ref_seq),
                "n_snps": len(snp_info),
                "snp_ids": "|".join(s["rsID"] for s in snp_info),
                "snp_positions": "|".join(str(s["BP"]) for s in snp_info),
            })
            win_refseqs[window_id] = (chrom_fa, win_start0, ref_seq)
            win_snp_info[window_id] = snp_info

    if skipped_chroms:
        print(f"  [warn] skipped chroms not in FASTA: {sorted(set(skipped_chroms))}")

    print(f"  Total windows: {win_id}")
    print(f"  Windows per chromosome:")
    win_by_chrom: dict[str, int] = defaultdict(int)
    for w in windows_meta:
        win_by_chrom[w["chrom"]] += 1
    for c in sorted(win_by_chrom, key=lambda x: int(x) if x.isdigit() else 99):
        n_snps_in_chrom = sum(w["n_snps"] for w in windows_meta if w["chrom"] == c)
        print(f"    chr{c:>3}: {win_by_chrom[c]:>2} windows  "
              f"({n_snps_in_chrom:>3} SNPs)")

    # Summary: compare to 8kb windows
    total_snps = sum(w["n_snps"] for w in windows_meta)
    avg_len = sum(w["length_bp"] for w in windows_meta) / max(len(windows_meta), 1)
    avg_snps = total_snps / max(len(windows_meta), 1)
    print(f"\n  Summary vs 8 kb windows:")
    print(f"    Total windows:     {win_id} (was ~25 at 8 kb)")
    print(f"    Avg window length: {avg_len:,.0f} bp (was ~8,192 bp)")
    print(f"    Avg SNPs/window:   {avg_snps:.1f}")
    print(f"    Total SNPs:        {total_snps}")

    if args.check:
        print(f"\n[check] would emit {len(dosage)} patients × {win_id} windows = "
              f"{len(dosage) * win_id} sequences")
        return

    args.outdir.mkdir(parents=True, exist_ok=True)
    windows_df = pd.DataFrame(windows_meta)
    windows_df.to_csv(args.outdir / "windows.tsv", sep="\t", index=False)
    print(f"\nWrote {args.outdir / 'windows.tsv'}")

    # ── Build per-patient sequences ───────────────────────────────────────────
    seqs_dir = args.outdir / "sequences"
    seqs_dir.mkdir(parents=True, exist_ok=True)
    out_path = seqs_dir / "all_patients.csv"

    print(f"\nBuilding sequences  →  {out_path}")
    rows_out = []
    for prow in dosage.itertuples(index=False):
        ptid = prow.PTID
        dosage_for_pt = {col: getattr(prow, col) for col in snp_cols}
        for window_id, (_, _, ref_seq) in win_refseqs.items():
            seq = list(ref_seq)
            n_carried = 0
            for s in win_snp_info[window_id]:
                d = dosage_for_pt.get(s["rsID"])
                if pd.notna(d) and int(d) >= 1:
                    if 0 <= s["offset"] < len(seq):
                        seq[s["offset"]] = s["A1"]
                        n_carried += 1
            rows_out.append({
                "PTID": ptid,
                "window_id": window_id,
                "sequence": "".join(seq),
                "n_carried": n_carried,
            })

    out_df = pd.DataFrame(rows_out)
    out_df.to_csv(out_path, index=False)
    print(f"  Rows: {len(out_df)}  ({out_df['n_carried'].sum()} A1 substitutions total)")
    print(f"  Avg carried per (patient, window): {out_df['n_carried'].mean():.2f}")
    print(f"  Avg sequence length: {out_df['sequence'].str.len().mean():,.0f} bp")


if __name__ == "__main__":
    main()
