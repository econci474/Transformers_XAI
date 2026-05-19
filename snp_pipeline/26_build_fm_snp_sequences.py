r"""
26_build_fm_snp_sequences.py   (LOCAL, env: snp)
================================================
Regenerate, per filtered SNP set, the **REF** and **ALT** 1 kb SNP-centred
sequences (+ reverse-complement) directly from the GRCh38 FASTA, driven by
the script-25 `<SET>_snps.tsv`. The legacy `bmfm_gwas_signed_regression_*`
files cover only the original 430 SNPs (and carry z-scores, EA/OA labels) —
the filtered sets contain re-extracted Wightman + novel Desikan SNPs absent
from them, and allele orientation must match the dosage/β (effect_allele =
ALT = A1 = dosage-counted), so we regenerate rather than edit.

Reuses verbatim: resolve_chrom + windowing math + _make_seq from
08c_prepare_bmfm_gwas_regression_inputs.py (L164-216); reverse_complement
from 08d_augment_bmfm_gwas_regression.py (L82-86). SNP ordering / allele
validation from fm_diff_lib (the single contract shared with 27/28).

Input  (per --set <S>, e.g. 9W27B / 9W27B5D / APOE_cluster):
  D:\ADNI_SNP_Omni2.5M_20140220\GWAS_filtered\<S>\<S>_snps.tsv
  ...\<S>\<S>_patient_dosage.tsv  +  ...\<S>\<S>_patient_dosage.bim
  D:\ADNI_SNP_Omni2.5M_20140220\Homo_sapiens.GRCh38.dna.primary_assembly.fa

Output → ...\GWAS_filtered\<S>\fm_sequences\ :
  <S>_snp_sequences.tsv   one row per (rsID, REF|ALT, fwd|rc)  + provenance
  <S>_snp_manifest.tsv    one row per SNP — qc_verdict (the single drop authority)
  <S>_seqbuild_qc.json    machine-readable QC report
  <S>_seqbuild_summary.txt

Usage:
  conda run -n snp python snp_pipeline/26_build_fm_snp_sequences.py --set 9W27B --check
  conda run -n snp python snp_pipeline/26_build_fm_snp_sequences.py --set all
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util as _il
import io
import json
import sys
from pathlib import Path

import pandas as pd

if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")

_LIB = Path(__file__).parent / "fm_diff_lib.py"
_spec = _il.spec_from_file_location("fm_diff_lib", _LIB)
fl = _il.module_from_spec(_spec)
_spec.loader.exec_module(fl)

BASE = fl.BASE
GWAS_FILTERED = fl.GWAS_FILTERED
FASTA = BASE / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# reverse_complement — verbatim 08d L82-86
_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", required=True,
                    help="9W27B | 9W27B5D | APOE_cluster | … | all")
    ap.add_argument("--gwas-root", type=Path, default=GWAS_FILTERED)
    ap.add_argument("--fasta", type=Path, default=FASTA)
    ap.add_argument("--flank", type=int, default=500,
                    help="bp each side → 1+2*flank window (default 1001)")
    ap.add_argument("--emit-rc", dest="emit_rc", action="store_true",
                    default=True)
    ap.add_argument("--no-emit-rc", dest="emit_rc", action="store_false")
    ap.add_argument("--check", action="store_true",
                    help="dry-run: QC + counts, write nothing")
    ap.add_argument("--seed", type=int, default=42, help="recorded only")
    args = ap.parse_args()

    sets = (sorted([p.name for p in args.gwas_root.iterdir()
                    if (p / f"{p.name}_snps.tsv").exists()])
            if args.set == "all" else [args.set])

    from pyfaidx import Fasta, FastaIndexingError
    print(f"[fasta] opening {args.fasta.name} …", flush=True)
    try:
        fa = Fasta(str(args.fasta), rebuild=False)
    except FastaIndexingError:
        fa = Fasta(str(args.fasta), rebuild=True)
    chrom_names = set(fa.keys())

    def resolve_chrom(chrom_val: str):                # verbatim 08c L164-169
        for cand in [str(chrom_val), f"chr{chrom_val}"]:
            if cand in chrom_names:
                return cand
        return None

    FLANK = args.flank
    seq_hashes_first: dict[str, str] = {}             # cross-set determinism

    for sname in sets:
        sdir = args.gwas_root / sname
        snps_tsv = sdir / f"{sname}_snps.tsv"
        dos_tsv = sdir / f"{sname}_patient_dosage.tsv"
        bim_tsv = sdir / f"{sname}_patient_dosage.bim"
        if not snps_tsv.exists():
            print(f"[skip] {sname}: no {snps_tsv.name}")
            continue
        print(f"\n=== {sname} ===", flush=True)

        snp = fl.canonical_snp_order(fl.load_snp_table(snps_tsv))
        dos_cols = (set(pd.read_csv(dos_tsv, sep="\t", nrows=0).columns)
                    - {"PTID"}) if dos_tsv.exists() else set()
        bim_rs = (set(pd.read_csv(bim_tsv, sep="\t")["rsID"].astype(str))
                  if bim_tsv.exists() else set())

        qc = {"set": sname, "n_snps_in": int(len(snp)), "flank": FLANK,
              "window_bp": 2 * FLANK + 1, "seed": args.seed,
              "drop_invalid_allele": [], "drop_nan_beta": [],
              "strand_flipped": [], "warn_palindromic": [],
              "warn_ref_unexplained": [], "warn_edge": [],
              "not_in_dosage": [], "not_in_bim": [], "rc_involution_ok": None,
              "centre_selfcheck_fail": []}
        seq_rows, man_rows = [], []

        for r in snp.itertuples(index=False):
            rsid = str(r.rsID)
            ea = str(r.effect_allele).strip().upper()   # ALT / A1 (β-referent)
            oa = str(r.other_allele).strip().upper()     # REF / A2
            verdict, reason = "ok", "ok"
            in_dos = rsid in dos_cols
            in_bim = rsid in bim_rs
            if not in_dos:
                qc["not_in_dosage"].append(rsid)
            if not in_bim:
                qc["not_in_bim"].append(rsid)

            ok, why = fl.validate_alleles(ea, oa)
            if not ok:
                verdict, reason = "drop_invalid_allele", why
                qc["drop_invalid_allele"].append(rsid)
            elif not pd.notna(r.beta_A1):
                verdict, reason = "drop_nan_beta", "nan_beta"
                qc["drop_nan_beta"].append(rsid)

            chrom_fa = resolve_chrom(str(r.CHR))
            if verdict == "ok" and chrom_fa is None:
                verdict, reason = "drop_invalid_allele", "chrom_not_in_fasta"
                qc["drop_invalid_allele"].append(rsid)

            man = {"snp_idx": int(r.snp_idx), "rsID": rsid, "CHR": str(r.CHR),
                   "BP_GRCh38": int(r.BP_GRCh38), "effect_allele": ea,
                   "other_allele": oa,
                   "beta_A1": (float(r.beta_A1) if pd.notna(r.beta_A1)
                               else float("nan")),
                   "qc_verdict": verdict, "qc_reason": reason,
                   "in_dosage": in_dos, "in_bim": in_bim,
                   "fasta_ref_base": "", "ref_matches_fasta": ""}
            if verdict != "ok":
                man_rows.append(man)
                print(f"  [drop] {rsid}: {reason}")
                continue

            pos1 = int(r.BP_GRCh38)                       # 1-based GRCh38
            pos0 = pos1 - 1
            clen = len(fa[chrom_fa])
            start = max(0, pos0 - FLANK)
            end = min(clen, pos0 + FLANK + 1)
            base_seq = str(fa[chrom_fa][start:end]).upper()
            snp_in_win = pos0 - start
            wlen = len(base_seq)
            fasta_ref = (base_seq[snp_in_win]
                         if 0 <= snp_in_win < wlen else "")
            # Strand-normalise the alleles to the FASTA(+) strand before
            # substitution (the script-25 export carries dosage-bim A1/A2,
            # which can be minus-strand vs GRCh38). Inserting a minus-strand
            # base into plus-strand flanks would be a chimeric sequence.
            # Mirrors the strand-aware harmoniser in script 24 (palindromic
            # A/T & C/G guard). β/dosage orientation is unchanged — only the
            # physical base representing the effect allele is corrected.
            _C = {"A": "T", "T": "A", "C": "G", "G": "C"}
            palindromic = ({ea, oa} == {"A", "T"} or {ea, oa} == {"C", "G"})
            ea_u, oa_u, strand = ea, oa, "plus"
            if fasta_ref in (ea, oa):
                strand = "plus"
            elif palindromic:
                strand = "palindromic"
                qc["warn_palindromic"].append(rsid)
                if man["qc_verdict"] == "ok":
                    man["qc_verdict"] = "warn_palindromic"
            elif fasta_ref in (_C.get(ea), _C.get(oa)):
                ea_u, oa_u, strand = _C[ea], _C[oa], "flip"
                qc["strand_flipped"].append(rsid)
            else:
                strand = "unexplained"
                qc["warn_ref_unexplained"].append(
                    {"rsID": rsid, "CHR": str(r.CHR), "BP": pos1,
                     "fasta": fasta_ref, "effect": ea, "other": oa})
                if man["qc_verdict"] == "ok":
                    man["qc_verdict"] = "warn_ref_unexplained"
            ref_match = (fasta_ref == oa_u)
            man["fasta_ref_base"] = fasta_ref
            man["ref_matches_fasta"] = bool(ref_match)
            man["strand"] = strand
            man["allele_alt_fasta"] = ea_u   # base actually substituted
            man["allele_ref_fasta"] = oa_u
            if wlen < 2 * FLANK + 1:
                qc["warn_edge"].append(rsid)
                if man["qc_verdict"] == "ok":
                    man["qc_verdict"] = "warn_edge"

            def _make_seq(allele: str) -> str:            # verbatim 08c L209
                s = list(base_seq)
                if 0 <= snp_in_win < len(s):
                    s[snp_in_win] = allele.upper()
                return "".join(s)

            ref_seq, alt_seq = _make_seq(oa_u), _make_seq(ea_u)
            # post-substitution centre self-check (asserts the index math)
            if not (ref_seq[snp_in_win] == oa_u
                    and alt_seq[snp_in_win] == ea_u):
                qc["centre_selfcheck_fail"].append(rsid)

            variants = [("REF", "fwd", oa_u, ref_seq),
                        ("ALT", "fwd", ea_u, alt_seq)]
            if args.emit_rc:
                variants += [("REF", "rc", oa_u, reverse_complement(ref_seq)),
                             ("ALT", "rc", ea_u, reverse_complement(alt_seq))]
            for atype, strand, allele, s in variants:
                seq_rows.append({
                    "seq_id": f"{rsid}__{atype}__{strand}",
                    "rsID": rsid, "snp_idx": int(r.snp_idx),
                    "allele_type": atype, "strand": strand, "allele": allele,
                    "CHR": str(r.CHR), "chrom_fa": chrom_fa,
                    "BP_GRCh38": pos1,
                    "pos_in_window": (snp_in_win if strand == "fwd"
                                      else wlen - 1 - snp_in_win),
                    "window_start_1based": start + 1, "window_len": wlen,
                    "fasta_ref_base": fasta_ref,
                    "ref_matches_fasta": bool(ref_match),
                    "allele_strand": strand,
                    "effect_allele_src": ea, "other_allele_src": oa,
                    "gene": r.gene, "lead_rsID": r.lead_rsID,
                    "source": r.source, "origin": r.origin,
                    "beta_A1": float(r.beta_A1), "sequence": s})
            man_rows.append(man)

        # RC involution check on one OK SNP
        ok_seqs = [x for x in seq_rows if x["strand"] == "fwd"]
        if ok_seqs:
            s0 = ok_seqs[0]["sequence"]
            qc["rc_involution_ok"] = bool(
                reverse_complement(reverse_complement(s0)) == s0)
        # cross-set determinism: hash the REF-fwd sequence per rsID
        for x in seq_rows:
            if x["allele_type"] == "REF" and x["strand"] == "fwd":
                h = hashlib.sha256(x["sequence"].encode()).hexdigest()[:16]
                key = x["rsID"]
                if key in seq_hashes_first and seq_hashes_first[key] != h:
                    qc.setdefault("cross_set_seq_mismatch", []).append(key)
                seq_hashes_first.setdefault(key, h)

        n_ok = sum(1 for m in man_rows
                   if m["qc_verdict"] in ("ok", "warn_ref_mismatch",
                                          "warn_edge"))
        qc["n_snps_usable"] = n_ok
        qc["n_seq_rows"] = len(seq_rows)
        summary = (
            f"{sname}: {len(snp)} SNPs in → {n_ok} usable "
            f"({len(qc['drop_invalid_allele'])} invalid-allele, "
            f"{len(qc['drop_nan_beta'])} nan-beta drop); "
            f"{len(qc['strand_flipped'])} strand-flipped(fixed), "
            f"{len(qc['warn_palindromic'])} palindromic, "
            f"{len(qc['warn_ref_unexplained'])} ref-unexplained, "
            f"{len(qc['warn_edge'])} edge, "
            f"{len(qc['not_in_dosage'])} not-in-dosage, "
            f"{len(qc['not_in_bim'])} not-in-bim; "
            f"{len(seq_rows)} seq rows "
            f"(rc_involution={qc['rc_involution_ok']}, "
            f"centre_fail={len(qc['centre_selfcheck_fail'])}, "
            f"xset_mismatch={len(qc.get('cross_set_seq_mismatch', []))})")
        print("  " + summary)
        if qc["centre_selfcheck_fail"] or qc.get("cross_set_seq_mismatch"):
            print("  [ERROR] centre/cross-set check failed — investigate")

        if args.check:
            print("  [check] no files written")
            continue
        outdir = sdir / "fm_sequences"
        outdir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(seq_rows).to_csv(
            outdir / f"{sname}_snp_sequences.tsv", sep="\t", index=False)
        pd.DataFrame(man_rows).to_csv(
            outdir / f"{sname}_snp_manifest.tsv", sep="\t", index=False)
        (outdir / f"{sname}_seqbuild_qc.json").write_text(
            json.dumps(qc, indent=2), encoding="utf-8")
        (outdir / f"{sname}_seqbuild_summary.txt").write_text(
            summary + "\n", encoding="utf-8")
        print(f"  [out] {outdir}")


if __name__ == "__main__":
    main()
