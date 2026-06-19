r"""
27b_assemble_set_embeddings.py   (LOCAL, env: snp)
==================================================
Smart-delta assembly: for each (model × new_set) in 5 × 7 = 35 pairs,
construct a per-set npz `<model>/<set>_snp_embeddings.npz` by
combining:
  - the existing legacy npz files (9W27B / 9W27B5D / APOE_cluster /
    APOE_causal) for rsIDs already embedded there;
  - the newly-extracted `delta_pool` npz for the 11 novel rsIDs.

The output schema is identical to script 27's:
  {rsids, snp_idx, ref_emb, alt_emb, model_id, H, pooling, strand_mode}
and the rsID order in the npz EXACTLY matches `<set>_snps.tsv`
(after manifest QC drop) so that `fm_diff_lib.load_diff` passes
its alignment assertion.

Numpy indexing only — no GPU, no model inference.

Usage:
  conda run -n snp python snp_pipeline/27b_assemble_set_embeddings.py
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
LEGACY_ROOT = BASE / "fm_embeddings_short_seq_1kb"  # bmfm_ref/ bmfm_snp/ ntv2/ caduceus_ph_d256/ caduceus_ps_d256/
DELTA_ROOT  = BASE / "diff_attn_drive_download"     # bmfm_ref/ bmfm_snp/ ntv2/ caduceus_ph/ caduceus_ps/
NEW_GWAS    = BASE / "GWAS_filtered_maf_resolved"

# Legacy set names ↔ filenames in <legacy_root>/<model_dir>/<set>_snp_embeddings.npz
LEGACY_SETS = ["9W27B5D", "9W27B", "APOE_cluster", "APOE_causal"]  # priority order
NEW_SETS    = ["5W26B", "5W26B14D", "5W26B14D2K",
                "sourceB", "sourceW", "sourceD", "sourceK"]

# Driver-side model key → directory name on disk for *legacy* embeddings.
# The bmfm/ntv2 dirs match by name; the caduceus dirs carry "_d256" suffix.
# The DELTA download uses the bare names (no suffix) — symmetric to the
# extraction.
MODELS = ["bmfm_ref", "bmfm_snp", "ntv2", "caduceus_ph", "caduceus_ps"]
LEGACY_DIR_ALIAS = {"caduceus_ph": "caduceus_ph_d256",
                     "caduceus_ps": "caduceus_ps_d256"}


def load_npz_index(npz_path: Path) -> tuple[dict[str, int], np.ndarray, np.ndarray, dict]:
    """Open an npz and return: {rsID -> row idx}, ref_emb, alt_emb,
    metadata dict (model_id, H, pooling, strand_mode)."""
    z = np.load(npz_path, allow_pickle=True)
    rsids = np.array([str(x) for x in z["rsids"]])
    idx = {r: i for i, r in enumerate(rsids)}
    meta = {k: z[k].item() if z[k].dtype.kind in "OU" and z[k].shape == ()
            else str(z[k]) for k in z.files
            if k not in {"rsids", "snp_idx", "ref_emb", "alt_emb"}}
    return idx, z["ref_emb"], z["alt_emb"], meta


def load_set_snps(set_name: str) -> pd.DataFrame:
    """Load the set's snps.tsv after applying manifest QC drops."""
    snps_path = NEW_GWAS / set_name / f"{set_name}_snps.tsv"
    man_path  = NEW_GWAS / set_name / "fm_sequences" / f"{set_name}_snp_manifest.tsv"
    snp = pd.read_csv(snps_path, sep="\t", dtype={"rsID": str})
    if man_path.exists():
        man = pd.read_csv(man_path, sep="\t", dtype={"rsID": str})
        # drop SNPs whose qc_verdict starts with 'drop_'
        if "qc_verdict" in man.columns:
            keep = set(man.loc[
                ~man["qc_verdict"].astype(str).str.startswith("drop_"),
                "rsID"].astype(str))
            snp = snp[snp["rsID"].astype(str).isin(keep)]
    return snp.reset_index(drop=True)


def assemble_set(model: str, set_name: str, legacy_indexes: list[tuple],
                  delta_idx: dict, delta_ref: np.ndarray, delta_alt: np.ndarray,
                  delta_meta: dict) -> tuple[np.ndarray, np.ndarray, list[str], list[int]]:
    """Build (ref_emb, alt_emb, rsids, snp_idx) for one set in this set's
    canonical order. Returns row counts. Per-rsID lookup: try delta first,
    then legacy in priority order."""
    snp = load_set_snps(set_name)
    rsids_in_order = snp["rsID"].astype(str).tolist()
    H = delta_ref.shape[1]
    ref_out = np.zeros((len(rsids_in_order), H), dtype=delta_ref.dtype)
    alt_out = np.zeros_like(ref_out)
    sources = []
    missing = []
    for i, r in enumerate(rsids_in_order):
        # delta first (newly-embedded)
        if r in delta_idx:
            j = delta_idx[r]
            ref_out[i] = delta_ref[j]; alt_out[i] = delta_alt[j]
            sources.append("delta"); continue
        # legacy npz priority order
        hit = None
        for ls_name, lidx, lref, lalt in legacy_indexes:
            if r in lidx:
                hit = (ls_name, lidx[r], lref, lalt); break
        if hit is None:
            missing.append(r); continue
        ls_name, j, lref, lalt = hit
        ref_out[i] = lref[j]; alt_out[i] = lalt[j]
        sources.append(ls_name)
    return ref_out, alt_out, rsids_in_order, sources, missing


def main() -> None:
    print(f"=== Assemble per-(model, new_set) npz for {len(MODELS)} models × {len(NEW_SETS)} sets ===\n")
    total_ok = 0; total_fail = 0; total_missing = 0
    for model in MODELS:
        legacy_dir_name = LEGACY_DIR_ALIAS.get(model, model)
        legacy_root_m   = LEGACY_ROOT / legacy_dir_name
        delta_root_m    = DELTA_ROOT / model
        delta_npz       = delta_root_m / "delta_pool_snp_embeddings.npz"

        # Load delta + legacy once per model
        delta_idx, delta_ref, delta_alt, delta_meta = load_npz_index(delta_npz)
        legacy_indexes = []
        for ls in LEGACY_SETS:
            lp = legacy_root_m / f"{ls}_snp_embeddings.npz"
            if not lp.exists():
                print(f"  [skip-legacy] {model} / {ls}: {lp.name} not found")
                continue
            lidx, lref, lalt, _ = load_npz_index(lp)
            legacy_indexes.append((ls, lidx, lref, lalt))
        H = delta_ref.shape[1]
        print(f"\n[{model}]  delta:{len(delta_idx)} rsIDs   legacy:"
              f"{[f'{ls}={len(lidx)}' for ls,lidx,*_ in legacy_indexes]}   H={H}")

        for set_name in NEW_SETS:
            try:
                ref_out, alt_out, rsids, sources, missing = assemble_set(
                    model, set_name, legacy_indexes,
                    delta_idx, delta_ref, delta_alt, delta_meta)
                if missing:
                    print(f"  [FAIL] {set_name}: {len(missing)} rsIDs missing: "
                          f"{missing[:5]}…")
                    total_fail += 1
                    total_missing += len(missing)
                    continue
                # source breakdown
                from collections import Counter
                src_ctr = Counter(sources)
                # QC: per-row REF≠ALT (most should differ on z; caduceus may not)
                d_norm = float(np.linalg.norm(alt_out - ref_out, axis=1).mean())

                # Write npz mirroring script-27 schema
                outdir = legacy_root_m   # WRITE BACK into the existing dir so
                                          #  script 28 finds it by --model name
                outpath = outdir / f"{set_name}_snp_embeddings.npz"
                np.savez(
                    outpath,
                    rsids=np.array(rsids, dtype=object),
                    snp_idx=np.arange(len(rsids), dtype=np.int64),
                    ref_emb=ref_out,
                    alt_emb=alt_out,
                    model_id=str(delta_meta.get("model_id", model)),
                    H=H,
                    pooling=str(delta_meta.get("pooling", "?")),
                    strand_mode=str(delta_meta.get("strand_mode", "fwd")),
                )
                src_summary = ", ".join(f"{k}={v}" for k, v in src_ctr.most_common())
                print(f"  [out] {set_name:<14} ({len(rsids):3d} SNPs, ‖diff‖={d_norm:.4f})"
                      f"  sources: {src_summary}")
                total_ok += 1
            except Exception as e:
                print(f"  [ERR] {set_name}: {type(e).__name__}: {e}")
                total_fail += 1

    print(f"\n=== done ===  ok={total_ok}  fail={total_fail}  total_missing_rsids={total_missing}")


if __name__ == "__main__":
    main()
