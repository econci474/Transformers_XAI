r"""
28_extract_alphagenome_scores.py   (COLAB or local snp env)
==========================================================
Zero-shot variant effect scoring via DeepMind AlphaGenome API for the
128-SNP recover_all_pool. Emits SIX strict-filter feature sets, slotting
into the existing per-model npz convention so the downstream
diff-attention loader treats them as six more models:

  alphagenome_eqtl             RNA-seq variant effect (GeneMaskLFCScorer)
                                 → max-abs LFC across cis genes per
                                   allow-listed brain tissue.
  alphagenome_accessibility    max_signed(DNase, ATAC) per allow-listed
                                 biosample (CenterMaskScorer).
  alphagenome_TSS              max_signed(CAGE, PRO-cap) per allow-listed
                                 biosample.
  alphagenome_splicing_merged  per-variant scalar:
                                 max|splice_sites| + max|splice_site_usage|
                                  + max|splice_junctions| / 5
                                 (literal Python precedence per user;
                                  raw per-component max-abs also saved).
  alphagenome_chip_histone     per (biosample × histone_mark) CHIP_HISTONE
                                 (CenterMaskScorer). 11 tracks (6 DLPFC +
                                 5 Layer-of-hippocampus, no H3K27me3 in
                                 Layer-of-hippo).
  alphagenome_chip_tf          per (biosample × TF) CHIP_TF
                                 (CenterMaskScorer). 1 track (DLPFC CTCF).

Strict allow-list filter (user 2026-06-07, replaces the 2026-05-26
brain-filter that had three independent bugs — wrong UBERON codes,
wrong GTEx string format, and a too-loose substring rule that pulled
in spinal-cord / cerebellar / motor-neuron / generic-astrocyte tracks):

  Region-anchored bulk tissue (adult or unspecified life stage only):
  - UBERON:0009834 dorsolateral prefrontal cortex (DLPFC, BA9)
  - UBERON:0001954 hippocampus (Ammon's horn)
  - UBERON:0002305 layer of hippocampus (adds hippocampal CHIP_HISTONE)

  Region-anchored cell-sorted (any life stage):
  - CL:0002604 astrocyte of the hippocampus (only cell-sorted
    region-anchored track in all of AlphaGenome for AD-relevant regions;
    1 DNase-seq track from ENCODE).

  Note: microglia, generic astrocyte/neuron, parietal-lobe-embryonic,
  motor neuron, neuronal stem cell, OPC, cerebellar/spinal-cord cells
  are NOT region-anchored to AD regions and are excluded by design.

Per-modality npz keys (mirrors 27's format where useful):
  rsids, snp_idx, effect_scores, track_biosamples, track_assays,
  model_tag, model_version, scorer, aggregation, interval_bp,
  ref_emb (zeros), alt_emb (= effect_scores).
For accessibility / TSS: extra source_assay (which assay won max_signed).
For splicing_merged: per-component max-abs sidecar arrays.

API key resolved in order: $ALPHA_GENOME_API_KEY → --api-key-file path.

Usage (Colab cell):
  !python snp_pipeline/28_extract_alphagenome_scores.py \
      --snps /content/drive/MyDrive/ADNI_SNP/diff_attn_drive_upload/recover_all_pool/recover_all_pool_snps.tsv \
      --set recover_all_pool \
      --outdir /content/drive/MyDrive/ADNI_SNP/fm_embeddings_short_seq_1kb
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ─────────────────────────────────────────────────────────────────────────
# Strict allow-list filter (user 2026-06-07).
# Replaces the buggy 2026-05-26 OR-filter — see module docstring.
# Matches AlphaGenome's `tidy_scores` columns: ontology_curie (UBERON/CL),
# biosample_type ("tissue" / "primary_cell" / "in_vitro_differentiated_cells"
# / "organoid"), biosample_life_stage ("adult" / "embryonic" / "" / …).
# ─────────────────────────────────────────────────────────────────────────
REGION_BULK_ALLOW = {
    # ontology_curie : human-readable label
    "UBERON:0009834": "DLPFC",                 # dorsolateral prefrontal cortex (BA9)
    "UBERON:0001954": "Hippocampus",           # Ammon's horn
    "UBERON:0002305": "Layer of hippocampus",  # 5 adult ENCODE CHIP_HISTONE tracks
}
CELLTYPE_ALLOW = {
    "CL:0002604": "astrocyte of the hippocampus",  # 1 ENCODE DNase track
}

SPLICE_SCORER_KEYS = ("SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS")


def _resolve_api_key(api_key_file: str | None) -> str:
    key = os.environ.get("ALPHA_GENOME_API_KEY", "").strip()
    if key:
        return key
    if api_key_file and Path(api_key_file).exists():
        key = Path(api_key_file).read_text(encoding="utf-8").strip()
        if key:
            return key
    sys.exit("[ERROR] no API key: set $ALPHA_GENOME_API_KEY or pass --api-key-file")


def _is_brain_track(row) -> bool:
    """Strict allow-list — True iff the track is one of:
      (a) bulk tissue from an allow-listed UBERON region, adult or no
          life-stage tag (excludes embryonic + organoid/cell tracks
          assigned to the same UBERON code), OR
      (b) one of the explicit allow-listed cell types (any life stage —
          currently only CL:0002604 hippocampal astrocyte qualifies).
    """
    onto = str(row.get("ontology_curie", "") or "")
    if onto in REGION_BULK_ALLOW:
        biosample_type = str(row.get("biosample_type", "") or "")
        life_stage = str(row.get("biosample_life_stage", "") or "")
        return biosample_type == "tissue" and life_stage in ("adult", "")
    return onto in CELLTYPE_ALLOW


def _max_signed(values: np.ndarray, axis: int = -1) -> np.ndarray:
    """Take the value with max abs magnitude along `axis`, preserving sign.
    `values` may have NaN — treated as missing (skipped)."""
    abs_v = np.abs(np.where(np.isnan(values), -np.inf, values))
    idx = abs_v.argmax(axis=axis)
    return np.take_along_axis(values, np.expand_dims(idx, axis), axis=axis).squeeze(axis)


def _save_modality_npz(out_dir: Path, set_name: str, modality: str,
                        rsids: list[str], scores: np.ndarray,
                        track_biosamples: list[str],
                        track_assays: list[str],
                        meta: dict,
                        extras: dict | None = None) -> Path:
    mod_dir = out_dir / f"alphagenome_{modality}"
    mod_dir.mkdir(parents=True, exist_ok=True)
    npz_path = mod_dir / f"{set_name}_snp_embeddings_{modality}.npz"
    payload = {
        "rsids":          np.array(rsids, dtype=object),
        "snp_idx":        np.arange(len(rsids), dtype=np.int64),
        "effect_scores":  scores.astype(np.float32, copy=False),
        "track_biosamples": np.array(track_biosamples, dtype=object),
        "track_assays":   np.array(track_assays, dtype=object),
        "model_tag":      f"alphagenome_{modality}",
        "model_version":  meta.get("model_version", ""),
        "scorer":         meta.get("scorer", ""),
        "aggregation":    meta.get("aggregation", ""),
        "interval_bp":    int(meta.get("interval_bp", 0)),
        # diff-attention compat shims
        "ref_emb":        np.zeros_like(scores, dtype=np.float32) if scores.ndim == 2
                            else np.zeros_like(scores.reshape(-1, 1), dtype=np.float32),
        "alt_emb":        (scores.astype(np.float32) if scores.ndim == 2
                            else scores.reshape(-1, 1).astype(np.float32)),
    }
    if extras:
        payload.update(extras)
    np.savez_compressed(npz_path, **payload)
    # Sidecar QC / run files
    (mod_dir / f"{set_name}_embed_qc_{modality}.json").write_text(
        json.dumps({**meta, "n_rsids": len(rsids), "n_tracks":
                     scores.shape[1] if scores.ndim == 2 else 1}, indent=2))
    (mod_dir / f"{set_name}_embed_run_{modality}.json").write_text(
        json.dumps({"timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
                    **meta}, indent=2))
    return npz_path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snps", required=True, type=Path,
                    help="recover_all_pool_snps.tsv (CHR + BP_GRCh38 + REF + ALT + rsID)")
    ap.add_argument("--set", default="recover_all_pool",
                    help="set name (used in output file names)")
    ap.add_argument("--outdir", required=True, type=Path,
                    help="parent dir; alphagenome_<modality> subdirs will be created")
    ap.add_argument("--api-key-file", default=None,
                    help="fallback path to a file containing the API key")
    ap.add_argument("--no-brain-filter", action="store_true",
                    help="dump all tracks instead of brain-filtered")
    ap.add_argument("--limit", type=int, default=None,
                    help="dev limit on # variants (default: all)")
    args = ap.parse_args()

    api_key = _resolve_api_key(args.api_key_file)
    print(f"[auth] API key resolved (len={len(api_key)})")

    # Lazy import — AlphaGenome must be pip-installed; do this AFTER auth so
    # a missing key fails fast.
    from alphagenome.models import dna_client, variant_scorers
    from alphagenome.data import genome

    print(f"[snps] reading {args.snps}")
    snps = pd.read_csv(args.snps, sep="\t", dtype=str, keep_default_na=False)
    if args.limit:
        snps = snps.head(args.limit)
    n = len(snps)
    print(f"[snps] {n} variants to score")

    model = dna_client.create(api_key)
    model_version = getattr(model, "version", "alphagenome_unknown")
    print(f"[model] dna_client created (version={model_version})")

    # Recommended scorers come from the module-level immutabledict
    # (get_recommended_scorers expects an Organism enum and returns a list;
    # the dict gives us name-keyed access we want here).
    rec = dict(variant_scorers.RECOMMENDED_VARIANT_SCORERS)
    print(f"[scorers] available recommended scorers: {sorted(rec.keys())}")

    # Pick the per-modality scorer keys defensively (skip if missing).
    wanted_keys = {
        "eqtl":           ("RNA_SEQ",),
        "accessibility":  ("ATAC", "DNASE"),
        "TSS":            ("CAGE", "PROCAP"),
        "splicing":       SPLICE_SCORER_KEYS,
        "chip_histone":   ("CHIP_HISTONE",),
        "chip_tf":        ("CHIP_TF",),
    }
    selected = {}
    for mod, keys in wanted_keys.items():
        selected[mod] = [rec[k] for k in keys if k in rec]
        present = [k for k in keys if k in rec]
        missing = [k for k in keys if k not in rec]
        print(f"  {mod:14s}: using {present} {'(missing: ' + ','.join(missing) + ')' if missing else ''}")

    # Flatten the union for the actual API call
    all_scorers = [s for slist in selected.values() for s in slist]

    # Per-variant accumulators (use Python lists then stack)
    rsids = []
    per_variant_tidy: list[pd.DataFrame] = []

    SEQUENCE_LENGTH_1MB = dna_client.SEQUENCE_LENGTH_1MB
    t0 = time.time()
    for i, row in snps.iterrows():
        rs = row["rsID"]
        chrom = f"chr{row['CHR']}".replace("chrchr", "chr")
        try:
            pos = int(row["BP_GRCh38"])
        except (KeyError, ValueError):
            pos = int(row.get("BP_pub", 0))
        ref = str(row.get("REF", row.get("other_allele", ""))).upper()
        alt = str(row.get("ALT", row.get("effect_allele", ""))).upper()
        if not ref or not alt or not pos:
            print(f"  [skip] {rs}: missing pos/ref/alt")
            continue

        variant = genome.Variant(chrom, pos, ref, alt)
        interval = genome.Interval(chrom, pos, pos + 1).resize(SEQUENCE_LENGTH_1MB)
        try:
            scores_list = model.score_variant(
                interval=interval, variant=variant,
                variant_scorers=all_scorers)
        except Exception as e:
            print(f"  [error] {rs}: {e}")
            continue

        # tidy_scores flattens AnnData list → long DataFrame
        try:
            tidy = variant_scorers.tidy_scores(scores_list)
        except Exception:
            # Fallback if tidy_scores isn't exposed: manual flatten
            tidy = _manual_tidy(scores_list)
        tidy["rsID"] = rs
        rsids.append(rs)
        per_variant_tidy.append(tidy)
        if (i + 1) % 10 == 0 or (i + 1) == n:
            dt = time.time() - t0
            print(f"  [{i+1}/{n}] {rs}  ({dt:.0f}s elapsed, "
                   f"{dt/(i+1):.1f}s/variant)", flush=True)

    if not per_variant_tidy:
        sys.exit("[ERROR] no variants scored successfully")

    big = pd.concat(per_variant_tidy, ignore_index=True)
    print(f"\n[tidy] long DataFrame: {len(big)} rows; "
          f"columns={list(big.columns)[:12]}…")
    print(f"[tidy] output_type values: "
          f"{big['output_type'].value_counts().to_dict() if 'output_type' in big.columns else 'no output_type col'}")

    # Apply brain filter at the row level (uses ontology_curie / biosample_type
    # / biosample_life_stage). Normalise NaN → "" first so the strict-filter
    # allow-list correctly accepts FANTOM CAGE rows (which AlphaGenome's
    # tidy_scores emits with `biosample_life_stage=NaN` instead of "").
    if not args.no_brain_filter:
        for c in ("ontology_curie", "biosample_type", "biosample_life_stage"):
            if c in big.columns:
                big[c] = big[c].fillna("")
        mask = big.apply(_is_brain_track, axis=1)
        print(f"[brain] kept {int(mask.sum())} / {len(big)} rows after brain filter")
        big = big.loc[mask].reset_index(drop=True)

    # Pivot helper: (rsID × (output_type, biosample_name)) → matrix
    def _pivot(df_sub: pd.DataFrame, score_col: str = "raw_score",
               agg: str = "max_abs") -> tuple[np.ndarray, list, list]:
        if df_sub.empty:
            return np.zeros((n, 0), dtype=np.float32), [], []
        if agg == "max_abs":
            df_sub = df_sub.copy()
            df_sub["_abs"] = df_sub[score_col].abs()
            df_sub = (df_sub.sort_values("_abs", ascending=False)
                       .drop_duplicates(subset=["rsID", "output_type",
                                                "biosample_name"]))
        pv = df_sub.pivot_table(index="rsID",
                                  columns=["output_type", "biosample_name"],
                                  values=score_col, aggfunc="first")
        pv = pv.reindex(snps["rsID"].tolist())
        scores = pv.fillna(0.0).values
        track_cols = list(pv.columns)
        track_assays = [str(c[0]) for c in track_cols]
        track_biosamples = [str(c[1]) for c in track_cols]
        return scores, track_biosamples, track_assays

    meta_base = {"model_version": str(model_version),
                  "interval_bp": int(SEQUENCE_LENGTH_1MB)}

    # ── 1) alphagenome_eqtl: RNA_SEQ rows; pivot by (track_name, track_strand)
    # so strand-specific ENCODE tracks stay as distinct columns. The 3 DLPFC
    # RNA_SEQ tracks split as: 2 ENCODE total RNA-seq (strand + and -) + 1
    # GTEx polyA plus RNA-seq (strand "."). The 1 Hippocampus track is GTEx
    # polyA plus RNA-seq (strand "."). Total = 4 columns.
    rna = big[big["output_type"] == "RNA_SEQ"] if "output_type" in big.columns else pd.DataFrame()
    if not rna.empty and "track_name" in rna.columns:
        rna = rna.assign(_track_id=rna["track_name"].astype(str) + "|strand=" +
                            rna.get("track_strand", pd.Series("", index=rna.index)).fillna("").astype(str))
        ag = (rna.assign(_abs=rna["raw_score"].abs())
                .sort_values("_abs", ascending=False)
                .drop_duplicates(subset=["rsID", "_track_id"]))
        pv = ag.pivot_table(index="rsID", columns="_track_id",
                              values="raw_score", aggfunc="first")
        pv = pv.reindex(snps["rsID"].tolist())
        eqtl_scores = pv.fillna(0.0).values.astype(np.float32)
        track_ids = list(pv.columns)
        tid_to_bs = dict(zip(rna["_track_id"].astype(str), rna["biosample_name"].astype(str)))
        eqtl_bio = [tid_to_bs.get(tid, "") for tid in track_ids]
        eqtl_asys = ["RNA_SEQ"] * len(track_ids)
    else:
        eqtl_scores, eqtl_bio, eqtl_asys = _pivot(rna, agg="max_abs")
    _save_modality_npz(args.outdir, args.set, "eqtl", snps["rsID"].tolist(),
                       eqtl_scores, eqtl_bio, eqtl_asys,
                       {**meta_base, "scorer": "RNA_SEQ_GeneMaskLFCScorer",
                        "aggregation": "max_abs_lfc_across_cis_genes_per_(track_name,track_strand)"})
    print(f"[eqtl] scores shape = {eqtl_scores.shape}")

    # ── 2) alphagenome_accessibility: DNase + ATAC; max_signed per biosample
    acc = big[big["output_type"].isin(["DNASE", "ATAC"])] if "output_type" in big.columns else pd.DataFrame()
    if not acc.empty:
        # Per (rsID, biosample) — take both DNase + ATAC, then max_signed
        ag = (acc.groupby(["rsID", "biosample_name", "output_type"])["raw_score"]
                .agg(lambda s: s.iloc[s.abs().argmax()])
                .reset_index())
        pv = ag.pivot_table(index="rsID", columns=["biosample_name", "output_type"],
                              values="raw_score", aggfunc="first")
        pv = pv.reindex(snps["rsID"].tolist())
        # For each (rsID × biosample), take max_signed across [DNase, ATAC]
        biosamples = sorted({b for b, _ in pv.columns})
        scores = np.full((len(snps), len(biosamples)), np.nan, dtype=np.float32)
        source_assay = []
        for j, b in enumerate(biosamples):
            cols = [(b, "DNASE"), (b, "ATAC")]
            present_cols = [c for c in cols if c in pv.columns]
            if not present_cols: continue
            sub = pv[present_cols].values
            scores[:, j] = _max_signed(sub, axis=1)
            # Which assay won?
            winners = []
            for r_idx in range(len(snps)):
                row_vals = sub[r_idx, :]
                if np.all(np.isnan(row_vals)):
                    winners.append("")
                else:
                    abs_vals = np.where(np.isnan(row_vals), -np.inf, np.abs(row_vals))
                    winners.append(present_cols[int(abs_vals.argmax())][1])
            source_assay.append(winners)
        source_assay_arr = (np.array(source_assay, dtype=object).T
                              if source_assay else np.zeros((len(snps), 0), dtype=object))
        scores = np.nan_to_num(scores)
        _save_modality_npz(args.outdir, args.set, "accessibility",
                           snps["rsID"].tolist(),
                           scores, biosamples, ["DNASE_or_ATAC"]*len(biosamples),
                           {**meta_base, "scorer": "DNASE+ATAC_CenterMaskScorer",
                            "aggregation": "max_signed(DNase,ATAC)_per_biosample"},
                           extras={"source_assay": source_assay_arr})
        print(f"[accessibility] scores shape = {scores.shape}")
    else:
        print(f"[accessibility] no DNase/ATAC rows after brain filter — skipping")

    # ── 3) alphagenome_TSS: CAGE + PRO-cap; pivot by (track_name, track_strand)
    # so strand-specific FANTOM CAGE tracks stay distinct. For the strict
    # allow-list, only Hippocampus FANTOM hCAGE survives (2 strand-specific
    # tracks); no PROCAP at allow-listed regions.
    tss = big[big["output_type"].isin(["CAGE", "PROCAP"])] if "output_type" in big.columns else pd.DataFrame()
    if not tss.empty:
        tss = tss.assign(_track_id=tss["track_name"].astype(str) + "|strand=" +
                            tss.get("track_strand", pd.Series("", index=tss.index)).fillna("").astype(str))
        ag = (tss.assign(_abs=tss["raw_score"].abs())
                .sort_values("_abs", ascending=False)
                .drop_duplicates(subset=["rsID", "_track_id"]))
        pv = ag.pivot_table(index="rsID", columns="_track_id",
                              values="raw_score", aggfunc="first")
        pv = pv.reindex(snps["rsID"].tolist())
        tss_scores = pv.fillna(0.0).values.astype(np.float32)
        track_ids = list(pv.columns)
        tid_to_bs = dict(zip(tss["_track_id"].astype(str), tss["biosample_name"].astype(str)))
        tid_to_ot = dict(zip(tss["_track_id"].astype(str), tss["output_type"].astype(str)))
        tss_bio = [tid_to_bs.get(tid, "") for tid in track_ids]
        tss_asys = [tid_to_ot.get(tid, "") for tid in track_ids]
        _save_modality_npz(args.outdir, args.set, "TSS",
                           snps["rsID"].tolist(),
                           tss_scores, tss_bio, tss_asys,
                           {**meta_base, "scorer": "CAGE+PROCAP_CenterMaskScorer",
                            "aggregation": "max_abs_per_(track_name,track_strand)"})
        print(f"[TSS] scores shape = {tss_scores.shape}")
    else:
        print(f"[TSS] no CAGE/PRO-cap rows after brain filter — skipping")

    # ── 4) alphagenome_splicing_merged: per-variant scalar formula
    spl = big[big["output_type"].isin(["SPLICE_SITES", "SPLICE_SITE_USAGE",
                                          "SPLICE_JUNCTIONS"])] if "output_type" in big.columns else pd.DataFrame()
    if not spl.empty:
        # Per (rsID, output_type), take max-abs across all rows (genes/sites)
        per_rs_per_kind = (spl.assign(_abs=spl["raw_score"].abs())
                             .groupby(["rsID", "output_type"])["raw_score"]
                             .agg(lambda s: s.iloc[s.abs().argmax()])
                             .unstack("output_type"))
        per_rs_per_kind = per_rs_per_kind.reindex(snps["rsID"].tolist())
        for col in ("SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"):
            if col not in per_rs_per_kind.columns:
                per_rs_per_kind[col] = 0.0
        s_sites = per_rs_per_kind["SPLICE_SITES"].abs().fillna(0).values
        s_usage = per_rs_per_kind["SPLICE_SITE_USAGE"].abs().fillna(0).values
        s_junc  = per_rs_per_kind["SPLICE_JUNCTIONS"].abs().fillna(0).values
        # Literal Python precedence: A + B + C/5
        merged = s_sites + s_usage + (s_junc / 5.0)
        # NB: also save (rsID × 3) component matrix for re-aggregation
        comp = np.stack([s_sites, s_usage, s_junc], axis=1).astype(np.float32)
        _save_modality_npz(args.outdir, args.set, "splicing_merged",
                           snps["rsID"].tolist(),
                           merged.reshape(-1, 1).astype(np.float32),
                           ["MERGED"], ["SPLICING"],
                           {**meta_base, "scorer": "SPLICE_SITES+SPLICE_SITE_USAGE+SPLICE_JUNCTIONS",
                            "aggregation": "max_sites + max_usage + max_junctions/5  "
                                          "(literal Python precedence)"},
                           extras={"components": comp,
                                    "component_names": np.array(
                                        ["SPLICE_SITES", "SPLICE_SITE_USAGE",
                                         "SPLICE_JUNCTIONS"], dtype=object)})
        print(f"[splicing_merged] saved (128,1) merged + (128,3) components")
    else:
        print(f"[splicing_merged] no splice rows after brain filter — skipping")

    # ── 5) alphagenome_chip_histone: per (biosample × histone_mark)
    chip_h = big[big["output_type"] == "CHIP_HISTONE"] if "output_type" in big.columns else pd.DataFrame()
    if not chip_h.empty:
        # Identify the histone_mark column (AlphaGenome tidy_scores carries
        # `histone_mark`; fall back to parsing `name` if absent).
        if "histone_mark" in chip_h.columns:
            chip_h = chip_h.assign(_mark=chip_h["histone_mark"].fillna(""))
        else:
            chip_h = chip_h.assign(_mark=chip_h["name"].fillna("").str.extract(
                r"Histone ChIP-seq (\S+)", expand=False).fillna(""))
        # Per (rsID, biosample, mark) — max-abs across any duplicates
        ag = (chip_h.assign(_abs=chip_h["raw_score"].abs())
                    .sort_values("_abs", ascending=False)
                    .drop_duplicates(subset=["rsID", "biosample_name", "_mark"]))
        pv = ag.pivot_table(index="rsID",
                              columns=["biosample_name", "_mark"],
                              values="raw_score", aggfunc="first")
        pv = pv.reindex(snps["rsID"].tolist())
        scores = pv.fillna(0.0).values.astype(np.float32)
        track_cols = list(pv.columns)
        track_biosamples = [str(c[0]) for c in track_cols]
        track_assays = [str(c[1]) for c in track_cols]  # histone marks
        _save_modality_npz(args.outdir, args.set, "chip_histone",
                           snps["rsID"].tolist(),
                           scores, track_biosamples, track_assays,
                           {**meta_base, "scorer": "CHIP_HISTONE_CenterMaskScorer",
                            "aggregation": "max_signed_per_(biosample,histone_mark)"})
        print(f"[chip_histone] scores shape = {scores.shape}")
    else:
        print(f"[chip_histone] no CHIP_HISTONE rows after strict filter — skipping")

    # ── 6) alphagenome_chip_tf: per (biosample × transcription_factor)
    chip_tf = big[big["output_type"] == "CHIP_TF"] if "output_type" in big.columns else pd.DataFrame()
    if not chip_tf.empty:
        if "transcription_factor" in chip_tf.columns:
            chip_tf = chip_tf.assign(_tf=chip_tf["transcription_factor"].fillna(""))
        else:
            chip_tf = chip_tf.assign(_tf=chip_tf["name"].fillna("").str.extract(
                r"TF ChIP-seq (\S+)", expand=False).fillna(""))
        ag = (chip_tf.assign(_abs=chip_tf["raw_score"].abs())
                     .sort_values("_abs", ascending=False)
                     .drop_duplicates(subset=["rsID", "biosample_name", "_tf"]))
        pv = ag.pivot_table(index="rsID",
                              columns=["biosample_name", "_tf"],
                              values="raw_score", aggfunc="first")
        pv = pv.reindex(snps["rsID"].tolist())
        scores = pv.fillna(0.0).values.astype(np.float32)
        track_cols = list(pv.columns)
        track_biosamples = [str(c[0]) for c in track_cols]
        track_assays = [str(c[1]) for c in track_cols]  # TF names
        _save_modality_npz(args.outdir, args.set, "chip_tf",
                           snps["rsID"].tolist(),
                           scores, track_biosamples, track_assays,
                           {**meta_base, "scorer": "CHIP_TF_CenterMaskScorer",
                            "aggregation": "max_signed_per_(biosample,transcription_factor)"})
        print(f"[chip_tf] scores shape = {scores.shape}")
    else:
        print(f"[chip_tf] no CHIP_TF rows after strict filter — skipping")

    # Long-form tidy backup for full provenance (gz to save space)
    big_path = args.outdir / f"alphagenome_eqtl/{args.set}_tidy_long.tsv.gz"
    big_path.parent.mkdir(parents=True, exist_ok=True)
    big.to_csv(big_path, sep="\t", index=False, compression="gzip")
    print(f"\n[provenance] full tidy long-form: {big_path}")

    print("\n=== AlphaGenome extraction complete ===")


def _manual_tidy(scores_list):
    """Fallback if variant_scorers.tidy_scores not available — manual flatten."""
    rows = []
    for ad in scores_list:
        try:
            X = ad.X
            obs = ad.obs.reset_index()
            var = ad.var.reset_index()
            for r_i in range(X.shape[0]):
                for c_i in range(X.shape[1]):
                    rows.append({
                        "gene_id": obs.iloc[r_i].get("gene_id", ""),
                        "biosample_name": var.iloc[c_i].get("biosample_name", ""),
                        "ontology_curie": var.iloc[c_i].get("ontology_curie", ""),
                        "output_type": var.iloc[c_i].get("output_type",
                                                          getattr(ad.uns, "get",
                                                                  dict().get)("scorer_name", "?")),
                        "raw_score": float(X[r_i, c_i]),
                    })
        except Exception:
            continue
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
