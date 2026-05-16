"""
03_extract_patient_embeddings.py
=================================
Extract per-SNP-window mean-pooled embeddings from frozen Caduceus models.

Caduceus uses single-nucleotide tokenisation (1 token ≈ 1 bp), so our
8,192 bp patient windows → ~8,192 tokens, well within the 131k context.

Two variants are supported:
  --variant ps  → caduceus-ps_seqlen-131k_d_model-256_n_layer-16
                  (parameter sharing; RC equivariant)
  --variant ph  → caduceus-ph_seqlen-131k_d_model-256_n_layer-16
                  (post-hoc conjoining; requires forward + RC at inference)

Output format matches the BMFM / NT v2 pipeline: one .npz per model variant
with ``mean__<PTID>`` arrays (B, N_windows, H=256) and ``window_ids``.

Usage:
  python caduceus/03_extract_patient_embeddings.py \\
      --sequences D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/patient_window_sequences/sequences/all_patients.csv \\
      --outdir    D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/patient_embeddings \\
      --variant ps \\
      --device cuda --batch-size 4
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm


VARIANTS = {
    "ps": "kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16",
    "ph": "kuleshov-group/caduceus-ph_seqlen-131k_d_model-256_n_layer-16",
}


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def mean_pool(hidden_states: torch.Tensor,
              attention_mask: torch.Tensor) -> torch.Tensor:
    """Mean-pool last hidden state, respecting the attention mask."""
    mask = attention_mask.unsqueeze(-1).float()            # (B, L, 1)
    return (hidden_states * mask).sum(dim=1) / mask.sum(dim=1).clamp_min(1)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--sequences", type=Path, required=True,
                    help="all_patients.csv from 13_build_patient_window_sequences.py")
    ap.add_argument("--outdir", type=Path, required=True,
                    help="Output root; subdir per variant is created.")
    ap.add_argument("--variant", choices=list(VARIANTS.keys()), required=True,
                    help="'ps' = parameter-sharing (RC equivariant); "
                         "'ph' = post-hoc conjoining (needs fwd + RC)")
    ap.add_argument("--batch-size", type=int, default=4)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--limit", type=int, default=None,
                    help="DEBUG: only process first N rows.")
    args = ap.parse_args()

    if not args.sequences.exists():
        sys.exit(f"[ERROR] sequences CSV not found: {args.sequences}")

    model_id = VARIANTS[args.variant]
    tag = f"caduceus_{args.variant}"
    out_dir = args.outdir / tag
    out_dir.mkdir(parents=True, exist_ok=True)
    out_npz = out_dir / "embeddings.npz"

    print(f"[model]    {model_id}")
    print(f"[variant]  {args.variant}")
    print(f"[output]   {out_npz}")

    # ── Load model + tokenizer ────────────────────────────────────────────────
    from transformers import AutoTokenizer, AutoModelForMaskedLM

    device = torch.device(
        args.device if torch.cuda.is_available() or args.device == "cpu"
        else "cpu"
    )
    print(f"\n  torch = {torch.__version__}; cuda = {torch.cuda.is_available()}")
    print(f"  device = {device}")

    print(f"\n  Loading tokenizer ...")
    tokenizer = AutoTokenizer.from_pretrained(model_id, trust_remote_code=True)
    print(f"  Vocab size: {tokenizer.vocab_size}")

    print(f"  Loading model ...")
    model = AutoModelForMaskedLM.from_pretrained(
        model_id, trust_remote_code=True,
    )
    model = model.to(device).eval()

    # Caduceus single-char tokenisation: no special max_length concern for 8k bp
    # (model supports 131k context). We pad to longest-in-batch for efficiency.

    # ── Load patient sequences ────────────────────────────────────────────────
    print(f"\nLoading sequences: {args.sequences}")
    df = pd.read_csv(args.sequences)
    if args.limit:
        df = df.head(args.limit)
    print(f"  {len(df)} rows  ({df['PTID'].nunique()} patients × "
          f"{df['window_id'].nunique()} windows)")

    window_ids_unique = sorted(df["window_id"].unique())
    win_idx = {w: i for i, w in enumerate(window_ids_unique)}
    n_windows = len(window_ids_unique)

    hidden_dim = None
    buffers: dict[str, np.ndarray] = {}

    # ── Forward in batches ────────────────────────────────────────────────────
    n_batches = (len(df) + args.batch_size - 1) // args.batch_size
    print(f"\nForward pass: {n_batches} batches of {args.batch_size}")

    is_ph = args.variant == "ph"
    if is_ph:
        print("  [info] ph variant: averaging fwd + RC embeddings")

    with torch.no_grad():
        for b_idx in tqdm(range(n_batches), desc="Embedding"):
            chunk = df.iloc[b_idx * args.batch_size:(b_idx + 1) * args.batch_size]
            seqs = chunk["sequence"].tolist()
            ptids = chunk["PTID"].tolist()
            wids = chunk["window_id"].tolist()

            # Tokenize forward sequences
            tokenized = tokenizer(
                seqs,
                return_tensors="pt",
                padding=True,
                truncation=False,
            )
            input_ids = tokenized["input_ids"].to(device)
            attention_mask = tokenized["attention_mask"].to(device)

            out = model(
                input_ids,
                attention_mask=attention_mask,
                output_hidden_states=True,
            )
            last_hidden = out.hidden_states[-1]
            emb_fwd = mean_pool(last_hidden, attention_mask)

            if is_ph:
                # Post-hoc variant: also run reverse complement and average
                rc_seqs = [reverse_complement(s) for s in seqs]
                rc_tok = tokenizer(
                    rc_seqs,
                    return_tensors="pt",
                    padding=True,
                    truncation=False,
                )
                rc_ids = rc_tok["input_ids"].to(device)
                rc_mask = rc_tok["attention_mask"].to(device)

                rc_out = model(
                    rc_ids,
                    attention_mask=rc_mask,
                    output_hidden_states=True,
                )
                rc_hidden = rc_out.hidden_states[-1]
                emb_rc = mean_pool(rc_hidden, rc_mask)

                # Average fwd + RC
                emb = (emb_fwd + emb_rc) / 2.0
            else:
                # PS variant: RC equivariant — forward pass only
                emb = emb_fwd

            emb_np = emb.detach().to(torch.float32).cpu().numpy()

            if hidden_dim is None:
                hidden_dim = emb_np.shape[1]
                print(f"  hidden_dim = {hidden_dim}")

            for ptid, wid, vec in zip(ptids, wids, emb_np):
                if ptid not in buffers:
                    buffers[ptid] = np.full((n_windows, hidden_dim),
                                            np.nan, dtype=np.float32)
                buffers[ptid][win_idx[wid]] = vec

    # ── Save (dual-pool format: mean__<PTID> keys) ───────────────────────────
    print(f"\nWriting {out_npz}")
    save_dict = {"window_ids": np.array(window_ids_unique),
                 "upstream": np.array(tag)}
    for ptid, arr in buffers.items():
        save_dict[f"mean__{ptid}"] = arr
    np.savez_compressed(out_npz, **save_dict)
    print(f"  {len(buffers)} patient arrays, each {n_windows} × {hidden_dim}")

    nan_pat = [pt for pt, arr in buffers.items() if np.isnan(arr).any()]
    if nan_pat:
        print(f"  [warn] {len(nan_pat)} patients have at least one NaN window-row")


if __name__ == "__main__":
    main()
"""Caduceus extraction script for frozen patient embeddings."""
