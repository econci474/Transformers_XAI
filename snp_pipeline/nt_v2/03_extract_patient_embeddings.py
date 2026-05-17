"""
03_extract_patient_embeddings.py
=================================
Extract per-SNP-window mean-pooled embeddings from frozen NT v2.

Usage:
  python nt_v2/03_extract_patient_embeddings.py \
      --sequences D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/patient_window_sequences/all_patients.csv \
      --outdir    D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/patient_embeddings/ntv2 \
      --device    cuda --batch-size 8
"""

import argparse
import pathlib

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm
from transformers import AutoTokenizer, AutoModelForMaskedLM


MODEL_ID = "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species"


def mean_pool(hidden_states, attention_mask):
    """Mean-pool last hidden state, respecting the attention mask."""
    mask = attention_mask.unsqueeze(-1).float()           # (B, L, 1)
    return (hidden_states * mask).sum(dim=1) / mask.sum(dim=1)  # (B, H)


def main():
    parser = argparse.ArgumentParser(description="Extract NT v2 embeddings from patient window sequences.")
    parser.add_argument("--sequences", required=True, help="CSV with columns: PTID, window_id, sequence")
    parser.add_argument("--outdir", required=True, help="Output directory for embeddings.npz")
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--max-length", type=int, default=2048,
                        help="Max token length (NT v2 context window)")
    parser.add_argument("--limit", type=int, default=None, help="Process only first N rows (for testing)")
    args = parser.parse_args()

    out_dir = pathlib.Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_npz = out_dir / "embeddings.npz"

    # ── Load model + tokenizer ────────────────────────────────────────────────
    print(f"Loading model: {MODEL_ID}")
    tokenizer = AutoTokenizer.from_pretrained(MODEL_ID, trust_remote_code=True)
    model = AutoModelForMaskedLM.from_pretrained(MODEL_ID, trust_remote_code=True)
    model = model.to(args.device).eval()
    # We only need the backbone, not the MLM head — but we extract hidden states
    print(f"  Device: {args.device}")
    print(f"  Max length: {args.max_length}")

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

    # Pre-allocate per-patient buffers
    hidden_dim = None
    buffers: dict[str, np.ndarray] = {}

    # ── Resume from checkpoint if available ───────────────────────────────────
    ckpt_path = out_dir / "checkpoint.npz"
    n_batches = (len(df) + args.batch_size - 1) // args.batch_size
    start_batch = 0
    if ckpt_path.exists():
        print(f"\n  Resuming from checkpoint: {ckpt_path}")
        ckpt = np.load(ckpt_path, allow_pickle=True)
        start_batch = int(ckpt["completed_batches"])
        hidden_dim = int(ckpt["hidden_dim"])
        for k in ckpt.files:
            if k.startswith("mean__"):
                ptid = k[len("mean__"):]
                buffers[ptid] = ckpt[k]
        print(f"  Loaded {len(buffers)} patients, resuming from batch {start_batch}/{n_batches}")

    CKPT_EVERY = 500

    # ── Forward in batches ────────────────────────────────────────────────────
    print(f"\nForward pass: {n_batches} batches of {args.batch_size}"
          f" (starting from batch {start_batch})")

    with torch.no_grad():
        for b_idx in tqdm(range(start_batch, n_batches), desc="Embedding",
                          initial=start_batch, total=n_batches):
            chunk = df.iloc[b_idx * args.batch_size:(b_idx + 1) * args.batch_size]
            seqs = chunk["sequence"].tolist()
            ptids = chunk["PTID"].tolist()
            wids = chunk["window_id"].tolist()

            tokenized = tokenizer.batch_encode_plus(
                seqs,
                return_tensors="pt",
                padding="max_length",
                max_length=args.max_length,
                truncation=True,
            )
            input_ids = tokenized["input_ids"].to(args.device)
            attention_mask = (input_ids != tokenizer.pad_token_id).long()

            out = model(
                input_ids,
                attention_mask=attention_mask,
                encoder_attention_mask=attention_mask,
                output_hidden_states=True,
            )
            # Mean-pool last hidden state
            last_hidden = out.hidden_states[-1]       # (B, L, H)
            emb = mean_pool(last_hidden, attention_mask)  # (B, H)
            emb_np = emb.detach().to(torch.float32).cpu().numpy()

            if hidden_dim is None:
                hidden_dim = emb_np.shape[1]
                print(f"  hidden_dim = {hidden_dim}")

            for ptid, wid, vec in zip(ptids, wids, emb_np):
                if ptid not in buffers:
                    buffers[ptid] = np.full((n_windows, hidden_dim),
                                            np.nan, dtype=np.float32)
                buffers[ptid][win_idx[wid]] = vec

            # Periodic checkpoint
            if (b_idx + 1) % CKPT_EVERY == 0:
                _save_ckpt = {"completed_batches": np.array(b_idx + 1),
                              "hidden_dim": np.array(hidden_dim),
                              "window_ids": np.array(window_ids_unique)}
                for pt, arr in buffers.items():
                    _save_ckpt[f"mean__{pt}"] = arr
                np.savez_compressed(ckpt_path, **_save_ckpt)
                tqdm.write(f"  [ckpt] saved at batch {b_idx + 1}/{n_batches}")

    # ── Save ──────────────────────────────────────────────────────────────────
    print(f"\nWriting {out_npz}")
    save_dict = {"window_ids": np.array(window_ids_unique)}
    for ptid, arr in buffers.items():
        save_dict[f"mean__{ptid}"] = arr
    np.savez_compressed(out_npz, **save_dict)
    print(f"  {len(buffers)} patient arrays, each {n_windows} × {hidden_dim}")

    # Clean up checkpoint
    if ckpt_path.exists():
        ckpt_path.unlink()
        print("  Removed checkpoint (final output saved).")

    nan_pat = [pt for pt, arr in buffers.items() if np.isnan(arr).any()]
    if nan_pat:
        print(f"  [warn] {len(nan_pat)} patients have at least one NaN window-row")


if __name__ == "__main__":
    main()
