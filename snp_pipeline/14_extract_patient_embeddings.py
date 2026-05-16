"""
14_extract_patient_embeddings.py
=================================
Run a frozen BMFM backbone over the per-patient × per-window DNA sequences
produced by 13_build_patient_window_sequences.py, and cache the CLS embeddings.

For each upstream model we want to compare, run this script once. The output
is a single .npz file per upstream tag containing one (n_windows, hidden_dim)
tensor per patient — a flat dict keyed by PTID, plus a parallel `window_ids`
vector so the column order is recoverable.

Supported upstream identifiers
------------------------------
  --upstream base_ref      → ibm-research/biomed.dna.ref.modernbert.113m.v1
                             (no fine-tune; HF download)
  --upstream base_snp      → ibm-research/biomed.dna.snp.modernbert.113m.v1
                             (no fine-tune; HF download)
  --upstream <ckpt_path>   → any local Lightning .ckpt (e.g. an L4 regression
                             fine-tune). The script auto-loads via BMFM main.

Compute
-------
Designed for the L4 dev-gpu-acs box (or any single CUDA device with ≥16 GB).
fp32 on the L4 is fine since this is a forward-only pass — no autograd graph,
no backward. ~15 k sequences × ~5 ms forward = ~75 s total per upstream at
batch_size=4 on the L4 GPU.

Output
------
  bmfm_inputs/patient_embeddings/<upstream_tag>/embeddings.npz
    keys:
      window_ids : str[n_windows]
      <PTID>     : float32[n_windows, hidden_dim]    (one per patient)

Run on dev-gpu-acs
------------------
  scp sequences/all_patients.csv ec474@dev-gpu-acs:~/ADNI_SNP/patient_seqs/  # one-off
  ssh dev-gpu-acs
  conda activate bmfm
  python ~/Transformers_XAI/snp_pipeline/14_extract_patient_embeddings.py \\
      --sequences ~/ADNI_SNP/patient_seqs/all_patients.csv \\
      --upstream base_ref \\
      --outdir ~/ADNI_SNP/patient_embeddings/

Then scp the .npz file back to local for downstream training.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ── Upstream alias map ───────────────────────────────────────────────────────
UPSTREAM_ALIASES = {
    "base_ref": "ibm-research/biomed.dna.ref.modernbert.113m.v1",
    "base_snp": "ibm-research/biomed.dna.snp.modernbert.113m.v1",
}


def resolve_upstream(arg: str) -> tuple[str, str]:
    """Return (tag, checkpoint_arg_for_bmfm)."""
    if arg in UPSTREAM_ALIASES:
        return arg, UPSTREAM_ALIASES[arg]
    p = Path(arg)
    if p.is_file() or p.is_dir():
        tag = p.stem if p.is_file() else p.name
        return tag, str(p)
    sys.exit(f"[ERROR] unknown --upstream alias and not a path: {arg}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sequences", type=Path, required=True,
                    help="all_patients.csv from 13_build_patient_window_sequences.py")
    ap.add_argument("--upstream", type=str, required=True,
                    help="One of base_ref / base_snp, or a path to a Lightning .ckpt")
    ap.add_argument("--outdir", type=Path, required=True,
                    help="Output directory; one subdir per upstream tag is created.")
    ap.add_argument("--batch-size", type=int, default=4)
    ap.add_argument("--device", default="cuda",
                    help="cuda | cpu | cuda:0")
    ap.add_argument("--max-length", type=int, default=2048,
                    help="Tokeniser max_length / padding (default: 2048; matches model pre-training).")
    ap.add_argument("--limit", type=int, default=None,
                    help="DEBUG: only process first N sequences.")
    args = ap.parse_args()

    if not args.sequences.exists():
        sys.exit(f"[ERROR] sequences CSV not found: {args.sequences}")

    tag, ckpt = resolve_upstream(args.upstream)
    out_dir = args.outdir / tag
    out_dir.mkdir(parents=True, exist_ok=True)
    out_npz = out_dir / "embeddings.npz"
    print(f"[upstream]   tag={tag}   ckpt={ckpt}")
    print(f"[output]     {out_npz}")

    # ── Load BMFM model + tokenizer via the package's helpers ────────────────
    print("\nImporting bmfm_targets …")
    import torch
    from bmfm_targets.config.main_config import load_checkpoint_configs
    from bmfm_targets.models.model_utils import load_bmfm_model
    from bmfm_targets.tokenization.multifield_tokenizer import MultiFieldTokenizer

    print(f"  torch = {torch.__version__}; cuda available = {torch.cuda.is_available()}")
    device = torch.device(args.device if torch.cuda.is_available() or args.device == "cpu" else "cpu")
    print(f"  device   = {device}")

    # The cleanest way to load both the base IBM checkpoints and our Lightning
    # .ckpt files is via BMFM's `load_checkpoint_configs` (handles both HF IDs
    # and local .ckpt paths) followed by `load_bmfm_model`. If those helpers
    # are missing or renamed in the installed package, fall back to a direct
    # HF AutoModel path.
    try:
        cfg = load_checkpoint_configs(ckpt)
        model = load_bmfm_model(cfg)
        # Force SDPA (avoid the FA2 multi-field bug — see force_sdpa_wrapper.py)
        if hasattr(model, "config"):
            model.config._attn_implementation = "sdpa"
        tokenizer = MultiFieldTokenizer.from_pretrained(ckpt) \
            if hasattr(MultiFieldTokenizer, "from_pretrained") else None
    except Exception as e:
        sys.exit(
            f"[ERROR] failed to load BMFM model via package helpers: {e!r}\n"
            f"  Check that ckpt='{ckpt}' is either a valid HF id or a local .ckpt path,\n"
            f"  and that the bmfm conda env has biomed-multi-omic installed."
        )

    if tokenizer is None:
        # Fall back: tokenizer files usually shipped alongside the checkpoint.
        from bmfm_targets.models.model_utils import load_tokenizer
        tokenizer = load_tokenizer(cfg)

    model.eval().to(device)

    # ── Load patient sequences ───────────────────────────────────────────────
    print(f"\nLoading sequences: {args.sequences}")
    df = pd.read_csv(args.sequences)
    if args.limit:
        df = df.head(args.limit)
    print(f"  {len(df)} rows  ({df['PTID'].nunique()} patients × "
          f"{df['window_id'].nunique()} windows)")

    window_ids_unique = sorted(df["window_id"].unique())
    win_idx = {w: i for i, w in enumerate(window_ids_unique)}
    n_windows = len(window_ids_unique)

    # Pre-allocate per-patient buffers.
    hidden_dim = None
    buffers: dict[str, np.ndarray] = {}

    # ── Forward in batches ───────────────────────────────────────────────────
    n_batches = (len(df) + args.batch_size - 1) // args.batch_size
    print(f"\nForward pass: {n_batches} batches of {args.batch_size}")

    with torch.no_grad():
        for b_idx in range(n_batches):
            chunk = df.iloc[b_idx * args.batch_size:(b_idx + 1) * args.batch_size]
            seqs = chunk["sequence"].tolist()
            ptids = chunk["PTID"].tolist()
            wids = chunk["window_id"].tolist()

            # BMFM tokenizer expects a list of dicts with the field name.
            tokenized = tokenizer(
                [{"dna_chunks": s} for s in seqs],
                padding="max_length",
                truncation=True,
                max_length=args.max_length,
                return_tensors="pt",
            )
            input_ids = tokenized["input_ids"].to(device)
            attention_mask = tokenized["attention_mask"].to(device)

            out = model(input_ids=input_ids, attention_mask=attention_mask)
            # The BMFM model returns a dict-like with "embeddings" (CLS-pooled).
            emb = out["embeddings"] if isinstance(out, dict) else out.embeddings
            emb_np = emb.detach().to(torch.float32).cpu().numpy()  # (B, H)

            if hidden_dim is None:
                hidden_dim = emb_np.shape[1]
                print(f"  hidden_dim = {hidden_dim}")

            for ptid, wid, vec in zip(ptids, wids, emb_np):
                if ptid not in buffers:
                    buffers[ptid] = np.full((n_windows, hidden_dim),
                                            np.nan, dtype=np.float32)
                buffers[ptid][win_idx[wid]] = vec

            if (b_idx + 1) % 20 == 0 or b_idx == n_batches - 1:
                print(f"  batch {b_idx + 1}/{n_batches}")

    # ── Save ─────────────────────────────────────────────────────────────────
    print(f"\nWriting {out_npz}")
    np.savez_compressed(
        out_npz,
        window_ids=np.array(window_ids_unique),
        **buffers,
    )
    print(f"  {len(buffers)} patient arrays, each {n_windows} × {hidden_dim}")

    # Sanity: any patient with a NaN row means a missing (PTID, window) pair.
    nan_pat = [pt for pt, arr in buffers.items() if np.isnan(arr).any()]
    if nan_pat:
        print(f"  [warn] {len(nan_pat)} patients have at least one NaN window-row")


if __name__ == "__main__":
    main()
