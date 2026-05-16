"""
01_test_ntv2_tokeniser_lengths.py
==================================
Verify that 8,192 bp DNA sequences tokenise within the NT v2 2,048-token limit.

NT v2 uses non-overlapping 6-mer tokenisation (vocab=4,105).
  8,192 bp ÷ 6 bp/token ≈ 1,365 tokens + 1 CLS = ~1,366 tokens → fits 2,048.

Usage:
  python nt_v2/01_test_ntv2_tokeniser_lengths.py
"""

import pathlib
import random
import sys

try:
    from transformers import AutoTokenizer
except ImportError:
    print("[ERROR] transformers not installed. Run:  pip install transformers")
    sys.exit(1)

MODEL_ID = "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species"
MAX_TOKEN_LIMIT = 2048

print(f"Loading tokenizer: {MODEL_ID}")
tok = AutoTokenizer.from_pretrained(MODEL_ID, trust_remote_code=True)
print(f"  Vocab size       : {tok.vocab_size:,}")
print(f"  Model max length : {getattr(tok, 'model_max_length', 'N/A')}")

# ── Synthetic random DNA ─────────────────────────────────────────────────────
random.seed(42)
BASES = "ACGT"

test_lengths = [500, 1000, 2000, 4000, 6000, 8000, 8192, 10000, 12000]

print(f"\n{'BP length':>10}  {'Tokens':>8}  {'Fits 2048?':>10}  {'Utilisation':>11}")
print("-" * 55)

for bp_len in test_lengths:
    seq = "".join(random.choice(BASES) for _ in range(bp_len))
    ids = tok.encode(seq, add_special_tokens=True)
    n_tokens = len(ids)
    fits = n_tokens <= MAX_TOKEN_LIMIT
    util = n_tokens / MAX_TOKEN_LIMIT * 100
    marker = "  ✓" if fits else "  ✗ EXCEEDS"
    print(f"{bp_len:>10,}  {n_tokens:>8,}  {marker:>10}  {util:>10.1f}%")

# ── Real FASTA sequence ──────────────────────────────────────────────────────
fasta_path = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
if fasta_path.exists():
    try:
        from pyfaidx import Fasta
        fa = Fasta(str(fasta_path), rebuild=False)
        real_seq = str(fa["1"][100000:108192]).upper()
        ids = tok.encode(real_seq, add_special_tokens=True)
        print(f"\n  Real chr1 [100000:108192] (8,192 bp):")
        print(f"    → {len(ids):,} tokens  "
              f"({'✓ fits' if len(ids) <= MAX_TOKEN_LIMIT else '✗ EXCEEDS'} 2048)")
    except Exception as e:
        print(f"\n  (Could not test real FASTA sequence: {e})")

# ── Actual by_chrom_combos dataset ───────────────────────────────────────────
combo_csv = pathlib.Path(
    "D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/"
    "bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom_combos/"
    "forward/train.csv"
)
if combo_csv.exists():
    import pandas as pd
    print(f"\n  Testing by_chrom_combos: forward/train.csv")

    cols = pd.read_csv(combo_csv, nrows=0).columns.tolist()
    use = ["sequence"] + (["win_len_bp"] if "win_len_bp" in cols else [])
    df = pd.read_csv(combo_csv, usecols=use)
    print(f"    Total rows: {len(df):,}")

    sample = df.sample(n=min(200, len(df)), random_state=42)
    token_counts = []
    exceeds = 0
    for _, row in sample.iterrows():
        ids = tok.encode(str(row["sequence"]), add_special_tokens=True)
        n_tok = len(ids)
        bp_len = int(row.get("win_len_bp", len(str(row["sequence"]))))
        token_counts.append((bp_len, n_tok))
        if n_tok > MAX_TOKEN_LIMIT:
            exceeds += 1

    bps = [bp for bp, _ in token_counts]
    toks = [nt for _, nt in token_counts]
    print(f"    Sampled {len(token_counts)} rows:")
    print(f"      bp range   : {min(bps):>5,}–{max(bps):>5,}")
    print(f"      token range: {min(toks):>5,}–{max(toks):>5,}  (avg {sum(toks)/len(toks):,.0f})")
    print(f"    Exceeds 2048? {'✗ YES — ' + str(exceeds) + ' rows!' if exceeds else '✓ none'}")
else:
    print(f"\n  (by_chrom_combos CSV not found — skipped)")

# ── Summary ──────────────────────────────────────────────────────────────────
seq_8k = "".join(random.choice(BASES) for _ in range(8192))
n_8k = len(tok.encode(seq_8k, add_special_tokens=True))
print(f"\n{'='*55}")
print(f"  8,192 bp → {n_8k:,} tokens")
print(f"  2,048 token limit → {'SAFE' if n_8k <= MAX_TOKEN_LIMIT else 'TRUNCATION NEEDED'}")
if n_8k <= MAX_TOKEN_LIMIT:
    print(f"  Headroom: {MAX_TOKEN_LIMIT - n_8k:,} tokens ({(MAX_TOKEN_LIMIT - n_8k) / MAX_TOKEN_LIMIT * 100:.1f}%)")
print(f"  Compression ratio: ~{8192 / n_8k:.1f} bp per token")
print(f"{'='*55}")
