"""
Quick test: how many BPE tokens does the BMFM-DNA-SNP tokeniser produce
for DNA sequences of various lengths (including the 8,192 bp max window)?

The BMFM tokenizer lives in a non-standard path on HuggingFace:
  ibm-research/biomed.dna.snp.modernbert.113m.v1 → tokenizers/dna_chunks/

Usage:
  python snp_pipeline/test_tokeniser_lengths.py
"""
import sys
import random
import pathlib

try:
    from transformers import PreTrainedTokenizerFast
except ImportError:
    print("[ERROR] transformers not installed. Run:  pip install transformers")
    sys.exit(1)

# ── Download & load the tokenizer from the correct subdirectory ───────────────
MODEL_ID = "ibm-research/biomed.dna.snp.modernbert.113m.v1"
MAX_TOKEN_LIMIT = 2048

print(f"Downloading tokenizer files from: {MODEL_ID}")

try:
    from huggingface_hub import hf_hub_download
except ImportError:
    print("[ERROR] huggingface_hub not installed. Run:  pip install huggingface_hub")
    sys.exit(1)

# Download the tokenizer.json from the tokenizers/dna_chunks/ subdirectory
tok_json = hf_hub_download(
    repo_id=MODEL_ID,
    filename="tokenizers/dna_chunks/tokenizer.json",
)
tok_dir = str(pathlib.Path(tok_json).parent)
print(f"  Tokenizer dir: {tok_dir}")

tok = PreTrainedTokenizerFast(tokenizer_file=tok_json)
print(f"  Vocab size    : {tok.vocab_size:,}")

# ── Generate test sequences ──────────────────────────────────────────────────
random.seed(42)
BASES = "ACGT"

test_lengths = [500, 1000, 2000, 4000, 6000, 8000, 8192, 10000, 12000]

print(f"\n{'BP length':>10}  {'Tokens':>8}  {'Fits 2048?':>10}  {'Utilisation':>11}")
print("-" * 55)

for bp_len in test_lengths:
    seq = "".join(random.choice(BASES) for _ in range(bp_len))
    ids = tok.encode(seq, add_special_tokens=False)
    n_tokens = len(ids)
    fits = n_tokens <= MAX_TOKEN_LIMIT
    util = n_tokens / MAX_TOKEN_LIMIT * 100
    marker = "  ✓" if fits else "  ✗ EXCEEDS"
    print(f"{bp_len:>10,}  {n_tokens:>8,}  {marker:>10}  {util:>10.1f}%")

# ── Also test with a real FASTA sequence if available ─────────────────────────
fasta_path = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
if fasta_path.exists():
    try:
        from pyfaidx import Fasta
        fa = Fasta(str(fasta_path), rebuild=False)
        real_seq = str(fa["1"][100000:108192]).upper()
        ids = tok.encode(real_seq, add_special_tokens=False)
        print(f"\n  Real chr1 [100000:108192] (8192 bp):")
        print(f"    → {len(ids):,} tokens  "
              f"({'✓ fits' if len(ids) <= MAX_TOKEN_LIMIT else '✗ EXCEEDS'} 2048)")
    except Exception as e:
        print(f"\n  (Could not test real FASTA sequence: {e})")

# ── Test the by_chrom_combos dataset (8192 bp windows) ────────────────────────
combo_csv = pathlib.Path(
    "D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/"
    "bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom_combos/"
    "forward/train.csv"
)
if combo_csv.exists():
    import pandas as pd
    print(f"\n  Testing by_chrom_combos: {combo_csv.parent.parent.name}/forward/train.csv")

    # Only load the columns we need
    cols_available = pd.read_csv(combo_csv, nrows=0).columns.tolist()
    use_cols = ["sequence"]
    if "win_len_bp" in cols_available:
        use_cols.append("win_len_bp")

    df = pd.read_csv(combo_csv, usecols=use_cols)
    print(f"    Total rows: {len(df):,}")
    if "win_len_bp" in df.columns:
        print(f"    BP lengths: {df['win_len_bp'].min():,}–{df['win_len_bp'].max():,} "
              f"(mean {df['win_len_bp'].mean():,.0f})")

    # Sample up to 200 rows
    sample = df.sample(n=min(200, len(df)), random_state=42)

    token_counts = []
    exceeds = 0
    for _, row in sample.iterrows():
        ids = tok.encode(str(row["sequence"]), add_special_tokens=False)
        n_tok = len(ids)
        bp_len = int(row.get("win_len_bp", len(str(row["sequence"]))))
        token_counts.append((bp_len, n_tok))
        if n_tok > MAX_TOKEN_LIMIT:
            exceeds += 1

    bps  = [bp for bp, _ in token_counts]
    toks = [nt for _, nt in token_counts]
    print(f"    Sampled {len(token_counts)} rows:")
    print(f"      bp range   : {min(bps):>5,}–{max(bps):>5,}")
    print(f"      token range: {min(toks):>5,}–{max(toks):>5,}  (avg {sum(toks)/len(toks):,.0f})")
    print(f"    Exceeds 2048? {'✗ YES — ' + str(exceeds) + ' rows!' if exceeds else '✓ none'}")
else:
    print(f"\n  (by_chrom_combos CSV not found — skipped)")

# ── Summary ───────────────────────────────────────────────────────────────────
# Recompute for clean summary
seq_8k = "".join(random.choice(BASES) for _ in range(8192))
n_8k = len(tok.encode(seq_8k, add_special_tokens=False))
print(f"\n{'='*55}")
print(f"  8,192 bp → {n_8k:,} tokens")
print(f"  2,048 token limit → {'SAFE' if n_8k <= 2048 else 'TRUNCATION NEEDED'}")
if n_8k <= 2048:
    print(f"  Headroom: {2048 - n_8k:,} tokens ({(2048 - n_8k) / 2048 * 100:.1f}%)")
print(f"  Compression ratio: ~{8192 / n_8k:.1f} bp per token")
print(f"{'='*55}")
