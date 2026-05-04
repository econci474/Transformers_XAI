"""Check Generated_Text context length statistics."""
import pandas as pd

PATH = r"D:\ADNI_BIDS_project\derivatives\clinical\verbose\longitudinal\master_clinical_verbose.csv"
df = pd.read_csv(PATH, low_memory=False)

texts = df["Generated_Text"].dropna().astype(str)

chars  = texts.str.len()
words  = texts.str.split().str.len()
approx_tokens = chars / 4   # GPT-style: ~4 chars per token

print(f"Rows with Generated_Text: {len(texts)}")
print()
print(f"{'Metric':<25} {'Min':>8} {'Median':>8} {'Mean':>8} {'Max':>8} {'P95':>8}")
print("-" * 65)
for label, series in [("Characters", chars), ("Words", words), ("Tokens (≈chars/4)", approx_tokens)]:
    print(f"{label:<25} {series.min():>8.0f} {series.median():>8.0f} "
          f"{series.mean():>8.0f} {series.max():>8.0f} {series.quantile(0.95):>8.0f}")

print(f"\nNote: Token count is approximate (chars/4). For exact GPT-4o / LLaMA tokenisation,")
print(f"install tiktoken:  pip install tiktoken")
print(f"  import tiktoken; enc = tiktoken.get_encoding('cl100k_base')")
print(f"  tokens = texts.apply(lambda t: len(enc.encode(t)))")
