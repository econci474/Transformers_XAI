"""
diagnose_encoding.py
====================
Inspects a sample verbose baseline CSV to find the exact form of the
plasma biomarker string and diagnose what encoding artefact is present.
"""

import pandas as pd
from pathlib import Path

# Pick one verbose baseline CSV that should have plasma data
SAMPLE = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\verbose\baseline\seed_0\train.csv")

print(f"Reading: {SAMPLE}\n")
df = pd.read_csv(SAMPLE, low_memory=False)

print(f"Columns: {list(df.columns)}\n")

# Find rows whose Generated_Text mentions plasma
col = "Generated_Text"
mask = df[col].str.contains("Plasma|plasma|beta|Beta|42|NfL", na=False, regex=True)
sample_rows = df.loc[mask, col].dropna().head(3)

for i, text in sample_rows.items():
    # Find the relevant snippet
    idx = text.lower().find("plasma")
    if idx == -1:
        idx = 0
    snippet = text[max(0, idx-5):idx+120]
    print(f"--- Row {i} (chars {idx}-{idx+120}) ---")
    print(repr(snippet))   # repr() shows exact Unicode code points
    print(snippet)
    print()

# Also scan every text column for any non-ASCII characters
print("\n=== Non-ASCII character census ===")
for col in df.select_dtypes(include="object").columns:
    all_text = " ".join(df[col].dropna().astype(str))
    non_ascii = [(c, hex(ord(c)), i) for i, c in enumerate(all_text) if ord(c) > 127]
    if non_ascii:
        # Show unique chars
        unique = sorted(set((c, h) for c, h, _ in non_ascii))
        print(f"  Col={col!r}: {len(non_ascii)} non-ASCII chars, unique: {unique[:20]}")
