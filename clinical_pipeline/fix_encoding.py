"""
fix_encoding.py
===============
Replaces encoding artefacts in all CSV files under
D:\\ADNI_BIDS_project\\derivatives\\clinical (recursive).

Tries multiple candidate forms of the bad string since the exact
Unicode representation depends on how the CSV was originally written.
"""

import pandas as pd
from pathlib import Path

ROOT = Path(r"D:\ADNI_BIDS_project\derivatives\clinical")

# Candidates to try: (bad_string, good_string, description)
# Covers: Latin-1 mojibake form, actual Greek beta, and numeric suffixes
REPLACEMENTS = [
    # Form 1: mojibake — Î (U+00CE) + ² (U+00B2)
    ("AÎ²",  "Abeta",  "mojibake AÎ²"),
    # Form 2: actual Greek lowercase beta (U+03B2) as pandas may decode it
    ("A\u03b2", "Abeta", "Greek beta Aβ"),
]

csv_files = list(ROOT.rglob("*.csv"))
print(f"Found {len(csv_files)} CSV files under {ROOT}\n")

total_cells = 0
total_files = 0

for fp in csv_files:
    try:
        df = pd.read_csv(fp, low_memory=False)
    except Exception as e:
        print(f"  [ERROR reading] {fp.relative_to(ROOT)}: {e}")
        continue

    changed = False
    for col in df.select_dtypes(include="object").columns:
        for bad, good, desc in REPLACEMENTS:
            mask = df[col].str.contains(bad, na=False, regex=False)
            if mask.any():
                n = int(mask.sum())
                df[col] = df[col].str.replace(bad, good, regex=False)
                print(f"  {fp.relative_to(ROOT)}  |  col={col!r}  |  {n} cell(s) [{desc}] → '{good}'")
                total_cells += n
                changed = True

    if changed:
        try:
            df.to_csv(fp, index=False, encoding="utf-8")
            total_files += 1
        except Exception as e:
            print(f"  [ERROR writing] {fp.relative_to(ROOT)}: {e}")

print(f"\n{'='*60}")
print(f"Done.  Files modified: {total_files}  |  Total cells fixed: {total_cells}")
