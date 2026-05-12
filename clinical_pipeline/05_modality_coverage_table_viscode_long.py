"""
05_modality_coverage_table.py
==============================
CLI flags
---------
  --aligned   Use the clinical-aligned master (bids_ses replaced by clinical
              VISCODE_long for matched scans). Output goes to a parallel
              folder so the two views can be compared side by side.


Cross-modality coverage audit per ADNI visit code, restricted to the cohort
with both SNP and MRI data (subjects_with_snp_and_mri.tsv).

Counts shown per session:
  Clinical (N)  : unique Patient_IDs with a clinical row for that VISCODE_long
  MRI (N)       : unique Patient_IDs whose MRI scan was matched to a clinical
                  visit with that VISCODE_long under the union matching policy
                  below (Strategy A OR Strategy B)
  SNP (N)       : unique Patient_IDs with SNP genotyping for that session.
                  ADNI collects SNP only once at baseline, so this column is
                  populated at 'bl' only and reads "—" elsewhere
  All 3 (N)     : unique Patient_IDs with all three modalities present at that
                  session, treating SNP as permanent (collected once at bl,
                  applies to every subsequent visit). The intersection that
                  matters for downstream multimodal modelling.

MRI-to-clinical matching policy (UNION of two strategies)
---------------------------------------------------------
Strategy A — direct ADNI VISCODE2 join: ADNI assigned each MRI scan and each
clinical visit a VISCODE2 from the same scheduling system. For matched cases
(both sides have a row at the same VISCODE for the same subject), the
pairing is authoritative.

Strategy B — nearest-by-date within +-14 days: for scans whose ADNI VISCODE
has no matching clinical row (e.g. subject missed that clinical assessment,
or ADNI3 patients whose first MRI was at the screening visit but tagged
'bl' by ADNI), match to the closest clinical visit by date within +-14 days.

The +-14d window comes from the ADNI-2 / ADNI-3 protocol:
  "Visits must be conducted within 2 weeks before or after the target date.
   Once the visit begins, all imaging studies, biofluid collection and
   cognitive and clinical assessments must take place within 2 weeks from
   the start of the visit."  (ADNI-2 Schedule of Events; ADNI-3 §6.0)

This is read directly from `master_mri_clinical_matched_labels.csv` (built
by `mri_pipeline/03b_match_mri_to_clinical.py`). MRI scans with neither
strategy matching are reported as 'unmatched' in the footnote.

Outputs
-------
  outputs/modality_coverage_table.csv  (machine-readable counts)
  outputs/modality_coverage_table.pdf
  outputs/modality_coverage_table.png

Run
---
  python clinical_pipeline/05_modality_coverage_table.py
"""

import argparse
import re
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"


# ── CLI ────────────────────────────────────────────────────────────────────────
ap = argparse.ArgumentParser(description="ADNI cohort coverage by modality and visit.")
ap.add_argument("--aligned", action="store_true",
                help="Use the clinical-aligned master (bids_ses == VISCODE_long for matched scans). "
                     "Output goes to a parallel 'clinical_aligned_intermediates' folder.")
args = ap.parse_args()

# ── Paths ──────────────────────────────────────────────────────────────────────
REPO_ROOT      = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
MASTER_CLIN    = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\verbose\longitudinal\master_clinical_verbose.csv")
SNP_TSV        = Path(r"D:\ADNI_BIDS_project\bids\genotype\subjects_with_snp_and_mri.tsv")

if args.aligned:
    MASTER_MRI_CLIN = Path(r"D:\ADNI_BIDS_project\derivatives\mri_clinical_matched\clinical_aligned_intermediates\master_mri_clinical_matched_labels.csv")
    OUT_DIR = REPO_ROOT / "clinical_pipeline" / "outputs" / "modality_coverage" / "clinical_aligned_intermediates"
    VARIANT_LABEL = " (clinical-aligned: bids_ses = VISCODE_long for matched scans)"
else:
    MASTER_MRI_CLIN = Path(r"D:\ADNI_BIDS_project\derivatives\mri_clinical_matched\master_mri_clinical_matched_labels.csv")
    OUT_DIR = REPO_ROOT / "clinical_pipeline" / "outputs" / "modality_coverage"
    VARIANT_LABEL = ""

TOLERANCE_DAYS = 14  # bidirectional window (ADNI-2/3 within-visit tolerance), used in subtitle/footnote
MATCHED_STATUSES = ("viscode_exact", "nearest_within_14d")  # union of A and B


# ── Helpers ────────────────────────────────────────────────────────────────────
def to_pid(participant_id: str):
    """'sub-002S0413' -> '002_S_0413'."""
    m = re.match(r"sub-(\d+)S(\d+)", str(participant_id))
    return f"{m.group(1)}_S_{m.group(2)}" if m else None


def session_order(s: str) -> int:
    """Sort key for visit codes: sc < bl < m03 < m06 < ... < m240."""
    s = str(s)
    if s == "sc":
        return -1
    if s == "bl":
        return 0
    m = re.match(r"^m(\d+)$", s)
    return int(m.group(1)) if m else 9999


# ── Load sources ──────────────────────────────────────────────────────────────
print(f"Loading clinical master           : {MASTER_CLIN}")
clin = pd.read_csv(MASTER_CLIN, usecols=["Patient_ID", "VISCODE_long", "Date"])
clin["Date"] = pd.to_datetime(clin["Date"], errors="coerce")
clin = clin.dropna(subset=["Date"]).copy()

print(f"Loading MRI<->clinical match table: {MASTER_MRI_CLIN}")
mri_match = pd.read_csv(MASTER_MRI_CLIN, low_memory=False)
# bids_sub format '002S0413' -> Patient_ID '002_S_0413'
mri_match["Patient_ID"] = mri_match["bids_sub"].apply(
    lambda s: f"{str(s)[:3]}_S_{str(s)[4:]}" if str(s)[3:4] == "S" else None)

print(f"Loading SNP+MRI TSV               : {SNP_TSV}")
snp_df  = pd.read_csv(SNP_TSV, sep="\t")
snp_df["Patient_ID"] = snp_df["participant_id"].apply(to_pid)
snp_pids = set(snp_df["Patient_ID"].dropna().unique())

# Restrict everything to the SNP+MRI cohort for an apples-to-apples table
clin      = clin[clin["Patient_ID"].isin(snp_pids)].copy()
mri_match = mri_match[mri_match["Patient_ID"].isin(snp_pids)].copy()

n_total_mri = len(mri_match)
n_matched   = int(mri_match["match_status"].isin(MATCHED_STATUSES).sum())
n_unmatched = n_total_mri - n_matched

print(f"  cohort size: {len(snp_pids)} subjects (SNP+MRI overlap)")
print(f"  clinical rows in cohort        : {len(clin)}")
print(f"  MRI scans   in cohort          : {n_total_mri}")
print(f"  MRI matched (A or B)           : {n_matched} ({100*n_matched/max(n_total_mri,1):.1f}%)")
print(f"  MRI unmatched                  : {n_unmatched} ({100*n_unmatched/max(n_total_mri,1):.1f}%)")
print(f"    by status:")
print(mri_match["match_status"].value_counts().to_string())


# ── Per-session counts (union of Strategy A and B from master CSV) ───────────
# Group MRIs by bids_ses (the session label as the master records it).
# - Standard master: bids_ses = ADNI's MRI VISCODE2 (always a yearly visit
#   present in the MRI-side spine).
# - Aligned master: bids_ses = matched clinical VISCODE_long for matched
#   scans, so intermediate clinical VISCODEs (m06, m18, m30, ...) appear.
matched_only = mri_match[mri_match["match_status"].isin(MATCHED_STATUSES)].copy()

clinical_per_ses = clin.groupby("VISCODE_long")["Patient_ID"].nunique()
mri_per_ses      = matched_only.groupby("bids_ses")["Patient_ID"].nunique()

sessions = sorted(set(clin["VISCODE_long"].dropna().astype(str).unique()), key=session_order)
sessions = [s for s in sessions if s != "unknown"]

rows = []
for s in sessions:
    clin_pids_s = set(clin.loc[clin["VISCODE_long"] == s, "Patient_ID"].unique())
    mri_pids_s  = set(matched_only.loc[matched_only["bids_ses"] == s, "Patient_ID"].dropna().unique())
    snp_n       = len(snp_pids) if s == "bl" else None    # ADNI collects SNP once at baseline only
    triple_n    = len(clin_pids_s & mri_pids_s & snp_pids)  # SNP carryover
    rows.append({
        "Session":  s,
        "Clinical": int(clinical_per_ses.get(s, 0)),
        "MRI":      int(mri_per_ses.get(s, 0)),
        "SNP":      snp_n,
        "All 3":    triple_n,
    })

df_out = pd.DataFrame(rows)
OUT_DIR.mkdir(parents=True, exist_ok=True)
df_out.to_csv(OUT_DIR / "modality_coverage_table.csv", index=False)
print("\n" + df_out.to_string(index=False))


# ── Render table styled like 02b_plot_baseline_table.py ───────────────────────
# Trim trailing empty MRI rows: keep all sessions up to the last one with MRI > 0,
# plus always-show sc and bl. Late sessions with no MRI in this cohort add no signal.
last_mri_idx = max((i for i, r in enumerate(rows) if r["MRI"] > 0), default=-1)
display_rows = [r for i, r in enumerate(rows) if i <= last_mri_idx or r["Session"] in ("sc", "bl")]

HEADERS = ["Session", "Clinical N", "MRI N", "SNP N", "All-3 N"]
DATA_ROWS = [[
    r["Session"],
    f"{r['Clinical']:,}",
    f"{r['MRI']:,}" if r["MRI"] else "0",
    "-" if r["SNP"] is None else f"{r['SNP']:,}",
    f"{r['All 3']:,}" if r["All 3"] else "0",
] for r in display_rows]

n_cols = len(HEADERS)
n_rows = len(DATA_ROWS)

# Wider columns + increased font sizes; subtitle and footnote constrained to table width
COL_W      = [1.40, 1.70, 1.40, 1.40, 1.40]  # inches per column (widened)
ROW_H      = 0.32
# Title/subtitle heights scale with the variant label (longer for --aligned)
_has_variant = bool(VARIANT_LABEL)
TITLE_H1   = 0.80 if _has_variant else 0.38  # title (2 lines for --aligned)
TITLE_H2   = 0.65                             # subtitle (2 lines)
HEADER_H   = 0.36
FOOTNOTE_H = 1.30   # taller to fit wrapped footnote
TOP_PAD    = 0.12
BOT_PAD    = 0.12

LEFT      = 0.35
RIGHT_PAD = 0.35
fig_w = LEFT + sum(COL_W) + RIGHT_PAD
fig_h = (TOP_PAD + TITLE_H1 + TITLE_H2 + HEADER_H +
         ROW_H * n_rows + FOOTNOTE_H + BOT_PAD)

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.axis("off")
ax.set_xlim(0, fig_w)
ax.set_ylim(0, fig_h)

col_left = [LEFT]
for w in COL_W:
    col_left.append(col_left[-1] + w)
RIGHT  = col_left[-1]
col_cx = [(col_left[i] + col_left[i + 1]) / 2 for i in range(n_cols)]

# Row y-positions (top to bottom).
# row 0 = title line 1, 1 = title line 2, 2 = headers, 3..3+n_rows-1 = data
TOP    = fig_h - TOP_PAD
y_top  = [TOP, TOP - TITLE_H1, TOP - TITLE_H1 - TITLE_H2,
          TOP - TITLE_H1 - TITLE_H2 - HEADER_H]
for _ in range(n_rows):
    y_top.append(y_top[-1] - ROW_H)
BOTTOM = y_top[-1]
def row_mid(i): return (y_top[i] + y_top[i + 1]) / 2

# Horizontal rules
def hline(y, lw=1.0, ls="-"):
    ax.plot([LEFT, RIGHT], [y, y], color="black", linewidth=lw, linestyle=ls,
            solid_capstyle="butt", zorder=3)

hline(y_top[0], lw=1.5)               # top border
hline(y_top[2], lw=0.8)               # below title block
hline(y_top[3], lw=1.2)               # below header (midrule)
hline(BOTTOM,   lw=1.5)               # bottom border

# Dotted vertical separator after Session column
xv = col_left[1]
ax.plot([xv, xv], [BOTTOM, y_top[2]],
        color="black", linewidth=0.6, linestyle=(0, (3, 3)), zorder=2)

# Title (line 1: main) — explicit line break before variant label if present
title_base = "ADNI Cohort Coverage by Modality and Visit"
if _has_variant:
    title_text = title_base + "\n" + VARIANT_LABEL.strip()
    title_fs = 11
else:
    title_text = title_base
    title_fs = 12
ax.text((LEFT + RIGHT) / 2, row_mid(0),
        title_text,
        ha="center", va="center", fontsize=title_fs, fontweight="bold",
        color="black", linespacing=1.3)

# Title (line 2: subtitle) — each sentence on its own line
subtitle_line1 = f"SNP+MRI cohort, n={len(snp_pids)} subjects."
subtitle_line2 = (f"MRI matched on ADNI VISCODE2 OR nearest clinical "
                  f"visit within +/-{TOLERANCE_DAYS} days.")
subtitle_text = subtitle_line1 + "\n" + subtitle_line2
ax.text((LEFT + RIGHT) / 2, row_mid(1),
        subtitle_text,
        ha="center", va="center", fontsize=9.5, fontstyle="italic",
        color="black", linespacing=1.4)

# Header row
for j, h in enumerate(HEADERS):
    weight = "bold" if j == 0 else "normal"
    style  = "italic" if j > 0 else "normal"
    ax.text(col_cx[j], row_mid(2), h,
            ha="center", va="center", fontsize=10,
            fontstyle=style, fontweight=weight, color="black")

# Data rows
for r_idx, row_vals in enumerate(DATA_ROWS):
    rmid = row_mid(3 + r_idx)
    # Session label centred + bold, numbers centred
    ax.text(col_cx[0], rmid, row_vals[0],
            ha="center", va="center", fontsize=9.5, fontweight="bold", color="black")
    for j in range(1, n_cols):
        ax.text(col_cx[j], rmid, row_vals[j],
                ha="center", va="center", fontsize=9.5, color="black")

# Outer box
rect = plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT, y_top[0] - BOTTOM,
                     facecolor="none", edgecolor="black", linewidth=1.5, zorder=5)
ax.add_patch(rect)

# Footnote — wrap text to stay within the table width
foot_y = BOTTOM - 0.20
n_dropped_late = len(rows) - len(display_rows)
footnote_raw = (
    f"SNP genotyping is collected once at baseline (column SNP N reads \u2018-\u2019 for "
    f"non-baseline visits). \u2018All-3 N\u2019 treats SNP as permanent and counts subjects "
    f"with all three modalities at that visit. "
    f"MRI N counts subjects whose scan matched a clinical visit by ADNI VISCODE2 "
    f"OR nearest clinical visit within +/-{TOLERANCE_DAYS} days. "
    f"Unmatched MRI scans: {n_unmatched} of "
    f"{n_total_mri} ({100*n_unmatched/max(n_total_mri,1):.1f}%), excluded from the "
    f"MRI column. Tolerance derived from ADNI-2 / ADNI-3 protocols: all imaging and "
    f"clinical components of one visit must complete within 2 weeks of the visit "
    f"start (ADNI-2 Schedule of Events; ADNI-3 Sec 6.0)."
    + (f"  {n_dropped_late} late session(s) with no MRI in cohort omitted." if n_dropped_late else "")
)
# Wrap footnote to span the full table width (~17 chars/inch at fontsize 8)
foot_wrap_chars = int((RIGHT - LEFT) * 17)
footnote_wrapped = "\n".join(
    textwrap.wrap(footnote_raw, width=foot_wrap_chars))
ax.text(LEFT, foot_y,
        footnote_wrapped,
        ha="left", va="top", fontsize=8, color="black")

plt.tight_layout(pad=0.1)
fig.savefig(OUT_DIR / "modality_coverage_table.png", bbox_inches="tight", dpi=300)
fig.savefig(OUT_DIR / "modality_coverage_table.pdf", bbox_inches="tight", dpi=300)
print(f"\nSaved CSV -> {OUT_DIR / 'modality_coverage_table.csv'}")
print(f"Saved PNG -> {OUT_DIR / 'modality_coverage_table.png'}")
print(f"Saved PDF -> {OUT_DIR / 'modality_coverage_table.pdf'}")
plt.close()
