"""
05b_modality_coverage_table_viscode2.py
========================================
Cross-modality coverage audit using VISCODE2 alignment.

Unlike 05 (which uses VISCODE_long), this script:
  - Counts clinical visits by VISCODE_2
  - Counts MRI scans by adni_viscode (which == VISCODE_2)
  - Uses the VISCODE2-aligned master from 03c

This shows exactly how many MRI bids_ses match clinical VISCODE_2 codes
for the SNP+MRI cohort.

Outputs (to outputs/modality_coverage/VISCODE_2_aligned/)
-------
  modality_coverage_viscode2.csv
  modality_coverage_viscode2.png
  modality_coverage_viscode2.pdf

Run
---
  python clinical_pipeline/05b_modality_coverage_table_viscode2.py
"""

import re
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"


# ── Paths ──────────────────────────────────────────────────────────────────────
REPO_ROOT      = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
MASTER_CLIN    = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified"
                      r"\tabular\longitudinal\master_clinical_tabular.csv")
MASTER_MRI_V2  = Path(r"D:\ADNI_BIDS_project\derivatives\mri_clinical_matched"
                      r"\viscode_2_aligned\master_mri_clinical_matched_viscode2.csv")
SNP_TSV        = Path(r"D:\ADNI_BIDS_project\bids\genotype\subjects_with_snp_and_mri.tsv")

OUT_DIR = REPO_ROOT / "clinical_pipeline" / "outputs" / "modality_coverage" / "viscode_2_aligned"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MATCHED_STATUSES = ("viscode2_exact", "nearest_within_14d")


# ── Helpers ────────────────────────────────────────────────────────────────────
def to_pid(participant_id: str):
    """'sub-002S0413' -> '002_S_0413'."""
    m = re.match(r"sub-(\d+)S(\d+)", str(participant_id))
    return f"{m.group(1)}_S_{m.group(2)}" if m else None


def session_order(s: str) -> int:
    s = str(s)
    if s == "sc":  return -1
    if s == "bl":  return 0
    m = re.match(r"^m(\d+)$", s)
    return int(m.group(1)) if m else 9999


def pid_from_bids(bids_sub: str) -> str:
    """002S0413 -> 002_S_0413."""
    s = str(bids_sub).strip()
    if len(s) >= 5 and s[3] == "S":
        return f"{s[:3]}_S_{s[4:]}"
    return s


# ── Load sources ──────────────────────────────────────────────────────────────
print(f"Loading clinical master     : {MASTER_CLIN}")
clin = pd.read_csv(MASTER_CLIN,
                    usecols=["Patient_ID", "VISCODE_2", "VISCODE_long", "Date"],
                    low_memory=False)
clin["Date"] = pd.to_datetime(clin["Date"], errors="coerce")
clin = clin.dropna(subset=["Date"]).copy()

print(f"Loading MRI VISCODE2 master : {MASTER_MRI_V2}")
mri = pd.read_csv(MASTER_MRI_V2, low_memory=False)
mri["Patient_ID"] = mri["bids_sub"].apply(pid_from_bids)

print(f"Loading SNP+MRI cohort      : {SNP_TSV}")
snp_df = pd.read_csv(SNP_TSV, sep="\t")
snp_df["Patient_ID"] = snp_df["participant_id"].apply(to_pid)
snp_pids = set(snp_df["Patient_ID"].dropna().unique())

# Filter to cohort
clin = clin[clin["Patient_ID"].isin(snp_pids)].copy()
mri  = mri[mri["Patient_ID"].isin(snp_pids)].copy()
matched_mri = mri[mri["match_status"].isin(MATCHED_STATUSES)].copy()

n_total_mri = len(mri)
n_matched   = len(matched_mri)
n_unmatched = n_total_mri - n_matched

print(f"  cohort: {len(snp_pids)} subjects")
print(f"  clinical rows: {len(clin)}")
print(f"  MRI scans: {n_total_mri} (matched: {n_matched}, unmatched: {n_unmatched})")
print(f"  match status breakdown:")
print(mri["match_status"].value_counts().to_string())
print()


# ── Build the visit list (union of clinical VISCODE_2 and MRI adni_viscode) ──
clinical_v2s = set(clin["VISCODE_2"].dropna().astype(str).unique())
mri_v2s      = set(mri["adni_viscode"].dropna().astype(str).unique())
all_sessions = sorted(clinical_v2s | mri_v2s, key=session_order)
all_sessions = [s for s in all_sessions if s not in ("unknown", "nan", "")]


# ── Per-session counts ────────────────────────────────────────────────────────
clinical_per_ses = clin.groupby("VISCODE_2")["Patient_ID"].nunique()
mri_all_per_ses  = mri.groupby("adni_viscode")["Patient_ID"].nunique()
mri_matched_per_ses = matched_mri.groupby("adni_viscode")["Patient_ID"].nunique()

# SNP is baseline-only
snp_bl_pids = snp_pids  # collected once at baseline, applies to all visits

result_rows = []
for ses in all_sessions:
    clin_n   = int(clinical_per_ses.get(ses, 0))
    mri_all  = int(mri_all_per_ses.get(ses, 0))
    mri_n    = int(mri_matched_per_ses.get(ses, 0))

    # SNP: present at bl only, but conceptually available at all visits
    if ses == "bl":
        snp_n = len(snp_bl_pids)
    else:
        snp_n = 0  # show 0 to indicate SNP not re-collected

    # All-3 intersection: subjects with clinical + matched MRI + SNP at this session
    clin_pids = set(clin[clin["VISCODE_2"] == ses]["Patient_ID"].unique())
    mri_pids  = set(matched_mri[matched_mri["adni_viscode"] == ses]["Patient_ID"].unique())
    all3_pids = clin_pids & mri_pids & snp_pids
    all3_n    = len(all3_pids)

    result_rows.append({
        "Session_VISCODE2": ses,
        "Clinical_N": clin_n,
        "MRI_total_N": mri_all,
        "MRI_matched_N": mri_n,
        "SNP_N": snp_n if ses == "bl" else "—",
        "All3_N": all3_n,
        "MRI_match_pct": f"{100 * mri_n / mri_all:.1f}%" if mri_all > 0 else "—",
    })

result_df = pd.DataFrame(result_rows)
result_df.to_csv(OUT_DIR / "modality_coverage_viscode2.csv", index=False)
print(f"Saved CSV -> {OUT_DIR / 'modality_coverage_viscode2.csv'}")
print(result_df.to_string(index=False))
print()


# ── Render publication-quality table ──────────────────────────────────────────
HEADERS = ["Session\n(VISCODE2)", "Clinical\nN", "MRI\n(total)", "MRI\n(matched)",
           "Match\n%", "SNP\nN", "All 3\nN"]
data_cols = ["Session_VISCODE2", "Clinical_N", "MRI_total_N", "MRI_matched_N",
             "MRI_match_pct", "SNP_N", "All3_N"]

# Trim to sessions with at least some data
display_df = result_df[(result_df["Clinical_N"] > 0) | (result_df["MRI_total_N"] > 0)].copy()

# Add TOTAL row
totals = {
    "Session_VISCODE2": "TOTAL",
    "Clinical_N": int(result_df["Clinical_N"].sum()),
    "MRI_total_N": int(result_df["MRI_total_N"].sum()),
    "MRI_matched_N": int(result_df["MRI_matched_N"].sum()),
    "MRI_match_pct": f"{100 * result_df['MRI_matched_N'].sum() / max(result_df['MRI_total_N'].sum(), 1):.1f}%",
    "SNP_N": len(snp_pids),
    "All3_N": "—",
}
display_df = pd.concat([display_df, pd.DataFrame([totals])], ignore_index=True)

table_rows = display_df[data_cols].astype(str).values.tolist()
n_rows = len(table_rows)
n_cols = len(HEADERS)

COL_W = [1.00, 0.85, 0.85, 0.85, 0.80, 0.70, 0.80]
ROW_H = 0.32
HDR_H = 0.55

fig_w = 0.30 + sum(COL_W) + 0.30
fig_h = 0.80 + HDR_H + ROW_H * n_rows + 0.60

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.axis("off")
ax.set_xlim(0, fig_w)
ax.set_ylim(0, fig_h)

LEFT = 0.30
col_left = [LEFT]
for w in COL_W:
    col_left.append(col_left[-1] + w)
col_cx = [(col_left[i] + col_left[i + 1]) / 2 for i in range(n_cols)]
RIGHT = col_left[-1]

TOP = fig_h - 0.15
y_cursor = TOP

def hline(y, lw=1.0, ls="-"):
    ax.plot([LEFT, RIGHT], [y, y], color="black", linewidth=lw, linestyle=ls,
            solid_capstyle="butt", zorder=3)

# Title
hline(y_cursor, lw=1.5)
y_cursor -= 0.40
ax.text((LEFT + RIGHT) / 2, y_cursor + 0.20,
        "Cross-Modality Coverage by ADNI Session (VISCODE2 alignment)",
        ha="center", va="center", fontsize=10, fontweight="bold")

# Subtitle
y_cursor -= 0.30
ax.text((LEFT + RIGHT) / 2, y_cursor + 0.15,
        f"SNP+MRI cohort: {len(snp_pids)} subjects  |  "
        f"MRI matched: {n_matched}/{n_total_mri} "
        f"({100*n_matched/max(n_total_mri,1):.1f}%)  |  "
        f"Unmatched: {n_unmatched}",
        ha="center", va="center", fontsize=8, fontstyle="italic", color="#333333")

hline(y_cursor, lw=1.2)

# Column headers
hdr_top = y_cursor
y_cursor -= HDR_H
for j, h in enumerate(HEADERS):
    ax.text(col_cx[j], (hdr_top + y_cursor) / 2, h,
            ha="center", va="center", fontsize=8, fontstyle="italic")

hline(y_cursor, lw=0.8)

# Data rows
for i, row in enumerate(table_rows):
    y_top = y_cursor
    y_cursor -= ROW_H
    y_mid = (y_top + y_cursor) / 2
    is_total = (row[0] == "TOTAL")

    if is_total:
        hline(y_top, lw=1.0)

    for j, val in enumerate(row):
        weight = "bold" if j == 0 or is_total else "normal"
        ax.text(col_cx[j], y_mid, val, ha="center", va="center",
                fontsize=8, fontweight=weight)

hline(y_cursor, lw=1.5)

# Vertical separators
for x_idx in [1, 4, 5, 6]:
    x = col_left[x_idx]
    ax.plot([x, x], [y_cursor, hdr_top],
            color="black", linewidth=0.5, linestyle=(0, (3, 3)), zorder=2)

# Outer box
rect = plt.Rectangle((LEFT, y_cursor), RIGHT - LEFT, TOP - y_cursor,
                      facecolor="none", edgecolor="black", linewidth=1.5, zorder=5)
ax.add_patch(rect)

# Footnote
foot_y = y_cursor - 0.12
footnote = (
    f"MRI (matched) = scans with VISCODE2-exact or ±14d nearest clinical match. "
    f"SNP collected once at baseline. "
    f"All 3 = subjects with clinical + matched MRI + SNP at that session."
)
foot_wrap = int((RIGHT - LEFT) * 18)
ax.text(LEFT, foot_y, "\n".join(textwrap.wrap(footnote, width=foot_wrap)),
        ha="left", va="top", fontsize=7, color="#444444")

plt.tight_layout(pad=0.1)
fig.savefig(OUT_DIR / "modality_coverage_viscode2.png", bbox_inches="tight", dpi=300)
fig.savefig(OUT_DIR / "modality_coverage_viscode2.pdf", bbox_inches="tight", dpi=300)
plt.close()
print(f"Saved PNG -> {OUT_DIR / 'modality_coverage_viscode2.png'}")
print(f"Saved PDF -> {OUT_DIR / 'modality_coverage_viscode2.pdf'}")
print("\nDone.")
