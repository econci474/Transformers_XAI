"""
prepare_bmfm_labels_from_external_gwas.py
==========================================
Generate BMFM-DNA-SNP fine-tuning labels (label=1/0) by intersecting
published AD GWAS summary statistics with the ADNI BIM file.

Sources
-------
  Bellenguez 2022 (Nat Genet, PMID 35379992)
    - gwas-association-downloaded_*.tsv   (GWAS Catalog format, 90 GW-sig SNPs)

  Wightman 2021 (Nat Genet, PMID 34385711)
    Two mutually exclusive summary-stat files are available:

    [DEFAULT]  ExcludingUKBand23andME  — excludes both UK Biobank and 23andMe
      PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis (1).txt.gz

    [--with-ukb]  ExcludingOnly23andMe  — includes UK Biobank, excludes only 23andMe
      PGCALZ2sumstatsExcluding23andMe.txt.gz

Labelling strategy
------------------
  label=1 : SNP is genome-wide significant (p < 5e-8) in at least one
             external study AND is present in the ADNI genotyped BIM.
  label=0 : SNP is in ADNI BIM, NOT in any external significant list,
             and has p-value >= 0.5 in our internal GWAS (if available).
             If no internal GWAS is available, all remaining ADNI SNPs
             are used as nulls.

Outputs (default)
-----------------
  D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/external_gwas_labels.tsv
  D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/external_gwas_summary.txt

Outputs (--with-ukb)
--------------------
  D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/external_gwas_labels_with_ukb.tsv
  D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/external_gwas_summary_with_ukb.txt

Usage
-----
  python prepare_bmfm_labels_from_external_gwas.py
  python prepare_bmfm_labels_from_external_gwas.py --with-ukb
  python prepare_bmfm_labels_from_external_gwas.py --null-p-min 0.5
  python prepare_bmfm_labels_from_external_gwas.py --no-internal-gwas
"""

import argparse
import gzip
import io
import math
import pathlib
import re
import sys
from datetime import datetime

import pandas as pd

if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR   = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220")
GWAS_DIR   = BASE_DIR / "GWAS"
BIM_FILE   = BASE_DIR / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr.bim"
OUT_DIR    = BASE_DIR / "bmfm_inputs"

BELLENGUEZ_TSV = next(
    (GWAS_DIR / "Bellenguez_2022" / f
     for f in ["gwas-association-downloaded_2026-04-28-pubmedId_35379992.tsv"]
     if (GWAS_DIR / "Bellenguez_2022" / f).exists()),
    None
)

# Pre-lifted GRCh38 TSVs — produced by liftover_hg19_to_hg38_wightman_2021.py
WIGHTMAN_EXCL_UKB_TSV = GWAS_DIR / "Wightman_2021" / "wightman_without_ukb_GRCh38.tsv"
WIGHTMAN_INCL_UKB_TSV = GWAS_DIR / "Wightman_2021" / "wightman_with_ukb_GRCh38.tsv"
# Resolved after CLI parsing (see below)

INTERNAL_GWAS = BASE_DIR / "results" / "gwas" / "gwas_CN_vs_AD_GRCh38.assoc.logistic"

# ── CLI ───────────────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser(
    description="Build BMFM-DNA-SNP labels from external AD GWAS summary stats.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
_p.add_argument("--bim",    default=str(BIM_FILE))
_p.add_argument("--out-dir", default=str(OUT_DIR), dest="out_dir")
_p.add_argument("--gw-thresh", type=float, default=5e-8, dest="gw_thresh",
                help="Genome-wide significance threshold (default: 5e-8)")
_p.add_argument("--null-p-min", type=float, default=0.5, dest="null_p_min",
                help="Min internal-GWAS p for a SNP to be labelled null (default: 0.5)")
_p.add_argument("--no-internal-gwas", action="store_true", dest="no_internal",
                help="Skip internal GWAS filtering; use all non-significant ADNI SNPs as nulls")
_p.add_argument("--with-ukb", action="store_true", dest="with_ukb",
                help=(
                    "Use the Wightman 2021 file that INCLUDES UK Biobank "
                    "(excludes only 23andMe). "
                    "Default uses the file excluding BOTH UKB and 23andMe. "
                    "Output files are suffixed '_with_ukb'."
                ))
args = _p.parse_args()

OUT_DIR  = pathlib.Path(args.out_dir)
BIM_FILE = pathlib.Path(args.bim)
OUT_DIR.mkdir(parents=True, exist_ok=True)

GW_THRESH  = args.gw_thresh

# ── Wightman file selection (mutually exclusive via --with-ukb) ───────────────
WIGHTMAN_TSV = WIGHTMAN_INCL_UKB_TSV if args.with_ukb else WIGHTMAN_EXCL_UKB_TSV

# Output file suffix
_out_suffix = "_with_ukb" if args.with_ukb else ""

if args.with_ukb:
    print("[INFO] --with-ukb selected: using Wightman 2021 file INCLUDING UK Biobank "
          "(excluding only 23andMe).")
else:
    print("[INFO] Default: using Wightman 2021 file EXCLUDING both UKB and 23andMe.")

# ── Load ADNI BIM ─────────────────────────────────────────────────────────────
print(f"Loading ADNI BIM: {BIM_FILE}")
bim = pd.read_csv(BIM_FILE, sep="\t", header=None,
                  names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
                  dtype={"CHR": str, "SNP": str, "BP": int, "A1": str, "A2": str})
adni_rsids = set(bim["SNP"])
print(f"  {len(adni_rsids):,} SNPs in ADNI BIM")

# Positional lookup: (str(CHR), int(BP)) -> bim row  [used for Wightman matching]
bim_by_pos = {(str(r.CHR), int(r.BP)): r for r in bim.itertuples(index=False)}

def _match_ea(ea, bim_row):
    """
    Strictly match the GWAS effect allele (EA) against BIM A1 / A2.
    Returns:
      allele_col : "A1" | "A2" | None   (which BIM allele the EA maps to)
      flip       : bool                 (True → negate z-score to orient to BIM A1)
    If EA maps to A1 (the ALT/counted allele) → no flip needed.
    If EA maps to A2 (the REF allele after refcorr) → flip = True.
    """
    if not ea:
        return None, False
    ea_u = ea.upper()
    a1_u = (str(bim_row.A1) if bim_row.A1 else "").upper()
    a2_u = (str(bim_row.A2) if bim_row.A2 else "").upper()
    if ea_u == a1_u and ea_u not in {"", "0", "NA", "."}:
        return "A1", False      # EA is ALT — effect correctly signed for A1
    if ea_u == a2_u and ea_u not in {"", "0", "NA", "."}:
        return "A2", True       # EA is REF — must negate z to express effect per A1 copy
    return None, False

def _z_from_or_ci(or_val, ci_text):
    """
    Calculate z-score from OR and 95% CI text (e.g. '[1.15-1.32]').
    Returns float or None.
    """
    try:
        or_val = float(or_val)
        if or_val <= 0:
            return None
        m = re.search(r"\[([\d.]+)[^\d]+([\d.]+)\]", str(ci_text or ""))
        if not m:
            return None
        ci_lo, ci_hi = float(m.group(1)), float(m.group(2))
        if ci_lo <= 0 or ci_hi <= 0 or ci_lo >= ci_hi:
            return None
        se = (math.log(ci_hi) - math.log(ci_lo)) / (2 * 1.96)
        return math.log(or_val) / se if se > 0 else None
    except (TypeError, ValueError, ZeroDivisionError):
        return None

# ── Parse Bellenguez 2022 (GWAS Catalog TSV) ──────────────────────────────────
bellenguez_sig  = set()          # rsIDs that are GW-sig
bellenguez_ea   = {}             # rsID → risk allele (EA after '-' in catalog field)
bellenguez_z_raw = {}            # rsID → raw z-score (unsigned; sign determined by EA vs BIM)
if BELLENGUEZ_TSV and BELLENGUEZ_TSV.exists():
    print(f"\nLoading Bellenguez 2022: {BELLENGUEZ_TSV.name}")
    b = pd.read_csv(BELLENGUEZ_TSV, sep="\t")
    b["P-VALUE"] = pd.to_numeric(b["P-VALUE"], errors="coerce")
    b_sig = b[b["P-VALUE"] < GW_THRESH].copy()
    for _, row in b_sig.iterrows():
        rsid = str(row.get("SNPS", "") or "").strip()
        if not rsid:
            continue
        bellenguez_sig.add(rsid)
        # Extract EA from e.g. "rs141749679-C"
        risk_field = str(row.get("STRONGEST SNP-RISK ALLELE", "") or "")
        if "-" in risk_field:
            ea = risk_field.split("-")[-1].strip().upper()
            if ea:
                bellenguez_ea[rsid] = ea
        # Calculate z from OR and 95% CI
        z = _z_from_or_ci(row.get("OR or BETA"), row.get("95% CI (TEXT)"))
        if z is not None:
            bellenguez_z_raw[rsid] = z
    print(f"  {len(b_sig)} genome-wide significant SNPs (p < {GW_THRESH:.0e})")
    print(f"  Risk alleles parsed : {len(bellenguez_ea):,} / {len(bellenguez_sig):,}")
    print(f"  z-scores calculated : {len(bellenguez_z_raw):,} / {len(bellenguez_sig):,}")
else:
    print("\n[WARNING] Bellenguez 2022 TSV not found — skipping.")

# ── Load Wightman 2021 (pre-lifted GRCh38 TSV, positional matching) ─────────
# wightman_matched: bim_rsid → dict with keys:
#   match_method, ea_gwas, allele_col, flip, z_gwas, z_adni
wightman_matched = {}
wightman_n_sig   = 0
wightman_gwas_rsids = set()   # all GW-sig rsIDs from Wightman TSV (for GWAS-level overlap)
if WIGHTMAN_TSV.exists():
    print(f"\nLoading Wightman 2021 GRCh38 TSV: {WIGHTMAN_TSV.name}")
    w = pd.read_csv(WIGHTMAN_TSV, sep="\t")
    w = w[w["lifted"] == True].copy()
    w["CHR_38"] = w["CHR_38"].astype(str)
    w["BP_38"]  = pd.to_numeric(w["BP_38"], errors="coerce").dropna().astype(int)
    w = w.dropna(subset=["CHR_38", "BP_38"])
    wightman_n_sig = len(w)
    # Collect Wightman rsIDs for GWAS-level intersection with Bellenguez
    for col in ["rsID", "rsID_canonical"]:
        if col in w.columns:
            wightman_gwas_rsids |= set(w[col].dropna())
    has_rsid_col = "rsID" in w.columns
    print(f"  {wightman_n_sig:,} GW-significant lifted SNPs  "
          f"({'rsID column present' if has_rsid_col else 'no rsID column — positional only'})")

    n_rsid_match, n_pos_match, n_allele_fail = 0, 0, 0

    for _, row in w.iterrows():
        ea = str(row.get("effect_allele") or "").upper().strip()
        z_gwas = row.get("z")
        try:
            z_gwas = float(z_gwas) if z_gwas is not None and str(z_gwas) not in ("", "nan", "NA") else None
        except (TypeError, ValueError):
            z_gwas = None
        # Fallback: compute z from beta/SE
        # (without-UKB Wightman has no 'z' column; with-UKB file has 'z' directly)
        if z_gwas is None:
            try:
                z_gwas = float(row.get("beta")) / float(row.get("standard_error"))
            except (TypeError, ValueError, ZeroDivisionError):
                pass

        def _store_wightman(rsid, bim_row, method):
            """Strict EA match to BIM; compute harmonized z; store in wightman_matched."""
            allele_col, flip = _match_ea(ea, bim_row)
            if allele_col is None and ea:
                return False   # EA not in BIM alleles — reject
            z_adni = (-z_gwas if flip else z_gwas) if z_gwas is not None else None
            if rsid not in wightman_matched:
                wightman_matched[rsid] = {
                    "match_method": method,
                    "ea_gwas"     : ea,
                    "allele_col"  : allele_col,
                    "flip"        : flip,
                    "z_gwas"      : z_gwas,
                    "z_adni"      : z_adni,
                }
            return True

        # ── Tier 1: rsID or rsID_canonical match ────────────────────────────
        # Wightman rsIDs are deprecated; ADNI BIM has canonical IDs.
        # Try rsID first, then rsID_canonical (from RsMergeArch lookup in liftover).
        if has_rsid_col:
            matched_t1 = False
            for rsid_key in ["rsID_canonical", "rsID"]:
                rsid_val = row.get(rsid_key)
                if pd.notna(rsid_val) and str(rsid_val) in adni_rsids:
                    rsid = str(rsid_val)
                    bim_rows = bim[bim["SNP"] == rsid]
                    if not bim_rows.empty:
                        if _store_wightman(rsid, bim_rows.iloc[0],
                                           f"rsID_canonical+EA" if rsid_key == "rsID_canonical" else "rsID+EA"):
                            n_rsid_match += 1
                        matched_t1 = True
                    break
            if matched_t1:
                continue

        # ── Tier 2: hg38 positional match ────────────────────────────────────
        key = (str(row["CHR_38"]), int(row["BP_38"]))
        bim_row = bim_by_pos.get(key)
        if bim_row is None:
            continue
        if _store_wightman(bim_row.SNP, bim_row, "positional+EA"):
            n_pos_match += 1
        else:
            n_allele_fail += 1

    print(f"  rsID matches              : {n_rsid_match:,}")
    print(f"  Positional matches (hg38) : {n_pos_match:,}")
    if n_allele_fail:
        print(f"  Rejected (EA not in BIM)  : {n_allele_fail:,}")
    print(f"  Total matched to ADNI     : {len(wightman_matched):,}")

else:
    print(f"\n[WARNING] Wightman 2021 GRCh38 TSV not found: {WIGHTMAN_TSV}")


# ── Merge matched SNPs from all sources ──────────────────────────────────────
# Bellenguez: rsID intersection + EA-to-BIM allele check for z-sign
bellenguez_matched_rsids = bellenguez_sig & adni_rsids
bellenguez_matched = {}   # rsid → {ea_gwas, z_adni, allele_col, flip}
for rsid in bellenguez_matched_rsids:
    ea = bellenguez_ea.get(rsid, "")
    z_raw = bellenguez_z_raw.get(rsid)
    bim_rows = bim[bim["SNP"] == rsid]
    allele_col, flip = (None, False)
    if not bim_rows.empty and ea:
        allele_col, flip = _match_ea(ea, bim_rows.iloc[0])
    z_adni = None
    if z_raw is not None:
        z_adni = -z_raw if flip else z_raw
    bellenguez_matched[rsid] = {
        "ea_gwas"   : ea,
        "allele_col": allele_col,
        "flip"      : flip,
        "z_adni"    : z_adni,
    }

in_adni_sig = set(bellenguez_matched) | set(wightman_matched)
overlap_both = set(bellenguez_matched) & set(wightman_matched)

print(f"\nMatched to ADNI BIM (label=1 candidates):")
print(f"  Bellenguez 2022 (rsID + EA-strict)  : {len(bellenguez_matched):,}  of {len(bellenguez_sig):,} GW-sig")
print(f"  Wightman 2021 (rsID/pos + EA-strict): {len(wightman_matched):,}  of {wightman_n_sig:,} GW-sig lifted")
print(f"  Union label=1                        : {len(in_adni_sig):,}")

gwas_level_overlap = wightman_gwas_rsids & bellenguez_sig
print(f"\nIntersection between Bellenguez 2022 and Wightman 2021:")
print(f"  GW-sig in both studies (GWAS level)  : {len(gwas_level_overlap):,}  "
      f"(Wightman rsIDs matched to Bellenguez GW-sig)")
print(f"  Both matched to ADNI BIM             : {len(overlap_both):,}")

# ── Load internal GWAS for null selection (optional) ─────────────────────────
null_p_filter = None
if not args.no_internal and INTERNAL_GWAS.exists():
    print(f"\nLoading internal GWAS for null selection: {INTERNAL_GWAS.name}")
    gwas_df = pd.read_csv(INTERNAL_GWAS, sep=r"\s+")
    gwas_df = gwas_df[gwas_df["TEST"] == "ADD"].copy()
    gwas_df["P"] = pd.to_numeric(gwas_df["P"], errors="coerce")
    null_p_filter = set(
        gwas_df[gwas_df["P"] >= args.null_p_min]["SNP"].dropna()
    )
    print(f"  {len(null_p_filter):,} SNPs with p >= {args.null_p_min} "
          f"(usable as nulls)")
elif not args.no_internal:
    print(f"\n[INFO] Internal GWAS not found at {INTERNAL_GWAS}")
    print("  Using all non-significant ADNI SNPs as null set.")

# ── Build label dataframe ─────────────────────────────────────────────────────
bim_sig = bim[bim["SNP"].isin(in_adni_sig)].copy()
bim_sig["label"]       = 1
bim_sig["label_basis"] = "external_GW_sig"
bim_sig["source"] = bim_sig["SNP"].apply(
    lambda r: "+".join(filter(None, [
        "Bellenguez2022" if r in bellenguez_matched else "",
        "Wightman2021"   if r in wightman_matched   else "",
    ]))
)
# ── Source preference + dual z columns ─────────────────────────────────────
# Priority for primary z/EA:
#   1. Both EA-matched → Bellenguez primary, Wightman secondary
#   2. Only one EA-matched → that source is primary, no secondary
#   3. Neither EA-matched → Bellenguez (rsID) primary if available, else Wightman
def _resolve(rsid):
    bell  = bellenguez_matched.get(rsid) or {}
    wight = wightman_matched.get(rsid)   or {}
    b_ok  = bool(bell.get("allele_col"))
    w_ok  = bool(wight.get("allele_col"))
    if b_ok and w_ok:
        return bell, wight, "Bellenguez2022_logOR", "Wightman2021"
    elif b_ok:
        return bell, None, "Bellenguez2022_logOR", None
    elif w_ok:
        return wight, None, "Wightman2021", None
    else:   # no EA match in either — use available data, flag for review
        if bell:
            return bell, (wight or None), "Bellenguez2022_logOR", ("Wightman2021" if wight else None)
        return wight, None, "Wightman2021", None

_res = {r: _resolve(r) for r in in_adni_sig}

bim_sig["REF"] = bim_sig["A2"]
bim_sig["EA"]  = bim_sig["SNP"].map(lambda r: _res[r][0].get("ea_gwas") or pd.NA)
bim_sig["OA"]  = bim_sig.apply(
    lambda row: (row["A2"] if str(row["EA"]).upper() == str(row["A1"]).upper()
                 else row["A1"]
                 if str(row["EA"]).upper() == str(row["A2"]).upper()
                 else pd.NA), axis=1
)
bim_sig["z_score"]             = bim_sig["SNP"].map(lambda r: _res[r][0].get("z_adni"))
bim_sig["z_source"]            = bim_sig["SNP"].map(lambda r: _res[r][2])
bim_sig["allele_flip"]         = bim_sig["SNP"].map(lambda r: _res[r][0].get("flip"))
bim_sig["z_score_additional"]  = bim_sig["SNP"].map(lambda r: (_res[r][1] or {}).get("z_adni"))
bim_sig["z_source_additional"] = bim_sig["SNP"].map(lambda r: _res[r][3])

n_z_primary = bim_sig["z_score"].notna().sum()
n_z_addl    = bim_sig["z_score_additional"].notna().sum()
print(f"  z_score (primary)  : {n_z_primary:,} / {len(bim_sig):,} label=1 SNPs have a z-score")
print(f"  z_score (secondary): {n_z_addl:,} SNPs have a second-source z-score")

# Null SNPs: not in external significant list
if null_p_filter is not None:
    null_candidates = (adni_rsids - in_adni_sig) & null_p_filter
else:
    null_candidates = adni_rsids - in_adni_sig

bim_null = bim[bim["SNP"].isin(null_candidates)].copy()
bim_null["label"]              = 0
bim_null["label_basis"]        = "not_GW_sig"
bim_null["source"]             = "null"
bim_null["REF"]                = bim_null["A2"]
bim_null["EA"]                 = pd.NA
bim_null["OA"]                 = pd.NA
bim_null["z_score"]            = pd.NA
bim_null["z_source"]           = pd.NA
bim_null["allele_flip"]        = pd.NA
bim_null["z_score_additional"] = pd.NA
bim_null["z_source_additional"]= pd.NA

# Combine
labels = pd.concat([bim_sig, bim_null], ignore_index=True)
labels_out = OUT_DIR / f"external_gwas_labels{_out_suffix}.tsv"
labels[[
    "SNP", "CHR", "BP",
    "REF", "EA", "OA",
    "label", "label_basis", "source",
    "z_score", "z_source", "allele_flip",
    "z_score_additional", "z_source_additional",
]].to_csv(labels_out, sep="\t", index=False)

# ── Summary report ────────────────────────────────────────────────────────────
SEP = "=" * 65
report_lines = [
    SEP,
    "  BMFM Label Generation — External AD GWAS",
    f"  Generated : {datetime.now().strftime('%Y-%m-%d %H:%M')}",
    SEP,
    "",
    f"  ADNI BIM SNPs total           : {len(adni_rsids):>10,}",
    f"  GW-sig Bellenguez (GWAS Cat.) : {len(bellenguez_sig):>10,}  → {len(bellenguez_matched):,} in ADNI (rsID + EA-strict)",
    f"  GW-sig Wightman (lifted TSV)  : {wightman_n_sig:>10,}  → {len(wightman_matched):,} in ADNI (rsID/pos + EA-strict)",
    f"  Replicated in both studies    : {len(overlap_both):>10,}  (ADNI-matched)",
    f"  GWAS-level overlap (rsID)     : {len(gwas_level_overlap):>10,}  (Wightman rsIDs ∩ Bellenguez)",
    f"  Union matched in ADNI (label=1): {len(in_adni_sig):>9,}",
    "",
    f"  Overlap with ADNI BIM (label=1) : {len(bim_sig):>8,}",
    f"  Null SNPs in ADNI BIM (label=0) : {len(bim_null):>8,}",
    f"  Total labelled                  : {len(labels):>8,}",
    "",
    "  Top label=1 SNPs (by source):",
]
for _, row in bim_sig.head(10).iterrows():
    report_lines.append(
        f"    {row['SNP']:<20}  chr{row['CHR']}:{row['BP']:,}  [{row['source']}]"
    )
report_lines += [
    "",
    f"  Output : {labels_out}",
    SEP,
]
report_text = "\n".join(report_lines)
print("\n" + report_text)

summary_path = OUT_DIR / f"external_gwas_summary{_out_suffix}.txt"
summary_path.write_text(report_text, encoding="utf-8")
print(f"Summary written: {summary_path}")
print(f"Labels written : {labels_out}")
