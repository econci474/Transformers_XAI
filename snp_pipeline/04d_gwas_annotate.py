"""
gwas_annotate.py – Annotate the top GWAS hits from ADNI CN vs AD analysis.

For each significant SNP (p < 1e-5) this script looks up:
  1. Nearest gene and functional annotation  → via MyVariant.info (free REST API)
  2. Prior GWAS evidence                     → via GWAS Catalog REST API
  3. Allele frequency (gnomAD EUR population) → via MyVariant.info

Inputs
------
  D:/ADNI_SNP_Omni2.5M_20140220/results/gwas_no_ld/gwas_CN_vs_AD_no_pruning.assoc.logistic

Outputs
-------
  D:/ADNI_SNP_Omni2.5M_20140220/results/gwas_no_ld/annotation/top_snps_annotated.tsv
  D:/ADNI_SNP_Omni2.5M_20140220/results/gwas_no_ld/annotation/top_snps_annotated.html

Usage
-----
  pip install pandas requests
  python gwas_annotate.py

Notes
-----
  - Requires internet access to query REST APIs.
  - Rate limit: MyVariant.info allows ~1000 free queries/day; we have ~6 SNPs.
  - GWAS Catalog queries association evidence for each rsID.
  - kgp* IDs (Illumina internal) cannot be resolved by rsID; the script will
    attempt a lift-over by chr:pos query instead.
"""

import json
import pathlib
import time
import warnings

import pandas as pd
import requests

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
GWAS_FILE = pathlib.Path(
    "D:/ADNI_SNP_Omni2.5M_20140220/results/gwas_no_ld/"
    "gwas_CN_vs_AD_no_pruning.assoc.logistic"
)
OUT_DIR = GWAS_FILE.parent / "annotation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

P_THRESH = 1e-5   # significance threshold

# ─────────────────────────────────────────────────────────────────────────────
# 1. Load top hits
# ─────────────────────────────────────────────────────────────────────────────
print("Loading GWAS results …")
df = pd.read_csv(GWAS_FILE, sep=r"\s+")
df.columns = df.columns.str.strip()
df = df[df["TEST"] == "ADD"].copy()
df["P"] = pd.to_numeric(df["P"], errors="coerce")
df = df.dropna(subset=["P"])
df = df[df["P"] < P_THRESH].sort_values("P")

print(f"  {len(df)} SNPs below p < {P_THRESH:.0e}")
if df.empty:
    print("No significant hits — exiting.")
    raise SystemExit

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
MYVARIANT_URL  = "https://myvariant.info/v1/variant"
MYVARIANT_GENE = "https://myvariant.info/v1/variant"
GWAS_CAT_URL   = "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms"
UNIPROT_URL    = "https://mygene.info/v3/gene"

HEADERS = {"User-Agent": "ADNI-SNP-Annotator/1.0 (academic research)"}


def mv_query_rsid(rsid: str) -> dict:
    """Query MyVariant.info by rsID."""
    fields = "dbsnp.gene,dbsnp.vartype,cadd.annotype,gnomad_exome.af,gnomad_genome.af," \
             "snpeff.ann"
    try:
        r = requests.get(
            f"{MYVARIANT_URL}/{rsid}",
            params={"fields": fields, "assembly": "hg19"},
            headers=HEADERS, timeout=15
        )
        if r.status_code == 200:
            return r.json()
    except Exception as e:
        print(f"    MyVariant rsID query failed: {e}")
    return {}


def mv_query_chrpos(chrom: int, bp: int) -> dict:
    """Query MyVariant.info by chr:pos (for kgp* IDs)."""
    variant_id = f"chr{chrom}:g.{bp}"
    fields = "dbsnp.gene,dbsnp.vartype,cadd.annotype,gnomad_exome.af,gnomad_genome.af,snpeff.ann"
    try:
        r = requests.get(
            f"{MYVARIANT_URL}/{variant_id}",
            params={"fields": fields, "assembly": "hg19"},
            headers=HEADERS, timeout=15
        )
        if r.status_code == 200:
            return r.json()
    except Exception as e:
        print(f"    MyVariant chrpos query failed: {e}")
    # Fallback: search by position
    try:
        r = requests.get(
            f"https://myvariant.info/v1/query",
            params={
                "q": f"chr{chrom}:{bp}",
                "fields": fields,
                "assembly": "hg19",
                "size": 1,
            },
            headers=HEADERS, timeout=15,
        )
        if r.status_code == 200:
            hits = r.json().get("hits", [])
            if hits:
                return hits[0]
    except Exception as e:
        print(f"    MyVariant fallback search failed: {e}")
    return {}


def extract_gene(mv: dict) -> str:
    """Extract nearest gene symbol from MyVariant response."""
    # Try dbsnp first
    gene_info = mv.get("dbsnp", {}).get("gene", None)
    if gene_info:
        if isinstance(gene_info, list):
            return ", ".join(g.get("symbol", "") for g in gene_info if g.get("symbol"))
        if isinstance(gene_info, dict):
            return gene_info.get("symbol", "N/A")
    # Try snpeff annotations
    ann = mv.get("snpeff", {}).get("ann", None)
    if ann:
        if isinstance(ann, list):
            return ann[0].get("gene_name", "N/A")
        if isinstance(ann, dict):
            return ann.get("gene_name", "N/A")
    return "N/A"


def extract_annotype(mv: dict) -> str:
    """Extract functional annotation (intronic/exonic/intergenic)."""
    # cadd.annotype
    annotype = mv.get("cadd", {}).get("annotype", None)
    if annotype:
        return str(annotype)
    # snpeff
    ann = mv.get("snpeff", {}).get("ann", None)
    if ann:
        if isinstance(ann, list):
            return ann[0].get("effect", "N/A")
        if isinstance(ann, dict):
            return ann.get("effect", "N/A")
    return "N/A"


def extract_af(mv: dict) -> str:
    """Extract allele frequency (gnomAD genome EUR preferred)."""
    # gnomad_genome has population-level AFs
    gaf = mv.get("gnomad_genome", {})
    if gaf:
        # Try af (overall)
        af = gaf.get("af", {})
        if isinstance(af, dict):
            eur = af.get("af_eur", af.get("af", None))
            if eur is not None:
                return f"{eur:.4f} (gnomAD EUR)"
        if isinstance(af, (float, int)):
            return f"{af:.4f} (gnomAD overall)"
    # Fallback: gnomad_exome
    gaf_ex = mv.get("gnomad_exome", {})
    if gaf_ex:
        af = gaf_ex.get("af", {})
        if isinstance(af, dict):
            eur = af.get("af_eur", af.get("af", None))
            if eur is not None:
                return f"{eur:.4f} (gnomAD exome EUR)"
    return "N/A"


def query_gwas_catalog(rsid: str) -> str:
    """Query GWAS Catalog for prior associations with this rsID."""
    if not rsid.startswith("rs"):
        return "kgp ID – no GWAS Catalog lookup"
    try:
        r = requests.get(
            f"{GWAS_CAT_URL}/{rsid}/associations",
            headers={**HEADERS, "Accept": "application/json"},
            timeout=15,
        )
        if r.status_code == 404:
            return "No prior GWAS associations"
        if r.status_code == 200:
            data = r.json()
            embedded = data.get("_embedded", {})
            assocs   = embedded.get("associations", [])
            if not assocs:
                return "No prior GWAS associations"
            # Summarise unique traits
            traits = []
            for a in assocs[:5]:  # cap at 5
                for trait in a.get("efoTraits", []):
                    traits.append(trait.get("trait", ""))
            traits = list(dict.fromkeys(t for t in traits if t))
            return "; ".join(traits) if traits else f"{len(assocs)} associations (see GWAS Catalog)"
    except Exception as e:
        return f"Query error: {e}"
    return "N/A"


# ─────────────────────────────────────────────────────────────────────────────
# 2. Annotate each top hit
# ─────────────────────────────────────────────────────────────────────────────
rows = []
for _, snp_row in df.iterrows():
    snp  = snp_row["SNP"]
    chrom = int(snp_row["CHR"]) if pd.notna(snp_row["CHR"]) else None
    bp   = int(snp_row["BP"])   if pd.notna(snp_row["BP"])   else None

    print(f"\nAnnotating {snp}  (chr{chrom}:{bp}  p={snp_row['P']:.2e})")

    # Query MyVariant
    if snp.startswith("rs"):
        mv = mv_query_rsid(snp)
    else:
        mv = mv_query_chrpos(chrom, bp)

    gene     = extract_gene(mv)
    annotype = extract_annotype(mv)
    af_str   = extract_af(mv)

    # Query GWAS Catalog
    gwas_evidence = query_gwas_catalog(snp)
    print(f"    Gene={gene}  Annotation={annotype}  AF={af_str}")
    print(f"    GWAS evidence: {gwas_evidence}")

    rows.append({
        "SNP":            snp,
        "CHR":            chrom,
        "BP":             bp,
        "A1":             snp_row.get("A1", ""),
        "OR":             snp_row.get("OR", ""),
        "SE":             snp_row.get("SE", ""),
        "P":              snp_row["P"],
        "Nearest_Gene":   gene,
        "Annotation":     annotype,
        "AF_gnomAD":      af_str,
        "GWAS_Evidence":  gwas_evidence,
    })

    time.sleep(0.3)   # polite rate limiting

# ─────────────────────────────────────────────────────────────────────────────
# 3. Save outputs
# ─────────────────────────────────────────────────────────────────────────────
results = pd.DataFrame(rows)

# TSV
tsv_path = OUT_DIR / "top_snps_annotated.tsv"
results.to_csv(tsv_path, sep="\t", index=False)
print(f"\n✓ TSV saved → {tsv_path}")

# HTML table (styled)
html_path = OUT_DIR / "top_snps_annotated.html"
html_style = """
<style>
body { font-family: 'Segoe UI', Arial, sans-serif; background: #f8f9fa; }
h2 { color: #2c3e50; }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
th { background: #2c3e50; color: white; padding: 8px 12px; text-align: left; }
td { padding: 6px 12px; border-bottom: 1px solid #dee2e6; }
tr:nth-child(even) { background: #f2f5f9; }
tr:hover { background: #d6e4f0; }
.hit { color: #c0392b; font-weight: bold; }
</style>
"""

with open(html_path, "w") as f:
    f.write(f"<html><head><title>GWAS Top Hits</title>{html_style}</head><body>")
    f.write(f"<h2>ADNI GWAS Top SNPs (p &lt; {P_THRESH:.0e})</h2>")
    f.write(results.to_html(index=False, classes="annotated-table", border=0,
                             float_format=lambda x: f"{x:.2e}" if abs(x) < 0.001 else f"{x:.4f}"
                             if isinstance(x, float) else x))
    f.write("</body></html>")

print(f"✓ HTML saved → {html_path}")

# ─────────────────────────────────────────────────────────────────────────────
# 4. Print summary
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'═'*80}")
print("ANNOTATION SUMMARY")
print(f"{'═'*80}")
with pd.option_context("display.max_colwidth", 50, "display.width", 200):
    print(results[["SNP", "CHR", "BP", "OR", "P",
                   "Nearest_Gene", "Annotation",
                   "AF_gnomAD", "GWAS_Evidence"]].to_string(index=False))
print(f"{'═'*80}")
