"""
update_snp_ids.py
Convert kgp* probe IDs in a PLINK BIM file to dbSNP rs-IDs.

Strategy (priority cascade)
----------------------------
Each step only processes kgp SNPs not yet resolved by the previous step.

  Step 1  --b138      (RECOMMENDED — fast, offline, unambiguous)
          Direct probe-name → rsID from the Illumina b138 companion file
          (HumanOmni2-5-8-v1-2-A-b138-rsIDs.txt, Name<TAB>RsID).
          No positional or allele check needed — the kgp probe name is a
          globally unique key that Illumina mapped to dbSNP build 138.
          This resolves the vast majority of kgp IDs in seconds.

  Step 2  --kgp-pvar  (RECOMMENDED fallback — allele-aware positional)
          PLINK 1000-Genomes phase-3 build-37 pvar file (all_phase3.pvar).
          Matches remaining kgp probes by (CHR, BP) and verifies that the
          BIM alleles are consistent with the reference REF/ALT alleles on
          the same strand or their complement strand (A<->T, C<->G).
          Palindromic SNPs (A/T or C/G) are accepted when position+alleles
          agree, and are flagged 'kgp_pvar_palindrome' in the lookup table.
          Download: all_phase3.pvar.zst from
            https://www.cog-genomics.org/plink/2.0/resources
          Decompress with:
            plink2 --zst-decompress all_phase3.pvar.zst all_phase3.pvar
          Only the pvar file is needed — pgen and psam are NOT required.
          
          Advice: Use all_phase3_noannot.pvar.zst for a smaller file/faster processing.
          No need to decompress, script handles compressed files automatically.

  Step 3  --manifest  (optional — only if Steps 1+2 leave many unresolved)
          Large Illumina manifest CSV (~1 GB).  Useful if it contains an
          explicit RsID column; otherwise falls back to positional + allele
          matching against BIM rsID SNPs (weaker than KGP pvar).
          NOTE: The --manifest step is largely redundant when --b138 and
          --kgp-pvar are both provided.

  Step 4  --use-api   (opt-in online fallback)
          MyVariant.info REST API, queried by chr:pos with allele check.
          Use only for the long tail of unresolved probes.

Unresolved kgp IDs are written to  <bim_dir>/unmatched_kgp.tsv.

Outputs
-------
  <bim_dir>/update_name.txt       PLINK --update-name input
  <bim_dir>/kgp_to_rs_lookup.tsv  full lookup table (includes method column)
  <bim_dir>/unmatched_kgp.tsv     probes that could not be resolved

Recommended usage (hg19/GRCh37 BIM)
-------------------------------------
  # 1. Decompress the KGP pvar (one-time, ~614 MB → ~3 GB uncompressed)
  plink2 --zst-decompress all_phase3.pvar.zst all_phase3.pvar

  # 2. Run conversion
  python update_snp_ids.py \\
      --bim <stem>.bim \\
      --b138 D:/.../HumanOmni2-5-8-v1-2-A-b138-rsIDs.txt \\
      --kgp-pvar D:/.../all_phase3.pvar

  # 3. Apply to PLINK dataset
  plink --bfile <stem> \\
        --update-name update_name.txt \\
        --make-bed --out <stem>_rsid
"""

import argparse
import pathlib
import re
import sys
import time

import pandas as pd

# ── CLI ──────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Convert kgp* probe IDs in a PLINK BIM to dbSNP rs-IDs (allele-aware).",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument("--bim", required=True,
    help="Path to PLINK .bim file (hg19/GRCh37).")
parser.add_argument("--b138", default=None,
    help="Path to HumanOmni2-5-8-v1-2-A-b138-rsIDs.txt  (Name<TAB>RsID).")
parser.add_argument("--extra-rsid-files", nargs="*", default=[], dest="extra_rsid_files",
    help=("Additional Illumina rsID companion files in the same Name<TAB>RsID format, "
          "tried in order on probes still unresolved after --b138. "
          "e.g. --extra-rsid-files b144_rsids.txt b150_rsids.txt b151_rsids.txt"))
parser.add_argument("--manifest", default=None,
    help="Path to Illumina manifest CSV or ZIP (optional supplement to b138).")
parser.add_argument("--kgp-pvar", default=None, dest="kgp_pvar",
    help=("Path to KGP phase-3 build-37 pvar file, e.g. all_phase3_ns.pvar. "
          "Download all_phase3_ns_noannot.pvar.zst from "
          "https://www.cog-genomics.org/plink/2.0/resources, decompress "
          "(plink2 --zst-decompress ...), rename to all_phase3_ns.pvar."))
parser.add_argument("--use-api", action="store_true", dest="use_api",
    help="Fall back to MyVariant.info REST API for still-unresolved probes.")
parser.add_argument("--validate-rsids", action="store_true", dest="validate_rsids",
    help=("After conversion, batch-check assigned rsIDs against MyVariant.info to detect "
          "deprecated/merged IDs (dbSNP build 138 rsIDs may have been superseded). "
          "Results written to deprecated_rsids.tsv. Requires internet access. "
          "Not needed for PLINK GWAS — only relevant before VEP or GWAS-Catalog annotation."))
parser.add_argument("--batch-size", type=int, default=1000, dest="batch_size",
    help="Batch size for MyVariant.info API requests (default: 1000).")
parser.add_argument("--assembly", default="hg19", choices=["hg19", "hg38"],
    help="Genome assembly of the BIM file (default: hg19).")
args = parser.parse_args()

BIM_PATH = pathlib.Path(args.bim)
if not BIM_PATH.exists():
    print(f"[ERROR] BIM file not found: {BIM_PATH}")
    sys.exit(1)

OUT_DIR           = BIM_PATH.parent
UPDATE_NAME_FILE  = OUT_DIR / "update_name.txt"
LOOKUP_TSV        = OUT_DIR / "kgp_to_rs_lookup.tsv"
UNMATCHED_TSV     = OUT_DIR / "unmatched_kgp.tsv"
DEPRECATED_TSV    = OUT_DIR / "deprecated_rsids.tsv"


# ── 1. Load BIM ───────────────────────────────────────────────────────────────
print(f"\nLoading BIM: {BIM_PATH}")
bim = pd.read_csv(BIM_PATH, sep="\t", header=None,
                  names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
                  dtype={"CHR": str, "SNP": str, "CM": str,
                         "BP": int, "A1": str, "A2": str})

kgp_mask = bim["SNP"].str.startswith("kgp", na=False)
kgp_snps = bim[kgp_mask].copy()

print(f"  Total SNPs : {len(bim):>10,}")
print(f"  rs-ID SNPs : {len(bim[~kgp_mask]):>10,}  (no change needed)")
print(f"  kgp SNPs   : {len(kgp_snps):>10,}  (to be converted)")

if kgp_snps.empty:
    print("\nNo kgp SNPs found — nothing to do.")
    sys.exit(0)


# ── 2. Allele helpers ─────────────────────────────────────────────────────────
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
_MISSING = {"0", ".", "N", ""}


def _comp_set(allele_set):
    """Return the complement-strand version of an allele set."""
    return {_COMP.get(a, a) for a in allele_set}


def _clean(a1, a2):
    """Return frozenset of non-missing alleles (uppercase)."""
    return frozenset({str(a1).upper(), str(a2).upper()} - _MISSING)


def _is_palindromic(allele_set):
    clean = allele_set - _MISSING
    return len(clean) == 2 and _comp_set(clean) == clean


def _alleles_match(bim_a1, bim_a2, ref_a, alt_a):
    """
    Return (match: bool, palindromic: bool).
    Accepts same-strand or complement-strand allele pairs.
    """
    b = _clean(bim_a1, bim_a2)
    r = _clean(ref_a, alt_a)
    if not b or not r:
        return False, False
    palindromic = _is_palindromic(b)
    match = (b == r) or (_comp_set(b) == r)
    return match, palindromic


def _parse_illumina_snp(snp_str):
    """Parse Illumina '[A/G]' SNP column → (allele1, allele2) or (None, None)."""
    m = re.match(r"\[([ACGTI])/([ACGTI])\]", str(snp_str).strip().upper())
    if m:
        a1, a2 = m.group(1), m.group(2)
        # I/D (insertion/deletion) — treat as indel, skip allele check
        if "I" in (a1, a2) or "D" in (a1, a2):
            return None, None
        return a1, a2
    return None, None


def _norm_chrom(c):
    """Normalise chromosome label to a plain numeric string (e.g. 'chr1' → '1', 'X' → '23')."""
    c = str(c).strip().upper()
    # Remove leading 'CHR' prefix robustly
    if c.startswith("CHR"):
        c = c[3:]
    mapping = {"X": "23", "Y": "24", "XY": "25", "M": "26", "MT": "26"}
    return mapping.get(c, c)


# ── Step 1: b138 direct probe → rsID ─────────────────────────────────────────
def method_b138(b138_path, snps_df):
    """
    Read Illumina b138 rsIDs companion file (Name TAB RsID).
    No positional or allele matching needed — probe name is the unique key.
    """
    b138_path = pathlib.Path(b138_path)
    print(f"\n[Step 1 – b138] {b138_path.name}")

    b138 = pd.read_csv(b138_path, sep="\t", usecols=["Name", "RsID"],
                       dtype=str, low_memory=False)
    b138 = b138.dropna(subset=["Name", "RsID"])
    b138["RsID"] = b138["RsID"].str.strip()
    b138 = b138[b138["RsID"].str.match(r"^rs\d+$", na=False)]

    probe_to_rs = dict(zip(b138["Name"], b138["RsID"]))

    rows = []
    for _, row in snps_df.iterrows():
        rsid = probe_to_rs.get(row["SNP"])
        if rsid:
            rows.append({"kgp_id": row["SNP"], "new_id": rsid, "method": "b138"})

    df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=["kgp_id","new_id","method"])
    print(f"  Resolved {len(df):,} / {len(snps_df):,} kgp probes")
    return df


# ── Step 2: Illumina manifest (allele-aware) ──────────────────────────────────
def _open_manifest(manifest_path):
    import io, zipfile as zf
    if manifest_path.suffix.lower() == ".zip":
        z = zf.ZipFile(manifest_path)
        csvs = [n for n in z.namelist() if n.lower().endswith(".csv")]
        if not csvs:
            print(f"[ERROR] No CSV in {manifest_path.name}"); sys.exit(1)
        print(f"  Extracting '{csvs[0]}' from zip …")
        return io.TextIOWrapper(z.open(csvs[0]), encoding="utf-8", errors="replace")
    return open(manifest_path, "r", encoding="utf-8", errors="replace")


def method_manifest(manifest_path, snps_df):
    """
    Parse Illumina manifest CSV.  Uses RsID column when available;
    otherwise falls back to positional + allele-verified matching.
    """
    manifest_path = pathlib.Path(manifest_path)
    print(f"\n[Step 2 – manifest] {manifest_path.name}")

    # Find header row
    header_row = None
    with _open_manifest(manifest_path) as f:
        for i, line in enumerate(f):
            if line.startswith("Name,") or line.startswith("Name\t"):
                header_row = i; break
    if header_row is None:
        print("  [WARNING] Header row not found — assuming row 0.")
        header_row = 0

    with _open_manifest(manifest_path) as f:
        mfst = pd.read_csv(f, sep=",", skiprows=header_row, low_memory=False)

    mfst = mfst[mfst["Name"].notna()]
    mfst = mfst[~mfst["Name"].astype(str).str.startswith("[")]
    mfst.columns = mfst.columns.str.strip()
    print(f"  Rows: {len(mfst):,}  |  Cols: {list(mfst.columns[:10])}")

    # Restrict to kgp probes in our remaining list
    want = set(snps_df["SNP"])
    mfst_kgp = mfst[mfst["Name"].isin(want)].copy()

    # ── Branch A: explicit RsID column ───────────────────────────────────────
    rs_col = next(
        (c for c in mfst.columns
         if c.lower().replace(" ", "").replace("_", "") in ("rsid","dbsnprsid","snprsid")),
        None)

    if rs_col:
        print(f"  Using RsID column: '{rs_col}'")
        mfst_kgp["_rsid"] = mfst_kgp[rs_col].astype(str).str.strip()
        rows = []
        for _, r in mfst_kgp.iterrows():
            rsid = r["_rsid"]
            if re.match(r"^rs\d+$", rsid):
                rows.append({"kgp_id": r["Name"], "new_id": rsid, "method": "manifest_rsid"})
        df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=["kgp_id","new_id","method"])
        print(f"  Resolved {len(df):,} / {len(snps_df):,} kgp probes")
        return df

    # ── Branch B: positional + allele matching ────────────────────────────────
    print("  No RsID column — positional + allele matching")
    pos_col = next((c for c in mfst.columns
                    if c.lower() in ("mapinfo","position","bp","chromosomecoordinate")), None)
    chr_col = next((c for c in mfst.columns
                    if c.lower() in ("chr","chromosome","chrom")), None)
    snp_col = next((c for c in mfst.columns if c.lower() == "snp"), None)

    if not pos_col or not chr_col:
        print(f"  [WARNING] Cannot find Chr/MapInfo columns — skipping manifest.")
        return pd.DataFrame(columns=["kgp_id","new_id","method"])

    # Build lookup from BIM rs-ID SNPs: (norm_chrom, pos) → (rsid, a1, a2)
    bim_rs = bim[~bim["SNP"].str.startswith("kgp", na=False)].copy()
    bim_rs["_chrom"] = bim_rs["CHR"].apply(_norm_chrom)
    bim_rs_lookup = {}
    for _, r in bim_rs.iterrows():
        bim_rs_lookup.setdefault((r["_chrom"], r["BP"]), []).append(
            (r["SNP"], str(r["A1"]).upper(), str(r["A2"]).upper()))

    rows = []
    for _, r in mfst_kgp.iterrows():
        chrom = _norm_chrom(r[chr_col])
        try:
            pos = int(r[pos_col])
        except (ValueError, TypeError):
            continue

        a1_m, a2_m = _parse_illumina_snp(r[snp_col]) if snp_col else (None, None)

        candidates = bim_rs_lookup.get((chrom, pos), [])
        for rsid, b_a1, b_a2 in candidates:
            if a1_m is None:
                # No allele info — accept position match only (weaker)
                rows.append({"kgp_id": r["Name"], "new_id": rsid,
                             "method": "manifest_pos"})
                break
            ok, palindromic = _alleles_match(a1_m, a2_m, b_a1, b_a2)
            if ok and not palindromic:
                rows.append({"kgp_id": r["Name"], "new_id": rsid,
                             "method": "manifest_pos_allele"})
                break
            elif ok and palindromic:
                # Palindromic: accept but flag
                rows.append({"kgp_id": r["Name"], "new_id": rsid,
                             "method": "manifest_pos_palindrome"})
                break

    df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=["kgp_id","new_id","method"])
    print(f"  Resolved {len(df):,} / {len(snps_df):,} kgp probes")
    return df


# ── Step 3: KGP phase-3 build-37 pvar (allele-aware) ─────────────────────────
import contextlib

@contextlib.contextmanager
def _open_pvar(path):
    """
    Open a pvar file for line-by-line text reading, handling .zst compression
    transparently.  No need to decompress with plink2 first.

    Requires 'zstandard' for .zst files:  pip install zstandard
    """
    import io
    path = pathlib.Path(path)
    if path.suffix.lower() == ".zst":
        try:
            import zstandard
        except ImportError:
            print("[ERROR] zstandard not installed.  Run:  pip install zstandard")
            print("        Or decompress first and pass the plain .pvar file.")
            sys.exit(1)
        raw = open(path, "rb")
        reader = zstandard.ZstdDecompressor().stream_reader(raw)
        try:
            yield io.TextIOWrapper(reader, encoding="utf-8")
        finally:
            reader.close()
            raw.close()
    else:
        fh = open(path, "r", encoding="utf-8")
        try:
            yield fh
        finally:
            fh.close()


def method_kgp_pvar(pvar_path, snps_df):
    """
    Read a PLINK pvar file (tab-separated, header line starts with #CHROM).
    Accepts both plain .pvar and compressed .pvar.zst directly
    (no plink2 needed — requires: pip install zstandard).
    Builds (norm_chrom, pos, frozenset{REF, ALT}) → rsID lookup.
    Only accepts entries where ID starts with 'rs'.
    Palindromic matches are accepted but flagged in the method column.
    """
    pvar_path = pathlib.Path(pvar_path)
    if not pvar_path.exists():
        print(f"[ERROR] KGP pvar not found: {pvar_path}"); sys.exit(1)

    print(f"\n[Step 3 – KGP pvar] {pvar_path.name}  (loading …)")
    if pvar_path.suffix.lower() == ".zst":
        print("  Streaming directly from .zst (no decompression needed)")

    # Count ## metadata lines to skip; the #CHROM line becomes the data header
    n_skip = 0
    with _open_pvar(pvar_path) as fh:
        for line in fh:
            if line.startswith("##"):
                n_skip += 1
            else:
                break  # first non-## line is the #CHROM header

    # Build lookup directly from file to conserve memory
    pos_lookup = {}   # (norm_chrom, pos) → [(allele_frozenset, rsid), …]
    n_loaded = 0
    with _open_pvar(pvar_path) as fh:
        for _ in range(n_skip):
            next(fh)
        # Next line is the #CHROM header — skip it
        next(fh)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom_raw, pos_raw, vid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            if not vid.startswith("rs"):
                continue
            try:
                pos = int(pos_raw)
            except ValueError:
                continue
            chrom = _norm_chrom(chrom_raw)
            key = (chrom, pos)
            alleles = frozenset({ref.upper(), alt.upper()} - _MISSING)
            pos_lookup.setdefault(key, []).append((alleles, vid))
            n_loaded += 1

    print(f"  Indexed {n_loaded:,} rsID-annotated variants from KGP pvar")

    rows = []
    n_palindrome = 0
    for _, row in snps_df.iterrows():
        chrom = _norm_chrom(row["CHR"])
        pos   = int(row["BP"])
        a1    = str(row["A1"]).upper()
        a2    = str(row["A2"]).upper()
        key   = (chrom, pos)

        for (allele_set, rsid) in pos_lookup.get(key, []):
            ok, palindromic = _alleles_match(a1, a2, *allele_set)
            if ok:
                method = "kgp_pvar_palindrome" if palindromic else "kgp_pvar"
                if palindromic:
                    n_palindrome += 1
                rows.append({"kgp_id": row["SNP"], "new_id": rsid, "method": method})
                break  # take first matching variant at this position

    df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=["kgp_id","new_id","method"])
    print(f"  Resolved {len(df):,} / {len(snps_df):,} kgp probes "
          f"({n_palindrome} palindromic)")
    return df


# ── Step 4: MyVariant.info API (opt-in) ───────────────────────────────────────
def method_api(snps_df, batch_size=1000, assembly="hg19"):
    """
    Query MyVariant.info in batches by chr:pos.
    Only used when --use-api is specified.
    """
    try:
        import requests
    except ImportError:
        print("[ERROR] requests not installed: pip install requests"); sys.exit(1)

    print(f"\n[Step 4 – API] MyVariant.info  ({len(snps_df):,} probes, {assembly})")
    url  = "https://myvariant.info/v1/query"
    rows = snps_df[["SNP","CHR","BP","A1","A2"]].reset_index(drop=True)
    n    = len(rows)
    results = []

    for start in range(0, n, batch_size):
        batch = rows.iloc[start:start+batch_size]
        queries = [f"chr{r.CHR}:g.{r.BP}" for _, r in batch.iterrows()]
        try:
            resp = requests.post(url, json={
                "q": queries,
                "fields": "dbsnp.rsid,dbsnp.ref,dbsnp.alt",
                "assembly": assembly,
                "size": 1,
                "dotfield": True,
            }, timeout=30, headers={"Content-Type": "application/json"})
            data = resp.json()
        except Exception as e:
            print(f"  [WARNING] Batch {start}–{start+batch_size} failed: {e}")
            for _, r in batch.iterrows():
                results.append({"kgp_id": r.SNP, "new_id": r.SNP, "method": "api_fail"})
            continue

        for i, hit in enumerate(data):
            r = batch.iloc[i]
            rsid = None
            if hit.get("notfound") is not True:
                rsid = hit.get("dbsnp.rsid")
                if isinstance(rsid, list):
                    rsid = rsid[0]
                if rsid and not str(rsid).startswith("rs"):
                    rsid = f"rs{rsid}"
                # Allele check
                if rsid:
                    ref_api = hit.get("dbsnp.ref", "")
                    alt_api = hit.get("dbsnp.alt", "")
                    ok, palindromic = _alleles_match(r.A1, r.A2, ref_api, alt_api)
                    if not ok:
                        rsid = None  # reject mismatched allele
            method = "api_palindrome" if (rsid and palindromic) else "api"
            results.append({"kgp_id": r.SNP,
                            "new_id": rsid if rsid else r.SNP,
                            "method": method if rsid else "api_no_match"})

        done = min(start + batch_size, n)
        print(f"  {done:>7,} / {n:,}  ({100*done/n:.0f}%)", flush=True)
        time.sleep(0.1)

    df = pd.DataFrame(results)
    resolved = df[df["kgp_id"] != df["new_id"]]
    print(f"  Resolved {len(resolved):,} / {len(snps_df):,} kgp probes via API")
    # Only return actually resolved rows
    return resolved.reset_index(drop=True)


# ── Cascade ───────────────────────────────────────────────────────────────────
print(f"\n{'─'*50}")
print(f"Running kgp → rsID cascade  ({len(kgp_snps):,} probes)")
print(f"{'─'*50}")

resolved_parts = []
remaining = kgp_snps.copy()

if args.b138 and not remaining.empty:
    part = method_b138(args.b138, remaining)
    if not part.empty:
        resolved_parts.append(part)
        remaining = remaining[~remaining["SNP"].isin(part["kgp_id"])]

# Step 1b-1d: extra rsID files (b144, b150, b151, …) — same format as b138
for extra_path in (args.extra_rsid_files or []):
    if remaining.empty:
        break
    part = method_b138(extra_path, remaining)
    if not part.empty:
        # re-label method to reflect which build the match came from
        part["method"] = pathlib.Path(extra_path).stem
        resolved_parts.append(part)
        remaining = remaining[~remaining["SNP"].isin(part["kgp_id"])]

if args.manifest and not remaining.empty:
    part = method_manifest(args.manifest, remaining)
    if not part.empty:
        resolved_parts.append(part)
        remaining = remaining[~remaining["SNP"].isin(part["kgp_id"])]

if args.kgp_pvar and not remaining.empty:
    part = method_kgp_pvar(args.kgp_pvar, remaining)
    if not part.empty:
        resolved_parts.append(part)
        remaining = remaining[~remaining["SNP"].isin(part["kgp_id"])]

if args.use_api and not remaining.empty:
    part = method_api(remaining, args.batch_size, args.assembly)
    if not part.empty:
        resolved_parts.append(part)
        remaining = remaining[~remaining["SNP"].isin(part["kgp_id"])]


# ── Optional: validate assigned rsIDs for deprecation/merges ─────────────────
def validate_rsids(lookup_df, batch_size=1000, assembly="hg19"):
    """
    Batch-check all assigned rsIDs against MyVariant.info.
    dbSNP periodically merges rsIDs (old rsID superseded by a newer one).
    This is orthogonal to coordinate liftover: rsIDs are assembly-agnostic,
    so lifting hg19→hg38 does NOT fix stale rsIDs.

    Returns a DataFrame of deprecated entries:
        kgp_id | old_rsid | current_rsid | status
    where status is 'merged' (replaced) or 'not_found' (removed from dbSNP).
    """
    try:
        import requests
    except ImportError:
        print("[ERROR] requests not installed: pip install requests"); sys.exit(1)

    rs_rows = lookup_df[lookup_df["new_id"].str.match(r"^rs\d+$", na=False)].copy()
    rsids   = rs_rows["new_id"].tolist()
    n       = len(rsids)
    print(f"\n[Validate rsIDs] Checking {n:,} assigned rsIDs against MyVariant.info …")
    print(f"  Note: rsID deprecation is independent of hg19/hg38 coordinate liftover.")

    url      = "https://myvariant.info/v1/query"
    deprecated = []  # rows: (kgp_id, old_rsid, current_rsid, status)

    # Build reverse map: new_id → kgp_id (for reporting)
    rsid_to_kgp = dict(zip(rs_rows["new_id"], rs_rows["kgp_id"]))

    for start in range(0, n, batch_size):
        batch = rsids[start : start + batch_size]
        try:
            resp = requests.post(
                url,
                json={"q": batch, "fields": "dbsnp.rsid",
                      "assembly": assembly, "dotfield": True},
                timeout=30,
                headers={"Content-Type": "application/json"},
            )
            data = resp.json()
        except Exception as e:
            print(f"  [WARNING] Validation batch {start}–{start+batch_size} failed: {e}")
            continue

        for i, hit in enumerate(data):
            queried = batch[i]
            kgp_id  = rsid_to_kgp.get(queried, "?")
            if hit.get("notfound"):
                deprecated.append({"kgp_id": kgp_id, "old_rsid": queried,
                                   "current_rsid": None, "status": "not_found"})
                continue
            current = hit.get("dbsnp.rsid")
            if isinstance(current, list):
                current = current[0]
            if current is None:
                continue
            current_str = str(current) if str(current).startswith("rs") else f"rs{current}"
            if current_str != queried:
                deprecated.append({"kgp_id": kgp_id, "old_rsid": queried,
                                   "current_rsid": current_str, "status": "merged"})

        done = min(start + batch_size, n)
        print(f"  {done:>7,} / {n:,}  ({100*done/n:.0f}%)", flush=True)
        time.sleep(0.1)

    df = pd.DataFrame(deprecated) if deprecated else pd.DataFrame(
        columns=["kgp_id", "old_rsid", "current_rsid", "status"])
    print(f"  Deprecated/merged rsIDs found: {len(df):,} "
          f"({len(df[df.status=='merged']) if not df.empty else 0} merged, "
          f"{len(df[df.status=='not_found']) if not df.empty else 0} not found)")
    return df


# ── Outputs ───────────────────────────────────────────────────────────────────
print(f"\n{'─'*50}")
print(f"Writing outputs to {OUT_DIR}")

# Full lookup table
if resolved_parts:
    lookup = pd.concat(resolved_parts, ignore_index=True)
else:
    lookup = pd.DataFrame(columns=["kgp_id", "new_id", "method"])

lookup.to_csv(LOOKUP_TSV, sep="\t", index=False)
print(f"  Lookup table    → {LOOKUP_TSV.name}  ({len(lookup):,} rows)")

# PLINK --update-name file  (old_id  new_id, only changed rows)
changed = lookup[lookup["kgp_id"] != lookup["new_id"]][["kgp_id", "new_id"]]
changed.to_csv(UPDATE_NAME_FILE, sep="\t", header=False, index=False)
print(f"  update_name.txt → {UPDATE_NAME_FILE.name}  ({len(changed):,} SNPs renamed)")

# Unmatched probes
unmatched = remaining[["SNP", "CHR", "BP", "A1", "A2"]].rename(columns={"SNP": "kgp_id"})
unmatched.to_csv(UNMATCHED_TSV, sep="\t", index=False)
print(f"  unmatched_kgp   → {UNMATCHED_TSV.name}  ({len(unmatched):,} probes)")

# Optional: rsID deprecation check
if args.validate_rsids and not lookup.empty:
    deprecated_df = validate_rsids(lookup, args.batch_size, args.assembly)
    deprecated_df.to_csv(DEPRECATED_TSV, sep="\t", index=False)
    print(f"  deprecated_rsids→ {DEPRECATED_TSV.name}  ({len(deprecated_df):,} rows)")
    if not deprecated_df.empty:
        print(f"  [NOTE] Review deprecated_rsids.tsv before downstream annotation "
              f"(VEP, GWAS Catalog). For PLINK GWAS analysis, stale rsIDs are "
              f"harmless — they are labels only and dbSNP maintains merge history.")


# ── Summary ───────────────────────────────────────────────────────────────────
total     = len(kgp_snps)
n_res     = len(changed)
n_unmatch = len(unmatched)

print(f"\n── Summary {'─'*40}")
print(f"  kgp SNPs in BIM      : {total:>10,}")
print(f"  Resolved to rsID     : {n_res:>10,}  ({100*n_res/total:.1f}%)")
print(f"  Unresolved (kgp kept): {n_unmatch:>10,}  ({100*n_unmatch/total:.1f}%)")

if resolved_parts:
    print("\n  Breakdown by method:")
    for method, grp in lookup.groupby("method"):
        print(f"    {method:<30s}: {len(grp):>8,}")

print(f"{'─'*50}")

print(f"\nTo apply to your PLINK dataset:")
stem = BIM_PATH.with_suffix("")
print(f"  plink --bfile {stem} \\")
print(f"        --update-name {UPDATE_NAME_FILE} \\")
print(f"        --make-bed \\")
print(f"        --out {stem}_rsid")
