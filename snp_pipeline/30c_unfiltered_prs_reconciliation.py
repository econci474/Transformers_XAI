r"""
30c_unfiltered_prs_reconciliation.py   (LOCAL, env: snp)
========================================================
Disambiguate the "not_on_array_or_below_maf" drop reason from 30b by
checking each dropped published-PRS SNP against the **original
unfiltered Omni2.5M genotype file** (`WGS_Omni25_BIN_wo_ConsentsIssues.bim`),
which represents the raw Illumina manifest before any QC (missingness,
HWE, MAF, LD-prune, strand/ref-correction).

The unfiltered BIM stores most SNPs under Illumina internal **kgp probe
names**; `05a_update_snp_ids.py` renames them to rsIDs downstream using
the Illumina dbSNP-build companion manifests at `D:\…\liftover\Illumina\`
(b138 primary, then b144/b150/b151 supplements). To correctly recognise a
published rsID on the chip we therefore need the reverse map (rsID → kgp),
not just direct rsID matching.

Resolution chain (per published SNP):
  1. direct rsID in unfiltered BIM
  2. RsMergeArch-canonical rsID in unfiltered BIM
  3. Illumina kgp name (from b138/b144/b150/b151 manifests) in unfiltered BIM
  4. positional match in GRCh37 (Wightman/Kunkle/Desikan only; the
     unfiltered BIM is GRCh37 — confirmed via spot-checks against
     rs2075650/rs1476679/rs3851179)

Reconciliation categories:
  - on_chip_passed_maf                       : in MAF BIM AND in unfilt BIM
  - on_chip_filtered_out                     : direct rsID in unfilt; lost to QC
  - on_chip_via_rsID_merge_filtered_out      : only via RsMergeArch canonical
  - on_chip_via_kgp_renamed_filtered_out     : only via Illumina kgp manifest
  - on_chip_via_position_filtered_out        : only via GRCh37 position
  - in_maf_but_not_on_chip_manifest          : in MAF BIM but absent from
                                               unfilt BIM by all four methods
                                               (genuine imputation / extract)
  - not_on_chip                              : absent from unfilt BIM and
                                               not in MAF either
  - apoe_haplotype_special                   : Desikan APOE haplotype (pass through)

Output → `source_prs/unfiltered_SNP_reconciliation/`:
  <source>_unfiltered_reconciliation.tsv   per-published-SNP augmented row
  reconciliation_summary.txt               per-source totals + rescue candidates
  recoverable_snps.tsv                     flat roster of on-chip-but-filtered

Usage:
  conda run -n snp python snp_pipeline/30c_unfiltered_prs_reconciliation.py
"""
from __future__ import annotations

import argparse
import importlib.util
import math
import re
from pathlib import Path

import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
SRC_PRS = BASE / "source_prs"
UNFILT_BIM = BASE / "WGS_Omni25_BIN_wo_ConsentsIssues.bim"
OUT = SRC_PRS / "unfiltered_SNP_reconciliation"

# QC report files (from 01_qc.sh runs; pre-filter on 812 subjects)
LMISS = BASE / "results" / "missingness" / "missing.lmiss"
FRQ = BASE / "results" / "maf" / "maf.frq"
HWE = BASE / "results" / "hwe" / "hwe.hwe"

# Intermediate BIMs for per-stage membership detection
STAGE_BIMS = {
    "geno_0.02":    BASE / "SNP_filtered.bim",
    "hwe_1e-7":     BASE / "SNP_filtered_hwe.bim",
    "refcorr_step": BASE / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr.bim",
    "maf_0.01":     BASE / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr_maf001.bim",
}

ILLUMINA = BASE / "liftover" / "Illumina"
ILLUMINA_MANIFESTS = (
    ("b138", ILLUMINA / "HumanOmni2-5-8-v1-2-A-b138-rsIDs.txt"),
    ("b144", ILLUMINA / "InfiniumOmni2-5-8v1-3_A1_b144_rsids.txt"),
    ("b150", ILLUMINA / "InfiniumOmni2-5-8v1-4_A1_b150_rsids.txt"),
    ("b151", ILLUMINA / "InfiniumOmni2-5-8v1-5_A1_b151_rsids.txt"),
)

# RsMergeArch reused via importlib (script 30b's file starts with a digit)
SCRIPT_30B = Path(__file__).with_name("30b_published_prs_resolution.py")
_spec = importlib.util.spec_from_file_location("_s30b", SCRIPT_30B)
_s30b = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_s30b)
_build_rsmerge_map = _s30b._build_rsmerge_map

SOURCES = ("Bellenguez", "Wightman", "Kunkle", "Desikan")
_RSID_RE = re.compile(r"^rs\d+$")


def _read_resolution_tsv(src: str) -> pd.DataFrame:
    return pd.read_csv(SRC_PRS / f"{src}_full_snps_resolution.tsv",
                       sep="\t", dtype=str, keep_default_na=False)


def _load_unfilt_bim():
    """Return (set_of_names, {CHR:BP -> {name,A1,A2}}, by_name DataFrame).
    The unfiltered BIM `name` column carries Illumina kgp probe IDs and
    some rsIDs; we don't pre-separate them."""
    print(f"[unfilt] loading {UNFILT_BIM.name} …")
    bim = pd.read_csv(UNFILT_BIM, sep="\t", header=None,
                      names=["CHR", "name", "cM", "BP", "A1", "A2"],
                      dtype=str)
    n_total = len(bim)
    name_set = set(bim["name"])
    mapped = bim[(bim["CHR"] != "0") & (bim["BP"] != "0")]
    n_unmapped = n_total - len(mapped)
    pos_idx: dict[str, dict] = {}
    for chrom, name, bp, a1, a2 in mapped[["CHR", "name", "BP", "A1",
                                            "A2"]].itertuples(index=False):
        pos_idx[f"{chrom}:{bp}"] = {"name": name, "A1": a1, "A2": a2,
                                    "CHR": chrom, "BP": bp}
    print(f"  unfiltered: {n_total:,} rows ({n_unmapped:,} chr=0/bp=0); "
          f"|name set|={len(name_set):,}; |mapped pos|={len(pos_idx):,}")
    return name_set, pos_idx, bim.set_index("name")


def _load_illumina_rsid_to_kgp() -> dict[str, list[tuple[str, str]]]:
    """rsID → [(kgp_name, manifest_tag), …] union across b138/b144/b150/b151.
    Mirrors the parse used by 05a_update_snp_ids.method_b138 (L196-218).
    """
    rsid2kgp: dict[str, list[tuple[str, str]]] = {}
    for tag, path in ILLUMINA_MANIFESTS:
        if not path.exists():
            print(f"  [WARN] missing {path}")
            continue
        df = pd.read_csv(path, sep="\t", usecols=["Name", "RsID"],
                         dtype=str).dropna(subset=["Name", "RsID"])
        df["RsID"] = df["RsID"].str.strip()
        df = df[df["RsID"].str.match(_RSID_RE, na=False)]
        # filter to kgp probe names only (skip Illumina internal rsXXX-named probes)
        df = df[df["Name"].str.startswith("kgp", na=False)]
        for rs, kgp in zip(df["RsID"], df["Name"]):
            rsid2kgp.setdefault(rs, []).append((kgp, tag))
        print(f"  [illumina] {tag}: {len(df):,} kgp→rsID rows")
    print(f"  [illumina] rsID → kgp map: {len(rsid2kgp):,} unique rsIDs")
    return rsid2kgp


def _load_qc_reports() -> tuple[dict[str, float], dict[str, float], dict[str, float]]:
    """Load per-SNP F_MISS, MAF, HWE p-value (812-subject pre-filter) from
    the QC reports under `results/`. plink output is whitespace-separated."""
    def _read(path: Path, key_col: str, val_col: str) -> dict[str, float]:
        print(f"  [qc] {path.name}")
        df = pd.read_csv(path, delim_whitespace=True, dtype=str,
                         usecols=[key_col, val_col])
        df[val_col] = pd.to_numeric(df[val_col], errors="coerce")
        return dict(zip(df[key_col], df[val_col]))

    print(f"[qc] loading 812-subject QC reports …")
    f_miss = _read(LMISS, "SNP", "F_MISS")
    maf = _read(FRQ, "SNP", "MAF")
    hwe = _read(HWE, "SNP", "P")
    print(f"  [qc] F_MISS: {len(f_miss):,} | MAF: {len(maf):,} | HWE: {len(hwe):,}")
    return f_miss, maf, hwe


def _load_stage_membership() -> dict[str, set[str]]:
    """Return {stage_key: set_of_names} from each intermediate BIM. The
    union of the kgp and rsID columns is just the BIM's `name` column —
    one probe is either kgpXXX or rsXXX at each stage (not both)."""
    print(f"[stages] loading intermediate BIM membership …")
    out: dict[str, set[str]] = {}
    for tag, path in STAGE_BIMS.items():
        if not path.exists():
            print(f"  [WARN] missing {path}")
            out[tag] = set()
            continue
        df = pd.read_csv(path, sep="\t", header=None, usecols=[1],
                         names=["name"], dtype=str)
        out[tag] = set(df["name"])
        print(f"  [stages] {tag:14s}: {len(out[tag]):>10,} probes — {path.name}")
    return out


def _determine_filter_reason(kgp: str, rsid: str,
                              f_miss: dict[str, float],
                              maf: dict[str, float],
                              hwe: dict[str, float],
                              stages: dict[str, set[str]]
                              ) -> tuple[str, str, float, float, float]:
    """Walk QC stages in order; return (filter_reason, drop_stage,
    F_MISS, MAF, HWE_p) for the first stage where this probe drops.
    `kgp` and/or `rsid` may be blank — try both at each stage since the
    probe is named kgpXXX in early stages and rsXXX after step 5 rename."""

    def _in(stage: str) -> bool:
        s = stages.get(stage, set())
        return (bool(kgp) and kgp in s) or (bool(rsid) and rsid in s)

    # statistics keyed by kgp name (the early stages' BIM name)
    fm = f_miss.get(kgp, float("nan"))
    mf = maf.get(kgp, float("nan"))
    hp = hwe.get(kgp, float("nan"))

    # step 2: --geno 0.02
    if not _in("geno_0.02"):
        reason = (f"missingness>0.02 (F_MISS={fm:.4f})"
                  if not math.isnan(fm) else "missingness>0.02 (F_MISS=NA)")
        return reason, "geno_0.02", fm, mf, hp
    # step 3: --hwe 1e-7
    if not _in("hwe_1e-7"):
        reason = (f"HWE<1e-7 (p={hp:.2e})"
                  if not math.isnan(hp) else "HWE<1e-7 (p=NA)")
        return reason, "hwe_1e-7", fm, mf, hp
    # step 5: rsID rename / liftover / ref-correction
    if not _in("refcorr_step"):
        return "rename_or_liftover_or_refcorr", "refcorr_step", fm, mf, hp
    # step 6: --maf 0.01
    if not _in("maf_0.01"):
        if not math.isnan(mf):
            mono = ", monomorphic" if mf == 0 else ""
            reason = f"MAF<0.01 (MAF={mf:.4f}{mono})"
        else:
            reason = "MAF<0.01 (MAF=NA)"
        return reason, "maf_0.01", fm, mf, hp
    return "unknown_still_in_maf_bim", "unknown", fm, mf, hp


def _parse_or_beta(s: str) -> float | None:
    if not s:
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    if math.isnan(v):
        return None
    if 0.5 < v < 3.0:
        return abs(math.log(v)) if v > 0 else None
    return abs(v)


def _kgp_lookup(rsid: str,
                rsid2kgp: dict[str, list[tuple[str, str]]],
                name_set: set[str]) -> tuple[str, str]:
    """Return (kgp_name, manifest_tag) for the first match found in the
    unfiltered BIM, or ("", "") if none. If multiple kgp predecessors
    exist, the kgp present in the BIM wins; if multiple are present
    (rare), return the first; the manifest tag is the union of all the
    manifests that list that winning kgp→rsID."""
    if not rsid:
        return "", ""
    cands = rsid2kgp.get(rsid, [])
    if not cands:
        return "", ""
    # find a kgp that exists in the unfiltered BIM
    winning_kgp = None
    for kgp, _ in cands:
        if kgp in name_set:
            winning_kgp = kgp
            break
    if winning_kgp is None:
        return "", ""
    # collect every manifest where the winning kgp→rsID is recorded
    tags = sorted({tag for kgp, tag in cands if kgp == winning_kgp})
    return winning_kgp, "+".join(tags)


def _base_row(r, src, assembly, in_maf_bim) -> dict:
    return {
        "source": src,
        "rsID_pub": str(r["rsID_pub"]),
        "rsID_canonical": str(r["rsID_canonical"]),
        "original_assembly": assembly,
        "CHR_pub": str(r.get("CHR_pub", "")),
        "BP_pub": str(r.get("BP_pub", "")),
        "effect_allele_pub": str(r.get("effect_allele_pub", "")),
        "OR_or_beta_pub": str(r.get("OR_or_beta_pub", "")),
        "P_pub": str(r.get("P_pub", "")),
        "locus_name": str(r.get("locus_name", "")),
        "drop_reason_30b": str(r["drop_reason"]),
        "in_maf_bim": in_maf_bim,
    }


def _reconcile_one(df: pd.DataFrame, src: str,
                   name_set: set[str],
                   pos_idx: dict[str, dict],
                   by_name: pd.DataFrame,
                   rsid2kgp: dict[str, list[tuple[str, str]]],
                   f_miss: dict[str, float],
                   maf_dict: dict[str, float],
                   hwe_dict: dict[str, float],
                   stages: dict[str, set[str]]) -> pd.DataFrame:
    rows = []
    src_assembly = str(df.iloc[0]["original_assembly"]) if len(df) else ""
    can_use_pos = src_assembly == "GRCh37"

    if "rsID_canonical" not in df.columns:
        rsmap = _build_rsmerge_map(set(df["rsID_pub"]))
        df = df.assign(rsID_canonical=df["rsID_pub"].map(
            lambda r: rsmap.get(r, r)))

    for _, r in df.iterrows():
        rs_pub = str(r["rsID_pub"])
        rs_can = str(r["rsID_canonical"])
        chrom = str(r["CHR_pub"]).strip()
        bp_str = str(r["BP_pub"]).strip()
        drop30b = str(r["drop_reason"])
        in_maf_bim = (drop30b == "ok")

        if drop30b == "apoe_haplotype_special":
            rows.append({**_base_row(r, src, src_assembly, in_maf_bim),
                         "in_unfiltered_bim": False,
                         "unfilt_match_method": "none",
                         "unfilt_name": "",
                         "unfilt_kgp_name": "",
                         "kgp_manifest_source": "",
                         "unfilt_CHR_grch37": "",
                         "unfilt_BP_grch37": "",
                         "unfilt_A1": "", "unfilt_A2": "",
                         "reconciliation_category": "apoe_haplotype_special",
                         "filter_reason": "",
                         "drop_stage": "",
                         "F_MISS_812": float("nan"),
                         "MAF_812": float("nan"),
                         "HWE_p_812": float("nan")})
            continue

        method = "none"
        ufname = ufkgp = uftag = ""
        ufchr = ufbp = ufa1 = ufa2 = ""

        # 1) direct rsID
        if rs_pub and rs_pub.startswith("rs") and rs_pub in name_set:
            method = "rsID"
            ufname = rs_pub
        # 2) RsMergeArch canonical rsID
        elif rs_can and rs_can.startswith("rs") and rs_can in name_set:
            method = "rsID_merged"
            ufname = rs_can
        # 3) Illumina kgp reverse-map (try rs_pub then rs_can)
        else:
            for candidate in (rs_pub, rs_can):
                kgp, tag = _kgp_lookup(candidate, rsid2kgp, name_set)
                if kgp:
                    method = "kgp_renamed"
                    ufname = kgp
                    ufkgp, uftag = kgp, tag
                    break
            # 4) positional fallback (GRCh37 sources only)
            if method == "none" and can_use_pos and chrom and bp_str.isdigit():
                key = f"{chrom}:{bp_str}"
                hit = pos_idx.get(key)
                if hit is not None:
                    method = "position"
                    ufname, ufchr, ufbp = hit["name"], hit["CHR"], hit["BP"]
                    ufa1, ufa2 = hit["A1"], hit["A2"]
                    # if this positional hit is itself a kgp, record the
                    # manifest tag where it maps to the published rsID
                    if ufname.startswith("kgp"):
                        ufkgp = ufname
                        # find which manifest(s) name this kgp → rs_pub
                        for kgp, tag in rsid2kgp.get(rs_pub, []):
                            if kgp == ufname:
                                uftag = (uftag + "+" + tag).lstrip("+")

        # fill CHR/BP/A1/A2 from BIM for non-positional hits
        if method in ("rsID", "rsID_merged", "kgp_renamed") and ufname in by_name.index:
            br = by_name.loc[ufname]
            if isinstance(br, pd.DataFrame):
                br = br.iloc[0]
            ufchr = ufchr or str(br["CHR"])
            ufbp = ufbp or str(br["BP"])
            ufa1 = ufa1 or str(br["A1"])
            ufa2 = ufa2 or str(br["A2"])

        in_unfilt = method != "none"
        if in_maf_bim:
            cat = ("on_chip_passed_maf" if in_unfilt
                   else "in_maf_but_not_on_chip_manifest")
        elif method == "rsID":
            cat = "on_chip_filtered_out"
        elif method == "rsID_merged":
            cat = "on_chip_via_rsID_merge_filtered_out"
        elif method == "kgp_renamed":
            cat = "on_chip_via_kgp_renamed_filtered_out"
        elif method == "position":
            cat = "on_chip_via_position_filtered_out"
        else:
            cat = "not_on_chip"

        # filter-reason attribution for on-chip-but-not-in-MAF probes
        filter_reason = ""
        drop_stage = ""
        fm = mf = hp = float("nan")
        if in_unfilt and not in_maf_bim:
            # the probe IS on the chip — figure out where it dropped
            probe_kgp = ufkgp if ufkgp else (ufname if ufname.startswith("kgp") else "")
            probe_rs = ufname if ufname.startswith("rs") else rs_pub
            filter_reason, drop_stage, fm, mf, hp = _determine_filter_reason(
                probe_kgp, probe_rs, f_miss, maf_dict, hwe_dict, stages)

        rows.append({**_base_row(r, src, src_assembly, in_maf_bim),
                     "in_unfiltered_bim": in_unfilt,
                     "unfilt_match_method": method,
                     "unfilt_name": ufname,
                     "unfilt_kgp_name": ufkgp,
                     "kgp_manifest_source": uftag,
                     "unfilt_CHR_grch37": ufchr,
                     "unfilt_BP_grch37": ufbp,
                     "unfilt_A1": ufa1, "unfilt_A2": ufa2,
                     "reconciliation_category": cat,
                     "filter_reason": filter_reason,
                     "drop_stage": drop_stage,
                     "F_MISS_812": fm,
                     "MAF_812": mf,
                     "HWE_p_812": hp})

    return pd.DataFrame(rows)


def _summarise(per_source: dict[str, pd.DataFrame], n_unfilt: int) -> str:
    lines = [
        "Unfiltered-BIM PRS-SNP reconciliation summary",
        f"Reconciliation target: {UNFILT_BIM.name} ({n_unfilt:,} rows; GRCh37)",
        "Methods: direct rsID → RsMergeArch canonical → Illumina kgp manifest "
        "→ position (GRCh37 sources only)",
        "kgp→rsID reverse map built from "
        "liftover/Illumina/{HumanOmni2-5-8-v1-2-A-b138-rsIDs.txt, "
        "InfiniumOmni2-5-8v1-{3,4,5}_A1_b{144,150,151}_rsids.txt} — "
        "the same manifests used by 05a_update_snp_ids.py to rename "
        "kgp probes → rsIDs in the active pipeline.",
        "Bellenguez positional fallback skipped (published in GRCh38; no "
        "reverse liftover chain available — rsID/canonical/kgp matching only).",
        "Categories:",
        "  on_chip_passed_maf                       : in MAF>0.01 BIM AND on chip",
        "  on_chip_filtered_out                     : direct rsID in unfilt; lost to QC",
        "  on_chip_via_rsID_merge_filtered_out      : only via canonical rsID",
        "  on_chip_via_kgp_renamed_filtered_out     : only via Illumina kgp manifest",
        "  on_chip_via_position_filtered_out        : only via GRCh37 position",
        "  in_maf_but_not_on_chip_manifest          : in MAF but NOT in unfilt by any",
        "                                             route (genuine extract / imputation)",
        "  not_on_chip                              : absent from chip and not in MAF",
        "  apoe_haplotype_special                   : Desikan APOE haplotype (pass-through)",
        "",
    ]
    for src, df in per_source.items():
        n_pub = len(df)
        cats = df["reconciliation_category"].value_counts().to_dict()
        n_maf = cats.get("on_chip_passed_maf", 0)
        n_maf_nochip = cats.get("in_maf_but_not_on_chip_manifest", 0)
        n_filt_keys = (
            "on_chip_filtered_out",
            "on_chip_via_rsID_merge_filtered_out",
            "on_chip_via_kgp_renamed_filtered_out",
            "on_chip_via_position_filtered_out",
        )
        n_filt = sum(cats.get(k, 0) for k in n_filt_keys)
        n_nochip = cats.get("not_on_chip", 0)
        n_hap = cats.get("apoe_haplotype_special", 0)

        # how many of the on-MAF SNPs were matched via kgp rename?
        on_maf_by_kgp = int(((df["reconciliation_category"] == "on_chip_passed_maf")
                             & (df["unfilt_match_method"] == "kgp_renamed")).sum())
        on_maf_direct = int(((df["reconciliation_category"] == "on_chip_passed_maf")
                             & (df["unfilt_match_method"] == "rsID")).sum())

        lines.append(f"## {src}")
        lines.append(f"  published entries                  : {n_pub}")
        lines.append(f"  on chip & in MAF (ok in 30b)       : {n_maf}")
        if on_maf_direct:
            lines.append(f"    matched by direct rsID                : {on_maf_direct}")
        if on_maf_by_kgp:
            lines.append(f"    matched by Illumina kgp rename        : {on_maf_by_kgp}")
        if n_maf_nochip:
            lines.append(f"  in MAF but NOT on Illumina manifest: {n_maf_nochip}")
            lines.append(f"     (genuine extract / imputation)")
        lines.append(f"  on chip but filtered out           : {n_filt}")
        for sub_key, sub_label in (
            ("on_chip_filtered_out", "by direct rsID"),
            ("on_chip_via_rsID_merge_filtered_out", "by canonical rsID"),
            ("on_chip_via_kgp_renamed_filtered_out", "by Illumina kgp rename"),
            ("on_chip_via_position_filtered_out", "by GRCh37 position"),
        ):
            if cats.get(sub_key, 0):
                lines.append(f"    {sub_label:33s}: {cats[sub_key]}")
        lines.append(f"  not on Illumina manifest           : {n_nochip}")
        if n_hap:
            lines.append(f"  APOE haplotype (special)           : {n_hap}")
        if n_maf + n_maf_nochip + n_filt + n_nochip + n_hap != n_pub:
            lines.append(f"  ! WARNING categories sum to "
                         f"{n_maf + n_maf_nochip + n_filt + n_nochip + n_hap} ≠ {n_pub}")

        rescue = df[df["reconciliation_category"].isin(n_filt_keys)].copy()
        if len(rescue):
            rescue["abs_effect"] = rescue["OR_or_beta_pub"].map(_parse_or_beta)
            rescue = rescue.dropna(subset=["abs_effect"]).sort_values(
                "abs_effect", ascending=False)
            if len(rescue):
                lines.append(f"  top rescue candidates (sorted by |effect|):")
                for _, r in rescue.head(3).iterrows():
                    eff = r["OR_or_beta_pub"]
                    p = r["P_pub"]
                    locus = r["locus_name"] or "-"
                    kgp_tag = (f" [{r['unfilt_kgp_name']}, {r['kgp_manifest_source']}]"
                               if r["unfilt_kgp_name"] else "")
                    filt = (f"  filter={r['filter_reason']}"
                            if r.get("filter_reason") else "")
                    lines.append(
                        f"    {r['rsID_pub']:14s} "
                        f"OR/β={eff:>8s}  P={p:>10s}  "
                        f"chr{r['unfilt_CHR_grch37']}:"
                        f"{r['unfilt_BP_grch37']}  {locus}{kgp_tag}{filt}")
        lines.append("")

    lines.extend([
        "Caveats:",
        "  - 'on chip but filtered out' SNPs are attributed to the earliest QC",
        "    stage that drops them. Stages tested (per 01_qc.sh): --geno 0.02",
        "    (step 2) → --hwe 1e-7 (step 3) → rsID rename / liftover / ref-corr",
        "    (step 5) → --maf 0.01 (step 6). Per-SNP statistics (F_MISS, MAF,",
        "    HWE p) computed on the 812-subject pre-filter cohort — see the",
        "    `filter_reason`, `drop_stage`, `F_MISS_812`, `MAF_812`, `HWE_p_812`",
        "    columns of the per-source and recoverable TSVs.",
        "  - 'on_chip_via_kgp_renamed_*' uses the SAME Illumina kgp→rsID",
        "    manifests that 05a_update_snp_ids.py applies to rename probes",
        "    in the active pipeline; matching here is therefore self-consistent",
        "    with how the MAF BIM was built. Manifest provenance is in the",
        "    'kgp_manifest_source' TSV column (b138 / b144 / b150 / b151 or a",
        "    '+'-joined union).",
        "  - APOE causal variants:",
        "      rs429358 (ε4 SNP, chr19:45,411,941 GRCh37)  IS on the Omni2.5M",
        "        chip as kgp9680313 (b138/b144 manifests); the probe is",
        "        monomorphic (A1=0) in the consent-cleared cohort → dropped",
        "        at MAF QC. Categorised 'on_chip_via_kgp_renamed_filtered_out'.",
        "      rs7412   (ε2 SNP, chr19:45,412,079 GRCh37)  NOT on the chip at",
        "        all (no probe at this position; no kgp in any manifest).",
        "        Genuinely 'not_on_chip'.",
        "  - Bellenguez positional fallback is unavailable (GRCh38→GRCh37 chain",
        "    not in liftover/); however the kgp manifests are assembly-",
        "    independent so Bellenguez kgp matches still work.",
    ])
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    name_set, pos_idx, by_name = _load_unfilt_bim()
    n_unfilt = by_name.shape[0]
    print(f"\n[illumina] loading kgp→rsID manifests …")
    rsid2kgp = _load_illumina_rsid_to_kgp()
    f_miss, maf_dict, hwe_dict = _load_qc_reports()
    stages = _load_stage_membership()

    per_source: dict[str, pd.DataFrame] = {}
    for src in SOURCES:
        print(f"\n[{src}] reconciling …")
        df = _read_resolution_tsv(src)
        out_df = _reconcile_one(df, src, name_set, pos_idx, by_name, rsid2kgp,
                                f_miss, maf_dict, hwe_dict, stages)
        per_source[src] = out_df
        out_path = args.out / f"{src}_unfiltered_reconciliation.tsv"
        out_df.to_csv(out_path, sep="\t", index=False)
        cats = out_df["reconciliation_category"].value_counts().to_dict()
        n_filt = sum(cats.get(k, 0) for k in (
            "on_chip_filtered_out",
            "on_chip_via_rsID_merge_filtered_out",
            "on_chip_via_kgp_renamed_filtered_out",
            "on_chip_via_position_filtered_out"))
        print(f"  published={len(out_df)}  on_MAF&chip={cats.get('on_chip_passed_maf', 0)}  "
              f"in_MAF_not_on_chip={cats.get('in_maf_but_not_on_chip_manifest', 0)}  "
              f"filtered_out={n_filt}  not_on_chip={cats.get('not_on_chip', 0)}  "
              f"haplotype={cats.get('apoe_haplotype_special', 0)}")
        print(f"  [out] {out_path}")

    recoverable = pd.concat([
        df[df["reconciliation_category"].isin([
            "on_chip_filtered_out",
            "on_chip_via_rsID_merge_filtered_out",
            "on_chip_via_kgp_renamed_filtered_out",
            "on_chip_via_position_filtered_out"])][[
            "source", "rsID_pub", "rsID_canonical", "unfilt_name",
            "unfilt_kgp_name", "kgp_manifest_source",
            "unfilt_CHR_grch37", "unfilt_BP_grch37", "unfilt_A1", "unfilt_A2",
            "effect_allele_pub", "OR_or_beta_pub", "P_pub", "locus_name",
            "unfilt_match_method", "reconciliation_category",
            "filter_reason", "drop_stage",
            "F_MISS_812", "MAF_812", "HWE_p_812"]]
        for df in per_source.values()
    ], ignore_index=True)
    rec_path = args.out / "recoverable_snps.tsv"
    recoverable.to_csv(rec_path, sep="\t", index=False)
    print(f"\n[out] {rec_path}  (n={len(recoverable)})")

    summ = _summarise(per_source, n_unfilt)
    summ_path = args.out / "reconciliation_summary.txt"
    summ_path.write_text(summ, "utf-8")
    print(f"[out] {summ_path}")


if __name__ == "__main__":
    main()
