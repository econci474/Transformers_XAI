"""
20_genetic_baselines.py
=======================
Classical genetic baselines for ADNI `ever_convert`, to contextualise the
exhaustive BMFM/NT/Caduceus negative (test AUC ~0.5). If APOE / a β-weighted
PRS clears chance on the IDENTICAL splits while the foundation models do not,
the FM negative is a representation/extraction failure, not absent signal.

Baselines (all = sklearn LogisticRegression(class_weight="balanced") on a
feature vector; fit train∪val, eval test; seeds 0/1/2; metrics.json schema
== script 15 `evaluate`):

  * betaPRS_z = Σ β_A1 · dosage, z-standardised on TRAIN before the logit,
    at 3 harmonised APOE tiers: withAPOE (all 430) / onlyAPOE (166, isolated
    chr19 APOE LD block) / noAPOE (263, APOE region removed, rs2072561 kept).
  * Logistic regression on the raw dosage matrix (same SNP subsets).
  * APOE-ε4 single-SNP reference (rs429358 dosage) — the positive control.
  * Covariate grid: PRS / APOE × {–, +AGE, +SEX, +AGE+SEX} + demographic-only.

β-harmonisation: re-parse β from the two raw GWAS files (Wightman `beta`;
Bellenguez `log(OR)`), pick the per-SNP primary study from
`external_gwas_labels.z_source`, cross-check the effect allele vs that file's
`EA`, and orient to the dosage-counted allele (BIM A1) by reusing the
already-validated `allele_flip` column:  β_A1 = −β_raw if allele_flip else β_raw
(dosage = #A1 copies; A2 = REF — script 12).

PRS standardisation (reference-population SD) → `prs_standardisation.tsv`:
empirical (full 616-cohort mean/SD) + theoretical discovery SD
√(Σ β_A1²·2p(1−p)) from GWAS effect-allele freq (ADNI-freq fallback when
the Catalog says NR); OR per +1 SD reported for 1-feature PRS models.

Usage
-----
  conda run -n snp python snp_pipeline/20_genetic_baselines.py --sanity
  conda run -n snp python snp_pipeline/20_genetic_baselines.py --seeds 0 1 2
"""
from __future__ import annotations

import argparse
import importlib.util as _il
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_auc_score, balanced_accuracy_score, f1_score,
    precision_score, recall_score, confusion_matrix,
)

warnings.filterwarnings("ignore")

# ── Reuse script 15 load_labels (digit-leading filename) ────────────────────
_S15 = Path(__file__).parent / "15_train_classification_head.py"
_spec = _il.spec_from_file_location("script15", _S15)
_s15 = _il.module_from_spec(_spec)
_spec.loader.exec_module(_s15)
load_labels = _s15.load_labels

# ── Paths ───────────────────────────────────────────────────────────────────
DATA       = Path("D:/ADNI_SNP_Omni2.5M_20140220")
DOSAGE_TSV = DATA / "bmfm_inputs/patient_genotypes/patient_dosage.tsv"
DOSAGE_BIM = DATA / "bmfm_inputs/patient_genotypes/patient_dosage.bim"
LABELS_TSV = DATA / "bmfm_inputs/external_gwas_labels.tsv"
TOPK25_TSV = DATA / "classification_finetune_topk/topk25_snps.tsv"
WIGHTMAN   = DATA / "GWAS/Wightman_2021/wightman_without_ukb_GRCh38.tsv"
BELLENGUEZ = DATA / "GWAS/Bellenguez_2022/gwas-association-downloaded_2026-04-28-pubmedId_35379992.tsv"
SPLITS = {
    "ever_convert": Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                          "no_cdr_stratified_ever_convert/tabular/baseline"),
    "bl_multi":     Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                          "no_cdr_stratified/tabular/baseline"),
}
OUT_ROOT_DEFAULT = DATA / "genetic_baselines"

E_SNPS_DEFAULT = ["rs429358", "rs7412"]          # APOE-ε defining SNPs
APOE_CLUSTER_DEFAULT = ("19", 44_400_000, 46_500_000)   # GRCh38 broad APOE region
# rs2072561 (chr19:46,022,609) is the single physically-isolated SNP in the
# broad window — a 241 kb gap from the next-nearest, ~600 kb past the bulk —
# so it is NOT part of the APOE LD block: it is EXCLUDED from `onlyAPOE`
# and RETAINED in `noAPOE` (treated as non-APOE polygenic background,
# alongside rs429358 which is monomorphic in the array).
APOE_ISOLATED_DEFAULT = ["rs2072561"]
# Harmonised PRS tiers (same names used by script 23 Cox):
#   withAPOE  all 430 GW-sig SNPs (rs429358 monomorphic & rs7412 absent are
#             inert, so this == "all minus the ε-SNPs" numerically)
#   onlyAPOE  166 = chr19:44.4–46.5 Mb window − rs2072561 − rs429358
#             (the isolated APOE LD block only)
#   noAPOE    263 = the 430 with the APOE region removed (rs2072561 kept as
#             non-APOE background)
APOE_TIERS = ["withAPOE", "onlyAPOE", "noAPOE"]

# Strict conversion endpoints from script 22 conversion_labels.tsv:
# mode → (at-risk baseline-dx set, label column). Active only when
# --label-mode is one of these AND --labels-tsv is given (else byte-identical
# to the legacy load_labels path). Split partitions are unchanged: membership
# comes from the existing ever_convert seed splits, joined by Patient_ID.
STRICT_MODES = {
    "ever_conversion_AD":        ({"CN", "MCI"}, "ever_conversion_AD"),
    "ever_conversion_AD_or_MCI": ({"CN", "MCI"}, "ever_conversion_AD_or_MCI"),
    "ever_conversion_MCI":       ({"CN"},        "ever_conversion_MCI"),
    "progressive_mci":           ({"MCI"},       "progressive_MCI"),
    # classification case/control: prevalent baseline-AD kept as positives
    # (at-risk = everyone) alongside the incident converters.
    "ad_case":        ({"CN", "MCI", "AD"}, "ad_case"),
    "ad_or_mci_case": ({"CN", "MCI", "AD"}, "ad_or_mci_case"),
}


# ── Metric dict — mirrors snp_pipeline/15 evaluate() exactly ────────────────
def _metrics(y_true: np.ndarray, proba: np.ndarray) -> dict:
    y = np.asarray(y_true).astype(int)
    pred = (proba > 0.5).astype(int)
    out = {
        "auc": float(roc_auc_score(y, proba)) if len(np.unique(y)) > 1 else float("nan"),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "f1": float(f1_score(y, pred, zero_division=0)),
        "precision": float(precision_score(y, pred, pos_label=1, zero_division=0)),
        "recall": float(recall_score(y, pred, pos_label=1, zero_division=0)),
        "sensitivity": float(recall_score(y, pred, pos_label=1, zero_division=0)),
        "specificity": float(recall_score(y, pred, pos_label=0, zero_division=0)),
    }
    cm = confusion_matrix(y, pred, labels=[0, 1])
    out["confusion_matrix"] = {"tn": int(cm[0, 0]), "fp": int(cm[0, 1]),
                               "fn": int(cm[1, 0]), "tp": int(cm[1, 1])}
    return out


# ── Parse raw GWAS effect sizes ─────────────────────────────────────────────
def _load_wightman_allowlist(path) -> set | None:
    """Resolved-lead allowlist from script 24 (cols pipeline_rsID, lead_rsID,
    lead_rsID_canonical). Returns a lowercased rsID set, or None (no filter)."""
    if not path:
        return None
    a = pd.read_csv(path, sep="\t")
    ids = set()
    for c in ("pipeline_rsID", "lead_rsID", "lead_rsID_canonical"):
        if c in a.columns:
            ids |= {str(x).strip().lower() for x in a[c].dropna()
                    if str(x).strip()}
    return ids


def parse_wightman(allowlist: set | None = None) -> dict:
    """rsID & rsID_canonical → (effect_allele, beta, eaf). If `allowlist`
    given, restrict to those rsIDs (script-24 resolved Wightman leads) —
    default None ⇒ byte-identical to the unfiltered behaviour."""
    w = pd.read_csv(WIGHTMAN, sep="\t", low_memory=False)
    w = w[w["lifted"].astype(str).str.lower().isin(["true", "1"])]
    out: dict[str, tuple] = {}
    for r in w.itertuples(index=False):
        try:
            rec = (str(r.effect_allele).upper().strip(),
                   float(r.beta),
                   float(r.effect_allele_frequency))
        except (ValueError, TypeError):
            continue
        for key in (r.rsID, r.rsID_canonical):
            k = str(key).strip()
            if k and k.lower() not in ("nan", "none"):
                if allowlist is not None and k.lower() not in allowlist:
                    continue
                out.setdefault(k, rec)          # first wins (canonical dups rare)
    return out


def parse_bellenguez() -> dict:
    """SNPS → (risk_allele, log_or, raf).  Vectorised: GWAS-Catalog column
    names contain spaces (itertuples attr access would mangle them)."""
    b = pd.read_csv(BELLENGUEZ, sep="\t", low_memory=False)
    out: dict[str, tuple] = {}
    sa = b["STRONGEST SNP-RISK ALLELE"].astype(str)
    risk = sa.str.rsplit("-", n=1).str[-1].str.upper().str.strip()
    orb = pd.to_numeric(b["OR or BETA"], errors="coerce")
    raf = pd.to_numeric(b["RISK ALLELE FREQUENCY"], errors="coerce")
    for rsid, ra, o, f in zip(b["SNPS"].astype(str).str.strip(), risk, orb, raf):
        if not rsid or rsid.lower() in ("nan", "none") or not np.isfinite(o) or o <= 0:
            continue
        out.setdefault(rsid, (ra, math.log(float(o)), float(f) if np.isfinite(f) else float("nan")))
    return out


# ── β-harmonisation table for the 430 GW-sig SNPs ───────────────────────────
def build_beta_table(verbose: bool = False,
                     wightman_allowlist: set | None = None) -> pd.DataFrame:
    lab = pd.read_csv(LABELS_TSV, sep="\t", low_memory=False)
    lab = lab[lab["label"] == 1].copy()
    lab["SNP"] = lab["SNP"].astype(str)
    lab = lab.drop_duplicates("SNP").set_index("SNP")
    bim = pd.read_csv(DOSAGE_BIM, sep="\t")
    bim["rsID"] = bim["rsID"].astype(str)
    bim = bim.drop_duplicates("rsID").set_index("rsID")        # CHR rsID BP A1 A2
    snps = [s for s in bim.index if s in lab.index]            # the ~430 genotyped

    wig = parse_wightman(wightman_allowlist)
    bel = parse_bellenguez()
    rows = []
    n_drop = 0
    for s in snps:
        L = lab.loc[s]
        zsrc = str(L.get("z_source", "")) or str(L.get("source", ""))
        ea_ext = str(L["EA"]).upper().strip()
        flip = str(L["allele_flip"]).strip().lower() in ("true", "1")
        z_ext = float(L["z_score"]) if pd.notna(L["z_score"]) else float("nan")
        if "Bellenguez" in zsrc:
            rec, src = bel.get(s), "Bellenguez2022_logOR"
        elif "Wightman" in zsrc:
            rec, src = wig.get(s), "Wightman2021"
        else:                                                  # defensive
            rec, src = (wig.get(s) or bel.get(s)), zsrc or "unknown"
        if rec is None:
            n_drop += 1
            if verbose:
                print(f"  [drop] {s}: no β in {src} raw file")
            continue
        ea_raw, eff, eaf = rec[0], float(rec[1]), float(rec[2])
        if ea_raw != ea_ext:
            n_drop += 1
            if verbose:
                print(f"  [drop] {s}: EA mismatch raw={ea_raw} vs external={ea_ext}")
            continue
        beta_a1 = -eff if flip else eff
        # GWAS EAF is freq of the GWAS effect allele (== external EA).
        # A1-allele freq: if EA==A1 (no flip) → eaf; if flip (EA==A2) → 1−eaf.
        if np.isfinite(eaf):
            eaf_a1, eaf_src = ((1.0 - eaf) if flip else eaf), src
        else:
            eaf_a1, eaf_src = float("nan"), "missing"
        rows.append({
            "SNP": s, "CHR": str(bim.loc[s, "CHR"]), "BP": int(bim.loc[s, "BP"]),
            "A1": str(bim.loc[s, "A1"]).upper(), "A2": str(bim.loc[s, "A2"]).upper(),
            "source": src, "ea_raw": ea_raw, "ea_external": ea_ext,
            "beta_raw": eff, "allele_flip": flip, "beta_A1": beta_a1,
            "z_A1": z_ext, "eaf_A1": eaf_a1, "eaf_source": eaf_src,
        })
    bt = pd.DataFrame(rows).set_index("SNP")
    bt.attrs["n_dropped"] = n_drop
    bt.attrs["n_candidate"] = len(snps)
    return bt


# ── APOE-tier SNP subsetting ────────────────────────────────────────────────
def snps_for_tier(bt: pd.DataFrame, base_snps: list[str], tier: str,
                  e_snps: list[str], cluster: tuple,
                  isolated: list[str] | None = None) -> list[str]:
    """Harmonised APOE tiers (same names in script 23 Cox):
      withAPOE   all SNPs (rs429358 monomorphic & rs7412 absent are inert)
      onlyAPOE   the isolated APOE LD block only: in the broad window, minus
                 the physically-isolated SNP(s) and the (monomorphic) e_snps
      noAPOE     drop the APOE region, but KEEP the physically-isolated
                 SNP(s) (`isolated`) as non-APOE polygenic background
    """
    s = [x for x in base_snps if x in bt.index]
    iso = set(APOE_ISOLATED_DEFAULT if isolated is None else isolated)
    if tier == "withAPOE":
        return s
    c, lo, hi = cluster

    def _in_apoe(x: str) -> bool:           # in window AND not the isolated SNP
        return (str(bt.loc[x, "CHR"]) == str(c)
                and lo <= int(bt.loc[x, "BP"]) <= hi
                and x not in iso)
    if tier == "noAPOE":                    # drop APOE region; keep isolated
        return [x for x in s if not _in_apoe(x)]
    if tier == "onlyAPOE":                  # isolated APOE LD block only
        return [x for x in s if _in_apoe(x) and x not in set(e_snps)]
    raise ValueError(tier)


# ── Feature builders ────────────────────────────────────────────────────────
def _covars(split_df: pd.DataFrame) -> pd.DataFrame:
    """Split-derived features, index = Patient_ID:
      AGE   = year(Date) − YOB
      SEX   = Female→1 / Male→0
      APOE4 = clinical ε4 allele count 0/1/2 (APOE4_Dosage, or '4' count
              in APOE_Alleles). Clinical source because the Omni2.5M array
              does NOT genotype rs429358 (monomorphic A1='.').
      APOE2 = clinical ε2 allele count 0/1/2 ('2' count in APOE_Alleles)
              — protective allele, modelled alongside ε4.
    """
    d = split_df.copy()
    d["Patient_ID"] = d["Patient_ID"].astype(str)
    yr = pd.to_datetime(d["Date"], errors="coerce").dt.year
    age = yr - pd.to_numeric(d["YOB"], errors="coerce")
    sex = (d["Sex"].astype(str).str.strip().str.lower()
           .map({"female": 1, "f": 1, "2": 1, "male": 0, "m": 0, "1": 0}))
    apoe4 = pd.to_numeric(d.get("APOE4_Dosage"), errors="coerce")
    al = d.get("APOE_Alleles", pd.Series([""] * len(d))).astype(str)
    apoe2 = al.str.count("2").where(al.str.contains(r"\d", na=False))
    apoe4 = apoe4.where(apoe4.notna(),
                        al.str.count("4").where(al.str.contains(r"\d",
                                                                na=False)))
    return pd.DataFrame({"AGE": age.values, "SEX": sex.values,
                         "APOE4": apoe4.values, "APOE2": apoe2.values},
                        index=d["Patient_ID"].values)


def _strict_labels(split_root: Path, conv_tsv: Path, mode: str,
                    seeds) -> dict:
    """Per-seed {train/val/test → DataFrame[Patient_ID, y]} drop-in for
    load_labels, sourced from script-22 conversion_labels.tsv. Partitions
    (train/val/test membership) are taken verbatim from the existing seed
    splits; only patients in the transition's at-risk baseline-dx set are
    kept (e.g. baseline-AD excluded from AD endpoints). y = strict label."""
    atrisk, ycol = STRICT_MODES[mode]
    cv = pd.read_csv(conv_tsv, sep="\t")
    cv["Patient_ID"] = cv["Patient_ID"].astype(str)
    cv = cv[cv["bl_dx"].astype(str).str.upper().isin(atrisk)]
    ymap = dict(zip(cv["Patient_ID"], cv[ycol].astype(int)))
    out: dict = {}
    for sd in seeds:
        out[sd] = {}
        for sp in ("train", "val", "test"):
            ids = (pd.read_csv(split_root / f"seed_{sd}" / f"{sp}.csv",
                               usecols=["Patient_ID"])["Patient_ID"]
                   .astype(str))
            keep = [p for p in ids if p in ymap]
            out[sd][sp] = pd.DataFrame({"Patient_ID": keep,
                                        "y": [ymap[p] for p in keep]})
    return out


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--label-mode", default="ever_convert",
                    choices=["ever_convert", "bl_multi",
                             *STRICT_MODES.keys()])
    ap.add_argument("--labels-tsv", type=Path, default=None,
                    help="script-22 conversion_labels.tsv. Required for the "
                         "strict ever_conversion_* / progressive_mci modes; "
                         "default None ⇒ legacy load_labels (byte-identical "
                         "for ever_convert / bl_multi).")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--out-root", type=Path, default=OUT_ROOT_DEFAULT)
    ap.add_argument("--apoe-esnps", nargs="+", default=E_SNPS_DEFAULT)
    ap.add_argument("--apoe-cluster-region", default="chr19:44400000-46500000",
                    help="GRCh38 APOE region (defines onlyAPOE / excluded by "
                         "noAPOE).")
    ap.add_argument("--skip-if-exists", action="store_true")
    ap.add_argument("--sanity", action="store_true",
                    help="Print β-harmonisation + APOE-ε4 assertion, then exit.")
    # Enrichment (all default None ⇒ byte-identical to the base run)
    ap.add_argument("--enriched-weights", type=Path, default=None,
                    help="enriched_snp_weights.tsv from script 21 (novel loci).")
    ap.add_argument("--enriched-dosage", type=Path, default=None,
                    help="patient_genotypes_enriched/patient_dosage.tsv (439).")
    ap.add_argument("--desikan-pgs", type=Path, default=None,
                    help="desikan_pgs_resolved.tsv from script 21.")
    ap.add_argument("--enriched-label", default="enriched",
                    help="names the injected-PRS feature set "
                         "betaPRS_z_<label>__<tier> (e.g. 9W27B / "
                         "9W27B+5Desikan) so the table is self-documenting.")
    ap.add_argument("--wightman-filtered", type=Path, default=None,
                    help="script-24 wightman_resolved_allowlist.tsv. When "
                         "set: restrict Wightman β to the resolved leads and "
                         "use only the noAPOE PRS tier. Default None ⇒ "
                         "byte-identical to the unfiltered run. (Desikan "
                         "coverage report is produced by script 24.)")
    args = ap.parse_args()

    wig_allow = _load_wightman_allowlist(args.wightman_filtered)

    _c = args.apoe_cluster_region.replace("chr", "").split(":")
    _lo, _hi = (int(x) for x in _c[1].split("-"))
    cluster = (_c[0], _lo, _hi)

    # Filtered run: Wightman→resolved leads only ⇒ onlyAPOE collapses to the
    # inert rs429358 and withAPOE≈noAPOE, so use the single noAPOE tier.
    tiers = ["noAPOE"] if wig_allow is not None else list(APOE_TIERS)

    # Strict modes reuse the ever_convert seed splits for membership (same
    # 616-patient universe / partitions); y comes from conversion_labels.tsv.
    split_root = SPLITS.get(args.label_mode) or SPLITS["ever_convert"]
    if args.label_mode in STRICT_MODES and args.labels_tsv is None:
        raise SystemExit(f"--label-mode {args.label_mode} requires "
                         f"--labels-tsv (script 22 conversion_labels.tsv)")

    print("[β] building harmonised β table …", flush=True)
    bt = build_beta_table(verbose=True, wightman_allowlist=wig_allow)
    n_drop, n_cand = bt.attrs["n_dropped"], bt.attrs["n_candidate"]
    all_snps = list(bt.index)
    print(f"[β] {len(all_snps)}/{n_cand} SNPs harmonised "
          f"({n_drop} dropped)", flush=True)

    # ── PRS reference standardisation provenance ────────────────────────────
    dos = pd.read_csv(DOSAGE_TSV, sep="\t")
    dos["PTID"] = dos["PTID"].astype(str)
    dos = dos.set_index("PTID")
    prov_rows = []
    for s in all_snps:
        b1, p = float(bt.loc[s, "beta_A1"]), bt.loc[s, "eaf_A1"]
        esrc = bt.loc[s, "eaf_source"]
        if not np.isfinite(p):                       # ADNI-freq fallback
            col = dos[s] if s in dos.columns else None
            p = (float(np.nanmean(col)) / 2.0) if col is not None else float("nan")
            esrc = "adni_fallback"
        prov_rows.append({"SNP": s, "beta_A1": b1, "eaf_A1": p,
                          "eaf_source": esrc, "CHR": bt.loc[s, "CHR"],
                          "BP": bt.loc[s, "BP"]})
    prov = pd.DataFrame(prov_rows)

    def _prs(snp_set, weight_col):
        w = bt.loc[snp_set, weight_col].astype(float)
        cols = [s for s in snp_set if s in dos.columns]
        D = dos[cols].copy()
        D = D.apply(lambda c: c.fillna(c.mean()))     # cohort-mean fill (ref)
        return D.mul(w.reindex(cols).values, axis=1).sum(axis=1)

    args.out_root.mkdir(parents=True, exist_ok=True)
    # Cohort-reference provenance only (NOT used for the reported OR — that
    # uses the per-seed empirical TRAIN PRS SD, see the fit loop). betaPRS_z
    # = the β-weighted PRS that is z-standardised before the logistic fit.
    std_rows = []
    for tier in tiers:
        ss = snps_for_tier(bt, all_snps, tier, args.apoe_esnps, cluster)
        if not ss:
            continue
        prs_all = _prs(ss, "beta_A1")
        emp_mu, emp_sd = float(prs_all.mean()), float(prs_all.std(ddof=0))
        pr = prov[prov["SNP"].isin(ss)]
        theo = float(np.sqrt(np.nansum(
            (pr["beta_A1"].values ** 2) * 2.0
            * pr["eaf_A1"].values * (1.0 - pr["eaf_A1"].values))))
        std_rows.append({"weight": "betaPRS_z", "snp_set": "all",
                         "apoe_tier": tier, "n_snps": len(ss),
                         "cohort_mean": emp_mu, "empirical_sd": emp_sd,
                         "theoretical_sd": theo})
    pd.DataFrame(std_rows).to_csv(args.out_root / "prs_standardisation.tsv",
                                  sep="\t", index=False)
    prov.to_csv(args.out_root / "prs_snp_weights.tsv", sep="\t", index=False)

    # ── Sanity gate (rs429358-independent) ──────────────────────────────────
    if args.sanity:
        print("\n=== SANITY: APOE-tier SNP counts ===")
        for tier in tiers:
            ss = snps_for_tier(bt, all_snps, tier, args.apoe_esnps, cluster)
            print(f"  {tier:13s} n={len(ss)}")
        iso = [x for x in APOE_ISOLATED_DEFAULT if x in bt.index]
        print(f"  isolated (kept in noAPOE, dropped from onlyAPOE): {iso}")
        print("\n=== SANITY: β-harmonisation sample (first 20) ===")
        cols = ["source", "ea_raw", "ea_external", "beta_raw",
                "allele_flip", "beta_A1", "z_A1", "A1", "A2"]
        with pd.option_context("display.width", 220,
                               "display.max_columns", None):
            print(bt[cols].head(20).to_string())
        mono = [s for s in all_snps
                if str(bt.loc[s, "A1"]).strip() in (".", "0")]
        print(f"\n[monomorphic in ADNI array] {len(mono)} SNP(s): {mono}")
        if "rs429358" in mono:
            print("  → rs429358 (APOE-ε4) NOT genotyped by Omni2.5M (dosage "
                  "constant). APOE control uses clinical APOE4_Dosage instead.")
        # Correctness gate: β_A1 and z_A1 are both oriented to BIM A1 via the
        # SAME allele_flip ⇒ they MUST share sign. Study-agnostic, all SNPs.
        m = bt[(bt["beta_A1"].abs() > 1e-12) & (bt["z_A1"].abs() > 1e-9)]
        disc = m[np.sign(m["beta_A1"]) != np.sign(m["z_A1"])]
        print(f"\n[gate] β_A1 vs z_A1 sign concordance: "
              f"{len(m)-len(disc)}/{len(m)} agree ({len(disc)} discordant)")
        if len(disc):
            print(disc[["source", "beta_A1", "z_A1"]].head(20).to_string())
        assert len(disc) == 0, ("β/z sign discordance ⇒ harmonisation bug "
                                 "(allele_flip or OR/logOR sign).")
        assert n_drop <= 12, f"too many SNPs dropped ({n_drop})"
        sd = pd.DataFrame(std_rows)
        assert (sd["theoretical_sd"] > 0).all() and (sd["empirical_sd"] > 0).all()
        print(f"\n[sanity] OK — harmonisation validated ({len(m)} informative "
              f"SNPs β/z-sign-concordant); {len(all_snps)} SNPs, {n_drop} "
              f"dropped; prs_standardisation.tsv written. (no models fitted)")
        return

    # ── Run registry: (model, fset, prs_spec|None, split_cols) ──────────────
    #   prs_spec = ("beta"|"dosage", tier) | ("precomp", name) | None
    #   split_cols ⊆ {"AGE","SEX","APOE4"} (pulled from the split CSV)
    #   betaPRS_z = β-weighted PRS, z-standardised on TRAIN before the logit.
    COVAR = {"+age": ["AGE"], "+sex": ["SEX"], "+age+sex": ["AGE", "SEX"]}
    # covariate-combo family applied to EACH filtered PRS set (6W27B /
    # 9W27B / 9W27B5D) so every set is reported in isolation and adjusted.
    PRS_COMBOS = [("", "prs", []),
                  ("+age", "prs_covar", ["AGE"]),
                  ("+sex", "prs_covar", ["SEX"]),
                  ("+age+sex", "prs_covar", ["AGE", "SEX"]),
                  ("+APOE4+age+sex", "prs_covar", ["APOE4", "AGE", "SEX"]),
                  ("+APOE2+APOE4+age+sex", "prs_covar",
                   ["APOE2", "APOE4", "AGE", "SEX"])]
    RUNS = []
    if wig_allow is None:
        # ── unfiltered: tier-based betaPRS_z / dosage / prs_covar ──────────
        for tier in tiers:
            RUNS.append(("prs", f"betaPRS_z_{tier}", ("beta", tier), []))
            RUNS.append(("logreg_dosage", f"dosage_{tier}",
                         ("dosage", tier), []))
            for cn, cc in COVAR.items():
                RUNS.append(("prs_covar", f"betaPRS_z_{tier}{cn}",
                             ("beta", tier), cc))
        RUNS.append(("prs_covar", "betaPRS_z_noAPOE+APOE4+age+sex",
                     ("beta", "noAPOE"), ["APOE4", "AGE", "SEX"]))
        RUNS.append(("prs_covar", "betaPRS_z_noAPOE+APOE2+APOE4+age+sex",
                     ("beta", "noAPOE"), ["APOE2", "APOE4", "AGE", "SEX"]))
        if "onlyAPOE" in tiers:
            RUNS.append(("prs_covar", "betaPRS_z_onlyAPOE+APOE4+age+sex",
                         ("beta", "onlyAPOE"), ["APOE4", "AGE", "SEX"]))
    # (filtered: betaPRS_z_{6W27B,9W27B,9W27B5D} families are registered in
    #  the enrichment block below.)
    RUNS.append(("apoe_e4", "APOE4_dosage", None, ["APOE4"]))
    RUNS.append(("apoe_e2", "APOE2_dosage", None, ["APOE2"]))
    RUNS.append(("apoe_e2e4", "APOE2+APOE4_dosage", None,
                 ["APOE2", "APOE4"]))
    for cn, cc in COVAR.items():
        RUNS.append(("apoe_e2", f"APOE2_dosage{cn}", None, ["APOE2"] + cc))
        RUNS.append(("apoe_e2e4", f"APOE2+APOE4_dosage{cn}", None,
                     ["APOE2", "APOE4"] + cc))
    for cn, cc in COVAR.items():
        RUNS.append(("apoe_covar", f"APOE4_dosage{cn}", None, ["APOE4"] + cc))
    for cn, cc in COVAR.items():
        RUNS.append(("demog", cn.lstrip("+"), None, cc))

    # ── Enrichment (script-21 artefacts): betaPRS_enriched + Desikan PHS ─────
    _precomp: dict = {}     # name → (PTID-indexed PRS Series, n_snps)
    if args.enriched_weights and args.enriched_dosage:
        ew = pd.read_csv(args.enriched_weights, sep="\t", low_memory=False)
        edos = pd.read_csv(args.enriched_dosage, sep="\t")
        edos["PTID"] = edos["PTID"].astype(str)
        edos = edos.set_index("PTID").apply(lambda c: c.fillna(c.mean()))
        # combined weight+pos table = 430 (from bt) ∪ novel (from script 21)
        comb = bt[["beta_A1", "CHR", "BP"]].copy()
        comb["CHR"] = comb["CHR"].astype(str)
        nv = ew[ew["novel"] == True].dropna(subset=["beta_A1"])  # noqa: E712
        for r in nv.itertuples(index=False):
            comb.loc[str(r.rsID)] = {"beta_A1": float(r.beta_A1),
                                     "CHR": str(r.CHR).split(".")[0],
                                     "BP": int(float(r.BP))}
        comb["BP"] = comb["BP"].astype(int)

        def _eprs(snp_list):
            cc2 = [s for s in snp_list if s in edos.columns]
            w = comb.loc[cc2, "beta_A1"].astype(float)
            return edos[cc2].mul(w.values, axis=1).sum(axis=1), len(cc2)

        # Named PRS sets, each reported in isolation AND with the full
        # covariate-combo family (PRS_COMBOS). base = the allowlist-filtered
        # β-table (6 in-430 Wightman + 27 Bellenguez); +W = the 3 re-extract
        # Wightman leads; +D = the 5 novel Desikan loci.
        base_rsids = list(bt.index)
        src = ew.get("source", pd.Series([""] * len(ew))).astype(str)
        w_novel = [str(r) for r, sc, nv in
                   zip(ew["rsID"], src, ew["novel"])
                   if nv and "wightman" in sc.lower()]
        d_novel = [str(r) for r, sc, nv in
                   zip(ew["rsID"], src, ew["novel"])
                   if nv and "desikan" in sc.lower()]
        prs_sets = {"6W27B": base_rsids,
                    "9W27B": base_rsids + w_novel}
        if d_novel:
            prs_sets["9W27B5D"] = base_rsids + w_novel + d_novel
        for S, rsids in prs_sets.items():
            for tier in tiers:                       # filtered ⇒ noAPOE only
                ss = snps_for_tier(comb, rsids, tier,
                                   args.apoe_esnps, cluster)
                s, n = _eprs(ss)
                key = f"betaPRS_z_{S}"
                _precomp[key] = (s, n)
                for suf, model, cc in PRS_COMBOS:
                    RUNS.append((model, f"{key}{suf}",
                                 ("precomp", key), cc))
                # raw-dosage LogReg on the same set's SNPs
                dcol = [x for x in ss if x in edos.columns]
                _precomp[f"dos_{S}"] = (edos[dcol], len(dcol))
                RUNS.append(("logreg_dosage", f"dosage_{S}",
                             ("precompmat", f"dos_{S}"), []))

        if args.desikan_pgs and Path(args.desikan_pgs).exists():
            dpg = pd.read_csv(args.desikan_pgs, sep="\t")
            dcols = [s for s in dpg["rsID"].astype(str) if s in edos.columns]
            dw = dpg.set_index("rsID").loc[dcols, "beta_A1"].astype(float)
            phs = edos[dcols].mul(dw.values, axis=1).sum(axis=1)
            _precomp["Desikan_PHS_noAPOE"] = (phs, len(dcols))
            RUNS.append(("desikan_phs", "Desikan_PHS_noAPOE",
                         ("precomp", "Desikan_PHS_noAPOE"), []))
            # + APOE-ε term from clinical APOE_Alleles (rs429358 not on
            # array). desikan_apoe_weights.tsv lives next to --desikan-pgs
            # (genetic_baselines/); fall back to --enriched-weights' dir.
            apw_f = Path(args.desikan_pgs).parent / "desikan_apoe_weights.tsv"
            if not apw_f.exists():
                apw_f = (Path(args.enriched_weights).parent
                         / "desikan_apoe_weights.tsv")
            w_e2 = w_e4 = 0.0
            if apw_f.exists():
                apw = pd.read_csv(apw_f, sep="\t", index_col=0)["weight"]
                w_e2 = float(apw.get("APOE_e2", 0.0))
                w_e4 = float(apw.get("APOE_e4", 0.0))
            al = pd.concat([pd.read_csv(split_root / f"seed_{args.seeds[0]}" / f"{sp}.csv",
                                        usecols=["Patient_ID", "APOE_Alleles"])
                            for sp in ("train", "val", "test")],
                           ignore_index=True).drop_duplicates("Patient_ID")
            al["Patient_ID"] = al["Patient_ID"].astype(str)

            def _apoe_term(a):
                a = str(a)
                return a.count("4") * w_e4 + a.count("2") * w_e2
            term = al.set_index("Patient_ID")["APOE_Alleles"].map(_apoe_term)
            phs_apoe = phs.add(term.reindex(phs.index).fillna(0.0), fill_value=0.0)
            _precomp["Desikan_PHS_withAPOE"] = (phs_apoe, len(dcols))
            RUNS.append(("desikan_phs", "Desikan_PHS_withAPOE",
                         ("precomp", "Desikan_PHS_withAPOE"), []))
        print(f"[enrich] PRS sets: {list(prs_sets)}  (each ×{len(PRS_COMBOS)} "
              f"covariate combos) + "
              f"{len([r for r in RUNS if r[0]=='desikan_phs'])} Desikan-PHS",
              flush=True)

    _cache: dict = {}
    def _feat(kind, key):
        ck = (kind, key)
        if ck in _cache:
            return _cache[ck]
        if kind == "precomp":
            s, n = _precomp[key]
            v = (pd.DataFrame({"PRS": s}), n)
        elif kind == "precompmat":             # raw multi-col dosage matrix
            df, n = _precomp[key]
            v = (df.copy(), n)
        else:
            ss = snps_for_tier(bt, all_snps, key, args.apoe_esnps, cluster)
            if kind == "dosage":
                cols = [s for s in ss if s in dos.columns]
                v = (dos[cols].copy(), len(cols))
            else:                                  # beta → z-scored on train
                v = (pd.DataFrame({"PRS": _prs(ss, "beta_A1")}), len(ss))
        v[0].index = v[0].index.astype(str)
        _cache[ck] = v
        return v

    if args.label_mode in STRICT_MODES:
        labels = _strict_labels(split_root, args.labels_tsv,
                                args.label_mode, args.seeds)
    else:
        labels = {sd: load_labels(split_root, sd,
                                  label_mode=args.label_mode)
                  for sd in args.seeds}
    raw_splits = {sd: {sp: pd.read_csv(split_root / f"seed_{sd}" / f"{sp}.csv")
                       for sp in ("train", "val", "test")}
                  for sd in args.seeds}

    done = skipped = failed = 0
    summary = []
    for model, fset, prs_spec, split_cols in RUNS:
        for sd in args.seeds:
            out_dir = args.out_root / args.label_mode / model / fset / f"seed_{sd}"
            mfile = out_dir / "metrics.json"
            if args.skip_if_exists and mfile.exists():
                skipped += 1
                continue
            try:
                feat, n_used = (None, 0)
                if prs_spec is not None:
                    feat, n_used = _feat(*prs_spec)
                parts = {}
                for sp in ("train", "val", "test"):
                    lab = labels[sd][sp].copy()
                    lab["Patient_ID"] = lab["Patient_ID"].astype(str)
                    cov = _covars(raw_splits[sd][sp])
                    cov.index = cov.index.astype(str)
                    rows = []
                    for pid, y in zip(lab["Patient_ID"], lab["y"]):
                        if feat is not None and pid not in feat.index:
                            continue
                        rec = {"y": int(y)}
                        if feat is not None:
                            rec.update(feat.loc[pid].to_dict())
                        for cc in split_cols:
                            rec[cc] = (float(cov.loc[pid, cc])
                                       if pid in cov.index
                                       and np.isfinite(cov.loc[pid, cc])
                                       else np.nan)
                        rows.append(rec)
                    parts[sp] = pd.DataFrame(rows)
                # TRAIN-ONLY fit: StandardScaler mean/std AND the logistic
                # coefficients use the train fold only (no val leakage); val
                # is evaluated for reference. Tables report TEST.
                tr = parts["train"].reset_index(drop=True)
                Xcols = [c for c in tr.columns if c != "y"]
                if not Xcols or tr.empty:
                    raise ValueError(f"no features ({Xcols}) / empty train")
                pipe = Pipeline([
                    ("imp", SimpleImputer(strategy="median")),
                    ("scl", StandardScaler()),
                    ("clf", LogisticRegression(max_iter=2000,
                                               class_weight="balanced",
                                               random_state=sd)),
                ])
                pipe.fit(tr[Xcols].astype(float), tr["y"].astype(int))
                te, vp = parts["test"], parts["val"]
                test_m = _metrics(te["y"].astype(int).values,
                                  pipe.predict_proba(te[Xcols].astype(float))[:, 1])
                val_m = _metrics(vp["y"].astype(int).values,
                                 pipe.predict_proba(vp[Xcols].astype(float))[:, 1])
                # OR per +1 SD of the PRS, using THIS seed's empirical TRAIN
                # PRS SD. Pooled across seeds (with 95% CI) by 20b.
                notes = {}
                if "PRS" in Xcols:
                    j = Xcols.index("PRS")
                    coef = float(pipe.named_steps["clf"].coef_[0][j])
                    scl = float(pipe.named_steps["scl"].scale_[j]) or 1.0
                    logor_per_unit = coef / scl
                    train_sd = float(tr["PRS"].astype(float).std(ddof=0))
                    notes = {"train_prs_sd": train_sd,
                             "logor_per_unit": logor_per_unit,
                             "or_per_train_sd":
                             math.exp(logor_per_unit * train_sd)}

                def _pn(p):
                    y = (p["y"].astype(int) if len(p)
                         else pd.Series([], dtype=int))
                    return int((y == 1).sum()), int((y == 0).sum())
                trp, vpn, tep = _pn(tr), _pn(vp), _pn(te)
                out_dir.mkdir(parents=True, exist_ok=True)
                json.dump({"label_mode": args.label_mode, "model": model,
                           "method": model, "aggregation": fset, "head": model,
                           "seed": sd, "n_snps": int(n_used),
                           "n_dropped_ea_mismatch": int(n_drop),
                           "n_train_pos": trp[0], "n_train_neg": trp[1],
                           "n_val_pos": vpn[0], "n_val_neg": vpn[1],
                           "n_test_pos": tep[0], "n_test_neg": tep[1],
                           "val": val_m, "test": test_m, "notes": notes},
                          open(mfile, "w"), indent=2)
                done += 1
                summary.append((model, fset, sd, test_m["auc"],
                                test_m["balanced_accuracy"]))
                print(f"[ok] {model}/{fset}/seed{sd}  "
                      f"test_auc={test_m['auc']:.3f} "
                      f"bACC={test_m['balanced_accuracy']:.3f} (n={n_used})",
                      flush=True)
            except Exception as e:
                failed += 1
                print(f"[FAIL] {model}/{fset}/seed{sd}: {type(e).__name__}: {e}",
                      flush=True)

    print(f"\n[done] {done} ok, {skipped} skipped, {failed} failed "
          f"→ {args.out_root}")
    if summary:
        s = pd.DataFrame(summary, columns=["model", "feature_set", "seed",
                                           "test_auc", "test_bACC"])
        agg = (s.groupby(["model", "feature_set"])
                 .agg(test_auc=("test_auc", "mean"),
                      auc_sd=("test_auc", "std"),
                      bACC=("test_bACC", "mean"))
                 .reset_index().sort_values("test_auc", ascending=False))
        with pd.option_context("display.width", 200, "display.max_rows", None):
            print("\n=== test AUC (mean over seeds), best first ===")
            print(agg.to_string(index=False))


if __name__ == "__main__":
    main()
