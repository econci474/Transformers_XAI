"""
_prs_quantile_plot.py
=====================
20-quantile PRS plot comparing three risk scores:
  1. FM-PRS — v3 diff-attention top classifier predicted probability
     (bmfm_snp / per_modality / mlp2 / 101bp / chrom_hier / w=256 / d=0.3 / lr=1e-3),
     pooled across seeds 0/1/2. Per-patient predictions from
     <PRED_DIR>/predictions_concat.tsv (concatenated from the per-seed
     predictions.tsv files written by --save-predictions).
  2. LD-top PRS_only — Kosteridis_shared_AD_CV β × dosage on the LD-pruned
     strict-QC pool at ld_1000kb_r2_0.8 (the prs_only top performer per
     the master leaderboard).
  3. DeRojas PRS-CS — DeRojas posterior β × dosage on the unpruned 115-SNP
     strict-QC pool (the prs_only top performer in the PRS-CS leaderboard).

Two panels per cohort split (val + test):
  A) AD case prevalence per quantile (line + scatter).
  B) Odds ratio per quantile vs reference (middle) quantile, with 95% CI.
Lee (2012) liability-scale pseudo-R² is reported in each PRS's legend
entry. Sample prevalence K_s = positives / total in the pooled cohort;
population prevalence K_p defaults to 0.10 (10% in 65+) but is a CLI flag.

Usage (env: snp):
  python snp_pipeline/_prs_quantile_plot.py
  python snp_pipeline/_prs_quantile_plot.py --split test
  python snp_pipeline/_prs_quantile_plot.py --K_p 0.05
"""
from __future__ import annotations
import argparse
import math
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, balanced_accuracy_score

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _prs_strict_qc_lib import (load_source_beta, compute_prs,
                                  load_ldpruned_dosage, load_strictqc_dosage)

ROOT = Path("D:/ADNI_SNP_Omni2.5M_20140220")
PRED_DIR = ROOT / "outputs/diff_attn_v3/fm_top_predictions"
SPLITS = Path("D:/ADNI_BIDS_project/derivatives/clinical/"
               "no_cdr_stratified_post_exclusion/tabular/baseline")
OUT_DIR = ROOT / "outputs/prs_quantile_plot"
N_BINS = 20

# Standard reference for AD population prevalence (~10% in 65+).
DEFAULT_K_POP = 0.10


def _load_fm_predictions() -> pd.DataFrame:
    """Concatenate per-seed predictions.tsv from --save-predictions runs."""
    paths = sorted(PRED_DIR.rglob("predictions.tsv"))
    if not paths:
        print(f"[warn] no predictions.tsv under {PRED_DIR}")
        return pd.DataFrame()
    frames = [pd.read_csv(p, sep="\t") for p in paths]
    df = pd.concat(frames, ignore_index=True)
    print(f"  loaded FM predictions from {len(paths)} files, "
          f"{len(df)} rows (across seeds + splits)")
    return df


def _load_labels(seed: int, split: str) -> pd.DataFrame:
    df = pd.read_csv(SPLITS / f"seed_{seed}/{split}.csv", dtype=str)
    df["Patient_ID"] = df["Patient_ID"].astype(str)
    pos = ((df["AD_bl"].astype(int) == 1)
            | (df["pMCI"].astype(int) == 1)
            | (df["CN_to_AD"].astype(int) == 1)).astype(int)
    neg = (df["sCN"].astype(int) == 1).astype(int)
    keep = (pos == 1) | (neg == 1)
    df = df[keep].copy()
    df["y"] = pos[keep].values
    df["seed"] = seed; df["split"] = split
    return df[["Patient_ID","seed","split","y"]]


def _load_baseline_labels_pooled(split: str) -> pd.DataFrame:
    return pd.concat([_load_labels(s, split) for s in (0, 1, 2)], ignore_index=True)


def _compute_kosteridis_prs() -> pd.DataFrame:
    """Per-patient Kosteridis (unified MTAG_AD ∪ shared_AD_CV) PRS on the
    LD-pruned strict-QC pool at 1000kb/r²<0.8."""
    bim, dos = load_ldpruned_dosage("ld_1000kb_r2_0.8")
    st = load_source_beta("Kosteridis", bim, beta_source="raw")
    if st.empty:
        return pd.DataFrame()
    prs = compute_prs(st, dos)
    out = prs.reset_index().rename(columns={"PTID":"Patient_ID", 0:"prs"})
    out.columns = ["Patient_ID", "prs"]
    return out


def _compute_derojas_prscs() -> pd.DataFrame:
    """Per-patient DeRojas PRS using PRS-CS posterior β on the unpruned strict-QC pool."""
    bim, dos = load_strictqc_dosage()
    st = load_source_beta("DeRojas", bim, beta_source="prscs")
    if st.empty:
        return pd.DataFrame()
    prs = compute_prs(st, dos)
    out = prs.reset_index()
    out.columns = ["Patient_ID", "prs"]
    return out


def _quantile_bins(prs_series: pd.Series, n_bins: int = N_BINS) -> pd.Series:
    """Return integer quantile bin (1..n_bins) per patient, using qcut on ranks
    so duplicates don't bias the binning."""
    ranks = prs_series.rank(method="first")
    q = pd.qcut(ranks, q=n_bins, labels=False, duplicates="drop") + 1
    return q.astype(int)


def _lee_R2(y: np.ndarray, prs: np.ndarray, K_p: float) -> float:
    """Lee et al. 2012 liability-scale pseudo-R² from logistic regression
    Nagelkerke R² on observed scale → liability scale.
    Returns NaN if degenerate."""
    y = np.asarray(y).astype(int); prs = np.asarray(prs).astype(float)
    if len(y) < 8 or len(np.unique(y)) < 2 or np.std(prs) < 1e-9:
        return float("nan")
    K_s = float(y.mean())
    # Standardise PRS
    prs_z = (prs - prs.mean()) / (prs.std(ddof=0) or 1.0)
    # Null + fitted LL
    pi_null = K_s
    ll_null = (y * math.log(pi_null + 1e-12)
                + (1 - y) * math.log(1 - pi_null + 1e-12)).sum()
    lr = LogisticRegression(max_iter=2000).fit(prs_z.reshape(-1, 1), y)
    ll_fit  = lr.predict_log_proba(prs_z.reshape(-1, 1))[np.arange(len(y)), y].sum()
    n = len(y)
    cs = 1.0 - math.exp(2 * (ll_null - ll_fit) / n)              # Cox–Snell
    nag = cs / (1.0 - math.exp(2 * ll_null / n))                  # Nagelkerke
    # Lee (2012) liability-scale conversion (PMID 22714019).
    t = norm.ppf(1 - K_p)
    phi_t = norm.pdf(t)
    if K_s in (0, 1) or phi_t < 1e-12: return float("nan")
    factor = (K_p * (1 - K_p))**2 / (K_s * (1 - K_s) * phi_t**2)
    return float(nag * factor)


def _per_quantile_stats(df: pd.DataFrame, ref_q: int) -> pd.DataFrame:
    """Per quantile: n, n_cases, prevalence, OR vs ref quantile + 95% CI."""
    rows = []
    ref = df[df["q"] == ref_q]
    ref_y = ref["y"].astype(int).values
    a_ref = int(ref_y.sum()); b_ref = int((ref_y == 0).sum())
    for q in sorted(df["q"].unique()):
        sub = df[df["q"] == q]
        y = sub["y"].astype(int).values
        n = len(y); n_pos = int(y.sum())
        prev = n_pos / n if n else float("nan")
        a, b = n_pos, n - n_pos
        # OR vs reference; Haldane–Anscombe 0.5 continuity to avoid zero/inf
        oa, ob = a + 0.5, b + 0.5
        ora, orb = a_ref + 0.5, b_ref + 0.5
        OR = (oa * orb) / (ob * ora)
        se_logOR = math.sqrt(1/oa + 1/ob + 1/ora + 1/orb)
        lo = math.exp(math.log(OR) - 1.96 * se_logOR)
        hi = math.exp(math.log(OR) + 1.96 * se_logOR)
        rows.append({"q": int(q), "n": int(n), "n_pos": int(n_pos),
                      "prevalence": prev, "OR": OR, "OR_lo": lo, "OR_hi": hi})
    return pd.DataFrame(rows)


def _safe(merge_args):
    df, prs, name = merge_args
    if prs is None or prs.empty: return None
    m = df.merge(prs, on="Patient_ID", how="inner").rename(columns={"prs": name})
    m = m[m[name].notna()].copy()
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", default="test", choices=["val","test","both"])
    ap.add_argument("--K_p", type=float, default=DEFAULT_K_POP,
                     help=f"AD population prevalence for Lee R² (default {DEFAULT_K_POP}).")
    ap.add_argument("--n-bins", type=int, default=N_BINS)
    args = ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"K_p (pop prevalence) = {args.K_p};  bins = {args.n_bins};  split = {args.split}")

    # Per-PRS per-patient scores
    print("\nLoading FM predictions...")
    fm = _load_fm_predictions()
    print("\nComputing Kosteridis (unified MTAG_AD ∪ shared_AD_CV) PRS (LD-pruned, raw β)...")
    kost = _compute_kosteridis_prs()
    print(f"  n_patients with Kosteridis PRS = {len(kost)}")
    print("\nComputing DeRojas PRS-CS (unpruned strict-QC)...")
    der = _compute_derojas_prscs()
    print(f"  n_patients with DeRojas PRS-CS = {len(der)}")

    splits_to_plot = ["val","test"] if args.split == "both" else [args.split]
    for split in splits_to_plot:
        print(f"\n========== split = {split} ==========")
        lab = _load_baseline_labels_pooled(split)
        print(f"  pooled label cohort size: {len(lab)}  "
              f"(POS={int(lab.y.sum())}, NEG={int((lab.y==0).sum())})")

        # Each PRS gets its own merge — n_patients may differ across PRS
        if not fm.empty:
            fm_split = fm[fm["split"] == split][["Patient_ID","seed","prob"]].copy()
            # Average prob across seeds for patients in multiple seed splits
            fm_avg = fm_split.groupby("Patient_ID")["prob"].mean().reset_index()
            fm_avg["Patient_ID"] = fm_avg["Patient_ID"].astype(str)
            fm_merged = lab.merge(fm_avg, on="Patient_ID", how="inner").rename(columns={"prob":"prs"})
        else:
            fm_merged = pd.DataFrame()
        kost_merged = lab.merge(kost.assign(Patient_ID=kost["Patient_ID"].astype(str)),
                                  on="Patient_ID", how="inner").rename(columns={"prs":"prs"})
        kost_merged = kost_merged[kost_merged["prs"].notna()].copy()
        der_merged = lab.merge(der.assign(Patient_ID=der["Patient_ID"].astype(str)),
                                 on="Patient_ID", how="inner").rename(columns={"prs":"prs"})
        der_merged = der_merged[der_merged["prs"].notna()].copy()

        prs_panels = [
            ("FM diff-attention (bmfm_snp/per_modality)", fm_merged),
            ("Kosteridis (LD-pruned, raw β, MTAG_AD ∪ shared_AD_CV)", kost_merged),
            ("DeRojas PRS-CS (unpruned 115-SNP pool)",      der_merged),
        ]

        ref_q = args.n_bins // 2  # middle quantile as the reference for OR

        # Compute per-PRS stats + Lee R²
        per_prs = []
        for name, m in prs_panels:
            if m.empty:
                print(f"  [skip] {name}: empty"); continue
            m = m.copy()
            m["q"] = _quantile_bins(m["prs"], args.n_bins)
            R2 = _lee_R2(m["y"].values, m["prs"].values, K_p=args.K_p)
            stats = _per_quantile_stats(m, ref_q)
            per_prs.append({"name": name, "data": m, "stats": stats, "R2": R2})
            print(f"  {name}: n={len(m)}  Lee R² = {R2:.4f}")
            stats_p = OUT_DIR / f"{split}_{name.replace(' ', '_').replace('/', '-')}_stats.tsv"
            stats.to_csv(stats_p, sep="\t", index=False)

        if not per_prs:
            print(f"  no PRS available for {split} — skipping plot"); continue

        # Plot: 2 panels (prevalence + OR), 3 curves
        fig, axs = plt.subplots(1, 2, figsize=(14, 5.5))
        colors = {"FM diff-attention (bmfm_snp/per_modality)": "tab:blue",
                  "Kosteridis (LD-pruned, raw β, MTAG_AD ∪ shared_AD_CV)": "tab:orange",
                  "DeRojas PRS-CS (unpruned 115-SNP pool)": "tab:green"}
        for d in per_prs:
            c = colors.get(d["name"], "k")
            label = f"{d['name']}\nLee R²={d['R2']:.3f}, n={len(d['data'])}"
            # Panel A: prevalence
            axs[0].plot(d["stats"]["q"], d["stats"]["prevalence"], "o-",
                         color=c, label=label, markersize=4, lw=1.5)
            # Panel B: OR vs ref with CI
            axs[1].errorbar(d["stats"]["q"], d["stats"]["OR"],
                             yerr=[d["stats"]["OR"] - d["stats"]["OR_lo"],
                                    d["stats"]["OR_hi"] - d["stats"]["OR"]],
                             fmt="o-", color=c, label=label, markersize=4,
                             capsize=2, lw=1.2)
        axs[0].set_xlabel(f"PRS quantile (1=lowest, {args.n_bins}=highest)")
        axs[0].set_ylabel("AD case prevalence")
        axs[0].set_title(f"AD prevalence by PRS quantile  ({split})", fontsize=10)
        axs[0].grid(alpha=0.3)
        axs[0].set_xticks(range(1, args.n_bins + 1))
        axs[0].legend(fontsize=7, loc="upper left")

        axs[1].axhline(1.0, color="grey", lw=0.6, ls="--")
        axs[1].set_yscale("log")
        axs[1].set_xlabel(f"PRS quantile (ref = q{ref_q})")
        axs[1].set_ylabel("OR vs reference quantile (log scale)")
        axs[1].set_title(f"Odds ratio by PRS quantile  ({split})", fontsize=10)
        axs[1].grid(alpha=0.3, which="both")
        axs[1].set_xticks(range(1, args.n_bins + 1))
        axs[1].legend(fontsize=7, loc="upper left")

        fig.suptitle(f"PRS quantile comparison — pooled seeds 0/1/2 — split={split}   "
                      f"(K_p={args.K_p})",
                     fontsize=11, fontweight="bold")
        fig.tight_layout()
        png = OUT_DIR / f"prs_quantile_plot_{split}.png"
        fig.savefig(png, dpi=200, bbox_inches="tight"); plt.close(fig)
        print(f"  wrote {png}")


if __name__ == "__main__":
    main()
