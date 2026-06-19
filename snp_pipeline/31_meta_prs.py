r"""
31_meta_prs.py   (LOCAL, env: snp)
==================================
Meta-PRS via L1/L2/elastic-net logistic regression on the 4 source-
specific PRS (sourceB, sourceW, sourceD, sourceK) + a 4-PRS correlation
heatmap (Pearson + Spearman) + label point-biserial correlation row.

Inputs (S2 outputs at `D:\…\GWAS_filtered_maf_resolved\<COMBO>\`):
  sourceB / sourceW / sourceD / sourceK
    <combo>_snps.tsv           rsID, beta_A1, ...
    <combo>_patient_dosage.tsv PTID + per-rsID dosage (0/1/2 A1 count)

For each patient, compute one weighted PRS per source:
  PRS_<src>(p) = Σ_i  beta_A1_i · dosage_{p,i}    over source's SNPs.

Then per label-mode + per seed (TRAIN-only fit, val reported):

  1. Correlation heatmap (5×5 = 4 PRS × 4 PRS + label point-biserial)
     Pearson + Spearman, cohort-wide.
  2. Unregularised LR on 4 PRS — coefficient α_B/α_W/α_D/α_K + AUC etc.
  3. Elastic-net LR (`penalty='elasticnet', solver='saga', l1_ratio ∈
     {0.1,0.3,0.5,0.7,0.9}` × `C ∈ logspace(-2,2,5)`) selected by
     val-AUC max — sparse coefficients + AUC.
  4. Same with +APOE4+age+sex (per [[evaluate-on-val-not-test]]).

Outputs → `D:\…\genetic_baselines_meta\<label_mode>\`:
  prs_correlation.csv            5×5 Pearson + Spearman + p-values
  prs_correlation_heatmap.{png,pdf}
  meta_prs_metrics.json          per-seed val/test AUC + bACC + CMs
  meta_prs_coefficients.{csv,png}  α,β,γ,ε weights (mean ± SD)
  meta_prs_table.{png,md}        unified table (LR / elastic-net /
                                  +APOE4+age+sex variants × 4 source-PRS
                                  features × 3 seeds → mean ± SD).

Usage:
  conda run -n snp python snp_pipeline/31_meta_prs.py \
    --label-mode ad_case \
    --labels-tsv "D:/ADNI_SNP_Omni2.5M_20140220/conversion_labels/conversion_labels.tsv"
"""
from __future__ import annotations

import argparse
import importlib.util as _il
import json
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, pointbiserialr
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (roc_auc_score, balanced_accuracy_score, f1_score,
                              recall_score, precision_score, confusion_matrix)

warnings.filterwarnings("ignore")

# ── importlib reuse of script 20 ────────────────────────────────────────────
_S20 = Path(__file__).parent / "20_genetic_baselines.py"
_spec = _il.spec_from_file_location("s20", _S20)
_s20 = _il.module_from_spec(_spec)
_spec.loader.exec_module(_s20)
_covars = _s20._covars
_strict_labels = _s20._strict_labels
load_labels = _s20.load_labels
SPLITS = _s20.SPLITS
STRICT_MODES = _s20.STRICT_MODES

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
COMBO_ROOT = BASE / "GWAS_filtered_maf_resolved"
SOURCES = ["sourceB", "sourceW", "sourceD", "sourceK"]
PRETTY = {"sourceB": "Bellenguez", "sourceW": "Wightman",
           "sourceD": "Desikan", "sourceK": "Kunkle"}
PALETTE = "RdBu_r"

plt.rcParams["font.family"] = "DejaVu Serif"


def compute_per_source_prs(combo_root: Path) -> pd.DataFrame:
    """Return DataFrame indexed by PTID with one PRS_<source> column."""
    out = None
    for src in SOURCES:
        d = combo_root / src
        snp = pd.read_csv(d / f"{src}_snps.tsv", sep="\t").set_index("rsID")
        dos = pd.read_csv(d / f"{src}_patient_dosage.tsv", sep="\t")
        dos["PTID"] = dos["PTID"].astype(str)
        dos = dos.set_index("PTID").apply(lambda c: c.fillna(c.mean()))
        cols = [s for s in snp.index if s in dos.columns]
        w = snp.loc[cols, "beta_A1"].astype(float).values
        prs = dos[cols].astype(float).mul(w, axis=1).sum(axis=1)
        prs.name = f"PRS_{PRETTY[src]}"
        out = prs.to_frame() if out is None else out.join(prs, how="outer")
    return out


def _metrics(y, proba):
    pred = (proba > 0.5).astype(int)
    y = np.asarray(y).astype(int)
    out = {
        "auc": float(roc_auc_score(y, proba)) if len(np.unique(y)) > 1 else float("nan"),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "f1": float(f1_score(y, pred, zero_division=0)),
        "sensitivity": float(recall_score(y, pred, pos_label=1, zero_division=0)),
        "specificity": float(recall_score(y, pred, pos_label=0, zero_division=0)),
        "precision": float(precision_score(y, pred, pos_label=1, zero_division=0)),
    }
    cm = confusion_matrix(y, pred, labels=[0, 1])
    out["confusion_matrix"] = {"tn": int(cm[0, 0]), "fp": int(cm[0, 1]),
                                "fn": int(cm[1, 0]), "tp": int(cm[1, 1])}
    return out


def correlation_table(prs: pd.DataFrame, y: pd.Series) -> pd.DataFrame:
    """Pairwise Pearson + Spearman among 4 PRS; point-biserial against y."""
    prs = prs[~prs.index.duplicated(keep="first")]
    cols = list(prs.columns)
    rows = []
    # PRS×PRS
    for i, a in enumerate(cols):
        for b in cols[i+1:]:
            x = prs[[a, b]].dropna()
            xa = np.asarray(x[a]).astype(float).ravel()
            xb = np.asarray(x[b]).astype(float).ravel()
            pr = pearsonr(xa, xb)
            sr = spearmanr(xa, xb)
            rows.append({"x": a, "y": b,
                          "pearson_r": float(np.asarray(pr.statistic).item()),
                          "pearson_p": float(np.asarray(pr.pvalue).item()),
                          "spearman_r": float(np.asarray(sr.statistic).item()),
                          "spearman_p": float(np.asarray(sr.pvalue).item()),
                          "n": int(len(x))})
    # PRS vs y (point-biserial)
    for a in cols:
        m = pd.concat([prs[a], y], axis=1).dropna()
        if len(m) < 4:
            continue
        yy = np.asarray(m.iloc[:, 1]).astype(int).ravel()
        xx = np.asarray(m.iloc[:, 0]).astype(float).ravel()
        pb = pointbiserialr(yy, xx)
        rows.append({"x": a, "y": "ad_label",
                      "pearson_r": float(np.asarray(pb.statistic).item()),
                      "pearson_p": float(np.asarray(pb.pvalue).item()),
                      "spearman_r": float("nan"),
                      "spearman_p": float("nan"), "n": int(len(m))})
    return pd.DataFrame(rows)


def render_heatmap(prs: pd.DataFrame, y: pd.Series, out_png: Path,
                   title: str) -> None:
    """5×5 heatmap: rows/cols = 4 PRS + ad_label. Pearson r in cells."""
    cols = list(prs.columns) + ["ad_label"]
    M = np.zeros((len(cols), len(cols))) * np.nan
    P = np.full_like(M, np.nan)
    base = prs.copy()
    base = base[~base.index.duplicated(keep="first")]
    base["ad_label"] = y.astype(float).reindex(base.index)
    for i, a in enumerate(cols):
        for j, b in enumerate(cols):
            x = base[[a, b]].dropna()
            if len(x) < 4:
                continue
            xa = np.asarray(x[a]).astype(float).ravel()
            xb = np.asarray(x[b]).astype(float).ravel()
            if a == "ad_label" or b == "ad_label":
                if a == b:
                    M[i, j] = 1.0; P[i, j] = 0.0; continue
                if a == "ad_label":
                    pb = pointbiserialr(xa.astype(int), xb)
                else:
                    pb = pointbiserialr(xb.astype(int), xa)
                M[i, j] = float(np.asarray(pb.statistic).item())
                P[i, j] = float(np.asarray(pb.pvalue).item())
            else:
                pr = pearsonr(xa, xb)
                M[i, j] = float(np.asarray(pr.statistic).item())
                P[i, j] = float(np.asarray(pr.pvalue).item())

    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    im = ax.imshow(M, cmap=PALETTE, vmin=-1, vmax=1)
    ax.set_xticks(range(len(cols)))
    ax.set_yticks(range(len(cols)))
    short = [c.replace("PRS_", "") for c in cols]
    ax.set_xticklabels(short, rotation=45, ha="right")
    ax.set_yticklabels(short)
    for i in range(len(cols)):
        for j in range(len(cols)):
            v = M[i, j]
            if not np.isfinite(v):
                continue
            star = " *" if P[i, j] < 0.05 else ""
            ax.text(j, i, f"{v:+.2f}{star}", ha="center", va="center",
                    fontsize=8, color="white" if abs(v) > 0.5 else "black")
    ax.set_title(title, fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Pearson r")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def fit_lr(Xtr, ytr, Xv, yv, kind: str, seed: int):
    """Fit LR or elastic-net LR; return (model, val_metrics, coefs).
    Elastic-net: grid over l1_ratio × C, selected by val AUC."""
    Xtr = Xtr.astype(float); Xv = Xv.astype(float)
    if kind == "lr":
        pipe = Pipeline([("imp", SimpleImputer(strategy="median")),
                          ("scl", StandardScaler()),
                          ("clf", LogisticRegression(max_iter=4000,
                                                      class_weight="balanced",
                                                      random_state=seed))])
        pipe.fit(Xtr, ytr.astype(int))
        proba_v = pipe.predict_proba(Xv)[:, 1]
        m = _metrics(yv, proba_v)
        coefs = dict(zip(Xtr.columns,
                          pipe.named_steps["clf"].coef_[0]))
        return pipe, m, coefs, {"kind": "lr"}
    if kind == "elasticnet":
        best = None
        for l1 in (0.1, 0.3, 0.5, 0.7, 0.9):
            for C in np.logspace(-2, 2, 5):
                p = Pipeline([("imp", SimpleImputer(strategy="median")),
                                ("scl", StandardScaler()),
                                ("clf", LogisticRegression(
                                    penalty="elasticnet", solver="saga",
                                    l1_ratio=l1, C=C, max_iter=8000,
                                    class_weight="balanced",
                                    random_state=seed))])
                p.fit(Xtr, ytr.astype(int))
                proba_v = p.predict_proba(Xv)[:, 1]
                if len(np.unique(yv)) < 2:
                    continue
                auc = float(roc_auc_score(yv, proba_v))
                if best is None or auc > best[0]:
                    best = (auc, p, {"l1_ratio": l1, "C": float(C)},
                             dict(zip(Xtr.columns,
                                       p.named_steps["clf"].coef_[0])))
        auc, p, hp, coefs = best
        proba_v = p.predict_proba(Xv)[:, 1]
        m = _metrics(yv, proba_v)
        return p, m, coefs, {"kind": "elasticnet", **hp}
    raise ValueError(kind)


def assemble_xy(prs: pd.DataFrame, labels: dict, raw_splits: dict, seed: int,
                add_covars: list[str], apoe_stratum: str | None = None
                ) -> dict:
    """Per-split (X, y) frames joining 4 PRS + optional covariates +
    label, restricted to label-mode at-risk patients."""
    out = {}
    for sp in ("train", "val", "test"):
        lab = labels[seed][sp].copy()
        lab["Patient_ID"] = lab["Patient_ID"].astype(str)
        cov = _covars(raw_splits[seed][sp])
        cov.index = cov.index.astype(str)
        rows = []
        for pid, y in zip(lab["Patient_ID"], lab["y"]):
            if pid not in prs.index:
                continue
            a4 = float(cov.loc[pid, "APOE4"]) if pid in cov.index and "APOE4" in cov.columns else float("nan")
            a2 = float(cov.loc[pid, "APOE2"]) if pid in cov.index and "APOE2" in cov.columns else float("nan")
            if apoe_stratum == "e33":
                if not (np.isfinite(a4) and np.isfinite(a2)
                         and a4 == 0 and a2 == 0):
                    continue
            rec = {"y": int(y)}
            for c in prs.columns:
                rec[c] = float(prs.loc[pid, c]) if pid in prs.index else float("nan")
            for cc in add_covars:
                v = float(cov.loc[pid, cc]) if pid in cov.index else float("nan")
                rec[cc] = v if np.isfinite(v) else float("nan")
            rows.append(rec)
        df = pd.DataFrame(rows).dropna(subset=list(prs.columns))
        out[sp] = df
    return out


def _render_meta_png(agg: pd.DataFrame, out_png: Path,
                      title: str, subtitle: list[str]) -> None:
    """20b-style PNG: DejaVu Serif, B&W, thin outline. One row per
    (variant, method). Columns mirror the genetic-baseline tables:
    AUC, BalAcc (bold > 0.5), F1, Sens, Spec, Prec, Recall, OR/+1 SD
    (blank — meta has 4 PRS coefficients, no single PRS β), TN, FP,
    FN, TP. All metrics are VAL fold."""
    # Order rows: PRSonly first, then PRS+covar variants; lr first, then elasticnet
    method_order = {"lr": 0, "elasticnet": 1}
    agg = agg.copy()
    agg["__m"] = agg["method"].map(method_order).fillna(99)
    agg["__v"] = agg["variant"].map(lambda v: 0 if v == "PRSonly" else 1)
    agg = agg.sort_values(["__v", "__m"]).reset_index(drop=True)

    # Columns: (key, label, width, alignment)
    cols = [
        ("variant",   "Variant",       0.18, "left"),
        ("method",    "Method",        0.10, "left"),
        ("val_auc",   "AUC",           None, "center"),
        ("val_bacc",  "BalAcc",        None, "center"),
        ("val_f1",    "F1",            None, "center"),
        ("val_sens",  "Sens",          None, "center"),
        ("val_spec",  "Spec",          None, "center"),
        ("val_prec",  "Prec",          None, "center"),
        ("val_recall", "Recall",       None, "center"),
        ("or_sd",     "OR /+1 SD (95% CI)", 0.105, "center"),
        ("val_tn",    "TN",            None, "center"),
        ("val_fp",    "FP",            None, "center"),
        ("val_fn",    "FN",            None, "center"),
        ("val_tp",    "TP",            None, "center"),
    ]
    n_num = sum(1 for c in cols if c[2] is None)
    num_w = (1.0 - 0.18 - 0.10 - 0.105) / max(1, n_num)
    widths = [c[2] if c[2] is not None else num_w for c in cols]
    xs = [0.0]
    for w in widths:
        xs.append(xs[-1] + w)

    n = len(agg)
    row_h, hdr_h, title_h, margin = 0.34, 0.5, 0.5, 0.22
    sub_h = 0.32 * len(subtitle)
    fig_w = 22.0
    fig_h = title_h + sub_h + hdr_h + n * row_h + 2 * margin
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fh = fig_h

    def Y(v): return v / fh

    def box(x, y, w, h, b=0.4):
        ax.add_patch(plt.Rectangle((x, Y(y)), w, h / fh,
                                     transform=ax.transAxes,
                                     facecolor="white", edgecolor="black",
                                     linewidth=b))

    def txt(x, y, s, fs, ha, weight="normal"):
        ax.text(x, Y(y), s, ha=ha, va="center", fontsize=fs,
                 fontweight=weight, transform=ax.transAxes)

    cur = fh - margin
    txt(0.5, cur - title_h * 0.5, title, 12.5, "center")
    cur -= title_h
    for sline in subtitle:
        txt(0.5, cur - 0.16, sline, 8.4, "center")
        cur -= 0.32
    # Header row
    for i, (_, lbl, _, _) in enumerate(cols):
        box(xs[i], cur - hdr_h, xs[i + 1] - xs[i], hdr_h, b=0.7)
        txt((xs[i] + xs[i + 1]) / 2, cur - hdr_h / 2, lbl, 8.6, "center")
    cur -= hdr_h

    def _ms(mean_col, sd_col, r):
        if pd.isna(r[mean_col]):
            return ""
        m = float(r[mean_col])
        s = float(r[sd_col]) if (sd_col in r.index and pd.notna(r[sd_col])) else 0.0
        return f"{m:.3f}±{s:.3f}"

    def _cnt(col, r):
        v = r.get(f"{col}_m")
        return f"{float(v):.1f}" if pd.notna(v) else ""

    prev_v = None
    for _, r in agg.iterrows():
        v = str(r["variant"])
        if prev_v is not None and v != prev_v:
            ax.plot([0, 1], [Y(cur), Y(cur)], transform=ax.transAxes,
                     linewidth=1.1, color="black")
        prev_v = v
        disp = {
            "variant":   v,
            "method":    str(r["method"]),
            "val_auc":   _ms("val_auc_m",  "val_auc_s",  r),
            "val_bacc":  _ms("val_bacc_m", "val_bacc_s", r),
            "val_f1":    _ms("val_f1_m",   "val_f1_s",   r),
            "val_sens":  _ms("val_sens_m", "val_sens_s", r),
            "val_spec":  _ms("val_spec_m", "val_spec_s", r),
            "val_prec":  _ms("val_prec_m", "val_prec_s", r),
            "val_recall": _ms("val_sens_m", "val_sens_s", r),  # Recall = Sens
            "or_sd":     "",
            "val_tn":    _cnt("val_tn", r),
            "val_fp":    _cnt("val_fp", r),
            "val_fn":    _cnt("val_fn", r),
            "val_tp":    _cnt("val_tp", r),
        }
        bacc = r.get("val_bacc_m")
        bacc_bold = pd.notna(bacc) and float(bacc) > 0.5
        for i, (key, _, _, al) in enumerate(cols):
            box(xs[i], cur - row_h, xs[i + 1] - xs[i], row_h)
            val = disp.get(key, "")
            w = "bold" if (key == "val_bacc" and bacc_bold) else "normal"
            if al == "left":
                txt(xs[i] + 0.005, cur - row_h / 2, val, 7.8, "left", w)
            else:
                txt((xs[i] + xs[i + 1]) / 2, cur - row_h / 2, val,
                     7.4, "center", w)
        cur -= row_h

    fig.savefig(out_png, dpi=150, facecolor="white",
                 bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"Wrote {out_png}  ({n} rows, {fig_w} x {fig_h:.1f} in)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--combo-root", type=Path, default=COMBO_ROOT)
    ap.add_argument("--out-root", type=Path,
                     default=BASE / "genetic_baselines_meta",
                     help="Output root for full-cohort runs. e33 stratum "
                          "outputs are auto-routed to a sibling folder "
                          "<out-root>_apoe33/ per the parallel-output-folders "
                          "convention.")
    ap.add_argument("--label-mode", default="ad_case",
                     choices=["ever_convert", *STRICT_MODES.keys()])
    ap.add_argument("--labels-tsv", type=Path,
                     default=BASE / "conversion_labels/conversion_labels.tsv")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--apoe-stratum", default=None,
                     choices=[None, "e33", "e4_carrier", "e2_carrier"])
    args = ap.parse_args()

    # Route e33 outputs to a parallel folder (matches the
    # `genetic_baselines_apoe33/` convention used by 20c).
    _STRATUM_SUFFIX = {"e33": "_apoe33",
                        "e4_carrier": "_apoe4carrier",
                        "e2_carrier": "_apoe2carrier"}
    if args.apoe_stratum:
        out_root = args.out_root.with_name(
            args.out_root.name + _STRATUM_SUFFIX[args.apoe_stratum])
    else:
        out_root = args.out_root
    out = out_root / args.label_mode
    out.mkdir(parents=True, exist_ok=True)
    print(f"[setup] out={out}  label_mode={args.label_mode}  "
          f"apoe_stratum={args.apoe_stratum}")

    # ── 1. Compute 4 per-source PRS ──────────────────────────────────────
    prs = compute_per_source_prs(args.combo_root)
    prs.index = prs.index.astype(str)
    print(f"[prs] {prs.shape[0]} patients × {prs.shape[1]} PRS")
    prs.to_csv(out / "patient_source_prs.tsv", sep="\t")

    # ── 2. Labels + splits ───────────────────────────────────────────────
    split_root = SPLITS.get(args.label_mode) or SPLITS["ever_convert"]
    if args.label_mode in STRICT_MODES:
        labels = _strict_labels(split_root, args.labels_tsv,
                                 args.label_mode, args.seeds)
    else:
        labels = {sd: load_labels(split_root, sd, label_mode=args.label_mode)
                   for sd in args.seeds}
    raw_splits = {sd: {sp: pd.read_csv(split_root / f"seed_{sd}" / f"{sp}.csv")
                       for sp in ("train", "val", "test")}
                   for sd in args.seeds}

    # ── 3. Cohort-wide correlation heatmap (label = AD case across all seeds) ──
    all_pids, all_y = [], []
    for sd in args.seeds:
        for sp in ("train", "val", "test"):
            lab = labels[sd][sp]
            for pid, y in zip(lab["Patient_ID"].astype(str), lab["y"]):
                if pid not in prs.index:
                    continue
                if args.apoe_stratum == "e33":
                    cov = _covars(raw_splits[sd][sp])
                    cov.index = cov.index.astype(str)
                    if pid not in cov.index:
                        continue
                    a4 = float(cov.loc[pid, "APOE4"])
                    a2 = float(cov.loc[pid, "APOE2"])
                    if not (np.isfinite(a4) and np.isfinite(a2)
                             and a4 == 0 and a2 == 0):
                        continue
                all_pids.append(pid); all_y.append(int(y))
    cohort_y = pd.Series(all_y, index=all_pids).groupby(level=0).first()
    cohort_prs = prs.loc[cohort_y.index]

    ct = correlation_table(cohort_prs, cohort_y)
    ct.to_csv(out / "prs_correlation.csv", index=False)
    n = len(cohort_y); pos = int(cohort_y.sum())
    render_heatmap(cohort_prs, cohort_y,
                    out / "prs_correlation_heatmap.png",
                    f"4-PRS correlation — {args.label_mode}"
                    + (f" / {args.apoe_stratum}" if args.apoe_stratum else "")
                    + f"  (n={n}, AD={pos})")
    print(f"[corr] {len(ct)} rows  → prs_correlation.csv + heatmap.png/pdf")

    # ── 4. Meta-PRS LR + elastic-net per seed (× 2 covariate variants) ───
    # In e33 (APOE4==0 ∧ APOE2==0 by stratum construction), APOE4 is a
    # constant zero column → degenerate feature. Swap the variant to
    # "PRS+age+sex" so the row label honestly reflects what the model
    # is using.
    if args.apoe_stratum == "e33":
        variants = [("PRSonly", []),
                      ("PRS+age+sex", ["AGE", "SEX"])]
    else:
        variants = [("PRSonly", []),
                      ("PRS+APOE4+age+sex", ["APOE4", "AGE", "SEX"])]
    methods = ["lr", "elasticnet"]
    long_rows = []
    coef_rows = []
    metrics_dict = {}
    for vname, vcov in variants:
        for kind in methods:
            for sd in args.seeds:
                parts = assemble_xy(prs, labels, raw_splits, sd, vcov,
                                     apoe_stratum=args.apoe_stratum)
                tr = parts["train"]; vp = parts["val"]; te = parts["test"]
                if len(tr) < 10 or len(vp) < 5:
                    print(f"  [skip] {vname}/{kind}/seed{sd}: "
                          f"n_train={len(tr)} n_val={len(vp)}")
                    continue
                Xcols = [c for c in tr.columns if c != "y"]
                try:
                    model, val_m, coefs, hp = fit_lr(
                        tr[Xcols], tr["y"], vp[Xcols], vp["y"], kind, sd)
                    # also evaluate on test for completeness (held out)
                    test_proba = model.predict_proba(te[Xcols].astype(float))[:, 1]
                    test_m = _metrics(te["y"], test_proba)
                    cm = val_m.get("confusion_matrix", {})
                    long_rows.append({"variant": vname, "method": kind,
                                       "seed": sd,
                                       "val_auc": val_m["auc"],
                                       "val_bacc": val_m["balanced_accuracy"],
                                       "val_f1": val_m["f1"],
                                       "val_sens": val_m["sensitivity"],
                                       "val_spec": val_m["specificity"],
                                       "val_prec": val_m["precision"],
                                       "val_tn": cm.get("tn", 0),
                                       "val_fp": cm.get("fp", 0),
                                       "val_fn": cm.get("fn", 0),
                                       "val_tp": cm.get("tp", 0),
                                       "n_train": len(tr), "n_val": len(vp),
                                       "n_test": len(te),
                                       **hp})
                    for k, v in coefs.items():
                        coef_rows.append({"variant": vname, "method": kind,
                                            "seed": sd, "feature": k,
                                            "coef": float(v)})
                    metrics_dict[f"{vname}/{kind}/seed_{sd}"] = {
                        "val": val_m, "test": test_m, **hp,
                        "n_train": len(tr), "n_val": len(vp), "n_test": len(te)}
                    print(f"  [{kind:11s}] {vname}/seed{sd}  "
                          f"val_auc={val_m['auc']:.3f}  "
                          f"bacc={val_m['balanced_accuracy']:.3f}  "
                          + (f"l1={hp.get('l1_ratio')} C={hp.get('C'):.2f}"
                             if kind == "elasticnet" else ""))
                except Exception as e:
                    print(f"  [FAIL] {vname}/{kind}/seed{sd}: "
                          f"{type(e).__name__}: {e}")

    long_df = pd.DataFrame(long_rows)
    coef_df = pd.DataFrame(coef_rows)
    long_df.to_csv(out / "meta_prs_long.csv", index=False)
    coef_df.to_csv(out / "meta_prs_coefficients.csv", index=False)
    json.dump(metrics_dict, open(out / "meta_prs_metrics.json", "w"),
               indent=2)

    # ── 5. Mean ± SD over seeds → meta_prs_table.{md,png} ────────────────
    if not long_df.empty:
        # Aggregate every metric per (variant, method) over seeds.
        # Counts are means (whole rounded for display).
        metric_cols = ["val_auc", "val_bacc", "val_f1", "val_sens",
                        "val_spec", "val_prec",
                        "val_tn", "val_fp", "val_fn", "val_tp"]
        agg_full = long_df.groupby(["variant", "method"]).agg(
            **{f"{c}_m": (c, "mean") for c in metric_cols},
            **{f"{c}_s": (c, "std") for c in metric_cols},
            n_train_m=("n_train", "mean"),
            n_val_m=("n_val", "mean"),
            n_test_m=("n_test", "mean"),
        ).reset_index()

        # Summary CSV (mean±std for all metrics)
        agg_disp = pd.DataFrame({"variant": agg_full["variant"],
                                   "method": agg_full["method"]})
        for c in metric_cols:
            if c.startswith("val_t") or c.startswith("val_f"):
                # counts → mean rounded to 1 decimal, no ± SD column
                if c in {"val_tn", "val_fp", "val_fn", "val_tp"}:
                    agg_disp[c] = agg_full[f"{c}_m"].map(lambda x: f"{x:.1f}")
                    continue
            agg_disp[c] = (agg_full[f"{c}_m"].map(lambda x: f"{x:.3f}") + "±"
                            + agg_full[f"{c}_s"].fillna(0).map(lambda x: f"{x:.3f}"))
        agg_disp.to_csv(out / "meta_prs_summary.csv", index=False)

        # MD table (compact form for quick scan)
        md_cols = ["variant", "method", "val_auc", "val_bacc", "val_f1",
                    "val_sens", "val_spec", "val_prec",
                    "val_tn", "val_fp", "val_fn", "val_tp"]
        md = ["| " + " | ".join(md_cols) + " |",
              "|" + "|".join(["---"] * len(md_cols)) + "|"]
        for _, r in agg_disp.iterrows():
            md.append("| " + " | ".join(str(r[c]) for c in md_cols) + " |")
        (out / "meta_prs_table.md").write_text("\n".join(md), "utf-8")

        # 20b-style PNG render (DejaVu Serif, B&W, thin outline)
        n_pos_v = int(round(long_df["n_val"].mean())) if "n_val" in long_df else 0
        title_mode = args.label_mode + (f" / APOE-ε3/ε3"
                                          if args.apoe_stratum == "e33"
                                          else "")
        subtitle = [
            f"N(VAL) ≈ {n_pos_v} subjects  ·  "
            f"meta-PRS = LR (sklearn L2, C=1.0) or elastic-net "
            f"(grid l1_ratio∈[0.1–0.9] × C∈logspace(-2,2,5), "
            f"selected by val AUC).",
            "Logistic regression on z-scored 4-source PRS [+ covariates]; "
            "mean ± std over seeds 0/1/2. **Bold BalAcc > 0.5**. "
            "TN/FP/FN/TP = mean VAL confusion-matrix counts. "
            "OR /+1 SD is undefined for a multi-feature meta-PRS (4 PRS "
            "coefficients) — shown blank.",
        ]
        _render_meta_png(agg_full, out / "meta_prs_table.png",
                          f"Meta-PRS — {title_mode} (ADNI, VAL fold, "
                          f"mean ± std seeds 0/1/2)",
                          subtitle)

        # Coefficient table
        coef_agg = (coef_df.groupby(["variant", "method", "feature"])
                     .agg(mean=("coef", "mean"), sd=("coef", "std"))
                     .reset_index())
        coef_agg.to_csv(out / "meta_prs_coefficients_meanstd.csv",
                         index=False)
        # Render coefficient bar plot (PRSonly variants only, focus on the
        # 4 source-PRS coefficients).
        prs_feats = list(prs.columns)
        variant_names = [v for v, _ in variants]
        for vname in variant_names:
            sub = coef_agg[(coef_agg["variant"] == vname)
                            & coef_agg["feature"].isin(prs_feats)]
            if sub.empty:
                continue
            fig, ax = plt.subplots(figsize=(5.5, 3.2))
            xs = np.arange(len(prs_feats))
            for i, m in enumerate(["lr", "elasticnet"]):
                d = (sub[sub["method"] == m]
                      .set_index("feature").reindex(prs_feats))
                ax.bar(xs + (i - 0.5) * 0.35, d["mean"], 0.35,
                        yerr=d["sd"].fillna(0), capsize=3,
                        label=m, edgecolor="black", linewidth=0.6,
                        color=("#cccccc" if m == "lr" else "#666666"))
            ax.axhline(0, color="black", lw=0.5)
            ax.set_xticks(xs)
            ax.set_xticklabels([f.replace("PRS_", "") for f in prs_feats],
                                rotation=20)
            ax.set_ylabel("Logistic coefficient (z-scored features)")
            ax.set_title(f"Meta-PRS coefficients — {vname} "
                         f"({args.label_mode}"
                         + (f"/{args.apoe_stratum}" if args.apoe_stratum
                            else "") + ")", fontsize=10)
            ax.legend()
            fig.tight_layout()
            fig.savefig(out / f"meta_prs_coefficients_{vname}.png",
                         dpi=180, bbox_inches="tight")
            plt.close(fig)

    print(f"\n[done] → {out}")


if __name__ == "__main__":
    main()
