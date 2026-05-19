"""
20b_render_genetic_baseline_table.py
=====================================
Aggregate the `genetic_baselines/` metrics.json tree (script 20) and render a
B&W publication table (house style: DejaVu Serif, white fill, thin outline,
no bold).

  rows = (group / model / feature_set)
  cols = test  AUC BalAcc F1 Precision Recall Specificity  +  TN FP FN TP
         (mean ± std over seeds; the 4 counts let a confusion matrix be
          reconstructed). Recall == Sensitivity (positive-class TPR).

Row groups are ordered: Demographics → APOE-ε4 → APOE+covar → β-PRS+covar →
β-PRS → z-PRS → LogReg-dosage → Foundation models (frozen).

`--fm-root` reads a prior `frozen_model_comparison/comparison_summary.csv`
(19b output) and appends the frozen-FM rows for context (their confusion
counts are not stored by 19b → shown blank).

Outputs (under --root): comparison_long.csv, comparison_summary.csv,
genetic_baseline_table.png, genetic_baseline_table.md
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.family"] = "DejaVu Serif"

EVAL_LABEL = "VAL"          # set by main() from --eval-split (VAL/TEST)
RAW = ["auc", "balanced_accuracy", "f1", "precision", "recall",
       "sensitivity", "specificity"]
CM_KEYS = ["tn", "fp", "fn", "tp"]
# displayed metric columns (mean ± std)
METRICS = [("test_auc", "AUC"), ("test_balanced_accuracy", "BalAcc"),
           ("test_f1", "F1"), ("test_sensitivity", "Sens"),
           ("test_specificity", "Spec"), ("test_precision", "Prec"),
           ("test_recall", "Recall")]
# displayed confusion-matrix columns (mean count over seeds)
CM = [("test_tn", "TN"), ("test_fp", "FP"), ("test_fn", "FN"),
      ("test_tp", "TP")]

GROUP_ORDER = ["demog", "apoe_e4", "apoe_e2", "apoe_e2e4", "apoe_covar",
               "desikan_phs", "prs_covar", "prs", "logreg_dosage"]
GROUP_PRETTY = {
    "demog": "Demographics only (reference)",
    "apoe_e4": "APOE-ε4 (clinical ε4 dosage)",
    "apoe_e2": "APOE-ε2 (clinical ε2 dosage)",
    "apoe_e2e4": "APOE-ε2 + ε4 (clinical dosages)",
    "apoe_covar": "APOE-ε4 + covariates",
    "desikan_phs": "Desikan PHS (published, external)",
    "prs_covar": "betaPRS_z + covariates",
    "prs": "betaPRS_z (GWAS-weighted, z-scored)",
    "logreg_dosage": "LogReg on raw dosages",
}
# OR per +1 SD of the PRS (pooled across seeds, 95% CI) — PRS rows only
EXTRA = [("or_sd", "OR /+1 SD (95% CI)")]

# One-line case/control criteria per label-mode (shown as a subtitle).
CRITERIA = {
    "ever_convert": "case = ANY non-CN visit ever (weak: prevalent MCI + "
                    "reverters counted); control = always-CN",
    "ever_conversion_AD": "case = strict incident CN/MCI→AD (reaches AD AND "
                          "AD at last visit); prevalent baseline-AD EXCLUDED",
    "ever_conversion_AD_or_MCI": "case = strict →AD OR sustained CN→MCI "
                                 "(confirmed at last visit); baseline-AD excl.",
    "ever_conversion_MCI": "case = strict CN→MCI sustained (MCI/AD at last "
                           "visit); at-risk = baseline-CN",
    "progressive_mci": "case = baseline-MCI reaching AD, AD at last visit",
    "ad_case": "case = prevalent baseline-AD OR strict ever_conversion_AD "
               "(CN/MCI reaching AD, AD-confirmed at last visit); control = "
               "all others (at-risk = everyone)",
    "ad_or_mci_case": "case = prevalent baseline-AD OR strict (→AD or "
                      "sustained CN→MCI); control = all others",
    "bl_multi": "baseline multiclass label (see split definition)",
}


def collect(root: Path, split: str = "val") -> pd.DataFrame:
    rows = []
    for mf in root.rglob("metrics.json"):
        try:
            m = json.load(open(mf))
        except Exception:
            continue
        blk = m.get(split) or m.get("test")
        if blk is None:
            continue
        r = {"group": m.get("model", m.get("method", "?")),
             "feature_set": m.get("aggregation", "?"),
             "seed": m.get("seed", -1),
             "or_per_train_sd": (m.get("notes") or {}).get("or_per_train_sd")}
        for k in RAW:
            r[f"test_{k}"] = blk.get(k)            # key kept; holds `split`
        cm = blk.get("confusion_matrix", {}) or {}
        for c in CM_KEYS:
            r[f"test_{c}"] = cm.get(c)
        for s in ("train", "val", "test"):
            r[f"n_{s}_pos"] = m.get(f"n_{s}_pos")
            r[f"n_{s}_neg"] = m.get(f"n_{s}_neg")
        r["n_snps"] = m.get("n_snps")
        rows.append(r)
    if not rows:
        raise SystemExit(f"No metrics.json under {root}")
    return pd.DataFrame(rows)


def fm_context(fm_root: Path) -> pd.DataFrame:
    csv = fm_root / "comparison_summary.csv"
    if not csv.exists():
        print(f"[fm] {csv} absent — skipping frozen-FM context rows")
        return pd.DataFrame()
    s = pd.read_csv(csv)
    out = []
    for _, x in s.iterrows():
        fs = "/".join(str(x.get(k, "")) for k in ("model", "aggregation", "head")
                      if k in s.columns and pd.notna(x.get(k)))
        row = {"group": "FM", "feature_set": fs or str(x.get("model", "FM"))}
        for mc, _ in METRICS:
            mu, sd = x.get(f"{mc}_mean"), x.get(f"{mc}_std")
            row[mc] = (f"{float(mu):.3f}±{float(sd):.3f}"
                       if pd.notna(mu) and pd.notna(sd)
                       else (f"{float(mu):.3f}" if pd.notna(mu) else ""))
            row[mc + "_mean"] = float(mu) if pd.notna(mu) else np.nan
        for mc, _ in CM:
            row[mc] = ""                       # 19b does not store CM counts
        out.append(row)
    return pd.DataFrame(out)


def summarise(long_df: pd.DataFrame) -> pd.DataFrame:
    out = []
    for (g, fs), grp in long_df.groupby(["group", "feature_set"]):
        row = {"group": g, "feature_set": fs, "n_seeds": len(grp)}
        for k in RAW:
            v = grp[f"test_{k}"].dropna().to_numpy(dtype=float)
            if len(v):
                row[f"test_{k}"] = f"{v.mean():.3f}±{v.std(ddof=0):.3f}"
                row[f"test_{k}_mean"] = float(v.mean())
            else:
                row[f"test_{k}"] = ""
        for c in CM_KEYS:
            v = grp[f"test_{c}"].dropna().to_numpy(dtype=float)
            row[f"test_{c}"] = f"{v.mean():.1f}" if len(v) else ""
        # OR per +1 SD: pool the per-seed (train-SD-based) ORs on the log
        # scale → exp(mean ± 1.96·SE), SE from the between-seed spread.
        ov = grp["or_per_train_sd"].dropna().to_numpy(dtype=float)
        ov = ov[ov > 0]
        if len(ov):
            ln = np.log(ov)
            mu = float(ln.mean())
            se = float(ln.std(ddof=1) / np.sqrt(len(ln))) if len(ln) > 1 else 0.0
            row["or_sd"] = (f"{np.exp(mu):.2f} "
                            f"({np.exp(mu-1.96*se):.2f}–{np.exp(mu+1.96*se):.2f})")
        else:
            row["or_sd"] = ""
        out.append(row)
    return pd.DataFrame(out)


def class_counts(long_df: pd.DataFrame, split: str):
    """Per-seed most-complete split (max total N — i.e. the rows with no
    covariate-NA dropping), then mean pos/neg across seeds. None if absent."""
    pc, nc = f"n_{split}_pos", f"n_{split}_neg"
    if pc not in long_df.columns:
        return None
    pos, neg = [], []
    for _, g in long_df.groupby("seed"):
        tot = (pd.to_numeric(g[pc], errors="coerce").fillna(-1)
               + pd.to_numeric(g[nc], errors="coerce").fillna(-1))
        if tot.empty or tot.max() < 0:
            continue
        r = g.loc[tot.idxmax()]
        pos.append(float(r[pc])); neg.append(float(r[nc]))
    if not pos:
        return None
    return int(round(np.mean(pos))), int(round(np.mean(neg)))


def _order(df: pd.DataFrame) -> pd.DataFrame:
    df = df[df["group"].isin(GROUP_ORDER)].copy()      # drop FM / stray groups
    df["__g"] = pd.Categorical(df["group"], categories=GROUP_ORDER, ordered=True)
    if "test_auc_mean" not in df.columns:
        df["test_auc_mean"] = pd.to_numeric(
            df["test_auc"].astype(str).str.split("±").str[0], errors="coerce")
    return df.sort_values(["__g", "test_auc_mean"], ascending=[True, False])


def render_png(df: pd.DataFrame, out_png: Path, title: str,
               subtitle: list[str] | None = None) -> None:
    df = _order(df)
    n = len(df)
    subtitle = subtitle or []
    cols = [("group", "Group", 0.15, "left"),
            ("feature_set", "Model / feature set", 0.20, "left")]
    cols += [(m[0], m[1], None, "center") for m in METRICS]
    cols += [(e[0], e[1], 0.105, "center") for e in EXTRA]
    cols += [(c[0], c[1], None, "center") for c in CM]
    n_num = len(METRICS) + len(CM)
    num_w = (1.0 - 0.15 - 0.20 - 0.105) / n_num
    widths = [c[2] if c[2] is not None else num_w for c in cols]
    xs = [0.0]
    for w in widths:
        xs.append(xs[-1] + w)

    row_h, hdr_h, title_h, margin = 0.34, 0.5, 0.5, 0.22
    sub_h = 0.32 * len(subtitle)
    fig_w = 25.0
    fig_h = title_h + sub_h + hdr_h + n * row_h + 2 * margin
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fh = fig_h

    def Y(v): return v / fh

    def box(x, y, w, h, b=0.4):
        ax.add_patch(plt.Rectangle((x, Y(y)), w, h / fh, transform=ax.transAxes,
                                   facecolor="white", edgecolor="black",
                                   linewidth=b))

    def txt(x, y, s, fs, ha, weight="normal"):
        ax.text(x, Y(y), s, ha=ha, va="center", fontsize=fs,
                fontweight=weight, transform=ax.transAxes)

    cur = fh - margin
    txt(0.5, cur - title_h * 0.5, title, 13, "center")
    cur -= title_h
    for sline in subtitle:
        txt(0.5, cur - 0.16, sline, 8.4, "center")
        cur -= 0.32
    for i, (_, lbl, _, _) in enumerate(cols):
        box(xs[i], cur - hdr_h, xs[i + 1] - xs[i], hdr_h, b=0.7)
        txt((xs[i] + xs[i + 1]) / 2, cur - hdr_h / 2, lbl, 8.6, "center")
    cur -= hdr_h
    prev = None
    for _, r in df.iterrows():
        g = str(r["group"])
        if prev is not None and g != prev:
            ax.plot([0, 1], [Y(cur), Y(cur)], transform=ax.transAxes,
                    linewidth=1.1, color="black")
        prev = g
        ba = r.get("test_balanced_accuracy_mean")
        ba_bold = pd.notna(ba) and float(ba) > 0.5
        for i, (key, _, _, al) in enumerate(cols):
            box(xs[i], cur - row_h, xs[i + 1] - xs[i], row_h)
            v = GROUP_PRETTY.get(str(r[key]), str(r[key])) if key == "group" else str(r[key])
            w = "bold" if (key == "test_balanced_accuracy" and ba_bold) else "normal"
            if al == "left":
                txt(xs[i] + 0.004, cur - row_h / 2, v, 7.6, "left", w)
            else:
                txt((xs[i] + xs[i + 1]) / 2, cur - row_h / 2, v, 7.4, "center", w)
        cur -= row_h
    fig.savefig(out_png, dpi=150, facecolor="white",
                bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"Wrote {out_png}  ({n} rows, {fig_w} x {fig_h:.1f} in)")


def render_md(df: pd.DataFrame, out_md: Path, title: str,
              subtitle: list[str] | None = None) -> None:
    df = _order(df)
    subtitle = subtitle or []
    hdr = (["Group", "Model / feature set"] + [m[1] for m in METRICS]
           + [e[1] for e in EXTRA] + [c[1] for c in CM])
    lines = [f"# {title}", ""]
    lines += [f"*{s}*" for s in subtitle] + [""]
    lines += [f"**Metrics = {EVAL_LABEL} set**; mean ± std over seeds 0/1/2. "
              "PRS is z-standardised and the model fit on **TRAIN only**, "
              f"evaluated on the held-out **{EVAL_LABEL}** fold. `OR /+1 SD` "
              "= odds ratio per +1 SD of the PRS using each seed's empirical "
              "train SD, pooled across seeds (95% CI, normal approx). "
              f"TN/FP/FN/TP = mean {EVAL_LABEL} confusion-matrix counts. "
              "Sens == Recall. **Bold "
              "BalAcc > 0.5**. PRS tiers (counts in the subtitle): "
              "`withAPOE` = all 430 GW-sig SNPs (rs429358 monomorphic & "
              "rs7412 absent ⇒ inert); `onlyAPOE` = the isolated chr19 APOE "
              "LD block (44.4–46.5 Mb window − rs2072561 − rs429358); "
              "`noAPOE` = APOE region removed (rs2072561 kept as non-APOE "
              "background).", "",
              "| " + " | ".join(hdr) + " |",
              "|" + "|".join(["---"] * len(hdr)) + "|"]
    for _, r in df.iterrows():
        cells = [GROUP_PRETTY.get(str(r["group"]), str(r["group"])),
                 str(r["feature_set"])]
        ba = r.get("test_balanced_accuracy_mean")
        ba_bold = pd.notna(ba) and float(ba) > 0.5
        for mc, _ in METRICS:
            val = str(r[mc])
            cells.append(f"**{val}**" if (mc == "test_balanced_accuracy"
                                          and ba_bold and val) else val)
        cells += [str(r.get(e[0], "")) for e in EXTRA]
        cells += [str(r[c[0]]) for c in CM]
        lines.append("| " + " | ".join(cells) + " |")
    out_md.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {out_md}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path,
                    default=Path("D:/ADNI_SNP_Omni2.5M_20140220/genetic_baselines"))
    ap.add_argument("--label-mode", default="ever_convert")
    ap.add_argument("--eval-split", default="val", choices=["val", "test"],
                    help="which held-out fold's metrics to report "
                         "(default val).")
    ap.add_argument("--fm-root", type=Path, default=None)
    args = ap.parse_args()
    global EVAL_LABEL
    SPL = args.eval_split.upper()
    EVAL_LABEL = SPL

    src = (args.root / args.label_mode
           if (args.root / args.label_mode).exists() else args.root)
    # Write outputs alongside the metrics they summarise (per-label-mode dir
    # when one exists) so multiple label-modes don't overwrite each other;
    # falls back to args.root (byte-identical) for the flat/legacy layout.
    dst = src
    long_df = collect(src, args.eval_split)
    long_df.to_csv(dst / "comparison_long.csv", index=False)
    print(f"Wrote {dst/'comparison_long.csv'} ({len(long_df)} rows)")
    summ = summarise(long_df)
    if args.fm_root:
        print("[fm] --fm-root ignored: foundation-model rows are excluded "
              "from the table (per request — kept it concise).")
    summ.to_csv(dst / "comparison_summary.csv", index=False)
    print(f"Wrote {dst/'comparison_summary.csv'} ({len(summ)} rows)")

    tc = class_counts(long_df, args.eval_split)
    subtitle = []
    if tc:
        subtitle.append(f"N({SPL}): control = {tc[1]}   vs   "
                        f"case = {tc[0]}")
    crit = CRITERIA.get(args.label_mode)
    if crit:
        subtitle.append(f"Criteria: {crit}")
    prs = long_df[long_df["group"] == "prs"]
    sn = []
    for t in ("withAPOE", "onlyAPOE", "noAPOE",
              "6W27B", "9W27B", "9W27B5D"):
        v = prs.loc[prs["feature_set"] == f"betaPRS_z_{t}", "n_snps"].dropna()
        if len(v):
            sn.append(f"{t} = {int(v.iloc[0])} SNPs")
    if sn:
        subtitle.append("betaPRS_z:   " + "   ·   ".join(sn))

    title = (f"Genetic baselines — {args.label_mode} "
             f"(ADNI; {SPL} set; mean ± std, seeds 0/1/2; "
             f"bold BalAcc > 0.5)")
    render_png(summ, dst / "genetic_baseline_table.png", title, subtitle)
    render_md(summ, dst / "genetic_baseline_table.md", title, subtitle)

    show = _order(summ)[["group", "feature_set", "test_auc",
                         "test_recall", "test_specificity"]]
    with pd.option_context("display.width", 200, "display.max_rows", None):
        print("\n=== ordered (group order; AUC desc within group) ===")
        print(show.to_string(index=False))


if __name__ == "__main__":
    main()
