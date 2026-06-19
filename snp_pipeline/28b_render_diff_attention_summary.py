r"""
28b_render_diff_attention_summary.py   (LOCAL, env: snp)
=======================================================
Aggregate the script-28 results tree and render the summary, reusing the
19b house style (DejaVu Serif, white fill, thin outline, mean ± std over
seeds, ddof=0). Outputs under --root:

  summary_long.csv      one row per run (test+val+train metrics)
  summary_mean_std.csv  grouped by (snp_set,model,agg,depth,delta),
                        mean ± std over seeds (+ n_seeds; flag <3)
  diff_attention_summary.csv  the flat spec table (TEST, per seed):
        model_name,aggregation_mode,mlp_depth,seed,balanced_accuracy,
        AUROC,F1,sensitivity,specificity,precision,recall
  diff_attention_table.{png,md}  mean ± std table (primary = BalAcc)

Usage:
  conda run -n snp python snp_pipeline/28b_render_diff_attention_summary.py \
      --root snp_pipeline/outputs/diff_attn
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

RAW = ["balanced_accuracy", "auc", "f1", "sensitivity", "specificity",
       "precision", "recall"]
KEYS = ["snp_set", "model", "aggregation", "mlp_depth", "delta"]


def collect(root: Path) -> pd.DataFrame:
    rows = []
    for mf in root.rglob("metrics.json"):
        try:
            m = json.load(open(mf))
        except Exception:
            continue
        if "test" not in m:
            continue
        r = {k: m.get(k) for k in KEYS + ["seed", "loss"]}
        r["best_val_auc"] = m.get("fit_info", {}).get("best_val_auc")
        r["gamma"] = m.get("learned", {}).get("gamma")
        r["delta_val"] = m.get("learned", {}).get("delta")
        for sp in ("train", "val", "test"):
            for k in RAW:
                r[f"{sp}_{k}"] = m.get(sp, {}).get(k)
        rows.append(r)
    if not rows:
        raise SystemExit(f"No metrics.json under {root}")
    return pd.DataFrame(rows)


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    out = []
    for key, g in df.groupby(KEYS, dropna=False):
        row = dict(zip(KEYS, key))
        row["n_seeds"] = len(g)
        for k in RAW:
            v = pd.to_numeric(g[f"test_{k}"], errors="coerce").dropna()
            if len(v):
                row[f"test_{k}"] = f"{v.mean():.3f} ± {v.std(ddof=0):.3f}"
                row[f"test_{k}_mean"] = float(v.mean())
            else:
                row[f"test_{k}"] = ""
        out.append(row)
    s = pd.DataFrame(out)
    return s.sort_values(["snp_set", "model", "aggregation", "mlp_depth",
                          "delta"]).reset_index(drop=True)


def render(summ: pd.DataFrame, out_png: Path, out_md: Path) -> None:
    metr = [("test_balanced_accuracy", "BalAcc"), ("test_auc", "AUROC"),
            ("test_f1", "F1"), ("test_sensitivity", "Sens"),
            ("test_specificity", "Spec"), ("test_precision", "Prec"),
            ("test_recall", "Rec")]
    left = KEYS + ["n_seeds"]
    cols = left + [m[0] for m in metr]
    hdr = [c.replace("_", " ") for c in KEYS] + ["n"] + [m[1] for m in metr]
    n = len(summ)
    # Layout: 1 title row + 3 subtitle rows + 1 header row + n data rows
    # All rows share height `rh = 1/(n+5)`. fig_h grows linearly with n.
    fig_w = 6 + 1.5 * len(metr)
    fig_h = 0.9 + 0.34 * (n + 4) + 0.4
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    nc = len(cols)
    lw = [0.07, 0.10, 0.10, 0.06, 0.06] + [0.045] + \
         [(1 - 0.45) / len(metr)] * len(metr)
    lw = [w / sum(lw) for w in lw]
    xs = [0.0]
    for w in lw:
        xs.append(xs[-1] + w)
    rh = 1.0 / (n + 5)

    def box(x, y, w, h):
        ax.add_patch(plt.Rectangle((x, y), w, h, transform=ax.transAxes,
                                   facecolor="white", edgecolor="black",
                                   linewidth=0.4))

    def txt(x, y, s, fs=7.5, ha="center"):
        ax.text(x, y, s, ha=ha, va="center", fontsize=fs,
                transform=ax.transAxes)

    # Title (slot 0)
    txt(0.5, 1 - rh * 0.5,
        "Frozen DNA-FM diff-embedding → β-biased attention — "
        "ad_case vs stable_CN  (TEST fold; mean ± std over seeds; "
        "primary BalAcc; bold > 0.5)",
        fs=10)
    # Subtitles (slots 1,2,3). Definitions match scripts 22/15:
    #   ad_case  = bl_dx=='AD' ∪ ever_conversion_AD (CN/MCI reaching AD,
    #              AD-confirmed at last visit) — baseline-AD plus strict
    #              converters
    #   stable_CN = bl_dx=='CN' AND never any non-CN visit (sustained CN)
    # The two are mutually exclusive; labeled cohort = ad_case ∪ stable_CN.
    sub1 = ("Labels: ad_case = baseline-AD ∪ ever-conversion-AD (CN/MCI "
            "reaching AD, AD at last visit)  ·  stable_CN = sustained CN "
            "(mutually exclusive)")
    sub2 = ("Cohort (SNP+MRI, n=616):  ad_case=182  ·  stable_CN=152  ·  "
            "labeled total=334  (the rest = MCI, reverters, "
            "non-AD-converters — excluded)")
    sub3 = ("N(test) per seed (mean): 36 subjects  ·  ad_case≈20  ·  "
            "stable_CN≈16  (val ≈ 36; train ≈ 262)")
    txt(0.5, 1 - rh * 1.5, sub1, fs=7.6)
    txt(0.5, 1 - rh * 2.5, sub2, fs=7.6)
    txt(0.5, 1 - rh * 3.5, sub3, fs=7.6)
    # Header row (slot 4)
    yb = 1 - 5 * rh
    for i, h in enumerate(hdr):
        box(xs[i], yb, xs[i + 1] - xs[i], rh)
        txt((xs[i] + xs[i + 1]) / 2, yb + rh / 2, h, fs=8)
    for ri in range(n):
        y0 = yb - (ri + 1) * rh
        row = summ.iloc[ri]
        # Bold the BalAcc cell when mean > 0.5 (matches 20b convention)
        bacc = pd.to_numeric(row.get("test_balanced_accuracy_mean"),
                              errors="coerce")
        bacc_bold = pd.notna(bacc) and float(bacc) > 0.5
        for ci, c in enumerate(cols):
            box(xs[ci], y0, xs[ci + 1] - xs[ci], rh)
            v = str(row.get(c, ""))
            weight = "bold" if (c == "test_balanced_accuracy"
                                  and bacc_bold) else "normal"
            ax.text((xs[ci] + xs[ci + 1]) / 2, y0 + rh / 2, v,
                     ha="center", va="center", transform=ax.transAxes,
                     fontsize=6.8 if ci >= len(left) else 7.0,
                     fontweight=weight)
        if ri < n - 1 and str(summ.iloc[ri]["model"]) != \
                str(summ.iloc[ri + 1]["model"]):
            ax.plot([0, 1], [y0, y0], transform=ax.transAxes,
                    linewidth=0.9, color="black")
    fig.savefig(out_png, dpi=150, facecolor="white", bbox_inches="tight",
                pad_inches=0.05)
    plt.close(fig)
    print(f"Wrote {out_png} ({n} rows)")

    lines = ["# Frozen DNA-FM diff-attention — ad_case vs stable_CN", "",
             "TEST metrics, mean ± std over seeds; primary = balanced "
             "accuracy. **Bold BalAcc > 0.5.** Rows with n<3 seeds are "
             "under-powered.", "",
             "| " + " | ".join(hdr) + " |",
             "|" + "|".join(["---"] * len(hdr)) + "|"]
    for ri in range(n):
        row = summ.iloc[ri]
        bacc = pd.to_numeric(row.get("test_balanced_accuracy_mean"),
                              errors="coerce")
        bacc_bold = pd.notna(bacc) and float(bacc) > 0.5
        cells = []
        for c in cols:
            v = str(row.get(c, ""))
            if c == "test_balanced_accuracy" and bacc_bold and v:
                v = f"**{v}**"
            cells.append(v)
        lines.append("| " + " | ".join(cells) + " |")
    out_md.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {out_md}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True)
    args = ap.parse_args()

    long_df = collect(args.root)
    long_df.to_csv(args.root / "summary_long.csv", index=False)
    print(f"Wrote {args.root/'summary_long.csv'} ({len(long_df)} rows)")

    flat = long_df.rename(columns={
        "model": "model_name", "aggregation": "aggregation_mode",
        "mlp_depth": "mlp_depth", "test_auc": "AUROC",
        "test_balanced_accuracy": "balanced_accuracy", "test_f1": "F1",
        "test_sensitivity": "sensitivity", "test_specificity": "specificity",
        "test_precision": "precision", "test_recall": "recall"})[
        ["snp_set", "model_name", "aggregation_mode", "mlp_depth", "delta",
         "seed", "balanced_accuracy", "AUROC", "F1", "sensitivity",
         "specificity", "precision", "recall"]]
    flat.to_csv(args.root / "diff_attention_summary.csv", index=False)
    print(f"Wrote {args.root/'diff_attention_summary.csv'}")

    summ = summarise(long_df)
    summ.to_csv(args.root / "summary_mean_std.csv", index=False)
    print(f"Wrote {args.root/'summary_mean_std.csv'} ({len(summ)} rows)")
    if (summ["n_seeds"] < 3).any():
        print("  [warn] some configs have <3 seeds (under-powered):")
        print(summ.loc[summ["n_seeds"] < 3, KEYS + ["n_seeds"]]
              .to_string(index=False))
    render(summ, args.root / "diff_attention_table.png",
           args.root / "diff_attention_table.md")


if __name__ == "__main__":
    main()
