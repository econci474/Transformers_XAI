"""
02_render_failure_table.py
==========================
Render a fusion `failure_table.csv` into colour-coded PNG(s) so per-patient
agreement / complementarity between modalities is readable at a glance.

Schema-agnostic: auto-detects modality columns as `<name>_pred` paired with a
`<name>_correct` (or `<first-token>_correct`) flag. Handles both:
  - MRI-arm  : clinical_pred/clinical_correct, mri_pred/mri_correct
  - integ-arm: clin_temporal_pred/clin_correct, mri_pred/mri_correct, fused_pred/fused_correct

Per CSV it writes, next to it:
  - failure_table_seed{n}.png  : one colour-coded table per seed (green=correct,
                                 red=wrong, grey=modality absent); cell text = predicted class.
  - failure_summary.png        : agreement-pattern counts per seed + overall.

Usage:
  python integration/T2_classification/02_render_failure_table.py \
      --csv integration/T2_classification/outputs/m12_only/failure_table.csv
  # or point at a dir / many dirs and it finds failure_table.csv under each:
  python integration/T2_classification/02_render_failure_table.py \
      --dir integration/T2_classification/outputs
"""
from __future__ import annotations

import argparse
import glob
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GREEN = "#cdebc5"   # correct
RED   = "#f4c7c3"   # wrong
GREY  = "#e6e6e6"   # modality absent / NaN
HEAD  = "#d9d9d9"


def detect_modalities(cols):
    """Return [(label, pred_col, correct_col)] auto-detected from columns."""
    preds = [c for c in cols if c.endswith("_pred")]
    correct = [c for c in cols if c.endswith("_correct")]
    out = []
    for pc in preds:
        pref = pc[:-5]                                  # strip "_pred"
        for cand in (pref + "_correct", pref.split("_")[0] + "_correct"):
            if cand in correct:
                out.append((pref.replace("_", " "), pc, cand))
                break
        else:
            out.append((pref.replace("_", " "), pc, None))
    # Column order: fused first, then clinical, then the rest (stable within rank).
    rank = lambda lbl: 0 if "fused" in lbl else (1 if lbl.startswith("clin") else 2)
    out.sort(key=lambda m: rank(m[0]))
    return out


def derive_agreement(df, mods):
    """If no 'agreement' column, derive one for the clinical+mri pair (2-modality)."""
    if "agreement" in df.columns:
        return df["agreement"].astype(str)
    cc = next((c for lbl, p, c in mods if lbl.startswith("clin") and c), None)
    mc = next((c for lbl, p, c in mods if lbl.startswith("mri") and c), None)
    if cc is None or mc is None:
        return pd.Series(["-"] * len(df), index=df.index)
    out = []
    for _, r in df.iterrows():
        c, m = r[cc], r[mc]
        if pd.isna(m):
            out.append("mri_absent")
        elif bool(c) and bool(m):
            out.append("both_correct")
        elif bool(c) and not bool(m):
            out.append("clinical_only")
        elif not bool(c) and bool(m):
            out.append("mri_only")
        else:
            out.append("both_wrong")
    return pd.Series(out, index=df.index)


FILTER_BANNER = ("FILTERED / WRONG-ENRICHED COHORT — only TEST patients where >=1 modality "
                 "is WRONG are shown (both-correct rows hidden). Do NOT read accuracy off these "
                 "cells; see the REAL full-cohort numbers below.")


def build_stats_footer(stats_path, seed=None):
    """Build the 'real numbers' footer from a sibling failure_stats.csv. If seed is
    None, pool across seeds (raw-acc = Σcorrect/Σn; bACC = mean over seeds)."""
    if not os.path.exists(stats_path):
        return ("(failure_stats.csv not found -- re-run 01_fuse_T2.py to emit real "
                "full-cohort numbers.)")
    st = pd.read_csv(stats_path)
    if seed is not None and "seed" in st.columns:
        st = st[st["seed"] == seed]
    if st.empty:
        return ""
    lines = ["REAL performance on the FULL test cohort (all patients"
             + ("" if seed is not None else ", pooled across seeds") + "):"]
    for mod in ("clinical", "mri"):
        m = st[st["modality"] == mod]
        if m.empty:
            continue
        nc, n = int(m["n_correct"].sum()), int(m["n"].sum())
        raw = nc / n if n else float("nan")
        bacc = float(m["balanced_acc"].mean())
        lines.append(f"  {mod:8s}: raw-acc {raw:.3f} ({nc}/{n} correct)   "
                     f"balanced-acc {bacc:.3f}")
    # cohort counts: one row per seed carries them; dedupe by seed.
    cc = st.drop_duplicates("seed") if "seed" in st.columns else st.head(1)
    n_test = int(cc["n_test"].sum())
    n_both = int(cc["n_both_present"].sum())
    n_hidden = int(cc["n_both_correct_hidden"].sum())
    lines.append(f"  cohort: {n_test} TEST total, {n_both} both-present, "
                 f"{n_hidden} both-correct (HIDDEN from this filtered view).")
    return "\n".join(lines)


def _cell_color(correct):
    if correct is None or (isinstance(correct, float) and np.isnan(correct)):
        return GREY
    return GREEN if bool(correct) else RED


def render_seed_table(sub, mods, agree, title, out_path, stats_footer=""):
    show_pred = "fused_pred" in sub.columns
    cols = ["Patient_ID", "y_true"] + [lbl for lbl, _, _ in mods] + ["agreement"]
    body, colors = [], []
    for i, (_, r) in enumerate(sub.iterrows()):
        row = [str(r["Patient_ID"]), str(r["y_true"])]
        crow = ["white", "white"]
        for lbl, pc, cc in mods:
            pred = r[pc]
            row.append("absent" if pd.isna(pred) else str(pred))
            crow.append(_cell_color(r[cc] if cc else None))
        row.append(str(agree.loc[r.name]))
        crow.append("white")
        body.append(row); colors.append(crow)

    n_rows, n_cols = len(body), len(cols)
    fig, ax = plt.subplots(figsize=(1.6 + 1.5 * n_cols, 0.5 + 0.27 * (n_rows + 2)))
    ax.axis("off")
    tab = ax.table(cellText=body, colLabels=cols, cellColours=colors,
                   loc="upper center", cellLoc="center", colLoc="center")
    tab.auto_set_font_size(False); tab.set_fontsize(8); tab.scale(1.0, 1.25)
    for j in range(n_cols):
        tab[0, j].set_facecolor(HEAD); tab[0, j].set_text_props(weight="bold")
    plt.title(title, pad=30, fontsize=11)
    # red filtered-cohort banner directly under the title
    ax.text(0.5, 1.005, FILTER_BANNER, transform=ax.transAxes, ha="center", va="bottom",
            fontsize=8, color="#b00000", weight="bold", wrap=True)
    legend = ("green = modality correct    red = modality wrong    grey = modality absent\n"
              "cell text = predicted class; rows = TEST patients in the FILTERED failure table\n"
              + (stats_footer or ""))
    ax.text(0.0, -0.04, legend, transform=ax.transAxes, ha="left", va="top",
            fontsize=8, family="monospace", linespacing=1.4)
    fig.savefig(out_path, dpi=170, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"  PNG: {out_path}")


def render_summary(df, agree, title, out_path, stats_footer=""):
    seed_col = "seed" if "seed" in df.columns else None
    patterns = sorted(agree.unique())
    rows = []
    seeds = sorted(df[seed_col].unique()) if seed_col else [None]
    for s in seeds:
        m = (df[seed_col] == s) if seed_col else np.ones(len(df), bool)
        vc = agree[m].value_counts()
        rows.append([f"seed {s}" if s is not None else "all"] + [int(vc.get(p, 0)) for p in patterns])
    vc_all = agree.value_counts()
    rows.append(["OVERALL"] + [int(vc_all.get(p, 0)) for p in patterns])

    cols = ["cohort"] + patterns
    fig, ax = plt.subplots(figsize=(2.0 + 1.7 * len(cols), 0.6 + 0.3 * (len(rows) + 2)))
    ax.axis("off")
    tab = ax.table(cellText=rows, colLabels=cols, loc="upper center",
                   cellLoc="center", colLoc="center")
    tab.auto_set_font_size(False); tab.set_fontsize(9); tab.scale(1.0, 1.4)
    for j in range(len(cols)):
        tab[0, j].set_facecolor(HEAD); tab[0, j].set_text_props(weight="bold")
    for i in range(1, len(rows) + 1):
        tab[i, 0].set_text_props(weight="bold", ha="left")
    plt.title(title, pad=30, fontsize=11)
    ax.text(0.5, 1.01, FILTER_BANNER, transform=ax.transAxes, ha="center", va="bottom",
            fontsize=8, color="#b00000", weight="bold", wrap=True)
    foot = ("Agreement patterns over the FILTERED failure-table rows (both-correct hidden). "
            "'clinical_only' = clinical\nright where MRI wrong; 'mri_only' = MRI right where "
            "clinical wrong (complementary lift); 'both_wrong' = neither.\n\n"
            + (stats_footer or ""))
    ax.text(0.0, -0.06, foot, transform=ax.transAxes, ha="left", va="top", fontsize=8,
            family="monospace", linespacing=1.4)
    fig.savefig(out_path, dpi=170, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"  PNG: {out_path}")


def render_csv(path):
    df = pd.read_csv(path)
    mods = detect_modalities(df.columns)
    if not mods:
        print(f"  [skip] no *_pred columns in {path}")
        return
    agree = derive_agreement(df, mods)
    out_dir = os.path.dirname(path)
    tag = os.path.relpath(out_dir, start=os.path.dirname(os.path.dirname(path)))
    stats_path = os.path.join(out_dir, "failure_stats.csv")
    print(f"[{path}]  modalities: {[m[0] for m in mods]}  rows={len(df)}")
    render_summary(df, agree, f"Failure-mode agreement summary — {tag}",
                   os.path.join(out_dir, "failure_summary.png"),
                   stats_footer=build_stats_footer(stats_path))
    seed_col = "seed" if "seed" in df.columns else None
    seeds = sorted(df[seed_col].unique()) if seed_col else [None]
    for s in seeds:
        sub = df[df[seed_col] == s] if seed_col else df
        render_seed_table(sub, mods, agree, f"Failure table — {tag}  (seed {s})",
                          os.path.join(out_dir, f"failure_table_seed{s}.png"),
                          stats_footer=build_stats_footer(stats_path, seed=s))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", help="A single failure_table.csv to render.")
    ap.add_argument("--dir", help="Render every failure_table.csv found under this dir.")
    args = ap.parse_args()
    targets = []
    if args.csv:
        targets.append(args.csv)
    if args.dir:
        targets += sorted(glob.glob(os.path.join(args.dir, "**", "failure_table.csv"),
                                    recursive=True))
    if not targets:
        ap.error("provide --csv or --dir")
    for t in targets:
        render_csv(t)


if __name__ == "__main__":
    main()
