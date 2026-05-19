"""
19d_head_hyperparams.py
========================
Hyperparameter reference for the four classifier heads used by
``19_compare_frozen_models.py``.

Two tables, B&W house style (DejaVu Serif, no bold), written next to the
other comparison deliverables:

  A. SPEC  — the FIXED hyperparameters per head (architecture + optimiser +
     regularisation + class-imbalance handling + the one tuned hyperparameter
     and its selection rule). Values mirror script 19 / script 15 defaults.
  B. SELECTED — the data-driven tuned value actually chosen per run, summarised
     across every metrics.json ``fit_info`` in the results tree (MLP → early-
     stop epoch; RF → n_estimators picked on val AUC; XGB → best_iteration
     from early stopping).

Outputs (under --root):
  head_hyperparameters.csv, head_hyperparameters.md, head_hyperparameters.png

Usage
-----
  python snp_pipeline/19d_head_hyperparams.py \\
      --root D:/ADNI_SNP_Omni2.5M_20140220/frozen_model_comparison
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.family"] = "DejaVu Serif"

HEAD_PRETTY = {"mlp2": "2-layer MLP", "mlp3": "3-layer MLP",
               "rf": "RandomForest", "xgb": "XGBoost"}
HEAD_ORDER = ["mlp2", "mlp3", "rf", "xgb"]

# ── A. Fixed spec — mirrors snp_pipeline/19_compare_frozen_models.py
#    (argparse defaults: epochs 200, patience 20, batch 64, lr 1e-3, wd 1e-4;
#     MLPHead hidden=256 dropout=0.3 from script 15; MLP3Head h1=256 h2=128
#     dropout=0.3; RF n_estimators grid {200,500}; XGB constructor in
#     run_cell) — single source of truth is the code; keep in sync if changed.
SPEC = [
    {"head": "mlp2", "family": "Neural net",
     "arch": "Linear(D→256) → GELU → Dropout(0.3) → Linear(256→1)",
     "tuned": "training length: ≤200 epochs, early-stop patience 20 on val AUC "
              "(best-val checkpoint restored)",
     "fixed": "AdamW, lr 1e-3, weight_decay 1e-4, batch_size 64",
     "imbalance": "weighted BCEWithLogits, pos_weight = n_neg/n_pos (train)"},
    {"head": "mlp3", "family": "Neural net",
     "arch": "Linear(D→256)→GELU→Drop(0.3)→Linear(256→128)→GELU→Drop(0.3)"
             "→Linear(128→1)",
     "tuned": "training length: ≤200 epochs, early-stop patience 20 on val AUC "
              "(best-val checkpoint restored)",
     "fixed": "AdamW, lr 1e-3, weight_decay 1e-4, batch_size 64",
     "imbalance": "weighted BCEWithLogits, pos_weight = n_neg/n_pos (train)"},
    {"head": "rf", "family": "Bagged trees",
     "arch": "RandomForestClassifier; trees default "
             "(max_depth=None, max_features=sqrt, bootstrap=True)",
     "tuned": "n_estimators ∈ {200, 500} → pick the value with higher val AUC",
     "fixed": "random_state = seed, n_jobs = −1",
     "imbalance": "class_weight = \"balanced\""},
    {"head": "xgb", "family": "Gradient boosting",
     "arch": "XGBClassifier; max_depth 4, learning_rate 0.05, "
             "subsample 0.9, colsample_bytree 0.9",
     "tuned": "n_estimators ≤ 400 with early_stopping_rounds = 30 on val "
              "logloss → best_iteration",
     "fixed": "eval_metric logloss, random_state = seed, n_jobs = −1",
     "imbalance": "scale_pos_weight = n_neg/n_pos (train)"},
]
SPEC_COLS = [("head", "Head"), ("family", "Family"),
             ("arch", "Architecture / key params"),
             ("tuned", "Tuned hyperparameter (search → selection)"),
             ("fixed", "Fixed training / optimiser"),
             ("imbalance", "Class-imbalance handling")]


def selected_summary(root: Path) -> pd.DataFrame:
    """Summarise the data-driven tuned value across every metrics.json."""
    by_head: dict[str, list] = {h: [] for h in HEAD_ORDER}
    n = 0
    for mf in root.rglob("metrics.json"):
        try:
            m = json.load(open(mf))
        except Exception:
            continue
        head = m.get("head") or m.get("clf")
        fi = m.get("fit_info", {})
        if head not in by_head:
            continue
        n += 1
        if head in ("mlp2", "mlp3"):
            by_head[head].append(("n_epochs", fi.get("n_epochs")))
        else:
            by_head[head].append(("n_estimators", fi.get("n_estimators")))
    rows = []
    for h in HEAD_ORDER:
        vals = [v for _, v in by_head[h] if v is not None]
        if not vals:
            rows.append({"head": h, "tuned_param": "—",
                         "n_runs": 0, "selected_values": "(no runs found)"})
            continue
        param = by_head[h][0][0]
        arr = np.array(vals, dtype=float)
        if h == "rf":
            cnt = Counter(int(v) for v in vals)
            desc = ", ".join(f"{k}: {cnt[k]}/{len(vals)} runs"
                             for k in sorted(cnt))
        else:
            desc = (f"min {int(arr.min())}, median {int(np.median(arr))}, "
                    f"mean {arr.mean():.0f}, max {int(arr.max())}")
        rows.append({"head": h, "tuned_param": param,
                     "n_runs": len(vals), "selected_values": desc})
    return pd.DataFrame(rows)


def _render_png(spec_df: pd.DataFrame, sel_df: pd.DataFrame, out_png: Path):
    """Hard-wrap every cell with textwrap to its column width and use
    variable per-row heights so nothing overlaps between columns."""
    import textwrap

    FS, CPI = 8.0, 6.3                       # font size; ~chars per inch
    LINE_IN, VPAD_IN = 0.18, 0.16            # text line pitch / cell v-padding
    HDR_IN, SEC_IN, TITLE_IN, MARGIN_IN = 0.34, 0.44, 0.62, 0.30
    fig_w = 19.0

    A_w = [0.06, 0.08, 0.30, 0.24, 0.16, 0.16]   # head/family/arch/tuned/fixed/imbalance
    B_w = [0.12, 0.17, 0.08, 0.63]               # head/param/n/values

    def wrap_cell(s, frac):
        width = max(6, int((frac * fig_w - 0.16) * CPI))
        out: list[str] = []
        for para in str(s).split("\n"):
            out += textwrap.wrap(para, width=width) or [""]
        return out

    A_hdr = [n for _, n in SPEC_COLS]
    A_rows, A_h = [], []
    for _, row in spec_df.iterrows():
        cells = []
        for (key, _), frac in zip(SPEC_COLS, A_w):
            v = (HEAD_PRETTY.get(row["head"], row["head"])
                 if key == "head" else str(row[key]))
            cells.append(wrap_cell(v, frac))
        A_rows.append(cells)
        A_h.append(max(len(c) for c in cells) * LINE_IN + VPAD_IN)

    B_hdr = ["Head", "Tuned parameter", "n runs", "Selected values"]
    B_rows, B_h = [], []
    for _, row in sel_df.iterrows():
        vals = [HEAD_PRETTY.get(row["head"], row["head"]),
                str(row["tuned_param"]), str(row["n_runs"]),
                str(row["selected_values"])]
        cells = [wrap_cell(v, f) for v, f in zip(vals, B_w)]
        B_rows.append(cells)
        B_h.append(max(len(c) for c in cells) * LINE_IN + VPAD_IN)

    fig_h = (TITLE_IN + SEC_IN + HDR_IN + sum(A_h)
             + SEC_IN + HDR_IN + sum(B_h) + 2 * MARGIN_IN)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    def Y(v):
        return v / fig_h

    def box(x, y_in, w, h_in, lw=0.4):
        ax.add_patch(plt.Rectangle((x, Y(y_in)), w, h_in / fig_h,
                     transform=ax.transAxes, facecolor="white",
                     edgecolor="black", linewidth=lw))

    def ctxt(x, y_in, s, fs, ha="center"):
        ax.text(x, Y(y_in), s, ha=ha, va="center", fontsize=fs,
                transform=ax.transAxes)

    def cell(x, y_top_in, lines, fs=FS):
        ax.text(x + 0.004, Y(y_top_in - 0.02), "\n".join(lines),
                ha="left", va="top", fontsize=fs, transform=ax.transAxes,
                linespacing=1.3)

    cur = fig_h - MARGIN_IN
    ctxt(0.5, cur - TITLE_IN * 0.5,
         "Classifier-head hyperparameters  "
         "(snp_pipeline/19_compare_frozen_models.py)", 14)
    cur -= TITLE_IN

    # ── Table A ─────────────────────────────────────────────────────
    ctxt(0.0, cur - SEC_IN * 0.55, "A.  Fixed specification per head",
         11, ha="left")
    cur -= SEC_IN
    xA = [0.0]
    for w in A_w:
        xA.append(xA[-1] + w)
    for i, name in enumerate(A_hdr):
        box(xA[i], cur - HDR_IN, A_w[i], HDR_IN, lw=0.7)
        ctxt(xA[i] + A_w[i] / 2, cur - HDR_IN / 2, name, 9)
    cur -= HDR_IN
    for cells, rh in zip(A_rows, A_h):
        for ci in range(len(A_w)):
            box(xA[ci], cur - rh, A_w[ci], rh)
            cell(xA[ci], cur, cells[ci])
        cur -= rh

    # ── Table B ─────────────────────────────────────────────────────
    ctxt(0.0, cur - SEC_IN * 0.55,
         "B.  Data-driven value actually selected (across all "
         f"{int(sel_df['n_runs'].sum())} runs)", 11, ha="left")
    cur -= SEC_IN
    xB = [0.0]
    for w in B_w:
        xB.append(xB[-1] + w)
    for i, name in enumerate(B_hdr):
        box(xB[i], cur - HDR_IN, B_w[i], HDR_IN, lw=0.7)
        ctxt(xB[i] + B_w[i] / 2, cur - HDR_IN / 2, name, 9)
    cur -= HDR_IN
    for cells, rh in zip(B_rows, B_h):
        for ci in range(len(B_w)):
            box(xB[ci], cur - rh, B_w[ci], rh)
            cell(xB[ci], cur, cells[ci])
        cur -= rh

    fig.savefig(out_png, dpi=150, facecolor="white",
                bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print(f"Wrote {out_png}  (figsize={fig_w} x {fig_h:.1f} in)")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path,
                    default=Path("D:/ADNI_SNP_Omni2.5M_20140220/frozen_model_comparison"))
    args = ap.parse_args()

    spec_df = pd.DataFrame(SPEC)[[c for c, _ in SPEC_COLS]]
    sel_df = selected_summary(args.root)

    csv = args.root / "head_hyperparameters.csv"
    spec_df.assign(head=spec_df["head"].map(HEAD_PRETTY)).to_csv(csv, index=False)
    print(f"Wrote {csv}  ({len(spec_df)} heads)")

    md = args.root / "head_hyperparameters.md"
    lines = ["# Classifier-head hyperparameters", "",
             "Heads used by `snp_pipeline/19_compare_frozen_models.py`. "
             "Per-window pool = mean; aggregation feature = "
             "`chrom_mean_then_mean` (all heads) and `chrom_mean_then_attn` "
             "(MLP heads only).", "", "## A. Fixed specification", ""]
    lines.append("| " + " | ".join(n for _, n in SPEC_COLS) + " |")
    lines.append("|" + "|".join(["---"] * len(SPEC_COLS)) + "|")
    for _, r in spec_df.iterrows():
        cells = [(HEAD_PRETTY.get(r["head"], r["head"]) if k == "head"
                  else str(r[k])) for k, _ in SPEC_COLS]
        lines.append("| " + " | ".join(cells) + " |")
    lines += ["", "## B. Data-driven value actually selected", "",
              "| Head | Tuned parameter | n runs | Selected values |",
              "|---|---|---|---|"]
    for _, r in sel_df.iterrows():
        lines.append(f"| {HEAD_PRETTY.get(r['head'], r['head'])} | "
                     f"{r['tuned_param']} | {r['n_runs']} | "
                     f"{r['selected_values']} |")
    md.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {md}")

    _render_png(spec_df, sel_df, args.root / "head_hyperparameters.png")

    print("\n=== Fixed spec ===")
    print(spec_df.to_string(index=False))
    print("\n=== Selected (data-driven) ===")
    print(sel_df.to_string(index=False))


if __name__ == "__main__":
    main()
