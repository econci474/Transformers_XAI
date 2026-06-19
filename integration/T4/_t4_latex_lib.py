"""
_t4_latex_lib.py
================
Shared LaTeX writer for the two T4 conversion late-fusion summary tables (04_summary_t4_cleaned.py and
05_fuse_T4_braindino_direct.py). Splits the packed `Method` column into Clinical | MRI | Aggregation and
drops the three per-bucket F1 columns; keeps Cohort grouping + TEST/pooled bACC, macro-AUC, macro-F1.
Portrait booktabs (\\multicolumn title, bold-best per cohort excl. MRI-only, bottom caption).
"""
import pandas as pd

_AGG = {"equal_0.5": "weighted-avg (equal)", "weighted_avg": "weighted-avg (val bACC)",
        "weighted_avg_J": "weighted-avg (Youden J)", "elastic_net": "elastic-net", "lr": "logistic-reg"}


def _esc(s):
    s = str(s)
    for a, b in [("≤", r"$\leq$"), ("≥", r"$\geq$"), ("±", r"$\pm$"), ("·", r"$\cdot$"),
                 ("⊕", r"$\oplus$"), ("−", r"$-$"), ("≈", r"$\approx$"), ("×", r"$\times$"),
                 ("–", "--"), ("—", "--"), ("<", r"$<$"), (">", r"$>$"),
                 ("&", r"\&"), ("%", r"\%"), ("_", r"\_"), ("#", r"\#")]:
        s = s.replace(a, b)
    return s


def _decompose(variant, clin_model, mri_model):
    if variant == "clinical_only":
        return clin_model, "—", "—"
    if variant == "mri_only":
        return "—", mri_model, "—"
    return clin_model, mri_model, _AGG.get(variant, variant)


def write_t4_latex(df, out_path, title1, caption_body, clin_model, mri_model, label):
    """df cols: ckey, Cohort, variant, bacc, bstd, pooled, macro_auc, macro_f1 (+ f1_h_* ignored)."""
    METRICS = [("bacc", "TEST bACC", True), ("pooled", "pooled bACC", False),
               ("macro_auc", "macro-AUC", False), ("macro_f1", "macro-F1", False)]
    # bold-best per cohort, excluding mri_only (different patient set)
    best = {}
    for ckey, g in df.groupby("ckey", sort=False):
        gf = g[g.variant != "mri_only"]
        for key, _, _ in METRICS:
            best[(ckey, key)] = gf[key].idxmax() if gf[key].notna().any() else -1

    def cell(val, std, bold):
        if pd.isna(val):
            return "--"
        txt = f"{val:.3f} \\pm {std:.3f}" if std is not None else f"{val:.3f}"
        return f"$\\mathbf{{{txt}}}$" if bold else f"${txt}$"

    head = [r"\textbf{Cohort}", r"\textbf{Clinical}", r"\textbf{MRI}", r"\textbf{Aggregation}"] \
        + [r"\textbf{" + lbl + "}" for _, lbl, _ in METRICS]
    L = [r"% " + title1,
         r"% Requires: \usepackage{booktabs, graphicx}",
         r"\begin{table}[ht]", r"\centering",
         r"\resizebox{\textwidth}{!}{%", r"\normalsize",
         r"\begin{tabular}{llll cccc}", r"\toprule",
         r"\multicolumn{8}{c}{\textbf{" + _esc(title1) + r"}} \\", r"\midrule",
         " & ".join(head) + r" \\", r"\midrule"]
    prev_c = None
    for idx, r in df.iterrows():
        if prev_c is not None and r.ckey != prev_c:
            L.append(r"\midrule")
        clin, mri, agg = _decompose(r.variant, clin_model, mri_model)
        coh = _esc(r.Cohort) if r.ckey != prev_c else ""
        cells = [coh, _esc(clin), _esc(mri), _esc(agg)]
        for key, _, has_std in METRICS:
            cells.append(cell(r[key], r["bstd"] if has_std else None, best[(r.ckey, key)] == idx))
        L.append(" & ".join(cells) + r" \\")
        prev_c = r.ckey
    L += [r"\bottomrule", r"\end{tabular}%", r"}",
          r"\caption{" + _esc(caption_body) + r"}", r"\label{" + label + r"}", r"\end{table}"]
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print(f"  TEX: {out_path}")
