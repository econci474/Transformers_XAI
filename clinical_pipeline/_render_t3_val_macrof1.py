"""
_render_t3_val_macrof1.py
=========================
Regenerate the clinical T3 VALIDATION table (table_val_t3at3bt3d.tex) WITH a true Macro-F1 column per
horizon (T3a ≤3y / T3b ≤5y / T3c ≤7y). No tracked generator existed for this .tex, and every stored
"F1" is positive-class — so macro-F1 is computed per seed here:
  - Encoders: exact, from the val confusion_matrix {tn,fp,fn,tp} in each metrics.json.
  - Tabular baselines: reconstructed per seed from (Accuracy, Sensitivity, Specificity) in
    results_per_seed.csv  ->  prevalence -> fractional confusion -> F1_pos / F1_neg -> macro-F1.
Bal.Acc and AUC are recomputed from the same per-seed sources (match the existing table). mean ± std
(ddof=1) across seeds 0/1/2, 2 dp.

Out: clinical_pipeline/outputs/encoder_only_post_exclusions/table_val_t3at3bt3d.tex  (regenerated)
Run: PYTHONIOENCODING=utf-8 python clinical_pipeline/_render_t3_val_macrof1.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE / "outputs" / "encoder_only_post_exclusions"
ENC_ROOT = Path("D:/ADNI_BIDS_project/derivatives/encoder_outputs_no_cdr_post_exclusion")
SEEDS = [0, 1, 2]

# horizons: (display title, encoder task, baseline task string)
HORIZONS = [
    ("T3a: AD conversion $\\leq$ 3 yrs", "T3a_conv3y", "Conversion to AD within 3 years"),
    ("T3b: AD conversion $\\leq$ 5 yrs", "T3b_conv5y", "Conversion to AD within 5 years"),
    ("T3c: AD conversion $\\leq$ 7 yrs", "T3d_conv7y", "Conversion to AD within 7 years"),
]

# tabular baselines (results_per_seed.csv Model names) -> display
BASELINES = [("LogReg", "Log.\\ Reg."), ("SVM", "SVM"), ("XGBoost", "XGBoost")]
# encoders: (metrics.json ModelDir, results display, strategy key, strategy display)
ENC_MODELS = [("ModernBERT-base", "ModernBERT-B"), ("ModernBERT-large", "ModernBERT-L"),
              ("BioClinical-ModernBERT-base", "BioClin-MBERT-B"),
              ("BioClinical-ModernBERT-large", "BioClin-MBERT-L")]
STRATS = [("frozen", "frozen"), ("full_ft", "full fine-tune")]


def _f1(a, bc):
    """F1 for a class with TP-equivalent `a` and the two error masses summed in `bc`."""
    denom = 2 * a + bc
    return (2 * a / denom) if denom > 0 else 0.0


def macro_f1_from_conf(tn, fp, fn, tp):
    f1_pos = _f1(tp, fp + fn)
    f1_neg = _f1(tn, fn + fp)
    return 0.5 * (f1_pos + f1_neg)


def macro_f1_from_rates(acc, sens, spec, f1_pos_stored):
    """Reconstruct per-seed macro-F1 for a binary baseline from accuracy/sens/spec (fractional)."""
    if abs(sens - spec) < 1e-6:                       # degenerate / near-constant predictor
        return 0.5 * (f1_pos_stored + 0.0) if np.isnan(acc) else _safe_macro(0.5, sens, spec, f1_pos_stored)
    p = (acc - spec) / (sens - spec)                  # prevalence (fraction positive)
    p = min(max(p, 0.0), 1.0)
    return _safe_macro(p, sens, spec, f1_pos_stored)


def _safe_macro(p, sens, spec, f1_pos_stored):
    tp, fn = sens * p, (1 - sens) * p
    tn, fp = spec * (1 - p), (1 - spec) * (1 - p)
    f1_pos = _f1(tp, fp + fn)
    f1_neg = _f1(tn, fn + fp)
    if f1_pos_stored is not None and not np.isnan(f1_pos_stored) and f1_pos > 0:
        # prefer the exact stored positive-class F1 (avoids prevalence rounding on the positive side)
        f1_pos = f1_pos_stored
    return 0.5 * (f1_pos + f1_neg)


def agg(vals):
    v = np.array([x for x in vals if x is not None and not np.isnan(x)], float)
    if len(v) == 0:
        return None
    return v.mean(), (v.std(ddof=1) if len(v) > 1 else 0.0)


# ── encoder rows: read val metrics.json per (model, strategy, horizon, seed) ──────
def encoder_metrics(model_dir, strat, enc_task):
    baccs, aucs, mf1s = [], [], []
    for s in SEEDS:
        jf = ENC_ROOT / model_dir / enc_task / f"seed_{s}" / strat / "metrics.json"
        if not jf.exists():
            continue
        vm = json.loads(jf.read_text()).get("val_metrics", {})
        cm = vm.get("confusion_matrix")
        baccs.append(vm.get("balanced_acc"))
        aucs.append(vm.get("auc_roc"))
        mf1s.append(macro_f1_from_conf(cm["tn"], cm["fp"], cm["fn"], cm["tp"]) if cm else None)
    return agg(baccs), agg(aucs), agg(mf1s)


# ── baseline rows: per-seed val from results_per_seed.csv ─────────────────────────
def baseline_metrics(ps, model, bl_task):
    g = ps[(ps.Model == model) & (ps.Task == bl_task) & (ps.Split == "val")]
    baccs, aucs, mf1s = [], [], []
    for _, r in g.iterrows():
        baccs.append(float(r["BalancedAcc"]))
        aucs.append(float(r["AUC-ROC"]))
        mf1s.append(macro_f1_from_rates(float(r["Accuracy"]), float(r["Sensitivity"]),
                                        float(r["Specificity"]), float(r["F1"])))
    return agg(baccs), agg(aucs), agg(mf1s)


def cell(stat, bold=False):
    if stat is None:
        return "--"
    txt = f"{stat[0]:.2f} \\pm {stat[1]:.2f}"
    return f"$\\mathbf{{{txt}}}$" if bold else f"${txt}$"


def main():
    ps = pd.read_csv(OUT_DIR / "results_per_seed.csv")

    # rows: list of (model_disp, strat_disp, [per-horizon (bacc, auc, mf1)]), grouped by family
    families = []   # (family_rows, ...)
    bl_rows = []
    for model, disp in BASELINES:
        cells = [baseline_metrics(ps, model, h[2]) for h in HORIZONS]
        bl_rows.append((disp, "tabular", cells))
    families.append(bl_rows)
    for strat, sdisp in STRATS:
        rows = []
        for mdir, disp in ENC_MODELS:
            cells = [encoder_metrics(mdir, strat, h[1]) for h in HORIZONS]
            rows.append((disp, sdisp, cells))
        families.append(rows)

    # bold = best Bal.Acc per window within each family
    bold_idx = {}
    for fi, fam in enumerate(families):
        for hj in range(len(HORIZONS)):
            best, bri = -1.0, None
            for ri, (_, _, cells) in enumerate(fam):
                st = cells[hj][0]
                if st is not None and st[0] > best:
                    best, bri = st[0], ri
            bold_idx[(fi, hj)] = bri

    L = [r"% Clinical T3 validation table with Macro-F1 (regenerated by _render_t3_val_macrof1.py)",
         r"% Macro-F1 = mean of per-class F1 (encoders: exact from val confusion matrix; "
         r"baselines: reconstructed from accuracy/sens/spec).",
         r"\begin{table}[ht]", r"\centering", r"\resizebox{\textwidth}{!}{%", r"\normalsize",
         r"\begin{tabular}{ll *{9}{c}}", r"\toprule",
         " & & " + " & ".join(r"\multicolumn{3}{c}{\textbf{" + h[0] + "}}" for h in HORIZONS) + r" \\",
         r"\cmidrule(lr){3-5}\cmidrule(lr){6-8}\cmidrule(lr){9-11}",
         r"\textbf{Model} & \textbf{Strategy} & "
         + " & ".join([r"\textbf{Bal.Acc} & \textbf{AUC} & \textbf{Macro-F1}"] * 3) + r" \\",
         r"\midrule"]
    for fi, fam in enumerate(families):
        if fi > 0:
            L.append(r"\midrule")
        for ri, (mdisp, sdisp, cells) in enumerate(fam):
            row = [mdisp, sdisp]
            for hj, (bacc, auc, mf1) in enumerate(cells):
                row.append(cell(bacc, bold=(bold_idx[(fi, hj)] == ri)))
                row.append(cell(auc))
                row.append(cell(mf1))
            L.append(" & ".join(row) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}%", r"}",
          r"\caption{Validation performance for AD conversion prediction across three time windows, "
          r"with \textbf{Macro-F1} (mean of the two per-class F1 scores). T3a ($\leq$3 yrs): NC=37, C=4, "
          r"N(val)=41. T3b ($\leq$5 yrs): NC=22, C=8, N(val)=30. T3c ($\leq$7 yrs): NC=18, C=9, "
          r"N(val)=27. Bold = best Bal.Acc per window within each model family. Results reported as "
          r"mean $\pm$ std across seeds. B=Base, L=Large.}",
          r"\label{tab:clinical_t3at3bt3c_val}", r"\end{table}"]

    out = OUT_DIR / "table_val_t3at3bt3d.tex"
    out.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"  TEX: {out}  ({len(out.read_text(encoding='utf-8').splitlines())} lines)")


if __name__ == "__main__":
    main()
