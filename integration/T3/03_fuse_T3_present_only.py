"""
03_fuse_T3_present_only.py — T3 per-horizon fusion, PRESENT-ONLY cohort
=======================================================================
The base-engine M12 framing (01_fuse_T3.py) is COMPLETE-CASE: it intersects to
the m12-MRI cohort and drops clinical-test patients without an m12 scan (too
restrictive). This evaluates the same per-horizon binary fusion on the FULL
clinical-test cohort: where the m12 MRI is absent, the single present modality
(clinical) proceeds — mirroring the T1d present-only approach.

Per patient, present-only blend: w·CL + (1-w)·MRI where both present, else CL
alone (w tuned on VAL). clinical-only/EN/LR are scored on the full cohort
(missing MRI uniform-imputed + present flag for EN/LR); mri-only on the
MRI-present subset (reference). Reuses the engine `01_fuse_T2.py` loaders +
fit + metrics.

Per-horizon clinical encoder (full_ft) + MRI source:
  3yr  BioClinical-ModernBERT-large / T3a_conv3y   MRI BrainMVP full_ft plus_original (T3a)
  5yr  BioClinical-ModernBERT-base  / T3b_conv5y   MRI BrainMVP frozen/none cached    (T3b; no full_ft)
  7yr  ModernBERT-large             / T3d_conv7y   MRI BrainMVP full_ft plus_original (T3c)

Run:  conda run -n snp python integration/T3/03_fuse_T3_present_only.py
Outputs: integration/T3/SUMMARY_T3_present_only.{png,csv}
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import balanced_accuracy_score as bacc

HERE = Path(__file__).resolve().parent
ENGINE = HERE.parent / "T2_classification" / "01_fuse_T2.py"
_spec = importlib.util.spec_from_file_location("_fuse_engine_t3", ENGINE)
h = importlib.util.module_from_spec(_spec); sys.modules["_fuse_engine_t3"] = h
_spec.loader.exec_module(h)
h._set_classes(["noConv", "conv"])                         # binary
K = len(h.CLASS_NAMES)

DERIV = Path("D:/ADNI_BIDS_project/derivatives")
CLIN_ROOT = DERIV / "encoder_outputs_no_cdr_post_exclusion"
MRI_ROOT = DERIV / "mri_embeddings"
BMVP_FT_ROOT = DERIV / "brainmvp_embeddings"
SEEDS = [0, 1, 2]

HORIZONS = {
    "3yr": dict(model="BioClinical-ModernBERT-large", clin_task="T3a_conv3y", mri_task="T3a_conv3y",
                mri="ft",     mri_label="BrainMVP full_ft plus_original", label="T3 conv<=3y"),
    "5yr": dict(model="BioClinical-ModernBERT-base",  clin_task="T3b_conv5y", mri_task="T3b_conv5y",
                mri="cached", mri_label="BrainMVP frozen/none",          label="T3 conv<=5y"),
    "7yr": dict(model="ModernBERT-large",             clin_task="T3d_conv7y", mri_task="T3c_conv7y",
                mri="ft",     mri_label="BrainMVP full_ft plus_original", label="T3 conv<=7y"),
}
ORDER = ["clinical_only", "mri_only", "equal_0.5", "weighted_avg", "weighted_avg_J", "elastic_net", "lr"]
VAR_LABEL = {"clinical_only": "clinical-only", "mri_only": "MRI-only (subset)",
             "equal_0.5": "equal 0.5", "weighted_avg": "weighted-avg bACC",
             "weighted_avg_J": "weighted-avg Youden", "elastic_net": "elastic-net", "lr": "logistic-reg"}


def render_horizon_present_only(per_seed, hz, cfg, out_path):
    """Per-horizon PRESENT-ONLY table in the engine's fusion_table_M12 layout
    (rows=variants, cols=bACC/macro-F1/macro-AUC/F1 noConv/F1 conv), saved as
    fusion_table_M12_present_only.png next to the complete-case engine table."""
    METR = [("bacc", "bACC"), ("macro_f1", "macro-F1"), ("macro_auc", "macro-AUC"),
            ("f1_noConv", "F1 noConv"), ("f1_conv", "F1 conv")]
    po = per_seed[(per_seed.horizon == hz) & (per_seed.cohort == "present_only")]
    mri = per_seed[(per_seed.horizon == hz) & (per_seed.cohort == "both_present")
                   & (per_seed.variant == "mri_only")]
    lab = {"clinical_only": f"clinical-only [{cfg['model']}/full_ft]",
           "mri_only": f"MRI-only [{cfg['mri_label']}] (MRI-present subset)",
           "equal_0.5": "equal-0.5", "weighted_avg": "weighted-avg (bACC-tuned)",
           "weighted_avg_J": "weighted-avg (Youden's J)",
           "elastic_net": "elastic-net (present-only)", "lr": "logistic regression"}
    order = ["clinical_only", "mri_only", "equal_0.5", "weighted_avg", "weighted_avg_J", "elastic_net", "lr"]
    body, baccs = [], []
    for v in order:
        src = mri if v == "mri_only" else po[po.variant == v]
        if src.empty:
            continue
        row = [lab[v]]
        for mk, _ in METR:
            m, s, n = src[mk].mean(), src[mk].std(), src["n"].mean()
            row.append(f"{m:.3f} ± {s:.3f}" + (f" ({n:.0f})" if mk == "bacc" else ""))
        body.append(row)
        baccs.append((v, src["bacc"].mean()))
    if not body:
        return
    cols = ["Variant (fusion method / model)"] + [lbl for _, lbl in METR]
    nC = len(cols)
    cw = [max(len(str(x)) for x in [cols[j]] + [b[j] for b in body]) + 2 for j in range(nC)]
    cw[0] += 6
    tot = sum(cw)
    fig, ax = plt.subplots(figsize=(0.115 * tot, 0.55 + 0.42 * (len(body) + 1)))
    ax.axis("off")
    tab = ax.table(cellText=body, colLabels=cols, loc="upper center", cellLoc="center",
                   colLoc="center", colWidths=[c / tot for c in cw])
    tab.auto_set_font_size(False); tab.set_fontsize(9); tab.scale(1.0, 1.5)
    for j in range(nC):
        tab[0, j].set_facecolor("#ECECEC"); tab[0, j].set_text_props(weight="bold")
    for i in range(1, len(body) + 1):
        tab[i, 0].set_text_props(ha="left")
    # bold best bACC among full-cohort variants (exclude MRI-only subset)
    cand = [(i, b) for i, (v, b) in enumerate(baccs) if v != "mri_only" and not np.isnan(b)]
    if cand:
        bi = max(cand, key=lambda t: t[1])[0]
        tab[bi + 1, 1].set_text_props(weight="bold")
    plt.title(f"T3 {cfg['label']} late fusion -- M12 (CL_bl + MRI_m12) — PRESENT-ONLY cohort\n"
              f"clinical={cfg['model']}/full_ft   MRI={cfg['mri_label']}   (mean ± std across seeds, n)",
              pad=8, fontsize=11)
    foot = ["PRESENT-ONLY: full clinical-test cohort; where the m12 MRI is absent the single present modality",
            "(clinical) proceeds (weighted-avg falls back to clinical; EN/LR uniform-impute MRI + present flag).",
            "MRI-only is scored on the MRI-present subset only (smaller n) — reference. (n) after bACC = mean #TEST",
            "patients across seeds. Compare with the complete-case fusion_table_M12.png (both-present only)."]
    ax.text(0.0, -0.04, "\n".join(foot), transform=ax.transAxes, ha="left", va="top",
            fontsize=8, family="monospace", linespacing=1.4)
    fig.savefig(out_path, dpi=180, bbox_inches="tight", pad_inches=0.06); plt.close(fig)
    print(f"  wrote {out_path}")


def _mri_template(cfg):
    t = cfg["mri_task"]
    if cfg["mri"] == "ft":
        return str(BMVP_FT_ROOT / t / "aug_plus_original" / "seed_{seed}" / "embeddings_seed_{seed}.csv")
    return str(MRI_ROOT / t / "brainmvp_frozen_none_cached" / "embeddings_seed_{seed}.csv")


def _probs(df, cols):
    return df[cols].to_numpy(float)


def main():
    rows, all_preds = [], []
    for hz, cfg in HORIZONS.items():
        tmpl = _mri_template(cfg)
        if not Path(tmpl.format(seed=0)).exists():
            print(f"[skip] {hz}: no MRI CSV at {tmpl.format(seed=0)}"); continue
        for seed in SEEDS:
            clin = h.load_clinical(seed, CLIN_ROOT, cfg["model"], "full_ft", cfg["clin_task"])
            mri = h.load_mri(seed, tmpl, "m12")
            f = clin.merge(mri, on="Patient_ID", how="left")     # keep ALL clinical patients
            cp = _probs(f, h.CP_COLS)
            mp = _probs(f, h.MP_COLS)                             # NaN where MRI absent
            pres = ~np.isnan(mp).any(axis=1)
            y = f["y_clin"].to_numpy(int)
            vm, tm = (f["split"] == "val").to_numpy(), (f["split"] == "test").to_numpy()
            yv, yt = y[vm], y[tm]
            cpv, cpt, mpv, mpt = cp[vm], cp[tm], mp[vm], mp[tm]
            prv, prt = pres[vm], pres[tm]
            pid_t = f["Patient_ID"].to_numpy()[tm]
            if len(yt) == 0 or len(yv) == 0:
                continue

            def fuse_po(cl, mr, p, w):
                out = cl.copy()
                if p.any():
                    out[p] = w * cl[p] + (1.0 - w) * mr[p]
                return out

            def tune_w(obj):
                scorer = h.youden_j if obj == "youden" else bacc
                bw, bs = 0.5, -2.0
                for w in np.linspace(0, 1, 21):
                    s = scorer(yv, fuse_po(cpv, mpv, prv, w).argmax(1))
                    if s > bs:
                        bs, bw = s, w
                return bw

            def Xpo(cl, mr, p):
                mrf = np.where(p[:, None], np.nan_to_num(mr), 1.0 / K)
                return np.hstack([cl, mrf, p.astype(float)[:, None]])

            def rec(variant, cohort, yy, pred, proba, vb, pids):
                r = h.metric_row(np.asarray(yy, int), np.asarray(pred, int), proba)
                r.update(horizon=hz, seed=seed, variant=variant, cohort=cohort, val_bacc=float(vb),
                         clin_model=cfg["model"], mri_model=cfg["mri_label"])
                rows.append(r)
                for i in range(len(yy)):
                    all_preds.append(dict(horizon=hz, seed=seed, variant=variant, cohort=cohort,
                                          Patient_ID=pids[i], y_true=int(yy[i]), pred=int(pred[i])))

            def add(variant, yy, pred, proba, vb, pids, pres):
                yy, pred, proba, pids = np.asarray(yy), np.asarray(pred), np.asarray(proba), np.asarray(pids)
                rec(variant, "present_only", yy, pred, proba, vb, pids)
                if pres is not None and pres.sum() > 0:
                    rec(variant, "both_present", yy[pres], pred[pres], proba[pres], vb, pids[pres])

            add("clinical_only", yt, cpt.argmax(1), cpt, bacc(yv, cpv.argmax(1)), pid_t, prt)
            if prt.any():
                vb = bacc(yv[prv], mpv[prv].argmax(1)) if prv.sum() > 1 else float("nan")
                rec("mri_only", "both_present", yt[prt], mpt[prt].argmax(1), mpt[prt], vb, pid_t[prt])
            et = fuse_po(cpt, mpt, prt, 0.5)
            add("equal_0.5", yt, et.argmax(1), et, bacc(yv, fuse_po(cpv, mpv, prv, 0.5).argmax(1)), pid_t, prt)
            wb = tune_w("bacc"); wt = fuse_po(cpt, mpt, prt, wb)
            add("weighted_avg", yt, wt.argmax(1), wt, bacc(yv, fuse_po(cpv, mpv, prv, wb).argmax(1)), pid_t, prt)
            wj = tune_w("youden"); wjt = fuse_po(cpt, mpt, prt, wj)
            add("weighted_avg_J", yt, wjt.argmax(1), wjt, bacc(yv, fuse_po(cpv, mpv, prv, wj).argmax(1)), pid_t, prt)
            Xv, Xt = Xpo(cpv, mpv, prv), Xpo(cpt, mpt, prt)
            for name, fn in [("elastic_net", h.fit_en), ("lr", h.fit_lr)]:
                res = fn(Xv, yv, Xt, seed)
                if res is not None:
                    add(name, yt, res[0], res[1], float("nan"), pid_t, prt)
            print(f"  {hz} seed{seed}: n_test={len(yt)} (MRI present {int(prt.sum())})  "
                  f"clin={bacc(yt, cpt.argmax(1)):.3f}  wavg={bacc(yt, wt.argmax(1)):.3f}")

    if not rows:
        sys.exit("no rows")
    per_seed = pd.DataFrame(rows)
    per_seed.to_csv(HERE / "SUMMARY_T3_present_only_per_seed.csv", index=False)
    ap = pd.DataFrame(all_preds)

    # per-horizon PRESENT-ONLY tables (engine layout) beside the complete-case fusion_table_M12.png
    for hz, cfg in HORIZONS.items():
        outd = HERE / hz / ("brainmvp_ft" if cfg["mri"] == "ft" else "brainmvp")
        if outd.is_dir():
            render_horizon_present_only(per_seed, hz, cfg, outd / "fusion_table_M12_present_only.png")

    # aggregate: seed-mean + pooled bACC per (horizon, variant)
    body = []
    for hz in HORIZONS:
        hd = per_seed[per_seed.horizon == hz]
        if hd.empty:
            continue
        for ckey, clabel in [("both_present", "both-present"), ("present_only", "present-only")]:
            cd = hd[hd.cohort == ckey]
            for v in ORDER:
                vd = cd[cd.variant == v]
                if vd.empty:
                    continue
                sm, ss, au, n = vd.bacc.mean(), vd.bacc.std(), vd.macro_auc.mean(), vd.n.mean()
                pe = ap[(ap.horizon == hz) & (ap.variant == v) & (ap.cohort == ckey)]
                pool = bacc(pe.y_true, pe.pred) if len(pe) and (pe.pred >= 0).all() else float("nan")
                body.append([HORIZONS[hz]["label"], clabel, VAR_LABEL[v],
                             f"{sm:.3f} ± {ss:.3f}", f"{pool:.3f}", f"{au:.3f}", f"{n:.0f}"])

    cols = ["Horizon", "Cohort", "Method", "TEST bACC (seed-mean)", "pooled bACC", "macro-AUC", "n"]
    nC = len(cols)
    cw = [max(len(str(x)) for x in [cols[j]] + [b[j] for b in body]) + 2 for j in range(nC)]
    cw[2] += 6
    tot = sum(cw)
    fig, ax = plt.subplots(figsize=(0.12 * tot, 0.6 + 0.38 * (len(body) + 1)))
    ax.axis("off")
    tab = ax.table(cellText=body, colLabels=cols, loc="upper center", cellLoc="center",
                   colLoc="center", colWidths=[c / tot for c in cw])
    tab.auto_set_font_size(False); tab.set_fontsize(8.5); tab.scale(1.0, 1.4)
    for j in range(nC):
        tab[0, j].set_text_props(weight="bold")
    for i in range(1, len(body) + 1):
        for j in (0, 1, 2):
            tab[i, j].set_text_props(ha="left")
    plt.title("T3 per-horizon late fusion — BOTH cohorts (clinical@bl ⊕ MRI@m12)\n"
              "both-present (complete-case) vs present-only (full clinical cohort; single modality proceeds "
              "where the other is absent)\nmean ± std across seeds (fit VAL / report TEST)", pad=10, fontsize=10.5)
    foot = ["COHORTS: both-present = only patients with an m12 MRI (smaller); present-only = FULL clinical cohort",
            "  (where MRI absent, clinical proceeds; weighted-avg falls back to clinical; EN/LR uniform-impute MRI + flag).",
            "MRI-only exists only on the MRI-present subset (both-present) — reference. pooled bACC = over all test",
            "  patients pooled across seeds. MRI: 3y/7y BrainMVP full_ft plus_original, 5y BrainMVP frozen/none (no",
            "  full_ft for T3b). Clinical: 3y BioClin-L, 5y BioClin-B, 7y ModernBERT-L (all full_ft)."]
    ax.text(0.0, -0.04, "\n".join(foot), transform=ax.transAxes, ha="left", va="top",
            fontsize=8, family="monospace", linespacing=1.4)
    png = HERE / "SUMMARY_T3_fusion_by_cohort.png"
    fig.savefig(png, dpi=180, bbox_inches="tight", pad_inches=0.06); plt.close(fig)
    print(f"wrote {png}")


if __name__ == "__main__":
    main()
