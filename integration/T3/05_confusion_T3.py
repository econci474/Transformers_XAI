"""
05_confusion_T3.py
==================
3-panel confusion-matrix figures (T2/T1d style: counts pooled across seeds, row-normalised %) for the
TOP T3 integration method per horizon — <=3y and <=5y (the <=7y horizon has no fusion that beats
clinical-only, so it is omitted). Panels: clinical-only (full cohort), MRI-only (MRI-present subset),
best fusion = weighted-avg (equal, w=0.5; CL@bl ⊕ MRI@m12, full cohort with CL-only fallback where MRI
absent). Reuses the 03 loaders/fusion (no metrics re-fit needed).

Out: integration/T3/confusion_T3_{3y,5y}.png (+ .pdf)
Run: python integration/T3/05_confusion_T3.py
"""
import importlib.util
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, balanced_accuracy_score

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("t3fuse", os.path.join(HERE, "03_fuse_T3_present_only.py"))
t3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(t3)
h = t3.h

CLASS_NAMES = ["noConv", "conv"]
CLIN_ABBR = {"3yr": "BioClin-L-ft", "5yr": "BioClin-B-ft"}
MRI_ABBR = {"3yr": "BrainMVP-ft (plus-orig)", "5yr": "BrainMVP-frozen (none)"}
HZLAB = {"3yr": "≤3y", "5yr": "≤5y"}


def panels_for(hk):
    cfg = t3.HORIZONS[hk]
    tmpl = t3._mri_template(cfg)
    pooled = {"clin": [[], []], "mri": [[], []], "fuse": [[], []]}  # each: [y_true, pred]
    for seed in t3.SEEDS:
        clin = h.load_clinical(seed, t3.CLIN_ROOT, cfg["model"], "full_ft", cfg["clin_task"])
        mri = h.load_mri(seed, tmpl, "m12")
        f = clin.merge(mri, on="Patient_ID", how="left")
        cp = f[h.CP_COLS].to_numpy(float)
        mp = f[h.MP_COLS].to_numpy(float)
        pres = ~np.isnan(mp).any(axis=1)
        y = f["y_clin"].to_numpy(int)
        tm = (f["split"] == "test").to_numpy()
        cpt, mpt, prt, yt = cp[tm], mp[tm], pres[tm], y[tm]
        # clinical-only (full cohort)
        pooled["clin"][0] += list(yt); pooled["clin"][1] += list(cpt.argmax(1))
        # MRI-only (MRI-present subset)
        pooled["mri"][0] += list(yt[prt]); pooled["mri"][1] += list(mpt[prt].argmax(1))
        # best fusion = equal-weight avg where both present, else CL (full cohort)
        fuse = cpt.copy()
        if prt.any():
            fuse[prt] = 0.5 * cpt[prt] + 0.5 * mpt[prt]
        pooled["fuse"][0] += list(yt); pooled["fuse"][1] += list(fuse.argmax(1))
    return [
        (f"clinical-only [{CLIN_ABBR[hk]}]", "clinical @bl · full cohort", *pooled["clin"]),
        (f"MRI-only [{MRI_ABBR[hk]}]", "MRI @m12 · MRI-present subset", *pooled["mri"]),
        ("best fusion [weighted-avg (equal)]", "CL @bl ⊕ MRI @m12 · full cohort", *pooled["fuse"]),
    ]


def render(hk):
    panels = panels_for(hk)
    fig, axes = plt.subplots(1, 3, figsize=(3.6 * 3, 4.3))
    for ax, (label, tp, y, pred) in zip(axes, panels):
        y = np.asarray(y, int); pred = np.asarray(pred, int)
        cm = confusion_matrix(y, pred, labels=[0, 1])
        row = cm.sum(1, keepdims=True)
        pct = np.divide(cm, np.where(row == 0, 1, row)) * 100
        bacc = balanced_accuracy_score(y, pred) if len(np.unique(y)) > 1 else float("nan")
        n_tot = len(y); n_ps = round(n_tot / len(t3.SEEDS))
        ax.imshow(pct, cmap="Blues", vmin=0, vmax=100)
        ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
        ax.set_xticklabels(CLASS_NAMES); ax.set_yticklabels(CLASS_NAMES)
        ax.tick_params(axis="both", length=0, pad=2)          # tick labels closer to the box
        ax.set_xlabel("predicted", labelpad=3); ax.set_ylabel("true", labelpad=1)
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f"{cm[i, j]}\n{pct[i, j]:.0f}%", ha="center", va="center",
                        fontsize=9, color="white" if pct[i, j] > 55 else "black")
        ax.set_title(f"{label}\n{tp}\n{n_ps}/seed TEST (Σ{n_tot})  bACC {bacc:.3f}", fontsize=9)
    fig.suptitle(f"T3 {HZLAB[hk]} Multimodal late integration: Confusion matrix",
                 fontsize=13, fontweight="bold", y=0.99)
    fig.text(0.5, 0.925, "counts pooled across seeds, row-normalised %",
             ha="center", va="center", fontsize=9.5, fontstyle="italic")
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    out = os.path.join(HERE, f"confusion_T3_{HZLAB[hk].replace('≤', '').strip()}.png")
    fig.savefig(out, dpi=170, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"wrote {out}")
    for label, tp, y, pred in panels:
        print(f"  {label:40s} n={len(y)} bACC={balanced_accuracy_score(y, pred):.3f}")


def main():
    for hk in ("3yr", "5yr"):
        render(hk)


if __name__ == "__main__":
    main()
