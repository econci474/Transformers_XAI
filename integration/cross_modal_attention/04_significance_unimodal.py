"""
04_significance_unimodal.py
===========================
Is the best 1-year (m12) cross-modal FUSION significantly different from the UNIMODAL 1-year
references — i.e. each model's OWN stored predictions (no retrained head)?

Unimodal references (on the M12 cohort, n=451, per-seed test fold):
  - Clinical-only 1-year : the T2_m12_multiclass model's own softmax (prob_0/1/2).
  - MRI-only 1-year      : BrainDINO-T2's own softmax (probs[:, :3]) at the m12 scan.
                           (BrainMVP-T1b is binary -> cannot be a 3-way T2 reference.)
Best fusion = the joined-table winner: M12 / BrainMVP-T1b / cross_attn / ce_patient_label / lam=1.0 / wu10.

Tests (fusion vs each reference):
  - per-seed paired t-test on test bACC (n=3 seeds, df=2 — underpowered, reported with caveat)
  - pooled McNemar exact on per-patient correctness across all 3 test folds (more powerful)

Writes outputs/unimodal_refs.csv (consumed by 03c) and prints the significance verdict.
Run:  python integration/cross_modal_attention/04_significance_unimodal.py
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _xmodal_lib as L  # noqa: E402

OUT = os.path.join(HERE, "outputs")
SEEDS = [0, 1, 2]
BEST = os.path.join(OUT, "sweep", "lam1.0_wu10", "M12", "B", "cross_attn", "ce_patient_label")


def _stored_clin_m12(seed):
    df = pd.read_parquet(L.CL_M12_TMPL.format(s=seed)).set_index("Patient_ID")
    return df[["prob_0", "prob_1", "prob_2"]]


def _stored_mri_dino(seed):
    z = np.load(L.MRI_TMPL["A"].format(s=seed), allow_pickle=True)
    pr = z["probs"]
    df = pd.DataFrame({"pid": z["Patient_ID"].astype(str), "vis": z["VISCODE_2"].astype(str),
                       "p0": pr[:, 0], "p1": pr[:, 1], "p2": pr[:, 2]})
    return df[df.vis == "m12"].drop_duplicates("pid").set_index("pid")[["p0", "p1", "p2"]]


def per_seed():
    """Return per-seed test (y, fusion_pred, clin_probs, mri_probs) on the M12 451 cohort."""
    rows = []
    for s in SEEDS:
        d = L.build_cohort(s, "M12", "B")              # variant B = the best-fusion cohort (451)
        test = d["split"] == "test"
        pid_test = [p for p, m in zip(d["pids"], test) if m]
        y = d["y"][test]
        # fusion per-sample preds (predictions.parquet is the test fold)
        fp = pd.read_parquet(os.path.join(BEST, f"seed_{s}", "predictions.parquet"))
        fp = fp.set_index("Patient_ID").loc[pid_test]
        f_pred = fp["pred"].to_numpy(int)
        # unimodal stored probs on the same test pids
        cprobs = _stored_clin_m12(s).loc[pid_test].to_numpy()
        mprobs = _stored_mri_dino(s).loc[pid_test].to_numpy()
        rows.append({"seed": s, "y": y, "f_pred": f_pred, "cprobs": cprobs, "mprobs": mprobs,
                     "f_bacc": L.compute_metrics(y, fp[["prob_0", "prob_1", "prob_2"]].to_numpy())["balanced_acc"],
                     "c": L.compute_metrics(y, cprobs), "m": L.compute_metrics(y, mprobs), "n": len(y)})
    return rows


def mcnemar(correct_a, correct_b):
    """Exact McNemar on paired per-sample correctness. Returns (b, c, p)."""
    b = int(np.sum(correct_a & ~correct_b))   # A right, B wrong
    c = int(np.sum(~correct_a & correct_b))   # A wrong, B right
    p = stats.binomtest(min(b, c), b + c, 0.5).pvalue if (b + c) > 0 else 1.0
    return b, c, p


def main():
    rows = per_seed()
    f = np.array([r["f_bacc"] for r in rows])
    c = np.array([r["c"]["balanced_acc"] for r in rows])
    m = np.array([r["m"]["balanced_acc"] for r in rows])
    print("Per-seed test bACC (M12 cohort, n=451):")
    for r in rows:
        print(f"  seed {r['seed']}: fusion {r['f_bacc']:.3f} | clinical-1yr {r['c']['balanced_acc']:.3f} "
              f"| MRI-1yr {r['m']['balanced_acc']:.3f}  (n_test={r['n']})")
    print(f"\nmean bACC  fusion {f.mean():.3f}±{f.std(ddof=1):.3f} | "
          f"clinical-1yr {c.mean():.3f}±{c.std(ddof=1):.3f} | MRI-1yr {m.mean():.3f}±{m.std(ddof=1):.3f}")

    # per-seed paired t-tests (df=2)
    tc = stats.ttest_rel(f, c); tm = stats.ttest_rel(f, m)
    print(f"\nPaired t-test across seeds (n=3, df=2):")
    print(f"  fusion vs clinical-1yr : diff {f.mean()-c.mean():+.3f}  t={tc.statistic:.3f}  p={tc.pvalue:.4f}")
    print(f"  fusion vs MRI-1yr      : diff {f.mean()-m.mean():+.3f}  t={tm.statistic:.3f}  p={tm.pvalue:.4f}")

    # pooled McNemar (per-patient correctness across all 3 test folds)
    fc = np.concatenate([r["f_pred"] == r["y"] for r in rows])
    cc = np.concatenate([r["cprobs"].argmax(1) == r["y"] for r in rows])
    mc = np.concatenate([r["mprobs"].argmax(1) == r["y"] for r in rows])
    bc, cc_, pc = mcnemar(fc, cc); bm, cm, pm = mcnemar(fc, mc)
    print(f"\nPooled McNemar exact ({len(fc)} test predictions across 3 folds):")
    print(f"  fusion vs clinical-1yr : fusion-only-right={bc}, clin-only-right={cc_}, p={pc:.4f}")
    print(f"  fusion vs MRI-1yr      : fusion-only-right={bm}, MRI-only-right={cm},  p={pm:.4f}")

    # significance marker (best 1-year fusion vs each reference, McNemar exact)
    def star(p):
        return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."

    # write unimodal reference rows for 03c (mean over seeds) + significance marker
    def agg(key):
        arr = {kk: np.array([r[key][kk] for r in rows]) for kk in ("balanced_acc", "auc_roc_ovr", "macro_f1")}
        return {f"{kk}_m": arr[kk].mean() for kk in arr} | {f"{kk}_s": arr[kk].std(ddof=1) for kk in arr}
    nbar = int(np.mean([r["n"] for r in rows]))
    ref = pd.DataFrame([
        {"ref": "clinical_only_1yr", "label": "Clinical-only 1-year (unimodal)", "n": nbar,
         "sig": star(pc), "mcnemar_p": pc, "ttest_p": tc.pvalue, **agg("c")},
        {"ref": "mri_only_1yr", "label": "MRI-only 1-year (unimodal)", "n": nbar,
         "sig": star(pm), "mcnemar_p": pm, "ttest_p": tm.pvalue, **agg("m")},
    ])
    ref.to_csv(os.path.join(OUT, "unimodal_refs.csv"), index=False)
    print(f"\nWrote {os.path.join(OUT, 'unimodal_refs.csv')}  "
          f"(sig vs best fusion: clinical {star(pc)}, MRI {star(pm)})")


if __name__ == "__main__":
    main()
