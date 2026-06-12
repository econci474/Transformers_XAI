"""
04_significance_unimodal.py
===========================
Is the best 1-year (m12) cross-modal FUSION significantly different from the UNIMODAL 1-year
references — i.e. each model's OWN stored predictions (no retrained head)?

Unimodal references (own stored predictions, no retrained head):
  - Clinical-only 1-year      : T2_m12_multiclass own softmax (prob_0/1/2), M12 cohort (n=451).
  - MRI-only 1-year           : BrainDINO-T2 own softmax (probs[:, :3]) at the m12 scan, M12 cohort.
                                (BrainMVP-T1b is binary -> cannot be a 3-way T2 reference.)
  - Clinical-only Longitudinal: T2_long_multiclass own softmax at the m12 row, LONG cohort (n=478),
                                compared against the best LONG fusion (auto-discovered).

Standardised statistics (fusion vs each reference, two markers — one per metric it tests):
  - ACCURACY : pooled McNemar exact on per-patient correctness across the 3 test folds.
  - AUC      : OVO-pairwise DeLong (Hand-Till). For each class pair (CN/MCI, CN/AD, MCI/AD) restrict
               to that pair, run paired DeLong on the same patients; report the mean pairwise AUC and
               a combined marker (mean pairwise z). Caveat: pooling across the 3 seed folds reuses
               some patients (mild non-independence); the OVO pairs also share patients, so the
               combined z is approximate. Per-seed bACC mean±std and the per-seed paired t-test are
               kept as DESCRIPTIVE context only (n=3 is underpowered).

Writes outputs/unimodal_refs.csv (consumed by 03c) + outputs/delong_pairwise.csv (appendix detail).
Run:  python integration/cross_modal_attention/04_significance_unimodal.py
"""
import glob
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
# BrainMVP T2 (3-way) own softmax. NB: fusion variant B uses BrainMVP-T1b (binary) embeddings, which
# have no 3-way own-softmax; the multiclass T2 head is the only valid 3-way BrainMVP reference.
MVP_T2_TMPL = os.path.join(L.MRI_ROOT, "brainmvp_embeddings", "T2_multiclass", "aug_stochastic",
                           "seed_{s}", "embeddings_seed_{s}.npz")


def _stored_clin_m12(seed):
    df = pd.read_parquet(L.CL_M12_TMPL.format(s=seed)).set_index("Patient_ID")
    return df[["prob_0", "prob_1", "prob_2"]]


def _stored_mri_3way(npz_path):
    z = np.load(npz_path, allow_pickle=True)
    pr = z["probs"]
    df = pd.DataFrame({"pid": z["Patient_ID"].astype(str), "vis": z["VISCODE_2"].astype(str),
                       "p0": pr[:, 0], "p1": pr[:, 1], "p2": pr[:, 2]})
    return df[df.vis == "m12"].drop_duplicates("pid").set_index("pid")[["p0", "p1", "p2"]]


def _stored_mri_dino(seed):
    return _stored_mri_3way(L.MRI_TMPL["A"].format(s=seed))


def _stored_mri_mvp(seed):
    return _stored_mri_3way(MVP_T2_TMPL.format(s=seed))


def _stored_clin_long(seed):
    """The longitudinal model (T2_long_multiclass): own softmax at the m12 row, per patient."""
    df = pd.read_parquet(L.CL_LONG_TMPL.format(s=seed))
    df = df[df["VISCODE_2"] == "m12"].drop_duplicates("Patient_ID").set_index("Patient_ID")
    return df[["prob_0", "prob_1", "prob_2"]]


def _best_long_dir():
    """Auto-discover the best LONG fusion: glob every LONG run, group by config dir (parent of
    seed_*), pick the max mean test bACC. Returns (dir, variant, arch, clin_pool, loss, lam, warmup)."""
    metas = glob.glob(os.path.join(OUT, "**", "LONG", "**", "metrics.json"), recursive=True)
    by_dir = {}
    for mf in metas:
        cfgdir = os.path.dirname(os.path.dirname(mf))          # strip seed_*/metrics.json
        by_dir.setdefault(cfgdir, []).append(
            json.load(open(mf))["test_metrics"]["balanced_acc"])
    best = max(by_dir, key=lambda d: float(np.mean(by_dir[d])))
    parts = best.replace("\\", "/").split("/")
    i = parts.index("LONG")
    variant, arch = parts[i + 1], parts[i + 2]
    pool = parts[i + 3] if arch == "cross_attn_x" else "mean"
    loss = parts[i + 4] if arch == "cross_attn_x" else parts[i + 3]
    lam, warmup = 0.0, 0
    for p in parts:
        if p.startswith("lam") and "_wu" in p:
            lam = float(p[3:p.index("_wu")]); warmup = int(p[p.index("_wu") + 3:])
    return best, variant, arch, pool, loss, lam, warmup


# ── DeLong's test (Sun & Xu 2014 fast midrank algorithm) ─────────────────────────
def _midrank(x):
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1.0       # 1-based midrank, ties averaged
        i = j
    out = np.empty(N, dtype=float)
    out[J] = T
    return out


def _fast_delong(preds_sorted, m):
    """preds_sorted: (k, N) predictor scores, columns ordered positives (m) then negatives (n).
    Returns (aucs[k], covariance[k,k]) of the k correlated AUCs (DeLong)."""
    k, N = preds_sorted.shape
    n = N - m
    pos, neg = preds_sorted[:, :m], preds_sorted[:, m:]
    tx = np.vstack([_midrank(pos[r]) for r in range(k)])
    ty = np.vstack([_midrank(neg[r]) for r in range(k)])
    tz = np.vstack([_midrank(preds_sorted[r]) for r in range(k)])
    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    cov = np.cov(v01) / m + np.cov(v10) / n
    return aucs, np.atleast_2d(cov)


def delong_binary(y_bin, sa, sb):
    """Paired DeLong on two correlated AUCs over the same binary labels. (sa=A, sb=B.)
    Returns (auc_a, auc_b, z, p_two_sided); z>0 means A separates better than B."""
    y_bin = np.asarray(y_bin, int)
    order = np.argsort(-y_bin, kind="stable")              # positives (1) first
    m = int(y_bin.sum())
    preds = np.vstack([np.asarray(sa, float), np.asarray(sb, float)])[:, order]
    aucs, cov = _fast_delong(preds, m)
    var = cov[0, 0] + cov[1, 1] - 2.0 * cov[0, 1]
    if var <= 0:
        z = 0.0
        p = 1.0
    else:
        z = (aucs[0] - aucs[1]) / np.sqrt(var)
        p = float(2.0 * stats.norm.sf(abs(z)))
    return float(aucs[0]), float(aucs[1]), float(z), float(p)


def ovo_delong(y, Pa, Pb, n_classes=3):
    """OVO (Hand-Till) DeLong: for each pair (i,j) restrict to y∈{i,j}, binary target (y==j),
    discriminant p_j-p_i; paired DeLong A vs B. A = first arg (fusion), B = second (reference).
    Returns dict: pairs[], macro_auc_a, macro_auc_b, z_combined, p_combined."""
    y = np.asarray(y, int)
    pairs, zs = [], []
    for i in range(n_classes):
        for j in range(i + 1, n_classes):
            mask = np.isin(y, [i, j])
            yb = (y[mask] == j).astype(int)
            if yb.sum() in (0, len(yb)):                  # one class absent in this fold -> skip
                continue
            sa = Pa[mask, j] - Pa[mask, i]
            sb = Pb[mask, j] - Pb[mask, i]
            auc_a, auc_b, z, p = delong_binary(yb, sa, sb)
            pairs.append({"pair": f"{i}v{j}", "n": int(mask.sum()),
                          "auc_a": auc_a, "auc_b": auc_b, "z": z, "p": p})
            zs.append(z)
    macro_a = float(np.mean([pp["auc_a"] for pp in pairs])) if pairs else float("nan")
    macro_b = float(np.mean([pp["auc_b"] for pp in pairs])) if pairs else float("nan")
    zc = float(np.mean(zs)) if zs else 0.0
    pc = float(2.0 * stats.norm.sf(abs(zc))) if zs else 1.0
    return {"pairs": pairs, "macro_auc_a": macro_a, "macro_auc_b": macro_b,
            "z_combined": zc, "p_combined": pc}


def per_seed():
    """Return per-seed test (y, fusion_pred/probs, clin_probs, mri_probs) on the M12 451 cohort."""
    rows = []
    for s in SEEDS:
        d = L.build_cohort(s, "M12", "B")              # variant B = the best-fusion cohort (451)
        test = d["split"] == "test"
        pid_test = [p for p, m in zip(d["pids"], test) if m]
        y = d["y"][test]
        # fusion per-sample preds + probs (predictions.parquet is the test fold)
        fp = pd.read_parquet(os.path.join(BEST, f"seed_{s}", "predictions.parquet"))
        fp = fp.set_index("Patient_ID").loc[pid_test]
        f_pred = fp["pred"].to_numpy(int)
        fprobs = fp[["prob_0", "prob_1", "prob_2"]].to_numpy()
        # unimodal stored probs on the same test pids
        cprobs = _stored_clin_m12(s).loc[pid_test].to_numpy()
        mprobs = _stored_mri_dino(s).loc[pid_test].to_numpy()
        vprobs = _stored_mri_mvp(s).loc[pid_test].to_numpy()
        rows.append({"seed": s, "y": y, "f_pred": f_pred, "fprobs": fprobs,
                     "cprobs": cprobs, "mprobs": mprobs, "vprobs": vprobs,
                     "f_bacc": L.compute_metrics(y, fprobs)["balanced_acc"],
                     "c": L.compute_metrics(y, cprobs), "m": L.compute_metrics(y, mprobs),
                     "v": L.compute_metrics(y, vprobs), "n": len(y)})
    return rows


def long_per_seed(best_dir, variant):
    """Per-seed test (y, fusion_pred/probs, clin_long_probs) on the LONG cohort for the best LONG fusion."""
    rows = []
    for s in SEEDS:
        d = L.build_cohort(s, "LONG", variant)
        test = d["split"] == "test"
        pid_test = [p for p, m in zip(d["pids"], test) if m]
        y = d["y"][test]
        fp = pd.read_parquet(os.path.join(best_dir, f"seed_{s}", "predictions.parquet"))
        fp = fp.set_index("Patient_ID").loc[pid_test]
        cprobs = _stored_clin_long(s).loc[pid_test].to_numpy()
        rows.append({"seed": s, "y": y, "f_pred": fp["pred"].to_numpy(int),
                     "fprobs": fp[["prob_0", "prob_1", "prob_2"]].to_numpy(),
                     "cprobs": cprobs, "c": L.compute_metrics(y, cprobs), "n": len(y)})
    return rows


def mcnemar(correct_a, correct_b):
    """Exact McNemar on paired per-sample correctness. Returns (b, c, p)."""
    b = int(np.sum(correct_a & ~correct_b))   # A right, B wrong
    c = int(np.sum(~correct_a & correct_b))   # A wrong, B right
    p = stats.binomtest(min(b, c), b + c, 0.5).pvalue if (b + c) > 0 else 1.0
    return b, c, p


def star(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."


def compare(rows, ref_key, fusion_label):
    """Pool fusion-vs-(ref_key) across seeds → McNemar (accuracy) + OVO-DeLong (AUC). Returns a dict
    of markers/p-values + the per-pair DeLong detail rows (tagged with fusion_label)."""
    y = np.concatenate([r["y"] for r in rows])
    fpred = np.concatenate([r["f_pred"] for r in rows])
    fprobs = np.concatenate([r["fprobs"] for r in rows])
    rprobs = np.concatenate([r[ref_key] for r in rows])
    fc = fpred == y
    rc = rprobs.argmax(1) == y
    _, _, p_acc = mcnemar(fc, rc)
    dl = ovo_delong(y, fprobs, rprobs)
    detail = [{**pp, "comparison": fusion_label} for pp in dl["pairs"]]
    return {"sig_acc": star(p_acc), "mcnemar_p": p_acc,
            "sig_auc": star(dl["p_combined"]), "delong_p": dl["p_combined"],
            "auc_ovo_fusion": dl["macro_auc_a"], "auc_ovo_ref": dl["macro_auc_b"]}, detail


def _agg(rows, key, nseeds_ddof=1):
    """Mean ± std (over seeds) of the displayed OVR metrics for reference `key` ('c' or 'm')."""
    arr = {kk: np.array([r[key][kk] for r in rows])
           for kk in ("balanced_acc", "auc_roc_ovr", "macro_f1")}
    return ({f"{kk}_m": float(arr[kk].mean()) for kk in arr}
            | {f"{kk}_s": float(arr[kk].std(ddof=nseeds_ddof)) for kk in arr})


def main():
    # ── 1-year (M12) references ──────────────────────────────────────────────────
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

    # descriptive only: per-seed paired t-tests (n=3, df=2 — underpowered)
    tc = stats.ttest_rel(f, c); tm = stats.ttest_rel(f, m)
    print(f"\n[descriptive] paired t-test across seeds (n=3, df=2):")
    print(f"  fusion vs clinical-1yr : diff {f.mean()-c.mean():+.3f}  t={tc.statistic:.3f}  p={tc.pvalue:.4f}")
    print(f"  fusion vs MRI-1yr      : diff {f.mean()-m.mean():+.3f}  t={tm.statistic:.3f}  p={tm.pvalue:.4f}")

    # primary tests: McNemar (accuracy) + OVO-DeLong (AUC), pooled across folds
    mk_c, det_c = compare(rows, "cprobs", "M12 fusion vs clinical-1yr")
    mk_m, det_m = compare(rows, "mprobs", "M12 fusion vs MRI-DINO-1yr")
    mk_v, det_v = compare(rows, "vprobs", "M12 fusion vs MRI-MVP-1yr")
    print(f"\nM12 fusion vs references (pooled over 3 folds):")
    print(f"  vs clinical-1yr     : acc {mk_c['sig_acc']} (McNemar p={mk_c['mcnemar_p']:.4f}) | "
          f"AUC {mk_c['sig_auc']} (OVO-DeLong p={mk_c['delong_p']:.4f})")
    print(f"  vs MRI BrainDINO-1yr: acc {mk_m['sig_acc']} (McNemar p={mk_m['mcnemar_p']:.4f}) | "
          f"AUC {mk_m['sig_auc']} (OVO-DeLong p={mk_m['delong_p']:.4f})")
    print(f"  vs MRI BrainMVP-1yr : acc {mk_v['sig_acc']} (McNemar p={mk_v['mcnemar_p']:.4f}) | "
          f"AUC {mk_v['sig_auc']} (OVO-DeLong p={mk_v['delong_p']:.4f})  "
          f"(BrainMVP-T2 multiclass own softmax)")

    # ── Longitudinal clinical reference (auto-discover best LONG fusion) ──────────
    best_dir, lvar, larch, lpool, lloss, llam, lwu = _best_long_dir()
    print(f"\nBest LONG fusion: {os.path.relpath(best_dir, OUT)}  (variant {lvar}, {larch}/{lpool}, "
          f"{lloss}, λ={llam:g} wu{lwu})")
    lrows = long_per_seed(best_dir, lvar)
    lf = np.array([L.compute_metrics(r["y"], r["fprobs"])["balanced_acc"] for r in lrows])
    lc = np.array([r["c"]["balanced_acc"] for r in lrows])
    print(f"  per-seed bACC  LONG-fusion {lf.mean():.3f}±{lf.std(ddof=1):.3f} | "
          f"clinical-long {lc.mean():.3f}±{lc.std(ddof=1):.3f}")
    mk_l, det_l = compare(lrows, "cprobs", "LONG fusion vs clinical-long")
    print(f"  vs clinical-long: acc {mk_l['sig_acc']} (McNemar p={mk_l['mcnemar_p']:.4f}) | "
          f"AUC {mk_l['sig_auc']} (OVO-DeLong p={mk_l['delong_p']:.4f})")

    # ── write reference rows for 03c (mean over seeds) + both markers ────────────
    nbar = int(np.mean([r["n"] for r in rows]))
    nbar_long = int(np.mean([r["n"] for r in lrows]))
    ref = pd.DataFrame([
        {"ref": "clinical_only_1yr", "label": "Clinical-only 1-year (unimodal)", "n": nbar,
         "sig": mk_c["sig_acc"], "sig_acc": mk_c["sig_acc"], "mcnemar_p": mk_c["mcnemar_p"],
         "sig_auc": mk_c["sig_auc"], "delong_p": mk_c["delong_p"],
         "auc_ovo_m": mk_c["auc_ovo_ref"], "ttest_p": float(tc.pvalue), **_agg(rows, "c")},
        {"ref": "mri_only_1yr", "label": "MRI-only 1-year — BrainDINO (unimodal)", "n": nbar,
         "sig": mk_m["sig_acc"], "sig_acc": mk_m["sig_acc"], "mcnemar_p": mk_m["mcnemar_p"],
         "sig_auc": mk_m["sig_auc"], "delong_p": mk_m["delong_p"],
         "auc_ovo_m": mk_m["auc_ovo_ref"], "ttest_p": float(tm.pvalue), **_agg(rows, "m")},
        {"ref": "mri_mvp_1yr", "label": "MRI-only 1-year — BrainMVP (unimodal)", "n": nbar,
         "sig": mk_v["sig_acc"], "sig_acc": mk_v["sig_acc"], "mcnemar_p": mk_v["mcnemar_p"],
         "sig_auc": mk_v["sig_auc"], "delong_p": mk_v["delong_p"],
         "auc_ovo_m": mk_v["auc_ovo_ref"], "ttest_p": float("nan"), **_agg(rows, "v")},
        {"ref": "clinical_only_long", "label": "Clinical-only Longitudinal (unimodal)", "n": nbar_long,
         "sig": mk_l["sig_acc"], "sig_acc": mk_l["sig_acc"], "mcnemar_p": mk_l["mcnemar_p"],
         "sig_auc": mk_l["sig_auc"], "delong_p": mk_l["delong_p"],
         "auc_ovo_m": mk_l["auc_ovo_ref"], "ttest_p": float("nan"), **_agg(lrows, "c")},
    ])
    ref.to_csv(os.path.join(OUT, "unimodal_refs.csv"), index=False)
    det = det_c + det_m + det_v + det_l
    pd.DataFrame(det).to_csv(os.path.join(OUT, "delong_pairwise.csv"), index=False)
    print(f"\nWrote {os.path.join(OUT, 'unimodal_refs.csv')} (4 refs) "
          f"+ delong_pairwise.csv ({len(det)} pair-rows)")


if __name__ == "__main__":
    main()
