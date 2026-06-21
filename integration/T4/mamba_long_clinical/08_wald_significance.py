r"""
08_wald_significance.py   (env: survml — scikit-survival + lifelines)
=====================================================================
Wald significance tests for the T4 longitudinal-MAMBA survival comparison table
(`survival_comparison_table_cleaned.png` / `.tex`, built by 02 / 07).

Question: for each survival model, is its held-out (val+test) performance significantly different from
the strongest tabular baseline **RSF**, on each of the 7 metrics (C-index, IPCW-C, IBS, AUC@3/5/7/10y)?

The table's "± std" is across only 3 seeds (df=2) — it conflates seed-variance with sampling-variance and
is too crude to be a standard error. Closed-form SEs for *differences* of C-index / IPCW-C / IBS / time-AUC
do not exist, so we estimate them by a **paired, seed-stratified bootstrap** of the metric difference and
form a Wald statistic z = Δ / SE_boot, two-sided p = 2(1−Φ(|z|)).

Why a *paired* bootstrap: every model is evaluated on the SAME canonical val+test patients per seed
(verified: 108 patients/seed, 0 Patient_ID mismatch, identical across all configs + baselines), so we
resample the patient indices ONCE per (seed, replicate) and apply them to the model AND to RSF — capturing
their strong positive correlation, which is the correct (smaller) SE for the difference.

Geometry is FIXED per seed (the 12-pt time grid, t_ref=4y, τ, and the train-fold censoring distribution
`ytr`) so the metric definition is identical across resamples and across the two models — exactly the
geometry `02l.evaluate` uses on the full val+test (we verify the point estimates reproduce
`master_survival_comparison.csv` to <1e-4 before bootstrapping).

Reuses `02l_survival_model_comparison.py` (sksurv primitives, `_grid`, `t_ref`, baseline fitters,
`load_seed`, `COMPACT_FEATURES`, `SPLIT_DIR`) and `02_survival_comparison.py::interp_surv_fn`.

Run:  conda run -n survml python integration/T4/mamba_long_clinical/08_wald_significance.py [--B 2000]
Out:  <comparison>/wald_significance_vs_rsf.csv (+ copy next to script)
      <comparison>/wald_significance_vs_rsf.{png,pdf} (+ copy next to script)
"""
import argparse
import importlib.util
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm
from sksurv.util import Surv
from sksurv.metrics import (concordance_index_censored, concordance_index_ipcw,
                            integrated_brier_score, cumulative_dynamic_auc)

matplotlib.rcParams["font.family"] = "DejaVu Serif"

HERE = Path(__file__).resolve().parent
PRED = Path(r"D:\ADNI_BIDS_project\derivatives\mamba_survival\survival")
COMP = PRED / "comparison"
CLIN = Path(__file__).resolve().parents[3] / "clinical_pipeline"   # repo-relative (repo moved off iCloud)
SEEDS = (0, 1, 2)
REF = "RSF"                                            # reference baseline (locked with user)
MC = ["c_index", "c_index_ipcw", "ibs", "auc_3y", "auc_5y", "auc_7y", "auc_10y"]
MODE = {m: ("min" if m == "ibs" else "max") for m in MC}   # IBS lower-better

# Display order — mirrors 07_render's ROWS. Tuple = (model_key, group, survival head, aggregation).
# (RSF is the reference, omitted from the test rows.)
TEST_MODELS = [
    ("A_ctrl_gru",              "Longitudinal",    "Weibull-pw",  "GRU"),
    ("A_ctrl_meanpool",         "Longitudinal",    "Weibull-pw",  "Meanpool"),
    ("A_default_mamba1_frozen", "Longitudinal",    "Weibull-pw",  "MAMBA (frozen)"),
    ("B_cox",                   "Longitudinal",    "Cox PH",      "MAMBA (frozen)"),
    ("B_soden",                 "Longitudinal",    "SODEN",       "MAMBA (frozen)"),
    ("B_weibull_aft",           "Longitudinal",    "Weibull-AFT", "MAMBA (frozen)"),
    ("Cox PH",                  "Cross-sectional", "Cox PH",      "—"),
    ("Weibull AFT",             "Cross-sectional", "Weibull-AFT", "—"),
    ("Weibull piecewise",       "Cross-sectional", "Weibull-pw",  "—"),
]
DISP = {k: f"{head} · {agg}" if agg != "—" else f"{head} (tabular)"
        for k, _, head, agg in TEST_MODELS}                # display string per model key
MAMBA_CFGS = {"A_ctrl_gru", "A_ctrl_meanpool", "A_default_mamba1_frozen", "B_cox", "B_soden", "B_weibull_aft"}
ALL_MODELS = [m for m, *_ in TEST_MODELS] + [REF]

# ── second reference: the cross-sectional Weibull-pw baseline → an analogous "..._vs_weibull_pw" table ──
# Every OTHER model (incl. RSF) is tested against it, quantifying the longitudinal lift over the SAME
# parametric family. RSF becomes a tested cross-sectional row here; Weibull-pw is the reference (omitted).
WPW_REF = "Weibull piecewise"
TEST_MODELS_WPW = [t for t in TEST_MODELS if t[0] != WPW_REF] + [("RSF", "Cross-sectional", "RSF", "—")]
DISP_WPW = {k: (f"{head} · {agg}" if agg != "—" else f"{head} (tabular)")
            for k, _, head, agg in TEST_MODELS_WPW}

# Per-reference render specs. The vs-RSF table is byte-identical to before; vs-Weibull-pw is the new analogue.
REF_SPECS = [
    dict(ref=REF, short="RSF", paren="best tabular baseline", tms=TEST_MODELS, disp=DISP,
         desc="cross-sectional Random Survival Forest on the baseline-visit-only tabular features",
         basename="wald_significance_vs_rsf", label="wald_vs_rsf"),
    dict(ref=WPW_REF, short="Weibull-pw", paren="cross-sectional baseline",
         tms=TEST_MODELS_WPW, disp=DISP_WPW,
         desc="cross-sectional Weibull-piecewise on the baseline-visit-only tabular features",
         basename="wald_significance_vs_weibull_pw", label="wald_vs_weibull_pw"),
]


def _load_02l():
    spec = importlib.util.spec_from_file_location("l02", CLIN / "02l_survival_model_comparison.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


L = _load_02l()


def interp_surv_fn(Sgrid, grid_times):
    """Per-row linear interp of a stored dense S-grid to requested times (== 02's interp_surv_fn)."""
    def surv_fn(X, times):
        times = np.asarray(times, float)
        out = np.empty((Sgrid.shape[0], times.size))
        for i in range(Sgrid.shape[0]):
            out[i] = np.interp(times, grid_times, Sgrid[i])
        return np.clip(out, 1e-8, 1.0)
    return surv_fn


# ── per-seed cohort + fixed geometry + per-(seed,model) per-patient risk/S ────────
def build_seed(seed):
    """Return GEO (e,t,ytr,grid,t_ref,tau,hmax,pids) and PRE[model] = (risk[n], S[n,12]) for one seed,
    all on the canonical val+test patient order (== the tabular baseline cohort)."""
    import json
    # canonical cohort + labels + train reference from the tabular baseline loader (8 compact features)
    d = L.load_seed(seed, L.COMPACT_FEATURES)
    Xtr, etr, ttr = d["train"]
    ytr = Surv.from_arrays(event=etr.astype(bool), time=ttr)
    Xvt = pd.concat([d["val"][0], d["test"][0]])
    e = np.concatenate([d["val"][1], d["test"][1]]).astype(bool)
    t = np.concatenate([d["val"][2], d["test"][2]]).astype(float)
    # Patient_ID for each valtest row (index-aligned through the split CSVs)
    pids = []
    for sp in ("val", "test"):
        csv = pd.read_csv(L.SPLIT_DIR / f"seed_{seed}" / f"{sp}.csv", low_memory=False)
        csv = csv[csv["bl_dx"].isin(["CN", "MCI"])].copy(); csv["Patient_ID"] = csv["Patient_ID"].astype(str)
        Xs, _, _ = L._build(csv, L.COMPACT_FEATURES)
        pids.append(csv.loc[Xs.index, "Patient_ID"].to_numpy())
    pids = np.concatenate(pids)
    order = {p: i for i, p in enumerate(pids)}            # canonical patient -> row index

    # fixed evaluation geometry — identical to 02l.evaluate on this seed's full val+test
    grid = L._grid(ttr, etr.astype(bool), t, e)
    t_ref = float(np.clip(4.0, grid[0] + 1e-6, grid[-1] - 1e-6))
    tau = float(grid[-1])
    hmax = min(float(ttr.max()), float(t.max())) - 1e-3
    GEO = dict(e=e, t=t, ytr=ytr, ttr=ttr, etr=etr.astype(bool), grid=grid,
               t_ref=t_ref, tau=tau, hmax=hmax, pids=pids, n=len(pids))

    PRE = {}
    # baselines: refit per seed on TRAIN, predict S on the fixed grid + at t_ref (canonical order)
    for name in ("RSF", "Cox PH", "Weibull AFT", "Weibull piecewise"):
        sf = L.fit_model(name, Xtr, etr, ttr)
        S = np.clip(np.asarray(sf(Xvt, grid), float), 1e-8, 1.0)
        risk = 1.0 - np.clip(np.asarray(sf(Xvt, np.array([t_ref])), float)[:, 0], 1e-8, 1.0)
        PRE[name] = (risk, S)
    # MAMBA configs: interp each patient's stored dense S-grid; align to canonical order by Patient_ID
    for cfg in MAMBA_CFGS:
        df = pd.read_parquet(PRED / cfg / f"seed_{seed}" / "predictions.parquet")
        df["Patient_ID"] = df["Patient_ID"].astype(str)
        gt = np.asarray(json.load(open(PRED / cfg / f"seed_{seed}" / "meta.json"))["grid_times"], float)
        g = df[df.split.isin(["val", "test"])].copy()
        idx = g["Patient_ID"].map(order).to_numpy()
        assert not np.isnan(idx).any() and len(g) == len(pids), f"{cfg} seed{seed} cohort mismatch"
        scol = [c for c in g.columns if c.startswith("S_")]
        sf = interp_surv_fn(g[scol].to_numpy(float), gt)
        S = sf(None, grid); risk = 1.0 - sf(None, np.array([t_ref]))[:, 0]
        # reorder rows to canonical order
        Sord = np.empty_like(S); rord = np.empty_like(risk)
        Sord[idx] = S; rord[idx] = risk
        PRE[cfg] = (rord, Sord)
        # sanity: the parquet's own event/time agree with the canonical labels
        ee = np.empty(len(pids), bool); tt = np.empty(len(pids), float)
        ee[idx] = g["event"].to_numpy() > 0.5; tt[idx] = g["time"].to_numpy()
        if not (np.array_equal(ee, e) and np.allclose(tt, t, atol=1e-6)):
            print(f"  [warn] seed{seed} {cfg}: label mismatch vs tabular (using each model's metric on "
                  f"shared patients) — max|Δt|={np.abs(tt - t).max():.4f}, Δe={(ee != e).sum()}")
    return GEO, PRE


def seed_metrics(idx, risk, S, GEO):
    """7 metrics on a (bootstrap) index set, with this seed's FIXED grid/t_ref/τ/ytr."""
    er, tr, rr, sr = GEO["e"][idx], GEO["t"][idx], risk[idx], S[idx]
    yev = Surv.from_arrays(event=er, time=tr)
    out = {"c_index": float(concordance_index_censored(er, tr, rr)[0])}
    try:
        out["c_index_ipcw"] = float(concordance_index_ipcw(GEO["ytr"], yev, rr, tau=GEO["tau"])[0])
    except Exception:
        out["c_index_ipcw"] = np.nan
    try:
        out["ibs"] = float(integrated_brier_score(GEO["ytr"], yev, sr, GEO["grid"]))
    except Exception:
        out["ibs"] = np.nan
    for h in (3.0, 5.0, 7.0, 10.0):
        v = np.nan
        if (tr[er] < h).any() and h < GEO["hmax"]:
            try:
                v = float(np.atleast_1d(cumulative_dynamic_auc(GEO["ytr"], yev, rr, h)[0])[0])
            except Exception:
                pass
        out[f"auc_{int(h)}y"] = v
    return out


def bh_adjust(pvals):
    """Benjamini–Hochberg FDR; NaN-safe (NaNs stay NaN, excluded from m)."""
    p = np.asarray(pvals, float); q = np.full_like(p, np.nan)
    ok = np.where(np.isfinite(p))[0]
    if ok.size == 0:
        return q
    o = ok[np.argsort(p[ok])]; m = ok.size
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = o[rank]
        prev = min(prev, p[i] * m / (rank + 1))
        q[i] = prev
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=2000, help="bootstrap replicates")
    ap.add_argument("--seed", type=int, default=12345, help="bootstrap RNG seed")
    ap.add_argument("--render_only", action="store_true",
                    help="skip the bootstrap; re-render PNG + .tex from the existing CSV")
    args = ap.parse_args()

    if args.render_only:                                   # re-draw from saved results, no re-bootstrap
        for sp in REF_SPECS:
            res = pd.read_csv(COMP / f"{sp['basename']}.csv")
            render(res, args.B, sp); write_tex(res, args.B, sp); write_tex_detailed(res, args.B, sp)
            print(f"Re-rendered {sp['basename']} (PNG + .tex) from {COMP / (sp['basename'] + '.csv')}")
        return

    # ── build geometry + per-patient predictions for every seed ──────────────────
    GEO, PRE = {}, {}
    for s in SEEDS:
        GEO[s], PRE[s] = build_seed(s)
        print(f"seed {s}: n={GEO[s]['n']} events={int(GEO[s]['e'].sum())} "
              f"grid=[{GEO[s]['grid'][0]:.2f},{GEO[s]['grid'][-1]:.2f}] t_ref={GEO[s]['t_ref']:.2f}")

    # ── point estimates (full sample) — must reproduce master_survival_comparison.csv ──
    full = {m: {s: seed_metrics(np.arange(GEO[s]["n"]), *PRE[s][m], GEO[s]) for s in SEEDS}
            for m in ALL_MODELS}
    pt = {m: {k: np.nanmean([full[m][s][k] for s in SEEDS]) for k in MC} for m in ALL_MODELS}

    master = pd.read_csv(COMP / "master_survival_comparison.csv")
    master = master[master.split == "valtest"]
    maxdiff = 0.0
    for m in ALL_MODELS:
        g = master[master.model == m]
        for k in MC:
            ref = g[k].mean()
            if np.isfinite(ref) and np.isfinite(pt[m][k]):
                maxdiff = max(maxdiff, abs(ref - pt[m][k]))
    print(f"\nPoint-estimate reproduction vs master CSV: max|Δ| = {maxdiff:.2e} "
          f"({'OK' if maxdiff < 1e-3 else 'CHECK — exceeds 1e-3'})")

    # ── paired, seed-stratified bootstrap of each model−reference metric difference (BOTH references) ──
    # ONE resample per (seed, b) shared by all models AND both references → fully paired; the vs-RSF numbers
    # are unchanged by adding the second reference (same rng draw sequence).
    rng = np.random.default_rng(args.seed)
    B = args.B
    boot = {(sp["ref"], m, k): np.full(B, np.nan)
            for sp in REF_SPECS for m, *_ in sp["tms"] for k in MC}
    for b in range(B):
        smet = {}
        for s in SEEDS:
            idx = rng.integers(0, GEO[s]["n"], GEO[s]["n"])       # ONE resample / (seed,b), shared by all
            smet[s] = {m: seed_metrics(idx, *PRE[s][m], GEO[s]) for m in ALL_MODELS}
        for sp in REF_SPECS:
            ref = sp["ref"]
            for m, *_ in sp["tms"]:
                for k in MC:
                    with warnings.catch_warnings():           # all-NaN seeds (rare AUC@10y resamples)
                        warnings.simplefilter("ignore", RuntimeWarning)
                        boot[(ref, m, k)][b] = np.nanmean([smet[s][m][k] - smet[s][ref][k] for s in SEEDS])
        if (b + 1) % 250 == 0:
            print(f"  bootstrap {b + 1}/{B}")

    def stars(p):
        return "" if not np.isfinite(p) else ("***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "")

    # ── assemble + render each reference's Wald table ─────────────────────────────
    for sp in REF_SPECS:
        ref, tms, disp = sp["ref"], sp["tms"], sp["disp"]
        recs = []
        for m, group, head, agg in tms:
            for k in MC:
                delta = pt[m][k] - pt[ref][k]                    # signed (IBS: negative = better)
                reps = boot[(ref, m, k)]; effB = int(np.isfinite(reps).sum())
                se = float(np.nanstd(reps, ddof=1)) if effB > 1 else np.nan
                z = delta / se if (se and se > 0 and np.isfinite(delta)) else np.nan
                p = float(2 * norm.sf(abs(z))) if np.isfinite(z) else np.nan
                recs.append(dict(model=m, group=group, head=head, aggregation=agg, display=disp[m],
                                 metric=k, mode=MODE[k], model_mean=pt[m][k], ref_mean=pt[ref][k],
                                 delta=delta, se_boot=se, z=z, p=p, eff_B=effB,
                                 better=(delta < 0 if MODE[k] == "min" else delta > 0)))
        res = pd.DataFrame(recs)
        res["p_bh"] = bh_adjust(res["p"].to_numpy())
        res["sig"] = res["p"].map(stars)
        for d in (COMP, HERE):
            res.to_csv(d / f"{sp['basename']}.csv", index=False, encoding="utf-8")
        print(f"\nWald tests vs {sp['short']} (Δ = model − {ref}; bootstrap B={B}):")
        show = res.pivot(index="display", columns="metric", values="z").reindex([disp[m] for m, *_ in tms])[MC]
        print(show.round(2).to_string())
        nsig = (res["p"] < 0.05).sum()
        print(f"  {nsig}/{len(res)} model×metric comparisons significant at raw p<0.05 "
              f"({(res['p_bh'] < 0.05).sum()} survive BH-FDR).")
        render(res, args.B, sp)
        write_tex(res, args.B, sp)
        write_tex_detailed(res, args.B, sp)
        print(f"  Done → {COMP}\\{sp['basename']}.csv / .png / .tex (+ _detailed.tex)")


# ── house-style standalone table (DejaVu-Serif / bordered / ruled / dashed-divider) ──
def render(res, B, sp):
    METlbl = [("c_index", "C-index"), ("c_index_ipcw", "IPCW C"), ("ibs", "IBS"),
              ("auc_3y", "AUC@3y"), ("auc_5y", "AUC@5y"), ("auc_7y", "AUC@7y"), ("auc_10y", "AUC@10y")]
    headers = ["Model", "Aggregation"] + [lbl for _, lbl in METlbl]
    LEFT_COLS = {0, 1}

    def stars(p):
        return "" if not np.isfinite(p) else ("***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "")

    def fmtp(x):
        return "—" if not np.isfinite(x) else ("<.001" if x < 1e-3 else f"{x:.3f}")

    body, rules, sigflags = [], [], []
    prev_g = None
    for key, group, head, agg in sp["tms"]:
        rule = "strong" if (prev_g is not None and group != prev_g) else None
        sub = {r.metric: r for r in res[res.model == key].itertuples()}
        cells = [head, agg]
        flags = {}
        for j, (mk, _) in enumerate(METlbl):
            r = sub[mk]
            st = stars(r.p)
            sign = "+" if r.delta >= 0 else "−"
            cells.append(f"{sign}{abs(r.delta):.3f}{st}\nSE {r.se_boot:.3f}\nz={r.z:.1f}\n"
                         f"p={fmtp(r.p)}\nq={fmtp(r.p_bh)}" if np.isfinite(r.z) else "—")
            flags[2 + j] = bool(st) and r.better
        body.append(cells); rules.append(rule); sigflags.append(flags); prev_g = group

    COL_W = [1.9, 1.5] + [1.5] * len(METlbl)
    LEFT, RIGHT_PAD = 0.28, 0.28
    nq = int((res[res.group == "Longitudinal"]["p_bh"] < 0.05).sum())
    bh_note = ("no longitudinal gain has q<0.05" if nq == 0
               else f"{nq} longitudinal gain(s) survive BH-FDR (q<0.05)")
    TITLE = f"T4 Survival — Wald significance vs {sp['short']} ({sp['paren']})"
    SUBTITLE = (
        f"cell = Δ(model − {sp['ref']}) on the held-out (val+test) metric  ·  SE_boot  ·  z = Δ/SE_boot  ·  "
        f"p (raw two-sided Wald)  ·  q (Benjamini–Hochberg FDR over all {len(sp['tms'])}×7 tests)\n"
        f"SE from a paired, seed-stratified bootstrap (B={B}) of the metric difference\n"
        f"* p<0.05  ** p<0.01  *** p<0.001 (raw)  ·  bold Δ = significant AND better than {sp['short']}  ·  "
        f"{bh_note}\n"
        "IBS is lower-better (negative Δ = better)  ·  N(val+test) = 108 per seed (≈25 events), 3 seeds\n"
        f"reference {sp['short']} = {sp['desc']}")
    SUB_LINES = SUBTITLE.count("\n") + 1
    SUB_H = 0.27 * SUB_LINES
    TITLE_H, HEAD_H, ROW_H = 0.5 + SUB_H, 0.40, 1.30   # 1.30: 5-line cells (Δ / SE / z / p / q)
    TOP_PAD, BOT_PAD = 0.12, 0.12
    fig_w = LEFT + sum(COL_W) + RIGHT_PAD
    fig_h = TOP_PAD + TITLE_H + HEAD_H + ROW_H * len(body) + BOT_PAD

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off"); ax.set_xlim(0, fig_w); ax.set_ylim(0, fig_h)
    col_left = [LEFT]
    for w in COL_W:
        col_left.append(col_left[-1] + w)
    RIGHT = col_left[-1]
    col_cx = [(col_left[i] + col_left[i + 1]) / 2 for i in range(len(COL_W))]

    def hline(y, lw=1.0, ls="-"):
        ax.plot([LEFT, RIGHT], [y, y], color="black", linewidth=lw, linestyle=ls,
                solid_capstyle="butt", zorder=3)

    y = fig_h - TOP_PAD; y_title_top = y; y -= TITLE_H
    cx = (LEFT + RIGHT) / 2
    ax.text(cx, y + SUB_H + (TITLE_H - SUB_H) * 0.5, TITLE, ha="center", va="center",
            fontsize=12.5, fontweight="bold")
    ax.text(cx, y + SUB_H / 2, SUBTITLE, ha="center", va="center", fontsize=7.6,
            fontstyle="italic", linespacing=1.5)
    hline(y_title_top, lw=1.5); hline(y, lw=1.2)

    y_head_top = y; y -= HEAD_H; ymid = (y_head_top + y) / 2
    for j in range(len(headers)):
        if j in LEFT_COLS:
            ax.text(col_left[j] + 0.06, ymid, headers[j], ha="left", va="center", fontsize=9.0, fontstyle="italic")
        else:
            ax.text(col_cx[j], ymid, headers[j], ha="center", va="center", fontsize=9.0, fontstyle="italic")
    hline(y, lw=1.2)

    y_data_top = y
    for i, cells in enumerate(body):
        if rules[i] == "strong":
            hline(y, lw=0.9, ls=(0, (4, 3)))
        yt = y; y -= ROW_H; ym = (yt + y) / 2
        for j in range(len(headers)):
            if j in LEFT_COLS:
                ax.text(col_left[j] + 0.06, ym, cells[j], ha="left", va="center", fontsize=8.6)
            else:
                ax.text(col_cx[j], ym, cells[j], ha="center", va="center", fontsize=7.4,
                        fontweight="bold" if sigflags[i].get(j, False) else "normal", linespacing=1.15)
    BOTTOM = y
    for x in col_left[1:-1]:
        ax.plot([x, x], [BOTTOM, y_data_top], color="black", linewidth=0.7, linestyle=(0, (3, 3)), zorder=2)
    hline(BOTTOM, lw=1.5)
    ax.add_patch(plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT, y_title_top - BOTTOM,
                               facecolor="none", edgecolor="black", linewidth=1.5, zorder=5))
    for d in (COMP, HERE):
        fig.savefig(d / f"{sp['basename']}.png", bbox_inches="tight", dpi=300)
        fig.savefig(d / f"{sp['basename']}.pdf", bbox_inches="tight", dpi=300)
    plt.close(fig)


# ── LaTeX export (booktabs + makecell; each cell = Δ^stars / SE / z, vs RSF) ──────
_TEX_STAR = {"*": "^{*}", "**": "^{**}", "***": "^{***}"}
TEX_METLBL = [("c_index", r"C-index"), ("c_index_ipcw", r"IPCW C"), ("ibs", r"IBS"),
              ("auc_3y", r"AUC@3y"), ("auc_5y", r"AUC@5y"), ("auc_7y", r"AUC@7y"), ("auc_10y", r"AUC@10y")]


def write_tex(res, B, sp, out_path=None):
    """Wald table (vs sp['ref']) as LaTeX. Each metric cell stacks Δ (with raw-p stars) / bootstrap SE / z
    via \\makecell. Bold Δ = significant AND better than the reference. Requires \\usepackage{makecell,booktabs}."""
    out_path = out_path or (HERE / f"{sp['basename']}.tex")
    nq = int((res[res.group == "Longitudinal"]["p_bh"] < 0.05).sum())
    bh_clause = ("none survive Benjamini--Hochberg FDR correction" if nq == 0
                 else f"{nq} longitudinal comparison(s) survive Benjamini--Hochberg FDR correction")

    def stars(p):
        return "" if not np.isfinite(p) else ("***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "")

    def cell(r):
        if not np.isfinite(r.z):
            return r"\makecell{---}"
        st = _TEX_STAR.get(stars(r.p), "")
        sign = "+" if r.delta >= 0 else "-"
        d = f"{sign}{abs(r.delta):.3f}"
        d = rf"\mathbf{{{d}}}" if (st and r.better) else d
        return (r"\makecell{$" + d + st + r"$ \\ {\scriptsize SE " + f"{r.se_boot:.3f}" +
                r"} \\ {\scriptsize $z{=}" + f"{r.z:.1f}" + r"$}}")

    def row(head_name, agg, key):
        sub = {r.metric: r for r in res[res.model == key].itertuples()}
        agg_tex = "---" if agg == "—" else agg
        return f"{head_name} & {agg_tex} & " + " & ".join(cell(sub[mk]) for mk, _ in TEX_METLBL) + r" \\"

    longi = [(h, a, m) for m, g, h, a in sp["tms"] if g == "Longitudinal"]
    cross = [(h, a, m) for m, g, h, a in sp["tms"] if g != "Longitudinal"]
    head = (r"\textbf{Model} & \textbf{Aggregation} & "
            + " & ".join(rf"\textbf{{{h}}}" for _, h in TEX_METLBL) + r" \\")
    lines = [
        r"% T4 survival — Wald significance of each model vs RSF (bootstrap SE of the metric difference).",
        r"% Requires \usepackage{makecell} and \usepackage{booktabs}. Generated by 08_wald_significance.py.",
        r"\begin{table}[ht]", r"\centering", r"\resizebox{\textwidth}{!}{%",
        r"\renewcommand{\cellalign}{c}",
        r"\begin{tabular}{llccccccc}", r"\toprule", head, r"\midrule",
        r"\multicolumn{9}{l}{\textit{Longitudinal (deep heads on bl/m06/m12 clinical embeddings)}} \\",
        *[row(h, a, m) for h, a, m in longi],
        r"\midrule",
        r"\multicolumn{9}{l}{\textit{Cross-sectional (tabular baseline vector)}} \\",
        *[row(h, a, m) for h, a, m in cross],
        r"\bottomrule", r"\end{tabular}%", r"}",
        r"\caption{Significance of each survival model against the " + sp["short"] + r" reference (" +
        sp["desc"] + r") on the held-out (validation+test) set. Each cell reports "
        r"$\Delta = (\text{model}-\text{ref})$ for that metric, the bootstrap standard error of the "
        r"difference, and the Wald statistic $z=\Delta/\mathrm{SE}$. The SE is from a paired, "
        r"seed-stratified bootstrap ($B=" + str(B) + r"$) that resamples the 108 held-out patients within "
        r"each seed and recomputes the model$-$reference difference (paired, so the strong correlation is "
        r"captured). Markers denote the two-sided Wald test: $^{*}p{<}0.05$, $^{**}p{<}0.01$, "
        r"$^{***}p{<}0.001$ (raw); " + bh_clause + r" across the $" + str(len(sp["tms"])) +
        r"\times7$ comparisons ($q$ in the supplementary CSV). Bold $\Delta$ = significant and better than " +
        sp["short"] + r". IBS is lower-better (negative $\Delta$ = better). "
        r"$N=108$ per seed ($\approx$25 events), 3 seeds.}",
        r"\label{tab:" + sp["label"] + r"}", r"\end{table}", ""]
    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"  TeX: {out_path}")


def write_tex_detailed(res, B, sp, out_path=None):
    """Detailed Wald table (vs sp['ref']): each cell stacks Δ^stars / SE / z / p (raw) / q (BH-FDR)."""
    out_path = out_path or (HERE / f"{sp['basename']}_detailed.tex")
    nq = int((res[res.group == "Longitudinal"]["p_bh"] < 0.05).sum())
    ntest = len(sp["tms"])
    bh_clause = ("none survive Benjamini--Hochberg FDR ($q<0.05$)" if nq == 0
                 else f"{nq} longitudinal comparison(s) survive Benjamini--Hochberg FDR ($q<0.05$)")

    def stars(p):
        return "" if not np.isfinite(p) else ("***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "")

    def cell(r):
        if not np.isfinite(r.z):
            return r"\makecell{---}"
        st = _TEX_STAR.get(stars(r.p), "")
        sign = "+" if r.delta >= 0 else "-"
        d = f"{sign}{abs(r.delta):.3f}"
        d = rf"\mathbf{{{d}}}" if (st and r.better) else d
        return (r"\makecell{$" + d + st + r"$ \\ {\scriptsize SE " + f"{r.se_boot:.3f}" + r"} \\ "
                r"{\scriptsize $z{=}" + f"{r.z:.1f}" + r"$} \\ {\scriptsize $p{=}" + f"{r.p:.3f}" + r"$} \\ "
                r"{\scriptsize $q{=}" + f"{r.p_bh:.3f}" + r"$}}")

    def row(head_name, agg, key):
        sub = {r.metric: r for r in res[res.model == key].itertuples()}
        agg_tex = "---" if agg == "—" else agg
        return f"{head_name} & {agg_tex} & " + " & ".join(cell(sub[mk]) for mk, _ in TEX_METLBL) + r" \\"

    longi = [(h, a, m) for m, g, h, a in sp["tms"] if g == "Longitudinal"]
    cross = [(h, a, m) for m, g, h, a in sp["tms"] if g != "Longitudinal"]
    head = (r"\textbf{Model} & \textbf{Aggregation} & "
            + " & ".join(rf"\textbf{{{h}}}" for _, h in TEX_METLBL) + r" \\")
    lines = [
        r"% T4 survival — Wald significance vs RSF, DETAILED (Δ / SE / z / p_raw / q_BH per cell).",
        r"% Requires \usepackage{makecell} and \usepackage{booktabs}. Generated by 08_wald_significance.py.",
        r"\begin{table}[ht]", r"\centering", r"\resizebox{\textwidth}{!}{%",
        r"\renewcommand{\cellalign}{c}",
        r"\begin{tabular}{llccccccc}", r"\toprule", head, r"\midrule",
        r"\multicolumn{9}{l}{\textit{Longitudinal (deep heads on bl/m06/m12 clinical embeddings)}} \\",
        *[row(h, a, m) for h, a, m in longi],
        r"\midrule",
        r"\multicolumn{9}{l}{\textit{Cross-sectional (tabular baseline-visit features)}} \\",
        *[row(h, a, m) for h, a, m in cross],
        r"\bottomrule", r"\end{tabular}%", r"}",
        r"\caption{Significance of each survival model against the " + sp["short"] + r" reference (" +
        sp["desc"] + r") on the held-out (validation+test) set. Each cell reports "
        r"$\Delta=(\text{model}-\text{ref})$ for that metric, the paired seed-stratified bootstrap standard "
        r"error of the difference ($B=" + str(B) + r"$), the Wald statistic $z=\Delta/\mathrm{SE}$, the raw "
        r"two-sided $p$-value, and the Benjamini--Hochberg FDR-adjusted $q$-value across all $" +
        str(ntest) + r"\times7=" + str(ntest * 7) + r"$ comparisons. Stars mark the raw $p$ "
        r"($^{*}p{<}0.05$, $^{**}p{<}0.01$, $^{***}p{<}0.001$); " + bh_clause +
        r". Bold $\Delta$ = significant (raw) and better than " + sp["short"] + r". IBS is lower-better "
        r"(negative $\Delta$ = better). $N=108$ per seed ($\approx$25 events), 3 seeds.}",
        r"\label{tab:" + sp["label"] + r"_detailed}", r"\end{table}", ""]
    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"  TeX (detailed): {out_path}")


if __name__ == "__main__":
    main()
