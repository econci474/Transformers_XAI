r"""
23_cox_survival_models.py
=========================
Cox proportional-hazards survival models for *time to clinical conversion*
on the ADNI SNP cohort, the survival counterpart to the binary genetic
baselines (script 20).  Answers: given covariates measured at the baseline
visit, how well can we predict *when* a patient converts?

Three at-risk transitions (event = first **sustained** transition: reaches
the target state AND is confirmed there at the last visit — transient blips
are censored, mirroring the strict labels from script 22):

  CN_to_MCI    at-risk bl_dx==CN     event ever_conversion_MCI
               time   first_conversion_MCI (→AD fallback) | censor followup
  MCI_to_AD    at-risk bl_dx==MCI    event progressive_MCI
               time   first_conversion_AD                 | censor followup
  CNMCI_to_AD  at-risk bl_dx∈{CN,MCI} event ever_conversion_AD
               time   first_conversion_AD                 | censor followup

Covariate grid (atoms: age, sex, APOE4[clinical ε4 dosage], and 3 harmonised
PRS tiers betaPRS_withAPOE[all 430] / betaPRS_onlyAPOE[166 isolated chr19
APOE LD block] / betaPRS_noAPOE[263, APOE region removed]): every non-empty
subset — singles, pairwise, joint — under the APOE-aware rules:
  • at most one PRS variant per model;
  • APOE4 + an APOE-region PRS (withAPOE/onlyAPOE) is excluded as a
    double-count; APOE4 + betaPRS_noAPOE is allowed.

Reuses (importlib, unmodified) from 20_genetic_baselines.py:
  build_beta_table, snps_for_tier, _covars, DOSAGE_TSV, SPLITS,
  E_SNPS_DEFAULT, APOE_CLUSTER_DEFAULT.  `_prs` is the 4-line cohort-mean
  weighted sum, replicated here (identical to script 20).

Per (transition, covset, seed): fit lifelines CoxPHFitter (penalizer=1e-3,
covariates z-standardised on train) on the seed's train∪val at-risk patients,
evaluate on the held-out test partition.  Metrics on test: Harrell's
C-index, Royston–Sauerbrei R²_D (explained variation from the prognostic
index), time-dependent AUC at 3/5/10 yr (Uno IPCW), calibration-in-the-
large observed/expected ratio at 3/5/10 yr, and the **hazard ratio per
1 SD** of every covariate with 95% CI (e.g. "each +1 SD non-APOE PRS →
HR 1.25 for AD conversion"; covariates are z-scored on train so exp(coef)
is already the per-SD HR). Seeds 0/1/2 → mean ± std rendered by 23b.

Output: …\cox_survival\<transition>\<covset>\seed_<n>\metrics.json

Usage:
  conda run -n snp python snp_pipeline/23_cox_survival_models.py --sanity
  conda run -n snp python snp_pipeline/23_cox_survival_models.py
"""
from __future__ import annotations

import argparse
import importlib.util as _il
import itertools
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

warnings.filterwarnings("ignore")

# ── Reuse script 20 (digit-leading filename → importlib) ────────────────────
_S20 = Path(__file__).parent / "20_genetic_baselines.py"
_spec = _il.spec_from_file_location("script20", _S20)
_s20 = _il.module_from_spec(_spec)
_spec.loader.exec_module(_s20)
build_beta_table = _s20.build_beta_table
snps_for_tier = _s20.snps_for_tier
_covars = _s20._covars
DOSAGE_TSV = _s20.DOSAGE_TSV
SPLITS = _s20.SPLITS
E_SNPS_DEFAULT = _s20.E_SNPS_DEFAULT
APOE_CLUSTER_DEFAULT = _s20.APOE_CLUSTER_DEFAULT
_load_wightman_allowlist = _s20._load_wightman_allowlist

LABELS_TSV = Path("D:/ADNI_SNP_Omni2.5M_20140220/conversion_labels/"
                  "conversion_labels.tsv")
OUT_DEFAULT = Path("D:/ADNI_SNP_Omni2.5M_20140220/cox_survival")
SPLIT_ROOT = SPLITS["ever_convert"]           # 616-patient universe (membership)
MIN_T = 1.0 / 365.25                          # floor for non-positive durations

# transition → (at-risk bl_dx set, event col, primary time col, fallback time)
TRANSITIONS = {
    "CN_to_MCI":   ({"CN"},       "ever_conversion_MCI",
                    "first_conversion_MCI", "first_conversion_AD"),
    "MCI_to_AD":   ({"MCI"},      "progressive_MCI",
                    "first_conversion_AD", None),
    "CNMCI_to_AD": ({"CN", "MCI"}, "ever_conversion_AD",
                    "first_conversion_AD", None),
}
ATOMS = ["age", "sex", "APOE4", "APOE2", "betaPRS_withAPOE",
         "betaPRS_onlyAPOE", "betaPRS_noAPOE"]
# APOE-region PRS variants (contain the chr19 APOE block) — excluded
# alongside a clinical ε-allele to avoid double-counting. The filtered
# PRS sets (6W27B/9W27B/9W27B5D) and noAPOE are APOE-region-free.
APOE_PRS = {"betaPRS_withAPOE", "betaPRS_onlyAPOE"}
APOE_GENE = {"APOE4", "APOE2"}                        # clinical ε-allele counts
KAPPA = np.sqrt(8.0 / np.pi)                   # Royston–Sauerbrei scale
PRS_PENALIZER = 1e-3                            # single-/few-covariate models
# Raw-dosage blocks are 29–37 LD-correlated SNP columns on a ~tens-of-events
# train fold. Structural rank-deficiency (monomorphic-in-fold SNPs →
# zero-variance columns; perfect-LD SNP pairs → exact-duplicate columns)
# is removed first (drop train-constant, collapse exact duplicates — the
# survival analog of L2-penalised LogReg / collapsing perfect-LD tags),
# then a moderate ridge handles residual near-collinearity. Applied ONLY
# to covsets containing a dosage_<SET> block; PRS/clinical models unchanged
# (those columns are continuous → pruning is a no-op → byte-identical).
DOSAGE_PENALIZER = 1.0


def covsets(atoms: list[str] | None = None) -> list[list[str]]:
    """Non-empty subsets of `atoms` (default ATOMS), APOE-aware filtered,
    ordered by size: ≤1 PRS/dosage signature (any betaPRS* OR dosage_* raw-
    dosage block) per model; a clinical APOE ε-allele (APOE4/APOE2) + an
    APOE-region PRS (withAPOE/onlyAPOE) is a double-count and excluded.
    APOE4/APOE2 + a noAPOE/6W27B/9W27B(5D) PRS or dosage block is allowed
    (those sets are APOE-region-free); APOE2+APOE4 together is allowed."""
    atoms = atoms if atoms is not None else ATOMS
    sig = {a for a in atoms
           if str(a).startswith("betaPRS") or str(a).startswith("dosage_")}
    out = []
    for k in range(1, len(atoms) + 1):
        for c in itertools.combinations(atoms, k):
            s = set(c)
            if len(s & sig) > 1:                  # ≤1 PRS/dosage signature
                continue
            if (s & APOE_GENE) and (s & APOE_PRS):         # APOE double-count
                continue
            out.append(list(c))
    return out


def _expand(atoms: list[str], cols) -> list[str]:
    """Map covset atoms → real design-matrix columns. A `dosage_<SET>` block
    atom expands to every raw per-SNP column `dos__<SET>__*` in `cols` (the
    survival analog of the script-20 `dosage_<set>` LogReg feature); every
    other atom is itself a single column. Order/uniqueness preserved."""
    out: list[str] = []
    for a in atoms:
        if str(a).startswith("dosage_"):
            S = str(a).split("dosage_", 1)[1]
            out += [c for c in cols if str(c).startswith(f"dos__{S}__")]
        elif a not in out:
            out.append(a)
    return out


def _prs(bt, dos, snp_set, weight_col="beta_A1") -> pd.Series:
    """Identical to 20_genetic_baselines.py main()._prs (cohort-mean fill)."""
    w = bt.loc[snp_set, weight_col].astype(float)
    cols = [s for s in snp_set if s in dos.columns]
    D = dos[cols].apply(lambda c: c.fillna(c.mean()))
    return D.mul(w.reindex(cols).values, axis=1).sum(axis=1)


def r2_d(time, event, pi) -> float:
    """Royston–Sauerbrei R²_D from the prognostic index `pi` (higher=worse)."""
    try:
        from lifelines import CoxPHFitter
        n = len(pi)
        if n < 8 or int(np.sum(event)) < 3 or np.std(pi) < 1e-9:
            return float("nan")
        order = np.argsort(np.argsort(pi))                 # 0..n-1 ranks
        rankit = norm.ppf((order + 1 - 3.0 / 8.0) / (n + 0.25))
        scaled = rankit / KAPPA
        d = pd.DataFrame({"t": time, "e": event, "z": scaled})
        cph = CoxPHFitter(penalizer=1e-3)
        cph.fit(d, duration_col="t", event_col="e")
        D = float(cph.params_["z"])
        v = (D * D) / (KAPPA * KAPPA)
        return float(v / (np.pi ** 2 / 6.0 + v))
    except Exception:
        return float("nan")


def calibration(cph, Xte, te_t, te_e, times=(3.0, 5.0, 10.0)) -> dict:
    """Calibration-in-the-large at fixed horizons: observed/expected event
    probability ratio. Observed = Kaplan–Meier event prob at t in the test
    fold; Expected = mean model-predicted event prob 1−Ŝ(t|x). O/E≈1 ⇒ well
    calibrated, >1 ⇒ risk under-predicted, <1 ⇒ over-predicted. Robust at
    small n (no binning). NaN if t beyond test follow-up / not estimable."""
    res = {}
    for t in times:
        res[f"calib_obs_{int(t)}y"] = float("nan")
        res[f"calib_exp_{int(t)}y"] = float("nan")
        res[f"calib_oe_{int(t)}y"] = float("nan")
    try:
        from lifelines import KaplanMeierFitter
        te_t = np.asarray(te_t, float); te_e = np.asarray(te_e, int)
        km = KaplanMeierFitter().fit(te_t, event_observed=te_e)
        tmax = float(te_t.max())
        sf = cph.predict_survival_function(Xte, times=list(times))
        for t in times:
            k = int(t)
            if t > tmax:                                  # outside follow-up
                continue
            obs = 1.0 - float(km.survival_function_at_times(t).iloc[0])
            exp = float((1.0 - sf.loc[t]).mean())
            res[f"calib_obs_{k}y"] = obs
            res[f"calib_exp_{k}y"] = exp
            if exp > 1e-9:
                res[f"calib_oe_{k}y"] = obs / exp
    except Exception:
        pass
    return res


def hr_per_sd(cph, atoms: list[str]) -> dict:
    """Hazard ratio per 1 SD of each covariate (covariates were z-scored on
    train ⇒ exp(coef) is already the per-1-SD HR), with 95% CI and p."""
    out = {}
    try:
        sm = cph.summary
        for a in atoms:
            out[a] = {"hr": float(sm.loc[a, "exp(coef)"]),
                      "lo": float(sm.loc[a, "exp(coef) lower 95%"]),
                      "hi": float(sm.loc[a, "exp(coef) upper 95%"]),
                      "p": float(sm.loc[a, "p"])}
    except Exception:
        pass
    return out


def td_auc(tr_t, tr_e, te_t, te_e, te_pi, times=(3.0, 5.0, 10.0)) -> dict:
    """Cumulative/dynamic time-dependent AUC at fixed horizons, Uno's IPCW
    estimator: cases = events by t (weighted 1/Ĝ(Tᵢ)), controls = event-free
    past t.  Ĝ = Kaplan–Meier of the censoring distribution estimated on the
    training fold (mirrors sksurv).  NaN if a horizon is infeasible."""
    res = {f"tdauc_{int(t)}y": float("nan") for t in times}
    try:
        from lifelines import KaplanMeierFitter
        tr_t = np.asarray(tr_t, float); tr_e = np.asarray(tr_e, int)
        te_t = np.asarray(te_t, float); te_e = np.asarray(te_e, int)
        te_pi = np.asarray(te_pi, float)
        kmc = KaplanMeierFitter().fit(tr_t, event_observed=1 - tr_e)
        ev_t = te_t[te_e == 1]
        if len(ev_t) < 2:
            return res
        lo, hi = float(ev_t.min()), float(te_t.max())
        for t in times:
            if not (lo < t < hi):
                continue
            case = (te_e == 1) & (te_t <= t)
            ctrl = te_t > t
            if case.sum() < 2 or ctrl.sum() < 1:
                continue
            g = kmc.survival_function_at_times(te_t[case]).to_numpy()
            w = 1.0 / np.clip(g, 1e-3, None)               # IPCW case weights
            mc, mk = te_pi[case], te_pi[ctrl]
            comp = (mc[:, None] > mk[None, :]).astype(float) \
                + 0.5 * (mc[:, None] == mk[None, :])
            num = float((w[:, None] * comp).sum())
            den = float(w.sum() * ctrl.sum())
            if den > 0:
                res[f"tdauc_{int(t)}y"] = num / den
    except Exception:
        pass
    return res


def c_index(time, pi, event) -> float:
    try:
        from lifelines.utils import concordance_index
        if int(np.sum(event)) < 2:
            return float("nan")
        return float(concordance_index(time, -pi, event))  # higher pi = worse
    except Exception:
        return float("nan")


def build_frame(lab: pd.DataFrame, atrisk: set, ev_col: str,
                t_col: str, t_fb: str | None) -> pd.DataFrame:
    """Per-PTID (event, time) for one transition; drops uninformative rows."""
    d = lab[lab["bl_dx"].isin(atrisk)].copy()
    e = d[ev_col].astype(int).to_numpy()
    tp = pd.to_numeric(d[t_col], errors="coerce")
    if t_fb is not None:
        tp = tp.fillna(pd.to_numeric(d[t_fb], errors="coerce"))
    fu = pd.to_numeric(d["followup_years"], errors="coerce")
    time = np.where(e == 1, tp, fu).astype(float)
    out = pd.DataFrame({"Patient_ID": d["Patient_ID"].astype(str).to_numpy(),
                        "event": e, "time": time})
    out = out[np.isfinite(out["time"])]
    out["time"] = out["time"].clip(lower=MIN_T)             # lifelines needs t>0
    return out.reset_index(drop=True)


def covar_matrix(bt, dos, seed: int, enr_w=None, enr_d=None) -> pd.DataFrame:
    """Per-PTID atom covariates for `seed` (covars constant across splits).
    If enr_w/enr_d given (script-24 filtered enriched bundle), betaPRS_noAPOE
    is computed over comb = bt ∪ novel using the enriched dosage — so the
    filtered Cox PRS matches the script-20 filtered baseline (9 W + 27 B
    [+ 5 Desikan for the withdesikan variant])."""
    parts = [pd.read_csv(SPLIT_ROOT / f"seed_{seed}" / f"{sp}.csv")
             for sp in ("train", "val", "test")]
    raw = pd.concat(parts, ignore_index=True).drop_duplicates("Patient_ID")
    cv = _covars(raw)                                  # AGE SEX APOE4 APOE2
    cv.index = cv.index.astype(str)
    base = list(bt.index)

    def _p(tier):
        s = _prs(bt, dos, snps_for_tier(bt, base, tier,
                                        E_SNPS_DEFAULT, APOE_CLUSTER_DEFAULT))
        s.index = s.index.astype(str)
        return s.reindex(cv.index)
    m = pd.DataFrame(index=cv.index)
    m["age"] = cv["AGE"]
    m["sex"] = cv["SEX"]
    m["APOE4"] = cv["APOE4"]
    m["APOE2"] = cv["APOE2"]
    if enr_w and enr_d:                                 # filtered PRS sets
        ew = pd.read_csv(enr_w, sep="\t", low_memory=False)
        edos = pd.read_csv(enr_d, sep="\t")
        edos["PTID"] = edos["PTID"].astype(str)
        edos = edos.set_index("PTID").apply(lambda c: c.fillna(c.mean()))
        comb = bt[["beta_A1", "CHR", "BP"]].copy()
        comb["CHR"] = comb["CHR"].astype(str)
        nv = ew[ew["novel"] == True].dropna(subset=["beta_A1"])  # noqa: E712
        for r in nv.itertuples(index=False):
            comb.loc[str(r.rsID)] = {"beta_A1": float(r.beta_A1),
                                     "CHR": str(r.CHR).split(".")[0],
                                     "BP": int(float(r.BP))}
        comb["BP"] = comb["BP"].astype(int)
        src = ew.get("source", pd.Series([""] * len(ew))).astype(str)
        w_novel = [str(r) for r, sc, nvl in zip(ew["rsID"], src, ew["novel"])
                   if nvl and "wightman" in sc.lower()]
        d_novel = [str(r) for r, sc, nvl in zip(ew["rsID"], src, ew["novel"])
                   if nvl and "desikan" in sc.lower()]

        def _setcols(rsids):             # noAPOE-filtered SNPs present in edos
            ss = snps_for_tier(comb, rsids, "noAPOE",
                               E_SNPS_DEFAULT, APOE_CLUSTER_DEFAULT)
            return [s for s in ss if s in edos.columns]

        def _eset(cc):                   # noAPOE-filtered weighted PRS
            w = comb.loc[cc, "beta_A1"].astype(float)
            p = edos[cc].mul(w.values, axis=1).sum(axis=1)
            p.index = p.index.astype(str)
            return p.reindex(cv.index)
        base = list(bt.index)
        sets = {"6W27B": base, "9W27B": base + w_novel}
        if d_novel:
            sets["9W27B5D"] = base + w_novel + d_novel
        blocks = []
        for S, rsids in sets.items():
            cc = _setcols(rsids)
            m[f"betaPRS_{S}"] = _eset(cc)            # weighted-PRS covariate
            # raw per-SNP dosage block — the survival analog of the
            # script-20 `dosage_<set>` LogReg model: one Cox covariate
            # signature, expanded to its SNP columns in fit_eval.
            dblk = edos[cc].copy()
            dblk.index = dblk.index.astype(str)
            dblk = dblk.reindex(cv.index)
            dblk.columns = [f"dos__{S}__{c}" for c in cc]
            blocks.append(dblk)
        if blocks:
            m = pd.concat([m] + blocks, axis=1)
    else:                                               # unfiltered tiers
        m["betaPRS_withAPOE"] = _p("withAPOE")
        m["betaPRS_onlyAPOE"] = _p("onlyAPOE")
        m["betaPRS_noAPOE"] = _p("noAPOE")
    return m


def split_ids(seed: int) -> dict:
    return {sp: set(pd.read_csv(SPLIT_ROOT / f"seed_{seed}" / f"{sp}.csv",
                                usecols=["Patient_ID"])["Patient_ID"]
                    .astype(str))
            for sp in ("train", "val", "test")}


def fit_eval(frame: pd.DataFrame, cov: pd.DataFrame, atoms: list[str],
             ids: dict) -> dict:
    from lifelines import CoxPHFitter
    cols = _expand(atoms, cov.columns)            # dosage_<SET> → SNP columns
    if not cols:
        return {"status": "insufficient", "n_train": 0, "n_test": 0,
                "n_train_events": 0, "n_test_events": 0}
    df = frame.merge(cov[cols], left_on="Patient_ID", right_index=True,
                     how="inner").dropna(subset=cols)
    # fit on TRAIN only, evaluate on the held-out VALIDATION fold
    # (consistent with the script-20 baselines).
    tr = df[df["Patient_ID"].isin(ids["train"])]
    te = df[df["Patient_ID"].isin(ids["val"])]
    if len(tr) < 10 or len(te) < 5 or tr["event"].sum() < 3 \
            or te["event"].sum() < 2:
        return {"status": "insufficient", "n_train": int(len(tr)),
                "n_test": int(len(te)),
                "n_train_events": int(tr["event"].sum()),
                "n_test_events": int(te["event"].sum())}
    is_dos = any(str(a).startswith("dosage_") for a in atoms)
    if is_dos:
        # remove structural rank-deficiency from the raw-dosage block:
        # drop train-constant (monomorphic-in-fold) columns, then collapse
        # exact-duplicate (perfect-LD) columns to one representative.
        cols = [c for c in cols if tr[c].std(ddof=0) > 1e-12]
        seen, ded = {}, []
        for c in cols:
            key = tuple(np.round(tr[c].to_numpy(float), 6))
            if key not in seen:
                seen[key] = c
                ded.append(c)
        cols = ded
        if not cols:
            return {"status": "insufficient", "n_train": int(len(tr)),
                    "n_test": int(len(te)),
                    "n_train_events": int(tr["event"].sum()),
                    "n_test_events": int(te["event"].sum())}
    mu, sd = tr[cols].mean(), tr[cols].std(ddof=0).replace(0, 1.0)
    Xtr, Xte = (tr[cols] - mu) / sd, (te[cols] - mu) / sd
    dtr = pd.concat([tr[["time", "event"]].reset_index(drop=True),
                     Xtr.reset_index(drop=True)], axis=1)
    pen = DOSAGE_PENALIZER if is_dos else PRS_PENALIZER
    cph = CoxPHFitter(penalizer=pen)
    cph.fit(dtr, duration_col="time", event_col="event")
    pi_te = cph.predict_log_partial_hazard(Xte).to_numpy()
    t_te = te["time"].to_numpy(float)
    e_te = te["event"].to_numpy(int)
    out = {"status": "ok", "n_train": int(len(tr)), "n_test": int(len(te)),
           "n_train_events": int(tr["event"].sum()),
           "n_test_events": int(te["event"].sum()),
           "c_index": c_index(t_te, pi_te, e_te),
           "r2_d": r2_d(t_te, e_te, pi_te)}
    out.update(td_auc(tr["time"].to_numpy(float), tr["event"].to_numpy(int),
                      t_te, e_te, pi_te))
    out.update(calibration(cph, Xte, t_te, e_te))
    # per-SD HR only for single-column atoms; a dosage_<SET> block is a
    # multivariate signature (no single HR), mirroring the baseline
    # dosage_<set> rows that carry no OR/1SD.
    out["hr_per_sd"] = hr_per_sd(
        cph, [a for a in atoms if not str(a).startswith("dosage_")])
    return out


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--labels-tsv", type=Path, default=LABELS_TSV)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--skip-if-exists", action="store_true")
    ap.add_argument("--sanity", action="store_true",
                    help="Per-transition n/events/median-fu + joint-model "
                         "asserts + PH check, then exit (no full grid).")
    ap.add_argument("--wightman-filtered", type=Path, default=None,
                    help="script-24 wightman_resolved_allowlist.tsv → "
                         "Wightman β restricted to resolved leads; ATOMS "
                         "reduced to age/sex/APOE4/APOE2/betaPRS_noAPOE.")
    ap.add_argument("--enriched-weights", type=Path, default=None,
                    help="script-24 wfilt_(desikan_)enriched_weights.tsv")
    ap.add_argument("--enriched-dosage", type=Path, default=None,
                    help="script-24 wfilt_dosage/patient_dosage.tsv")
    args = ap.parse_args()

    wig_allow = _load_wightman_allowlist(args.wightman_filtered)
    lab = pd.read_csv(args.labels_tsv, sep="\t")
    lab["Patient_ID"] = lab["Patient_ID"].astype(str)
    print("[β] building harmonised β table …", flush=True)
    bt = build_beta_table(wightman_allowlist=wig_allow)
    dos = pd.read_csv(DOSAGE_TSV, sep="\t")
    dos["PTID"] = dos["PTID"].astype(str)
    dos = dos.set_index("PTID")

    cov = {sd: covar_matrix(bt, dos, sd, args.enriched_weights,
                            args.enriched_dosage) for sd in args.seeds}
    if wig_allow is not None:
        # filtered: the named PRS sets that covar_matrix actually built,
        # plus one raw-dosage block atom per set (survival analog of the
        # script-20 dosage_<set> LogReg rows; expanded in fit_eval).
        c0 = cov[args.seeds[0]].columns
        prs_cols = [c for c in c0 if c.startswith("betaPRS_")]
        dos_sets = sorted({str(c).split("dos__", 1)[1].split("__", 1)[0]
                           for c in c0 if str(c).startswith("dos__")})
        dos_atoms = [f"dosage_{S}" for S in dos_sets]
        fa = ["age", "sex", "APOE4", "APOE2"] + prs_cols + dos_atoms
        cs = covsets(fa)
        # per-set SNP count (betaPRS_<S> and dosage_<S> share the same
        # noAPOE-filtered SNP set → one count each) for the 23b subtitle.
        prs_set_sizes = {S: len([c for c in c0
                                 if str(c).startswith(f"dos__{S}__")])
                         for S in dos_sets}
    else:
        cs = covsets()
        prs_set_sizes = {}
    print(f"[grid] {len(TRANSITIONS)} transitions × {len(cs)} covsets × "
          f"{len(args.seeds)} seeds = "
          f"{len(TRANSITIONS)*len(cs)*len(args.seeds)} fits", flush=True)

    if args.sanity:
        from lifelines import CoxPHFitter
        for tn, (ar, ev, tc, fb) in TRANSITIONS.items():
            fr = build_frame(lab, ar, ev, tc, fb)
            med = float(np.median(fr["time"]))
            print(f"\n[{tn}] at-risk={len(fr)}  events={int(fr['event'].sum())}"
                  f"  median follow-up/time={med:.2f} yr")
            assert fr["event"].sum() > 0, f"{tn}: zero events"
            j = next(c for c in cs if len(c) == max(len(x) for x in cs))
            r = fit_eval(fr, cov[args.seeds[0]], j, split_ids(args.seeds[0]))
            print(f"   joint covset {j} → {r}")
            if r["status"] == "ok":
                assert 0.0 < r["c_index"] < 1.0, f"{tn}: C-index out of range"
                assert np.isnan(r["r2_d"]) or 0.0 <= r["r2_d"] <= 1.0, \
                    f"{tn}: R²_D out of [0,1]"
                # PH-assumption diagnostic (reported, not gated)
                jc = _expand(j, cov[args.seeds[0]].columns)
                df = fr.merge(cov[args.seeds[0]][jc], left_on="Patient_ID",
                              right_index=True, how="inner").dropna(subset=jc)
                mu = df[jc].mean(); sd2 = df[jc].std(ddof=0).replace(0, 1.0)
                dd = pd.concat([df[["time", "event"]].reset_index(drop=True),
                                ((df[jc] - mu) / sd2).reset_index(drop=True)],
                               axis=1)
                c = CoxPHFitter(penalizer=1e-3)
                c.fit(dd, "time", "event")
                print(f"   [PH check {tn}]")
                try:
                    c.check_assumptions(dd, p_value_threshold=0.05,
                                        show_plots=False)
                except Exception as e:
                    print(f"   (check_assumptions skipped: {e})")
        print("\n[sanity] OK — events>0, C-index∈(0,1), R²_D∈[0,1] "
              "(no full grid fitted)")
        return

    done = skipped = failed = 0
    for tn, (ar, ev, tc, fb) in TRANSITIONS.items():
        fr = build_frame(lab, ar, ev, tc, fb)
        for atoms in cs:
            name = "+".join(atoms)
            for sd in args.seeds:
                od = args.out / tn / name / f"seed_{sd}"
                mf = od / "metrics.json"
                if args.skip_if_exists and mf.exists():
                    skipped += 1
                    continue
                try:
                    r = fit_eval(fr, cov[sd], atoms, split_ids(sd))
                    od.mkdir(parents=True, exist_ok=True)
                    json.dump({"transition": tn, "covset": name,
                               "atoms": atoms, "seed": sd,
                               "n_atrisk": int(len(fr)),
                               "n_events": int(fr["event"].sum()),
                               "prs_set_sizes": prs_set_sizes, **r},
                              open(mf, "w"), indent=2)
                    if r["status"] == "ok":
                        done += 1
                        print(f"[ok] {tn}/{name}/seed{sd}  "
                              f"C={r['c_index']:.3f} R2D={r['r2_d']:.3f} "
                              f"AUC3={r['tdauc_3y']:.3f} AUC5={r['tdauc_5y']:.3f}"
                              f" (ntr={r['n_train']},ev={r['n_test_events']})",
                              flush=True)
                    else:
                        skipped += 1
                        print(f"[skip] {tn}/{name}/seed{sd}: {r['status']} "
                              f"(ntr={r['n_train']},evtr={r['n_train_events']})",
                              flush=True)
                except Exception as e:
                    failed += 1
                    print(f"[FAIL] {tn}/{name}/seed{sd}: "
                          f"{type(e).__name__}: {e}", flush=True)
    print(f"\n[done] {done} ok, {skipped} skipped, {failed} failed → {args.out}")


if __name__ == "__main__":
    main()
