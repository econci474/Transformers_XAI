r"""
20c_genetic_baselines_resolved.py   (LOCAL, env: snp)
=====================================================
Classical genetic baselines for the **8 MAF-resolved PRS combinations**
exported by `25d_export_maf_resolved_sets.py` to
`D:\…\GWAS_filtered_maf_resolved\<COMBO>\`. Parallels script 20 but
consumes the new self-contained `(snps.tsv + patient_dosage.tsv)` trio
per combo — no `build_beta_table` / no `external_gwas_labels.tsv` /
no APOE-tier partitioning (each combo is already the noAPOE-equivalent
resolved set; published APOE-region SNPs are inert post-MAF-filter).

Per combo (8): sourceB(26) sourceW(8) sourceD(15) sourceK(3) 5W26B(31)
5W26B14D(45) 5W26B13D(44) 5W26B14D2K(47).

Per combo, the runs are:
  betaPRS_<combo>                  Σ β·dosage, z-scored on TRAIN
  betaPRS_<combo>+age              + AGE covariate
  betaPRS_<combo>+sex              + SEX
  betaPRS_<combo>+age+sex          + AGE,SEX
  betaPRS_<combo>+APOE4+age+sex    + clinical APOE4_Dosage
  betaPRS_<combo>+APOE2+APOE4+age+sex  + clinical APOE2,APOE4
  dosage_<combo>                   LR on raw dosage matrix (no β)

Plus the global APOE/demographic reference rows (single set, not
per-combo, to keep the table compact):
  APOE4_dosage                     clinical APOE4 only
  APOE4_dosage+age, +sex, +age+sex
  APOE2_dosage, +covars
  APOE2+APOE4_dosage, +covars
  age, sex, age+sex (demographics only)

Output schema matches script 20 so script 20b can render it:
  <out-root>/<label_mode>/<model>/<feature_set>/seed_<n>/metrics.json

Usage:
  conda run -n snp python snp_pipeline/20c_genetic_baselines_resolved.py \
    --label-mode ad_case \
    --labels-tsv "D:/ADNI_SNP_Omni2.5M_20140220/conversion_labels.tsv" \
    --out-root  "D:/ADNI_SNP_Omni2.5M_20140220/genetic_baselines_resolved"
"""
from __future__ import annotations

import argparse
import importlib.util as _il
import json
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

# ── importlib reuse of script 20 (digit-leading filename) ───────────────────
_S20 = Path(__file__).parent / "20_genetic_baselines.py"
_spec = _il.spec_from_file_location("s20", _S20)
_s20 = _il.module_from_spec(_spec)
_spec.loader.exec_module(_s20)
_covars = _s20._covars
_metrics = _s20._metrics
_strict_labels = _s20._strict_labels
load_labels = _s20.load_labels
SPLITS = _s20.SPLITS
STRICT_MODES = _s20.STRICT_MODES

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
COMBO_ROOT_DEFAULT = BASE / "GWAS_filtered_maf_resolved"
COMBOS = ["sourceB", "sourceW", "sourceD", "sourceK",
          "5W26B", "5W26B14D", "5W26B13D", "5W26B14D2K"]

COVAR = {"+age": ["AGE"], "+sex": ["SEX"], "+age+sex": ["AGE", "SEX"]}
PRS_COMBOS = [("", []), ("+age", ["AGE"]), ("+sex", ["SEX"]),
              ("+age+sex", ["AGE", "SEX"]),
              ("+APOE4+age+sex", ["APOE4", "AGE", "SEX"]),
              ("+APOE2+APOE4+age+sex", ["APOE2", "APOE4", "AGE", "SEX"])]


def load_combo(combo_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (snps_df: rsID→{beta_A1,CHR,BP_GRCh38,A1,A2,...},
                dosage_df: PTID-indexed with one col per rsID, 0/1/2)."""
    snp = pd.read_csv(combo_dir / f"{combo_dir.name}_snps.tsv", sep="\t")
    dos = pd.read_csv(combo_dir / f"{combo_dir.name}_patient_dosage.tsv",
                       sep="\t")
    dos["PTID"] = dos["PTID"].astype(str)
    dos = dos.set_index("PTID")
    # cohort-mean impute for the PRS reference standardisation (TRAIN-fold
    # impute is re-done inside each fit; this is for the global PRS series).
    dos = dos.apply(lambda c: c.fillna(c.mean()))
    snp = snp.set_index("rsID")
    return snp, dos


def compute_prs(snp: pd.DataFrame, dos: pd.DataFrame) -> pd.Series:
    """Σ_i β_A1_i · dosage_{p,i} over the combo's SNPs (cohort-mean imputed)."""
    cols = [s for s in snp.index if s in dos.columns]
    w = snp.loc[cols, "beta_A1"].astype(float).values
    return dos[cols].astype(float).mul(w, axis=1).sum(axis=1)


def fit_eval(parts: dict[str, pd.DataFrame], Xcols: list[str], seed: int) -> dict:
    """Fit LR on train, evaluate on val/test. Same pipeline as script 20."""
    tr = parts["train"].reset_index(drop=True)
    if not Xcols or tr.empty:
        raise ValueError(f"empty train / no Xcols={Xcols}")
    pipe = Pipeline([("imp", SimpleImputer(strategy="median")),
                     ("scl", StandardScaler()),
                     ("clf", LogisticRegression(max_iter=2000,
                                                 class_weight="balanced",
                                                 random_state=seed))])
    pipe.fit(tr[Xcols].astype(float), tr["y"].astype(int))
    te, vp = parts["test"], parts["val"]
    test_m = _metrics(te["y"].astype(int).values,
                       pipe.predict_proba(te[Xcols].astype(float))[:, 1])
    val_m = _metrics(vp["y"].astype(int).values,
                      pipe.predict_proba(vp[Xcols].astype(float))[:, 1])
    notes = {}
    if "PRS" in Xcols:
        j = Xcols.index("PRS")
        coef = float(pipe.named_steps["clf"].coef_[0][j])
        scl = float(pipe.named_steps["scl"].scale_[j]) or 1.0
        train_sd = float(tr["PRS"].astype(float).std(ddof=0))
        logor_per_unit = coef / scl
        notes = {"train_prs_sd": train_sd,
                  "logor_per_unit": logor_per_unit,
                  "or_per_train_sd": math.exp(logor_per_unit * train_sd)}
    return val_m, test_m, notes


def build_parts(feat: pd.DataFrame | None, split_cols: list[str],
                lab: pd.DataFrame, cov: pd.DataFrame,
                apoe_stratum: str | None = None) -> pd.DataFrame:
    """Assemble (y, [PRS columns from feat], [covariate columns from cov]).
    `apoe_stratum`: None (no filter) | 'e33' (APOE4==0 AND APOE2==0,
    ie ε3/ε3) | 'e4_carrier' (APOE4>=1) | 'e2_carrier' (APOE2>=1).
    Filtering happens on the cov frame BEFORE the per-patient assembly,
    so partitions shrink consistently across splits."""
    rows = []
    for pid, y in zip(lab["Patient_ID"].astype(str), lab["y"]):
        if feat is not None and pid not in feat.index:
            continue
        if apoe_stratum and pid in cov.index:
            a4 = float(cov.loc[pid, "APOE4"]) if "APOE4" in cov.columns else float("nan")
            a2 = float(cov.loc[pid, "APOE2"]) if "APOE2" in cov.columns else float("nan")
            if apoe_stratum == "e33":
                if not (np.isfinite(a4) and np.isfinite(a2)
                         and a4 == 0 and a2 == 0):
                    continue
            elif apoe_stratum == "e4_carrier":
                if not (np.isfinite(a4) and a4 >= 1):
                    continue
            elif apoe_stratum == "e2_carrier":
                if not (np.isfinite(a2) and a2 >= 1):
                    continue
        rec = {"y": int(y)}
        if feat is not None:
            rec.update(feat.loc[pid].to_dict())
        for cc in split_cols:
            v = float(cov.loc[pid, cc]) if pid in cov.index else float("nan")
            rec[cc] = v if np.isfinite(v) else float("nan")
        rows.append(rec)
    return pd.DataFrame(rows)


def _pn(p: pd.DataFrame) -> tuple[int, int]:
    y = p["y"].astype(int) if len(p) else pd.Series([], dtype=int)
    return int((y == 1).sum()), int((y == 0).sum())


def run_one(model: str, fset: str, feat: pd.DataFrame | None,
            n_used: int, split_cols: list[str],
            labels: dict, raw_splits: dict, args, out_root: Path) -> dict | None:
    out = {}
    for sd in args.seeds:
        out_dir = out_root / args.label_mode / model / fset / f"seed_{sd}"
        mfile = out_dir / "metrics.json"
        if args.skip_if_exists and mfile.exists():
            continue
        try:
            parts = {}
            for sp in ("train", "val", "test"):
                lab = labels[sd][sp].copy()
                cov = _covars(raw_splits[sd][sp])
                cov.index = cov.index.astype(str)
                parts[sp] = build_parts(feat, split_cols, lab, cov,
                                         apoe_stratum=args.apoe_stratum)
            Xcols = [c for c in parts["train"].columns if c != "y"]
            val_m, test_m, notes = fit_eval(parts, Xcols, sd)
            trp, vpn, tep = _pn(parts["train"]), _pn(parts["val"]), _pn(parts["test"])
            out_dir.mkdir(parents=True, exist_ok=True)
            json.dump({"label_mode": args.label_mode, "model": model,
                       "method": model, "aggregation": fset, "head": model,
                       "seed": sd, "n_snps": int(n_used),
                       "n_train_pos": trp[0], "n_train_neg": trp[1],
                       "n_val_pos": vpn[0], "n_val_neg": vpn[1],
                       "n_test_pos": tep[0], "n_test_neg": tep[1],
                       "val": val_m, "test": test_m, "notes": notes},
                       open(mfile, "w"), indent=2)
            out[sd] = (test_m["auc"], test_m["balanced_accuracy"])
            print(f"  [ok] {model}/{fset}/seed{sd}  "
                  f"test_auc={test_m['auc']:.3f} "
                  f"bACC={test_m['balanced_accuracy']:.3f} (n={n_used})",
                  flush=True)
        except Exception as e:
            print(f"  [FAIL] {model}/{fset}/seed{sd}: {type(e).__name__}: {e}",
                  flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--combo-root", type=Path, default=COMBO_ROOT_DEFAULT)
    ap.add_argument("--combos", nargs="+", default=COMBOS)
    ap.add_argument("--out-root", type=Path,
                     default=BASE / "genetic_baselines_resolved")
    ap.add_argument("--label-mode", default="ad_case",
                     choices=["ever_convert", *STRICT_MODES.keys()])
    ap.add_argument("--labels-tsv", type=Path,
                     default=BASE / "conversion_labels.tsv")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--skip-if-exists", action="store_true")
    ap.add_argument("--apoe-stratum", default=None,
                     choices=[None, "e33", "e4_carrier", "e2_carrier"],
                     help="Restrict cohort to APOE stratum. e33 = clinical "
                          "APOE4==0 AND APOE2==0 (ε3/ε3, ~55%% of Europeans, "
                          "no APOE confounding); e4_carrier = APOE4>=1; "
                          "e2_carrier = APOE2>=1. Default None (full cohort).")
    args = ap.parse_args()

    split_root = SPLITS.get(args.label_mode) or SPLITS["ever_convert"]
    if args.label_mode in STRICT_MODES and args.labels_tsv is None:
        raise SystemExit(f"--label-mode {args.label_mode} requires --labels-tsv")

    args.out_root.mkdir(parents=True, exist_ok=True)
    print(f"[setup] label_mode={args.label_mode}  out_root={args.out_root}")
    print(f"[setup] combos: {args.combos}")
    print(f"[setup] seeds: {args.seeds}")

    # Labels + raw splits (shared across combos)
    if args.label_mode in STRICT_MODES:
        labels = _strict_labels(split_root, args.labels_tsv,
                                args.label_mode, args.seeds)
    else:
        labels = {sd: load_labels(split_root, sd,
                                   label_mode=args.label_mode)
                  for sd in args.seeds}
    raw_splits = {sd: {sp: pd.read_csv(split_root / f"seed_{sd}" / f"{sp}.csv")
                       for sp in ("train", "val", "test")}
                   for sd in args.seeds}

    # ── Stratum-aware RUN filtering ──────────────────────────────────────
    # In e33 (APOE4==0 AND APOE2==0), any feature that references APOE4
    # or APOE2 is a constant zero column → uninformative. We skip those
    # runs entirely so the rendered table only contains non-degenerate
    # rows (demog, β-PRS alone, β-PRS + age/sex/age+sex, raw-dosage LR).
    def _is_degenerate(split_cols, model_kind=None):
        if args.apoe_stratum != "e33":
            return False
        if model_kind in {"apoe_e4", "apoe_e2", "apoe_e2e4", "apoe_covar"}:
            return True
        return any(c in {"APOE4", "APOE2"} for c in (split_cols or []))

    # PRS_COMBOS pruned by stratum (e33 drops +APOE4 / +APOE2+APOE4 variants)
    prs_combos = [(suf, cc) for suf, cc in PRS_COMBOS
                   if not _is_degenerate(cc)]
    if args.apoe_stratum == "e33":
        skipped_prs = [s for s, _ in PRS_COMBOS
                        if (s, _) not in prs_combos]
        print(f"[e33] dropping degenerate PRS_COMBOS suffixes: "
              f"{[s for s, _ in PRS_COMBOS if any(c in ('APOE4','APOE2') for c in _)]}")

    # ── 1) Per-combo betaPRS_z + (×covars) + raw-dosage LR ────────────────
    summary = []
    for combo in args.combos:
        combo_dir = args.combo_root / combo
        if not combo_dir.exists():
            print(f"[skip] {combo}: {combo_dir} not found")
            continue
        print(f"\n=== Combo: {combo} ===")
        snp, dos = load_combo(combo_dir)
        prs = compute_prs(snp, dos)
        prs_df = pd.DataFrame({"PRS": prs})
        n_used = int(len(snp))
        # PRS_COMBOS family (alone + 5 covariate combos; e33-filtered)
        for suffix, cc in prs_combos:
            fset = f"betaPRS_{combo}{suffix}"
            model = "prs" if not suffix else "prs_covar"
            r = run_one(model, fset, prs_df, n_used, cc,
                         labels, raw_splits, args, args.out_root)
            if r:
                for sd, (auc, bacc) in r.items():
                    summary.append({"combo": combo, "model": model,
                                     "feature_set": fset, "seed": sd,
                                     "test_auc": auc, "test_bACC": bacc})
        # Raw-dosage LR (all combo SNPs as features)
        dos_used = [s for s in snp.index if s in dos.columns]
        if dos_used:
            r = run_one("logreg_dosage", f"dosage_{combo}",
                         dos[dos_used], len(dos_used), [],
                         labels, raw_splits, args, args.out_root)
            if r:
                for sd, (auc, bacc) in r.items():
                    summary.append({"combo": combo, "model": "logreg_dosage",
                                     "feature_set": f"dosage_{combo}",
                                     "seed": sd, "test_auc": auc,
                                     "test_bACC": bacc})

    # ── 2) APOE + demographic reference rows (single set, all combos share) ──
    # e33 stratum: skip every APOE-touching row; keep only demog.
    print("\n=== APOE + demographic reference rows ===")
    apoe_runs = []
    apoe_runs.append(("apoe_e4",   "APOE4_dosage",         ["APOE4"]))
    apoe_runs.append(("apoe_e2",   "APOE2_dosage",         ["APOE2"]))
    apoe_runs.append(("apoe_e2e4", "APOE2+APOE4_dosage",   ["APOE2", "APOE4"]))
    for cn, cc in COVAR.items():
        apoe_runs.append(("apoe_covar", f"APOE4_dosage{cn}", ["APOE4"] + cc))
        apoe_runs.append(("apoe_e2",   f"APOE2_dosage{cn}", ["APOE2"] + cc))
        apoe_runs.append(("apoe_e2e4", f"APOE2+APOE4_dosage{cn}",
                          ["APOE2", "APOE4"] + cc))
    for cn, cc in COVAR.items():
        apoe_runs.append(("demog", cn.lstrip("+"), cc))
    # Filter degenerate runs out in e33
    apoe_runs = [(m, f, c) for m, f, c in apoe_runs
                  if not _is_degenerate(c, model_kind=m)]
    for model, fset, cc in apoe_runs:
        r = run_one(model, fset, None, 0, cc,
                     labels, raw_splits, args, args.out_root)
        if r:
            for sd, (auc, bacc) in r.items():
                summary.append({"combo": "_global", "model": model,
                                "feature_set": fset, "seed": sd,
                                "test_auc": auc, "test_bACC": bacc})

    # ── Summary across runs ────────────────────────────────────────────────
    if summary:
        s = pd.DataFrame(summary)
        s.to_csv(args.out_root / "summary_long.csv", index=False)
        agg = (s.groupby(["combo", "model", "feature_set"])
                 .agg(test_auc=("test_auc", "mean"),
                       auc_sd=("test_auc", "std"),
                       bACC=("test_bACC", "mean"))
                 .reset_index().sort_values("test_auc", ascending=False))
        agg.to_csv(args.out_root / "summary_meanstd.csv", index=False)
        with pd.option_context("display.width", 220, "display.max_rows", None):
            print("\n=== test AUC (mean over seeds), best first ===")
            print(agg.head(30).to_string(index=False))

    print(f"\n[done] → {args.out_root}")


if __name__ == "__main__":
    main()
