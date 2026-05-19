"""
23b_render_cox_table.py
=======================
Aggregate the `cox_survival/` metrics.json tree (script 23) and render, per
transition, two B&W publication tables (house style: DejaVu Serif, white
fill, thin outline):

  cox_table_<transition>.{png,md}  — discrimination/calibration
    rows = covariate set (single → pairwise → joint, then C-index desc)
    cols = C-index  R²_D  tdAUC 3/5/10y  Calib O/E 3/5/10y  n_events  n_test
           (mean ± std over seeds; **bold C-index > 0.5**)

  cox_hr_<transition>.{png,md}  — hazard ratio per 1 SD
    rows = covariate set;  cols = age, sex, APOE4, βPRS(+APOE),
           βPRS(noAPOE);  cell = HR (95% CI) per +1 SD when that covariate
           is in the model (covariates z-scored on train ⇒ exp(coef) is the
           per-SD HR), mean over seeds; **bold HR CI excluding 1**.

Calib O/E = observed (KM) / expected (model) event prob at the horizon;
≈1 well-calibrated, >1 risk under-predicted. Outputs under --root, plus
cox_summary.csv and cox_hr_summary.csv.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.family"] = "DejaVu Serif"

OUT_DEFAULT_ROOT = Path("D:/ADNI_SNP_Omni2.5M_20140220/cox_survival")
METRICS = [("c_index", "C-index"), ("r2_d", "R²_D"),
           ("tdauc_3y", "tdAUC 3y"), ("tdauc_5y", "tdAUC 5y"),
           ("tdauc_10y", "tdAUC 10y"),
           ("calib_oe_3y", "Calib O/E 3y"), ("calib_oe_5y", "Calib O/E 5y"),
           ("calib_oe_10y", "Calib O/E 10y")]
ATOMS = ["age", "sex", "APOE4", "APOE2", "betaPRS_withAPOE",
         "betaPRS_onlyAPOE", "betaPRS_noAPOE",
         "betaPRS_6W27B", "betaPRS_9W27B", "betaPRS_9W27B5D"]
ATOM_PRETTY = {"age": "age", "sex": "sex", "APOE4": "APOE4",
               "APOE2": "APOE2",
               "betaPRS_withAPOE": "βPRS_withAPOE",
               "betaPRS_onlyAPOE": "βPRS_onlyAPOE",
               "betaPRS_noAPOE": "βPRS_noAPOE",
               "betaPRS_6W27B": "βPRS_6W27B",
               "betaPRS_9W27B": "βPRS_9W27B",
               "betaPRS_9W27B5D": "βPRS_9W27B5D"}
TRANS_PRETTY = {
    "CN_to_MCI":   "CN → MCI  (sustained, sex/age/APOE/PRS at baseline)",
    "MCI_to_AD":   "MCI → AD  (progressive MCI; baseline covariates)",
    "CNMCI_to_AD": "CN + MCI → AD  (confirmed AD; baseline covariates)",
}


def _is_dosage(covset: str) -> bool:
    """A covset is a raw per-SNP dosage model iff any atom is a dosage_<S>
    block (the survival analog of the baseline raw-dosage LogReg)."""
    return any(tok.startswith("dosage_") for tok in str(covset).split("+"))


def collect(root: Path):
    rows, hrs, prs_sizes = [], [], {}
    for mf in root.rglob("metrics.json"):
        try:
            m = json.load(open(mf))
        except Exception:
            continue
        if m.get("status") != "ok":
            continue
        if not prs_sizes:
            ps = m.get("prs_set_sizes") or {}
            if ps:
                prs_sizes = dict(ps)
        rows.append({"transition": m["transition"], "covset": m["covset"],
                     "n_atoms": len(m.get("atoms", m["covset"].split("+"))),
                     "seed": m.get("seed", -1),
                     "n_events": m.get("n_events"),
                     "n_atrisk": m.get("n_atrisk"),
                     "n_test": m.get("n_test"),
                     "n_test_events": m.get("n_test_events"),
                     **{k: m.get(k) for k, _ in METRICS}})
        for atom, d in (m.get("hr_per_sd") or {}).items():
            hrs.append({"transition": m["transition"], "covset": m["covset"],
                        "n_atoms": len(m.get("atoms",
                                             m["covset"].split("+"))),
                        "seed": m.get("seed", -1), "atom": atom,
                        "hr": d.get("hr"), "lo": d.get("lo"),
                        "hi": d.get("hi"), "p": d.get("p")})
    if not rows:
        raise SystemExit(f"No usable metrics.json under {root}")
    return pd.DataFrame(rows), pd.DataFrame(hrs), prs_sizes


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    out = []
    for (tn, cs), g in df.groupby(["transition", "covset"]):
        r = {"transition": tn, "covset": cs,
             "n_atoms": int(g["n_atoms"].iloc[0]), "n_seeds": len(g),
             "n_events": (f"{g['n_events'].dropna().mean():.0f}"
                          if g["n_events"].notna().any() else ""),
             "n_test": (f"{g['n_test'].dropna().mean():.0f}"
                        if g["n_test"].notna().any() else "")}
        for k, _ in METRICS:
            v = g[k].dropna().to_numpy(dtype=float)
            r[k] = f"{v.mean():.3f}±{v.std(ddof=0):.3f}" if len(v) else ""
            r[f"{k}_mean"] = float(v.mean()) if len(v) else np.nan
        out.append(r)
    return pd.DataFrame(out)


def summarise_hr(hr: pd.DataFrame) -> pd.DataFrame:
    """Mean HR / mean 95% CI over seeds per (transition, covset, atom)."""
    out = []
    for (tn, cs, at), g in hr.groupby(["transition", "covset", "atom"]):
        v = g[["hr", "lo", "hi"]].apply(pd.to_numeric, errors="coerce").dropna()
        if v.empty:
            continue
        h, lo, hi = v["hr"].mean(), v["lo"].mean(), v["hi"].mean()
        out.append({"transition": tn, "covset": cs,
                    "n_atoms": int(g["n_atoms"].iloc[0]), "atom": at,
                    "cell": f"{h:.2f} ({lo:.2f}–{hi:.2f})",
                    "hr_mean": float(h), "excl1": bool(lo > 1 or hi < 1)})
    return pd.DataFrame(out)


def _order(df: pd.DataFrame) -> pd.DataFrame:
    key = "c_index_mean" if "c_index_mean" in df.columns else "covset"
    asc = [True, False] if key == "c_index_mean" else [True, True]
    return df.sort_values(["n_atoms", key], ascending=asc)


def _grid(df, cols, out_png, title, fig_w, bold_fn=None, subtitle=None,
          heavy_break_col=None):
    """Generic B&W table. cols = [(key,label,width|None,align)]; bold_fn(row,
    key)→bool selects bold cells; None-width cols share the remaining width.
    `subtitle` = list of centred lines drawn under the title. When
    `heavy_break_col` is given, a thick rule is drawn wherever that column's
    value changes between consecutive rows (section separator)."""
    subtitle = subtitle or []
    n = len(df)
    fixed = sum(c[2] for c in cols if c[2] is not None)
    nfree = sum(1 for c in cols if c[2] is None)
    free = (1.0 - fixed) / nfree if nfree else 0.0
    widths = [c[2] if c[2] is not None else free for c in cols]
    xs = [0.0]
    for w in widths:
        xs.append(xs[-1] + w)
    row_h, hdr_h, title_h, margin = 0.34, 0.5, 0.5, 0.22
    sub_h = 0.30 * len(subtitle)
    fig_h = title_h + sub_h + hdr_h + n * row_h + 2 * margin
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fh = fig_h

    def Y(v): return v / fh

    def box(x, y, w, h, b=0.4):
        ax.add_patch(plt.Rectangle((x, Y(y)), w, h / fh,
                                   transform=ax.transAxes, facecolor="white",
                                   edgecolor="black", linewidth=b))

    def txt(x, y, s, fs, ha, wt="normal"):
        ax.text(x, Y(y), s, ha=ha, va="center", fontsize=fs, fontweight=wt,
                transform=ax.transAxes)

    cur = fh - margin
    txt(0.5, cur - title_h * 0.5, title, 12, "center")
    cur -= title_h
    for sline in subtitle:
        txt(0.5, cur - 0.15, sline, 8.6, "center")
        cur -= 0.30
    for i, (_, lbl, _, _) in enumerate(cols):
        box(xs[i], cur - hdr_h, xs[i + 1] - xs[i], hdr_h, b=0.7)
        txt((xs[i] + xs[i + 1]) / 2, cur - hdr_h / 2, lbl, 8.2, "center")
    cur -= hdr_h
    prev = prev_blk = None
    for _, r in df.iterrows():
        blk = r.get(heavy_break_col) if heavy_break_col else None
        if heavy_break_col is not None and prev_blk is not None \
                and blk != prev_blk:
            ax.plot([0, 1], [Y(cur), Y(cur)], transform=ax.transAxes,
                    linewidth=2.4, color="black")          # section break
        elif prev is not None and int(r["n_atoms"]) != prev:
            ax.plot([0, 1], [Y(cur), Y(cur)], transform=ax.transAxes,
                    linewidth=1.1, color="black")
        prev = int(r["n_atoms"]); prev_blk = blk
        for i, (key, _, _, al) in enumerate(cols):
            box(xs[i], cur - row_h, xs[i + 1] - xs[i], row_h)
            v = str(r.get(key, ""))
            wt = "bold" if (bold_fn and bold_fn(r, key)) else "normal"
            if al == "left":
                txt(xs[i] + 0.005, cur - row_h / 2, v, 7.0, "left", wt)
            else:
                txt((xs[i] + xs[i + 1]) / 2, cur - row_h / 2, v, 7.0,
                    "center", wt)
        cur -= row_h
    fig.savefig(out_png, dpi=150, facecolor="white", bbox_inches="tight",
                pad_inches=0.05)
    plt.close(fig)
    print(f"Wrote {out_png}  ({n} rows)")


def render_metrics(df, root, tn, title, subtitle):
    df = df.copy()
    df["block"] = df["covset"].map(lambda c: 1 if _is_dosage(c) else 0)
    df = df.sort_values(["block", "n_atoms", "c_index_mean"],
                        ascending=[True, True, False]).reset_index(drop=True)
    cols = [("covset", "Covariate set", 0.26, "left")]
    cols += [(m[0], m[1], None, "center") for m in METRICS]

    def bold(r, k):
        return (k == "c_index" and pd.notna(r.get("c_index_mean"))
                and float(r["c_index_mean"]) > 0.5)

    _grid(df, cols, root / f"cox_table_{tn}.png", title, 22.0, bold,
          subtitle=subtitle, heavy_break_col="block")
    hdr = ["Covariate set"] + [m[1] for m in METRICS]
    lines = [f"# {title}", ""]
    lines += [f"*{s}*" for s in subtitle] + [""]
    lines += ["Mean ± std over seeds 0/1/2; fit on TRAIN, evaluated on the "
              "held-out VALIDATION fold. The Cox model is **per-subject** "
              "(baseline covariates → time-to-AD): at-risk / events / "
              "censored count SUBJECTS (one baseline record each), NOT "
              "subject-visits. **Bold C-index > 0.5**. R²_D = "
              "Royston–Sauerbrei explained variation; tdAUC = Uno IPCW "
              "time-dependent AUC; Calib O/E = observed(KM)/expected(model) "
              "event prob at the horizon (≈1 = well-calibrated, >1 = risk "
              "under-predicted, blank = horizon beyond val follow-up).", "",
              "**`betaPRS_<set>` vs `dosage_<set>` (the two rows below the "
              "thick rule are the dosage models).** Both are Cox PH models "
              "on the SAME SNP set; they differ only in how those SNPs enter "
              "the model:", "",
              "- **`betaPRS_<set>`** — the SNPs are collapsed into **one** "
              "covariate, the GWAS-β-weighted polygenic score "
              "Σᵢ βᵢ·dosageᵢ using the *external published* effect sizes "
              "(βᵢ fixed a-priori, not estimated here), z-scored. "
              "Parsimonious (1 df).",
              "- **`dosage_<set>`** — the *same* SNPs entered as **raw "
              "per-SNP allele-dosage covariates** (each SNP its own "
              "covariate, 0/1/2, z-scored); the Cox fit estimates every "
              "SNP's coefficient **de-novo from this ADNI fold** instead of "
              "imposing the external β (ridge-penalised; SNPs monomorphic in "
              "the train fold dropped, perfect-LD duplicates collapsed). It "
              "does **not** use the GWAS betas — it asks whether the cohort "
              "can *learn* the weighting itself. Many more parameters "
              "(≈SNP-count df) ⇒ higher variance / overfitting on the small "
              "fold, so it typically underperforms the 1-df β-PRS — the "
              "survival analog of the script-20 raw-dosage LogReg.",
              "", "PRS tiers (unfiltered family only): `betaPRS_withAPOE` = "
              "all 430; `betaPRS_onlyAPOE` = 166 isolated chr19 APOE LD "
              "block; `betaPRS_noAPOE` = 263 APOE region removed. APOE4 / "
              "APOE2 = clinical ε4 / ε2 allele dosage.", "",
              "| " + " | ".join(hdr) + " |",
              "|" + "|".join(["---"] * len(hdr)) + "|"]
    for _, r in df.iterrows():
        cb = bold(r, "c_index")
        cells = [str(r["covset"])]
        for k, _ in METRICS:
            val = str(r[k])
            cells.append(f"**{val}**" if (k == "c_index" and cb and val)
                         else val)
        lines.append("| " + " | ".join(cells) + " |")
    (root / f"cox_table_{tn}.md").write_text("\n".join(lines),
                                             encoding="utf-8")
    print(f"Wrote {root/f'cox_table_{tn}.md'}")


def render_hr(hsum, root, tn, title, subtitle):
    sub = hsum[hsum["transition"] == tn]
    if sub.empty:
        return
    wide = sub.pivot_table(index=["covset", "n_atoms"], columns="atom",
                           values="cell", aggfunc="first")
    bold = sub.pivot_table(index=["covset", "n_atoms"], columns="atom",
                           values="excl1", aggfunc="first")
    wide = wide.reset_index()
    bold = bold.reset_index()
    wide["block"] = wide["covset"].map(lambda c: 1 if _is_dosage(c) else 0)
    wide = wide.sort_values(["block", "n_atoms", "covset"]).reset_index(
        drop=True)
    bset = {(rw["covset"], a): bool(rw.get(a) is True)
            for _, rw in bold.iterrows() for a in ATOMS if a in bold.columns}
    cols = [("covset", "Covariate set", 0.30, "left")]
    cols += [(a, ATOM_PRETTY[a], None, "center") for a in ATOMS]
    for a in ATOMS:
        if a not in wide.columns:
            wide[a] = ""
    wide[ATOMS] = wide[ATOMS].fillna("")

    def bfn(r, k):
        return k in ATOMS and bset.get((r["covset"], k), False)

    _grid(wide, cols, root / f"cox_hr_{tn}.png", title, 19.0, bfn,
          subtitle=subtitle, heavy_break_col="block")
    hdr = ["Covariate set"] + [ATOM_PRETTY[a] for a in ATOMS]
    lines = [f"# {title}", ""]
    lines += [f"*{s}*" for s in subtitle] + [""]
    lines += ["Hazard ratio per **+1 SD** of the covariate (z-scored on "
              "train), mean over seeds 0/1/2, (mean 95% CI). Blank = "
              "covariate not in that model. **Bold = CI excludes 1**. Rows "
              "below the thick rule are the raw per-SNP `dosage_<set>` "
              "models; a `dosage_<set>` block is a multivariable signature "
              "with no single per-SD HR, so only its accompanying clinical "
              "covariates (age/sex/APOE2/APOE4) get an HR column here — see "
              "the discrimination table for those models' C-index. PRS "
              "tiers (unfiltered family only): `βPRS_withAPOE` = all 430; "
              "`βPRS_onlyAPOE` = 166 isolated chr19 APOE LD block; "
              "`βPRS_noAPOE` = 263 APOE region removed. APOE4/APOE2 = "
              "clinical ε4/ε2 dosage.", "",
              "| " + " | ".join(hdr) + " |",
              "|" + "|".join(["---"] * len(hdr)) + "|"]
    for _, r in wide.iterrows():
        cells = [str(r["covset"])]
        for a in ATOMS:
            val = str(r.get(a, "") or "")
            cells.append(f"**{val}**" if (val and bset.get((r["covset"], a)))
                         else val)
        lines.append("| " + " | ".join(cells) + " |")
    (root / f"cox_hr_{tn}.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {root/f'cox_hr_{tn}.md'}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=OUT_DEFAULT_ROOT)
    args = ap.parse_args()

    long_df, hr_df, prs_sizes = collect(args.root)
    summ = summarise(long_df)
    summ.to_csv(args.root / "cox_summary.csv", index=False)
    print(f"Wrote {args.root/'cox_summary.csv'} ({len(summ)} rows)")
    hsum = summarise_hr(hr_df) if not hr_df.empty else pd.DataFrame()
    if not hsum.empty:
        hsum.to_csv(args.root / "cox_hr_summary.csv", index=False)
        print(f"Wrote {args.root/'cox_hr_summary.csv'} ({len(hsum)} rows)")

    def _tsub(tn):
        g = long_df[long_df["transition"] == tn]

        def _f(col):
            v = pd.to_numeric(g[col], errors="coerce").dropna()
            return v.iloc[0] if len(v) else None

        def _m(col):
            v = pd.to_numeric(g[col], errors="coerce").dropna()
            return v.mean() if len(v) else None
        na, ne, nt = _f("n_atrisk"), _f("n_events"), _m("n_test")
        if na is None or ne is None:
            return []
        return [f"N(val) ≈ {nt:.0f} subjects     ·     at-risk = {int(na)} "
                f"subjects     ·     events = {int(ne)}     ·     censored "
                f"= {int(na - ne)}"]

    def _setlines():
        """SNP counts per PRS set + a one-line betaPRS-vs-dosage key,
        shown only for the filtered families (prs_set_sizes recorded)."""
        if not prs_sizes:
            return []
        order = ["6W27B", "9W27B", "9W27B5D"]
        seq = [s for s in order if s in prs_sizes] + \
              [s for s in prs_sizes if s not in order]
        sizes = "  ·  ".join(f"{s} = {prs_sizes[s]} SNPs" for s in seq)
        return [
            f"PRS sets (noAPOE-filtered, GWAS-β-weighted):   {sizes}",
            "betaPRS_<set> = ONE covariate = Σ(published GWAS β × dosage), "
            "z-scored.   dosage_<set> = the SAME SNPs as raw per-SNP "
            "dosage covariates (β learned in-cohort, ridge; mono/LD-dup "
            "pruned) — these are the rows in the lower section."]

    setlines = _setlines()
    for tn in summ["transition"].unique():
        sub = summ[summ["transition"] == tn].copy()
        tp = TRANS_PRETTY.get(tn, tn)
        st = _tsub(tn) + setlines
        render_metrics(sub, args.root, tn,
                       f"Cox survival — {tp}  (ADNI; VALIDATION fold; "
                       f"mean ± std, seeds 0/1/2; bold C-index > 0.5)", st)
        if not hsum.empty:
            render_hr(hsum, args.root, tn,
                      f"Cox HR per +1 SD — {tp}  (ADNI; VALIDATION fold; "
                      f"mean over seeds 0/1/2; bold CI excludes 1)", st)
        show = _order(sub)[["covset", "c_index", "tdauc_5y", "tdauc_10y",
                            "calib_oe_5y"]]
        with pd.option_context("display.width", 220, "display.max_rows", None):
            print(f"\n=== {tn} (size order; C-index desc) ===")
            print(show.to_string(index=False))


if __name__ == "__main__":
    main()
