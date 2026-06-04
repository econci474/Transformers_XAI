"""
02_failure_mode_T1e.py
======================
Per-modality FAILURE-MODE contingency for the T1e late-fusion arm. Reads the
8 (CL, SNP) combos, re-derives per-patient clinical and SNP predictions on
TEST (both modalities present), and tabulates:

  4x4 contingency       clin_outcome  in  {TP, FP, FN, TN}
                      x  snp_outcome   in  {TP, FP, FN, TN}

aggregated over the 3 seeds. The diagonal ({TP,TP}, {TN,TN}) is the both-agree
mass; off-diagonal cells are complementary (one modality catches what the
other misses), which is what motivates fusion. We also flag the two
*complementary-correct* cells -- (clin=TP & snp=FN) and (clin=FN & snp=TP) --
because those are the cases where late fusion can plausibly recover a hit
that one modality alone would miss.

REUSE, don't duplicate: imports `01_fuse_T1e.py` as a module (load_clinical,
load_snp, build_frame, CLIN_VARIANTS, SNP_VARIANTS, SEEDS_DEFAULT, CLASS_NAMES)
via importlib.util.spec_from_file_location. No edits to 01.

Inputs (auto-resolved):
  outputs/baseline/coverage.csv        (just for sanity / cross-check)
  on-disk CL + SNP parquets (same as 01)

Outputs (outputs/baseline/):
  failure_mode_contingency.csv         long-form 4x4 cells per combo
  failure_mode_heatmap.png             2x4 panel of 4x4 heatmaps

Usage:
  python integration/T1e_classification/02_failure_mode_T1e.py
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Import the 01_fuse_T1e.py harness as a module -- reuse, don't edit.
# --------------------------------------------------------------------------- #
_HARNESS = Path(__file__).resolve().parent / "01_fuse_T1e.py"
_spec = importlib.util.spec_from_file_location("fuse_t1e_harness", _HARNESS)
h = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(h)

CLIN_VARIANTS = h.CLIN_VARIANTS
SNP_VARIANTS = h.SNP_VARIANTS
SEEDS_DEFAULT = h.SEEDS_DEFAULT
CLASS_NAMES = h.CLASS_NAMES        # ["sCN", "pCN"]
OUT_DIR_DEFAULT = h.OUT_DIR_DEFAULT

OUTCOMES = ["TP", "FN", "FP", "TN"]   # row/col order: positive class first


# --------------------------------------------------------------------------- #
# Outcome classification (binary: pCN=positive, sCN=negative)
# --------------------------------------------------------------------------- #
def _outcome(y: int, pred: int) -> str:
    if y == 1 and pred == 1: return "TP"
    if y == 1 and pred == 0: return "FN"
    if y == 0 and pred == 1: return "FP"
    return "TN"


def collect_contingencies(seeds=SEEDS_DEFAULT):
    """For each (clin_variant, snp_variant), build the aggregated 4x4 contingency
    of (clin_outcome, snp_outcome) counts over TEST patients (both modalities
    present) summed across the requested seeds.

    Returns:
      contingency: dict[(clin_key, snp_key)] -> 4x4 np.ndarray (int counts)
      per_patient: list of dicts (one row per (combo, seed, test patient))
    """
    contingency = {}
    per_patient = []
    for clin_v in CLIN_VARIANTS:
        for snp_v in SNP_VARIANTS:
            grid = np.zeros((len(OUTCOMES), len(OUTCOMES)), dtype=int)
            for seed in seeds:
                try:
                    clin = h.load_clinical(seed, clin_v["model"], clin_v["strat"])
                    snp = h.load_snp(seed, snp_v["parquet_tmpl"])
                except FileNotFoundError as e:
                    print(f"  [WARN] {e}; skipping (cl='{clin_v['key']}', "
                          f"snp='{snp_v['key']}', seed={seed}).")
                    continue
                fr = h.build_frame(clin, snp, fill_missing_snp=False)
                test_b = fr[(fr["split"] == "test")
                            & fr["clin_present"] & fr["snp_present"]]
                if len(test_b) == 0:
                    continue
                cp = test_b[["cp0", "cp1"]].to_numpy(float).argmax(1)
                sp = test_b[["sp0", "sp1"]].to_numpy(float).argmax(1)
                yt = test_b["y_true"].to_numpy(int)
                for i, (_, row) in enumerate(test_b.reset_index(drop=True).iterrows()):
                    co = _outcome(int(yt[i]), int(cp[i]))
                    so = _outcome(int(yt[i]), int(sp[i]))
                    grid[OUTCOMES.index(co), OUTCOMES.index(so)] += 1
                    per_patient.append({
                        "clin_variant": clin_v["key"],
                        "snp_variant":  snp_v["key"],
                        "seed": seed,
                        "Patient_ID": row["Patient_ID"],
                        "y_true": CLASS_NAMES[int(yt[i])],
                        "clin_pred": CLASS_NAMES[int(cp[i])],
                        "snp_pred":  CLASS_NAMES[int(sp[i])],
                        "clin_outcome": co,
                        "snp_outcome":  so,
                    })
            contingency[(clin_v["key"], snp_v["key"])] = grid
    return contingency, per_patient


# --------------------------------------------------------------------------- #
# Rendering -- 2x4 panel of 4x4 heatmaps
# --------------------------------------------------------------------------- #
def render_panel_png(contingency: dict, out_path: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARN] matplotlib not installed -- skipping PNG.")
        return

    cl_short = {v["key"]: v["short"] for v in CLIN_VARIANTS}
    sn_short = {v["key"]: v["short"] for v in SNP_VARIANTS}
    n_cl, n_sn = len(CLIN_VARIANTS), len(SNP_VARIANTS)

    fig, axes = plt.subplots(n_cl, n_sn, figsize=(3.2 * n_sn, 3.9 * n_cl),
                              squeeze=False)
    cmap = plt.cm.Blues
    # Shared vmax across panels for visual comparability
    vmax = max((g.max() for g in contingency.values() if g.size), default=1)

    for r, clin_v in enumerate(CLIN_VARIANTS):
        for c, snp_v in enumerate(SNP_VARIANTS):
            ax = axes[r][c]
            grid = contingency.get((clin_v["key"], snp_v["key"]),
                                    np.zeros((4, 4), dtype=int))
            ax.imshow(grid, cmap=cmap, vmin=0, vmax=max(vmax, 1))
            ax.set_xticks(range(len(OUTCOMES)))
            ax.set_yticks(range(len(OUTCOMES)))
            ax.set_xticklabels(OUTCOMES, fontsize=9)
            ax.set_yticklabels(OUTCOMES, fontsize=9)
            ax.set_xlabel("SNP outcome", fontsize=9)
            if c == 0:
                ax.set_ylabel(f"{cl_short[clin_v['key']]}\nCL outcome",
                               fontsize=9)
            # Title every panel so each row's SNP variant is self-labelled
            # (previously only row 0 had titles -- row 1 then overlapped the
            # bottom of row 0's heatmap content).
            ax.set_title(sn_short[snp_v["key"]], fontsize=10)
            # cell counts -- bold when complementary-correct
            for i in range(len(OUTCOMES)):
                for j in range(len(OUTCOMES)):
                    val = int(grid[i, j])
                    is_complement = ((OUTCOMES[i] == "TP" and OUTCOMES[j] == "FN")
                                     or (OUTCOMES[i] == "FN" and OUTCOMES[j] == "TP")
                                     or (OUTCOMES[i] == "TN" and OUTCOMES[j] == "FP")
                                     or (OUTCOMES[i] == "FP" and OUTCOMES[j] == "TN"))
                    weight = "bold" if is_complement and val > 0 else "normal"
                    # contrast: dark text on light cells, white text on dark cells
                    color = "white" if grid[i, j] > 0.55 * vmax else "black"
                    ax.text(j, i, str(val), ha="center", va="center",
                             fontsize=10, color=color, weight=weight)
            ax.tick_params(length=0)

    fig.suptitle("T1e late fusion -- per-modality failure-mode contingency on TEST "
                  "(counts summed over seeds 0,1,2)",
                  fontsize=12, y=1.005)
    foot_lines = [
        "Each panel: rows = CL prediction outcome, cols = SNP prediction outcome.",
        "Outcomes (binary, pCN = positive): TP = correct positive (predict pCN, is pCN);",
        "  FN = missed positive (predict sCN, is pCN); FP = false alarm (predict pCN, is sCN);",
        "  TN = correct negative (predict sCN, is sCN).",
        "BOLD cells = COMPLEMENTARY-CORRECT (one modality right, the other wrong).",
        "  -- TP/FN and FN/TP for the positive class",
        "  -- TN/FP and FP/TN for the negative class",
        "  These are the cases where late-fusion can plausibly recover by averaging.",
        "Diagonal cells (TP/TP, FN/FN, FP/FP, TN/TN) = both modalities agree.",
    ]
    fig.text(0.01, -0.02, "\n".join(foot_lines), ha="left", va="top",
              fontsize=7.5, family="monospace", linespacing=1.4)
    fig.tight_layout(rect=[0, 0.04, 1, 0.99])
    fig.subplots_adjust(hspace=0.55, wspace=0.30)
    fig.savefig(out_path, dpi=180, bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)
    print(f"  PNG  -> {out_path}")


def contingency_to_df(contingency: dict) -> pd.DataFrame:
    rows = []
    for (cl_key, sn_key), grid in contingency.items():
        for i, co in enumerate(OUTCOMES):
            for j, so in enumerate(OUTCOMES):
                rows.append({
                    "clin_variant": cl_key,
                    "snp_variant":  sn_key,
                    "clin_outcome": co,
                    "snp_outcome":  so,
                    "n_patients":   int(grid[i, j]),
                })
    return pd.DataFrame(rows)


def print_complementarity(contingency: dict) -> None:
    print("\n[complementarity per combo]")
    print("  off_diag_correct = n(CL=TP & SNP=FN) + n(CL=FN & SNP=TP) "
          "+ n(CL=TN & SNP=FP) + n(CL=FP & SNP=TN)")
    print("  both_correct     = n(CL=TP & SNP=TP) + n(CL=TN & SNP=TN)")
    print("  both_wrong       = n(CL=FN & SNP=FN) + n(CL=FP & SNP=FP) "
           "+ disagreeing-wrongs (CL=FN & SNP=FP, etc.)")
    idx = {o: i for i, o in enumerate(OUTCOMES)}
    for (cl_key, sn_key), grid in contingency.items():
        total = int(grid.sum())
        if total == 0:
            continue
        both_correct = int(grid[idx["TP"], idx["TP"]] + grid[idx["TN"], idx["TN"]])
        comp_correct = int(grid[idx["TP"], idx["FN"]] + grid[idx["FN"], idx["TP"]]
                            + grid[idx["TN"], idx["FP"]] + grid[idx["FP"], idx["TN"]])
        both_wrong = total - both_correct - comp_correct
        pct = (comp_correct / total) * 100 if total else 0
        print(f"  {cl_key:>26s} x {sn_key:<22s} "
              f"n_test={total:3d}  both_correct={both_correct:3d}  "
              f"complement={comp_correct:3d} ({pct:4.1f}%)  "
              f"both_wrong={both_wrong:3d}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS_DEFAULT))
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR_DEFAULT)
    args = ap.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    contingency, per_patient = collect_contingencies(seeds=args.seeds)

    df_cells = contingency_to_df(contingency)
    df_cells.to_csv(out_dir / "failure_mode_contingency.csv", index=False)
    df_pp = pd.DataFrame(per_patient)
    df_pp.to_csv(out_dir / "failure_mode_per_patient.csv", index=False)
    print(f"  wrote -> {out_dir / 'failure_mode_contingency.csv'} "
           f"({len(df_cells)} cells)")
    print(f"  wrote -> {out_dir / 'failure_mode_per_patient.csv'} "
           f"({len(df_pp)} rows)")

    render_panel_png(contingency, out_dir / "failure_mode_heatmap.png")
    print_complementarity(contingency)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
