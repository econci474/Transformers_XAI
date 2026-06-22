r"""
02m_soden_xsec.py   (env: snp [torch-cpu + torchdiffeq] locally, or clinical/torch on GPU)
=========================================================================================
Cross-sectional SODEN baseline: train the neural-ODE SODEN survival head DIRECTLY on the 8-feature
cross-sectional clinical baseline vector (the same features RSF / Cox PH / Weibull AFT / Weibull-piecewise
use), per seed. This is the cross-sectional counterpart of the deep "SODEN / MAMBA (frozen)" head — SODEN
has no classical tabular analogue (it is a neural ODE), so it must be trained rather than fit in closed form.

Reuses the SODENHead from _mamba_seq_lib (z = the 8 standardised features, d_out=8; no MAMBA encoder) and
exports, per seed, a per-patient predicted survival curve on the SAME dense grid as 01 — so the existing
02_survival_comparison.py auto-discovers it (PRED_ROOT/<config>/ scan) and 08_wald_significance.py loads it
through the MAMBA_CFGS parquet path. Downstream we place it in the CROSS-SECTIONAL baseline block.

The 8-feature build / median-impute / standardise (fit on each seed's TRAIN) mirrors
clinical_pipeline/02l_survival_model_comparison.py verbatim (02l can't be imported here — it pulls
sksurv/lifelines at import; snp has neither).

Outputs: PRED_ROOT/soden_xsec/seed_{s}/predictions.parquet (Patient_ID,split,event,time,S_000..S_060)
         + meta.json (grid_times). NO weights kept.

Run (local CPU):
  conda run -n snp python integration/T4/mamba_long_clinical/02m_soden_xsec.py
Then locally: re-run 02 (master) -> 08 (wald) -> 07 (tables).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _mamba_seq_lib as lib

SPLIT_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion"
                 r"\tabular\baseline")
PRED_ROOT = Path(r"D:\ADNI_BIDS_project\derivatives\mamba_survival\survival")
CONFIG = "soden_xsec"
SEEDS = (0, 1, 2)
EPS = 0.01
COMPACT_FEATURES = ["MMSE_Total", "ADAS_Total13", "MoCA_Total", "FAQ_Total",
                    "age_at_baseline", "sex_M", "Education_Years", "APOE4_Dosage"]
GRID_TIMES = np.round(np.linspace(0.0, 15.0, 61), 4)        # dense grid for S export (matches 01)


# ── 8-feature build + per-seed TRAIN-fitted impute/scale (verbatim from 02l) ─────────
def _build(df, features):
    df = df[df["bl_dx"].isin(["CN", "MCI"])].copy()
    ev = df["conversion_group"].isin(["pCN_to_AD", "pMCI"]).astype(bool).values
    t = np.where(ev, pd.to_numeric(df["years_to_AD"], errors="coerce"),
                 pd.to_numeric(df["FU_years"], errors="coerce"))
    t = pd.Series(t).fillna(EPS).clip(lower=EPS).values
    df["Date"] = pd.to_datetime(df["Date"], errors="coerce")
    birth = pd.to_datetime(pd.DataFrame(
        {"year": df["YOB"].astype(int), "month": 7, "day": 1}).set_index(df.index))
    df["age_at_baseline"] = (df["Date"] - birth) / np.timedelta64(1, "D") / 365.25
    df["sex_M"] = (df["Sex"].astype(str) == "Male").astype(float)
    X = df[features].apply(pd.to_numeric, errors="coerce")
    pid = df["Patient_ID"].astype(str).to_numpy()
    return X, ev, t, pid


def load_seed(seed, split_dir):
    parts = {sp: pd.read_csv(split_dir / f"seed_{seed}" / f"{sp}.csv", low_memory=False)
             for sp in ("train", "val", "test")}
    Xtr, _, _, _ = _build(parts["train"], COMPACT_FEATURES)
    imp = SimpleImputer(strategy="median").fit(Xtr)
    scl = StandardScaler().fit(imp.transform(Xtr))
    out = {}
    for sp in ("train", "val", "test"):
        X, e, t, pid = _build(parts[sp], COMPACT_FEATURES)
        Xs = scl.transform(imp.transform(X)).astype(np.float32)
        out[sp] = (Xs, e.astype(np.float32), t.astype(np.float32), pid)
    return out


# ── train the SODEN head on the cross-sectional vector; early-stop on val NLL ────────
def train_soden(dd, device, lr=1e-2, max_epochs=3000, patience=150):
    import torch
    head = lib.make_head("soden", len(COMPACT_FEATURES)).to(device)

    def tens(sp):
        X, e, t, _ = dd[sp]
        return (torch.tensor(X, device=device),
                torch.tensor(t, device=device), torch.tensor(e, device=device))

    ztr, ttr, etr = tens("train"); zva, tva, eva = tens("val")
    opt = torch.optim.Adam(head.parameters(), lr=lr, weight_decay=1e-4)
    best, best_state, bad = float("inf"), None, 0
    for ep in range(max_epochs):
        head.train(); opt.zero_grad()
        loss = head.nll(ztr, ttr, etr)
        if not torch.isfinite(loss):
            break
        loss.backward()
        torch.nn.utils.clip_grad_norm_(head.parameters(), 5.0)
        opt.step()
        head.eval()
        with torch.no_grad():
            vnll = float(head.nll(zva, tva, eva))
        if vnll < best - 1e-5:
            best, best_state, bad = vnll, {k: v.clone() for k, v in head.state_dict().items()}, 0
        else:
            bad += 1
            if bad >= patience:
                break
    if best_state is not None:
        head.load_state_dict(best_state)
    return head, best


def run_seed(seed, split_dir, out_root, device):
    import torch
    lib.seed_everything(seed)
    dd = load_seed(seed, split_dir)
    head, val_nll = train_soden(dd, device)

    grid = torch.tensor(GRID_TIMES.astype(np.float32), device=device)
    frames = []
    head.eval()
    with torch.no_grad():
        for sp in ("train", "val", "test"):
            X, e, t, pid = dd[sp]
            z = torch.tensor(X, device=device)
            S = head.survival(z, grid).cpu().numpy()                # [N, 61]
            df = pd.DataFrame({"Patient_ID": pid, "split": sp,
                               "event": e.astype(int), "time": t.astype(float)})
            for j in range(len(GRID_TIMES)):
                df[f"S_{j:03d}"] = S[:, j]
            frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    od = out_root / CONFIG / f"seed_{seed}"
    od.mkdir(parents=True, exist_ok=True)
    out.to_parquet(od / "predictions.parquet", index=False)
    json.dump({"config": CONFIG, "seed": seed, "head": "soden", "lr": 1e-2,
               "val_nll": float(val_nll), "grid_times": GRID_TIMES.tolist(), "n": int(len(out))},
              open(od / "meta.json", "w"), indent=2)
    n_ev = int(out["event"].sum())
    print(f"  [soden_xsec seed{seed}] val_nll={val_nll:.4f}  N={len(out)} (events={n_ev}) -> {od}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS))
    ap.add_argument("--split_dir", type=str, default=str(SPLIT_DIR),
                    help="dir with seed_*/{train,val,test}.csv (override on HPC)")
    ap.add_argument("--out_dir", type=str, default=str(PRED_ROOT),
                    help="PRED_ROOT that 02/08 scan; soden_xsec/<seed>/ is written under it")
    ap.add_argument("--device", type=str, default=None)
    args = ap.parse_args()
    import torch
    device = args.device or ("cuda" if torch.cuda.is_available() else "cpu")
    split_dir, out_root = Path(args.split_dir), Path(args.out_dir)
    print(f"device={device}  split_dir={split_dir}  out_root={out_root}")
    for seed in args.seeds:
        run_seed(seed, split_dir, out_root, device)
    print("Done. Now re-run 02 (master) -> 08 (wald) -> 07 (tables).")


if __name__ == "__main__":
    main()
