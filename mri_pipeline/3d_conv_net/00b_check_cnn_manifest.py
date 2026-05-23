"""
00b_check_cnn_manifest.py — fast sanity check on cnn_manifest.csv

Reads the manifest written by 00_prepare_CNN_inputs.py and verifies that the
z-scored volumes have znorm_mean ~ 0 and znorm_std ~ 1, status='ok' for all
rows, and flags offenders. No volume I/O — relies on the stats already
recorded at write time by the preprocessing step.

Healthy run expectations:
  - status histogram: all rows 'ok'
  - znorm_mean : median within +/- 0.01, max |.| < 0.05
  - znorm_std  : median within +/- 0.02 of 1.0, range [0.95, 1.05]
  - 0 offenders in both tolerance bins below

If any offender appears, look at its brain_voxels: the mean/std are computed
over voxels where |val| > 1e-6, so a tiny/sparse brain mask is the usual
cause. Healthy adult T1 ~ 1-2 M brain voxels after z-scoring; < 500 k is
suspect.

Usage (on CSD3, env: mri or any env with pandas):
  python mri_pipeline/3d_conv_net/00b_check_cnn_manifest.py
  python mri_pipeline/3d_conv_net/00b_check_cnn_manifest.py /path/to/cnn_manifest.csv
"""

from __future__ import annotations

import sys
import pandas as pd

DEFAULT_MANIFEST = (
    "/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/cnn_inputs/cnn_manifest.csv"
)
TOL_MEAN = 0.05        # tolerated |znorm_mean|
TOL_STD  = 0.05        # tolerated |znorm_std - 1|


def main(manifest_path: str = DEFAULT_MANIFEST) -> int:
    df = pd.read_csv(manifest_path)
    ok = df[df["status"] == "ok"].copy()
    ok["znorm_mean"] = pd.to_numeric(ok["znorm_mean"], errors="coerce")
    ok["znorm_std"]  = pd.to_numeric(ok["znorm_std"],  errors="coerce")
    m, s = ok["znorm_mean"], ok["znorm_std"]

    print(f"manifest            : {manifest_path}")
    print(f"rows total          : {len(df)}")
    print(f"rows ok             : {len(ok)}")
    print(f"status histogram    : {df['status'].value_counts().to_dict()}")
    if len(ok) == 0:
        print("[ERROR] no 'ok' rows in the manifest.")
        return 1

    print(f"znorm_mean  [min / median / max] : "
          f"{m.min():+.4f} / {m.median():+.4f} / {m.max():+.4f}")
    print(f"znorm_std   [min / median / max] : "
          f"{s.min():.4f} / {s.median():.4f} / {s.max():.4f}")

    bad_m = ok[m.abs() > TOL_MEAN]
    bad_s = ok[(s - 1).abs() > TOL_STD]
    print(f"rows with |mean|  > {TOL_MEAN}    : {len(bad_m)}")
    print(f"rows with |std-1| > {TOL_STD}    : {len(bad_s)}")

    for label, bad, ref in [("mean", bad_m, 0.0), ("std", bad_s, 1.0)]:
        if len(bad):
            col = "znorm_mean" if label == "mean" else "znorm_std"
            worst = bad.reindex(
                (bad[col] - ref).abs().sort_values(ascending=False).index
            ).head(5)
            print(f"  worst 5 {label} offenders:")
            for _, r in worst.iterrows():
                print(f"    {r['bids_sub']} {r['bids_ses']}: "
                      f"mean={r['znorm_mean']:+.4f}  std={r['znorm_std']:.4f}  "
                      f"brain_voxels={r.get('brain_voxels', '?')}")

    n_bad_status = int((df["status"] != "ok").sum())
    n_bad_stats  = len(bad_m) + len(bad_s)
    if n_bad_status == 0 and n_bad_stats == 0:
        print("RESULT: PASS — all rows ok, znorm stats within tolerance.")
        return 0
    print(f"RESULT: REVIEW — {n_bad_status} non-ok status, "
          f"{n_bad_stats} out-of-tolerance stat rows.")
    return 1


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_MANIFEST
    sys.exit(main(path))
