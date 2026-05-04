"""
check_class_distribution.py
============================
Prints class distributions for all 4 tasks across seeds 0, 1, 2
to assess whether class-weighted CE loss is warranted.
"""

import pandas as pd
import numpy as np
from pathlib import Path

DATA_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\verbose\baseline")

TASKS = {
    "T1_binary (CN vs MCI+AD)": {
        "label_col": "Label_bl_multi",
        "label_map": {"CN": "CN (0)", "MCI": "MCI+AD (1)", "AD": "MCI+AD (1)"},
        "filter_non_ad": False,
    },
    "T2_multiclass (CN/MCI/AD)": {
        "label_col": "Label_bl_multi",
        "label_map": {"CN": "CN (0)", "MCI": "MCI (1)", "AD": "AD (2)"},
        "filter_non_ad": False,
    },
    "T3a_conv3y (non-AD baseline)": {
        "label_col": "Label_3y",
        "label_map": None,
        "filter_non_ad": True,
    },
    "T3b_conv5y (non-AD baseline)": {
        "label_col": "Label_5y",
        "label_map": None,
        "filter_non_ad": True,
    },
}

for seed in [0, 1, 2]:
    seed_dir = DATA_DIR / f"seed_{seed}"
    splits = {}
    for split in ["train", "val", "test"]:
        df = pd.read_csv(seed_dir / f"{split}.csv", low_memory=False)
        splits[split] = df

    print(f"\n{'#'*65}")
    print(f"  SEED {seed}")
    print(f"  Total — Train: {len(splits['train'])}  Val: {len(splits['val'])}  Test: {len(splits['test'])}")
    print(f"{'#'*65}")

    for task_name, cfg in TASKS.items():
        print(f"\n  [{task_name}]")
        label_col = cfg["label_col"]

        rows = []
        for split_name, df in splits.items():
            d = df.copy()
            if cfg["filter_non_ad"] and "Label_bl_multi" in d.columns:
                d = d[d["Label_bl_multi"].isin(["CN", "MCI"])]
            if cfg["label_map"]:
                d[label_col] = d[label_col].map(cfg["label_map"])
            else:
                d[label_col] = pd.to_numeric(d[label_col], errors="coerce")
            d = d.dropna(subset=[label_col])

            counts = d[label_col].value_counts().sort_index()
            n_total = len(d)
            row = {"Split": split_name, "N": n_total}
            for cls, cnt in counts.items():
                row[str(cls)] = f"{cnt} ({100*cnt/n_total:.0f}%)"
            rows.append(row)

        print(pd.DataFrame(rows).to_string(index=False))
