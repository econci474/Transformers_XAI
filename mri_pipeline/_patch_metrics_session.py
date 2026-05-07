"""
One-off patcher for legacy metrics.json files written by 04_supervised_finetuning_ViT.py
before the schema captured longitudinal cohort info.

Sets the `session` field to a HUMAN-READABLE cohort label (so the JSON is
self-describing in the table at a glance) and adds machine-readable
`long_mode` / `max_months` / `session_policy` fields for replay/aggregation.

The 04 sweep on CSD3 was launched with --long 1 for all 4 tasks, but:
- T1/T2 used session_policy='current' so the cohort really was bl + m12
  (~422 train scans).
- T3a/T3b used the legacy session_policy='baseline_only', which silently
  dropped non-bl scans -> effective cohort is bl-only (~21 train scans).

Heuristics (by top-level dir under outputs/ViT_B_mae75/):
  smoke_test/*  -> session='bl', long_mode=None, session_policy='current'
  T1_binary/*   -> session='bl + m12', long_mode='cutoff', max_months=12, sp='current'
  T2_multiclass -> session='bl + m12', long_mode='cutoff', max_months=12, sp='current'
  T3a_conv3y/*  -> session='bl', long_mode='cutoff', max_months=12, sp='baseline_only'
  T3b_conv5y/*  -> session='bl', long_mode='cutoff', max_months=12, sp='baseline_only'

Idempotent: re-running overwrites the same fields with the same values.
"""

import json
from pathlib import Path

OUT_ROOT = Path(__file__).resolve().parent / "outputs" / "ViT_B_mae75"

# Top-level dir -> (session_label, long_mode, max_months, session_policy)
LEGACY_RUN_PROFILE = {
    "smoke_test":    ("bl",         None,     None, "current"),
    "T1_binary":     ("bl + m12",   "cutoff", 12,   "current"),
    "T2_multiclass": ("bl + m12",   "cutoff", 12,   "current"),
    "T3a_conv3y":    ("bl",         "cutoff", 12,   "baseline_only"),
    "T3b_conv5y":    ("bl",         "cutoff", 12,   "baseline_only"),
}

n_patched = 0
for jf in OUT_ROOT.rglob("metrics.json"):
    rel = jf.relative_to(OUT_ROOT)
    top = rel.parts[0]
    if top not in LEGACY_RUN_PROFILE:
        print(f"  [skip] unrecognized layout: {rel}")
        continue
    session, long_mode, max_months, session_policy = LEGACY_RUN_PROFILE[top]

    with open(jf, "r") as f:
        data = json.load(f)
    cfg = data.get("config", {})
    cfg["session"]        = session
    cfg["long_mode"]      = long_mode
    cfg["max_months"]     = max_months
    cfg["session_policy"] = session_policy
    data["config"] = cfg

    with open(jf, "w") as f:
        json.dump(data, f, indent=2)

    print(f"  patched {rel}: "
          f"session='{session}' long_mode={long_mode} max_months={max_months} "
          f"session_policy={session_policy}")
    n_patched += 1

print(f"\nPatched {n_patched} metrics.json files under {OUT_ROOT}.")
