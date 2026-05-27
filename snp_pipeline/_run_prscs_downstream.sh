#!/usr/bin/env bash
# Run the 3 downstream tasks ONCE with --beta-source prscs against the
# unpruned 115-SNP strict-QC pool. PRS-CS already does LD-aware shrinkage
# internally, so we do NOT loop over LD configs (that would mean applying
# PLINK's pairwise-r² LD pruning on top of PRS-CS posteriors, which discards
# information PRS-CS deliberately preserved).
#
# Then render the single-block PRS-CS master leaderboard.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PY="/c/Users/elena/miniforge3/envs/snp/python.exe"

for task in 45_prs_classification_strictQC 46_prs_cox_strictQC 47_prs_aao_regression_strictQC; do
  echo "=== $task (prscs, no LD pruning) ==="
  "$PY" "$HERE/${task}.py" --beta-source prscs 2>&1 | tail -2
done

echo "=== rendering PRS-CS master leaderboard ==="
"$PY" "$HERE/48_render_strictQC_master_leaderboard.py" --beta-source prscs 2>&1 | tail -5
