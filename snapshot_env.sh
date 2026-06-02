#!/usr/bin/env bash
# =============================================================================
# snapshot_env.sh — version-controlled, dated snapshot of a conda environment
# =============================================================================
# Writes the EXACT installed package versions of a conda env to a dated file
# under env_logs/<env>/, so we can (a) diff against a known-good state and
# (b) reproduce / roll back. This exists because a SILENT torch 2.3.1 -> 2.4.1
# drift in the `clinical` env (it diverged from clinical_pipeline/environment.yml,
# which pins 2.3.1) broke ModernBERT full fine-tune with step-1 NaNs and cost
# real compute to diagnose. A dated log makes such drift visible and reversible.
#
# Usage (run on the machine whose env you want to record):
#   conda activate clinical && bash snapshot_env.sh        # infers the active env
#   bash snapshot_env.sh clinical                          # or name it explicitly
#
# Then commit from wherever you ran it (HPC or local) — see env_logs/README.md:
#   git add env_logs/ && git commit -m "env(clinical): snapshot $(date +%F)" && git push
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- resolve which env to snapshot -------------------------------------------
ENV_NAME="${1:-${CONDA_DEFAULT_ENV:-}}"
if [[ -z "${ENV_NAME}" ]] || { [[ "${ENV_NAME}" == "base" ]] && [[ -z "${1:-}" ]]; }; then
  echo "ERROR: no env specified and no non-base env active." >&2
  echo "  Activate it first  (conda activate <env>)  or pass it  (bash snapshot_env.sh <env>)." >&2
  exit 1
fi

# Introspect inside the named env via `conda run` (works even from another env).
RUN=(conda run -n "${ENV_NAME}")

DATE="$(date +%Y-%m-%d)"
HOST="$(hostname -s 2>/dev/null || echo unknown)"

OUTDIR="${SCRIPT_DIR}/env_logs/${ENV_NAME}"
mkdir -p "${OUTDIR}"
OUT="${OUTDIR}/${ENV_NAME}_${DATE}_${HOST}.txt"

{
  echo "# ============================================================"
  echo "# conda env : ${ENV_NAME}"
  echo "# date      : ${DATE}"
  echo "# host      : ${HOST}"
  echo "# python    : $("${RUN[@]}" python --version 2>&1)"
  echo "# ============================================================"
  echo "#"
  echo "# --- key versions ---"
  "${RUN[@]}" python - <<'PY' 2>/dev/null || true
import importlib
# (display name, import name)
mods = [
    ("torch", "torch"), ("torchvision", "torchvision"),
    ("transformers", "transformers"), ("tokenizers", "tokenizers"),
    ("huggingface_hub", "huggingface_hub"), ("accelerate", "accelerate"),
    ("flash-attn", "flash_attn"), ("numpy", "numpy"), ("scipy", "scipy"),
    ("pandas", "pandas"), ("scikit-learn", "sklearn"),
    ("lightning", "lightning"), ("pytorch-lightning", "pytorch_lightning"),
    ("lifelines", "lifelines"), ("scikit-survival", "sksurv"),
    ("xgboost", "xgboost"),
]
for disp, imp in mods:
    try:
        m = importlib.import_module(imp)
        print(f"{disp:20s} {getattr(m, '__version__', '?')}")
    except Exception:
        pass
PY
  echo "#"
  echo "# --- pip freeze ---"
  "${RUN[@]}" python -m pip freeze 2>/dev/null
  echo "#"
  echo "# --- conda list ---"
  conda list -n "${ENV_NAME}"
} > "${OUT}"

echo "Wrote ${OUT}"
echo "Commit it:  git add env_logs/ && git commit -m \"env(${ENV_NAME}): snapshot ${DATE}\" && git push"
