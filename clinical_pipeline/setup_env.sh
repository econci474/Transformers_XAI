#!/bin/bash
# =============================================================================
# setup_env.sh
# =============================================================================
# Creates the 'clinical' conda environment on the department cluster or CSD3.
# Run this ONCE from a login node (or interactive session with internet access).
#
# Usage:
#   bash setup_env.sh
#
# On CSD3, load CUDA module first:
#   module load cuda/12.1  (or: module load cuda/12.1.1-iccifort-2023.1.0)
#   bash setup_env.sh
# =============================================================================

set -e

ENV_NAME="clinical"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "============================================================"
echo "  Setting up conda environment: ${ENV_NAME}"
echo "============================================================"

# --- 1. Create base environment from environment.yml ---------------------------
echo "[1/4] Creating conda env from environment.yml (this may take ~5 min)..."
conda env create -f "${SCRIPT_DIR}/environment.yml" --name "${ENV_NAME}" -y \
    || conda env update -f "${SCRIPT_DIR}/environment.yml" --name "${ENV_NAME}"

# --- 2. Activate and verify torch + CUDA ---------------------------------------
echo "[2/4] Verifying torch + CUDA..."
conda run -n "${ENV_NAME}" python -c "
import torch
print(f'  torch version : {torch.__version__}')
print(f'  CUDA available: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'  GPU           : {torch.cuda.get_device_name(0)}')
    print(f'  CUDA version  : {torch.version.cuda}')
else:
    print('  WARNING: CUDA not available — check your CUDA module is loaded')
"

# --- 3. Verify transformers version (must be >= 4.47 for ModernBERT) ----------
echo "[3/4] Verifying transformers version..."
conda run -n "${ENV_NAME}" python -c "
import transformers
v = transformers.__version__
major, minor, *_ = v.split('.')
assert int(major) > 4 or (int(major) == 4 and int(minor) >= 47), \
    f'transformers {v} is too old — ModernBERT requires >= 4.47.0'
print(f'  transformers  : {v}  OK')
"

# --- 4. Flash-attention check (optional) --------------------------------------
echo "[4/4] Checking flash-attn..."
conda run -n "${ENV_NAME}" python -c "
try:
    import flash_attn
    print(f'  flash-attn    : {flash_attn.__version__}  OK')
except ImportError:
    print('  flash-attn    : NOT installed (training will still work, just slower)')
"

echo ""
echo "============================================================"
echo "  Environment '${ENV_NAME}' is ready."
echo "  Activate with:  conda activate ${ENV_NAME}"
echo "============================================================"
