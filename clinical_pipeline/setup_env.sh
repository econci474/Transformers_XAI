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
# torch is PINNED to 2.3.x: torch 2.4.x NaNs ModernBERT full fine-tune (a compiled/
# SDPA backward regression — frozen runs are spared, full_ft NaNs from step 1). The
# env once silently drifted 2.3.1 -> 2.4.1 and broke the sweep; ASSERT, don't just print.
echo "[2/4] Verifying torch + CUDA (torch must be 2.3.x)..."
conda run -n "${ENV_NAME}" python -c "
import torch
v = torch.__version__
print(f'  torch version : {v}')
print(f'  CUDA available: {torch.cuda.is_available()}')
assert v.startswith('2.3.'), (
    f'torch {v} != pinned 2.3.x — torch 2.4+ NaNs ModernBERT full fine-tune. '
    'Reinstall: pip install torch==2.3.1 torchvision==0.18.1 '
    '--index-url https://download.pytorch.org/whl/cu121')
print('  torch pin     : 2.3.x  OK')
if torch.cuda.is_available():
    print(f'  GPU           : {torch.cuda.get_device_name(0)}')
    print(f'  CUDA version  : {torch.version.cuda}')
else:
    print('  NOTE: CUDA not available here (expected on a login node).')
"

# --- 3. Verify transformers version (must be >= 4.47 AND < 5 for ModernBERT) ---
# transformers 5.x is a breaking major release (new Trainer / checkpoint-save API)
# that crashes 03_encoder_finetune.py. ModernBERT needs >= 4.47. Enforce 4.47 <= v < 5.
echo "[3/4] Verifying transformers version (4.47 <= v < 5)..."
conda run -n "${ENV_NAME}" python -c "
import transformers
v = transformers.__version__
major, minor, *_ = v.split('.')
major, minor = int(major), int(minor)
assert major == 4 and minor >= 47, \
    f'transformers {v} is incompatible — need >= 4.47.0 AND < 5.0 ' \
    f'(5.x breaks the Trainer/checkpoint API used by 03_encoder_finetune.py)'
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
