#!/bin/bash
# =============================================================================
# download_models.sh
# =============================================================================
# Pre-downloads all 4 encoder-only model weights to a local cache directory
# using your HuggingFace token.
#
# Run this ONCE from a login node (or any node with internet access)
# BEFORE submitting training jobs — SLURM compute nodes may not have internet.
#
# Usage:
#   export HF_TOKEN="hf_xxxxxxxxxxxxxxxxxxxx"   # your HF read token
#   bash download_models.sh
#
# Optional: override cache directory (default: ~/.cache/huggingface)
# On CSD3 your home quota may be tight — point to RDS instead:
#   export HF_HOME="/rds/user/YOUR_USERNAME/hpc-work/hf_cache"
#   bash download_models.sh
# =============================================================================

set -e
conda activate clinical 2>/dev/null || true   # activate if not already active

# --- Config -------------------------------------------------------------------
MODELS=(
    "answerdotai/ModernBERT-base"
    "answerdotai/ModernBERT-large"
    "thomas-sounack/BioClinical-ModernBERT-base"
    "thomas-sounack/BioClinical-ModernBERT-large"
)

# HF cache — override with HF_HOME env var if needed
HF_CACHE="${HF_HOME:-$HOME/.cache/huggingface}"
echo "HuggingFace cache directory: ${HF_CACHE}"
echo ""

# --- Login check --------------------------------------------------------------
if [ -z "${HF_TOKEN}" ]; then
    echo "ERROR: HF_TOKEN is not set."
    echo "Get your read token from: https://huggingface.co/settings/tokens"
    echo "Then run:  export HF_TOKEN=\"hf_xxxx\""
    exit 1
fi

echo "Logging in to HuggingFace..."
huggingface-cli login --token "${HF_TOKEN}" --add-to-git-credential 2>/dev/null || true

# --- Download -----------------------------------------------------------------
for MODEL_ID in "${MODELS[@]}"; do
    echo "------------------------------------------------------------"
    echo "Downloading: ${MODEL_ID}"
    echo "------------------------------------------------------------"
    huggingface-cli download "${MODEL_ID}" \
        --include "*.json" "*.safetensors" "tokenizer*" "config*" "vocab*" \
        --cache-dir "${HF_CACHE}"
    echo "  Done: ${MODEL_ID}"
    echo ""
done

# --- Verify -------------------------------------------------------------------
echo "============================================================"
echo "  Verifying downloads..."
echo "============================================================"
python - <<'EOF'
from transformers import AutoTokenizer, AutoConfig
import os

HF_CACHE = os.environ.get("HF_HOME", os.path.expanduser("~/.cache/huggingface"))

models = [
    "answerdotai/ModernBERT-base",
    "answerdotai/ModernBERT-large",
    "thomas-sounack/BioClinical-ModernBERT-base",
    "thomas-sounack/BioClinical-ModernBERT-large",
]

all_ok = True
for m in models:
    try:
        cfg = AutoConfig.from_pretrained(m, cache_dir=HF_CACHE)
        tok = AutoTokenizer.from_pretrained(m, cache_dir=HF_CACHE)
        print(f"  OK  {m}  (hidden_size={cfg.hidden_size}, vocab={tok.vocab_size})")
    except Exception as e:
        print(f"  FAIL {m}: {e}")
        all_ok = False

if all_ok:
    print("\nAll models downloaded successfully.")
else:
    print("\nSome models failed — check errors above.")
    exit(1)
EOF

echo ""
echo "============================================================"
echo "  All models ready. You can now submit training jobs."
echo "  Cache: ${HF_CACHE}"
echo "============================================================"
