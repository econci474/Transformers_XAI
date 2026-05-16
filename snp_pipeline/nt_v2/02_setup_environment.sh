#!/bin/bash
# 02_setup_environment.sh
# ========================
# Set up the ntv2 conda environment and download the NT v2 model.
#
# Usage:
#   bash nt_v2/02_setup_environment.sh

set -e

echo "============================================================"
echo " Nucleotide Transformer v2 — Environment Setup"
echo "============================================================"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="$SCRIPT_DIR/environment.yml"

# ── 1) Create conda environment ──────────────────────────────────────────────
echo ""
echo "Step 1: Creating conda environment from $ENV_FILE"
if conda env list | grep -q "ntv2"; then
    echo "  Environment 'ntv2' already exists. Updating..."
    conda env update -f "$ENV_FILE" --prune
else
    conda env create -f "$ENV_FILE"
fi

echo ""
echo "Step 2: Activate and verify"
echo "  Run:  conda activate ntv2"

# ── 2) Download the model ────────────────────────────────────────────────────
echo ""
echo "Step 3: Downloading NT v2 model (this will cache to ~/.cache/huggingface)"
echo "  Model: InstaDeepAI/nucleotide-transformer-v2-500m-multi-species"

# Activate env in subshell for download
eval "$(conda shell.bash hook)"
conda activate ntv2

python -c "
from transformers import AutoTokenizer, AutoModelForMaskedLM
import torch

MODEL_ID = 'InstaDeepAI/nucleotide-transformer-v2-500m-multi-species'

print('Downloading tokenizer...')
tok = AutoTokenizer.from_pretrained(MODEL_ID, trust_remote_code=True)
print(f'  Vocab size: {tok.vocab_size}')
print(f'  Model max length: {tok.model_max_length}')

print('Downloading model...')
model = AutoModelForMaskedLM.from_pretrained(MODEL_ID, trust_remote_code=True)
n_params = sum(p.numel() for p in model.parameters())
print(f'  Parameters: {n_params:,}')

# Quick forward pass sanity check
print('Running sanity check...')
seq = 'ACGTACGTACGTACGTACGT'
ids = tok(seq, return_tensors='pt')['input_ids']
with torch.no_grad():
    out = model(ids, output_hidden_states=True)
emb = out.hidden_states[-1]
print(f'  Input: {seq[:20]}...')
print(f'  Tokens: {ids.shape[1]}')
print(f'  Hidden states shape: {emb.shape}')
print(f'  Hidden dim: {emb.shape[-1]}')
print('✓ Model downloaded and verified successfully!')
"

echo ""
echo "============================================================"
echo " Setup complete!"
echo " To use:  conda activate ntv2"
echo "============================================================"
