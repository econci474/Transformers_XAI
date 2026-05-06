#!/usr/bin/env python
"""
Wrapper around bmfm-targets-run that forces SDPA attention.

The BMFM SCModernBert flash_attention_2 code path has several incompatibilities
with the installed flash-attn version and the multi-field embedding architecture:
  1) rope_theta UnboundLocalError on global-attention layers.
  2) pos_idx_in_fp32 API mismatch with certain flash-attn versions.
  3) FA2 unpadding flattens the (batch, num_fields, seq_len) input_ids to 2D,
     breaking the multi-field SCEmbeddingsLayer which expects 3D input.

Rather than patching around each issue, we force SDPA (PyTorch's native scaled
dot-product attention) which is memory-efficient, correctly handles the padded
multi-field tensor layout, and avoids all flash-attn API compatibility issues.

Combined with bf16-mixed precision and gradient accumulation, SDPA is more than
sufficient for the 113M-parameter model on A100.

Usage: python force_sdpa_wrapper.py [same args as bmfm-targets-run]
"""
import importlib
import sys

# ── Import the BMFM module ───────────────────────────────────────────────────
mod = importlib.import_module(
    "bmfm_targets.models.predictive.scmodernbert.modeling_scmodernbert"
)

# ── Patch: force SDPA on the top-level model ─────────────────────────────────
_OrigCls = mod.SCModernBertForMultiTaskModeling
_orig_init = _OrigCls.__init__

def _sdpa_init(self, config, *args, **kwargs):
    config._attn_implementation = "sdpa"
    print("[force_sdpa_wrapper] Forced _attn_implementation → sdpa")
    return _orig_init(self, config, *args, **kwargs)

_OrigCls.__init__ = _sdpa_init

# ── Delegate to the normal BMFM entry point ──────────────────────────────────
from bmfm_targets.tasks.scbert.scbert_main import main  # noqa: E402
main()
