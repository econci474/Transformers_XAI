#!/usr/bin/env python
"""
Wrapper around bmfm-targets-run that forces flash_attention_2.

Fixes two issues with the BMFM SCModernBert implementation:
  1) Forces _attn_implementation = "flash_attention_2" on the model config.
  2) Patches SCModernBertAttention.__init__ to initialise rope_theta for
     *global* attention layers (the upstream code only sets rope_theta inside
     the `if self.local_attention != (-1, -1)` branch, causing an
     UnboundLocalError for global-attention layers when FA2 is enabled).

Usage: python force_fa2_wrapper.py [same args as bmfm-targets-run]
"""
import importlib
import sys

# ── Import the BMFM module ───────────────────────────────────────────────────
mod = importlib.import_module(
    "bmfm_targets.models.predictive.scmodernbert.modeling_scmodernbert"
)

# ── Patch 1: fix rope_theta UnboundLocalError in SCModernBertAttention ────────
# The upstream __init__ only assigns rope_theta when self.local_attention !=
# (-1, -1).  For global-attention layers (local_attention == (-1, -1)) the
# variable is never set, but the flash_attention_2 branch still references it.
#
# Fix: replace __init__ with a corrected copy that defaults
#      rope_theta = config.global_rope_theta before the conditional.
_OrigAttnCls = mod.SCModernBertAttention

def _safe_attn_init(self, config, layer_id=None):
    """SCModernBertAttention.__init__ with rope_theta fix for global layers."""
    from torch import nn

    super(_OrigAttnCls, self).__init__()
    self.config = config
    self.layer_id = layer_id

    if config.hidden_size % config.num_attention_heads != 0:
        raise ValueError(
            f"The hidden size ({config.hidden_size}) is not a multiple of "
            f"the number of attention heads ({config.num_attention_heads})"
        )

    self.attention_dropout = config.attention_dropout
    self.deterministic_flash_attn = config.deterministic_flash_attn
    self.num_heads = config.num_attention_heads
    self.head_dim = config.hidden_size // config.num_attention_heads
    self.all_head_size = self.head_dim * self.num_heads
    self.Wqkv = nn.Linear(
        config.hidden_size, 3 * self.all_head_size, bias=config.attention_bias
    )

    if layer_id % config.global_attn_every_n_layers != 0:
        self.local_attention = (
            config.local_attention // 2,
            config.local_attention // 2,
        )
    else:
        self.local_attention = (-1, -1)

    max_position_embeddings = config.max_position_embeddings

    # ── FIX: always initialise rope_theta (default = global_rope_theta) ──
    rope_theta = config.global_rope_theta

    if self.local_attention != (-1, -1):
        rope_theta = (
            config.global_rope_theta
            if config.local_rope_theta is None
            else config.local_rope_theta
        )
        max_position_embeddings = config.local_attention

    if config._attn_implementation == "flash_attention_2":
        self.rotary_emb = mod.SCModernBertUnpaddedRotaryEmbedding(
            dim=self.head_dim, max_seqlen=max_position_embeddings, base=rope_theta
        )
    else:
        self.rotary_emb = mod.SCModernBertRotaryEmbedding(config=config)

    self.Wo = nn.Linear(
        config.hidden_size, config.hidden_size, bias=config.attention_bias
    )
    self.out_drop = (
        nn.Dropout(config.attention_dropout)
        if config.attention_dropout > 0.0
        else nn.Identity()
    )
    self.pruned_heads = set()

_OrigAttnCls.__init__ = _safe_attn_init
print("[force_fa2_wrapper] Patched SCModernBertAttention.__init__ "
      "(rope_theta fix for global-attention layers)")

# ── Patch 2: force flash_attention_2 on the top-level model ──────────────────
_OrigCls = mod.SCModernBertForMultiTaskModeling
_orig_init = _OrigCls.__init__

def _fa2_init(self, config, *args, **kwargs):
    config._attn_implementation = "flash_attention_2"
    print("[force_fa2_wrapper] Forced _attn_implementation → flash_attention_2")
    return _orig_init(self, config, *args, **kwargs)

_OrigCls.__init__ = _fa2_init

# ── Delegate to the normal BMFM entry point ──────────────────────────────────
from bmfm_targets.tasks.scbert.scbert_main import main  # noqa: E402
main()
