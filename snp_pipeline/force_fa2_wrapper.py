#!/usr/bin/env python
"""
Wrapper around bmfm-targets-run that forces flash_attention_2.
Usage: python force_fa2_wrapper.py [same args as bmfm-targets-run]
"""
import importlib
import sys

# ── Monkey-patch: force flash_attention_2 on SCModernBert models ──────────────
mod = importlib.import_module(
    "bmfm_targets.models.predictive.scmodernbert.modeling_scmodernbert"
)

# Patch the top-level model so config._attn_implementation is set before any
# encoder layer is constructed.
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
