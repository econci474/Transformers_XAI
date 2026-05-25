# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This software may be used and distributed in accordance with
# the terms of the DINOv3 License Agreement.
#
# Vendored copy: imports rewritten as package-relative; the
# `convert_linears_to_fp8` symbol from upstream `fp8_linear.py` is
# stubbed as a no-op (we don't use FP8 inference downstream).

from .attention import CausalSelfAttention, LinearKMaskedBias, SelfAttention
from .block import CausalSelfAttentionBlock, SelfAttentionBlock
from .ffn_layers import Mlp, SwiGLUFFN
from .layer_scale import LayerScale
from .patch_embed import PatchEmbed
from .rms_norm import RMSNorm
from .rope_position_encoding import RopePositionEmbedding


def convert_linears_to_fp8(*args, **kwargs):
    """No-op stub. Upstream upgrades nn.Linear -> Fp8Linear for H100
    inference; we run downstream in fp16/bf16 AMP and never call this."""
    return None


__all__ = [
    "CausalSelfAttention", "LinearKMaskedBias", "SelfAttention",
    "CausalSelfAttentionBlock", "SelfAttentionBlock",
    "Mlp", "SwiGLUFFN",
    "LayerScale", "PatchEmbed", "RMSNorm", "RopePositionEmbedding",
    "convert_linears_to_fp8",
]
