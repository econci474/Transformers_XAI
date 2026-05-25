"""
AGMS3DCNN.py
============
Attention-Gated Multi-Scale 3D CNN classifier for AD diagnosis on ADNI.

Adapted from:
    Nature Sci. Rep. 2025, "AG-MS3D-CNN ..." (PDF in this folder)
    https://www.nature.com/articles/s41598-025-09351-x

The original paper is a multi-task brain-tumor model on BraTS (4-channel
multi-modal MRI; segmentation + classification + volume estimation heads).
This file ports the *backbone + attention gate + classification head* to the
ADNI diagnostic harness as instructed in
``mri_pipeline/3d_cnn_vit/3d_cnn_vit.txt``:

    "Adapt this architecture from Nature paper: AG-MS3D-CNN for the
     classification branch only."

Adaptations vs. the paper
-------------------------
- Single-channel input (T1 only)         <- paper uses 4-channel (T1/T1ce/T2/FLAIR)
- 193 x 229 x 193, sMRIprep MNI grid     <- paper uses BraTS 240x240x155
- Depthwise-separable Conv3D backbone    <- paper uses standard Conv3D
- Single classification head             <- paper has seg+cls+vol heads
- Weighted CrossEntropyLoss              <- paper uses Dice/CE/MSE compound loss
- No GRL adversarial training            <- not applicable; single domain
- No MC dropout inference                <- not part of this evaluation

Architecture
------------
Three parallel pathways at three input scales (fine / medium / coarse),
each a 5-block depthwise-separable CNN. Per the paper, the filter count
doubles only at the two downsampling steps (stride-2 convs at B2 and B4),
starting from 32 in the fine-scale block (default base=32):

    B1: in    -> 32          (no DS)
    B2: 32    -> 64          (stride=2 DS, doubles)
    B3: 64    -> 64          (no DS)
    B4: 64    -> 128         (stride=2 DS, doubles)
    B5: 128   -> 128         (no DS)

The medium and coarse pathways downsample the input by avg-pool 2x and 4x
respectively before entering the pathway; their outputs are trilinearly
resampled back to the fine pathway's spatial size and fused as a
weighted sum with learnable per-pathway weights (softmaxed at forward).

A CBAM-style 3D Spatial Attention module (avg- and max-pool over the
channel axis, concat, 7x7x7 conv, sigmoid) gates the fused feature map.
Global average pooling collapses spatial dims, then a small MLP head
(256 -> 128 -> n_outputs) emits raw logits.

Note on the "attention" in the paper: it is **not** Transformer
self-attention (Q/K/V softmax). It is a one-conv spatial gating map -
O(C*D*H*W), not O(N^2) - which is why it is feasible at 193^3.

Input
-----
Tensor of shape [B, in_channels, 193, 229, 193] (z-scored on non-zero
voxels by the preprocessing step 00_prepare_CNN_inputs.py). The model
expects to consume the same cnn_inputs/ tree as ``3d_conv_net``.

Output
------
``n_outputs`` raw logits per sample. Train with ``nn.CrossEntropyLoss``
(n_outputs=2 for binary, n_outputs=3 for multiclass) - same convention
as the rest of the MRI pipeline (ViT, Spasov, BrainMVP).

Scope
-----
MODEL DEFINITION ONLY - no data loading, no training loop. The training
script ``train_agms3d.py`` (sibling) loads this via ``importlib`` and
materialises the LazyLinear head with one dummy forward before the
optimiser/checkpoint code touches it.
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


DEFAULT_DROPOUT = 0.1


# =============================================================================
#                       Depthwise-separable 3D convolution
# =============================================================================
class SeparableConv3d(nn.Module):
    """Depthwise-separable 3D conv = depthwise (per-channel spatial) + 1x1x1
    pointwise (channel mixing). For C_in x C_out x k^3 weights of a standard
    Conv3d, the separable factorisation costs C_in*k^3 + C_in*C_out -- in our
    3-pathway design that is roughly 5-8x fewer params overall.

    Duplicated from ``3d_conv_net/3DSCNN.py:SeparableConv3d`` (the plan keeps
    ``3d_cnn_vit/`` self-contained to avoid the cross-folder importlib
    coupling that has bitten BrainMVP <-> ViT this session).
    """

    def __init__(self, in_channels, out_channels, kernel_size,
                 stride=1, padding=0, bias=False):
        super().__init__()
        self.depthwise = nn.Conv3d(
            in_channels, in_channels, kernel_size,
            stride=stride, padding=padding,
            groups=in_channels,           # <-- the key line that makes it depthwise
            bias=bias)
        self.pointwise = nn.Conv3d(
            in_channels, out_channels, kernel_size=1,
            stride=1, padding=0, bias=bias)

    def forward(self, x):
        return self.pointwise(self.depthwise(x))


# =============================================================================
#                                Pathway block
# =============================================================================
class SepBlock3D(nn.Module):
    """One pathway block, implementing the paper's residual formulation.

        Eq. (1):  F_block = ReLU( BN( Conv( F_in ) ) )      -- block transform
        Eq. (3):  F_l = F_{l-1} + Conv3D( F_{l-1} )         -- residual sum

    Combined (and using ELU in place of ReLU per user preference, matching
    the Spasov 3DSCNN template):

        F_out = Dropout3d( F_in_proj + ELU( BN( SepConv( F_in ) ) ) )

    Every block has a residual (the paper's "residual connections within
    each pathway" applies pathway-wide). When the block changes shape
    (channel count differs from stride>1 spatial downsampling), the
    shortcut uses a 1x1x1 projection conv with matching stride to align
    the input to the output's shape -- the standard ResNet trick. When
    shape is preserved, the shortcut is the identity."""

    def __init__(self, in_channels, out_channels, kernel_size=3,
                 stride=1, dropout=DEFAULT_DROPOUT, backbone: str = "separable"):
        super().__init__()
        padding = kernel_size // 2
        if backbone == "separable":
            self.conv = SeparableConv3d(
                in_channels, out_channels, kernel_size,
                stride=stride, padding=padding, bias=False)
        elif backbone == "vanilla":
            self.conv = nn.Conv3d(
                in_channels, out_channels, kernel_size,
                stride=stride, padding=padding, bias=False)
        else:
            raise ValueError(
                f"backbone must be 'separable' or 'vanilla', got {backbone!r}")
        self.norm = nn.BatchNorm3d(out_channels)
        self.act = nn.ELU(inplace=True)
        self.drop = nn.Dropout3d(p=dropout)

        # Residual shortcut: identity if shape preserved, else 1x1x1
        # projection with matching stride (paper Eq. 3 implies every block).
        if in_channels == out_channels and stride == 1:
            self.shortcut = nn.Identity()
        else:
            self.shortcut = nn.Sequential(
                nn.Conv3d(in_channels, out_channels, kernel_size=1,
                          stride=stride, bias=False),
                nn.BatchNorm3d(out_channels),
            )

    def forward(self, x):
        y = self.act(self.norm(self.conv(x)))
        y = y + self.shortcut(x)              # paper Eq. (3): residual sum
        y = self.drop(y)                      # dropout after the sum
        return y


class Pathway(nn.Module):
    """A 5-block separable-CNN pathway.

    Channel progression (per paper: "doubling the number of filters after every
    downsampling, starting from 32 in fine-scale"):
        B1: in -> base       (no DS)
        B2: base -> 2*base   (stride=2, DS, channels double)
        B3: 2*base -> 2*base (no DS)
        B4: 2*base -> 4*base (stride=2, DS, channels double)
        B5: 4*base -> 4*base (no DS)

    With base=32 (paper default): 32 -> 64 -> 64 -> 128 -> 128.
    Stride pattern:               1  -> 2  -> 1  -> 2   -> 1    (~4x intra-pathway DS)
    """

    def __init__(self, in_channels=1, base=32, dropout=DEFAULT_DROPOUT,
                 backbone: str = "separable"):
        super().__init__()
        c1, c2, c3, c4, c5 = base, base * 2, base * 2, base * 4, base * 4
        kw = dict(kernel_size=3, dropout=dropout, backbone=backbone)
        self.b1 = SepBlock3D(in_channels, c1, stride=1, **kw)
        self.b2 = SepBlock3D(c1,          c2, stride=2, **kw)
        self.b3 = SepBlock3D(c2,          c3, stride=1, **kw)
        self.b4 = SepBlock3D(c3,          c4, stride=2, **kw)
        self.b5 = SepBlock3D(c4,          c5, stride=1, **kw)
        self.out_channels = c5

    def forward(self, x):
        x = self.b1(x)
        x = self.b2(x)
        x = self.b3(x)
        x = self.b4(x)
        x = self.b5(x)
        return x


# =============================================================================
#                       CBAM-style 3D Spatial Attention
# =============================================================================
class SpatialAttention3D(nn.Module):
    """CBAM 3D spatial gating.

        avg, mx = x.mean(C, keepdim=True), x.amax(C, keepdim=True)
        attn    = sigmoid(Conv3d_1ch(k=7)(cat([avg, mx], dim=1)))
        return x * attn

    Cost: one 7x7x7 conv on a 2-channel input. NOT Transformer self-attention.
    """

    def __init__(self, kernel_size: int = 7):
        super().__init__()
        self.conv = nn.Conv3d(2, 1, kernel_size=kernel_size,
                              padding=kernel_size // 2, bias=False)

    def forward(self, x):
        avg = x.mean(dim=1, keepdim=True)
        mx  = x.amax(dim=1, keepdim=True)
        attn = torch.sigmoid(self.conv(torch.cat([avg, mx], dim=1)))
        return x * attn


# =============================================================================
#                                AG-MS3D-CNN
# =============================================================================
class AGMS3DCNN(nn.Module):
    """Classification-only adaptation of AG-MS3D-CNN for ADNI AD diagnosis.

    Args:
        in_channels: 1 for single-channel T1 (default).
        n_outputs:   number of classes (2 for binary tasks, 3 for T2_multiclass).
                     Train with ``nn.CrossEntropyLoss`` (same convention as
                     04_supervised_finetuning_ViT.py and train_3dcnn.py).
        dropout:     per-block and head dropout probability.
        base:        starting channel count (paper default: 32). Channels
                     double only at the two downsampling steps -> reaches
                     4*base in the deepest blocks.
    """

    def __init__(self, in_channels: int = 1, n_outputs: int = 2,
                 dropout: float = DEFAULT_DROPOUT, base: int = 32,
                 backbone: str = "separable"):
        super().__init__()
        kw = dict(base=base, dropout=dropout, backbone=backbone)
        self.fine   = Pathway(in_channels, **kw)
        self.medium = Pathway(in_channels, **kw)
        self.coarse = Pathway(in_channels, **kw)
        self.backbone = backbone

        # Learnable pathway fusion weights (softmaxed in forward so they sum to 1).
        self.alpha = nn.Parameter(torch.ones(3) / 3.0)

        self.attn = SpatialAttention3D(kernel_size=7)
        self.gap = nn.AdaptiveAvgPool3d(1)

        # Classification head per paper: GAP -> FC(1024) -> FC(256) -> FC(n_outputs).
        # LazyLinear adapts to whatever channel count the pathway emits
        # (initialised on the first forward; the training script does one
        # eval+no_grad() dummy forward before the optimiser is built).
        # NOTE on capacity: the paper's 1024 -> 256 head adds ~393k params vs
        # a 256 -> 128 head. With ~600 train subjects this may overfit; if
        # validation balanced-acc plateaus low / train-val gap blows up,
        # consider reducing the head to LazyLinear(256) -> Linear(128) ->
        # Linear(n_outputs) and noting the deviation in the writeup.
        self.head = nn.Sequential(
            nn.Flatten(),
            nn.LazyLinear(1024), nn.ELU(inplace=True), nn.Dropout(p=dropout),
            nn.Linear(1024, 256), nn.ELU(inplace=True), nn.Dropout(p=dropout),
            nn.Linear(256, n_outputs),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Args:
            x: [B, in_channels, D, H, W] -- expects D, H, W ~ 193, 229, 193.
        Returns:
            [B, n_outputs] raw logits.
        """
        # Pre-pool the inputs to medium and coarse pathways (paper: 2x, 4x).
        x_m = F.avg_pool3d(x, kernel_size=2)
        x_c = F.avg_pool3d(x, kernel_size=4)

        f = self.fine(x)        # [B, 4b, ~D/4, ~H/4, ~W/4]
        m = self.medium(x_m)    # [B, 4b, ~D/8, ~H/8, ~W/8]
        c = self.coarse(x_c)    # [B, 4b, ~D/16, ~H/16, ~W/16]

        # Trilinear-resample medium and coarse back to the fine pathway's grid
        # so the weighted sum lines up voxel-wise.
        target = f.shape[2:]
        if m.shape[2:] != target:
            m = F.interpolate(m, size=target, mode="trilinear", align_corners=False)
        if c.shape[2:] != target:
            c = F.interpolate(c, size=target, mode="trilinear", align_corners=False)

        w = torch.softmax(self.alpha, dim=0)
        fused = w[0] * f + w[1] * m + w[2] * c    # [B, 4b, ~D/4, ~H/4, ~W/4]
        gated = self.attn(fused)                  # CBAM-style spatial gating
        pooled = self.gap(gated)                  # [B, 4b, 1, 1, 1]
        return self.head(pooled)                  # [B, n_outputs]
