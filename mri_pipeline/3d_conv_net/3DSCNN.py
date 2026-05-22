"""
3DSCNN.py
=========
Separable-convolution single-stream 3D CNN for MRI-only Alzheimer's-disease
classification on the ADNI cohort — the *parameter-efficient* variant.

This reproduces the headline idea of:

    Spasov et al. 2019, NeuroImage 189:276-287,
    "A parameter-efficient deep learning approach to predict conversion
     from mild cognitive impairment to Alzheimer's disease."
    Code: https://github.com/simeon-spasov/MCI  (utils/sepconv3D.py)

The network is INTENTIONALLY identical in topology to the vanilla baseline
in the sibling file ``3DCNN.py`` — same number of blocks, same channel
counts (1->24->48->96->48->8), same kernel sizes, same pooling, same
classifier head. The ONLY difference is that every standard ``nn.Conv3d``
is replaced by a depthwise-separable ``SeparableConv3d``. Keeping every
other choice fixed means any performance difference between the two
networks can be attributed to the separable convolution alone — a clean,
controlled comparison.

Same two thesis-specific deviations from the published network as the
vanilla file: MRI ONLY (no Jacobian-determinant stream, no clinical
features) and a SINGLE stream (no redundant two-stream topology).

Input
-----
A single-channel 3D volume of shape ``[B, 1, 193, 229, 193]`` — the
MNI152NLin2009cAsym res-1 (1 mm) voxel grid from sMRIprep, z-scored on the
non-zero (brain) voxels by the preprocessing step. No normalisation is
done inside this module.

Output
------
``n_outputs`` raw logits.
  * Binary task (AD-visit vs stable-CN): ``n_outputs=1`` + BCEWithLogitsLoss.
  * K-class task: ``n_outputs=K`` + CrossEntropyLoss.

Scope
-----
MODEL DEFINITION ONLY — no data loading, no training loop. A separate
shared training script imports ``MRI3DSeparableCNN`` from here via
``importlib`` (the file name starts with a digit, so a plain
``import 3DSCNN`` is not valid Python).

Approximate size: ~0.03 M trainable parameters — roughly 20x smaller than
the ~0.5 M of the vanilla 3DCNN.py, because the separable factorisation
removes most of the weights from the spatial convolutions.
"""

import torch
import torch.nn as nn


# Per-block dropout probability — see the note in 3DCNN.py. Spasov's public
# training script uses 0.1; kept as the tunable default here too.
DEFAULT_DROPOUT = 0.1


class SeparableConv3d(nn.Module):
    """
    Depthwise-separable 3D convolution = a depthwise conv + a pointwise conv.

    WHY THIS EXISTS
    ---------------
    A standard ``Conv3d`` with C_in input channels, C_out output channels
    and a k x k x k kernel holds

        C_in * C_out * k^3

    weights. It does TWO jobs at once: it mixes information across SPACE
    (the k^3 kernel) and across CHANNELS (the C_in x C_out mixing matrix).

    A separable convolution factorises those two jobs into two cheaper
    layers:

      1. DEPTHWISE conv — one independent k^3 spatial filter PER input
         channel. Implemented with ``groups=C_in`` so channel i is
         convolved only with its own filter. It mixes SPACE but not
         channels. Cost: C_in * k^3 weights.

      2. POINTWISE conv — a 1 x 1 x 1 convolution. It mixes CHANNELS but
         not space. Cost: C_in * C_out weights.

    Total: ``C_in*k^3 + C_in*C_out`` instead of ``C_in*C_out*k^3``. For the
    block2 layer of this network (C_in=24, C_out=48, k=5x6x5 -> 150 voxels)
    that is 24*150 + 24*48 = 4,752 weights versus 24*48*150 = 172,800 — a
    ~36x reduction. Summed over the whole network the separable variant has
    ~0.03 M parameters versus ~0.5 M for the vanilla one. This is exactly
    the "parameter efficiency" Spasov et al. report.

    FIRST-LAYER NOTE
    ----------------
    When C_in == 1 the depthwise conv has ``groups=1``, i.e. it is just an
    ordinary single-channel spatial convolution; the factorisation only
    starts saving weights once C_in > 1. The first block is therefore
    effectively a normal spatial conv (1->1) followed by a 1x1x1 channel
    expansion (1->24). This is harmless and keeps the code uniform.
    """

    def __init__(self, in_channels, out_channels, kernel_size,
                 stride=1, padding=0, bias=False):
        super().__init__()

        # ── Depthwise convolution ────────────────────────────────────────
        # in_channels == out_channels here, and groups == in_channels, so
        # each input channel gets exactly one k^3 spatial filter and is
        # processed completely independently of the others. The stride and
        # padding are applied HERE (the spatial step happens in this conv).
        self.depthwise = nn.Conv3d(
            in_channels, in_channels, kernel_size,
            stride=stride, padding=padding,
            groups=in_channels,   # <-- the key line that makes it depthwise
            bias=bias)

        # ── Pointwise convolution ────────────────────────────────────────
        # A 1x1x1 conv. It does not look at spatial neighbours at all; its
        # only job is to mix the C_in depthwise outputs into C_out channels.
        self.pointwise = nn.Conv3d(
            in_channels, out_channels, kernel_size=1,
            stride=1, padding=0, bias=bias)

    def forward(self, x):
        # Spatial filtering first, then channel mixing.
        return self.pointwise(self.depthwise(x))


class SeparableConvBlock3D(nn.Module):
    """
    One convolutional stage of the separable network:

        SeparableConv3d -> BatchNorm3d -> ELU -> (MaxPool3d) -> Dropout3d

    This is byte-for-byte the same block as ``ConvBlock3D`` in 3DCNN.py
    EXCEPT that the leading ``nn.Conv3d`` is swapped for a
    ``SeparableConv3d``. Keeping BatchNorm / ELU / MaxPool / Dropout
    identical is deliberate: it guarantees the vanilla and separable
    networks differ ONLY in the convolution type.

    See ConvBlock3D in 3DCNN.py for the rationale behind each component
    (bias=False because BatchNorm follows, ELU activation, MaxPool with
    receptive field 3 / stride 2, channel-wise Dropout3d).
    """

    def __init__(self, in_channels, out_channels, kernel_size,
                 stride=1, dropout=DEFAULT_DROPOUT, use_pool=True):
        super().__init__()

        # Allow an int kernel/stride to mean "the same value on all 3 axes".
        if isinstance(kernel_size, int):
            kernel_size = (kernel_size,) * 3
        if isinstance(stride, int):
            stride = (stride,) * 3

        # "Same-style" padding: half the kernel width on each side (the
        # paper does not specify padding; this keeps spatial sizes
        # predictable). The padding is handed to the depthwise conv inside
        # SeparableConv3d, which is where the spatial step happens.
        padding = tuple(k // 2 for k in kernel_size)

        self.conv = SeparableConv3d(in_channels, out_channels, kernel_size,
                                    stride=stride, padding=padding,
                                    bias=False)
        self.norm = nn.BatchNorm3d(out_channels)
        self.act = nn.ELU(inplace=True)

        self.use_pool = use_pool
        self.pool = nn.MaxPool3d(kernel_size=3, stride=2, padding=1)
        self.drop = nn.Dropout3d(p=dropout)

    def forward(self, x):
        x = self.conv(x)        # depthwise + pointwise convolution
        x = self.norm(x)        # normalise per channel
        x = self.act(x)         # non-linearity
        if self.use_pool:
            x = self.pool(x)    # downsample by 2x
        x = self.drop(x)        # channel-wise dropout
        return x


class MRI3DSeparableCNN(nn.Module):
    """
    Separable-convolution single-stream 3D CNN.

    Topology is IDENTICAL to ``MRI3DCNN`` in 3DCNN.py — only the
    convolution type differs. For a [B, 1, 193, 229, 193] input the spatial
    size shrinks exactly as in the vanilla network (D x H x W after each
    block):

        input                                          1 x 193 x 229 x 193
        block1  SepConv 1 ->24  k(11,13,11) s4 + pool  24 x  25 x  29 x  25
        block2  SepConv 24->48  k(5,6,5)    s1 + pool  48 x  13 x  15 x  13
        block3  SepConv 48->96  k(3,3,3)    s1 + pool  96 x   7 x   8 x   7
        block4  SepConv 96->48  k(3,4,3)    s1 + pool  48 x   4 x   5 x   4
        block5  SepConv 48-> 8  k(3,4,3)    s1 + pool   8 x   2 x   3 x   2
        flatten                                               96 features
        fc1     Linear 96 -> 32  (+ BN + ELU + drop)
        head    Linear 32 -> n_outputs                        logits
    """

    def __init__(self, in_channels=1, n_outputs=1, dropout=DEFAULT_DROPOUT):
        """
        Args:
            in_channels: number of input image channels. 1 for a single
                T1-weighted MRI volume.
            n_outputs:   number of output logits. 1 for the binary
                AD-visit-vs-stable-CN task (train with BCEWithLogitsLoss);
                set to K for a K-class task (train with CrossEntropyLoss).
            dropout:     dropout probability used in every conv block and
                in the classifier head.
        """
        super().__init__()

        # ── Convolutional feature extractor ──────────────────────────────
        # Same channel progression and kernels as the vanilla net; every
        # block uses a SeparableConvBlock3D instead of a plain ConvBlock3D.
        self.block1 = SeparableConvBlock3D(in_channels, 24, (11, 13, 11),
                                           stride=4, dropout=dropout)
        self.block2 = SeparableConvBlock3D(24, 48, (5, 6, 5),
                                           stride=1, dropout=dropout)
        self.block3 = SeparableConvBlock3D(48, 96, (3, 3, 3),
                                           stride=1, dropout=dropout)
        self.block4 = SeparableConvBlock3D(96, 48, (3, 4, 3),
                                           stride=1, dropout=dropout)
        self.block5 = SeparableConvBlock3D(48, 8, (3, 4, 3),
                                           stride=1, dropout=dropout)

        # ── Classifier head ──────────────────────────────────────────────
        self.flatten = nn.Flatten()
        # LazyLinear infers the flattened feature length on the first
        # forward pass (see the note in 3DCNN.py).
        self.fc1 = nn.LazyLinear(32)
        self.fc_norm = nn.BatchNorm1d(32)
        self.fc_act = nn.ELU(inplace=True)
        self.fc_drop = nn.Dropout(p=dropout)

        # Final linear layer producing raw logits — no activation here, the
        # loss function applies sigmoid/softmax internally.
        self.head = nn.Linear(32, n_outputs)

    def forward(self, x):
        """
        Args:
            x: input tensor of shape [B, in_channels, 193, 229, 193].
        Returns:
            logits of shape [B, n_outputs].
        """
        # Convolutional feature extraction.
        x = self.block1(x)
        x = self.block2(x)
        x = self.block3(x)
        x = self.block4(x)
        x = self.block5(x)

        # Flatten the [B, 8, 2, 3, 2] feature map to [B, 96].
        x = self.flatten(x)

        # Classifier head.
        x = self.fc1(x)
        x = self.fc_norm(x)
        x = self.fc_act(x)
        x = self.fc_drop(x)

        # Raw logits — apply sigmoid/softmax inside the loss, not here.
        return self.head(x)
