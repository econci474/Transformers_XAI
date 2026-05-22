"""
3DCNN.py
========
Vanilla (standard-convolution) single-stream 3D CNN for MRI-only
Alzheimer's-disease classification on the ADNI cohort.

This is a faithful-in-spirit reproduction of the convolutional backbone in:

    Spasov et al. 2019, NeuroImage 189:276-287,
    "A parameter-efficient deep learning approach to predict conversion
     from mild cognitive impairment to Alzheimer's disease."
    Code: https://github.com/simeon-spasov/MCI

Two deliberate deviations from the published network (both requested for
this thesis):

  1. MRI ONLY. The original Spasov network is multi-modal: it fuses a
     structural-MRI stream, a Jacobian-determinant-map stream, and a
     tabular clinical/demographic stream. Here we keep a SINGLE
     structural-MRI stream — no Jacobian maps, no clinical features.
     The earlier draft of this file carried a redundant two-stream
     ("f" and "p") topology in which both streams received the SAME MRI
     volume; with no second modality that is just a wider network for no
     architectural reason, so it has been collapsed to one stream.

  2. STANDARD convolutions. The headline contribution of the 2019 paper is
     the 3D *separable* convolution (a depthwise + pointwise factorisation)
     that makes the network parameter-efficient. THIS file uses ordinary
     ``nn.Conv3d`` so it can act as the "vanilla" baseline. The separable
     variant lives in the sibling file ``3DSCNN.py`` and shares an
     identical topology, so the two can be compared cleanly.

Input
-----
A single-channel 3D volume of shape ``[B, 1, 182, 218, 182]`` — the
canonical MNI152NLin2009cAsym 1 mm voxel grid produced by sMRIprep.
Intensities are expected to be z-scored on the non-zero (brain) voxels by
the preprocessing step; this module performs NO normalisation of its own.

Output
------
``n_outputs`` raw, un-activated logits.
  * Binary task (AD-visit vs stable-CN): use ``n_outputs=1`` and train with
    ``nn.BCEWithLogitsLoss`` (this is what Spasov et al. do).
  * K-class task: set ``n_outputs=K`` and train with
    ``nn.CrossEntropyLoss``.

Scope
-----
This file is a MODEL DEFINITION ONLY — there is no data loading and no
training loop. A separate shared training script imports ``MRI3DCNN`` from
here. (Because the file name starts with a digit it cannot be imported with
a plain ``import 3DCNN``; the training script loads it via ``importlib``.)

Approximate size: ~0.5 M trainable parameters (compare ~0.03 M for the
separable variant in 3DSCNN.py).
"""

import torch
import torch.nn as nn


# Per-block dropout probability. Spasov's public training script
# (mci_train.py) uses 0.1; an earlier draft of this file quoted 0.02 from
# the paper text. It is exposed as a constructor argument so it can be
# tuned; 0.1 is the default because that is the value in the reference
# implementation.
DEFAULT_DROPOUT = 0.1


class ConvBlock3D(nn.Module):
    """
    One convolutional stage of the network:

        Conv3d -> BatchNorm3d -> ELU -> (MaxPool3d) -> Dropout3d

    Rationale for each component:

      * Conv3d(bias=False) — the convolution learns the spatial filters.
        Bias is disabled because the BatchNorm3d immediately afterwards has
        its own learnable shift parameter (beta); a separate conv bias
        would be mathematically redundant.

      * BatchNorm3d — normalises each channel's activations across the
        batch, which stabilises and speeds up training and acts as a mild
        regulariser. IMPORTANT: BatchNorm needs a batch size > 1 at train
        time, so the training DataLoader should use ``drop_last=True``.

      * ELU (Exponential Linear Unit) — the activation used by Spasov et
        al. Unlike ReLU it has smooth, non-zero gradients for negative
        inputs, which helps gradient flow.

      * MaxPool3d(kernel_size=3, stride=2) — the pooling Spasov specify
        ("receptive field 3, stride 2"). It halves each spatial dimension
        and gives a little translation invariance. ``padding=1`` keeps the
        output-size arithmetic clean for both odd and even inputs.

      * Dropout3d — zeroes ENTIRE feature-map channels at random (not
        individual voxels). For convolutional feature maps, whose voxels
        within a channel are spatially correlated, dropping whole channels
        is the correct form of dropout regularisation.
    """

    def __init__(self, in_channels, out_channels, kernel_size,
                 stride=1, dropout=DEFAULT_DROPOUT, use_pool=True):
        super().__init__()

        # Allow an int kernel/stride to mean "the same value on all 3 axes".
        if isinstance(kernel_size, int):
            kernel_size = (kernel_size,) * 3
        if isinstance(stride, int):
            stride = (stride,) * 3

        # "Same-style" padding: half the kernel width on each side. The
        # paper does not state the padding it used, so this is a pragmatic
        # choice that keeps the spatial sizes predictable through the net.
        padding = tuple(k // 2 for k in kernel_size)

        self.conv = nn.Conv3d(in_channels, out_channels, kernel_size,
                              stride=stride, padding=padding, bias=False)
        self.norm = nn.BatchNorm3d(out_channels)
        self.act = nn.ELU(inplace=True)

        self.use_pool = use_pool
        self.pool = nn.MaxPool3d(kernel_size=3, stride=2, padding=1)
        self.drop = nn.Dropout3d(p=dropout)

    def forward(self, x):
        x = self.conv(x)        # learn spatial filters
        x = self.norm(x)        # normalise per channel
        x = self.act(x)         # non-linearity
        if self.use_pool:
            x = self.pool(x)    # downsample by 2x
        x = self.drop(x)        # channel-wise dropout
        return x


class MRI3DCNN(nn.Module):
    """
    Vanilla single-stream 3D CNN.

    Topology (channel counts and kernel sizes follow the Spasov MRI
    pathway). For a [B, 1, 182, 218, 182] input the spatial size shrinks as
    annotated on the right (D x H x W after each block):

        input                                         1 x 182 x 218 x 182
        block1  Conv 1 ->24  k(11,13,11) s4 + pool   24 x  23 x  28 x  23
        block2  Conv 24->48  k(5,6,5)    s1 + pool   48 x  12 x  15 x  12
        block3  Conv 48->96  k(3,3,3)    s1 + pool   96 x   6 x   8 x   6
        block4  Conv 96->48  k(3,4,3)    s1 + pool   48 x   3 x   5 x   3
        block5  Conv 48-> 8  k(3,4,3)    s1 + pool    8 x   2 x   3 x   2
        flatten                                              96 features
        fc1     Linear 96 -> 32  (+ BN + ELU + drop)
        head    Linear 32 -> n_outputs                       logits

    The first block does the heavy lifting on memory: a large anisotropic
    11x13x11 kernel with stride 4 aggressively downsamples the full-
    resolution 182x218x182 volume straight away, which is what keeps the
    rest of the network small enough to train on a single GPU.
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
        # Channel progression 1 -> 24 -> 48 -> 96 -> 48 -> 8. The "96 then
        # back down to 48 then 8" funnel mirrors Spasov's design: widen to
        # build a rich representation, then squeeze to a compact descriptor
        # before flattening.
        self.block1 = ConvBlock3D(in_channels, 24, (11, 13, 11),
                                  stride=4, dropout=dropout)
        self.block2 = ConvBlock3D(24, 48, (5, 6, 5),
                                  stride=1, dropout=dropout)
        self.block3 = ConvBlock3D(48, 96, (3, 3, 3),
                                  stride=1, dropout=dropout)
        self.block4 = ConvBlock3D(96, 48, (3, 4, 3),
                                  stride=1, dropout=dropout)
        self.block5 = ConvBlock3D(48, 8, (3, 4, 3),
                                  stride=1, dropout=dropout)

        # ── Classifier head ──────────────────────────────────────────────
        self.flatten = nn.Flatten()
        # LazyLinear infers the flattened feature length on the FIRST
        # forward pass, so we do not have to hard-code the post-convolution
        # spatial size (it depends on the exact input shape and the padding
        # arithmetic). After one forward pass it behaves like a normal
        # nn.Linear.
        self.fc1 = nn.LazyLinear(32)
        self.fc_norm = nn.BatchNorm1d(32)
        self.fc_act = nn.ELU(inplace=True)
        self.fc_drop = nn.Dropout(p=dropout)

        # Final linear layer producing the raw logits. No activation here:
        # the loss function (BCEWithLogitsLoss / CrossEntropyLoss) applies
        # the sigmoid / softmax internally in a numerically stable way.
        self.head = nn.Linear(32, n_outputs)

    def forward(self, x):
        """
        Args:
            x: input tensor of shape [B, in_channels, 182, 218, 182].
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
