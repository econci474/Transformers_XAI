"""
3D Vision Transformer for ADNI MRI fine-tuning.

Vendored from qasymjomart/ViT_recipe_for_AD/models/vit3d.py (MIT License).
Trimmed to the supervised-classification subset:
  - DeiT_Transformer3D removed (we don't use distillation)
  - get_3d_sincos_pos_embed dependency removed (only used when pos_embed_type='abs';
    we use the default 'learnable')
  - torchvision import removed (was unused upstream)
"""

import math

import numpy as np
import torch
import torch.nn as nn

from .weight_init import (
    trunc_normal_,
    init_weights_vit_timm,
    get_init_weights_vit,
    named_apply,
)


def drop_path(x, drop_prob: float = 0., training: bool = False):
    if drop_prob == 0. or not training:
        return x
    keep_prob = 1 - drop_prob
    shape = (x.shape[0],) + (1,) * (x.ndim - 1)
    random_tensor = keep_prob + torch.rand(shape, dtype=x.dtype, device=x.device)
    random_tensor.floor_()
    output = x.div(keep_prob) * random_tensor
    return output


class DropPath(nn.Module):
    def __init__(self, drop_prob=None):
        super().__init__()
        self.drop_prob = drop_prob

    def forward(self, x):
        return drop_path(x, self.drop_prob, self.training)


class PatchEmbed3D(nn.Module):
    def __init__(self, img_size, patch_size, embed_dim=768, patch_embed_fun='conv3d'):
        super().__init__()
        self.img_size = img_size
        self.patch_size = patch_size
        self.n_patches = (img_size[0] // patch_size) * (img_size[1] // patch_size) * (img_size[2] // patch_size)

        sample_torch = torch.rand((1, 1, *self.img_size))

        if patch_embed_fun == 'conv3d':
            self.proj = nn.Conv3d(
                in_channels=1,
                out_channels=embed_dim,
                kernel_size=patch_size,
                stride=patch_size,
            )

        out = self.proj(sample_torch)
        self.n_patches = out.flatten(2).shape[2]

    def forward(self, x):
        x = self.proj(x)
        x = x.flatten(2)
        x = x.transpose(1, 2)
        return x


class Attention(nn.Module):
    def __init__(self, dim, n_heads=12, qkv_bias=True, attn_p=0., proj_p=0.):
        super().__init__()
        self.dim = dim
        self.n_heads = n_heads
        self.head_dim = dim // n_heads
        self.scale = self.head_dim ** -0.5

        self.qkv = nn.Linear(dim, dim * 3, bias=qkv_bias)
        self.attn_drop = nn.Dropout(attn_p)
        self.proj = nn.Linear(dim, dim)
        self.proj_drop = nn.Dropout(proj_p)

    def forward(self, x):
        n_samples, n_tokens, dim = x.shape
        if dim != self.dim:
            raise ValueError

        qkv = self.qkv(x)
        qkv = qkv.reshape(n_samples, n_tokens, 3, self.n_heads, self.head_dim)
        qkv = qkv.permute(2, 0, 3, 1, 4)
        q, k, v = qkv[0], qkv[1], qkv[2]

        k_t = k.transpose(-2, -1)
        dp = (q @ k_t) * self.scale
        attn = dp.softmax(dim=-1)
        attn = self.attn_drop(attn)

        weighted_avg = attn @ v
        weighted_avg = weighted_avg.transpose(1, 2).flatten(2)
        x = self.proj(weighted_avg)
        x = self.proj_drop(x)
        return x


class MLP(nn.Module):
    def __init__(self, in_features, hidden_features, out_features, p=0.):
        super().__init__()
        self.fc1 = nn.Linear(in_features, hidden_features)
        self.act = nn.GELU()
        self.fc2 = nn.Linear(hidden_features, out_features)
        self.drop = nn.Dropout(p)

    def forward(self, x):
        x = self.fc1(x)
        x = self.act(x)
        x = self.drop(x)
        x = self.fc2(x)
        x = self.drop(x)
        return x


class Block(nn.Module):
    def __init__(self, dim, n_heads, mlp_ratio=4.0, qkv_bias=True, drop_path=0., p=0., attn_p=0.):
        super().__init__()
        self.norm1 = nn.LayerNorm(dim, eps=1e-6)
        self.attn = Attention(
            dim,
            n_heads=n_heads,
            qkv_bias=qkv_bias,
            attn_p=attn_p,
            proj_p=p,
        )
        self.norm2 = nn.LayerNorm(dim, eps=1e-6)
        hidden_features = int(dim * mlp_ratio)
        self.mlp = MLP(
            in_features=dim,
            hidden_features=hidden_features,
            out_features=dim,
        )
        self.drop_path = DropPath(drop_prob=drop_path) if drop_path > 0. else nn.Identity()

    def forward(self, x):
        x = x + self.drop_path(self.attn(self.norm1(x)))
        x = x + self.drop_path(self.mlp(self.norm2(x)))
        return x


class Vision_Transformer3D(nn.Module):
    """
    3D Vision Transformer matching the upstream ViT_recipe_for_AD architecture.

    For supervised fine-tuning from the MAE-pretrained checkpoint, use:
        img_size=(128,128,128), patch_size=16, in_chans=1, embed_dim=768,
        depth=12, n_heads=12, mlp_ratio=4.0, qkv_bias=True, drop_path_rate=0.1,
        p=0.0, attn_p=0.1, global_avg_pool=False, pos_embed_type='learnable',
        patch_embed_fun='conv3d'.
    """

    def __init__(self,
                 img_size=(128, 128, 128),
                 patch_size=16,
                 in_chans=1,
                 n_classes=2,
                 embed_dim=768,
                 depth=12,
                 n_heads=12,
                 mlp_ratio=4.,
                 qkv_bias=True,
                 drop_path_rate=0.,
                 p=0.,
                 attn_p=0.,
                 patch_embed_fun='conv3d',
                 weight_init='',
                 global_avg_pool=False,
                 pos_embed_type='learnable'):
        super().__init__()

        if patch_embed_fun != 'conv3d':
            raise ValueError(f"Vendored vit3d.py supports patch_embed_fun='conv3d' only "
                             f"(got {patch_embed_fun!r}).")
        if pos_embed_type != 'learnable':
            raise ValueError(f"Vendored vit3d.py supports pos_embed_type='learnable' only "
                             f"(got {pos_embed_type!r} — 'abs' would require get_3d_sincos_pos_embed).")

        self.patch_embed = PatchEmbed3D(
            img_size=img_size,
            patch_size=patch_size,
            embed_dim=embed_dim,
            patch_embed_fun=patch_embed_fun,
        )

        self.cls_token = nn.Parameter(torch.rand(1, 1, embed_dim)) if not global_avg_pool else None
        embed_len = self.patch_embed.n_patches if global_avg_pool else 1 + self.patch_embed.n_patches
        self.pos_embed = nn.Parameter(torch.rand(1, embed_len, embed_dim), requires_grad=True)

        self.pos_drop = nn.Dropout(p=p)
        self.dpr = [x.item() for x in torch.linspace(0, drop_path_rate, depth)]

        self.blocks = nn.ModuleList([
            Block(
                dim=embed_dim,
                n_heads=n_heads,
                mlp_ratio=mlp_ratio,
                qkv_bias=qkv_bias,
                drop_path=self.dpr[ii],
                p=p,
                attn_p=attn_p,
            )
            for ii in range(depth)
        ])

        self.norm = nn.LayerNorm(embed_dim, eps=1e-6)
        self.head = nn.Linear(embed_dim, n_classes)

        if self.cls_token is not None:
            trunc_normal_(self.cls_token, std=.02)

        self.init_weights(weight_init)

    def init_weights(self, mode=''):
        assert mode in ('jax', 'jax_nlhb', 'moco', '')
        head_bias = 0.
        trunc_normal_(self.pos_embed, std=.02)
        if self.cls_token is not None:
            nn.init.normal_(self.cls_token, std=1e-6)
        named_apply(get_init_weights_vit(mode, head_bias), self)

    @torch.jit.ignore
    def no_weight_decay(self):
        return {'pos_embed', 'cls_token'}

    @torch.jit.ignore
    def get_classifier(self):
        return self.head

    def forward(self, x):
        n_samples = x.shape[0]
        x = self.patch_embed(x)

        if self.cls_token is not None:
            cls_token = self.cls_token.expand(n_samples, -1, -1)
            x = torch.cat((cls_token, x), dim=1)
        x = x + self.pos_embed
        x = self.pos_drop(x)

        for block in self.blocks:
            x = block(x)

        x = self.norm(x)
        cls_token_final = x[:, 0] if self.cls_token is not None else x.mean(dim=1)
        x = self.head(cls_token_final)
        return x
