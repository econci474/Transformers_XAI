"""
Vendored minimal subset of qasymjomart/ViT_recipe_for_AD (MIT License)
for ADNI supervised fine-tuning. Source: https://github.com/qasymjomart/ViT_recipe_for_AD

Modules:
  - vit3d:      Vision_Transformer3D + supporting blocks
  - weight_init: trunc_normal_, named_apply, init_weights_vit_timm helpers
  - checkpoint:  load_pretrained_checkpoint (reads ckpt['net'], strict=False)
"""
