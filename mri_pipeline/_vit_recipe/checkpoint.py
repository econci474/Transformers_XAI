"""
Checkpoint loading for the MAE-pretrained ViT-B/3D from ViT_recipe_for_AD.

Vendored from qasymjomart/ViT_recipe_for_AD/utils/utils.py:load_pretrained_checkpoint
(MIT License). The 'imagenet_weights/' branch is preserved but requires timm to be
installed if used; the standard .pth branch (our case) does not need timm.
"""

import torch


def load_pretrained_checkpoint(model, pre_trained_model_path, checkpoint_type=None):
    """
    Load MAE-pretrained encoder weights into a Vision_Transformer3D supervised model.

    Reads `ckpt['net']` from the upstream save format, drops keys with shape
    mismatches against the target model (head, pos_embed, patch_embed.proj),
    and uses strict=False so MAE decoder keys in the checkpoint are silently
    ignored.

    Returns the model. Prints the missing-keys list — for a fresh supervised head
    it should be ['head.weight', 'head.bias'] (and 'pos_embed' if shape mismatch).
    """
    if pre_trained_model_path == 'imagenet_weights/':
        import timm  # only required if user asks for ImageNet init
        keys_to_remove = ['head.weight', 'head.bias', 'pos_embed',
                          'patch_embed.proj.weight', 'patch_embed.proj.bias']
        checkpoint_model = timm.create_model('vit_base_patch16_224').state_dict()
        print('Loaded ImageNet pre-trained checkpoint')
    else:
        checkpoint = torch.load(pre_trained_model_path, map_location='cpu')
        print(f"Loaded pre-trained checkpoint from: {pre_trained_model_path}")
        if 'net' not in checkpoint:
            raise KeyError(
                f"Expected key 'net' in checkpoint. Found keys: {list(checkpoint.keys())}"
            )
        checkpoint_model = checkpoint['net']
        keys_to_remove = ['head.weight', 'head.bias', 'pos_embed',
                          'patch_embed.proj.weight', 'patch_embed.proj.bias']

    state_dict = model.state_dict()
    for k in keys_to_remove:
        if k in checkpoint_model and k in state_dict and checkpoint_model[k].shape != state_dict[k].shape:
            print(f"  Removing key {k} from pretrained checkpoint (shape mismatch)")
            del checkpoint_model[k]

    msg = model.load_state_dict(checkpoint_model, strict=False)
    print(f"  Missing keys (will be randomly initialised): {msg.missing_keys}")
    if msg.unexpected_keys:
        print(f"  Unexpected keys (ignored — likely MAE decoder): {len(msg.unexpected_keys)} keys")

    return model
