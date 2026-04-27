import os
import torch
from monai.networks.nets import ViT
from safetensors.torch import load_file
import nibabel as nib
import numpy as np

# ── Paths ─────────────────────────────────────────────────────────────────────
# Weights downloaded from: https://huggingface.co/eugenehp/brainiac/tree/main
WEIGHTS_DIR      = r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\hugging_face\brainiac"
BACKBONE_WEIGHTS = os.path.join(WEIGHTS_DIR, "backbone.safetensors")  # pretrained ViT encoder only
MCI_WEIGHTS      = os.path.join(WEIGHTS_DIR, "mci.safetensors")       # backbone + MCI classification head

# Local path to preprocessed MRI (space-MNI96, 96x96x96 voxels, single channel)
PREPROCESSED_MRI = r'D:\ADNI_BIDS_project\derivatives\brainiac_inputs\sub-002S0413\sub-002S0413_space-MNI96_desc-brainiac_T1w.nii.gz'

# MCI classification labels (binary: 0 = CN/other, 1 = MCI)
LABELS = ["CN", "MCI"]


# ── Load NIfTI → tensor ───────────────────────────────────────────────────────
def load_nifti_as_tensor(path: str) -> torch.Tensor:
    """
    Load a preprocessed NIfTI file and return a (1, 1, 96, 96, 96) float32 tensor.
    The image is expected to already be registered to MNI space and sized 96³.
    """
    img = nib.load(path)
    data = img.get_fdata(dtype=np.float32)            # (96, 96, 96)
    tensor = torch.from_numpy(data)
    tensor = tensor.unsqueeze(0).unsqueeze(0)         # → (1, 1, 96, 96, 96)
    return tensor


# ── Shared ViT config (BRAINIAC backbone) ─────────────────────────────────────
VIT_KWARGS = dict(
    in_channels=1,
    img_size=(96, 96, 96),
    patch_size=(16, 16, 16),
    hidden_size=768,
    mlp_dim=3072,
    num_layers=12,
    num_heads=12,
)


# ── Part 1: Feature extraction with backbone only ─────────────────────────────
# features[0][:, 0] → 768-dim embedding of the first patch token.
# BRAINIAC does not use a CLS token; it uses features[:, 0] as the global rep.
print("=== BRAINIAC Feature Extraction (backbone) ===")
backbone = ViT(**VIT_KWARGS)
backbone_weights = load_file(BACKBONE_WEIGHTS)
backbone.load_state_dict(backbone_weights, strict=False)
backbone.eval()

volume = load_nifti_as_tensor(PREPROCESSED_MRI)
print(f"Input tensor shape: {tuple(volume.shape)}")   # (1, 1, 96, 96, 96)

with torch.no_grad():
    features, _ = backbone(volume)          # features: (1, 216, 768)
    embedding = features[:, 0]              # first patch token → (1, 768)

print(f"Feature vector shape:          {tuple(embedding.shape)}")
print(f"Feature vector (first 8 dims): {embedding[0, :8].numpy()}")


# ── Part 2: MCI binary classification ─────────────────────────────────────────
# BRAINIAC's mci.safetensors stores:
#   - All ViT backbone weights (same keys as backbone.safetensors)
#   - head.fc.weight / head.fc.bias  →  a Linear(768, 2) applied to features[:, 0]
# It does NOT use MONAI's built-in classification=True head, so we reconstruct
# the architecture manually by separating the two parts of the checkpoint.
print("\n=== BRAINIAC MCI Classification ===")
mci_weights = load_file(MCI_WEIGHTS)

# Separate backbone weights from the classification head weights
backbone_w = {k: v for k, v in mci_weights.items() if not k.startswith("head.")}
head_w     = {k: v for k, v in mci_weights.items() if k.startswith("head.")}

# Load backbone
clf_backbone = ViT(**VIT_KWARGS)
clf_backbone.load_state_dict(backbone_w, strict=False)
clf_backbone.eval()

# Build classification head from the actual weight shape in the checkpoint.
# head.fc.weight has shape (num_outputs, 768).
# BRAINIAC uses num_outputs=1 → single sigmoid logit (binary: 0=CN, 1=MCI).
num_outputs = head_w["head.fc.weight"].shape[0]
head_fc = torch.nn.Linear(768, num_outputs)
head_fc.weight = torch.nn.Parameter(head_w["head.fc.weight"])
head_fc.bias   = torch.nn.Parameter(head_w["head.fc.bias"])
head_fc.eval()

with torch.no_grad():
    features, _ = clf_backbone(volume)      # (1, 216, 768)
    embedding   = features[:, 0]            # first patch token → (1, 768)
    logit       = head_fc(embedding)        # (1, num_outputs)

if num_outputs == 1:
    # Single sigmoid logit: P(MCI) = sigmoid(logit)
    prob_mci = torch.sigmoid(logit[0, 0]).item()
    prob_cn  = 1.0 - prob_mci
    pred_label = "MCI" if prob_mci >= 0.5 else "CN"
else:
    # Two-class softmax (fallback)
    probs    = torch.softmax(logit, dim=1)
    prob_cn  = probs[0, 0].item()
    prob_mci = probs[0, 1].item()
    pred_label = LABELS[probs.argmax(dim=1).item()]

print(f"Raw logit:        {logit[0].numpy()}")
print(f"P(CN):            {prob_cn:.4f}")
print(f"P(MCI):           {prob_mci:.4f}")
print(f"Predicted label:  {pred_label}")