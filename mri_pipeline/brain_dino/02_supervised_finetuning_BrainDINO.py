"""
02_supervised_finetuning_BrainDINO.py
======================================
Supervised fine-tuning of the BrainDINO ViT-B/16 encoder (arXiv:2604.27277)
on ADNI T1w MRIs preprocessed by 01_prepare_braindino_inputs.py
((256, 256, 128) in RAS axis order = 128 axial x 256 cor x 256 sag).

Architecture (paper Section 4.3):
    For each volume of shape (R, A, S=axial):
        permute to (axial, A, R)              -> (128, 256, 256)
        for each axial slice s in [0, 128):
            replicate 1ch -> 3ch              -> (3, 256, 256)
            encoder(s)                        -> CLS token (768-d)
        mean-pool CLS tokens across slices    -> (768,)
        head(pooled)                          -> (num_classes,)

Head + training config reconciled against upstream
github.com/mclwu22/BrainDINO/blob/main/runners/classification.py:
    head            LayerNorm -> Linear(768->768) -> GELU -> Dropout(0.2)
                    -> Linear(768->n_classes)
    optimizer       Adam (not AdamW), weight_decay=1e-5
    scheduler       ReduceLROnPlateau (mode=max, factor=0.5, patience=10),
                    stepped on val balanced-accuracy each epoch -- replaces
                    the cosine+warmup schedule used elsewhere in the
                    pipeline. No warmup.
    grad_accum      1 (paper)
    drop_path_rate  0.2 (encoder kwarg; head dropout = 0.2 = teacher_dropout_prob)
    label_smoothing 0.0 (paper default; bump to 0.1 if T1b/T2 collapse)

Deviation from paper (documented): use_amp=True with bf16 (paper uses
fp32). Needed for A100-80GB memory at batch=2, 128 axial slices, 256x256,
ViT-B; running fp32 would force batch=1 + grad_accum>=8 or OOM. Effect on
results is expected to be sub-percentage-point.

Strategies:
    --strategy frozen      Encoder requires_grad=False, only head trains.
                           LR=1e-3, epochs=100, patience=50.
    --strategy full_ft     Encoder trains end-to-end. LR=1e-4, epochs=50.
    --strategy lora        Paper-canonical: backbone frozen + LoRA-rank-8
                           adapters injected into every attention layer's
                           qkv + proj. Uses peft. LR=1e-4, epochs=8.

Augmentation arms (--augment):
    none           No augmentation (deterministic).
    stochastic     RandFlip (R/L axis) + RandAffine + intensity jitter.
    plus_original  ConcatDataset(orig, K x always-on stochastic).

Output schema: same metrics.json as BrainMVP (config, test_metrics,
test_diagnostics) so the aggregator picks it up unchanged.
"""

import argparse
import importlib.util
import json
import os
import sys
import warnings
from datetime import datetime
from pathlib import Path

# torch FIRST on Windows to dodge MKL DLL collision (see prep script).
import torch  # noqa: F401
import numpy as np
import pandas as pd
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import balanced_accuracy_score, confusion_matrix
from sklearn.utils.class_weight import compute_class_weight
from torch.utils.data import DataLoader, ConcatDataset

import monai.transforms as mt
from monai.data import Dataset as MonaiDataset

# ── Path setup ────────────────────────────────────────────────────────────────
THIS_DIR = Path(__file__).resolve().parent
MRI_DIR = THIS_DIR.parent
REPO_ROOT = MRI_DIR.parent
sys.path.insert(0, str(MRI_DIR))     # for the importlib import of ViT pipeline
sys.path.insert(0, str(REPO_ROOT))   # for bidsification.exclusions
sys.path.insert(0, str(THIS_DIR))    # for the vendored dinov3 package

# Vendored ViT (loads paper-fidelity ViT-B with storage tokens + RoPE).
from dinov3 import vit_base

# Shared MRI-pipeline utilities (TASK_CONFIG, splits loader, metrics, etc).
# Mirrors the BrainMVP pattern so the BrainDINO trainer can drop into the
# same harness without re-implementing label resolution, NaN-aware early
# stop, balanced-accuracy checkpoint selection, etc.
_vit_spec = importlib.util.spec_from_file_location(
    "_vit_pipeline", MRI_DIR / "04_supervised_finetuning_ViT.py")
_vit = importlib.util.module_from_spec(_vit_spec)
sys.modules["_vit_pipeline"] = _vit
_vit_spec.loader.exec_module(_vit)

TASK_CONFIG          = _vit.TASK_CONFIG
cohort_label         = _vit.cohort_label
class_names_for      = _vit.class_names_for
load_matched_labels  = _vit.load_matched_labels
session_to_months    = _vit.session_to_months
patient_id_to_bids_sub = _vit.patient_id_to_bids_sub
compute_test_metrics = _vit.compute_test_metrics
compute_diagnostics  = _vit.compute_diagnostics
run_one_epoch        = _vit.run_one_epoch
lr_multiplier_at     = _vit.lr_multiplier_at

warnings.filterwarnings("ignore")

MODEL_SLUG = "BrainDINO_vitb16"

# Paper-fidelity ViT-B/16 kwargs (must match the upstream checkpoint).
VIT_KWARGS = dict(
    drop_path_rate=0.2,
    layerscale_init=1.0e-5,
    n_storage_tokens=4,
    qkv_bias=False,
    mask_k_bias=True,
)
N_AXIAL = 128            # paper Section 4.3 explicit
EMBED_DIM = 768          # vit_base
PATCH = 16
SLICE_SIZE = 256         # 256x256 axial slice fed to ViT (no resize needed)


# ── Init wandb (local copy: no ViT-only args) ─────────────────────────────────
def init_wandb(args, task_cfg, extra=None):
    if not args.wandb:
        return None
    try:
        import wandb
    except ImportError:
        print("  [WARN] --wandb set but wandb not installed; continuing local-only.")
        return None
    name = args.wandb_run_name or (
        f"{args.task}-s{args.seed}-{args.strategy}-{args.augment}")
    cfg = {
        "task":            args.task,
        "seed":            args.seed,
        "strategy":        args.strategy,
        "augment":         args.augment,
        "aug_copies":      args.aug_copies,
        "epochs":          args.epochs,
        "lr":              args.lr,
        "weight_decay":    args.weight_decay,
        "batch_size":      args.batch_size,
        "warmup_epochs":   args.warmup_epochs,
        "patience":        args.patience,
        "drop_rate":       args.drop_rate,
        "label_smoothing": args.label_smoothing,
        "lora_rank":       args.lora_rank if args.strategy == "lora" else None,
        "grad_accum":      args.grad_accum_steps,
        "long_mode":       args.long_mode,
        "max_months":      args.max_months,
        "session_policy":  task_cfg["session_policy"],
        "num_labels":      task_cfg["num_labels"],
        "pretrained_ckpt": args.pretrained_ckpt,
    }
    if extra:
        cfg.update(extra)
    return wandb.init(project=args.wandb_project, entity=args.wandb_entity,
                      name=name, config=cfg, reinit=True)


# ── Data transforms (volume-level only; slice extraction is in forward) ───────
def _bdino_eval_transform() -> mt.Compose:
    return mt.Compose([
        mt.LoadImaged(keys=["image"]),
        mt.EnsureChannelFirstd(keys=["image"]),
        mt.ToTensord(keys=["image"]),
    ])


def _bdino_stochastic_transform(p=0.5) -> mt.Compose:
    """Augment the (R, A, S) volume BEFORE slicing -- spatial aug applies to
    all slices uniformly, intensity jitter is per-volume.

    Affine rotation/scale ranges match the user's established preference
    for paper-faithful pretrained-model finetuning (small magnitudes that
    won't distort hippocampal landmarks)."""
    return mt.Compose([
        mt.LoadImaged(keys=["image"]),
        mt.EnsureChannelFirstd(keys=["image"]),
        # Volume-level spatial aug. Axis 0 = R, so axis=0 flip is the
        # left/right anatomical flip (the only one used in MRI aug).
        mt.RandFlipd(keys=["image"], prob=p, spatial_axis=0),
        mt.RandAffined(keys=["image"], prob=p,
                       rotate_range=(0.087,) * 3,   # ~+/- 5 degrees
                       scale_range=(0.05,) * 3,
                       padding_mode="zeros"),
        mt.RandScaleIntensityd(keys=["image"], factors=0.1, prob=0.5),
        mt.RandShiftIntensityd(keys=["image"], offsets=0.1, prob=0.5),
        mt.ToTensord(keys=["image"]),
    ])


def build_bdino_dataset(df, train: bool,
                        augment: str = "stochastic", aug_copies: int = 1):
    items = [{"image": row.image_path, "label": int(row.label)}
             for row in df.itertuples()]

    if not train or augment == "none":
        return MonaiDataset(data=items, transform=_bdino_eval_transform())
    if augment == "stochastic":
        return MonaiDataset(data=items, transform=_bdino_stochastic_transform())
    if augment == "plus_original":
        k = max(1, int(aug_copies))
        orig = MonaiDataset(data=items, transform=_bdino_eval_transform())
        aug = [MonaiDataset(data=items,
                            transform=_bdino_stochastic_transform(p=1.0))
               for _ in range(k)]
        return ConcatDataset([orig] + aug)
    raise ValueError(f"Unknown augment mode: {augment!r}")


def _resolve_bdino_path(patient_id: str, viscode: str, inputs_dir: Path) -> str:
    sub = patient_id_to_bids_sub(patient_id)
    ses_label = f"ses-{viscode}"
    return str(inputs_dir / sub / ses_label /
               f"{sub}_{ses_label}_space-MNI128x256_desc-braindino_T1w.nii.gz")


def load_split_bdino(baseline_csv, matched_df, task_cfg, inputs_dir,
                     session_filter=None, max_months=None):
    """Same load_split logic as BrainMVP but with BrainDINO file paths."""
    from bidsification.exclusions import is_excluded_subject

    subjects = pd.read_csv(baseline_csv, usecols=["Patient_ID"])["Patient_ID"].unique().tolist()
    df = matched_df[matched_df["Patient_ID"].isin(subjects)].copy()

    if session_filter is not None:
        df = df[df["bids_ses"] == session_filter].copy()
    elif max_months is not None:
        df["_months"] = df["bids_ses"].map(session_to_months)
        df = df.dropna(subset=["_months"])
        df = df[df["_months"] <= max_months].drop(columns=["_months"])

    if task_cfg["session_policy"] == "baseline_anchored":
        cap = int(task_cfg.get("label_anchor_max_months", 0))
        df["_months"] = df["bids_ses"].map(session_to_months)
        df = df.dropna(subset=["_months"])
        df = df[df["_months"] <= cap].drop(columns=["_months"]).copy()
    elif task_cfg["session_policy"] == "baseline_only":
        df = df[df["bids_ses"] == "bl"].copy()

    if task_cfg["filter_non_ad"]:
        df = df[df["Label_bl_multi"] != "AD"].copy()

    df = df[~df["Patient_ID"].apply(is_excluded_subject)].copy()
    df = _vit._resolve_labels(df, task_cfg)

    df["image_path"] = df.apply(
        lambda r: _resolve_bdino_path(r["Patient_ID"], r["bids_ses"], inputs_dir),
        axis=1)
    have = df["image_path"].apply(os.path.isfile)
    n_missing = int((~have).sum())
    df = df[have].copy()

    out_cols = ["Patient_ID", "bids_ses", "label", "image_path"]
    return df.reset_index(drop=True)[out_cols], n_missing


# ── BrainDINO classifier wrapper ──────────────────────────────────────────────
class BrainDINOClassifier(nn.Module):
    """Slice-encode (BrainDINO ViT-B/16) -> mean-pool CLS tokens -> 2-layer
    MLP head matching the paper's downstream classification head:

        head = LayerNorm(768) -> Linear(768, 768) -> GELU
               -> Dropout(drop_rate) -> Linear(768, n_classes)

    drop_rate=0.2 corresponds to the paper's `teacher_dropout_prob=0.2`
    in runners/classification.py.

    Input volume shape: (B, 1, R=256, A=256, S=128) -- the C.2 preprocessing
    output, in RAS axis order. The forward permutes axial first so we can
    iterate over slices along dim 1.
    """

    def __init__(self, n_classes: int, drop_rate: float = 0.2):
        super().__init__()
        # img_size=256 sets the PatchEmbed grid to 256/16=16x16 patches. The
        # ViT uses RoPE position embedding which is shape-agnostic, so even
        # an img_size mismatch wouldn't break it at runtime; we pass the
        # actual slice size for clarity.
        self.encoder = vit_base(img_size=SLICE_SIZE, **VIT_KWARGS)
        # 2-layer MLP head (paper runners/classification.py + slice_student).
        # Dropout sits between the inner activation and the classifier --
        # i.e. teacher_dropout_prob is applied to the pooled representation
        # after the inner non-linearity.
        self.head = nn.Sequential(
            nn.LayerNorm(EMBED_DIM),
            nn.Linear(EMBED_DIM, EMBED_DIM),
            nn.GELU(),
            nn.Dropout(drop_rate),
            nn.Linear(EMBED_DIM, n_classes),
        )

    def encode_slices(self, slices: torch.Tensor) -> torch.Tensor:
        """slices: (M, 1, 256, 256) -> CLS tokens (M, 768)."""
        # 3-ch fan-out: ViT-B/16 was pretrained with in_chans=3.
        slices_3 = slices.repeat(1, 3, 1, 1)
        feats = self.encoder(slices_3, is_training=True)
        return feats["x_norm_clstoken"]            # (M, 768)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # x: (B, 1, 256, 256, 128) -- (B, C, R, A, S=axial). Permute axial first.
        if x.shape[-1] != N_AXIAL:
            raise RuntimeError(
                f"Expected (B, 1, 256, 256, {N_AXIAL}) input, got {tuple(x.shape)}. "
                f"Was the volume preprocessed with 01_prepare_braindino_inputs.py?")
        B = x.shape[0]
        # (B, 1, R, A, S) -> (B, S, A, R) -> (B*S, 1, A, R)
        slices = x.squeeze(1).permute(0, 3, 2, 1).contiguous()  # (B, S, A, R)
        slices = slices.reshape(B * N_AXIAL, 1, SLICE_SIZE, SLICE_SIZE)
        cls = self.encode_slices(slices)                        # (B*S, 768)
        cls = cls.view(B, N_AXIAL, EMBED_DIM)                    # (B, S, 768)
        pooled = cls.mean(dim=1)                                 # (B, 768)
        return self.head(pooled)                                 # (B, n_classes)


def load_braindino_pretrained(model: BrainDINOClassifier, ckpt_path: str):
    """Load BrainDINO teacher weights into the ViT encoder.

    The checkpoint is the upstream DINO teacher dict; we strip the
    `backbone.` prefix and drop any ibot / dino_head keys -- matches
    upstream braindino_trainer.py verbatim.
    """
    print(f"  Loading BrainDINO checkpoint: {ckpt_path}")
    ck = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    if "teacher" not in ck:
        raise RuntimeError(
            f"Expected 'teacher' key in BrainDINO checkpoint; "
            f"got keys: {list(ck.keys())[:5]}")
    sd = {k.replace("backbone.", ""): v
          for k, v in ck["teacher"].items()
          if "ibot" not in k and "dino_head" not in k}
    missing, unexpected = model.encoder.load_state_dict(sd, strict=False)
    print(f"  Loaded {len(sd) - len(unexpected)} encoder tensors  "
          f"(missing={len(missing)}, unexpected={len(unexpected)})")
    if missing:
        print(f"    first 3 missing: {missing[:3]}")
    if unexpected:
        print(f"    first 3 unexpected: {unexpected[:3]}")
    return model


# ── LoRA wrapping (paper-canonical strategy) ──────────────────────────────────
def apply_lora_to_encoder(model: BrainDINOClassifier, rank: int = 8):
    """Inject LoRA-rank adapters into every attention layer (qkv + proj)
    of the encoder. Uses peft. Freezes the backbone; only LoRA adapters
    + head are trainable."""
    try:
        from peft import LoraConfig, get_peft_model
    except ImportError as e:
        raise RuntimeError(
            "--strategy lora requires the `peft` package. Install with:\n"
            "    pip install peft>=0.7") from e

    # Freeze all encoder params first; LoRA adapters get added inside the
    # attention layers and will be the only trainable parameters in the
    # encoder. The head stays unfrozen because it's outside the encoder.
    for p in model.encoder.parameters():
        p.requires_grad = False

    lora_cfg = LoraConfig(
        r=rank, lora_alpha=rank * 2, lora_dropout=0.05,
        # peft's target_modules accepts substrings; the dinov3 attention
        # modules have qkv + proj nn.Linear submodules.
        target_modules=["qkv", "proj"],
        bias="none",
        # The encoder is a plain nn.Module (not a HF task model) so we
        # use the generic FEATURE_EXTRACTION task type which doesn't add a
        # task-specific head.
        task_type=None,
    )
    model.encoder = get_peft_model(model.encoder, lora_cfg)
    n_train = sum(p.numel() for p in model.encoder.parameters() if p.requires_grad)
    print(f"  LoRA-{rank} adapters injected; encoder trainable params: {n_train:,}")
    return model


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="BrainDINO ViT-B/16 fine-tuning on ADNI MRI")
    p.add_argument("--task", type=str, required=True,
                   choices=list(TASK_CONFIG.keys()))
    p.add_argument("--seed", type=int, required=True, choices=[0, 1, 2])
    p.add_argument("--strategy", type=str, default="lora",
                   choices=["frozen", "full_ft", "lora"])
    p.add_argument("--lora_rank", type=int, default=8,
                   help="LoRA rank for --strategy lora (paper default: 8).")
    p.add_argument("--pretrained_ckpt", type=str, required=True,
                   help="Path to brain_dino_model.pth (DINO teacher dict).")
    session_group = p.add_mutually_exclusive_group()
    session_group.add_argument("--session", type=str, default="bl")
    session_group.add_argument("--long", type=str, default=None)
    p.add_argument("--matched_labels_csv", type=str,
                   default=r"D:\ADNI_BIDS_project\derivatives\mri_clinical_matched"
                           r"\master_mri_clinical_matched_labels.csv")
    p.add_argument("--data_dir", type=str,
                   default=r"D:\ADNI_BIDS_project\derivatives\clinical\verbose\baseline")
    p.add_argument("--braindino_inputs_dir", type=str,
                   default=r"D:\ADNI_BIDS_project\derivatives\braindino_inputs")
    p.add_argument("--out_dir", type=str, required=True)
    p.add_argument("--epochs", type=int, default=None)
    p.add_argument("--batch_size", type=int, default=2,
                   help="Volumes per batch. 128 axial slices x 256^2 x ViT-B "
                        "is ~9 GB AMP / 18 GB fp32 per volume on A100-80GB; "
                        "batch=2 fits comfortably with AMP.")
    p.add_argument("--grad_accum_steps", type=int, default=1,
                   help="Paper default. Effective batch = batch_size * grad_accum_steps.")
    p.add_argument("--lr", type=float, default=None)
    p.add_argument("--weight_decay", type=float, default=1e-5,
                   help="Paper default (runners/classification.py): 1e-5.")
    p.add_argument("--warmup_epochs", type=int, default=0,
                   help="Paper has no warmup; lr_multiplier_at returns 1.0.")
    p.add_argument("--scheduler_patience", type=int, default=10,
                   help="ReduceLROnPlateau patience (paper default: 10).")
    p.add_argument("--scheduler_factor", type=float, default=0.5,
                   help="ReduceLROnPlateau factor (paper default: 0.5).")
    p.add_argument("--patience", type=int, default=25,
                   help="Early-stopping patience (epochs without val_bacc "
                        "improvement). Distinct from --scheduler_patience.")
    p.add_argument("--num_workers", type=int, default=0)
    p.add_argument("--max_subjects", type=int, default=None)
    p.add_argument("--augment", type=str, default="stochastic",
                   choices=["none", "stochastic", "plus_original"])
    p.add_argument("--aug_copies", type=int, default=1)
    p.add_argument("--drop_rate", type=float, default=0.2,
                   help="Head dropout = paper's teacher_dropout_prob (0.2).")
    p.add_argument("--label_smoothing", type=float, default=0.0,
                   help="Paper default (runners/classification.py): 0.0. "
                        "Bump to 0.1 if T1b/T2 class collapse appears.")
    p.add_argument("--no_resume", action="store_true")
    p.add_argument("--wandb", action="store_true")
    p.add_argument("--wandb_project", type=str, default="braindino-mri-finetune")
    p.add_argument("--wandb_entity", type=str, default=None)
    p.add_argument("--wandb_run_name", type=str, default=None)
    args = p.parse_args()

    # Strategy-aware defaults
    if args.epochs is None:
        args.epochs = {"frozen": 100, "full_ft": 50, "lora": 8}[args.strategy]
    if args.lr is None:
        # Paper (runners/classification.py) sets lr=1e-4 for the frozen
        # variant -- inherited here for the head-only run too (no LoRA).
        args.lr = {"frozen": 1e-4, "full_ft": 1e-4, "lora": 1e-4}[args.strategy]

    # Resolve session selection mode
    if args.long is None:
        args.long_mode = None
        args.max_months = None
    else:
        if args.long.lower() == "all" or args.long == "0":
            args.long_mode = "all"; args.max_months = None
        else:
            n = int(args.long); args.long_mode = "cutoff"; args.max_months = n * 12
        args.session = None

    # Compat attrs for shared utilities
    args.drop_path_rate = VIT_KWARGS["drop_path_rate"]
    args.attn_dropout = 0.0
    return args


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    task_cfg = TASK_CONFIG[args.task]
    torch.manual_seed(args.seed); np.random.seed(args.seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    out_dir = (Path(args.out_dir) / f"aug_{args.augment}"
               / f"{MODEL_SLUG}_{args.strategy}" / args.task
               / f"seed_{args.seed}")
    out_dir.mkdir(parents=True, exist_ok=True)

    if (out_dir / "metrics.json").exists():
        print(f"  metrics.json already exists at {out_dir} -- skipping.")
        return

    print("=" * 70)
    print(f"  BrainDINO ViT-B/16 -- {args.task} | seed={args.seed} | "
          f"{args.strategy} | aug={args.augment}")
    print(f"  Device: {device}")
    print(f"  Output: {out_dir}")
    print("=" * 70)

    # ── Splits ────────────────────────────────────────────────────────────────
    seed_dir = Path(args.data_dir) / f"seed_{args.seed}"
    dino_dir = Path(args.braindino_inputs_dir)
    matched_df = load_matched_labels(args.matched_labels_csv)

    common = dict(matched_df=matched_df, task_cfg=task_cfg, inputs_dir=dino_dir)
    if args.long_mode is None:
        sf = args.session
        train_df, n_miss_tr = load_split_bdino(seed_dir / "train.csv", session_filter=sf, **common)
        val_df,   n_miss_va = load_split_bdino(seed_dir / "val.csv",   session_filter=sf, **common)
        test_df,  n_miss_te = load_split_bdino(seed_dir / "test.csv",  session_filter=sf, **common)
    else:
        mm = args.max_months
        train_df, n_miss_tr = load_split_bdino(seed_dir / "train.csv", max_months=mm, **common)
        val_df,   n_miss_va = load_split_bdino(seed_dir / "val.csv",   max_months=mm, **common)
        test_df,  n_miss_te = load_split_bdino(seed_dir / "test.csv",  max_months=mm, **common)

    if args.max_subjects:
        train_df = train_df.head(args.max_subjects)
        val_df   = val_df.head(min(args.max_subjects, len(val_df)))
        test_df  = test_df.head(min(args.max_subjects, len(test_df)))

    print(f"  Splits -- train: {len(train_df)}  val: {len(val_df)}  test: {len(test_df)}")
    print(f"  Missing files -- train: {n_miss_tr}  val: {n_miss_va}  test: {n_miss_te}")
    manifest = pd.concat([train_df.assign(split="train"),
                          val_df.assign(split="val"),
                          test_df.assign(split="test")], ignore_index=True)
    manifest.to_csv(out_dir / "dataset_manifest.csv", index=False)

    # ── Class weights ─────────────────────────────────────────────────────────
    num_labels = task_cfg["num_labels"]
    classes_present = np.unique(train_df["label"].values).astype(int)
    cw_present = compute_class_weight("balanced", classes=classes_present,
                                      y=train_df["label"].values)
    cw = np.ones(num_labels, dtype=np.float32)
    for c, w in zip(classes_present, cw_present):
        cw[int(c)] = float(w)
    class_weights = torch.tensor(cw, dtype=torch.float, device=device)
    print(f"  Class weights: {dict(enumerate(cw.round(3).tolist()))}")

    # ── Dataloaders ───────────────────────────────────────────────────────────
    train_ds = build_bdino_dataset(train_df, train=True,
                                   augment=args.augment, aug_copies=args.aug_copies)
    print(f"  Train dataset size: {len(train_ds)} "
          f"(augment={args.augment}, aug_copies={args.aug_copies})")
    train_loader = DataLoader(
        train_ds, batch_size=args.batch_size, shuffle=True,
        num_workers=args.num_workers, pin_memory=True, drop_last=True)
    val_loader = DataLoader(
        build_bdino_dataset(val_df, train=False),
        batch_size=1, shuffle=False,
        num_workers=args.num_workers, pin_memory=True)
    test_loader = DataLoader(
        build_bdino_dataset(test_df, train=False),
        batch_size=1, shuffle=False,
        num_workers=args.num_workers, pin_memory=True)

    # ── Model ─────────────────────────────────────────────────────────────────
    model = BrainDINOClassifier(n_classes=num_labels, drop_rate=args.drop_rate)
    model = load_braindino_pretrained(model, args.pretrained_ckpt)

    if args.strategy == "frozen":
        for n, p in model.named_parameters():
            if not n.startswith("head"):
                p.requires_grad = False
        n_train = sum(p.numel() for p in model.parameters() if p.requires_grad)
        print(f"  Frozen: {n_train:,} trainable params (head only)")
    elif args.strategy == "lora":
        model = apply_lora_to_encoder(model, rank=args.lora_rank)
        n_train = sum(p.numel() for p in model.parameters() if p.requires_grad)
        print(f"  LoRA: {n_train:,} trainable params (adapters + head)")

    model = model.to(device)

    # ── Optimizer / scheduler / criterion ─────────────────────────────────────
    # Paper: Adam (not AdamW), weight_decay=1e-5, ReduceLROnPlateau on
    # val_bacc (mode=max, factor=0.5, patience=10). No warmup.
    trainable = [p for p in model.parameters() if p.requires_grad]
    optimizer = torch.optim.Adam(trainable, lr=args.lr,
                                 weight_decay=args.weight_decay)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="max",
        factor=args.scheduler_factor, patience=args.scheduler_patience)
    scaler = torch.cuda.amp.GradScaler(enabled=(device.type == "cuda"))
    criterion = nn.CrossEntropyLoss(weight=class_weights,
                                    label_smoothing=args.label_smoothing)

    wb = init_wandb(args, task_cfg, extra={
        "n_train": len(train_df), "n_val": len(val_df),
        "n_test": len(test_df), "model": MODEL_SLUG})

    # ── Resume ────────────────────────────────────────────────────────────────
    ckpt_path = out_dir / "last_checkpoint.pt"
    best_metric = -1.0
    best_val_loss_at_best = float("inf")
    best_epoch = -1
    epochs_since_improve = 0
    log_rows = []
    start_epoch = 0

    if (not args.no_resume) and ckpt_path.exists():
        ck = torch.load(ckpt_path, map_location=device)
        model.load_state_dict(ck["model"])
        optimizer.load_state_dict(ck["optimizer"])
        if "scheduler" in ck:
            scheduler.load_state_dict(ck["scheduler"])
        scaler.load_state_dict(ck["scaler"])
        best_metric = ck.get("best_metric", -1.0)
        best_val_loss_at_best = ck.get("best_val_loss_at_best", float("inf"))
        best_epoch = ck.get("best_epoch", -1)
        epochs_since_improve = ck.get("epochs_since_improve", 0)
        log_rows = ck.get("log_rows", [])
        start_epoch = ck["epoch"] + 1
        print(f"  Resuming from epoch {start_epoch + 1}")

    def _save_checkpoint(ep):
        state = {"model": model.state_dict(), "optimizer": optimizer.state_dict(),
                 "scheduler": scheduler.state_dict(),
                 "scaler": scaler.state_dict(), "epoch": ep,
                 "best_metric": best_metric,
                 "best_val_loss_at_best": best_val_loss_at_best,
                 "best_epoch": best_epoch,
                 "epochs_since_improve": epochs_since_improve,
                 "log_rows": log_rows}
        tmp = ckpt_path.with_suffix(".pt.tmp")
        torch.save(state, tmp)
        os.replace(tmp, ckpt_path)

    try:
        for epoch in range(start_epoch, args.epochs):
            # Adam keeps its own LR; ReduceLROnPlateau (stepped below on
            # val_bacc) is the only LR controller. Read the current LR
            # from the optimizer rather than computing it from a schedule.
            cur_lr = optimizer.param_groups[0]["lr"]

            tr_loss, tr_acc, _, _ = run_one_epoch(
                model, train_loader, criterion, optimizer, scaler, device, True,
                grad_accum_steps=args.grad_accum_steps)
            va_loss, va_acc, va_logits, va_labels = run_one_epoch(
                model, val_loader, criterion, optimizer, scaler, device, False)

            if not (np.isfinite(tr_loss) and np.isfinite(va_loss)):
                print(f"  [ABORT] non-finite loss at epoch {epoch+1} "
                      f"(tr_loss={tr_loss}, va_loss={va_loss}) -- stopping.")
                break

            # Full val-metric set (matches compute_test_metrics so W&B has
            # the same numbers we record for test): bacc, accuracy, AUC,
            # F1, precision, recall + binary-only sensitivity/specificity
            # and confusion-matrix entries (TP/FP/TN/FN). va_acc/val_loss
            # came from run_one_epoch; everything else is derived here.
            _val_full, va_preds, _ = compute_test_metrics(
                va_labels.astype(int), va_logits, task_cfg["task_type"])
            va_bacc = float(_val_full["balanced_acc"])
            val_extra = {
                "val_auc":    float(_val_full.get(
                    "auc_roc", _val_full.get("auc_roc_ovr", float("nan")))),
                "val_auc_pr": float(_val_full.get(
                    "auc_pr",  _val_full.get("auc_pr_macro", float("nan")))),
                "val_f1":     float(_val_full.get(
                    "f1",      _val_full.get("macro_f1", float("nan")))),
                "val_prec":   float(_val_full.get(
                    "precision", _val_full.get("precision_macro", float("nan")))),
                "val_recall": float(_val_full.get(
                    "recall",  _val_full.get("recall_macro", float("nan")))),
            }
            if task_cfg["task_type"] == "binary":
                val_extra["val_sens"] = float(_val_full["sensitivity"])
                val_extra["val_spec"] = float(_val_full["specificity"])
                _cm = confusion_matrix(
                    va_labels.astype(int), va_preds, labels=[0, 1])
                val_extra["val_TN"] = int(_cm[0, 0])
                val_extra["val_FP"] = int(_cm[0, 1])
                val_extra["val_FN"] = int(_cm[1, 0])
                val_extra["val_TP"] = int(_cm[1, 1])

            # Paper's ReduceLROnPlateau: factor=0.5, patience=10, on val_bacc.
            scheduler.step(va_bacc)

            log_rows.append({"epoch": epoch + 1, "lr": cur_lr,
                             "train_loss": tr_loss, "train_acc": tr_acc,
                             "val_loss": va_loss, "val_acc": va_acc,
                             "val_bacc": va_bacc, **val_extra})
            print(f"  [epoch {epoch+1:>3}/{args.epochs}] "
                  f"lr={cur_lr:.2e}  tr_loss={tr_loss:.4f}  "
                  f"va_loss={va_loss:.4f}  va_acc={va_acc:.4f}  "
                  f"va_bacc={va_bacc:.4f}  "
                  f"va_auc={val_extra['val_auc']:.3f}  "
                  f"va_f1={val_extra['val_f1']:.3f}")

            improved = (va_bacc > best_metric + 1e-6) or (
                abs(va_bacc - best_metric) <= 1e-6
                and va_loss < best_val_loss_at_best)
            if improved:
                best_metric = va_bacc
                best_val_loss_at_best = va_loss
                best_epoch = epoch + 1
                epochs_since_improve = 0
                torch.save({"net": model.state_dict(), "epoch": best_epoch},
                           out_dir / "best_model.pt")
            else:
                epochs_since_improve += 1

            pd.DataFrame(log_rows).to_csv(out_dir / "train_log.csv", index=False)
            _save_checkpoint(epoch)

            if wb is not None:
                import wandb
                wandb.log({"epoch": epoch + 1, "lr": cur_lr,
                           "train_loss": tr_loss, "val_loss": va_loss,
                           "val_acc": va_acc, "val_bacc": va_bacc,
                           "best_val_balanced_acc": best_metric,
                           **val_extra})

            if epochs_since_improve >= args.patience:
                print(f"  Early stopping at epoch {epoch+1}")
                break

        pd.DataFrame(log_rows).to_csv(out_dir / "train_log.csv", index=False)

        # ── Test ──────────────────────────────────────────────────────────────
        if not (out_dir / "best_model.pt").exists():
            print("  [ERROR] no best_model.pt -- training produced no usable epoch. "
                  "No metrics.json written.")
            return
        best_state = torch.load(out_dir / "best_model.pt", map_location=device)
        model.load_state_dict(best_state["net"])
        _, _, test_logits, test_labels = run_one_epoch(
            model, test_loader, criterion, optimizer, scaler, device, False)
        test_metrics, test_preds, test_probs = compute_test_metrics(
            test_labels, test_logits, task_cfg["task_type"])
        print(f"  Test metrics: {test_metrics}")

        label_names = class_names_for(task_cfg)
        test_diagnostics = compute_diagnostics(
            test_labels, test_preds, num_labels, label_names)

        pred_df = pd.DataFrame({"y_true": test_labels.astype(int),
                                "y_pred": test_preds.astype(int)})
        for c in range(test_probs.shape[1]):
            pred_df[f"prob_{c}"] = test_probs[:, c]
        pred_df.to_csv(out_dir / "test_predictions.csv", index=False)

        config = {
            "model_id":              MODEL_SLUG,
            "pretrained_ckpt":       args.pretrained_ckpt,
            "task":                  args.task,
            "task_description":      task_cfg["description"],
            "session":               cohort_label(args, task_cfg),
            "long_mode":             args.long_mode,
            "max_months":            args.max_months,
            "session_policy":        task_cfg["session_policy"],
            "seed":                  args.seed,
            "strategy":              args.strategy,
            "lora_rank":             args.lora_rank if args.strategy == "lora" else None,
            "augment":               args.augment,
            "aug_copies":            args.aug_copies,
            "epochs":                args.epochs,
            "best_epoch":            best_epoch,
            "best_val_balanced_acc": round(float(best_metric), 4),
            "lr":                    args.lr,
            "weight_decay":          args.weight_decay,
            "batch_size":            args.batch_size,
            "grad_accum_steps":      args.grad_accum_steps,
            "warmup_epochs":         args.warmup_epochs,
            "patience":              args.patience,
            "scheduler":             "ReduceLROnPlateau",
            "scheduler_factor":      args.scheduler_factor,
            "scheduler_patience":    args.scheduler_patience,
            "optimizer":             "Adam",
            "drop_rate":             args.drop_rate,
            "label_smoothing":       args.label_smoothing,
            "n_train":               len(train_df),
            "n_val":                 len(val_df),
            "n_test":                len(test_df),
            "class_weights":         {i: round(float(w), 4) for i, w in enumerate(cw)},
            "timestamp":             datetime.now().isoformat(),
        }
        with open(out_dir / "metrics.json", "w") as f:
            json.dump({"config": config, "test_metrics": test_metrics,
                       "test_diagnostics": test_diagnostics}, f, indent=2)
        print(f"  Saved: {out_dir}/metrics.json")
        print("=" * 70)
    finally:
        if wb is not None:
            import wandb
            wandb.finish()


if __name__ == "__main__":
    main()
