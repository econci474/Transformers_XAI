#!/usr/bin/env python
"""
Diagnostic probe: dump the first train batch that BMFM's dataloader hands to
the model, then exit before training starts.

Question we're answering: are the float z-score labels arriving at the model as
real continuous values, or have they been silently collapsed to zero (which
would explain the val_loss=0 / R²=1 trivial-prediction artefact on the combos
and by_chrom runs)?

Usage (CSD3 login node, env bmfm; no GPU needed):
    python snp_pipeline/11_debug_bmfm_dataloader.py \
        --config-path /home/ec474/rds/hpc-work/Transformers_XAI/snp_pipeline \
        -cn 09_bmfm_gwas_regression_finetuning \
        input_directory=/home/ec474/rds/hpc-work/ADNI_SNP/bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom_combos \
        working_dir=/tmp/bmfm_probe \
        "checkpoint='ibm-research/biomed.dna.ref.modernbert.113m.v1'"

What it does:
    1. Forces SDPA (same as force_sdpa_wrapper.py) so model init doesn't crash.
    2. Monkeypatches pytorch_lightning.Trainer.fit to:
         - reach into the datamodule, run setup("fit"),
         - pull one batch from train_dataloader(),
         - print the keys/shapes/dtypes and the first ~8 label values,
         - cross-reference against the raw z_score column in train.csv,
         - sys.exit(0).
    3. Delegates to bmfm_targets.tasks.scbert.scbert_main:main so all the
       Hydra/instantiate plumbing runs exactly as in production.
"""
import importlib
import sys
from pathlib import Path

# ── 1) Force SDPA on the BMFM model (same patch as force_sdpa_wrapper.py) ────
_mod = importlib.import_module(
    "bmfm_targets.models.predictive.scmodernbert.modeling_scmodernbert"
)
_OrigCls = _mod.SCModernBertForMultiTaskModeling
_orig_init = _OrigCls.__init__


def _sdpa_init(self, config, *args, **kwargs):
    config._attn_implementation = "sdpa"
    return _orig_init(self, config, *args, **kwargs)


_OrigCls.__init__ = _sdpa_init


# ── 1b) Wrap model.forward to dump output structure on first call ────────────
_orig_forward = _OrigCls.forward
_forward_call_count = {"n": 0}


def _instrumented_forward(self, *args, **kwargs):
    out = _orig_forward(self, *args, **kwargs)
    if _forward_call_count["n"] == 0:
        _forward_call_count["n"] = 1
        print("\n  >> [forward intercept] first forward pass — output structure:")
        if hasattr(out, "_fields"):
            keys = out._fields
            for k in keys:
                v = getattr(out, k)
                _print_any(f"    out.{k}", v)
        elif isinstance(out, dict):
            for k, v in out.items():
                _print_any(f"    out[{k!r}]", v)
        else:
            _print_any("    out", out)
    return out


def _print_any(label, v):
    import torch

    if isinstance(v, torch.Tensor):
        flat = v.detach().cpu().float().reshape(-1)
        head = flat[: min(8, flat.numel())].tolist()
        print(f"{label}: Tensor shape={tuple(v.shape)} dtype={v.dtype}  "
              f"mean={flat.mean().item():.4g} std={flat.std().item():.4g} "
              f"min={flat.min().item():.4g} max={flat.max().item():.4g}")
        print(f"{label}    head[:8]={head}")
    elif isinstance(v, dict):
        for k, sub in v.items():
            _print_any(f"{label}[{k!r}]", sub)
    elif isinstance(v, (list, tuple)):
        for i, sub in enumerate(v):
            _print_any(f"{label}[{i}]", sub)
    elif v is None:
        print(f"{label}: None")
    else:
        print(f"{label}: {type(v).__name__} = {v!r}")


_OrigCls.forward = _instrumented_forward


# ── 2) Intercept Trainer.fit to dump one batch and run one training step ────
import pytorch_lightning as pl  # noqa: E402

_orig_fit = pl.Trainer.fit


def _summarize_tensor(name, t):
    import torch

    if not isinstance(t, torch.Tensor):
        print(f"    {name}: {type(t).__name__} = {t!r}")
        return
    flat = t.detach().cpu().reshape(-1)
    head = flat[: min(16, flat.numel())].tolist()
    extras = ""
    if t.is_floating_point():
        extras = f"  mean={flat.float().mean().item():.6g}  std={flat.float().std().item():.6g}  min={flat.float().min().item():.6g}  max={flat.float().max().item():.6g}"
    else:
        extras = f"  unique={len(set(flat.tolist()))}  min={int(flat.min().item())}  max={int(flat.max().item())}"
    print(f"    {name}: shape={tuple(t.shape)}  dtype={t.dtype}{extras}")
    print(f"      head[:16] = {head}")


def _probe_fit(self, model=None, train_dataloaders=None, val_dataloaders=None, datamodule=None, **kwargs):
    print("\n" + "=" * 72)
    print("=== BMFM dataloader probe — first train batch ===")
    print("=" * 72)
    dm = datamodule if datamodule is not None else getattr(self, "datamodule", None)
    if dm is None and hasattr(model, "trainer") and getattr(model.trainer, "datamodule", None) is not None:
        dm = model.trainer.datamodule
    if dm is not None:
        try:
            dm.setup("fit")
        except Exception as e:
            print(f"  [warn] dm.setup raised: {e!r}")
        loader = dm.train_dataloader()
        # Try to find the underlying csv path for cross-reference.
        cand = None
        for attr in ("processed_data_source", "input_directory", "data_dir"):
            v = getattr(dm, attr, None)
            if v:
                cand = Path(str(v)) / "train.csv"
                if cand.exists():
                    break
                cand = None
        if cand is None and hasattr(dm, "hparams"):
            for k in ("processed_data_source", "input_directory", "data_dir"):
                v = dm.hparams.get(k) if hasattr(dm.hparams, "get") else None
                if v:
                    cand = Path(str(v)) / "train.csv"
                    if cand.exists():
                        break
                    cand = None
        if cand is not None:
            print(f"\n  Cross-reference: raw train.csv → {cand}")
            try:
                import pandas as pd

                head = pd.read_csv(cand, usecols=["z_score"], nrows=16)["z_score"].tolist()
                print(f"  raw z_score[:16] = {head}")
            except Exception as e:
                print(f"  [warn] failed to read train.csv: {e!r}")
        else:
            print("\n  Cross-reference: could not locate train.csv from datamodule attrs.")
    else:
        loader = train_dataloaders
        print("  [warn] datamodule not found; falling back to train_dataloaders arg.")

    if loader is None:
        print("  [err] no dataloader available.")
        sys.exit(2)

    print("\n  Pulling first batch from train_dataloader …")
    batch = next(iter(loader))
    print(f"  Batch type: {type(batch).__name__}")
    if isinstance(batch, dict):
        for k, v in batch.items():
            _summarize_tensor(k, v)
    elif isinstance(batch, (list, tuple)):
        for i, v in enumerate(batch):
            _summarize_tensor(f"[{i}]", v)
    else:
        print(f"  unknown batch shape: {batch!r}")

    # ── Run forward+loss TWICE: fp32 then bf16-autocast, compare ──
    print("\n" + "-" * 72)
    print("  fp32 vs bf16-autocast forward+loss comparison")
    print("-" * 72)
    import torch  # noqa: E402
    from contextlib import nullcontext  # noqa: E402

    module = model  # In Lightning, .fit(model, …) — `model` is the LightningModule
    if module is None:
        print("  [err] LightningModule not available; can't run training_step.")
        sys.exit(3)

    dev = next(module.parameters()).device
    print(f"  Module device: {dev}; first param dtype: {next(module.parameters()).dtype}")

    def _to_dev(x):
        if isinstance(x, torch.Tensor):
            return x.to(dev)
        if isinstance(x, dict):
            return {k: _to_dev(v) for k, v in x.items()}
        if isinstance(x, list):
            return [_to_dev(v) for v in x]
        return x

    batch_dev = _to_dev(batch)

    # Replicate BMFM's calculate_losses without retaining autograd.
    from bmfm_targets.training.losses.utils import calculate_losses as _calc  # noqa: E402

    def _summary(label, t):
        if not torch.is_tensor(t):
            return f"{label}={t!r}"
        f = t.detach().float().reshape(-1)
        has_nan = bool(torch.isnan(t).any().item())
        has_inf = bool(torch.isinf(t).any().item())
        return (f"{label}: shape={tuple(t.shape)} dtype={t.dtype} "
                f"mean={f.mean().item():.4g} std={f.std().item() if f.numel() > 1 else float('nan'):.4g} "
                f"min={f.min().item():.4g} max={f.max().item():.4g} "
                f"has_nan={has_nan} has_inf={has_inf}")

    def _check(name, ctx):
        print(f"\n>>> {name} <<<")
        _forward_call_count["n"] = 0
        module.eval()  # disable dropout to avoid extra allocations
        with torch.no_grad(), ctx:
            # Forward pass only — no backward graph retained.
            try:
                out = module.model(
                    input_ids=batch_dev["input_ids"],
                    attention_mask=batch_dev["attention_mask"],
                )
            except Exception as e:
                print(f"  [err] model.forward raised: {type(e).__name__}: {e}")
                import traceback
                traceback.print_exc()
                return None
            # Inspect logits
            logits = out["logits"] if isinstance(out, dict) else getattr(out, "logits", None)
            if logits is None:
                print(f"  [err] no .logits on model output: keys={list(out.keys()) if hasattr(out, 'keys') else type(out)}")
                return None
            for k, v in logits.items():
                print(f"    {_summary(f'logits[{k!r}]', v)}")
            # Now compute MSE via BMFM's own pipeline
            try:
                all_losses = _calc(module.loss_tasks, logits, batch_dev["labels"])
            except Exception as e:
                print(f"  [err] calculate_losses raised: {type(e).__name__}: {e}")
                import traceback
                traceback.print_exc()
                return None
            for k, v in all_losses.items():
                print(f"    {_summary(f'all_losses[{k!r}]', v) if torch.is_tensor(v) else f'all_losses[{k!r}]={v!r}'}")
        return all_losses.get("loss")

    # If CUDA is available, also try GPU + cuda-autocast(bfloat16) — that's
    # what bf16-mixed actually runs on CSD3.
    if torch.cuda.is_available():
        print(f"  CUDA available: {torch.cuda.get_device_name(0)} — moving module to GPU.")
        module.to("cuda")
        batch_dev = _to_dev(batch)  # re-do with cuda as the new module device

        print("\n[1/4] CPU fp32 (sanity)")
        module.to("cpu")
        batch_dev_cpu = _to_dev(batch)
        fp32_cpu = _check("CPU FP32", nullcontext())
        # rebind dev for later prints
        module.to("cuda")
        batch_dev = _to_dev(batch)

        print("\n[2/4] CUDA fp32")
        cuda_fp32 = _check("CUDA FP32", nullcontext())

        print("\n[3/4] CUDA bf16 autocast")
        cuda_bf16 = _check("CUDA BF16", torch.autocast(device_type="cuda", dtype=torch.bfloat16))

        print("\n[4/4] CUDA fp16 autocast (for comparison)")
        cuda_fp16 = _check("CUDA FP16", torch.autocast(device_type="cuda", dtype=torch.float16))

        print("\n" + "=" * 72)
        print(f"  CPU fp32     = {fp32_cpu}")
        print(f"  CUDA fp32    = {cuda_fp32}")
        print(f"  CUDA bf16    = {cuda_bf16}   <-- this should be 0 / NaN if the CSD3 bug reproduces")
        print(f"  CUDA fp16    = {cuda_fp16}")
        print("=" * 72)
    else:
        print("\n[1/2] fp32 (no autocast)")
        fp32_loss = _check("FP32", nullcontext())

        print("\n[2/2] bf16 autocast on cpu")
        try:
            bf16_ctx = torch.autocast(device_type="cpu", dtype=torch.bfloat16)
        except Exception as e:
            print(f"  [warn] torch.autocast(cpu, bfloat16) failed: {e}; falling back to no-op")
            bf16_ctx = nullcontext()
        bf16_loss = _check("BF16 autocast", bf16_ctx)

        print("\n" + "=" * 72)
        print(f"  SUMMARY:  fp32_loss = {fp32_loss}   bf16_loss = {bf16_loss}")
        print("=" * 72)
    sys.exit(0)


pl.Trainer.fit = _probe_fit


# ── 3) Delegate to BMFM main ─────────────────────────────────────────────────
from bmfm_targets.tasks.scbert.scbert_main import main  # noqa: E402

if __name__ == "__main__":
    main()
