"""
colab_setup_env.py
===================
Codify the (hard-won) Colab A100 environment that ``bmfm_targets`` + script 18
(``18_patient_classification_finetune.py``) need to import and run.

A fresh Colab VM ships packages that break BMFM in non-obvious ways. This
module makes the working environment **reproducible and idempotent** so the
dependency archaeology never has to be repeated.

What it does
------------
1. ``set_env()``  — export the exact env vars BMFM needs (TF off, torch-only
   transformers, protobuf python impl, HF cache on Drive, expandable CUDA
   allocator, TorchDynamo/compile disabled) and prepend the persisted libs
   dir to ``sys.path`` + ``PYTHONPATH``.
2. ``repair_libs()`` — one-time, persisted-on-Drive package repairs (each
   guarded by an import probe so re-runs are no-ops):
     * remove the *TensorFlow* ``focal-loss`` pkg, install the *pytorch*
       ``focal_loss_torch`` (bmfm wants ``from focal_loss.focal_loss import
       FocalLoss``);
     * pin ``transformers>=4.54,<5.0`` (base Colab ships 5.0 which removed
       ``find_pruneable_heads_and_indices`` that bmfm needs);
     * uninstall ``torchao`` (peft 0.19 hard-raises on Colab's old 0.10);
     * ensure ``numcodecs``/``zarr``; strip stray ``tensorflow``/``keras``
       *inside the persisted libs dir only* (never base site-packages).
3. ``verify()`` — import ``bmfm_targets`` + ``SCModernBertForMultiTaskModeling``,
   print torch/transformers/peft versions + the CUDA device, and emit the
   sentinel ``[colab_setup_env] OK``.

Usage on Colab
--------------
Run once per fresh VM for the persistent repairs + verification::

    !python /content/Transformers_XAI/snp_pipeline/colab_setup_env.py

Then, in the *kernel cell that launches training*, make the env vars reach the
``!python …18….py`` subprocess by setting them in the kernel first::

    import sys; sys.path.insert(0, "/content/Transformers_XAI/snp_pipeline")
    import colab_setup_env; colab_setup_env.set_env()
    # …now any `!python …18….py` inherits USE_TF=0 / PYTHONPATH / HF_HOME / …

``set_env()`` has no heavy side effects and is safe to call from the kernel;
the expensive ``repair_libs()`` only needs the one ``!python`` invocation
(its installs persist on Drive).

Notes
-----
- ``LIBS`` is on Drive so pip ``--target`` installs survive VM resets.
- ``HF_HOME`` on Drive keeps the BMFM-DNA-REF ckpt + ref2vec tokenizer cached.
- This is generic to any BMFM Colab session, not just script 18.
"""
from __future__ import annotations

import importlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

# ── Persisted locations (on Google Drive → survive VM resets) ────────────────
DRIVE_ROOT = Path("/content/drive/MyDrive/ADNI_SNP")
LIBS    = DRIVE_ROOT / "python_libs"
HF_HOME = DRIVE_ROOT / "hf_cache"

# ── The exact env BMFM needs (see module docstring for the "why" of each) ────
ENV: dict[str, str] = {
    "USE_TF": "0",                 # transformers torch-only; base Colab TF is
    "USE_TORCH": "1",              #   broken (ABSL symbol) + drags a protobuf
                                   #   gencode/runtime conflict via TF.
    "PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION": "python",
    "TRANSFORMERS_NO_ADVISORY_WARNINGS": "1",
    "HF_HOME": str(HF_HOME),
    "PYTORCH_CUDA_ALLOC_CONF": "expandable_segments:True",
    "TORCHDYNAMO_DISABLE": "1",    # Drive-mounted triton ptxas is non-exec
    "TORCH_COMPILE_DISABLE": "1",  #   (FUSE strips +x) → inductor crash.
}


def _log(msg: str) -> None:
    print(f"[colab_setup_env] {msg}", flush=True)


def ensure_dirs() -> None:
    for d in (LIBS, HF_HOME):
        d.mkdir(parents=True, exist_ok=True)


def set_env() -> None:
    """Export the BMFM env vars and put the persisted libs dir first on the
    import path. Cheap + idempotent — safe to call from the Colab kernel so
    that subsequent ``!python`` subprocesses inherit it."""
    ensure_dirs()
    for k, v in ENV.items():
        os.environ[k] = v
    libs = str(LIBS)
    if libs not in sys.path:
        sys.path.insert(0, libs)
    cur = os.environ.get("PYTHONPATH", "")
    parts = [p for p in cur.split(os.pathsep) if p]
    if libs not in parts:
        os.environ["PYTHONPATH"] = os.pathsep.join([libs, *parts])
    _log(f"env set ({len(ENV)} vars); PYTHONPATH[0]={libs}")


def _pip(*args: str) -> int:
    """Run pip in this interpreter; return its exit code (non-raising)."""
    cmd = [sys.executable, "-m", "pip", *args]
    _log("pip " + " ".join(args))
    return subprocess.run(cmd).returncode


def _importable(modpath: str) -> bool:
    """True if ``modpath`` imports cleanly in a *fresh* subprocess (so a
    failed/cached import here does not poison the probe)."""
    code = f"import importlib,sys; importlib.import_module('{modpath}')"
    return subprocess.run([sys.executable, "-c", code]).returncode == 0


def _rm_libs_pkg(*names: str) -> None:
    """Delete stray top-level package dirs *inside LIBS only* (never base
    site-packages)."""
    for name in names:
        for p in (LIBS / name, LIBS / f"{name}.py"):
            if p.exists():
                _log(f"removing stray {p} from persisted libs")
                shutil.rmtree(p, ignore_errors=True) if p.is_dir() else p.unlink()
        for dist in LIBS.glob(f"{name}-*.dist-info"):
            shutil.rmtree(dist, ignore_errors=True)


def repair_libs() -> None:
    """One-time persisted package repairs. Each step is guarded by an import
    probe so re-running on an already-fixed Drive is a no-op."""
    ensure_dirs()
    target = ["--target", str(LIBS), "--no-warn-script-location"]

    # (1) focal loss: bmfm wants the *pytorch* one
    #     (`from focal_loss.focal_loss import FocalLoss`). The TF `focal-loss`
    #     pkg ships an incompatible top-level `focal_loss` → uninstall it,
    #     install `focal_loss_torch` into the persisted libs (so it shadows).
    if not _importable("focal_loss.focal_loss"):
        _pip("uninstall", "-y", "focal-loss", "focal_loss")
        _rm_libs_pkg("focal_loss")
        _pip("install", *target, "focal_loss_torch")
    else:
        _log("focal_loss_torch OK (skip)")

    # (2) transformers: base Colab ships 5.0 which removed
    #     `find_pruneable_heads_and_indices` that bmfm imports. Pin <5.0.
    def _tf_ok() -> bool:
        code = ("import transformers as t; from packaging.version import parse "
                "as v; import sys; "
                "sys.exit(0 if v('4.54') <= v(t.__version__) < v('5.0') else 1)")
        return subprocess.run([sys.executable, "-c", code]).returncode == 0
    if not _tf_ok():
        _pip("install", *target, "transformers>=4.54,<5.0")
    else:
        _log("transformers in [4.54,5.0) (skip)")

    # (3) torchao: peft 0.19 hard-raises against Colab's old base 0.10. We do
    #     not use quantized LoRA, so removing it → graceful skip in peft.
    if _importable("torchao"):
        _pip("uninstall", "-y", "torchao")
        _rm_libs_pkg("torchao")
    else:
        _log("torchao absent (skip)")

    # (4) zarr stack bmfm_targets imports transitively.
    for mod, pkg in (("numcodecs", "numcodecs"), ("zarr", "zarr")):
        if not _importable(mod):
            _pip("install", *target, pkg)
        else:
            _log(f"{mod} OK (skip)")

    # (5) never let a stray persisted TF/keras shadow the torch-only path.
    _rm_libs_pkg("tensorflow", "keras")


def verify() -> int:
    """Import the critical stack and print versions + device. Returns 0 on
    success, 1 on failure (so ``main`` can exit non-zero)."""
    set_env()
    try:
        import torch  # noqa: F401
        import transformers  # noqa: F401
        import peft  # noqa: F401
        import bmfm_targets  # noqa: F401
        from bmfm_targets.models.predictive.scmodernbert.modeling_scmodernbert import (  # noqa: F401,E501
            SCModernBertForMultiTaskModeling,
        )
    except Exception as e:  # pragma: no cover - environment-dependent
        _log(f"VERIFY FAILED: {type(e).__name__}: {e}")
        return 1
    import torch
    import transformers
    import peft
    dev = (torch.cuda.get_device_name(0)
           if torch.cuda.is_available() else "CPU (no CUDA!)")
    _log(f"torch={torch.__version__} transformers={transformers.__version__} "
         f"peft={peft.__version__}")
    _log(f"cuda_available={torch.cuda.is_available()} device={dev}")
    if not torch.cuda.is_available():
        _log("WARNING: no CUDA — switch the Colab runtime to an A100 GPU.")
    return 0


def main() -> None:
    _log(f"LIBS={LIBS}")
    _log(f"HF_HOME={HF_HOME}")
    set_env()
    repair_libs()
    rc = verify()
    if rc == 0:
        _log("OK")  # sentinel grepped by the orchestration
    else:
        _log("FAILED")
    sys.exit(rc)


if __name__ == "__main__":
    main()
