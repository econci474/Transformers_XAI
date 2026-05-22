"""
test_vit_preprocessing.py
=========================
Unit tests for the ViT preprocessing pipeline and the training-time input gate:

  * 03_prepare_ViT.py               -- build_transform() correctness
  * 04_supervised_finetuning_ViT.py -- preflight_check_inputs() accepts good
                                       volumes and rejects malformed ones;
                                       TASK_CONFIG label maps are consistent.

All tests use SYNTHETIC data only -- no real preprocessed scans are needed, so
they run anywhere the `mri` conda env is installed (no GPU required).

Run:  conda run -n mri pytest mri_pipeline/tests/test_vit_preprocessing.py -v
"""
import importlib.util
from pathlib import Path

import numpy as np
import nibabel as nib
import pandas as pd
import pytest

MRI_DIR = Path(__file__).resolve().parents[1]   # .../mri_pipeline


def _load_script(filename: str, modname: str):
    """Import a digit-prefixed pipeline script (e.g. 03_prepare_ViT.py) by path
    -- a plain `import 03_prepare_ViT` is not valid Python."""
    spec = importlib.util.spec_from_file_location(modname, MRI_DIR / filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


prep = _load_script("03_prepare_ViT.py", "_prep_vit")
train = _load_script("04_supervised_finetuning_ViT.py", "_train_vit")


# ── helpers ───────────────────────────────────────────────────────────────────
def _save_nifti(path: Path, arr: np.ndarray, spacing=(1.0, 1.0, 1.0)) -> Path:
    affine = np.diag([*spacing, 1.0])
    nib.save(nib.Nifti1Image(arr.astype(np.float32), affine), str(path))
    return path


def _synthetic_brain(shape=(150, 150, 150), seed=0) -> np.ndarray:
    """A brain-like volume: a positive-intensity ellipsoid in zero padding."""
    rng = np.random.default_rng(seed)
    vol = np.zeros(shape, np.float32)
    cx, cy, cz = (s // 2 for s in shape)
    rx, ry, rz = (s // 3 for s in shape)
    zz, yy, xx = np.ogrid[:shape[0], :shape[1], :shape[2]]
    mask = (((xx - cx) / rx) ** 2 + ((yy - cy) / ry) ** 2
            + ((zz - cz) / rz) ** 2) <= 1.0
    vol[mask] = rng.normal(600.0, 120.0, size=int(mask.sum())).astype(np.float32)
    return vol


def _zscored_volume(seed=0) -> np.ndarray:
    """A 128^3 already-standardised volume (mean ~0, std ~1, all nonzero)."""
    rng = np.random.default_rng(seed)
    return rng.standard_normal((128, 128, 128)).astype(np.float32)


def _make_split(tmp_path: Path, vols, labels, tag: str) -> pd.DataFrame:
    rows = []
    for i, (v, lab) in enumerate(zip(vols, labels)):
        p = _save_nifti(tmp_path / f"{tag}_{i}.nii.gz", v)
        rows.append({"label": lab, "image_path": str(p)})
    return pd.DataFrame(rows)


# ── 03_prepare_ViT.py : build_transform() ─────────────────────────────────────
def test_build_transform_shape_and_zscore(tmp_path):
    """The transform must yield a 128^3 volume, finite, roughly standardised."""
    src = _save_nifti(tmp_path / "syn.nii.gz", _synthetic_brain())
    out = prep.build_transform()({"image": str(src)})["image"]
    arr = np.asarray(out).squeeze()

    assert arr.shape == (128, 128, 128), f"shape {arr.shape}"
    assert np.isfinite(arr).all(), "output has NaN/Inf"

    nz = arr[arr != 0]
    assert nz.size > 800_000, f"brain too small: {nz.size} nonzero voxels"
    m, s = float(nz.mean()), float(nz.std())
    # Loose bands -- the final resize perturbs exact 0/1 (see build_transform).
    assert -0.2 <= m <= 0.2, f"nonzero mean {m:.3f} outside band"
    assert 0.75 <= s <= 1.15, f"nonzero std {s:.3f} outside band"


# ── 04 : TASK_CONFIG consistency ──────────────────────────────────────────────
def test_task_config_label_maps_consistent():
    """Every label_map value must be a valid class index in [0, num_labels)."""
    for name, cfg in train.TASK_CONFIG.items():
        n = cfg["num_labels"]
        lm = cfg.get("label_map")
        if lm is not None:
            assert set(lm.values()) <= set(range(n)), (
                f"{name}: label_map values {set(lm.values())} exceed "
                f"num_labels={n}")


def test_t1c_binary_is_cn_vs_ad():
    """T1c_binary must be CN-vs-AD with MCI absent (so MCI scans are dropped)."""
    cfg = train.TASK_CONFIG["T1c_binary"]
    assert cfg["label_map"] == {"CN": 0, "AD": 1}
    assert cfg["num_labels"] == 2
    assert cfg["task_type"] == "binary"
    assert "MCI" not in cfg["label_map"]


# ── 04 : preflight_check_inputs() ─────────────────────────────────────────────
def test_preflight_accepts_good_inputs(tmp_path):
    """Well-formed standardised 128^3 volumes must pass without raising."""
    vols = [_zscored_volume(s) for s in range(6)]
    df = _make_split(tmp_path, vols, [0, 1, 0, 1, 0, 1], "good")
    train.preflight_check_inputs({"train": df, "val": df, "test": df},
                                 sample_k=6)


def test_preflight_rejects_unnormalised(tmp_path):
    """Raw MRI-scale intensities (NormalizeIntensity not applied) must fail."""
    rng = np.random.default_rng(1)
    raw = (rng.standard_normal((128, 128, 128)) * 200 + 800).astype(np.float32)
    df = _make_split(tmp_path, [raw] * 3, [0, 1, 0], "raw")
    with pytest.raises(RuntimeError, match="z-scored"):
        train.preflight_check_inputs({"train": df}, sample_k=3)


def test_preflight_rejects_wrong_shape(tmp_path):
    bad = _zscored_volume()[:64, :64, :64]
    df = _make_split(tmp_path, [bad] * 3, [0, 1, 0], "shape")
    with pytest.raises(RuntimeError, match="shape"):
        train.preflight_check_inputs({"train": df}, sample_k=3)


def test_preflight_rejects_nan(tmp_path):
    v = _zscored_volume()
    v[0, 0, 0] = np.nan
    df = _make_split(tmp_path, [v] * 3, [0, 1, 0], "nan")
    with pytest.raises(RuntimeError, match="NaN"):
        train.preflight_check_inputs({"train": df}, sample_k=3)


def test_preflight_rejects_single_class(tmp_path):
    """A train/val split with only one class cannot be trained/validated."""
    vols = [_zscored_volume(s) for s in range(3)]
    df = _make_split(tmp_path, vols, [0, 0, 0], "oneclass")
    with pytest.raises(RuntimeError, match="<2 classes"):
        train.preflight_check_inputs({"train": df}, sample_k=3)
