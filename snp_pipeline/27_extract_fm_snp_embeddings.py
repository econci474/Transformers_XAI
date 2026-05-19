r"""
27_extract_fm_snp_embeddings.py   (COLAB A100, env: see colab_setup_env)
=======================================================================
Run the 4 FROZEN DNA foundation models over the script-26 per-SNP REF/ALT
1 kb sequences and cache, per (snp_set, model), the REF and ALT embeddings.
The signed contrast diff = ALT_emb − REF_emb and the per-patient
d = diff·β·dosage are computed downstream by fm_diff_lib (one place for the
sign/orientation logic).

Models (frozen, no grad; record HF id):
  bmfm_ref     ibm-research/biomed.dna.ref.modernbert.113m.v1     CLS  768
  bmfm_snp     ibm-research/biomed.dna.snp.modernbert.113m.v1     CLS  768
  ntv2         InstaDeepAI/nucleotide-transformer-v2-500m-multi…  mean 1024
  caduceus_ph  kuleshov-group/caduceus-ph_seqlen-1k_d_model-256_… mean 256
                                            (fwd+RC averaged; RC in-script)

Reuses verbatim: BMFM load/tokenise/CLS — 14_extract_patient_embeddings.py
L103-182; NT-v2 mean_pool/forward — nt_v2/03 L26-120; Caduceus load/forward/
fwd+RC-avg — caduceus/03 L43-211; Colab env — colab_setup_env.set_env().

Input  : <Drive>/ADNI_SNP/fm/fm_sequences/<SET>_snp_sequences.tsv  (script 26)
Output : <outdir>/<model>/<SET>_snp_embeddings.npz  (+ _embed_qc.json,
         _embed_run.json). Keys: rsids, snp_idx, ref_emb, alt_emb,
         model_tag, model_id, hidden_dim, pooling, strand_mode.

Usage (on Colab, after colab_setup_env in the kernel):
  !python snp_pipeline/27_extract_fm_snp_embeddings.py \
      --sequences /content/drive/MyDrive/ADNI_SNP/fm/fm_sequences/9W27B_snp_sequences.tsv \
      --set 9W27B --model bmfm_ref \
      --outdir /content/drive/MyDrive/ADNI_SNP/fm/fm_embeddings --limit 8
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# Defensive Colab env (idempotent; no-op off-Colab / if module absent)
try:
    sys.path.insert(0, str(Path(__file__).parent))
    import colab_setup_env  # noqa: E402
    colab_setup_env.set_env()
except Exception as _e:                                   # pragma: no cover
    print(f"[colab_setup_env] skipped ({_e!r})")

import torch  # noqa: E402

MODELS = {
    "bmfm_ref":  {"id": "ibm-research/biomed.dna.ref.modernbert.113m.v1",
                  "kind": "bmfm", "pool": "cls", "H": 768, "max_len": 2048},
    "bmfm_snp":  {"id": "ibm-research/biomed.dna.snp.modernbert.113m.v1",
                  "kind": "bmfm", "pool": "cls", "H": 768, "max_len": 2048},
    "ntv2":      {"id": "InstaDeepAI/nucleotide-transformer-v2-500m-"
                        "multi-species",
                  "kind": "hf_mlm", "pool": "mean", "H": 1024,
                  "max_len": 2048, "rc": False},
    "caduceus_ph": {"id": "kuleshov-group/caduceus-ph_seqlen-1k_"
                          "d_model-256_n_layer-4_lr-8e-3",
                    "kind": "hf_mlm", "pool": "mean", "H": 256,
                    "max_len": 1024, "rc": True},
}

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:                  # caduceus/03 L43-46
    return seq.translate(_COMP)[::-1]


def mean_pool(hidden, attn):                              # nt_v2/03 L26-29
    m = attn.unsqueeze(-1).float()
    return (hidden * m).sum(1) / m.sum(1).clamp_min(1)


def _load_bmfm(model_id, device):                         # 14 L103-139
    from bmfm_targets.config.main_config import load_checkpoint_configs
    from bmfm_targets.models.model_utils import load_bmfm_model
    from bmfm_targets.tokenization.multifield_tokenizer import (
        MultiFieldTokenizer)
    cfg = load_checkpoint_configs(model_id)
    model = load_bmfm_model(cfg)
    if hasattr(model, "config"):
        model.config._attn_implementation = "sdpa"
    tok = (MultiFieldTokenizer.from_pretrained(model_id)
           if hasattr(MultiFieldTokenizer, "from_pretrained") else None)
    if tok is None:
        from bmfm_targets.models.model_utils import load_tokenizer
        tok = load_tokenizer(cfg)
    return model.eval().to(device), tok


def _embed_bmfm(model, tok, seqs, device, max_len):       # 14 L169-182
    t = tok([{"dna_chunks": s} for s in seqs], padding="max_length",
            truncation=True, max_length=max_len, return_tensors="pt")
    out = model(input_ids=t["input_ids"].to(device),
                attention_mask=t["attention_mask"].to(device))
    emb = out["embeddings"] if isinstance(out, dict) else out.embeddings
    return emb.detach().to(torch.float32).cpu().numpy()


def _load_hf(model_id, device):                           # nt_v2/caduceus
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    tok = AutoTokenizer.from_pretrained(model_id, trust_remote_code=True)
    model = AutoModelForMaskedLM.from_pretrained(
        model_id, trust_remote_code=True).to(device).eval()
    return model, tok


def _embed_hf(model, tok, seqs, device, max_len):
    t = tok(seqs, return_tensors="pt", padding="max_length",
            truncation=True, max_length=max_len)
    ids = t["input_ids"].to(device)
    am = (t["attention_mask"].to(device) if "attention_mask" in t
          else (ids != (tok.pad_token_id or 0)).long())
    out = model(ids, attention_mask=am, output_hidden_states=True)
    return mean_pool(out.hidden_states[-1], am).detach().to(
        torch.float32).cpu().numpy()


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sequences", type=Path, required=True)
    ap.add_argument("--set", required=True)
    ap.add_argument("--model", required=True, choices=list(MODELS))
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--device",
                    default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--limit", type=int, default=None,
                    help="smoke: only first N SNPs")
    args = ap.parse_args()

    spec = MODELS[args.model]
    device = torch.device(args.device)
    print(f"[27] set={args.set} model={args.model} id={spec['id']} "
          f"device={device} torch={torch.__version__}", flush=True)

    df = pd.read_csv(args.sequences, sep="\t")
    df = df[df["strand"] == "fwd"].copy()                 # rc done in-script
    piv = df.pivot_table(index=["snp_idx", "rsID"], columns="allele_type",
                         values="sequence", aggfunc="first").reset_index()
    piv = piv.sort_values("snp_idx").reset_index(drop=True)
    if args.limit:
        piv = piv.head(args.limit)
    if not {"REF", "ALT"}.issubset(piv.columns):
        sys.exit("[ERROR] sequence file lacks REF/ALT rows")
    rsids = piv["rsID"].astype(str).tolist()
    snp_idx = piv["snp_idx"].astype(int).to_numpy()
    ref_seqs = piv["REF"].tolist()
    alt_seqs = piv["ALT"].tolist()
    n = len(rsids)
    print(f"  {n} SNPs (REF+ALT)")

    if spec["kind"] == "bmfm":
        model, tok = _load_bmfm(spec["id"], device)
        emb_fn = lambda S: _embed_bmfm(model, tok, S, device, spec["max_len"])
        strand_mode = "fwd"
    else:
        model, tok = _load_hf(spec["id"], device)
        emb_fn = lambda S: _embed_hf(model, tok, S, device, spec["max_len"])
        strand_mode = "fwd_rc_avg" if spec.get("rc") else "fwd"

    def embed_all(seqs):
        out = np.zeros((len(seqs), spec["H"]), np.float32)
        bs = args.batch_size
        with torch.no_grad():
            for i in range(0, len(seqs), bs):
                chunk = seqs[i:i + bs]
                e = emb_fn(chunk)
                if spec.get("rc"):                        # caduceus_ph
                    e = (e + emb_fn([reverse_complement(s)
                                     for s in chunk])) / 2.0
                out[i:i + bs] = e
                print(f"  {min(i + bs, len(seqs))}/{len(seqs)}", flush=True)
        return out

    ref_emb = embed_all(ref_seqs)
    alt_emb = embed_all(alt_seqs)

    # ── QC ───────────────────────────────────────────────────────────────
    qc = {"set": args.set, "model": args.model, "n_snps": n,
          "hidden_dim": int(ref_emb.shape[1]),
          "expected_H": spec["H"], "strand_mode": strand_mode}
    assert ref_emb.shape == alt_emb.shape == (n, spec["H"]), \
        f"shape {ref_emb.shape}/{alt_emb.shape} != ({n},{spec['H']})"
    fin = bool(np.isfinite(ref_emb).all() and np.isfinite(alt_emb).all())
    qc["all_finite"] = fin
    if not fin:
        sys.exit("[ERROR] non-finite embeddings — aborting (not saved)")
    rn = np.linalg.norm(ref_emb, axis=1)
    an = np.linalg.norm(alt_emb, axis=1)
    qc["min_norm"] = float(min(rn.min(), an.min()))
    if qc["min_norm"] < 1e-6:
        sys.exit("[ERROR] zero-norm embedding — tokenise/forward failure")
    # allele-sensitivity headline: per-SNP cos(ref,alt) & ‖alt-ref‖
    cos = ((ref_emb * alt_emb).sum(1)
           / (rn * an + 1e-12)).astype(float)
    dn = np.linalg.norm(alt_emb - ref_emb, axis=1).astype(float)
    qc["cos_ref_alt"] = {"mean": float(cos.mean()), "min": float(cos.min()),
                         "max": float(cos.max())}
    qc["diff_norm"] = {"mean": float(dn.mean()), "min": float(dn.min()),
                       "max": float(dn.max())}
    near = [rsids[i] for i in range(n)
            if cos[i] > 0.99999 or dn[i] < 1e-5]
    qc["near_identical_ref_alt"] = near
    print(f"  [qc] cos(ref,alt) mean={cos.mean():.5f} "
          f"min={cos.min():.5f}; ‖diff‖ mean={dn.mean():.4f}; "
          f"near_identical={len(near)}")

    # determinism re-embed (first ≤4)
    k = min(4, n)
    re_ref = embed_all(ref_seqs[:k])
    qc["determinism_allclose"] = bool(
        np.allclose(re_ref, ref_emb[:k], atol=1e-4))

    outdir = args.outdir / args.model
    outdir.mkdir(parents=True, exist_ok=True)
    npz = outdir / f"{args.set}_snp_embeddings.npz"
    tmp = outdir / f".{args.set}_snp_embeddings.tmp.npz"   # Drive-safe write
    np.savez_compressed(
        tmp, rsids=np.array(rsids), snp_idx=snp_idx,
        ref_emb=ref_emb, alt_emb=alt_emb,
        model_tag=args.model, model_id=spec["id"],
        hidden_dim=spec["H"], pooling=spec["pool"],
        strand_mode=strand_mode)
    shutil.move(str(tmp), str(npz))
    (outdir / f"{args.set}_embed_qc.json").write_text(
        json.dumps(qc, indent=2), encoding="utf-8")
    seq_sha = hashlib.sha256(
        Path(args.sequences).read_bytes()).hexdigest()[:16]
    (outdir / f"{args.set}_embed_run.json").write_text(json.dumps({
        "set": args.set, "model": args.model, "model_id": spec["id"],
        "torch": torch.__version__,
        "gpu": (torch.cuda.get_device_name(0)
                if torch.cuda.is_available() else "cpu"),
        "sequences": str(args.sequences), "seq_sha256_16": seq_sha,
        "n_snps": n, "limit": args.limit,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")},
        indent=2), encoding="utf-8")
    print(f"[out] {npz}  ({n}×{spec['H']}, strand={strand_mode}, "
          f"determinism={qc['determinism_allclose']})")


if __name__ == "__main__":
    main()
