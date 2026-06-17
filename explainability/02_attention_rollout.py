r"""
02_attention_rollout.py   (env: clinical; GPU or CPU)
=====================================================
Attention-rollout explainability for the BEST fine-tuned LLM in each of T1d / T1e / T2, overlaid on the
clinical free-text so one can see which words / test-result blocks drive the prediction.

Best LLM per task (post-exclusion leaderboard, val balanced accuracy; full fine-tune):
  T1d -> answerdotai/ModernBERT-large            (T1d_pmci_smci)
  T1e -> answerdotai/ModernBERT-base             (T1e_scn_pcn)
  T2  -> thomas-sounack/BioClinical-ModernBERT-large (T2_multiclass)

Requires the saved checkpoints from 01_refit_best_llms.sh (the sweep did NOT save weights). We load with
attn_implementation="eager" FORCED — ModernBERT under FlashAttention-2 returns NO attentions, and SDPA is
numerically unstable for this model; eager is mathematically identical and returns attentions.

Method: Abnar-Zuidema attention rollout — mean attention over heads per layer, add identity (residual),
row-normalise, matmul across layers; take the [CLS] row as per-token importance. Aggregate sub-tokens ->
words (offset mapping) and words -> clinical blocks (section anchors in the generated summary).

Outputs (per task) -> explainability/<T1d|T1e|T2>/:
  attention_<class>_<rank>_<PID>.png   - clinical text highlighted by word importance + per-block bar
  attention_block_summary.png          - mean per-block importance by class (correct test predictions)
  attention_words_<PID>.csv / attention_block_summary.csv

Run:  conda run -n clinical python explainability/02_attention_rollout.py \
        --ckpt_root D:\ADNI_BIDS_project\derivatives\clinical\encoder_explain_checkpoints
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors

HERE = Path(__file__).resolve().parent
DATA_DIR = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified_post_exclusion"
                r"\verbose\baseline")
CKPT_ROOT_DEFAULT = Path(r"D:\ADNI_BIDS_project\derivatives\clinical\encoder_explain_checkpoints")
TEXT_COL = "Generated_Text"
MAX_LEN = 1024

# task subdir -> (model_id, task_key, label_col, label_map, class names). Mirrors 03_encoder_finetune.
TASKS = {
    "T1d": dict(model_id="answerdotai/ModernBERT-large", task_key="T1d_pmci_smci",
                label_col="conversion_group", label_map={"sMCI": 0, "pMCI": 1},
                classes=["sMCI", "pMCI"], title="T1d - pMCI vs sMCI"),
    "T1e": dict(model_id="answerdotai/ModernBERT-base", task_key="T1e_scn_pcn",
                label_col="conversion_group", label_map={"sCN": 0, "pCN_to_MCI": 1, "pCN_to_AD": 1},
                classes=["sCN", "pCN"], title="T1e - sCN vs pCN"),
    "T2":  dict(model_id="thomas-sounack/BioClinical-ModernBERT-large", task_key="T2_multiclass",
                label_col="Label_bl_multi", label_map={"CN": 0, "MCI": 1, "AD": 2},
                classes=["CN", "MCI", "AD"], title="T2 - CN / MCI / AD"),
}

# Clinical-section anchors (case-insensitive substrings) in the order they appear in the generated
# summary (format_summarised_text). The char span from one anchor to the next is one "block".
# Text before the first anchor is Demographics. Anchors absent in a given summary are skipped.
BLOCK_ANCHORS = [
    ("APOE",            "APOE4 allele status"),
    ("Plasma",          "Plasma "),
    ("Depression/GDS",  "Geriatric Depression"),
    ("Hachinski",       "Hachinski"),
    ("ADAS-Cog13",      "Alzheimer's Disease Assessment"),
    ("CDR",             "Clinical Dementia Rating"),
    ("MMSE",            "Mini-Mental"),
    ("MoCA",            "Montreal Cognitive"),
    ("Clock Drawing",   "Clock Drawing"),
    ("Category Fluency", "Category Fluency"),
    ("Logical Memory",  "Logical Memory"),
    ("AVLT",            "Auditory Verbal Learning"),
    ("Trails",          "Trail Making"),
    ("FAQ",             "Functional Activities"),
]


# ── model loading (FORCE eager so attentions are returned) ────────────────────
def load_eager(ckpt_path):
    import torch
    from transformers import AutoModelForSequenceClassification, AutoTokenizer
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    tok = AutoTokenizer.from_pretrained(str(ckpt_path))
    model = AutoModelForSequenceClassification.from_pretrained(
        str(ckpt_path), attn_implementation="eager", reference_compile=False)
    model.to(dev).eval()
    return model, tok, dev


# ── attention rollout ─────────────────────────────────────────────────────────
def rollout_cls(attentions):
    """attentions: tuple of L tensors [1,H,T,T] -> per-token importance to CLS [T] (numpy)."""
    import torch
    R = None
    for A in attentions:
        a = A[0].mean(0)                          # [T,T] mean over heads
        a = a + torch.eye(a.shape[0], device=a.device)
        a = a / a.sum(-1, keepdim=True)
        R = a if R is None else a @ R
    cls = R[0].detach().float().cpu().numpy()     # CLS row: attention paid to each token
    cls[0] = 0.0                                  # drop CLS self
    return cls


def words_from_offsets(text, offsets, special_mask, tok_imp):
    """Group sub-token importances into whitespace-delimited words via char offsets. Each sub-token's
    importance is added to the word whose char-span contains the token's start. Returns a list of
    [word, start, end, importance] in text order."""
    spans = [(m.start(), m.end(), m.group()) for m in re.finditer(r"\S+", text)]
    imp = np.zeros(len(spans), float)
    for (s, e), is_sp, ti in zip(offsets, special_mask, tok_imp):
        if is_sp or e <= s:
            continue
        for wi, (ws, we, _) in enumerate(spans):
            if ws <= s < we:
                imp[wi] += float(ti)
                break
    return [[spans[i][2], spans[i][0], spans[i][1], float(imp[i])] for i in range(len(spans))]


def block_of(char_start, boundaries):
    """boundaries: sorted list of (start_char, name). Return the block name for a char position."""
    name = "Demographics"
    for st, nm in boundaries:
        if char_start >= st:
            name = nm
        else:
            break
    return name


def find_boundaries(text):
    low = text.lower()
    bs = []
    for name, anchor in BLOCK_ANCHORS:
        i = low.find(anchor.lower())
        if i >= 0:
            bs.append((i, name))
    bs.sort()
    return bs


# ── rendering ────────────────────────────────────────────────────────────────
def render_text_overlay(text, words, block_imp, title, out_png):
    """Highlight each word by its (normalised) importance; show a per-block importance bar beside it."""
    imp = np.array([w[3] for w in words], float)
    if imp.max() > 0:
        norm = imp / imp.max()
    else:
        norm = imp
    cmap = plt.get_cmap("OrRd")

    fig = plt.figure(figsize=(16, 9.5))
    # 2x2: text (top-left) + a thin horizontal colorbar in its OWN row beneath the text (bottom-left),
    # with the per-block bar spanning the full height on the right. This keeps the colorbar out from
    # between the panels (its old placement overlapped the block bar) and gives the bar labels room.
    gs = fig.add_gridspec(2, 2, width_ratios=[3.0, 1.05], height_ratios=[24, 1],
                          wspace=0.26, hspace=0.08)
    axt = fig.add_subplot(gs[0, 0]); axt.axis("off"); axt.set_xlim(0, 1); axt.set_ylim(0, 1)
    axt.set_title(title, fontsize=11, fontweight="bold", loc="left")
    cax = fig.add_subplot(gs[1, 0])                 # dedicated colorbar axis (horizontal, under the text)
    axb = fig.add_subplot(gs[:, 1])                 # per-block bar, full height, own column

    x, y = 0.01, 0.97
    line_h, max_x = 0.022, 0.99
    fig.canvas.draw()
    for w, nrm in zip(words, norm):
        word = w[0]
        chunk = word + " "
        wlen = 0.0085 * len(chunk)
        if x + wlen > max_x:
            x = 0.01; y -= line_h
        if y < 0.02:
            break
        bg = cmap(0.15 + 0.85 * float(nrm)) if nrm > 0.04 else (1, 1, 1, 0)
        axt.text(x, y, word, fontsize=8.0, ha="left", va="top",
                 bbox=dict(boxstyle="round,pad=0.05", fc=bg, ec="none"))
        x += wlen

    # block bar (own column, full height)
    bi = pd.Series(block_imp).sort_values()
    axb.barh(bi.index, bi.values, color="#c0392b")
    axb.set_title("per-block attention", fontsize=9)
    axb.tick_params(axis="y", labelsize=7.5); axb.tick_params(axis="x", labelsize=7)

    # colorbar in its dedicated axis (no longer overlapping the block bar)
    sm = cm.ScalarMappable(cmap=cmap, norm=colors.Normalize(0, 1))
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("word attention (rollout, normalised)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.savefig(out_png, dpi=200, bbox_inches="tight"); plt.close(fig)


def main():
    import torch
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt_root", type=str, default=str(CKPT_ROOT_DEFAULT))
    ap.add_argument("--data_dir", type=str, default=str(DATA_DIR),
                    help="dir CONTAINING seed_*/test.csv (verbose baseline); override on HPC")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n_per_class", type=int, default=3, help="correct exemplars per class")
    ap.add_argument("--n_wrong", type=int, default=2, help="high-confidence misclassifications")
    args = ap.parse_args()
    ckpt_root = Path(args.ckpt_root)

    for sub, cfg in TASKS.items():
        out = HERE / sub; out.mkdir(parents=True, exist_ok=True)
        slug = cfg["model_id"].split("/")[-1]
        ckpt = ckpt_root / slug / cfg["task_key"] / f"seed_{args.seed}" / "full_ft" / "best_checkpoint"
        if not ckpt.exists():
            print(f"[skip] {sub}: checkpoint not found at {ckpt}\n  -> run 01_refit_best_llms.sh first.")
            continue
        print(f"\n{'='*70}\n  {cfg['title']}  ({slug})\n{'='*70}")
        model, tok, dev = load_eager(ckpt)

        df = pd.read_csv(Path(args.data_dir) / f"seed_{args.seed}" / "test.csv", low_memory=False)
        df[cfg["label_col"]] = df[cfg["label_col"]].map(cfg["label_map"])
        df = df.dropna(subset=[cfg["label_col"], TEXT_COL]).reset_index(drop=True)
        df[cfg["label_col"]] = df[cfg["label_col"]].astype(int)
        pid_col = "Patient_ID" if "Patient_ID" in df.columns else df.columns[0]

        # forward all rows, collect prob/pred + rollout per row
        recs, block_rows = [], []
        for i, row in df.iterrows():
            text = str(row[TEXT_COL])
            enc = tok(text, truncation=True, max_length=MAX_LEN, return_tensors="pt",
                      return_offsets_mapping=True)
            offsets = enc.pop("offset_mapping")[0].tolist()
            special = tok.get_special_tokens_mask(enc["input_ids"][0].tolist(),
                                                  already_has_special_tokens=True)
            enc = {k: v.to(dev) for k, v in enc.items()}
            with torch.no_grad():
                o = model(**enc, output_attentions=True)
            prob = torch.softmax(o.logits, -1)[0].detach().float().cpu().numpy()
            pred = int(prob.argmax()); y = int(row[cfg["label_col"]])
            tok_imp = rollout_cls(o.attentions)
            words = words_from_offsets(text, offsets, special, tok_imp)
            bnds = find_boundaries(text)
            bimp = {}
            for w in words:
                b = block_of(w[1], bnds); bimp[b] = bimp.get(b, 0.0) + w[3]
            tot = sum(bimp.values()) or 1.0
            bimp = {k: v / tot for k, v in bimp.items()}                 # fraction of attention per block
            recs.append(dict(idx=i, pid=str(row[pid_col]), y=y, pred=pred,
                             conf=float(prob[pred]), correct=(pred == y),
                             words=words, text=text, bimp=bimp))
            r = {"pid": str(row[pid_col]), "y": cfg["classes"][y], "pred": cfg["classes"][pred],
                 "correct": pred == y}
            r.update(bimp); block_rows.append(r)
        print(f"  n_test={len(recs)}  acc={np.mean([r['correct'] for r in recs]):.3f}")

        # exemplars: top-n correct by confidence per class + n_wrong most-confident errors
        for ci, cname in enumerate(cfg["classes"]):
            corr = sorted([r for r in recs if r["correct"] and r["y"] == ci],
                          key=lambda r: -r["conf"])[:args.n_per_class]
            for rank, r in enumerate(corr):
                ttl = (f"{cfg['title']}  |  {cname} (correct, p={r['conf']:.2f})  |  {r['pid']}  "
                       f"|  attention rollout - {slug}")
                render_text_overlay(r["text"], r["words"], r["bimp"], ttl,
                                    out / f"attention_{cname}_correct{rank}_{r['pid']}.png")
                pd.DataFrame([(w[0], w[3]) for w in r["words"]], columns=["word", "rollout"]) \
                    .to_csv(out / f"attention_words_{cname}_correct{rank}_{r['pid']}.csv", index=False)
        wrong = sorted([r for r in recs if not r["correct"]], key=lambda r: -r["conf"])[:args.n_wrong]
        for rank, r in enumerate(wrong):
            ttl = (f"{cfg['title']}  |  pred {cfg['classes'][r['pred']]} but true "
                   f"{cfg['classes'][r['y']]} (p={r['conf']:.2f})  |  {r['pid']}  |  rollout - {slug}")
            render_text_overlay(r["text"], r["words"], r["bimp"], ttl,
                                out / f"attention_wrong{rank}_{r['pid']}.png")

        # per-task block summary (correct predictions, mean per-block fraction by true class)
        bdf = pd.DataFrame(block_rows)
        bdf.to_csv(out / "attention_block_summary.csv", index=False)
        block_cols = [c for c in bdf.columns if c not in ("pid", "y", "pred", "correct")]
        summ = bdf[bdf["correct"]].groupby("y")[block_cols].mean()
        if len(summ):
            ax = summ.T.plot(kind="barh", figsize=(8.5, 6.5), width=0.8)
            ax.set_xlabel("mean attention fraction (correct test predictions)")
            ax.set_title(f"{cfg['title']} - per-block attention by class ({slug})", fontsize=10)
            ax.legend(title="true class", fontsize=8)
            plt.tight_layout(); plt.savefig(out / "attention_block_summary.png", dpi=200); plt.close()
        print(f"  -> {out}")

    print("\nDone. Attention overlays written to explainability/<task>/.")


if __name__ == "__main__":
    main()
