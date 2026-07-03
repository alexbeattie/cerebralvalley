"""regmodel CLI: train, ism, variant, compare, and a fully-offline `demo`.

    python -m regmodel.cli demo                       # whole synthetic pipeline -> artifacts/
    python -m regmodel.cli train --out artifacts       # train + save model + sidecar
    python -m regmodel.cli ism --seq ACGT...           # ISM map for one sequence
    python -m regmodel.cli variant --seq ACGT... --pos 100 --alt G
    python -m regmodel.cli compare --model artifacts/model.pt   # AlphaGenome cross-check

Every path that shows a number also writes a JSON sidecar recording how it was produced.
`demo`, `ism`, and `variant` need no network and no API key; `compare` skips the AlphaGenome
half cleanly when the key/SDK is absent.
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import asdict

import numpy as np

from .compare import VariantSpec, compare_variants
from .data import ACTIVATOR_MOTIF, make_labeled_sequence, make_synthetic_mpra
from .encoding import BASES
from .ism import (
    importance_track,
    ism_matrix,
    predict_activity,
    span_contains_top,
    top_positions,
    variant_effect,
)
from .model import ModelConfig
from .plots import plot_comparison_scatter, plot_importance_track, plot_ism_heatmap
from .train import TrainConfig, load_model, now_iso, save_model, train_model

DEFAULT_ARTIFACTS = "artifacts"
# Length of the fallback synthetic model trained when no --model is given. The CNN is
# length-agnostic (padded convs + global pooling), so it scores sequences of ANY length at
# ISM/variant time. This is deliberately NOT tied to the user's --seq length: a short query
# sequence must neither shrink the training task nor (when shorter than a motif) break it.
DEFAULT_TRAIN_LENGTH = 200


# --- helpers ---------------------------------------------------------------------------

def _train_synthetic(n: int, length: int, seed: int, fast: bool):
    dataset = make_synthetic_mpra(n=n, length=length, seed=seed)
    tc = TrainConfig.fast(seed=seed) if fast else TrainConfig(seed=seed)
    return train_model(dataset, tc, ModelConfig(seq_length=length)), dataset


def _load_or_train(model_path: str | None, seed: int, *, train_length: int = DEFAULT_TRAIN_LENGTH):
    """Load a saved model, or train a quick synthetic one so single-shot commands work.

    The fallback is trained at `train_length` (a sensible fixed size), independent of whatever
    sequence ISM/variant will later score -- the model's global pooling makes it length-
    agnostic at inference, so the query sequence length never feeds back into training.
    """
    if model_path and os.path.exists(model_path):
        return load_model(model_path), f"loaded {model_path}"
    result, _ = _train_synthetic(n=1500, length=train_length, seed=seed, fast=True)
    return result.model, f"trained a quick synthetic model at L={train_length} (no --model given)"


def _demo_variants(length: int) -> list[VariantSpec]:
    """Illustrative variants for the offline cross-check.

    Genomic coords are real disease-gene loci (so the AlphaGenome call is well-formed when a
    key exists), but the local windows are *synthetic* -- an activator motif planted at the
    center with the variant breaking it -- because we cannot fetch the genome offline. Flagged
    as illustrative in the row note so nothing here is mistaken for a real measurement.
    """
    specs = [
        ("11", 47337421, "MYBPC3"),
        ("7", 117559593, "CFTR"),
        ("3", 38550350, "SCN5A"),
    ]
    out: list[VariantSpec] = []
    for i, (chrom, pos, gene) in enumerate(specs):
        seq, (start, end) = make_labeled_sequence(length=length, seed=1000 + i, activator_motif=ACTIVATOR_MOTIF)
        center = (start + end) // 2  # a base inside the planted motif
        ref = seq[center]
        alt = next(b for b in BASES if b != ref)  # any different base breaks the motif
        out.append(
            VariantSpec(
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                gene_symbol=gene,
                local_window_seq=seq,
                window_pos=center,
                note="illustrative synthetic window (offline)",
            )
        )
    return out


def _write_json(path: str, payload: dict) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, default=str)


# --- subcommands -----------------------------------------------------------------------

def cmd_train(args) -> None:
    result, dataset = _train_synthetic(args.n, args.length, args.seed, args.fast)
    paths = save_model(result, args.out)
    print(f"trained on {dataset.provenance['source']}  (n={args.n}, L={args.length}, seed={args.seed})")
    m = result.metrics
    print(f"held-out test:  Pearson r={m['pearson']:.3f}  Spearman r={m['spearman']:.3f}  MSE={m['mse']:.3f}")
    print(f"wrote model    -> {paths['model']}")
    print(f"wrote sidecar  -> {paths['sidecar']}")


def cmd_ism(args) -> None:
    seq = args.seq.strip().upper()
    model, how = _load_or_train(args.model, args.seed)
    matrix = ism_matrix(model, seq)
    track = importance_track(matrix)
    tops = top_positions(track, args.top)
    print(f"model: {how}")
    print(f"sequence length: {len(seq)}")
    print(f"reference activity: {predict_activity(model, seq):.3f}")
    print(f"top-{args.top} importance positions (0-based): {tops}")
    for i in tops:
        best_alt = BASES[int(np.argmin(matrix[i]))]  # base with the largest activity drop
        print(f"  pos {i:4d} (ref {seq[i]}): mean|Δ|={track[i]:.3f}, strongest loss -> {best_alt} (Δ={matrix[i].min():+.3f})")

    if args.out:
        os.makedirs(args.out, exist_ok=True)
        heat = plot_ism_heatmap(matrix, seq, os.path.join(args.out, "ism_heatmap.png"))
        line = plot_importance_track(track, os.path.join(args.out, "ism_importance.png"))
        _write_json(
            os.path.join(args.out, "ism.json"),
            {
                "created_at": now_iso(),
                "model": how,
                "sequence": seq,
                "reference_activity": predict_activity(model, seq),
                "top_positions": tops,
                "importance_track": track.tolist(),
                "research_use_only": True,
            },
        )
        print(f"wrote {heat}\nwrote {line}")


def cmd_variant(args) -> None:
    seq = args.seq.strip().upper()
    model, how = _load_or_train(args.model, args.seed)
    delta = variant_effect(model, seq, args.pos, args.alt)
    direction = "increase" if delta > 0 else ("decrease" if delta < 0 else "no change")
    print(f"model: {how}")
    print(f"variant: pos {args.pos} {seq[args.pos]}>{args.alt.upper()}")
    print(f"predicted activity delta: {delta:+.4f}  ({direction})")


def cmd_compare(args) -> None:
    length = args.length
    model, how = _load_or_train(args.model, args.seed, train_length=length)
    variants = _demo_variants(length)
    result = compare_variants(model, variants, api_key=args.api_key)

    print(f"model: {how}")
    print(f"AlphaGenome available: {result.alphagenome_available}")
    print(result.note)
    print(result.table.to_string(index=False))
    if result.correlation is not None:
        print(f"Pearson(our |Δ|, AG magnitude) = {result.correlation:.3f}")

    if args.out:
        os.makedirs(args.out, exist_ok=True)
        csv_path = os.path.join(args.out, "comparison.csv")
        result.table.to_csv(csv_path, index=False)
        _write_json(
            os.path.join(args.out, "comparison.json"),
            {
                "created_at": now_iso(),
                "model": how,
                "alphagenome_available": result.alphagenome_available,
                "correlation": result.correlation,
                "note": result.note,
                "rows": [asdict(r) for r in result.rows],
                "research_use_only": True,
            },
        )
        print(f"wrote {csv_path}")
        if result.alphagenome_available:
            ours = np.array([r.our_magnitude for r in result.rows])
            ag = np.array([r.alphagenome_magnitude or 0.0 for r in result.rows])
            scatter = plot_comparison_scatter(
                ours, ag, os.path.join(args.out, "comparison_scatter.png"),
                labels=[r.gene_symbol for r in result.rows], correlation=result.correlation,
            )
            print(f"wrote {scatter}")


def cmd_demo(args) -> None:
    out = args.out
    os.makedirs(out, exist_ok=True)
    length = args.length
    seed = args.seed
    print("=" * 72)
    print("  regmodel offline demo  (synthetic MPRA -> train -> ISM -> cross-check)")
    print("=" * 72)

    # 1) train on synthetic data
    result, dataset = _train_synthetic(args.n, length, seed, fast=args.fast)
    paths = save_model(result, out)
    m = result.metrics
    print(f"\n[1] trained on {dataset.provenance['source']}  (n={args.n}, L={length}, seed={seed})")
    print(f"    held-out test:  Pearson r={m['pearson']:.3f}  Spearman r={m['spearman']:.3f}  MSE={m['mse']:.3f}")
    print(f"    model   -> {paths['model']}")
    print(f"    sidecar -> {paths['sidecar']}")

    # 2) ISM on a held-out sequence with a known planted motif
    seq, span = make_labeled_sequence(length=length, seed=999, activator_motif=ACTIVATOR_MOTIF)
    matrix = ism_matrix(result.model, seq)
    track = importance_track(matrix)
    tops = top_positions(track, 5)
    hit = span_contains_top(track, span, k=3)
    print(f"\n[2] ISM on a held-out sequence; activator motif planted at {span}")
    print(f"    top-5 importance positions: {tops}")
    print(f"    motif recovered in top-3? {hit}")
    heat = plot_ism_heatmap(matrix, seq, os.path.join(out, "ism_heatmap.png"), highlight=span)
    line = plot_importance_track(track, os.path.join(out, "ism_importance.png"), highlight=span)
    print(f"    heatmap    -> {heat}")
    print(f"    importance -> {line}")

    # 3) metrics.json (top-level provenance for the whole demo run)
    metrics_path = os.path.join(out, "metrics.json")
    _write_json(
        metrics_path,
        {
            "created_at": now_iso(),
            "seed": seed,
            "data_provenance": dataset.provenance,
            "heldout_metrics": m,
            "ism_probe": {"planted_span": list(span), "top_positions": tops, "recovered_in_top3": hit},
            "research_use_only": True,
        },
    )
    print(f"    metrics    -> {metrics_path}")

    # 4) AlphaGenome cross-check (skips cleanly offline)
    variants = _demo_variants(length)
    comp = compare_variants(result.model, variants, api_key=args.api_key)
    print("\n[3] AlphaGenome cross-check")
    print(f"    available: {comp.alphagenome_available}")
    print(f"    {comp.note}")
    print(comp.table.to_string(index=False).replace("\n", "\n    "))
    csv_path = os.path.join(out, "comparison.csv")
    comp.table.to_csv(csv_path, index=False)
    _write_json(
        os.path.join(out, "comparison.json"),
        {
            "created_at": now_iso(),
            "alphagenome_available": comp.alphagenome_available,
            "correlation": comp.correlation,
            "note": comp.note,
            "rows": [asdict(r) for r in comp.rows],
            "research_use_only": True,
        },
    )
    print(f"    table      -> {csv_path}")
    if comp.alphagenome_available:
        ours = np.array([r.our_magnitude for r in comp.rows])
        ag = np.array([r.alphagenome_magnitude or 0.0 for r in comp.rows])
        scatter = plot_comparison_scatter(
            ours, ag, os.path.join(out, "comparison_scatter.png"),
            labels=[r.gene_symbol for r in comp.rows], correlation=comp.correlation,
        )
        print(f"    scatter    -> {scatter}")

    print(f"\nDone. All artifacts in ./{out}/  (research use only; synthetic data by default).")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="regmodel: MPRA sequence->activity model with ISM + AlphaGenome cross-check.")
    sub = parser.add_subparsers(dest="command")

    p_train = sub.add_parser("train", help="train on synthetic MPRA and save model + sidecar")
    p_train.add_argument("--n", type=int, default=4000)
    p_train.add_argument("--length", type=int, default=200)
    p_train.add_argument("--seed", type=int, default=0)
    p_train.add_argument("--out", default=DEFAULT_ARTIFACTS)
    p_train.add_argument("--fast", action="store_true", help="few-epoch path (quick smoke run)")
    p_train.set_defaults(func=cmd_train)

    p_ism = sub.add_parser("ism", help="single-base ISM map for a sequence")
    p_ism.add_argument("--seq", required=True)
    p_ism.add_argument("--model", default=None, help="path to a saved model.pt (else trains a quick one)")
    p_ism.add_argument("--seed", type=int, default=0)
    p_ism.add_argument("--top", type=int, default=5)
    p_ism.add_argument("--out", default=None, help="dir to write heatmap/track PNGs + json")
    p_ism.set_defaults(func=cmd_ism)

    p_var = sub.add_parser("variant", help="predicted activity delta for one substitution")
    p_var.add_argument("--seq", required=True)
    p_var.add_argument("--pos", type=int, required=True, help="0-based position")
    p_var.add_argument("--alt", required=True, help="alternate base (A/C/G/T)")
    p_var.add_argument("--model", default=None)
    p_var.add_argument("--seed", type=int, default=0)
    p_var.set_defaults(func=cmd_variant)

    p_cmp = sub.add_parser("compare", help="AlphaGenome cross-check (skips cleanly w/o key)")
    p_cmp.add_argument("--model", default=None)
    p_cmp.add_argument("--length", type=int, default=200)
    p_cmp.add_argument("--seed", type=int, default=0)
    p_cmp.add_argument("--api-key", default=None, help="AlphaGenome key (else uses ALPHAGENOME_API_KEY / skips)")
    p_cmp.add_argument("--out", default=None)
    p_cmp.set_defaults(func=cmd_compare)

    p_demo = sub.add_parser("demo", help="run the whole synthetic pipeline offline -> artifacts/")
    p_demo.add_argument("--n", type=int, default=4000)
    p_demo.add_argument("--length", type=int, default=200)
    p_demo.add_argument("--seed", type=int, default=0)
    p_demo.add_argument("--out", default=DEFAULT_ARTIFACTS)
    p_demo.add_argument("--api-key", default=None)
    p_demo.add_argument("--fast", action="store_true", help="few-epoch path for a quicker demo")
    p_demo.set_defaults(func=cmd_demo)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "command", None):
        # No subcommand: run the offline demo, matching variant_curator's "just works" default.
        print("(no subcommand given; running the offline demo)\n")
        demo_args = parser.parse_args(["demo"])
        demo_args.func(demo_args)
        return
    args.func(args)


if __name__ == "__main__":
    main()
