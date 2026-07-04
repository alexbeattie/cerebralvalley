"""CLI: rank a trio's variants, integrating non-coding candidates with the rest.

    python -m trio_prioritizer.cli demo
        Fully offline synthetic trio -> ranked candidates. No network, no key, no torch.
        Recovers the planted de novo + comp-het pair at the top; writes candidates.json.

    python -m trio_prioritizer.cli run --table trio.tsv
        Rank a real TSV of trio variants. Uses the live AlphaGenome scorer when
        ALPHAGENOME_API_KEY is set, otherwise degrades cleanly (ranks on inheritance).

Columns for `run` (tab-separated, one header row):
    chrom  pos  ref  alt  gene  consequence  proband_gt  mother_gt  father_gt
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from typing import Optional

from variant_curator.clients.alphagenome import fetch_alphagenome
from variant_curator.models import VariantConsequence

from . import RESEARCH_USE_CAVEAT
from .data import build_synthetic_trio, causal_variant_ids, offline_fake_scorer
from .models import Genotype, PrioritizedCandidate, TrioVariant
from .prioritize import prioritize
from .scoring import Scorer

TSV_COLUMNS = ("chrom", "pos", "ref", "alt", "gene", "consequence",
               "proband_gt", "mother_gt", "father_gt")


def _parse_consequence(text: str) -> VariantConsequence:
    """Map a VEP/SO term to VariantConsequence; unknown terms fall back to OTHER."""
    try:
        return VariantConsequence(text.strip())
    except ValueError:
        # Accept the enum name too (e.g. "INTRON"), else abstain to OTHER.
        try:
            return VariantConsequence[text.strip().upper()]
        except KeyError:
            return VariantConsequence.OTHER


def load_table(path: str) -> list[TrioVariant]:
    """Read a TSV of trio variants into TrioVariant records."""
    variants: list[TrioVariant] = []
    with open(path, encoding="utf-8") as fh:
        header: Optional[list[str]] = None
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = [c.strip().lower() for c in fields]
                missing = [c for c in TSV_COLUMNS if c not in header]
                if missing:
                    raise ValueError(f"{path}: missing columns {missing}; expected {TSV_COLUMNS}")
                continue
            row = dict(zip(header, fields))
            variants.append(
                TrioVariant(
                    chrom=row["chrom"].strip(),
                    pos=int(row["pos"]),
                    ref=row["ref"].strip(),
                    alt=row["alt"].strip(),
                    gene=row["gene"].strip(),
                    consequence=_parse_consequence(row["consequence"]),
                    proband=Genotype.parse(row["proband_gt"]),
                    mother=Genotype.parse(row["mother_gt"]),
                    father=Genotype.parse(row["father_gt"]),
                )
            )
    if not variants:
        raise ValueError(f"{path}: no variant rows found")
    return variants


def _candidate_to_dict(rank: int, cand: PrioritizedCandidate) -> dict:
    return {
        "rank": rank,
        "gene": cand.gene,
        "inheritance_mode": cand.inheritance.mode.value,
        "parent_of_origin": cand.inheritance.parent_of_origin,
        "confidence": cand.inheritance.confidence,
        "combined_score": cand.combined_score,
        "rationale": cand.rationale,
        "variants": [
            {
                "variant_id": v.variant_id,
                "gene": v.gene,
                "consequence": v.consequence.value,
                "is_coding": v.is_coding,
                "genotypes": {
                    "proband": v.proband.value,
                    "mother": v.mother.value,
                    "father": v.father.value,
                },
                "noncoding_score": {
                    "available": s.available,
                    "scorer": s.scorer,
                    "magnitude": s.magnitude,
                    "top_effect": s.top_effect,
                    "reason": s.reason,
                    "secondary_scorer": s.secondary_scorer,
                    "secondary_magnitude": s.secondary_magnitude,
                },
            }
            for v, s in zip(cand.variants, cand.scores)
        ],
        "sources": [asdict(src) for src in cand.sources()],
        "research_use_only": True,
    }


def _print_candidate(rank: int, cand: PrioritizedCandidate) -> None:
    call = cand.inheritance
    phase = f" — {call.parent_of_origin}" if call.parent_of_origin else ""
    kind = "COMPOUND-HET PAIR" if cand.is_pair else call.mode.value.replace("_", " ").upper()
    print(f"\n#{rank}  [{cand.combined_score:.3f}]  {cand.gene}  {kind}{phase}")
    print(f"    inheritance : {call.mode.value}  ({call.confidence})")
    print(f"    rationale   : {call.rationale}")
    for v, s in zip(cand.variants, cand.scores):
        origin = ""  # phase per allele for the pair
        if cand.is_pair:
            origin = "  <paternal>" if v is cand.variants[0] else "  <maternal>"
        print(f"    variant     : {v.label}  "
              f"GT(p/m/f)={v.proband.value}/{v.mother.value}/{v.father.value}{origin}")
        if s.available:
            print(f"        non-coding: {s.scorer} magnitude={s.magnitude:.2f} — {s.top_effect}")
            if s.source:
                print(f"        source    : {s.source.name} | {s.source.url} | {s.source.detail}")
        elif s.scorer == "n/a":
            print("        non-coding: skipped (coding allele; protein path handles it)")
        else:
            print(f"        non-coding: unavailable — {s.reason}")
            if s.source:
                print(f"        source    : {s.source.name} | {s.source.url}")
        if s.secondary_magnitude is not None:
            print(f"        secondary : {s.secondary_scorer} |Δ|={s.secondary_magnitude:.3f}")


def _run_pipeline(variants: list[TrioVariant], scorer: Scorer, *, out_path: Optional[str],
                  title: str) -> list[PrioritizedCandidate]:
    print("=" * 78)
    print(f"  {title}")
    print("=" * 78)
    print(f"  {RESEARCH_USE_CAVEAT}")
    ranked = prioritize(variants, scorer)
    for i, cand in enumerate(ranked, start=1):
        _print_candidate(i, cand)
    if out_path:
        payload = {
            "research_use_only": True,
            "caveat": RESEARCH_USE_CAVEAT,
            "candidates": [_candidate_to_dict(i, c) for i, c in enumerate(ranked, start=1)],
        }
        with open(out_path, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)
        print(f"\n  wrote {out_path} ({len(ranked)} candidates, with provenance)")
    return ranked


def cmd_demo(args: argparse.Namespace) -> int:
    """Offline synthetic-trio demo: prove the planted answer is recovered at the top."""
    variants = build_synthetic_trio()
    ranked = _run_pipeline(
        variants, offline_fake_scorer, out_path=args.out,
        title="trio_prioritizer demo — synthetic trio (OFFLINE, illustrative scores)",
    )

    # Self-check: the causal de novo and both comp-het alleles must be at the top.
    causal = causal_variant_ids()
    top_ids = {v.variant_id for cand in ranked[:2] for v in cand.variants}
    recovered = causal.issubset(top_ids)
    print("\n" + "-" * 78)
    print("  ANSWER-KEY CHECK")
    print(f"    causal variants : {sorted(causal)}")
    print(f"    top-2 candidates: {sorted(top_ids)}")
    print(f"    recovered at top: {'YES ✓' if recovered else 'NO ✗'}")
    print("-" * 78)
    return 0 if recovered else 1


def cmd_run(args: argparse.Namespace) -> int:
    """Rank a real TSV; live AlphaGenome when keyed, else a clean degrade."""
    variants = load_table(args.table)
    import os

    keyed = bool(os.environ.get("ALPHAGENOME_API_KEY"))
    note = ("live AlphaGenome scorer" if keyed
            else "no ALPHAGENOME_API_KEY — non-coding scores degrade; ranking uses inheritance")
    _run_pipeline(
        variants, fetch_alphagenome, out_path=args.out,
        title=f"trio_prioritizer run — {args.table} ({note})",
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="trio_prioritizer.cli",
        description="Inheritance-aware prioritization of non-coding variants in trios.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    d = sub.add_parser("demo", help="offline synthetic-trio demo (no network/key)")
    d.add_argument("--out", default="candidates.json", help="write ranked candidates + provenance")
    d.set_defaults(func=cmd_demo)

    r = sub.add_parser("run", help="rank a TSV of trio variants")
    r.add_argument("--table", required=True, help="TSV: " + ",".join(TSV_COLUMNS))
    r.add_argument("--out", default="candidates.json", help="write ranked candidates + provenance")
    r.set_defaults(func=cmd_run)
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
