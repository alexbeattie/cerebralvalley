"""Integrate inheritance + non-coding score into one ranked candidate list.

This is the "integrate the non-coding candidates with the other potentially-positive
variants" step. Compound-het pairs, de novo singletons, and homozygous-recessive
candidates are interleaved into a single ranking so a deep-intronic second hit is weighed
directly against everything else on the table.

The combined score is deliberately simple, transparent, and documented: a genotype-derived
inheritance **prior** blended with the non-coding **magnitude**. The weights are provisional
constants (named below) and need calibration against real cases before any clinical read.
"""

from __future__ import annotations

from typing import Optional

from variant_curator.clients.alphagenome import fetch_alphagenome

from .inheritance import classify_single, find_compound_het_pairs
from .models import InheritanceMode, PrioritizedCandidate, TrioVariant
from .scoring import Scorer, score_variant

SCORER_LABEL = "AlphaGenome"

# --- Provisional weights (need calibration; do not read clinically) ---------------------
# Inheritance dominates ranking because genotype pattern is the strongest signal for a
# severe-disease trio; the non-coding magnitude then orders candidates within a mode.
W_INHERITANCE = 0.6
W_NONCODING = 0.4

# Inheritance priors in [0, 1]. A de novo in a severe-disease proband and a completed
# comp-het pair rank high; a lone inherited het (second hit still missing) ranks low.
INHERITANCE_PRIOR: dict[InheritanceMode, float] = {
    InheritanceMode.DE_NOVO: 1.0,
    InheritanceMode.COMPOUND_HET: 0.9,
    InheritanceMode.HOMOZYGOUS_RECESSIVE: 0.85,
    InheritanceMode.INHERITED_DOMINANT: 0.3,
    InheritanceMode.UNKNOWN: 0.1,
}


def _clamp01(x: float) -> float:
    return 0.0 if x < 0 else 1.0 if x > 1 else x


def combined_score(mode: InheritanceMode, magnitude: Optional[float]) -> float:
    """Blend the inheritance prior with the non-coding magnitude.

    When the magnitude is unavailable (scorer degraded, or a coding allele) the non-coding
    term contributes 0 — ranking still runs on the inheritance prior, so we never fabricate
    a score to fill the gap.
    """
    prior = INHERITANCE_PRIOR.get(mode, 0.1)
    mag = _clamp01(magnitude) if magnitude is not None else 0.0
    return round(W_INHERITANCE * prior + W_NONCODING * mag, 4)


def _pair_key(v: TrioVariant) -> tuple:
    return (v.chrom, v.pos, v.ref, v.alt)


def prioritize(
    variants: list[TrioVariant], scorer: Scorer = fetch_alphagenome
) -> list[PrioritizedCandidate]:
    """Rank a trio's variants into one integrated candidate list, highest combined first.

    Compound-het pairs are formed first and their member alleles are consumed, so a pair
    member is never also listed as a lone inherited-dominant candidate. Remaining variants
    are classified individually. Every scored candidate carries its AlphaGenome `Source`.
    """
    candidates: list[PrioritizedCandidate] = []
    consumed: set[tuple] = set()

    # 1) Compound-het pairs (paternal allele, maternal allele).
    for paternal, maternal, call in find_compound_het_pairs(variants):
        consumed.add(_pair_key(paternal))
        consumed.add(_pair_key(maternal))
        scores = [score_variant(paternal, scorer), score_variant(maternal, scorer)]
        cand = PrioritizedCandidate(
            variants=[paternal, maternal],
            inheritance=call,
            scores=scores,
            combined_score=0.0,
        )
        cand.combined_score = combined_score(call.mode, cand.top_magnitude)
        cand.rationale = _rationale(cand)
        candidates.append(cand)

    # 2) Everything else, classified individually.
    for v in variants:
        if _pair_key(v) in consumed:
            continue
        call = classify_single(v)
        score = score_variant(v, scorer)
        cand = PrioritizedCandidate(
            variants=[v],
            inheritance=call,
            scores=[score],
            combined_score=0.0,
        )
        cand.combined_score = combined_score(call.mode, cand.top_magnitude)
        cand.rationale = _rationale(cand)
        candidates.append(cand)

    # Stable, deterministic order: score desc, then gene, then label.
    candidates.sort(key=lambda c: (-c.combined_score, c.gene, c.label))
    return candidates


def _rationale(cand: PrioritizedCandidate) -> str:
    """A one-line human summary: mode + phase + the non-coding evidence, if any."""
    call = cand.inheritance
    mode = call.mode.value.replace("_", " ")
    phase = f" ({call.parent_of_origin})" if call.parent_of_origin else ""
    mag = cand.top_magnitude
    if mag is not None:
        effect = next(
            (s.top_effect for s in cand.scores if s.available and s.top_effect),
            "non-coding signal",
        )
        ev = f"; {SCORER_LABEL} magnitude {mag:.2f} — {effect}"
    elif any(s.available is False and s.scorer != "n/a" for s in cand.scores):
        reason = next((s.reason for s in cand.scores if not s.available), "")
        ev = f"; non-coding score unavailable ({reason})"
    else:
        ev = "; coding allele scored on the protein path elsewhere"
    return f"{cand.gene} {mode}{phase}{ev}"
