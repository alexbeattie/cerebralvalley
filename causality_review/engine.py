"""Causality-scoring engine.

Given a patient's HPO phenotype set and the variants the lab reported, score how
well each variant's gene explains THIS patient, rank them, and surface the two
traps a clinician's intuition falls into:

  1. Overfitting -- we score against an explicit feature set, so a partial match
     reads as "partial", not "close enough".
  2. Partial explanations -- features left unexplained by every reported variant
     are called out, because they may mean an unreported second cause. That is
     the trigger to re-analyze the genome, not a failure of the review.
"""

from __future__ import annotations

import httpx

from .clients.hpo import fetch_gene_phenotypes
from .models import (
    CausalityReport,
    FitTier,
    GenePhenotypeKnowledge,
    PatientProfile,
    ReportedVariant,
    VariantFit,
)

# Fraction of the patient's features a gene must explain to reach each tier.
_TIER_THRESHOLDS = [
    (1.0, FitTier.BEST_FIT),
    (0.6, FitTier.POSSIBLE),
    (0.3, FitTier.PARTIAL),
    (0.0001, FitTier.WEAK),
]


def _tier_for(score: float) -> FitTier:
    for threshold, tier in _TIER_THRESHOLDS:
        if score >= threshold:
            return tier
    return FitTier.UNLIKELY


def _score_one(variant: ReportedVariant, patient: PatientProfile, knowledge: GenePhenotypeKnowledge) -> VariantFit:
    if not knowledge.found:
        return VariantFit(
            variant=variant,
            tier=FitTier.UNLIKELY,
            score=0.0,
            unexplained=list(patient.phenotypes),
            rationale=f"No HPO gene-phenotype knowledge found for {variant.gene}; cannot assess fit.",
            source=knowledge.source,
            knowledge_found=False,
        )

    explained = [p for p in patient.phenotypes if p.hpo_id in knowledge.phenotype_ids]
    unexplained = [p for p in patient.phenotypes if p.hpo_id not in knowledge.phenotype_ids]
    total = len(patient.phenotypes) or 1
    score = len(explained) / total
    tier = _tier_for(score)

    disease_hint = f" ({knowledge.diseases[0]})" if knowledge.diseases else ""
    rationale = (
        f"{variant.gene}{disease_hint} explains {len(explained)}/{total} of the "
        f"patient's features"
    )
    if unexplained:
        rationale += f"; leaves unexplained: {', '.join(p.label for p in unexplained)}."
    else:
        rationale += "; accounts for the full presented picture."

    return VariantFit(
        variant=variant,
        tier=tier,
        score=score,
        explained=explained,
        unexplained=unexplained,
        diseases=knowledge.diseases,
        rationale=rationale,
        source=knowledge.source,
    )


def review_causality(
    client: httpx.Client, patient: PatientProfile, variants: list[ReportedVariant]
) -> CausalityReport:
    fits: list[VariantFit] = []
    for v in variants:
        knowledge = fetch_gene_phenotypes(client, v.gene)
        fits.append(_score_one(v, patient, knowledge))

    # Rank: higher score first, then by lab-classification prominence is left to
    # the clinician -- we sort purely on fit so the ranking stays objective.
    fits.sort(key=lambda f: (f.score, len(f.explained)), reverse=True)

    # Features explained by NO reported variant at all.
    explained_ids = {p.hpo_id for f in fits for p in f.explained}
    residual = [p for p in patient.phenotypes if p.hpo_id not in explained_ids]

    report = CausalityReport(patient=patient, fits=fits, residual_unexplained=residual)
    report.flags = _build_flags(fits, residual, len(patient.phenotypes))
    return report


def _build_flags(fits: list[VariantFit], residual: list, n_features: int) -> list[str]:
    flags: list[str] = []
    top = fits[0] if fits else None

    if top and top.tier == FitTier.BEST_FIT:
        flags.append(
            f"Strong single-variant fit: {top.variant.label} explains the full presented picture."
        )
    elif top and top.explained:
        flags.append(
            f"No single variant explains everything; best candidate ({top.variant.label}) "
            f"explains {len(top.explained)}/{n_features}."
        )

    if residual:
        labels = ", ".join(p.label for p in residual)
        flags.append(
            f"{len(residual)} feature(s) explained by NO reported variant ({labels}). "
            "Consider an unreported extension of a listed gene's phenotype, a second "
            "independent cause (~5% of solved cases), or genome re-analysis to chase them down."
        )

    for f in fits:
        if not f.knowledge_found:
            flags.append(f"Could not retrieve HPO knowledge for {f.variant.gene}; its fit is unscored.")

    return flags
