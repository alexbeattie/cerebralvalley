"""The planted causal de novo AND the comp-het pair must rank at the top, offline."""

from __future__ import annotations

from trio_prioritizer.data import build_synthetic_trio, causal_variant_ids, offline_fake_scorer
from trio_prioritizer.models import Genotype, InheritanceMode, TrioVariant
from trio_prioritizer.prioritize import combined_score, prioritize
from variant_curator.models import VariantConsequence as VC


def test_planted_answer_ranks_at_top():
    ranked = prioritize(build_synthetic_trio(), offline_fake_scorer)

    # Top candidate is the de novo SCN1A variant.
    assert ranked[0].inheritance.mode is InheritanceMode.DE_NOVO
    assert ranked[0].gene == "SCN1A"

    # Second is the completed PAH compound-het pair (a coding + a non-coding allele).
    assert ranked[1].inheritance.mode is InheritanceMode.COMPOUND_HET
    assert ranked[1].gene == "PAH"
    assert ranked[1].is_pair
    kinds = {v.is_coding for v in ranked[1].variants}
    assert kinds == {True, False}  # the classic coding + non-coding second-hit pair

    # Both causal answers occupy the top two candidates.
    top_ids = {v.variant_id for cand in ranked[:2] for v in cand.variants}
    assert causal_variant_ids().issubset(top_ids)


def test_lone_het_ranks_below_completed_comphet():
    ranked = prioritize(build_synthetic_trio(), offline_fake_scorer)
    order = [c.gene for c in ranked]
    # GAA is a lone het (one hit, no partner) — must sit below the comp-het pair.
    assert order.index("PAH") < order.index("GAA")


def test_comphet_members_not_double_listed():
    ranked = prioritize(build_synthetic_trio(), offline_fake_scorer)
    pah_candidates = [c for c in ranked if c.gene == "PAH"]
    # Exactly one PAH candidate (the pair); its alleles are not re-listed individually.
    assert len(pah_candidates) == 1
    assert pah_candidates[0].is_pair


def test_ranking_still_runs_when_scores_unavailable(monkeypatch):
    # With every non-coding score unavailable, ranking falls back to inheritance priors:
    # de novo still beats the inherited-dominant distractor. No fabricated scores.
    def degraded(chrom, pos, ref, alt, gene, **kw):
        from variant_curator.models import AlphaGenomeEvidence
        return AlphaGenomeEvidence(found=False, variant_id=f"{chrom}:{pos}", reason="no key")

    ranked = prioritize(build_synthetic_trio(), degraded)
    assert ranked[0].inheritance.mode in (
        InheritanceMode.DE_NOVO, InheritanceMode.COMPOUND_HET, InheritanceMode.HOMOZYGOUS_RECESSIVE
    )
    assert all(c.top_magnitude is None for c in ranked)  # nothing invented


def test_combined_score_blends_prior_and_magnitude():
    de_novo_high = combined_score(InheritanceMode.DE_NOVO, 0.93)
    inherited_high = combined_score(InheritanceMode.INHERITED_DOMINANT, 0.93)
    # Same magnitude, different prior -> de novo must outrank inherited.
    assert de_novo_high > inherited_high
    # Missing magnitude contributes zero, not a guess.
    assert combined_score(InheritanceMode.DE_NOVO, None) == round(0.6 * 1.0, 4)
