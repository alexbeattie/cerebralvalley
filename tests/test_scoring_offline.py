"""Scoring degrades cleanly with no key; Source is attached when a variant is scored."""

from __future__ import annotations

from trio_prioritizer.data import offline_fake_scorer
from trio_prioritizer.models import Genotype, TrioVariant
from trio_prioritizer.scoring import score_variant
from variant_curator.clients import alphagenome
from variant_curator.models import VariantConsequence as VC


def _v(gene, cons, pos=100, chrom="2", ref="G", alt="A") -> TrioVariant:
    return TrioVariant(
        chrom=chrom, pos=pos, ref=ref, alt=alt, gene=gene, consequence=cons,
        proband=Genotype.HET, mother=Genotype.HOM_REF, father=Genotype.HOM_REF,
    )


def test_degrades_without_key_and_never_touches_network(monkeypatch):
    monkeypatch.delenv(alphagenome.API_KEY_ENV, raising=False)
    # Real scorer with no key must not call the SDK (_score) at all.
    called = {"n": 0}
    monkeypatch.setattr(alphagenome, "_score", lambda *a, **k: called.__setitem__("n", 1))
    score = score_variant(_v("SCN1A", VC.INTRON))
    assert called["n"] == 0
    assert score.available is False
    assert alphagenome.API_KEY_ENV in score.reason
    assert score.source is not None  # provenance carried even on abstention


def test_coding_variant_skips_alphagenome():
    score = score_variant(_v("PAH", VC.MISSENSE))
    assert score.available is False
    assert score.scorer == "n/a"
    assert "coding" in score.reason


def test_offline_scorer_attaches_source_and_magnitude():
    # Use the SCN1A de novo coordinates the synthetic scenario plants a strong score for.
    v = TrioVariant(
        chrom="2", pos=165_991_200, ref="G", alt="A", gene="SCN1A", consequence=VC.INTRON,
        proband=Genotype.HET, mother=Genotype.HOM_REF, father=Genotype.HOM_REF,
    )
    score = score_variant(v, offline_fake_scorer)
    assert score.available is True
    assert score.scorer == "AlphaGenome"
    assert score.magnitude == 0.93
    assert "SCN1A" in score.top_effect
    assert score.source is not None
    assert "ILLUSTRATIVE" in score.source.name  # synthetic origin is honest end to end


def test_offline_scorer_abstains_for_unknown_variant():
    score = score_variant(_v("SCN1A", VC.INTRON, pos=999), offline_fake_scorer)
    assert score.available is False
    assert score.magnitude is None
