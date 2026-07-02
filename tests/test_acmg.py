"""ACMG mapper: supporting-only PP3/BP4, each carrying an AlphaGenome Source."""

from __future__ import annotations

from variant_curator.acmg import BP4_QUANTILE, PP3_QUANTILE, map_alphagenome_to_acmg
from variant_curator.models import (
    AlphaGenomeEvidence,
    RegulatorySignal,
    Source,
    SplicingSignal,
)

SOURCE = Source(name="AlphaGenome (Google DeepMind)", url="https://deepmind.google.com/science/alphagenome", retrieved_at="2026-07-02T00:00:00+00:00", detail="chr7:1:C>T")


def _splice(q):
    return SplicingSignal(
        scorer="SPLICE_SITES", output_type="SPLICE_SITES", top_gene="CFTR",
        max_quantile=q, max_raw=q * 3, direction=None, interpretation="x", source=SOURCE,
    )


def _reg(q):
    return RegulatorySignal(
        modality="RNA_SEQ", output_type="RNA_SEQ", top_tissue="Lung", ontology_curie="UBERON:0002048",
        quantile_score=q, raw_score=q * 2, direction="down", interpretation="x", source=SOURCE,
    )


def test_high_splicing_emits_supporting_pp3_with_source():
    ev = AlphaGenomeEvidence(found=True, splicing=[_splice(PP3_QUANTILE + 0.05)], source=SOURCE)
    hits = map_alphagenome_to_acmg(ev)
    assert [h.code for h in hits] == ["PP3"]
    pp3 = hits[0]
    assert pp3.strength == "supporting"
    assert pp3.direction == "pathogenic"
    assert pp3.source is SOURCE
    assert "not clinically validated" in pp3.caveat


def test_high_regulatory_emits_pp3():
    ev = AlphaGenomeEvidence(found=True, regulatory=[_reg(-(PP3_QUANTILE + 0.02))], source=SOURCE)
    hits = map_alphagenome_to_acmg(ev)
    assert [h.code for h in hits] == ["PP3"]
    assert hits[0].source is SOURCE


def test_all_low_emits_supporting_bp4_with_source():
    ev = AlphaGenomeEvidence(
        found=True,
        splicing=[_splice(BP4_QUANTILE - 0.4)],
        regulatory=[_reg(BP4_QUANTILE - 0.3)],
        source=SOURCE,
    )
    hits = map_alphagenome_to_acmg(ev)
    assert [h.code for h in hits] == ["BP4"]
    bp4 = hits[0]
    assert bp4.direction == "benign"
    assert bp4.source is SOURCE


def test_grey_zone_abstains():
    # Between BP4 and PP3 thresholds: no criterion, don't guess.
    mid = (PP3_QUANTILE + BP4_QUANTILE) / 2
    ev = AlphaGenomeEvidence(found=True, splicing=[_splice(mid)], source=SOURCE)
    assert map_alphagenome_to_acmg(ev) == []


def test_not_found_abstains():
    ev = AlphaGenomeEvidence(found=False, reason="key unset", source=SOURCE)
    assert map_alphagenome_to_acmg(ev) == []
