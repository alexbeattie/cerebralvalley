"""AlphaGenome client: graceful degradation and tidy_scores reduction, all offline."""

from __future__ import annotations

from unittest import mock

from variant_curator.clients import alphagenome
from variant_curator.models import AlphaGenomeEvidence

from . import fixtures


def test_import_variant_curator_needs_no_sdk_or_key():
    # Importing the package must never require the alphagenome SDK or a key.
    import importlib

    import variant_curator

    importlib.reload(variant_curator)


def test_missing_key_returns_found_false_with_reason(monkeypatch):
    monkeypatch.delenv(alphagenome.API_KEY_ENV, raising=False)
    with mock.patch.object(alphagenome, "_score") as scored:
        ev = alphagenome.fetch_alphagenome("7", 117559593, "C", "T", "CFTR")
    scored.assert_not_called()  # no network / SDK contact
    assert isinstance(ev, AlphaGenomeEvidence)
    assert ev.found is False
    assert alphagenome.API_KEY_ENV in ev.reason
    assert ev.source is not None and ev.source.url.startswith("https://")
    assert ev.source.detail  # variant id recorded on the source


def test_reduces_tidy_scores_into_signals_and_top_effect():
    with mock.patch.object(
        alphagenome,
        "_score",
        side_effect=[fixtures.splicing_df("CFTR"), fixtures.regulatory_df("CFTR")],
    ):
        ev = alphagenome.fetch_alphagenome("7", 117559593, "C", "T", "CFTR", api_key="fake")

    assert ev.found is True
    assert ev.research_use_only is True

    # Splicing: one signal per scorer, neighbouring gene filtered, ranked by magnitude.
    scorers = {s.scorer for s in ev.splicing}
    assert scorers == {"SPLICE_SITES", "SPLICE_JUNCTIONS", "SPLICE_SITE_USAGE"}
    top_splice = ev.splicing[0]
    assert top_splice.scorer == "SPLICE_SITES"
    assert top_splice.max_quantile == 0.95  # 0.95 beat the 0.40 track and the 0.99 neighbour
    assert top_splice.top_gene == "CFTR"
    assert top_splice.source is not None

    # Regulatory: disease-relevant tissue (lung) wins over the larger liver score.
    rna = next(r for r in ev.regulatory if r.modality == "RNA_SEQ")
    assert rna.top_tissue == "Lung"
    assert rna.quantile_score == -0.70
    assert rna.direction == "down"
    assert rna.ontology_curie == "UBERON:0002048"

    # Positional accessibility track survived (no gene_name).
    assert any(r.modality == "DNASE" for r in ev.regulatory)

    # Top effect is the strongest signal overall — the 0.95 splice-site alteration.
    assert "splice-site" in ev.top_effect
    assert "CFTR" in ev.top_effect


def test_api_error_degrades_to_found_false():
    with mock.patch.object(alphagenome, "_score", side_effect=RuntimeError("boom")):
        ev = alphagenome.fetch_alphagenome("7", 117559593, "C", "T", "CFTR", api_key="fake")
    assert ev.found is False
    assert "boom" in ev.reason
    assert ev.source is not None


def test_import_error_degrades_to_found_false():
    with mock.patch.object(alphagenome, "_score", side_effect=ImportError("no alphagenome")):
        ev = alphagenome.fetch_alphagenome("7", 117559593, "C", "T", "CFTR", api_key="fake")
    assert ev.found is False
    assert "not installed" in ev.reason
