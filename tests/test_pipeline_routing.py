"""Pipeline routing: non-coding -> AlphaGenome; missense/nonsense path unchanged."""

from __future__ import annotations

from unittest import mock

from variant_curator import pipeline
from variant_curator.models import (
    AlphaGenomeEvidence,
    ClinvarEvidence,
    GnomadEvidence,
    VariantConsequence,
    VariantInput,
    VepEvidence,
)


def _vep(consequence: VariantConsequence, raw: str) -> VepEvidence:
    return VepEvidence(
        found=True,
        consequence=consequence,
        consequence_raw=raw,
        gene_symbol="CFTR",
        protein_hgvs=None,
        chrom="7",
        pos=117559593,
        ref="C",
        alt="T",
    )


def _patched(vep: VepEvidence):
    """Patch the pipeline's network/SDK boundary; return the mock namespace."""
    ctx = {
        "get_client": mock.patch.object(pipeline, "get_client", return_value=mock.MagicMock()),
        "fetch_vep": mock.patch.object(pipeline, "fetch_vep", return_value=vep),
        "fetch_gnomad": mock.patch.object(
            pipeline, "fetch_gnomad",
            return_value=GnomadEvidence(found=True, variant_id="7-117559593-C-T", dataset="gnomad_r4"),
        ),
        "fetch_clinvar": mock.patch.object(
            pipeline, "fetch_clinvar", return_value=ClinvarEvidence(found=False),
        ),
        "fetch_alphagenome": mock.patch.object(
            pipeline, "fetch_alphagenome",
            return_value=AlphaGenomeEvidence(found=True, variant_id="chr7:117559593:C>T"),
        ),
    }
    return ctx


def test_noncoding_routes_to_alphagenome_and_still_fetches_gnomad_clinvar():
    patches = _patched(_vep(VariantConsequence.INTRON, "intron_variant"))
    with patches["get_client"], patches["fetch_vep"], patches["fetch_gnomad"] as gnomad, patches[
        "fetch_clinvar"
    ] as clinvar, patches["fetch_alphagenome"] as alphagenome:
        bundle = pipeline.assemble_evidence(VariantInput(gene="CFTR", hgvs_c="c.3717+12191C>T"))

    alphagenome.assert_called_once()
    assert bundle.alphagenome is not None and bundle.alphagenome.found is True
    gnomad.assert_called_once()  # frequency still relevant for non-coding
    clinvar.assert_called_once()  # priors still relevant for non-coding
    assert bundle.gnomad is not None
    assert bundle.clinvar is not None


def test_missense_path_unchanged_and_never_calls_alphagenome():
    patches = _patched(_vep(VariantConsequence.MISSENSE, "missense_variant"))
    with patches["get_client"], patches["fetch_vep"], patches["fetch_gnomad"] as gnomad, patches[
        "fetch_clinvar"
    ] as clinvar, patches["fetch_alphagenome"] as alphagenome:
        bundle = pipeline.assemble_evidence(VariantInput(gene="CFTR", hgvs_c="c.1521_1523del"))

    alphagenome.assert_not_called()
    assert bundle.alphagenome is None
    gnomad.assert_called_once()
    clinvar.assert_called_once()


def test_genuinely_out_of_scope_still_warns_and_bails():
    patches = _patched(_vep(VariantConsequence.OTHER, "synonymous_variant"))
    with patches["get_client"], patches["fetch_vep"], patches["fetch_gnomad"] as gnomad, patches[
        "fetch_clinvar"
    ] as clinvar, patches["fetch_alphagenome"] as alphagenome:
        bundle = pipeline.assemble_evidence(VariantInput(gene="CFTR", hgvs_c="c.1408G>A"))

    alphagenome.assert_not_called()
    gnomad.assert_not_called()
    clinvar.assert_not_called()
    assert bundle.alphagenome is None
    assert any("Out of scope" in w for w in bundle.warnings)
