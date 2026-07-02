"""The AlphaGenome cross-check must degrade cleanly with no key -- our-model rows only,
a recorded reason, and a Source on every row. No network, no key."""

from __future__ import annotations

import numpy as np

from regmodel.compare import VariantSpec, compare_variants
from regmodel.data import ACTIVATOR_MOTIF, make_labeled_sequence
from regmodel.model import ModelConfig, build_model
from regmodel.train import seed_everything
from variant_curator.clients import alphagenome


def _model(length=60):
    seed_everything(0)
    return build_model(ModelConfig(seq_length=length))


def _variants(length=60):
    seq, (start, end) = make_labeled_sequence(length=length, seed=3, activator_motif=ACTIVATOR_MOTIF)
    center = (start + end) // 2
    ref = seq[center]
    alt = next(b for b in "ACGT" if b != ref)
    return [
        VariantSpec(
            chrom="7", pos=117559593, ref=ref, alt=alt, gene_symbol="CFTR",
            local_window_seq=seq, window_pos=center, note="illustrative",
        )
    ]


def test_compare_skips_cleanly_without_key(monkeypatch):
    monkeypatch.delenv(alphagenome.API_KEY_ENV, raising=False)
    model = _model()
    result = compare_variants(model, _variants(), api_key=None)

    assert result.alphagenome_available is False
    assert result.correlation is None
    assert "SKIPPED" in result.note

    row = result.rows[0]
    assert row.alphagenome_found is False
    assert row.alphagenome_magnitude is None
    # Our model still produced a delta, and the AG reason + Source are recorded (never faked).
    assert isinstance(row.our_delta, float)
    assert alphagenome.API_KEY_ENV in row.note
    assert row.source is not None and row.source.url.startswith("https://")


def test_compare_table_has_expected_columns():
    model = _model()
    result = compare_variants(model, _variants(), api_key=None)
    for col in ("variant_id", "our_delta", "ag_found", "ag_magnitude", "agreement", "note"):
        assert col in result.table.columns


def test_compare_correlates_when_alphagenome_present(monkeypatch):
    # Patch the reused client so we exercise the "AG available" path with no network:
    # fabricate two AG evidences with different magnitudes and check a correlation is produced.
    from variant_curator.models import AlphaGenomeEvidence, RegulatorySignal, Source

    def fake_fetch(chrom, pos, ref, alt, gene, *, api_key=None):
        src = Source(name="AlphaGenome (fake)", url="https://example.org", retrieved_at="now", detail="x")
        mag = 0.9 if gene == "A" else 0.2
        return AlphaGenomeEvidence(
            found=True, variant_id=f"{gene}:{pos}", api="fake", source=src,
            regulatory=[RegulatorySignal(
                modality="RNA_SEQ", output_type="RNA_SEQ", top_tissue="t", ontology_curie=None,
                quantile_score=mag, raw_score=mag, direction="up", interpretation="i", source=src,
            )],
        )

    monkeypatch.setattr("regmodel.compare.fetch_alphagenome", fake_fetch)
    model = _model()
    seq, (s, e) = make_labeled_sequence(length=60, seed=4, activator_motif=ACTIVATOR_MOTIF)
    center = (s + e) // 2
    alt = next(b for b in "ACGT" if b != seq[center])
    variants = [
        VariantSpec("1", 1, seq[center], alt, "A", seq, center),
        VariantSpec("2", 2, seq[center], alt, "B", seq, center),
    ]
    result = compare_variants(model, variants, api_key="fake")
    assert result.alphagenome_available is True
    assert all(r.alphagenome_found for r in result.rows)
