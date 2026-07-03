"""Synthetic MPRA: determinism, planted-signal correlation, provenance -- all offline."""

from __future__ import annotations

import numpy as np
import pytest

from regmodel.data import (
    ACTIVATOR_MOTIF,
    make_labeled_sequence,
    make_synthetic_mpra,
)


def test_deterministic_for_a_seed():
    a = make_synthetic_mpra(n=100, length=80, seed=7)
    b = make_synthetic_mpra(n=100, length=80, seed=7)
    assert a.sequences == b.sequences
    assert np.allclose(a.activities, b.activities)


def test_different_seeds_differ():
    a = make_synthetic_mpra(n=100, length=80, seed=1)
    b = make_synthetic_mpra(n=100, length=80, seed=2)
    assert a.sequences != b.sequences


def test_activator_presence_raises_activity():
    ds = make_synthetic_mpra(n=800, length=120, seed=3)
    with_motif = [ds.activities[i] for i, s in enumerate(ds.activator_spans) if s]
    without = [ds.activities[i] for i, s in enumerate(ds.activator_spans) if not s]
    # The planted activator must make the label correlate with its presence.
    assert np.mean(with_motif) > np.mean(without) + 1.0


def test_spans_point_at_the_planted_motif():
    ds = make_synthetic_mpra(n=200, length=120, seed=5)
    for i, spans in enumerate(ds.activator_spans):
        for start, end in spans:
            assert ds.sequences[i][start:end] == ACTIVATOR_MOTIF


def test_provenance_is_complete():
    ds = make_synthetic_mpra(n=10, length=50, seed=0)
    prov = ds.provenance
    for key in ("source", "seed", "n", "seq_length", "activator_motif", "noise_sd"):
        assert key in prov
    assert prov["research_use_only"] is True


def test_labeled_sequence_has_motif_at_reported_span():
    seq, (start, end) = make_labeled_sequence(length=100, seed=42)
    assert seq[start:end] == ACTIVATOR_MOTIF


def test_too_short_length_raises_clear_error_not_numpy_high():
    # Regression: length shorter than a motif used to surface a cryptic numpy
    # "high <= 0" from rng.integers. It must now raise a clear, guiding ValueError.
    with pytest.raises(ValueError, match="too short"):
        make_synthetic_mpra(n=4, length=4, seed=0)
    with pytest.raises(ValueError, match="too short"):
        make_labeled_sequence(length=4, seed=0)
