"""ISM matrix shape/invariants and single-variant query -- offline, on a tiny random model."""

from __future__ import annotations

import numpy as np
import pytest

from regmodel.encoding import BASE_TO_IDX
from regmodel.ism import (
    importance_track,
    ism_matrix,
    predict_activity,
    top_positions,
    variant_effect,
)
from regmodel.model import ModelConfig, build_model
from regmodel.train import seed_everything


def _model(length=40):
    seed_everything(0)
    return build_model(ModelConfig(seq_length=length))


def test_ism_matrix_shape_and_ref_zeros():
    model = _model(40)
    seq = "ACGT" * 10
    matrix = ism_matrix(model, seq)
    assert matrix.shape == (40, 4)
    # The reference base's own entry must be exactly 0 at every position.
    for i, base in enumerate(seq):
        assert matrix[i, BASE_TO_IDX[base]] == 0.0


def test_importance_track_shape_and_modes():
    model = _model(24)
    matrix = ism_matrix(model, "ACGTACGTACGTACGTACGTACGT")
    for mode in ("mean_abs", "min", "mean"):
        track = importance_track(matrix, mode=mode)
        assert track.shape == (24,)
    assert (importance_track(matrix, "mean_abs") >= 0).all()
    with pytest.raises(ValueError):
        importance_track(matrix, "bogus")


def test_variant_effect_matches_manual_delta():
    model = _model(40)
    seq = "ACGT" * 10
    pos, alt = 5, "T"
    ref_base = seq[pos]
    assert ref_base != alt
    manual = predict_activity(model, seq[:pos] + alt + seq[pos + 1 :]) - predict_activity(model, seq)
    assert variant_effect(model, seq, pos, alt) == pytest.approx(manual, abs=1e-5)
    # It also equals the corresponding ISM matrix entry.
    matrix = ism_matrix(model, seq)
    assert variant_effect(model, seq, pos, alt) == pytest.approx(matrix[pos, BASE_TO_IDX[alt]], abs=1e-5)


def test_variant_effect_noop_is_zero_and_validates():
    model = _model(40)
    seq = "ACGT" * 10
    assert variant_effect(model, seq, 0, seq[0]) == 0.0  # same base -> no change
    with pytest.raises(ValueError):
        variant_effect(model, seq, 0, "Z")
    with pytest.raises(ValueError):
        variant_effect(model, seq, 999, "A")


def test_top_positions_orders_by_importance():
    track = np.array([0.1, 0.9, 0.2, 0.5])
    assert top_positions(track, 2) == [1, 3]


def test_ism_runs_on_sequence_shorter_than_training_length():
    # Regression: `ism --seq ACGT` crashed because the fallback trained at len(seq).
    # The CNN is length-agnostic (padded convs + global pooling), so a model built for
    # L=200 must still ISM a 4 bp query without error.
    seed_everything(0)
    model = build_model(ModelConfig(seq_length=200))
    matrix = ism_matrix(model, "ACGT")
    assert matrix.shape == (4, 4)
