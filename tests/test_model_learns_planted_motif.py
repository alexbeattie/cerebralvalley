"""The core scientific claim, checked offline: the model learns the planted signal AND its
ISM lights up the planted motif. Trains briefly on synthetic data on CPU with a fixed seed."""

from __future__ import annotations

import numpy as np

from regmodel.data import ACTIVATOR_MOTIF, make_labeled_sequence, make_synthetic_mpra
from regmodel.ism import importance_track, ism_matrix, span_contains_top, top_positions
from regmodel.model import ModelConfig
from regmodel.train import TrainConfig, train_model


def _trained():
    ds = make_synthetic_mpra(n=1800, length=120, seed=11)
    result = train_model(ds, TrainConfig.fast(seed=11), ModelConfig(seq_length=120))
    return result


def test_model_recovers_planted_signal_on_heldout():
    result = _trained()
    # A modest but real generalization bar: the planted motif is strong, so a working model
    # clears this comfortably; a broken pipeline (e.g. shuffled labels) would not.
    assert result.metrics["pearson"] > 0.5
    assert result.metrics["spearman"] > 0.5


def test_ism_peaks_at_the_planted_motif():
    result = _trained()
    seq, span = make_labeled_sequence(length=120, seed=777, activator_motif=ACTIVATOR_MOTIF)
    matrix = ism_matrix(result.model, seq)
    assert matrix.shape == (120, 4)

    track = importance_track(matrix)
    # The single most-sensitive base should sit inside the planted motif, and the motif's
    # mean importance should tower over the genomic background.
    assert span_contains_top(track, span, k=3)
    start, end = span
    inside = track[start:end].mean()
    outside_mask = np.ones(track.shape[0], dtype=bool)
    outside_mask[start:end] = False
    outside = track[outside_mask].mean()
    assert inside > 3 * outside
