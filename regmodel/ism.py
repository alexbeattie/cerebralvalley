"""In-silico mutagenesis: read a per-base mutation-effect map out of the trained model.

ISM is the interpretability workhorse the Pollard lab uses to ask "what does this single
base do?". We predict the reference activity, then for every position substitute each of the
three alternative bases, predict again, and record the delta. Where the model has learned a
motif, mutating inside it swings the prediction; elsewhere deltas are ~0. That signed map is
both the explanation and, on the synthetic task, a checkable ground truth.
"""

from __future__ import annotations

import numpy as np
import torch

from .encoding import BASE_TO_IDX, BASES, one_hot_encode, validate_sequence
from .model import ActivityCNN


def predict_activity(model: ActivityCNN, seq: str) -> float:
    """Predict scalar activity for a single sequence."""
    model.eval()
    x = torch.from_numpy(one_hot_encode(seq)).unsqueeze(0)
    with torch.no_grad():
        return float(model(x).item())


def _predict_batch(model: ActivityCNN, encoded: np.ndarray) -> np.ndarray:
    model.eval()
    x = torch.from_numpy(encoded)
    with torch.no_grad():
        return model(x).cpu().numpy()


def ism_matrix(model: ActivityCNN, seq: str) -> np.ndarray:
    """Return an (L, 4) array of predicted-activity deltas vs the reference.

    Entry [i, b] = activity(seq with position i set to base b) - activity(seq). The
    reference base's own entry is exactly 0. Channel order is ACGT (matches encoding). All
    L*3 mutant predictions run in one batched forward pass for speed.
    """
    s = validate_sequence(seq)
    length = len(s)
    ref_encoded = one_hot_encode(s)
    ref_pred = _predict_batch(model, ref_encoded[None, :, :])[0]

    # Build every single-base mutant; skip the ref base (its delta is 0 by definition).
    mutants = []
    coords = []  # (position, base_idx) each mutant corresponds to
    for i in range(length):
        ref_idx = BASE_TO_IDX.get(s[i])
        for b in range(4):
            if b == ref_idx:
                continue
            m = ref_encoded.copy()
            m[:, i] = 0.0
            m[b, i] = 1.0
            mutants.append(m)
            coords.append((i, b))

    matrix = np.zeros((length, 4), dtype=np.float32)
    if mutants:
        preds = _predict_batch(model, np.stack(mutants, axis=0))
        for (i, b), p in zip(coords, preds):
            matrix[i, b] = p - ref_pred
    return matrix


def importance_track(matrix: np.ndarray, mode: str = "mean_abs") -> np.ndarray:
    """Collapse the (L, 4) ISM matrix to a per-position scalar for plotting/ranking.

    - "mean_abs" (default): mean |delta| over the 3 alternative bases -- overall sensitivity,
      which is what peaks at a functional motif regardless of effect direction.
    - "min": most negative delta (largest activity *loss* on mutation).
    - "mean": signed mean delta.
    """
    if matrix.ndim != 2 or matrix.shape[1] != 4:
        raise ValueError(f"expected (L, 4) ISM matrix, got {matrix.shape}")
    if mode == "mean_abs":
        # Divide by 3 (the non-ref bases); the ref entry is 0 and doesn't bias the sum.
        return np.abs(matrix).sum(axis=1) / 3.0
    if mode == "min":
        return matrix.min(axis=1)
    if mode == "mean":
        return matrix.sum(axis=1) / 3.0
    raise ValueError(f"unknown importance mode: {mode!r}")


def top_positions(track: np.ndarray, k: int = 5) -> list[int]:
    """Indices of the `k` highest-importance positions, most important first."""
    k = min(k, track.shape[0])
    return [int(i) for i in np.argsort(track)[::-1][:k]]


def variant_effect(model: ActivityCNN, ref_seq: str, pos: int, alt: str) -> float:
    """Predicted activity delta for a single substitution: activity(alt) - activity(ref).

    `pos` is 0-based into `ref_seq`; `alt` is a single ACGT base. A positive value means the
    model predicts the change *increases* activity, negative means it decreases it.
    """
    s = validate_sequence(ref_seq)
    alt = alt.strip().upper()
    if len(alt) != 1 or alt not in BASES:
        raise ValueError(f"alt must be a single ACGT base, got {alt!r}")
    if not (0 <= pos < len(s)):
        raise ValueError(f"pos {pos} out of range for sequence length {len(s)}")
    if s[pos] == alt:
        return 0.0  # no-op substitution
    ref_pred = predict_activity(model, s)
    mutant = s[:pos] + alt + s[pos + 1 :]
    return predict_activity(model, mutant) - ref_pred


def span_contains_top(track: np.ndarray, span: tuple[int, int], k: int = 3) -> bool:
    """True if any of the top-`k` importance positions falls inside `span` (start, end excl).

    A tiny helper so the 'ISM recovers the planted motif' test reads as one assertion.
    """
    start, end = span
    return any(start <= p < end for p in top_positions(track, k))
