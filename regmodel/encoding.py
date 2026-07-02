"""One-hot DNA encoding and small sequence helpers.

The model reads DNA as a (4, L) one-hot tensor over the fixed alphabet ACGT, channel
order A,C,G,T. We keep this in one place because encoding, decoding, and single-base
mutation (for ISM) must all agree on that channel order -- a silent mismatch there would
quietly corrupt every prediction and every ISM delta.
"""

from __future__ import annotations

import numpy as np

# Channel order is fixed and shared by every consumer (model, ISM, plots).
BASES = "ACGT"
BASE_TO_IDX = {b: i for i, b in enumerate(BASES)}
IDX_TO_BASE = {i: b for i, b in enumerate(BASES)}
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def validate_sequence(seq: str) -> str:
    """Uppercase and check the alphabet. Raises on anything but ACGTN.

    We validate rather than silently coerce because a stray character usually means the
    caller passed the wrong window; failing loudly is cheaper than a wrong prediction.
    """
    s = seq.strip().upper()
    if not s:
        raise ValueError("empty sequence")
    bad = set(s) - set("ACGTN")
    if bad:
        raise ValueError(f"sequence has non-ACGTN characters: {sorted(bad)}")
    return s


def one_hot_encode(seq: str) -> np.ndarray:
    """Encode a DNA string as a float32 array of shape (4, L), channel order ACGT.

    'N' (and any zero-signal position) maps to an all-zero column rather than a guessed
    base, so ambiguity contributes nothing instead of a fabricated one-hot.
    """
    s = validate_sequence(seq)
    arr = np.zeros((4, len(s)), dtype=np.float32)
    for j, base in enumerate(s):
        idx = BASE_TO_IDX.get(base)
        if idx is not None:  # 'N' stays all-zero
            arr[idx, j] = 1.0
    return arr


def one_hot_decode(arr: np.ndarray) -> str:
    """Inverse of `one_hot_encode`; an all-zero column decodes back to 'N'."""
    if arr.ndim != 2 or arr.shape[0] != 4:
        raise ValueError(f"expected (4, L) array, got {arr.shape}")
    out = []
    for j in range(arr.shape[1]):
        col = arr[:, j]
        out.append("N" if not col.any() else IDX_TO_BASE[int(np.argmax(col))])
    return "".join(out)


def reverse_complement(seq: str) -> str:
    """Reverse-complement a DNA string (kept here so augmentation/ISM share one impl)."""
    s = validate_sequence(seq)
    return "".join(COMPLEMENT[b] for b in reversed(s))


def encode_batch(seqs: list[str]) -> np.ndarray:
    """Encode a list of equal-length sequences to (N, 4, L). Raises on ragged input."""
    if not seqs:
        raise ValueError("no sequences to encode")
    length = len(seqs[0])
    if any(len(s) != length for s in seqs):
        raise ValueError("encode_batch requires equal-length sequences")
    return np.stack([one_hot_encode(s) for s in seqs], axis=0)


def random_sequence(rng: np.random.Generator, length: int) -> str:
    """A uniform-random ACGT sequence; the synthetic background is built from these."""
    idx = rng.integers(0, 4, size=length)
    return "".join(BASES[i] for i in idx)
