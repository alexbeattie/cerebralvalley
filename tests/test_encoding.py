"""One-hot encoding round-trips, channel order, and reverse-complement -- all offline."""

from __future__ import annotations

import numpy as np
import pytest

from regmodel.encoding import (
    BASES,
    encode_batch,
    one_hot_decode,
    one_hot_encode,
    reverse_complement,
    validate_sequence,
)


def test_encode_shape_and_channel_order():
    arr = one_hot_encode("ACGT")
    assert arr.shape == (4, 4)
    # channel order is ACGT: the diagonal is hot.
    assert np.array_equal(arr, np.eye(4, dtype=np.float32))


def test_encode_decode_round_trips():
    seq = "ACGTACGTTTGGCCAA"
    assert one_hot_decode(one_hot_encode(seq)) == seq


def test_n_maps_to_zero_column_and_decodes_back():
    arr = one_hot_encode("ANG")
    assert arr[:, 1].sum() == 0.0  # 'N' is an all-zero column, not a guessed base
    assert one_hot_decode(arr) == "ANG"


def test_lowercase_is_accepted():
    assert one_hot_decode(one_hot_encode("acgt")) == "ACGT"


def test_reverse_complement():
    assert reverse_complement("AACCGGTT") == "AACCGGTT"[::-1].translate(str.maketrans("ACGT", "TGCA"))
    assert reverse_complement("ATGC") == "GCAT"


def test_validate_rejects_bad_alphabet():
    with pytest.raises(ValueError):
        validate_sequence("ACGTX")
    with pytest.raises(ValueError):
        validate_sequence("")


def test_encode_batch_shape_and_ragged_guard():
    batch = encode_batch(["ACGT", "TTTT"])
    assert batch.shape == (2, 4, 4)
    with pytest.raises(ValueError):
        encode_batch(["ACGT", "TTT"])


def test_bases_constant():
    assert BASES == "ACGT"
