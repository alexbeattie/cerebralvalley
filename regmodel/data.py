"""MPRA dataset abstraction: a built-in synthetic task plus a network-gated real loader.

Why synthetic-by-default: a hackathon demo must run on a plane. The synthetic MPRA below
plants a known activator motif (and a weaker repressor) into random background, so activity
is a *known function of sequence*. That gives us ground truth twice over: a well-trained
model should (a) predict held-out activity well, and (b) have its ISM light up exactly the
planted motif positions -- both of which the tests assert. Real MPRA loaders are optional
and gated behind an explicit opt-in because they need downloads/network.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .encoding import BASES, random_sequence

# The planted "activator" is an AP-1-like site (TGACTCA); the "repressor" is an arbitrary
# distinct 6-mer. Real TF motifs are degenerate PWMs -- we use exact k-mers on purpose so
# the ground-truth ISM peak is unambiguous and a test can assert on it.
ACTIVATOR_MOTIF = "TGACTCA"
REPRESSOR_MOTIF = "GCGATCGC"


@dataclass
class MPRADataset:
    """A batch of (sequence, activity) pairs plus provenance.

    `activator_spans` is the ground-truth annotation the synthetic generator emits (list of
    (start, end_exclusive) per sequence, empty when none planted). Real loaders leave it
    None -- there is no ground-truth motif to compare ISM against outside the synthetic task.
    """

    sequences: list[str]
    activities: np.ndarray  # shape (N,), float32
    seq_length: int
    provenance: dict
    activator_spans: list[list[tuple[int, int]]] | None = None
    repressor_spans: list[list[tuple[int, int]]] | None = None

    def __post_init__(self) -> None:
        if len(self.sequences) != len(self.activities):
            raise ValueError("sequences and activities length mismatch")


def _plant(chars: list[str], motif: str, start: int) -> None:
    """Overwrite `chars[start:start+len(motif)]` with `motif` in place."""
    for k, base in enumerate(motif):
        chars[start + k] = base


def make_synthetic_mpra(
    n: int = 4000,
    length: int = 200,
    seed: int = 0,
    *,
    activator_motif: str = ACTIVATOR_MOTIF,
    repressor_motif: str = REPRESSOR_MOTIF,
    activator_weight: float = 3.0,
    repressor_weight: float = 1.5,
    p_activator: float = 0.5,
    p_repressor: float = 0.3,
    noise_sd: float = 0.5,
    intercept: float = 0.0,
) -> MPRADataset:
    """Generate a learnable, interpretable synthetic MPRA.

    Each sequence is uniform-random background of length `length`. With probability
    `p_activator` we plant the activator motif at a random position (raising activity by
    `activator_weight` per copy); independently with `p_repressor` we plant the repressor
    (lowering it by `repressor_weight`). Activity is that signal plus Gaussian noise, so the
    label genuinely depends on motif presence *and position content* -- ISM can recover it.

    Seeded end to end: same `seed` -> identical sequences, activities, and planted spans.
    """
    rng = np.random.default_rng(seed)
    a_len, r_len = len(activator_motif), len(repressor_motif)

    sequences: list[str] = []
    activities = np.empty(n, dtype=np.float32)
    activator_spans: list[list[tuple[int, int]]] = []
    repressor_spans: list[list[tuple[int, int]]] = []

    for i in range(n):
        chars = list(random_sequence(rng, length))
        a_here: list[tuple[int, int]] = []
        r_here: list[tuple[int, int]] = []

        if rng.random() < p_activator:
            start = int(rng.integers(0, length - a_len + 1))
            _plant(chars, activator_motif, start)
            a_here.append((start, start + a_len))

        if rng.random() < p_repressor:
            start = int(rng.integers(0, length - r_len + 1))
            # Avoid clobbering a just-planted activator so the label stays interpretable.
            if not _overlaps(start, start + r_len, a_here):
                _plant(chars, repressor_motif, start)
                r_here.append((start, start + r_len))

        signal = (
            intercept
            + activator_weight * len(a_here)
            - repressor_weight * len(r_here)
            + rng.normal(0.0, noise_sd)
        )
        sequences.append("".join(chars))
        activities[i] = signal
        activator_spans.append(a_here)
        repressor_spans.append(r_here)

    provenance = {
        "source": "SyntheticMPRA (regmodel.data.make_synthetic_mpra)",
        "kind": "synthetic",
        "description": (
            "Random ACGT background with a planted AP-1-like activator motif and a weaker "
            "repressor; activity = weighted motif counts + Gaussian noise."
        ),
        "n": n,
        "seq_length": length,
        "seed": seed,
        "alphabet": BASES,
        "activator_motif": activator_motif,
        "repressor_motif": repressor_motif,
        "activator_weight": activator_weight,
        "repressor_weight": repressor_weight,
        "p_activator": p_activator,
        "p_repressor": p_repressor,
        "noise_sd": noise_sd,
        "research_use_only": True,
    }

    return MPRADataset(
        sequences=sequences,
        activities=activities,
        seq_length=length,
        provenance=provenance,
        activator_spans=activator_spans,
        repressor_spans=repressor_spans,
    )


def make_labeled_sequence(
    length: int = 200,
    seed: int = 12345,
    *,
    activator_motif: str = ACTIVATOR_MOTIF,
    position: int | None = None,
) -> tuple[str, tuple[int, int]]:
    """Build one background sequence with a single activator planted at a known span.

    Used by the demo and tests to probe ISM: we know exactly where the motif is, so we can
    check the model's ISM importance peaks there. Returns (sequence, (start, end_exclusive)).
    """
    rng = np.random.default_rng(seed)
    a_len = len(activator_motif)
    chars = list(random_sequence(rng, length))
    start = position if position is not None else int(rng.integers(0, length - a_len + 1))
    _plant(chars, activator_motif, start)
    return "".join(chars), (start, start + a_len)


def _overlaps(start: int, end: int, spans: list[tuple[int, int]]) -> bool:
    return any(start < s_end and s_start < end for s_start, s_end in spans)


def load_real_mpra(name: str, cache_dir: str | None = None) -> MPRADataset:  # pragma: no cover - network-gated stub
    """Network-gated loader for a real published MPRA. Not exercised offline.

    Candidate public datasets to wire up here (all sequence -> activity, MPRA-style):

      - **lentiMPRA, developing human brain** (Pollard / Ahituv labs) -- the directly
        on-theme dataset; enhancer activity of human sequences incl. human accelerated
        regions.
      - **Sharpr-MPRA** (Ernst et al., 2016) -- large tiled promoter/enhancer MPRA in
        K562/HepG2, a common sequence-model benchmark.

    Implementing this means: download the supplementary table, parse (sequence, activity)
    columns, and return an `MPRADataset` whose provenance records the accession, URL, and
    retrieval timestamp -- matching the offline path's provenance contract. Deliberately a
    stub: the whole demo/test suite must run with no network, so we abstain rather than
    ship a half-wired downloader.
    """
    raise NotImplementedError(
        f"Real MPRA loader for {name!r} is not wired up (offline-first build). "
        "Use make_synthetic_mpra() for the offline demo, or implement this against one of "
        "the datasets named in the docstring (lentiMPRA developing brain, Sharpr-MPRA)."
    )
