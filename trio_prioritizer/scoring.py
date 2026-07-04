"""Attach a non-coding pathogenicity score to a variant by reusing the AlphaGenome engine.

We do not re-implement scoring. `variant_curator.clients.alphagenome.fetch_alphagenome`
is the engine; here we (a) call it through an *injectable* `scorer` callable so tests and
the offline demo pass a fake and never touch the network, (b) reduce its splicing +
regulatory signals to one calibrated magnitude, and (c) degrade honestly to
`available=False` with a reason when it can't run — carrying the `Source` either way.

AlphaGenome is coordinate-based, which is what trio inputs are, so it is the primary
scorer. `regmodel` (sequence->activity + ISM) is a strictly optional secondary signal:
it needs a local sequence window, so it only runs when one is supplied, and it is imported
lazily so the core demo never requires torch.
"""

from __future__ import annotations

from typing import Callable, Optional

from variant_curator.clients.alphagenome import fetch_alphagenome
from variant_curator.models import AlphaGenomeEvidence, Source

from .models import NoncodingScore, TrioVariant

# The scorer signature we depend on: the same one fetch_alphagenome exposes.
Scorer = Callable[..., AlphaGenomeEvidence]

SCORER_NAME = "AlphaGenome"


def _evidence_magnitude(ev: AlphaGenomeEvidence) -> Optional[float]:
    """Reduce all splicing + regulatory signals to the single strongest magnitude.

    Matches how `AlphaGenomeEvidence.top_effect` picks its headline signal, so the number
    and the phrase always describe the same effect.
    """
    best = None
    for sig in (*ev.splicing, *ev.regulatory):
        mag = sig.magnitude
        if mag is not None and (best is None or mag > best):
            best = mag
    return best


def score_variant(
    v: TrioVariant,
    scorer: Scorer = fetch_alphagenome,
    *,
    window: Optional[str] = None,
    regmodel_predictor=None,
) -> NoncodingScore:
    """Produce a `NoncodingScore` for one variant.

    Coding variants skip AlphaGenome entirely (their evidence path lives elsewhere in the
    project); we still return an explicit, provenance-free abstention so the candidate can
    be ranked on inheritance. Non-coding variants are scored via the injected `scorer`.

    When `window` (a local reference sequence) and `regmodel_predictor` are both supplied,
    an ISM-based delta is attached as a secondary signal — never required for the core demo.
    """
    if v.is_coding:
        return NoncodingScore(
            available=False,
            scorer="n/a",
            reason="coding variant — non-coding scorer not applicable (coding path handled elsewhere)",
        )

    ev = scorer(v.chrom, v.pos, v.ref, v.alt, v.gene)

    if not ev.found:
        score = NoncodingScore(
            available=False,
            scorer=SCORER_NAME,
            reason=ev.reason or "scorer returned no result",
            source=ev.source,
        )
    else:
        score = NoncodingScore(
            available=True,
            scorer=SCORER_NAME,
            magnitude=_evidence_magnitude(ev),
            top_effect=ev.top_effect,
            source=ev.source,
        )

    _maybe_attach_regmodel(score, v, window, regmodel_predictor)
    return score


def _maybe_attach_regmodel(
    score: NoncodingScore,
    v: TrioVariant,
    window: Optional[str],
    regmodel_predictor,
) -> None:
    """Optionally add a regmodel ISM delta as a secondary signal. No-op without a window.

    `window` is a (ref_seq, pos0, alt) tuple locating the variant inside a local sequence.
    `regmodel_predictor` is a loaded `regmodel.model.ActivityCNN` (or any callable-compatible
    predictor). Imports are lazy so torch is never required for the AlphaGenome-only path.
    """
    if window is None or regmodel_predictor is None:
        return
    try:
        ref_seq, pos0, alt = window
        from regmodel.ism import variant_effect  # lazy: only when a window is supplied

        delta = variant_effect(regmodel_predictor, ref_seq, pos0, alt)
    except Exception:
        return  # optional signal: never let it break the primary score
    score.secondary_scorer = "regmodel ISM"
    score.secondary_magnitude = abs(float(delta))
    score.secondary_source = Source(
        name="regmodel (MPRA sequence->activity ISM)",
        url="local:regmodel",
        retrieved_at="",
        detail=f"|Δactivity|={abs(float(delta)):.3f} at {v.variant_id}",
    )


def score_variants(
    variants: list[TrioVariant], scorer: Scorer = fetch_alphagenome
) -> dict[str, NoncodingScore]:
    """Score a list of variants; returns a map keyed by `TrioVariant.variant_id`."""
    return {v.variant_id: score_variant(v, scorer) for v in variants}
