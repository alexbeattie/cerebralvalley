"""Cross-check: where does our small MPRA model agree or disagree with AlphaGenome?

This is the differentiated result. For each variant we have two independent read-outs of
"how regulatory is this change": our model's ISM delta on a local sequence window, and
AlphaGenome's calibrated splicing/regulatory magnitude on the 1 Mb genomic context. We line
them up into a tidy table, sign the agreement, and (with >=2 comparable rows) correlate them.

Reuses `variant_curator.clients.alphagenome.fetch_alphagenome` verbatim -- same graceful
degradation: no key / no SDK / API error -> the external half is skipped with a recorded
reason and Source, and we still emit our-model-only rows. We never fabricate an AG score.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from variant_curator.clients.alphagenome import fetch_alphagenome
from variant_curator.models import AlphaGenomeEvidence, Source

from .ism import variant_effect
from .model import ActivityCNN


@dataclass
class VariantSpec:
    """One variant to cross-check: genomic coords for AlphaGenome + a local window for us.

    `local_window_seq` is the DNA our model scores; `window_pos` is the 0-based index of the
    variant base within that window (defaults to the center). In the offline demo the windows
    are synthetic/illustrative -- flagged as such -- because we can't fetch the genome without
    network.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    gene_symbol: str
    local_window_seq: str
    window_pos: int | None = None
    note: str = ""

    def center_pos(self) -> int:
        return self.window_pos if self.window_pos is not None else len(self.local_window_seq) // 2


@dataclass
class ComparisonRow:
    variant_id: str
    gene_symbol: str
    our_delta: float  # signed ISM delta from our model
    our_magnitude: float  # |our_delta|, for correlation against AG's magnitude
    alphagenome_found: bool
    alphagenome_magnitude: float | None  # top |quantile|/|raw| across AG signals
    alphagenome_signed: float | None  # signed regulatory quantile if AG gave direction
    alphagenome_top_effect: str
    agreement_sign: str  # "agree" / "disagree" / "n/a" (needs both signs)
    source: Source | None  # AG Source (present even when found=False)
    note: str = ""


@dataclass
class ComparisonResult:
    rows: list[ComparisonRow]
    table: pd.DataFrame
    alphagenome_available: bool
    correlation: float | None  # Pearson(our_magnitude, ag_magnitude) over comparable rows
    note: str
    sources: list[Source] = field(default_factory=list)


def _ag_magnitude(ev: AlphaGenomeEvidence) -> tuple[float | None, float | None, str]:
    """Reduce AG evidence to (magnitude, signed_regulatory, top_effect_text).

    magnitude is the strongest |signal| across splicing+regulatory (splicing is
    non-directional); signed is the strongest *regulatory* signed quantile when present, so
    we can compare direction where AG offers one.
    """
    if not ev.found:
        return None, None, "not evaluated"
    best_mag = None
    for sig in (*ev.splicing, *ev.regulatory):
        mag = sig.magnitude
        if mag is not None and (best_mag is None or mag > best_mag):
            best_mag = mag
    signed = None
    best_reg = None
    for sig in ev.regulatory:
        q = sig.quantile_score if sig.quantile_score is not None else sig.raw_score
        if q is not None and (best_reg is None or abs(q) > best_reg):
            best_reg, signed = abs(q), q
    return best_mag, signed, ev.top_effect


def _agreement(our_delta: float, ag_signed: float | None) -> str:
    if ag_signed is None or our_delta == 0 or ag_signed == 0:
        return "n/a"
    return "agree" if (our_delta > 0) == (ag_signed > 0) else "disagree"


def compare_variants(
    model: ActivityCNN,
    variants: list[VariantSpec],
    *,
    api_key: str | None = None,
) -> ComparisonResult:
    """Run the head-to-head over `variants`. Degrades cleanly with no AlphaGenome access."""
    rows: list[ComparisonRow] = []
    sources: list[Source] = []
    any_ag = False

    for v in variants:
        our_delta = variant_effect(model, v.local_window_seq, v.center_pos(), v.alt)
        ev = fetch_alphagenome(v.chrom, v.pos, v.ref, v.alt, v.gene_symbol, api_key=api_key)
        if ev.source is not None:
            sources.append(ev.source)
        mag, signed, top = _ag_magnitude(ev)
        any_ag = any_ag or ev.found

        rows.append(
            ComparisonRow(
                variant_id=ev.variant_id or f"chr{v.chrom}:{v.pos}:{v.ref}>{v.alt}",
                gene_symbol=v.gene_symbol,
                our_delta=round(our_delta, 5),
                our_magnitude=round(abs(our_delta), 5),
                alphagenome_found=ev.found,
                alphagenome_magnitude=None if mag is None else round(mag, 5),
                alphagenome_signed=None if signed is None else round(signed, 5),
                alphagenome_top_effect=top,
                agreement_sign=_agreement(our_delta, signed),
                source=ev.source,
                note=v.note if ev.found else (v.note + f" | AG skipped: {ev.reason}").strip(" |"),
            )
        )

    table = pd.DataFrame(
        [
            {
                "variant_id": r.variant_id,
                "gene": r.gene_symbol,
                "our_delta": r.our_delta,
                "our_magnitude": r.our_magnitude,
                "ag_found": r.alphagenome_found,
                "ag_magnitude": r.alphagenome_magnitude,
                "ag_signed": r.alphagenome_signed,
                "ag_top_effect": r.alphagenome_top_effect,
                "agreement": r.agreement_sign,
                "note": r.note,
            }
            for r in rows
        ]
    )

    correlation = _correlation(rows) if any_ag else None
    if any_ag:
        note = "AlphaGenome cross-check ran; correlation computed over comparable rows."
    else:
        reason = rows[0].note if rows else "no variants"
        note = (
            "AlphaGenome comparison SKIPPED (no key / SDK / API) -- our-model deltas only. "
            f"Reason on first row: {reason}. Local windows are illustrative in offline mode."
        )

    return ComparisonResult(
        rows=rows,
        table=table,
        alphagenome_available=any_ag,
        correlation=correlation,
        note=note,
        sources=sources,
    )


def _correlation(rows: list[ComparisonRow]) -> float | None:
    pairs = [
        (r.our_magnitude, r.alphagenome_magnitude)
        for r in rows
        if r.alphagenome_found and r.alphagenome_magnitude is not None
    ]
    if len(pairs) < 2:
        return None
    ours = np.array([p[0] for p in pairs])
    ag = np.array([p[1] for p in pairs])
    if ours.std() == 0 or ag.std() == 0:
        return None
    from scipy.stats import pearsonr

    return float(pearsonr(ours, ag)[0])
