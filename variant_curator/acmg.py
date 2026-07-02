"""Map AlphaGenome signals to ACMG-style supporting evidence.

Provenance-first and deliberately conservative: an AlphaGenome prediction is
research-model output, not clinically validated, so it can contribute *supporting*
strength only — `PP3` (computational evidence toward deleterious) or `BP4` (toward
benign) — and never a standalone Pathogenic/Benign assertion. Every emitted
criterion carries the AlphaGenome `Source` and the research-use caveat, so a curator
can see exactly what drove it and weigh it accordingly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from .models import AlphaGenomeEvidence, Source

# Provisional calibrated-quantile thresholds. These are placeholders and MUST be
# calibrated against a labelled set (e.g. ClinVar splicing/regulatory variants with
# known effect) before any real use. quantile_score is calibrated genome-wide, so a
# value near 1.0 is an extreme effect.
PP3_QUANTILE = 0.90  # at/above this |quantile| we call the predicted effect strong
BP4_QUANTILE = 0.50  # every signal below this -> predicted benign / no effect

CAVEAT = (
    "AlphaGenome is a research model, not clinically validated. Treat as supporting "
    "computational evidence only; it cannot on its own drive a Pathogenic/Benign call."
)


@dataclass
class AcmgCriterion:
    """One ACMG/AMP criterion hit, with the source that justifies it."""

    code: str  # "PP3" / "BP4"
    strength: str  # always "supporting" here
    direction: str  # "pathogenic" / "benign"
    rationale: str
    caveat: str
    source: Optional[Source] = None


def _max_quantile(magnitudes: list[float | None]) -> float | None:
    present = [m for m in magnitudes if m is not None]
    return max(present) if present else None


def map_alphagenome_to_acmg(evidence: AlphaGenomeEvidence) -> list[AcmgCriterion]:
    """Return the supporting ACMG criteria implied by an AlphaGenome result.

    Abstains (empty list) when AlphaGenome was not run, when no calibrated quantile
    is available, or when the strongest effect sits between the benign and
    pathogenic thresholds — we don't guess in the grey zone.
    """
    if not evidence.found:
        return []

    splice_mag = _max_quantile([s.max_quantile for s in evidence.splicing])
    reg_mag = _max_quantile([abs(r.quantile_score) if r.quantile_score is not None else None for r in evidence.regulatory])

    hits: list[AcmgCriterion] = []

    strong_splice = splice_mag is not None and splice_mag >= PP3_QUANTILE
    strong_reg = reg_mag is not None and reg_mag >= PP3_QUANTILE

    if strong_splice or strong_reg:
        reasons = []
        if strong_splice:
            reasons.append(f"strong predicted splicing disruption (|quantile|={splice_mag:.2f})")
        if strong_reg:
            reasons.append(f"strong predicted regulatory effect (|quantile|={reg_mag:.2f})")
        hits.append(
            AcmgCriterion(
                code="PP3",
                strength="supporting",
                direction="pathogenic",
                rationale="; ".join(reasons) + f" — thresholds provisional (|quantile|>={PP3_QUANTILE}).",
                caveat=CAVEAT,
                source=evidence.source,
            )
        )
        return hits

    # No strong effect: apply BP4 only when we actually evaluated calibrated scores
    # and every one of them is below the benign threshold.
    evaluated = [m for m in (splice_mag, reg_mag) if m is not None]
    if evaluated and all(m < BP4_QUANTILE for m in evaluated):
        hits.append(
            AcmgCriterion(
                code="BP4",
                strength="supporting",
                direction="benign",
                rationale=(
                    f"no predicted splicing or regulatory effect (max |quantile|={max(evaluated):.2f} "
                    f"< {BP4_QUANTILE}) — thresholds provisional."
                ),
                caveat=CAVEAT,
                source=evidence.source,
            )
        )
    return hits
