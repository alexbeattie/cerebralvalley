"""A seeded synthetic trio scenario with a known 'answer', plus an offline fake scorer.

The demo has to *prove itself* with no network and no API key, so this module plants a
ground truth the pipeline must recover:

  - a **de novo** deep-intronic variant in SCN1A (a plausible severe-epilepsy gene), with a
    strong predicted splice effect,
  - a **compound-het pair** in PAH (classic recessive, phenylketonuria): a coding missense
    inherited from the father + a deep-intronic variant inherited from the mother — the
    textbook "the second hit is non-coding" case,
  - several benign distractors: het alleles inherited from an unaffected parent with weak
    predicted effect, plus a lone het (one hit, no partner) that must rank low.

The offline scorer returns deterministic magnitudes keyed to these planted variants so the
demo runs with no key. Its `Source` is explicitly labelled ILLUSTRATIVE so a synthetic
number is never mistaken for a real AlphaGenome prediction.

No PHI: coordinates are illustrative (roughly the genes' GRCh38 loci) and carry no sample.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from variant_curator.http import now_iso
from variant_curator.models import (
    AlphaGenomeEvidence,
    RegulatorySignal,
    Source,
    SplicingSignal,
    VariantConsequence,
)

from .models import Genotype, TrioVariant

OFFLINE_SCORER_NAME = "AlphaGenome (offline synthetic scorer — ILLUSTRATIVE)"


@dataclass(frozen=True)
class _Planted:
    """A synthetic variant plus the illustrative effect the offline scorer should report."""

    variant: TrioVariant
    kind: str  # "splicing" or "regulatory"
    magnitude: float  # 0..1 illustrative calibrated magnitude
    phrase: str  # human effect description
    note: str  # why it's in the scenario (the 'answer key')


def _v(chrom, pos, ref, alt, gene, cons, proband, mother, father) -> TrioVariant:
    return TrioVariant(
        chrom=chrom, pos=pos, ref=ref, alt=alt, gene=gene, consequence=cons,
        proband=Genotype.parse(proband), mother=Genotype.parse(mother), father=Genotype.parse(father),
    )


# --- The planted answer + distractors ---------------------------------------------------
# Kept as module data so both the trio builder and the offline scorer read the same table.
_SCENARIO: list[_Planted] = [
    # ANSWER 1 — de novo deep-intronic in SCN1A: present in child, 0/0 in both parents.
    _Planted(
        variant=_v("2", 165_991_200, "G", "A", "SCN1A", VariantConsequence.INTRON,
                   "0/1", "0/0", "0/0"),
        kind="splicing", magnitude=0.93,
        phrase="predicted cryptic splice-site activation in SCN1A",
        note="causal de novo (should rank #1)",
    ),
    # ANSWER 2a — comp-het PAH, PATERNAL coding missense (coding path; not AG-scored).
    _Planted(
        variant=_v("12", 102_852_100, "C", "T", "PAH", VariantConsequence.MISSENSE,
                   "0/1", "0/0", "0/1"),
        kind="none", magnitude=0.0,
        phrase="coding missense (scored on the protein path)",
        note="causal comp-het allele 1 (paternal, coding)",
    ),
    # ANSWER 2b — comp-het PAH, MATERNAL deep-intronic 'second hit' (the non-coding one).
    _Planted(
        variant=_v("12", 102_912_400, "A", "G", "PAH", VariantConsequence.INTRON,
                   "0/1", "0/1", "0/0"),
        kind="splicing", magnitude=0.88,
        phrase="predicted branch-point / deep-intronic splice defect in PAH",
        note="causal comp-het allele 2 (maternal, non-coding second hit)",
    ),
    # DISTRACTOR — inherited het from mother, weak effect (benign background).
    _Planted(
        variant=_v("2", 179_400_000, "T", "C", "TTN", VariantConsequence.INTRON,
                   "0/1", "0/1", "0/0"),
        kind="splicing", magnitude=0.08,
        phrase="no meaningful splice effect predicted in TTN",
        note="benign inherited distractor",
    ),
    # DISTRACTOR — inherited het from father, weak regulatory effect.
    _Planted(
        variant=_v("14", 23_430_000, "G", "C", "MYH7", VariantConsequence.FIVE_PRIME_UTR,
                   "0/1", "0/0", "0/1"),
        kind="regulatory", magnitude=0.12,
        phrase="negligible expression change predicted for MYH7",
        note="benign inherited distractor",
    ),
    # DISTRACTOR — het carried by BOTH parents (biparental, inherited), weak effect.
    _Planted(
        variant=_v("2", 21_000_500, "A", "T", "APOB", VariantConsequence.INTRON,
                   "0/1", "0/1", "0/1"),
        kind="splicing", magnitude=0.15,
        phrase="no significant splice effect predicted in APOB",
        note="benign inherited distractor",
    ),
    # DISTRACTOR — lone het in a recessive gene: ONE hit, no partner. Must rank low even
    # though its non-coding effect is moderate — this is the 'second hit still missing' case.
    _Planted(
        variant=_v("17", 80_110_000, "C", "G", "GAA", VariantConsequence.INTRON,
                   "0/1", "0/1", "0/0"),
        kind="splicing", magnitude=0.55,
        phrase="moderate predicted splice effect in GAA (but no second hit found)",
        note="lone het — deliberately weak candidate without a partner allele",
    ),
]

# Index the illustrative effects by variant coordinates for the offline scorer.
_BY_ID: dict[str, _Planted] = {p.variant.variant_id: p for p in _SCENARIO}


def build_synthetic_trio() -> list[TrioVariant]:
    """Return the planted trio candidate set (deterministic; order as authored)."""
    return [p.variant for p in _SCENARIO]


def answer_key() -> dict[str, str]:
    """Map variant_id -> why it's in the scenario, for the demo to self-check against."""
    return {p.variant.variant_id: p.note for p in _SCENARIO}


def causal_variant_ids() -> set[str]:
    """The variant ids the tool must surface at the top (de novo + the comp-het pair)."""
    return {vid for vid, note in answer_key().items() if note.startswith("causal")}


def offline_fake_scorer(
    chrom: str, pos: int, ref: str, alt: str, gene_symbol: str, *, api_key: Optional[str] = None
) -> AlphaGenomeEvidence:
    """A deterministic, offline stand-in for `fetch_alphagenome`.

    Same signature as the real scorer, so it is injected identically. Returns the planted,
    ILLUSTRATIVE magnitude for known synthetic variants and abstains for anything else —
    never a fabricated real score. The `Source` name marks it synthetic so the caveat is
    honest end to end.
    """
    variant_id = f"{chrom}:{pos}:{ref}>{alt}"
    source = Source(
        name=OFFLINE_SCORER_NAME,
        url="offline:synthetic-trio",
        retrieved_at=now_iso(),
        detail=f"{variant_id} (illustrative; not a real AlphaGenome call)",
    )

    planted = _BY_ID.get(variant_id)
    if planted is None or planted.kind == "none":
        return AlphaGenomeEvidence(
            found=False,
            variant_id=variant_id,
            reason="no offline synthetic score for this variant",
            source=source,
        )

    interp = f"{planted.phrase} (|magnitude|={planted.magnitude:.2f}, illustrative)"
    splicing: list[SplicingSignal] = []
    regulatory: list[RegulatorySignal] = []
    if planted.kind == "splicing":
        splicing.append(
            SplicingSignal(
                scorer="SPLICE_SITES", output_type="SPLICE_SITES", top_gene=gene_symbol,
                max_quantile=planted.magnitude, max_raw=None, direction=None,
                interpretation=interp, source=source,
            )
        )
    else:  # regulatory
        regulatory.append(
            RegulatorySignal(
                modality="RNA_SEQ", output_type="RNA_SEQ", top_tissue="illustrative tissue",
                ontology_curie=None, quantile_score=planted.magnitude, raw_score=None,
                direction="up", interpretation=interp, source=source,
            )
        )

    return AlphaGenomeEvidence(
        found=True,
        variant_id=variant_id,
        api="offline synthetic scorer",
        splicing=splicing,
        regulatory=regulatory,
        source=source,
    )
