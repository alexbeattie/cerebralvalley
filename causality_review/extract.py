"""Clinical-note -> HPO extraction (AI at the edge, ontology in control).

Pipeline:
  1. An LLM reads the free-text note and PROPOSES short phenotype phrases.
  2. Each phrase is resolved against the HPO ontology by the deterministic
     term search -- the model never emits an HPO id, so it cannot invent one.

The result keeps provenance: which note phrase produced which HPO term, and
which phrases could not be grounded. The scoring engine only ever sees real,
validated HPO terms.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import httpx

from .clients.hpo import resolve_term
from .clients.llm import extract_phenotype_phrases, extract_variants_raw
from .models import Phenotype


@dataclass
class ExtractedPhenotype:
    phrase: str          # what the LLM proposed from the note
    phenotype: Phenotype  # the HPO term it grounded to


@dataclass
class ExtractionResult:
    phenotypes: list[ExtractedPhenotype] = field(default_factory=list)
    ungrounded: list[str] = field(default_factory=list)  # phrases with no HPO match

    def deduped_phenotypes(self) -> list[Phenotype]:
        seen: set[str] = set()
        out: list[Phenotype] = []
        for e in self.phenotypes:
            if e.phenotype.hpo_id not in seen:
                seen.add(e.phenotype.hpo_id)
                out.append(e.phenotype)
        return out


def extract_from_notes(client: httpx.Client, notes: str, *, model: str | None = None) -> ExtractionResult:
    """Free-text clinical notes -> validated HPO phenotypes (LLM proposes, HPO validates)."""

    phrases = extract_phenotype_phrases(client, notes, model=model)
    result = ExtractionResult()
    for phrase in phrases:
        term = resolve_term(client, phrase)
        if term is None:
            result.ungrounded.append(phrase)
        else:
            result.phenotypes.append(ExtractedPhenotype(phrase=phrase, phenotype=term))
    return result


@dataclass
class ReportIngest:
    """Everything pulled from a dropped lab-report PDF."""

    variants: list[dict] = field(default_factory=list)  # {gene, hgvs, classification}
    phenotypes: ExtractionResult = field(default_factory=ExtractionResult)


def ingest_report(client: httpx.Client, report_text: str, *, model: str | None = None) -> ReportIngest:
    """Lab-report text -> reported variants + validated HPO phenotypes.

    The LLM only proposes text (variants it read off the page, phenotype phrases);
    HPO grounding stays deterministic so the scoring input is always real terms.
    """

    variants = extract_variants_raw(client, report_text, model=model)
    phenotypes = extract_from_notes(client, report_text, model=model)
    return ReportIngest(variants=variants, phenotypes=phenotypes)
