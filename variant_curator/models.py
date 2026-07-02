"""Core data models.

Design principle: every piece of evidence carries the source it came from, so the
UI and the ACMG mapper can link each claim back to a URL. Nothing is asserted
without provenance. This is what a curator needs to trust an automated draft.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class VariantConsequence(str, Enum):
    """The subset of molecular consequences the MVP handles."""

    MISSENSE = "missense_variant"
    NONSENSE = "stop_gained"
    OTHER = "other"  # in-scope check will reject these


@dataclass
class VariantInput:
    """What the curator types in: a gene symbol plus an HGVS coding change."""

    gene: str
    hgvs_c: str  # e.g. "c.181T>G"

    @property
    def label(self) -> str:
        return f"{self.gene} {self.hgvs_c}"


@dataclass
class Source:
    """Provenance for a single evidence item. `url` should be directly openable."""

    name: str  # e.g. "gnomAD v4", "ClinVar", "Ensembl VEP"
    url: str
    retrieved_at: str  # ISO timestamp
    detail: str = ""


@dataclass
class PopulationFrequency:
    """Allele frequency in one gnomAD genetic ancestry group."""

    ancestry_id: str  # gnomAD population id, e.g. "nfe", "afr", "amr", "eas", "sas"
    ancestry_label: str
    ac: int  # allele count
    an: int  # allele number (total observed)

    @property
    def af(self) -> Optional[float]:
        return (self.ac / self.an) if self.an else None


@dataclass
class GnomadEvidence:
    """Ancestry-stratified allele frequency from gnomAD, or explicit absence."""

    found: bool
    variant_id: Optional[str]
    dataset: str
    global_ac: int = 0
    global_an: int = 0
    global_af: Optional[float] = None
    faf95_popmax: Optional[float] = None  # filtering AF (95% CI) across pops
    faf95_popmax_ancestry: Optional[str] = None
    populations: list[PopulationFrequency] = field(default_factory=list)
    source: Optional[Source] = None

    def ancestry(self, ancestry_id: str) -> Optional[PopulationFrequency]:
        for p in self.populations:
            if p.ancestry_id == ancestry_id:
                return p
        return None


@dataclass
class ClinvarSubmission:
    """One ClinVar record for the variant."""

    accession: str  # VCV / RCV
    clinical_significance: str
    review_status: str
    star_rating: int
    condition: str
    last_evaluated: str


@dataclass
class ClinvarEvidence:
    found: bool
    submissions: list[ClinvarSubmission] = field(default_factory=list)
    has_conflict: bool = False
    aggregate_significance: str = ""
    source: Optional[Source] = None


@dataclass
class VepEvidence:
    """Consequence and coordinates resolved from HGVS by Ensembl VEP."""

    found: bool
    consequence: VariantConsequence
    consequence_raw: str  # the exact VEP term
    gene_symbol: str
    protein_hgvs: Optional[str]  # e.g. "p.Cys61Gly"
    # GRCh38 minimal representation used to build a gnomAD variant id.
    chrom: Optional[str] = None
    pos: Optional[int] = None
    ref: Optional[str] = None
    alt: Optional[str] = None
    source: Optional[Source] = None

    @property
    def gnomad_variant_id(self) -> Optional[str]:
        if None in (self.chrom, self.pos, self.ref, self.alt):
            return None
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alt}"


@dataclass
class EvidenceBundle:
    """Everything the ACMG mapper gets to reason over for a single variant."""

    variant: VariantInput
    vep: Optional[VepEvidence] = None
    gnomad: Optional[GnomadEvidence] = None
    clinvar: Optional[ClinvarEvidence] = None
    warnings: list[str] = field(default_factory=list)

    def sources(self) -> list[Source]:
        out: list[Source] = []
        for ev in (self.vep, self.gnomad, self.clinvar):
            if ev is not None and getattr(ev, "source", None) is not None:
                out.append(ev.source)
        return out
