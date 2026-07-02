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
    """The molecular consequences the tool handles.

    Missense/nonsense go down the protein-evidence path (gnomAD + ClinVar). The
    non-coding classes below are routed to AlphaGenome, which can see splicing and
    regulatory effects that consequence-based tools and protein predictors miss.
    """

    MISSENSE = "missense_variant"
    NONSENSE = "stop_gained"
    # Non-coding classes routed to AlphaGenome.
    SPLICE_DONOR = "splice_donor_variant"
    SPLICE_ACCEPTOR = "splice_acceptor_variant"
    SPLICE_REGION = "splice_region_variant"
    INTRON = "intron_variant"
    FIVE_PRIME_UTR = "5_prime_UTR_variant"
    THREE_PRIME_UTR = "3_prime_UTR_variant"
    UPSTREAM = "upstream_gene_variant"
    DOWNSTREAM = "downstream_gene_variant"
    REGULATORY = "regulatory_region_variant"
    INTERGENIC = "intergenic_variant"
    OTHER = "other"  # the true fallback; in-scope check will reject these


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
class SplicingSignal:
    """One AlphaGenome splicing scorer reduced to its strongest effect on a gene.

    Splicing scores are non-directional (always positive): magnitude is the whole
    signal, so we keep the max calibrated |quantile| and raw score. `direction` is
    populated only where a scorer distinguishes gain vs loss.
    """

    scorer: str  # e.g. "SPLICE_SITES"
    output_type: str  # AlphaGenome output modality the scorer reads
    top_gene: str
    max_quantile: Optional[float]  # |quantile_score|, calibrated genome-wide
    max_raw: Optional[float]  # |raw_score|, fallback magnitude
    direction: Optional[str]  # usually None; splicing scores are non-directional
    interpretation: str
    source: Optional[Source] = None

    @property
    def magnitude(self) -> Optional[float]:
        return self.max_quantile if self.max_quantile is not None else self.max_raw


@dataclass
class RegulatorySignal:
    """One AlphaGenome regulatory modality reduced to its top tissue/biosample.

    Regulatory scores are signed: for expression (`RNA_SEQ`) the sign is up- vs
    down-regulation; for accessibility it is gain vs loss. We keep the signed
    quantile (primary) and raw score.
    """

    modality: str  # e.g. "RNA_SEQ", "DNASE"
    output_type: str
    top_tissue: str  # gtex_tissue or biosample_name
    ontology_curie: Optional[str]
    quantile_score: Optional[float]  # signed, calibrated genome-wide
    raw_score: Optional[float]  # signed, fallback magnitude
    direction: Optional[str]  # "up" / "down"
    interpretation: str
    source: Optional[Source] = None

    @property
    def magnitude(self) -> Optional[float]:
        primary = self.quantile_score if self.quantile_score is not None else self.raw_score
        return abs(primary) if primary is not None else None


@dataclass
class AlphaGenomeEvidence:
    """AlphaGenome interpretation of a non-coding variant, or explicit absence.

    `found=False` with a `reason` means AlphaGenome could not be run (no API key,
    package not installed, or an API error) — we abstain rather than fabricate a
    score. `research_use_only` stays true: these are research-model predictions,
    not clinically validated, and map to supporting ACMG evidence at most.
    """

    found: bool
    variant_id: Optional[str] = None
    api: str = ""  # model/API version string
    splicing: list[SplicingSignal] = field(default_factory=list)
    regulatory: list[RegulatorySignal] = field(default_factory=list)
    reason: str = ""  # why found is False, surfaced to the curator
    source: Optional[Source] = None
    research_use_only: bool = True

    @property
    def top_effect(self) -> str:
        """The single strongest predicted effect across splicing + regulatory."""
        if not self.found:
            return "not evaluated"
        best = None
        best_mag = -1.0
        for sig in (*self.splicing, *self.regulatory):
            mag = sig.magnitude
            if mag is not None and mag > best_mag:
                best_mag, best = mag, sig
        if best is None:
            return "no signal returned"
        return best.interpretation


@dataclass
class EvidenceBundle:
    """Everything the ACMG mapper gets to reason over for a single variant."""

    variant: VariantInput
    vep: Optional[VepEvidence] = None
    gnomad: Optional[GnomadEvidence] = None
    clinvar: Optional[ClinvarEvidence] = None
    alphagenome: Optional[AlphaGenomeEvidence] = None
    warnings: list[str] = field(default_factory=list)

    def sources(self) -> list[Source]:
        out: list[Source] = []
        for ev in (self.vep, self.gnomad, self.clinvar, self.alphagenome):
            if ev is not None and getattr(ev, "source", None) is not None:
                out.append(ev.source)
        return out
