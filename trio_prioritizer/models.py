"""Core data models for trio prioritization.

Mirrors `variant_curator.models`: dataclasses, string-valued enums, and a `Source` on
every scored claim so each pathogenicity number links back to where it came from. We
reuse `variant_curator.models.Source` and `VariantConsequence` directly rather than
redefining them — the whole point is one shared provenance/consequence vocabulary.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from variant_curator.models import Source, VariantConsequence

# Consequences that go down the *coding* path (handled elsewhere in the project). Every
# other VariantConsequence is non-coding and is what this tool scores via AlphaGenome.
CODING_CONSEQUENCES = frozenset({VariantConsequence.MISSENSE, VariantConsequence.NONSENSE})


class Genotype(str, Enum):
    """A single-sample genotype at one biallelic site."""

    HOM_REF = "hom_ref"  # 0/0
    HET = "het"  # 0/1
    HOM_ALT = "hom_alt"  # 1/1
    MISSING = "missing"  # ./.

    @property
    def carries_alt(self) -> bool:
        """True when at least one alternate allele is present (HET or HOM_ALT)."""
        return self in (Genotype.HET, Genotype.HOM_ALT)

    @classmethod
    def parse(cls, text: str) -> "Genotype":
        """Parse a VCF GT ('0/1', '1|1', './.') or a word ('het', 'hom_alt')."""
        t = text.strip().lower()
        aliases = {
            "hom_ref": cls.HOM_REF, "homref": cls.HOM_REF, "ref": cls.HOM_REF,
            "het": cls.HET, "heterozygous": cls.HET,
            "hom_alt": cls.HOM_ALT, "homalt": cls.HOM_ALT, "hom": cls.HOM_ALT, "alt": cls.HOM_ALT,
            "missing": cls.MISSING, "none": cls.MISSING, "nocall": cls.MISSING,
        }
        if t in aliases:
            return aliases[t]
        # VCF-style: split on / or | and count alt alleles; any '.' means missing.
        alleles = t.replace("|", "/").split("/")
        if any(a == "." for a in alleles) or not alleles:
            return cls.MISSING
        try:
            alt_count = sum(1 for a in alleles if int(a) > 0)
        except ValueError:
            return cls.MISSING
        if alt_count == 0:
            return cls.HOM_REF
        if alt_count >= len(alleles):
            return cls.HOM_ALT
        return cls.HET


class InheritanceMode(str, Enum):
    """How a candidate is transmitted in the trio (autosomal MVP scope)."""

    DE_NOVO = "de_novo"
    HOMOZYGOUS_RECESSIVE = "homozygous_recessive"
    COMPOUND_HET = "compound_het"
    INHERITED_DOMINANT = "inherited_dominant"
    UNKNOWN = "unknown"


@dataclass
class TrioVariant:
    """One variant with the trio's genotypes.

    `consequence` reuses `variant_curator.models.VariantConsequence`, so `is_coding`
    keys off the same vocabulary the rest of the repo routes on.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    consequence: VariantConsequence
    proband: Genotype
    mother: Genotype
    father: Genotype

    @property
    def variant_id(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"

    @property
    def label(self) -> str:
        return f"{self.gene} {self.variant_id} [{self.consequence.value}]"

    @property
    def is_coding(self) -> bool:
        return self.consequence in CODING_CONSEQUENCES

    @property
    def is_noncoding(self) -> bool:
        return not self.is_coding


@dataclass
class InheritanceCall:
    """The inheritance interpretation of a candidate (single variant or comp-het pair)."""

    mode: InheritanceMode
    parent_of_origin: Optional[str] = None  # "maternal" / "paternal" / None (biparental/unknown)
    confidence: str = ""  # qualitative flag: "genotype-consistent" vs "needs-QC..."
    rationale: str = ""


@dataclass
class NoncodingScore:
    """A non-coding pathogenicity score attached to one variant, or explicit absence.

    `available=False` with a `reason` means the scorer could not run (no API key, coding
    variant, package missing, API error). We abstain rather than invent a number. The
    `source` carries the scorer's provenance either way.
    """

    available: bool
    scorer: str  # e.g. "AlphaGenome"
    magnitude: Optional[float] = None  # calibrated 0..~1 where available
    top_effect: str = ""
    reason: str = ""  # why unavailable
    source: Optional[Source] = None
    # Optional secondary signal from regmodel ISM (only when a local window is supplied).
    secondary_scorer: Optional[str] = None
    secondary_magnitude: Optional[float] = None
    secondary_source: Optional[Source] = None


@dataclass
class PrioritizedCandidate:
    """One ranked candidate: a single variant or a compound-het *pair*.

    `variants` holds one entry for de novo / homozygous / inherited candidates and two
    (paternal + maternal allele) for a compound-het pair. `scores` is aligned positionally
    with `variants`.
    """

    variants: list[TrioVariant]
    inheritance: InheritanceCall
    scores: list[NoncodingScore]
    combined_score: float
    rationale: str = ""

    @property
    def gene(self) -> str:
        return self.variants[0].gene

    @property
    def is_pair(self) -> bool:
        return len(self.variants) == 2

    @property
    def top_magnitude(self) -> Optional[float]:
        """Strongest available non-coding magnitude across the candidate's variants."""
        mags = [s.magnitude for s in self.scores if s.available and s.magnitude is not None]
        return max(mags) if mags else None

    @property
    def label(self) -> str:
        return " + ".join(v.variant_id for v in self.variants)

    def sources(self) -> list[Source]:
        """Every provenance carried by this candidate's scores (deduplicated by url+detail)."""
        out: list[Source] = []
        seen: set[tuple[str, str]] = set()
        for s in self.scores:
            for src in (s.source, s.secondary_source):
                if src is None:
                    continue
                key = (src.url, src.detail)
                if key not in seen:
                    seen.add(key)
                    out.append(src)
        return out
