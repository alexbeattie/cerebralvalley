"""Trio genotypes -> inheritance mode, plus compound-het pairing by parent-of-origin.

This is the new core the brief asks for: the per-variant non-coding *scoring* already
exists (`variant_curator.clients.alphagenome`); here we add the inheritance reasoning that
decides *which* variants matter and *how* they are transmitted.

Scope (MVP): autosomal de novo and recessive (homozygous + compound het). X-linked and
imprinting-aware inheritance are documented future extensions — deliberately not
implemented, because guessing chromosome/sex handling here would be worse than abstaining.
"""

from __future__ import annotations

from typing import Optional

from .models import Genotype, InheritanceCall, InheritanceMode, TrioVariant

# Qualitative confidence flags (the brief asks for a simple genotype-consistent vs
# needs-QC signal, not a probability). De novo carries an explicit QC caveat because a
# de novo *call* is only as good as coverage/quality at that site.
CONF_CONSISTENT = "genotype-consistent"
CONF_DE_NOVO = "genotype-consistent (de novo calling needs QC/coverage confirmation)"
CONF_NEEDS_QC = "needs-QC (genotype missing or ambiguous)"


def transmitting_parent(v: TrioVariant) -> Optional[str]:
    """Which parent transmitted the alt allele to the proband, from genotypes alone.

    Returns "maternal" (mother carries, father hom-ref), "paternal" (mirror), or None
    when both or neither parent carries it (biparental / can't phase from genotypes).
    """
    mom = v.mother.carries_alt
    dad = v.father.carries_alt
    if mom and not dad:
        return "maternal"
    if dad and not mom:
        return "paternal"
    return None


def _has_missing(v: TrioVariant) -> bool:
    return Genotype.MISSING in (v.proband, v.mother, v.father)


def classify_single(v: TrioVariant) -> InheritanceCall:
    """Classify a single variant's inheritance mode from the trio genotypes.

    Compound-het is *not* decided here — it is a property of a pair within a gene, handled
    by `find_compound_het_pairs`. A lone het transmitted from one parent falls through to
    INHERITED_DOMINANT (and is deprioritized downstream), which is exactly the "one hit,
    second hit still missing" state the tool should rank low.
    """
    p, mom, dad = v.proband, v.mother, v.father

    # A variant absent from the proband is not a candidate for the proband's disease.
    if not p.carries_alt:
        return InheritanceCall(
            mode=InheritanceMode.UNKNOWN,
            confidence=CONF_NEEDS_QC if _has_missing(v) else CONF_CONSISTENT,
            rationale="not present in proband; not a candidate",
        )

    if _has_missing(v):
        return InheritanceCall(
            mode=InheritanceMode.UNKNOWN,
            confidence=CONF_NEEDS_QC,
            rationale="a trio genotype is missing; inheritance cannot be assigned",
        )

    # De novo: present in proband, absent in both parents.
    if mom == Genotype.HOM_REF and dad == Genotype.HOM_REF:
        return InheritanceCall(
            mode=InheritanceMode.DE_NOVO,
            parent_of_origin=None,
            confidence=CONF_DE_NOVO,
            rationale="present in proband, absent in both parents (0/0) — de novo",
        )

    # Homozygous recessive: proband hom-alt, both parents carriers.
    if p == Genotype.HOM_ALT and mom == Genotype.HET and dad == Genotype.HET:
        return InheritanceCall(
            mode=InheritanceMode.HOMOZYGOUS_RECESSIVE,
            parent_of_origin=None,
            confidence=CONF_CONSISTENT,
            rationale="proband homozygous alt; both parents heterozygous carriers",
        )

    # Otherwise it is inherited from at least one parent -> dominant-model candidate.
    origin = transmitting_parent(v)
    who = origin or "both parents"
    return InheritanceCall(
        mode=InheritanceMode.INHERITED_DOMINANT,
        parent_of_origin=origin,
        confidence=CONF_CONSISTENT,
        rationale=(
            f"inherited from {who}; for a severe-disease proband a single inherited allele "
            "is a low prior (lone het = second hit still missing)"
        ),
    )


def _is_phased_het(v: TrioVariant) -> Optional[str]:
    """A proband het transmitted cleanly from exactly one parent returns that parent."""
    if v.proband != Genotype.HET:
        return None
    return transmitting_parent(v)  # "maternal"/"paternal"/None


def find_compound_het_pairs(
    variants: list[TrioVariant],
) -> list[tuple[TrioVariant, TrioVariant, InheritanceCall]]:
    """Find compound-het pairs: two *different* het variants in the same gene, one from the
    mother and one from the father.

    This is the "find the second hit" case. It must work when one allele is coding and the
    other is a deep-intronic/regulatory variant — nothing here inspects consequence, only
    parent-of-origin, so a coding+non-coding pair is found exactly like any other. The
    emitted call phases the pair (paternal allele, maternal allele).
    """
    by_gene: dict[str, list[TrioVariant]] = {}
    for v in variants:
        by_gene.setdefault(v.gene, []).append(v)

    pairs: list[tuple[TrioVariant, TrioVariant, InheritanceCall]] = []
    for gene, group in by_gene.items():
        # Only proband-hets phased to a single parent can form a comp-het pair.
        maternal = [v for v in group if _is_phased_het(v) == "maternal"]
        paternal = [v for v in group if _is_phased_het(v) == "paternal"]
        for mv, pv in _best_pairs(maternal, paternal):
            call = InheritanceCall(
                mode=InheritanceMode.COMPOUND_HET,
                parent_of_origin=None,  # biparental by definition; per-allele phase below
                confidence=CONF_CONSISTENT,
                rationale=(
                    f"compound het in {gene}: paternal {pv.variant_id} "
                    f"[{'coding' if pv.is_coding else 'non-coding'}] + maternal "
                    f"{mv.variant_id} [{'coding' if mv.is_coding else 'non-coding'}] "
                    "— trans-configured second hit"
                ),
            )
            # Emit as (paternal, maternal) so the pair reads allele-1/allele-2 consistently.
            pairs.append((pv, mv, call))
    return pairs


def _best_pairs(
    maternal: list[TrioVariant], paternal: list[TrioVariant]
) -> list[tuple[TrioVariant, TrioVariant]]:
    """Pair maternal x paternal het alleles, requiring two *distinct* sites.

    For the MVP we emit at most one pair per gene (the biologically relevant "one from each
    parent" pattern); with multiple candidates on a side we take the first distinct combo so
    the demo/tests are deterministic. Real phasing would rank by score — noted for later.
    """
    for mv, pv in ((m, p) for m in maternal for p in paternal):
        if (mv.chrom, mv.pos, mv.ref, mv.alt) != (pv.chrom, pv.pos, pv.ref, pv.alt):
            return [(mv, pv)]
    return []


# Kept for symmetry / potential future use: enumerate *all* distinct comp-het combos in a
# gene (not just the first). Unused by the MVP ranking but handy for exploration/tests.
def all_compound_het_combos(
    maternal: list[TrioVariant], paternal: list[TrioVariant]
) -> list[tuple[TrioVariant, TrioVariant]]:
    combos: list[tuple[TrioVariant, TrioVariant]] = []
    for mv, pv in ((m, p) for m in maternal for p in paternal):
        if (mv.chrom, mv.pos, mv.ref, mv.alt) != (pv.chrom, pv.pos, pv.ref, pv.alt):
            combos.append((mv, pv))
    # Distinct as unordered pairs.
    seen: set[frozenset] = set()
    out = []
    for mv, pv in combos:
        key = frozenset({mv.variant_id, pv.variant_id})
        if key not in seen:
            seen.add(key)
            out.append((mv, pv))
    return out
