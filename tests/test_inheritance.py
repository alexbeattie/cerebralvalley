"""Inheritance classification + compound-het pairing, all offline."""

from __future__ import annotations

from trio_prioritizer.inheritance import (
    classify_single,
    find_compound_het_pairs,
    transmitting_parent,
)
from trio_prioritizer.models import Genotype, InheritanceMode, TrioVariant
from variant_curator.models import VariantConsequence as VC


def _v(gene, pos, cons, p, m, f, chrom="1", ref="A", alt="G") -> TrioVariant:
    return TrioVariant(
        chrom=chrom, pos=pos, ref=ref, alt=alt, gene=gene, consequence=cons,
        proband=Genotype.parse(p), mother=Genotype.parse(m), father=Genotype.parse(f),
    )


def test_de_novo_detected():
    v = _v("SCN1A", 100, VC.INTRON, "0/1", "0/0", "0/0")
    call = classify_single(v)
    assert call.mode is InheritanceMode.DE_NOVO
    assert call.parent_of_origin is None
    assert "de novo calling needs QC" in call.confidence  # QC caveat is surfaced


def test_de_novo_hom_alt_in_proband_still_de_novo():
    v = _v("SCN1A", 100, VC.INTRON, "1/1", "0/0", "0/0")
    assert classify_single(v).mode is InheritanceMode.DE_NOVO


def test_homozygous_recessive_detected():
    v = _v("PAH", 200, VC.INTRON, "1/1", "0/1", "0/1")
    call = classify_single(v)
    assert call.mode is InheritanceMode.HOMOZYGOUS_RECESSIVE
    assert call.parent_of_origin is None


def test_inherited_dominant_deprioritized_but_listed():
    v = _v("MYH7", 300, VC.FIVE_PRIME_UTR, "0/1", "0/1", "0/0")
    call = classify_single(v)
    assert call.mode is InheritanceMode.INHERITED_DOMINANT
    assert call.parent_of_origin == "maternal"


def test_not_in_proband_is_unknown():
    v = _v("MYH7", 300, VC.INTRON, "0/0", "0/1", "0/0")
    assert classify_single(v).mode is InheritanceMode.UNKNOWN


def test_missing_genotype_needs_qc():
    v = _v("MYH7", 300, VC.INTRON, "0/1", "./.", "0/0")
    call = classify_single(v)
    assert call.mode is InheritanceMode.UNKNOWN
    assert "needs-QC" in call.confidence


def test_transmitting_parent_phasing():
    assert transmitting_parent(_v("X", 1, VC.INTRON, "0/1", "0/1", "0/0")) == "maternal"
    assert transmitting_parent(_v("X", 1, VC.INTRON, "0/1", "0/0", "0/1")) == "paternal"
    assert transmitting_parent(_v("X", 1, VC.INTRON, "0/1", "0/1", "0/1")) is None  # biparental


def test_genotype_parse_vcf_and_words():
    assert Genotype.parse("0/1") is Genotype.HET
    assert Genotype.parse("1|1") is Genotype.HOM_ALT
    assert Genotype.parse("0/0") is Genotype.HOM_REF
    assert Genotype.parse("./.") is Genotype.MISSING
    assert Genotype.parse("het") is Genotype.HET
    assert Genotype.parse("hom_alt") is Genotype.HOM_ALT


def test_compound_het_requires_one_from_each_parent():
    # Coding allele from father + non-coding allele from mother -> a comp-het pair.
    coding = _v("PAH", 200, VC.MISSENSE, "0/1", "0/0", "0/1")
    intronic = _v("PAH", 900, VC.INTRON, "0/1", "0/1", "0/0")
    pairs = find_compound_het_pairs([coding, intronic])
    assert len(pairs) == 1
    paternal, maternal, call = pairs[0]
    assert call.mode is InheritanceMode.COMPOUND_HET
    assert paternal is coding and maternal is intronic  # emitted (paternal, maternal)
    assert paternal.is_coding and maternal.is_noncoding  # coding + non-coding pair works


def test_compound_het_not_formed_from_same_parent():
    # Two hets both inherited from the mother is NOT a comp-het (cis, not trans).
    a = _v("PAH", 200, VC.INTRON, "0/1", "0/1", "0/0")
    b = _v("PAH", 900, VC.INTRON, "0/1", "0/1", "0/0")
    assert find_compound_het_pairs([a, b]) == []


def test_compound_het_needs_two_distinct_sites():
    # The "same" variant appearing once cannot pair with itself.
    only = _v("PAH", 200, VC.INTRON, "0/1", "0/1", "0/0")
    assert find_compound_het_pairs([only]) == []


def test_compound_het_scoped_within_gene():
    # A maternal het in gene A and a paternal het in gene B do not pair across genes.
    a = _v("GENEA", 200, VC.INTRON, "0/1", "0/1", "0/0")
    b = _v("GENEB", 900, VC.INTRON, "0/1", "0/0", "0/1")
    assert find_compound_het_pairs([a, b]) == []
