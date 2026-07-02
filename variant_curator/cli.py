"""CLI: assemble and print the evidence bundle for one variant.

Day 1-2 deliverable: one variant flowing end to end through VEP, gnomAD, and
ClinVar, proving the data path before anything smart sits on top of it.

Usage:
    python -m variant_curator.cli --gene BRCA1 --hgvs "c.181T>G"
    python -m variant_curator.cli            # runs the hardcoded demo variant
"""

from __future__ import annotations

import argparse

from .models import EvidenceBundle, VariantInput
from .pipeline import assemble_evidence

# Hardcoded demo variant: BRCA1 c.181T>G (p.Cys61Gly), a well-known pathogenic
# missense in the RING domain. Good end-to-end smoke test.
DEMO_VARIANT = VariantInput(gene="BRCA1", hgvs_c="c.181T>G")


def _print_bundle(bundle: EvidenceBundle) -> None:
    v = bundle.variant
    print("=" * 72)
    print(f"  Evidence bundle for {v.label}")
    print("=" * 72)

    vep = bundle.vep
    print("\n[VEP] consequence")
    if vep and vep.found:
        print(f"  consequence : {vep.consequence_raw} ({vep.consequence.value})")
        print(f"  protein     : {vep.protein_hgvs}")
        print(f"  GRCh38       : {vep.chrom}:{vep.pos} {vep.ref}>{vep.alt}")
        print(f"  gnomAD id   : {vep.gnomad_variant_id}")
        print(f"  source      : {vep.source.url}")
    else:
        print("  not resolved")

    g = bundle.gnomad
    print("\n[gnomAD] ancestry-stratified allele frequency")
    if g and g.found:
        af = f"{g.global_af:.2e}" if g.global_af else "0"
        print(f"  global AF   : {af}  (AC={g.global_ac}, AN={g.global_an})")
        print(f"  FAF95 popmax: {g.faf95_popmax}  in {g.faf95_popmax_ancestry}")
        for p in sorted(g.populations, key=lambda x: -(x.af or 0)):
            af_p = f"{p.af:.2e}" if p.af else "0"
            print(f"    {p.ancestry_id:10s} {p.ancestry_label:32s} AF={af_p}  (AC={p.ac}, AN={p.an})")
        print(f"  source      : {g.source.url}")
    elif g:
        print(f"  NOT FOUND in gnomAD {g.dataset} (absence is itself evidence: supports PM2)")
        print(f"  source      : {g.source.url}")
    else:
        print("  not queried")

    c = bundle.clinvar
    print("\n[ClinVar] prior classifications")
    if c and c.found:
        print(f"  aggregate   : {c.aggregate_significance}  (conflict={c.has_conflict})")
        for s in c.submissions:
            stars = "*" * s.star_rating + "." * (4 - s.star_rating)
            print(f"    [{stars}] {s.clinical_significance:32s} {s.accession}  ({s.condition})")
        print(f"  source      : {c.source.url}")
    elif c:
        print("  no matching ClinVar record")
        print(f"  source      : {c.source.url}")
    else:
        print("  not queried")

    if bundle.warnings:
        print("\n[warnings]")
        for w in bundle.warnings:
            print(f"  - {w}")
    print()


def main() -> None:
    parser = argparse.ArgumentParser(description="Assemble ACMG evidence for one variant.")
    parser.add_argument("--gene", help="Gene symbol (must be in the allowlist)")
    parser.add_argument("--hgvs", help="HGVS coding change, e.g. c.181T>G")
    args = parser.parse_args()

    if args.gene and args.hgvs:
        variant = VariantInput(gene=args.gene, hgvs_c=args.hgvs)
    else:
        print("(no --gene/--hgvs given; running hardcoded demo variant)\n")
        variant = DEMO_VARIANT

    bundle = assemble_evidence(variant)
    _print_bundle(bundle)


if __name__ == "__main__":
    main()
