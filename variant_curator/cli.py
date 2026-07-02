"""CLI: assemble and print the evidence bundle for one variant.

Day 1-2 deliverable: one variant flowing end to end through VEP, gnomAD, and
ClinVar, proving the data path before anything smart sits on top of it.

Usage:
    python -m variant_curator.cli --gene BRCA1 --hgvs "c.181T>G"
    python -m variant_curator.cli            # runs the hardcoded demo variant
"""

from __future__ import annotations

import argparse

from .acmg import map_alphagenome_to_acmg
from .models import EvidenceBundle, VariantInput
from .pipeline import assemble_evidence

# Hardcoded demo variant: BRCA1 c.181T>G (p.Cys61Gly), a well-known pathogenic
# missense in the RING domain. Good end-to-end smoke test.
DEMO_VARIANT = VariantInput(gene="BRCA1", hgvs_c="c.181T>G")

# Illustrative non-coding demos for the AlphaGenome path. HGVS is indicative;
# without a key these print a clean "skipped" line rather than crashing.
DEMO_NONCODING = [
    # CFTR c.3717+12191C>T (legacy 3849+10kbC>T): classic deep intronic variant
    # that creates a cryptic splice site — invisible to protein predictors.
    VariantInput(gene="CFTR", hgvs_c="c.3717+12191C>T"),
    # LDLR 5' regulatory-region example (familial hypercholesterolemia, liver).
    VariantInput(gene="LDLR", hgvs_c="c.-135C>G"),
]


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

    a = bundle.alphagenome
    if a is not None:
        print("\n[AlphaGenome] non-coding interpretation (research use only)")
        if a.found:
            print(f"  variant     : {a.variant_id}  ({a.api})")
            print(f"  top effect  : {a.top_effect}")
            if a.splicing:
                print("  splicing:")
                for s in a.splicing:
                    q = f"|q|={s.max_quantile:.2f}" if s.max_quantile is not None else "|q|=n/a"
                    print(f"    - {s.scorer:18s} {q}  {s.interpretation}")
            if a.regulatory:
                print("  regulatory:")
                for r in a.regulatory:
                    q = f"q={r.quantile_score:+.2f}" if r.quantile_score is not None else "q=n/a"
                    curie = f" [{r.ontology_curie}]" if r.ontology_curie else ""
                    print(f"    - {r.modality:12s} {q}  {r.top_tissue}{curie}  ({r.direction or 'change'})")
            print(f"  source      : {a.source.url}")
            print("  CAVEAT      : research model, not clinically validated — supporting evidence only.")

            criteria = map_alphagenome_to_acmg(a)
            print("\n[ACMG] supporting criteria from AlphaGenome")
            if criteria:
                for crit in criteria:
                    print(f"  {crit.code} ({crit.strength}, {crit.direction}): {crit.rationale}")
                    if crit.source:
                        print(f"      source: {crit.source.url}")
                    print(f"      caveat: {crit.caveat}")
            else:
                print("  none — AlphaGenome effect sits in the grey zone; abstaining.")
        else:
            print(f"  skipped: {a.reason}")
            if a.source:
                print(f"  source      : {a.source.url}")

    if bundle.warnings:
        print("\n[warnings]")
        for w in bundle.warnings:
            print(f"  - {w}")
    print()


def main() -> None:
    parser = argparse.ArgumentParser(description="Assemble ACMG evidence for one variant.")
    parser.add_argument("--gene", help="Gene symbol (must be in the allowlist)")
    parser.add_argument("--hgvs", help="HGVS change, e.g. c.181T>G or a deep intronic c.")
    parser.add_argument(
        "--noncoding-demo",
        action="store_true",
        help="Run the illustrative non-coding (AlphaGenome) demo variants",
    )
    args = parser.parse_args()

    if args.gene and args.hgvs:
        variants = [VariantInput(gene=args.gene, hgvs_c=args.hgvs)]
    elif args.noncoding_demo:
        print("(running illustrative non-coding demo variants)\n")
        variants = DEMO_NONCODING
    else:
        print("(no --gene/--hgvs given; running hardcoded demo variant)\n")
        variants = [DEMO_VARIANT]

    for variant in variants:
        bundle = assemble_evidence(variant)
        _print_bundle(bundle)


if __name__ == "__main__":
    main()
