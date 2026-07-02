"""Evidence-assembly pipeline.

Given a gene + HGVS coding change, this:
  1. validates the gene against the MVP allowlist,
  2. resolves consequence + coordinates via VEP,
  3. enforces the missense/nonsense scope,
  4. fetches gnomAD ancestry frequencies and ClinVar priors.

The output is an EvidenceBundle with a linked source on every item. No ACMG
reasoning happens here yet (that is the Day 3-4 layer); this is the data path.
"""

from __future__ import annotations

from .clients.clinvar import fetch_clinvar
from .clients.gnomad import fetch_gnomad
from .clients.vep import fetch_vep
from .genes import get_gene
from .http import get_client
from .models import EvidenceBundle, VariantConsequence, VariantInput

IN_SCOPE_CONSEQUENCES = {VariantConsequence.MISSENSE, VariantConsequence.NONSENSE}


def assemble_evidence(variant: VariantInput) -> EvidenceBundle:
    gene = get_gene(variant.gene)  # raises for out-of-allowlist genes
    bundle = EvidenceBundle(variant=variant)

    with get_client() as client:
        # --- VEP: consequence + coordinates -------------------------------
        vep = fetch_vep(client, gene.refseq_transcript, variant.hgvs_c, gene.symbol)
        bundle.vep = vep

        if not vep.found:
            bundle.warnings.append(
                f"VEP could not resolve {gene.refseq_transcript}:{variant.hgvs_c}. "
                "Check the HGVS notation and transcript."
            )
            return bundle

        if vep.consequence not in IN_SCOPE_CONSEQUENCES:
            bundle.warnings.append(
                f"Out of MVP scope: consequence is '{vep.consequence_raw}'. "
                "This prototype only handles missense and nonsense SNVs."
            )
            return bundle

        # --- gnomAD: ancestry-stratified frequency ------------------------
        variant_id = vep.gnomad_variant_id
        if variant_id:
            bundle.gnomad = fetch_gnomad(client, variant_id)
        else:
            bundle.warnings.append("Could not build a gnomAD variant id from VEP coordinates.")

        # --- ClinVar: prior classifications -------------------------------
        bundle.clinvar = fetch_clinvar(client, gene.symbol, variant.hgvs_c)

    return bundle
