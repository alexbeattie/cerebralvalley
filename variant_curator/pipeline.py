"""Evidence-assembly pipeline.

Given a gene + HGVS change, this:
  1. validates the gene against the MVP allowlist,
  2. resolves consequence + coordinates via VEP,
  3. routes by consequence:
       - missense / nonsense -> gnomAD + ClinVar (the protein-evidence path),
       - non-coding (intronic / splice / UTR / regulatory) -> AlphaGenome, plus
         gnomAD + ClinVar (frequency and priors matter just as much here),
       - anything else -> refused as out of scope.

The output is an EvidenceBundle with a linked source on every item.
"""

from __future__ import annotations

from .clients.alphagenome import fetch_alphagenome
from .clients.clinvar import fetch_clinvar
from .clients.gnomad import fetch_gnomad
from .clients.vep import fetch_vep
from .genes import get_gene
from .http import get_client
from .models import EvidenceBundle, VariantConsequence, VariantInput

IN_SCOPE_CONSEQUENCES = {VariantConsequence.MISSENSE, VariantConsequence.NONSENSE}

# Non-coding consequences routed to AlphaGenome instead of being rejected.
NONCODING_CONSEQUENCES = {
    VariantConsequence.SPLICE_DONOR,
    VariantConsequence.SPLICE_ACCEPTOR,
    VariantConsequence.SPLICE_REGION,
    VariantConsequence.INTRON,
    VariantConsequence.FIVE_PRIME_UTR,
    VariantConsequence.THREE_PRIME_UTR,
    VariantConsequence.UPSTREAM,
    VariantConsequence.DOWNSTREAM,
    VariantConsequence.REGULATORY,
    VariantConsequence.INTERGENIC,
}


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

        if vep.consequence in NONCODING_CONSEQUENCES:
            # --- AlphaGenome: splicing + regulatory interpretation --------
            # Uses SDK/env, not the httpx client. Degrades gracefully offline.
            ag = fetch_alphagenome(vep.chrom, vep.pos, vep.ref, vep.alt, gene.symbol)
            bundle.alphagenome = ag
            if not ag.found:
                bundle.warnings.append(f"AlphaGenome skipped: {ag.reason}")
        elif vep.consequence not in IN_SCOPE_CONSEQUENCES:
            bundle.warnings.append(
                f"Out of scope: consequence is '{vep.consequence_raw}'. "
                "This tool handles missense/nonsense and non-coding (intronic, "
                "splice, UTR, regulatory) variants."
            )
            return bundle

        # gnomAD frequency and ClinVar priors are relevant on both paths.
        # --- gnomAD: ancestry-stratified frequency ------------------------
        variant_id = vep.gnomad_variant_id
        if variant_id:
            bundle.gnomad = fetch_gnomad(client, variant_id)
        else:
            bundle.warnings.append("Could not build a gnomAD variant id from VEP coordinates.")

        # --- ClinVar: prior classifications -------------------------------
        bundle.clinvar = fetch_clinvar(client, gene.symbol, variant.hgvs_c)

    return bundle
