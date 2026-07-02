"""Ensembl VEP client.

Resolves an HGVS coding change (against a RefSeq transcript) into:
  - the molecular consequence (we only proceed for missense / nonsense),
  - the protein-level HGVS,
  - GRCh38 minimal coordinates, from which we build a gnomAD variant id.

Public REST endpoint, no key required. GRCh38 assembly.
"""

from __future__ import annotations

import httpx

from ..http import get_json, now_iso
from ..models import Source, VariantConsequence, VepEvidence

VEP_BASE = "https://rest.ensembl.org"

# VEP returns Sequence Ontology terms. Missense/nonsense go down the protein path;
# the non-coding terms are mapped to the classes the pipeline routes to AlphaGenome.
# Several finer-grained splice terms collapse onto SPLICE_REGION on purpose — the
# AlphaGenome splicing scorers care about the effect, not the exact SO label.
_CONSEQUENCE_MAP = {
    "missense_variant": VariantConsequence.MISSENSE,
    "stop_gained": VariantConsequence.NONSENSE,
    "splice_donor_variant": VariantConsequence.SPLICE_DONOR,
    "splice_acceptor_variant": VariantConsequence.SPLICE_ACCEPTOR,
    "splice_region_variant": VariantConsequence.SPLICE_REGION,
    "splice_donor_5th_base_variant": VariantConsequence.SPLICE_REGION,
    "splice_donor_region_variant": VariantConsequence.SPLICE_REGION,
    "splice_polypyrimidine_tract_variant": VariantConsequence.SPLICE_REGION,
    "intron_variant": VariantConsequence.INTRON,
    "5_prime_UTR_variant": VariantConsequence.FIVE_PRIME_UTR,
    "3_prime_UTR_variant": VariantConsequence.THREE_PRIME_UTR,
    "upstream_gene_variant": VariantConsequence.UPSTREAM,
    "downstream_gene_variant": VariantConsequence.DOWNSTREAM,
    "regulatory_region_variant": VariantConsequence.REGULATORY,
    "TF_binding_site_variant": VariantConsequence.REGULATORY,
    "intergenic_variant": VariantConsequence.INTERGENIC,
}


def _classify(term: str) -> VariantConsequence:
    return _CONSEQUENCE_MAP.get(term, VariantConsequence.OTHER)


def fetch_vep(client: httpx.Client, transcript: str, hgvs_c: str, gene_symbol: str) -> VepEvidence:
    """Query VEP by HGVS notation, e.g. NM_007294.4:c.181T>G."""

    hgvs = f"{transcript}:{hgvs_c}"
    url = f"{VEP_BASE}/vep/human/hgvs/{hgvs}"
    # vcf_string=1 returns the forward-strand VCF representation
    # (CHROM-POS-REF-ALT). This is essential: for minus-strand genes (e.g. BRCA1)
    # the top-level allele_string is in coding orientation and would build the
    # wrong gnomAD id, producing a false "not found".
    params = {"content-type": "application/json", "hgvs": "1", "canonical": "1", "vcf_string": "1"}
    web_url = f"https://www.ensembl.org/Homo_sapiens/Tools/VEP?hgvs={hgvs}"
    source = Source(name="Ensembl VEP (GRCh38)", url=web_url, retrieved_at=now_iso(), detail=hgvs)

    try:
        data = get_json(client, url, params=params)
    except httpx.HTTPStatusError as exc:
        # VEP returns 400 for HGVS it can't validate against the transcript.
        if exc.response.status_code == 400:
            source.detail = f"{hgvs} (VEP rejected as invalid HGVS)"
            return VepEvidence(
                found=False,
                consequence=VariantConsequence.OTHER,
                consequence_raw="",
                gene_symbol=gene_symbol,
                protein_hgvs=None,
                source=source,
            )
        raise

    if not isinstance(data, list) or not data:
        return VepEvidence(
            found=False,
            consequence=VariantConsequence.OTHER,
            consequence_raw="",
            gene_symbol=gene_symbol,
            protein_hgvs=None,
            source=source,
        )

    record = data[0]
    most_severe = record.get("most_severe_consequence", "")

    # Prefer vcf_string (forward-strand CHROM-POS-REF-ALT) for coordinates so the
    # gnomAD id is correct regardless of gene strand. Fall back to the top-level
    # allele_string only if vcf_string is absent.
    chrom = pos = ref = alt = None
    vcf_string = record.get("vcf_string")
    if vcf_string:
        first = vcf_string.split(",")[0]  # multiallelic sites are comma-joined
        parts = first.split("-")
        if len(parts) == 4:
            chrom, pos_s, ref, alt = parts
            pos = int(pos_s)
    if chrom is None:
        chrom = str(record.get("seq_region_name")) if record.get("seq_region_name") is not None else None
        pos = record.get("start")
        allele_string = record.get("allele_string", "")
        if "/" in allele_string:
            ref, alt = allele_string.split("/", 1)

    # Pick the transcript consequence that carries the most-severe term, preferring
    # the canonical transcript. Fall back to the first available.
    chosen = None
    for tc in record.get("transcript_consequences", []):
        terms = tc.get("consequence_terms", [])
        if most_severe in terms:
            chosen = tc
            if tc.get("canonical") == 1:
                break
    if chosen is None and record.get("transcript_consequences"):
        chosen = record["transcript_consequences"][0]

    protein_hgvs = None
    resolved_gene = gene_symbol
    if chosen is not None:
        raw_hgvsp = chosen.get("hgvsp")
        if raw_hgvsp and ":" in raw_hgvsp:
            protein_hgvs = raw_hgvsp.split(":", 1)[1]
        resolved_gene = chosen.get("gene_symbol", gene_symbol)

    return VepEvidence(
        found=True,
        consequence=_classify(most_severe),
        consequence_raw=most_severe,
        gene_symbol=resolved_gene,
        protein_hgvs=protein_hgvs,
        chrom=chrom,
        pos=int(pos) if pos is not None else None,
        ref=ref,
        alt=alt,
        source=source,
    )
