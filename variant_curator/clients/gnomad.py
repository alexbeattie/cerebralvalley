"""gnomAD GraphQL client.

Fetches ancestry-stratified allele frequency for a variant id (CHROM-POS-REF-ALT,
GRCh38). We use the gnomAD v4 `joint` dataset (combined exomes + genomes) and its
FAF95 (filtering allele frequency, 95% CI) which is the value ACMG BA1/BS1 should
key off, not the raw AF.

A "Variant not found" result is meaningful evidence in its own right (absence from a
large reference population supports PM2_supporting), so we return found=False rather
than raising.
"""

from __future__ import annotations

import httpx

from ..http import now_iso, post_json
from ..models import GnomadEvidence, PopulationFrequency, Source

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"

# gnomAD genetic-ancestry group ids -> human labels. Sex-split ids (…_XX/_XY)
# are filtered out so we report one row per ancestry group.
ANCESTRY_LABELS = {
    "afr": "African / African-American",
    "amr": "Admixed American / Latino",
    "asj": "Ashkenazi Jewish",
    "eas": "East Asian",
    "fin": "Finnish (European)",
    "nfe": "Non-Finnish European",
    "mid": "Middle Eastern",
    "sas": "South Asian",
    "ami": "Amish",
    "remaining": "Remaining (unassigned)",
}

_QUERY = """
query VariantAncestry($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id
    joint {
      ac
      an
      faf95 { popmax popmax_population }
      populations { id ac an }
    }
    genome { ac an af }
    exome { ac an af }
  }
}
"""


def fetch_gnomad(client: httpx.Client, variant_id: str) -> GnomadEvidence:
    web_url = f"https://gnomad.broadinstitute.org/variant/{variant_id}?dataset={DATASET}"
    source = Source(name="gnomAD v4 (joint)", url=web_url, retrieved_at=now_iso(), detail=variant_id)

    payload = {"query": _QUERY, "variables": {"variantId": variant_id, "dataset": DATASET}}
    data = post_json(client, GNOMAD_API, payload)
    variant = (data.get("data") or {}).get("variant")

    if not variant:
        return GnomadEvidence(found=False, variant_id=variant_id, dataset=DATASET, source=source)

    joint = variant.get("joint")
    genome = variant.get("genome")
    exome = variant.get("exome")

    # Prefer joint; fall back to genome then exome for the global figures.
    if joint:
        global_ac = joint.get("ac", 0)
        global_an = joint.get("an", 0)
        pop_source = joint.get("populations", [])
        faf = joint.get("faf95") or {}
    elif genome:
        global_ac = genome.get("ac", 0)
        global_an = genome.get("an", 0)
        pop_source = []
        faf = {}
    elif exome:
        global_ac = exome.get("ac", 0)
        global_an = exome.get("an", 0)
        pop_source = []
        faf = {}
    else:
        return GnomadEvidence(found=False, variant_id=variant_id, dataset=DATASET, source=source)

    populations: list[PopulationFrequency] = []
    for p in pop_source:
        pid = p.get("id", "")
        # Keep only top-level genetic-ancestry groups. This drops sex totals
        # ("XX"/"XY"), sex-split rows ("nfe_XX"), and any subcontinental splits,
        # leaving one clean row per ancestry group.
        if pid not in ANCESTRY_LABELS:
            continue
        populations.append(
            PopulationFrequency(
                ancestry_id=pid,
                ancestry_label=ANCESTRY_LABELS.get(pid, pid),
                ac=p.get("ac", 0),
                an=p.get("an", 0),
            )
        )

    global_af = (global_ac / global_an) if global_an else None

    return GnomadEvidence(
        found=True,
        variant_id=variant.get("variant_id", variant_id),
        dataset=DATASET,
        global_ac=global_ac,
        global_an=global_an,
        global_af=global_af,
        faf95_popmax=faf.get("popmax"),
        faf95_popmax_ancestry=faf.get("popmax_population"),
        populations=populations,
        source=source,
    )
