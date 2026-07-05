"""HPO (Human Phenotype Ontology) client via the Jax ontology API.

Two jobs:
  1. gene symbol -> the set of HPO phenotype ids known for that gene's diseases
     (this is the knowledge we match the patient against),
  2. free text -> HPO term (so a curator can type "seizures" instead of looking
     up HP:0001250).

Public API (ontology.jax.org), no key required. Every result carries an
openable source URL so a clinician can verify the gene-disease link.
"""

from __future__ import annotations

import re

import httpx

from ..http import get_json, now_iso
from ..models import GenePhenotypeKnowledge, Phenotype, Source

JAX_API = "https://ontology.jax.org/api"

# Ontology roots carry no clinical meaning.
_ROOT_NAMES = {"All", "Phenotypic abnormality"}

# Organ-system "container" nodes (e.g. "Abnormality of the cardiovascular
# system", "Abnormal nervous system physiology"). Matching a patient feature to
# a gene THROUGH one of these is spurious -- it only means the gene has some
# annotation in that whole organ system. We exclude them from ancestor matching,
# while keeping real (if broad) phenotypes like "Seizure" or "Cardiomyopathy".
_SYSTEM_CONTAINER_RE = re.compile(
    r"^(abnormality of (the )?.+ system|abnormal .+ system (morphology|physiology))$",
    re.IGNORECASE,
)


def resolve_term(client: httpx.Client, query: str) -> Phenotype | None:
    """Resolve free text (or an HP:xxxxx id) to a single HPO term."""

    q = query.strip()
    if q.upper().startswith("HP:"):
        data = get_json(client, f"{JAX_API}/hp/terms/{q.upper()}")
        if isinstance(data, dict) and data.get("id"):
            return Phenotype(hpo_id=data["id"], label=data.get("name", q))
        return None

    data = get_json(client, f"{JAX_API}/hp/search", params={"q": q, "max": "1"})
    terms = data.get("terms") or data.get("results") or [] if isinstance(data, dict) else []
    if not terms:
        return None
    top = terms[0]
    return Phenotype(hpo_id=top["id"], label=top.get("name", q))


def _gene_id(client: httpx.Client, gene: str) -> str | None:
    """Map a gene symbol to its NCBIGene id via the network search endpoint."""

    data = get_json(client, f"{JAX_API}/network/search/GENE", params={"q": gene, "limit": "10"})
    results = data.get("results", []) if isinstance(data, dict) else []
    for r in results:
        if r.get("name", "").upper() == gene.upper():
            return r.get("id")
    return None


# Ancestors are shared across genes/patients within a run, so cache them to keep
# the is_a graph walk to one API call per distinct term.
_ANCESTOR_CACHE: dict[str, set[str]] = {}


def ancestor_ids(client: httpx.Client, hpo_id: str) -> set[str]:
    """The term itself plus all its HPO is_a ancestors (excluding the roots).

    Used so a gene annotated to a broad term (e.g. Seizure) still explains a
    patient's more specific feature (e.g. Focal-onset seizure).
    """

    if hpo_id in _ANCESTOR_CACHE:
        return _ANCESTOR_CACHE[hpo_id]

    ids = {hpo_id}
    try:
        data = get_json(client, f"{JAX_API}/hp/terms/{hpo_id}/ancestors")
        if isinstance(data, list):
            # Drop ontology roots and organ-system container nodes: matching
            # through "Abnormality of the cardiovascular system" would let any
            # gene with any cardiac annotation "explain" a specific cardiac
            # feature, which is exactly the overfit this tool exists to prevent.
            for t in data:
                name = t.get("name", "")
                if not t.get("id") or name in _ROOT_NAMES or _SYSTEM_CONTAINER_RE.match(name):
                    continue
                ids.add(t["id"])
    except Exception:
        pass  # fall back to exact-only matching for this term
    _ANCESTOR_CACHE[hpo_id] = ids
    return ids


# Gene knowledge is stable within a run and often queried repeatedly (a gene
# reported twice, or the same decoy panel across many eval cases), so cache it.
_GENE_CACHE: dict[str, GenePhenotypeKnowledge] = {}


def fetch_gene_phenotypes(client: httpx.Client, gene: str) -> GenePhenotypeKnowledge:
    """All HPO phenotypes and diseases associated with a gene."""

    key = gene.upper()
    if key in _GENE_CACHE:
        return _GENE_CACHE[key]

    web_url = f"https://hpo.jax.org/browse/gene/{gene}"
    source = Source(name="HPO (Jax)", url=web_url, retrieved_at=now_iso(), detail=gene)

    gene_id = _gene_id(client, gene)
    if not gene_id:
        result = GenePhenotypeKnowledge(gene=gene, found=False, source=source)
        _GENE_CACHE[key] = result
        return result

    source.detail = f"{gene} ({gene_id})"
    source.url = f"https://ontology.jax.org/api/network/annotation/{gene_id}"

    data = get_json(client, f"{JAX_API}/network/annotation/{gene_id}")
    if not isinstance(data, dict) or "phenotypes" not in data:
        result = GenePhenotypeKnowledge(gene=gene, found=False, source=source)
        _GENE_CACHE[key] = result
        return result

    phenotype_labels = {p["id"]: p.get("name", p["id"]) for p in data.get("phenotypes", []) if p.get("id")}
    phenotype_ids = set(phenotype_labels)
    diseases = [d.get("name", "") for d in data.get("diseases", []) if d.get("name")]

    result = GenePhenotypeKnowledge(
        gene=gene,
        found=bool(phenotype_ids),
        diseases=diseases,
        phenotype_ids=phenotype_ids,
        phenotype_labels=phenotype_labels,
        source=source,
    )
    _GENE_CACHE[key] = result
    return result
