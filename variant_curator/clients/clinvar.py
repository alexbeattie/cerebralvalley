"""ClinVar client via NCBI E-utilities.

Returns prior germline classifications for a variant, the ClinVar review status
(mapped to the familiar star rating), and whether submitters conflict. This is the
evidence behind ACMG PP5 / BP6 (used cautiously) and, more importantly, the sanity
check a curator wants: "has anyone called this before, and did they agree?"

Public E-utilities, no key required (rate-limited to ~3 req/s without a key).
"""

from __future__ import annotations

from urllib.parse import quote_plus

import httpx

from ..http import get_json, now_iso
from ..models import ClinvarEvidence, ClinvarSubmission, Source

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

_STAR_BY_REVIEW_STATUS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
}


def _stars(review_status: str) -> int:
    return _STAR_BY_REVIEW_STATUS.get(review_status.strip().lower(), 0)


def fetch_clinvar(client: httpx.Client, gene: str, hgvs_c: str) -> ClinvarEvidence:
    term = f'{gene}[gene] AND {hgvs_c}'
    web_url = "https://www.ncbi.nlm.nih.gov/clinvar/?term=" + quote_plus(term)
    source = Source(name="ClinVar", url=web_url, retrieved_at=now_iso(), detail=term)

    search = get_json(
        client,
        f"{EUTILS}/esearch.fcgi",
        params={"db": "clinvar", "term": term, "retmode": "json", "retmax": "20"},
    )
    id_list = (search.get("esearchresult") or {}).get("idlist", [])
    if not id_list:
        return ClinvarEvidence(found=False, source=source)

    summary = get_json(
        client,
        f"{EUTILS}/esummary.fcgi",
        params={"db": "clinvar", "id": ",".join(id_list), "retmode": "json"},
    )
    result = summary.get("result", {})
    uids = result.get("uids", [])

    submissions: list[ClinvarSubmission] = []
    hgvs_lower = hgvs_c.lower().replace(" ", "")
    for uid in uids:
        rec = result.get(uid, {})
        title = rec.get("title", "")
        # Keep only records whose title actually contains this coding change, so a
        # broad gene search doesn't drag in neighbouring variants.
        if hgvs_lower not in title.lower().replace(" ", ""):
            continue
        gc = rec.get("germline_classification", {}) or {}
        desc = gc.get("description", "")
        if not desc:
            continue
        review_status = gc.get("review_status", "")
        traits = gc.get("trait_set", []) or []
        condition = next(
            (t.get("trait_name") for t in traits if t.get("trait_name") and t["trait_name"] != "not provided"),
            "not specified",
        )
        submissions.append(
            ClinvarSubmission(
                accession=rec.get("accession", uid),
                clinical_significance=desc,
                review_status=review_status,
                star_rating=_stars(review_status),
                condition=condition,
                last_evaluated=gc.get("last_evaluated", ""),
            )
        )

    if not submissions:
        return ClinvarEvidence(found=False, source=source)

    # Prefer the highest-reviewed record for the aggregate call.
    best = max(submissions, key=lambda s: s.star_rating)
    has_conflict = "conflict" in best.review_status.lower() or "conflicting" in best.clinical_significance.lower()

    return ClinvarEvidence(
        found=True,
        submissions=submissions,
        has_conflict=has_conflict,
        aggregate_significance=best.clinical_significance,
        source=source,
    )
