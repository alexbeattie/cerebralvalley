"""trio_prioritizer: inheritance-aware prioritization of non-coding variants in trios.

The clinical need (from Matt Deardorff, clinical geneticist): in rare-disease trios,
two patterns dominate. A **de novo** variant present in the affected child but neither
parent is the usual mechanism for severe dominant disease. A **recessive / compound
heterozygous** case has an affected child and two unaffected carrier parents — and the
elusive "second hit" is increasingly a deep-intronic, UTR, or regulatory variant that a
coding-only pipeline never surfaces. Finding and scoring that non-coding allele, in the
context of inheritance, is the unmet need.

This package supplies the inheritance logic, compound-het pairing by parent-of-origin,
and an integrated ranking on top of the repo's existing non-coding *scoring* engine
(`variant_curator.clients.alphagenome`). It interleaves de novo, compound-het pairs, and
homozygous-recessive candidates into one ranked list, so a non-coding candidate is
weighed alongside the other potentially-causal variants already on the table.

Offline-first: the whole `demo` and the entire test suite run on a CPU laptop with **no
network and no API key** on a built-in synthetic trio. The AlphaGenome scorer is
injectable, so tests never touch the network; when it is unavailable we mark the score
unavailable with a reason and still rank on inheritance — we never fabricate a score.

No PHI (synthetic / coordinate data only). Designed so the same engine could run behind a
clinical firewall unchanged. Research use only: this is decision-support, not a diagnosis.
"""

from __future__ import annotations

__version__ = "0.1.0"

# Surfaced wherever a pathogenicity score is shown (CLI, JSON, rationale strings).
RESEARCH_USE_CAVEAT = (
    "RESEARCH USE ONLY — inheritance-aware prioritization with research-model "
    "(AlphaGenome) non-coding predictions. Decision-support, not a diagnosis; "
    "de novo and phasing calls require QC/coverage and orthogonal confirmation."
)
