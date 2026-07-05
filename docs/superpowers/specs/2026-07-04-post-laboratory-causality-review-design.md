# Post-Laboratory Causality Review — Design Spec

**Date:** 2026-07-04
**Status:** Approved for implementation
**Scope target:** Demoable end-to-end slice, hours not days (live hackathon).

## Problem

A genomics laboratory has already analyzed a genome and produced a report with
~2–5 candidate variants. The clinician's job is the *second* stage of
interpretation: given those variants and the full patient story, determine
**"Does this reported variant actually explain this specific patient?"**

This is a **reverse matching** problem, and the inverse of what `variant_curator`
(the existing lab-facing MVP) does. Lab pipelines ask "which variants might
explain this phenotype?" This tool asks "given these variants, which one actually
explains this patient?"

It must protect the clinician from two cognitive traps:

1. **Overfitting** — convincing yourself a disease matches when the overlap is
   only partial. Counter: always compare phenotype objectively; show supporting,
   conflicting, and missing features; never assert from impression.
2. **Partial explanations** — a variant explains most but not all features. The
   remainder implies either a broader-than-reported disease spectrum *or* a
   second independent diagnosis (~5% of solved cases get two molecular
   diagnoses). Counter: explicitly detect unexplained phenotype clusters and
   recommend genome re-analysis when appropriate.

## Non-goals (explicit scope cuts for the hours-scale slice)

- **No real PDF binary parsing.** Inputs are pasted text or `.txt`/`.md` upload.
  A "Load example case" button preloads the demo case.
- **No live PubMed literature synthesis.** Use ClinVar citations + link-outs.
- **No directive management recommendations.** Only "suggested topics for
  clinician review," clearly framed as non-directive.
- **No genes outside the existing 30-gene allowlist** (`variant_curator/genes.py`).
  Anything else is refused gracefully with a clear message.
- **Evaluator is a scoped grounding-check pass**, not an iterate-until-perfect
  agent loop (see below).

## Users

Experienced pediatric geneticists performing post-laboratory causality review.
Target audience for tone, depth, and the honesty/uncertainty bar.

## The demo case (drives every design decision)

A synthetic, PHI-free pediatric report lists 3 candidate variants across genes in
the existing allowlist. The patient's notes describe two phenotype clusters:

- **Cluster 1 — neurodevelopmental/epilepsy:** cleanly explained by an *SCN1A*
  variant (Dravet/epilepsy).
- **Cluster 2 — connective tissue:** tall stature, arachnodactyly, lens
  dislocation — **not** explained by *SCN1A*, matching *FBN1* (Marfan) instead.

Expected system behavior: rank the *SCN1A* variant as the strong explanation for
cluster 1, refuse to overfit it onto cluster 2, surface cluster 2 as an
**unexplained phenotype cluster**, and recommend re-analysis / flag a possible
**second molecular diagnosis**. This single case exercises both failure modes and
the full ranked causality output.

## Architecture

New package `causality_review/` alongside the existing `variant_curator/`, reusing
its plumbing. Nothing in `variant_curator/` needs to change; we import it.

```
causality_review/
  models.py          data models (see below); every item carries provenance
  extract.py         Claude calls: report -> candidate variants; notes -> HPO phenotype
  hpo.py             HPO/Jax ontology API client: NCBI gene id -> disease HPO terms (cited)
  match.py           deterministic phenotype matching + confidence rubric (no LLM)
  narrate.py         Claude call: grounded per-variant narrative, constrained to retrieved evidence
  evaluate.py        Claude call: house-rules grounding check over the assembled review
  pipeline.py        orchestrates extract -> retrieve -> match -> narrate -> evaluate
  server.py          FastAPI app: POST /review -> JSON; serves the static UI
  web/index.html     single 3-pane page, vanilla JS/CSS, no build step
  examples/          the synthetic demo case (report.txt, notes.txt)
```

### Data flow

1. **Ingest** — UI posts `{ report_text, notes_text }` to `POST /review`.
2. **Extract** (`extract.py`, Claude, grounded):
   - Report -> `CandidateVariant[]`: gene, HGVS, zygosity, lab-reported call, and
     the **verbatim report snippet** each was drawn from.
   - Notes -> `PatientPhenotype`: a list of `PhenotypeObservation` (HPO term id +
     label + the **verbatim note passage** it came from). Claude maps free text to
     HPO ids; the passage is retained so every term is traceable.
3. **Retrieve** (deterministic, cited):
   - `hpo.py`: for each candidate variant's gene, resolve NCBI gene id (already in
     `genes.py`) -> associated diseases + their HPO terms from the Jax HPO API
     (`https://ontology.jax.org/api/network/annotation/NCBIGene:{id}`). Keep the
     source URL. This is the **objective gene->phenotype source** — never Claude's
     memory.
   - Reuse `variant_curator.pipeline.assemble_evidence` for VEP/gnomAD/ClinVar
     (the molecular axis: consequence, rarity, prior classification).
4. **Match** (`match.py`, deterministic, no LLM):
   - For each variant: `supporting = patient_terms ∩ disease_terms`,
     `missing = disease_terms − patient_terms`.
   - Across all variants: `unexplained = patient_terms − ∪ disease_terms`,
     grouped into clusters (by shared HPO ancestry / organ system where feasible;
     otherwise a flat list).
   - Confidence tier from a **fixed, visible rubric** combining phenotype overlap
     fraction and molecular evidence (ClinVar significance, rarity). The rubric is
     data, not prose, so every score change is explainable by pointing at the
     input that moved it.
5. **Narrate** (`narrate.py`, Claude, grounded): per variant, write the causality
   narrative constrained to only the retrieved terms/snippets/molecular evidence.
   Every statement tagged one of **Known / Likely / Possible / Speculative /
   Unknown**. Model is instructed to abstain rather than assert beyond evidence.
6. **Evaluate** (`evaluate.py`, Claude, fresh context): a grounding-check pass
   that reads the assembled review + the evidence it is allowed to cite, and flags
   any statement not traceable to patient records, the lab report, or a retrieved
   source. Flags are surfaced in the UI, not silently dropped.
7. **Output** — ranked `CausalityReview` JSON -> 3-pane UI.

### Confidence rubric (deterministic, shown in UI)

Tiers map to a scored combination; exact thresholds finalized in the plan, but the
inputs are fixed:

- **phenotype overlap fraction** = |supporting| / |disease_terms| (capped, with a
  floor on |disease_terms| to avoid tiny-denominator noise).
- **molecular support** = ClinVar significance + star rating, gnomAD rarity, VEP
  consequence in-scope.
- **unexplained penalty** = presence of a large unexplained cluster lowers the
  claim that this variant explains *the patient* (vs. explaining *a* phenotype).

Every variant card shows the rubric inputs so the score is never opaque.

## Data models (sketch)

```python
@dataclass
class ReportSnippet:            # provenance for an extracted claim
    text: str                  # verbatim passage
    source: str                # "laboratory report" | "clinical notes"

@dataclass
class CandidateVariant:
    gene: str
    hgvs_c: str
    zygosity: str | None
    reported_significance: str | None
    snippet: ReportSnippet

@dataclass
class PhenotypeObservation:
    hpo_id: str                # e.g. "HP:0001250"
    label: str                 # "Seizure"
    passage: ReportSnippet     # note text it was drawn from

@dataclass
class DiseasePhenotype:        # retrieved, cited
    disease: str
    hpo_terms: list[tuple[str, str]]   # (id, label)
    source: Source             # HPO/Jax URL, reused from variant_curator.models

@dataclass
class VariantAssessment:
    variant: CandidateVariant
    disease: DiseasePhenotype
    molecular: EvidenceBundle          # from variant_curator
    supporting: list[PhenotypeObservation]
    missing: list[tuple[str, str]]
    confidence_tier: str               # Known/Likely/Possible/Speculative/Unknown
    rubric: dict                       # the inputs that produced the tier
    narrative: str                     # grounded, statement-tagged
    grounding_flags: list[str]         # from evaluate.py

@dataclass
class UnexplainedCluster:
    observations: list[PhenotypeObservation]
    recommendation: str        # e.g. "consider genome re-analysis / second diagnosis"

@dataclass
class CausalityReview:
    assessments: list[VariantAssessment]   # ranked, best explanation first
    unexplained: list[UnexplainedCluster]
    sources: list[Source]
```

## Interface

Single page, three panes:

- **Left** — laboratory report (the pasted/loaded text; extracted variants
  highlighted with their snippets).
- **Middle** — clinical notes (extracted phenotype passages highlighted).
- **Right** — **ranked causality review**. One card per candidate variant, ranked
  best-explanation-first. Each card expands to an evidence page:
  - summary + confidence tier (with rubric inputs visible)
  - supporting findings (patient ∩ disease), each linking to its note passage
  - conflicting / missing findings (disease − patient)
  - phenotype overlap fraction
  - molecular evidence (VEP/gnomAD/ClinVar with source links)
  - grounded narrative with Known/Likely/Possible/Speculative/Unknown tags
  - grounding-check flags (if any)
  - sources
  - A prominent banner above the cards when an unexplained cluster exists:
    "N documented features remain unexplained — consider genome re-analysis /
    possible second molecular diagnosis," listing the features.

Aesthetic: calm, organized, trustworthy. Evidence-first. No hidden scoring.

## House rules (enforced structurally, not just prompted)

- **Never invent evidence.** Extraction retains verbatim snippets; matching is
  deterministic over retrieved term sets; narrative is constrained to those sets
  and abstains otherwise; the evaluator flags anything untraceable.
- **Every conclusion traces to** patient records, the lab report, or a retrieved
  source (HPO/Jax, ClinVar, gnomAD, VEP) — each carries a `Source`/snippet.
- **Always distinguish** Known / Likely / Possible / Speculative / Unknown.
- **Never hide uncertainty; never overstate confidence.** Missing features and
  unexplained clusters are first-class output, not footnotes.
- **Always explain why a score is what it is** — the rubric inputs are shown.

## Error handling

- Gene outside allowlist -> variant card shows a clear "out of scope" state, still
  listed (not dropped), phenotype match skipped.
- HPO/Jax API failure -> variant retains molecular evidence; phenotype match marked
  "unavailable (source unreachable)," never fabricated.
- Claude extraction returning nothing / malformed -> surfaced as an ingest error,
  not a silent empty review.
- Missing `ANTHROPIC_API_KEY` -> server refuses to start reasoning endpoints with a
  clear message; static UI + example still viewable.

## Testing

- **Deterministic core (`match.py`, `hpo.py` parsing)** — unit tests with fixed
  term sets and a recorded HPO/Jax response fixture. This is the part that must be
  correct and is cheap to test without the LLM.
- **Rubric** — table tests: given overlap fraction + molecular inputs, assert the
  tier and that the rubric dict explains it.
- **End-to-end demo case** — a smoke test that runs the example case through the
  pipeline (LLM live) and asserts structural invariants: 3 variants extracted, the
  SCN1A variant ranks highest for the epilepsy cluster, and a connective-tissue
  unexplained cluster is surfaced. Kept as a manual/opt-in test given it hits the
  API.

## Stack / dependencies

- Python 3.12 (existing `.venv`), `httpx` (present), add `anthropic` SDK and
  `fastapi` + `uvicorn`. No frontend build step (vanilla HTML/CSS/JS).
- `ANTHROPIC_API_KEY` read from env / `.env`. Claude model: latest available
  (`claude-fable-5` if permitted, else current default) — confirm at build time.
- Reuses `variant_curator` package unchanged.

## What "good enough for the demo" means

An experienced pediatric geneticist watching the demo sees the tool: rank the
variants objectively, show supporting *and* conflicting evidence side by side,
refuse to overfit the SCN1A variant to the connective-tissue features, surface the
unexplained cluster, and recommend re-analysis — with every claim linking back to a
source. That is the bar.
