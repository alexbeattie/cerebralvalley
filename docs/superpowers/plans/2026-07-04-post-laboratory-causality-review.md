# Post-Laboratory Causality Review — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a demoable clinician-facing causality-review web app: upload a genetics report + clinical notes, get a ranked, evidence-traceable analysis of which reported variant actually explains this patient — surfacing supporting *and* conflicting evidence, unexplained phenotype clusters, and possible dual diagnoses.

**Architecture:** A new `causality_review/` package alongside the existing `variant_curator/` (reused unchanged for the VEP→gnomAD→ClinVar molecular axis). A deterministic phenotype-matching core (`match.py`) does the reverse-matching math; Claude (`claude-fable-5`) handles only grounded extraction and narrative, constrained to retrieved evidence. Gene→phenotype ground truth comes from the HPO/Jax ontology API (real, keyless, cited). A FastAPI server serves one vanilla three-pane HTML page.

**Tech Stack:** Python 3.12, `httpx` (present), `anthropic` SDK, `fastapi` + `uvicorn`, vanilla HTML/CSS/JS (no build step). Spec: `docs/superpowers/specs/2026-07-04-post-laboratory-causality-review-design.md`.

## Global Constraints

- **Model:** `claude-fable-5` by default, read from env `CAUSALITY_MODEL` (fallback to that literal). Never hardcode elsewhere.
- **Fable 5 API rules (all LLM calls go through `causality_review/llm.py`):** omit the `thinking` parameter entirely (thinking is always on; `{type:"disabled"}` 400s); never send `budget_tokens`, `temperature`, `top_p`, `top_k`; no assistant-turn prefill. Use `output_config={"effort": ..., "format": {...}}` for structured JSON. Always pass `betas=["server-side-fallback-2026-06-01"]` + `fallbacks=[{"model":"claude-opus-4-8"}]` and check `stop_reason == "refusal"` before reading content — genomics/clinical text can trip Fable's bio/cyber classifiers.
- **Structured-output JSON Schema limits:** every object must set `additionalProperties: false` and list all keys in `required`; no `minLength`/`maxLength`/`minimum`/`maximum`, no recursive `$ref`.
- **House rules (enforced in code, not just prompts):** every extracted claim keeps its verbatim source snippet; matching is deterministic and its inputs are returned in the response; the narrative is constrained to retrieved terms and abstains otherwise; the evaluator flags anything untraceable. Confidence tiers are exactly `Known | Likely | Possible | Speculative | Unknown`.
- **Scope:** only genes in `variant_curator/genes.py` (30-gene allowlist); anything else is surfaced as an explicit "out of scope" state, never dropped or guessed. No real PDF binary parsing (text/paste only). No PHI — synthetic example case only.
- **Env:** the app reads `ANTHROPIC_API_KEY` from the environment (or an `ant auth login` profile). If neither the key nor a profile is present, reasoning endpoints return a clear 503; the static UI + example still load.
- **Provenance model reuse:** use `variant_curator.models.Source(name, url, retrieved_at, detail)` for all retrieved-source citations.

---

### Task 1: Project scaffolding + dependencies + the Claude wrapper

Everything else calls `causality_review/llm.py`. Build it first with a live smoke test so the Fable 5 wiring (structured output + refusal fallback) is proven before layering logic on top.

**Files:**
- Create: `causality_review/__init__.py`
- Create: `causality_review/llm.py`
- Modify: `requirements.txt`
- Test: `tests/__init__.py`, `tests/test_llm_smoke.py`

**Interfaces:**
- Produces: `causality_review.llm.call_json(system: str, user: str, schema: dict, *, effort: str = "high", max_tokens: int = 8000) -> dict` — one grounded structured-output call; returns the parsed JSON object. Raises `llm.LLMUnavailable` (no credentials) or `llm.LLMRefused` (whole fallback chain refused).
- Produces: `causality_review.llm.credentials_available() -> bool`.
- Produces: `causality_review.llm.MODEL` (str, resolved from env).

- [ ] **Step 1: Add dependencies**

Overwrite `requirements.txt`:

```
httpx==0.27.2
anthropic>=0.69
fastapi>=0.115
uvicorn>=0.30
```

Install:

```bash
./.venv/bin/pip install -r requirements.txt
```

- [ ] **Step 2: Write `causality_review/__init__.py`**

```python
"""Post-laboratory causality review: reverse-matching reported variants to a patient."""
```

- [ ] **Step 3: Write the failing smoke test**

`tests/__init__.py`: empty file.

`tests/test_llm_smoke.py`:

```python
"""Live smoke test for the Claude wrapper. Skipped without credentials."""

import pytest

from causality_review import llm


@pytest.mark.skipif(not llm.credentials_available(), reason="no Anthropic credentials")
def test_call_json_returns_validated_object():
    schema = {
        "type": "object",
        "additionalProperties": False,
        "properties": {"capital": {"type": "string"}},
        "required": ["capital"],
    }
    out = llm.call_json(
        system="You answer with facts only.",
        user="What is the capital of France? Reply as JSON.",
        schema=schema,
        effort="low",
        max_tokens=200,
    )
    assert isinstance(out, dict)
    assert "paris" in out["capital"].lower()
```

- [ ] **Step 4: Run it to confirm it fails**

Run: `./.venv/bin/python -m pytest tests/test_llm_smoke.py -v`
Expected: FAIL — `ModuleNotFoundError: causality_review.llm` (or collection error).

- [ ] **Step 5: Write `causality_review/llm.py`**

```python
"""Single choke point for all Claude calls.

Every LLM call in this app goes through call_json(). It centralizes the Fable 5
request rules (no thinking/sampling params, structured output via
output_config.format, refusal fallback to Opus 4.8) and the credential check, so
no other module has to know them.
"""

from __future__ import annotations

import json
import os

MODEL = os.environ.get("CAUSALITY_MODEL", "claude-fable-5")
FALLBACK_MODEL = "claude-opus-4-8"


class LLMUnavailable(RuntimeError):
    """No Anthropic credentials are configured."""


class LLMRefused(RuntimeError):
    """The model (and its fallback) declined the request."""


def _client():
    import anthropic

    # Generous timeout: Fable 5 turns can run minutes on hard tasks.
    return anthropic.Anthropic(timeout=600.0)


def credentials_available() -> bool:
    """True if an API key env var is set. (Profile-only auth still works at call
    time; this is a best-effort pre-check for the UI's 503 path.)"""
    return bool(os.environ.get("ANTHROPIC_API_KEY") or os.environ.get("ANTHROPIC_AUTH_TOKEN"))


def call_json(
    system: str,
    user: str,
    schema: dict,
    *,
    effort: str = "high",
    max_tokens: int = 8000,
) -> dict:
    """One grounded structured-output call. Returns the parsed JSON object.

    Raises LLMUnavailable if the SDK cannot authenticate, LLMRefused if the
    request (including the Opus 4.8 fallback) was declined by safety classifiers.
    """
    import anthropic

    try:
        client = _client()
        resp = client.beta.messages.create(
            model=MODEL,
            max_tokens=max_tokens,
            betas=["server-side-fallback-2026-06-01"],
            fallbacks=[{"model": FALLBACK_MODEL}],
            output_config={"effort": effort, "format": {"type": "json_schema", "schema": schema}},
            system=system,
            messages=[{"role": "user", "content": user}],
        )
    except anthropic.AuthenticationError as e:
        raise LLMUnavailable(str(e)) from e

    if resp.stop_reason == "refusal":
        detail = getattr(resp.stop_details, "explanation", "") or ""
        raise LLMRefused(f"Request declined by safety classifiers. {detail}".strip())

    text = next((b.text for b in resp.content if b.type == "text"), "")
    return json.loads(text)
```

- [ ] **Step 6: Run the smoke test**

Run: `./.venv/bin/python -m pytest tests/test_llm_smoke.py -v`
Expected: PASS if credentials are set, else SKIPPED. Both are acceptable. (If you have a key, `export ANTHROPIC_API_KEY=...` first to actually exercise the path.)

- [ ] **Step 7: Commit**

```bash
git add requirements.txt causality_review/__init__.py causality_review/llm.py tests/__init__.py tests/test_llm_smoke.py
git commit -m "feat(causality): Claude wrapper (Fable 5 structured output + refusal fallback)"
```

---

### Task 2: Data models

Plain dataclasses shared across the pipeline. No logic, no tests of their own — they're exercised by later tasks.

**Files:**
- Create: `causality_review/models.py`

**Interfaces:**
- Produces: `ReportSnippet`, `CandidateVariant`, `PhenotypeObservation`, `DiseasePhenotype`, `VariantAssessment`, `UnexplainedCluster`, `CausalityReview` — and `to_dict()` on `CausalityReview` for JSON serialization.
- Consumes: `variant_curator.models.Source`, `variant_curator.models.EvidenceBundle`.

- [ ] **Step 1: Write `causality_review/models.py`**

```python
"""Data models for the causality review. Every extracted claim carries its
verbatim source snippet; every retrieved fact carries a Source. Serialization to
plain dicts (for the JSON API) lives here so the server stays thin.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import Optional

from variant_curator.models import EvidenceBundle, Source

# Confidence tiers — the only allowed values, per the house rules.
TIERS = ("Known", "Likely", "Possible", "Speculative", "Unknown")


@dataclass
class ReportSnippet:
    text: str            # verbatim passage
    source: str          # "laboratory report" | "clinical notes"


@dataclass
class CandidateVariant:
    gene: str
    hgvs_c: str
    zygosity: Optional[str]
    reported_significance: Optional[str]
    snippet: ReportSnippet


@dataclass
class PhenotypeObservation:
    hpo_id: str          # e.g. "HP:0001250"; "" if the model could not map one
    label: str           # e.g. "Seizure"
    passage: ReportSnippet


@dataclass
class DiseasePhenotype:
    genes_disease_names: list[str]          # associated disease names for the gene
    hpo_terms: list[tuple[str, str]]        # (hpo_id, label) — the gene's known spectrum
    source: Optional[Source]
    available: bool = True                  # False if the HPO source was unreachable


@dataclass
class VariantAssessment:
    variant: CandidateVariant
    in_scope: bool
    disease: Optional[DiseasePhenotype]
    molecular: Optional[EvidenceBundle]
    supporting: list[PhenotypeObservation]              # patient ∩ disease
    missing: list[tuple[str, str]]                      # disease − patient (id, label)
    confidence_tier: str
    rubric: dict                                        # the inputs that produced the tier
    rank_score: float
    narrative: str = ""
    grounding_flags: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)      # e.g. "out of scope", errors


@dataclass
class UnexplainedCluster:
    observations: list[PhenotypeObservation]
    recommendation: str


@dataclass
class CausalityReview:
    assessments: list[VariantAssessment]     # ranked, best explanation first
    unexplained: list[UnexplainedCluster]
    sources: list[Source]
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return _to_jsonable(self)


def _to_jsonable(obj):
    if dataclasses.is_dataclass(obj):
        return {k: _to_jsonable(v) for k, v in dataclasses.asdict(obj).items()}
    if isinstance(obj, (list, tuple)):
        return [_to_jsonable(v) for v in obj]
    if isinstance(obj, dict):
        return {k: _to_jsonable(v) for k, v in obj.items()}
    return obj
```

> Note: `dataclasses.asdict` deep-converts nested dataclasses (including the reused `EvidenceBundle`/`Source`), so `to_dict()` yields a fully JSON-serializable tree.

- [ ] **Step 2: Sanity check import**

Run: `./.venv/bin/python -c "from causality_review.models import CausalityReview, TIERS; print(TIERS)"`
Expected: `('Known', 'Likely', 'Possible', 'Speculative', 'Unknown')`

- [ ] **Step 3: Commit**

```bash
git add causality_review/models.py
git commit -m "feat(causality): data models with provenance + JSON serialization"
```

---

### Task 3: HPO/Jax gene→phenotype client

The objective gene→phenotype source that makes "never invent evidence" true. Deterministic parsing; tested against a recorded fixture so it needs no network.

**Files:**
- Create: `causality_review/hpo.py`
- Test: `tests/test_hpo.py`, `tests/fixtures/hpo_scn1a.json`

**Interfaces:**
- Consumes: `variant_curator.http.get_client`, `variant_curator.models.Source`.
- Produces: `causality_review.hpo.parse_gene_annotation(payload: dict, *, retrieved_at: str, url: str) -> DiseasePhenotype`.
- Produces: `causality_review.hpo.fetch_gene_phenotype(client, ncbi_gene_id: int) -> DiseasePhenotype`.

- [ ] **Step 1: Record the fixture**

```bash
mkdir -p tests/fixtures
curl -s "https://ontology.jax.org/api/network/annotation/NCBIGene:6323" -o tests/fixtures/hpo_scn1a.json
./.venv/bin/python -c "import json;d=json.load(open('tests/fixtures/hpo_scn1a.json'));print(len(d['diseases']),'diseases',len(d['phenotypes']),'phenotypes')"
```

Expected: prints nonzero counts (e.g. `11 diseases 100+ phenotypes`). The API shape is `{"diseases":[{"id","name",...}], "phenotypes":[{"id","name"}]}`.

- [ ] **Step 2: Write the failing test**

`tests/test_hpo.py`:

```python
import json
from pathlib import Path

from causality_review import hpo

FIX = Path(__file__).parent / "fixtures" / "hpo_scn1a.json"


def test_parse_gene_annotation_extracts_terms_and_disease_names():
    payload = json.loads(FIX.read_text())
    dp = hpo.parse_gene_annotation(
        payload, retrieved_at="2026-07-04T00:00:00Z",
        url="https://ontology.jax.org/api/network/annotation/NCBIGene:6323",
    )
    assert dp.available is True
    # SCN1A → Dravet syndrome is in the associated disease list.
    assert any("dravet" in name.lower() for name in dp.genes_disease_names)
    # Phenotype terms are (HP:id, label) pairs; seizure-family term present.
    ids = {t[0] for t in dp.hpo_terms}
    assert any(i.startswith("HP:") for i in ids)
    assert dp.source is not None
    assert dp.source.name == "HPO (Jax)"


def test_parse_gene_annotation_handles_empty_payload():
    dp = hpo.parse_gene_annotation(
        {"diseases": [], "phenotypes": []},
        retrieved_at="2026-07-04T00:00:00Z", url="http://x",
    )
    assert dp.available is True   # reachable, just empty
    assert dp.hpo_terms == []
```

- [ ] **Step 3: Run it to confirm it fails**

Run: `./.venv/bin/python -m pytest tests/test_hpo.py -v`
Expected: FAIL — `ModuleNotFoundError: causality_review.hpo`.

- [ ] **Step 4: Write `causality_review/hpo.py`**

```python
"""Gene -> known phenotype spectrum from the HPO / Jax ontology API.

This is the objective, citable source for what a gene's disease *should* look
like. We never use the model's memory for this — the reverse-matching math runs
only over terms retrieved here.

Endpoint: GET https://ontology.jax.org/api/network/annotation/NCBIGene:{id}
Shape:    {"diseases":[{"id","name",...}], "phenotypes":[{"id","name"}]}
The phenotypes list is the gene-level union across its associated diseases.
"""

from __future__ import annotations

import datetime

from variant_curator.models import Source
from causality_review.models import DiseasePhenotype

BASE = "https://ontology.jax.org/api/network/annotation/NCBIGene:{gene_id}"


def parse_gene_annotation(payload: dict, *, retrieved_at: str, url: str) -> DiseasePhenotype:
    diseases = [d.get("name", "") for d in payload.get("diseases", []) if d.get("name")]
    terms: list[tuple[str, str]] = []
    seen: set[str] = set()
    for p in payload.get("phenotypes", []):
        hid, label = p.get("id", ""), p.get("name", "")
        if hid and hid not in seen:
            seen.add(hid)
            terms.append((hid, label))
    return DiseasePhenotype(
        genes_disease_names=diseases,
        hpo_terms=terms,
        source=Source(name="HPO (Jax)", url=url, retrieved_at=retrieved_at,
                      detail=f"{len(terms)} phenotype terms, {len(diseases)} diseases"),
        available=True,
    )


def fetch_gene_phenotype(client, ncbi_gene_id: int) -> DiseasePhenotype:
    """Fetch and parse. On any network/parse failure, return an unavailable
    DiseasePhenotype (never fabricate terms)."""
    url = BASE.format(gene_id=ncbi_gene_id)
    now = datetime.datetime.now(datetime.timezone.utc).isoformat()
    try:
        resp = client.get(url)
        resp.raise_for_status()
        return parse_gene_annotation(resp.json(), retrieved_at=now, url=url)
    except Exception as e:  # noqa: BLE001 — degrade gracefully, don't guess
        return DiseasePhenotype(
            genes_disease_names=[], hpo_terms=[],
            source=Source(name="HPO (Jax)", url=url, retrieved_at=now, detail=f"unreachable: {e}"),
            available=False,
        )
```

- [ ] **Step 5: Run the test to confirm it passes**

Run: `./.venv/bin/python -m pytest tests/test_hpo.py -v`
Expected: PASS (both tests).

- [ ] **Step 6: Commit**

```bash
git add causality_review/hpo.py tests/test_hpo.py tests/fixtures/hpo_scn1a.json
git commit -m "feat(causality): HPO/Jax gene->phenotype client (cited, offline-tested)"
```

---

### Task 4: Deterministic matching + confidence rubric (the core)

The part that must be correct and is cheapest to test without the LLM. No Claude, no network — pure functions over term sets.

**Files:**
- Create: `causality_review/match.py`
- Test: `tests/test_match.py`

**Interfaces:**
- Consumes: `PhenotypeObservation`, `DiseasePhenotype`, `CandidateVariant`, `VariantAssessment`, `UnexplainedCluster` from `causality_review.models`; `EvidenceBundle` from `variant_curator.models`.
- Produces: `causality_review.match.assess_variant(variant, disease, molecular, patient) -> VariantAssessment`.
- Produces: `causality_review.match.find_unexplained(patient, assessments) -> list[UnexplainedCluster]`.
- Produces: `causality_review.match.rank(assessments) -> list[VariantAssessment]` (returns a new sorted list, best first).

- [ ] **Step 1: Write the failing tests**

`tests/test_match.py`:

```python
from causality_review import match
from causality_review.models import (
    CandidateVariant, DiseasePhenotype, PhenotypeObservation, ReportSnippet,
)


def _obs(hpo_id, label):
    return PhenotypeObservation(hpo_id=hpo_id, label=label,
                               passage=ReportSnippet(text=label, source="clinical notes"))


def _variant(gene="SCN1A"):
    return CandidateVariant(gene=gene, hgvs_c="c.1A>T", zygosity="heterozygous",
                            reported_significance="Pathogenic",
                            snippet=ReportSnippet(text="x", source="laboratory report"))


def _disease(terms, names=("Dravet syndrome",)):
    return DiseasePhenotype(genes_disease_names=list(names), hpo_terms=terms, source=None)


def _clinvar(sig="Pathogenic", stars=3):
    class _C:
        found = True
        aggregate_significance = sig
        submissions = [type("S", (), {"star_rating": stars})()]
        has_conflict = False
    class _B:
        clinvar = _C()
        gnomad = type("G", (), {"found": True, "global_af": 1e-6})()
    return _B()


def test_supporting_and_missing_partition():
    patient = [_obs("HP:0001250", "Seizure"), _obs("HP:0001263", "Global developmental delay")]
    disease = _disease([("HP:0001250", "Seizure"), ("HP:0002133", "Status epilepticus")])
    a = match.assess_variant(_variant(), disease, _clinvar(), patient)
    assert [o.hpo_id for o in a.supporting] == ["HP:0001250"]
    assert ("HP:0002133", "Status epilepticus") in a.missing


def test_strong_molecular_plus_overlap_is_known():
    patient = [_obs("HP:0001250", "Seizure"), _obs("HP:0001263", "Delay")]
    disease = _disease([("HP:0001250", "Seizure"), ("HP:0001263", "Delay")])
    a = match.assess_variant(_variant(), disease, _clinvar("Pathogenic"), patient)
    assert a.confidence_tier == "Known"
    # Rubric exposes the inputs that produced the tier.
    assert a.rubric["n_supporting"] == 2
    assert a.rubric["clinvar_significance"] == "Pathogenic"


def test_no_overlap_is_speculative_not_known_even_if_pathogenic():
    patient = [_obs("HP:0004322", "Short stature")]
    disease = _disease([("HP:0001250", "Seizure")])
    a = match.assess_variant(_variant(), disease, _clinvar("Pathogenic"), patient)
    assert a.confidence_tier in ("Speculative", "Possible")
    assert a.rubric["n_supporting"] == 0


def test_out_of_scope_variant_is_marked_not_dropped():
    patient = [_obs("HP:0001250", "Seizure")]
    a = match.assess_variant(_variant("NOTAGENE"), None, None, patient, in_scope=False)
    assert a.in_scope is False
    assert a.confidence_tier == "Unknown"
    assert any("scope" in n.lower() for n in a.notes)


def test_unexplained_cluster_surfaces_leftover_patient_features():
    patient = [_obs("HP:0001250", "Seizure"), _obs("HP:0032045", "Ectopia lentis")]
    disease = _disease([("HP:0001250", "Seizure")])
    a = match.assess_variant(_variant(), disease, _clinvar(), patient)
    clusters = match.find_unexplained(patient, [a])
    leftover = {o.hpo_id for c in clusters for o in c.observations}
    assert "HP:0032045" in leftover
    assert "HP:0001250" not in leftover
    assert clusters and clusters[0].recommendation  # non-empty recommendation


def test_rank_orders_best_explanation_first():
    patient = [_obs("HP:0001250", "Seizure"), _obs("HP:0001263", "Delay")]
    strong = match.assess_variant(
        _variant("SCN1A"), _disease([("HP:0001250", "Seizure"), ("HP:0001263", "Delay")]),
        _clinvar("Pathogenic"), patient)
    weak = match.assess_variant(
        _variant("FBN1"), _disease([("HP:0004322", "Short stature")]),
        _clinvar("Uncertain significance", stars=1), patient)
    ranked = match.rank([weak, strong])
    assert ranked[0].variant.gene == "SCN1A"
```

- [ ] **Step 2: Run to confirm failure**

Run: `./.venv/bin/python -m pytest tests/test_match.py -v`
Expected: FAIL — `ModuleNotFoundError: causality_review.match`.

- [ ] **Step 3: Write `causality_review/match.py`**

```python
"""Deterministic phenotype matching + confidence rubric.

No LLM, no network. This is the reverse-matching core: given a patient's HPO
terms and a gene's known HPO spectrum, compute what supports, what conflicts,
what's unexplained, and a confidence tier from a fixed, inspectable rubric. The
rubric dict is returned so the UI can always show *why* a score is what it is.
"""

from __future__ import annotations

from typing import Optional

from causality_review.models import (
    CandidateVariant, DiseasePhenotype, PhenotypeObservation, UnexplainedCluster,
    VariantAssessment,
)

_PATHOGENIC = ("pathogenic", "likely pathogenic")


def _molecular_summary(molecular) -> dict:
    """Flatten the reused EvidenceBundle into rubric-ready primitives."""
    sig, stars, rare = "unknown", 0, None
    if molecular is not None:
        cv = getattr(molecular, "clinvar", None)
        if cv is not None and getattr(cv, "found", False):
            sig = getattr(cv, "aggregate_significance", "") or "unknown"
            subs = getattr(cv, "submissions", []) or []
            stars = max((getattr(s, "star_rating", 0) for s in subs), default=0)
        gn = getattr(molecular, "gnomad", None)
        if gn is not None and getattr(gn, "found", False):
            af = getattr(gn, "global_af", None)
            rare = (af is None) or (af < 1e-4)
        elif gn is not None:
            rare = True  # absent from gnomAD == rare
    return {"clinvar_significance": sig, "clinvar_stars": stars, "gnomad_rare": rare}


def _tier(n_supporting: int, patient_frac: float, mol: dict) -> str:
    sig = (mol["clinvar_significance"] or "").lower()
    pathogenic = any(p in sig for p in _PATHOGENIC)
    if pathogenic and n_supporting >= 2 and patient_frac >= 0.4:
        return "Known"
    if (pathogenic or n_supporting >= 2) and n_supporting >= 1:
        return "Likely"
    if n_supporting >= 1:
        return "Possible"
    if pathogenic or mol["clinvar_significance"] != "unknown":
        return "Speculative"
    return "Unknown"


_TIER_RANK = {"Known": 4, "Likely": 3, "Possible": 2, "Speculative": 1, "Unknown": 0}


def assess_variant(
    variant: CandidateVariant,
    disease: Optional[DiseasePhenotype],
    molecular,
    patient: list[PhenotypeObservation],
    *,
    in_scope: bool = True,
) -> VariantAssessment:
    notes: list[str] = []
    if not in_scope:
        notes.append("Gene is outside the supported allowlist — out of scope; phenotype match skipped.")
        return VariantAssessment(
            variant=variant, in_scope=False, disease=disease, molecular=molecular,
            supporting=[], missing=[], confidence_tier="Unknown",
            rubric={"n_supporting": 0, "patient_explained_fraction": 0.0,
                    "clinvar_significance": "unknown", "clinvar_stars": 0, "gnomad_rare": None,
                    "phenotype_available": False},
            rank_score=-1.0, notes=notes,
        )

    disease_ids = {t[0] for t in (disease.hpo_terms if disease and disease.available else [])}
    supporting = [o for o in patient if o.hpo_id and o.hpo_id in disease_ids]
    supporting_ids = {o.hpo_id for o in supporting}
    missing = [(hid, label) for (hid, label) in (disease.hpo_terms if disease else [])
               if hid not in supporting_ids]

    n_supporting = len(supporting)
    n_patient = len([o for o in patient if o.hpo_id]) or 1
    patient_frac = n_supporting / n_patient
    mol = _molecular_summary(molecular)
    pheno_available = bool(disease and disease.available)
    if not pheno_available:
        notes.append("Gene phenotype source unavailable — phenotype match could not be computed.")

    tier = _tier(n_supporting, patient_frac, mol) if pheno_available else "Unknown"

    rubric = {
        "n_supporting": n_supporting,
        "patient_explained_fraction": round(patient_frac, 3),
        "clinvar_significance": mol["clinvar_significance"],
        "clinvar_stars": mol["clinvar_stars"],
        "gnomad_rare": mol["gnomad_rare"],
        "phenotype_available": pheno_available,
    }
    # Rank primarily by tier, then by how much of the patient it explains, then stars.
    rank_score = _TIER_RANK[tier] * 100 + patient_frac * 10 + mol["clinvar_stars"]

    return VariantAssessment(
        variant=variant, in_scope=True, disease=disease, molecular=molecular,
        supporting=supporting, missing=missing, confidence_tier=tier,
        rubric=rubric, rank_score=rank_score, notes=notes,
    )


def find_unexplained(
    patient: list[PhenotypeObservation],
    assessments: list[VariantAssessment],
) -> list[UnexplainedCluster]:
    """Patient features not in ANY assessed gene's spectrum. One flat cluster for
    the demo (grouping by organ system is out of scope for the hours-slice)."""
    explained: set[str] = set()
    for a in assessments:
        for o in a.supporting:
            explained.add(o.hpo_id)
    leftover = [o for o in patient if o.hpo_id and o.hpo_id not in explained]
    if not leftover:
        return []
    rec = (
        "These documented features are not explained by any reported variant. "
        "Consider whether the reported disease spectrum is broader than currently "
        "annotated, or whether a second, independent molecular diagnosis is present — "
        "genome/exome re-analysis may be warranted."
    )
    return [UnexplainedCluster(observations=leftover, recommendation=rec)]


def rank(assessments: list[VariantAssessment]) -> list[VariantAssessment]:
    return sorted(assessments, key=lambda a: a.rank_score, reverse=True)
```

- [ ] **Step 4: Run the tests**

Run: `./.venv/bin/python -m pytest tests/test_match.py -v`
Expected: PASS (all six tests).

- [ ] **Step 5: Commit**

```bash
git add causality_review/match.py tests/test_match.py
git commit -m "feat(causality): deterministic matching + confidence rubric (core, tested)"
```

---

### Task 5: Grounded extraction (report→variants, notes→phenotype)

Two Claude calls, each constrained to return only what's in the source text, with the verbatim snippet retained for every item.

**Files:**
- Create: `causality_review/extract.py`
- Test: `tests/test_extract.py`

**Interfaces:**
- Consumes: `causality_review.llm.call_json`; `CandidateVariant`, `PhenotypeObservation`, `ReportSnippet` from models.
- Produces: `causality_review.extract.extract_variants(report_text: str) -> list[CandidateVariant]`.
- Produces: `causality_review.extract.extract_phenotype(notes_text: str) -> list[PhenotypeObservation]`.
- Produces: `causality_review.extract.VARIANTS_SCHEMA`, `causality_review.extract.PHENOTYPE_SCHEMA` (module-level dicts, so tests can assert schema validity without the network).

- [ ] **Step 1: Write the failing test**

`tests/test_extract.py` (schema-shape tests run offline; the extraction call is smoke-tested only with credentials):

```python
import pytest

from causality_review import extract, llm


def _valid_schema(s):
    # Structured-output rules: objects declare additionalProperties False and
    # list every property in required.
    if s.get("type") == "object":
        assert s.get("additionalProperties") is False
        assert set(s.get("required", [])) == set(s.get("properties", {}))
        for v in s["properties"].values():
            _valid_schema(v)
    if s.get("type") == "array":
        _valid_schema(s["items"])


def test_variants_schema_is_structured_output_safe():
    _valid_schema(extract.VARIANTS_SCHEMA)


def test_phenotype_schema_is_structured_output_safe():
    _valid_schema(extract.PHENOTYPE_SCHEMA)


@pytest.mark.skipif(not llm.credentials_available(), reason="no Anthropic credentials")
def test_extract_variants_finds_reported_gene():
    report = ("Reportable findings: SCN1A c.4933C>T (p.Arg1645Ter), heterozygous, "
              "classified Pathogenic. FBN1 c.1A>G, heterozygous, Uncertain significance.")
    variants = extract.extract_variants(report)
    genes = {v.gene.upper() for v in variants}
    assert "SCN1A" in genes
    assert all(v.snippet.text for v in variants)  # every variant keeps a snippet
```

- [ ] **Step 2: Run to confirm failure**

Run: `./.venv/bin/python -m pytest tests/test_extract.py -v`
Expected: FAIL — `ModuleNotFoundError: causality_review.extract`.

- [ ] **Step 3: Write `causality_review/extract.py`**

```python
"""Grounded extraction with Claude.

Two calls: the report -> candidate variants, the notes -> patient phenotype (HPO
terms). Each returned item keeps the verbatim source passage it was drawn from,
so nothing downstream is asserted without provenance. The model is instructed to
extract only what is literally present and to leave a field empty rather than
guess.
"""

from __future__ import annotations

from causality_review import llm
from causality_review.models import CandidateVariant, PhenotypeObservation, ReportSnippet

VARIANTS_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "properties": {
        "variants": {
            "type": "array",
            "items": {
                "type": "object",
                "additionalProperties": False,
                "properties": {
                    "gene": {"type": "string"},
                    "hgvs_c": {"type": "string"},
                    "zygosity": {"type": "string"},
                    "reported_significance": {"type": "string"},
                    "snippet": {"type": "string"},
                },
                "required": ["gene", "hgvs_c", "zygosity", "reported_significance", "snippet"],
            },
        }
    },
    "required": ["variants"],
}

PHENOTYPE_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "properties": {
        "observations": {
            "type": "array",
            "items": {
                "type": "object",
                "additionalProperties": False,
                "properties": {
                    "hpo_id": {"type": "string"},
                    "label": {"type": "string"},
                    "passage": {"type": "string"},
                },
                "required": ["hpo_id", "label", "passage"],
            },
        }
    },
    "required": ["observations"],
}

_VARIANTS_SYSTEM = (
    "You extract reported genetic variants from a laboratory report. Return ONLY "
    "variants literally present in the text. For each: the gene symbol, the HGVS "
    "coding change (c. notation), zygosity if stated (else empty string), the "
    "lab's reported classification if stated (else empty string), and 'snippet' = "
    "the verbatim sentence or phrase you drew it from. Never invent a variant or a "
    "classification not in the text."
)

_PHENOTYPE_SYSTEM = (
    "You extract the patient's phenotype from clinical notes as HPO terms. For each "
    "distinct clinical feature actually documented: 'hpo_id' = the best-matching "
    "HPO id (format 'HP:0000000'; use an empty string if you are not confident of "
    "the exact id), 'label' = the standard HPO term name, and 'passage' = the "
    "verbatim phrase from the notes that documents it. Extract only documented "
    "findings — not negated findings, not family history, not possibilities. Never "
    "invent a feature that is not in the notes."
)


def extract_variants(report_text: str) -> list[CandidateVariant]:
    data = llm.call_json(_VARIANTS_SYSTEM, report_text, VARIANTS_SCHEMA, effort="medium")
    out: list[CandidateVariant] = []
    for v in data.get("variants", []):
        out.append(CandidateVariant(
            gene=(v.get("gene") or "").strip(),
            hgvs_c=(v.get("hgvs_c") or "").strip(),
            zygosity=(v.get("zygosity") or "").strip() or None,
            reported_significance=(v.get("reported_significance") or "").strip() or None,
            snippet=ReportSnippet(text=(v.get("snippet") or "").strip(), source="laboratory report"),
        ))
    return out


def extract_phenotype(notes_text: str) -> list[PhenotypeObservation]:
    data = llm.call_json(_PHENOTYPE_SYSTEM, notes_text, PHENOTYPE_SCHEMA, effort="medium")
    out: list[PhenotypeObservation] = []
    for o in data.get("observations", []):
        out.append(PhenotypeObservation(
            hpo_id=(o.get("hpo_id") or "").strip(),
            label=(o.get("label") or "").strip(),
            passage=ReportSnippet(text=(o.get("passage") or "").strip(), source="clinical notes"),
        ))
    return out
```

- [ ] **Step 4: Run the tests**

Run: `./.venv/bin/python -m pytest tests/test_extract.py -v`
Expected: schema tests PASS; the extraction test PASSES with credentials or SKIPS without.

- [ ] **Step 5: Commit**

```bash
git add causality_review/extract.py tests/test_extract.py
git commit -m "feat(causality): grounded variant + phenotype extraction"
```

---

### Task 6: Grounded narrative + grounding-check evaluator

Two more Claude calls, both constrained to the retrieved evidence. The narrator writes the per-variant story with Known/Likely/Possible/Speculative/Unknown tags; the evaluator (fresh context) flags any statement not traceable to the evidence.

**Files:**
- Create: `causality_review/narrate.py`
- Create: `causality_review/evaluate.py`
- Test: `tests/test_narrate.py`

**Interfaces:**
- Consumes: `causality_review.llm.call_json`; `VariantAssessment`, `PhenotypeObservation`.
- Produces: `causality_review.narrate.write_narrative(assessment: VariantAssessment) -> str`.
- Produces: `causality_review.evaluate.grounding_check(assessment: VariantAssessment) -> list[str]` (returns a list of flag strings; empty means clean).
- Produces: `causality_review.narrate.NARRATIVE_SCHEMA`, `causality_review.evaluate.EVAL_SCHEMA`.

- [ ] **Step 1: Write the failing test**

`tests/test_narrate.py`:

```python
import pytest

from causality_review import evaluate, llm, narrate
from causality_review.match import assess_variant
from causality_review.models import (
    CandidateVariant, DiseasePhenotype, PhenotypeObservation, ReportSnippet,
)


def _fixture_assessment():
    patient = [PhenotypeObservation("HP:0001250", "Seizure",
                                    ReportSnippet("recurrent seizures since infancy", "clinical notes"))]
    disease = DiseasePhenotype(["Dravet syndrome"], [("HP:0001250", "Seizure")], source=None)
    variant = CandidateVariant("SCN1A", "c.4933C>T", "heterozygous", "Pathogenic",
                               ReportSnippet("SCN1A c.4933C>T Pathogenic", "laboratory report"))
    class _C:
        found = True; aggregate_significance = "Pathogenic"; has_conflict = False
        submissions = [type("S", (), {"star_rating": 3})()]
    class _B:
        clinvar = _C(); gnomad = type("G", (), {"found": True, "global_af": 1e-6})()
    return assess_variant(variant, disease, _B(), patient)


def test_narrative_schema_present():
    assert narrate.NARRATIVE_SCHEMA["type"] == "object"
    assert evaluate.EVAL_SCHEMA["type"] == "object"


@pytest.mark.skipif(not llm.credentials_available(), reason="no Anthropic credentials")
def test_write_narrative_mentions_supporting_feature():
    text = narrate.write_narrative(_fixture_assessment())
    assert isinstance(text, str) and len(text) > 0


@pytest.mark.skipif(not llm.credentials_available(), reason="no Anthropic credentials")
def test_grounding_check_returns_list():
    flags = evaluate.grounding_check(_fixture_assessment())
    assert isinstance(flags, list)
```

- [ ] **Step 2: Run to confirm failure**

Run: `./.venv/bin/python -m pytest tests/test_narrate.py -v`
Expected: FAIL — `ModuleNotFoundError`.

- [ ] **Step 3: Write `causality_review/narrate.py`**

```python
"""Grounded per-variant narrative.

The narrator is given ONLY the retrieved evidence for one variant (its supporting
and missing terms, molecular summary, and confidence tier) and must write a
causality assessment constrained to that evidence, tagging each statement with one
of Known / Likely / Possible / Speculative / Unknown. It is told to abstain rather
than assert beyond the evidence.
"""

from __future__ import annotations

import json

from causality_review import llm
from causality_review.match import _molecular_summary
from causality_review.models import VariantAssessment

NARRATIVE_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "properties": {"narrative": {"type": "string"}},
    "required": ["narrative"],
}

_SYSTEM = (
    "You are assisting a clinical geneticist with post-laboratory causality review. "
    "Write a concise assessment of whether ONE reported variant explains THIS "
    "patient, using ONLY the evidence provided in the user message. Address: which "
    "documented patient features the gene's phenotype supports, which expected "
    "features are missing, and what the molecular evidence adds. Tag each statement "
    "with its certainty in brackets — one of [Known], [Likely], [Possible], "
    "[Speculative], [Unknown]. Do not introduce any fact, gene function, or "
    "phenotype not present in the evidence. If the evidence is thin, say so plainly "
    "rather than overstating. Do not give directive medical advice."
)


def _evidence_payload(a: VariantAssessment) -> str:
    return json.dumps({
        "variant": {"gene": a.variant.gene, "hgvs_c": a.variant.hgvs_c,
                    "zygosity": a.variant.zygosity,
                    "reported_significance": a.variant.reported_significance},
        "gene_diseases": a.disease.genes_disease_names if a.disease else [],
        "supporting_patient_features": [{"hpo_id": o.hpo_id, "label": o.label,
                                         "note": o.passage.text} for o in a.supporting],
        "missing_expected_features": [{"hpo_id": h, "label": l} for (h, l) in a.missing],
        "molecular": _molecular_summary(a.molecular),
        "confidence_tier": a.confidence_tier,
        "rubric_inputs": a.rubric,
    }, indent=2)


def write_narrative(a: VariantAssessment) -> str:
    if not a.in_scope:
        return ("This gene is outside the tool's supported allowlist, so no phenotype "
                "matching was performed. [Unknown]")
    data = llm.call_json(_SYSTEM, _evidence_payload(a), NARRATIVE_SCHEMA, effort="high", max_tokens=2000)
    return (data.get("narrative") or "").strip()
```

- [ ] **Step 4: Write `causality_review/evaluate.py`**

```python
"""Grounding-check evaluator — the vision's independent reviewer, scoped to one
pass. Given the narrative AND the exact evidence it was allowed to cite, it flags
any statement not traceable to that evidence. Flags surface in the UI; they are
never silently dropped.
"""

from __future__ import annotations

import json

from causality_review import llm
from causality_review.match import _molecular_summary
from causality_review.models import VariantAssessment

EVAL_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "properties": {
        "flags": {
            "type": "array",
            "items": {
                "type": "object",
                "additionalProperties": False,
                "properties": {"statement": {"type": "string"}, "issue": {"type": "string"}},
                "required": ["statement", "issue"],
            },
        }
    },
    "required": ["flags"],
}

_SYSTEM = (
    "You are an independent grounding checker. You are given a causality narrative "
    "and the ONLY evidence its author was permitted to use. Identify statements in "
    "the narrative that are NOT traceable to that evidence — fabricated phenotype "
    "links, gene facts not provided, or confidence not justified by the rubric. "
    "Return one flag per problem with the offending statement and a short reason. "
    "If everything is traceable, return an empty list. Do not flag correctly-hedged "
    "or correctly-cited statements."
)


def grounding_check(a: VariantAssessment) -> list[str]:
    if not a.in_scope or not a.narrative:
        return []
    payload = json.dumps({
        "narrative": a.narrative,
        "allowed_evidence": {
            "gene_diseases": a.disease.genes_disease_names if a.disease else [],
            "supporting": [o.label for o in a.supporting],
            "missing": [l for (_h, l) in a.missing],
            "molecular": _molecular_summary(a.molecular),
            "rubric": a.rubric,
        },
    }, indent=2)
    data = llm.call_json(_SYSTEM, payload, EVAL_SCHEMA, effort="high", max_tokens=1500)
    return [f"{f.get('statement','').strip()} — {f.get('issue','').strip()}"
            for f in data.get("flags", [])]
```

- [ ] **Step 5: Run the tests**

Run: `./.venv/bin/python -m pytest tests/test_narrate.py -v`
Expected: schema test PASSES; the two live tests PASS with credentials or SKIP without.

- [ ] **Step 6: Commit**

```bash
git add causality_review/narrate.py causality_review/evaluate.py tests/test_narrate.py
git commit -m "feat(causality): grounded narrative + grounding-check evaluator"
```

---

### Task 7: Pipeline orchestration

Wire the pieces: extract → retrieve (HPO + reuse variant_curator molecular) → match → narrate → evaluate → rank → assemble.

**Files:**
- Create: `causality_review/pipeline.py`
- Test: `tests/test_pipeline.py`

**Interfaces:**
- Consumes: everything above; `variant_curator.pipeline.assemble_evidence`, `variant_curator.models.VariantInput`, `variant_curator.genes.is_supported_gene`, `variant_curator.genes.get_gene`, `variant_curator.http.get_client`.
- Produces: `causality_review.pipeline.run_review(report_text: str, notes_text: str) -> CausalityReview`.

- [ ] **Step 1: Write the failing test** (uses monkeypatch to keep it offline)

`tests/test_pipeline.py`:

```python
from causality_review import pipeline
from causality_review.models import CandidateVariant, PhenotypeObservation, ReportSnippet, DiseasePhenotype


def test_run_review_ranks_and_finds_unexplained(monkeypatch):
    # Stub the LLM-backed and network-backed pieces so this runs offline.
    scn1a = CandidateVariant("SCN1A", "c.4933C>T", "heterozygous", "Pathogenic",
                             ReportSnippet("SCN1A c.4933C>T Pathogenic", "laboratory report"))
    patient = [
        PhenotypeObservation("HP:0001250", "Seizure", ReportSnippet("seizures", "clinical notes")),
        PhenotypeObservation("HP:0032045", "Ectopia lentis", ReportSnippet("lens dislocation", "clinical notes")),
    ]
    monkeypatch.setattr(pipeline.extract, "extract_variants", lambda t: [scn1a])
    monkeypatch.setattr(pipeline.extract, "extract_phenotype", lambda t: patient)
    monkeypatch.setattr(pipeline, "_molecular_for", lambda client, v: None)
    monkeypatch.setattr(pipeline.hpo, "fetch_gene_phenotype",
                        lambda client, gid: DiseasePhenotype(["Dravet syndrome"],
                                                             [("HP:0001250", "Seizure")], source=None))
    monkeypatch.setattr(pipeline.narrate, "write_narrative", lambda a: "stub narrative [Likely]")
    monkeypatch.setattr(pipeline.evaluate, "grounding_check", lambda a: [])

    review = pipeline.run_review("report", "notes")
    assert review.assessments[0].variant.gene == "SCN1A"
    # The connective-tissue feature is unexplained by SCN1A.
    leftover = {o.hpo_id for c in review.unexplained for o in c.observations}
    assert "HP:0032045" in leftover
    d = review.to_dict()
    assert d["assessments"][0]["confidence_tier"] in ("Known", "Likely", "Possible")
```

- [ ] **Step 2: Run to confirm failure**

Run: `./.venv/bin/python -m pytest tests/test_pipeline.py -v`
Expected: FAIL — `ModuleNotFoundError: causality_review.pipeline`.

- [ ] **Step 3: Write `causality_review/pipeline.py`**

```python
"""Orchestrates one causality review end to end.

extract (Claude) -> retrieve (HPO/Jax + reuse variant_curator VEP/gnomAD/ClinVar)
-> match (deterministic) -> narrate (Claude) -> evaluate (Claude) -> rank -> assemble.
Each stage degrades safely: a gene outside the allowlist is marked out-of-scope
(not dropped); an unreachable source is recorded, never fabricated.
"""

from __future__ import annotations

from causality_review import evaluate, extract, hpo, narrate
from causality_review.match import assess_variant, find_unexplained, rank
from causality_review.models import CausalityReview, VariantAssessment
from variant_curator.genes import get_gene, is_supported_gene
from variant_curator.http import get_client
from variant_curator.models import VariantInput


def _molecular_for(client, variant):
    """Reuse the variant_curator molecular pipeline. Returns an EvidenceBundle or
    None on failure (kept separate so tests can stub it)."""
    try:
        from variant_curator.pipeline import assemble_evidence
        return assemble_evidence(VariantInput(gene=variant.gene, hgvs_c=variant.hgvs_c))
    except Exception:  # noqa: BLE001 — molecular axis is best-effort
        return None


def run_review(report_text: str, notes_text: str) -> CausalityReview:
    warnings: list[str] = []
    variants = extract.extract_variants(report_text)
    patient = extract.extract_phenotype(notes_text)
    if not variants:
        warnings.append("No reported variants could be extracted from the laboratory report.")
    if not patient:
        warnings.append("No phenotype could be extracted from the clinical notes.")

    assessments: list[VariantAssessment] = []
    sources = []
    with get_client() as client:
        for v in variants:
            if not v.gene or not is_supported_gene(v.gene):
                assessments.append(assess_variant(v, None, None, patient, in_scope=False))
                continue
            gene = get_gene(v.gene)
            disease = hpo.fetch_gene_phenotype(client, gene.ncbi_gene_id)
            if disease.source is not None:
                sources.append(disease.source)
            molecular = _molecular_for(client, v)
            if molecular is not None:
                sources.extend(molecular.sources())
            assessments.append(assess_variant(v, disease, molecular, patient))

    # Narrative + grounding check per in-scope assessment.
    for a in assessments:
        a.narrative = narrate.write_narrative(a)
        a.grounding_flags = evaluate.grounding_check(a)

    ranked = rank(assessments)
    unexplained = find_unexplained(patient, ranked)
    return CausalityReview(assessments=ranked, unexplained=unexplained,
                           sources=sources, warnings=warnings)
```

- [ ] **Step 4: Run the tests**

Run: `./.venv/bin/python -m pytest tests/test_pipeline.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add causality_review/pipeline.py tests/test_pipeline.py
git commit -m "feat(causality): pipeline orchestration (extract->retrieve->match->narrate->eval->rank)"
```

---

### Task 8: FastAPI server + three-pane UI + example case

The demo surface. `POST /review` runs the pipeline; `GET /` serves one page; `GET /example` returns the baked case text.

**Files:**
- Create: `causality_review/server.py`
- Create: `causality_review/web/index.html`
- Create: `causality_review/examples/report.txt`
- Create: `causality_review/examples/notes.txt`
- Test: `tests/test_server.py`

**Interfaces:**
- Consumes: `causality_review.pipeline.run_review`, `causality_review.llm` (credential check + error types).
- Produces: FastAPI `app`; routes `GET /`, `GET /example`, `POST /review`.

- [ ] **Step 1: Write the example case files**

`causality_review/examples/report.txt`:

```
GENOMIC SEQUENCING REPORT — REPORTABLE FINDINGS (synthetic, no PHI)

Patient: pediatric proband. Test: trio exome sequencing.

Finding 1: SCN1A (NM_001165963.4) c.4933C>T (p.Arg1645Ter), heterozygous,
de novo. Classification: Pathogenic. Associated with Dravet syndrome /
developmental and epileptic encephalopathy.

Finding 2: FBN1 (NM_000138.5) c.4082G>A (p.Cys1361Tyr), heterozygous,
maternally inherited. Classification: Uncertain significance. FBN1 is
associated with Marfan syndrome.

Finding 3: MYBPC3 (NM_000256.3) c.1090G>A (p.Ala364Thr), heterozygous.
Classification: Uncertain significance.
```

`causality_review/examples/notes.txt`:

```
CLINICAL SUMMARY (synthetic, no PHI)

3-year-old with onset of prolonged febrile and afebrile seizures beginning at
6 months of age, including episodes of status epilepticus. Global developmental
delay with regression noted after seizure onset. EEG with multifocal epileptiform
discharges.

On examination the child is notably tall for age with long, slender fingers
(arachnodactyly) and hyperextensible joints. Ophthalmology documented bilateral
lens dislocation (ectopia lentis). Echocardiogram: mild aortic root dilatation.

No family history of sudden cardiac death. Cardiac exam otherwise unremarkable.
```

> This case is the whole demo: SCN1A cleanly explains the epilepsy cluster; the connective-tissue cluster (tall stature, arachnodactyly, ectopia lentis, aortic root dilatation) is left unexplained by SCN1A and points at FBN1/Marfan — exercising overfitting avoidance + the dual-diagnosis surface.

- [ ] **Step 2: Write the failing test**

`tests/test_server.py`:

```python
from fastapi.testclient import TestClient

from causality_review.server import app

client = TestClient(app)


def test_index_serves_html():
    r = client.get("/")
    assert r.status_code == 200
    assert "causality" in r.text.lower()


def test_example_returns_report_and_notes():
    r = client.get("/example")
    assert r.status_code == 200
    body = r.json()
    assert "SCN1A" in body["report"]
    assert "seizure" in body["notes"].lower()


def test_review_without_credentials_returns_503(monkeypatch):
    import causality_review.server as srv
    monkeypatch.setattr(srv.llm, "credentials_available", lambda: False)
    r = client.post("/review", json={"report": "x", "notes": "y"})
    assert r.status_code == 503
```

- [ ] **Step 3: Run to confirm failure**

Run: `./.venv/bin/python -m pytest tests/test_server.py -v`
Expected: FAIL — `ModuleNotFoundError: causality_review.server`.

- [ ] **Step 4: Write `causality_review/server.py`**

```python
"""FastAPI app: serves the three-pane UI and runs the causality pipeline.

The static page and the /example case load with no credentials; only /review
needs Claude, and it returns a clear 503 when credentials are absent.
"""

from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse
from pydantic import BaseModel

from causality_review import llm
from causality_review.pipeline import run_review

app = FastAPI(title="Post-Laboratory Causality Review")

_WEB = Path(__file__).parent / "web"
_EXAMPLES = Path(__file__).parent / "examples"


class ReviewRequest(BaseModel):
    report: str
    notes: str


@app.get("/", response_class=HTMLResponse)
def index() -> str:
    return (_WEB / "index.html").read_text()


@app.get("/example")
def example() -> dict:
    return {
        "report": (_EXAMPLES / "report.txt").read_text(),
        "notes": (_EXAMPLES / "notes.txt").read_text(),
    }


@app.post("/review")
def review(req: ReviewRequest) -> dict:
    if not llm.credentials_available():
        raise HTTPException(status_code=503,
                            detail="No Anthropic credentials configured. Set ANTHROPIC_API_KEY.")
    if not req.report.strip() or not req.notes.strip():
        raise HTTPException(status_code=400, detail="Both a report and clinical notes are required.")
    try:
        return run_review(req.report, req.notes).to_dict()
    except llm.LLMUnavailable as e:
        raise HTTPException(status_code=503, detail=str(e)) from e
    except llm.LLMRefused as e:
        raise HTTPException(status_code=422,
                            detail=f"The model declined this request: {e}") from e
```

- [ ] **Step 5: Write `causality_review/web/index.html`**

A single self-contained page. Left = report, middle = notes, right = ranked review; each variant card expands to its evidence page; a banner surfaces unexplained clusters.

```html
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Post-Laboratory Causality Review</title>
  <style>
    :root { --bg:#0f1417; --panel:#161d22; --line:#26313a; --ink:#e7edf1; --muted:#94a6b2;
            --known:#3fb27f; --likely:#5aa9e6; --possible:#d8a24a; --spec:#c77dcb; --unknown:#7a8a95;
            --flag:#e0685f; }
    * { box-sizing:border-box; }
    body { margin:0; font:14px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;
           background:var(--bg); color:var(--ink); }
    header { padding:14px 20px; border-bottom:1px solid var(--line); display:flex; gap:16px; align-items:center; }
    header h1 { font-size:16px; margin:0; font-weight:600; }
    header button { margin-left:auto; background:#22303a; color:var(--ink); border:1px solid var(--line);
                    padding:8px 14px; border-radius:8px; cursor:pointer; font-size:13px; }
    header button.primary { background:#2d6cdf; border-color:#2d6cdf; }
    header button:disabled { opacity:.5; cursor:default; }
    main { display:grid; grid-template-columns:1fr 1fr 1.3fr; gap:1px; background:var(--line);
            height:calc(100vh - 57px); }
    .pane { background:var(--panel); overflow:auto; padding:16px; }
    .pane h2 { font-size:12px; text-transform:uppercase; letter-spacing:.06em; color:var(--muted); margin:0 0 10px; }
    textarea { width:100%; height:calc(100% - 30px); background:#0e1418; color:var(--ink);
               border:1px solid var(--line); border-radius:8px; padding:10px; resize:none; font:13px/1.5 ui-monospace,monospace; }
    .banner { background:#3a2320; border:1px solid #6e3c34; color:#f2c7bf; padding:12px 14px;
              border-radius:10px; margin-bottom:14px; }
    .banner b { color:#ffd9d1; }
    .card { border:1px solid var(--line); border-radius:12px; margin-bottom:12px; overflow:hidden; }
    .card > summary { list-style:none; cursor:pointer; padding:12px 14px; display:flex; gap:10px; align-items:center; }
    .card > summary::-webkit-details-marker { display:none; }
    .tier { font-size:11px; font-weight:700; padding:3px 9px; border-radius:20px; color:#0c1013; }
    .tier.Known{background:var(--known)} .tier.Likely{background:var(--likely)}
    .tier.Possible{background:var(--possible)} .tier.Speculative{background:var(--spec)}
    .tier.Unknown{background:var(--unknown)}
    .gene { font-weight:600; } .hgvs { color:var(--muted); font-family:ui-monospace,monospace; font-size:12px; }
    .body { padding:0 14px 14px; border-top:1px solid var(--line); }
    .sec { margin-top:12px; }
    .sec h3 { font-size:11px; text-transform:uppercase; letter-spacing:.05em; color:var(--muted); margin:0 0 6px; }
    .chip { display:inline-block; background:#1d262c; border:1px solid var(--line); border-radius:6px;
            padding:2px 8px; margin:0 6px 6px 0; font-size:12px; }
    .chip.miss { opacity:.65; }
    .quote { border-left:3px solid var(--line); padding:2px 0 2px 10px; color:var(--muted); margin:4px 0; font-size:12px; }
    .narr { white-space:pre-wrap; background:#0e1418; border:1px solid var(--line); border-radius:8px; padding:10px; }
    .flag { color:var(--flag); font-size:12px; margin:4px 0; }
    .rubric { font-size:12px; color:var(--muted); }
    .rubric code { color:var(--ink); }
    a { color:#6db3ff; } .muted { color:var(--muted); }
    .spin { color:var(--muted); }
  </style>
</head>
<body>
  <header>
    <h1>Post-Laboratory Causality Review</h1>
    <span class="muted" id="status"></span>
    <button id="load">Load example case</button>
    <button id="run" class="primary">Run causality review</button>
  </header>
  <main>
    <section class="pane">
      <h2>Laboratory report</h2>
      <textarea id="report" placeholder="Paste the genetic testing report (reportable variants)…"></textarea>
    </section>
    <section class="pane">
      <h2>Clinical notes</h2>
      <textarea id="notes" placeholder="Paste clinical notes, H&amp;P, consults, imaging…"></textarea>
    </section>
    <section class="pane" id="out">
      <h2>Ranked causality review</h2>
      <p class="muted">Load the example case (or paste your own), then run the review.</p>
    </section>
  </main>
  <script>
    const $ = (id) => document.getElementById(id);
    const esc = (s) => (s||"").replace(/[&<>]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;'}[c]));

    $("load").onclick = async () => {
      const r = await fetch("/example"); const d = await r.json();
      $("report").value = d.report; $("notes").value = d.notes;
    };

    $("run").onclick = async () => {
      $("run").disabled = true; $("status").textContent = "Analyzing… (this can take a minute)";
      $("out").innerHTML = '<h2>Ranked causality review</h2><p class="spin">Extracting variants and phenotype, retrieving HPO / ClinVar / gnomAD, matching…</p>';
      try {
        const r = await fetch("/review", {method:"POST", headers:{"Content-Type":"application/json"},
          body: JSON.stringify({report:$("report").value, notes:$("notes").value})});
        if (!r.ok) { const e = await r.json(); throw new Error(e.detail || r.statusText); }
        render(await r.json());
      } catch (e) {
        $("out").innerHTML = '<h2>Ranked causality review</h2><p class="flag">'+esc(e.message)+'</p>';
      } finally { $("run").disabled = false; $("status").textContent = ""; }
    };

    function render(review) {
      let h = '<h2>Ranked causality review</h2>';
      (review.warnings||[]).forEach(w => h += '<div class="banner">'+esc(w)+'</div>');
      const unexplained = review.unexplained||[];
      if (unexplained.length) {
        unexplained.forEach(c => {
          const feats = c.observations.map(o => esc(o.label||o.hpo_id)).join(", ");
          h += '<div class="banner"><b>Unexplained phenotype cluster.</b> '
             + esc(c.recommendation) + '<br><span class="muted">Features: '+feats+'</span></div>';
        });
      }
      (review.assessments||[]).forEach((a,i) => h += card(a, i===0));
      h += sources(review.sources||[]);
      $("out").innerHTML = h;
    }

    function card(a, open) {
      const v = a.variant;
      let h = '<details class="card"'+(open?' open':'')+'><summary>'
        + '<span class="tier '+esc(a.confidence_tier)+'">'+esc(a.confidence_tier)+'</span>'
        + '<span class="gene">'+esc(v.gene)+'</span> <span class="hgvs">'+esc(v.hgvs_c)+'</span>'
        + (v.reported_significance?'<span class="muted"> · '+esc(v.reported_significance)+'</span>':'')
        + '</summary><div class="body">';
      if (a.in_scope === false) {
        h += '<p class="muted">'+esc((a.notes||[]).join(" "))+'</p></div></details>';
        return h;
      }
      if (a.narrative) h += '<div class="sec"><h3>Assessment</h3><div class="narr">'+esc(a.narrative)+'</div></div>';
      h += '<div class="sec"><h3>Supporting patient features</h3>'
         + (a.supporting.length ? a.supporting.map(o =>
             '<div class="chip">'+esc(o.label||o.hpo_id)+'</div>'
           + '<div class="quote">'+esc(o.passage.text)+'</div>').join("")
           : '<span class="muted">none</span>') + '</div>';
      h += '<div class="sec"><h3>Expected but not documented</h3>'
         + (a.missing.length ? a.missing.slice(0,12).map(m => '<span class="chip miss">'+esc(m[1]||m[0])+'</span>').join("")
           : '<span class="muted">none</span>')
         + (a.missing.length>12?'<span class="muted"> +'+(a.missing.length-12)+' more</span>':'') + '</div>';
      const R = a.rubric||{};
      h += '<div class="sec"><h3>Why this tier</h3><div class="rubric">'
         + 'supporting features: <code>'+R.n_supporting+'</code> · '
         + 'patient explained: <code>'+Math.round((R.patient_explained_fraction||0)*100)+'%</code> · '
         + 'ClinVar: <code>'+esc(R.clinvar_significance)+'</code> ('+R.clinvar_stars+'★) · '
         + 'gnomAD rare: <code>'+String(R.gnomad_rare)+'</code></div></div>';
      if ((a.grounding_flags||[]).length) {
        h += '<div class="sec"><h3>Grounding check</h3>'
           + a.grounding_flags.map(f => '<div class="flag">⚑ '+esc(f)+'</div>').join("") + '</div>';
      }
      h += '</div></details>';
      return h;
    }

    function sources(srcs) {
      const seen = {}; const uniq = srcs.filter(s => s && s.url && !seen[s.url] && (seen[s.url]=1));
      if (!uniq.length) return '';
      return '<div class="sec"><h3>Sources</h3>'
        + uniq.map(s => '<div><a href="'+esc(s.url)+'" target="_blank" rel="noopener">'
          + esc(s.name)+'</a> <span class="muted">'+esc(s.detail||"")+'</span></div>').join("") + '</div>';
    }
  </script>
</body>
</html>
```

- [ ] **Step 6: Run the server tests**

Run: `./.venv/bin/python -m pytest tests/test_server.py -v`
Expected: PASS (all three).

- [ ] **Step 7: Commit**

```bash
git add causality_review/server.py causality_review/web/index.html causality_review/examples/report.txt causality_review/examples/notes.txt tests/test_server.py
git commit -m "feat(causality): FastAPI server, three-pane UI, synthetic example case"
```

---

### Task 9: Run script, README, and end-to-end verification

Make it launchable and document it; verify the demo path end to end (live if credentials are present).

**Files:**
- Create: `causality_review/README.md`
- Create: `run_causality.sh`
- Modify: `tests/test_pipeline.py` (add an opt-in live end-to-end test)

**Interfaces:** none new — this task ties off the deliverable.

- [ ] **Step 1: Write `run_causality.sh`**

```bash
#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [ -z "${ANTHROPIC_API_KEY:-}" ] && [ -z "${ANTHROPIC_AUTH_TOKEN:-}" ]; then
  echo "warning: no ANTHROPIC_API_KEY set — the UI and example load, but /review will 503." >&2
fi
exec ./.venv/bin/python -m uvicorn causality_review.server:app --reload --port 8000
```

```bash
chmod +x run_causality.sh
```

- [ ] **Step 2: Add an opt-in live end-to-end test**

Append to `tests/test_pipeline.py`:

```python
import os
import pytest
from causality_review import llm


@pytest.mark.skipif(not (llm.credentials_available() and os.environ.get("RUN_LIVE")),
                    reason="set RUN_LIVE=1 with credentials for the live end-to-end test")
def test_live_example_case_surfaces_scn1a_and_unexplained_cluster():
    from pathlib import Path
    ex = Path("causality_review/examples")
    review = pipeline.run_review((ex / "report.txt").read_text(), (ex / "notes.txt").read_text())
    genes = [a.variant.gene.upper() for a in review.assessments]
    assert "SCN1A" in genes
    # SCN1A should rank at or near the top for the epilepsy cluster.
    assert review.assessments[0].variant.gene.upper() in ("SCN1A", "FBN1")
    # The connective-tissue features should not be fully explained by SCN1A alone.
    assert review.unexplained, "expected an unexplained phenotype cluster for the dual-diagnosis case"
```

- [ ] **Step 3: Write `causality_review/README.md`**

```markdown
# Post-Laboratory Causality Review

A clinician-facing copilot for the *second* stage of genomic interpretation. The
lab has already reported ~2–5 candidate variants; this tool answers the different
question: **given these variants and the full patient story, which one actually
explains this patient?**

It is the reverse of `variant_curator/` (which classifies one variant into an
ACMG tier) and reuses that package's VEP→gnomAD→ClinVar plumbing for the
molecular axis.

## What it does

1. **Extract** (Claude, grounded): reported variants from the report; HPO
   phenotype from the notes — each item keeps its verbatim source snippet.
2. **Retrieve** (cited): each gene's known HPO spectrum from the HPO/Jax ontology
   API; molecular evidence via `variant_curator`.
3. **Match** (deterministic): supporting (patient ∩ disease), missing
   (disease − patient), and unexplained clusters (patient − all genes), with a
   fixed confidence rubric whose inputs are shown in the UI.
4. **Narrate + grounding-check** (Claude): a per-variant assessment tagged
   Known/Likely/Possible/Speculative/Unknown, then an independent pass that flags
   any statement not traceable to the retrieved evidence.
5. **Rank** best-explanation-first and surface unexplained clusters + a
   re-analysis / possible-second-diagnosis recommendation.

## Run

```bash
export ANTHROPIC_API_KEY=...        # or `ant auth login`
./run_causality.sh                  # http://localhost:8000
```

Click **Load example case**, then **Run causality review**. The synthetic case
shows the two failure modes the tool guards against: SCN1A cleanly explains the
epilepsy cluster, while the connective-tissue features (FBN1/Marfan-shaped)
remain unexplained — so the tool refuses to overfit SCN1A and flags a possible
second molecular diagnosis.

## Model

Uses `claude-fable-5` by default (override with `CAUSALITY_MODEL`), with an
automatic Opus 4.8 refusal-fallback since genomics text can trip Fable's safety
classifiers.

## House rules (enforced in code)

Never invent evidence · every claim traces to the report, the notes, or a
retrieved source · distinguish Known/Likely/Possible/Speculative/Unknown · never
hide uncertainty · always show why a tier is what it is. Scope is the 30-gene
allowlist in `variant_curator/genes.py`; other genes are marked out-of-scope, not
guessed. Synthetic data only — no PHI.

## Tests

```bash
./.venv/bin/python -m pytest tests/ -v          # offline: deterministic core + schemas + server
RUN_LIVE=1 ./.venv/bin/python -m pytest tests/ -v   # + live LLM path (needs credentials)
```
```

- [ ] **Step 4: Full offline test run**

Run: `./.venv/bin/python -m pytest tests/ -v`
Expected: all deterministic/schema/server tests PASS; live-only tests SKIP.

- [ ] **Step 5: Live end-to-end verification (if credentials available)**

Run:
```bash
export ANTHROPIC_API_KEY=...   # your Fable 5 key
RUN_LIVE=1 ./.venv/bin/python -m pytest tests/test_pipeline.py -v -k live
```
Then launch and eyeball the demo:
```bash
./run_causality.sh
# open http://localhost:8000, Load example case, Run causality review
```
Expected: SCN1A ranks at/near top with supporting seizure/developmental features;
an unexplained-cluster banner lists the connective-tissue features with a
re-analysis recommendation; each card shows its rubric inputs and sources; no
grounding flags on well-formed narratives.

- [ ] **Step 6: Commit**

```bash
git add causality_review/README.md run_causality.sh tests/test_pipeline.py
git commit -m "feat(causality): run script, README, live end-to-end verification"
```

---

## Self-Review

**Spec coverage:** ingest (Task 8 UI + `/example`), grounded extraction with snippets (Task 5), HPO/Jax retrieval + molecular reuse (Tasks 3, 7), deterministic matching + rubric (Task 4), unexplained clusters + dual-diagnosis recommendation (Task 4 `find_unexplained`, surfaced in Task 8 banner), narrative with certainty tags (Task 6), grounding-check evaluator (Task 6), ranked three-pane output with per-variant evidence pages + sources (Task 8), house rules enforced structurally (Global Constraints + Tasks 4–6), error handling for out-of-scope genes / unreachable HPO / missing key / refusal (Tasks 4, 3, 8, 1), tests for the deterministic core and an opt-in live end-to-end (Tasks 4, 9). All spec sections map to a task.

**Placeholder scan:** no TBD/TODO; every code step shows complete code; the one non-network approximation (flat unexplained cluster instead of organ-system grouping) is called out explicitly in the spec's non-goals and in Task 4's docstring, not left vague.

**Type consistency:** `call_json(system, user, schema, *, effort, max_tokens)` is defined in Task 1 and called with that exact signature in Tasks 5/6; `assess_variant(..., in_scope=)`, `find_unexplained`, `rank` signatures in Task 4 match their calls in Task 7; `CausalityReview.to_dict()` (Task 2) is what the server returns (Task 8) and what the UI consumes (field names `assessments`, `unexplained`, `sources`, `warnings`, and per-assessment `variant`, `confidence_tier`, `supporting`, `missing`, `rubric`, `grounding_flags`, `in_scope`, `notes` all align between Task 2 dataclasses and the Task 8 renderer). `_molecular_summary` is defined once in Task 4 and reused in Tasks 6. HPO `ncbi_gene_id` matches the field name in `variant_curator/genes.py:GeneSpec`.
