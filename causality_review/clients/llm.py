"""Anthropic (Claude) client for clinical-note phenotype extraction.

Used only at the ingestion edge: turn free-text clinical notes into candidate
phenotype phrases. The LLM only *proposes* phrases -- every phrase is then
validated and normalized to a real HPO term by the deterministic HPO search
(see `extract.py`), so the model never invents an HPO id and the scoring core
stays AI-free and sourced.

No new dependency: calls the Anthropic Messages REST API over httpx. The API
key is read from the ANTHROPIC_API_KEY env var (never hard-coded, never logged);
the model is overridable via ANTHROPIC_MODEL.
"""

from __future__ import annotations

import json
import os
import re

import httpx

ANTHROPIC_URL = "https://api.anthropic.com/v1/messages"
DEFAULT_MODEL = "claude-sonnet-5"  # override with ANTHROPIC_MODEL
ANTHROPIC_VERSION = "2023-06-01"

_SYSTEM = (
    "You are a clinical phenotyping assistant. Given a free-text clinical note, "
    "extract the patient's OBSERVED phenotypic abnormalities as short, canonical "
    "clinical terms suitable for looking up in the Human Phenotype Ontology (HPO). "
    "Rules: include only positive findings actually present in THIS patient; exclude "
    "explicitly negated findings, normal findings, family history, and procedures/"
    "medications. Prefer standard terminology (e.g. 'seizure', 'hypertrophic "
    "cardiomyopathy', 'global developmental delay'). "
    "Respond with ONLY a JSON array of strings, nothing else."
)

_VARIANT_SYSTEM = (
    "You extract reported genetic variants from a diagnostic lab report. For each "
    "variant the report lists, return an object with: 'gene' (HGNC symbol), 'hgvs' "
    "(the coding-level HGVS such as 'c.3637C>T', empty string if absent), and "
    "'classification' (the lab's call, e.g. 'Pathogenic', 'Likely pathogenic', 'VUS', "
    "empty string if absent). Include only variants the lab actually reports for this "
    "patient; ignore methodology, references, and genes listed only as panel coverage. "
    "Respond with ONLY a JSON array of objects, nothing else."
)


class LLMError(RuntimeError):
    """Raised when the LLM call cannot be made or its output can't be parsed."""


def is_configured() -> bool:
    return bool(os.environ.get("ANTHROPIC_API_KEY"))


def _parse_json_array(text: str) -> list:
    """Pull a JSON array out of the model's reply, tolerating surrounding prose."""

    match = re.search(r"\[.*\]", text, re.DOTALL)
    if not match:
        raise LLMError(f"Model did not return a JSON array. Got: {text[:200]!r}")
    try:
        return json.loads(match.group(0))
    except json.JSONDecodeError as exc:
        raise LLMError(f"Could not parse model JSON: {exc}") from exc


def _parse_array(text: str) -> list[str]:
    """A JSON array of strings (phenotype phrases)."""

    return [str(x).strip() for x in _parse_json_array(text) if str(x).strip()]


def _call(client: httpx.Client, system: str, user: str, *, model: str | None = None) -> str:
    """One Anthropic Messages call; returns the concatenated text output."""

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        raise LLMError(
            "ANTHROPIC_API_KEY is not set. Export it (export ANTHROPIC_API_KEY=sk-ant-...) "
            "or add it to a .env file to enable AI extraction, or enter data manually."
        )

    payload = {
        "model": model or os.environ.get("ANTHROPIC_MODEL", DEFAULT_MODEL),
        "max_tokens": 1024,
        "system": system,
        "messages": [{"role": "user", "content": user}],
    }
    headers = {
        "x-api-key": api_key,
        "anthropic-version": ANTHROPIC_VERSION,
        "content-type": "application/json",
    }

    try:
        resp = client.post(ANTHROPIC_URL, json=payload, headers=headers, timeout=60.0)
        resp.raise_for_status()
    except httpx.HTTPStatusError as exc:
        # Surface the API's message (e.g. auth/model errors) without the key.
        detail = exc.response.text[:300]
        raise LLMError(f"Anthropic API error {exc.response.status_code}: {detail}") from exc
    except httpx.HTTPError as exc:
        raise LLMError(f"Could not reach Anthropic API: {exc}") from exc

    data = resp.json()
    blocks = data.get("content", [])
    text = "".join(b.get("text", "") for b in blocks if b.get("type") == "text")
    if not text:
        raise LLMError("Empty response from the model.")
    return text


def extract_phenotype_phrases(client: httpx.Client, notes: str, *, model: str | None = None) -> list[str]:
    """Ask Claude for candidate phenotype phrases from a clinical note."""

    return _parse_array(_call(client, _SYSTEM, notes, model=model))


def extract_variants_raw(client: httpx.Client, report_text: str, *, model: str | None = None) -> list[dict]:
    """Ask Claude for the reported variants (gene/hgvs/classification) from lab-report text."""

    out: list[dict] = []
    for obj in _parse_json_array(_call(client, _VARIANT_SYSTEM, report_text, model=model)):
        if isinstance(obj, dict) and obj.get("gene"):
            out.append({
                "gene": str(obj.get("gene", "")).strip(),
                "hgvs": str(obj.get("hgvs", "")).strip(),
                "classification": str(obj.get("classification", "")).strip(),
            })
    return out
