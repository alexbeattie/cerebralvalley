"""Shared HTTP helpers with polite defaults for public APIs."""

from __future__ import annotations

import time
from datetime import datetime, timezone

import httpx

USER_AGENT = "variant-curator/0.1 (research prototype; public data only)"
DEFAULT_TIMEOUT = 30.0


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def get_client() -> httpx.Client:
    return httpx.Client(
        timeout=DEFAULT_TIMEOUT,
        headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
        follow_redirects=True,
    )


def _retry_after(resp: httpx.Response, attempt: int) -> float:
    """Honor the server's Retry-After header (Ensembl sends it on 429/503)."""
    header = resp.headers.get("Retry-After")
    if header:
        try:
            return max(float(header), 1.0)
        except ValueError:
            pass
    return min(2.0 * (2**attempt), 30.0)  # exponential backoff, capped


def get_json(client: httpx.Client, url: str, *, params: dict | None = None, retries: int = 5) -> dict | list:
    last_exc: Exception | None = None
    for attempt in range(retries):
        try:
            resp = client.get(url, params=params)
            resp.raise_for_status()
            return resp.json()
        except httpx.HTTPStatusError as exc:
            last_exc = exc
            if exc.response.status_code in (429, 500, 502, 503, 504):
                time.sleep(_retry_after(exc.response, attempt))
                continue
            raise
        except (httpx.TransportError, httpx.DecodingError) as exc:
            last_exc = exc
            time.sleep(min(2.0 * (2**attempt), 30.0))
    raise RuntimeError(f"GET {url} failed after {retries} attempts: {last_exc}")


def post_json(client: httpx.Client, url: str, payload: dict, *, retries: int = 3) -> dict:
    last_exc: Exception | None = None
    for attempt in range(retries):
        try:
            resp = client.post(url, json=payload)
            resp.raise_for_status()
            return resp.json()
        except httpx.HTTPStatusError as exc:
            last_exc = exc
            if exc.response.status_code in (429, 500, 502, 503, 504):
                time.sleep(1.5 * (attempt + 1))
                continue
            raise
        except (httpx.TransportError, httpx.DecodingError) as exc:
            last_exc = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"POST {url} failed after {retries} attempts: {last_exc}")
