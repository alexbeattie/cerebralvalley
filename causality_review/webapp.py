"""Minimal web UI for the causality review — Matt's sketch, made real.

Zero third-party dependencies: Python's stdlib HTTP server serves one static
page and a single JSON endpoint that runs the real scoring engine. The point is
that the UI calls the same `review_causality` the CLI does — no mock data.

    python -m causality_review.webapp          # then open http://localhost:8000
    python -m causality_review.webapp --port 8080

Input: the variants the lab reported + the patient's features (free text or
HP:xxxxxxx). Output: the ranked causality list, with explained vs unexplained
features, source links, and the second-cause / dual-diagnosis flags.
"""

from __future__ import annotations

import argparse
import json
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from .clients.hpo import resolve_term
from .clients.llm import LLMError, is_configured
from .engine import review_causality
from .extract import extract_from_notes
from .http import get_client
from .models import PatientProfile, Phenotype, ReportedVariant

STATIC_DIR = Path(__file__).parent / "static"


def _parse_variants(raw: list[str]) -> list[ReportedVariant]:
    variants: list[ReportedVariant] = []
    for line in raw:
        line = line.strip()
        if not line:
            continue
        gene, _, hgvs = line.partition(":")
        variants.append(ReportedVariant(gene=gene.strip(), hgvs_c=hgvs.strip()))
    return variants


def _run_review(payload: dict) -> dict:
    """Resolve features to HPO, run the engine, and serialize the report."""

    variants = _parse_variants(payload.get("variants", []))
    feature_lines = [f.strip() for f in payload.get("features", []) if f.strip()]

    with get_client() as client:
        phenotypes: list[Phenotype] = []
        unresolved: list[str] = []
        for raw in feature_lines:
            term = resolve_term(client, raw)
            if term is None:
                unresolved.append(raw)
            else:
                phenotypes.append(term)

        if not variants:
            return {"error": "No reported variants given."}
        if not phenotypes:
            return {"error": "No patient features could be resolved to HPO terms.", "unresolved": unresolved}

        patient = PatientProfile(phenotypes=phenotypes)
        report = review_causality(client, patient, variants)

    return {
        "patient": [{"hpo_id": p.hpo_id, "label": p.label} for p in report.patient.phenotypes],
        "unresolved": unresolved,
        "fits": [
            {
                "rank": i,
                "gene": f.variant.gene,
                "hgvs": f.variant.hgvs_c,
                "label": f.variant.label,
                "stars": f.tier.stars,
                "tier": f.tier.display,
                "tier_level": int(f.tier),
                "score": round(f.score, 3),
                "explained_count": len(f.explained),
                "total": len(report.patient.phenotypes),
                "explained": [
                    {"label": m.phenotype.label, "exact": m.exact, "via": m.via}
                    for m in f.explained
                ],
                "unexplained": [p.label for p in f.unexplained],
                "diseases": f.diseases[:3],
                "knowledge_found": f.knowledge_found,
                "source": f.source.url if f.source else "",
            }
            for i, f in enumerate(report.fits, 1)
        ],
        "residual": [p.label for p in report.residual_unexplained],
        "flags": report.flags,
    }


def _run_extract(payload: dict) -> dict:
    """Turn free-text clinical notes into validated HPO phenotypes via the LLM edge."""

    notes = (payload.get("notes") or "").strip()
    if not notes:
        return {"error": "No clinical notes provided."}
    with get_client() as client:
        try:
            result = extract_from_notes(client, notes)
        except LLMError as exc:
            return {"error": str(exc)}
    return {
        "phenotypes": [
            {"phrase": e.phrase, "hpo_id": e.phenotype.hpo_id, "label": e.phenotype.label}
            for e in result.phenotypes
        ],
        "ungrounded": result.ungrounded,
    }


class _Server(ThreadingHTTPServer):
    # Lets us rebind immediately after a restart instead of hitting "Address
    # already in use" while the old socket lingers in TIME_WAIT.
    allow_reuse_address = True


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *args):  # keep the console clean
        pass

    def _send(self, code: int, body: bytes, content_type: str) -> None:
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        if self.path in ("/", "/index.html"):
            html = (STATIC_DIR / "index.html").read_bytes()
            self._send(200, html, "text/html; charset=utf-8")
        elif self.path == "/api/config":
            # Lets the UI show/hide the AI panel based on whether a key is set.
            self._send(200, json.dumps({"ai_enabled": is_configured()}).encode(), "application/json")
        else:
            self._send(404, b"not found", "text/plain")

    def do_POST(self) -> None:
        routes = {"/api/review": _run_review, "/api/extract": _run_extract}
        handler = routes.get(self.path)
        if handler is None:
            self._send(404, b"not found", "text/plain")
            return
        length = int(self.headers.get("Content-Length", 0))
        try:
            payload = json.loads(self.rfile.read(length) or b"{}")
            result = handler(payload)
        except Exception as exc:  # pylint: disable=broad-except  # surface to UI, don't 500 silently
            result = {"error": f"{type(exc).__name__}: {exc}"}
        self._send(200, json.dumps(result).encode(), "application/json")


def _bind(port: int, tries: int = 20) -> _Server:
    """Bind to `port`, or the next free port if it's already in use."""

    last: OSError | None = None
    for candidate in range(port, port + tries):
        try:
            return _Server(("127.0.0.1", candidate), Handler)
        except OSError as exc:  # EADDRINUSE (48/98): try the next port
            last = exc
            if candidate == port:
                print(f"Port {port} is busy; trying {port + 1}…")
    raise SystemExit(f"No free port in {port}-{port + tries - 1}: {last}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Causality review web UI.")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--no-open", action="store_true", help="don't auto-open a browser")
    args = parser.parse_args()

    server = _bind(args.port)
    actual_port = server.server_address[1]
    url = f"http://localhost:{actual_port}"
    print(f"Causality review UI running at {url}  (Ctrl-C to stop)")
    if not args.no_open:
        try:
            webbrowser.open(url)
        except Exception:  # pylint: disable=broad-except  # opening a browser is best-effort
            pass
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped.")
        server.shutdown()


if __name__ == "__main__":
    main()
