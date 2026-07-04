"""The `demo` command runs end to end offline and recovers the planted answer."""

from __future__ import annotations

import json

from trio_prioritizer import cli
from trio_prioritizer.data import causal_variant_ids


def test_demo_runs_offline_and_recovers_answer(tmp_path, capsys, monkeypatch):
    # Guarantee no key is present, so nothing can reach the network.
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    out = tmp_path / "candidates.json"

    rc = cli.main(["demo", "--out", str(out)])
    assert rc == 0  # the demo self-check passed: causal variants at the top

    captured = capsys.readouterr().out
    assert "RESEARCH USE ONLY" in captured  # caveat surfaced wherever scores are shown
    assert "recovered at top: YES" in captured

    # candidates.json written with provenance and the research flag.
    payload = json.loads(out.read_text())
    assert payload["research_use_only"] is True
    assert payload["candidates"], "expected ranked candidates"
    top_two = payload["candidates"][:2]
    top_ids = {v["variant_id"] for c in top_two for v in c["variants"]}
    assert causal_variant_ids().issubset(top_ids)

    # Every scored candidate carries at least one Source (provenance on every claim).
    scored = [c for c in payload["candidates"]
              if any(v["noncoding_score"]["available"] for v in c["variants"])]
    assert scored and all(c["sources"] for c in scored)


def test_run_table_offline_degrades(tmp_path, capsys, monkeypatch):
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    tsv = tmp_path / "trio.tsv"
    tsv.write_text(
        "chrom\tpos\tref\talt\tgene\tconsequence\tproband_gt\tmother_gt\tfather_gt\n"
        "2\t165991200\tG\tA\tSCN1A\tintron_variant\t0/1\t0/0\t0/0\n"
        "12\t102852100\tC\tT\tPAH\tmissense_variant\t0/1\t0/0\t0/1\n"
        "12\t102912400\tA\tG\tPAH\tintron_variant\t0/1\t0/1\t0/0\n"
    )
    out = tmp_path / "out.json"
    rc = cli.main(["run", "--table", str(tsv), "--out", str(out)])
    assert rc == 0
    captured = capsys.readouterr().out
    assert "no ALPHAGENOME_API_KEY" in captured  # honest degrade note
    payload = json.loads(out.read_text())
    # De novo SCN1A ranks first even with no non-coding score available.
    assert payload["candidates"][0]["gene"] == "SCN1A"
    assert payload["candidates"][0]["inheritance_mode"] == "de_novo"
