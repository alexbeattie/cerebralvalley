# Variant-Curation Assistant (MVP)

Drafts an [ACMG/AMP](https://www.ncbi.nlm.nih.gov/pubmed/25741868) classification for a
**single missense or nonsense SNV**, with every piece of evidence linked to its source,
and **abstains where the evidence isn't there**.

This is the daily bottleneck in a clinical variant-curation lab: a curator opens gnomAD,
ClinVar, VEP, a prediction tool or two, and PubMed, gathers evidence for one variant, maps
it to ACMG criteria, and writes a defensible call — 30 to 60 minutes per variant, and the
queue never empties. This tool hands them a solid first draft to edit instead of a blank
page.

## Two hard rules

1. **No PHI, ever.** This is built entirely on public data (Ensembl VEP, gnomAD, ClinVar).
   It is designed so a clinical site (e.g. CHLA) could later point the same engine at real
   cases behind their own firewall without changing the code — the input is just a gene and
   an HGVS change.
2. **Guard scope.** A narrow slice that works and is measured beats a broad platform that
   half-runs. The MVP handles only **missense and nonsense SNVs** in a fixed allowlist of
   **30 well-characterized disease genes** (`variant_curator/genes.py`). Anything else is
   refused with a clear message, not guessed.

## Status: Day 1–2 (plumbing) — DONE

One variant flows end to end through **VEP → gnomAD → ClinVar**, proving the data path
before anything smart sits on top of it. Every evidence item carries a linked, openable
source URL.

```
python -m variant_curator.cli --gene BRCA1 --hgvs "c.181T>G"
```

Example output (abridged) for the pathogenic RING-domain missense *BRCA1* c.181T>G
(p.Cys61Gly):

```
[VEP]     missense_variant  p.Cys61Gly  GRCh38 17:43106487 A>C  -> gnomAD id 17-43106487-A-C
[gnomAD]  global AF 1.87e-05; FAF95 popmax 1.66e-05 (nfe); per-ancestry AC/AN listed,
          0 observations in afr/eas/sas/amr/asj/mid
[ClinVar] Pathogenic (3 stars, no conflict)  VCV000017661
```

Note the ancestry angle already visible: this variant is seen almost only in Non-Finnish
Europeans and is absent from African, East/South Asian, Admixed American, Ashkenazi, and
Middle Eastern samples. When a patient's ancestry is thinly represented at a locus, the
frequency-based criteria must be flagged low-confidence rather than asserted — that layer
lands on Day 5.

## Non-coding variants (AlphaGenome)

Deep intronic and regulatory variants are exactly what a consequence-based tool
rejects — they aren't missense/nonsense — yet they can be disease-causing by creating
cryptic splice sites, disrupting branch points, or altering transcription,
accessibility, or expression. Instead of refusing them, the pipeline now routes
intronic / splice / UTR / regulatory consequences to **[AlphaGenome](https://deepmind.google.com/science/alphagenome)**
(Google DeepMind) for splicing and regulatory interpretation, while still pulling
gnomAD frequency and ClinVar priors (both still matter on the non-coding path).

```
export ALPHAGENOME_API_KEY=...        # non-commercial research key
python -m variant_curator.cli --gene CFTR --hgvs "c.3717+12191C>T"
python -m variant_curator.cli --noncoding-demo   # illustrative deep-intronic + regulatory demos
```

What it reports, each with a linked `Source`:

- **Splicing signals** (`SPLICE_SITES`, `SPLICE_SITE_USAGE`, `SPLICE_JUNCTIONS`) — the
  primary signal for deep intronic variants, filtered to the gene of interest.
- **Regulatory signals** (`RNA_SEQ`, `CAGE`, `PROCAP`, `DNASE`, `ATAC`, `CHIP_TF`,
  `CHIP_HISTONE`) — reported tissue-first using a small, provisional per-gene disease-tissue
  hint (e.g. LDLR→liver, MYH7→heart, SCN1A→brain, CFTR→lung), falling back to the top
  tissue across all tracks.
- **ACMG mapping** — a strong calibrated splicing/regulatory effect maps to **`PP3`**, no
  predicted effect to **`BP4`**, and the grey zone abstains. Thresholds are provisional and
  need calibration against a labelled set.

Two hard rules carry over. **Abstain, don't guess:** with no key, no package, or an API
error, the tool prints `AlphaGenome skipped: ALPHAGENOME_API_KEY not set` (or the reason) and
never fabricates a score — `import variant_curator` and `pytest` both run fully offline.
**Supporting evidence only:** AlphaGenome is a research model, **not clinically validated**,
so a prediction maps to `PP3`/`BP4` at most and can never on its own drive a Pathogenic/Benign
call. The research-use caveat is surfaced wherever scores are shown.

## Architecture

```
variant_curator/
  genes.py            30-gene allowlist with MANE Select transcripts (scope control)
  models.py           evidence data models; every item carries a Source(url, retrieved_at)
  http.py             polite HTTP with Retry-After + exponential backoff
  clients/
    vep.py            HGVS -> consequence, protein change, forward-strand GRCh38 coords
    gnomad.py         variant id -> ancestry-stratified AF + FAF95 (gnomAD v4 joint)
    clinvar.py        gene + HGVS -> prior classifications, star rating, conflict flag
    alphagenome.py    non-coding path: splicing + regulatory scoring (lazy SDK, key-gated)
  acmg.py             maps AlphaGenome signals to supporting-only PP3/BP4 with a Source
  pipeline.py         orchestrates the above into one EvidenceBundle; routes by consequence
  cli.py              prints the bundle for one variant (hardcoded demo or --gene/--hgvs)
```

Design decision worth calling out: gnomAD variant ids are built from VEP's `vcf_string`
(forward-strand VCF), **not** the coding-orientation `allele_string`. For minus-strand
genes like *BRCA1* the latter would build the wrong id and produce a false "not found".

## Setup

```
python3 -m venv .venv
./.venv/bin/pip install -r requirements.txt
./.venv/bin/python -m variant_curator.cli            # hardcoded demo variant
./.venv/bin/python -m variant_curator.cli --gene SCN5A --hgvs "c.1673A>G"
```

Requires network access. VEP, gnomAD, and ClinVar are public and keyless; NCBI E-utilities
are rate-limited to ~3 req/s without a key and Ensembl VEP occasionally returns 503 under
load — the HTTP layer retries with backoff. The non-coding path additionally needs
`ALPHAGENOME_API_KEY` (see below); without it that step is skipped cleanly, not fatal.

## Roadmap

- **Day 1–2 — plumbing.** One variant end to end through gnomAD, ClinVar, VEP. ✅ done
- **Day 3–4 — the agent.** ACMG code mapping, grounded narrative, and abstention: refuse to
  assert any criterion the evidence can't support (the faithfulness-judge pattern, made
  visible in the UI, not buried).
- **Day 5 — ancestry-confidence layer + ClinVar eval harness.** Mark frequency criteria
  low-confidence when the patient's ancestry is underrepresented at the locus. Pull N
  ClinVar variants with known calls, hide the labels, run the tool, and measure agreement on
  both the ACMG criteria and the final 5-tier call. Track regressions.
- **Day 6 — thin UI + curator override loop.** Draft call, editable narrative, sources
  panel; override any criterion and watch the call update.
- **Day 7 — polish, demo, reproducibility notes** so a judge can rerun it.

Planned additional evidence sources for the agent phase: AlphaMissense (in-silico missense
prediction), ClinGen (gene–disease validity), PubMed/LitVar (variant-linked literature). The
non-coding path is handled by **AlphaGenome** (see above), which supersedes the originally
planned SpliceAI — it covers cryptic-splice effects *and* regulatory (expression /
accessibility / TSS) signal in one model, calibrated genome-wide.

## Blocked on user input before Day 3+

The build plan is gated on a short interview with the two named users (Matt, Bridget). Their
answers become the grounding spec; guessing here would waste the week. See
[`docs/interview.md`](docs/interview.md) for the exact questions.
