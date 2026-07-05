# Variant-Curation Assistant (MVP)

> **This branch (`feat/causality-review`) adds a second tool: `causality_review`.**
> `variant_curator` is **step 1** (the lab's variant classification, described below).
> `causality_review` is **step 2** — the clinician's causality review: given the variants
> the lab already reported plus the patient's phenotype, rank which one actually explains
> THIS patient. See **[`docs/causality-review.md`](docs/causality-review.md)** and the
> [Causality review](#causality-review-step-2) section below.

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
  pipeline.py         orchestrates the above into one EvidenceBundle; enforces scope
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

Requires network access. All APIs are public and keyless. NCBI E-utilities are rate-limited
to ~3 req/s without a key; Ensembl VEP occasionally returns 503 under load — the HTTP layer
retries with backoff.

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

Planned additional evidence sources for the agent phase: AlphaMissense and SpliceAI
(in-silico prediction), ClinGen (gene–disease validity), PubMed/LitVar (variant-linked
literature).

## Blocked on user input before Day 3+

The build plan is gated on a short interview with the two named users (Matt, Bridget). Their
answers become the grounding spec; guessing here would waste the week. See
[`docs/interview.md`](docs/interview.md) for the exact questions.

---

# Causality review (step 2)

The clinician's question, downstream of classification: the lab reports 2–5 suspicious
variants; the clinician holds the patient's full phenotype and asks **"does this variant
actually explain THIS patient?"** `causality_review` runs the match in reverse — for each
reported variant, how well does its gene's known disease phenotype explain the patient's
HPO features? — and ranks them, separating explained from unexplained features and flagging
a possible second cause / genome re-analysis. Full design: [`docs/causality-review.md`](docs/causality-review.md).

### Web UI (Matt's sketch)

A point-and-click version of the review: enter the reported variants and the patient's
features, get the ranked causality list back. Zero third-party dependencies — Python's
stdlib server serves one page and calls the *same* engine the CLI does (no mock data).

```
./run_ui.sh                    # opens http://localhost:8000 (or the next free port)
./run_ui.sh --port 8080        # pick a port
```

If the port is busy it automatically falls back to the next free one, so "Address already
in use" won't stop you. (Equivalent long form: `./.venv/bin/python -m causality_review.webapp`.)

**Optional AI phenotype extraction.** Set an Anthropic key and the UI gains a "Clinical
notes" box: paste free text, click *Extract phenotypes with AI*, and it fills the features
list. The LLM only **proposes** phrases; each is then grounded in a real HPO term by the
deterministic ontology search (the model never emits an HPO id), and the scoring stays
AI-free and sourced. Without a key the app runs exactly as before (the panel stays hidden).

```
export ANTHROPIC_API_KEY=sk-ant-...        # enables the notes -> HPO extraction
export ANTHROPIC_MODEL=claude-3-5-sonnet-latest   # optional override
./run_ui.sh

# CLI equivalent:
./.venv/bin/python -m causality_review.cli --variant SCN1A:c.3637C>T \
    --notes "4yo with recurrent febrile seizures, developmental delay, gait ataxia"
```

### CLI

```
# hardcoded demo case
./.venv/bin/python -m causality_review.cli

# your own case: reported variants + patient features (free text or HP:xxxxxxx)
./.venv/bin/python -m causality_review.cli \
    --variant SCN1A:c.3637C>T --variant MYH7:c.1063G>A \
    --hpo "seizures" --hpo "global developmental delay" --hpo "ataxia"

./demo.sh                 # four guided cases (best fit, partial+orphan, two partials, dual dx)
```

Key properties:

- **Ontology-aware matching** — walks the HPO `is_a` graph so a gene annotated to a broad
  term (*Seizure*) explains a specific feature (*Focal-onset seizure*), while excluding
  organ-system container nodes so a cardiac gene can't "explain" any cardiac feature.
- **Honest partials** — a partial match reads as `Partial (explains 3/5)`, not "close enough".
- **Second-cause + dual-diagnosis flags** — surfaces features no reported variant explains,
  and the ~5% case where two variants are each partly responsible.
- **Sourced + abstaining** — every gene–phenotype link carries an openable HPO URL; a gene
  with no retrievable knowledge is marked unscored, not guessed.

## Eval

A ranking sanity/regression harness: 9 solved cases (real classic HPO features), each
scored against a panel of every fixture gene, measuring how often the true gene ranks first.

```
./.venv/bin/python -m causality_review.eval
```

Current result: **top-1 89%, top-3 100%, MRR 0.926** (9 cases, 9-gene panel). This tests
ranking and discrimination against live HPO data; it is not a held-out clinical validation
(phenotypes are curated classic features matched against the same HPO annotations), so a real
deployment must still be measured on real solved cases. The one miss (RYR1 malignant
hyperthermia, ranked 3rd) is honest: MH is an anesthesia-triggered reaction whose baseline
HPO features overlap other myopathies.

Offline unit tests (no network, fake ontology):

```
./.venv/bin/python -m unittest discover -s tests
```
