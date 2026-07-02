# Build brief: AlphaGenome interpretation of deep intronic & regulatory variants

You (Claude Code) are extending an existing Python project, `variant_curator`, in this
worktree. Implement the feature described below end to end. This brief is self-contained —
read the existing code first, then build.

## Context: what the project is

`variant_curator` is a scoped ACMG/AMP variant-curation assistant. Today it drafts an
ACMG classification for a **single missense or nonsense SNV** in a fixed allowlist of 30
disease genes, assembling evidence from **Ensembl VEP, gnomAD, and ClinVar**, with a linked
`Source` (URL + retrieved-at timestamp) on **every** piece of evidence, and it **abstains**
where evidence is missing.

Read these files before writing anything:

- `variant_curator/__init__.py` — scope statement.
- `variant_curator/models.py` — dataclasses. Note `Source`, `VariantInput`, `VepEvidence`,
  `EvidenceBundle`, and the `VariantConsequence` enum (only MISSENSE/NONSENSE/OTHER today).
- `variant_curator/genes.py` — the 30-gene allowlist with MANE Select transcripts.
- `variant_curator/http.py` — polite HTTP helpers with retry/backoff.
- `variant_curator/clients/vep.py`, `clients/gnomad.py`, `clients/clinvar.py` — the existing
  evidence clients. Match their style exactly.
- `variant_curator/pipeline.py` — `assemble_evidence()`; note `IN_SCOPE_CONSEQUENCES` and how
  out-of-scope consequences are currently **rejected** with a warning.
- `variant_curator/cli.py` — how a bundle is printed.
- `README.md` — tone, and the roadmap (AlphaMissense/SpliceAI were already planned; this
  feature supersedes SpliceAI with AlphaGenome for the non-coding path).

## The goal

Deep intronic and regulatory (non-coding) variants are exactly what the current tool
**rejects**: they are not missense/nonsense, so `pipeline.assemble_evidence()` bails with an
"Out of MVP scope" warning. These variants can still be disease-causing — they act by
creating cryptic splice sites, disrupting branch points / splicing regulatory elements, or
altering transcription, accessibility, or expression. Consequence-based tools and protein
missense predictors are blind to them.

Add a new evidence source, **AlphaGenome** (Google DeepMind), that interprets these variants,
and route intronic + regulatory variants to it instead of rejecting them.

## AlphaGenome SDK — exact API to use

Package: `alphagenome` (PyPI). Install: `pip install alphagenome`. Requires an API key
(non-commercial research use). **Read the key from the environment variable
`ALPHAGENOME_API_KEY`.** If the key is missing, do NOT crash — return an evidence object with
`found=False` and a warning explaining the key is unset, so the rest of the pipeline and all
imports still work offline.

Import lazily *inside* the client function (not at module top level) so that importing
`variant_curator` never fails when `alphagenome` isn't installed. Catch `ImportError` and
degrade gracefully with a warning.

Core usage (this is the real, current API — use it as-is):

```python
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

model = dna_client.create(API_KEY)

variant = genome.Variant(
    chromosome='chr17',            # note the 'chr' prefix
    position=43106487,             # 1-based GRCh38
    reference_bases='A',
    alternate_bases='C',
)
interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

scores = model.score_variant(
    interval=interval,
    variant=variant,
    variant_scorers=SELECTED_SCORERS,
    organism=dna_client.Organism.HOMO_SAPIENS,
)
df = variant_scorers.tidy_scores([scores])   # tidy long-format pandas DataFrame
```

`variant_scorers.RECOMMENDED_VARIANT_SCORERS` is a dict keyed by string. Use these keys:

- **Splicing (primary signal for deep intronic):** `'SPLICE_SITES'`, `'SPLICE_SITE_USAGE'`,
  `'SPLICE_JUNCTIONS'`.
- **Regulatory (primary signal for promoter/enhancer/UTR):** `'RNA_SEQ'` (expression
  log-fold-change), `'CAGE'`, `'PROCAP'` (TSS activity), `'DNASE'`, `'ATAC'` (accessibility),
  `'CHIP_TF'`, `'CHIP_HISTONE'`.

`tidy_scores` returns a DataFrame with columns including: `variant_id`, `gene_name`,
`gene_id`, `output_type`, `variant_scorer`, `track_name`, `ontology_curie`, `gtex_tissue`,
`biosample_name`, `raw_score`, and `quantile_score` when available. Splicing scores are only
returned for genes overlapping the variant, so **filter the tidy DataFrame to the gene of
interest** (`gene_name == <symbol>`) before summarizing.

Use `quantile_score` (calibrated against the genome-wide distribution) as the primary
magnitude for interpretation when present; fall back to `raw_score`. For splicing scorers,
scores are non-directional (always positive). For `RNA_SEQ` LFC, sign matters.

Notes / gotchas:
- Chromosome must be `'chr'`-prefixed (`chr1`…`chr22`, `chrX`, `chrY`).
- Coordinates are 1-based GRCh38, matching what VEP already returns in this project.
- Keep it to a handful of `score_variant` calls per variant (the API suits 1000s of
  predictions, not millions) — one call with the splicing scorer set and one with the
  regulatory scorer set is fine.

## What to build (deliverables)

1. **`variant_curator/models.py` additions:**
   - Extend `VariantConsequence` with the non-coding classes routed to AlphaGenome, e.g.
     `INTRON` (`intron_variant`), `SPLICE_REGION`, `SPLICE_DONOR`, `SPLICE_ACCEPTOR`,
     `FIVE_PRIME_UTR`, `THREE_PRIME_UTR`, `UPSTREAM` (promoter-proximal), `DOWNSTREAM`,
     `REGULATORY`/`INTERGENIC`. Keep `OTHER` as the true fallback.
   - New dataclasses (mirror the existing evidence-object style, each with an optional
     `Source`):
     - `SplicingSignal` — per-scorer summary: scorer name, output type, top affected
       gene, max |quantile_score|, max |raw_score|, direction where applicable, and a short
       human interpretation (e.g. "predicted cryptic acceptor gain 212 bp into intron 4").
     - `RegulatorySignal` — per-modality summary: modality (RNA_SEQ/DNASE/…), top tissue /
       biosample (`gtex_tissue` or `biosample_name` + `ontology_curie`), signed
       quantile/raw score, direction (up/down for expression).
     - `AlphaGenomeEvidence` — `found: bool`, `variant_id`, `model_version`/`api` string,
       `splicing: list[SplicingSignal]`, `regulatory: list[RegulatorySignal]`, a computed
       `top_effect` summary, `source: Source`, and a `research_use_only: bool = True` flag.
   - Add `alphagenome: Optional[AlphaGenomeEvidence]` to `EvidenceBundle` and include it in
     `EvidenceBundle.sources()`.

2. **`variant_curator/clients/alphagenome.py`:** new client, same shape/idioms as the other
   clients (module docstring explaining *why*, `from __future__ import annotations`, a single
   `fetch_alphagenome(...)` entry point). It takes GRCh38 coordinates (chrom, pos, ref, alt —
   the same fields VEP already resolves) plus the gene symbol, calls the splicing and
   regulatory scorer sets, filters `tidy_scores` output to the gene, and reduces it to the
   `SplicingSignal` / `RegulatorySignal` summaries. Read `ALPHAGENOME_API_KEY` from env;
   lazy import; graceful degradation (return `found=False` + reason) if the key or package is
   missing, or on API error. Set the `Source` URL to
   `https://deepmind.google.com/science/alphagenome` with a retrieved-at timestamp and the
   variant id in `detail`.
   - Provide an optional disease-relevant tissue hint: add a per-gene tissue/ontology mapping
     (UBERON/GTEx term) so regulatory scoring can be reported tissue-first for the gene's
     disease context (e.g. LDLR→liver, MYH7/MYBPC3→heart, SCN1A/MECP2→brain, CFTR→lung/
     epithelium). Keep it small and clearly marked as a heuristic; fall back to reporting the
     top tissue across all tracks when no mapping exists.

3. **`variant_curator/pipeline.py`:** define a set of non-coding consequences
   (`NONCODING_CONSEQUENCES`). Change the flow so that:
   - missense/nonsense → existing gnomAD + ClinVar path (unchanged);
   - a non-coding consequence → call `fetch_alphagenome(...)` using the GRCh38 coords VEP
     already returns, attach `bundle.alphagenome`, and also still fetch gnomAD + ClinVar
     (frequency and prior classifications are just as relevant for non-coding variants);
   - anything genuinely unhandled → keep the existing "out of scope" warning.
   Do NOT break the existing missense/nonsense behaviour.

4. **ACMG mapping (lightweight, provenance-first):** add a small module
   `variant_curator/acmg.py` (or a function) that maps AlphaGenome signals to ACMG-style
   supporting evidence and returns structured, sourced criterion hits — **supporting strength
   only**, never a standalone pathogenic assertion:
   - strong predicted splicing disruption (high calibrated splicing score) → `PP3`
     (computational evidence of a deleterious splicing effect);
   - strong predicted regulatory effect (large expression/accessibility shift in the disease-
     relevant tissue) → `PP3`;
   - no predicted effect across splicing + regulatory → `BP4`.
   Use conservative quantile thresholds, name them as constants with a comment that they are
   provisional and need calibration against a labelled set, and attach the AlphaGenome
   `Source` to each emitted criterion. Include a clear caveat string that AlphaGenome is a
   research model, **not** clinically validated, so a curator treats it as supporting only.

5. **`variant_curator/cli.py`:** print the AlphaGenome block when present (splicing signals,
   regulatory signals with tissue, top effect, the source URL, and the research-use caveat),
   plus any ACMG PP3/BP4 hits with their source. Add a couple of demo non-coding variants
   (with a comment noting they are illustrative), e.g. a deep intronic splicing example.

6. **`requirements.txt`:** add `alphagenome` (and `pandas` if the summarization needs it —
   the SDK already depends on pandas, but list it if you import it directly). Pin nothing you
   can't verify; a lower bound is fine.

7. **Tests** under `tests/` using `pytest` + `unittest.mock`:
   - `fetch_alphagenome` returns `found=False` with a clear reason when the API key is unset
     (no network).
   - Given a mocked `tidy_scores` DataFrame (construct a small pandas DataFrame fixture with
     the documented columns), the client reduces it correctly into `SplicingSignal` /
     `RegulatorySignal` and picks the right top effect.
   - The ACMG mapper emits `PP3` for a high-score fixture and `BP4` for an all-low fixture,
     each carrying a `Source`.
   - The pipeline routes a non-coding consequence to the AlphaGenome path (mock the client)
     and still leaves the missense/nonsense path unchanged.
   Mock the `alphagenome` SDK so tests run offline with no key and no package installed
   (e.g. patch `variant_curator.clients.alphagenome._score` or inject a fake).

8. **`README.md`:** add a short section documenting the new non-coding capability, the
   `ALPHAGENOME_API_KEY` env var, an example invocation, and the research-use-only caveat.
   Update the roadmap bullet that mentioned SpliceAI.

## Hard constraints (do not violate)

- **Provenance on every claim.** Every evidence item and every emitted ACMG criterion must
  carry a `Source` with an openable URL and a retrieved-at timestamp. This is the core
  design principle of the project.
- **Abstain, don't guess.** If AlphaGenome can't be reached (no key, no package, API error),
  say so explicitly via `found=False` + a warning; never fabricate scores.
- **Supporting evidence only.** AlphaGenome predictions map at most to `PP3`/`BP4`. Never let
  a prediction alone drive a Pathogenic/Benign call. Surface the "research model, not
  clinically validated" caveat wherever scores are shown.
- **No PHI.** Public data + coordinates only, consistent with the rest of the project.
- **Don't break the existing missense/nonsense flow** or the existing tests/CLI behaviour.
- Match the existing code style: `from __future__ import annotations`, dataclasses, module
  docstrings that explain *why*, small focused functions, no needless comments.

## Definition of done

- `python -c "import variant_curator"` works with **no** `alphagenome` installed and **no**
  API key set (lazy import + graceful degradation).
- `pytest` passes offline (SDK mocked).
- `python -m variant_curator.cli --gene <g> --hgvs <deep-intronic-c.>` runs; without a key it
  prints a clean "AlphaGenome skipped: ALPHAGENOME_API_KEY not set" style message rather than
  crashing, and with a key it prints splicing + regulatory signals, top effect, PP3/BP4, and
  sources.
- New code committed on this branch (`feat/alphagenome-noncoding`) with a clear message. Do
  not push.

When done, print a concise summary of the files you created/changed and how to run it.
