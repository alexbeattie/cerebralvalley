# Build brief: `trio_prioritizer` — inheritance-aware prioritization of non-coding variants in trios

You (Claude Code) are building a **new, self-contained Python package `trio_prioritizer/`** in this
worktree, alongside the existing `variant_curator/` and `regmodel/` packages. Read this brief
fully, read the existing code for style and to reuse the scoring engine, then build end to end.

## Where this comes from (the clinical need, in the user's words)

A clinical geneticist (Matt Deardorff) described the core unmet need for interpreting the
non-coding genome in rare-disease diagnostics. Two trio patterns dominate:

- **De novo**: a variant present in the affected child but in *neither* parent — the usual
  mechanism for severe dominant disease. With whole-genome sequencing these can now fall in
  introns / UTRs / regulatory regions, not just coding exons.
- **Recessive / compound heterozygous**: the child is affected, each parent is an unaffected
  carrier. Classically you find one (often coding) variant but *cannot find the second hit* —
  and increasingly that second hit is a deep-intronic, UTR, or regulatory variant that creates
  an aberrant splice site or disrupts an enhancer/promoter. Finding and scoring that second
  allele is the key pain point.

The ask: **identify intronic/intergenic/regulatory variants suggestive of disease in a trio
(de novo or recessive), score their pathogenicity, and prioritize them integrated with the
other candidate variants already on the table.** This package is that workflow. The per-variant
non-coding *scoring* already exists in this repo — you are building the inheritance logic,
the compound-het pairing, and the integrated ranking on top of it.

## Read first (reuse, don't reinvent)

- `variant_curator/models.py` — dataclass + `Source(name,url,retrieved_at,detail)` style;
  every evidence item carries provenance. Match this exactly.
- `variant_curator/clients/alphagenome.py` — **the scoring engine you will call**:
  `fetch_alphagenome(chrom, pos, ref, alt, gene_symbol) -> AlphaGenomeEvidence` with
  `.found`, `.splicing`, `.regulatory`, `.top_effect`, `.source`, and each signal's
  `.magnitude`. It degrades to `found=False` + `reason` with no `ALPHAGENOME_API_KEY`. Do NOT
  modify it.
- `variant_curator/http.py` — `now_iso()` helper.
- `regmodel/` — optional secondary scorer (sequence→activity + ISM). It needs a local sequence
  window, so treat it as *optional*: only use it when a window is supplied. AlphaGenome
  (coordinate-based) is the primary scorer because trio inputs are genomic coordinates.
- The three existing feature branches are merged into this one's base; `variant_curator` and
  `regmodel` tests must keep passing.

## Non-negotiable engineering constraints (must RUN offline)

1. **Offline-first.** The full pipeline (`demo`) and the entire test suite MUST run on a CPU
   laptop with **no network and no API key**, using a built-in **synthetic trio scenario**.
   The AlphaGenome scorer must be injectable/mockable so tests never touch the network.
2. **Abstain, don't fabricate.** When AlphaGenome is unavailable, scores are marked
   unavailable with a reason; ranking still runs on inheritance + whatever evidence exists.
   Never invent a score.
3. **Provenance on every scored claim** (carry the AlphaGenome `Source`), matching repo ethos.
4. **Research-use-only** caveat surfaced wherever pathogenicity scores are shown. This is
   decision-support, not a diagnosis.
5. **No PHI.** Synthetic/coordinate data only. (The real clinician tool will be on-prem; this
   engine must be designed so it *could* run behind a firewall unchanged.)
6. Small, focused, honest. `from __future__ import annotations`, dataclasses, why-comments.

## Package layout to build

```
trio_prioritizer/
  __init__.py         scope statement (trio inheritance + non-coding prioritization; research use)
  models.py           core dataclasses (below)
  inheritance.py      trio genotype -> inheritance mode; compound-het pairing by parent-of-origin
  scoring.py          attach a non-coding pathogenicity score to a variant via AlphaGenome
                      (injectable scorer fn; graceful degrade); optional regmodel hook
  prioritize.py       integrate inheritance + score into one ranked candidate list
  data.py             synthetic offline trio scenario generator (plants a known "answer")
  cli.py              `demo` (offline synthetic trio -> ranked candidates) + `run` (table input)
tests/
  test_inheritance.py       de novo / hom-recessive / comp-het detection; comp-het requires
                            one variant from EACH parent; genotype edge cases
  test_prioritize.py        the planted causal de novo AND the planted comp-het pair rank at top
  test_scoring_offline.py   scoring degrades cleanly with no key; Source attached when scored
  test_cli_demo_offline.py  `demo` runs end to end offline and recovers the planted answer
```

## Data model (mirror variant_curator's style)

- `Genotype` enum: `HOM_REF`, `HET`, `HOM_ALT`, `MISSING`.
- `TrioVariant`: `chrom, pos, ref, alt, gene, consequence` (reuse
  `variant_curator.models.VariantConsequence`), and `proband/mother/father: Genotype`.
  Include an `is_coding`/`is_noncoding` helper keyed off the consequence.
- `InheritanceCall`: `mode` (`DE_NOVO`, `HOMOZYGOUS_RECESSIVE`, `COMPOUND_HET`,
  `INHERITED_DOMINANT`, `UNKNOWN`), `parent_of_origin` (`maternal`/`paternal`/`None`),
  `confidence` (a simple qualitative flag — genotype-consistent vs. needs-QC), `rationale`.
- `NoncodingScore`: `available: bool`, `scorer` (e.g. "AlphaGenome"), `magnitude` (calibrated,
  0..~1 where available), `top_effect: str`, `reason` (why unavailable), `source: Optional[Source]`.
- `PrioritizedCandidate`: the variant (or the *pair*, for comp-het), its `InheritanceCall`,
  its `NoncodingScore`(s), a `combined_score: float`, a human `rationale`, and `sources()`.

## Inheritance logic (the new core)

- **De novo**: proband `HET`/`HOM_ALT`, mother `HOM_REF`, father `HOM_REF`. Flag confidence as
  "genotype-consistent (de novo calling needs QC/coverage confirmation)".
- **Homozygous recessive**: proband `HOM_ALT`, both parents `HET`.
- **Compound heterozygous**: within a single gene, **two different** `HET` variants in the
  proband where one is transmitted from the mother (mother `HET`, father `HOM_REF`) and the
  other from the father (father `HET`, mother `HOM_REF`). This is the "find the second hit"
  case — and it must work when one member of the pair is coding and the other is a
  deep-intronic/regulatory variant. Emit the *pair* as one candidate with parent-of-origin phase.
- **Inherited dominant**: present in a parent (parent `HET`/`HOM_ALT`) — for a severe-disease
  proband, deprioritize (lower prior) but still list.
- Scope: autosomal de novo + recessive (hom + comp-het) for the MVP. Note X-linked as a
  documented future extension; don't implement it now.

## Scoring (reuse the engine)

- `scoring.py` takes an injectable `scorer` callable defaulting to
  `variant_curator.clients.alphagenome.fetch_alphagenome`, so tests pass a fake. For each
  variant (prioritize the non-coding ones) it produces a `NoncodingScore` from the AlphaGenome
  splicing/regulatory magnitude + `top_effect`, carrying the `Source`. Coding variants may skip
  AlphaGenome (note that the coding path is handled elsewhere in the project).
- Optional: if a local sequence window is available, also compute a `regmodel` ISM-based delta
  as a secondary signal. Keep this strictly optional so nothing requires torch at runtime for
  the core demo.

## Prioritization / integration

- Combine an inheritance **prior** (de novo in a severe-disease proband and a completed
  comp-het pair rank high; a lone het with no second hit ranks low; inherited-dominant lower)
  with the **non-coding score magnitude** into `combined_score`. Keep the formula simple,
  transparent, and documented (name the weights as constants; note they're provisional and
  need calibration). Produce ONE ranked list that interleaves de novo, comp-het pairs, and
  homozygous candidates — i.e. it "integrates the non-coding candidates with the other
  potentially-positive variants," which is exactly the ask.
- Output each candidate with: gene, inheritance mode (+ phase), the variant(s), the non-coding
  evidence (score + top effect + source), combined score, and a one-line rationale.

## Synthetic scenario (make the demo prove itself)

`data.py` builds a seeded synthetic trio candidate set that contains a known "answer" the tool
must recover:
- a **de novo** deep-intronic/regulatory variant in a plausible gene (high scorer magnitude),
- a **compound-het pair** in a recessive gene: one coding `HET` from the father + one
  deep-intronic `HET` from the mother (the classic "second hit is non-coding" case),
- several benign inherited het variants as distractors (low magnitude / inherited from an
  unaffected parent),
so the demo/test can assert the causal de novo and the comp-het pair land at the top of the
ranking. Provide a fake scorer for offline runs that returns deterministic magnitudes keyed to
the planted variants (so `demo` works with no key), and clearly mark scores as illustrative
when the real AlphaGenome scorer is absent.

## CLI / definition of done

- `python -m trio_prioritizer.cli demo` runs fully offline: builds the synthetic trio, scores
  (fake/degraded scorer), prints the ranked candidate list with inheritance mode, non-coding
  evidence, sources, and the research-use caveat; recovers the planted de novo + comp-het pair
  at the top. Writes a `candidates.json` with provenance.
- `python -m trio_prioritizer.cli run --table trio.tsv` reads a simple TSV of trio variants
  (chrom,pos,ref,alt,gene,consequence,proband_gt,mother_gt,father_gt) and ranks them; uses the
  real AlphaGenome scorer when `ALPHAGENOME_API_KEY` is set, otherwise degrades with a clear note.
- `python -c "import trio_prioritizer"` works with no network / no torch / no key.
- `pytest` passes offline (AlphaGenome mocked); existing `variant_curator` + `regmodel` tests
  still pass.
- Update `requirements.txt` only if needed. Add a `trio_prioritizer` section to `README.md`
  (the trio patterns, the offline demo, integration with the AlphaGenome scorer, the
  research-use caveat).
- Commit on the current branch `feat/trio-noncoding-prioritizer` with a clear message. Do NOT push.

When done, print a concise summary of files created and exactly how to run the offline demo.
