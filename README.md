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
python -m variant_curator.cli --gene MYBPC3 --hgvs "c.1224-52G>A"   # real intronic, ClinVar P/LP
python -m variant_curator.cli --noncoding-demo   # illustrative intronic + splice-site demos
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

## `regmodel` — a miniature MPRA sequence→activity model with single-base ISM

A second, self-contained package in this repo (`regmodel/`). Where `variant_curator` curates
one variant from public databases, `regmodel` reproduces in miniature what **Katie Pollard's
lab** does with deep learning: train a compact CNN on a **massively parallel reporter assay
(MPRA)** that maps short DNA sequences to regulatory activity, then run **in-silico
mutagenesis (ISM)** — mutate every base and read off a per-position mutation-effect map (the
approach behind PARM, Nature 2025, and SuPreMo/Akita). The differentiated step is a
**head-to-head cross-check against [AlphaGenome](https://deepmind.google.com/science/alphagenome)**,
a 1 Mb genome foundation model, on the same variants: *where does a small task-specific MPRA
model agree or disagree with a foundation model?*

### Offline demo (no network, no GPU, no API key)

```
python -m regmodel.cli demo          # or: python -m regmodel.cli demo --fast
```

This trains a Basset/DeepSTARR-like CNN on a **built-in synthetic MPRA** (random sequences
with a planted AP-1-like activator motif `TGACTCA` and a weaker repressor; activity = weighted
motif content + Gaussian noise, fully seeded), prints held-out **Pearson/Spearman/MSE**, runs
ISM on a held-out sequence, and writes everything to `artifacts/`:

- `model.pt` + `model.json` — weights and a sidecar recording arch, hyperparameters, data
  provenance, seed, timestamp, and metrics (the full recipe for the numbers).
- `ism_heatmap.png`, `ism_importance.png` — the (L×4) mutation-effect map and per-base
  importance track; on the synthetic task the ISM peak lands on the planted motif (the tests
  assert this).
- `metrics.json`, `comparison.csv`, `comparison.json` — run-level provenance and the
  cross-check table.

Other subcommands:

```
python -m regmodel.cli train --out artifacts            # train + save model + sidecar
python -m regmodel.cli ism --seq TGACTCAGCTAGCTAGCTAG   # ISM map for one sequence
python -m regmodel.cli variant --seq ACGT... --pos 100 --alt G   # predicted activity delta
python -m regmodel.cli compare --model artifacts/model.pt        # AlphaGenome cross-check
```

### AlphaGenome cross-check

`regmodel.compare` **reuses `variant_curator`'s AlphaGenome client verbatim**
(`fetch_alphagenome`). For each variant it lines up our model's ISM delta on a local sequence
window against AlphaGenome's calibrated splicing/regulatory magnitude, signs the agreement,
and correlates them. The same graceful-degradation rule carries over: with no
`ALPHAGENOME_API_KEY`, no SDK, or an API error, the external half is **skipped with a recorded
reason and `Source`** — we still emit our-model rows and **never fabricate** an AlphaGenome
score. Every AlphaGenome-derived row carries its linked `Source`.

```
export ALPHAGENOME_API_KEY=...        # non-commercial research key; without it, skips cleanly
python -m regmodel.cli compare --out artifacts
```

### Research use only

By default this is a **toy model trained on synthetic data** with a planted motif — enough to
demonstrate the train→ISM→cross-check pattern and to test it deterministically, **not** to draw
biological conclusions. Real conclusions require loading a real MPRA: `regmodel.data` ships a
network-gated `load_real_mpra` stub naming candidate public datasets (lentiMPRA developing
human brain [Pollard/Ahituv]; Sharpr-MPRA). Nothing here is clinically validated.

`torch` is the one heavy dependency; CPU-only wheels are sufficient (`pip install -r
requirements.txt`). `python -c "import regmodel"`, `python -m regmodel.cli demo`, and `pytest`
all run fully offline.

## `trio_prioritizer` — inheritance-aware prioritization of non-coding variants in trios

A third self-contained package (`trio_prioritizer/`). Where `variant_curator` scores one
variant and `regmodel` reads a per-base effect map, `trio_prioritizer` answers the trio
question a clinical geneticist actually asks (per Matt Deardorff): *given an affected child
and two parents, which intronic/regulatory variants look causal, and how do they rank
against everything else on the table?* Two inheritance patterns dominate rare disease:

- **De novo** — present in the child, absent in both parents (`0/1` vs `0/0`,`0/0`): the
  usual mechanism for severe dominant disease. With WGS these now land in introns, UTRs,
  and regulatory regions, not just exons.
- **Compound heterozygous** — the child is affected and each parent an unaffected carrier.
  You often find one (coding) hit but *cannot find the second*, and increasingly that
  second allele is a deep-intronic / regulatory variant creating a cryptic splice site or
  disrupting an enhancer. Finding and scoring that second hit is the pain point.

This package supplies the **new** logic — trio genotype → inheritance mode, compound-het
pairing by parent-of-origin (it pairs a coding allele from one parent with a non-coding
allele from the other), and an integrated ranking — on top of the **existing** non-coding
scorer. It reuses `variant_curator.clients.alphagenome.fetch_alphagenome` as the scoring
engine (coordinate-based, which is what trio inputs are) via an **injectable** `scorer`, so
tests and the demo never touch the network. Scores blend a genotype-derived inheritance
prior with the AlphaGenome magnitude (`combined = 0.6·prior + 0.4·magnitude`; weights are
provisional constants needing calibration). When AlphaGenome is unavailable the score is
marked unavailable with a reason and ranking still runs on inheritance — it never fabricates
a number. `regmodel` ISM is an optional secondary signal, used only when a local sequence
window is supplied (so the core path needs no torch).

### Offline demo (no network, no API key)

```
python -m trio_prioritizer.cli demo            # synthetic trio -> ranked candidates + candidates.json
python -m trio_prioritizer.cli run --table trio.tsv   # rank a real TSV (live AlphaGenome if keyed)
```

The `demo` builds a **seeded synthetic trio** with a planted answer — a de novo deep-intronic
variant in *SCN1A*, a compound-het pair in *PAH* (a paternal coding missense + a maternal
deep-intronic "second hit"), and benign inherited distractors — and asserts the causal de
novo and the comp-het pair land in the top two, using a deterministic offline scorer whose
`Source` is explicitly labelled **ILLUSTRATIVE**. `run` reads a TSV
(`chrom,pos,ref,alt,gene,consequence,proband_gt,mother_gt,father_gt`) and uses the live
AlphaGenome scorer when `ALPHAGENOME_API_KEY` is set, otherwise degrades with a clear note.

Both commands write `candidates.json` with full provenance and a research-use caveat on
every scored claim. Scope for the MVP is autosomal de novo + recessive (homozygous and
compound het); X-linked/imprinting-aware inheritance is a documented future extension.
`python -c "import trio_prioritizer"`, the demo, and `pytest` all run fully offline with no
key and no torch. **Research use only** — decision-support, not a diagnosis.

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

regmodel/
  encoding.py         one-hot DNA (4×L, ACGT order), decode, reverse-complement
  data.py             SyntheticMPRA (offline default) + network-gated real-loader stub
  model.py            compact Conv1d CNN (motif detectors -> global pool -> MLP -> scalar)
  train.py            seeded split/train/early-stop; held-out Pearson/Spearman/MSE; sidecar
  ism.py              (L×4) ISM delta matrix, importance track, single-variant delta
  compare.py          AlphaGenome cross-check reusing variant_curator's client (skips w/o key)
  plots.py            matplotlib (Agg) ISM heatmap, importance track, cross-check scatter
  cli.py              train / ism / variant / compare / demo (offline)

trio_prioritizer/
  models.py           Genotype/TrioVariant/InheritanceCall/NoncodingScore/PrioritizedCandidate
                      (reuses variant_curator's Source + VariantConsequence)
  inheritance.py      trio genotypes -> inheritance mode; compound-het pairing by parent-of-origin
  scoring.py          non-coding score via injectable AlphaGenome scorer; optional regmodel ISM hook
  prioritize.py       blends inheritance prior + non-coding magnitude into one ranked list
  data.py             seeded synthetic trio (planted de novo + comp-het) + offline fake scorer
  cli.py              demo (offline synthetic) / run --table trio.tsv; writes candidates.json
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
