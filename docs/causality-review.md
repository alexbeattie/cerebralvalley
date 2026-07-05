# Causality Review (Step 2) — does this variant actually explain THIS patient?

## The reframe

Diagnostic variant work has three stages. Matt's real ask is the **second** one.

1. **Lab — variant classification.** Prioritize candidate variants off relatively
   thin clinical detail. (This is what `variant_curator` on `main` does.)
2. **Clinical care — causality review.** The lab hands the clinician 2–5 reported
   variants, each somewhat suspicious. The clinician holds what the lab did not:
   the patient's full phenotype. For each reported variant, *does it actually
   explain THIS patient?*
3. **Hard frontier — intronic/regulatory variants** the exome missed. (Separate
   experiments on the other feature branches.)

This tool is **step 2**.

## Matching, run in reverse

- **Lab (typical):** find variants that *might* explain a patient.
- **Clinician (this tool):** does *this* variant actually explain *this* patient?

Same knowledge bases, opposite direction — and the clinician supplies the input
the lab never had: the fully worked-up phenotype.

## The two traps this exists to counter

1. **The mind overfits.** It is easy to talk yourself into a match when the
   overlap is only partial. We score against an explicit HPO feature set, so a
   partial match reads as *"Partial (explains 3/5)"*, not *"close enough"*.
2. **Partial explanations.** A variant may account for most features but not all.
   Leftover features are either (a) an unreported extension of a listed gene's
   disease, or (b) a **second independent cause** (~5% of solved cases). Because
   the lab may never have prioritized a variant for those features, the report
   can contain *no* candidate for them — which is the trigger to **re-analyze the
   genome**, not a failure of the review.

## What it does today (working MVP)

Inputs:
- **Reported variants** — `GENE` or `GENE:hgvs` (the lab's call is taken as given).
- **Patient profile** — HPO terms, entered as `HP:xxxxxxx` or free text that the
  tool resolves against the HPO ontology (e.g. `"seizures"` → `HP:0001250`).

Engine (`causality_review/engine.py`):
1. For each reported variant, pull its gene's known disease phenotypes from HPO/Jax.
2. Split the patient's features into **explained** vs **unexplained** by that gene.
   Matching is **ontology-aware**: a gene annotated to a broad term (e.g. *Seizure*)
   explains a patient's more specific feature (e.g. *Focal-onset seizure*) by walking
   the HPO `is_a` graph. Such matches are shown as "via broader: <term>".
3. Score = fraction of the patient's features the gene explains; map to a tier:
   `Best fit (****) ≥ all · Possible (***) ≥ 0.6 · Partial (**) ≥ 0.3 · Weak (*) > 0 · Unlikely (.)`.
4. Rank variants by fit (objective, not by the lab's ordering).
5. Compute **residual** features explained by *no* reported variant → second-cause /
   re-analysis flag.
6. Detect a **dual diagnosis**: when no single variant covers the picture but two
   reported variants are complementary (each explains features the other does not)
   and together account for it, flag the ~5% "two independent causes" case.

Every gene–phenotype claim carries an openable HPO source URL, and the tool
abstains (marks a variant unscored) when it cannot retrieve knowledge for a gene —
the same "cite everything, don't assert what the data can't support" discipline as
step 1.

Run it:

```
# web UI (Matt's sketch): enter variants + features, get the ranked list back
python -m causality_review.webapp                   # http://localhost:8000

# or the CLI
python -m causality_review.cli                      # hardcoded demo case
python -m causality_review.cli \
    --variant SCN1A:c.3637C>T --variant MYH7:c.1063G>A \
    --hpo "seizures" --hpo "global developmental delay" --hpo "ataxia"
```

The web UI (`causality_review/webapp.py` + `static/index.html`) is a stdlib HTTP server
with a single JSON endpoint that runs the real `review_causality` — same engine as the CLI.
It covers steps 2 and 3 of Matt's interface sketch (enter reported variants, enter/paste
patient features, get the ranked causality list with explained/unexplained features, source
links, and flags). Step 1 of his sketch — **drag-and-drop the lab-report PDF** — is still on
the "not built yet" list below; today variants and features are entered as text.

## Evaluation

`causality_review/eval.py` is a ranking sanity/regression harness: 9 solved cases
(each a set of real classic HPO features + the true causal gene) scored against a
panel of every fixture gene, measuring how often the true gene ranks first.

```
python -m causality_review.eval
```

Current result: **top-1 89%, top-3 100%, MRR 0.926** (9 cases, 9-gene panel).

Honest caveat: this tests ranking and discrimination against live HPO data, not a
held-out clinical validation — the phenotypes are curated classic features and the
tool matches against the same HPO annotations, so a real deployment must still be
measured on real solved cases. The single miss (RYR1 malignant hyperthermia, 3rd) is
genuine signal: MH is an anesthesia-triggered reaction whose baseline HPO features
overlap other myopathies.

Deriving the eval also drove a real engine fix: naive ancestor matching let genes
"explain" a feature through organ-system container nodes (e.g. FBN1 matching *Prolonged
QT interval* via *Abnormality of the cardiovascular system*), inflating cross-system
scores. Those container terms are now excluded, and exact matches tie-break above
ancestor matches. Top-1 went 56% → 89%.

Offline unit tests (no network, fake ontology) cover matching, tiering, the
explained/unexplained split, ranking, residual, and the dual-diagnosis flag:

```
python -m unittest discover -s tests
```

## Data sources

- **HPO / Jax ontology API** (`ontology.jax.org`) — gene → diseases → phenotype
  ids, term ancestors (is_a graph), and free-text → HPO term resolution. Public, keyless.

## AI at the edge (optional)

Free-text clinical notes → HPO terms is the one place an LLM earns its keep. When
`ANTHROPIC_API_KEY` is set, the UI's "Clinical notes" box and the CLI's `--notes` /
`--notes-file` extract phenotypes via Claude (`causality_review/clients/llm.py`,
`extract.py`). The design keeps AI strictly at the ingestion boundary:

1. the LLM **proposes** short phenotype phrases from the note;
2. each phrase is **grounded** in a real HPO term by the deterministic ontology search —
   the model never emits an HPO id, so it can't invent one;
3. the scoring engine only ever sees validated HPO terms, and stays AI-free and sourced.

Ungrounded phrases are reported, not silently dropped, and the whole feature degrades
gracefully (panel hidden / clear message) when no key is set. No new dependency — the
Anthropic Messages API is called over the existing httpx layer; the key is read from the
environment and never logged.

## Not built yet (next)

- **PDF lab-report ingestion** — step 1 of the sketch (drag-and-drop the lab PDF to
  auto-extract the reported variants). Notes → HPO extraction now exists; the PDF/variant
  side does not, so variants are still entered as text.
- **Penetrance / inheritance / zygosity weighting**, and a management layer
  (natural history, what to test next, referrals).
- **Information-content weighting** — a rare, specific feature (e.g. *Ectopia lentis*)
  should count more than a common one (*Developmental delay*); current scoring weights
  every feature equally.

## Two hard rules (unchanged)

1. **No PHI.** Phenotypes are ontology ids; gene knowledge is public. The engine
   can be pointed at real EMR-derived HPO profiles behind a firewall without code
   changes.
2. **Guard scope.** A narrow slice that works and is measured beats a broad
   platform that half-runs.
