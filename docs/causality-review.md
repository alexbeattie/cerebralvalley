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
3. Score = fraction of the patient's features the gene explains; map to a tier:
   `Best fit (****) ≥ all · Possible (***) ≥ 0.6 · Partial (**) ≥ 0.3 · Weak (*) > 0 · Unlikely (.)`.
4. Rank variants by fit (objective, not by the lab's ordering).
5. Compute **residual** features explained by *no* reported variant → second-cause /
   re-analysis flag.

Every gene–phenotype claim carries an openable HPO source URL, and the tool
abstains (marks a variant unscored) when it cannot retrieve knowledge for a gene —
the same "cite everything, don't assert what the data can't support" discipline as
step 1.

Run it:

```
python -m causality_review.cli                      # hardcoded demo case
python -m causality_review.cli \
    --variant SCN1A:c.3637C>T --variant MYH7:c.1063G>A \
    --hpo "seizures" --hpo "global developmental delay" --hpo "ataxia"
```

## Data sources

- **HPO / Jax ontology API** (`ontology.jax.org`) — gene → diseases → phenotype
  ids, and free-text → HPO term resolution. Public, keyless.

## Not built yet (next)

- **PDF lab-report ingestion** — today variants are passed on the CLI; step 2's
  sketch has the clinician drag-and-drop the lab PDF.
- **EMR-note → HPO extraction** — today the curator supplies HPO terms; production
  extracts and normalizes them from clinical notes behind the site's firewall.
- **Ontology-aware matching** — currently exact HPO-id membership. Should credit a
  patient's specific term (e.g. *focal-onset seizure*) when the gene is annotated
  to an ancestor (*seizure*), via HPO's `is_a` graph.
- **Penetrance / inheritance / zygosity weighting**, and a management layer
  (natural history, what to test next, referrals).
- **Eval harness** — solved cases with known causal genes, labels hidden, measure
  how often the true gene ranks first.

## Two hard rules (unchanged)

1. **No PHI.** Phenotypes are ontology ids; gene knowledge is public. The engine
   can be pointed at real EMR-derived HPO profiles behind a firewall without code
   changes.
2. **Guard scope.** A narrow slice that works and is measured beats a broad
   platform that half-runs.
