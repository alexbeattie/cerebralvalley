# Build brief: `regmodel` — an MPRA-trained sequence→activity model with single-base ISM, cross-checked against AlphaGenome

You (Claude Code) are building a **new, self-contained Python package `regmodel/`** in this
worktree, alongside the existing `variant_curator/` package. This is a hackathon "Research
track" project for a Gladstone Institutes event (Claude Science). Read this brief fully, look
at the existing repo for style, then build it end to end.

## The scientific idea (one paragraph)

Katie Pollard's lab trains deep-learning models on **massively parallel reporter assays
(MPRAs)** that map short DNA sequences to regulatory activity, then uses **in-silico
mutagenesis (ISM)** to ask what a single-base change does to predicted activity (this is how
they study human accelerated regions and disease variants; see PARM, Nature 2025, and
SuPreMo/Akita). We reproduce that pattern in miniature: train a compact sequence→activity
CNN on an MPRA, run single-base ISM to produce a per-base mutation-effect map, and then
**cross-check our lightweight model against AlphaGenome** (a 1 Mb foundation model) on the
same variants. The head-to-head — "where does a small task-specific MPRA model agree or
disagree with a genome foundation model?" — is the differentiated result.

## Non-negotiable engineering constraints (this is a hackathon; it must RUN)

1. **Offline-first, no GPU, no downloads required to demo or test.** The entire pipeline
   (train → evaluate → ISM → variant query) MUST run to completion on a CPU laptop with **no
   network access** using a built-in **synthetic MPRA dataset** (described below). Real-data
   loaders are optional and network-gated.
2. **Deterministic where it matters.** Seed numpy/torch. The synthetic task must be learnable
   well enough that tests can assert the model recovers the planted signal.
3. **Graceful degradation for AlphaGenome.** Reuse the existing
   `variant_curator/clients/alphagenome.py` (already on this branch). If `ALPHAGENOME_API_KEY`
   is unset or the package/API is unavailable, the cross-check step must skip cleanly with a
   message — never crash, never fabricate.
4. **Provenance & reproducibility**, matching the repo ethos: every saved artifact (model,
   metrics, ISM map, comparison) is written with a small JSON sidecar recording the config,
   data source, seed, timestamp, and (for AlphaGenome) the `Source`. A curator/judge must be
   able to see exactly how each number was produced.
5. **Keep it small and honest.** A compact CNN a laptop can train in a couple of minutes.
   Don't over-engineer. No half-finished features. Comments explain *why*, not *what*.

## Repository context to match

- `variant_curator/` is the existing package: dataclasses with a `Source` on every evidence
  item, module docstrings that explain *why*, `from __future__ import annotations`, small
  focused functions, graceful degradation, abstain-don't-guess. Match this style.
- Reuse `variant_curator.clients.alphagenome.fetch_alphagenome(chrom, pos, ref, alt,
  gene_symbol)` and `variant_curator.http` helpers where useful. Do not modify
  `variant_curator/` except tiny, clearly-justified changes if truly needed.
- `requirements.txt` currently has `httpx`, `alphagenome`, `pandas`. Add what you need.

## Package layout to build

```
regmodel/
  __init__.py          scope statement (what this is, offline-first, research use only)
  encoding.py          one-hot DNA encoding/decoding; reverse-complement; helpers
  data.py              dataset abstraction: SyntheticMPRA (default, offline) + a network-
                       gated real MPRA loader stub with a clear docstring naming candidate
                       public datasets (lentiMPRA developing human brain [Pollard/Ahituv],
                       Sharpr-MPRA). Returns (sequences, activities) + a provenance dict.
  model.py             a compact CNN (Basset/DeepSTARR-like: 2-3 conv blocks -> global pool
                       -> MLP head) in PyTorch for scalar activity regression. Small.
  train.py             training loop: split, train, early-ish stop, return metrics
                       (Pearson/Spearman r, MSE) on held-out test set; save model + sidecar.
  ism.py               in-silico mutagenesis: for a sequence, mutate every position to all 3
                       alternative bases, predict, return a (L x 4) delta matrix and a
                       per-position importance track; a single-variant query (ref seq, 0-based
                       pos, alt base) -> predicted activity delta.
  plots.py             matplotlib: ISM heatmap / per-base importance track; a scatter of
                       our-model vs AlphaGenome for the cross-check. Save PNGs.
  compare.py           cross-check: given variants (genomic coords + a local sequence window),
                       get AlphaGenome regulatory score via variant_curator's client and our
                       model's ISM delta, align them, compute correlation; degrade gracefully
                       with no key.
  cli.py               subcommands: `train`, `ism`, `variant`, `compare`, plus a `demo` that
                       runs the whole synthetic pipeline offline and writes artifacts.
tests/
  test_encoding.py
  test_synthetic_data.py
  test_model_learns_planted_motif.py   (trains briefly on synthetic; asserts test r above a
                                        modest threshold AND ISM importance peaks at the
                                        planted motif location)
  test_ism.py                          (variant query sign/shape; ISM matrix shape = L x 4)
  test_compare_offline.py              (compare skips cleanly with no ALPHAGENOME_API_KEY)
artifacts/            (gitignored) where demo writes model + plots + sidecars
```

## The synthetic MPRA task (make it genuinely learnable and interpretable)

- Generate N random sequences of fixed length L (e.g. N=4000, L=200) over {A,C,G,T}.
- Define a short "activator" motif (e.g. `TGACTCA`, an AP-1-like site) and optionally a
  weaker repressor motif. Activity = base noise + (weight if activator motif present, scaled
  by count/position) − (repressor effect) + Gaussian noise. Plant motifs in a random subset
  of sequences so the label correlates with motif presence.
- This gives a ground truth: a well-trained model's ISM should light up the planted motif
  positions, which the test can assert. Keep the generator seeded and documented.

## Model / training specifics

- PyTorch, CPU-friendly. One-hot input shape (batch, 4, L). 2–3 Conv1d blocks (e.g. 64→128
  filters, kernel ~7–11, ReLU, maxpool), global average+max pool, small MLP to a scalar.
- Regression to activity. Adam, MSE loss, a modest number of epochs with a fixed seed;
  training on the synthetic set should finish in well under a minute on CPU for the test
  (allow a `--fast`/small-epoch path used by tests).
- Report Pearson and Spearman on a held-out split. Save `model.pt` + `model.json` sidecar
  (arch, hyperparams, data provenance, seed, metrics, timestamp).

## ISM specifics

- `ism_matrix(model, seq) -> np.ndarray` of shape (L, 4): predicted activity for each
  single-base substitution minus the reference prediction (ref base entries = 0).
- `importance_track(ism_matrix)` -> per-position scalar (e.g. mean or min delta) for plotting.
- `variant_effect(model, ref_seq, pos, alt) -> float`: the predicted delta for one change.
- A helper to locate the strongest ISM positions so the test can compare against the planted
  motif coordinates.

## AlphaGenome cross-check specifics

- In `compare.py`, accept a small list of variants as (chrom, pos, ref, alt, gene_symbol,
  local_window_seq). For each: call `fetch_alphagenome(...)` (reuse the branch's client) to
  get a regulatory magnitude, and compute our model's `variant_effect` on the local window.
- Produce a tidy table (our_delta, alphagenome_quantile/raw, agreement sign) + a scatter PNG
  and a correlation. If AlphaGenome is unavailable (no key/package/API), write the table with
  our-model values only and a clear note that the external comparison was skipped, and why.
- Every AlphaGenome-derived row carries the `Source` from the client. Never fabricate.

## CLI / definition of done

- `python -m regmodel.cli demo` runs fully offline: trains on synthetic data, prints metrics,
  runs ISM on a held-out sequence, writes `artifacts/` (model + sidecar, ISM heatmap PNG,
  importance track PNG, a metrics.json). Prints where everything was written.
- `python -m regmodel.cli ism --seq ACGT...` prints/saves an ISM map for a given sequence.
- `python -m regmodel.cli variant --seq ACGT... --pos 100 --alt G` prints the predicted delta.
- `python -m regmodel.cli compare ...` runs the AlphaGenome cross-check (skips cleanly w/o key).
- `python -c "import regmodel"` works with no network.
- `pytest` passes offline with no GPU and no API key.
- Update `requirements.txt` (torch CPU, numpy, matplotlib, pandas as needed; note in a
  comment that torch is the heavy dep and CPU wheels are fine).
- Add a `regmodel` section to `README.md`: the idea, the offline `demo`, the Pollard/MPRA
  framing, the AlphaGenome cross-check, the `ALPHAGENOME_API_KEY` note, and a **research-use-
  only** caveat (this is a toy model on synthetic data by default; real MPRA data must be
  loaded for real conclusions).
- Commit on the current branch `feat/mpra-regulatory-model` with a clear message. Do NOT push.

## Hard rules (carry over from the project)

- Provenance on every artifact; abstain/skip rather than fabricate; research-use-only caveats
  where results are shown; match the existing code style; do not break `variant_curator` or
  its tests.

When done, print a concise summary of files created and exactly how to run the offline demo.
