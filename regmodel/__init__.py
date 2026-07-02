"""regmodel: a miniature MPRA-trained sequence -> activity model with single-base ISM.

What this is: a compact, CPU-trainable CNN that learns to map short DNA sequences to a
scalar regulatory activity (the massively-parallel-reporter-assay setting), plus
in-silico mutagenesis (ISM) that mutates every base to read off a per-position
mutation-effect map. This reproduces, in miniature, the pattern Katie Pollard's lab uses
(PARM, SuPreMo/Akita): train a sequence model on an MPRA, then ISM to ask what one base
change does. The differentiated step is a head-to-head **cross-check against AlphaGenome**
(a 1 Mb genome foundation model) on the same variants -- where does a small task-specific
model agree or disagree with a foundation model?

Offline-first: the whole train -> evaluate -> ISM -> variant pipeline runs on a CPU laptop
with **no network and no API key** using a built-in synthetic MPRA. The AlphaGenome
cross-check degrades cleanly (skips with a reason) when `ALPHAGENOME_API_KEY` is unset or
the SDK is unavailable -- it never fabricates a score.

Research use only: by default this is a toy model trained on *synthetic* data with a
planted motif. Real conclusions require loading a real MPRA (see `regmodel.data`). Nothing
here is clinically validated.
"""

from __future__ import annotations

__version__ = "0.1.0"
