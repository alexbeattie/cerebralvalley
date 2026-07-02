"""AlphaGenome client (Google DeepMind).

Deep intronic and regulatory variants are exactly what a consequence-based tool
rejects: they are not missense/nonsense, yet they can be disease-causing by
creating cryptic splice sites, disrupting branch points, or altering transcription,
accessibility, or expression. AlphaGenome scores those effects directly, so this is
the evidence source that lets the curator reason about the non-coding path at all.

Design notes that keep the rest of the project honest:
  - The API needs a key (`ALPHAGENOME_API_KEY`, non-commercial research use). If it
    is unset we return `found=False` with a reason rather than crashing, so imports
    and the whole offline pipeline still work.
  - The SDK is imported lazily *inside* `_score`, so `import variant_curator` never
    depends on `alphagenome` being installed. Any ImportError / API error degrades
    to `found=False` with a reason. We never fabricate a score.
  - Predictions are research-model output, not clinically validated: they map to
    supporting ACMG evidence only (see `variant_curator.acmg`).
"""

from __future__ import annotations

import os

import pandas as pd

from ..http import now_iso
from ..models import AlphaGenomeEvidence, RegulatorySignal, Source, SplicingSignal

API_KEY_ENV = "ALPHAGENOME_API_KEY"
ALPHAGENOME_URL = "https://deepmind.google.com/science/alphagenome"

# Scorer keys into variant_scorers.RECOMMENDED_VARIANT_SCORERS. Splicing is the
# primary signal for deep intronic variants; the regulatory set covers
# promoter / enhancer / UTR effects.
SPLICING_SCORERS = ("SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS")
REGULATORY_SCORERS = ("RNA_SEQ", "CAGE", "PROCAP", "DNASE", "ATAC", "CHIP_TF", "CHIP_HISTONE")

_SPLICING_OUTPUT_TYPES = set(SPLICING_SCORERS)
_REGULATORY_OUTPUT_TYPES = set(REGULATORY_SCORERS)

# Human phrasing per splicing scorer, so the interpretation reads like a note a
# curator would write rather than a column name.
_SPLICING_PHRASE = {
    "SPLICE_SITES": "predicted splice-site alteration",
    "SPLICE_SITE_USAGE": "predicted change in splice-site usage",
    "SPLICE_JUNCTIONS": "predicted altered splice junction",
}

# Heuristic disease-relevant tissue per gene so regulatory scoring can be reported
# tissue-first for the gene's disease context. Deliberately small and provisional;
# when a gene is absent we fall back to the top tissue across all tracks. Curies are
# UBERON terms for the primary affected tissue.
GENE_TISSUE_HINT: dict[str, tuple[str, str]] = {
    "LDLR": ("liver", "UBERON:0002107"),
    "APOB": ("liver", "UBERON:0002107"),
    "PCSK9": ("liver", "UBERON:0002107"),
    "PAH": ("liver", "UBERON:0002107"),
    "GAA": ("skeletal muscle", "UBERON:0001134"),
    "MYH7": ("heart", "UBERON:0000948"),
    "MYBPC3": ("heart", "UBERON:0000948"),
    "KCNQ1": ("heart", "UBERON:0000948"),
    "KCNH2": ("heart", "UBERON:0000948"),
    "SCN5A": ("heart", "UBERON:0000948"),
    "LMNA": ("heart", "UBERON:0000948"),
    "RYR1": ("skeletal muscle", "UBERON:0001134"),
    "SCN1A": ("brain", "UBERON:0000955"),
    "MECP2": ("brain", "UBERON:0000955"),
    "CFTR": ("lung", "UBERON:0002048"),
}


def _chrom(chrom: str) -> str:
    """AlphaGenome wants a 'chr'-prefixed contig; VEP hands us the bare name."""
    c = str(chrom)
    return c if c.startswith("chr") else f"chr{c}"


def _score(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    scorer_keys: tuple[str, ...],
    *,
    api_key: str,
) -> pd.DataFrame:
    """Run one scorer set against AlphaGenome and return the tidy long DataFrame.

    All SDK contact lives here so tests can patch this single function and run with
    no key and no package installed.
    """
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers

    model = dna_client.create(api_key)
    variant = genome.Variant(
        chromosome=_chrom(chrom),
        position=int(pos),
        reference_bases=ref,
        alternate_bases=alt,
    )
    interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    selected = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[k] for k in scorer_keys]
    scores = model.score_variant(
        interval=interval,
        variant=variant,
        variant_scorers=selected,
        organism=dna_client.Organism.HOMO_SAPIENS,
    )
    return variant_scorers.tidy_scores([scores])


def _for_gene(df: pd.DataFrame, gene_symbol: str) -> pd.DataFrame:
    """Splicing scores are returned per overlapping gene; keep only ours."""
    if df is None or df.empty or "gene_name" not in df.columns:
        return df
    return df[df["gene_name"] == gene_symbol]


def _regulatory_rows(df: pd.DataFrame, gene_symbol: str) -> pd.DataFrame:
    """Keep gene-scoped rows for our gene plus positional tracks (accessibility /
    ChIP) that carry no gene, which a strict gene filter would wrongly drop."""
    if df is None or df.empty or "gene_name" not in df.columns:
        return df
    names = df["gene_name"]
    mask = (names == gene_symbol) | names.isna() | (names.astype(str).str.strip() == "")
    return df[mask]


def _num(value) -> float | None:
    if value is None:
        return None
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return None if pd.isna(f) else f


def _summarize_splicing(df: pd.DataFrame, gene_symbol: str, source: Source) -> list[SplicingSignal]:
    df = _for_gene(df, gene_symbol)
    if df is None or df.empty:
        return []
    signals: list[SplicingSignal] = []
    for output_type, group in df.groupby("output_type"):
        if output_type not in _SPLICING_OUTPUT_TYPES:
            continue
        # Non-directional: rank by |quantile| (calibrated) then |raw|.
        ranked = group.assign(_mag=group.apply(_row_magnitude, axis=1))
        top = ranked.sort_values("_mag", ascending=False).iloc[0]
        q = _num(top.get("quantile_score"))
        r = _num(top.get("raw_score"))
        phrase = _SPLICING_PHRASE.get(str(output_type), "predicted splicing effect")
        mag = abs(q) if q is not None else (abs(r) if r is not None else None)
        mag_str = f"|quantile|={abs(q):.2f}" if q is not None else (f"|raw|={abs(r):.2f}" if r is not None else "no score")
        signals.append(
            SplicingSignal(
                scorer=str(top.get("variant_scorer", output_type)),
                output_type=str(output_type),
                top_gene=str(top.get("gene_name", gene_symbol)),
                max_quantile=abs(q) if q is not None else None,
                max_raw=abs(r) if r is not None else None,
                direction=None,
                interpretation=f"{phrase} in {gene_symbol} ({mag_str})",
                source=source,
            )
        )
    return sorted(signals, key=lambda s: (s.magnitude or 0.0), reverse=True)


def _summarize_regulatory(
    df: pd.DataFrame, gene_symbol: str, source: Source
) -> list[RegulatorySignal]:
    df = _regulatory_rows(df, gene_symbol)
    if df is None or df.empty:
        return []
    hint = GENE_TISSUE_HINT.get(gene_symbol.upper())
    signals: list[RegulatorySignal] = []
    for output_type, group in df.groupby("output_type"):
        if output_type not in _REGULATORY_OUTPUT_TYPES:
            continue
        top = _pick_regulatory_row(group, hint)
        q = _num(top.get("quantile_score"))
        r = _num(top.get("raw_score"))
        tissue = _tissue_of(top)
        signed = q if q is not None else r
        direction = None if signed in (None, 0) else ("up" if signed > 0 else "down")
        mag_str = f"quantile={q:+.2f}" if q is not None else (f"raw={r:+.2f}" if r is not None else "no score")
        signals.append(
            RegulatorySignal(
                modality=str(output_type),
                output_type=str(output_type),
                top_tissue=tissue,
                ontology_curie=_str_or_none(top.get("ontology_curie")),
                quantile_score=q,
                raw_score=r,
                direction=direction,
                interpretation=f"{output_type} {direction or 'change'} in {tissue} ({mag_str})",
                source=source,
            )
        )
    return sorted(signals, key=lambda s: (s.magnitude or 0.0), reverse=True)


def _row_magnitude(row) -> float:
    q = _num(row.get("quantile_score"))
    if q is not None:
        return abs(q)
    r = _num(row.get("raw_score"))
    return abs(r) if r is not None else 0.0


def _pick_regulatory_row(group: pd.DataFrame, hint: tuple[str, str] | None):
    """Prefer the disease-relevant tissue when the hint matches a track; else the
    strongest track across all tissues."""
    if hint is not None:
        tissue_term, curie = hint
        for _, row in group.iterrows():
            tissue = _tissue_of(row).lower()
            row_curie = (_str_or_none(row.get("ontology_curie")) or "").upper()
            if tissue_term.lower() in tissue or curie.upper() == row_curie:
                return row
    ranked = group.assign(_mag=group.apply(_row_magnitude, axis=1))
    return ranked.sort_values("_mag", ascending=False).iloc[0]


def _tissue_of(row) -> str:
    for key in ("gtex_tissue", "biosample_name", "track_name"):
        val = _str_or_none(row.get(key))
        if val:
            return val
    return "unspecified tissue"


def _str_or_none(value) -> str | None:
    if value is None:
        return None
    if isinstance(value, float) and pd.isna(value):
        return None
    text = str(value).strip()
    return text or None


def fetch_alphagenome(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    gene_symbol: str,
    *,
    api_key: str | None = None,
) -> AlphaGenomeEvidence:
    """Score a GRCh38 variant for splicing + regulatory effects on `gene_symbol`.

    Takes the same coordinates VEP already resolves. Reads the API key from
    `ALPHAGENOME_API_KEY` when not passed. Degrades to `found=False` + a reason on a
    missing key, a missing package, or any API error — never a fabricated score.
    """
    key = api_key if api_key is not None else os.environ.get(API_KEY_ENV)
    variant_id = f"{_chrom(chrom)}:{pos}:{ref}>{alt}"
    source = Source(
        name="AlphaGenome (Google DeepMind)",
        url=ALPHAGENOME_URL,
        retrieved_at=now_iso(),
        detail=variant_id,
    )

    if not key:
        return AlphaGenomeEvidence(
            found=False,
            variant_id=variant_id,
            reason=f"{API_KEY_ENV} not set",
            source=source,
        )

    try:
        splice_df = _score(chrom, pos, ref, alt, SPLICING_SCORERS, api_key=key)
        reg_df = _score(chrom, pos, ref, alt, REGULATORY_SCORERS, api_key=key)
    except ImportError:
        return AlphaGenomeEvidence(
            found=False,
            variant_id=variant_id,
            reason="alphagenome package not installed (pip install alphagenome)",
            source=source,
        )
    except Exception as exc:  # any API / network / SDK failure: abstain, don't guess
        return AlphaGenomeEvidence(
            found=False,
            variant_id=variant_id,
            reason=f"AlphaGenome API error: {exc}",
            source=source,
        )

    splicing = _summarize_splicing(splice_df, gene_symbol, source)
    regulatory = _summarize_regulatory(reg_df, gene_symbol, source)

    return AlphaGenomeEvidence(
        found=True,
        variant_id=variant_id,
        api="alphagenome dna_client",
        splicing=splicing,
        regulatory=regulatory,
        source=source,
    )
