"""Shared tidy_scores-shaped fixtures for the AlphaGenome tests.

These mirror the documented columns of `variant_scorers.tidy_scores()` so the
client's reduction can be exercised offline with no SDK and no API key.
"""

from __future__ import annotations

import pandas as pd

COLUMNS = [
    "variant_id",
    "gene_name",
    "gene_id",
    "output_type",
    "variant_scorer",
    "track_name",
    "ontology_curie",
    "gtex_tissue",
    "biosample_name",
    "raw_score",
    "quantile_score",
]


def _row(**over):
    base = {c: None for c in COLUMNS}
    base.update(variant_id="chr7:117559593:C:T", gene_id="ENSG00000001626")
    base.update(over)
    return base


def splicing_df(gene: str = "CFTR") -> pd.DataFrame:
    rows = [
        # Strong splice-site alteration on our gene: this is the top effect.
        _row(gene_name=gene, output_type="SPLICE_SITES", variant_scorer="SPLICE_SITES",
             raw_score=3.1, quantile_score=0.95),
        # A weaker SPLICE_SITES track for the same gene — must be dominated.
        _row(gene_name=gene, output_type="SPLICE_SITES", variant_scorer="SPLICE_SITES",
             raw_score=1.0, quantile_score=0.40),
        _row(gene_name=gene, output_type="SPLICE_JUNCTIONS", variant_scorer="SPLICE_JUNCTIONS",
             raw_score=2.2, quantile_score=0.88),
        _row(gene_name=gene, output_type="SPLICE_SITE_USAGE", variant_scorer="SPLICE_SITE_USAGE",
             raw_score=0.5, quantile_score=0.30),
        # A neighbouring gene in the 1MB window — must be filtered out.
        _row(gene_name="NEIGHBOR", output_type="SPLICE_SITES", variant_scorer="SPLICE_SITES",
             raw_score=9.9, quantile_score=0.99),
    ]
    return pd.DataFrame(rows, columns=COLUMNS)


def regulatory_df(gene: str = "CFTR") -> pd.DataFrame:
    rows = [
        # RNA_SEQ in a non-disease tissue with a larger score...
        _row(gene_name=gene, output_type="RNA_SEQ", variant_scorer="RNA_SEQ",
             gtex_tissue="Liver", ontology_curie="UBERON:0002107",
             raw_score=2.0, quantile_score=0.90),
        # ...and in lung, the CFTR disease-relevant tissue: hint must prefer this.
        _row(gene_name=gene, output_type="RNA_SEQ", variant_scorer="RNA_SEQ",
             gtex_tissue="Lung", ontology_curie="UBERON:0002048",
             raw_score=-2.0, quantile_score=-0.70),
        # Positional accessibility track: no gene_name, must NOT be filtered out.
        _row(output_type="DNASE", variant_scorer="DNASE", biosample_name="lung fibroblast",
             ontology_curie="UBERON:0002048", raw_score=1.1, quantile_score=0.50),
    ]
    return pd.DataFrame(rows, columns=COLUMNS)


def all_low_regulatory_df(gene: str = "CFTR") -> pd.DataFrame:
    rows = [
        _row(gene_name=gene, output_type="RNA_SEQ", variant_scorer="RNA_SEQ",
             gtex_tissue="Lung", ontology_curie="UBERON:0002048",
             raw_score=0.1, quantile_score=0.05),
    ]
    return pd.DataFrame(rows, columns=COLUMNS)


def all_low_splicing_df(gene: str = "CFTR") -> pd.DataFrame:
    rows = [
        _row(gene_name=gene, output_type="SPLICE_SITES", variant_scorer="SPLICE_SITES",
             raw_score=0.2, quantile_score=0.10),
    ]
    return pd.DataFrame(rows, columns=COLUMNS)
