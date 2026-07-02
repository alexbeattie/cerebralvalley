"""Allowlist of well-characterized disease genes for the MVP.

We deliberately restrict to a small set so that behavior is predictable and
reviewable. Each entry carries the MANE Select / canonical RefSeq transcript we
resolve HGVS coding changes against, so consequence calls are reproducible.

Restricting scope is a hard rule: a narrow slice that works and is measured beats
a broad platform that half-runs.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class GeneSpec:
    symbol: str
    refseq_transcript: str  # MANE Select where available
    ncbi_gene_id: int
    disease_context: str


# 30 well-characterized disease genes spanning cancer predisposition, cardiac,
# metabolic, and neurodevelopmental contexts. Transcripts are MANE Select.
GENE_ALLOWLIST: dict[str, GeneSpec] = {
    g.symbol: g
    for g in [
        GeneSpec("BRCA1", "NM_007294.4", 672, "Hereditary breast/ovarian cancer"),
        GeneSpec("BRCA2", "NM_000059.4", 675, "Hereditary breast/ovarian cancer"),
        GeneSpec("TP53", "NM_000546.6", 7157, "Li-Fraumeni syndrome"),
        GeneSpec("MLH1", "NM_000249.4", 4292, "Lynch syndrome"),
        GeneSpec("MSH2", "NM_000251.3", 4436, "Lynch syndrome"),
        GeneSpec("MSH6", "NM_000179.3", 2956, "Lynch syndrome"),
        GeneSpec("PMS2", "NM_000535.7", 5395, "Lynch syndrome"),
        GeneSpec("APC", "NM_000038.6", 324, "Familial adenomatous polyposis"),
        GeneSpec("PTEN", "NM_000314.8", 5728, "PTEN hamartoma tumor syndrome"),
        GeneSpec("STK11", "NM_000455.5", 6794, "Peutz-Jeghers syndrome"),
        GeneSpec("CDH1", "NM_004360.5", 999, "Hereditary diffuse gastric cancer"),
        GeneSpec("VHL", "NM_000551.4", 7428, "Von Hippel-Lindau disease"),
        GeneSpec("RET", "NM_020975.6", 5979, "Multiple endocrine neoplasia 2"),
        GeneSpec("MEN1", "NM_130799.3", 4221, "Multiple endocrine neoplasia 1"),
        GeneSpec("LDLR", "NM_000527.5", 3949, "Familial hypercholesterolemia"),
        GeneSpec("APOB", "NM_000384.3", 338, "Familial hypercholesterolemia"),
        GeneSpec("PCSK9", "NM_174936.4", 255738, "Familial hypercholesterolemia"),
        GeneSpec("MYH7", "NM_000257.4", 4625, "Hypertrophic cardiomyopathy"),
        GeneSpec("MYBPC3", "NM_000256.3", 4607, "Hypertrophic cardiomyopathy"),
        GeneSpec("KCNQ1", "NM_000218.3", 3784, "Long QT syndrome"),
        GeneSpec("KCNH2", "NM_000238.4", 3757, "Long QT syndrome"),
        GeneSpec("SCN5A", "NM_198056.3", 6331, "Long QT / Brugada syndrome"),
        GeneSpec("LMNA", "NM_170707.4", 4000, "Dilated cardiomyopathy / laminopathy"),
        GeneSpec("FBN1", "NM_000138.5", 2200, "Marfan syndrome"),
        GeneSpec("RYR1", "NM_000540.3", 6261, "Malignant hyperthermia"),
        GeneSpec("CFTR", "NM_000492.4", 1080, "Cystic fibrosis"),
        GeneSpec("PAH", "NM_000277.3", 5053, "Phenylketonuria"),
        GeneSpec("GAA", "NM_000152.5", 2548, "Pompe disease"),
        GeneSpec("MECP2", "NM_004992.4", 4204, "Rett syndrome"),
        GeneSpec("SCN1A", "NM_001165963.4", 6323, "Dravet syndrome / epilepsy"),
    ]
}


def is_supported_gene(symbol: str) -> bool:
    return symbol.upper() in GENE_ALLOWLIST


def get_gene(symbol: str) -> GeneSpec:
    key = symbol.upper()
    if key not in GENE_ALLOWLIST:
        raise ValueError(
            f"Gene {symbol!r} is not in the MVP allowlist "
            f"({len(GENE_ALLOWLIST)} supported genes). This is intentional scope control."
        )
    return GENE_ALLOWLIST[key]
