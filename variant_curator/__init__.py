"""variant_curator: a scoped ACMG variant-curation assistant.

Scope (MVP): missense and nonsense SNVs in a fixed allowlist of well-characterized
disease genes. Assembles evidence from public sources (Ensembl VEP, gnomAD, ClinVar),
maps it to ACMG criteria with a linked source for every code, and abstains where the
evidence isn't there.

No PHI. Built entirely on public data. Designed so a clinical site could later point
the same engine at real cases behind its own firewall without code changes.
"""

__version__ = "0.1.0"
