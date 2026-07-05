"""CLI: run a causality review for one patient against the reported variants.

Step 2 deliverable: take the variants the lab reported plus the patient's HPO
phenotype profile, and rank the variants by how well each explains THIS patient,
separating explained from unexplained features and flagging when leftover
features imply a second cause or a genome re-analysis.

Usage:
    # hardcoded demo patient (SCN1A best fit; a mismatched MYH7 candidate)
    python -m causality_review.cli

    # your own case: reported variants (gene[:hgvs]) + patient features
    python -m causality_review.cli \
        --variant SCN1A:c.3637C>T --variant MYH7:c.1063G>A \
        --hpo "seizures" --hpo "global developmental delay" --hpo "ataxia"
"""

from __future__ import annotations

import argparse

from .clients.hpo import resolve_term
from .engine import review_causality
from .http import get_client
from .models import CausalityReport, PatientProfile, Phenotype, ReportedVariant

# Hardcoded demo: an infant with an epileptic-encephalopathy picture. The lab
# reported two candidates -- a neuronal sodium-channel gene (SCN1A) and a
# cardiomyopathy gene (MYH7). The review should rank SCN1A far above MYH7 and
# flag the one feature neither gene explains.
DEMO_VARIANTS = [
    ReportedVariant("SCN1A", "c.3637C>T", "Likely pathogenic"),
    ReportedVariant("MYH7", "c.1063G>A", "VUS"),
]
DEMO_HPO = [
    Phenotype("HP:0001250", "Seizure"),
    Phenotype("HP:0001263", "Global developmental delay"),
    Phenotype("HP:0002373", "Febrile seizure"),
    Phenotype("HP:0001251", "Ataxia"),
    Phenotype("HP:0000175", "Cleft palate"),  # explained by neither gene
]


def _print_report(report: CausalityReport) -> None:
    p = report.patient
    print("=" * 72)
    print("  Causality review — does a reported variant explain THIS patient?")
    print("=" * 72)

    print("\n[patient] presented features (HPO)")
    for ph in p.phenotypes:
        print(f"  - {ph.hpo_id:12s} {ph.label}")

    print("\n[causality ranking] reported variants, best fit first")
    for rank, f in enumerate(report.fits, 1):
        print(
            f"\n  {rank}. [{f.tier.stars}] {f.variant.label:22s} "
            f"{f.tier.display:9s} (explains {len(f.explained)}/{len(p.phenotypes)})"
        )
        if f.variant.lab_classification:
            print(f"       lab call  : {f.variant.lab_classification}")
        if f.explained:
            print(f"       explains  : {', '.join(m.display for m in f.explained)}")
        if f.unexplained:
            print(f"       leaves    : {', '.join(ph.label for ph in f.unexplained)}")
        if f.diseases:
            print(f"       disease   : {'; '.join(f.diseases[:3])}")
        if f.source:
            print(f"       source    : {f.source.url}")

    if report.flags:
        print("\n[flags] what the clinician should weigh")
        for fl in report.flags:
            print(f"  - {fl}")
    print()


def _build_patient(client, hpo_args: list[str]) -> PatientProfile:
    phenotypes: list[Phenotype] = []
    for raw in hpo_args:
        term = resolve_term(client, raw)
        if term is None:
            print(f"  (could not resolve HPO term for {raw!r}; skipping)")
            continue
        phenotypes.append(term)
    return PatientProfile(phenotypes=phenotypes)


def main() -> None:
    parser = argparse.ArgumentParser(description="Rank reported variants by fit to a patient.")
    parser.add_argument(
        "--variant", action="append", default=[],
        help="Reported variant as GENE or GENE:hgvs, e.g. SCN1A:c.3637C>T (repeatable)",
    )
    parser.add_argument(
        "--hpo", action="append", default=[],
        help="Patient feature as free text or HP:xxxxxxx (repeatable)",
    )
    args = parser.parse_args()

    with get_client() as client:
        if args.variant and args.hpo:
            variants = []
            for spec in args.variant:
                gene, _, hgvs = spec.partition(":")
                variants.append(ReportedVariant(gene=gene.strip(), hgvs_c=hgvs.strip()))
            patient = _build_patient(client, args.hpo)
        else:
            print("(no --variant/--hpo given; running hardcoded demo case)\n")
            variants = DEMO_VARIANTS
            patient = PatientProfile(phenotypes=DEMO_HPO)

        if not patient.phenotypes:
            print("No resolvable patient features; nothing to score.")
            raise SystemExit(2)

        report = review_causality(client, patient, variants)
        _print_report(report)


if __name__ == "__main__":
    main()
