"""causality_review: does this reported variant actually explain THIS patient?

Step 2 of the diagnostic pipeline (the clinician's causality review), as opposed
to step 1 (the lab's variant classification, handled by `variant_curator`).

The lab hands back 2-5 reported variants scored on thin clinical detail. The
clinician holds what the lab did not: the patient's full phenotype. This tool
runs the match in reverse -- for each reported variant, how well does its
gene's known disease phenotype explain THIS patient's features? -- and ranks
them, honestly separating explained from unexplained features and flagging when
leftover features imply a second cause or a genome re-analysis.

No PHI. Phenotypes are HPO term ids; gene-phenotype knowledge is public
(HPO / Jax). Designed so a site could point the same engine at real EMR-derived
HPO profiles behind their own firewall without code changes.
"""

__version__ = "0.1.0"
