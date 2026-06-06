# v1.2 Phase 1 PPP Mapping Audit

## Scope

This phase expands the CHOmpact interpretation dictionary only. It does not
change iCHO3K, qMet constraints, pFBA, or FVA calculations.

The iCHO3K production model contains 39 reactions assigned to the
`PENTOSE PHOSPHATE PATHWAY` subsystem. The v1.1 dictionary mapped 6 canonical
cytosolic PPP reactions.

## Added High-Confidence Reactions

| Reaction | CHOmpact subpathway | Reason for inclusion |
| --- | --- | --- |
| `G6PDH2rer` | Oxidative PPP | ER glucose-6-phosphate dehydrogenase using NADP/NADPH |
| `PGLer` | Oxidative PPP | ER 6-phosphogluconolactonase |
| `GNDer` | Oxidative PPP | ER phosphogluconate dehydrogenase using NADP/NADPH |
| `RPE` | Non-oxidative PPP | Canonical ribulose-5-phosphate epimerase |
| `RPI` | Non-oxidative PPP | Canonical ribose-5-phosphate isomerase |
| `r0249` | Non-oxidative PPP | ER ribose-5-phosphate/ribulose-5-phosphate isomerase |
| `PRPPS` | PPP-linked nucleotide precursor bridge | R5P to PRPP handoff into nucleotide biosynthesis; not core oxidative or non-oxidative PPP |

Expected dictionary coverage changes from 6/39 (15.38%) to 13/39 (33.33%).
This remains `partially_covered`.

`PRPPS` is included to preserve the biologically important
ribose-5-phosphate -> PRPP -> nucleotide biosynthesis handoff. It is a
PPP-linked bridge and must not be interpreted as a core oxidative or
non-oxidative PPP reaction.

## Confidence Guardrails

Mapping coverage is an interpretation/QC metric, not affirmative evidence that
a pathway is active or robust. Increasing mapping coverage must not
automatically increase `confidence_score`.

- Coverage may cap or limit confidence when ontology support is incomplete.
- FVA robustness, feasible-range separation, and reproducibility should
  dominate confidence.
- A newly mapped reaction that is available only in full FVA remains a
  model-emergent hypothesis until robust FVA support is demonstrated.
- Loop-prone, very broad, or under-constrained FVA ranges must not inflate
  confidence merely because the reaction now has a pathway label.

No confidence or ranking calculation was changed in this phase.

## Intentionally Excluded

The following groups remain unmapped in Phase 1:

- NAD-dependent alternatives such as `G6PDH1rer` and `GAUGE-R10221`.
- Gap-fill reactions without sufficiently specific biological annotation.
- Pentose, deoxyribose, xylitol, arabinitol, and glucuronate salvage reactions.
- UDP-glucuronate hydrolysis reactions, which are better reviewed with the
  nucleotide-sugar/glycosylation ontology.
- Reactions assigned to the PPP subsystem whose reaction chemistry is unrelated
  to canonical PPP, such as peptide hydrolysis.

Their exclusion does not imply zero flux or biological irrelevance. It means
that Phase 1 does not classify them as high-confidence PPP interpretation
reactions.

## Validation

Validation dataset: `practice_20aa`

```text
python run_pipeline.py --dataset practice_20aa --input_format tsv \
  --steps 10,11,12,13,14 --fva_source auto
```

Result:

```text
Steps 10-14: PASS
model_reaction_count: 39
mapped_reaction_count: 13
mapping_coverage_percent: 33.3333
analysis_status: partially_covered
Fig14: generated
Fig14 companion CSV: generated
```

The validation worktree contained focused FVA but no full all/internal FVA
table, so `auto` selected `focused_fva / central_mab_panel`. PPP dictionary
mapping increased from 6 to 13 reactions, but focused-FVA scoring and ranking
still covered only the original 6 reactions. The seven newly mapped reactions
are outside that focused panel.

Full-FVA workstation validation is required to confirm scoring behavior for
the seven newly mapped reactions, including whether their feasible ranges are
robust, separable, or loop-prone. The 33.33% mapping coverage result alone must
not be interpreted as higher pathway confidence.
