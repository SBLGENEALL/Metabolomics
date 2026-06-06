# Metabolomics v1.1 PR1 — Industrial Decision Layer Plan

Branch:

```text
v1.1-pr1-release-polish
```

Parent branch:

```text
v1.1
```

Purpose:

```text
Turn v1.1 from a CHOmpact pathway-reporting workflow into a defensible industrial CHO cell-line development decision-support layer.
```

This PR must remain a v1.1 polish/fix release. It should not add new v1.2-scale capabilities such as flux sampling, ecFBA, transcriptomics-GPR integration, or secretion-pathway modeling.

---

## PR1 implementation status

Implemented in steps 10-14:

```text
1. Reaction-level evidence provenance labels
2. Measured High-Low qMet effect scoring
3. Reaction-level FBA/FVA separation before pathway aggregation
4. Separate measured, model-emergent, and demand-conditioned ranking tracks
5. Separate priority and confidence scores
6. NA-preserving missing-evidence handling
7. Predictive-ranking exclusion for product-demand-driven IgG reactions
8. Evidence-aware Figures 12-14
9. Three-layer industrial executive summary
```

New decision-support tables:

```text
results/<dataset>/tables/chompact/measured_screening_markers.csv
results/<dataset>/tables/chompact/model_emergent_pathway_hypotheses.csv
results/<dataset>/tables/chompact/demand_conditioned_explanations.csv
results/<dataset>/tables/chompact/reaction_level_candidates.csv
results/<dataset>/tables/chompact/pathway_level_candidates.csv
```

Supporting evidence tables:

```text
chompact_measured_high_low_effects.csv
chompact_reaction_high_low_separation.csv
chompact_industrial_ranking_qc.csv
```

Validation command:

```text
python run_pipeline.py --dataset practice_20aa --steps 10,11,12,13,14
```

PR1 retains the existing v1.1 filenames for compatibility, but their content
uses candidate-signature and hypothesis-ranking terminology. Steps 00-09 are
unchanged.

---

## 1. Consensus from Claude and Codex review

Both biological and analytics reviews agree on the same core issue:

```text
The current v1.1 ranking mixes measured evidence, model-predicted evidence, exchange-constraint-driven outputs, and product-demand-driven outputs into one score.
```

For industrial CLD, this is dangerous because it can make an input constraint look like an independent discovery.

Therefore, PR1 must split evidence by provenance and decision value before ranking or reporting.

---

## 2. Evidence framework

Use two independent axes plus one reconciliation flag.

### Axis 1 — Observability

```text
[M] measured
  Experimental or directly derived from experimental data.

[P] predicted-only
  Internal model-predicted value, not directly measured.
```

### Axis 2 — Source of High-Low difference

```text
[C] exchange-constraint-driven
  Difference is directly driven by measured exchange bounds.

[D] product-demand-driven
  Difference is directly driven by imposed IgG/qP demand.

[E] model-emergent
  Difference is not directly imposed by a boundary or product-demand constraint.
```

### Reconciliation flag

```text
reconciled = yes / no / untested
```

This indicates whether the model-predicted signal agrees with independent measured phenotype or external validation.

### Examples

| Reaction/pathway | Observability | Difference source | Interpretation |
|---|---:|---:|---|
| Glucose uptake | M | C | Screening/feed marker, not emergent mechanism |
| Lactate exchange | M | C | Screening/feed marker; lactate-shift dynamics may be emergent if not directly imposed |
| Glutamine uptake | M | C | Feed/nutrient-risk marker |
| ATP synthase | P | E | Model-emergent oxidative hypothesis |
| ETC Complex I/III | P | E | Model-emergent oxidative hypothesis |
| Citrate synthase | P | E | Model-emergent TCA hypothesis |
| Pyruvate dehydrogenase | P | E | Model-emergent oxidative-entry hypothesis |
| DM_igg_g | P | D | Product-demand reconstruction, not independent biomarker |
| IgG assembly / HC / LC | P | D | By construction; exclude from predictive biomarker ranking |

---

## 3. Industrial interpretation hierarchy

Present results in this order.

### Level 1 — Measured phenotype

Highest evidence strength.

Examples:

- qIgG / titer
- VCD / viability / growth
- glucose uptake
- lactate secretion
- Lac/Glc ratio
- ammonia burden
- amino-acid uptake/depletion risk

Use for:

- clone screening
- feed/media triage
- sanity check of CHO physiology

### Level 2 — Model-emergent metabolic signatures

Mechanistic hypotheses supported by iCHO3K under measured constraints.

Examples:

- OxPhos / ETC
- ATP synthase
- citrate synthase
- TCA activity
- PDH/pyruvate oxidative entry
- glutamine catabolism/anaplerosis

Use for:

- mechanistic interpretation
- engineering hypotheses
- validation experiment prioritization

### Level 3 — Demand-conditioned explanations

Outputs that explain imposed measured IgG demand.

Examples:

- IgG demand
- IgG assembly
- HC synthesis
- LC synthesis

Use for:

- production-burden reconstruction
- explanatory context

Do not use for:

- independent productivity biomarker discovery
- clone selection claims
- engineering target claims without independent validation

---

## 4. Required v1.1 PR1 fixes

### 4.1 Evidence labels

Add evidence provenance columns to key v1.1 tables.

Required columns:

```text
evidence_observability
high_low_difference_source
reconciliation_status
decision_role
interpretation_guardrail
```

Recommended values:

```text
evidence_observability: measured | predicted_only
high_low_difference_source: exchange_constraint_driven | product_demand_driven | model_emergent | mixed | unknown
reconciliation_status: reconciled | not_reconciled | untested | not_applicable
decision_role: clone_selection_marker | feed_media_marker | engineering_hypothesis | explanatory_context | do_not_rank_as_predictive
```

### 4.2 Ranking separation

Separate outputs into at least three tracks:

```text
measured_screening_markers.csv
model_emergent_pathway_hypotheses.csv
demand_conditioned_explanations.csv
```

Do not allow product-demand-driven IgG reactions to appear as independent predictive biomarkers.

### 4.3 Rename composite score

Current language such as:

```text
biomarker score
validated biomarker
```

must be replaced with:

```text
hypothesis ranking score
candidate pathway signature
constraint-informed pathway hypothesis
```

### 4.4 Priority and confidence separation

Industrial decision support must separate:

```text
priority_score
confidence_score
```

Priority answers:

```text
How important/actionable is this signal?
```

Confidence answers:

```text
How reliable is this signal?
```

### 4.5 Missing evidence handling

Do not treat missing evidence as zero.

Required distinction:

```text
not measured
measured but no difference
measured with strong difference
```

### 4.6 Reaction-level vs pathway-level outputs

Create separate tables for:

```text
reaction_level_candidates.csv
pathway_level_candidates.csv
```

Reaction-level ranking is for engineering target hypotheses.

Pathway-level ranking is for biological interpretation and CLD decision context.

### 4.7 Product-demand warning

IgG demand, IgG assembly, HC synthesis, and LC synthesis should be labeled:

```text
product-demand-driven / by construction / explanatory only
```

They may be reported, but should not be ranked as model-emergent discoveries.

---

## 5. Industrial decision-support outputs

PR1 should start moving the report toward this structure.

### 5.1 Executive decision section

Must answer:

1. Which clone or clone group should advance?
2. What is the main measured evidence?
3. What is the main model-emergent hypothesis?
4. What are the key nutrient/feed risks?
5. What confirmation experiment is recommended?
6. What is the confidence level?

### 5.2 Clone advancement status

Use these categories:

```text
ADVANCE
ADVANCE_WITH_CONDITIONS
HOLD_FOR_CONFIRMATION
DO_NOT_ADVANCE
INSUFFICIENT_DATA
```

### 5.3 Clone advancement score components

For future multi-clone datasets, use separate component scores:

```text
product_performance_score
culture_robustness_score
metabolic_efficiency_score
model_supported_capacity_score
operational_confidence_score
```

A consensus score may be shown, but component scores must remain visible.

### 5.4 Nutrient limitation risk

Use cautious language:

```text
confirmed_depletion
high_limitation_risk
potential_limitation
no_current_evidence
insufficient_data
```

Do not claim true limitation without concentration/depletion and perturbation evidence.

### 5.5 Engineering hypothesis table

Use the term:

```text
engineering target hypothesis
```

not:

```text
validated engineering target
```

Required fields:

```text
target_reaction
associated_gene_or_gpr
predicted_intervention_direction
supporting_evidence
confidence_level
potential_risk
recommended_validation_experiment
```

---

## 6. Figure direction for PR1

Current figures 12-14 can remain in v1.1, but their interpretation must change.

### Fig12

Current:

```text
CHOmpact pathway activity
```

PR1 direction:

```text
Signed pathway effect / evidence-aware pathway activity
```

Must distinguish measured, model-emergent, and demand-conditioned evidence.

### Fig13

Current:

```text
FVA robustness bar plot
```

PR1 direction:

```text
FVA robustness by evidence type
```

Product-demand-driven IgG reactions must be marked separately or excluded from predictive ranking.

### Fig14

Current:

```text
single composite biomarker ranking
```

PR1 direction:

```text
priority-confidence matrix or hypothesis ranking with evidence tags
```

Do not present the bar chart as validated biomarker discovery.

---

## 7. Recommended validation experiments

Top three biological validation experiments for CLD relevance:

1. Central carbon/TCA 13C-MFA in High vs Low clones.
2. Respiratory/OXPHOS functional assay such as Seahorse OCR/ECAR or OUR.
3. Causal feed or metabolic perturbation test, such as glucose-limited feed, lactate-management test, or PDH/lactate-axis perturbation.

Rationale:

```text
These validate whether the model-emergent oxidative signature is real, functional, and actionable.
```

---

## 8. v1.2 boundaries

Do not include these in v1.1 PR1 unless explicitly approved later.

These justify v1.2 or later:

- flux sampling
- ecFBA
- transcriptomics-GPR integration
- clone-specific context models
- dynamic culture-phase modeling
- glycosylation model integration
- ER folding/secretion model integration
- secretion-pathway model such as iCHO2048s-style integration
- 13C-MFA validation module
- hundreds-of-clones distributed computation
- DuckDB/Parquet backend
- machine-learning productivity classifier
- prospective clone-selection validation

---

## 9. Industrial wording standard

Use:

```text
candidate pathway signature
constraint-informed pathway hypothesis
model-emergent hypothesis
screening marker
feed/media marker
engineering target hypothesis
```

Avoid:

```text
validated biomarker
proven intracellular flux
model-validated productivity driver
independent prediction from measured-demand IgG flux
```

---

## 10. PR1 success criteria

PR1 is successful if a CLD scientist can read the output and answer:

```text
Which clone should I advance?
Which evidence is measured?
Which evidence is model-predicted?
Which finding is constraint-driven?
Which pathway is an actionable hypothesis?
What should I test next?
```

PR1 is not successful if it only produces a more attractive pathway figure without improving decision traceability.
