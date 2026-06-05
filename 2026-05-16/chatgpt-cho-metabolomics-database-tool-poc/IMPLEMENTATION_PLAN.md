# Metabolomics v1.1 Implementation Plan

## 1. Scope

Metabolomics v1.1 adds a biological interpretation layer on top of the validated v1.0 iCHO3K workflow.

Core architecture:

```text
iCHO3K = calculation engine
CHOmpact = interpretation / pathway aggregation layer
```

v1.1 must not modify the iCHO3K model, COBRApy FBA/FVA logic, or v1.0 calculation outputs. It consumes existing v1.0 output files and produces CHOmpact-style pathway summaries, biomarker rankings, demand sensitivity summaries, publication figures, and executive reports.

## 2. Existing v1.0 Inputs

The v1.1 layer consumes these existing outputs:

```text
results/<run>/spent_media/04_interval_rates.csv
results/<run>/icho3k_inputs/icho_exchange_constraints.csv
results/<run>/fba/fba_objective_results_by_scenario.csv
results/<run>/fba/fba_selected_fluxes_by_scenario.csv
results/<run>/fba/fba_all_reaction_fluxes_by_scenario.csv
results/<run>/fba/fva_selected_reactions_by_scenario.csv
results/<run>/fba/fva_active_reactions_by_scenario.csv
results/<run>/fba/fva_all_reactions_by_scenario.csv
results/<run>/fba/fva_reaction_range_summary.csv
```

Required scenario columns:

```text
passage_or_clone
producer_group
replicate
day_start
day_end
scenario
```

Comparison groups:

```text
High
Mother
Low
```

The implementation must allow arbitrary labels, but these three labels should be first-class supported groups for reports and defaults.

## 3. New v1.1 File Layout

Recommended new files:

```text
data/chompact_pathway_mapping.csv
pipeline/10_map_to_chompact.py
pipeline/11_score_chompact_pathways.py
pipeline/12_demand_sensitivity.py
pipeline/13_rank_pathway_biomarkers.py
pipeline/14_make_chompact_figures.py
```

Recommended output folder:

```text
results/<run>/chompact/
```

## 4. `data/chompact_pathway_mapping.csv`

### Purpose

Map iCHO3K reaction-level outputs to CHOmpact-style pathway categories without replacing iCHO3K calculations.

### Required Columns

```text
reaction_id
reaction_name
subsystem
chompact_pathway
chompact_subpathway
interpretation_note
```

### Recommended Additional Columns

```text
reaction_class
direction_hint
priority_for_figures
is_measured_exchange
is_product_related
is_central_carbon
```

These optional columns make scoring and figure filtering more reproducible.

### Pathway Vocabulary

Initial controlled vocabulary:

```text
Central carbon metabolism
Glycolysis
TCA cycle
Pentose phosphate pathway
Lactate metabolism
Glutamine and glutamate metabolism
Amino acid metabolism
Nucleotide metabolism
Lipid metabolism
Energy and redox metabolism
Transport and exchange
Product synthesis and IgG assembly
Other / unmapped
```

### Example Rows

```csv
reaction_id,reaction_name,subsystem,chompact_pathway,chompact_subpathway,interpretation_note
EX_glc_e,D-Glucose exchange,Exchange,Central carbon metabolism,Glucose uptake,Measured extracellular glucose uptake constraint
EX_lac_L_e,L-Lactate exchange,Exchange,Lactate metabolism,Lactate exchange,Measured lactate secretion or consumption
EX_nh4_e,Ammonium exchange,Exchange,Glutamine and glutamate metabolism,Ammonia burden,Measured nitrogen byproduct secretion
EX_gln_L_e,L-Glutamine exchange,Exchange,Glutamine and glutamate metabolism,Glutamine uptake,Measured glutamine uptake constraint
EX_asn_L_e,L-Asparagine exchange,Exchange,Amino acid metabolism,Asparagine uptake,Potential antibody-production amino acid limitation
DM_igg_g,IgG demand,Demand,Product synthesis and IgG assembly,IgG product demand,Product objective / secretion proxy
igg_formation,IgG formation,Product,Product synthesis and IgG assembly,IgG assembly,IgG assembly objective reaction
```

### Mapping Strategy

1. Start from all reaction IDs appearing in v1.0 FBA/FVA outputs.
2. Map measured exchange reactions first.
3. Map product-related reactions next.
4. Map high-flux and high-FVA-priority internal reactions.
5. Leave uncertain reactions as `Other / unmapped`.
6. Do not force weak CHOmpact assignments. Unmapped is preferable to misleading categories.

### Validation Rules

`10_map_to_chompact.py` should fail or warn if:

```text
reaction_id is missing
chompact_pathway is missing
duplicate reaction_id exists with conflicting pathway labels
>30% of high-priority reactions are unmapped
required measured exchange reactions are unmapped
```

Measured exchange reactions that must be mapped:

```text
EX_glc_e
EX_lac_L_e
EX_nh4_e
EX_gln_L_e
EX_glu_L_e
EX_ala_L_e
EX_arg_L_e
EX_asn_L_e
EX_asp_L_e
EX_cys_L_e
EX_gly_e
EX_his_L_e
EX_ile_L_e
EX_leu_L_e
EX_lys_L_e
EX_met_L_e
EX_phe_L_e
EX_pro_L_e
EX_ser_L_e
EX_thr_L_e
EX_trp_L_e
EX_tyr_L_e
EX_val_L_e
```

## 5. `pipeline/10_map_to_chompact.py`

### Purpose

Attach CHOmpact pathway labels to v1.0 FBA/FVA outputs.

### Inputs

```text
--fba-fluxes results/<run>/fba/fba_selected_fluxes_by_scenario.csv
--fba-all-fluxes results/<run>/fba/fba_all_reaction_fluxes_by_scenario.csv
--fva results/<run>/fba/fva_selected_reactions_by_scenario.csv
--fva-active results/<run>/fba/fva_active_reactions_by_scenario.csv
--mapping data/chompact_pathway_mapping.csv
--outdir results/<run>/chompact
```

All FBA/FVA inputs should be optional except at least one FBA flux table.

### Outputs

```text
mapped_fba_fluxes.csv
mapped_fva_ranges.csv
mapping_qc_summary.csv
unmapped_reactions.csv
```

### Core Logic

1. Load mapping table.
2. Load FBA and FVA outputs.
3. Left-join by `reaction_id`.
4. Fill missing `chompact_pathway` with `Other / unmapped`.
5. Write mapped outputs and QC.

### QC Metrics

```text
n_total_reactions
n_mapped_reactions
mapping_rate
n_unmapped_high_flux_reactions
n_unmapped_high_fva_reactions
n_unmapped_measured_exchange_reactions
```

### CLI Example

```bash
python pipeline/10_map_to_chompact.py \
  --fba-fluxes results/interactive_run/fba/fba_selected_fluxes_by_scenario.csv \
  --fva results/interactive_run/fba/fva_selected_reactions_by_scenario.csv \
  --mapping data/chompact_pathway_mapping.csv \
  --outdir results/interactive_run/chompact
```

## 6. `pipeline/11_score_chompact_pathways.py`

### Purpose

Aggregate reaction-level FBA/FVA outputs and measured qMet into pathway-level scores.

### Inputs

```text
--mapped-fba results/<run>/chompact/mapped_fba_fluxes.csv
--mapped-fva results/<run>/chompact/mapped_fva_ranges.csv
--rates results/<run>/spent_media/04_interval_rates.csv
--constraints results/<run>/icho3k_inputs/icho_exchange_constraints.csv
--reference-group Mother
--compare-groups Low High
--outdir results/<run>/chompact
```

### Outputs

```text
pathway_activity_scores.csv
pathway_high_low_separation.csv
pathway_measured_qmet_scores.csv
pathway_fva_robustness_scores.csv
pathway_score_qc.csv
```

### Score Definitions

#### 6.1 Measured qMet Pathway Score

Use measured exchange rates only.

```text
measured_qmet_activity =
median(abs(qmet_per_1e6_cells_day))
```

By group/pathway:

```text
group_qmet_score =
median(abs(qmet_per_1e6_cells_day)) grouped by producer_group + chompact_pathway
```

This is the strongest biological evidence because it comes from measured spent-media data.

#### 6.2 FBA Flux Activity Score

```text
pathway_flux_activity =
median(abs(flux)) across mapped reactions
```

Use median rather than mean to reduce dominance by one high-flux reaction.

Recommended grouping:

```text
producer_group
passage_or_clone
day_start
day_end
objective
chompact_pathway
```

#### 6.3 High-Low Separation Score

For each pathway:

```text
pathway_separation =
median(pathway_score_high) - median(pathway_score_low)
```

For Mother comparison:

```text
high_vs_mother =
median(pathway_score_high) - median(pathway_score_mother)

low_vs_mother =
median(pathway_score_low) - median(pathway_score_mother)
```

Standardized version:

```text
standardized_separation =
(mean_high - mean_low) / pooled_sd
```

If `pooled_sd = 0`, return null and flag as non-estimable.

#### 6.4 FVA Tightness Score

For each reaction:

```text
fva_range = maximum - minimum
fva_tightness = 1 / (1 + abs(fva_range))
```

For each pathway:

```text
pathway_fva_tightness =
median(fva_tightness)
```

Interpretation:

```text
High tightness = less alternate flux flexibility
Low tightness = flux can vary widely under the same objective
```

#### 6.5 Robust Pathway Score

```text
robust_pathway_score =
abs(standardized_high_low_separation) * pathway_fva_tightness
```

Optional composite:

```text
pathway_priority_score =
0.35 * measured_qMet_separation_z
+ 0.30 * FBA_flux_separation_z
+ 0.25 * FVA_tightness_z
+ 0.10 * time_consistency_z
```

The composite score is a candidate ranking score, not causal proof.

### Minimum QC

Flag pathways with:

```text
<3 mapped reactions
no measured qMet support
missing High/Mother/Low group
FVA unavailable
all-zero FBA flux
```

## 7. `pipeline/12_demand_sensitivity.py`

### Purpose

Review demand scaling assumptions and identify pathway rankings that remain stable across demand strategies.

This script should not modify the iCHO3K engine. It should either:

1. consume separate v1.0 FBA/FVA outputs generated with different demand settings, or
2. generate a command manifest for rerunning v1.0 FBA under defined demand settings.

### Demand Strategies

Compare:

```text
fixed_demand_scale
titer_normalized_demand
qP_normalized_demand
IVCD_normalized_demand
measured_IgG_secretion_rate_constraint
```

Biological defensibility ranking:

```text
1. measured_IgG_secretion_rate_constraint
2. qP_normalized_demand
3. IVCD_normalized_demand
4. titer_normalized_demand
5. fixed_demand_scale
```

### Inputs

```text
--runsheet data/demand_sensitivity_runsheet.csv
--score-files results/<run>/chompact/pathway_priority_scores.csv
--outdir results/<run>/chompact/demand_sensitivity
```

### Runsheet Columns

```text
demand_strategy
demand_scale
fba_result_dir
notes
```

### Outputs

```text
demand_scale_sensitivity.csv
robust_rank_across_demand_scales.csv
demand_strategy_qc.csv
demand_sensitivity_report.html
```

### Robustness Metrics

For each pathway:

```text
rank_median
rank_iqr
score_median
score_iqr
sign_consistency
top_n_frequency
```

Definitions:

```text
sign_consistency =
fraction of demand settings where High-Low separation has the same sign

top_n_frequency =
fraction of demand settings where pathway appears in top N
```

Recommended robust biomarker criteria:

```text
top_10_frequency >= 0.60
sign_consistency >= 0.80
rank_iqr <= 10
```

### Reporting Language

Use:

```text
Demand scaling sensitivity identified pathway rankings stable across biologically plausible product-demand assumptions.
```

Avoid:

```text
Demand scale proves the exact IgG flux.
```

## 8. `pipeline/13_rank_pathway_biomarkers.py`

### Purpose

Rank CHOmpact pathways by biological interpretability, High-Low separation, FVA robustness, and time/clone consistency.

### Inputs

```text
--pathway-scores results/<run>/chompact/pathway_activity_scores.csv
--separation results/<run>/chompact/pathway_high_low_separation.csv
--fva-scores results/<run>/chompact/pathway_fva_robustness_scores.csv
--demand-sensitivity results/<run>/chompact/demand_sensitivity/robust_rank_across_demand_scales.csv
--outdir results/<run>/chompact
```

### Outputs

```text
ranked_pathway_biomarkers.csv
top_pathway_biomarkers_for_report.csv
biomarker_evidence_matrix.csv
```

### Ranking Columns

```text
chompact_pathway
chompact_subpathway
measured_qmet_separation
fba_flux_separation
fva_tightness
time_consistency
clone_consistency
demand_rank_stability
pathway_priority_score
evidence_level
interpretation_note
```

### Evidence Level Rules

```text
Strong:
measured qMet support + FBA separation + FVA tightness + demand robustness

Moderate:
two or three evidence types agree

Hypothesis:
model-predicted only or weak measured support
```

### Classification Potential

For each pathway, compute simple clone classification metrics:

```text
effect_size_high_vs_low
auc_high_vs_low
threshold_direction
misclassification_count
```

This should be treated as exploratory because sample count may be low.

## 9. `pipeline/14_make_chompact_figures.py`

### Purpose

Generate publication-quality CHOmpact interpretation figures and an automated report.

### Inputs

```text
--pathway-scores results/<run>/chompact/pathway_activity_scores.csv
--biomarkers results/<run>/chompact/ranked_pathway_biomarkers.csv
--objectives results/<run>/fba/fba_objective_results_by_scenario.csv
--rates results/<run>/spent_media/04_interval_rates.csv
--outdir results/<run>/chompact/figures
```

### Outputs

```text
figure_01_culture_phenotype.svg
figure_02_measured_qmet_heatmap.svg
figure_03_icho3k_objective_by_clone.svg
figure_04_chompact_pathway_activity.svg
figure_05_fva_robustness_vs_separation.svg
figure_06_pathway_biomarker_summary.svg
chompact_report.html
executive_summary.md
```

### Figure Requirements

General style:

```text
white background
no 3D effects
consistent pathway colors
readable labels
explicit measured vs predicted labels
SVG and HTML outputs
```

Figure 1:

```text
Culture phenotype:
VCD, viability, titer, qP if available
```

Figure 2:

```text
Measured qMet heatmap:
Glucose, lactate, ammonia, glutamine, glutamate, 20 amino acids
```

Figure 3:

```text
iCHO3K objective prediction:
biomass_cho_prod, DM_igg_g, igg_formation
```

Figure 4:

```text
CHOmpact pathway activity:
High vs Mother vs Low pathway score
```

Figure 5:

```text
FVA robustness plot:
x-axis = High-Low separation
y-axis = FVA tightness
point size = measured qMet support
color = CHOmpact pathway
```

Figure 6:

```text
Biomarker summary:
Top pathways with evidence-level labels
```

### Report Structure

```text
1. Executive Summary
2. Data and Run QC
3. Measured Metabolomics Results
4. iCHO3K FBA/FVA Predictions
5. CHOmpact Pathway Interpretation
6. Demand Sensitivity
7. Biomarker Candidates
8. Limitations and Recommended Follow-up
```

## 10. Main Scientific Interpretation

Recommended wording:

```text
Measured spent-media rates were converted to feed-corrected qMet values and used as exchange constraints for iCHO3K FBA/FVA.

iCHO3K was retained as the calculation model, while CHOmpact-style pathway categories were used to aggregate reaction-level FBA/FVA outputs into biologically interpretable pathway scores.

High-productivity clone signatures were prioritized when measured qMet separation, iCHO3K-predicted flux separation, FVA robustness, and demand-sensitivity stability were concordant.
```

Avoid:

```text
Unmeasured intracellular fluxes were experimentally confirmed.
CHOmpact replaced iCHO3K.
FBA proves causal bottlenecks.
```

Use:

```text
model-predicted
constraint-based
candidate pathway
robust under FVA
stable across demand assumptions
```

## 11. Implementation Milestones

### Milestone 1: Mapping Layer

Deliver:

```text
data/chompact_pathway_mapping.csv
pipeline/10_map_to_chompact.py
mapped_fba_fluxes.csv
mapped_fva_ranges.csv
mapping_qc_summary.csv
```

Acceptance:

```text
All measured exchange reactions mapped
Product objective reactions mapped
Mapping QC generated
Unmapped reactions listed
```

### Milestone 2: Pathway Scoring

Deliver:

```text
pipeline/11_score_chompact_pathways.py
pathway_activity_scores.csv
pathway_high_low_separation.csv
pathway_fva_robustness_scores.csv
```

Acceptance:

```text
High/Mother/Low comparisons supported
Measured and predicted scores separated
FVA tightness included
```

### Milestone 3: Demand Sensitivity

Deliver:

```text
pipeline/12_demand_sensitivity.py
demand_scale_sensitivity.csv
robust_rank_across_demand_scales.csv
```

Acceptance:

```text
Multiple demand strategies compared
Rank stability quantified
Unstable pathway rankings flagged
```

### Milestone 4: Biomarker Ranking

Deliver:

```text
pipeline/13_rank_pathway_biomarkers.py
ranked_pathway_biomarkers.csv
biomarker_evidence_matrix.csv
```

Acceptance:

```text
Evidence levels assigned
Classification potential estimated
Pathway candidates ranked
```

### Milestone 5: Figures and Report

Deliver:

```text
pipeline/14_make_chompact_figures.py
chompact_report.html
executive_summary.md
publication SVGs
```

Acceptance:

```text
Measured vs predicted results visually separated
Pathway-level interpretation clear
Executive summary generated automatically
```

## 12. Risks and Controls

| Risk | Control |
|---|---|
| CHOmpact mapping is incomplete | Keep `Other / unmapped`; report mapping rate |
| Reaction-level FBA alternate optima overinterpreted | Use FVA tightness and demand sensitivity |
| Demand scaling biases pathway rankings | Compare multiple demand strategies |
| Measured and predicted values mixed in report | Label all figures as measured or model-predicted |
| Too many reactions in figures | Aggregate by pathway and show top-ranked candidates |
| Low sample count for classification | Mark classification as exploratory |

## 13. Final v1.1 Deliverable Definition

v1.1 is complete when the pipeline can:

```text
1. consume v1.0 iCHO3K FBA/FVA outputs
2. map reactions to CHOmpact pathway categories
3. calculate pathway-level qMet, FBA, and FVA scores
4. compare High / Mother / Low clone groups
5. evaluate demand-scale robustness
6. rank pathway biomarkers
7. generate publication-ready figures
8. generate an executive summary report
```

The primary v1.1 result should be:

```text
FVA-robust CHOmpact pathway biomarkers that separate High, Mother, and Low productivity phenotypes, supported by measured qMet and iCHO3K-constrained FBA/FVA predictions.
```
