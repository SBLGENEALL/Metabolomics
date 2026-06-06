# CHO Metabolomics v1.1 Implementation Plan

## Architecture

- `iCHO3K` remains the calculation engine.
- `CHOmpact` is used only as a pathway category, interpretation, and figure layer.
- v1.1 reads existing v1.0 tables from `results/<dataset>/tables`.
- v1.1 does not recompute flux using CHOmpact.
- v1.1 does not sum iCHO3K reactions to reconstruct CHOmpact lumped fluxes.
- Measured qMet/titer/VCD, model-predicted pFBA flux, and FVA feasible intervals are kept separate.
- `demand_scale` is treated as a sensitivity setting, not a fixed biological truth.

## File-by-File Plan

### `data/chompact_pathway_mapping.csv`

Purpose: reaction-level interpretation schema for iCHO3K results.

Required columns:

- `reaction_id`
- `reaction_name`
- `subsystem`
- `chompact_pathway`
- `chompact_subpathway`
- `interpretation_note`

Additional columns:

- `reaction_class`
- `is_measured_exchange`
- `is_product_related`
- `is_constraint_reaction`
- `priority_for_figures`

### `scripts/steps/10_map_to_chompact.py`

Purpose: attach CHOmpact pathway labels to v1.0 iCHO3K result tables.

Inputs:

- `results/<dataset>/tables/central_mab_flux_state.csv`
- `results/<dataset>/tables/central_mab_fva_report.csv`
- optional `results/<dataset>/tables/full_fva/full_fva_all_combined_report.csv`
- optional `results/<dataset>/tables/full_fva/full_fva_all_high_low_overlap.csv`
- `results/<dataset>/tables/exchange_rates.csv`
- `data/chompact_pathway_mapping.csv`

Outputs:

- `results/<dataset>/tables/chompact/chompact_mapped_flux.csv`
- `results/<dataset>/tables/chompact/chompact_mapped_fva.csv`
- `results/<dataset>/tables/chompact/chompact_mapped_fva_overlap.csv`
- `results/<dataset>/tables/chompact/chompact_mapped_measured_rates.csv`
- `results/<dataset>/tables/chompact/chompact_mapping_qc.csv`
- `results/<dataset>/tables/chompact/chompact_unmapped_reactions.csv`

Pseudocode:

```text
load mapping CSV
load v1.0 flux/FVA/rate tables when present
normalize rxn_id -> reaction_id
left-join mapping by reaction_id
fill unmapped reactions as "Other / unmapped"
add source_type and mapping_status
write mapped tables and QC
```

### `scripts/steps/11_score_chompact_pathways.py`

Purpose: generate pathway-level scores without reconstructing pathway flux.

Inputs:

- `chompact_mapped_flux.csv`
- `chompact_mapped_fva.csv`
- `chompact_mapped_measured_rates.csv`

Outputs:

- `chompact_pathway_fba_activity_scores.csv`
- `chompact_pathway_fva_robustness_scores.csv`
- `chompact_pathway_measured_qmet_scores.csv`
- `chompact_pathway_high_low_separation.csv`
- `chompact_pathway_score_qc.csv`

Key schema:

- `producer_group`: `High`, `Mother`, `Moderate`, `Mid`, `Low`, or inferred clone group
- `chompact_pathway`
- `chompact_subpathway`
- `source_type`
- FBA metrics: `median_flux`, `median_abs_flux`, `mean_abs_flux`
- FVA metrics: `median_fva_min`, `median_fva_max`, `median_fva_range`, `median_fva_tightness`
- separation metrics: `comparison`, `metric_family`, `delta_a_minus_b`, `abs_delta`, `fva_overlap_ratio`, `fva_non_overlap_score`

Pseudocode:

```text
load mapped flux, FVA, and measured-rate tables
infer producer_group from condition/clone labels
score FBA with median absolute flux by pathway and group
score FVA with median min/max/range/tightness by pathway and group
score measured qMet separately from exchange-rate table
compute pairwise High/Mother/Low separation
compute FVA interval overlap and non-overlap
write pathway score tables
```

### `scripts/steps/12_demand_sensitivity.py`

Purpose: summarize pathway robustness across demand_scale settings.

Inputs:

- default current run: `chompact_pathway_high_low_separation.csv`
- optional runsheet: `data/demand_sensitivity_runsheet.csv`

Optional runsheet schema:

- `run_label`
- `dataset`
- `demand_scale`
- `chompact_dir`

Outputs:

- `chompact_demand_scale_sensitivity.csv`
- `chompact_robust_rank_across_demand_scales.csv`
- `chompact_demand_strategy_qc.csv`

Pseudocode:

```text
if runsheet exists:
    load each run's CHOmpact separation table
else:
    use current run as single-setting baseline
rank pathways within each demand setting
calculate median rank, rank IQR, sign consistency, top10 frequency
label pathways as stable/sensitive/single_setting_not_assessed
write sensitivity tables
```

### `scripts/steps/13_rank_pathway_biomarkers.py`

Purpose: prioritize pathways for clone-productivity interpretation.

Inputs:

- `chompact_pathway_high_low_separation.csv`
- `chompact_robust_rank_across_demand_scales.csv`
- `chompact_pathway_measured_qmet_scores.csv`
- `chompact_pathway_fba_activity_scores.csv`

Outputs:

- `chompact_ranked_pathway_biomarkers.csv`
- `chompact_top_pathway_biomarkers_for_report.csv`
- `chompact_biomarker_evidence_matrix.csv`
- `chompact_clone_classification_potential.csv`

Pseudocode:

```text
aggregate High-vs-Low FBA separation by pathway
aggregate FVA non-overlap by pathway
aggregate measured qMet evidence by pathway
add demand-scale stability if multiple settings are available
calculate composite score from measured qMet, FBA separation, FVA non-overlap, demand stability
label evidence type and robustness
write ranked biomarker tables
```

### `scripts/steps/14_make_chompact_figures.py`

Purpose: create publication-style PNG figures and an executive summary.

Inputs:

- `chompact_pathway_fba_activity_scores.csv`
- `chompact_pathway_high_low_separation.csv`
- `chompact_ranked_pathway_biomarkers.csv`

Outputs:

- `results/<dataset>/figures/Fig12_measured_screening_markers.png`
- `results/<dataset>/figures/Fig13_model_emergent_fva_robustness.png`
- `results/<dataset>/figures/Fig14_candidate_pathway_priority_confidence.png`
- `results/<dataset>/tables/chompact/CHOmpact_v1_1_executive_summary.md`

Pseudocode:

```text
load pathway activity, FVA separation, and biomarker ranking tables
plot pathway activity grouped by High/Mother/Low
plot FVA non-overlap for robust High-vs-Low pathway separation
plot composite pathway biomarker ranking
write executive summary with interpretation boundaries
```

## CLI Commands

Run v1.0 first:

```bash
python run_pipeline.py --dataset practice_20aa --input_format tsv --steps 00,01,02,03,04,05,06,07,08,09 --analysis_mode both --fva_scope all
```

Run v1.1 extension only:

```bash
python run_pipeline.py --dataset practice_20aa --steps 10,11,12,13,14
```

Run v1.0 plus v1.1 together:

```bash
python run_pipeline.py --dataset practice_20aa --input_format tsv --steps 00,01,02,03,04,05,06,07,08,09,10,11,12,13,14 --analysis_mode both --fva_scope all
```

Run individual v1.1 steps:

```bash
python scripts/steps/10_map_to_chompact.py --dataset practice_20aa
python scripts/steps/11_score_chompact_pathways.py --dataset practice_20aa
python scripts/steps/12_demand_sensitivity.py --dataset practice_20aa
python scripts/steps/13_rank_pathway_biomarkers.py --dataset practice_20aa
python scripts/steps/14_make_chompact_figures.py --dataset practice_20aa
```

## Risks and Assumptions

- CHOmpact-compatible public map JSON is not assumed to exist; v1.1 uses a conservative CSV category schema.
- FVA non-overlap is a robustness indicator, not proof of causal biology.
- If only one `demand_scale` setting is present, sensitivity output is marked `single_setting_not_assessed`.
- Pathway scores summarize reaction-level iCHO3K outputs; they are not CHOmpact lumped fluxes.
- Measured qMet evidence depends on exchange-rate quality, feed correction, and correct units.
- Clone classification outputs are exploratory and require independent validation before operational use.

## Proposed Commits

1. `Add CHOmpact pathway mapping schema`
2. `Add optional CHOmpact interpretation steps 10-14`
3. `Wire CHOmpact v1.1 extension into pipeline runner`
4. `Document CHOmpact v1.1 implementation and CLI`
