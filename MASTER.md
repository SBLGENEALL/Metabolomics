# Metabolomics MASTER

This document is the single source of truth for the CHO Metabolomics project. Older chat history and legacy step numbers are superseded by this file.

---

## 0. Version policy

Current stable release:

```text
Metabolomics_v1.0
```

Stable branch:

```text
main
```

Current development branch:

```text
v1.1
```

Current development target:

```text
Metabolomics_v1.1
```

Version meaning:

```text
v1.0
  = validated iCHO3K production workflow
  = cleanup-step-order work merged into main
  = fixed by tag Metabolomics_v1.0

v1.1
  = CHOmpact interpretation integration release
  = optional steps 10-14 added after the validated 00-09 workflow
  = current branch for bugfixes, figure refinement, and interpretation polishing

v1.1.x / v1.1-pr*
  = v1.1-level fixes and polish only
  = examples: logging cleanup, figure title/style fixes, score normalization, report wording
  = merge back into v1.1 after local validation

v1.2+
  = future larger feature releases
  = examples: expanded pathway ontology, glycosylation module, nucleotide/lipid modules, ER folding/secretion module, multi-omics integration, CHOmpact network map redesign
```

Branch interpretation:

```text
main
  = latest validated stable branch
  = currently equivalent to Metabolomics_v1.0

v1.1
  = active development/release-candidate branch for CHOmpact integration

backup_before_archive_reorg
  = historical backup branch
```

Historical branches:

```text
feature/cleanup-step-order
  = completed v1.0 development branch
  = merged into main
  = deleted after v1.0 tag was created

feature/chompact-interpretation
  = original v1.1 development branch name
  = renamed to v1.1
```

---

## 1. Project purpose

The CHO Metabolomics pipeline converts fed-batch culture and spent-media metabolomics data into exchange-rate-constrained iCHO3K FBA/FVA outputs.

The pipeline is designed for:

- clone/group comparison
- productivity-group interpretation
- measured exchange-rate constraint generation
- pFBA flux-state comparison
- focused central/mAb pathway FVA
- full/internal genome-scale FVA
- pathway-level interpretation
- report and figure generation

This pipeline is **not** a stand-alone IgG/mAb titer predictor.

Actual clone grouping is intentionally deferred until real clone data are available. The code should support two-group, three-group, or multi-group comparisons rather than hard-coding High/Mother/Low assumptions.

---

## 2. Core modeling strategy

### iCHO3K

Role:

- calculation engine
- COBRApy FBA/pFBA/FVA
- exchange-rate-constrained clone comparison
- full/internal genome-scale FVA

Interpretation:

- iCHO3K reaction-level outputs are the quantitative source of truth.

### CHOmpact

Role in v1.1:

- pathway vocabulary
- biological grouping
- figure/report interpretation layer
- CHO-specific mechanistic framing

Important principle:

- CHOmpact must not replace iCHO3K as the calculation model.
- CHOmpact must not be used to recompute fluxes.
- CHOmpact lumped fluxes must not be reconstructed by summing iCHO3K reactions.
- CHOmpact should summarize, label, and visualize iCHO3K outputs only.

---

## 3. Environment

Validated v1.0 environment:

- Linux workstation
- Python 3.13 environment: `metabolomics_env_py313`
- COBRApy 0.30.0
- GLPK solver
- offline-compatible TSV input

The v1.0 report generator does not require the optional `tabulate` package.

v1.1 steps 10-14 are interpretation-layer scripts that read v1.0 CSV outputs and should not require COBRApy unless the full 00-09 calculation pipeline is rerun.

---

## 4. v1.0 validated workflow

The v1.0 pipeline uses clean sequential steps `00` through `09`.

```text
00_data_qc.py
01_load_and_rate.py
02_map_metabolites.py
03_run_pfba.py
04_export_escher_flux.py
05_focused_fva.py
06_full_fva.py
07_pathway_scores.py
08_make_figures.py
09_generate_report.py
```

Legacy step numbers such as 17, 10, 14, 21, 19, 20, and 6 are deprecated.

Validated full-model FVA command:

```bash
python run_pipeline.py \
  --dataset practice_20aa \
  --input_format tsv \
  --rate_days 7,10 \
  --analysis_mode both \
  --constraint_policy production_relaxed \
  --feed_volume_mode interval \
  --demand_scale auto \
  --steps 00,01,02,03,04,05,06,07,08,09 \
  --fva_scope all \
  --fva_targets group_avg \
  --fva_processes 16
```

Fast route without full FVA:

```bash
python run_pipeline.py \
  --dataset practice_20aa \
  --input_format tsv \
  --steps 00,01,02,03,04,05,07,08,09
```

---

## 5. v1.1 workflow

The v1.1 branch adds optional interpretation steps after v1.0.

```text
10_map_to_chompact.py
11_score_chompact_pathways.py
12_demand_sensitivity.py
13_rank_pathway_biomarkers.py
14_make_chompact_figures.py
```

Run v1.1 extension after v1.0 outputs already exist:

```bash
python run_pipeline.py --dataset practice_20aa --steps 10,11,12,13,14
```

Run v1.0 plus v1.1 together:

```bash
python run_pipeline.py \
  --dataset practice_20aa \
  --input_format tsv \
  --steps 00,01,02,03,04,05,06,07,08,09,10,11,12,13,14 \
  --analysis_mode both \
  --fva_scope all
```

Note: `--fva_scope all` in step 06 is computationally expensive. For quick v1.1 interpretation-layer smoke tests, run steps 10-14 only after existing v1.0 tables are available.

---

## 6. Input data structure

Practice TSV input:

```text
data/raw/practice_20aa_tsv/raw_timeseries.tsv
data/raw/practice_20aa_tsv/metabolite_map.tsv
data/raw/practice_20aa_tsv/feed_composition.tsv
```

Own experiment TSV input:

```text
data/raw/own_experiment/raw_timeseries.tsv
data/raw/own_experiment/metabolite_map.tsv
data/raw/own_experiment/feed_composition.tsv
```

---

## 7. v1.0 outputs

Main figures:

```text
Fig1_IgG_timecourse.png
Fig2_rate_heatmap.png
Fig3_lac_glc_ratio.png
Fig4_high_vs_low_rates.png
Fig5_central_mab_flux_heatmap.png
Fig6_central_mab_flux_zscore.png
Fig7_central_mab_flux_delta.png
Fig8_pathway_scores_fva_overlap.png
Fig9_full_fva_high_low_separation.png
Fig10_full_fva_range_delta.png
Fig11_summary_panel.png
```

Supplementary figure:

```text
SuppFig1_data_qc_overview.png
```

Main report:

```text
results/<dataset>/REPORT_SUMMARY.md
```

---

## 8. v1.1 outputs

CHOmpact interpretation tables:

```text
results/<dataset>/tables/chompact/
```

Expected v1.1 outputs include:

```text
chompact_mapping_qc.csv
chompact_mapped_flux.csv
chompact_mapped_fva.csv
chompact_mapped_fva_overlap.csv
chompact_mapped_measured_rates.csv
chompact_pathway_fba_activity_scores.csv
chompact_pathway_fva_robustness_scores.csv
chompact_pathway_measured_qmet_scores.csv
chompact_pathway_high_low_separation.csv
chompact_demand_scale_sensitivity.csv
chompact_robust_rank_across_demand_scales.csv
chompact_ranked_pathway_biomarkers.csv
chompact_top_pathway_biomarkers_for_report.csv
CHOmpact_v1_1_executive_summary.md
```

v1.1 figures:

```text
Fig12_chompact_pathway_activity.png
Fig13_chompact_fva_robustness.png
Fig14_chompact_biomarker_ranking.png
```

---

## 9. Interpretation rules

### Measured values

Measured or directly derived from experimental data:

- IgG titer
- VCD
- viability
- feed-corrected qMet / uptake / secretion rates
- glucose, lactate, ammonia, amino acid exchange rates

### Model-constrained values

Measured rates imposed as model constraints:

- exchange bounds derived from measured qMet
- product demand constraints derived from measured IgG/qP

### Model-predicted values

Model-derived:

- pFBA flux values
- FVA minimum/maximum ranges
- internal reaction fluxes
- model objective values
- pathway scores derived from iCHO3K outputs

Reports must clearly distinguish measured evidence from constrained model inputs and predicted model interpretation.

### FVA caution

FVA ranges are model-based feasible ranges, not direct intracellular flux measurements and not statistical error bars. Very broad FVA intervals may reflect under-constrained reactions, alternate optima, or thermodynamic loops.

### Constraint-driven vs emergent interpretation

If a High-vs-Low difference is directly imposed as a boundary constraint, it is constraint-driven and should be reported descriptively. Interior differences that persist after boundary equalization or demand-attribution checks may be treated as emergent model hypotheses.

---

## 10. v1.1 release criteria

Metabolomics v1.1 should not be merged to `main` or tagged until:

1. v1.0 steps 00-09 still pass on `practice_20aa`.
2. v1.1 steps 10-14 pass after valid v1.0 outputs are present.
3. CHOmpact outputs are non-empty and interpretable.
4. `chompact_mapping_qc.csv` confirms acceptable mapping coverage.
5. Figures 12-14 and the executive summary are generated correctly.
6. FVA failure/retry logging is clarified enough not to confuse prior failed attempts with final success.
7. `README.md`, `MASTER.md`, and `CHANGELOG.md` are updated.

The intended v1.1 tag will be:

```text
Metabolomics_v1.1
```

---

## 11. Future v1.2 candidates

The following should be deferred to v1.2 or later unless the v1.1 review explicitly decides otherwise:

- expanded CHOmpact/iCHO3K pathway ontology
- glycosylation and nucleotide-sugar donor module expansion
- nucleotide and lipid metabolism interpretation modules
- ER folding/secretion burden module
- publication-grade CHOmpact network map
- multi-omics integration
