# CHO_METABOLOMICS v1.0 — Offline CHO FBA/FVA Pipeline

This repository contains an offline-friendly CHO metabolomics FBA/FVA workflow for clone/group comparison using the iCHO3K production model.

## Scientific scope

This pipeline is **not** a stand-alone mAb/IgG titer predictor.

It uses measured extracellular exchange-rate constraints in the iCHO3K production model to compare:

- pFBA flux states
- FVA feasible flux ranges
- High/Mother/Low clone-group differences
- central metabolism and mAb-related pathway behavior
- genome-scale internal/full FVA separation when requested

Use the following language in reports:

- Step 05: **focused central/mAb pathway FVA**
- Step 06 with `--fva_scope internal`: **internal genome-scale FVA**
- Step 06 with `--fva_scope all`: **full model FVA**

Do not call focused FVA alone “full FVA.”

## Input format

For Linux/workstation use, TSV input is recommended.

Required TSV files for the practice dataset:

```text
data/raw/practice_20aa_tsv/raw_timeseries.tsv
data/raw/practice_20aa_tsv/metabolite_map.tsv
data/raw/practice_20aa_tsv/feed_composition.tsv
```

For your own data, place the same three files here:

```text
data/raw/own_experiment/raw_timeseries.tsv
data/raw/own_experiment/metabolite_map.tsv
data/raw/own_experiment/feed_composition.tsv
```

## v1.0 step order

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

## Fast standard workflow

Use this for routine clone comparison, pathway interpretation, figures, and report generation without full/genome-scale FVA.

```bash
python run_pipeline.py \
  --dataset practice_20aa \
  --input_format tsv \
  --rate_days 7,10 \
  --analysis_mode both \
  --constraint_policy production_relaxed \
  --feed_volume_mode interval \
  --demand_scale auto \
  --steps 00,01,02,03,04,05,07,08,09
```

## Workstation full-model FVA workflow

Use this on a Linux/workstation when full model FVA is required. This is the validated v1.0 workstation test route for `practice_20aa`.

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

If full model FVA is too slow, use internal FVA instead:

```bash
--fva_scope internal
```

If the run is still too slow, reduce the target set or process count:

```bash
--fva_targets representative
--fva_processes 8
```

## FVA scope definitions

```text
exchange  = uptake/secretion reactions only
focused   = curated central metabolism + mAb pathway panel
internal  = all non-exchange internal reactions
all       = every reaction in iCHO3K
```

## v1.0 figure numbering

```text
Fig1   IgG time course
Fig2   Exchange-rate heatmap
Fig3   Lactate/glucose phenotype
Fig4   High vs Low exchange-rate comparison
Fig5   Central/mAb flux heatmap
Fig6   Central/mAb flux z-score heatmap
Fig7   Central/mAb High-Low flux delta
Fig8   Pathway scores + FVA overlap
Fig9   Full/internal/all FVA High-Low separation
Fig10  Full/internal/all FVA range delta
Fig11  Summary panel
```

## Key output files

Focused analysis:

```text
results/<dataset>/tables/central_mab_reaction_panel.csv
results/<dataset>/tables/central_mab_flux_state.csv
results/<dataset>/tables/central_mab_fva_report.csv
results/<dataset>/figures/Fig5_central_mab_flux_heatmap.png
results/<dataset>/figures/Fig6_central_mab_flux_zscore.png
results/<dataset>/figures/Fig7_central_mab_flux_delta.png
```

Full model FVA:

```text
results/<dataset>/tables/full_fva/full_fva_all_summary.csv
results/<dataset>/tables/full_fva/full_fva_all_combined_report.csv
results/<dataset>/tables/full_fva/full_fva_all_high_low_overlap.csv
results/<dataset>/tables/full_fva/full_fva_all_central_mab_subset.csv
results/<dataset>/figures/Fig9_full_fva_high_low_separation.png
results/<dataset>/figures/Fig10_full_fva_range_delta.png
```

Run summary:

```text
results/<dataset>/REPORT_SUMMARY.md
results/<dataset>/figures/Fig11_summary_panel.png
```

## Legacy step-number mapping

The previous development branch used non-sequential step numbers. v1.0 uses only `00` through `09`.

```text
Legacy → v1.0
17     → 00_data_qc
01     → 01_load_and_rate
02     → 02_map_metabolites
03     → 03_run_pfba
10     → 04_export_escher_flux
14/16/18 → 05_focused_fva
21     → 06_full_fva
19     → 07_pathway_scores
06     → 08_make_figures
20     → 09_generate_report
```

Legacy numbering is retained only for historical reference.

## Interpretation caution

FVA ranges are model-based feasible ranges, not directly measured intracellular fluxes. Very large ranges can reflect under-constrained reactions or thermodynamic loops. Prioritize focused pathway reactions, robust High/Low separation, and consistency with measured exchange phenotypes.
