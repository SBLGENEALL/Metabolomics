# CHO_METABOLOMICS — Workstation Full/Internal FVA + TSV Input

This is the final offline-friendly CHO FBA/FVA pipeline.

## Scientific scope

This pipeline is **not** a stand-alone mAb/IgG titer predictor. It uses the iCHO3K production model as a mechanistic layer to compare clone/group flux states and feasible flux ranges under measured exchange-rate constraints.

Use the following language in reports:

- Step 14: **focused central/mAb pathway FVA**
- Step 21 with `--fva_scope internal`: **internal genome-scale FVA**
- Step 21 with `--fva_scope all`: **full model FVA**

Do not call step 14 alone “full FVA.”

## Offline/Linux input format

For workstations where Excel is inconvenient or DRM-wrapped, use TSV files instead of XLSX.

Required TSV files:

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

Then run with:

```bash
--input_format tsv
```

## Fast focused FBA/FVA workflow

Use this for quick development, figures, focused Escher maps, and presentation outputs:

```bash
python run_pipeline.py \
  --dataset practice_20aa \
  --input_format tsv \
  --rate_days 7,10 \
  --analysis_mode both \
  --constraint_policy production_relaxed \
  --feed_volume_mode interval \
  --demand_scale auto \
  --steps 17,1,2,3,10,14,16,18,19,20,6
```

## Workstation-scale internal FVA workflow

Use this on Linux/workstation. This runs internal FVA across all non-exchange reactions for group averages first.

```bash
python run_pipeline.py \
  --dataset practice_20aa \
  --input_format tsv \
  --rate_days 7,10 \
  --analysis_mode both \
  --constraint_policy production_relaxed \
  --feed_volume_mode interval \
  --demand_scale auto \
  --steps 17,1,2,3,10,14,16,21,19,20,6 \
  --fva_scope internal \
  --fva_targets group_avg \
  --fva_processes 16
```

If the workstation is very strong, you may try full model FVA:

```bash
--fva_scope all
```

If it is too slow, reduce the target set:

```bash
--fva_targets representative
```

or reduce processes:

```bash
--fva_processes 8
```

## FVA scope definitions

```text
exchange  = uptake/secretion reactions only
focused   = curated central metabolism + mAb pathway panel
internal  = all non-exchange internal reactions
all       = every reaction in iCHO3K
```

## Key output files

Focused analysis:

```text
results/<dataset>/tables/central_mab_fva_report.csv
results/<dataset>/figures/Fig10_central_mab_flux_heatmap.png
results/<dataset>/figures/Fig10B_central_mab_flux_zscore.png
results/<dataset>/figures/Fig12_focused_core_fva_range.png
results/<dataset>/escher_maps/focused/CHO_focus_core_carbon_map.json
```

Internal/full FVA:

```text
results/<dataset>/tables/full_fva/full_fva_internal_summary.csv
results/<dataset>/tables/full_fva/full_fva_internal_combined_report.csv
results/<dataset>/tables/full_fva/full_fva_internal_high_low_overlap.csv
results/<dataset>/tables/full_fva/full_fva_internal_central_mab_subset.csv
results/<dataset>/figures/Fig21_full_fva_high_low_separation.png
results/<dataset>/figures/Fig21B_full_fva_range_delta.png
```

## Interpretation caution

FVA ranges are model-based feasible ranges, not directly measured intracellular fluxes. Very large ranges can reflect under-constrained reactions or thermodynamic loops. Prioritize focused pathway reactions, robust High/Low separation, and consistency with measured exchange phenotypes.
