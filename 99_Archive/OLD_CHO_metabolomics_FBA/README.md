# CHO Metabolomics FBA/FVA Pipeline

Windows-friendly proof-of-concept pipeline for CHO fed-batch metabolomics analysis.

## What this pipeline does

- Converts CHO culture/process raw data into analysis-ready tables
- Calculates IVCD-normalized extracellular exchange rates
- Performs feed-corrected exchange-rate recalculation
- Runs interval-wise exchange-constrained FBA/pFBA
- Produces clone-level decision summaries
- Adds pathway-level pFBA/FVA summaries for model-inferred metabolic phenotype comparison
- Includes Escher/static visualization helper scripts

## Project layout

```text
scripts/
  00_setup_project.py
  01_model_check.py
  03_run_fba_fva.py
  03_run_fba_fva_relaxed.py
  08_interactive_runner.py
  09_convert_existing_excel_to_template.py
  10_debug_infeasible_constraints.py
  11_batch_run_intervals.py
  12_plot_batch_results.py
  13_recalculate_feed_corrected_rates.py
  14_make_clone_decision_summary.py
  15_model_pathway_scout.py
  16_batch_pathway_pfba_fva.py
  17_plot_pathway_results.py

templates/
  CHO_FBA_final_input_template.xlsx
  CHO_raw_data_with_feeding_columns_template.xlsx

docs/
  README_KR.md
  COMMANDS_QUICK_START_KR.md
  README_clone_decision_summary_KR.md
  README_pathway_FBA_FVA_upgrade_KR.md
```

## Quick start

```bat
pip install -r requirements_windows.txt
python scripts\00_setup_project.py --base C:\CHO_POC
```

Put raw data here:

```text
C:\CHO_POC\data\metabolomics\raw\CHO_raw_data.xlsx
```

Convert raw Excel:

```bat
python scripts\09_convert_existing_excel_to_template.py --base C:\CHO_POC --input "C:\CHO_POC\data\metabolomics\raw\CHO_raw_data.xlsx" --template "C:\CHO_POC\templates\CHO_FBA_final_input_template.xlsx" --volume 30 --igg-unit mg/L
```

Recalculate feed-corrected rates:

```bat
python scripts\13_recalculate_feed_corrected_rates.py --base C:\CHO_POC --raw-excel "C:\CHO_POC\data\metabolomics\raw\CHO_raw_data.xlsx" --volume 30
```

Run interval-wise batch FBA/pFBA:

```bat
python scripts\11_batch_run_intervals.py --base C:\CHO_POC --rates-file data\metabolomics\processed\exchange_rates_feed_corrected.csv --mode strict_then_relaxed --fva-mode none
```

Plot batch results:

```bat
python scripts\12_plot_batch_results.py --base C:\CHO_POC
```

Make clone decision summary:

```bat
python scripts\14_make_clone_decision_summary.py --base C:\CHO_POC
```

Pathway-level model scan:

```bat
python scripts\15_model_pathway_scout.py --base C:\CHO_POC
```

Pathway-level pFBA/FVA:

```bat
python scripts\16_batch_pathway_pfba_fva.py --base C:\CHO_POC --rates-file data\metabolomics\processed\exchange_rates_feed_corrected.csv --pathway-file results\tables\pathway_reaction_sets.csv --mode relaxed --run-fva
```

Pathway plots:

```bat
python scripts\17_plot_pathway_results.py --base C:\CHO_POC
```

## Interpretation

Recommended wording:

> Feed-corrected extracellular exchange rates were used to constrain the iCHO3K model. Interval-wise pFBA and pathway-focused FVA were performed to compare model-inferred metabolic phenotypes across CHO clones.

Avoid claiming exact intracellular flux reconstruction unless objective, pathway reaction sets, and constraints have been manually curated and validated.
