# Current CHO iCHO3K Workflow

Use these files for the real experiment workflow.

## Main input

- Copy `templates/CHO_actual_experiment_input_template.xlsx` to `inputs/CHO_actual_experiment_input.xlsx`.
  - Fill `Experiment_Input` with spent media concentrations, VCD, viability, titer, and optional 20 amino acids.
  - `inputs/` is git-ignored so real company data is not committed.
- The pipeline also accepts the raw feeding workbook format with a `Paste_Raw_Data`,
  `CHO_raw_data`, or similarly structured raw-data sheet.
  - `DAY` values like `Day5` are converted to numeric day values.
  - `Sample ID` is used as the clone/condition name.
  - `Gln`, `Glu`, `Gluc`, `Lac`, and `NH4+` are mapped to iCHO3K-ready metabolite names.
  - Feeding columns are used to calculate `feed_corrected_rate_per_day` and
    `rate_for_model_per_day`, so added glucose/feed is not mistaken for cellular production.
- Optional transcriptomics can be copied from `templates/transcriptomics_input_template.csv` to `inputs/transcriptomics_input.csv`.

## Main runner

```powershell
python pipeline/interactive_pipeline.py
```

To run directly on the raw feeding workbook:

```powershell
python pipeline/interactive_pipeline.py --input "C:\Users\j3das\Downloads\CHO_raw_data_with_feeding_columns_template.xlsx"
```

To create a POC workbook with simulated 20 amino acid spent-media columns from
an existing raw workbook:

```powershell
python pipeline/make_poc_raw_data_with_amino_acids.py --input "C:\path\to\CHO_raw_data.xlsx" --output inputs/CHO_raw_data_with_20AA_POC.xlsx
python pipeline/interactive_pipeline.py --input inputs/CHO_raw_data_with_20AA_POC.xlsx --skip-fba
```

The raw parser recognizes the core spent-media columns plus 3-letter amino acid
columns such as `Ala`, `Arg`, `Asn`, `Asp`, `Cys`, `Gly`, `His`, `Ile`, `Leu`,
`Lys`, `Met`, `Phe`, `Pro`, `Ser`, `Thr`, `Trp`, `Tyr`, and `Val`.

For real fed-batch workbooks with many clones/timepoints, the interactive runner
uses compact FBA flux output by default. To save every iCHO3K reaction flux, add:

```powershell
python pipeline/interactive_pipeline.py --input "C:\path\to\CHO_raw_data.xlsx" --fba-flux-scope all
```

For explicit group comparison, for example Mother versus High:

```powershell
python pipeline/interactive_pipeline.py --reference-group Mother --compare-group High
```

If FVA is too slow on a large company dataset, run FBA first and skip FVA:

```powershell
python pipeline/interactive_pipeline.py --skip-fva
```

For broader model-level FVA on the top active reactions that carry nonzero FBA flux, use:

```powershell
python pipeline/interactive_pipeline.py --reference-group Mother --compare-group High --fva-scope active
```

To calculate full-model FVA for every iCHO3K reaction, use the slower all-reaction mode:

```powershell
python pipeline/interactive_pipeline.py --reference-group Mother --compare-group High --fva-scope all
```

You can use any `producer_group` labels in the input sheet, such as `Mother`, `Low`, `Moderate`, and `High`.

The runner creates:

- `results/interactive_run/spent_media`
- `results/interactive_run/figures`
- `results/interactive_run/icho3k_inputs`
- `results/interactive_run/fba` if COBRApy is available

FBA/FVA output includes:

- `fba_objective_results_by_scenario.csv`: predicted growth and IgG objective values for each clone/group/time interval
- `fba_selected_fluxes_by_scenario.csv`: model-inferred fluxes for core exchange, biomass, and IgG reactions
- `fba_all_reaction_fluxes_by_scenario.csv`: full FBA solution for every iCHO3K reaction and objective
- `fba_flux_differences_by_group.csv`: model-predicted flux differences between the requested reference and comparison groups
- `fba_flux_design_candidates.csv`: high-priority flux-difference candidates for media/feed or pathway follow-up
- `fva_selected_reactions_by_scenario.csv`: flux variability ranges for selected reactions, showing which fluxes are tightly constrained versus flexible
- `fva_active_reactions_by_scenario.csv`: FVA for reactions with nonzero FBA flux when `--fva-scope active` is used
- `fva_all_reactions_by_scenario.csv`: full-model FVA output when `--fva-scope all` or `--full-fva` is used
- `fva_reaction_range_summary.csv`: average FVA flexibility/tightness by reaction
- `fba_fva_report.html`: browser-friendly summary tables

## Figures

The pipeline writes publication-style SVG files to `results/interactive_run/figures`:

- `figure_A_time_course.svg`: spent-media time-course profiles
- `figure_B_endpoint_contrast.svg`: selected group contrast
- `figure_C_log2_change_heatmap.svg`: mean log2 change heatmap
- `figure_D_qmet_heatmap.svg`: clone-level uptake/secretion heatmap
- `figure_E_dominant_exchange_fluxes.svg`: dominant qMet bar plot
- `figure_F_culture_profile.svg`: VCD, viability, and titer profiles
- `figure_report.html`: scrollable HTML report that is easier to inspect when
  SVG labels are dense

SVG files can be opened in a browser, PowerPoint, Illustrator, or Inkscape.

## Escher Overlay

The model-input step writes Escher-compatible reaction data to
`results/interactive_run/icho3k_inputs`:

- `escher_reaction_data_mean_flux.json`
- `escher_reaction_data_mean_flux.csv`
- `escher_overlay_instructions.html`
- `escher_flux_overlay.html` if COBRApy and Escher are installed
- `escher_by_clone/*_reaction_data.json`: clone-specific Escher reaction data

Use the JSON as reaction data in Escher Builder with an iCHO3K-compatible map.
The sign convention follows COBRA exchange flux: uptake is negative and
secretion is positive.

If you have an iCHO3K Escher map JSON, pass it directly:

```powershell
python pipeline/build_escher_overlay.py --map-json path/to/icho3k_map.json --reaction-data results/interactive_run/icho3k_inputs/escher_by_clone/CloneA_reaction_data.json --output results/interactive_run/icho3k_inputs/CloneA_escher.html
```

## iCHO3K model

The workflow uses the local production model:

```text
models\iCHO3K\Model\iCHO3K_cho_prod_generic_unblocked.json
```

Verified objectives:

- Growth: `biomass_cho_prod`
- Product: `DM_igg_g`
- IgG assembly: `igg_formation`

## Needed data

Minimum:

- VCD, viability, titer
- glucose, lactate, ammonia
- glutamine, glutamate, alanine

Recommended:

- 20 amino acids
- RNA-seq TPM table for high vs low producer or phase comparison
