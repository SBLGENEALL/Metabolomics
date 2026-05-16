# Current CHO iCHO3K Workflow

Use these files for the real experiment workflow.

## Main input

- Copy `templates/CHO_actual_experiment_input_template.xlsx` to `inputs/CHO_actual_experiment_input.xlsx`.
  - Fill `Experiment_Input` with spent media concentrations, VCD, viability, titer, and optional 20 amino acids.
  - `inputs/` is git-ignored so real company data is not committed.
- Optional transcriptomics can be copied from `templates/transcriptomics_input_template.csv` to `inputs/transcriptomics_input.csv`.

## Main runner

```powershell
python pipeline/interactive_pipeline.py
```

For explicit group comparison, for example Mother versus High:

```powershell
python pipeline/interactive_pipeline.py --reference-group Mother --compare-group High
```

You can use any `producer_group` labels in the input sheet, such as `Mother`, `Low`, `Moderate`, and `High`.

The runner creates:

- `results/interactive_run/spent_media`
- `results/interactive_run/figures`
- `results/interactive_run/icho3k_inputs`
- `results/interactive_run/fba` if COBRApy is available

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
