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

The runner creates:

- `results/interactive_run/spent_media`
- `results/interactive_run/figures`
- `results/interactive_run/icho3k_inputs`
- `results/interactive_run/fba` if COBRApy is available

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
