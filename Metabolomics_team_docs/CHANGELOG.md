# CHANGELOG

## v1.0.0 — Workstation-ready FBA/FVA pipeline

### Added
- TSV input mode for offline Linux workstation.
- Data QC step for required columns, missing values, negative values, duplicated clone-day rows, feed pattern, and rate sanity checks.
- Focused central/mAb pathway FVA.
- Internal/full FVA option with `--fva_scope exchange|focused|internal|all`.
- Group-average FVA target mode for HighAvg / MotherAvg / LowAvg.
- FVA overlap/separation summary for High vs Low producer groups.
- Pathway-level scores for glycolysis, lactate, PPP, TCA, glutamine/nitrogen, energy, and mAb-related pathway panels.
- Focused Escher map export.
- Automated report generation.

### Removed from main workflow
- KO screening as a main output.
- Product sequence requirement as a default step.
- RECON1 large-map annotation as a default workflow.
- Fig4 direct objective diagnostic from main figure set.

### Interpretation change
- The project is positioned as a mechanistic clone-comparison framework, not a direct titer-prediction model.
