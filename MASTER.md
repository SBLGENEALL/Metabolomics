# Metabolomics MASTER

Current project baseline: **Metabolomics v1.0**  
Validated branch: `feature/cleanup-step-order`  
Validated dataset: `practice_20aa`  
Validation status: passed on Linux/workstation with full FVA (`--fva_scope all`).

This document is the single source of truth for the current Metabolomics pipeline. Older chat history and legacy step numbers are superseded by this file.

---

## 1. Project purpose

The CHO Metabolomics pipeline converts fed-batch culture and spent-media metabolomics data into exchange-rate-constrained iCHO3K FBA/FVA outputs.

The pipeline is designed for:

- clone/group comparison
- High/Mother/Low productivity interpretation
- measured exchange-rate constraint generation
- pFBA flux-state comparison
- focused central/mAb pathway FVA
- full/internal genome-scale FVA
- report and figure generation

This pipeline is **not** a stand-alone IgG/mAb titer predictor.

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

Current v1.0 status:

- not yet implemented as code
- planned for v1.1 as an interpretation/visualization layer

Planned v1.1 role:

- pathway vocabulary
- biological grouping
- figure layout logic
- CHO-specific pathway interpretation

Important principle:

- CHOmpact should not replace iCHO3K as the calculation model.
- CHOmpact should be used to summarize and visualize iCHO3K outputs.

---

## 3. Environment

Validated environment:

- Linux workstation
- Python 3.13 environment: `metabolomics_env_py313`
- COBRApy 0.30.0
- GLPK solver
- offline-compatible TSV input

The v1.0 report generator no longer requires the optional `tabulate` package.

---

## 4. v1.0 step structure

The v1.0 pipeline uses only clean sequential steps `00` through `09`.

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

---

## 5. Validated v1.0 command

Validated full-model FVA route:

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

## 7. Current outputs

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

## 8. Interpretation rules

### Measured values

Measured or directly derived from experimental data:

- IgG titer
- VCD
- viability
- feed-corrected qMet / uptake / secretion rates
- glucose, lactate, ammonia, amino acid exchange rates

### Model-predicted values

Model-derived:

- pFBA flux values
- FVA minimum/maximum ranges
- internal reaction fluxes
- model objective values
- pathway scores derived from iCHO3K outputs

Reports must clearly distinguish measured evidence from model-predicted interpretation.

### FVA caution

FVA ranges are model-based feasible ranges, not direct intracellular flux measurements. Very broad FVA intervals may reflect under-constrained reactions, alternate optima, or thermodynamic loops.

---

## 9. v1.1 direction

The next development branch should be:

```text
feature/chompact-interpretation
```

v1.1 goals:

1. build `data/chompact_pathway_mapping.csv`
2. map iCHO3K reaction-level outputs to CHOmpact-style pathway categories
3. score pathways using measured qMet separation, pFBA flux separation, and FVA robustness
4. run demand-scale sensitivity analysis
5. redesign figures and report around biological pathway interpretation

v1.1 must preserve the v1.0 calculation engine and should not replace iCHO3K with CHOmpact.
