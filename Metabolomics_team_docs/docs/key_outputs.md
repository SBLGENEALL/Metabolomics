# Key outputs

## 1. QC

```text
results/<dataset>/tables/data_qc_warnings.csv
results/<dataset>/REPORT_SUMMARY.md
```

## 2. Exchange phenotype

```text
results/<dataset>/tables/exchange_rates.csv
results/<dataset>/figures/Fig2_rate_heatmap.png
results/<dataset>/figures/Fig3_lac_glc_ratio.png
```

## 3. pFBA flux difference

```text
results/<dataset>/figures/Fig10B_central_mab_flux_heatmap_zscore.png
results/<dataset>/figures/Fig11_central_mab_flux_delta.png
```

## 4. Full/internal FVA

```text
results/<dataset>/tables/full_fva/full_fva_all_high_low_overlap.csv
results/<dataset>/tables/full_fva/full_fva_all_central_mab_subset.csv
results/<dataset>/figures/Fig21_full_fva_high_low_separation.png
results/<dataset>/figures/Fig21B_full_fva_range_delta.png
```

Sort `full_fva_*_high_low_overlap.csv` by:

```text
overlap_fraction ascending
```

Then inspect central metabolism, lactate, PPP, TCA, glutamine/nitrogen, energy, and mAb-related reactions.

## 5. Pathway score

```text
results/<dataset>/tables/pathway_scores.csv
results/<dataset>/tables/pathway_score_high_low_delta.csv
results/<dataset>/figures/Fig15_pathway_scores_fva_overlap.png
```

## 6. Escher

Map:

```text
results/<dataset>/escher_maps/focused/CHO_focus_core_carbon_map.json
```

Data:

```text
results/<dataset>/tables/focused_escher/focused_escher_flux_high_minus_low_cho.json
```

In Escher:
- Red / positive = higher in HighAvg
- Blue / negative = higher in LowAvg
