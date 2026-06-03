# Metabolomics

CHO producer clone의 spent-media / culture data를 iCHO3K production model에 적용하여 **FBA / pFBA / focused FVA / internal or full FVA** 분석을 수행하는 팀 배포용 파이프라인입니다.

이 repository의 목적은 FBA objective 값으로 IgG titer를 직접 예측하는 것이 아니라, measured exchange phenotype을 기반으로 High / Mother / Low producer의 **central metabolism, lactate handling, PPP/TCA, glutamine/nitrogen metabolism, mAb-related flux, FVA flexibility** 차이를 해석하는 것입니다.

## 핵심 메시지

> FBA/FVA는 clone titer를 직접 맞히는 black-box predictor가 아니라, clone-specific measured exchange phenotype을 iCHO3K production model 위에 올려 feasible metabolic flux state와 pathway flexibility를 해석하는 mechanistic analysis framework입니다.

## Recommended workstation run

인터넷이 안 되는 Linux workstation에서는 TSV input을 권장합니다.

```bash
cd CHO_METABOLOMICS

rm -rf data/processed results logs

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

워크스테이션 성능이 충분하면 full model FVA도 가능합니다.

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
  --fva_scope all \
  --fva_targets group_avg \
  --fva_processes 16
```

## Repository structure

```text
Metabolomics/
├── README.md
├── CHANGELOG.md
├── .gitignore
├── environment/
├── model/
├── practice_data/
├── pipeline/
├── docs/
├── examples/
└── results_template/
```

## Key outputs to check first

1. QC: `results/<dataset>/tables/data_qc_warnings.csv`, `results/<dataset>/REPORT_SUMMARY.md`
2. Input phenotype: `exchange_rates.csv`, `Fig2_rate_heatmap.png`, `Fig3_lac_glc_ratio.png`
3. pFBA flux difference: `Fig10B_central_mab_flux_heatmap_zscore.png`, `Fig11_central_mab_flux_delta.png`
4. FVA range separation: `full_fva_all_high_low_overlap.csv`, `full_fva_all_central_mab_subset.csv`, `Fig21_full_fva_high_low_separation.png`
5. Pathway summary: `pathway_scores.csv`, `pathway_score_high_low_delta.csv`, `Fig15_pathway_scores_fva_overlap.png`
6. Escher: `CHO_focus_core_carbon_map.json` + `focused_escher_flux_high_minus_low_cho.json`

## FVA overlap interpretation

`overlap_fraction`은 High와 Low의 FVA range가 얼마나 겹치는지를 나타냅니다.

- `0`에 가까움: High/Low feasible flux range가 거의 안 겹침 → 강한 metabolic separation
- `1`에 가까움: 가능한 flux range가 거의 동일함 → 해당 reaction으로는 구분 어려움

## Interpretation caution

- `DM_igg_g` objective max가 flat할 수 있으며, 이 값을 titer prediction으로 과해석하지 않습니다.
- measured-demand mode는 실측 qIgG를 demand constraint로 고정한 explanation mode입니다.
- FBA/FVA 결과는 measured intracellular flux가 아니라 constraint-based feasible flux state입니다.
- FVA range가 매우 큰 reaction은 under-constrained cycle 또는 model loop 영향을 받을 수 있으므로 focused central/mAb subset과 overlap analysis를 함께 봅니다.

## Version

Current distribution target: `v1.0.0` — workstation-ready TSV + focused/internal/full FVA pipeline.
