
# Pathway-level FBA/FVA upgrade add-on

이 add-on은 기존 pipeline을 다음 수준으로 올립니다.

기존:
- extracellular metabolite 증감
- feed-corrected exchange rate
- FBA feasibility / pFBA QC

추가:
- model reaction scan
- objective candidate scan
- pathway-level internal reaction set 생성
- pathway-level pFBA flux summary
- pathway-level FVA flexibility summary
- clone별 pathway phenotype 비교

## 파일 위치

아래 파일들을 넣으세요.

```text
C:\CHO_POC_ChatGPT_260519\scripts\15_model_pathway_scout.py
C:\CHO_POC_ChatGPT_260519\scripts\16_batch_pathway_pfba_fva.py
C:\CHO_POC_ChatGPT_260519\scripts\17_plot_pathway_results.py
```

## 1. 모델 내부 reaction scan

```bat
cd C:\CHO_POC_ChatGPT_260519

python scripts\15_model_pathway_scout.py --base C:\CHO_POC_ChatGPT_260519
```

모델이 기본 위치가 아니면:

```bat
python scripts\15_model_pathway_scout.py --base C:\CHO_POC_ChatGPT_260519 --model-file "C:\경로\iCHO3K_cho_prod_generic_unblocked.json"
```

생성:

```text
results\tables\model_reaction_scout.csv
results\tables\objective_candidates.csv
results\tables\pathway_reaction_sets.csv
```

중요:
- `objective_candidates.csv`를 열어서 biomass/product/IgG 관련 objective 후보를 확인하세요.
- `pathway_reaction_sets.csv`는 자동 keyword 기반이므로, 필요하면 include 컬럼을 TRUE/FALSE로 손봐도 됩니다.

## 2. Pathway-level pFBA/FVA 실행

기본은 objective가 불확실하므로 `relaxed` + `fraction 0.0`으로 시작합니다.

```bat
python scripts\16_batch_pathway_pfba_fva.py --base C:\CHO_POC_ChatGPT_260519 --rates-file data\metabolomics\processed\exchange_rates_feed_corrected.csv --pathway-file results\tables\pathway_reaction_sets.csv --mode relaxed --run-fva
```

objective 후보를 골랐다면 예를 들어:

```bat
python scripts\16_batch_pathway_pfba_fva.py --base C:\CHO_POC_ChatGPT_260519 --rates-file data\metabolomics\processed\exchange_rates_feed_corrected.csv --pathway-file results\tables\pathway_reaction_sets.csv --mode strict_then_relaxed --objective biomass_cho_prod --fraction 0.9 --run-fva
```

## 3. Pathway 결과 plot

```bat
python scripts\17_plot_pathway_results.py --base C:\CHO_POC_ChatGPT_260519
```

생성:

```text
results\tables\pathway\pathway_interval_summary.csv
results\tables\pathway\pathway_reaction_fluxes.csv
results\tables\pathway\pathway_reaction_fva.csv
results\tables\pathway\pathway_fva_summary.csv
results\tables\pathway\pathway_clone_summary.csv

results\figures\pathway_pfba_heatmap.png
results\figures\pathway_fva_flexibility_heatmap.png
results\figures\pathway_<pathway>_by_interval.png
```

## 해석 포인트

이제 결과를 세 층으로 나눠서 말할 수 있습니다.

### A. Observed / feed-corrected exchange phenotype
- glucose uptake
- lactate secretion/reuptake
- glutamine uptake
- ammonia secretion
- titer/qP/viability

### B. Model feasibility
- strict optimal 여부
- relaxed 필요 여부
- warning 수

### C. Model-inferred pathway phenotype
- glycolysis pFBA activity
- TCA pFBA activity
- lactate/pyruvate pathway activity
- glutamine/glutamate pathway activity
- ammonia/nitrogen pathway activity
- FVA flexibility range

보고서 표현 예:

```text
Feed-corrected extracellular exchange rates were used to constrain the iCHO3K model.
Interval-wise pFBA and pathway-focused FVA were then used to compare model-inferred metabolic phenotypes across clones.
```

주의:
자동 reaction set은 keyword 기반이므로, 논문화하려면 `pathway_reaction_sets.csv`를 열어서 reaction 이름을 검토하고 include 컬럼을 수동 curate하는 것을 추천합니다.
