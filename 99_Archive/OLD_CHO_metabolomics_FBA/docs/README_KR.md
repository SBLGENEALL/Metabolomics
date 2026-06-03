# CHO metabolomics FBA/FVA final all-in-one package

## 0. 설치

```bat
pip install -r requirements_windows.txt
```

또는 conda 환경에서 이미 cobra, pandas, matplotlib, openpyxl, escher가 있으면 그대로 사용하면 됩니다.

## 1. 프로젝트 폴더 준비

```bat
cd C:\CHO_POC
python scripts\00_setup_project.py --base C:\CHO_POC
```

## 2. 원본 데이터 입력 양식

`templates\CHO_raw_data_with_feeding_columns_template.xlsx`의 `Paste_Raw_Data` 시트에 붙여넣으면 됩니다.

앞쪽 컬럼은 사용자 원본 형식입니다.

```text
DAY
Sample ID
Gln
Glu
Gluc
Lac
NH4+
Na+
K+
Ca++
pH
PO2
PCO2
Osm
Total Desity
Viable Density
Viability
Average Live Diameter
IgG
```

뒤쪽 feed correction 컬럼:

```text
Culture Volume mL
Sample Removed mL
Glucose Feed mL
Glucose Feed Conc mM
Feed4 mL
CellBoost mL
Feed4 Glucose mM
Feed4 Gln mM
Feed4 Glu mM
Feed4 Lac mM
Feed4 NH4 mM
CellBoost Glucose mM
CellBoost Gln mM
CellBoost Glu mM
CellBoost Lac mM
CellBoost NH4 mM
Feed Note
```

각 Day row의 feed 값은 **이전 sampling 이후 현재 Day까지 넣은 양**으로 입력합니다.

예: Day7 row의 `Glucose Feed mL = 0.3`, `Glucose Feed Conc mM = 1000`이면 Day5~Day7 사이에 0.3 mL 1000 mM glucose를 넣었다는 뜻입니다.

## 3. 기존 raw Excel 변환

```bat
python scripts\09_convert_existing_excel_to_template.py --base C:\CHO_POC --input "C:\CHO_POC\data\metabolomics\raw\CHO_raw_data.xlsx" --template "C:\CHO_POC\templates\CHO_FBA_final_input_template.xlsx" --volume 30 --igg-unit mg/L
```

생성:

```text
data\metabolomics\processed\converted_process_data.csv
data\metabolomics\processed\converted_extracellular_metabolomics.csv
data\metabolomics\processed\exchange_rates_from_user_format.csv
```

## 4. feed correction rate 생성

```bat
python scripts\13_recalculate_feed_corrected_rates.py --base C:\CHO_POC --raw-excel "C:\CHO_POC\data\metabolomics\raw\CHO_raw_data.xlsx" --volume 30
```

생성:

```text
data\metabolomics\processed\exchange_rates_feed_corrected.csv
data\metabolomics\processed\feed_events_parsed.csv
data\metabolomics\processed\feed_column_detection_debug.csv
```

## 5. 전체 clone × 전체 interval batch FBA

feed correction 전 apparent rate:

```bat
python scripts\11_batch_run_intervals.py --base C:\CHO_POC --rates-file data\metabolomics\processed\exchange_rates_from_user_format.csv --mode strict_then_relaxed --fva-mode none
```

feed correction 후:

```bat
python scripts\11_batch_run_intervals.py --base C:\CHO_POC --rates-file data\metabolomics\processed\exchange_rates_feed_corrected.csv --mode strict_then_relaxed --fva-mode none
```

## 6. 결과 plot

```bat
python scripts\12_plot_batch_results.py --base C:\CHO_POC
```

결과:

```text
results\tables\batch\batch_objective_summary.csv
results\tables\batch\batch_key_exchange_summary.csv
results\tables\batch\batch_clone_level_summary.csv

results\figures\batch_objective_by_interval.png
results\figures\batch_EX_glc_e_by_interval.png
results\figures\batch_EX_lac_L_e_by_interval.png
results\figures\batch_EX_gln_L_e_by_interval.png
results\figures\batch_EX_glu_L_e_by_interval.png
results\figures\batch_EX_nh4_e_by_interval.png
```

## 7. interactive 실행

```bat
python scripts\08_interactive_runner.py
```

또는 `run_interactive.bat` 더블클릭.

## 8. 해석 주의

- 14일 배양 전체를 하나의 FBA로 보지 말고, interval별 FBA를 batch로 돌린 뒤 clone별 trend로 요약합니다.
- feed correction을 하지 않은 rate는 true uptake/secretion이 아니라 apparent net rate입니다.
- Feed4/CellBoost 조성을 모르면 해당 feed 내 amino acid correction은 불완전합니다.
- strict가 infeasible이고 relaxed가 optimal이면 입력 constraint가 너무 좁거나 feed/dilution 보정이 필요하다는 뜻입니다.
