# 자주 쓰는 명령어

## 폴더 준비
```bat
python scripts\00_setup_project.py --base C:\CHO_POC
```

## 모델 체크
```bat
python scripts\01_model_check.py --base C:\CHO_POC --model-file C:\CHO_POC\model\iCHO3K-main\iCHO3K\Model\iCHO3K_cho_prod_generic_unblocked.json
```

## raw Excel 변환
```bat
python scripts\09_convert_existing_excel_to_template.py --base C:\CHO_POC --input "C:\CHO_POC\data\metabolomics\raw\CHO_raw_data.xlsx" --template "C:\CHO_POC\templates\CHO_FBA_final_input_template.xlsx" --volume 30 --igg-unit mg/L
```

## feed corrected rates 만들기
```bat
python scripts\13_recalculate_feed_corrected_rates.py --base C:\CHO_POC --raw-excel "C:\CHO_POC\data\metabolomics\raw\CHO_raw_data.xlsx" --volume 30
```

## batch FBA: feed correction 전
```bat
python scripts\11_batch_run_intervals.py --base C:\CHO_POC --rates-file data\metabolomics\processed\exchange_rates_from_user_format.csv --mode strict_then_relaxed --fva-mode none
```

## batch FBA: feed correction 후
```bat
python scripts\11_batch_run_intervals.py --base C:\CHO_POC --rates-file data\metabolomics\processed\exchange_rates_feed_corrected.csv --mode strict_then_relaxed --fva-mode none
```

## plot
```bat
python scripts\12_plot_batch_results.py --base C:\CHO_POC
```

## 한 interval만
```bat
python scripts\03_run_fba_fva.py --base C:\CHO_POC --condition Clone01 --rates-file data\metabolomics\processed\exchange_rates_feed_corrected.csv --day-start 5 --day-end 7 --skip-fva
```

## relaxed single interval
```bat
python scripts\03_run_fba_fva_relaxed.py --base C:\CHO_POC --condition Clone01 --rates-file data\metabolomics\processed\exchange_rates_feed_corrected.csv --day-start 5 --day-end 7 --relax-multiplier 10
```
