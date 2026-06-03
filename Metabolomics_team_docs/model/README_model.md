# Model files

이 pipeline은 iCHO3K production model을 사용합니다.

권장 위치:

```text
model/iCHO3K/iCHO3K_cho_prod_generic_unblocked.json
```

GitHub에는 대용량/라이선스 이슈를 피하기 위해 model JSON을 직접 올리지 않을 수 있습니다. 워크스테이션 실행 전 모델 파일을 위 경로에 배치하세요.

## Model interpretation

`iCHO3K_cho_prod_generic_unblocked`는 항체 생산 CHO 세포용 genome-scale metabolic model입니다. 하지만 이 모델 자체가 clone별 titer를 자동 예측하는 것은 아닙니다. Clone 차이는 measured exchange rate, growth/product constraint, optional omics constraints를 통해 부여됩니다.
