
# Clone decision summary add-on

## 넣는 위치

```text
C:\CHO_POC_ChatGPT_260519\scripts\14_make_clone_decision_summary.py
```

## 실행

```bat
cd C:\CHO_POC_ChatGPT_260519
python scripts\14_make_clone_decision_summary.py --base C:\CHO_POC_ChatGPT_260519
```

## 출력

```text
results\tables\clone_decision_summary.csv
results\figures\clone_decision_score.png
results\figures\final_titer_vs_lactate.png
```

## 해석

이 스크립트는 단순 FBA objective가 아니라 아래를 같이 봅니다.

- final titer
- mean qP
- final viability
- mean lactate secretion
- mean ammonia secretion
- glucose uptake burden

주의: decision score는 선별 보조용 composite score입니다. 실제 clone 선정은 raw titer/qP/viability/feed-corrected flux를 함께 보고 결정해야 합니다.
