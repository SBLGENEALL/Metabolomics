# Practice data

Practice dataset은 팀 교육과 pipeline sanity check를 위한 synthetic/demo data입니다.

권장 TSV 구조:

```text
data/raw/practice_20aa_tsv/raw_timeseries.tsv
data/raw/practice_20aa_tsv/metabolite_map.tsv
data/raw/practice_20aa_tsv/feed_composition.tsv
```

실제 실험 데이터는 같은 schema로 아래 위치에 넣습니다.

```text
data/raw/own_experiment/raw_timeseries.tsv
data/raw/own_experiment/metabolite_map.tsv
data/raw/own_experiment/feed_composition.tsv
```

필수 데이터:
- Clone or Sample ID
- Group, if available
- Day
- VCD / viable density / viability
- culture volume
- feed volume and sample removal
- IgG titer
- glucose, lactate, NH4
- amino acids, when available

Spent media가 없으면 exchange-constrained FBA/FVA 해석력이 크게 떨어집니다.
