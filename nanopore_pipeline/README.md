# Nanopore (MinION) Reference-Mapping Pipeline

Oxford Nanopore MinION으로 시퀀싱한 결과(barcode별 fastq)를 reference FASTA에
매핑하고, BAM/consensus FASTA/VCF/리포트를 생성하는 파이프라인입니다.

## 핵심 아이디어: 이름 자동화

여러 샘플을 분석하면 결과 파일이 `barcode01`, `barcode02`... 처럼 시퀀서가 붙인
이름으로 나오는데, 실제로는 각각 다른 reference(샘플)에 대응됩니다. 이 파이프
라인은 **`samplesheet.csv`에 barcode ↔ reference ↔ 샘플명을 한 번만 정의**해
두면, 처음부터 모든 출력 파일을 샘플명으로 생성하므로 분석 후에 일일이 이름을
바꿀 필요가 없습니다.

이미 다른 도구(EPI2ME 등)로 분석을 끝낸 결과가 barcode 이름으로 남아있다면,
`scripts/rename_outputs.py`로 한 번에 일괄 이름변경/정리할 수 있습니다.

## 디렉터리 구조

```
nanopore_pipeline/
├── samplesheet.csv      # barcode_id, reference_fasta, sample_name 매핑표
├── config.yaml          # 경로/파라미터 설정
├── run_pipeline.sh       # 메인 파이프라인 (minimap2 -> samtools -> consensus -> variant call)
├── scripts/
│   └── rename_outputs.py # 기존 결과 일괄 이름변경 유틸리티
├── references/           # 샘플별 reference FASTA (예: refA.fasta)
└── data/raw/fastq_pass/  # MinKNOW/Guppy/Dorado의 barcode별 fastq.gz 폴더
    ├── barcode01/
    ├── barcode02/
    └── ...
```

## 사전 준비

1. 필요한 도구 설치 (conda 권장):
   ```bash
   conda install -c bioconda -c conda-forge minimap2 samtools bcftools nanofilt
   ```
2. `data/raw/fastq_pass/`에 barcode별 fastq.gz 폴더를 둔다.
3. `references/`에 샘플별 reference FASTA를 넣는다.
4. `samplesheet.csv`를 실제 매핑에 맞게 수정한다:

   ```csv
   barcode_id,reference_fasta,sample_name
   barcode01,refA.fasta,SampleA
   barcode02,refB.fasta,SampleB
   barcode03,refC.fasta,SampleC
   ```

5. 필요하면 `config.yaml`에서 minimap2 preset, threads, QC 필터링 기준,
   variant caller(`bcftools` 또는 `medaka`)를 조정한다.

## 실행

```bash
cd nanopore_pipeline
./run_pipeline.sh
```

각 샘플마다 `results/<sample_name>/` 폴더에 다음 파일들이 생성됩니다
(전부 샘플명 기준, 추가 이름변경 불필요):

- `<sample_name>.fastq.gz` — 병합된 raw reads (필요시 QC 필터링본 추가)
- `<sample_name>.sorted.bam` / `.bai` — reference 정렬 결과
- `<sample_name>.flagstat.txt` — 매핑 통계
- `<sample_name>.depth.txt` — position별 coverage
- `<sample_name>.consensus.fasta` — consensus 서열
- `<sample_name>.vcf.gz` — variant calling 결과
- `<sample_name>_report.md` — 샘플별 요약 리포트

## 이미 분석을 끝낸 결과의 이름만 바꾸고 싶을 때

```bash
python scripts/rename_outputs.py \
    --samplesheet samplesheet.csv \
    --results-dir /path/to/existing_results \
    --dry-run            # 먼저 변경 내역만 확인

python scripts/rename_outputs.py \
    --samplesheet samplesheet.csv \
    --results-dir /path/to/existing_results \
    --rename-content     # 파일/폴더명 + 파일 내부 헤더(fasta/vcf 등)까지 변경
```

`--dry-run`을 빼면 실제로 파일/폴더 이름이 변경됩니다 (`--rename-content`를
주면 fasta 헤더, vcf 샘플명 등 텍스트 내용 안의 barcode ID도 함께 치환).
