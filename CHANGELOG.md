# CHANGELOG

All notable changes to the CHO Metabolomics project are documented in this file.

---

## [v1.0] - 2026-06

### Added

- Clean sequential workflow structure (`00`–`09`)
- Linux/workstation TSV workflow
- Focused central/mAb pathway FVA
- Full/internal genome-scale FVA
- Pathway score generation
- Automated figure generation
- Automated Markdown report generation
- Offline-compatible report generation without `tabulate`
- MASTER.md project baseline document

### Changed

- Legacy workflow numbering removed from primary workflow
- README updated for v1.0 execution
- Figure naming standardized
- Report generation standardized
- Linux workstation validation route documented

### Validated

Dataset:

```text
practice_20aa
```

Workflow:

```text
00_data_qc
01_load_and_rate
02_map_metabolites
03_run_pfba
04_export_escher_flux
05_focused_fva
06_full_fva
07_pathway_scores
08_make_figures
09_generate_report
```

FVA mode:

```text
fva_scope = all
```

Validation result:

```text
PASS
```

### Deprecated

Legacy step identifiers:

```text
17
10
14
16
18
21
19
20
6
```

These are retained only for historical reference.

---

## Planned [v1.1]

### CHOmpact interpretation layer

- chompact_pathway_mapping.csv
- pathway-level interpretation
- CHOmpact-style visualization

### Robust biomarker ranking

- FVA robustness ranking
- pathway prioritization
- demand sensitivity analysis

### Reporting

- executive summary report
- publication-ready figures
- High vs Low pathway interpretation
