# CHANGELOG

All notable changes to the CHO Metabolomics project are documented in this file.

---

## Version lineage

```text
main
  = latest validated stable branch
  = currently Metabolomics_v1.0

Metabolomics_v1.0 tag
  = fixed v1.0 release snapshot
  = cleanup-step-order work merged into main

v1.1 branch
  = active development branch
  = CHOmpact interpretation layer and pathway-biomarker workflow

v1.1-pr* branches
  = future subtask branches for v1.1 work
  = merge back into v1.1 after local validation
```

Historical branch notes:

```text
feature/cleanup-step-order
  = completed v1.0 development branch
  = merged into main
  = deleted after v1.0 tag was created

feature/chompact-interpretation
  = original v1.1 branch name
  = renamed to v1.1
```

---

## [v1.1 - Unreleased]

Branch:

```text
v1.1
```

Base:

```text
main / Metabolomics_v1.0
```

Target:

```text
Metabolomics_v1.1
```

### Added

- Optional CHOmpact interpretation workflow after v1.0 steps
- Step 10: `10_map_to_chompact.py`
- Step 11: `11_score_chompact_pathways.py`
- Step 12: `12_demand_sensitivity.py`
- Step 13: `13_rank_pathway_biomarkers.py`
- Step 14: `14_make_chompact_figures.py`
- `data/chompact_pathway_mapping.csv`
- `IMPLEMENTATION_PLAN.md`
- `docs/biological_interpretation_guide.md`
- CHOmpact executive summary output concept

### Changed

- `run_pipeline.py` can now dispatch optional v1.1 steps 10-14
- CHOmpact is formally treated as an interpretation and visualization layer, not a calculation model
- Development branch naming policy changed from descriptive feature names to version branches such as `v1.1`, `v1.2`, and subtask branches such as `v1.1-pr1`

### Current validation status

- Step 10-14 smoke test: PASS when run without crashing
- Full biological/data validation: pending
- Current issue: v1.1 outputs can be empty if v1.0 tables are missing or mapping coverage is low
- Current issue: FVA console/log output can show repeated earlier `!! FVA failed` messages before final success, which is confusing and should be cleaned before v1.1 release

### v1.1 release criteria

- v1.0 steps 00-09 still pass on `practice_20aa`
- v1.1 steps 10-14 pass after valid v1.0 outputs exist
- `results/<dataset>/tables/chompact/` outputs are non-empty and interpretable
- `chompact_mapping_qc.csv` confirms acceptable mapping coverage
- v1.1 figures are generated correctly
- `CHOmpact_v1_1_executive_summary.md` is informative
- FVA failure/retry logging is clarified
- `README.md`, `MASTER.md`, and `CHANGELOG.md` are updated before merge to main

---

## [v1.0] - 2026-06

Branch/source lineage:

```text
feature/cleanup-step-order
  -> main
  -> tag: Metabolomics_v1.0
```

Status:

```text
stable release
```

### Added

- Clean sequential workflow structure (`00`-`09`)
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
