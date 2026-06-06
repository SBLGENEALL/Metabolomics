# Metabolomics MASTER

This document is the single source of truth for the CHO Metabolomics project. Older chat history and legacy step numbers are superseded by this file.

---

## 0. Version policy

Current stable release:

```text
Metabolomics_v1.0
```

Stable branch:

```text
main
```

Current development branch:

```text
v1.1
```

Current development target:

```text
Metabolomics_v1.1
```

Version meaning:

```text
v1.0
  = validated iCHO3K production workflow
  = cleanup-step-order work merged into main
  = fixed by tag Metabolomics_v1.0

v1.1
  = CHOmpact interpretation integration release
  = optional steps 10-14 added after the validated 00-09 workflow
  = current branch for bugfixes, figure refinement, and interpretation polishing

v1.1.x / v1.1-pr*
  = v1.1-level fixes and polish only
  = examples: logging cleanup, figure title/style fixes, score normalization, report wording
  = merge back into v1.1 after local validation

v1.2+
  = future larger feature releases
  = examples: expanded pathway ontology, glycosylation module, nucleotide/lipid modules, ER folding/secretion module, multi-omics integration, CHOmpact network map redesign
```

Branch interpretation:

```text
main
  = latest validated stable branch
  = currently equivalent to Metabolomics_v1.0

v1.1
  = active development/release-candidate branch for CHOmpact integration

backup_before_archive_reorg
  = historical backup branch
```

Historical branches:

```text
feature/cleanup-step-order
  = completed v1.0 development branch
  = merged into main
  = deleted after v1.0 tag was created

feature/chompact-interpretation
  = original v1.1 development branch name
  = renamed to v1.1
```

---

## 1. Project purpose

The CHO Metabolomics pipeline converts fed-batch culture and spent-media metabolomics data into exchange-rate-constrained iCHO3K FBA/FVA outputs.

The pipeline is designed for:

- clone/group comparison
- productivity-group interpretation
- measured exchange-rate constraint generation
- pFBA flux-state comparison
- focused central/mAb pathway FVA
- full/internal genome-scale FVA
- pathway-level interpretation
- report and figure generation

This pipeline is **not** a stand-alone IgG/mAb titer predictor.

Actual clone grouping is intentionally deferred until real clone data are available. The code should support two-group, three-group, or multi-group comparisons rather than hard-coding High/Mother/Low assumptions.

---

## 2. Core modeling strategy

### iCHO3K

Role:

- calculation engine
- COBRApy FBA/pFBA/FVA
- exchange-rate-constrained clone comparison
- full/internal genome-scale FVA

Interpretation:

- iCHO3K reaction-level outputs are the quantitative source of truth.

### CHOmpact

Role in v1.1:

- pathway vocabulary
- biological grouping
- figure/report interpretation layer
- CHO-specific mechanistic framing

Important principle:

- CHOmpact must not replace iCHO3K as the calculation model.
- CHOmpact must not be used to recompute fluxes.
- CHOmpact lumped fluxes must not be reconstructed by summing iCHO3K reactions.
- CHOmpact should summarize, label, and visualize iCHO3K outputs only.

---

## 3. Environment

Validated v1.0 environment:

- Linux workstation
- Python 3.13 environment: `metabolomics_env_py313`
- COBRApy 0.30.0
- GLPK solver
- offline-compatible TSV input

The v1.0 report generator does not require the optional `tabulate` package.

v1.1 steps 10-14 are interpretation-layer scripts that read v1.0 CSV outputs and should not require COBRApy unless the full 00-09 calculation pipeline is rerun.

---

## 4. v1.0 validated workflow

The v1.0 pipeline uses clean sequential steps `00` through `09`.

```text
00_data_qc.py
01_load_and_rate.py
02_map_metabolites.py
03_run_pfba.py
04_export_escher_flux.py
05_focused_fva.py
06_full_fva.py
07_pathway_scores.py
08_make_figures.py
09_generate_report.py
```

Legacy step numbers such as 17, 10, 14, 21, 19, 20, and 6 are deprecated.

Validated full-model FVA command:

```bash
python run_pipeline.py \
  --dataset practice_20aa \
  --input_format tsv \
  --rate_days 7,10 \
  --analysis_mode both \
  --constraint_policy production_relaxed \
  --feed_volume_mode interval \
  --demand_scale auto \
  --steps 00,01,02,03,04,05,06,07,08,09 \
  --fva_scope all \
  --fva_targets group_avg \
  --fva_processes 16
```

Fast route without full FVA:

```bash
python run_pipeline.py \
  --dataset practice_20aa \
  --input_format tsv \
  --steps 00,01,02,03,04,05,07,08,09
```

---

## 5. v1.1 workflow

The v1.1 branch adds optional interpretation steps after v1.0.

```text
10_map_to_chompact.py
11_score_chompact_pathways.py
12_demand_sensitivity.py
13_rank_pathway_biomarkers.py
14_make_chompact_figures.py
```

Run v1.1 extension after v1.0 outputs already exist:

```bash
python run_pipeline.py \
  --dataset practice_20aa \
  --steps 10,11,12,13,14 \
  --fva_source auto
```

CHOmpact FVA source policy:

```text
--fva_source auto
  = prefer non-empty full FVA from all/internal scope
  = fall back to focused central/mAb FVA

--fva_source full
  = require full all/internal FVA
  = fail clearly when unavailable

--fva_source focused
  = use central/mAb focused FVA only
  = label results as focused confirmation with limited discovery coverage
```

All CHOmpact ranking outputs carry:

```text
fva_source
fva_scope
discovery_role
```

Run v1.0 plus v1.1 together:

```bash
python run_pipeline.py \
  --dataset practice_20aa \
  --input_format tsv \
  --steps 00,01,02,03,04,05,06,07,08,09,10,11,12,13,14 \
  --analysis_mode both \
  --fva_scope all
```

Note: `--fva_scope all` in step 06 is computationally expensive. For quick v1.1 interpretation-layer smoke tests, run steps 10-14 only after existing v1.0 tables are available.

---

## 6. Input data structure

Practice TSV input:

```text
data/raw/practice_20aa_tsv/raw_timeseries.tsv
data/raw/practice_20aa_tsv/metabolite_map.tsv
data/raw/practice_20aa_tsv/feed_composition.tsv
```

Own experiment TSV input:

```text
data/raw/own_experiment/raw_timeseries.tsv
data/raw/own_experiment/metabolite_map.tsv
data/raw/own_experiment/feed_composition.tsv
```

---

## 7. v1.0 outputs

Main figures:

```text
Fig1_IgG_timecourse.png
Fig2_exchange_rate_heatmap.png
Fig3_lactate_glucose_phenotype.png
Fig4_high_low_exchange_rates.png
Fig5_central_mab_flux_heatmap.png
Fig6_central_mab_flux_zscore.png
Fig7_central_mab_flux_delta.png
Fig8_pathway_scores_fva_overlap.png
Fig9_full_fva_high_low_separation.png
Fig10_full_fva_range_delta.png
Fig11_summary_panel.png
```

Supplementary figure:

```text
SuppFig1_data_qc_overview.png
```

Main report:

```text
results/<dataset>/REPORT_SUMMARY.md
```

---

## 8. v1.1 outputs

CHOmpact interpretation tables:

```text
results/<dataset>/tables/chompact/
```

Expected v1.1 outputs include:

```text
chompact_mapping_qc.csv
chompact_mapped_flux.csv
chompact_mapped_fva.csv
chompact_mapped_fva_overlap.csv
chompact_mapped_measured_rates.csv
chompact_pathway_fba_activity_scores.csv
chompact_pathway_fva_robustness_scores.csv
chompact_pathway_measured_qmet_scores.csv
chompact_pathway_high_low_separation.csv
chompact_demand_scale_sensitivity.csv
chompact_robust_rank_across_demand_scales.csv
chompact_ranked_pathway_biomarkers.csv
chompact_top_pathway_biomarkers_for_report.csv
chompact_domain_coverage_audit.csv
Fig14_candidate_pathway_priority_confidence_data.csv
CHOmpact_v1_1_executive_summary.md
```

v1.1 figures:

```text
Fig12_measured_screening_markers.png
Fig13_model_emergent_fva_robustness.png
Fig14_candidate_pathway_priority_confidence.png
```

---

## 9. Glossary

This glossary is intended to prevent future ambiguity between measured data, model constraints, model outputs, and engineering interpretation.

| Term | Meaning in this project | Interpretation caution |
|---|---|---|
| qMet | Feed-corrected metabolite-specific exchange rate calculated from spent-media time course data. | Measured/derived extracellular rate, not an intracellular flux. |
| pFBA | Parsimonious flux balance analysis. A model solution that satisfies constraints while minimizing total flux usage. | One optimal model solution, not a direct measurement. |
| FVA | Flux variability analysis. For each reaction, reports feasible minimum and maximum flux under model constraints. | Feasible range, not statistical SD/error. |
| FVA overlap | A condition comparison where two groups' feasible FVA intervals overlap. | Overlap suggests the model cannot strongly separate the groups for that reaction/pathway. |
| FVA non-overlap | A condition comparison where feasible FVA intervals do not overlap. | Non-overlap is stronger model evidence for separation, but still model-based. |
| Robustness score | A pathway or candidate score summarizing whether the signal remains stable across FVA ranges and/or demand-scale checks. | High robustness increases confidence but does not prove causality. |
| Measured evidence | Evidence directly measured or directly derived from measurements, such as IgG, VCD, viability, qMet, uptake, and secretion. | Highest evidence class for experimental observation. |
| Model-emergent evidence | Internal pFBA/FVA signal emerging from model constraints rather than directly measured outside the cell. | Hypothesis-generating unless independently validated. |
| Constraint-driven | A result that follows directly from imposed model constraints, such as measured exchange bounds or product demand. | Should be reported descriptively, not as an independently discovered mechanism. |
| Product-demand-driven | A model effect primarily caused by IgG/product demand constraints. | Important for production burden interpretation; do not overinterpret as measured metabolic regulation. |
| Screening marker | A measured or model-derived feature useful for ranking or triaging clones. | Can be useful even when mechanism is incomplete. |
| Engineering target | A pathway, reaction, or process proposed for perturbation or host/vector/process engineering. | Requires stronger biological plausibility and ideally validation beyond model output. |
| Priority score | Ranking score used to order candidate pathways/reactions by effect size, robustness, and relevance. | Ranking aid, not a probability of success. |
| Confidence score | Score representing how reliable the candidate is given evidence type, coverage, and robustness. | Should be interpreted together with evidence class and coverage. |

---

## 10. Figure dictionary

All figure filenames, internal figure numbers, and `REPORT_SUMMARY.md` references should remain synchronized.

| Figure | File | What it shows | How to interpret |
|---|---|---|---|
| Fig1 | `Fig1_IgG_timecourse.png` | IgG titer or production time course across culture days/groups. | Confirms productivity pattern and group separation. This is measured process output, not model prediction. |
| Fig2 | `Fig2_exchange_rate_heatmap.png` | Exchange-rate heatmap for measured metabolites. | Shows extracellular uptake/secretion patterns that become model exchange constraints. |
| Fig3 | `Fig3_lactate_glucose_phenotype.png` | Lactate/glucose relationship or ratio. | Interprets overflow metabolism and glucose-lactate phenotype. |
| Fig4 | `Fig4_high_low_exchange_rates.png` | High-vs-Low or group-level exchange-rate differences. | Highlights measured qMet features that may explain or constrain downstream model results. |
| Fig5 | `Fig5_central_mab_flux_heatmap.png` | Central/mAb-related pFBA flux state. | Model-predicted intracellular flux state under measured constraints. |
| Fig6 | `Fig6_central_mab_flux_zscore.png` | Z-scored central/mAb flux pattern. | Useful for pattern comparison; z-score magnitude is relative within the plotted matrix. |
| Fig7 | `Fig7_central_mab_flux_delta.png` | Group delta in central/mAb fluxes. | Highlights predicted reaction-level shifts between groups. |
| Fig8 | `Fig8_pathway_scores_fva_overlap.png` | Pathway score with FVA overlap information. | Separates robust non-overlap candidates from ambiguous overlapping FVA candidates. |
| Fig9 | `Fig9_full_fva_high_low_separation.png` | Full/internal genome-scale FVA group separation. | Broad model scan for reactions with potential feasible-range separation. |
| Fig10 | `Fig10_full_fva_range_delta.png` | Full FVA range-width or range-delta summary. | Identifies reactions/pathways with changed flexibility or constraint tightness. |
| Fig11 | `Fig11_summary_panel.png` | Summary panel of v1.0 outputs. | High-level dashboard for report-level interpretation. |
| Fig12 | `Fig12_measured_screening_markers.png` | Measured screening and feed-media markers. | Separates observed qMet evidence from model-emergent hypotheses. |
| Fig13 | `Fig13_model_emergent_fva_robustness.png` | Model-emergent pathway robustness from mapped FVA outputs. | Prioritizes pathways with stable model support across feasible ranges. |
| Fig14 | `Fig14_candidate_pathway_priority_confidence.png` | Candidate pathway priority-confidence view. | Displays priority and confidence separately with evidence provenance. |
| SuppFig1 | `SuppFig1_data_qc_overview.png` | Data QC overview. | Use before interpreting model output; bad input quality invalidates downstream conclusions. |

---

## 11. Output dictionary

Primary output files should be documented by biological meaning, modeling role, and evidence class.

| Output file | Meaning | Primary evidence class | Typical use |
|---|---|---|---|
| `exchange_rates.csv` | Feed-corrected uptake/secretion rates calculated from time-course extracellular data. | Measured/derived qMet evidence. | Input constraints, extracellular phenotype comparison, screening marker discovery. |
| `central_mab_flux_state.csv` | pFBA flux solution for selected central metabolism and mAb-related reactions. | Model-emergent pFBA evidence. | Compare predicted intracellular flux states between groups. |
| `full_fva_all_combined_report.csv` | Genome-scale/internal FVA report combining feasible min/max ranges and group comparison metrics. | Model-emergent FVA evidence. | Broad reaction-level candidate discovery and coverage audit. |
| `reaction_level_candidates.csv` | Ranked reaction-level candidate list derived from pFBA/FVA/pathway scoring. | Mixed; usually model-emergent unless tied to measured qMet. | Prioritize screening markers or engineering targets. |
| `pathway_scores.csv` | v1.0 pathway-level scores summarized from reaction-level outputs. | Model-derived pathway evidence. | Pathway interpretation and report summary. |
| `chompact_mapping_qc.csv` | QC table describing how iCHO3K reactions/metabolites map into CHOmpact interpretation groups. | Documentation/QC evidence. | Check pathway coverage before trusting CHOmpact summaries. |
| `chompact_mapped_flux.csv` | iCHO3K pFBA flux values relabeled into CHOmpact pathway vocabulary. | Model-emergent pFBA evidence. | CHOmpact-based activity visualization. |
| `chompact_mapped_fva.csv` | iCHO3K FVA ranges relabeled into CHOmpact pathway vocabulary. | Model-emergent FVA evidence. | CHOmpact robustness scoring. |
| `chompact_mapped_fva_overlap.csv` | FVA overlap/non-overlap summary after CHOmpact mapping. | Model-emergent FVA evidence. | Identify robust vs ambiguous pathway separation. |
| `chompact_mapped_measured_rates.csv` | Measured qMet rates mapped to CHOmpact pathway labels where possible. | Measured/derived evidence. | Distinguish measured extracellular support from internal model hypotheses. |
| `chompact_pathway_fba_activity_scores.csv` | CHOmpact pathway activity scores from mapped pFBA outputs. | Model-emergent pFBA evidence. | Pathway activity ranking. |
| `chompact_pathway_fva_robustness_scores.csv` | CHOmpact pathway robustness scores from mapped FVA intervals. | Model-emergent FVA evidence. | Confidence support for pathway ranking. |
| `chompact_pathway_measured_qmet_scores.csv` | CHOmpact pathway scores based on measured exchange-rate evidence. | Measured/derived qMet evidence. | Measured evidence layer for candidate ranking. |
| `chompact_pathway_high_low_separation.csv` | CHOmpact pathway-level High-vs-Low or group separation summary. | Mixed evidence. | Report-level group comparison. |
| `chompact_demand_scale_sensitivity.csv` | Sensitivity of pathway/candidate scores to product-demand scaling. | Product-demand-driven model evidence. | Detect whether a candidate is mainly demand-driven. |
| `chompact_robust_rank_across_demand_scales.csv` | Candidate/pathway ranking stability across demand scales. | Model robustness evidence. | Identify candidates that persist across product-demand assumptions. |
| `chompact_ranked_pathway_biomarkers.csv` | Main ranked CHOmpact candidate table with priority/confidence/evidence annotations. | Mixed evidence, explicitly classified. | Final screening-marker and engineering-target triage. |
| `chompact_top_pathway_biomarkers_for_report.csv` | Short report-facing subset of ranked CHOmpact candidates. | Mixed evidence, explicitly classified. | Executive summary and Fig14 input. |
| `chompact_domain_coverage_audit.csv` | Domain-level audit from iCHO3K model coverage through focused/full FVA, mapping, scoring, and ranking. | Coverage/QC evidence. | Distinguish absent activity from unavailable or insufficiently mapped biology. |
| `Fig14_candidate_pathway_priority_confidence_data.csv` | Exact plotted data and visual encodings for Fig14. | Mixed evidence, explicitly classified. | Audit priority, confidence, coverage, evidence type, pathway family, and FVA source. |
| `CHOmpact_v1_1_executive_summary.md` | Human-readable v1.1 CHOmpact interpretation summary. | Narrative summary. | Review, communication, and release validation. |

---

## 12. Interpretation rules

### Measured values

Measured or directly derived from experimental data:

- IgG titer
- VCD
- viability
- feed-corrected qMet / uptake / secretion rates
- glucose, lactate, ammonia, amino acid exchange rates

### Model-constrained values

Measured rates imposed as model constraints:

- exchange bounds derived from measured qMet
- product demand constraints derived from measured IgG/qP

### Model-predicted values

Model-derived:

- pFBA flux values
- FVA minimum/maximum ranges
- internal reaction fluxes
- model objective values
- pathway scores derived from iCHO3K outputs

Reports must clearly distinguish measured evidence from constrained model inputs and predicted model interpretation.

### FVA caution

FVA ranges are model-based feasible ranges, not direct intracellular flux measurements and not statistical error bars. Very broad FVA intervals may reflect under-constrained reactions, alternate optima, or thermodynamic loops.

### Constraint-driven vs emergent interpretation

If a High-vs-Low difference is directly imposed as a boundary constraint, it is constraint-driven and should be reported descriptively. Interior differences that persist after boundary equalization or demand-attribution checks may be treated as emergent model hypotheses.

### Confidence and coverage rules

The five decision metrics remain separate:

```text
priority_score
confidence_score
evidence_coverage_score
mapping_coverage_score
robustness_score
```

`evidence_coverage_score` is calculated as available evidence divided by expected
evidence. Missing evidence remains `NA`; it is not converted to zero. Missing
evidence lowers coverage and may limit confidence.

Mapping coverage is an ontology/QC measure, not positive biological evidence.
Improved mapping coverage must not automatically raise confidence. Coverage can
cap or limit confidence, while FVA robustness, feasible-range separation, and
reproducibility should dominate confidence. Newly mapped full-FVA-only,
under-constrained, or loop-prone reactions remain model-emergent hypotheses and
must not inflate confidence without robust FVA support.

Current conservative confidence caps:

```text
focused FVA only          -> confidence_score <= 70
mapping coverage < 20%    -> confidence_score <= 50
no FVA evidence           -> confidence_score <= 60
```

Model-only candidates remain classified as model-emergent hypotheses regardless
of score. These caps prevent focused or poorly mapped evidence from appearing as
broad, high-confidence pathway discovery.

### Coverage audit interpretation

The domain coverage audit uses only these status values:

```text
adequately_covered
partially_covered
mapping_missing
full_fva_only
not_evaluated
unavailable
```

`mapping_missing`, `full_fva_only`, `not_evaluated`, and `unavailable` must never
be interpreted as zero pathway activity.

---

## 13. Release-polish checklist for v1.1-pr1

### Documentation

- [ ] Keep `MASTER.md`, `README.md`, `CHANGELOG.md`, and report wording synchronized.
- [ ] Keep figure filenames, internal figure numbers, and `REPORT_SUMMARY.md` references synchronized.
- [ ] Ensure all top-ranked candidates carry an evidence class: measured, model-emergent, constraint-driven, or product-demand-driven.
- [ ] Ensure screening markers and engineering targets are not mixed without labeling.

### Fig14 visualization polish

Required Fig14 encoding:

```text
shape:
  circle = model-emergent
  triangle = measured

color:
  TCA
  PPP
  OXPHOS
  Glutamine
  Exchange
  Other / unmapped

hover / annotation fields:
  pathway
  priority
  confidence
  coverage
  evidence type
```

Fig14 should not show only priority and confidence. Evidence type and pathway coverage must be visible either directly in the figure or in hover/label/exported table fields.

Implemented Fig14 encoding:

```text
shape:
  circle = model-emergent
  triangle = measured

color:
  TCA
  PPP
  OXPHOS
  Glutamine
  Exchange
  Other

size:
  evidence_coverage_score
```

### Coverage audit

Priority biological coverage audit targets:

```text
PPP
nucleotide metabolism
lipid metabolism
glycosylation / nucleotide-sugar donors
```

For each target, audit whether the signal is absent because biology is absent, model coverage is weak, mapping is missing, or the pathway is buried in full FVA outputs.

### FVA source audit

Confirm and document whether CHOmpact interpretation uses:

```text
focused FVA
full FVA
both, with explicit labels
```

v1.1 interpretation should prefer full FVA for broad discovery when available, but may use focused FVA for fast smoke tests only if the report clearly labels this limitation.

---

## 14. v1.1 release criteria

Metabolomics v1.1 should not be merged to `main` or tagged until:

1. v1.0 steps 00-09 still pass on `practice_20aa`.
2. v1.1 steps 10-14 pass after valid v1.0 outputs are present.
3. CHOmpact outputs are non-empty and interpretable.
4. `chompact_mapping_qc.csv` confirms acceptable mapping coverage.
5. Figures 12-14 and the executive summary are generated correctly.
6. FVA failure/retry logging is clarified enough not to confuse prior failed attempts with final success.
7. `README.md`, `MASTER.md`, and `CHANGELOG.md` are updated.
8. `MASTER.md` glossary, figure dictionary, and output dictionary are complete enough for future users to interpret v1.0/v1.1 outputs without old chat history.
9. Fig14 shows priority, confidence, evidence type, pathway class, and coverage.
10. PPP, nucleotide, lipid, and glycosylation coverage limitations are documented or fixed.

The intended v1.1 tag will be:

```text
Metabolomics_v1.1
```

---

## 15. Future v1.2 candidates

The following should be deferred to v1.2 or later unless the v1.1 review explicitly decides otherwise:

- flux sampling
- OCR/ECAR integration
- 13C-MFA integration
- iCHO2048s secretory model
- transcriptomics constraints
- expanded CHOmpact/iCHO3K pathway ontology
- glycosylation and nucleotide-sugar donor module expansion
- nucleotide and lipid metabolism interpretation modules
- ER folding/secretion burden module
- publication-grade CHOmpact network map
- multi-omics integration
