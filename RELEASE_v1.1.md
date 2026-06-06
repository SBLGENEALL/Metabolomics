# Metabolomics v1.1 Release Notes

Release tag:

```text
Metabolomics_v1.1
```

Release branch:

```text
v1.1
```

Base release:

```text
Metabolomics_v1.0
```

---

## Release summary

Metabolomics v1.1 adds the CHOmpact interpretation layer on top of the validated iCHO3K v1.0 workflow.

The v1.0 workflow remains the quantitative calculation engine:

- feed-corrected qMet generation
- exchange-rate-constrained iCHO3K FBA/pFBA
- focused and full FVA
- reaction-level and pathway-level outputs
- report and figure generation

CHOmpact is used for pathway labeling, evidence classification, visualization, and interpretation only. It does not recalculate fluxes and does not replace iCHO3K as the source of quantitative model output.

---

## Major additions in v1.1

### CHOmpact interpretation workflow

Optional steps 10-14 were added after the validated v1.0 steps 00-09:

```text
10_map_to_chompact.py
11_score_chompact_pathways.py
12_demand_sensitivity.py
13_rank_pathway_biomarkers.py
14_make_chompact_figures.py
```

These steps read existing iCHO3K/v1.0 outputs and convert them into CHOmpact-labeled interpretation tables, candidate rankings, executive summaries, and figures.

### Evidence classification

v1.1 separates evidence provenance into interpretation-safe classes:

```text
measured evidence
model-emergent evidence
constraint-driven evidence
product-demand-driven evidence
```

This distinction is required because measured exchange rates, imposed model constraints, product-demand effects, and emergent internal pFBA/FVA hypotheses should not be interpreted as equivalent evidence.

### FVA provenance tracking

v1.1 adds explicit FVA source selection:

```text
--fva_source auto
--fva_source full
--fva_source focused
```

`--fva_source auto` uses full all/internal FVA when available and otherwise falls back to focused FVA. The selected source is recorded in CHOmpact ranking outputs using:

```text
fva_source
fva_scope
discovery_role
```

Full FVA with `--fva_scope all` is the preferred route for broad-discovery interpretation. Focused FVA should be interpreted as central/mAb confirmation or smoke-test support.

### Conservative confidence scoring

v1.1 keeps these metrics separate:

```text
priority_score
confidence_score
evidence_coverage_score
mapping_coverage_score
robustness_score
```

Missing evidence remains unavailable/NA rather than being converted to zero. Missing evidence lowers evidence coverage and may cap confidence.

Current conservative confidence caps include:

```text
focused FVA only          -> confidence_score <= 70
mapping coverage < 20%    -> confidence_score <= 50
no FVA evidence           -> confidence_score <= 60
```

Product-demand-driven IgG reactions are excluded from predictive biomarker ranking by design because they reflect the imposed product-demand constraint rather than emergent metabolism.

### Domain coverage audit

v1.1 adds:

```text
results/<dataset>/tables/chompact/chompact_domain_coverage_audit.csv
```

This file distinguishes absent biological activity from limited interpretation coverage.

Coverage audit status values include:

```text
adequately_covered
partially_covered
mapping_missing
full_fva_only
not_evaluated
unavailable
```

`mapping_missing`, `full_fva_only`, `not_evaluated`, and `unavailable` must not be interpreted as zero pathway activity.

### Fig14 candidate ranking visualization

v1.1 adds a release-polished Fig14:

```text
results/<dataset>/figures/Fig14_candidate_pathway_priority_confidence.png
```

Fig14 shows:

- priority score
- confidence score
- evidence type
- pathway family
- evidence coverage

A companion table is also generated:

```text
results/<dataset>/tables/chompact/Fig14_candidate_pathway_priority_confidence_data.csv
```

This table records the plotted candidate values and FVA provenance for auditability.

---

## Validation summary

The release candidate was validated on `practice_20aa`.

Confirmed outputs:

```text
Fig1-Fig14
SuppFig1
chompact_domain_coverage_audit.csv
Fig14_candidate_pathway_priority_confidence_data.csv
chompact_ranked_pathway_biomarkers.csv
CHOmpact_v1_1_executive_summary.md
```

Windows release-polish validation confirmed:

```text
fva_source = full_fva
fva_scope = all
discovery_role = broad_discovery
```

Fig14 generation, companion CSV generation, coverage audit generation, and FVA provenance propagation were confirmed.

---

## Known limitations

Absence from CHOmpact ranking must not be interpreted as zero biological activity.

Current domain interpretation status:

```text
PPP
  = partially covered

nucleotide metabolism
  = present in iCHO3K/full-FVA layer
  = insufficiently mapped in current CHOmpact interpretation layer
  = flagged as full_fva_only or not fully evaluated

lipid metabolism
  = present in iCHO3K/full-FVA layer
  = insufficiently mapped in current CHOmpact interpretation layer
  = flagged as full_fva_only or not fully evaluated

glycosylation
  = present in iCHO3K/full-FVA layer
  = insufficiently mapped in current CHOmpact interpretation layer
  = flagged as full_fva_only or not fully evaluated

nucleotide-sugar donor metabolism
  = present in iCHO3K/full-FVA layer
  = insufficiently mapped in current CHOmpact interpretation layer
  = flagged as full_fva_only or not fully evaluated
```

These domains are deferred to v1.2 ontology and mapping expansion.

---

## Deferred to v1.2+

The following features are intentionally not part of the v1.1 release:

- expanded CHOmpact/iCHO3K ontology
- PPP, nucleotide, lipid, glycosylation, and nucleotide-sugar donor module expansion
- flux sampling
- OCR/ECAR integration
- 13C-MFA integration
- iCHO2048s secretory model comparison
- transcriptomics constraints
- publication-grade CHOmpact network map

---

## Release recommendation

Metabolomics v1.1 is approved for tagging once this release note, `CHANGELOG.md`, `MASTER.md`, and `README.md` are synchronized on the `v1.1` branch.
