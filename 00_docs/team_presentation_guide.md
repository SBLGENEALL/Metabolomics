# Team presentation guide

## Recommended title

`iCHO3K 기반 FBA/FVA를 이용한 CHO Producer Clone의 Metabolic Flux Phenotype 분석`

## Slide structure

1. Why: high vs low producer metabolic phenotype
2. Why spent-media alone is insufficient
3. FBA/FVA concept with toy pathway
4. Why direct `DM_igg_g` objective is not enough
5. Final pipeline workflow
6. Data QC and exchange-rate calculation
7. Exchange phenotype result
8. pFBA flux heatmap and delta
9. FVA overlap/separation result
10. Pathway-level summary
11. Focused Escher map
12. Limitations and next steps

## Key message

> We use FBA/FVA not as a direct titer predictor but as a mechanistic framework to compare feasible metabolic flux states and pathway flexibility between High and Low CHO producer clones.

## Important caveats

- pFBA is a representative solution.
- FVA is a feasible range, not a measured flux.
- Measured-demand mode is explanation, not prediction.
- Full/internal FVA is more complete than exchange-only FVA but can still contain under-constrained loops.

## Expected questions

### Why did direct IgG objective become flat?

Because measured exchange constraints alone may not capture clone-specific expression, translation, folding, secretion, or stress limitations.

### Is this actual intracellular flux?

No. It is a model-derived feasible flux state. Direct validation would require 13C-MFA or equivalent intracellular flux evidence.

### Why FVA?

Genome-scale models have many alternative feasible solutions. FVA reveals whether a reaction is tightly constrained or flexible under each clone condition.
