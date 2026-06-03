# Workflow

## Analysis philosophy

This pipeline is designed for mechanistic interpretation of clone-specific metabolic phenotypes. It should not be presented as a direct titer prediction model.

## Workflow

```text
Raw culture / spent-media data
    ↓
Data QC
    ↓
Feed-corrected exchange rate calculation
    ↓
Clone-specific constraints on iCHO3K production model
    ↓
pFBA / focused FVA / internal or full FVA
    ↓
Flux heatmap, FVA overlap, pathway score, Escher map
    ↓
Mechanistic interpretation of High vs Low producer metabolic phenotype
```

## Recommended interpretation order

1. QC warnings and report summary
2. Exchange phenotype: glucose, lactate, NH4, amino acids
3. pFBA flux differences
4. FVA overlap/separation
5. Pathway-level scores
6. Focused Escher maps

## What not to claim

- Do not claim that FBA directly measured intracellular flux.
- Do not claim that `DM_igg_g` objective value predicts titer.
- Do not claim that measured-demand mode is prediction.

## What to claim

- Measured exchange phenotype was used to constrain feasible metabolic states.
- pFBA provides representative flux solutions.
- FVA provides feasible flux ranges.
- Low overlap between High and Low FVA ranges indicates distinct feasible metabolic spaces.
