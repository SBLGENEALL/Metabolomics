# FBA and FVA concepts

## Toy pathway

```text
Glucose --v1--> Pyruvate
Pyruvate --v2--> Lactate
Pyruvate --v3--> TCA
TCA --v4--> ATP
ATP + AA --v5--> IgG
Pyruvate --v6--> Biomass
```

## Mass balance

At steady state, each internal metabolite has no net accumulation.

```text
S · v = 0
```

Example balances:

```text
Pyruvate: v1 - v2 - v3 - v6 = 0
TCA:      v3 - v4 = 0
ATP:      v4 - v5 = 0
```

## Constraint

Measured glucose uptake can be represented as:

```text
0 ≤ v1 ≤ 10
```

## Objective

Canonical product objective:

```text
maximize v5
```

This asks: “How much IgG can the model theoretically produce under the given constraints?”

## FBA

FBA returns one optimal or representative flux state. In genome-scale models, many alternative optimal solutions can exist.

## FVA

FVA calculates the minimum and maximum possible flux for each reaction while satisfying the constraints and a chosen optimum fraction.

Example:

```text
High producer TCA range: 7–9
Low producer TCA range: 0–2
```

This indicates strong separation.

## Overlap fraction

`overlap_fraction` measures how much the High and Low FVA ranges overlap.

```text
0 = no overlap; distinct feasible flux ranges
1 = almost complete overlap; similar feasible flux ranges
```

Interpretation:
- Low overlap: stronger evidence for different feasible metabolic states.
- High overlap: that reaction is less useful for separating High vs Low.

## Measured-demand mode

Measured-demand mode fixes the IgG demand reaction to observed qIgG.

It asks:

```text
If this clone is actually producing this much IgG, what metabolic state is required to sustain it?
```

This is an explanation mode, not prediction.
