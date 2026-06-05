# CHO Productivity — Biological Interpretation Guide (Metabolomics v1.1)
**Calculation model:** iCHO3K (COBRApy 0.30.0) · **Interpretation layer:** CHOmpact-style categories · **Goal:** High vs Low productivity interpretation

## How to use this guide
iCHO3K produces all fluxes; CHOmpact supplies only the category labels and mechanistic framing. Every number you report should carry a **provenance tag**:
`[M]` measured · `[C]` constrained model input (measured rate imposed as a bound) · `[P]` predicted interior flux.
Interpretation strength is `[M] > [C] > [P]`. Interior `[P]` fluxes are hypotheses, never measurements.

> **Reaction IDs below are templates.** iCHO3K reconciles several databases, so confirm each reaction/exchange ID and compartment suffix (`_c` cytosol, `_m` mitochondria, `_e` extracellular) against your model build before wiring the biomarkers.

---

## 1. Final CHOmpact pathway category list

Five biological interpretation categories plus one methodological bin:

1. **Glycolysis / TCA / Oxidative phosphorylation** — central carbon & energy
2. **Amino acid & nucleotide metabolism** — anaplerosis & biosynthetic precursors
3. **Aspartate–malate shuttle** — cytosol↔mitochondria redox transfer
4. **Urea cycle / ammonia handling** — nitrogen disposal & detoxification
5. **Nucleotide-sugar-donor (NSD) metabolism** — glycosylation burden
6. **Constrained / cyclic fluxes** — identifiability scaffolding (**not a biomarker category**)

---

## 2. Pathway-by-pathway interpretation

### 1. Glycolysis / TCA / OxPhos
- **Biological meaning:** Glucose → pyruvate, then either lactate (fermentative overflow) or acetyl-CoA → TCA → electron transport chain (oxidative). This category sets the cell's fermentative-vs-oxidative balance and ATP source.
- **Expected relationship to productivity:** The primary axis. Oxidative/TCA-leaning, low-lactate metabolism — and especially the lactate-consumption "metabolic shift" — tracks the high-performer phenotype (longevity, higher yield). Sustained lactate overflow with glycolytic ATP tracks low performers.
- **Measured evidence to check `[M]`:** glucose and lactate specific rates; lactate/glucose yield; timing/depth of the lactate shift; oxygen uptake rate (OUR); optionally ¹³C-MFA TCA fluxes.
- **Model-derived evidence to check `[P]`:** PDH flux (`PDHm`), core TCA fluxes (`CSm`, `ICDHxm/ICDHym`, `AKGDm`), net LDH (`LDH_L`), and the oxidative-vs-glycolytic ATP fraction.
- **Caveats:** The lactate phenotype is strongly medium-dependent, not purely intrinsic to the clone. FBA tends to over-predict oxidative efficiency unless boundary-constrained. The pyruvate node is a classic alternate-optima site — use pFBA or sampling.

### 2. Amino acid & nucleotide metabolism
- **Biological meaning:** Amino-acid catabolism feeding the TCA cycle (anaplerosis), nitrogen flow, and nucleotide synthesis for growth.
- **Expected relationship to productivity:** High performers show stronger anaplerosis (glutamate → α-ketoglutarate via glutamate dehydrogenase; asparagine/aspartate routes), raising global TCA throughput. Nucleotide/biosynthetic demand tracks **growth** more than specific productivity — separate the two.
- **Measured evidence `[M]`:** amino-acid uptake/secretion rates (esp. Gln, Glu, Asn, Asp); ammonia; growth rate.
- **Model-derived evidence `[P]`:** glutamate dehydrogenase (`GLUDxm/GLUDym`), asparaginase/Asn handling, anaplerotic inflows into TCA, nucleotide-synthesis fluxes.
- **Caveats:** High constraint-driven risk — AA uptake is often imposed directly as bounds, so interior "differences" may echo inputs. Clonal drift confound. Disentangle growth-linked from product-linked demand.

### 3. Aspartate–malate shuttle
- **Biological meaning:** Moves cytosolic NADH reducing equivalents into the mitochondrion for oxidation while regenerating cytosolic NAD⁺ for glycolysis — no net carbon transfer.
- **Expected relationship to productivity:** Higher shuttle capacity lets glycolysis run without lactate overflow (NADH reoxidized oxidatively, not via LDH), supporting the oxidative high-performer phenotype.
- **Measured evidence `[M]`:** Largely indirect — infer from low lactate yield and, if available, NADH/NAD⁺ pools from metabolomics.
- **Model-derived evidence `[P]`:** malate–aspartate carrier fluxes (`AKGMALtm`, `ASPGLUm`), cytosolic vs mitochondrial malate dehydrogenase (`MDH`, `MDHm`), aspartate aminotransferase (`ASPTA`, `ASPTAm`).
- **Caveats:** Mostly inferred, rarely measured. GEMs carry parallel redox routes (e.g., glycerol-3-phosphate shuttle) that create alternate optima and make direction sensitive to constraints — sample, don't trust a single solution.

### 4. Urea cycle / ammonia handling
- **Biological meaning:** Nitrogen disposal and ammonia detoxification (CHOmpact specifically surfaced a mitochondrial ammonia-detox mechanism).
- **Expected relationship to productivity:** Secondary. High ammonia burden is growth- and glycosylation-toxic → associated with low performers and worse product quality; efficient detox is favorable. Not a direct titer driver.
- **Measured evidence `[M]`:** ammonia specific rate and accumulation; glutamine/glutamate dynamics.
- **Model-derived evidence `[P]`:** net NH₄ exchange (`EX_nh4_e`), urea-cycle reactions (`CBPS`, `OCBT`, `ARGSS`, `ARGSL`, `ARGN`), glutamine synthetase (`GLNS`).
- **Caveats:** Strongly process/medium-dependent. The CHO urea cycle is often incomplete, so the model may route detox through alternate reactions — read as a burden index, not a literal urea-cycle claim.

### 5. Nucleotide-sugar-donor (NSD) metabolism
- **Biological meaning:** Synthesis of UDP-/GDP-/CMP-sugars that donate to N- and O-glycosylation of host and recombinant glycoproteins.
- **Expected relationship to productivity:** Primarily a **product-quality** and secretory-burden signal (glycan occupancy/profile), not a titer driver. Scales with glycoprotein output.
- **Measured evidence `[M]`:** product glycan profile; glycoprotein qP; UDP-sugar pools if measured.
- **Model-derived evidence `[P]`:** NSD synthesis fluxes (`UDPG`, UDP-GlcNAc, GDP-Fuc, CMP-Neu5Ac routes); demand pulled toward glycosylation.
- **Caveats:** Treat as a quality axis. Titer correlation is weak/indirect. Demand is usually imposed via biomass/product composition → constraint-driven by construction.

### 6. Constrained / cyclic fluxes — **non-biomarker**
- **Biological meaning:** Numerical identifiability scaffolding (bounded futile/cyclic fluxes).
- **Expected relationship to productivity:** None. Do not use as a biomarker.
- **Evidence to check:** Monitor only that they stay within imposed bounds (sanity check).
- **Caveats:** Any apparent High–Low difference here is almost certainly numerical, not biological.

---

## 3. Biomarker definitions

All expressed as **ratios/yields** (scale-free, robust to demand-scale choice). State sign convention: in COBRApy exchanges, uptake is negative and secretion positive; formulas use magnitudes unless noted.

### 3.1 Lactate/glucose index — `Y_lac/glc`
```
Y_lac/glc = v(EX_lac__L_e, secretion) / |v(EX_glc__D_e, uptake)|
```
- **Direction:** ≈2 → homolactic Warburg overflow (worst); 0–1 → partial oxidative; **< 0 (net lactate consumption)** → metabolic shift (best).
- **Provenance:** `[M]` or `[C]` (both are exchange fluxes). This is the strongest single discriminator.

### 3.2 Oxidative entry fraction — `f_ox`
```
f_ox = v(PDHm) / [ v(PDHm) + max(0, v(EX_lac__L_e, secretion)) ]
```
- **Meaning:** fraction of pyruvate routed oxidatively into the TCA cycle vs spilled as lactate. Range 0–1; higher = more oxidative.
- **Provenance:** mixes `[P]` (PDH) and `[M]/[C]` (lactate) — flag accordingly.

### 3.3 Anaplerosis/cataplerosis ratio — `R_ana`
```
R_ana = [ v(PCm) + v(GLUDxm, glu→αKG) + Σ v(AA anaplerotic inflows) ]
        / [ v(ME, malate→pyruvate) + v(PEPCK, OAA→PEP) + Σ v(cataplerotic outflows) ]
```
- **Direction:** > 1 → net anaplerotic (high-producer-like TCA replenishment); < 1 → net cataplerotic.
- **Provenance:** `[P]`. Watch GLUD directionality (anaplerotic = glutamate → α-ketoglutarate).

### 3.4 Redox shuttle score — `S_redox`
```
S_redox = v(malate–aspartate carrier, cyto→mito; e.g. AKGMALtm) / v(GAPDH)
```
- **Meaning:** share of cytosolic (GAPDH-generated) NADH reoxidized via the shuttle/OxPhos rather than via LDH. Higher = oxidative redox handling. Roughly 0–1.
- **Provenance:** `[P]` (interior + inferred). Most fragile biomarker — always report alongside sampling spread.

### 3.5 Ammonia burden score — `S_NH4`
```
S_NH4 = v(EX_nh4_e, secretion) / |v(EX_gln__L_e, uptake)|
```
- **Direction:** lower = better nitrogen utilization/detox. Optionally append the fraction of nitrogen routed through urea-cycle detox reactions.
- **Provenance:** `[M]`/`[C]`. Quality/robustness marker, secondary to titer.

---

## 4. Interpretation rules

### 4.1 Measured vs predicted
- **Boundary fluxes:** trust `[M]`; use them as `[C]` constraints; predicted boundary values should *match* measured within error — this is validation, not discovery.
- **Interior fluxes:** `[P]` only — treat as hypotheses; never present an interior predicted flux as if measured.
- **Disagreement at the boundary** is a model/constraint problem to fix *before* interpreting any interior flux.
- Carry the `[M]/[C]/[P]` tag on every reported number.

### 4.2 Constraint-driven vs emergent (run before calling anything a biomarker)
1. **Source check.** If the discriminating flux *is* a measured exchange imposed as a bound, the separation exists by construction → label **constraint-driven**, report as descriptive only.
2. **Equalization test.** Set common boundary bounds for High and Low and re-solve. Separation persists → candidate **emergent**; separation collapses → **constraint-driven**.
3. **Demand-attribution test.** Does the interior flux change exceed what biomass/product/maintenance demand requires? Excess beyond demand → **emergent** (the canonical example: a pyruvate-carboxylase increase not explained by anabolic demand).
- **Pass criterion:** call a biomarker *emergent* only if it survives **both** the equalization and demand-attribution tests.

### 4.3 FVA range caution
- An FVA range is the **feasible span under the objective**, not biological variance and not a confidence interval — **never draw FVA ranges as error bars.**
- A wide range means an under-determined reaction (add constraints), not "high variability."
- Use High-vs-Low FVA **interval non-overlap** as a robustness flag — but a tightly bounded non-overlapping interval may be constraint-driven, so route it through Rule 4.2.
- Always report the objective fraction at which FVA was run (e.g., 99–100% of optimum) and keep it fixed across comparisons.

### 4.4 Pseudoreplication caution
- **The clone is the experimental unit, not the flux sample or the reaction.** Flux-sampling points from one clone, or many reactions within one clone, are *not* independent replicates — treating them as such inflates significance artificially.
- **Aggregate before testing.** Reduce each clone to one value per biomarker (e.g., the sampling mean/median), then compare across clones with clone-level *n*. Your statistical *n* is the number of clones (and biological/process replicates), never the number of sampled solutions.
- **Cross-validate at the clone level.** Any classifier separating High vs Low must be validated with clones held out as whole units (leave-clone-out), not with sampled points split across train/test from the same clone.
- **Correct for multiple testing across reactions.** When scanning the ~11,000 iCHO3K reactions, apply FDR control before declaring per-reaction significance.
- **Replicates vs samples — keep them distinct in reporting.** State explicitly whether an error bar reflects biological/process replicates `[M]` or flux-sampling spread `[P]`; the two must never be pooled.

---

## 5. Recommended figure captions for v1.1

**Figure 1. Central-carbon flux reallocation between high- and low-productivity CHO clones.**
iCHO3K parsimonious-FBA flux solutions mapped onto a CHOmpact-style central-carbon layout (rendered in Escher). Edge color encodes the direction of the high−low flux difference (oxidative TCA routing vs lactate overflow) and edge width its magnitude. Boundary (exchange) fluxes are constrained to measured rates `[C]`; interior fluxes are model-predicted `[P]`.

**Figure 2. Pathway-level productivity biomarkers ranked by high−low separation.**
Forest plot of the five interpretation biomarkers (lactate/glucose index, oxidative entry fraction, anaplerosis/cataplerosis ratio, redox shuttle score, ammonia burden score). Points show the high−low difference; filled markers denote non-overlapping FVA intervals (robust discriminators) and open markers overlapping intervals. Sampling spread shown as the distribution underlay.

**Figure 3. Robustness of biomarker separation across the demand-scale sensitivity sweep.**
Heatmap of biomarker high−low separation (color) across demand-scale settings (columns) for each biomarker (rows). Stable color across the row marks demand-scale-invariant — and therefore defensible — biomarkers; sign flips or fade-outs mark settings-dependent artifacts.

**Figure 4. Biomarker–productivity association across the clone panel.**
Scatter of each biomarker value against measured specific productivity (qP) across the full clone panel, with a fitted trend. Demonstrates that the biomarker tracks productivity across clones rather than merely separating a single high/low pair (guarding against the clonal-drift confound).

---

*Biological rationale and primary literature for the productivity links above are documented in the companion `pathway_review.md`. Confirm all reaction/exchange IDs against your iCHO3K build before implementation.*
