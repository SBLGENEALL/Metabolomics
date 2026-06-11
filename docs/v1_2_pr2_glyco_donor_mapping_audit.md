# v1.2 PR2 Glycosylation Donor-Supply Mapping Audit

## Scope

This PR expands CHOmpact interpretation coverage for nucleotide-sugar donor
supply reactions already represented in iCHO3K. It is a mapping and
documentation update only.

iCHO3K remains the calculation engine. CHOmpact remains an interpretation
layer. This PR does not modify the iCHO3K model, exchange-rate calculation,
qMet constraints, pFBA/FVA calculation logic, or pathway scoring logic.

## Scientific boundary

The mapped biology is nucleotide-sugar donor-supply capacity and
glycosylation precursor supply only.

These mappings must not be interpreted as predictions of:

- glycan structure
- glycan quality
- glycoform
- fucosylation level
- sialylation level
- site occupancy
- Neu5Gc/Neu5Ac ratio
- glycosyltransferase kinetics
- Golgi/ER donor availability
- SLC35 transporter capacity
- Golgi/ER transport or localization

Flux through donor-supply reactions can support hypotheses about precursor
availability. It cannot determine glycan quality without additional glycomics,
enzyme, transport, compartment, and process data.

Donor-supply flux does not directly represent nucleotide-sugar concentration,
intracellular pool size, compartment-specific availability, or Golgi-accessible
donor abundance. UDP-GlcA donor-supply flux is not an mAb-specific
glycosylation marker.

## Added donor-supply mappings

Before PR2, `data/chompact_pathway_mapping.csv` contained zero rows assigned
to `Glycosylation donor supply`.

PR2 adds 16 conservative donor-supply mappings: 12 independent core mappings
and 4 model-alternative, cofactor-alternative, or compartment-duplicate
mappings. The non-independent mappings remain visible for model audit but must
not be counted as separate biological evidence.

| reaction_id | reaction_name | donor class | CHOmpact pathway | CHOmpact subpathway | inclusion rationale |
|---|---|---|---|---|---|
| UAGDP | UDP-N-acetylglucosamine diphosphorylase | UDP-GlcNAc | Glycosylation donor supply | Hexosamine / UDP-GlcNAc / UDP-GalNAc donor supply | Direct activation of GlcNAc-1-phosphate to UDP-GlcNAc. |
| UAG4E | UDP-N-acetylglucosamine 4-epimerase | UDP-GalNAc bridge | Glycosylation donor supply | Hexosamine / UDP-GlcNAc / UDP-GalNAc donor supply | Interconversion from UDP-GlcNAc to UDP-GalNAc donor pool. |
| UAGALDP | UDP-N-acetylgalactosamine diphosphorylase | UDP-GalNAc | Glycosylation donor supply | Hexosamine / UDP-GlcNAc / UDP-GalNAc donor supply | Direct activation of GalNAc-1-phosphate to UDP-GalNAc. |
| GALU | UTP-glucose-1-phosphate uridylyltransferase | UDP-Glucose | Glycosylation donor supply | UDP-Hexose / UDP-Gal donor supply | Direct activation of glucose-1-phosphate to UDP-glucose. |
| GALT | BiGGRxn35 | UDP-Galactose | Glycosylation donor supply | UDP-Hexose / UDP-Gal donor supply | Model alternative to UGLT; retained as non-independent evidence. |
| UGLT | UDPglucose--hexose-1-phosphate uridylyltransferase | UDP-Hexose bridge | Glycosylation donor supply | UDP-Hexose / UDP-Gal donor supply | UDP-glucose and UDP-galactose donor-pool interconversion bridge. |
| UDPG4E | UDPglucose 4-epimerase | UDP-Galactose bridge | Glycosylation donor supply | UDP-Hexose / UDP-Gal donor supply | Interconversion from UDP-glucose to UDP-galactose. |
| UDPGD | UDPglucose 6-dehydrogenase | UDP-GlcA | Glycosylation donor supply | UDP-GlcA donor supply | Direct formation of UDP-glucuronate from UDP-glucose. |
| MAN1PT2 | mannose-1-phosphate guanylyltransferase (GDP) | GDP-Mannose | Glycosylation donor supply | GDP-Mannose donor supply | GDP-dependent model alternative to r0208; retained as non-independent evidence. |
| r0208 | GTP:alpha-D-mannose-1-phosphate guanylyltransferase | GDP-Mannose | Glycosylation donor supply | GDP-Mannose donor supply | Direct GTP-dependent GDP-mannose activation route. |
| GMAND | GDP-D-mannose dehydratase | GDP-Fucose precursor | Glycosylation donor supply | GDP-Fucose donor supply | Handoff from GDP-mannose toward GDP-fucose donor synthesis. |
| GFUCS | GDP-L-fucose synthase | GDP-Fucose | Glycosylation donor supply | GDP-Fucose donor supply | Direct GDP-fucose donor-supply formation. |
| r0782 | GDP-L-fucose:NADP+ 4-oxidoreductase (3,5-epimerizing) | GDP-Fucose | Glycosylation donor supply | GDP-Fucose donor supply | NADH-dependent cofactor alternative to GFUCS; retained as non-independent evidence. |
| F1PGT | fucose-1-phosphate guanylyltransferase | GDP-Fucose salvage | Glycosylation donor supply | GDP-Fucose donor supply | Salvage activation of fucose-1-phosphate to GDP-fucose. |
| CMPSAS | CMP sialic acid synthase | CMP-Sialic acid | Glycosylation donor supply | CMP-Sialic acid donor supply | Cytosolic model compartment duplicate of CMPSASn; retained as non-independent evidence. |
| CMPSASn | CMP sialic acid synthase, nuclear | CMP-Sialic acid | Glycosylation donor supply | CMP-Sialic acid donor supply | Model compartment-specific CMP-sialic acid activation; no localization claim is made. |

## Intentionally excluded reactions

Ambiguous or out-of-scope reactions were intentionally excluded.

### Glycosyltransferases and glycan assembly

Reactions that consume activated donors on glycan acceptors were excluded.
These include representative `ABO*`, `B3GNT*`, `FUT*`, `ST*`, `GLCNACPT`,
`GLCNACT`, `BDMT`, `DOLPMT`, `HAS*`, keratan, chondroitin, and other glycan
assembly reactions.

Reason for exclusion: these reactions may reflect glycan assembly or donor
consumption, not donor-supply generation. Mapping them as donor supply would
overstate what v1.2 PR2 can interpret.

### Donor transport and localization

Transport/localization reactions were excluded, including representative
`UDPGLCAtg`, `UGALNACtg`, `UGLCNACtg`, `UDPACGALtl`, `r0842`, `r0845`, and
`r2519`.

Reason for exclusion: PR2 does not model Golgi transport/localization capacity.

### Degradation, detoxification, exchange, and demand

UDP-sugar degradation, glucuronidation detoxification, exchange, and demand
reactions were excluded, including representative `UDPG1P`, `UDPGNP`, `UGT*`,
bile acid/steroid/vitamin A glucuronosyltransferases, `EX_*`, and `DM_*`.

Reason for exclusion: these reactions do not represent donor-supply capacity
for CHO recombinant antibody glycosylation interpretation.

### Upstream non-activated precursor reactions

Upstream hexosamine and monosaccharide precursor reactions were deferred unless
they directly generate an activated nucleotide-sugar donor.

Reason for exclusion: PR2 is intentionally limited to activated donor-supply
mapping to avoid inflating interpretation coverage.

## Confidence and coverage safeguard

Mapping coverage is an ontology/QC measure, not positive biological evidence.
Adding donor-supply mappings must not automatically increase confidence.
Coverage can cap or limit confidence, but FVA robustness, feasible-range
separation, and reproducibility must dominate confidence.

Newly mapped full-FVA-only reactions remain model-emergent hypotheses until
they show robust FVA support. Missing mapping remains unavailable or not
evaluated; it must not be reported as zero activity.

## Validation status

Full-FVA validation was completed on branch
`v1.2-pr2-glyco-donor-mapping`.

CHOmpact validation command:

```bash
python run_pipeline.py --dataset practice_20aa --input_format tsv --steps 10,11,12,13,14 --fva_source full
```

Validation results:

- Steps 10-14: PASS.
- FVA provenance: `full_fva / all / broad_discovery`.
- Full FVA reactions: 8,368.
- Conditions: HighAvg, MotherAvg, and LowAvg in no-IgG-input and
  measured-demand modes.
- PPP audit: 39 model reactions, 13 mapped reactions, 33.33% coverage,
  `partially_covered`.
- Nucleotide-sugar donor audit: 16 model donor-supply reactions, 16 raw mapped
  reactions, 12 independent core mappings, and 4 model-alternative/duplicate
  mappings.
- Raw donor mapping coverage: 100%.
- Effective independent core coverage used as the confidence cap: 75%.
- Six donor-supply subpathways were scored and ranked.
- Fig14 and `Fig14_candidate_pathway_priority_confidence_data.csv`: generated.
- Fig14 companion provenance: `full_fva / all / broad_discovery`.

The four non-independent reactions are:

- `GALT`: `model_alternative` to `UGLT`
- `MAN1PT2`: `model_alternative` to `r0208`
- `r0782`: `cofactor_alternative` to `GFUCS`
- `CMPSAS`: `compartment_duplicate` of `CMPSASn`

They remain mapped for model transparency but are marked
`predictive_ranking_eligible = False` at reaction level and do not increase
effective core mapping coverage.
