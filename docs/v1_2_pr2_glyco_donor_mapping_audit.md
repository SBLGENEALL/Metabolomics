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
- fucosylation level
- sialylation level
- site occupancy
- glycosyltransferase kinetics
- Golgi transport or localization

Flux through donor-supply reactions can support hypotheses about precursor
availability. It cannot determine glycan quality without additional glycomics,
enzyme, transport, compartment, and process data.

## Added donor-supply mappings

Before PR2, `data/chompact_pathway_mapping.csv` contained zero rows assigned
to `Glycosylation donor supply`.

PR2 adds 16 conservative, high-confidence donor-supply mappings.

| reaction_id | reaction_name | donor class | CHOmpact pathway | CHOmpact subpathway | inclusion rationale |
|---|---|---|---|---|---|
| UAGDP | UDP-N-acetylglucosamine diphosphorylase | UDP-GlcNAc | Glycosylation donor supply | Hexosamine / UDP-GlcNAc donor supply | Direct activation of GlcNAc-1-phosphate to UDP-GlcNAc. |
| UAG4E | UDP-N-acetylglucosamine 4-epimerase | UDP-GalNAc bridge | Glycosylation donor supply | Hexosamine / UDP-GlcNAc donor supply | Interconversion from UDP-GlcNAc to UDP-GalNAc donor pool. |
| UAGALDP | UDP-N-acetylgalactosamine diphosphorylase | UDP-GalNAc | Glycosylation donor supply | Hexosamine / UDP-GlcNAc donor supply | Direct activation of GalNAc-1-phosphate to UDP-GalNAc. |
| GALU | UTP-glucose-1-phosphate uridylyltransferase | UDP-Glucose | Glycosylation donor supply | UDP-Hexose / UDP-Gal donor supply | Direct activation of glucose-1-phosphate to UDP-glucose. |
| GALT | BiGGRxn35 | UDP-Galactose | Glycosylation donor supply | UDP-Hexose / UDP-Gal donor supply | Direct activation route from galactose-1-phosphate and UTP to UDP-galactose. |
| UGLT | UDPglucose--hexose-1-phosphate uridylyltransferase | UDP-Hexose bridge | Glycosylation donor supply | UDP-Hexose / UDP-Gal donor supply | UDP-glucose and UDP-galactose donor-pool interconversion bridge. |
| UDPG4E | UDPglucose 4-epimerase | UDP-Galactose bridge | Glycosylation donor supply | UDP-Hexose / UDP-Gal donor supply | Interconversion from UDP-glucose to UDP-galactose. |
| UDPGD | UDPglucose 6-dehydrogenase | UDP-GlcA | Glycosylation donor supply | UDP-GlcA donor supply | Direct formation of UDP-glucuronate from UDP-glucose. |
| MAN1PT2 | mannose-1-phosphate guanylyltransferase (GDP) | GDP-Mannose | Glycosylation donor supply | GDP-Mannose donor supply | Direct activation of mannose-1-phosphate to GDP-mannose. |
| r0208 | GTP:alpha-D-mannose-1-phosphate guanylyltransferase | GDP-Mannose | Glycosylation donor supply | GDP-Mannose donor supply | Direct GTP-dependent GDP-mannose activation route. |
| GMAND | GDP-D-mannose dehydratase | GDP-Fucose precursor | Glycosylation donor supply | GDP-Fucose donor supply | Handoff from GDP-mannose toward GDP-fucose donor synthesis. |
| GFUCS | GDP-L-fucose synthase | GDP-Fucose | Glycosylation donor supply | GDP-Fucose donor supply | Direct GDP-fucose donor-supply formation. |
| r0782 | GDP-L-fucose:NADP+ 4-oxidoreductase (3,5-epimerizing) | GDP-Fucose | Glycosylation donor supply | GDP-Fucose donor supply | Alternate redox route to GDP-fucose donor supply. |
| F1PGT | fucose-1-phosphate guanylyltransferase | GDP-Fucose salvage | Glycosylation donor supply | GDP-Fucose donor supply | Salvage activation of fucose-1-phosphate to GDP-fucose. |
| CMPSAS | CMP sialic acid synthase | CMP-Sialic acid | Glycosylation donor supply | CMP-Sialic acid donor supply | Direct activation of N-acetylneuraminate to CMP-sialic acid. |
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

Requested full-FVA validation command:

```bash
python run_pipeline.py --dataset practice_20aa --input_format tsv --steps 10,11,12,13,14 --fva_source full
```

Full-FVA outputs were not available in the accessible local workspace during
this implementation session, so workstation full-FVA validation was not run
here.

Mapping integrity was checked against iCHO3K reaction identifiers, reaction
names, subsystem annotations, and metabolite-level donor-supply roles. The
expected post-validation behavior is:

- `Glycosylation donor supply` mappings increase from 0 to 16.
- `chompact_domain_coverage_audit.csv` should show improved nucleotide-sugar
  donor/glycosylation mapping coverage when full FVA is supplied.
- These reactions should remain guardrailed as model-emergent and full-FVA
  dependent unless measured evidence is added separately.

Full-FVA workstation validation is required before merging this PR as a
validated release branch.
