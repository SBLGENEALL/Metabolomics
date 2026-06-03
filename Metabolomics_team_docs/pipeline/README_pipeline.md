# Pipeline steps

| Step | Purpose |
|---:|---|
| 17 | Data QC |
| 1 | Load raw data and calculate feed-corrected exchange rates |
| 2 | Map metabolite names to model reactions |
| 3 | Run FBA/pFBA diagnostic modes |
| 10 | Export Escher reaction data JSON |
| 14 | Focused central/mAb pathway FVA |
| 16 | Focused Escher map and focused flux/FVA JSON export |
| 21 | Internal/full FVA for group averages or clone-level targets |
| 19 | Pathway score, FVA overlap/separation, constraint sensitivity |
| 20 | Automated report generation |
| 6 | Figure generation |

Recommended default:

```bash
--steps 17,1,2,3,10,14,16,21,19,20,6
```

FVA scope:

```text
exchange  = exchange reactions only
focused   = curated central/mAb pathway panel
internal  = internal non-boundary reactions
all       = all model reactions
```

Recommended for workstation:

```bash
--fva_scope internal --fva_targets group_avg --fva_processes 16
```

Use `--fva_scope all` only when workstation performance is sufficient.
