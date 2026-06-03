# Repository cleanup plan

## Final repository

Keep as final team distribution repository:

```text
SBLGENEALL/Metabolomics
```

## Repositories to archive manually

The current ChatGPT GitHub connector can create or update files but cannot toggle the GitHub archive setting. Archive these manually in GitHub settings:

```text
SBLGENEALL/Metabolomics_final
SBLGENEALL/CHO_metabolomics_FBA
SBLGENEALL/metabolomocs_fba
SBLGENEALL/Claude_Metabolomics
SBLGENEALL/Claude_ChatGPT_collaboration_Metabolomics
```

## Manual archive steps

For each old repository:

1. Open repository on GitHub.
2. Go to `Settings`.
3. Scroll to `Danger Zone`.
4. Click `Archive this repository`.
5. Confirm archive.

## Naming policy

Only one active repository should remain for this project:

```text
Metabolomics
```

Future major changes should be tracked by release tags:

```text
v1.0.0  workstation-ready FBA/FVA pipeline
v1.1.0  transcriptomics integration
v1.2.0  genomics integration
v2.0.0  multi-omics/ML extension
```
