# CHO Metabolomics MASTER

Last Updated: 2026-06-04
Source of Truth: GitHub `main` branch + this `MASTER.md`
Status: ACTIVE

---

## Project Objective

Identify biologically meaningful metabolic differences between CHO clones/groups using experimentally constrained genome-scale metabolic modeling.

This repository is the final offline-friendly CHO FBA/FVA pipeline.

---

## Official Environment

Environment: `metabolomics_env_py313`
Python: 3.13
COBRApy: 0.30.0
Solver: GLPK

---

## Scientific Strategy

Calculation Layer:

```text
iCHO3K
```

Interpretation Layer:

```text
CHOmpact
```

Project consensus from recent development discussions:

- Use iCHO3K for all calculations.
- Use CHOmpact for biological interpretation and figure generation.

---

## Official Workflow

```text
Metabolite Data
→ Exchange Rate Calculation
→ iCHO3K Constraints
→ FBA
→ pFBA
→ FVA
→ Flux Comparison
→ CHOmpact Aggregation
→ Figure Generation
```

---

## FVA Scope Definitions

```text
exchange  = uptake/secretion reactions only
focused   = curated central metabolism + mAb pathway panel
internal  = all non-exchange internal reactions
all       = every reaction in iCHO3K
```

---

## Current Recommended Usage

Fast focused workflow:

- Development
- Presentation figures
- Focused Escher maps

Workstation workflow:

- Internal genome-scale FVA
- High vs Low producer comparison
- Pathway interpretation

---

## Important Interpretation Rule

FVA ranges are model-feasible ranges.
They are not directly measured intracellular fluxes.

Prioritize:

- focused pathway reactions
- robust High/Low separation
- consistency with measured exchange phenotypes

---

## Workspace Rules

- GitHub main = official project state.
- MASTER.md = lightweight consensus.
- Workspace chats are disposable.
- Use workspace rotation when chats become slow.

Example:

```text
METABOLOMICS_WORKSPACE_01
METABOLOMICS_WORKSPACE_02
METABOLOMICS_WORKSPACE_03
```

---

## New Chat Start Template

```text
Metabolomics 프로젝트.
GitHub MASTER.md 기준으로 진행.
오늘 작업:
- ...
```
