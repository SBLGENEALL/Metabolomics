# iCHO3K Model Files

This folder contains the iCHO3K production generic unblocked model files copied
into the project so the workflow can run without depending on an external local
path.

Included files:

- `Model/iCHO3K_cho_prod_generic_unblocked.json`
- `Model/iCHO3K_cho_prod_generic_unblocked.xml`
- `Model/iCHO3K_cho_prod_generic_unblocked.mat`
- `env/environment.yml`
- `env/requirements.txt`

The workflow uses the JSON model by default. Use `env/environment.yml` or
`env/requirements.txt` to create a Python environment with COBRApy, libSBML, and
GLPK support for FBA/pFBA/FVA.

