#!/usr/bin/env python
"""Run a minimal COBRApy FBA workflow on local iCHO3K.

Run this inside the iCHO3K environment, not the bundled Codex runtime:

conda env create -f C:\CHO_FBA\model\iCHO3K-main\environment.yml
conda activate icho3k
python pipeline/run_icho3k_cobra_fba.py --constraints outputs/multiomics_model_inputs_icho3k_verified/icho_exchange_constraints.csv --outdir outputs/icho3k_fba

The constraints file is produced by prepare_multiomics_model_inputs.py.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODEL = ROOT / "models" / "iCHO3K" / "Model" / "iCHO3K_cho_prod_generic_unblocked.json"


def load_model(model_path: Path):
    try:
        import cobra
    except ImportError as exc:
        raise SystemExit(
            "COBRApy is not installed in this Python environment. "
            "Use the iCHO3K conda environment from models\\iCHO3K\\env\\environment.yml."
        ) from exc
    return cobra.io.load_json_model(str(model_path))


def apply_constraints(model, constraints: pd.DataFrame, condition: str | None = None):
    sub = constraints.copy()
    if condition and "passage_or_clone" in sub.columns:
        sub = sub[sub["passage_or_clone"].astype(str) == condition]
    applied = []
    for row in sub.itertuples(index=False):
        rid = getattr(row, "model_exchange_reaction_id")
        if pd.isna(rid) or rid not in model.reactions:
            continue
        rxn = model.reactions.get_by_id(rid)
        lb = getattr(row, "lower_bound")
        ub = getattr(row, "upper_bound")
        if pd.notna(lb):
            rxn.lower_bound = float(lb)
        if pd.notna(ub):
            rxn.upper_bound = float(ub)
        applied.append(rid)
    return sorted(set(applied))


def optimize_objective(model, objective_id: str):
    if objective_id not in model.reactions:
        return {"objective": objective_id, "status": "missing", "value": None}
    model.objective = objective_id
    sol = model.optimize()
    return {"objective": objective_id, "status": sol.status, "value": sol.objective_value}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--constraints", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--condition", help="Optional passage_or_clone value to constrain one condition")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    model = load_model(args.model)
    constraints = pd.read_csv(args.constraints)
    applied = apply_constraints(model, constraints, args.condition)

    objective_results = [
        optimize_objective(model.copy(), "biomass_cho_prod"),
        optimize_objective(model.copy(), "DM_igg_g"),
        optimize_objective(model.copy(), "igg_formation"),
    ]
    pd.DataFrame(objective_results).to_csv(args.outdir / "fba_objective_results.csv", index=False)
    pd.DataFrame({"applied_exchange_reaction_id": applied}).to_csv(args.outdir / "applied_constraints.csv", index=False)
    print(pd.DataFrame(objective_results).to_string(index=False))
    print(f"Applied {len(applied)} exchange constraints. Wrote results to {args.outdir}")


if __name__ == "__main__":
    main()
