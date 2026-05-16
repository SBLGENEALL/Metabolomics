#!/usr/bin/env python
"""Run scenario-level COBRApy FBA/FVA workflow on local iCHO3K.

Run this inside the iCHO3K environment, not the bundled Codex runtime.
Use forward slashes in example paths so Python docstrings do not interpret
Windows backslashes as escape sequences:

conda env create -f models/iCHO3K/env/environment.yml
conda activate icho3k
python pipeline/run_icho3k_cobra_fba.py --constraints results/interactive_run/icho3k_inputs/icho_exchange_constraints.csv --outdir results/interactive_run/fba

The constraints file is produced by prepare_multiomics_model_inputs.py.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODEL = ROOT / "models" / "iCHO3K" / "Model" / "iCHO3K_cho_prod_generic_unblocked.json"
OBJECTIVES = ["biomass_cho_prod", "DM_igg_g", "igg_formation"]
CORE_REACTIONS = [
    "biomass_cho_prod",
    "DM_igg_g",
    "igg_formation",
    "igg_hc",
    "igg_lc",
    "EX_glc_e",
    "EX_lac_L_e",
    "EX_nh4_e",
    "EX_gln_L_e",
    "EX_glu_L_e",
    "EX_ala_L_e",
    "EX_asn_L_e",
    "EX_asp_L_e",
    "EX_arg_L_e",
    "EX_ser_L_e",
    "EX_gly_e",
    "EX_cys_L_e",
    "EX_met_L_e",
    "EX_leu_L_e",
    "EX_ile_L_e",
    "EX_val_L_e",
    "EX_lys_L_e",
    "EX_his_L_e",
    "EX_phe_L_e",
    "EX_tyr_L_e",
    "EX_trp_L_e",
    "EX_pro_L_e",
    "EX_thr_L_e",
    "EX_pyr_e",
    "EX_cit_e",
    "EX_ac_e",
]


def load_model(model_path: Path):
    try:
        import appdirs
        cache_dir = ROOT / ".cache" / "cobrapy"
        cache_dir.mkdir(parents=True, exist_ok=True)
        appdirs.user_cache_dir = lambda appname=None, appauthor=None, **kwargs: str(cache_dir)
    except ImportError:
        pass
    try:
        import cobra
    except ImportError as exc:
        raise SystemExit(
            "COBRApy is not installed in this Python environment. "
            "Use the iCHO3K conda environment from models\\iCHO3K\\env\\environment.yml."
        ) from exc
    try:
        cobra.Configuration().processes = 1
    except Exception:
        pass
    return cobra.io.load_json_model(str(model_path))


def scenario_columns(constraints: pd.DataFrame) -> list[str]:
    candidates = ["passage_or_clone", "producer_group", "replicate", "day_start", "day_end"]
    return [c for c in candidates if c in constraints.columns]


def apply_constraints(model, constraints: pd.DataFrame) -> list[str]:
    applied = []
    for row in constraints.itertuples(index=False):
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


def scenario_label(row: pd.Series | dict) -> str:
    parts = []
    for key in ["passage_or_clone", "producer_group", "replicate", "day_start", "day_end"]:
        if key in row and pd.notna(row[key]):
            parts.append(f"{key}={row[key]}")
    return "|".join(parts) if parts else "all_constraints"


def optimize_objective(model, objective_id: str) -> tuple[str, float | None, object | None]:
    if objective_id not in model.reactions:
        return "missing", None, None
    model.objective = objective_id
    sol = model.optimize()
    value = sol.objective_value if sol.status == "optimal" else None
    return sol.status, value, sol


def selected_fluxes(model, solution, scenario: dict, objective: str, reaction_ids: list[str]) -> list[dict]:
    rows = []
    if solution is None or solution.status != "optimal":
        return rows
    for rid in reaction_ids:
        if rid not in model.reactions:
            continue
        rows.append(
            {
                **scenario,
                "objective": objective,
                "reaction_id": rid,
                "reaction_name": model.reactions.get_by_id(rid).name,
                "flux": float(solution.fluxes.get(rid, math.nan)),
            }
        )
    return rows


def run_fva(model, scenario: dict, reaction_ids: list[str]) -> list[dict]:
    try:
        from cobra.flux_analysis import flux_variability_analysis
    except ImportError:
        return []
    valid = [rid for rid in reaction_ids if rid in model.reactions]
    if not valid:
        return []
    fva = flux_variability_analysis(model, reaction_list=valid, fraction_of_optimum=0.9, processes=1)
    rows = []
    for rid, row in fva.iterrows():
        rows.append(
            {
                **scenario,
                "reaction_id": rid,
                "reaction_name": model.reactions.get_by_id(rid).name,
                "minimum": row["minimum"],
                "maximum": row["maximum"],
                "range": row["maximum"] - row["minimum"],
            }
        )
    return rows


def write_fba_figures(outdir: Path, objective_df: pd.DataFrame, flux_df: pd.DataFrame, fva_df: pd.DataFrame) -> None:
    html_parts = [
        "<!doctype html><html><head><meta charset='utf-8'><title>iCHO3K FBA/FVA Report</title>",
        "<style>body{font-family:Arial,sans-serif;margin:28px;background:#f8fafc;color:#111827}"
        "section{background:white;border:1px solid #cbd5e1;margin:18px 0;padding:16px}"
        "table{border-collapse:collapse;width:100%;font-size:12px}td,th{border:1px solid #e5e7eb;padding:6px;text-align:right}"
        "td:first-child,th:first-child{text-align:left}.bar{height:16px;background:#2563eb;display:inline-block}.neg{background:#dc2626}</style></head><body>",
        "<h1>iCHO3K FBA/FVA Report</h1>",
    ]
    if not objective_df.empty:
        pivot = objective_df.pivot_table(index="scenario", columns="objective", values="value", aggfunc="mean").reset_index()
        html_parts.append("<section><h2>Objective values by scenario</h2>")
        html_parts.append(pivot.to_html(index=False, float_format=lambda x: f"{x:.4g}"))
        html_parts.append("</section>")
    if not flux_df.empty:
        core = flux_df[flux_df["reaction_id"].isin(["EX_glc_e", "EX_lac_L_e", "EX_nh4_e", "EX_gln_L_e", "DM_igg_g", "igg_formation"])]
        html_parts.append("<section><h2>Selected fluxes</h2>")
        html_parts.append(core.head(200).to_html(index=False, float_format=lambda x: f"{x:.4g}"))
        html_parts.append("</section>")
    if not fva_df.empty:
        html_parts.append("<section><h2>FVA ranges</h2>")
        html_parts.append(fva_df.head(200).to_html(index=False, float_format=lambda x: f"{x:.4g}"))
        html_parts.append("</section>")
    html_parts.append("</body></html>")
    (outdir / "fba_fva_report.html").write_text("\n".join(html_parts), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--constraints", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--condition", help="Optional passage_or_clone value to constrain one condition")
    parser.add_argument("--skip-fva", action="store_true", help="Skip FVA to make the run faster")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    base_model = load_model(args.model)
    constraints = pd.read_csv(args.constraints)
    if args.condition and "passage_or_clone" in constraints.columns:
        constraints = constraints[constraints["passage_or_clone"].astype(str) == args.condition]

    scenario_cols = scenario_columns(constraints)
    grouped = [((), constraints)] if not scenario_cols else list(constraints.groupby(scenario_cols, dropna=False))
    objective_rows = []
    flux_rows = []
    fva_rows = []
    applied_rows = []

    for keys, sub in grouped:
        if not isinstance(keys, tuple):
            keys = (keys,)
        scenario = dict(zip(scenario_cols, keys))
        scenario["scenario"] = scenario_label(scenario)
        model = base_model.copy()
        applied = apply_constraints(model, sub)
        for rid in applied:
            applied_rows.append({**scenario, "applied_exchange_reaction_id": rid})

        for objective in OBJECTIVES:
            objective_model = model.copy()
            status, value, solution = optimize_objective(objective_model, objective)
            objective_rows.append({**scenario, "objective": objective, "status": status, "value": value})
            flux_rows.extend(selected_fluxes(objective_model, solution, scenario, objective, CORE_REACTIONS))

        if not args.skip_fva:
            fva_model = model.copy()
            if "DM_igg_g" in fva_model.reactions:
                fva_model.objective = "DM_igg_g"
            fva_rows.extend(run_fva(fva_model, scenario, CORE_REACTIONS))

    objective_df = pd.DataFrame(objective_rows)
    flux_df = pd.DataFrame(flux_rows)
    fva_df = pd.DataFrame(fva_rows)
    applied_df = pd.DataFrame(applied_rows)

    objective_df.to_csv(args.outdir / "fba_objective_results_by_scenario.csv", index=False)
    flux_df.to_csv(args.outdir / "fba_selected_fluxes_by_scenario.csv", index=False)
    fva_df.to_csv(args.outdir / "fva_selected_reactions_by_scenario.csv", index=False)
    applied_df.to_csv(args.outdir / "applied_constraints_by_scenario.csv", index=False)
    write_fba_figures(args.outdir, objective_df, flux_df, fva_df)

    # Keep backwards-compatible filename for older instructions.
    objective_df.to_csv(args.outdir / "fba_objective_results.csv", index=False)
    print(objective_df.head(20).to_string(index=False))
    print(f"Ran {len(grouped)} scenario(s). Wrote FBA/FVA results to {args.outdir}")


if __name__ == "__main__":
    main()
