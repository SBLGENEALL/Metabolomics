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


def reaction_class(reaction_id: str) -> str:
    if reaction_id.startswith("EX_"):
        return "exchange_media"
    if reaction_id.startswith("DM_"):
        return "demand_sink"
    if "igg" in reaction_id.lower():
        return "antibody_product"
    return "internal_network"


def reaction_metadata(model, reaction_id: str) -> dict:
    rxn = model.reactions.get_by_id(reaction_id)
    return {
        "reaction_id": reaction_id,
        "reaction_name": rxn.name,
        "reaction_class": reaction_class(reaction_id),
        "subsystem": getattr(rxn, "subsystem", ""),
        "lower_bound": rxn.lower_bound,
        "upper_bound": rxn.upper_bound,
        "reaction_equation": rxn.reaction,
    }


def reaction_fluxes(model, solution, scenario: dict, objective: str, reaction_ids: list[str]) -> list[dict]:
    rows = []
    if solution is None or solution.status != "optimal":
        return rows
    for rid in reaction_ids:
        if rid not in model.reactions:
            continue
        reduced_cost = math.nan
        if hasattr(solution, "reduced_costs"):
            reduced_cost = float(solution.reduced_costs.get(rid, math.nan))
        rows.append(
            {
                **scenario,
                "objective": objective,
                "objective_value": solution.objective_value,
                **reaction_metadata(model, rid),
                "flux": float(solution.fluxes.get(rid, math.nan)),
                "abs_flux": abs(float(solution.fluxes.get(rid, math.nan))),
                "reduced_cost": reduced_cost,
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
                **reaction_metadata(model, rid),
                "minimum": row["minimum"],
                "maximum": row["maximum"],
                "range": row["maximum"] - row["minimum"],
                "abs_range": abs(row["maximum"] - row["minimum"]),
                "is_fixed_or_tight": abs(row["maximum"] - row["minimum"]) < 1e-9,
            }
        )
    return rows


def write_flux_differences(
    outdir: Path,
    all_flux_df: pd.DataFrame,
    reference_group: str | None,
    compare_group: str | None,
) -> pd.DataFrame:
    if (
        all_flux_df.empty
        or not reference_group
        or not compare_group
        or "producer_group" not in all_flux_df.columns
    ):
        return pd.DataFrame()
    sub = all_flux_df[all_flux_df["producer_group"].astype(str).isin([reference_group, compare_group])].copy()
    if sub.empty:
        return pd.DataFrame()
    keys = ["objective", "reaction_id", "reaction_name", "reaction_class", "subsystem", "reaction_equation"]
    means = sub.groupby(["producer_group", *keys], dropna=False)["flux"].mean().reset_index()
    pivot = means.pivot_table(index=keys, columns="producer_group", values="flux", aggfunc="mean").reset_index()
    if reference_group not in pivot.columns or compare_group not in pivot.columns:
        return pd.DataFrame()
    pivot["reference_group"] = reference_group
    pivot["compare_group"] = compare_group
    pivot["flux_reference"] = pivot[reference_group]
    pivot["flux_compare"] = pivot[compare_group]
    pivot["flux_delta_compare_minus_reference"] = pivot["flux_compare"] - pivot["flux_reference"]
    pivot["abs_flux_delta"] = pivot["flux_delta_compare_minus_reference"].abs()
    result = pivot.drop(columns=[reference_group, compare_group]).sort_values("abs_flux_delta", ascending=False)
    result.to_csv(outdir / "fba_flux_differences_by_group.csv", index=False)

    design = result[result["objective"].isin(["DM_igg_g", "igg_formation"])].copy()
    design["design_hint"] = design["reaction_class"].map(
        {
            "exchange_media": "media/feed candidate",
            "demand_sink": "product sink/objective check",
            "antibody_product": "IgG synthesis/assembly candidate",
            "internal_network": "intracellular pathway candidate",
        }
    ).fillna("model candidate")
    design.head(300).to_csv(outdir / "fba_flux_design_candidates.csv", index=False)
    return result


def write_fva_summary(outdir: Path, fva_df: pd.DataFrame) -> None:
    if fva_df.empty:
        return
    summary = (
        fva_df.groupby(["reaction_id", "reaction_name", "reaction_class", "subsystem"], dropna=False)
        .agg(
            mean_minimum=("minimum", "mean"),
            mean_maximum=("maximum", "mean"),
            mean_range=("range", "mean"),
            max_abs_range=("abs_range", "max"),
            tight_scenario_fraction=("is_fixed_or_tight", "mean"),
        )
        .reset_index()
        .sort_values(["tight_scenario_fraction", "max_abs_range"], ascending=[False, True])
    )
    summary.to_csv(outdir / "fva_reaction_range_summary.csv", index=False)


def choose_fva_reactions(model, scenario_flux_rows: list[dict], scope: str, threshold: float, max_reactions: int | None) -> list[str]:
    if scope == "all":
        selected = [rxn.id for rxn in model.reactions]
    elif scope == "active":
        selected = sorted(
            {
                row["reaction_id"]
                for row in scenario_flux_rows
                if abs(float(row.get("flux", 0.0))) > threshold
            }.union(CORE_REACTIONS)
        )
    else:
        selected = list(CORE_REACTIONS)

    if scope == "active" and max_reactions is None:
        max_reactions = 500

    if max_reactions and len(selected) > max_reactions:
        ranked = (
            pd.DataFrame(scenario_flux_rows)
            .assign(abs_flux=lambda df: df["flux"].abs())
            .groupby("reaction_id", dropna=False)["abs_flux"]
            .max()
            .sort_values(ascending=False)
        )
        keep = list(dict.fromkeys([*CORE_REACTIONS, *ranked.index.tolist()]))
        selected = [rid for rid in keep if rid in selected][:max_reactions]
    return selected


def write_fba_figures(
    outdir: Path,
    objective_df: pd.DataFrame,
    selected_flux_df: pd.DataFrame,
    fva_df: pd.DataFrame,
    diff_df: pd.DataFrame,
) -> None:
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
    if not selected_flux_df.empty:
        core = selected_flux_df[selected_flux_df["reaction_id"].isin(["EX_glc_e", "EX_lac_L_e", "EX_nh4_e", "EX_gln_L_e", "DM_igg_g", "igg_formation"])]
        html_parts.append("<section><h2>Selected fluxes</h2>")
        html_parts.append(core.head(200).to_html(index=False, float_format=lambda x: f"{x:.4g}"))
        html_parts.append("</section>")
    if not fva_df.empty:
        html_parts.append("<section><h2>FVA ranges</h2>")
        html_parts.append(fva_df.head(200).to_html(index=False, float_format=lambda x: f"{x:.4g}"))
        html_parts.append("</section>")
    if not diff_df.empty:
        html_parts.append("<section><h2>Largest model-predicted flux differences</h2>")
        html_parts.append(diff_df.head(80).to_html(index=False, float_format=lambda x: f"{x:.4g}"))
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
    parser.add_argument("--full-fva", action="store_true", help="Run FVA on every model reaction instead of the core reaction panel")
    parser.add_argument("--fva-scope", choices=["core", "active", "all"], default="core", help="FVA reaction set: core is fastest, active follows nonzero FBA reactions, all is full-model FVA")
    parser.add_argument("--active-flux-threshold", type=float, default=1e-9, help="Flux threshold for --fva-scope active")
    parser.add_argument("--fva-max-reactions", type=int, help="Optional cap for FVA reactions after ranking by absolute FBA flux; active scope defaults to 500")
    parser.add_argument("--reference-group", help="Reference producer_group for model-predicted flux difference calculations")
    parser.add_argument("--compare-group", help="Comparison producer_group for model-predicted flux difference calculations")
    args = parser.parse_args()
    if args.full_fva:
        args.fva_scope = "all"

    args.outdir.mkdir(parents=True, exist_ok=True)
    base_model = load_model(args.model)
    constraints = pd.read_csv(args.constraints)
    if args.condition and "passage_or_clone" in constraints.columns:
        constraints = constraints[constraints["passage_or_clone"].astype(str) == args.condition]

    scenario_cols = scenario_columns(constraints)
    grouped = [((), constraints)] if not scenario_cols else list(constraints.groupby(scenario_cols, dropna=False))
    objective_rows = []
    all_flux_rows = []
    fva_rows = []
    applied_rows = []

    for keys, sub in grouped:
        if not isinstance(keys, tuple):
            keys = (keys,)
        scenario = dict(zip(scenario_cols, keys))
        scenario["scenario"] = scenario_label(scenario)
        model = base_model.copy()
        applied = apply_constraints(model, sub)
        scenario_flux_rows = []
        for rid in applied:
            applied_rows.append({**scenario, "applied_exchange_reaction_id": rid})

        for objective in OBJECTIVES:
            objective_model = model.copy()
            status, value, solution = optimize_objective(objective_model, objective)
            objective_rows.append({**scenario, "objective": objective, "status": status, "value": value})
            rows = reaction_fluxes(
                objective_model,
                solution,
                scenario,
                objective,
                [rxn.id for rxn in objective_model.reactions],
            )
            scenario_flux_rows.extend(rows)
            all_flux_rows.extend(rows)

        if not args.skip_fva:
            fva_model = model.copy()
            if "DM_igg_g" in fva_model.reactions:
                fva_model.objective = "DM_igg_g"
            fva_reactions = choose_fva_reactions(
                fva_model,
                scenario_flux_rows,
                args.fva_scope,
                args.active_flux_threshold,
                args.fva_max_reactions,
            )
            fva_rows.extend(run_fva(fva_model, scenario, fva_reactions))

    objective_df = pd.DataFrame(objective_rows)
    all_flux_df = pd.DataFrame(all_flux_rows)
    selected_flux_df = (
        all_flux_df[all_flux_df["reaction_id"].isin(CORE_REACTIONS)].copy()
        if not all_flux_df.empty
        else pd.DataFrame()
    )
    fva_df = pd.DataFrame(fva_rows)
    selected_fva_df = (
        fva_df[fva_df["reaction_id"].isin(CORE_REACTIONS)].copy()
        if not fva_df.empty
        else pd.DataFrame()
    )
    applied_df = pd.DataFrame(applied_rows)
    diff_df = write_flux_differences(args.outdir, all_flux_df, args.reference_group, args.compare_group)
    write_fva_summary(args.outdir, fva_df)

    objective_df.to_csv(args.outdir / "fba_objective_results_by_scenario.csv", index=False)
    selected_flux_df.to_csv(args.outdir / "fba_selected_fluxes_by_scenario.csv", index=False)
    all_flux_df.to_csv(args.outdir / "fba_all_reaction_fluxes_by_scenario.csv", index=False)
    selected_fva_df.to_csv(args.outdir / "fva_selected_reactions_by_scenario.csv", index=False)
    if args.fva_scope == "active":
        fva_df.to_csv(args.outdir / "fva_active_reactions_by_scenario.csv", index=False)
    if args.fva_scope == "all":
        fva_df.to_csv(args.outdir / "fva_all_reactions_by_scenario.csv", index=False)
    applied_df.to_csv(args.outdir / "applied_constraints_by_scenario.csv", index=False)
    write_fba_figures(args.outdir, objective_df, selected_flux_df, selected_fva_df, diff_df)

    # Keep backwards-compatible filename for older instructions.
    objective_df.to_csv(args.outdir / "fba_objective_results.csv", index=False)
    print(objective_df.head(20).to_string(index=False))
    print(f"Ran {len(grouped)} scenario(s). Wrote full FBA flux results to {args.outdir}")
    if not args.skip_fva:
        print(f"FVA scope: {args.fva_scope}. Use --fva-scope active for broad practical FVA or --fva-scope all for full-model FVA.")


if __name__ == "__main__":
    main()
