"""
fba_utils.py — common utilities for CHO fed-batch FBA/FVA

2026-05 objective correction
---------------------------
Main analyses now maximize IgG/mAb production (default DM_igg_g), not a
correlation-selected arbitrary reaction. Optional biomass coupling is supported
through BIOMASS_MIN_FRACTION but is disabled by default.
"""
import os
import sys
import json
from typing import Dict, Iterable, Tuple, Optional, Any


def _find_root():
    cur = os.path.dirname(os.path.abspath(__file__))
    for _ in range(6):
        if os.path.exists(os.path.join(cur, "src", "config.py")):
            return cur
        cur = os.path.dirname(cur)
    return cur


ROOT = _find_root()
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from src.config import *  # noqa
import cobra
import pandas as pd
import numpy as np
from cobra.flux_analysis import flux_variability_analysis, pfba


def load_model(verbose: bool = True):
    """Load the iCHO3K production model and set the default production objective."""
    if not MODEL_PATH or not os.path.exists(MODEL_PATH):
        raise FileNotFoundError(
            f"Model not found: {MODEL_PATH}\n"
            "Expected iCHO3K files under model/iCHO3K/."
        )
    m = cobra.io.load_json_model(MODEL_PATH) if MODEL_PATH.endswith(".json") else cobra.io.read_sbml_model(MODEL_PATH)
    obj = resolve_objective(m)
    m.objective = obj
    if verbose:
        print(f"  Model : {os.path.basename(MODEL_PATH)}")
        print(f"  Rxns  : {len(m.reactions)} | Genes: {len(m.genes)}")
        print(f"  Obj   : {obj} (IgG/mAb production)")
    return m


def resolve_objective(model, preferred: Optional[str] = None) -> str:
    """Resolve a valid mAb production objective reaction id in the loaded model."""
    rxn_ids = {r.id for r in model.reactions}
    candidates = []
    if preferred:
        candidates.append(preferred)
    candidates += [OBJ_RXN, PRODUCTION_OBJECTIVE]
    candidates += list(PRODUCTION_OBJECTIVE_CANDIDATES)
    # de-duplicate while preserving order
    seen = set()
    ordered = []
    for x in candidates:
        if x and x not in seen:
            ordered.append(x)
            seen.add(x)
    for rid in ordered:
        if rid in rxn_ids:
            return rid
    raise ValueError(
        "No valid production objective found. Tried: " + ", ".join(ordered)
    )


def is_old_bad_objective(rid: Optional[str]) -> bool:
    """Return True for previous auto-screened non-IgG objective ids that should be ignored."""
    if not rid:
        return False
    if rid in PRODUCTION_OBJECTIVE_CANDIDATES:
        return False
    if rid == BIOMASS_RXN:
        return False
    # BiGGRxn05 was the problematic auto-picked reaction in the previous version.
    return rid.startswith("BiGGRxn") or rid.startswith("BiGGEx")


def choose_objective_from_data(model, data: dict, requested: Optional[str] = None) -> str:
    """Choose objective, ignoring stale selected_objective values from old auto-screening."""
    if requested:
        return resolve_objective(model, requested)
    saved = data.get("selected_objective") or data.get("objective")
    if saved and not is_old_bad_objective(saved):
        try:
            return resolve_objective(model, saved)
        except Exception:
            pass
    return resolve_objective(model)


def apply_bounds(model, constraints: dict) -> int:
    """Apply exchange constraints safely, avoiding cobra bound-order conflicts."""
    rxn_ids = {r.id for r in model.reactions}
    n = 0
    for key, (lb, ub) in (constraints or {}).items():
        rxn_id = EXCHANGE_IDS.get(key, key)
        if rxn_id not in rxn_ids:
            continue
        rxn = model.reactions.get_by_id(rxn_id)
        lb_f, ub_f = float(lb), float(ub)
        if lb_f > ub_f:
            lb_f, ub_f = ub_f, lb_f
        if lb_f > rxn.upper_bound:
            rxn.upper_bound = ub_f
            rxn.lower_bound = lb_f
        elif ub_f < rxn.lower_bound:
            rxn.lower_bound = lb_f
            rxn.upper_bound = ub_f
        else:
            rxn.lower_bound = lb_f
            rxn.upper_bound = ub_f
        n += 1
    return n


def fix_one_sided(constraints: dict) -> dict:
    """Convert same-sign exchange bounds to directional uptake/secretion bounds."""
    fixed = {}
    for rxn_id, (lb, ub) in (constraints or {}).items():
        if lb < 0 and ub < 0:
            fixed[rxn_id] = (lb, 0.0)
        elif lb > 0 and ub > 0:
            fixed[rxn_id] = (0.0, ub)
        else:
            fixed[rxn_id] = (lb, ub)
    return fixed


def rates_to_constraints(rates: dict, buffer: float = 0.20) -> dict:
    """Convert exchange rates to directional FBA bounds."""
    cst = {}
    for key, q in (rates or {}).items():
        rxn_id = EXCHANGE_IDS.get(key, key)
        if abs(q) < 1e-9:
            continue
        if q < 0:
            lb, ub = q * (1 + buffer), 0.0
        else:
            lb, ub = 0.0, q * (1 + buffer)
        cst[rxn_id] = (min(lb, ub), max(lb, ub))
    return fix_one_sided(cst)


def avg_constraints(clone_list: list, all_constraints: dict) -> dict:
    """Average constraints across a list of clones."""
    cst_list = [all_constraints[c] for c in clone_list if c in all_constraints]
    if not cst_list:
        return {}
    all_ids = set()
    for c in cst_list:
        all_ids.update(c.keys())
    return {
        rid: (
            float(np.mean([c[rid][0] for c in cst_list if rid in c])),
            float(np.mean([c[rid][1] for c in cst_list if rid in c])),
        )
        for rid in all_ids
    }


def set_biomass_minimum(model, constraints: dict, fraction: float) -> Tuple[float, float, str]:
    """
    Optionally impose biomass >= fraction * max biomass under the same constraints.
    Returns (biomass_capacity, imposed_lower_bound, status).
    """
    if not fraction or fraction <= 0 or BIOMASS_RXN not in {r.id for r in model.reactions}:
        return 0.0, 0.0, "not_used"

    with model:
        apply_bounds(model, constraints)
        model.objective = BIOMASS_RXN
        sol = model.optimize()
        if sol.status != "optimal" or sol.objective_value is None or sol.objective_value <= 1e-12:
            return 0.0, 0.0, sol.status
        cap = float(sol.objective_value)

    lb = cap * float(fraction)
    biomass_rxn = model.reactions.get_by_id(BIOMASS_RXN)
    biomass_rxn.lower_bound = max(float(biomass_rxn.lower_bound), lb)
    return cap, lb, "imposed"


def optimize_production(
    model,
    constraints: dict,
    objective: Optional[str] = None,
    biomass_fraction: Optional[float] = None,
    label: str = "",
    return_fluxes: bool = False,
) -> dict:
    """Run FBA with IgG/mAb production objective under optional biomass coupling."""
    obj_rxn = resolve_objective(model, objective)
    if biomass_fraction is None:
        biomass_fraction = BIOMASS_MIN_FRACTION
    with model:
        n = apply_bounds(model, constraints)
        bio_cap, bio_lb, bio_status = set_biomass_minimum(model, constraints, biomass_fraction)
        model.objective = obj_rxn
        sol = model.optimize()
        ok = sol.status == "optimal"
        obj_val = float(sol.objective_value) if ok and sol.objective_value is not None else 0.0
        fluxes = sol.fluxes.to_dict() if (ok and return_fluxes) else {}
        ex_flux = {}
        if ok:
            ex_flux = {
                r.id: float(sol.fluxes.get(r.id, 0.0))
                for r in model.exchanges
                if abs(float(sol.fluxes.get(r.id, 0.0))) > 1e-9
            }
    if label:
        print(
            f"  {label:32s} n={n:2d} {sol.status:10s} "
            f"{obj_rxn}={obj_val:.6g} bio_lb={bio_lb:.4g}"
        )
    return {
        "label": label,
        "objective": obj_rxn,
        "obj": obj_val,
        "status": sol.status,
        "n_constraints": n,
        "biomass_capacity": bio_cap,
        "biomass_lb": bio_lb,
        "biomass_status": bio_status,
        "ex_flux": ex_flux,
        "fluxes": fluxes,
    }


def run_fba(model, constraints: dict, label: str = "") -> dict:
    """Backward-compatible wrapper around optimize_production."""
    return optimize_production(model, constraints, label=label)


def get_fva_reaction_list(model, scope: str = "focused", focused_ids: list = None):
    """Return a COBRA reaction list for FVA.

    scope
    -----
    exchange : exchange/boundary uptake-secretion reactions only
    focused  : curated pathway IDs supplied by caller; if none, falls back to exchange
    internal : all non-boundary internal reactions
    all      : every reaction in the model
    """
    scope = (scope or "focused").lower()
    rxn_ids_set = {r.id for r in model.reactions}
    if scope == "exchange":
        return list(model.exchanges)
    if scope == "focused":
        if focused_ids:
            return [model.reactions.get_by_id(r) for r in focused_ids if r in rxn_ids_set]
        return list(model.exchanges)
    if scope == "internal":
        ex_ids = {r.id for r in model.exchanges}
        return [r for r in model.reactions if r.id not in ex_ids and not getattr(r, "boundary", False)]
    if scope == "all":
        return list(model.reactions)
    raise ValueError(f"Unknown FVA scope: {scope}")


def run_fva(model, constraints: dict, label: str = "", exchange_ids: list = None,
            scope: str = None, focused_ids: list = None, fraction: float = None,
            processes: int = None) -> pd.DataFrame:
    """Run FVA under the configured mAb production objective.

    Backward-compatible behavior: if exchange_ids are supplied, those are used as
    a focused list. If neither scope nor focused_ids is supplied, the default is
    the configured FVA_SCOPE_DEFAULT, which is 'focused' but falls back to
    exchange unless focused_ids are provided.
    """
    with model:
        apply_bounds(model, constraints)
        obj_rxn = resolve_objective(model)
        model.objective = obj_rxn
        pre = model.optimize()
        if pre.status != "optimal":
            print(f"  FVA {label}: infeasible before FVA")
            return pd.DataFrame()
        if exchange_ids and focused_ids is None:
            focused_ids = exchange_ids
            scope = "focused"
        if scope is None:
            scope = FVA_SCOPE_DEFAULT
        rxn_list = get_fva_reaction_list(model, scope=scope, focused_ids=focused_ids)
        try:
            fva = flux_variability_analysis(
                model,
                fraction_of_optimum=FVA_FRACTION if fraction is None else float(fraction),
                processes=FVA_PROCESSES if processes is None else int(processes),
                reaction_list=rxn_list,
            )
            fva["mean"] = (fva["minimum"] + fva["maximum"]) / 2
            fva["range"] = fva["maximum"] - fva["minimum"]
            fva["fva_scope"] = scope
            if label:
                print(f"  FVA {label}: OK ({len(fva)} reactions, scope={scope}), {obj_rxn}={pre.objective_value:.6g}")
            return fva
        except Exception as e:
            print(f"  FVA {label}: {e}")
            return pd.DataFrame()


def save_figure(fig, name: str, dataset: str, subdir: str = "figures") -> str:
    import matplotlib.pyplot as plt
    out = results_dir(dataset, subdir)
    path = os.path.join(out, name)
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"  [saved] results/{dataset}/{subdir}/{name}")
    return path


def save_table(df: pd.DataFrame, name: str, dataset: str) -> str:
    out = results_dir(dataset, "tables")
    path = os.path.join(out, name)
    if name.endswith(".csv"):
        df.to_csv(path, index=False)
    elif name.endswith(".xlsx"):
        df.to_excel(path, index=False)
    else:
        df.to_csv(path, index=False)
    print(f"  [saved] results/{dataset}/tables/{name}")
    return path


def ax_style(ax):
    ax.xaxis.grid(True, ls=":", color="0.88", zorder=0)
    ax.yaxis.grid(True, ls=":", color="0.88", zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)


def estimate_igg_interval_rates(data: dict, mw_g_per_mol: float = None) -> dict:
    """Estimate observed IgG production rate for the selected interval.

    Returns clone -> qIgG in mmol/gDCW/h using IgG titer increase, culture volume,
    and IVCD calculation consistent with 01_load_data.py.
    Assumes IgG/mAb MW ≈ 150 kDa unless configured otherwise.
    """
    mw = float(mw_g_per_mol or globals().get("IGG_MW_G_PER_MOL", 150000.0))
    df = data.get("df_raw")
    interval = data.get("day_interval")
    clones = data.get("clones", [])
    if df is None or interval is None:
        return {}
    d1, d2 = [int(x) for x in interval]
    out = {}
    for clone in clones:
        sub = df[df["Sample ID"] == clone].sort_values("DAY")
        rows = {int(r["DAY"]): r for _, r in sub.iterrows()}
        if d1 not in rows or d2 not in rows:
            continue
        r1, r2 = rows[d1], rows[d2]
        try:
            vcd1 = float(r1["Viable Density"])
            vcd2 = float(r2["Viable Density"])
            v1 = float(r1.get("Culture Volume mL", 50))
            v2 = float(r2.get("Culture Volume mL", 50))
            v_avg = (v1 + v2) / 2.0
            dt_h = (d2 - d1) * 24.0
            ivcd = ((vcd1 + vcd2) / 2.0) * 1e6 * v_avg * dt_h * 8e-12
            if ivcd <= 0:
                continue
            igg1 = float(r1.get("IgG", 0.0) if pd.notna(r1.get("IgG", 0.0)) else 0.0)
            igg2 = float(r2.get("IgG", 0.0) if pd.notna(r2.get("IgG", 0.0)) else 0.0)
            delta_mg = max(0.0, (igg2 - igg1) * (v_avg / 1000.0))  # mg in culture
            q_mmol_gdcwh = (delta_mg / mw) / ivcd if ivcd > 0 else 0.0
            out[clone] = q_mmol_gdcwh
        except Exception:
            continue
    return out


def optimize_with_measured_igg_demand(
    model,
    constraints: dict,
    measured_igg_rate: float,
    demand_fraction: float = None,
    demand_rxn: str = "DM_igg_g",
    secondary_objective: str = None,
    biomass_fraction: Optional[float] = None,
    label: str = "",
    return_fluxes: bool = False,
) -> dict:
    """Check feasibility while forcing observed IgG production as a demand.

    This is the recommended mode when maximizing DM_igg_g gives zero. It does
    not pretend to predict IgG production; it asks whether the measured qIgG can
    be supported by the model and then optimizes a secondary objective, usually
    biomass_cho_prod, under the measured production demand.
    """
    if demand_fraction is None:
        demand_fraction = float(globals().get("MEASURED_IGG_DEMAND_FRACTION", 0.80))
    if secondary_objective is None:
        secondary_objective = BIOMASS_RXN if BIOMASS_RXN in {r.id for r in model.reactions} else demand_rxn
    target = max(0.0, float(measured_igg_rate or 0.0) * float(demand_fraction))
    with model:
        n = apply_bounds(model, constraints)
        rxn_ids = {r.id for r in model.reactions}
        if demand_rxn not in rxn_ids:
            return {
                "label": label, "mode": "measured_igg_demand_missing", "objective": secondary_objective,
                "production_rxn": demand_rxn, "obj": 0.0, "production_flux": 0.0,
                "measured_igg_rate": measured_igg_rate, "enforced_igg_lb": target,
                "status": "missing_demand_rxn", "n_constraints": n, "fluxes": {}, "ex_flux": {},
                "biomass_capacity": 0.0, "biomass_lb": 0.0, "biomass_status": "not_used",
            }
        dr = model.reactions.get_by_id(demand_rxn)
        dr.lower_bound = max(float(dr.lower_bound), target)
        if dr.upper_bound < dr.lower_bound:
            dr.upper_bound = max(1000.0, dr.lower_bound * 10.0)
        bio_cap, bio_lb, bio_status = set_biomass_minimum(model, constraints, biomass_fraction or 0.0)
        if secondary_objective in rxn_ids:
            model.objective = secondary_objective
        else:
            model.objective = demand_rxn
        sol = model.optimize()
        ok = sol.status == "optimal"
        prod_flux = float(sol.fluxes.get(demand_rxn, 0.0)) if ok else 0.0
        obj_val = float(sol.objective_value) if ok and sol.objective_value is not None else 0.0
        fluxes = sol.fluxes.to_dict() if (ok and return_fluxes) else {}
        ex_flux = {}
        if ok:
            ex_flux = {
                r.id: float(sol.fluxes.get(r.id, 0.0))
                for r in model.exchanges
                if abs(float(sol.fluxes.get(r.id, 0.0))) > 1e-12
            }
    if label:
        print(
            f"  {label:32s} n={n:2d} {sol.status:10s} "
            f"forced {demand_rxn}>={target:.3e}; flux={prod_flux:.3e}"
        )
    return {
        "label": label,
        "mode": "measured_igg_demand",
        "objective": secondary_objective,
        "production_rxn": demand_rxn,
        "obj": prod_flux,
        "secondary_obj": obj_val,
        "production_flux": prod_flux,
        "measured_igg_rate": float(measured_igg_rate or 0.0),
        "enforced_igg_lb": target,
        "status": sol.status,
        "n_constraints": n,
        "biomass_capacity": bio_cap,
        "biomass_lb": bio_lb,
        "biomass_status": bio_status,
        "ex_flux": ex_flux,
        "fluxes": fluxes,
    }



# ─────────────────────────────────────────────────────────────────────────────
# Measured IgG-demand constrained pFBA utilities
# ─────────────────────────────────────────────────────────────────────────────

def compute_measured_demand_scale(q_rates: dict, strict_caps: dict, requested: Any = "auto", safety: float = None) -> float:
    """Return one global scale factor for observed qIgG demands.

    Rationale
    ---------
    Some practice/synthetic datasets can report qIgG values that are larger than
    the maximum model DM_igg_g capacity because of unit assumptions or an upper
    bound in the model. Maximizing DM_igg_g then makes every clone hit the same
    cap, which is not interpretable.

    We instead preserve the *relative* measured qIgG differences by applying one
    global scale factor to all clones. In auto mode, the factor is chosen so the
    largest scaled qIgG is safely below the smallest clone-specific strict
    capacity. This makes fixed-demand pFBA feasible without pretending that FBA
    predicted the production rate.
    """
    if safety is None:
        safety = float(globals().get("MEASURED_IGG_DEMAND_SCALE_SAFETY", 0.90))
    if requested is None:
        requested = globals().get("MEASURED_IGG_DEMAND_SCALE", "auto")
    if isinstance(requested, str) and requested.lower() != "auto":
        try:
            return max(0.0, float(requested))
        except Exception:
            requested = "auto"
    if not isinstance(requested, str):
        try:
            return max(0.0, float(requested))
        except Exception:
            pass

    ratios = []
    for c, q in (q_rates or {}).items():
        q = float(q or 0.0)
        cap = float((strict_caps or {}).get(c, 0.0) or 0.0)
        if q > PRODUCTION_ZERO_EPS and cap > PRODUCTION_ZERO_EPS:
            ratios.append(cap / q)
    if not ratios:
        return 0.0
    return min(1.0, max(0.0, float(safety) * min(ratios)))


def optimize_fixed_igg_demand_pfba(
    model,
    constraints: dict,
    measured_igg_rate: float,
    demand_scale: float,
    demand_rxn: str = "DM_igg_g",
    biomass_fraction: Optional[float] = None,
    label: str = "",
    return_fluxes: bool = False,
) -> dict:
    """Fix scaled measured IgG demand and run parsimonious FBA.

    This is a *data-enforced flux-state* calculation, not a production prediction.
    The IgG flux is fixed to measured_igg_rate × demand_scale so that downstream
    FVA/Escher maps compare flux distributions under measured production demand.
    """
    target = max(0.0, float(measured_igg_rate or 0.0) * float(demand_scale or 0.0))
    rxn_ids = {r.id for r in model.reactions}
    if demand_rxn not in rxn_ids:
        demand_rxn = resolve_objective(model, demand_rxn)

    with model:
        n = apply_bounds(model, constraints)
        bio_cap, bio_lb, bio_status = set_biomass_minimum(model, constraints, biomass_fraction or 0.0)
        dr = model.reactions.get_by_id(demand_rxn)
        # Fix, not maximize. This avoids the model's DM_igg_g cap dominating all clones.
        dr.lower_bound = max(float(dr.lower_bound), target)
        dr.upper_bound = min(float(dr.upper_bound), target) if float(dr.upper_bound) >= target else target
        # If numerical/bound-order issues remain, force exactly target.
        dr.lower_bound = target
        dr.upper_bound = target
        model.objective = demand_rxn
        check = model.optimize()
        ok = check.status == "optimal"
        sol = check
        total_flux = np.nan
        if ok:
            try:
                sol = pfba(model, fraction_of_optimum=1.0)
                total_flux = float(np.sum(np.abs(sol.fluxes.values))) if sol.status == "optimal" else np.nan
            except Exception:
                total_flux = float(np.sum(np.abs(check.fluxes.values))) if ok else np.nan
                sol = check
        prod_flux = float(sol.fluxes.get(demand_rxn, 0.0)) if sol.status == "optimal" else 0.0
        fluxes = sol.fluxes.to_dict() if (sol.status == "optimal" and return_fluxes) else {}
        ex_flux = {}
        if sol.status == "optimal":
            ex_flux = {
                r.id: float(sol.fluxes.get(r.id, 0.0))
                for r in model.exchanges
                if abs(float(sol.fluxes.get(r.id, 0.0))) > 1e-12
            }
    if label:
        print(
            f"  {label:32s} n={n:2d} {sol.status:10s} "
            f"fixed {demand_rxn}={target:.3e}; total_flux={total_flux:.3g}"
        )
    return {
        "label": label,
        "mode": "measured_scaled_demand_pfba",
        "objective": "pFBA_total_flux_minimization",
        "production_rxn": demand_rxn,
        "obj": prod_flux,
        "secondary_obj": total_flux,
        "pfba_total_flux": total_flux,
        "production_flux": prod_flux,
        "measured_igg_rate": float(measured_igg_rate or 0.0),
        "demand_scale": float(demand_scale or 0.0),
        "scaled_igg_demand": target,
        "enforced_igg_lb": target,
        "enforced_igg_ub": target,
        "status": sol.status,
        "n_constraints": n,
        "biomass_capacity": bio_cap,
        "biomass_lb": bio_lb,
        "biomass_status": bio_status,
        "ex_flux": ex_flux,
        "fluxes": fluxes,
    }



# ─────────────────────────────────────────────────────────────────────────────
# No-IgG-input mechanistic feature utilities
# ─────────────────────────────────────────────────────────────────────────────


def load_product_aa_requirement() -> dict:
    """Load product-specific amino-acid requirement created by step 15.

    Returns a dict keyed by model metabolite id, e.g. "leu_L_c" -> residue count
    per product molecule. If no product sequence has been provided, returns {}.
    """
    try:
        path = os.path.join(DATA_PROCESSED, "product_aa_requirement.json")
        if not os.path.exists(path):
            return {}
        with open(path, "r", encoding="utf-8") as f:
            payload = json.load(f)
        req = payload.get("requirement_by_model_metabolite", {}) or {}
        return {str(k): float(v) for k, v in req.items() if float(v) > 0}
    except Exception:
        return {}

def compute_igg_resource_proxy(model, constraints: dict, objective: str = None) -> dict:
    """Estimate a no-IgG-input nutrient support proxy for IgG synthesis.

    Why this exists
    ---------------
    In many genome-scale CHO models, maximizing a demand reaction such as
    DM_igg_g can hit a model-internal cap or remain insensitive to extracellular
    measurements. That FBA maximum is still a useful *diagnostic*, but it should
    not be presented as a productivity prediction when it is flat.

    This function creates an ML-ready mechanistic feature from the measured
    exchange constraints only. It asks: based on measured amino-acid uptake
    capacities, how much IgG stoichiometry could be supported? It does NOT use
    measured IgG titer/qIgG and it is not a substitute for a trained ML model.
    """
    rxn_ids = {r.id for r in model.reactions}
    aa_met_to_ex = {
        "ala_L_c": "EX_ala_L_e", "arg_L_c": "EX_arg_L_e", "asn_L_c": "EX_asn_L_e",
        "asp_L_c": "EX_asp_L_e", "cys_L_c": "EX_cys_L_e", "gln_L_c": "EX_gln_L_e",
        "glu_L_c": "EX_glu_L_e", "gly_c": "EX_gly_e", "his_L_c": "EX_his_L_e",
        "ile_L_c": "EX_ile_L_e", "leu_L_c": "EX_leu_L_e", "lys_L_c": "EX_lys_L_e",
        "met_L_c": "EX_met_L_e", "phe_L_c": "EX_phe_L_e", "pro_L_c": "EX_pro_L_e",
        "ser_L_c": "EX_ser_L_e", "thr_L_c": "EX_thr_L_e", "trp_L_c": "EX_trp_L_e",
        "tyr_L_c": "EX_tyr_L_e", "val_L_c": "EX_val_L_e",
    }

    # Prefer product-specific amino-acid requirement from step 15. If no sequence
    # was supplied, fall back to the generic iCHO3K IgG HC/LC stoichiometry.
    aa_req = load_product_aa_requirement()
    requirement_source = "product_sequence_step15" if aa_req else "generic_iCHO3K_igg_hc_lc"
    if not aa_req:
        try:
            if "igg_hc" in rxn_ids:
                for met, coeff in model.reactions.get_by_id("igg_hc").metabolites.items():
                    mid = met.id
                    if mid in aa_met_to_ex and coeff < 0:
                        aa_req[mid] = aa_req.get(mid, 0.0) + 2.0 * abs(float(coeff))
            if "igg_lc" in rxn_ids:
                for met, coeff in model.reactions.get_by_id("igg_lc").metabolites.items():
                    mid = met.id
                    if mid in aa_met_to_ex and coeff < 0:
                        aa_req[mid] = aa_req.get(mid, 0.0) + 2.0 * abs(float(coeff))
        except Exception:
            aa_req = {}

    capacities = {}
    uptake_sum = 0.0
    essential_uptake_sum = 0.0
    essential = {"arg_L_c", "his_L_c", "ile_L_c", "leu_L_c", "lys_L_c", "met_L_c", "phe_L_c", "thr_L_c", "trp_L_c", "val_L_c"}
    for met_id, ex_id in aa_met_to_ex.items():
        lb, ub = (constraints or {}).get(ex_id, (0.0, 0.0))
        uptake = max(0.0, -float(lb))
        uptake_sum += uptake
        if met_id in essential:
            essential_uptake_sum += uptake
        req = float(aa_req.get(met_id, 0.0) or 0.0)
        if req > 0:
            capacities[met_id] = uptake / req

    if capacities:
        bottleneck_met, min_cap = min(capacities.items(), key=lambda kv: kv[1])
        median_cap = float(np.median(list(capacities.values())))
        mean_cap = float(np.mean(list(capacities.values())))
    else:
        bottleneck_met, min_cap, median_cap, mean_cap = "", 0.0, 0.0, 0.0

    glc_lb, glc_ub = (constraints or {}).get("EX_glc_e", (0.0, 0.0))
    lac_lb, lac_ub = (constraints or {}).get("EX_lac_L_e", (0.0, 0.0))
    glc_uptake = max(0.0, -float(glc_lb))
    lac_rate = float(lac_ub) if float(lac_ub) > 0 else float(lac_lb)
    lac_glc_ratio = lac_rate / glc_uptake if glc_uptake > 1e-12 else np.nan

    # This is intentionally a feature/proxy, not a predicted titer. It tends to
    # increase with amino-acid support and improves when lactate burden is lower.
    lactate_penalty = max(0.0, lac_rate)  # positive lactate secretion/accumulation
    carbon_efficiency_proxy = glc_uptake / (1.0 + lactate_penalty)

    return {
        "resource_limited_igg_proxy": float(min_cap),
        "aa_support_median_proxy": float(median_cap),
        "aa_support_mean_proxy": float(mean_cap),
        "aa_total_uptake": float(uptake_sum),
        "essential_aa_total_uptake": float(essential_uptake_sum),
        "aa_bottleneck_metabolite": bottleneck_met,
        "glucose_uptake_capacity": float(glc_uptake),
        "lactate_exchange_rate": float(lac_rate),
        "lac_glc_ratio_feature": float(lac_glc_ratio) if np.isfinite(lac_glc_ratio) else np.nan,
        "carbon_efficiency_proxy": float(carbon_efficiency_proxy),
        "product_requirement_source": requirement_source,
        "product_aa_requirement_total": float(sum(aa_req.values())) if aa_req else 0.0,
        "resource_proxy_note": "No measured IgG used. ML-ready feature only; not a stand-alone validated prediction.",
    }


def try_production_capacity_strategies(model, constraints: dict, objective: str = None, biomass_fraction: float = None) -> pd.DataFrame:
    """Try strict and relaxed capacity checks for the true IgG objective.

    This is diagnostic only. It helps explain why the true production objective
    can be zero under measured exchange bounds.
    """
    obj = resolve_objective(model, objective)

    def widen(base, buf):
        out = {}
        for k, (lb, ub) in (base or {}).items():
            if lb < 0:
                out[k] = (lb * (1 + buf), 0.0)
            elif ub > 0:
                out[k] = (0.0, ub * (1 + buf))
            else:
                out[k] = (lb, ub)
        return out

    core_ids = {EXCHANGE_IDS.get(k, k) for k in ["Glucose", "Lactate", "Glutamine", "Glutamate", "NH4+"]}
    strategies = [
        ("strict_measured", constraints or {}),
        ("relaxed_50pct", widen(constraints or {}, 0.5)),
        ("relaxed_100pct", widen(constraints or {}, 1.0)),
        ("core_exchange_only", {k: v for k, v in (constraints or {}).items() if k in core_ids}),
        ("no_measured_exchange", {}),
    ]
    rows = []
    for name, cst in strategies:
        res = optimize_production(
            model, cst, objective=obj, biomass_fraction=biomass_fraction or 0.0, label="", return_fluxes=False
        )
        rows.append({
            "strategy": name,
            "objective": obj,
            "obj": res.get("obj", 0.0),
            "status": res.get("status", ""),
            "n_constraints": res.get("n_constraints", 0),
            "biomass_lb": res.get("biomass_lb", 0.0),
        })
    return pd.DataFrame(rows)
