"""
03_run_fba.py — split prediction and explanation modes

predict mode
------------
IgG titer/qIgG is NOT used as a model input. The fixed mAb objective
(default DM_igg_g) is maximized under measured exchange constraints. This is the
only mode that can be discussed as a production-capacity prediction.

explain mode
------------
Observed qIgG is scaled with one global factor, fixed as a DM_igg_g demand, and
pFBA is used to reconstruct a feasible flux state. This is explanation/flux
reconstruction, not prediction.
"""
import argparse, os, pickle, sys, warnings
warnings.filterwarnings("ignore")


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
from src.fba_utils import (
    load_model, save_table, optimize_production, choose_objective_from_data,
    estimate_igg_interval_rates, try_production_capacity_strategies,
    compute_measured_demand_scale, optimize_fixed_igg_demand_pfba,
    compute_igg_resource_proxy,
)
import pandas as pd
import numpy as np
from scipy import stats

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa", choices=["sowa2020", "practice_20aa", "own_experiment"])
parser.add_argument("--objective", default=os.environ.get("CHO_OBJECTIVE"), help="Fixed true production reaction id; default DM_igg_g")
parser.add_argument("--analysis_mode", default=os.environ.get("CHO_ANALYSIS_MODE", "both"), choices=["predict", "explain", "both"])
parser.add_argument("--biomass_fraction", type=float, default=None)
parser.add_argument("--demand_scale", default=os.environ.get("CHO_DEMAND_SCALE", MEASURED_IGG_DEMAND_SCALE), help="auto or numeric global scale for measured qIgG demand")
args = parser.parse_args()
DATASET = args.dataset
if args.biomass_fraction is None:
    args.biomass_fraction = float(os.environ.get("CHO_BIOMASS_FRACTION", BIOMASS_MIN_FRACTION))

print("=" * 65)
print(f"  03_run_fba.py — {DATASET}")
print("=" * 65)

pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
if not os.path.exists(pkl_path):
    print(f"  !! {pkl_path} 없음 — 02_map_metabolites.py 먼저 실행")
    sys.exit(1)
with open(pkl_path, "rb") as f:
    data = pickle.load(f)

all_constraints = data.get("all_constraints", {})
clones = data.get("clones", [])
igG_day14 = data.get("igG_day14", {})
high_clones = data.get("high_clones", [])
low_clones = data.get("low_clones", [])
if not all_constraints:
    print("  !! constraints 없음 — 02_map_metabolites.py 먼저 실행")
    sys.exit(1)

model = load_model(verbose=True)
OBJ = choose_objective_from_data(model, data, args.objective)
DEMAND_RXN = "DM_igg_g" if "DM_igg_g" in {r.id for r in model.reactions} else OBJ
igg_rates = estimate_igg_interval_rates(data)

print(f"\n  Fixed production objective : {OBJ}")
print(f"  Analysis mode              : {args.analysis_mode}")
print(f"  Biomass minimum fraction   : {args.biomass_fraction}")
print("  Rule: objective is never selected from correlation with IgG.")

prediction_results = {}
prediction_rows = []
capacity_rows = []
strict_by_clone = {}

if args.analysis_mode in ["predict", "both"]:
    print("\n  [Prediction mode] no-IgG-input diagnostics + ML-ready mechanistic features")
    print(f"  {'':2s} {'Clone':12s} {'status':10s} {'FBA_max_DM':>14s} {'IgG_day14':>10s}")
    print("  " + "─" * 60)
    for clone in clones:
        cst = all_constraints.get(clone, {})
        res = optimize_production(model, cst, objective=OBJ, biomass_fraction=args.biomass_fraction, return_fluxes=True)
        cap = float(res.get("obj", 0.0) or 0.0)
        strict_by_clone[clone] = cap
        proxy = compute_igg_resource_proxy(model, cst, objective=OBJ)
        tag = "★" if clone in high_clones else ("▼" if clone in low_clones else " ")
        print(f"  {tag} {clone:12s} {res['status']:10s} {cap:14.3e} {igG_day14.get(clone,0):10.0f}")
        row = {
            "analysis_type": "prediction",
            "clone": clone,
            "objective": OBJ,
            "production_rxn": OBJ,
            "obj": cap,
            # Backward-compatible column name. This is a diagnostic FBA capacity,
            # not necessarily a useful titer prediction if flat/capped.
            "predicted_capacity": cap,
            "fba_max_capacity_diagnostic": cap,
            "production_flux": cap,
            "measured_igg_rate": float(igg_rates.get(clone, 0.0) or 0.0),
            "uses_igg_as_input": False,
            "mode": "prediction_diagnostic_no_igg_input",
            "status": res.get("status"),
            "n_constraints": res.get("n_constraints", 0),
            "igG_day14": igG_day14.get(clone, 0.0),
            "biomass_capacity": res.get("biomass_capacity", 0.0),
            "biomass_lb": res.get("biomass_lb", 0.0),
            "biomass_status": res.get("biomass_status", ""),
            "group": "High" if clone in high_clones else ("Low" if clone in low_clones else "Mid"),
            **proxy,
        }
        prediction_rows.append(row)
        prediction_results[clone] = {**row, "fluxes": res.get("fluxes", {}), "ex_flux": res.get("ex_flux", {}), "constraints": cst}
        cap_df = try_production_capacity_strategies(model, cst, objective=OBJ, biomass_fraction=args.biomass_fraction)
        cap_df.insert(0, "clone", clone)
        capacity_rows.append(cap_df)
else:
    # Still compute strict caps because explain mode needs feasible global scaling.
    for clone in clones:
        res = optimize_production(model, all_constraints.get(clone, {}), objective=OBJ, biomass_fraction=args.biomass_fraction, return_fluxes=False)
        strict_by_clone[clone] = float(res.get("obj", 0.0) or 0.0)

prediction_df = pd.DataFrame(prediction_rows)
capacity_df = pd.concat(capacity_rows, ignore_index=True) if capacity_rows else pd.DataFrame()

# Predictability diagnostics: IgG appears here only as validation, never as input.
diag_rows = []
if len(prediction_df) > 0:
    ok = prediction_df[(prediction_df["status"] == "optimal") & prediction_df["predicted_capacity"].notna()]
    variance = float(ok["predicted_capacity"].var()) if len(ok) > 1 else 0.0
    can_corr = len(ok) > 2 and variance > 1e-18 and ok["igG_day14"].var() > 0
    if can_corr:
        r, p = stats.pearsonr(ok["igG_day14"], ok["predicted_capacity"])
    else:
        r, p = np.nan, np.nan
    proxy_col = "resource_limited_igg_proxy" if "resource_limited_igg_proxy" in ok.columns else None
    proxy_r, proxy_p = np.nan, np.nan
    if proxy_col and len(ok) > 2 and ok[proxy_col].var() > 0 and ok["igG_day14"].var() > 0:
        proxy_r, proxy_p = stats.pearsonr(ok["igG_day14"], ok[proxy_col])
    diag_rows.append({
        "analysis_type": "prediction",
        "objective": OBJ,
        "uses_igg_as_input": False,
        "n_clones": len(ok),
        "predicted_capacity_variance": variance,
        "pearson_r_vs_observed_igg": r,
        "pearson_p": p,
        "resource_proxy_r_vs_observed_igg": proxy_r,
        "resource_proxy_p": proxy_p,
        "is_fba_capacity_predictive_in_this_dataset": bool(can_corr and abs(r) >= 0.5),
        "warning": "fba_demand_capacity_flat_or_capped_use_ml_features_not_capacity_as_prediction" if not can_corr else "interpret_with_validation_only",
    })
    if can_corr:
        print(f"\n  Prediction-vs-observed validation: r={r:.3f}, p={p:.3g}")
    else:
        print("\n  Prediction-vs-observed validation: not meaningful (flat/low-variance predicted capacity).")

diagnostics_df = pd.DataFrame(diag_rows)

explanation_results = {}
explanation_rows = []
if args.analysis_mode in ["explain", "both"]:
    scale = compute_measured_demand_scale(
        igg_rates,
        strict_by_clone,
        requested=args.demand_scale,
        safety=MEASURED_IGG_DEMAND_SCALE_SAFETY,
    )
    print("\n  [Explanation mode] fix scaled measured qIgG demand and run pFBA")
    print(f"  Global measured-qIgG demand scale: {scale:.6g}")
    print("  NOTE: this is data-enforced flux reconstruction, not production prediction.")
    print(f"\n  {'':2s} {'Clone':12s} {'status':10s} {'strict_cap':>12s} {'qIgG_obs':>12s} {'scaled_dem':>12s} {'total_flux':>12s}")
    print("  " + "─" * 86)
    for clone in clones:
        cst = all_constraints.get(clone, {})
        q_obs = float(igg_rates.get(clone, 0.0) or 0.0)
        res = optimize_fixed_igg_demand_pfba(
            model, cst, measured_igg_rate=q_obs, demand_scale=scale,
            demand_rxn=DEMAND_RXN, biomass_fraction=args.biomass_fraction,
            return_fluxes=True,
        )
        cap = float(strict_by_clone.get(clone, 0.0) or 0.0)
        util = res.get("scaled_igg_demand", 0.0) / cap if cap > PRODUCTION_ZERO_EPS else np.nan
        tag = "★" if clone in high_clones else ("▼" if clone in low_clones else " ")
        print(f"  {tag} {clone:12s} {res['status']:10s} {cap:12.3e} {q_obs:12.3e} {res.get('scaled_igg_demand',0):12.3e} {res.get('pfba_total_flux',np.nan):12.3g}")
        row = {
            "analysis_type": "explanation",
            "clone": clone,
            "objective": "fixed_measured_qIgG_pFBA",
            "production_rxn": DEMAND_RXN,
            "obj": float(res.get("production_flux", 0.0) or 0.0),
            "production_flux": float(res.get("production_flux", 0.0) or 0.0),
            "strict_obj": cap,
            "strict_capacity": cap,
            "measured_igg_rate": q_obs,
            "demand_scale": scale,
            "scaled_igg_demand": float(res.get("scaled_igg_demand", 0.0) or 0.0),
            "capacity_utilization": util,
            "uses_igg_as_input": True,
            "mode": res.get("mode", "measured_scaled_demand_pfba"),
            "status": res.get("status"),
            "n_constraints": res.get("n_constraints", 0),
            "igG_day14": igG_day14.get(clone, 0.0),
            "biomass_capacity": res.get("biomass_capacity", 0.0),
            "biomass_lb": res.get("biomass_lb", 0.0),
            "biomass_status": res.get("biomass_status", ""),
            "pfba_total_flux": res.get("pfba_total_flux", np.nan),
            "secondary_obj": res.get("secondary_obj", np.nan),
            "group": "High" if clone in high_clones else ("Low" if clone in low_clones else "Mid"),
        }
        explanation_rows.append(row)
        explanation_results[clone] = {**row, "fluxes": res.get("fluxes", {}), "ex_flux": res.get("ex_flux", {}), "constraints": cst}
else:
    scale = np.nan

explanation_df = pd.DataFrame(explanation_rows)

# Lactate/glucose apparent ratio from measured constraints.
G, L = "EX_glc_e", "EX_lac_L_e"
lac_glc_rows = []
for clone, cst in all_constraints.items():
    glc_uptake = abs(min(cst.get(G, (0, 0))))
    lac_lb, lac_ub = cst.get(L, (0, 0))
    lac_rate = lac_ub if lac_ub > 0 else lac_lb
    ratio = lac_rate / glc_uptake if glc_uptake > 0 else 0.0
    lac_glc_rows.append({"clone": clone, "lac_glc_ratio": ratio, "igG_day14": igG_day14.get(clone, 0)})
lac_glc_df = pd.DataFrame(lac_glc_rows).sort_values("igG_day14", ascending=False)

key_ex = [G, L, "EX_gln_L_e", "EX_nh4_e", "EX_ala_L_e", "EX_glu_L_e"]
primary_results = explanation_results if args.analysis_mode in ["explain", "both"] else prediction_results
ex_compare = pd.DataFrame({clone: {ex: res.get("ex_flux", {}).get(ex, 0) for ex in key_ex} for clone, res in primary_results.items()}).T
ex_compare.index.name = "clone"

# Combined output table.
combined = pd.concat([d for d in [prediction_df, explanation_df] if len(d) > 0], ignore_index=True) if (len(prediction_df) or len(explanation_df)) else pd.DataFrame()
if len(combined) > 0:
    combined = combined.sort_values(["analysis_type", "igG_day14"], ascending=[True, False])

# Save for downstream.
data["analysis_mode"] = args.analysis_mode
data["selected_objective"] = DEMAND_RXN
# Compatibility: downstream Escher/FVA defaults to explanation when available, otherwise prediction.
data["objective"] = DEMAND_RXN
if args.analysis_mode in ["explain", "both"]:
    data["objective_mode"] = "fixed_scaled_measured_qIgG_pfba_explanation_not_prediction"
    data["fba_results"] = explanation_results
    data["fba_df"] = explanation_df
else:
    data["objective_mode"] = "prediction_capacity_no_igg_constraint"
    data["fba_results"] = prediction_results
    data["fba_df"] = prediction_df

data["prediction_fba_results"] = prediction_results
data["prediction_fba_df"] = prediction_df
data["explanation_fba_results"] = explanation_results
data["explanation_fba_df"] = explanation_df
data["model_predictability_diagnostics"] = diagnostics_df
data["production_capacity_diagnostics"] = capacity_df
data["biomass_fraction"] = args.biomass_fraction
data["demand_scale"] = scale
data["demand_scale_method"] = str(args.demand_scale)
data["igg_interval_rates_mmol_gDCWh"] = igg_rates
data["lac_glc_df"] = lac_glc_df
data["ex_compare"] = ex_compare
with open(pkl_path, "wb") as f:
    pickle.dump(data, f)

if len(prediction_df): save_table(prediction_df, "prediction_fba_results.csv", DATASET)
if len(explanation_df): save_table(explanation_df, "explanation_fba_results.csv", DATASET)
if len(combined): save_table(combined, "fba_results.csv", DATASET)
if len(diagnostics_df): save_table(diagnostics_df, "model_predictability_diagnostics.csv", DATASET)
if len(capacity_df): save_table(capacity_df, "production_capacity_diagnostics.csv", DATASET)
save_table(lac_glc_df, "lac_glc_ratio.csv", DATASET)
save_table(ex_compare.reset_index(), "fba_exchange_compare.csv", DATASET)

print("\n  OK  03_run_fba 완료")
print(f"  다음: python scripts/steps/04_run_fva.py --dataset {DATASET} --analysis_mode {args.analysis_mode}")
