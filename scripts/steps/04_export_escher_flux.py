"""
export_escher_flux.py — Escher-compatible flux JSON export

Corrected behavior
------------------
Fluxes are exported from FBA solutions that maximize the fixed IgG/mAb
production objective (default DM_igg_g). The script writes both original iCHO3K
reaction IDs and RECON/BiGG-friendly alias IDs.
"""
import argparse
import json
import os
import pickle
import sys
import warnings
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
from src.fba_utils import load_model, apply_bounds, avg_constraints, choose_objective_from_data, set_biomass_minimum, estimate_igg_interval_rates
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa")
parser.add_argument("--objective", default=os.environ.get("CHO_OBJECTIVE"), help="Production objective reaction id; default DM_igg_g")
parser.add_argument("--biomass_fraction", type=float, default=None, help="Optional biomass lower-bound fraction")
parser.add_argument(
    "--policy",
    default="auto",
    choices=["strict", "visualization", "none", "auto"],
    help="strict=use measured bounds; visualization=relaxed flux-capacity bounds; none=no exchange constraints; auto=strict then visualization.",
)
parser.add_argument("--tol", type=float, default=1e-9)
parser.add_argument("--analysis_mode", default=os.environ.get("CHO_ANALYSIS_MODE", "both"), choices=["predict", "explain", "both"], help="predict exports max-objective flux; explain/both exports fixed measured-demand flux for visualization.")
args = parser.parse_args()
DATASET = args.dataset
if args.biomass_fraction is None:
    args.biomass_fraction = float(os.environ.get("CHO_BIOMASS_FRACTION", BIOMASS_MIN_FRACTION))
TOL = args.tol

print("=" * 65)
print(f"  export_escher_flux.py — {DATASET}")
print("=" * 65)

pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
if not os.path.exists(pkl_path):
    raise FileNotFoundError(f"{pkl_path} 없음 — 01/02/03 먼저 실행")
with open(pkl_path, "rb") as f:
    data = pickle.load(f)

all_constraints = data.get("all_constraints", {})
clones = data.get("clones", [])
clone_groups = data.get("clone_groups", {}) or {}
high_clones = data.get("high_clones", [])
low_clones = data.get("low_clones", [])
igg_rates = data.get("igg_interval_rates_mmol_gDCWh") or estimate_igg_interval_rates(data)
analysis_mode = args.analysis_mode or data.get("analysis_mode", "both")
try:
    demand_scale = float(data.get("demand_scale", 1.0))
    if not np.isfinite(demand_scale):
        demand_scale = 0.0
except Exception:
    demand_scale = 0.0
if not all_constraints or not clones:
    raise ValueError("Missing constraints or clone information")

model = load_model(verbose=True)
OBJ = choose_objective_from_data(model, data, args.objective)
print(f"\n  Objective : {OBJ} (IgG/mAb production)")
print(f"  Biomass minimum fraction: {args.biomass_fraction}")
print(f"  Policy    : {args.policy}")
print(f"  Analysis mode for export: {analysis_mode}")

out_dir = results_dir(DATASET, "tables")

ICHOK_TO_RECON = {
    # Exchange / BiGG naming compatibility
    "EX_glc_e": "EX_glc__D_e",
    "EX_lac_L_e": "EX_lac__L_e",
    "EX_gln_L_e": "EX_gln__L_e",
    "EX_glu_L_e": "EX_glu__L_e",
    "EX_ala_L_e": "EX_ala__L_e",
    "EX_arg_L_e": "EX_arg__L_e",
    "EX_asn_L_e": "EX_asn__L_e",
    "EX_asp_L_e": "EX_asp__L_e",
    "EX_cys_L_e": "EX_cys__L_e",
    "EX_his_L_e": "EX_his__L_e",
    "EX_ile_L_e": "EX_ile__L_e",
    "EX_leu_L_e": "EX_leu__L_e",
    "EX_lys_L_e": "EX_lys__L_e",
    "EX_met_L_e": "EX_met__L_e",
    "EX_phe_L_e": "EX_phe__L_e",
    "EX_pro_L_e": "EX_pro__L_e",
    "EX_ser_L_e": "EX_ser__L_e",
    "EX_thr_L_e": "EX_thr__L_e",
    "EX_trp_L_e": "EX_trp__L_e",
    "EX_tyr_L_e": "EX_tyr__L_e",
    "EX_val_L_e": "EX_val__L_e",
    # Core central metabolism aliases where official RECON1 Escher maps use one combined id.
    "ACONTam": "ACONTm",
    "ACONTbm": "ACONTm",
}


def to_recon_id(rid):
    return ICHOK_TO_RECON.get(rid, rid)


def apply_visualization_bounds(model, constraints):
    """Relaxed display mode for Escher only; not for quantitative interpretation."""
    rxn_ids = {r.id for r in model.reactions}
    n = 0
    for key, (lb, ub) in (constraints or {}).items():
        rid = EXCHANGE_IDS.get(key, key)
        if rid not in rxn_ids:
            continue
        rxn = model.reactions.get_by_id(rid)
        lo, hi = min(float(lb), float(ub)), max(float(lb), float(ub))
        if hi <= 0:  # observed uptake
            rxn.lower_bound = min(rxn.lower_bound, lo)
            rxn.upper_bound = max(rxn.upper_bound, 1000.0)
        elif lo >= 0:  # observed secretion
            rxn.lower_bound = min(rxn.lower_bound, -1000.0)
            rxn.upper_bound = max(rxn.upper_bound, hi)
        else:
            rxn.lower_bound = lo
            rxn.upper_bound = hi
        n += 1
    return n


def solve_group(label, clone_list, policy):
    cst = avg_constraints(clone_list, all_constraints)
    qigg = float(np.mean([igg_rates.get(c, 0.0) for c in clone_list if c in igg_rates]) or 0.0)
    with model:
        if policy == "strict":
            n = apply_bounds(model, cst)
        elif policy == "visualization":
            n = apply_visualization_bounds(model, cst)
        elif policy == "none":
            n = 0
        else:
            raise ValueError(policy)
        bio_cap, bio_lb, bio_status = set_biomass_minimum(model, cst if policy != "none" else {}, args.biomass_fraction)
        model.objective = OBJ
        sol = model.optimize()
        mode = "strict_predicted_capacity"
        if sol.status == "optimal":
            obj = float(sol.objective_value) if sol.objective_value is not None else 0.0
        else:
            obj = 0.0

        # If true IgG objective is zero, do not export arbitrary zero-objective fluxes.
        # Instead, force observed qIgG demand and optimize biomass to obtain a feasible
        # flux state for visualization. This is labeled measured-demand mode.
        if analysis_mode in ["explain", "both"] and demand_scale > 0:  # export fixed measured-demand flux state when explanation mode is available
            demand_rxn = "DM_igg_g" if "DM_igg_g" in {r.id for r in model.reactions} else OBJ
            if demand_rxn in {r.id for r in model.reactions}:
                dr = model.reactions.get_by_id(demand_rxn)
                target = max(0.0, qigg * demand_scale)
                dr.lower_bound = target
                dr.upper_bound = target
                model.objective = demand_rxn
                sol = model.optimize()
                mode = "fixed_scaled_measured_qIgG_pfba_for_visualization"
                if sol.status == "optimal":
                    obj = float(sol.fluxes.get(demand_rxn, 0.0))
                else:
                    obj = 0.0
        if sol.status == "optimal":
            fluxes = {str(k): float(v) for k, v in sol.fluxes.to_dict().items()}
            nonzero = int(np.sum(np.abs(sol.fluxes.values) > TOL))
        else:
            fluxes, obj, nonzero = {}, 0.0, 0
    return {
        "label": label,
        "policy": policy,
        "mode": mode,
        "status": sol.status,
        "objective": OBJ,
        "objective_value": obj,
        "observed_qigg": qigg,
        "nonzero_fluxes": nonzero,
        "n_constraints": n,
        "biomass_lb": bio_lb,
        "biomass_capacity": bio_cap,
        "n_clones": len(clone_list),
        "clones": ",".join(clone_list),
        "fluxes": fluxes,
    }


def build_targets():
    targets = []
    # individual clone exports for demo/small datasets
    if DATASET == "practice_20aa" or len(clones) <= 12:
        for c in clones:
            targets.append((c, [c]))
    # biological group averages when available
    group_to_clones = {}
    for c in clones:
        g = str(clone_groups.get(c, "Mid") or "Mid")
        group_to_clones.setdefault(g, []).append(c)
    for g in ["High", "Mother", "Low", "Mid"] + sorted([x for x in group_to_clones if x not in {"High","Mother","Low","Mid"}]):
        if g in group_to_clones:
            targets.append((f"{g}Avg", group_to_clones[g]))
    # fallback top/bottom averages
    if high_clones and not any(lbl == "HighAvg" for lbl, _ in targets):
        targets.append(("HighAvg", high_clones))
    if low_clones and not any(lbl == "LowAvg" for lbl, _ in targets):
        targets.append(("LowAvg", low_clones))
    # de-duplicate while preserving order
    out=[]; seen=set()
    for lbl, clist in targets:
        if lbl in seen or not clist:
            continue
        out.append((lbl, clist))
        seen.add(lbl)
    return out

TARGETS = build_targets()

def run_policy(policy):
    return {label: solve_group(label, clist, policy) for label, clist in TARGETS}

def acceptable(res):
    return all(v["status"] == "optimal" and v["nonzero_fluxes"] > 0 and abs(v.get("objective_value", 0.0)) > TOL for v in res.values())


if args.policy == "auto":
    results = None
    chosen_policy = None
    for pol in ["strict", "visualization", "none"]:
        res = run_policy(pol)
        for label, r in res.items():
            print(f"  {label:6s}/{pol:13s}: {r['status']:10s} {OBJ}={r['objective_value']:.6g} nonzero={r['nonzero_fluxes']} mode={r.get('mode','')} n_cst={r['n_constraints']}")
        if acceptable(res):
            results, chosen_policy = res, pol
            break
    if results is None:
        results, chosen_policy = res, pol
else:
    chosen_policy = args.policy
    results = run_policy(chosen_policy)
    for label, r in results.items():
        print(f"  {label:6s}/{chosen_policy:13s}: {r['status']:10s} {OBJ}={r['objective_value']:.6g} nonzero={r['nonzero_fluxes']} mode={r.get('mode','')} n_cst={r['n_constraints']}")


def dict_cho(fluxes):
    return {rid: round(float(v), 6) for rid, v in fluxes.items()}


def dict_recon(fluxes):
    out = {}
    for rid, val in fluxes.items():
        tid = to_recon_id(rid)
        out[tid] = out.get(tid, 0.0) + float(val)
    return {rid: round(float(v), 6) for rid, v in out.items()}


def save_pair(label, data_dict, suffix):
    paths = []
    full = os.path.join(out_dir, f"escher_flux_{label}{suffix}_full.json")
    with open(full, "w", encoding="utf-8") as f:
        json.dump(data_dict, f, indent=2, ensure_ascii=False, sort_keys=True)
    paths.append(full)
    nz = {k: v for k, v in data_dict.items() if abs(v) > TOL}
    nonzero = os.path.join(out_dir, f"escher_flux_{label}{suffix}_nonzero.json")
    with open(nonzero, "w", encoding="utf-8") as f:
        json.dump(nz, f, indent=2, ensure_ascii=False, sort_keys=True)
    paths.append(nonzero)
    print(f"  [saved] {os.path.relpath(full, ROOT)} ({len(data_dict)} reactions)")
    print(f"  [saved] {os.path.relpath(nonzero, ROOT)} ({len(nz)} nonzero reactions)")
    return paths


saved = []
for label, r in results.items():
    flux = r["fluxes"]
    saved += save_pair(label, dict_recon(flux), "")
    saved += save_pair(label, dict_cho(flux), "_cho")

# group delta exports for common comparisons
def has_label(x):
    return x in results

comparison_pairs = []
if has_label("HighAvg") and has_label("LowAvg"):
    comparison_pairs.append(("HighAvg", "LowAvg", "high_minus_low"))
if has_label("HighAvg") and has_label("MotherAvg"):
    comparison_pairs.append(("HighAvg", "MotherAvg", "high_minus_mother"))
if has_label("MotherAvg") and has_label("LowAvg"):
    comparison_pairs.append(("MotherAvg", "LowAvg", "mother_minus_low"))

def save_delta(high_dict, low_dict, filename):
    ids = set(high_dict) | set(low_dict)
    delta = {
        rid: round(float(high_dict.get(rid, 0.0) - low_dict.get(rid, 0.0)), 6)
        for rid in ids
        if abs(high_dict.get(rid, 0.0)) > TOL or abs(low_dict.get(rid, 0.0)) > TOL
    }
    path = os.path.join(out_dir, filename)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(delta, f, indent=2, ensure_ascii=False, sort_keys=True)
    print(f"  [saved] {os.path.relpath(path, ROOT)} ({len(delta)} reactions)")
    saved.append(path)

for hi, lo, stem in comparison_pairs:
    save_delta(dict_recon(results[hi]["fluxes"]), dict_recon(results[lo]["fluxes"]), f"escher_flux_{stem}.json")
    save_delta(dict_cho(results[hi]["fluxes"]), dict_cho(results[lo]["fluxes"]), f"escher_flux_{stem}_cho.json")

status_rows = []
for label, r in results.items():
    row = {k: v for k, v in r.items() if k != "fluxes"}
    row["chosen_policy"] = chosen_policy
    status_rows.append(row)
status_df = pd.DataFrame(status_rows)
status_path = os.path.join(out_dir, "escher_flux_status.csv")
status_df.to_csv(status_path, index=False)
print(f"  [saved] {os.path.relpath(status_path, ROOT)}")

guide = f"""
Escher Flux Map Guide
=====================
Dataset      : {DATASET}
Objective    : {OBJ} (IgG/mAb production objective / diagnostic)
Policy       : {chosen_policy}
Demand scale : {demand_scale:.6g}
Biomass frac : {args.biomass_fraction}
Targets      : {', '.join(results.keys())}

Recommended workflow
--------------------
1) Use group-average maps first: `HighAvg`, `MotherAvg`, `LowAvg` when available.
2) Then inspect clone-specific maps to understand within-group diversity.
3) Use `*_minus_*` delta JSON files on the same Escher map to highlight pathway changes.

Useful files for CHO/custom maps (iCHO3K IDs):
- `escher_flux_HighAvg_cho_nonzero.json`, `escher_flux_MotherAvg_cho_nonzero.json`, `escher_flux_LowAvg_cho_nonzero.json`
- clone-level files such as `escher_flux_High01_cho_nonzero.json`
- delta files such as `escher_flux_high_minus_low_cho.json`, `escher_flux_high_minus_mother_cho.json`

Useful files for RECON/Escher maps:
- `recon1_pretty_map/` outputs from step 12
- corresponding `escher_flux_*_recon1_all_maps.json` exports if generated

Interpretation notes
--------------------
- Use these maps to compare central carbon metabolism, lactate handling, PPP, TCA, glutamine/nitrogen metabolism, energy metabolism, and mAb synthesis reactions.
- If mode indicates measured-demand visualization, treat the map as an explanation of an observed production phenotype, not a titer prediction.
- The most persuasive demonstration for teams is usually `HighAvg` vs `LowAvg` plus one representative clone from each group.
"""
with open(os.path.join(out_dir, "escher_guide.txt"), "w", encoding="utf-8") as f:
    f.write(guide)
print(f"  [saved] results/{DATASET}/tables/escher_guide.txt")

# Update pkl metadata.
data["selected_objective"] = OBJ
data["objective"] = OBJ
data["objective_mode"] = "zero_aware_true_igg_objective"
data["biomass_fraction"] = args.biomass_fraction
data["escher_flux_status"] = status_df
with open(pkl_path, "wb") as f:
    pickle.dump(data, f)

print("\n  완료: Escher JSON export")
