"""
21_full_internal_fva.py — workstation-scale internal/full FVA

Purpose
-------
The focused FVA steps (14/18) are fast and presentation-friendly, but they only
cover curated central/mAb pathway panels. This step runs a larger FVA scope for
Linux/workstation use:

  --fva_scope internal : all non-exchange internal reactions
  --fva_scope all      : all model reactions, including exchange/boundary
  --fva_scope exchange : exchange only, for QC/reference
  --fva_scope focused  : focused list if central_mab_reaction_panel.csv exists

Recommended workstation command:
  python run_pipeline.py --dataset practice_20aa --input_format tsv \
    --steps 17,1,2,3,10,14,16,21,19,20,6 \
    --fva_scope internal --fva_targets group_avg --fva_processes 16

Outputs
-------
results/<dataset>/tables/full_fva/
  full_fva_<scope>_<mode>_<condition>.csv
  full_fva_<scope>_summary.csv
  full_fva_<scope>_combined_report.csv
  full_fva_<scope>_high_low_overlap.csv
  full_fva_<scope>_central_mab_subset.csv
results/<dataset>/figures/
  Fig21_full_fva_high_low_separation.png
  Fig21B_full_fva_range_delta.png
"""
import argparse
import os
import pickle
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cobra.flux_analysis import flux_variability_analysis


def _find_root():
    cur = os.path.dirname(os.path.abspath(__file__))
    for _ in range(7):
        if os.path.exists(os.path.join(cur, "src", "config.py")):
            return cur
        cur = os.path.dirname(cur)
    return cur

ROOT = _find_root()
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from src.config import *  # noqa
from src.fba_utils import (
    load_model, results_dir, save_table, save_figure, apply_bounds,
    avg_constraints, choose_objective_from_data, estimate_igg_interval_rates,
    compute_measured_demand_scale, set_biomass_minimum, get_fva_reaction_list,
)

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa", choices=["sowa2020", "practice_20aa", "own_experiment"])
parser.add_argument("--objective", default=os.environ.get("CHO_OBJECTIVE"), help="Production objective; default DM_igg_g")
parser.add_argument("--analysis_mode", default=os.environ.get("CHO_ANALYSIS_MODE", "both"), choices=["predict", "explain", "both"])
parser.add_argument("--fva_scope", default=os.environ.get("CHO_FVA_SCOPE", "internal"), choices=["exchange", "focused", "internal", "all"], help="FVA reaction scope")
parser.add_argument("--fva_targets", default=os.environ.get("CHO_FVA_TARGETS", "group_avg"), choices=["group_avg", "representative", "all_clones"], help="Which conditions to run for expensive full/internal FVA")
parser.add_argument("--fraction", type=float, default=FVA_FRACTION, help="FVA fraction of optimum")
parser.add_argument("--fva_processes", type=int, default=int(os.environ.get("CHO_FVA_PROCESSES", FVA_PROCESSES)), help="Parallel FVA processes; use >1 on Linux workstation")
parser.add_argument("--biomass_fraction", type=float, default=None)
parser.add_argument("--demand_scale", default=os.environ.get("CHO_DEMAND_SCALE", MEASURED_IGG_DEMAND_SCALE))
parser.add_argument("--resume", action="store_true", help="Skip condition/mode CSVs that already exist")
args = parser.parse_args()
DATASET = args.dataset
if args.biomass_fraction is None:
    args.biomass_fraction = float(os.environ.get("CHO_BIOMASS_FRACTION", BIOMASS_MIN_FRACTION))

print("=" * 72)
print(f"  21_full_internal_fva.py — {DATASET}")
print("=" * 72)
print(f"  FVA scope     : {args.fva_scope}")
print(f"  FVA targets   : {args.fva_targets}")
print(f"  FVA fraction  : {args.fraction}")
print(f"  FVA processes : {args.fva_processes}")
print("  Note: internal/all FVA may take a long time. Run on Linux/workstation.")

pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
if not os.path.exists(pkl_path):
    print(f"  !! {pkl_path} 없음 — steps 1,2 먼저 실행")
    sys.exit(1)
with open(pkl_path, "rb") as f:
    data = pickle.load(f)

all_constraints = data.get("all_constraints", {})
clones = data.get("clones", [])
high_clones = data.get("high_clones", [])
low_clones = data.get("low_clones", [])
clone_groups = data.get("clone_groups", {}) or {}
if not all_constraints or not clones:
    print("  !! constraints/clones 없음")
    sys.exit(1)

model = load_model(verbose=True)
OBJ = choose_objective_from_data(model, data, args.objective)
rxn_ids = {r.id for r in model.reactions}
DEMAND_RXN = "DM_igg_g" if "DM_igg_g" in rxn_ids else OBJ
igg_rates = estimate_igg_interval_rates(data)

# Reaction list for scope. For focused scope, prefer central panel from step 14.
focused_ids = None
panel_path = os.path.join(results_dir(DATASET, "tables"), "central_mab_reaction_panel.csv")
if args.fva_scope == "focused" and os.path.exists(panel_path):
    focused_ids = pd.read_csv(panel_path)["rxn_id"].dropna().astype(str).tolist()
rxn_list = get_fva_reaction_list(model, scope=args.fva_scope, focused_ids=focused_ids)
print(f"  Reaction list : {len(rxn_list)} reactions")
if args.fva_scope in {"internal", "all"} and len(rxn_list) < 1000:
    print("  ⚠ Internal/all FVA reaction count looks unexpectedly small; check model/scope.")

# Build target conditions.
conditions = []
seen = set()
def add_condition(name, clist_or_name, cst, qobs):
    if name in seen or not cst:
        return
    conditions.append((name, clist_or_name, cst, float(qobs or 0.0)))
    seen.add(name)

group_to_clones = {}
for c in clones:
    g = str(clone_groups.get(c, "Mid") or "Mid")
    group_to_clones.setdefault(g, []).append(c)

if args.fva_targets in ["group_avg", "all_clones"]:
    for gname in ["High", "Mother", "Low", "Mid"] + sorted([g for g in group_to_clones if g not in {"High", "Mother", "Low", "Mid"}]):
        clist = group_to_clones.get(gname, [])
        if clist:
            qobs = float(np.mean([igg_rates.get(c, 0.0) for c in clist if c in igg_rates]) or 0.0)
            add_condition(f"{gname}Avg", ",".join(clist), avg_constraints(clist, all_constraints), qobs)
    if high_clones and "HighAvg" not in seen:
        add_condition("HighAvg", ",".join(high_clones), avg_constraints(high_clones, all_constraints), float(np.mean([igg_rates.get(c, 0.0) for c in high_clones if c in igg_rates]) or 0.0))
    if low_clones and "LowAvg" not in seen:
        add_condition("LowAvg", ",".join(low_clones), avg_constraints(low_clones, all_constraints), float(np.mean([igg_rates.get(c, 0.0) for c in low_clones if c in igg_rates]) or 0.0))

if args.fva_targets == "representative":
    # Use top high, median/mother, and bottom low representative clones.
    sorted_clones = data.get("sorted_clones", clones)
    reps = []
    if high_clones:
        reps.append(high_clones[0])
    mother = [c for c in clones if str(clone_groups.get(c, "")).lower().startswith("mother")]
    if mother:
        reps.append(mother[len(mother)//2])
    elif sorted_clones:
        reps.append(sorted_clones[len(sorted_clones)//2])
    if low_clones:
        reps.append(low_clones[-1])
    for c in dict.fromkeys(reps):
        add_condition(c, c, all_constraints.get(c, {}), igg_rates.get(c, 0.0))

if args.fva_targets == "all_clones":
    for c in clones:
        add_condition(c, c, all_constraints.get(c, {}), igg_rates.get(c, 0.0))

if not conditions:
    print("  !! No FVA conditions could be built")
    sys.exit(1)
print("  Conditions    : " + ", ".join([c[0] for c in conditions]))

modes = []
if args.analysis_mode in ["predict", "both"]:
    modes.append("no_igg_input")
if args.analysis_mode in ["explain", "both"]:
    modes.append("measured_demand")
if not modes:
    modes = ["no_igg_input"]

# Compute global demand scaling for measured-demand mode.
strict_caps = {}
for cond_name, source, cst, qobs in conditions:
    with model:
        apply_bounds(model, cst)
        set_biomass_minimum(model, cst, args.biomass_fraction)
        model.objective = OBJ
        sol = model.optimize()
        strict_caps[cond_name] = float(sol.objective_value) if sol.status == "optimal" and sol.objective_value is not None else 0.0
scale = compute_measured_demand_scale(
    {cond_name: qobs for cond_name, source, cst, qobs in conditions},
    strict_caps,
    requested=args.demand_scale,
    safety=MEASURED_IGG_DEMAND_SCALE_SAFETY,
)

out_dir = os.path.join(results_dir(DATASET, "tables"), "full_fva")
os.makedirs(out_dir, exist_ok=True)
fig_dir = results_dir(DATASET, "figures")
summary_rows = []
combined_paths = []

# Reaction metadata, used for annotation.
ex_ids = {r.id for r in model.exchanges}
meta = []
for r in rxn_list:
    meta.append({
        "rxn_id": r.id,
        "reaction_name": r.name,
        "subsystem": getattr(r, "subsystem", "") or "",
        "is_exchange": r.id in ex_ids,
        "is_boundary": bool(getattr(r, "boundary", False)),
        "reaction": r.reaction,
    })
meta_df = pd.DataFrame(meta).drop_duplicates("rxn_id")
meta_df.to_csv(os.path.join(out_dir, f"full_fva_{args.fva_scope}_reaction_metadata.csv"), index=False)

for mode in modes:
    for cond_name, source, cst, qobs in conditions:
        safe_cond = str(cond_name).replace("/", "_").replace(" ", "_")
        path = os.path.join(out_dir, f"full_fva_{args.fva_scope}_{mode}_{safe_cond}.csv")
        if args.resume and os.path.exists(path):
            print(f"  [resume] {mode:15s} {cond_name:12s} existing")
            combined_paths.append(path)
            continue
        with model:
            n = apply_bounds(model, cst)
            set_biomass_minimum(model, cst, args.biomass_fraction)
            target = 0.0
            if mode == "measured_demand":
                target = max(0.0, float(qobs or 0.0) * scale)
                if DEMAND_RXN in rxn_ids:
                    dr = model.reactions.get_by_id(DEMAND_RXN)
                    dr.lower_bound = max(float(dr.lower_bound), target)
                    dr.upper_bound = max(float(dr.lower_bound), target)
            model.objective = OBJ if mode == "no_igg_input" else DEMAND_RXN
            sol = model.optimize()
            if sol.status != "optimal":
                print(f"  !! {mode}/{cond_name}: infeasible ({sol.status})")
                summary_rows.append({"mode": mode, "condition": cond_name, "status": sol.status, "n_reactions": 0})
                continue
            obj_val = float(sol.objective_value) if sol.objective_value is not None else 0.0
            print(f"  FVA start {mode:15s} {cond_name:12s} rxns={len(rxn_list)} obj={obj_val:.4g}")
            try:
                fva = flux_variability_analysis(
                    model,
                    reaction_list=rxn_list,
                    fraction_of_optimum=args.fraction,
                    processes=args.fva_processes,
                )
                fva = fva.reset_index().rename(columns={"index": "rxn_id"})
                if "reaction" in fva.columns and "rxn_id" not in fva.columns:
                    fva = fva.rename(columns={"reaction": "rxn_id"})
                fva["mean"] = (fva["minimum"] + fva["maximum"]) / 2.0
                fva["range"] = fva["maximum"] - fva["minimum"]
                fva["mode"] = mode
                fva["condition"] = cond_name
                fva["source_clones"] = source
                fva["fva_scope"] = args.fva_scope
                fva["objective"] = OBJ
                fva["objective_value"] = obj_val
                fva["scaled_demand"] = target
                fva["n_constraints"] = n
                fva = fva.merge(meta_df, on="rxn_id", how="left")
                fva.to_csv(path, index=False)
                combined_paths.append(path)
                summary_rows.append({
                    "mode": mode,
                    "condition": cond_name,
                    "status": "optimal",
                    "n_reactions": len(fva),
                    "objective_value": obj_val,
                    "scaled_demand": target,
                    "median_range": float(fva["range"].median()),
                    "mean_range": float(fva["range"].mean()),
                    "n_nonzero_range": int((fva["range"].abs() > 1e-12).sum()),
                    "output_csv": os.path.relpath(path, ROOT),
                })
                print(f"  [saved] {os.path.relpath(path, ROOT)}")
            except Exception as e:
                print(f"  !! FVA failed {mode}/{cond_name}: {e}")
                summary_rows.append({"mode": mode, "condition": cond_name, "status": f"failed: {e}", "n_reactions": 0})

summary = pd.DataFrame(summary_rows)
summary_path = os.path.join(out_dir, f"full_fva_{args.fva_scope}_summary.csv")
summary.to_csv(summary_path, index=False)
print(f"  [saved] {os.path.relpath(summary_path, ROOT)}")

# Combine manageable CSVs into one report. This is useful for downstream filters,
# but the per-condition CSV files are the true checkpoint outputs.
frames = []
for path in combined_paths:
    try:
        frames.append(pd.read_csv(path))
    except Exception:
        pass
combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
if not combined.empty:
    combined_path = os.path.join(out_dir, f"full_fva_{args.fva_scope}_combined_report.csv")
    combined.to_csv(combined_path, index=False)
    print(f"  [saved] {os.path.relpath(combined_path, ROOT)}")

# High vs Low FVA overlap/separation. Prefer measured_demand if present because
# it explains observed production phenotype. Otherwise use no_igg_input.
def _pick_mode(df):
    modes_available = set(df["mode"].dropna().unique()) if "mode" in df else set()
    return "measured_demand" if "measured_demand" in modes_available else ("no_igg_input" if "no_igg_input" in modes_available else None)

if not combined.empty and {"HighAvg", "LowAvg"}.issubset(set(combined["condition"])):
    comp_mode = _pick_mode(combined)
    sub = combined[combined["mode"] == comp_mode].copy()
    h = sub[sub["condition"] == "HighAvg"].set_index("rxn_id")
    l = sub[sub["condition"] == "LowAvg"].set_index("rxn_id")
    common = h.index.intersection(l.index)
    rows = []
    for rid in common:
        hmin, hmax = float(h.loc[rid, "minimum"]), float(h.loc[rid, "maximum"])
        lmin, lmax = float(l.loc[rid, "minimum"]), float(l.loc[rid, "maximum"])
        overlap = max(0.0, min(hmax, lmax) - max(hmin, lmin))
        h_range = hmax - hmin
        l_range = lmax - lmin
        union = max(hmax, lmax) - min(hmin, lmin)
        separation = max(0.0, max(hmin, lmin) - min(hmax, lmax))
        rows.append({
            "mode": comp_mode,
            "rxn_id": rid,
            "reaction_name": h.loc[rid].get("reaction_name", ""),
            "subsystem": h.loc[rid].get("subsystem", ""),
            "is_exchange": h.loc[rid].get("is_exchange", False),
            "high_min": hmin,
            "high_max": hmax,
            "high_range": h_range,
            "low_min": lmin,
            "low_max": lmax,
            "low_range": l_range,
            "range_delta_high_minus_low": h_range - l_range,
            "mean_delta_high_minus_low": ((hmin + hmax) / 2.0) - ((lmin + lmax) / 2.0),
            "overlap_width": overlap,
            "union_width": union,
            "overlap_fraction": overlap / union if union > 1e-12 else 1.0,
            "separation_width": separation,
            "reaction": h.loc[rid].get("reaction", ""),
        })
    ov = pd.DataFrame(rows)
    ov["abs_mean_delta"] = ov["mean_delta_high_minus_low"].abs()
    ov["abs_range_delta"] = ov["range_delta_high_minus_low"].abs()
    ov = ov.sort_values(["separation_width", "abs_mean_delta", "abs_range_delta"], ascending=[False, False, False])
    ov_path = os.path.join(out_dir, f"full_fva_{args.fva_scope}_high_low_overlap.csv")
    ov.to_csv(ov_path, index=False)
    print(f"  [saved] {os.path.relpath(ov_path, ROOT)}")

    # Central/mAb subset if panel exists.
    if os.path.exists(panel_path):
        panel = pd.read_csv(panel_path)
        subset = ov[ov["rxn_id"].isin(panel["rxn_id"].astype(str))].merge(panel[["rxn_id", "pathway", "label"]], on="rxn_id", how="left")
        subset_path = os.path.join(out_dir, f"full_fva_{args.fva_scope}_central_mab_subset.csv")
        subset.to_csv(subset_path, index=False)
        print(f"  [saved] {os.path.relpath(subset_path, ROOT)}")

    # Figures
    try:
        top = ov.sort_values(["separation_width", "abs_mean_delta"], ascending=[False, False]).head(30).iloc[::-1]
        if not top.empty:
            labels = top["rxn_id"] + " | " + top["reaction_name"].fillna("").str.slice(0, 35)
            fig, ax = plt.subplots(figsize=(11, max(6, 0.28 * len(top) + 1.5)))
            ax.barh(labels, top["separation_width"])
            ax.set_xlabel("High/Low FVA non-overlap width")
            ax.set_title(f"Figure 21. Genome-scale/internal FVA High-vs-Low separation ({args.fva_scope}, {comp_mode})", weight="bold")
            ax.grid(True, axis="x", ls=":", alpha=0.5)
            fig.tight_layout()
            fig.savefig(os.path.join(fig_dir, "Fig21_full_fva_high_low_separation.png"), dpi=250, bbox_inches="tight")
            plt.close(fig)
        top2 = ov.sort_values("abs_range_delta", ascending=False).head(30).iloc[::-1]
        if not top2.empty:
            labels = top2["rxn_id"] + " | " + top2["reaction_name"].fillna("").str.slice(0, 35)
            fig, ax = plt.subplots(figsize=(11, max(6, 0.28 * len(top2) + 1.5)))
            ax.barh(labels, top2["range_delta_high_minus_low"])
            ax.axvline(0, color="black", lw=1)
            ax.set_xlabel("HighAvg FVA range - LowAvg FVA range")
            ax.set_title(f"Figure 21B. Genome-scale/internal FVA flexibility delta ({args.fva_scope}, {comp_mode})", weight="bold")
            ax.grid(True, axis="x", ls=":", alpha=0.5)
            fig.tight_layout()
            fig.savefig(os.path.join(fig_dir, "Fig21B_full_fva_range_delta.png"), dpi=250, bbox_inches="tight")
            plt.close(fig)
            print("  [saved] results/%s/figures/Fig21*.png" % DATASET)
    except Exception as e:
        print(f"  !! Fig21 skipped: {e}")

readme = f"""# Full/Internal FVA report

Scope: `{args.fva_scope}`  
Targets: `{args.fva_targets}`  
FVA fraction: `{args.fraction}`  
Processes: `{args.fva_processes}`

This step is the workstation-scale FVA layer. It is different from the focused
central/mAb FVA in step 14:

- `focused`: curated central metabolism / mAb panel.
- `internal`: all non-exchange internal reactions.
- `all`: every reaction in the model.

For scientific reporting, say **focused FVA** when using step 14 only, and say
**internal/genome-scale FVA** only when this step was run with `internal` or `all`.

Recommended interpretation:

1. Use `full_fva_{args.fva_scope}_summary.csv` to confirm which conditions ran.
2. Use `full_fva_{args.fva_scope}_high_low_overlap.csv` to find reactions whose
   feasible ranges differ between HighAvg and LowAvg.
3. Use `full_fva_{args.fva_scope}_central_mab_subset.csv` to connect the genome-scale
   FVA back to the focused pathway panel.
4. Very large FVA ranges may reflect under-constrained reactions or thermodynamic
   loops; prioritize focused pathway reactions and robust high/low separation.
"""
with open(os.path.join(out_dir, "README_FULL_INTERNAL_FVA.md"), "w", encoding="utf-8") as f:
    f.write(readme)
print("  OK full/internal FVA complete")
