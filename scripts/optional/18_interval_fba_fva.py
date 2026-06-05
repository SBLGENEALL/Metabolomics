"""
18_interval_fba_fva.py — Interval-wise focused pFBA/FVA

Runs the existing step 1 rate calculation logic for every consecutive day pair,
then performs focused pFBA/FVA on the central+mAb panel for HighAvg/MotherAvg/LowAvg.

This is useful for fed-batch interpretation because a single Day 7→10 interval
can miss lactate shifts and production-phase metabolic transitions.
"""
import argparse
import os
import pickle
import subprocess
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd


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
from src.config import DATA_PROCESSED, results_dir, FVA_FRACTION, FVA_PROCESSES, BIOMASS_MIN_FRACTION, MEASURED_IGG_DEMAND_SCALE, MEASURED_IGG_DEMAND_SCALE_SAFETY  # noqa
from src.fba_utils import load_model, apply_bounds, avg_constraints, choose_objective_from_data, estimate_igg_interval_rates, compute_measured_demand_scale, set_biomass_minimum  # noqa
from cobra.flux_analysis import pfba, flux_variability_analysis

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa", choices=["sowa2020", "practice_20aa", "own_experiment"])
parser.add_argument("--analysis_mode", default=os.environ.get("CHO_ANALYSIS_MODE", "both"), choices=["predict", "explain", "both"])
parser.add_argument("--objective", default=os.environ.get("CHO_OBJECTIVE"))
parser.add_argument("--biomass_fraction", type=float, default=float(os.environ.get("CHO_BIOMASS_FRACTION", BIOMASS_MIN_FRACTION)))
parser.add_argument("--demand_scale", default=os.environ.get("CHO_DEMAND_SCALE", MEASURED_IGG_DEMAND_SCALE))
parser.add_argument("--feed_volume_mode", default=os.environ.get("CHO_FEED_VOLUME_MODE", "interval"), choices=["interval", "cumulative"])
parser.add_argument("--intervals", default="auto", help="auto or semicolon list such as 3,5;5,7;7,10")
parser.add_argument("--fraction", type=float, default=FVA_FRACTION)
args = parser.parse_args()
DATASET = args.dataset
print("="*65)
print(f"  18_interval_fba_fva.py — {DATASET}")
print("="*65)

# Load current processed data to discover available days and use current clone grouping.
pkl_current = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
if not os.path.exists(pkl_current):
    raise FileNotFoundError("Run step 1 first so interval script can discover raw data/grouping.")
with open(pkl_current, "rb") as f:
    current = pickle.load(f)
raw = current.get("df_raw")
if raw is None or raw.empty:
    raise ValueError("Processed data does not include df_raw.")

if args.intervals == "auto":
    days = sorted([int(x) for x in pd.to_numeric(raw["DAY"], errors="coerce").dropna().unique()])
    intervals = [(days[i], days[i+1]) for i in range(len(days)-1)]
else:
    intervals = []
    for part in args.intervals.split(";"):
        a, b = [int(x) for x in part.split(",")]
        intervals.append((a,b))

# Use step 1 as the canonical rate calculator for each interval.
step1 = os.path.join(ROOT, "scripts", "steps", "01_load_data.py")
original_bytes = None
if os.path.exists(pkl_current):
    with open(pkl_current, "rb") as f:
        original_bytes = f.read()

interval_data = {}
for d1, d2 in intervals:
    cmd = [sys.executable, step1, "--dataset", DATASET, "--rate_days", f"{d1},{d2}", "--feed_volume_mode", args.feed_volume_mode]
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, encoding="utf-8", errors="replace")
    if r.returncode != 0:
        print(f"  !! rate calculation failed for {d1}->{d2}: {r.stderr[-500:]}")
        continue
    with open(pkl_current, "rb") as f:
        interval_data[(d1,d2)] = pickle.load(f)

# Restore original rates_{dataset}.pkl so downstream single-interval steps are not surprised.
if original_bytes is not None:
    with open(pkl_current, "wb") as f:
        f.write(original_bytes)

# Resolve panel from central_mab_reaction_panel if it exists; otherwise a compact default.
model = load_model(verbose=True)
rxn_ids = {r.id for r in model.reactions}
OBJ = choose_objective_from_data(model, current, args.objective)
DEMAND_RXN = "DM_igg_g" if "DM_igg_g" in rxn_ids else OBJ
panel_file = os.path.join(results_dir(DATASET, "tables"), "central_mab_reaction_panel.csv")
if os.path.exists(panel_file):
    panel = pd.read_csv(panel_file)
    panel = panel[panel["rxn_id"].isin(rxn_ids)].copy()
else:
    default_ids = ["EX_glc_e","EX_lac_L_e","EX_gln_L_e","EX_nh4_e","HEX1","PGI","PFK","PYK","LDH_L","PCm","CSm","ICDHyrm","AKGDm","SUCD1m","G6PDH2r","GND","GLUDym","ATPM","ATPS4m","igg_formation","DM_igg_g"]
    rows=[]
    for rid in default_ids:
        if rid in rxn_ids:
            r=model.reactions.get_by_id(rid)
            rows.append({"pathway":"Focused", "label":r.name or rid, "rxn_id":rid})
    panel=pd.DataFrame(rows)
rxn_list = [model.reactions.get_by_id(rid) for rid in panel["rxn_id"].tolist() if rid in rxn_ids]

def _conditions(data):
    clones = data.get("clones", [])
    clone_groups = data.get("clone_groups", {}) or {}
    all_constraints = data.get("all_constraints", {})
    igg_rates = estimate_igg_interval_rates(data)
    group_to_clones = {}
    for c in clones:
        group_to_clones.setdefault(str(clone_groups.get(c, "Mid") or "Mid"), []).append(c)
    out=[]
    for g in ["High","Mother","Low","Mid"] + sorted([x for x in group_to_clones if x not in {"High","Mother","Low","Mid"}]):
        if g in group_to_clones:
            clist=group_to_clones[g]
            out.append((f"{g}Avg", avg_constraints(clist, all_constraints), float(np.mean([igg_rates.get(c,0.0) for c in clist]) or 0.0)))
    return out

flux_rows=[]
fva_rows=[]
for (d1,d2), data in interval_data.items():
    conds=_conditions(data)
    # strict caps for demand scaling in this interval
    caps={}
    for name,cst,qobs in conds:
        with model:
            apply_bounds(model,cst)
            set_biomass_minimum(model,cst,args.biomass_fraction)
            model.objective=OBJ
            sol=model.optimize()
            caps[name]=float(sol.objective_value) if sol.status=="optimal" and sol.objective_value is not None else 0.0
    scale=compute_measured_demand_scale({name:q for name,cst,q in conds}, caps, requested=args.demand_scale, safety=MEASURED_IGG_DEMAND_SCALE_SAFETY)
    modes=[]
    if args.analysis_mode in ["predict","both"]: modes.append("no_igg_input")
    if args.analysis_mode in ["explain","both"]: modes.append("measured_demand")
    for mode in modes:
        for name,cst,qobs in conds:
            with model:
                apply_bounds(model,cst)
                set_biomass_minimum(model,cst,args.biomass_fraction)
                target=0.0
                if mode=="measured_demand" and DEMAND_RXN in rxn_ids:
                    target=max(0.0,float(qobs or 0.0)*scale)
                    dr=model.reactions.get_by_id(DEMAND_RXN)
                    dr.lower_bound=max(float(dr.lower_bound),target)
                    dr.upper_bound=max(float(dr.lower_bound),target)
                model.objective=OBJ if mode=="no_igg_input" else DEMAND_RXN
                sol=model.optimize()
                if sol.status!="optimal":
                    print(f"  !! {d1}->{d2} {mode} {name}: {sol.status}")
                    continue
                try:
                    psol=pfba(model); fluxes=psol.fluxes
                except Exception:
                    fluxes=sol.fluxes
                for _,p in panel.iterrows():
                    rid=p["rxn_id"]
                    flux_rows.append({"interval":f"{d1}->{d2}","day_start":d1,"day_end":d2,"mode":mode,"condition":name,"pathway":p.get("pathway",""),"label":p.get("label",rid),"rxn_id":rid,"flux":float(fluxes.get(rid,0.0)),"abs_flux":abs(float(fluxes.get(rid,0.0))),"scaled_demand":target})
                try:
                    fva=flux_variability_analysis(model,reaction_list=rxn_list,fraction_of_optimum=args.fraction,processes=FVA_PROCESSES)
                    for rid,row in fva.iterrows():
                        p=panel[panel["rxn_id"]==rid].iloc[0]
                        mn,mx=float(row["minimum"]),float(row["maximum"])
                        fva_rows.append({"interval":f"{d1}->{d2}","day_start":d1,"day_end":d2,"mode":mode,"condition":name,"pathway":p.get("pathway",""),"label":p.get("label",rid),"rxn_id":rid,"minimum":mn,"maximum":mx,"mean":(mn+mx)/2,"range":mx-mn,"scaled_demand":target})
                except Exception as e:
                    print(f"  !! FVA failed {d1}->{d2} {mode} {name}: {e}")
    print(f"  OK interval {d1}->{d2} conditions={len(conds)}")

out=results_dir(DATASET,"tables")
flux_df=pd.DataFrame(flux_rows); fva_df=pd.DataFrame(fva_rows)
flux_df.to_csv(os.path.join(out,"interval_flux_state.csv"),index=False)
fva_df.to_csv(os.path.join(out,"interval_fva_report.csv"),index=False)
print(f"  [saved] results/{DATASET}/tables/interval_flux_state.csv")
print(f"  [saved] results/{DATASET}/tables/interval_fva_report.csv")

# Plot pathway scores across intervals if available.
try:
    import matplotlib.pyplot as plt
    if not flux_df.empty:
        sub=flux_df[flux_df["mode"]==( "measured_demand" if "measured_demand" in set(flux_df["mode"]) else flux_df["mode"].iloc[0])].copy()
        scores=sub.groupby(["interval","condition","pathway"],as_index=False)["abs_flux"].mean()
        keep=["Glycolysis","Pyruvate/Lactate","TCA","PPP/Redox","Glutamine/Nitrogen","Energy","mAb synthesis","Exchange"]
        scores=scores[scores["pathway"].isin(keep)]
        piv=scores.pivot_table(index="pathway",columns=["interval","condition"],values="abs_flux",aggfunc="first").fillna(0)
        fig,ax=plt.subplots(figsize=(max(10,0.35*piv.shape[1]+3), max(5,0.35*len(piv)+2)))
        im=ax.imshow(piv.values,aspect="auto")
        ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index)
        ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels([f"{a}\n{b}" for a,b in piv.columns],rotation=45,ha="right",fontsize=8)
        ax.set_title("Figure 14. Interval-wise pathway flux scores")
        fig.colorbar(im,ax=ax,fraction=0.025,pad=0.02).set_label("Mean abs pFBA flux")
        fig.tight_layout()
        fig.savefig(os.path.join(results_dir(DATASET,"figures"),"Fig14_interval_pathway_scores.png"),dpi=200)
        plt.close(fig)
        print(f"  [saved] results/{DATASET}/figures/Fig14_interval_pathway_scores.png")
except Exception as e:
    print(f"  !! interval figure failed: {e}")
print("  OK interval-wise FBA/FVA complete")
