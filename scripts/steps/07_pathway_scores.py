"""
19_pathway_scores_sensitivity.py — Pathway scores, FVA separation/overlap, and constraint sensitivity

Runs after step 14. It turns reaction-level pFBA/FVA outputs into more interpretable
pathway-level summaries and checks whether High-Low flux deltas are robust to
small exchange-bound perturbations.
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
from src.config import DATA_PROCESSED, results_dir, BIOMASS_MIN_FRACTION, MEASURED_IGG_DEMAND_SCALE  # noqa
from src.fba_utils import load_model, apply_bounds, avg_constraints, choose_objective_from_data, estimate_igg_interval_rates, compute_measured_demand_scale, set_biomass_minimum  # noqa
from cobra.flux_analysis import pfba

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa", choices=["sowa2020", "practice_20aa", "own_experiment"])
parser.add_argument("--analysis_mode", default=os.environ.get("CHO_ANALYSIS_MODE", "both"), choices=["predict", "explain", "both"])
parser.add_argument("--objective", default=os.environ.get("CHO_OBJECTIVE"))
parser.add_argument("--biomass_fraction", type=float, default=float(os.environ.get("CHO_BIOMASS_FRACTION", BIOMASS_MIN_FRACTION)))
parser.add_argument("--demand_scale", default=os.environ.get("CHO_DEMAND_SCALE", MEASURED_IGG_DEMAND_SCALE))
parser.add_argument("--sensitivity_n", type=int, default=30)
parser.add_argument("--sensitivity_noise", type=float, default=0.10)
args = parser.parse_args()
DATASET=args.dataset
print("="*65)
print(f"  19_pathway_scores_sensitivity.py — {DATASET}")
print("="*65)

tables=results_dir(DATASET,"tables"); figs=results_dir(DATASET,"figures")
flux_file=os.path.join(tables,"central_mab_flux_state.csv")
fva_file=os.path.join(tables,"central_mab_fva_report.csv")
if not os.path.exists(flux_file):
    raise FileNotFoundError("Run step 14 first: central_mab_flux_state.csv not found")
flux=pd.read_csv(flux_file)
fva=pd.read_csv(fva_file) if os.path.exists(fva_file) else pd.DataFrame()
mode="measured_demand" if "measured_demand" in set(flux.get("mode",[])) else flux["mode"].iloc[0]
sub=flux[flux["mode"]==mode].copy()

# Pathway scores
scores=sub.groupby(["condition","pathway"],as_index=False).agg(
    mean_abs_flux=("abs_flux","mean"),
    sum_abs_flux=("abs_flux","sum"),
    mean_flux=("flux","mean"),
    n_reactions=("rxn_id","nunique"),
)
if {"HighAvg","LowAvg"}.issubset(set(scores["condition"])):
    h=scores[scores["condition"]=="HighAvg"].set_index("pathway")
    l=scores[scores["condition"]=="LowAvg"].set_index("pathway")
    common=h.index.intersection(l.index)
    rows=[]
    for p in common:
        rows.append({"pathway":p,"high_mean_abs_flux":h.loc[p,"mean_abs_flux"],"low_mean_abs_flux":l.loc[p,"mean_abs_flux"],"delta_high_minus_low":h.loc[p,"mean_abs_flux"]-l.loc[p,"mean_abs_flux"],"fold_high_over_low":(h.loc[p,"mean_abs_flux"]+1e-12)/(l.loc[p,"mean_abs_flux"]+1e-12)})
    pd.DataFrame(rows).sort_values("delta_high_minus_low",key=lambda s:s.abs(),ascending=False).to_csv(os.path.join(tables,"pathway_score_high_low_delta.csv"),index=False)
scores.to_csv(os.path.join(tables,"pathway_scores.csv"),index=False)
print(f"  [saved] results/{DATASET}/tables/pathway_scores.csv")

# FVA overlap/separation
if not fva.empty:
    fmode="measured_demand" if "measured_demand" in set(fva["mode"]) else fva["mode"].iloc[0]
    fs=fva[fva["mode"]==fmode].copy()
    if {"HighAvg","LowAvg"}.issubset(set(fs["condition"])):
        h=fs[fs["condition"]=="HighAvg"].set_index("rxn_id")
        l=fs[fs["condition"]=="LowAvg"].set_index("rxn_id")
        common=h.index.intersection(l.index)
        rows=[]
        for rid in common:
            hmin,hmax=float(h.loc[rid,"minimum"]),float(h.loc[rid,"maximum"])
            lmin,lmax=float(l.loc[rid,"minimum"]),float(l.loc[rid,"maximum"])
            overlap=max(0.0, min(hmax,lmax)-max(hmin,lmin))
            union=max(hmax,lmax)-min(hmin,lmin)
            sep=0.0
            direction="overlap"
            if hmin > lmax:
                sep=hmin-lmax; direction="High_range_above_Low"
            elif lmin > hmax:
                sep=lmin-hmax; direction="Low_range_above_High"
            rows.append({
                "rxn_id":rid,"pathway":h.loc[rid,"pathway"],"label":h.loc[rid,"label"],
                "high_min":hmin,"high_max":hmax,"high_range":hmax-hmin,
                "low_min":lmin,"low_max":lmax,"low_range":lmax-lmin,
                "overlap_width":overlap,"union_width":union,"overlap_fraction":overlap/union if union>1e-12 else np.nan,
                "separation_gap":sep,"separation_direction":direction,
                "range_delta_high_minus_low":(hmax-hmin)-(lmax-lmin),
            })
        odf=pd.DataFrame(rows).sort_values(["separation_gap","range_delta_high_minus_low"],ascending=False)
        odf.to_csv(os.path.join(tables,"fva_high_low_overlap_separation.csv"),index=False)
        path_sep=odf.groupby("pathway",as_index=False).agg(mean_overlap_fraction=("overlap_fraction","mean"),mean_abs_range_delta=("range_delta_high_minus_low",lambda x: float(np.mean(np.abs(x)))),n_nonoverlap=("separation_gap",lambda x:int((x>1e-12).sum())),n_reactions=("rxn_id","count"))
        path_sep.to_csv(os.path.join(tables,"fva_pathway_flexibility_summary.csv"),index=False)
        print(f"  [saved] results/{DATASET}/tables/fva_high_low_overlap_separation.csv")

        # Figure 8: pathway score + overlap summary
        fig,axes=plt.subplots(1,2,figsize=(14,5))
        ps=pd.read_csv(os.path.join(tables,"pathway_score_high_low_delta.csv")) if os.path.exists(os.path.join(tables,"pathway_score_high_low_delta.csv")) else pd.DataFrame()
        if not ps.empty:
            ps=ps.sort_values("delta_high_minus_low")
            axes[0].barh(ps["pathway"],ps["delta_high_minus_low"])
            axes[0].axvline(0,lw=0.7)
            axes[0].set_title("Pathway pFBA score delta\nHighAvg - LowAvg")
            axes[0].set_xlabel("Delta mean abs flux")
        if not path_sep.empty:
            p2=path_sep.sort_values("mean_abs_range_delta")
            axes[1].barh(p2["pathway"],p2["mean_abs_range_delta"])
            axes[1].set_title("FVA flexibility difference by pathway")
            axes[1].set_xlabel("Mean |range delta|")
        fig.suptitle("Figure 8. Pathway pFBA/FVA Summary", fontweight="bold")
        fig.tight_layout()
        fig.savefig(os.path.join(figs,"Fig8_pathway_scores_fva_overlap.png"),dpi=200)
        plt.close(fig)
        print(f"  [saved] results/{DATASET}/figures/Fig8_pathway_scores_fva_overlap.png")

# Sensitivity: perturb HighAvg and LowAvg constraints and see whether panel deltas keep their sign.
try:
    pkl=os.path.join(DATA_PROCESSED,f"rates_{DATASET}.pkl")
    with open(pkl,"rb") as f: data=pickle.load(f)
    model=load_model(verbose=False); rxn_ids={r.id for r in model.reactions}
    OBJ=choose_objective_from_data(model,data,args.objective); DEMAND_RXN="DM_igg_g" if "DM_igg_g" in rxn_ids else OBJ
    panel=pd.read_csv(os.path.join(tables,"central_mab_reaction_panel.csv"))
    panel=panel[panel["rxn_id"].isin(rxn_ids)].copy()
    clones=data.get("clones",[]); clone_groups=data.get("clone_groups",{}) or {}; all_c=data.get("all_constraints",{}); igg_rates=estimate_igg_interval_rates(data)
    g2c={}
    for c in clones: g2c.setdefault(str(clone_groups.get(c,"Mid") or "Mid"),[]).append(c)
    if "High" in g2c and "Low" in g2c:
        high_c=avg_constraints(g2c["High"],all_c); low_c=avg_constraints(g2c["Low"],all_c)
        high_q=float(np.mean([igg_rates.get(c,0.0) for c in g2c["High"]]) or 0.0)
        low_q=float(np.mean([igg_rates.get(c,0.0) for c in g2c["Low"]]) or 0.0)
        def perturb(cst, rng):
            out={}
            for rid,(lb,ub) in cst.items():
                factor=1.0+rng.uniform(-args.sensitivity_noise,args.sensitivity_noise)
                out[rid]=(float(lb)*factor,float(ub)*factor)
            return out
        def solve(cst,q):
            with model:
                apply_bounds(model,cst); set_biomass_minimum(model,cst,args.biomass_fraction)
                if DEMAND_RXN in rxn_ids and mode=="measured_demand":
                    dr=model.reactions.get_by_id(DEMAND_RXN)
                    # use a conservative fixed scale; robust enough for sign checks
                    target=max(0.0,float(q or 0.0)*0.10)
                    dr.lower_bound=max(float(dr.lower_bound),target); dr.upper_bound=max(float(dr.lower_bound),target)
                model.objective=DEMAND_RXN if mode=="measured_demand" else OBJ
                sol=model.optimize()
                if sol.status!="optimal": return None
                try: psol=pfba(model); fl=psol.fluxes
                except Exception: fl=sol.fluxes
                return {rid:float(fl.get(rid,0.0)) for rid in panel["rxn_id"]}
        rng=np.random.default_rng(42)
        rows=[]
        for i in range(args.sensitivity_n):
            h=solve(perturb(high_c,rng),high_q); l=solve(perturb(low_c,rng),low_q)
            if h is None or l is None: continue
            for rid in panel["rxn_id"]:
                rows.append({"iteration":i,"rxn_id":rid,"delta_high_minus_low":h.get(rid,0.0)-l.get(rid,0.0)})
        sraw=pd.DataFrame(rows)
        sraw.to_csv(os.path.join(tables,"constraint_sensitivity_raw.csv"),index=False)
        if not sraw.empty:
            summ=sraw.groupby("rxn_id",as_index=False).agg(delta_mean=("delta_high_minus_low","mean"),delta_sd=("delta_high_minus_low","std"),positive_fraction=("delta_high_minus_low",lambda x:float((x>0).mean())),negative_fraction=("delta_high_minus_low",lambda x:float((x<0).mean())),n=("delta_high_minus_low","count"))
            pmap=panel.set_index("rxn_id")[["pathway","label"]].to_dict("index")
            summ["pathway"]=summ["rxn_id"].map(lambda r:pmap.get(r,{}).get("pathway",""))
            summ["label"]=summ["rxn_id"].map(lambda r:pmap.get(r,{}).get("label",r))
            summ["robust_direction"] = np.where(summ["positive_fraction"]>=0.8,"High>Low",np.where(summ["negative_fraction"]>=0.8,"Low>High","unstable"))
            summ.sort_values("delta_mean",key=lambda s:s.abs(),ascending=False).to_csv(os.path.join(tables,"constraint_sensitivity_summary.csv"),index=False)
            print(f"  [saved] results/{DATASET}/tables/constraint_sensitivity_summary.csv")
except Exception as e:
    print(f"  !! sensitivity skipped/failed: {e}")

print("  OK pathway scores / FVA overlap / sensitivity complete")
