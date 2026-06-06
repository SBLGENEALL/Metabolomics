"""
14_central_mab_fva_report.py — Central metabolism + mAb-related pathway flux/FVA report

Purpose
-------
This step is for the question: "If direct DM_igg_g maximization is flat, can we
still use the FBA model to compare pathway-level flux states between clones?"

It does NOT invent a new titer-fitting objective. Instead it exports a curated
central metabolism / mAb-production panel and runs FBA/pFBA/FVA on those
reactions under two clearly separated modes:

1) no_igg_input
   IgG titer/qIgG is not used as an input. The normal production objective is
   optimized and the central pathway flux state/FVA is reported as a diagnostic.
   If DM_igg_g max is flat, report it as flat; do not over-interpret.

2) measured_demand
   Observed qIgG is fixed as a scaled demand. This is explanation mode, not
   prediction. It is useful for comparing which central and mAb-related flux
   states can support the observed High/Low/Mother production levels.

Outputs
-------
results/<dataset>/tables/central_mab_reaction_panel.csv
results/<dataset>/tables/central_mab_flux_state.csv
results/<dataset>/tables/central_mab_fva_report.csv
results/<dataset>/figures/Fig5_central_mab_flux_heatmap.png
results/<dataset>/figures/Fig6_central_mab_flux_zscore.png
results/<dataset>/figures/Fig7_central_mab_flux_delta.png
"""
import argparse
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
from src.fba_utils import (
    load_model, save_table, save_figure, apply_bounds, avg_constraints,
    choose_objective_from_data, estimate_igg_interval_rates,
    compute_measured_demand_scale, optimize_fixed_igg_demand_pfba,
    set_biomass_minimum,
)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cobra.flux_analysis import flux_variability_analysis, pfba

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa", choices=["sowa2020", "practice_20aa", "own_experiment"])
parser.add_argument("--objective", default=os.environ.get("CHO_OBJECTIVE"), help="Production reaction id; default DM_igg_g")
parser.add_argument("--analysis_mode", default=os.environ.get("CHO_ANALYSIS_MODE", "both"), choices=["predict", "explain", "both"])
parser.add_argument("--fva_mode", default="both", choices=["no_igg_input", "measured_demand", "both"], help="Which central pathway FVA mode to export")
parser.add_argument("--fraction", type=float, default=FVA_FRACTION, help="FVA fraction of optimum")
parser.add_argument("--biomass_fraction", type=float, default=None)
parser.add_argument("--demand_scale", default=os.environ.get("CHO_DEMAND_SCALE", MEASURED_IGG_DEMAND_SCALE))
args = parser.parse_args()
DATASET = args.dataset
if args.biomass_fraction is None:
    args.biomass_fraction = float(os.environ.get("CHO_BIOMASS_FRACTION", BIOMASS_MIN_FRACTION))

print("=" * 65)
print(f"  14_central_mab_fva_report.py — {DATASET}")
print("=" * 65)

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
igG_day14 = data.get("igG_day14", {})
if not all_constraints or not clones:
    print("  !! constraints/clones 없음")
    sys.exit(1)

model = load_model(verbose=True)
OBJ = choose_objective_from_data(model, data, args.objective)
rxn_ids = {r.id for r in model.reactions}
DEMAND_RXN = "DM_igg_g" if "DM_igg_g" in rxn_ids else OBJ
igg_rates = estimate_igg_interval_rates(data)

# Curated panel. Each row contains multiple possible iCHO3K/Recon-style IDs.
# The script keeps whichever IDs exist in the loaded model.
PANEL = [
    # Exchange / measured phenotype
    ("Exchange", "Glucose uptake", ["EX_glc_e", "EX_glc__D_e"]),
    ("Exchange", "Lactate exchange", ["EX_lac_L_e", "EX_lac__L_e"]),
    ("Exchange", "Glutamine uptake", ["EX_gln_L_e", "EX_gln__L_e"]),
    ("Exchange", "Glutamate exchange", ["EX_glu_L_e", "EX_glu__L_e"]),
    ("Exchange", "Ammonia exchange", ["EX_nh4_e", "EX_nh4__e"]),
    ("Exchange", "Alanine exchange", ["EX_ala_L_e", "EX_ala__L_e"]),
    ("Exchange", "Leucine exchange", ["EX_leu_L_e", "EX_leu__L_e"]),
    ("Exchange", "Lysine exchange", ["EX_lys_L_e", "EX_lys__L_e"]),
    ("Exchange", "Valine exchange", ["EX_val_L_e", "EX_val__L_e"]),
    ("Exchange", "Serine exchange", ["EX_ser_L_e", "EX_ser__L_e"]),

    # Glycolysis / pyruvate / lactate
    ("Glycolysis", "Hexokinase", ["HEX1", "HEX1_c"]),
    ("Glycolysis", "Glucose-6-phosphate isomerase", ["PGI"]),
    ("Glycolysis", "PFK", ["PFK", "PFK_2"]),
    ("Glycolysis", "Aldolase", ["FBA", "FBA2"]),
    ("Glycolysis", "TPI", ["TPI"]),
    ("Glycolysis", "GAPDH", ["GAPD", "GAPDH"]),
    ("Glycolysis", "PGK", ["PGK"]),
    ("Glycolysis", "PGM", ["PGM", "PGM1"]),
    ("Glycolysis", "ENO", ["ENO"]),
    ("Glycolysis", "Pyruvate kinase", ["PYK", "PYK_c"]),
    ("Pyruvate/Lactate", "LDH", ["LDH_L", "LDH_Lm", "LDH_D"]),
    ("Pyruvate/Lactate", "Pyruvate dehydrogenase", ["PDHm", "PDH"]),
    ("Pyruvate/Lactate", "Pyruvate carboxylase", ["PCm", "PC"]),
    ("Pyruvate/Lactate", "Malic enzyme", ["ME1m", "ME2m", "ME1", "ME2"]),

    # TCA
    ("TCA", "Citrate synthase", ["CSm", "CS"]),
    ("TCA", "Aconitase A", ["ACONTam", "ACONTa", "ACONTm"]),
    ("TCA", "Aconitase B", ["ACONTbm", "ACONTb", "ACONTm"]),
    ("TCA", "Isocitrate DH", ["ICDHyrm", "ICDHxm", "ICDHyr"]),
    ("TCA", "Alpha-KG DH", ["AKGDm", "AKGDH"]),
    ("TCA", "Succinyl-CoA synthetase", ["SUCOASm", "SUCOAS1m", "SUCOAS"]),
    ("TCA", "Succinate DH", ["SUCD1m", "SUCD1", "SUCDi"]),
    ("TCA", "Fumarase", ["FUMm", "FUM"]),
    ("TCA", "Malate DH", ["MDHm", "MDH"]),

    # PPP / redox
    ("PPP/Redox", "G6P DH", ["G6PDH2r", "G6PDH2", "G6PDH"]),
    ("PPP/Redox", "6PGL", ["PGL", "PGLYCP"]),
    ("PPP/Redox", "6PG DH", ["GND", "GNDc"]),
    ("PPP/Redox", "Transketolase 1", ["TKT1", "TKT"]),
    ("PPP/Redox", "Transketolase 2", ["TKT2"]),
    ("PPP/Redox", "Transaldolase", ["TALA"]),
    ("PPP/Redox", "NADPH oxidoreductase", ["NADPHtru", "NADPHtxu"]),

    # Glutamine / nitrogen
    ("Glutamine/Nitrogen", "Glutamine synthetase", ["GLNS", "GLNS_c"]),
    ("Glutamine/Nitrogen", "Glutaminase", ["GLUNm", "GLUN", "GLUN_c"]),
    ("Glutamine/Nitrogen", "Glutamate DH", ["GLUDym", "GLUDy", "GLUDxm", "GLUDx"]),
    ("Glutamine/Nitrogen", "Transaminase ALT", ["ALATA_L", "ALATA_Lm", "ALATA"]),
    ("Glutamine/Nitrogen", "Transaminase AST", ["ASPTA", "ASPTAm"]),

    # Energy / mitochondria
    ("Energy", "ATP maintenance", ["ATPM"]),
    ("Energy", "ATP synthase", ["ATPS4m", "ATPS4mi"]),
    ("Energy", "Complex I", ["NADH2_u10m", "NADH2_u10mi"]),
    ("Energy", "Complex III", ["CYOOm3", "CYOR_u10m"]),
    ("Energy", "Complex IV", ["CYOOm", "CYOOm2"]),

    # mAb production reactions in iCHO3K prod
    ("mAb synthesis", "IgG heavy chain", ["igg_hc"]),
    ("mAb synthesis", "IgG light chain", ["igg_lc"]),
    ("mAb synthesis", "IgG formation", ["igg_formation"]),
    ("mAb synthesis", "IgG demand/export", ["DM_igg_g", "DM_for_igg"]),
]


def resolve_panel(model):
    rows = []
    seen = set()
    ids = {r.id for r in model.reactions}
    for pathway, label, candidates in PANEL:
        matched = False
        for rid in candidates:
            if rid in ids and rid not in seen:
                r = model.reactions.get_by_id(rid)
                rows.append({
                    "pathway": pathway,
                    "label": label,
                    "rxn_id": rid,
                    "reaction_name": r.name,
                    "reaction": r.reaction,
                    "match_type": "curated_exact",
                })
                seen.add(rid)
                matched = True
                break
        # If no exact match, try a conservative ID substring match for common labels.
        if not matched:
            token = candidates[0].replace("m", "").replace("_c", "")
            # only short curated tokens; avoid broad accidental matches
            if len(token) >= 4:
                hits = [r for r in model.reactions if token.lower() in r.id.lower()]
                if hits:
                    r = hits[0]
                    if r.id not in seen:
                        rows.append({
                            "pathway": pathway,
                            "label": label,
                            "rxn_id": r.id,
                            "reaction_name": r.name,
                            "reaction": r.reaction,
                            "match_type": f"fallback_id_contains:{token}",
                        })
                        seen.add(r.id)
    return pd.DataFrame(rows)

panel = resolve_panel(model)
if panel.empty:
    print("  !! central/mAb panel reactions not found")
    sys.exit(1)
save_table(panel, "central_mab_reaction_panel.csv", DATASET)
rxn_list = [model.reactions.get_by_id(rid) for rid in panel["rxn_id"].tolist() if rid in rxn_ids]

# Conditions: emphasize clone-level comparison for demo datasets while also exporting
# group averages (High/Low/Mother) when biological group labels are available.
clone_groups = data.get("clone_groups", {}) or {}
conditions = []
seen_cond = set()

def _add_condition(cond_name, clone_name, cst, qobs):
    if cond_name in seen_cond:
        return
    conditions.append((cond_name, clone_name, cst, qobs))
    seen_cond.add(cond_name)

# For practice/demo datasets and smaller experiments, include every clone explicitly.
if DATASET == "practice_20aa" or len(clones) <= 12:
    for c in clones:
        _add_condition(c, c, all_constraints.get(c, {}), float(igg_rates.get(c, 0.0) or 0.0))

# Add biological group averages when available.
group_to_clones = {}
for c in clones:
    g = str(clone_groups.get(c, "Mid") or "Mid")
    group_to_clones.setdefault(g, []).append(c)
for gname in ["High", "Mother", "Low", "Mid"] + sorted([g for g in group_to_clones if g not in {"High","Mother","Low","Mid"}]):
    clist = group_to_clones.get(gname, [])
    if clist:
        qobs = float(np.mean([igg_rates.get(c, 0.0) for c in clist if c in igg_rates]) or 0.0)
        _add_condition(f"{gname}Avg", f"{gname}Avg", avg_constraints(clist, all_constraints), qobs)

# Fallback to top/bottom aggregates if explicit groups were not available.
if high_clones and "HighAvg" not in seen_cond:
    _add_condition("HighAvg", "HighAvg", avg_constraints(high_clones, all_constraints), float(np.mean([igg_rates.get(c, 0.0) for c in high_clones if c in igg_rates]) or 0.0))
if low_clones and "LowAvg" not in seen_cond:
    _add_condition("LowAvg", "LowAvg", avg_constraints(low_clones, all_constraints), float(np.mean([igg_rates.get(c, 0.0) for c in low_clones if c in igg_rates]) or 0.0))

modes_to_run = []
if args.fva_mode in ["no_igg_input", "both"] and args.analysis_mode in ["predict", "both"]:
    modes_to_run.append("no_igg_input")
if args.fva_mode in ["measured_demand", "both"] and args.analysis_mode in ["explain", "both"]:
    modes_to_run.append("measured_demand")
if not modes_to_run:
    modes_to_run = ["no_igg_input"]

# Scale measured qIgG demands only for explanation mode.
strict_caps = {}
for cond_name, clone_name, cst, qobs in conditions:
    with model:
        apply_bounds(model, cst)
        set_biomass_minimum(model, cst, args.biomass_fraction)
        model.objective = OBJ
        sol = model.optimize()
        strict_caps[cond_name] = float(sol.objective_value) if sol.status == "optimal" and sol.objective_value is not None else 0.0
scale = compute_measured_demand_scale(
    {cond_name: qobs for cond_name, clone_name, cst, qobs in conditions},
    strict_caps,
    requested=args.demand_scale,
    safety=MEASURED_IGG_DEMAND_SCALE_SAFETY,
)

flux_rows = []
fva_rows = []
print(f"  Objective diagnostic : {OBJ}")
print(f"  Demand reaction      : {DEMAND_RXN}")
print(f"  Conditions           : {', '.join([c[0] for c in conditions])}")
print(f"  Modes                : {', '.join(modes_to_run)}")

for mode in modes_to_run:
    for cond_name, clone_name, cst, qobs in conditions:
        with model:
            n = apply_bounds(model, cst)
            set_biomass_minimum(model, cst, args.biomass_fraction)
            target = 0.0
            if mode == "measured_demand":
                target = max(0.0, float(qobs or 0.0) * scale)
                if DEMAND_RXN in {r.id for r in model.reactions}:
                    dr = model.reactions.get_by_id(DEMAND_RXN)
                    dr.lower_bound = max(float(dr.lower_bound), target)
                    dr.upper_bound = max(float(dr.lower_bound), target)
            model.objective = OBJ if mode == "no_igg_input" else DEMAND_RXN
            sol = model.optimize()
            if sol.status != "optimal":
                print(f"  !! {mode}/{cond_name}: {sol.status}")
                continue
            # pFBA gives one representative flux distribution for the curated panel.
            try:
                psol = pfba(model)
                flux_source = psol.fluxes
                total_flux = float(psol.objective_value) if psol.objective_value is not None else np.nan
            except Exception:
                flux_source = sol.fluxes
                total_flux = np.nan
            obj_val = float(sol.objective_value) if sol.objective_value is not None else 0.0
            for _, prow in panel.iterrows():
                rid = prow["rxn_id"]
                val = float(flux_source.get(rid, 0.0)) if rid in flux_source.index else 0.0
                flux_rows.append({
                    "mode": mode,
                    "condition": cond_name,
                    "clone": clone_name,
                    "pathway": prow["pathway"],
                    "label": prow["label"],
                    "rxn_id": rid,
                    "flux": val,
                    "abs_flux": abs(val),
                    "objective": OBJ,
                    "objective_value": obj_val,
                    "measured_qIgG": qobs,
                    "scaled_demand": target,
                    "pfba_total_flux": total_flux,
                    "n_constraints": n,
                })
            try:
                fva = flux_variability_analysis(
                    model,
                    reaction_list=rxn_list,
                    fraction_of_optimum=args.fraction,
                    processes=FVA_PROCESSES,
                )
                for rid, row in fva.iterrows():
                    pinfo = panel[panel["rxn_id"] == rid].iloc[0].to_dict()
                    mn, mx = float(row["minimum"]), float(row["maximum"])
                    fva_rows.append({
                        "mode": mode,
                        "condition": cond_name,
                        "clone": clone_name,
                        "pathway": pinfo["pathway"],
                        "label": pinfo["label"],
                        "rxn_id": rid,
                        "minimum": mn,
                        "maximum": mx,
                        "mean": (mn + mx) / 2.0,
                        "range": mx - mn,
                        "objective": OBJ,
                        "objective_value": obj_val,
                        "measured_qIgG": qobs,
                        "scaled_demand": target,
                    })
                print(f"  OK {mode:15s} {cond_name:14s} obj={obj_val:.4g} panel_rxns={len(fva)}")
            except Exception as e:
                print(f"  !! FVA failed {mode}/{cond_name}: {e}")

flux_df = pd.DataFrame(flux_rows)
fva_df = pd.DataFrame(fva_rows)
save_table(flux_df, "central_mab_flux_state.csv", DATASET)
save_table(fva_df, "central_mab_fva_report.csv", DATASET)

# Also create concise high-vs-low delta table when possible.
if not flux_df.empty:
    for mode in sorted(flux_df["mode"].unique()):
        sub = flux_df[flux_df["mode"] == mode]
        conds = set(sub["condition"])
        hname = "HighAvg" if "HighAvg" in conds else ("HighProducer" if "HighProducer" in conds else None)
        lname = "LowAvg" if "LowAvg" in conds else ("LowProducer" if "LowProducer" in conds else None)
        if hname and lname:
            h = sub[sub["condition"] == hname].set_index("rxn_id")
            l = sub[sub["condition"] == lname].set_index("rxn_id")
            common = h.index.intersection(l.index)
            rows = []
            for rid in common:
                rows.append({
                    "mode": mode,
                    "rxn_id": rid,
                    "pathway": h.loc[rid, "pathway"],
                    "label": h.loc[rid, "label"],
                    "high_flux": float(h.loc[rid, "flux"]),
                    "low_flux": float(l.loc[rid, "flux"]),
                    "delta_high_minus_low": float(h.loc[rid, "flux"] - l.loc[rid, "flux"]),
                    "abs_delta": abs(float(h.loc[rid, "flux"] - l.loc[rid, "flux"])),
                })
            delta_df = pd.DataFrame(rows).sort_values("abs_delta", ascending=False)
            save_table(delta_df, f"central_mab_flux_delta_high_low_{mode}.csv", DATASET)

# Figure 5: central/mAb pFBA flux heatmap.
if not flux_df.empty:
    plot_mode = "measured_demand" if "measured_demand" in set(flux_df["mode"]) else sorted(flux_df["mode"].unique())[0]
    sub = flux_df[flux_df["mode"] == plot_mode].copy()
    sub["row_label"] = sub["pathway"] + " | " + sub["label"]
    # Keep the highest-activity reactions to avoid an unreadable figure.
    core_keep = sub[sub["pathway"].isin(["Glycolysis","Pyruvate/Lactate","TCA","PPP/Redox","Glutamine/Nitrogen","Energy","mAb synthesis","Exchange"])].copy()
    ranked = (core_keep.groupby("row_label")["abs_flux"].max().sort_values(ascending=False).index.tolist())
    # keep up to 50 rows but always preserve all mAb synthesis reactions
    mandatory = core_keep[core_keep["pathway"]=="mAb synthesis"]["row_label"].unique().tolist()
    order_rows = []
    for x in ranked + mandatory:
        if x not in order_rows:
            order_rows.append(x)
    order_rows = order_rows[:50]
    sub = core_keep[core_keep["row_label"].isin(order_rows)]
    piv = sub.pivot_table(index="row_label", columns="condition", values="flux", aggfunc="first").fillna(0.0)
    piv = piv.loc[[r for r in order_rows if r in piv.index]]
    fig_h = max(7, 0.28 * len(piv) + 2)
    fig, ax = plt.subplots(figsize=(11, fig_h))
    im = ax.imshow(piv.values, aspect="auto")
    ax.set_yticks(np.arange(len(piv.index)))
    ax.set_yticklabels(piv.index, fontsize=8)
    ax.set_xticks(np.arange(len(piv.columns)))
    ax.set_xticklabels(piv.columns, rotation=35, ha="right")
    ax.set_title("Figure 5. iCHO3K Central/mAb Flux Heatmap", weight="bold")
    ax.set_xlabel("Condition")
    ax.set_ylabel("Pathway | Reaction")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Flux")
    fig.tight_layout()
    save_figure(fig, "Fig5_central_mab_flux_heatmap.png", DATASET)

    # Figure 6: row-wise z-score heatmap. Absolute flux heatmaps are often
    # dominated by a few large energy/TCA reactions; z-scoring each reaction
    # exposes clone/group-specific relative differences. This is a visualization
    # aid, not a different FBA calculation.
    try:
        row_mean = piv.mean(axis=1)
        row_std = piv.std(axis=1).replace(0, np.nan)
        zpiv = piv.sub(row_mean, axis=0).div(row_std, axis=0).fillna(0.0)
        save_table(zpiv.reset_index().rename(columns={"row_label": "reaction"}), "central_mab_flux_state_zscore_matrix.csv", DATASET)
        fig, ax = plt.subplots(figsize=(11, fig_h))
        vmax = max(1.0, float(np.nanmax(np.abs(zpiv.values)))) if zpiv.size else 1.0
        im = ax.imshow(zpiv.values, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
        ax.set_yticks(np.arange(len(zpiv.index)))
        ax.set_yticklabels(zpiv.index, fontsize=8)
        ax.set_xticks(np.arange(len(zpiv.columns)))
        ax.set_xticklabels(zpiv.columns, rotation=35, ha="right")
        ax.set_title("Figure 6. Row-wise Relative Central/mAb Flux Pattern", weight="bold")
        ax.set_xlabel("Condition")
        ax.set_ylabel("Pathway | Reaction")
        cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
        cbar.set_label("Row-wise z-score of flux")
        fig.tight_layout()
        save_figure(fig, "Fig6_central_mab_flux_zscore.png", DATASET)
    except Exception as e:
        print(f"  !! Fig6 skipped: {e}")

# Figure 7: High-vs-low delta for central/mAb fluxes.
try:
    plot_mode = "measured_demand" if os.path.exists(os.path.join(results_dir(DATASET, "tables"), "central_mab_flux_delta_high_low_measured_demand.csv")) else "no_igg_input"
    delta_path = os.path.join(results_dir(DATASET, "tables"), f"central_mab_flux_delta_high_low_{plot_mode}.csv")
    if os.path.exists(delta_path):
        d = pd.read_csv(delta_path).sort_values("abs_delta", ascending=False).head(25).iloc[::-1]
        fig, ax = plt.subplots(figsize=(10, max(6, 0.32 * len(d) + 1.5)))
        ax.barh(d["pathway"] + " | " + d["label"], d["delta_high_minus_low"])
        ax.axvline(0, color="black", lw=1)
        ax.set_xlabel("High - Low flux")
        ax.set_title("Figure 7. High-vs-Low Central/mAb Flux Difference", weight="bold")
        ax.grid(True, axis="x", ls=":", alpha=0.5)
        fig.tight_layout()
        save_figure(fig, "Fig7_central_mab_flux_delta.png", DATASET)
except Exception as e:
    print(f"  !! Fig7 skipped: {e}")

# Markdown interpretation guide.
guide = f"""# Central metabolism + mAb FVA report

This report is designed to avoid the common over-interpretation that `DM_igg_g`\nmaximization alone predicts titer. It exports central pathway and mAb-related\nflux states/ranges for direct clone comparison.\n\n## Key files\n\n- `central_mab_reaction_panel.csv`: curated reaction list found in the loaded iCHO3K model.\n- `central_mab_flux_state.csv`: representative pFBA flux state for each condition.\n- `central_mab_fva_report.csv`: minimum/maximum feasible flux range for each curated reaction.\n- `central_mab_flux_delta_high_low_*.csv`: High-Low differences for the pathway panel.\n\n## Interpretation\n\n- `no_igg_input` mode does not use IgG titer as an input. If `DM_igg_g` maximum\n  is flat, report it as a diagnostic failure of stand-alone FBA capacity ranking.\n- `measured_demand` mode fixes observed qIgG as a scaled demand and reconstructs\n  pathway-level flux states. This is explanation mode, not prediction.\n- Use this report to compare glycolysis, lactate, TCA, PPP, glutamine/nitrogen,\n  energy, and mAb synthesis fluxes between High/Low/Mother.\n\nGenerated objective: `{OBJ}`; demand reaction: `{DEMAND_RXN}`.\n"""
with open(os.path.join(results_dir(DATASET, "tables"), "README_CENTRAL_MAB_FVA.md"), "w", encoding="utf-8") as f:
    f.write(guide)
print("  [saved] results/%s/tables/README_CENTRAL_MAB_FVA.md" % DATASET)
print("  OK 14_central_mab_fva_report 완료")
