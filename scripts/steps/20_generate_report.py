"""
20_generate_report.py — Generate a compact Markdown report for the CHO FBA/FVA run.
"""
import argparse
import os
import sys
from datetime import datetime
import pandas as pd


def _find_root():
    cur = os.path.dirname(os.path.abspath(__file__))
    for _ in range(7):
        if os.path.exists(os.path.join(cur, "src", "config.py")):
            return cur
        cur = os.path.dirname(cur)
    return cur

ROOT=_find_root()
if ROOT not in sys.path: sys.path.insert(0,ROOT)
from src.config import results_dir  # noqa

parser=argparse.ArgumentParser()
parser.add_argument("--dataset",default="practice_20aa")
parser.add_argument("--analysis_mode",default=os.environ.get("CHO_ANALYSIS_MODE","both"))
args=parser.parse_args(); DATASET=args.dataset
print("="*65); print(f"  20_generate_report.py — {DATASET}"); print("="*65)

tables=results_dir(DATASET,"tables"); figs=results_dir(DATASET,"figures")

def read_csv(name):
    p=os.path.join(tables,name)
    return pd.read_csv(p) if os.path.exists(p) else pd.DataFrame()

def top_rows(df,n=8):
    if df.empty: return "_Not available._\n"
    return df.head(n).to_markdown(index=False)

qc=read_csv("data_qc_warnings.csv")
path_delta=read_csv("pathway_score_high_low_delta.csv")
fva_sep=read_csv("fva_high_low_overlap_separation.csv")
sens=read_csv("constraint_sensitivity_summary.csv")
rate=read_csv("exchange_rates.csv")
figs_list=sorted([f for f in os.listdir(figs) if f.lower().endswith(".png")]) if os.path.isdir(figs) else []

if not path_delta.empty:
    path_delta=path_delta.reindex(path_delta["delta_high_minus_low"].abs().sort_values(ascending=False).index)
if not fva_sep.empty:
    fva_sep=fva_sep.sort_values(["separation_gap","range_delta_high_minus_low"],ascending=False)
if not sens.empty:
    sens=sens[sens["robust_direction"]!="unstable"].copy().sort_values("delta_mean",key=lambda s:s.abs(),ascending=False)

lines=[]
lines.append(f"# CHO FBA/FVA Run Summary — {DATASET}\n")
lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
lines.append("## Scientific scope\n")
lines.append("This pipeline is a mechanistic clone-comparison workflow, not a stand-alone IgG titer predictor. It uses measured exchange-rate constraints in the iCHO3K production model to compare pFBA flux states, FVA feasible ranges, pathway scores, and focused Escher maps.\n")
lines.append("## Key outputs to inspect\n")
for f in ["Fig1_IgG_timecourse.png","Fig2_rate_heatmap.png","Fig3_lac_glc_ratio.png","Fig5_high_vs_low_rates.png","Fig10_central_mab_flux_heatmap.png","Fig10B_central_mab_flux_zscore.png","Fig12_focused_core_fva_range.png","Fig12B_focused_core_fva_range_zscore.png","Fig12C_focused_fva_high_low_delta.png","Fig14_interval_pathway_scores.png","Fig15_pathway_scores_fva_overlap.png"]:
    if f in figs_list: lines.append(f"- `results/{DATASET}/figures/{f}`")
lines.append("\n## Data QC\n")
if qc.empty:
    lines.append("No QC warnings file found. Run step 17 for data QC.\n")
else:
    lines.append(f"QC warnings: {len(qc)} rows. ERROR={int((qc.get('level')=='ERROR').sum()) if 'level' in qc else 'NA'}, WARN={int((qc.get('level')=='WARN').sum()) if 'level' in qc else 'NA'}\n")
    lines.append(top_rows(qc,8)+"\n")
lines.append("## Pathway score summary: HighAvg vs LowAvg\n")
lines.append(top_rows(path_delta[[c for c in ["pathway","high_mean_abs_flux","low_mean_abs_flux","delta_high_minus_low","fold_high_over_low"] if c in path_delta.columns]],10)+"\n")
lines.append("## FVA overlap/separation summary\n")
cols=[c for c in ["pathway","label","rxn_id","overlap_fraction","separation_gap","separation_direction","range_delta_high_minus_low"] if c in fva_sep.columns]
lines.append(top_rows(fva_sep[cols],10)+"\n")
lines.append("## Robust sensitivity hits\n")
cols=[c for c in ["pathway","label","rxn_id","delta_mean","delta_sd","positive_fraction","negative_fraction","robust_direction"] if c in sens.columns]
lines.append(top_rows(sens[cols],10)+"\n")
lines.append("## Focused Escher map usage\n")
lines.append(f"Map: `results/{DATASET}/escher_maps/focused/CHO_focus_core_carbon_map.json`\n")
lines.append(f"Flux data: `results/{DATASET}/tables/focused_escher/focused_escher_flux_high_minus_low_cho.json`\n")
lines.append(f"FVA range delta data: `results/{DATASET}/tables/focused_escher/focused_fva_measured_demand_high_minus_low_range_cho.json`\n")
lines.append("## Interpretation cautions\n")
lines.append("- pFBA is one representative solution, not a unique intracellular flux measurement.\n- FVA ranges can be broad in under-constrained or loop-prone reactions; prioritize focused pathway panels and overlap/separation metrics.\n- Without spent-media exchange rates, FBA/FVA becomes weakly constrained and should not be over-interpreted.\n")

out=os.path.join(results_dir(DATASET,"."),"REPORT_SUMMARY.md")
with open(out,"w",encoding="utf-8") as f: f.write("\n".join(lines))
print(f"  [saved] results/{DATASET}/REPORT_SUMMARY.md")
print("  OK report complete")
