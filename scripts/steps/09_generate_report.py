"""
09_generate_report.py — Generate a compact Markdown report for the CHO FBA/FVA run.

This step intentionally avoids optional pandas markdown dependencies such as
`tabulate`, so it can run in offline workstation environments with minimal
packages.
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


ROOT = _find_root()
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from src.config import results_dir  # noqa

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa")
parser.add_argument("--analysis_mode", default=os.environ.get("CHO_ANALYSIS_MODE", "both"))
args = parser.parse_args()
DATASET = args.dataset
print("=" * 65)
print(f"  09_generate_report.py — {DATASET}")
print("=" * 65)

tables = results_dir(DATASET, "tables")
figs = results_dir(DATASET, "figures")


def read_csv(name):
    p = os.path.join(tables, name)
    return pd.read_csv(p) if os.path.exists(p) else pd.DataFrame()


def _markdown_table(df):
    """Small dependency-free markdown table writer for report summaries."""
    if df.empty:
        return "_Not available._\n"
    d = df.copy()
    d = d.fillna("")
    d = d.astype(str)
    cols = list(d.columns)
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join(["---"] * len(cols)) + " |"
    rows = []
    for _, row in d.iterrows():
        vals = [str(row[c]).replace("\n", " ").replace("|", "/") for c in cols]
        rows.append("| " + " | ".join(vals) + " |")
    return "\n".join([header, sep] + rows) + "\n"


def top_rows(df, n=8):
    if df.empty:
        return "_Not available._\n"
    return _markdown_table(df.head(n))


qc = read_csv("data_qc_warnings.csv")
path_delta = read_csv("pathway_score_high_low_delta.csv")
fva_sep = read_csv("fva_high_low_overlap_separation.csv")
sens = read_csv("constraint_sensitivity_summary.csv")
figs_list = sorted([f for f in os.listdir(figs) if f.lower().endswith(".png")]) if os.path.isdir(figs) else []

if not path_delta.empty and "delta_high_minus_low" in path_delta.columns:
    path_delta = path_delta.reindex(path_delta["delta_high_minus_low"].abs().sort_values(ascending=False).index)
if not fva_sep.empty and {"separation_gap", "range_delta_high_minus_low"}.issubset(fva_sep.columns):
    fva_sep = fva_sep.sort_values(["separation_gap", "range_delta_high_minus_low"], ascending=False)
if not sens.empty and "robust_direction" in sens.columns:
    sens = sens[sens["robust_direction"] != "unstable"].copy()
    if "delta_mean" in sens.columns:
        sens = sens.sort_values("delta_mean", key=lambda s: s.abs(), ascending=False)

figure_order = [
    "Fig1_IgG_timecourse.png",
    "Fig2_rate_heatmap.png",
    "Fig3_lac_glc_ratio.png",
    "Fig4_high_vs_low_rates.png",
    "Fig5_central_mab_flux_heatmap.png",
    "Fig6_central_mab_flux_zscore.png",
    "Fig7_central_mab_flux_delta.png",
    "Fig8_pathway_scores_fva_overlap.png",
    "Fig9_full_fva_high_low_separation.png",
    "Fig10_full_fva_range_delta.png",
    "Fig11_summary_panel.png",
    # Legacy names are kept as fallback if output renaming has not run yet.
    "Fig5_high_vs_low_rates.png",
    "Fig10_central_mab_flux_heatmap.png",
    "Fig10B_central_mab_flux_zscore.png",
    "Fig11_central_mab_flux_delta.png",
    "Fig15_pathway_scores_fva_overlap.png",
    "Fig21_full_fva_high_low_separation.png",
    "Fig21B_full_fva_range_delta.png",
    "Fig8_summary_panel.png",
]

lines = []
lines.append(f"# CHO FBA/FVA Run Summary — {DATASET}\n")
lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
lines.append("## Scientific scope\n")
lines.append(
    "This pipeline is a mechanistic clone-comparison workflow, not a stand-alone IgG titer predictor. "
    "It uses measured exchange-rate constraints in the iCHO3K production model to compare pFBA flux states, "
    "FVA feasible ranges, pathway scores, and focused/full FVA outputs.\n"
)
lines.append("## Key outputs to inspect\n")
for f in figure_order:
    if f in figs_list:
        lines.append(f"- `results/{DATASET}/figures/{f}`")
lines.append("\n## Data QC\n")
if qc.empty:
    lines.append("No QC warnings file found. Run step 00 for data QC.\n")
else:
    n_error = int((qc.get("level") == "ERROR").sum()) if "level" in qc else "NA"
    n_warn = int((qc.get("level") == "WARN").sum()) if "level" in qc else "NA"
    lines.append(f"QC warnings: {len(qc)} rows. ERROR={n_error}, WARN={n_warn}\n")
    lines.append(top_rows(qc, 8) + "\n")

lines.append("## Pathway score summary: HighAvg vs LowAvg\n")
cols = [c for c in ["pathway", "high_mean_abs_flux", "low_mean_abs_flux", "delta_high_minus_low", "fold_high_over_low"] if c in path_delta.columns]
lines.append(top_rows(path_delta[cols] if cols else pd.DataFrame(), 10) + "\n")

lines.append("## FVA overlap/separation summary\n")
cols = [c for c in ["pathway", "label", "rxn_id", "overlap_fraction", "separation_gap", "separation_direction", "range_delta_high_minus_low"] if c in fva_sep.columns]
lines.append(top_rows(fva_sep[cols] if cols else pd.DataFrame(), 10) + "\n")

lines.append("## Robust sensitivity hits\n")
cols = [c for c in ["pathway", "label", "rxn_id", "delta_mean", "delta_sd", "positive_fraction", "negative_fraction", "robust_direction"] if c in sens.columns]
lines.append(top_rows(sens[cols] if cols else pd.DataFrame(), 10) + "\n")

lines.append("## Full/internal FVA outputs\n")
lines.append(f"- `results/{DATASET}/tables/full_fva/`\n")
lines.append("## Interpretation cautions\n")
lines.append(
    "- pFBA is one representative solution, not a unique intracellular flux measurement.\n"
    "- FVA ranges can be broad in under-constrained or loop-prone reactions; prioritize focused pathway panels and overlap/separation metrics.\n"
    "- Without spent-media exchange rates, FBA/FVA becomes weakly constrained and should not be over-interpreted.\n"
)

out = os.path.join(results_dir(DATASET, "."), "REPORT_SUMMARY.md")
with open(out, "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
print(f"  [saved] results/{DATASET}/REPORT_SUMMARY.md")
print("  OK report complete")
