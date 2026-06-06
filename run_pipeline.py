import argparse
import glob
import os
import shutil
import subprocess
import sys
import time

parser = argparse.ArgumentParser(description="CHO Metabolomics clean FBA/FVA workflow")
parser.add_argument("--dataset", default="practice_20aa", choices=["sowa2020", "practice_20aa", "own_experiment"])
parser.add_argument("--steps", default="00,01,02,03,04,05,06,07,08,09")
parser.add_argument("--input_format", default="auto", choices=["auto", "xlsx", "tsv"])
parser.add_argument("--rate_days", default="7,10")
parser.add_argument("--top_n", type=int, default=3)
parser.add_argument("--objective", default=None)
parser.add_argument("--analysis_mode", default="both", choices=["predict", "explain", "both"])
parser.add_argument("--biomass_fraction", type=float, default=None)
parser.add_argument("--demand_fraction", type=float, default=None)
parser.add_argument("--demand_scale", default="auto")
parser.add_argument("--constraint_policy", default="production_relaxed", choices=["strict", "production_relaxed"])
parser.add_argument("--feed_volume_mode", default="interval", choices=["interval", "cumulative"])
parser.add_argument("--fva_scope", default="all", choices=["exchange", "focused", "internal", "all"])
parser.add_argument("--fva_targets", default="group_avg", choices=["group_avg", "representative", "all_clones"])
parser.add_argument("--fva_processes", type=int, default=16)
parser.add_argument("--fva_fraction", type=float, default=None)
parser.add_argument(
    "--fva_source",
    default="auto",
    choices=["auto", "full", "focused"],
    help="FVA source for CHOmpact steps 10-14; auto prefers full all/internal FVA.",
)
args = parser.parse_args()

ROOT = os.path.dirname(os.path.abspath(__file__))
STEPS_DIR = os.path.join(ROOT, "scripts", "steps")
LOGS_DIR = os.path.join(ROOT, "logs")
os.makedirs(LOGS_DIR, exist_ok=True)

STEP_MAP = {
    "00": ("00_data_qc.py", "Data QC and sanity checks"),
    "01": ("01_load_and_rate.py", "Data load + feed-corrected exchange-rate calculation"),
    "02": ("02_map_metabolites.py", "Metabolite-to-exchange mapping"),
    "03": ("03_run_pfba.py", "pFBA prediction diagnostic + measured-demand explanation"),
    "04": ("04_export_escher_flux.py", "Escher flux reaction-data JSON"),
    "05": ("05_focused_fva.py", "Focused central/mAb pathway FVA"),
    "06": ("06_full_fva.py", "Workstation-scale internal/full FVA"),
    "07": ("07_pathway_scores.py", "Pathway scores + FVA overlap + sensitivity"),
    "08": ("08_make_figures.py", "Final figures"),
    "09": ("09_generate_report.py", "Auto-generated Markdown report"),
    "10": ("10_map_to_chompact.py", "Optional v1.1 CHOmpact pathway mapping"),
    "11": ("11_score_chompact_pathways.py", "Optional v1.1 CHOmpact pathway scoring"),
    "12": ("12_demand_sensitivity.py", "Optional v1.1 demand-scale sensitivity review"),
    "13": ("13_rank_pathway_biomarkers.py", "Optional v1.1 pathway biomarker ranking"),
    "14": ("14_make_chompact_figures.py", "Optional v1.1 CHOmpact figures + executive summary"),
}

FIGURE_RENAME_MAP = {
    "Fig2_rate_heatmap.png": "Fig2_exchange_rate_heatmap.png",
    "Fig3_lac_glc_ratio.png": "Fig3_lactate_glucose_phenotype.png",
    "Fig4_high_vs_low_rates.png": "Fig4_high_low_exchange_rates.png",
    "Fig5_high_vs_low_rates.png": "Fig4_high_low_exchange_rates.png",
    "Fig10_central_mab_flux_heatmap.png": "Fig5_central_mab_flux_heatmap.png",
    "Fig10B_central_mab_flux_zscore.png": "Fig6_central_mab_flux_zscore.png",
    "Fig11_central_mab_flux_delta.png": "Fig7_central_mab_flux_delta.png",
    "Fig15_pathway_scores_fva_overlap.png": "Fig8_pathway_scores_fva_overlap.png",
    "Fig21_full_fva_high_low_separation.png": "Fig9_full_fva_high_low_separation.png",
    "Fig21B_full_fva_range_delta.png": "Fig10_full_fva_range_delta.png",
    "Fig8_summary_panel.png": "Fig11_summary_panel.png",
    "Fig13_data_qc_overview.png": "SuppFig1_data_qc_overview.png",
    "Fig12_chompact_pathway_activity.png": "Fig12_measured_screening_markers.png",
    "Fig13_chompact_fva_robustness.png": "Fig13_model_emergent_fva_robustness.png",
    "Fig14_chompact_biomarker_ranking.png": "Fig14_candidate_pathway_priority_confidence.png",
}

# Legacy or retired files that should not remain in release output folders.
RETIRED_FIGURES = {
    "Fig4_fba_results.png",
    "Fig7_ko_screen.png",
    "Fig10_central_mab_flux_heatmap.png",
    "Fig10B_central_mab_flux_zscore.png",
    "Fig11_central_mab_flux_delta.png",
    "Fig12_focused_core_fva_range.png",
    "Fig12B_focused_core_fva_range_zscore.png",
    "Fig12C_focused_fva_high_low_delta.png",
    "Fig13_data_qc_overview.png",
    "Fig14_interval_pathway_scores.png",
    "Fig15_pathway_scores_fva_overlap.png",
    "Fig21_full_fva_high_low_separation.png",
    "Fig21B_full_fva_range_delta.png",
    "Fig8_summary_panel.png",
    "Fig5_high_vs_low_rates.png",
    "Fig2_rate_heatmap.png",
    "Fig3_lac_glc_ratio.png",
    "Fig4_high_vs_low_rates.png",
    "Fig12_chompact_pathway_activity.png",
    "Fig13_chompact_fva_robustness.png",
    "Fig14_chompact_biomarker_ranking.png",
}


def _standardize_figures(dataset):
    """Rename legacy figure outputs to the release figure naming scheme.

    The computational step scripts intentionally remain close to the validated
    development versions. This post-processing layer makes the user-facing output
    folder clean and deterministic.
    """
    fig_dir = os.path.join(ROOT, "results", dataset, "figures")
    if not os.path.isdir(fig_dir):
        return

    renamed = []
    for old, new in FIGURE_RENAME_MAP.items():
        old_path = os.path.join(fig_dir, old)
        new_path = os.path.join(fig_dir, new)
        if os.path.exists(old_path):
            if os.path.exists(new_path):
                os.remove(new_path)
            shutil.move(old_path, new_path)
            renamed.append((old, new))

    # Remove stale legacy files left from previous runs in the same results folder.
    for name in RETIRED_FIGURES:
        path = os.path.join(fig_dir, name)
        if os.path.exists(path):
            try:
                os.remove(path)
            except Exception:
                pass

    if renamed:
        print("  [figures] release filenames standardized:")
        for old, new in renamed:
            print(f"    {old} -> {new}")


to_run = [x.strip().zfill(2) for x in args.steps.split(",") if x.strip()]

env = os.environ.copy()
env["PYTHONIOENCODING"] = "utf-8"
env["CHO_ANALYSIS_MODE"] = args.analysis_mode
env["CHO_DEMAND_SCALE"] = str(args.demand_scale)
env["CHO_CONSTRAINT_POLICY"] = args.constraint_policy
env["CHO_FEED_VOLUME_MODE"] = args.feed_volume_mode
env["CHO_FVA_SCOPE"] = args.fva_scope
env["CHO_FVA_TARGETS"] = args.fva_targets
env["CHO_FVA_PROCESSES"] = str(args.fva_processes)
env["CHO_FVA_SOURCE"] = args.fva_source

if args.objective:
    env["CHO_OBJECTIVE"] = args.objective
if args.biomass_fraction is not None:
    env["CHO_BIOMASS_FRACTION"] = str(args.biomass_fraction)
if args.demand_fraction is not None:
    env["CHO_DEMAND_FRACTION"] = str(args.demand_fraction)

print("=" * 76)
print(f"  CHO METABOLOMICS CLEAN WORKFLOW - {args.dataset}")
print(f"  Rate Day {args.rate_days.replace(',', '->')} | Steps {','.join(to_run)} | TopN {args.top_n}")
print(f"  Objective: {args.objective or 'DM_igg_g'} | analysis_mode: {args.analysis_mode}")
print(f"  Constraint policy: {args.constraint_policy} | Feed volume mode: {args.feed_volume_mode} | Input: {args.input_format}")
print(f"  demand_scale: {args.demand_scale} | biomass_fraction: {args.biomass_fraction if args.biomass_fraction is not None else 'config default'}")
print(
    f"  FVA: scope={args.fva_scope}, targets={args.fva_targets}, "
    f"processes={args.fva_processes}, CHOmpact source={args.fva_source}"
)
print("  Product sequence module: excluded from main workflow")
print("  KO screening / direct objective-ranking figures: excluded from main workflow")
print("=" * 76)

t_total = time.time()

for step_id in to_run:
    if step_id not in STEP_MAP:
        print(f"\n  [SKIP] step {step_id}: not included in clean workflow")
        continue

    script, desc = STEP_MAP[step_id]
    path = os.path.join(STEPS_DIR, script)

    if not os.path.exists(path):
        print(f"\n  [SKIP] {script} not found")
        continue

    print(f"\n  [{step_id}] {desc}")

    cmd = [sys.executable, path, "--dataset", args.dataset]

    if step_id == "00":
        cmd += ["--input_format", args.input_format, "--analysis_mode", args.analysis_mode]

    if step_id == "01":
        cmd += [
            "--rate_days", args.rate_days,
            "--top_n", str(args.top_n),
            "--feed_volume_mode", args.feed_volume_mode,
            "--input_format", args.input_format,
        ]

    if step_id == "02":
        cmd += ["--constraint_policy", args.constraint_policy]

    if step_id in ["03", "04", "05", "06", "07", "08", "09"]:
        cmd += ["--analysis_mode", args.analysis_mode]

    if step_id in ["03", "04", "05", "06", "07"] and args.objective:
        cmd += ["--objective", args.objective]

    if step_id in ["03", "04", "05", "06", "07"] and args.biomass_fraction is not None:
        cmd += ["--biomass_fraction", str(args.biomass_fraction)]

    if step_id in ["03", "06", "07"]:
        cmd += ["--demand_scale", str(args.demand_scale)]

    if step_id in ["03", "04", "05", "07"] and args.demand_fraction is not None:
        cmd += ["--demand_fraction", str(args.demand_fraction)]

    if step_id == "06":
        cmd += [
            "--fva_scope", args.fva_scope,
            "--fva_targets", args.fva_targets,
            "--fva_processes", str(args.fva_processes),
        ]
        if args.fva_fraction is not None:
            cmd += ["--fraction", str(args.fva_fraction)]

    if step_id == "10":
        cmd += [
            "--fva_source", args.fva_source,
            "--full_fva_scope", args.fva_scope,
        ]

    log = os.path.join(LOGS_DIR, f"{step_id}_{script.replace('.py', '')}_{args.dataset}.log")

    r = subprocess.run(
        cmd,
        cwd=ROOT,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=env,
    )

    with open(log, "w", encoding="utf-8") as lf:
        lf.write(r.stdout + "\n" + r.stderr)

    keywords = [
        "saved", "obj=", "optimal", "OK", "완료",
        "Error", "Objective", "Prediction", "Explain", "diagnostic", "FVA",
        "feature", "QC", "Pathway", "sensitivity", "report", "interval",
    ]

    for line in r.stdout.split("\n"):
        if any(k in line for k in keywords):
            if line.strip():
                print(f"  {line.rstrip()}")

    if r.returncode == 0:
        print("  OK")
        if step_id in {"00", "05", "06", "07", "08"}:
            _standardize_figures(args.dataset)
    else:
        print(f"  ERROR -> {log}")
        for l in (r.stdout + r.stderr).strip().split("\n")[-12:]:
            if l.strip():
                print(f"    {l}")
        break

_standardize_figures(args.dataset)

elapsed = time.time() - t_total
print(f"\n{'=' * 76}\n  Done. {elapsed:.1f} sec")

pngs = sorted(glob.glob(os.path.join(ROOT, "results", args.dataset, "figures", "*.png")))
if pngs:
    print(f"\n  Figures ({len(pngs)} files):")
    for p in pngs:
        print(f"    {os.path.basename(p)}  ({os.path.getsize(p)//1024} KB)")
