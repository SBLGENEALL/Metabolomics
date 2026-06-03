"""
run_pipeline.py — CHO_METABOLOMICS final clean workflow

Final scientific scope
----------------------
This pipeline is not a stand-alone IgG titer predictor. It uses the iCHO3K
production model as a mechanistic layer to compare clone/group flux states and
FVA ranges from measured exchange-rate constraints.

Main workflow:
 17  Raw data and rate-calculation QC
  1  Load culture data and calculate feed-corrected exchange rates
  2  Map measured metabolites to iCHO3K exchange reactions
  3  Run FBA/pFBA in prediction-diagnostic and measured-demand explanation modes
 10  Export Escher reaction-data JSON
 14  Export central-metabolism + mAb-pathway flux/FVA reports and figures
 16  Generate focused Escher maps and focused flux/FVA JSON
 18  Interval-wise focused FBA/FVA
 19  Pathway scores, FVA overlap/separation, and constraint sensitivity
 20  Auto-generate run summary report
 21  Workstation-scale internal/full FVA
  6  Generate summary figures without direct objective-ranking or KO panels

Recommended workstation command:
  python run_pipeline.py --dataset practice_20aa --input_format tsv --rate_days 7,10 \
    --analysis_mode both --constraint_policy production_relaxed \
    --feed_volume_mode interval --demand_scale auto \
    --steps 17,1,2,3,10,14,16,21,19,20,6 \
    --fva_scope internal --fva_targets group_avg --fva_processes 16
"""
import argparse
import glob
import os
import subprocess
import sys
import time

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa", choices=["sowa2020", "practice_20aa", "own_experiment"])
parser.add_argument("--steps", default="17,1,2,3,10,14,16,18,19,20,6", help="Default fast workflow. Add 21 for workstation-scale internal/full FVA.")
parser.add_argument("--input_format", default="auto", choices=["auto", "xlsx", "tsv"], help="Use TSV on offline Linux/workstation to avoid Excel/DRM issues")
parser.add_argument("--rate_days", default="7,10", help="Rate interval, e.g. 7,10 or 10,14")
parser.add_argument("--top_n", type=int, default=3, help="High/Low producer clone count if no Group column exists")
parser.add_argument("--objective", default=None, help="Production objective reaction id; default DM_igg_g")
parser.add_argument("--analysis_mode", default="both", choices=["predict", "explain", "both"], help="predict=no IgG input; explain=measured IgG-demand pFBA; both=run both")
parser.add_argument("--biomass_fraction", type=float, default=None, help="Optional biomass lower bound fraction. Default uses config.")
parser.add_argument("--demand_fraction", type=float, default=None, help="Legacy measured-demand fraction; normally unused when demand_scale is used.")
parser.add_argument("--demand_scale", default="auto", help="For explain mode: auto or numeric global scale for measured qIgG demand.")
parser.add_argument("--constraint_policy", default="production_relaxed", choices=["strict", "production_relaxed"], help="production_relaxed avoids forcing secretion-only bounds from noisy positive nutrient rates.")
parser.add_argument("--feed_volume_mode", default="interval", choices=["interval", "cumulative"], help="interval=feed volume added during interval ending at d2; cumulative=use d2-d1.")
parser.add_argument("--fva_scope", default="internal", choices=["exchange", "focused", "internal", "all"], help="Step 21 FVA scope: exchange/focused/internal/all")
parser.add_argument("--fva_targets", default="group_avg", choices=["group_avg", "representative", "all_clones"], help="Step 21 expensive FVA targets")
parser.add_argument("--fva_processes", type=int, default=1, help="Step 21 FVA parallel processes. Use 8-32 on Linux workstation.")
parser.add_argument("--fva_fraction", type=float, default=None, help="Optional override for FVA fraction_of_optimum")
args = parser.parse_args()

ROOT = os.path.dirname(os.path.abspath(__file__))
STEPS_DIR = os.path.join(ROOT, "scripts", "steps")
LOGS_DIR = os.path.join(ROOT, "logs")
os.makedirs(LOGS_DIR, exist_ok=True)

STEP_MAP = {
    1: ("01_load_data.py", "Data load + feed-corrected exchange-rate calculation"),
    2: ("02_map_metabolites.py", "Metabolite-to-exchange mapping"),
    3: ("03_run_fba.py", "FBA/pFBA prediction diagnostic + measured-demand explanation"),
    6: ("06_figures.py", "Main summary figures without objective-ranking or KO panels"),
    10: ("export_escher_flux.py", "Escher flux reaction-data JSON"),
    14: ("14_central_mab_fva_report.py", "Central metabolism + mAb pathway flux/FVA report"),
    16: ("16_focused_escher_maps.py", "Focused Escher maps + focused FVA JSON"),
    17: ("17_data_qc.py", "Raw data / rate QC and sanity checks"),
    18: ("18_interval_fba_fva.py", "Interval-wise focused FBA/FVA"),
    19: ("19_pathway_scores_sensitivity.py", "Pathway scores + FVA overlap + constraint sensitivity"),
    20: ("20_generate_report.py", "Auto-generated Markdown report"),
    21: ("21_full_internal_fva.py", "Workstation-scale internal/all-reaction FVA"),
}

to_run = [int(x.strip()) for x in args.steps.split(",") if x.strip()]
env = os.environ.copy()
env["PYTHONIOENCODING"] = "utf-8"
env["CHO_ANALYSIS_MODE"] = args.analysis_mode
env["CHO_DEMAND_SCALE"] = str(args.demand_scale)
env["CHO_CONSTRAINT_POLICY"] = args.constraint_policy
env["CHO_FEED_VOLUME_MODE"] = args.feed_volume_mode
env["CHO_FVA_SCOPE"] = args.fva_scope
env["CHO_FVA_TARGETS"] = args.fva_targets
env["CHO_FVA_PROCESSES"] = str(args.fva_processes)
if args.objective:
    env["CHO_OBJECTIVE"] = args.objective
if args.biomass_fraction is not None:
    env["CHO_BIOMASS_FRACTION"] = str(args.biomass_fraction)
if args.demand_fraction is not None:
    env["CHO_DEMAND_FRACTION"] = str(args.demand_fraction)

print("=" * 76)
print(f"  CHO_METABOLOMICS FINAL — {args.dataset}")
print(f"  Rate Day {args.rate_days.replace(',', '→')} | Steps {args.steps} | TopN {args.top_n}")
print(f"  Objective: {args.objective or 'DM_igg_g'} | analysis_mode: {args.analysis_mode}")
print(f"  Constraint policy: {args.constraint_policy} | Feed volume mode: {args.feed_volume_mode} | Input: {args.input_format}")
print(f"  demand_scale: {args.demand_scale} | biomass_fraction: {args.biomass_fraction if args.biomass_fraction is not None else 'config default'}")
print(f"  Workstation FVA: scope={args.fva_scope}, targets={args.fva_targets}, processes={args.fva_processes}")
print("  Product sequence module: excluded from final main workflow")
print("  KO screening / direct objective-ranking figures: excluded from final main workflow")
print("=" * 76)

t_total = time.time()
for num in to_run:
    if num not in STEP_MAP:
        print(f"\n  [SKIP] step {num}: not included in final clean workflow")
        continue
    script, desc = STEP_MAP[num]
    path = os.path.join(STEPS_DIR, script)
    if not os.path.exists(path):
        print(f"\n  [SKIP] {script} not found")
        continue
    print(f"\n  [{num}] {desc}")

    cmd = [sys.executable, path, "--dataset", args.dataset]
    if num == 1:
        cmd += ["--rate_days", args.rate_days, "--top_n", str(args.top_n), "--feed_volume_mode", args.feed_volume_mode, "--input_format", args.input_format]
    if num == 18:
        cmd += ["--feed_volume_mode", args.feed_volume_mode]
    if num == 17:
        cmd += ["--input_format", args.input_format]
    if num == 2:
        cmd += ["--constraint_policy", args.constraint_policy]
    if num in [3, 6, 10, 14, 16, 17, 18, 19, 20, 21]:
        cmd += ["--analysis_mode", args.analysis_mode]
    if num in [3, 10, 14, 18, 19, 21] and args.objective:
        cmd += ["--objective", args.objective]
    if num in [3, 10, 14, 18, 19, 21] and args.biomass_fraction is not None:
        cmd += ["--biomass_fraction", str(args.biomass_fraction)]
    if num in [3, 18, 19, 21]:
        cmd += ["--demand_scale", str(args.demand_scale)]
    if num in [3, 10, 14, 18, 19] and args.demand_fraction is not None:
        cmd += ["--demand_fraction", str(args.demand_fraction)]
    if num == 21:
        cmd += ["--fva_scope", args.fva_scope, "--fva_targets", args.fva_targets, "--fva_processes", str(args.fva_processes)]
        if args.fva_fraction is not None:
            cmd += ["--fraction", str(args.fva_fraction)]

    log = os.path.join(LOGS_DIR, f"{script.replace('.py', '')}_{args.dataset}.log")
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, encoding="utf-8", errors="replace", env=env)
    with open(log, "w", encoding="utf-8") as lf:
        lf.write(r.stdout + "\n" + r.stderr)

    for line in r.stdout.split("\n"):
        if any(k in line for k in ["saved", "obj=", "optimal", "OK", "✔", "✘", "★", "완료", "Error", "Objective", "Prediction", "Explain", "diagnostic", "FVA", "feature", "QC", "Pathway", "sensitivity", "report", "interval"]):
            if line.strip():
                print(f"  {line.rstrip()}")

    if r.returncode == 0:
        print("  ✔ 완료")
    else:
        print(f"  ✘ 오류 → {log}")
        for l in (r.stdout + r.stderr).strip().split("\n")[-12:]:
            if l.strip():
                print(f"    {l}")
        break

elapsed = time.time() - t_total
print(f"\n{'=' * 76}\n  완료! {elapsed:.1f}초")
pngs = sorted(glob.glob(os.path.join(ROOT, "results", args.dataset, "figures", "*.png")))
if pngs:
    print(f"\n  Figures ({len(pngs)} files):")
    for p in pngs:
        print(f"    {os.path.basename(p)}  ({os.path.getsize(p)//1024} KB)")
