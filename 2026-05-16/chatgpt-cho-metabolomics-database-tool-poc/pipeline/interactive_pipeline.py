#!/usr/bin/env python
"""Interactive CHO multi-omics pipeline runner.

This is the main script for real experiments.

Recommended command:
python pipeline/interactive_pipeline.py

It runs step-by-step:
1. Validate and summarize experiment input
2. Calculate spent-media deltas and uptake/secretion rates
3. Generate publication-style SVG figures
4. Prepare iCHO3K exchange constraints and transcript reaction scores
5. Optionally run COBRApy FBA if cobra/libsbml/glpk are installed
"""

from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "inputs" / "CHO_actual_experiment_input.xlsx"
DEFAULT_OUTDIR = ROOT / "results" / "interactive_run"
DEFAULT_MAPPING = ROOT / "data" / "icho3k_prod_exchange_mapping_verified.csv"
DEFAULT_GPR = ROOT / "data" / "gpr_mapping_template.csv"
DEFAULT_TRANSCRIPT = ROOT / "inputs" / "transcriptomics_input.csv"


def run_cmd(args: list[str]) -> None:
    print("\n> " + " ".join(str(a) for a in args))
    result = subprocess.run(args, cwd=ROOT)
    if result.returncode != 0:
        raise SystemExit(result.returncode)


def wait(enabled: bool) -> None:
    if enabled:
        input("\nEnter를 누르면 다음 단계로 진행합니다...")


def print_csv_preview(path: Path, n: int = 8) -> None:
    if not path.exists():
        print(f"  - not found: {path}")
        return
    df = pd.read_csv(path)
    print(f"\n[{path.name}] rows={len(df)}, columns={len(df.columns)}")
    if len(df):
        print(df.head(n).to_string(index=False))


def validate_input(input_path: Path) -> None:
    if not input_path.exists():
        template = ROOT / "templates" / "CHO_actual_experiment_input_template.xlsx"
        raise SystemExit(
            f"Input file not found: {input_path}\n"
            f"Copy the template first:\n"
            f"  copy {template} {input_path}\n"
            "The inputs/ folder is git-ignored so real experiment data stays local."
        )
    if input_path.suffix.lower() in {".xlsx", ".xls"}:
        xl = pd.ExcelFile(input_path)
        sheet = "Experiment_Input" if "Experiment_Input" in xl.sheet_names else xl.sheet_names[0]
        df = pd.read_excel(input_path, sheet_name=sheet)
    else:
        df = pd.read_csv(input_path)
    required = {"sample_id", "passage_or_clone", "producer_group", "day", "replicate"}
    missing = required.difference(df.columns)
    if missing:
        raise SystemExit(f"Missing required input columns: {sorted(missing)}")
    metabolite_cols = [
        c for c in df.columns
        if c not in {
            "sample_id", "passage_or_clone", "producer_group", "day", "replicate",
            "viable_cell_density_1e6_mL", "viability_pct", "titer_mg_L",
            "sample_type", "batch_id", "notes",
        }
    ]
    print(f"Input: {input_path}")
    print(f"Samples: {len(df)}")
    print(f"Conditions: {', '.join(map(str, sorted(df['producer_group'].dropna().unique())))}")
    print(f"Metabolite columns detected: {len(metabolite_cols)}")
    print(", ".join(metabolite_cols[:40]))


def cobra_available() -> bool:
    return importlib.util.find_spec("cobra") is not None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--non-interactive", action="store_true")
    parser.add_argument("--skip-fba", action="store_true")
    args = parser.parse_args()

    interactive = not args.non_interactive
    args.outdir.mkdir(parents=True, exist_ok=True)

    print("\nCHO iCHO3K multi-omics interactive pipeline")
    print("=" * 48)

    print("\nSTEP 1. Validate experiment input")
    validate_input(args.input)
    wait(interactive)

    print("\nSTEP 2. Calculate spent-media deltas and uptake/secretion rates")
    run_cmd([
        sys.executable,
        "pipeline/spent_media_minimal_pipeline.py",
        "--input",
        str(args.input),
        "--outdir",
        str(args.outdir / "spent_media"),
    ])
    print_csv_preview(args.outdir / "spent_media" / "04_interval_rates.csv")
    wait(interactive)

    print("\nSTEP 3. Generate publication-style figures")
    run_cmd([
        sys.executable,
        "pipeline/make_spent_media_figures.py",
        "--outdir",
        str(args.outdir / "spent_media"),
        "--figdir",
        str(args.outdir / "figures"),
    ])
    print(f"Figures written to: {args.outdir / 'figures'}")
    wait(interactive)

    print("\nSTEP 4. Prepare iCHO3K model input tables")
    transcript_args = []
    if DEFAULT_TRANSCRIPT.exists():
        transcript_args = ["--transcriptomics", str(DEFAULT_TRANSCRIPT), "--gpr-map", str(DEFAULT_GPR)]
    run_cmd([
        sys.executable,
        "pipeline/prepare_multiomics_model_inputs.py",
        "--rates",
        str(args.outdir / "spent_media" / "04_interval_rates.csv"),
        "--exchange-map",
        str(DEFAULT_MAPPING),
        "--outdir",
        str(args.outdir / "icho3k_inputs"),
        *transcript_args,
    ])
    print_csv_preview(args.outdir / "icho3k_inputs" / "icho_exchange_constraints.csv")
    wait(interactive)

    print("\nSTEP 5. Optional COBRApy FBA")
    if args.skip_fba:
        print("Skipped by --skip-fba.")
    elif not cobra_available():
        print("COBRApy is not available in this Python environment.")
        print("After packaging cobra/libsbml/glpk, rerun this script or run:")
        print("python pipeline/run_icho3k_cobra_fba.py --constraints results/interactive_run/icho3k_inputs/icho_exchange_constraints.csv --outdir results/interactive_run/fba")
    else:
        run_cmd([
            sys.executable,
            "pipeline/run_icho3k_cobra_fba.py",
            "--constraints",
            str(args.outdir / "icho3k_inputs" / "icho_exchange_constraints.csv"),
            "--outdir",
            str(args.outdir / "fba"),
        ])
        print_csv_preview(args.outdir / "fba" / "fba_objective_results.csv")

    print("\nDone.")
    print(f"Results folder: {args.outdir}")


if __name__ == "__main__":
    main()
