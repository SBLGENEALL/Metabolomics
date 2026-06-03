from __future__ import annotations

import subprocess
import sys
from pathlib import Path

def ask(prompt, default=""):
    msg = f"{prompt} [{default}]: " if default else f"{prompt}: "
    ans = input(msg).strip()
    return ans if ans else default

def ask_choice(prompt, choices, default):
    while True:
        ans = ask(f"{prompt} ({'/'.join(choices)})", default)
        if ans in choices:
            return ans
        print("Choose one of:", ", ".join(choices))

def run_cmd(args, cwd):
    print("\n[RUN]")
    print(" ".join([f'"{x}"' if " " in str(x) else str(x) for x in args]))
    print()
    p = subprocess.run(args, cwd=str(cwd))
    if p.returncode != 0:
        print(f"[ERROR] return code {p.returncode}")
    else:
        print("[OK] done")

def default_base():
    here = Path(__file__).resolve()
    if here.parent.name.lower() == "scripts":
        return here.parent.parent
    return Path.cwd()

def main():
    base = Path(ask("Project base folder", str(default_base()))).resolve()
    scripts = base / "scripts"
    default_rates = base / "data" / "metabolomics" / "processed" / "exchange_rates_from_user_format.csv"

    while True:
        print("\n=== CHO FBA interactive runner ===")
        print("1) Convert existing raw Excel")
        print("2) Recalculate feed-corrected rates")
        print("3) Batch run all clones × intervals")
        print("4) Plot batch results")
        print("5) Run one interval FBA")
        print("q) Quit")
        task = ask("Task", "3").lower()

        if task == "q":
            return

        if task == "1":
            src = Path(ask("Raw Excel", str(base / "data" / "metabolomics" / "raw" / "CHO_raw_data.xlsx"))).resolve()
            volume = ask("Culture volume mL", "30")
            igg_unit = ask_choice("IgG/titer unit", ["mg/L", "ug/mL", "g/L", "ug/L"], "mg/L")
            cmd = [
                sys.executable, str(scripts / "09_convert_existing_excel_to_template.py"),
                "--base", str(base),
                "--input", str(src),
                "--template", str(base / "templates" / "CHO_FBA_final_input_template.xlsx"),
                "--volume", volume,
                "--igg-unit", igg_unit,
            ]
            run_cmd(cmd, base)

        elif task == "2":
            raw = Path(ask("Raw Excel with feed columns", str(base / "data" / "metabolomics" / "raw" / "CHO_raw_data.xlsx"))).resolve()
            volume = ask("Default volume mL", "30")
            cmd = [
                sys.executable, str(scripts / "13_recalculate_feed_corrected_rates.py"),
                "--base", str(base),
                "--raw-excel", str(raw),
                "--volume", volume,
            ]
            run_cmd(cmd, base)

        elif task == "3":
            rates = Path(ask("Rates CSV", str(default_rates))).resolve()
            mode = ask_choice("Mode", ["strict", "relaxed", "strict_then_relaxed"], "strict_then_relaxed")
            fva_mode = ask_choice("FVA mode", ["none", "exchange", "all"], "none")
            relax = ask("Relax multiplier", "10")
            cmd = [
                sys.executable, str(scripts / "11_batch_run_intervals.py"),
                "--base", str(base),
                "--rates-file", str(rates),
                "--mode", mode,
                "--fva-mode", fva_mode,
                "--relax-multiplier", relax,
            ]
            run_cmd(cmd, base)

        elif task == "4":
            run_cmd([sys.executable, str(scripts / "12_plot_batch_results.py"), "--base", str(base)], base)

        elif task == "5":
            condition = ask("Condition / Sample ID", "")
            d0 = ask("Day start", "5")
            d1 = ask("Day end", "7")
            rates = Path(ask("Rates CSV", str(default_rates))).resolve()
            mode = ask_choice("Mode", ["strict", "relaxed"], "strict")
            script = "03_run_fba_fva.py" if mode == "strict" else "03_run_fba_fva_relaxed.py"
            cmd = [
                sys.executable, str(scripts / script),
                "--base", str(base),
                "--condition", condition,
                "--rates-file", str(rates),
                "--day-start", d0,
                "--day-end", d1,
            ]
            if mode == "strict":
                cmd += ["--skip-fva"]
            run_cmd(cmd, base)

        else:
            print("Unknown task")

if __name__ == "__main__":
    main()
