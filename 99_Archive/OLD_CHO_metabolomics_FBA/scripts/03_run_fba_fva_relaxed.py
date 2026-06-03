from __future__ import annotations

import argparse
import math
from pathlib import Path

import cobra
from cobra.flux_analysis import pfba
import pandas as pd

def load_model(model_file):
    p = Path(model_file)
    if p.suffix.lower() == ".json":
        return cobra.io.load_json_model(str(p))
    return cobra.io.read_sbml_model(str(p))

def pick_model(base, model_file):
    if model_file:
        p = Path(model_file)
        return p if p.is_absolute() else base / p
    for p in [
        base / "model" / "iCHO3K-main" / "iCHO3K" / "Model" / "iCHO3K_cho_prod_generic_unblocked.json",
        base / "model" / "iCHO3K-main" / "iCHO3K" / "Model" / "iCHO3K_cho_prod_generic_unblocked.xml",
    ]:
        if p.exists():
            return p
    raise FileNotFoundError("Model file not found. Use --model-file")

def finite(x):
    try:
        return math.isfinite(float(x))
    except Exception:
        return False

def relaxed_bounds_from_rate(rate, multiplier):
    rate = float(rate)
    m = float(multiplier)
    if rate < 0:
        return rate * m, 0.0
    if rate > 0:
        return 0.0, rate * m
    return 0.0, 0.0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--condition", required=True)
    ap.add_argument("--rates-file", required=True)
    ap.add_argument("--model-file", default=None)
    ap.add_argument("--relax-multiplier", type=float, default=10.0)
    ap.add_argument("--day-start", type=float, default=None)
    ap.add_argument("--day-end", type=float, default=None)
    args = ap.parse_args()

    base = Path(args.base)
    model_file = pick_model(base, args.model_file)
    rates_path = Path(args.rates_file)
    if not rates_path.is_absolute():
        rates_path = base / rates_path

    outdir = base / "results" / "tables"
    outdir.mkdir(parents=True, exist_ok=True)

    model = load_model(model_file)
    rates = pd.read_csv(rates_path)

    if "condition" in rates.columns:
        rates = rates[rates["condition"].astype(str) == str(args.condition)].copy()
    if args.day_start is not None and "day_start" in rates.columns:
        rates = rates[pd.to_numeric(rates["day_start"], errors="coerce") == float(args.day_start)].copy()
    if args.day_end is not None and "day_end" in rates.columns:
        rates = rates[pd.to_numeric(rates["day_end"], errors="coerce") == float(args.day_end)].copy()
    if "exchange_rxn" in rates.columns:
        rates = rates.drop_duplicates(subset=["exchange_rxn"], keep="last").copy()

    if rates.empty:
        raise ValueError("No rate rows left after condition/day interval filtering.")

    print("[INFO] filtered rate rows:", len(rates))
    print(rates[[c for c in ["metabolite", "exchange_rxn", "rate", "day_start", "day_end"] if c in rates.columns]].to_string(index=False))

    applied = []
    for _, row in rates.iterrows():
        rid = str(row.get("exchange_rxn", "")).strip()
        if rid not in model.reactions:
            continue
        rate = row.get("rate", None)
        if not finite(rate):
            continue
        lb, ub = relaxed_bounds_from_rate(float(rate), args.relax_multiplier)
        if lb > ub:
            lb, ub = ub, lb
        lb, ub = max(lb, -1000), min(ub, 1000)
        model.reactions.get_by_id(rid).bounds = (lb, ub)
        applied.append({"reaction": rid, "metabolite": row.get("metabolite", ""), "rate": rate, "lower_bound": lb, "upper_bound": ub})

    suffix = f"{args.condition}"
    if args.day_start is not None and args.day_end is not None:
        suffix += f"_D{args.day_start:g}_{args.day_end:g}"

    pd.DataFrame(applied).to_csv(outdir / f"{suffix}_relaxed_applied_constraints.csv", index=False, encoding="utf-8-sig")

    sol = model.optimize()
    print(f"[RELAXED FBA] {suffix}: status={sol.status}, objective={sol.objective_value}")
    if sol.status == "optimal":
        df = sol.fluxes.reset_index()
        df.columns = ["reaction", "flux"]
        df.to_csv(outdir / f"{suffix}_relaxed_fba_fluxes.csv", index=False, encoding="utf-8-sig")
        try:
            psol = pfba(model)
            pdf = psol.fluxes.reset_index()
            pdf.columns = ["reaction", "flux"]
            pdf.to_csv(outdir / f"{suffix}_relaxed_pfba_fluxes.csv", index=False, encoding="utf-8-sig")
            print(f"[RELAXED pFBA] objective={psol.objective_value}")
        except Exception as e:
            print("[WARN] relaxed pFBA failed:", e)

if __name__ == "__main__":
    main()
