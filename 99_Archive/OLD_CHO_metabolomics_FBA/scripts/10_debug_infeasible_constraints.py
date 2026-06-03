from __future__ import annotations

import argparse
import math
from pathlib import Path

import cobra
import pandas as pd

DEFAULT_OPEN = {
    "EX_h_e": (-1000, 1000),
    "EX_h2o_e": (-1000, 1000),
    "EX_o2_e": (-1000, 1000),
    "EX_co2_e": (-1000, 1000),
    "EX_hco3_e": (-1000, 1000),
    "EX_pi_e": (-1000, 1000),
    "EX_so4_e": (-1000, 1000),
}

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

def finite_number(x):
    try:
        return math.isfinite(float(x))
    except Exception:
        return False

def make_bounds(row, buffer=0.20):
    lb, ub = row.get("lower_bound"), row.get("upper_bound")
    if finite_number(lb) and finite_number(ub):
        lb, ub = float(lb), float(ub)
    else:
        rate = row.get("rate")
        if not finite_number(rate):
            return None
        rate = float(rate)
        lb, ub = (rate * 1.2, rate * 0.8) if rate < 0 else (rate * 0.8, rate * 1.2)
    if lb > ub:
        lb, ub = ub, lb
    return max(lb, -1000.0), min(ub, 1000.0)

def relax_bounds(lb, ub, multiplier):
    if ub <= 0:
        return lb * multiplier, 0.0
    if lb >= 0:
        return 0.0, ub * multiplier
    return lb * multiplier, ub * multiplier

def set_bounds(r, lb, ub):
    if lb > ub:
        lb, ub = ub, lb
    r.bounds = (float(lb), float(ub))

def open_defaults(model):
    for rid, bounds in DEFAULT_OPEN.items():
        if rid in model.reactions:
            try:
                set_bounds(model.reactions.get_by_id(rid), *bounds)
            except Exception:
                pass

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--condition", required=True)
    ap.add_argument("--rates-file", required=True)
    ap.add_argument("--model-file", default=None)
    ap.add_argument("--bound-buffer", type=float, default=0.20)
    ap.add_argument("--relax-multiplier", type=float, default=10.0)
    args = ap.parse_args()

    base = Path(args.base)
    model_file = pick_model(base, args.model_file)
    rates_file = Path(args.rates_file)
    if not rates_file.is_absolute():
        rates_file = base / rates_file

    outdir = base / "results" / "tables"
    outdir.mkdir(parents=True, exist_ok=True)

    model = load_model(model_file)
    rates = pd.read_csv(rates_file)
    if "condition" in rates.columns:
        rates = rates[rates["condition"].astype(str) == str(args.condition)].copy()

    print("[INFO] rows:", len(rates))
    open_defaults(model)

    ok, problem, missing = [], [], []
    for idx, row in rates.reset_index(drop=True).iterrows():
        rid = str(row.get("exchange_rxn", "")).strip()
        met = row.get("metabolite", "")
        if not rid or rid not in model.reactions:
            missing.append({"row": idx, "reaction": rid, "metabolite": met, "issue": "missing reaction in model or empty"})
            continue
        b = make_bounds(row, args.bound_buffer)
        if b is None:
            problem.append({"row": idx, "reaction": rid, "metabolite": met, "issue": "invalid bounds/rate"})
            continue

        lb, ub = b
        old = model.reactions.get_by_id(rid).bounds
        try:
            set_bounds(model.reactions.get_by_id(rid), lb, ub)
            sol = model.optimize()
            status, detail = sol.status, ""
        except Exception as e:
            status, detail = "error", str(e)

        if status == "optimal":
            ok.append({"row": idx, "reaction": rid, "metabolite": met, "lb": lb, "ub": ub, "status": status, "objective": sol.objective_value})
            continue

        try:
            set_bounds(model.reactions.get_by_id(rid), *old)
        except Exception:
            pass

        rlb, rub = relax_bounds(lb, ub, args.relax_multiplier)
        try:
            set_bounds(model.reactions.get_by_id(rid), rlb, rub)
            rsol = model.optimize()
            relaxed_status, relaxed_obj = rsol.status, rsol.objective_value
        except Exception as e:
            relaxed_status, relaxed_obj = "error", None
            detail += " relaxed_error=" + str(e)

        problem.append({
            "row": idx,
            "reaction": rid,
            "metabolite": met,
            "issue": "strict constraint caused infeasible/error",
            "lb": lb,
            "ub": ub,
            "relaxed_lb": rlb,
            "relaxed_ub": rub,
            "strict_status": status,
            "relaxed_status": relaxed_status,
            "relaxed_objective": relaxed_obj,
            "rate": row.get("rate"),
            "rate_type": row.get("rate_type"),
            "day_start": row.get("day_start"),
            "day_end": row.get("day_end"),
            "detail": detail,
        })
        if relaxed_status != "optimal":
            try:
                set_bounds(model.reactions.get_by_id(rid), *old)
            except Exception:
                pass

    pd.DataFrame(ok).to_csv(outdir / f"{args.condition}_debug_constraints_ok.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(problem).to_csv(outdir / f"{args.condition}_debug_infeasible_constraints.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(missing).to_csv(outdir / f"{args.condition}_debug_missing_reactions.csv", index=False, encoding="utf-8-sig")
    print("[OK] debug files saved in:", outdir)

if __name__ == "__main__":
    main()
