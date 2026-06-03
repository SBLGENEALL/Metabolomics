from __future__ import annotations

import argparse
from pathlib import Path
import math
import cobra
from cobra.flux_analysis import pfba, flux_variability_analysis
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

def make_bounds_from_rate(rate, buffer=0.20):
    rate = float(rate)
    if rate < 0:
        lb, ub = rate * (1 + buffer), rate * (1 - buffer)
    else:
        lb, ub = rate * (1 - buffer), rate * (1 + buffer)
    if lb > ub:
        lb, ub = ub, lb
    return lb, ub

def sanitize_bounds(lb, ub, min_bound=-1000.0, max_bound=1000.0):
    if not finite_number(lb) or not finite_number(ub):
        return None
    lb, ub = float(lb), float(ub)
    if lb > ub:
        lb, ub = ub, lb
    lb = max(lb, min_bound)
    ub = min(ub, max_bound)
    if lb > ub:
        return None
    return lb, ub

def set_bounds_safe(reaction, lb, ub):
    reaction.bounds = (float(lb), float(ub))

def filter_rates(rates, condition, day_start=None, day_end=None, deduplicate="last"):
    df = rates.copy()
    if "condition" in df.columns:
        df = df[df["condition"].astype(str) == str(condition)].copy()
    if day_start is not None and "day_start" in df.columns:
        df = df[pd.to_numeric(df["day_start"], errors="coerce") == float(day_start)].copy()
    if day_end is not None and "day_end" in df.columns:
        df = df[pd.to_numeric(df["day_end"], errors="coerce") == float(day_end)].copy()

    if deduplicate != "none" and "exchange_rxn" in df.columns:
        if deduplicate == "last":
            df = df.drop_duplicates(subset=["exchange_rxn"], keep="last").copy()
        elif deduplicate == "mean":
            numeric_cols = [c for c in ["rate", "lower_bound", "upper_bound", "day_start", "day_end"] if c in df.columns]
            meta_cols = [c for c in df.columns if c not in numeric_cols]
            grouped_num = df.groupby("exchange_rxn", as_index=False)[numeric_cols].mean(numeric_only=True)
            grouped_meta = df.drop_duplicates(subset=["exchange_rxn"], keep="first")[meta_cols]
            df = grouped_meta.merge(grouped_num, on="exchange_rxn", how="left")
    return df

def apply_rates(model, rates, buffer=0.20):
    applied, warnings = [], []

    for rid, (lb, ub) in DEFAULT_OPEN.items():
        if rid in model.reactions:
            try:
                set_bounds_safe(model.reactions.get_by_id(rid), lb, ub)
            except Exception as e:
                warnings.append({"reaction": rid, "issue": "failed_default_open", "detail": str(e)})

    for _, row in rates.iterrows():
        rid = str(row.get("exchange_rxn", "")).strip()
        if not rid:
            warnings.append({"reaction": rid, "issue": "empty_reaction_id", "detail": ""})
            continue
        if rid not in model.reactions:
            warnings.append({"reaction": rid, "issue": "reaction_not_in_model", "detail": ""})
            continue

        lb, ub = row.get("lower_bound", None), row.get("upper_bound", None)
        if not finite_number(lb) or not finite_number(ub):
            rate = row.get("rate", None)
            if not finite_number(rate):
                warnings.append({"reaction": rid, "issue": "missing_or_invalid_rate_and_bounds", "detail": f"rate={rate}, lb={lb}, ub={ub}"})
                continue
            lb, ub = make_bounds_from_rate(float(rate), buffer=buffer)

        clean = sanitize_bounds(lb, ub)
        if clean is None:
            warnings.append({"reaction": rid, "issue": "invalid_bounds_after_sanitize", "detail": f"lb={lb}, ub={ub}"})
            continue

        lb, ub = clean
        try:
            r = model.reactions.get_by_id(rid)
            set_bounds_safe(r, lb, ub)
            applied.append({
                "reaction": rid,
                "lower_bound": r.lower_bound,
                "upper_bound": r.upper_bound,
                "metabolite": row.get("metabolite", ""),
                "rate": row.get("rate", ""),
                "rate_type": row.get("rate_type", ""),
                "day_start": row.get("day_start", ""),
                "day_end": row.get("day_end", ""),
            })
        except Exception as e:
            warnings.append({"reaction": rid, "issue": "set_bounds_failed", "detail": str(e), "requested_lb": lb, "requested_ub": ub})

    return pd.DataFrame(applied), pd.DataFrame(warnings)

def save_fluxes(sol, path):
    df = sol.fluxes.reset_index()
    df.columns = ["reaction", "flux"]
    df.to_csv(path, index=False, encoding="utf-8-sig")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--condition", required=True)
    ap.add_argument("--rates-file", required=True)
    ap.add_argument("--model-file", default=None)
    ap.add_argument("--objective", default=None)
    ap.add_argument("--skip-fva", action="store_true")
    ap.add_argument("--fva-mode", choices=["exchange", "all"], default="exchange")
    ap.add_argument("--fraction", type=float, default=0.90)
    ap.add_argument("--bound-buffer", type=float, default=0.20)
    ap.add_argument("--day-start", type=float, default=None)
    ap.add_argument("--day-end", type=float, default=None)
    ap.add_argument("--deduplicate", choices=["last", "mean", "none"], default="last")
    args = ap.parse_args()

    base = Path(args.base)
    model_file = pick_model(base, args.model_file)
    table_dir = base / "results" / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)

    rates_path = Path(args.rates_file)
    if not rates_path.is_absolute():
        rates_path = base / rates_path

    rates = pd.read_csv(rates_path)
    rates = filter_rates(rates, args.condition, args.day_start, args.day_end, args.deduplicate)
    if rates.empty:
        raise ValueError("No rate rows left after condition/day interval filtering.")

    print("[INFO] filtered rate rows:", len(rates))
    if "day_start" in rates.columns and "day_end" in rates.columns:
        print("[INFO] intervals used:")
        print(rates[["day_start", "day_end"]].drop_duplicates().to_string(index=False))

    model = load_model(model_file)
    if args.objective:
        if args.objective not in model.reactions:
            raise ValueError(f"Objective not found: {args.objective}")
        model.objective = args.objective

    applied, warnings = apply_rates(model, rates, buffer=args.bound_buffer)

    suffix = str(args.condition)
    if args.day_start is not None and args.day_end is not None:
        suffix += f"_D{args.day_start:g}_{args.day_end:g}"

    applied.to_csv(table_dir / f"{suffix}_applied_constraints.csv", index=False, encoding="utf-8-sig")
    warnings.to_csv(table_dir / f"{suffix}_constraint_warnings.csv", index=False, encoding="utf-8-sig")

    print(f"[INFO] Applied constraints: {len(applied)}")
    if len(warnings):
        print(f"[WARN] Constraint warnings: {len(warnings)}")
        print(warnings.head(20).to_string(index=False))

    sol = model.optimize()
    print(f"[FBA] {suffix}: status={sol.status}, objective={sol.objective_value}")
    if sol.status == "optimal":
        save_fluxes(sol, table_dir / f"{suffix}_fba_fluxes.csv")

    try:
        psol = pfba(model)
        print(f"[pFBA] {suffix}: objective={psol.objective_value}")
        save_fluxes(psol, table_dir / f"{suffix}_pfba_fluxes.csv")
    except Exception as e:
        print("[WARN] pFBA failed:", e)

    if not args.skip_fva:
        try:
            rxn_list = list(model.exchanges) if args.fva_mode == "exchange" else None
            fva = flux_variability_analysis(model, reaction_list=rxn_list, fraction_of_optimum=args.fraction, processes=1)
            fva.to_csv(table_dir / f"{suffix}_fva.csv", encoding="utf-8-sig")
            print(f"[FVA] {suffix}: done")
        except Exception as e:
            print("[WARN] FVA failed:", e)

if __name__ == "__main__":
    main()
