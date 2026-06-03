from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import cobra
import pandas as pd
from cobra.flux_analysis import pfba, flux_variability_analysis

KEY_REACTIONS = ["EX_glc_e", "EX_lac_L_e", "EX_gln_L_e", "EX_glu_L_e", "EX_nh4_e", "biomass_cho_prod"]

DEFAULT_OPEN = {
    "EX_h_e": (-1000, 1000),
    "EX_h2o_e": (-1000, 1000),
    "EX_o2_e": (-1000, 1000),
    "EX_co2_e": (-1000, 1000),
    "EX_hco3_e": (-1000, 1000),
    "EX_pi_e": (-1000, 1000),
    "EX_so4_e": (-1000, 1000),
}

def safe_name(x):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x)).strip("_") or "condition"

def interval_tag(d0, d1):
    return f"D{float(d0):g}_{float(d1):g}"

def finite_number(x):
    try:
        return math.isfinite(float(x))
    except Exception:
        return False

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

def load_model(path):
    if path.suffix.lower() == ".json":
        return cobra.io.load_json_model(str(path))
    return cobra.io.read_sbml_model(str(path))

def set_bounds(rxn, lb, ub):
    if lb > ub:
        lb, ub = ub, lb
    lb = max(float(lb), -1000.0)
    ub = min(float(ub), 1000.0)
    if lb > ub:
        raise ValueError(f"Invalid bounds: {lb}, {ub}")
    rxn.bounds = (lb, ub)

def strict_bounds(row, buffer):
    lb, ub = row.get("lower_bound"), row.get("upper_bound")
    if finite_number(lb) and finite_number(ub):
        lb, ub = float(lb), float(ub)
    else:
        rate = row.get("rate")
        if not finite_number(rate):
            return None
        rate = float(rate)
        lb, ub = (rate * (1 + buffer), rate * (1 - buffer)) if rate < 0 else (rate * (1 - buffer), rate * (1 + buffer))
    if lb > ub:
        lb, ub = ub, lb
    return lb, ub

def relaxed_bounds(row, multiplier):
    rate = row.get("rate")
    if not finite_number(rate):
        return None
    rate = float(rate)
    if rate < 0:
        return rate * multiplier, 0.0
    if rate > 0:
        return 0.0, rate * multiplier
    return 0.0, 0.0

def apply_constraints(model, df, mode, bound_buffer, relax_multiplier):
    applied, warnings = [], []
    for rid, (lb, ub) in DEFAULT_OPEN.items():
        if rid in model.reactions:
            try:
                set_bounds(model.reactions.get_by_id(rid), lb, ub)
            except Exception as e:
                warnings.append({"reaction": rid, "metabolite": "", "issue": "default_open_failed", "detail": str(e)})

    df = df.drop_duplicates(subset=["exchange_rxn"], keep="last").copy()
    for _, row in df.iterrows():
        rid = str(row.get("exchange_rxn", "")).strip()
        met = row.get("metabolite", "")
        if not rid:
            warnings.append({"reaction": rid, "metabolite": met, "issue": "empty_reaction_id", "detail": ""})
            continue
        if rid not in model.reactions:
            warnings.append({"reaction": rid, "metabolite": met, "issue": "reaction_not_in_model", "detail": ""})
            continue

        b = strict_bounds(row, bound_buffer) if mode == "strict" else relaxed_bounds(row, relax_multiplier)
        if b is None:
            warnings.append({"reaction": rid, "metabolite": met, "issue": "invalid_rate_or_bounds", "detail": ""})
            continue

        lb, ub = b
        try:
            rxn = model.reactions.get_by_id(rid)
            set_bounds(rxn, lb, ub)
            applied.append({
                "reaction": rid,
                "metabolite": met,
                "rate": row.get("rate", ""),
                "rate_type": row.get("rate_type", ""),
                "day_start": row.get("day_start", ""),
                "day_end": row.get("day_end", ""),
                "lower_bound": rxn.lower_bound,
                "upper_bound": rxn.upper_bound,
                "mode": mode,
            })
        except Exception as e:
            warnings.append({"reaction": rid, "metabolite": met, "issue": "set_bounds_failed", "detail": str(e), "requested_lb": lb, "requested_ub": ub})

    return pd.DataFrame(applied), pd.DataFrame(warnings)

def save_solution(sol, path):
    df = sol.fluxes.reset_index()
    df.columns = ["reaction", "flux"]
    df.to_csv(path, index=False, encoding="utf-8-sig")
    return dict(zip(df["reaction"].astype(str), df["flux"]))

def run_interval(model_file, sub, outdir, condition, d0, d1, args, mode):
    model = load_model(model_file)
    if args.objective:
        if args.objective not in model.reactions:
            raise ValueError(f"Objective not found: {args.objective}")
        model.objective = args.objective

    applied, warnings = apply_constraints(model, sub, mode, args.bound_buffer, args.relax_multiplier)
    tag = f"{safe_name(condition)}_{interval_tag(d0, d1)}_{mode}"
    applied.to_csv(outdir / f"{tag}_applied_constraints.csv", index=False, encoding="utf-8-sig")
    warnings.to_csv(outdir / f"{tag}_constraint_warnings.csv", index=False, encoding="utf-8-sig")

    flux = {}
    sol = model.optimize()
    fba_status = sol.status
    fba_obj = sol.objective_value
    if fba_status == "optimal":
        flux = save_solution(sol, outdir / f"{tag}_fba_fluxes.csv")

    pfba_status, pfba_obj = "not_run", None
    try:
        psol = pfba(model)
        pfba_status = psol.status
        pfba_obj = psol.objective_value
        flux = save_solution(psol, outdir / f"{tag}_pfba_fluxes.csv")
    except Exception as e:
        pfba_status = f"failed: {e}"

    if args.fva_mode != "none" and fba_status == "optimal":
        try:
            rxn_list = list(model.exchanges) if args.fva_mode == "exchange" else None
            fva = flux_variability_analysis(model, reaction_list=rxn_list, fraction_of_optimum=args.fraction, processes=1)
            fva.to_csv(outdir / f"{tag}_fva.csv", encoding="utf-8-sig")
        except Exception as e:
            pd.DataFrame([{"issue": "fva_failed", "detail": str(e)}]).to_csv(outdir / f"{tag}_fva_warnings.csv", index=False, encoding="utf-8-sig")

    row = {
        "condition": condition,
        "day_start": d0,
        "day_end": d1,
        "interval": f"{float(d0):g}-{float(d1):g}",
        "mode": mode,
        "tag": tag,
        "fba_status": fba_status,
        "fba_objective": fba_obj,
        "pfba_status": pfba_status,
        "pfba_objective": pfba_obj,
        "n_constraints": len(applied),
        "n_warnings": len(warnings),
    }
    for rid in KEY_REACTIONS:
        row[rid] = flux.get(rid)

    key_rows = []
    for rid in KEY_REACTIONS:
        key_rows.append({"condition": condition, "day_start": d0, "day_end": d1, "interval": f"{float(d0):g}-{float(d1):g}", "mode": mode, "reaction": rid, "flux": flux.get(rid), "tag": tag})
    return row, key_rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--rates-file", required=True)
    ap.add_argument("--model-file", default=None)
    ap.add_argument("--objective", default=None)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--mode", choices=["strict", "relaxed", "strict_then_relaxed"], default="strict_then_relaxed")
    ap.add_argument("--bound-buffer", type=float, default=0.20)
    ap.add_argument("--relax-multiplier", type=float, default=10.0)
    ap.add_argument("--fva-mode", choices=["none", "exchange", "all"], default="none")
    ap.add_argument("--fraction", type=float, default=0.90)
    args = ap.parse_args()

    base = Path(args.base)
    model_file = pick_model(base, args.model_file)
    rates_path = Path(args.rates_file)
    if not rates_path.is_absolute():
        rates_path = base / rates_path

    outdir = base / "results" / "tables" / "batch"
    outdir.mkdir(parents=True, exist_ok=True)

    rates = pd.read_csv(rates_path)
    required = {"condition", "day_start", "day_end", "exchange_rxn"}
    missing = required - set(rates.columns)
    if missing:
        raise ValueError(f"rates file missing columns: {missing}")

    rates["day_start"] = pd.to_numeric(rates["day_start"], errors="coerce")
    rates["day_end"] = pd.to_numeric(rates["day_end"], errors="coerce")
    rates = rates.dropna(subset=["condition", "day_start", "day_end", "exchange_rxn"]).copy()

    if args.conditions:
        keep = set(map(str, args.conditions))
        rates = rates[rates["condition"].astype(str).isin(keep)].copy()

    intervals = rates[["condition", "day_start", "day_end"]].drop_duplicates().sort_values(["condition", "day_start", "day_end"]).reset_index(drop=True)
    print("[INFO] model:", model_file)
    print("[INFO] rates:", rates_path)
    print("[INFO] intervals to run:", len(intervals))
    print(intervals.head(80).to_string(index=False))

    summary_rows, key_rows_all = [], []
    for _, iv in intervals.iterrows():
        cond = str(iv["condition"])
        d0 = float(iv["day_start"])
        d1 = float(iv["day_end"])
        sub = rates[(rates["condition"].astype(str) == cond) & (rates["day_start"] == d0) & (rates["day_end"] == d1)].copy()

        print(f"\n[RUN] {cond} Day {d0:g}-{d1:g} rows={len(sub)}")
        if args.mode in ["strict", "relaxed"]:
            row, key_rows = run_interval(model_file, sub, outdir, cond, d0, d1, args, args.mode)
            summary_rows.append(row)
            key_rows_all.extend(key_rows)
            print(f"  -> {args.mode} {row['fba_status']}, objective={row['fba_objective']}")
        else:
            row, key_rows = run_interval(model_file, sub, outdir, cond, d0, d1, args, "strict")
            if row["fba_status"] == "optimal":
                summary_rows.append(row)
                key_rows_all.extend(key_rows)
                print(f"  -> strict optimal, objective={row['fba_objective']}")
            else:
                print(f"  -> strict {row['fba_status']}, retry relaxed")
                row, key_rows = run_interval(model_file, sub, outdir, cond, d0, d1, args, "relaxed")
                summary_rows.append(row)
                key_rows_all.extend(key_rows)
                print(f"  -> relaxed {row['fba_status']}, objective={row['fba_objective']}")

    summary = pd.DataFrame(summary_rows)
    key_summary = pd.DataFrame(key_rows_all)
    summary_path = outdir / "batch_objective_summary.csv"
    key_path = outdir / "batch_key_exchange_summary.csv"
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")
    key_summary.to_csv(key_path, index=False, encoding="utf-8-sig")

    print("\n[OK] saved:", summary_path)
    print("[OK] saved:", key_path)
    if not summary.empty:
        print(summary[["condition", "interval", "mode", "fba_status", "fba_objective", "pfba_objective", "n_constraints", "n_warnings"]].head(100).to_string(index=False))

if __name__ == "__main__":
    main()
