r"""
09_convert_existing_excel_to_template.py

Robust converter for user CHO raw Excel.

Supports columns:
DAY / DAY1, Sample ID, Gln, Glu, Gluc, Lac, NH4+, Na+, K+, Ca++, pH, PO2, PCO2, Osm,
Total Desity, Viable Density, Viability, Average Live Diameter, IgG,
and optional feeding columns.

Output:
data/metabolomics/processed/converted_from_user_format.xlsx
data/metabolomics/processed/converted_process_data.csv
data/metabolomics/processed/converted_extracellular_metabolomics.csv
data/metabolomics/processed/exchange_rates_from_user_format.csv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook, Workbook

COLUMN_ALIASES = {
    "day": ["day", "culture day", "sampling day", "sample day", "timepoint", "time point", "day1"],
    "flask_id": ["sample id", "sample_id", "flask id", "flask_id", "flaskid", "flask", "condition", "clone", "clone id", "id"],
    "gln": ["gln", "glutamine", "l-glutamine"],
    "glu": ["glu", "glutamate", "l-glutamate", "glutamic acid"],
    "gluc": ["gluc", "glucose", "glc"],
    "lac": ["lac", "lactate", "lactic acid"],
    "nh4": ["nh4+", "nh4", "ammonia", "ammonium"],
    "na": ["na+", "na", "sodium"],
    "k": ["k+", "k", "potassium"],
    "ca": ["ca++", "ca2+", "ca", "calcium"],
    "ph": ["ph", "pH"],
    "po2": ["po2", "pO2", "do", "dissolved oxygen"],
    "pco2": ["pco2", "pCO2", "co2"],
    "osm": ["osm", "osmolality", "osmolarity", "osmo"],
    "total_density": ["total density", "total desity", "total_density", "total_desity", "total cell density", "tcd"],
    "viable_density": ["viable density", "viable_density", "vcd", "viable cell density", "viable cells"],
    "viability": ["viability", "viability %", "viability_pct", "viab", "viab%"],
    "avg_live_diameter": ["average liv diameter", "average live diameter", "avg live diameter", "cell diameter"],
    "igg": ["igg", "titer", "titre", "mab", "antibody", "product"],
}

METABOLITE_MAP = {
    "gln": ("Glutamine", "mM", "EX_gln_L_e", True),
    "glu": ("Glutamate", "mM", "EX_glu_L_e", True),
    "gluc": ("Glucose", "mM", "EX_glc_e", True),
    "lac": ("Lactate", "mM", "EX_lac_L_e", True),
    "nh4": ("Ammonia", "mM", "EX_nh4_e", True),
    "na": ("Sodium", "mM", "", False),
    "k": ("Potassium", "mM", "", False),
    "ca": ("Calcium", "mM", "", False),
}

def norm_col(x):
    s = str(x).strip().lower()
    s = s.replace("\n", " ").replace("_", " ").replace("-", " ")
    s = re.sub(r"\s+", " ", s)
    return s.strip()

def compact_col(x):
    return re.sub(r"[^a-z0-9]+", "", norm_col(x))

def parse_day_value(x):
    if pd.isna(x):
        return None
    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x)
    s = str(x).strip()
    if not s:
        return None
    try:
        return float(s)
    except Exception:
        pass
    m = re.search(r"(?:day|d)\s*[_-]?\s*([0-9]+(?:\.[0-9]+)?)", s, flags=re.I)
    if m:
        return float(m.group(1))
    nums = re.findall(r"[0-9]+(?:\.[0-9]+)?", s)
    if len(nums) == 1 and re.search(r"day|^d", s, flags=re.I):
        return float(nums[0])
    return None

def canonicalize_columns(columns):
    norm_to_actual = {norm_col(c): c for c in columns}
    comp_to_actual = {compact_col(c): c for c in columns}
    found = {}
    for canon, aliases in COLUMN_ALIASES.items():
        for alias in aliases:
            if norm_col(alias) in norm_to_actual:
                found[canon] = norm_to_actual[norm_col(alias)]
                break
            if compact_col(alias) in comp_to_actual:
                found[canon] = comp_to_actual[compact_col(alias)]
                break
    if "day" not in found:
        for col in columns:
            if parse_day_value(col) is not None and re.search(r"day|^d\s*\d", str(col), flags=re.I):
                found["day"] = col
                break
    return found

def count_known_headers(values):
    vals = {compact_col(v) for v in values}
    count = 0
    for aliases in COLUMN_ALIASES.values():
        if any(compact_col(a) in vals for a in aliases):
            count += 1
    return count

def read_sheet_with_header(xlsx, sheet_name):
    raw = pd.read_excel(xlsx, sheet_name=sheet_name, header=None)
    header_row, best_score = 0, -1
    for i in range(min(80, len(raw))):
        score = count_known_headers(raw.iloc[i].tolist())
        if score > best_score:
            best_score, header_row = score, i
        if score >= 4:
            break
    df = pd.read_excel(xlsx, sheet_name=sheet_name, header=header_row)
    df = df.dropna(axis=0, how="all").dropna(axis=1, how="all")
    df.columns = [str(c).strip() for c in df.columns]
    return df, header_row

def guess_day_column(df):
    best_col, best_hits = None, 0
    for col in df.columns:
        hits = df[col].dropna().head(80).apply(lambda v: parse_day_value(v) is not None).sum()
        if hits > best_hits:
            best_col, best_hits = col, hits
    return best_col if best_hits >= 1 else None

def fallback_flask_column(df, cmap):
    known = set(cmap.values())
    for col in df.columns:
        if col not in known:
            non_na = df[col].dropna().astype(str).str.strip()
            if len(non_na[non_na != ""]) > 0:
                return col
    return None

def convert_igg_to_g_l(values, unit):
    unit = unit.lower().replace(" ", "")
    v = pd.to_numeric(values, errors="coerce")
    if unit in ["mg/l", "ug/ml", "µg/ml"]:
        return v / 1000.0
    if unit in ["g/l", "gpl"]:
        return v
    if unit in ["ug/l", "µg/l"]:
        return v / 1_000_000.0
    return v / 1000.0

def get_cell(row, cmap, key):
    col = cmap.get(key)
    if not col:
        return np.nan
    return row.get(col, np.nan)

def build_long_tables(source_xlsx, default_day, volume, igg_unit, platform, outdir):
    xl = pd.ExcelFile(source_xlsx)
    process_rows, metab_rows, debug_rows, preview_rows = [], [], [], []

    for sheet in xl.sheet_names:
        df, header_row = read_sheet_with_header(source_xlsx, sheet)
        if df.empty:
            continue
        cmap = canonicalize_columns(df.columns)

        for c in df.columns:
            debug_rows.append({"sheet": sheet, "header_row_0based": header_row, "column": c, "detected_as": ";".join([k for k, v in cmap.items() if v == c])})
        prev = df.head(10).copy()
        prev.insert(0, "_sheet", sheet)
        preview_rows.append(prev)

        day_col = cmap.get("day")
        if day_col is None:
            guessed = guess_day_column(df)
            if guessed:
                cmap["day"] = guessed
                day_col = guessed
                print(f"[INFO] guessed day column in sheet '{sheet}': {day_col}")

        if "flask_id" not in cmap:
            fallback = fallback_flask_column(df, cmap)
            if fallback:
                cmap["flask_id"] = fallback
                print(f"[WARN] Using '{fallback}' as Sample ID / condition.")
            else:
                continue

        if day_col is not None:
            parsed_days = df[day_col].apply(parse_day_value).ffill()
        else:
            sheet_day = parse_day_value(sheet)
            if sheet_day is None:
                sheet_day = default_day
            if sheet_day is None:
                raise ValueError(f"Cannot determine day for sheet '{sheet}'.")
            parsed_days = pd.Series([sheet_day] * len(df), index=df.index)

        sid_col = cmap["flask_id"]
        for idx, row in df.iterrows():
            sid = str(row.get(sid_col, "")).strip()
            if not sid or sid.lower() == "nan":
                continue
            day = parsed_days.loc[idx]
            if day is None or pd.isna(day):
                continue
            day = float(day)

            igg_raw = get_cell(row, cmap, "igg")
            igg_g_l = convert_igg_to_g_l(pd.Series([igg_raw]), igg_unit).iloc[0]

            process_rows.append({
                "condition": sid,
                "day": day,
                "time_h": day * 24.0,
                "vcd_10e6_cells_mL": pd.to_numeric(get_cell(row, cmap, "viable_density"), errors="coerce"),
                "viability_pct": pd.to_numeric(get_cell(row, cmap, "viability"), errors="coerce"),
                "culture_volume_mL": volume,
                "sample_removed_mL": np.nan,
                "ivcd_10e6_cell_day_mL": np.nan,
                "titer_g_L": igg_g_l,
                "qP_pg_cell_day": np.nan,
                "note": f"converted from sheet={sheet}; header_row={header_row}",
                "total_density_10e6_cells_mL": pd.to_numeric(get_cell(row, cmap, "total_density"), errors="coerce"),
                "average_live_diameter": pd.to_numeric(get_cell(row, cmap, "avg_live_diameter"), errors="coerce"),
                "pH": pd.to_numeric(get_cell(row, cmap, "ph"), errors="coerce"),
                "pO2": pd.to_numeric(get_cell(row, cmap, "po2"), errors="coerce"),
                "pCO2": pd.to_numeric(get_cell(row, cmap, "pco2"), errors="coerce"),
                "osmolality_mOsm_kg": pd.to_numeric(get_cell(row, cmap, "osm"), errors="coerce"),
                "Na_mM": pd.to_numeric(get_cell(row, cmap, "na"), errors="coerce"),
                "K_mM": pd.to_numeric(get_cell(row, cmap, "k"), errors="coerce"),
                "Ca_mM": pd.to_numeric(get_cell(row, cmap, "ca"), errors="coerce"),
            })

            for key, (met_name, unit, ex_rxn, fba_used) in METABOLITE_MAP.items():
                if key not in cmap:
                    continue
                val = pd.to_numeric(get_cell(row, cmap, key), errors="coerce")
                if pd.isna(val):
                    continue
                metab_rows.append({
                    "condition": sid,
                    "day": day,
                    "time_h": day * 24.0,
                    "replicate": "R1",
                    "metabolite": met_name,
                    "concentration": float(val),
                    "unit": unit,
                    "platform": platform,
                    "sample_type": "spent_media",
                    "qc_flag": "OK",
                    "note": f"converted from {sheet}; fba_used={fba_used}",
                })

    if debug_rows:
        pd.DataFrame(debug_rows).to_csv(outdir / "converter_debug_detected_columns.csv", index=False, encoding="utf-8-sig")
    if preview_rows:
        pd.concat(preview_rows, ignore_index=True, sort=False).to_csv(outdir / "converter_debug_sheet_preview.csv", index=False, encoding="utf-8-sig")

    process = pd.DataFrame(process_rows)
    metab = pd.DataFrame(metab_rows)

    if process.empty:
        raise ValueError("No process rows were converted. Check converter_debug files.")

    process = process.sort_values(["condition", "day"]).reset_index(drop=True)
    metab = metab.sort_values(["condition", "day", "metabolite"]).reset_index(drop=True)

    process["ivcd_10e6_cell_day_mL"] = np.nan
    process["qP_pg_cell_day"] = np.nan

    for cond, sub in process.groupby("condition", sort=False):
        idxs = sorted(list(sub.index), key=lambda i: process.loc[i, "day"])
        cum_ivcd = 0.0
        prev_i = None
        for i in idxs:
            if prev_i is None:
                process.loc[i, "ivcd_10e6_cell_day_mL"] = 0.0
                prev_i = i
                continue
            d0, d1 = process.loc[prev_i, "day"], process.loc[i, "day"]
            v0, v1 = process.loc[prev_i, "vcd_10e6_cells_mL"], process.loc[i, "vcd_10e6_cells_mL"]
            interval_ivcd = np.nan
            if pd.notna(v0) and pd.notna(v1) and d1 > d0:
                interval_ivcd = ((float(v0) + float(v1)) / 2.0) * (float(d1) - float(d0))
                cum_ivcd += interval_ivcd
            process.loc[i, "ivcd_10e6_cell_day_mL"] = cum_ivcd
            t0, t1 = process.loc[prev_i, "titer_g_L"], process.loc[i, "titer_g_L"]
            if pd.notna(t0) and pd.notna(t1) and pd.notna(interval_ivcd) and interval_ivcd > 0:
                process.loc[i, "qP_pg_cell_day"] = (float(t1) - float(t0)) * 1000.0 / interval_ivcd
            prev_i = i

    return process, metab

def build_exchange_rates(process, metab, rate_mode, bound_buffer):
    met_to_ex = {v[0]: v[2] for k, v in METABOLITE_MAP.items() if v[3]}
    rows = []

    avg = metab[metab["metabolite"].isin(met_to_ex)].copy()
    avg["concentration"] = pd.to_numeric(avg["concentration"], errors="coerce")
    avg = avg.dropna(subset=["concentration"])
    avg = avg.groupby(["condition", "day", "metabolite"], as_index=False)["concentration"].mean()
    proc_key = process.set_index(["condition", "day"])

    for cond in sorted(avg["condition"].astype(str).unique()):
        sub = avg[avg["condition"].astype(str) == cond]
        days = sorted(sub["day"].unique())
        for d0, d1 in zip(days[:-1], days[1:]):
            s0 = sub[sub["day"] == d0].set_index("metabolite")["concentration"]
            s1 = sub[sub["day"] == d1].set_index("metabolite")["concentration"]
            try:
                v0 = proc_key.loc[(cond, d0), "vcd_10e6_cells_mL"]
                v1 = proc_key.loc[(cond, d1), "vcd_10e6_cells_mL"]
                interval_ivcd = ((float(v0) + float(v1)) / 2.0) * (float(d1) - float(d0))
            except Exception:
                interval_ivcd = np.nan

            for met in sorted(set(s0.index) & set(s1.index)):
                ex = met_to_ex.get(met)
                delta = float(s1.loc[met] - s0.loc[met])
                delta_day = float(d1 - d0)
                if delta_day <= 0:
                    continue
                if rate_mode == "ivcd":
                    if pd.isna(interval_ivcd) or interval_ivcd <= 0:
                        continue
                    rate = delta / interval_ivcd
                    unit = "mM_per_10e6_cell_day_mL"
                else:
                    rate = delta / delta_day
                    unit = "mM_per_day"
                if rate < 0:
                    lb, ub = rate * (1 + bound_buffer), rate * (1 - bound_buffer)
                else:
                    lb, ub = rate * (1 - bound_buffer), rate * (1 + bound_buffer)
                if lb > ub:
                    lb, ub = ub, lb
                rows.append({
                    "condition": cond,
                    "metabolite": met,
                    "exchange_rxn": ex,
                    "rate": rate,
                    "rate_type": "uptake" if rate < 0 else "secretion",
                    "unit": unit,
                    "day_start": d0,
                    "day_end": d1,
                    "lower_bound": lb,
                    "upper_bound": ub,
                    "source": "converted_user_excel",
                    "note": f"rate_mode={rate_mode}; interval_ivcd={interval_ivcd}",
                })
    return pd.DataFrame(rows)

def clear_and_write(ws, df):
    ws.delete_rows(1, ws.max_row)
    ws.append(list(df.columns))
    for row in df.replace({np.nan: None}).itertuples(index=False, name=None):
        ws.append(list(row))

def write_to_template(template, output, process, metab):
    if template.exists():
        wb = load_workbook(template)
    else:
        wb = Workbook()
    for sheet_name in ["Process_Data", "Extracellular_Metabolomics"]:
        if sheet_name not in wb.sheetnames:
            wb.create_sheet(sheet_name)
    clear_and_write(wb["Process_Data"], process)
    clear_and_write(wb["Extracellular_Metabolomics"], metab)

    if "Experiment_Design" not in wb.sheetnames:
        wb.create_sheet("Experiment_Design")
    conditions = sorted(process["condition"].dropna().astype(str).unique())
    design = pd.DataFrame({
        "condition": conditions,
        "clone_id": conditions,
        "producer_class": "",
        "run_id": "converted",
        "vessel": "",
        "media": "",
        "feed_strategy": "",
        "objective_note": "",
        "include_in_fba": True,
    })
    clear_and_write(wb["Experiment_Design"], design)
    output.parent.mkdir(parents=True, exist_ok=True)
    wb.save(output)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--template", default="templates/CHO_FBA_final_input_template.xlsx")
    ap.add_argument("--output", default="data/metabolomics/processed/converted_from_user_format.xlsx")
    ap.add_argument("--volume", type=float, default=30.0)
    ap.add_argument("--day", type=float, default=None)
    ap.add_argument("--igg-unit", default="mg/L", choices=["mg/L", "ug/mL", "g/L", "ug/L"])
    ap.add_argument("--platform", default="BioProfile")
    ap.add_argument("--rate-mode", choices=["ivcd", "concentration"], default="ivcd")
    ap.add_argument("--bound-buffer", type=float, default=0.20)
    args = ap.parse_args()

    base = Path(args.base)
    src = Path(args.input)
    if not src.is_absolute():
        src = base / src
    template = Path(args.template)
    if not template.is_absolute():
        template = base / template
    output = Path(args.output)
    if not output.is_absolute():
        output = base / output

    outdir = base / "data" / "metabolomics" / "processed"
    outdir.mkdir(parents=True, exist_ok=True)

    process, metab = build_long_tables(src, args.day, args.volume, args.igg_unit, args.platform, outdir)
    rates = build_exchange_rates(process, metab, args.rate_mode, args.bound_buffer)

    process.to_csv(outdir / "converted_process_data.csv", index=False, encoding="utf-8-sig")
    metab.to_csv(outdir / "converted_extracellular_metabolomics.csv", index=False, encoding="utf-8-sig")
    rates.to_csv(outdir / "exchange_rates_from_user_format.csv", index=False, encoding="utf-8-sig")
    write_to_template(template, output, process, metab)

    print("[OK] saved converted workbook:", output)
    print("[OK] saved:", outdir / "exchange_rates_from_user_format.csv")
    print(process.head(10).to_string(index=False))
    print(rates.head(20).to_string(index=False))

if __name__ == "__main__":
    main()
