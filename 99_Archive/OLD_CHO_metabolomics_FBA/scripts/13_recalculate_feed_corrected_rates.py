r"""
13_recalculate_feed_corrected_rates.py

Recalculate exchange rates with feeding correction.

Feed columns in raw Excel:
Culture Volume mL, Sample Removed mL,
Glucose Feed mL, Glucose Feed Conc mM,
Feed4 mL, CellBoost mL,
Feed4 Glucose mM, Feed4 Gln mM, Feed4 Glu mM, Feed4 Lac mM, Feed4 NH4 mM,
CellBoost Glucose mM, CellBoost Gln mM, CellBoost Glu mM, CellBoost Lac mM, CellBoost NH4 mM

Meaning:
Values in a Day row are additions that occurred since the previous sampling
for that same Sample ID, up to that row's day.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

MET_TO_EX = {
    "Glucose": "EX_glc_e",
    "Lactate": "EX_lac_L_e",
    "Glutamine": "EX_gln_L_e",
    "Glutamate": "EX_glu_L_e",
    "Ammonia": "EX_nh4_e",
}

DAY_ALIASES = ["day", "culture day", "sample day", "sampling day", "day1"]
COND_ALIASES = ["sample id", "sample_id", "flask id", "flask_id", "condition", "clone", "id"]

FEED_ALIASES = {
    "culture_volume_mL": ["culture volume ml", "culture_volume_ml", "volume ml", "working volume ml", "volume"],
    "sample_removed_mL": ["sample removed ml", "sample_removed_ml", "sampling volume ml", "sample volume ml"],
    "glucose_feed_mL": ["glucose feed ml", "glucose_feed_ml", "gluc feed ml", "glucose shot ml"],
    "glucose_feed_conc_mM": ["glucose feed conc mm", "glucose_feed_conc_mm", "glucose stock mm", "gluc feed conc mm"],
    "feed4_mL": ["feed4 ml", "feed 4 ml", "feed4 volume ml", "feed4"],
    "cellboost_mL": ["cellboost ml", "cell boost ml", "cellboost volume ml", "cell boost"],
    "feed4_Glucose_mM": ["feed4 glucose mm", "feed4 gluc mm", "feed 4 glucose mm"],
    "feed4_Glutamine_mM": ["feed4 gln mm", "feed4 glutamine mm", "feed 4 gln mm"],
    "feed4_Glutamate_mM": ["feed4 glu mm", "feed4 glutamate mm", "feed 4 glu mm"],
    "feed4_Lactate_mM": ["feed4 lac mm", "feed4 lactate mm", "feed 4 lac mm"],
    "feed4_Ammonia_mM": ["feed4 nh4 mm", "feed4 ammonia mm", "feed4 ammonium mm"],
    "cellboost_Glucose_mM": ["cellboost glucose mm", "cell boost glucose mm", "cellboost gluc mm"],
    "cellboost_Glutamine_mM": ["cellboost gln mm", "cell boost gln mm", "cellboost glutamine mm"],
    "cellboost_Glutamate_mM": ["cellboost glu mm", "cell boost glu mm", "cellboost glutamate mm"],
    "cellboost_Lactate_mM": ["cellboost lac mm", "cell boost lac mm", "cellboost lactate mm"],
    "cellboost_Ammonia_mM": ["cellboost nh4 mm", "cell boost nh4 mm", "cellboost ammonia mm"],
}

def norm(x):
    s = str(x).strip().lower().replace("_", " ").replace("-", " ")
    return re.sub(r"\s+", " ", s).strip()

def compact(x):
    return re.sub(r"[^a-z0-9]+", "", norm(x))

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
    return None

def find_col(columns, aliases):
    nmap = {norm(c): c for c in columns}
    cmap = {compact(c): c for c in columns}
    for a in aliases:
        if norm(a) in nmap:
            return nmap[norm(a)]
        if compact(a) in cmap:
            return cmap[compact(a)]
    return None

def detect_header_row(xlsx, sheet_name):
    raw = pd.read_excel(xlsx, sheet_name=sheet_name, header=None)
    aliases_all = DAY_ALIASES + COND_ALIASES + ["gln", "glu", "gluc", "lac", "nh4+", "viable density", "igg"]
    for aliases in FEED_ALIASES.values():
        aliases_all += aliases
    aliases_comp = {compact(a) for a in aliases_all}
    best_i, best_score = 0, -1
    for i in range(min(80, len(raw))):
        vals = [compact(v) for v in raw.iloc[i].tolist()]
        score = sum(1 for v in vals if v in aliases_comp)
        if score > best_score:
            best_i, best_score = i, score
        if score >= 4:
            return i
    return best_i

def read_raw_feed_events(raw_excel, default_volume):
    xl = pd.ExcelFile(raw_excel)
    rows, debug = [], []
    for sheet in xl.sheet_names:
        header = detect_header_row(raw_excel, sheet)
        df = pd.read_excel(raw_excel, sheet_name=sheet, header=header)
        df = df.dropna(axis=0, how="all").dropna(axis=1, how="all")
        df.columns = [str(c).strip() for c in df.columns]

        day_col = find_col(df.columns, DAY_ALIASES)
        if day_col is None:
            for c in df.columns:
                if parse_day_value(c) is not None and re.search(r"day|^d", str(c), flags=re.I):
                    day_col = c
                    break
        cond_col = find_col(df.columns, COND_ALIASES)
        if cond_col is None:
            cols = list(df.columns)
            if day_col in cols:
                cols.remove(day_col)
            cond_col = cols[0] if cols else None
        if cond_col is None:
            continue

        feed_cols = {}
        for key, aliases in FEED_ALIASES.items():
            col = find_col(df.columns, aliases)
            if col is not None:
                feed_cols[key] = col

        debug.append({"sheet": sheet, "header_row_0based": header, "day_col": day_col, "condition_col": cond_col, "detected_feed_cols": "; ".join([f"{k}={v}" for k, v in feed_cols.items()])})

        if day_col is not None:
            days = df[day_col].apply(parse_day_value).ffill()
        else:
            sheet_day = parse_day_value(sheet)
            days = pd.Series([sheet_day] * len(df), index=df.index)

        for idx, r in df.iterrows():
            cond = str(r.get(cond_col, "")).strip()
            if not cond or cond.lower() == "nan":
                continue
            day = days.loc[idx]
            if day is None or pd.isna(day):
                continue
            row = {"condition": cond, "day": float(day), "culture_volume_mL": default_volume, "sample_removed_mL": 0.0, "source_sheet": sheet}
            for key, col in feed_cols.items():
                val = pd.to_numeric(r.get(col, np.nan), errors="coerce")
                if pd.notna(val):
                    row[key] = float(val)
            row["culture_volume_mL"] = row.get("culture_volume_mL", default_volume)
            row["sample_removed_mL"] = row.get("sample_removed_mL", 0.0)
            rows.append(row)

    events = pd.DataFrame(rows)
    debug_df = pd.DataFrame(debug)
    if events.empty:
        raise ValueError("No feed/process rows parsed from raw Excel.")
    for k in FEED_ALIASES:
        if k not in events.columns:
            events[k] = np.nan
    return events.sort_values(["condition", "day"]).reset_index(drop=True), debug_df

def feed_amount_umol(events_interval, metabolite):
    amount = 0.0
    notes = []
    if metabolite == "Glucose":
        a = (pd.to_numeric(events_interval.get("glucose_feed_mL", 0), errors="coerce").fillna(0) * pd.to_numeric(events_interval.get("glucose_feed_conc_mM", 0), errors="coerce").fillna(0)).sum()
        if a:
            amount += float(a)
            notes.append(f"glucose_bolus_umol={a:g}")

    for feed_name, vol_col in [("feed4", "feed4_mL"), ("cellboost", "cellboost_mL")]:
        comp_col = f"{feed_name}_{metabolite}_mM"
        if vol_col in events_interval.columns and comp_col in events_interval.columns:
            a = (pd.to_numeric(events_interval[vol_col], errors="coerce").fillna(0) * pd.to_numeric(events_interval[comp_col], errors="coerce").fillna(0)).sum()
            if a:
                amount += float(a)
                notes.append(f"{feed_name}_{metabolite}_umol={a:g}")
    return amount, notes

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--raw-excel", required=True)
    ap.add_argument("--process-csv", default="data/metabolomics/processed/converted_process_data.csv")
    ap.add_argument("--metabolomics-csv", default="data/metabolomics/processed/converted_extracellular_metabolomics.csv")
    ap.add_argument("--volume", type=float, default=30.0)
    ap.add_argument("--bound-buffer", type=float, default=0.20)
    args = ap.parse_args()

    base = Path(args.base)
    raw_excel = Path(args.raw_excel)
    if not raw_excel.is_absolute():
        raw_excel = base / raw_excel
    process_csv = Path(args.process_csv)
    if not process_csv.is_absolute():
        process_csv = base / process_csv
    metabolomics_csv = Path(args.metabolomics_csv)
    if not metabolomics_csv.is_absolute():
        metabolomics_csv = base / metabolomics_csv

    outdir = base / "data" / "metabolomics" / "processed"
    outdir.mkdir(parents=True, exist_ok=True)

    process = pd.read_csv(process_csv)
    metab = pd.read_csv(metabolomics_csv)
    events, debug = read_raw_feed_events(raw_excel, default_volume=args.volume)
    events.to_csv(outdir / "feed_events_parsed.csv", index=False, encoding="utf-8-sig")
    debug.to_csv(outdir / "feed_column_detection_debug.csv", index=False, encoding="utf-8-sig")

    process["day"] = pd.to_numeric(process["day"], errors="coerce")
    process["vcd_10e6_cells_mL"] = pd.to_numeric(process["vcd_10e6_cells_mL"], errors="coerce")
    process["culture_volume_mL"] = pd.to_numeric(process.get("culture_volume_mL", args.volume), errors="coerce").fillna(args.volume)

    metab["day"] = pd.to_numeric(metab["day"], errors="coerce")
    metab["concentration"] = pd.to_numeric(metab["concentration"], errors="coerce")

    avg_metab = metab[metab["metabolite"].isin(MET_TO_EX.keys())].dropna(subset=["condition", "day", "metabolite", "concentration"]).groupby(["condition", "day", "metabolite"], as_index=False)["concentration"].mean()
    proc_key = process.set_index(["condition", "day"])
    rows = []

    for cond in sorted(avg_metab["condition"].astype(str).unique()):
        cond_metab = avg_metab[avg_metab["condition"].astype(str) == cond]
        days = sorted(cond_metab["day"].dropna().unique())
        for d0, d1 in zip(days[:-1], days[1:]):
            m0 = cond_metab[cond_metab["day"] == d0].set_index("metabolite")["concentration"]
            m1 = cond_metab[cond_metab["day"] == d1].set_index("metabolite")["concentration"]
            try:
                p0 = proc_key.loc[(cond, d0)]
                p1 = proc_key.loc[(cond, d1)]
            except KeyError:
                continue
            vcd0, vcd1 = float(p0["vcd_10e6_cells_mL"]), float(p1["vcd_10e6_cells_mL"])
            vol0, vol1 = float(p0.get("culture_volume_mL", args.volume)), float(p1.get("culture_volume_mL", args.volume))
            avg_vol = (vol0 + vol1) / 2.0
            delta_day = float(d1 - d0)
            if delta_day <= 0:
                continue
            interval_ivcd_per_ml = ((vcd0 + vcd1) / 2.0) * delta_day
            total_ivcd = interval_ivcd_per_ml * avg_vol
            ev = events[(events["condition"].astype(str) == cond) & (events["day"] > d0) & (events["day"] <= d1)].copy()

            for met, ex in MET_TO_EX.items():
                if met not in m0.index or met not in m1.index:
                    continue
                c0, c1 = float(m0.loc[met]), float(m1.loc[met])
                start_umol = c0 * vol0
                end_umol = c1 * vol1
                feed_umol, feed_notes = feed_amount_umol(ev, met)
                cellular_net_umol = end_umol - start_umol - feed_umol
                if total_ivcd <= 0:
                    continue
                rate = cellular_net_umol / total_ivcd
                if rate < 0:
                    lb, ub = rate * (1 + args.bound_buffer), rate * (1 - args.bound_buffer)
                else:
                    lb, ub = rate * (1 - args.bound_buffer), rate * (1 + args.bound_buffer)
                if lb > ub:
                    lb, ub = ub, lb
                rows.append({
                    "condition": cond,
                    "metabolite": met,
                    "exchange_rxn": ex,
                    "rate": rate,
                    "rate_type": "uptake" if rate < 0 else "secretion",
                    "unit": "mmol_per_1e9_cells_day_feed_corrected",
                    "day_start": d0,
                    "day_end": d1,
                    "lower_bound": lb,
                    "upper_bound": ub,
                    "source": "feed_corrected_mass_balance",
                    "start_conc_mM": c0,
                    "end_conc_mM": c1,
                    "start_volume_mL": vol0,
                    "end_volume_mL": vol1,
                    "feed_amount_umol": feed_umol,
                    "cellular_net_umol": cellular_net_umol,
                    "interval_ivcd_10e6_cell_day_mL": interval_ivcd_per_ml,
                    "total_ivcd_10e6_cell_day": total_ivcd,
                    "note": "; ".join(feed_notes) if feed_notes else "no_feed_composition_used_for_this_metabolite",
                })

    out = pd.DataFrame(rows)
    out_path = outdir / "exchange_rates_feed_corrected.csv"
    out.to_csv(out_path, index=False, encoding="utf-8-sig")
    print("[OK] saved:", out_path)
    print(out.head(30).to_string(index=False))

if __name__ == "__main__":
    main()
