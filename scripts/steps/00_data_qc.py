"""
17_data_qc.py — Raw data / rate calculation QC for CHO FBA/FVA pipeline

This step does not run COBRA. It checks whether the culture/metabolite table is
suitable for exchange-rate-constrained FBA/FVA.

Outputs
-------
results/<dataset>/tables/data_qc_required_columns.csv
results/<dataset>/tables/data_qc_warnings.csv
results/<dataset>/tables/data_qc_timeseries_summary.csv
results/<dataset>/tables/data_qc_rate_flags.csv          # if step 1 has run
results/<dataset>/tables/data_qc_carbon_nitrogen_proxy.csv # if step 1 has run
results/<dataset>/figures/Fig13_data_qc_overview.png
"""
import argparse
import os
import pickle
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _find_root():
    cur = os.path.dirname(os.path.abspath(__file__))
    for _ in range(7):
        if os.path.exists(os.path.join(cur, "src", "config.py")):
            return cur
        cur = os.path.dirname(cur)
    return cur

ROOT = _find_root()
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from src.config import *  # noqa

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa", choices=["sowa2020", "practice_20aa", "own_experiment"])
parser.add_argument("--analysis_mode", default=os.environ.get("CHO_ANALYSIS_MODE", "both"))
parser.add_argument("--input_format", default="auto", choices=["auto", "xlsx", "tsv"], help="auto/xlsx/tsv raw data input")
args = parser.parse_args()
DATASET = args.dataset
print("="*65)
print(f"  17_data_qc.py — {DATASET}")
print("="*65)


def _excel_path():
    if DATASET == "practice_20aa":
        return AA20_EXCEL
    if DATASET == "sowa2020":
        return SOWA_MMC1
    if os.path.isdir(OWN_DIR):
        files = [os.path.join(OWN_DIR, f) for f in os.listdir(OWN_DIR) if f.lower().endswith(".xlsx")]
        return sorted(files)[0] if files else None
    return None


def _tsv_dir():
    cands = []
    if DATASET == "practice_20aa":
        cands = [AA20_TSV_DIR, os.path.join(os.path.dirname(AA20_EXCEL), "tsv")]
    elif DATASET == "own_experiment":
        cands = [OWN_DIR, os.path.join(OWN_DIR, "tsv")]
    else:
        cands = [os.path.join(DATA_RAW, DATASET + "_tsv")]
    for d in cands:
        if d and os.path.isdir(d) and os.path.exists(os.path.join(d, RAW_TIMESERIES_TSV)):
            return d
    return cands[0] if cands else None


def _read_first_or_named(path):
    xls = pd.ExcelFile(path)
    sheet = "CHO_raw_data_practice_20AA" if "CHO_raw_data_practice_20AA" in xls.sheet_names else xls.sheet_names[0]
    return pd.read_excel(path, sheet_name=sheet)


def _read_raw_table():
    xlsx_path = _excel_path()
    tsv_dir = _tsv_dir()
    use_tsv = False
    if args.input_format == "tsv":
        use_tsv = True
    elif args.input_format == "xlsx":
        use_tsv = False
    else:
        use_tsv = not (xlsx_path and os.path.exists(xlsx_path))
    if use_tsv:
        raw_tsv = os.path.join(tsv_dir, RAW_TIMESERIES_TSV) if tsv_dir else None
        if not raw_tsv or not os.path.exists(raw_tsv):
            raise FileNotFoundError(f"Raw TSV not found for {DATASET}: {raw_tsv}")
        print(f"  Input: TSV {raw_tsv}")
        return pd.read_csv(raw_tsv, sep="	")
    if not xlsx_path or not os.path.exists(xlsx_path):
        raise FileNotFoundError(f"Raw Excel not found for {DATASET}: {xlsx_path}")
    print(f"  Input: Excel {xlsx_path}")
    return _read_first_or_named(xlsx_path)

raw = _read_raw_table()
out_tables = results_dir(DATASET, "tables")
out_fig = results_dir(DATASET, "figures")

required = [
    "Sample ID", "DAY", "IgG", "Viable Density", "Total Density", "Viability",
    "Culture Volume mL", "Gluc", "Lac", "NH4+",
]
aa_cols = ["Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val"]
optional = ["Sample Removed mL", "FunctionMax mL", "CellBoost 7A mL", "CellBoost 7B mL", "Osm", "pH", "PO2"] + aa_cols

req_rows = []
for c in required + optional:
    req_rows.append({
        "column": c,
        "required": c in required,
        "present": c in raw.columns,
        "missing_fraction": float(raw[c].isna().mean()) if c in raw.columns else 1.0,
    })
req_df = pd.DataFrame(req_rows)
req_df.to_csv(os.path.join(out_tables, "data_qc_required_columns.csv"), index=False)

warnings_rows = []
def warn(level, check, message, clone="", day=""):
    warnings_rows.append({"level": level, "check": check, "message": message, "clone": clone, "day": day})

for c in required:
    if c not in raw.columns:
        warn("ERROR", "required_column", f"Missing required column: {c}")

if "Sample ID" in raw.columns and "DAY" in raw.columns:
    dup = raw.duplicated(["Sample ID", "DAY"], keep=False)
    if dup.any():
        for _, r in raw.loc[dup, ["Sample ID", "DAY"]].drop_duplicates().iterrows():
            warn("ERROR", "duplicate_clone_day", "Duplicate Sample ID + DAY rows", str(r["Sample ID"]), str(r["DAY"]))

for c in [x for x in ["IgG", "Viable Density", "Total Density", "Viability", "Culture Volume mL", "Gluc", "Lac", "NH4+"] + aa_cols if x in raw.columns]:
    vals = pd.to_numeric(raw[c], errors="coerce")
    if vals.isna().mean() > 0.1:
        warn("WARN", "numeric_parse", f"Column {c} has >10% non-numeric/missing values")
    if c not in ["pH"] and (vals < 0).any():
        warn("WARN", "negative_value", f"Column {c} contains negative values")

if "Viability" in raw.columns:
    v = pd.to_numeric(raw["Viability"], errors="coerce")
    if (v > 100).any() or (v < 0).any():
        warn("WARN", "viability_range", "Viability outside 0-100 detected")

feed_cols = [c for c in ["FunctionMax mL", "CellBoost 7A mL", "CellBoost 7B mL"] if c in raw.columns]
if feed_cols and "Sample ID" in raw.columns and "DAY" in raw.columns:
    for clone, g in raw.sort_values("DAY").groupby("Sample ID"):
        for fc in feed_cols:
            vals = pd.to_numeric(g[fc], errors="coerce").fillna(0).values
            if len(vals) > 1 and np.any(np.diff(vals) < -1e-9):
                warn("INFO", "feed_volume_pattern", f"{fc} decreases over time; likely interval feed volumes, not cumulative", clone)

summary_rows = []
if "Sample ID" in raw.columns:
    for clone, g in raw.groupby("Sample ID"):
        d = sorted(pd.to_numeric(g.get("DAY"), errors="coerce").dropna().unique()) if "DAY" in g else []
        row = {
            "clone": clone,
            "n_rows": len(g),
            "days": ",".join(str(int(x)) if float(x).is_integer() else str(x) for x in d),
            "day_min": min(d) if d else np.nan,
            "day_max": max(d) if d else np.nan,
            "IgG_last": float(g.sort_values("DAY")["IgG"].iloc[-1]) if "DAY" in g and "IgG" in g else np.nan,
            "VCD_last": float(g.sort_values("DAY")["Viable Density"].iloc[-1]) if "DAY" in g and "Viable Density" in g else np.nan,
            "Viability_last": float(g.sort_values("DAY")["Viability"].iloc[-1]) if "DAY" in g and "Viability" in g else np.nan,
        }
        summary_rows.append(row)
summary = pd.DataFrame(summary_rows)
summary.to_csv(os.path.join(out_tables, "data_qc_timeseries_summary.csv"), index=False)

# Optional rate-level QC if step 1 was already run.
pkl = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
rate_flags = []
carbon_n_proxy = []
if os.path.exists(pkl):
    with open(pkl, "rb") as f:
        data = pickle.load(f)
    rate_df = data.get("rate_df", pd.DataFrame())
    if not rate_df.empty:
        piv = rate_df.pivot_table(index="clone", columns="column", values="rate_mmol_gDCWh", aggfunc="first")
        for clone, row in piv.iterrows():
            glc = float(row.get("Gluc", 0) or 0)
            lac = float(row.get("Lac", 0) or 0)
            nh4 = float(row.get("NH4+", 0) or 0)
            gln = float(row.get("Gln", 0) or 0)
            lac_glc = lac / abs(glc) if abs(glc) > 1e-12 else np.nan
            aa_uptake = sum(abs(float(row.get(c, 0) or 0)) for c in aa_cols if float(row.get(c, 0) or 0) < 0)
            carbon_n_proxy.append({"clone": clone, "glc_q": glc, "lac_q": lac, "nh4_q": nh4, "gln_q": gln, "lac_glc_ratio": lac_glc, "total_aa_uptake_abs": aa_uptake})
            if glc > 1e-9:
                rate_flags.append({"level":"WARN", "clone":clone, "column":"Gluc", "message":"Glucose rate is secretion-positive; check sign/units/feed correction"})
            if abs(lac_glc) > 2 and np.isfinite(lac_glc):
                rate_flags.append({"level":"WARN", "clone":clone, "column":"Lac/Gluc", "message":f"Large lactate/glucose ratio: {lac_glc:.2f}"})
            if nh4 < -1e-9:
                rate_flags.append({"level":"INFO", "clone":clone, "column":"NH4+", "message":"NH4 rate is uptake-negative; check sign if unexpected"})
        pd.DataFrame(rate_flags).to_csv(os.path.join(out_tables, "data_qc_rate_flags.csv"), index=False)
        pd.DataFrame(carbon_n_proxy).to_csv(os.path.join(out_tables, "data_qc_carbon_nitrogen_proxy.csv"), index=False)

warnings_df = pd.DataFrame(warnings_rows)
warnings_df.to_csv(os.path.join(out_tables, "data_qc_warnings.csv"), index=False)

# Simple QC figure
try:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    ax = axes[0,0]
    if not summary.empty:
        ax.bar(summary["clone"], summary["IgG_last"])
        ax.set_title("Final IgG by clone")
        ax.tick_params(axis="x", rotation=45)
    ax = axes[0,1]
    if "DAY" in raw.columns and "IgG" in raw.columns:
        for clone, g in raw.groupby("Sample ID"):
            g = g.sort_values("DAY")
            ax.plot(g["DAY"], g["IgG"], marker="o", lw=1, label=clone)
        ax.set_title("IgG time course")
        ax.set_xlabel("Day")
        ax.set_ylabel("IgG")
    ax = axes[1,0]
    if carbon_n_proxy:
        cn = pd.DataFrame(carbon_n_proxy)
        ax.scatter(cn["glc_q"], cn["lac_q"])
        for _, r in cn.iterrows():
            ax.text(r["glc_q"], r["lac_q"], r["clone"], fontsize=7)
        ax.axhline(0, lw=0.7); ax.axvline(0, lw=0.7)
        ax.set_title("Rate QC: glucose vs lactate")
        ax.set_xlabel("Glucose q")
        ax.set_ylabel("Lactate q")
    ax = axes[1,1]
    counts = req_df.groupby(["required", "present"]).size().reset_index(name="n")
    labels = [f"{'Req' if r.required else 'Opt'} / {'Present' if r.present else 'Missing'}" for _, r in counts.iterrows()]
    ax.bar(labels, counts["n"])
    ax.set_title("Column availability")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(os.path.join(out_fig, "Fig13_data_qc_overview.png"), dpi=200)
    plt.close(fig)
except Exception as e:
    warn("WARN", "qc_figure", f"Could not generate QC figure: {e}")

n_err = sum(1 for w in warnings_rows if w["level"] == "ERROR")
n_warn = sum(1 for w in warnings_rows if w["level"] == "WARN")
print(f"  [saved] results/{DATASET}/tables/data_qc_required_columns.csv")
print(f"  [saved] results/{DATASET}/tables/data_qc_warnings.csv")
print(f"  QC summary: {n_err} ERROR, {n_warn} WARN")
print("  OK data QC complete")
