"""
01_load_data.py
Feed-Corrected Specific Rate 계산 (완전 수정판)
Ref: Goudar et al. (2005) Biotechnol Prog 21:1193
사용: python scripts/steps/01_load_data.py --dataset practice_20aa --rate_days 7,10
"""
import sys, os, argparse, warnings, pickle
warnings.filterwarnings("ignore")

def _find_root():
    cur = os.path.dirname(os.path.abspath(__file__))
    for _ in range(6):
        if os.path.exists(os.path.join(cur, "src", "config.py")): return cur
        cur = os.path.dirname(cur)
    return cur
ROOT = _find_root()
if ROOT not in sys.path: sys.path.insert(0, ROOT)
from src.config import *
import pandas as pd, numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--dataset",   default="practice_20aa",
                    choices=["sowa2020","practice_20aa","own_experiment"])
parser.add_argument("--rate_days", default="7,10")
parser.add_argument("--top_n",     type=int, default=3)
parser.add_argument("--feed_volume_mode", default="interval", choices=["interval", "cumulative"],
                    help="How to interpret feed volume columns. interval=volume added during interval ending at d2; cumulative=use d2-d1.")
parser.add_argument("--input_format", default="auto", choices=["auto", "xlsx", "tsv"],
                    help="auto=use xlsx if available, otherwise TSV; tsv=read raw_timeseries.tsv/metabolite_map.tsv/feed_composition.tsv")
args = parser.parse_args()
DATASET  = args.dataset
DAY_S, DAY_E = [int(x) for x in args.rate_days.split(",")]
TOP_N = args.top_n
FEED_VOLUME_MODE = args.feed_volume_mode
INPUT_FORMAT = args.input_format

print("="*65)
print(f"  01_load_data.py — {DATASET}  Day{DAY_S}→{DAY_E}")
print("="*65)
print("  Ref: Goudar et al. (2005) Biotechnol Prog 21:1193")
print(f"  Feed volume mode: {FEED_VOLUME_MODE}")
print(f"  Input format    : {INPUT_FORMAT}")

# -------------------------------------------------------------------------
# Input source resolution: Excel or TSV
# -------------------------------------------------------------------------
def _tsv_dir_for_dataset():
    if DATASET == "practice_20aa":
        # Prefer the dedicated TSV folder; fallback to data/raw/practice_20aa/tsv.
        cands = [AA20_TSV_DIR, os.path.join(os.path.dirname(AA20_EXCEL), "tsv")]
    elif DATASET == "own_experiment":
        cands = [OWN_DIR, os.path.join(OWN_DIR, "tsv")]
    else:
        cands = [os.path.join(DATA_RAW, DATASET + "_tsv")]
    for d in cands:
        if d and os.path.isdir(d) and os.path.exists(os.path.join(d, RAW_TIMESERIES_TSV)):
            return d
    return cands[0] if cands else None

EXCEL = {"practice_20aa": AA20_EXCEL, "sowa2020": SOWA_MMC1}.get(DATASET)
if DATASET == "own_experiment":
    files = [f for f in os.listdir(OWN_DIR) if f.lower().endswith(".xlsx")] if os.path.exists(OWN_DIR) else []
    EXCEL = os.path.join(OWN_DIR, sorted(files)[0]) if files else None
TSV_DIR = _tsv_dir_for_dataset()

USE_TSV = False
if INPUT_FORMAT == "tsv":
    USE_TSV = True
elif INPUT_FORMAT == "xlsx":
    USE_TSV = False
else:
    USE_TSV = not (EXCEL and os.path.exists(EXCEL))

if USE_TSV:
    raw_tsv = os.path.join(TSV_DIR, RAW_TIMESERIES_TSV) if TSV_DIR else None
    if not raw_tsv or not os.path.exists(raw_tsv):
        print(f"  !! TSV 입력 파일 없음: {raw_tsv}"); sys.exit(1)
    print(f"  입력: TSV folder = {TSV_DIR}")
else:
    if not EXCEL or not os.path.exists(EXCEL):
        print(f"  !! Excel 파일 없음: {EXCEL}"); sys.exit(1)
    print(f"  입력: Excel = {os.path.basename(EXCEL)}")

# -------------------------------------------------------------------------
# Excel reader helpers
# -------------------------------------------------------------------------
# Some generated practice workbooks may contain the raw time-course sheet as
# Excel's default "Sheet1" instead of "CHO_raw_data_practice_20AA". The raw
# demo workbook can also be regenerated without carrying Metabolite_Map and
# Feed_Composition sheets. To make the pipeline robust, read the raw sheet
# flexibly and fall back to bundled template workbooks for map/feed sheets.
def _read_tsv_required(tsv_dir, filename):
    path = os.path.join(tsv_dir, filename)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Required TSV not found: {path}")
    return pd.read_csv(path, sep="	")


def _sheet_names(path):
    try:
        return pd.ExcelFile(path).sheet_names
    except Exception:
        return []


def _read_excel_flexible(path, preferred_sheet=None, fallback_first=False):
    sheets = _sheet_names(path)
    if preferred_sheet and preferred_sheet in sheets:
        return pd.read_excel(path, sheet_name=preferred_sheet)
    if fallback_first and sheets:
        print(f"  ⚠ Sheet '{preferred_sheet}' not found in {os.path.basename(path)}; using first sheet '{sheets[0]}'")
        return pd.read_excel(path, sheet_name=sheets[0])
    raise ValueError(f"Worksheet named '{preferred_sheet}' not found in {path}; available={sheets}")


def _candidate_template_workbooks():
    here_raw = os.path.dirname(EXCEL) if EXCEL else ""
    candidates = [
        EXCEL,
        os.path.join(here_raw, "CHO_raw_data_practice_20AA_original_3clone_demo.xlsx"),
        os.path.join(os.path.dirname(here_raw), "practice_20aa_original", "CHO_raw_data_practice_20AA_original.xlsx"),
        os.path.join(os.path.dirname(here_raw), "practice_20aa_original", "CHO_raw_data_practice_20AA_previous_9clone.xlsx"),
    ]
    out = []
    for c in candidates:
        if c and os.path.exists(c) and c not in out:
            out.append(c)
    return out


def _read_required_sheet_from_candidates(sheet_name):
    for path in _candidate_template_workbooks():
        if sheet_name in _sheet_names(path):
            if path != EXCEL:
                print(f"  ⚠ Sheet '{sheet_name}' loaded from template: {os.path.basename(path)}")
            return pd.read_excel(path, sheet_name=sheet_name)
    raise ValueError(f"Worksheet named '{sheet_name}' not found in workbook or bundled templates")

if DATASET in ["practice_20aa", "own_experiment"]:
    if USE_TSV:
        df_raw = _read_tsv_required(TSV_DIR, RAW_TIMESERIES_TSV)
        df_map = _read_tsv_required(TSV_DIR, METABOLITE_MAP_TSV)
        try:
            fc_df = _read_tsv_required(TSV_DIR, FEED_COMPOSITION_TSV).set_index("Feed")
            has_feed = True
            print(f"  Feed TSV: {fc_df.shape[0]}종 × {fc_df.shape[1]}대사물질")
        except Exception as e:
            fc_df = pd.DataFrame(); has_feed = False
            print(f"  ⚠ Feed_Composition TSV 없음 또는 읽기 실패: {e}")
    else:
        df_raw = _read_excel_flexible(EXCEL, "CHO_raw_data_practice_20AA", fallback_first=True)
        df_map = _read_required_sheet_from_candidates("Metabolite_Map")
        try:
            fc_df = _read_required_sheet_from_candidates("Feed_Composition").set_index("Feed")
            has_feed = True
            print(f"  Feed: {fc_df.shape[0]}종 × {fc_df.shape[1]}대사물질")
        except Exception as e:
            fc_df = pd.DataFrame(); has_feed = False
            print(f"  ⚠ Feed_Composition 없음 또는 읽기 실패: {e}")

    MET_COLS = df_map["Column"].tolist()
    COL2EX   = dict(zip(df_map["Column"], df_map["Exchange reaction expected by v4"]))
    COL2NAME = dict(zip(df_map["Column"], df_map["Canonical metabolite"]))
    FEED_COLS = {
        "FunctionMax":  "FunctionMax mL",
        "CellBoost 7A": "CellBoost 7A mL",
        "CellBoost 7B": "CellBoost 7B mL",
    }

    clones   = sorted(df_raw["Sample ID"].unique())
    days_all = sorted(df_raw["DAY"].unique())
    igG_last = df_raw[df_raw["DAY"]==max(days_all)].set_index("Sample ID")["IgG"].to_dict()
    sorted_c = sorted(igG_last, key=igG_last.get, reverse=True)

    # Optional explicit biological grouping for demo/real datasets.
    # If a Group column exists, use it instead of forcing top/bottom ranking.
    # Accepted labels include High/HighProducer, Low/LowProducer, Mother/Parental.
    clone_groups = {}
    if "Group" in df_raw.columns:
        clone_groups = df_raw.groupby("Sample ID")["Group"].first().fillna("Unassigned").to_dict()
        def _norm_group(x):
            return str(x).strip().lower().replace(" ", "").replace("_", "")
        HIGH_CLONES = [c for c in clones if _norm_group(clone_groups.get(c, "")) in {"high", "highproducer", "highproducing", "highproducerclone"}]
        LOW_CLONES  = [c for c in clones if _norm_group(clone_groups.get(c, "")) in {"low", "lowproducer", "lowproducing", "lowproducerclone"}]
        if not HIGH_CLONES or not LOW_CLONES:
            n_eff = min(TOP_N, max(1, len(sorted_c)//3))
            HIGH_CLONES = sorted_c[:n_eff]
            LOW_CLONES  = sorted_c[-n_eff:]
    else:
        n_eff = min(TOP_N, max(1, len(sorted_c)//3))
        HIGH_CLONES = sorted_c[:n_eff]
        LOW_CLONES  = sorted_c[-n_eff:]
        clone_groups = {c: ("High" if c in HIGH_CLONES else ("Low" if c in LOW_CLONES else "Mid")) for c in clones}

    print(f"\n  Clones: {clones}")
    print(f"  Days  : {days_all}")
    print(f"  Groups: {clone_groups}")
    print(f"  High  : {HIGH_CLONES}")
    print(f"  Low   : {LOW_CLONES}")

    def compute_rate(clone, d1, d2):
        sub  = df_raw[df_raw["Sample ID"]==clone].sort_values("DAY")
        rows = {int(r["DAY"]): r for _, r in sub.iterrows()}
        if d1 not in rows or d2 not in rows: return {}
        r1, r2 = rows[d1], rows[d2]

        vcd1  = float(r1["Viable Density"])
        vcd2  = float(r2["Viable Density"])
        v1    = float(r1.get("Culture Volume mL", 50))
        v2    = float(r2.get("Culture Volume mL", 50))
        v_avg = (v1 + v2) / 2
        dt_h  = (d2 - d1) * 24

        # IVCD [gDCW·h]
        ivcd = ((vcd1+vcd2)/2) * 1e6 * v_avg * dt_h * 8e-12
        if ivcd <= 0: return {}

        # Feed volume interpretation
        # Most practice_20aa files store feed columns as interval additions for the
        # sampling interval ending at the current day, not as cumulative totals.
        # For Day 7→10, use the feed volumes recorded on Day 10 by default.
        # If your feed columns are truly cumulative, run with --feed_volume_mode cumulative.
        feed_vols = {}
        for fname, fcol in FEED_COLS.items():
            if fcol not in df_raw.columns:
                continue
            fv1 = float(r1.get(fcol, 0)) if pd.notna(r1.get(fcol)) else 0.0
            fv2 = float(r2.get(fcol, 0)) if pd.notna(r2.get(fcol)) else 0.0
            if FEED_VOLUME_MODE == "cumulative":
                feed_vols[fname] = max(0.0, fv2 - fv1)
            else:
                feed_vols[fname] = max(0.0, fv2)

        rates = {}
        for col in MET_COLS:
            if col not in df_raw.columns: continue
            c1 = float(r1[col]) if pd.notna(r1.get(col)) else 0.0
            c2 = float(r2[col]) if pd.notna(r2.get(col)) else 0.0

            # Observed amount change in the vessel. Using C*V is more robust than
            # (C2-C1)*Vavg when volume changes due to sampling/feed.
            delta_obs = (c2 * v2) - (c1 * v1)  # µmol

            feed_add = 0.0
            if has_feed:
                for fname, vol in feed_vols.items():
                    if fname in fc_df.index and col in fc_df.columns:
                        feed_add += vol * float(fc_df.loc[fname, col])

            net_mmol = (delta_obs - feed_add) / 1000
            q        = net_mmol / ivcd

            rates[col] = {
                "exchange_id":      COL2EX.get(col, ""),
                "metabolite_name":  COL2NAME.get(col, col),
                "c_start_mM":       c1,
                "c_end_mM":         c2,
                "delta_obs_umol":   delta_obs,
                "feed_add_umol":    feed_add,
                "net_mmol":         net_mmol,
                "ivcd_gDCWh":       ivcd,
                "rate_mmol_gDCWh":  q,
                "feed_volume_mode": FEED_VOLUME_MODE,
            }
        return rates

    print(f"\n  {'Clone':12s} {'Glc_q':>9s} {'Lac_q':>9s} {'Lac/Glc':>9s}  IgG")
    print("  " + "─"*52)

    all_rates    = {}
    rate_records = []
    for clone in clones:
        r = compute_rate(clone, DAY_S, DAY_E)
        if not r: continue
        all_rates[clone] = r
        gq = r.get("Gluc", {}).get("rate_mmol_gDCWh", 0)
        lq = r.get("Lac",  {}).get("rate_mmol_gDCWh", 0)
        ratio = lq/abs(gq) if abs(gq) > 1e-6 else 0
        tag = "★" if clone in HIGH_CLONES else ("▼" if clone in LOW_CLONES else " ")
        print(f"  {tag} {clone:10s} {gq:9.4f} {lq:9.4f} {ratio:9.3f}  {igG_last.get(clone,0):.0f}")
        for col, info in r.items():
            rate_records.append({"clone": clone, "column": col, **info})

    rate_df = pd.DataFrame(rate_records)

    # FBA Constraints (단방향 보정)
    all_constraints = {}
    for clone, rates in all_rates.items():
        cst = {}
        for col, info in rates.items():
            ex_id = info["exchange_id"]
            if not ex_id: continue
            q = info["rate_mmol_gDCWh"]
            if abs(q) < 1e-9: continue
            if q < 0:   cst[ex_id] = (q * 1.20, 0.0)    # uptake
            else:       cst[ex_id] = (0.0, q * 1.20)    # secretion
        all_constraints[clone] = cst

    n_total = sum(len(v) for v in all_constraints.values())
    print(f"\n  Constraints: {n_total} total ({n_total//len(all_constraints) if all_constraints else 0}/clone avg)")

    os.makedirs(DATA_PROCESSED, exist_ok=True)
    pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
    save_data = {
        "df_raw":           df_raw,
        "all_rates":        all_rates,
        "rate_df":          rate_df,
        "all_constraints":  all_constraints,
        "igG_day14":        igG_last,
        "high_clones":      HIGH_CLONES,
        "low_clones":       LOW_CLONES,
        "sorted_clones":    sorted_c,
        "clone_groups":     clone_groups,
        "clones":           clones,
        "MET_COLS":         MET_COLS,
        "COL2EX":           COL2EX,
        "COL2NAME":         COL2NAME,
        "day_interval":     (DAY_S, DAY_E),
        "dataset":          DATASET,
    }
    with open(pkl_path, "wb") as f:
        pickle.dump(save_data, f)

    out = results_dir(DATASET, "tables")
    rate_df.to_csv(os.path.join(out, "exchange_rates.csv"), index=False)
    df_raw.to_csv(os.path.join(out, "raw_timeseries.csv"),  index=False)
    print(f"  [saved] data/processed/rates_{DATASET}.pkl")
    print(f"  [saved] results/{DATASET}/tables/exchange_rates.csv")

print(f"\n  OK  완료")
print(f"  다음: python scripts/steps/02_map_metabolites.py --dataset {DATASET}")
