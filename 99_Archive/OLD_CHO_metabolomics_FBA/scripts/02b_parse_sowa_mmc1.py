import argparse
from pathlib import Path
import pandas as pd

DIRECT_MAP = {
    "glucose": "EX_glc_e",
    "lactate": "EX_lac_L_e",
    "glutamine": "EX_gln_L_e",
    "glutamate": "EX_glu_L_e",
    "alanine": "EX_ala_L_e",
    "asparagine": "EX_asn_L_e",
    "aspartate": "EX_asp_L_e",
    "histidine": "EX_his_L_e",
    "isoleucine": "EX_ile_L_e",
    "leucine": "EX_leu_L_e",
    "lysine": "EX_lys_L_e",
    "methionine": "EX_met_L_e",
    "phenylalanine": "EX_phe_L_e",
    "proline": "EX_pro_L_e",
    "serine": "EX_ser_L_e",
    "threonine": "EX_thr_L_e",
    "tryptophane": "EX_trp_L_e",
    "tryptophan": "EX_trp_L_e",
    "tyrosine": "EX_tyr_L_e",
    "valine": "EX_val_L_e",
    "glycine": "EX_gly_e",
}

def find_sheet(xl, kws):
    for s in xl.sheet_names:
        if any(k.lower() in s.lower() for k in kws):
            return s
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--input", default="data/metabolomics/raw/mmc1.xlsx")
    ap.add_argument("--day-start", type=int, default=1)
    ap.add_argument("--day-end", type=int, default=14)
    ap.add_argument("--scale", type=float, default=1.0)
    args = ap.parse_args()

    base = Path(args.base)
    xlsx = Path(args.input)
    if not xlsx.is_absolute():
        xlsx = base / xlsx

    outdir = base / "data" / "metabolomics" / "processed"
    outdir.mkdir(parents=True, exist_ok=True)

    xl = pd.ExcelFile(xlsx)
    nmr_sheet = find_sheet(xl, ["NMR", "Table3"])
    if nmr_sheet is None:
        raise ValueError("NMR sheet not found")

    df = pd.read_excel(xlsx, sheet_name=nmr_sheet, header=0)
    df.columns = [str(c).strip() for c in df.columns]
    sample_col, day_col, pass_col = df.columns[:3]
    df = df[~df[pass_col].isin(["Class", "class", "Separation Method"])].copy()
    df = df[df[pass_col].astype(str).str.match(r"^P\d+", na=False)].copy()
    df[day_col] = pd.to_numeric(df[day_col], errors="coerce")

    met_cols = list(df.columns[3:])
    for c in met_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    ts = df.groupby([pass_col, day_col])[met_cols].mean()
    ts.to_csv(outdir / "sowa_nmr_timeseries.csv", encoding="utf-8-sig")

    rows = []
    for passage in sorted(df[pass_col].dropna().unique()):
        try:
            start = ts.loc[(passage, args.day_start)]
            end = ts.loc[(passage, args.day_end)]
        except KeyError:
            continue

        for met in met_cols:
            key = str(met).strip().lower()
            ex = None
            for k, v in DIRECT_MAP.items():
                if k in key or key in k:
                    ex = v
                    break
            if ex is None:
                continue

            delta = end[met] - start[met]
            rate = (delta / (args.day_end - args.day_start)) * args.scale
            rate_type = "uptake" if rate < 0 else "secretion"
            lb, ub = (rate * 1.2, rate * 0.8) if rate < 0 else (rate * 0.8, rate * 1.2)
            if lb > ub:
                lb, ub = ub, lb

            rows.append({
                "condition": passage,
                "metabolite": met,
                "exchange_rxn": ex,
                "rate": rate,
                "rate_type": rate_type,
                "unit": "NMR_AU_per_day_scaled",
                "day_start": args.day_start,
                "day_end": args.day_end,
                "lower_bound": lb,
                "upper_bound": ub,
                "source": "Sowa_mmc1_NMR_Table3",
                "note": "POC pseudo-rate from NMR intensity; not absolute flux",
            })

    out = pd.DataFrame(rows)
    out_path = outdir / "exchange_rates_sowa_poc.csv"
    out.to_csv(out_path, index=False, encoding="utf-8-sig")
    print("[OK] saved:", out_path)

if __name__ == "__main__":
    main()
