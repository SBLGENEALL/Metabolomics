import argparse
from pathlib import Path
import pandas as pd

EXCHANGE_MAP_DEFAULT = {
    "glucose": "EX_glc_e",
    "lactate": "EX_lac_L_e",
    "glutamine": "EX_gln_L_e",
    "glutamate": "EX_glu_L_e",
    "ammonia": "EX_nh4_e",
    "ammonium": "EX_nh4_e",
    "alanine": "EX_ala_L_e",
    "asparagine": "EX_asn_L_e",
    "aspartate": "EX_asp_L_e",
    "serine": "EX_ser_L_e",
    "glycine": "EX_gly_e",
    "proline": "EX_pro_L_e",
    "leucine": "EX_leu_L_e",
    "isoleucine": "EX_ile_L_e",
    "valine": "EX_val_L_e",
    "lysine": "EX_lys_L_e",
    "arginine": "EX_arg_L_e",
    "histidine": "EX_his_L_e",
    "threonine": "EX_thr_L_e",
    "phenylalanine": "EX_phe_L_e",
    "tyrosine": "EX_tyr_L_e",
    "methionine": "EX_met_L_e",
    "tryptophan": "EX_trp_L_e",
    "pyruvate": "EX_pyr_e",
    "citrate": "EX_cit_e",
    "succinate": "EX_succ_e",
    "fumarate": "EX_fum_e",
    "malate": "EX_mal_L_e",
}

def norm(x):
    return str(x).strip().lower()

def find_exchange(met, mapping):
    k = norm(met)
    if k in mapping:
        return mapping[k]
    for name, ex in mapping.items():
        if name in k or k in name:
            return ex
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--interval-start", type=float, default=None)
    ap.add_argument("--interval-end", type=float, default=None)
    ap.add_argument("--bound-buffer", type=float, default=0.20)
    args = ap.parse_args()

    base = Path(args.base)
    xlsx = Path(args.input)
    if not xlsx.is_absolute():
        xlsx = base / xlsx

    outdir = base / "data" / "metabolomics" / "processed"
    outdir.mkdir(parents=True, exist_ok=True)

    metab = pd.read_excel(xlsx, sheet_name="Extracellular_Metabolomics")
    exmap = pd.read_excel(xlsx, sheet_name="Exchange_Map_iCHO3K")

    mapping = dict(EXCHANGE_MAP_DEFAULT)
    for _, r in exmap.dropna(subset=["metabolite", "exchange_rxn"]).iterrows():
        mapping[norm(r["metabolite"])] = str(r["exchange_rxn"]).strip()

    required = {"condition", "day", "metabolite", "concentration"}
    missing = required - set(metab.columns)
    if missing:
        raise ValueError(f"Extracellular_Metabolomics missing columns: {missing}")

    metab = metab.dropna(subset=["condition", "day", "metabolite", "concentration"]).copy()
    metab["day"] = pd.to_numeric(metab["day"], errors="coerce")
    metab["concentration"] = pd.to_numeric(metab["concentration"], errors="coerce")
    metab = metab.dropna(subset=["day", "concentration"])

    avg = metab.groupby(["condition", "day", "metabolite"], as_index=False)["concentration"].mean()
    rows = []

    for cond in sorted(avg["condition"].astype(str).unique()):
        sub = avg[avg["condition"].astype(str) == cond]
        days = sorted(sub["day"].unique())
        if len(days) < 2:
            continue

        intervals = [(args.interval_start, args.interval_end)] if args.interval_start is not None and args.interval_end is not None else list(zip(days[:-1], days[1:]))

        for d0, d1 in intervals:
            s0 = sub[sub["day"] == d0].set_index("metabolite")["concentration"]
            s1 = sub[sub["day"] == d1].set_index("metabolite")["concentration"]
            common = sorted(set(s0.index) & set(s1.index))

            for met in common:
                ex = find_exchange(met, mapping)
                if ex is None:
                    continue

                delta = float(s1.loc[met] - s0.loc[met])
                delta_day = float(d1 - d0)
                if delta_day <= 0:
                    continue

                rate = delta / delta_day
                rate_type = "uptake" if rate < 0 else "secretion"
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
                    "rate_type": rate_type,
                    "unit": "mM_per_day_pseudo",
                    "day_start": d0,
                    "day_end": d1,
                    "lower_bound": lb,
                    "upper_bound": ub,
                    "source": xlsx.name,
                    "note": "POC concentration-derived pseudo-rate. Use feed-corrected/IVCD-normalized rate for final interpretation.",
                })

    out = pd.DataFrame(rows)
    out_path = outdir / "exchange_rates_from_template.csv"
    out.to_csv(out_path, index=False, encoding="utf-8-sig")
    print("[OK] saved:", out_path)
    print(out.head(20).to_string(index=False))

if __name__ == "__main__":
    main()
