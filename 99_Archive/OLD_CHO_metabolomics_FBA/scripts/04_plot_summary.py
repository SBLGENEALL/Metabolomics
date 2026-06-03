import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

KEY_EX = ["EX_glc_e","EX_lac_L_e","EX_gln_L_e","EX_glu_L_e","EX_nh4_e","EX_ala_L_e","EX_asn_L_e","EX_leu_L_e","EX_ile_L_e","EX_val_L_e"]

def load_flux(path):
    df = pd.read_csv(path)
    id_col = "reaction" if "reaction" in df.columns else df.columns[0]
    fl_col = "flux" if "flux" in df.columns else df.columns[1]
    return dict(zip(df[id_col].astype(str), pd.to_numeric(df[fl_col], errors="coerce").fillna(0)))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--conditions", nargs="+", required=True)
    args = ap.parse_args()

    base = Path(args.base)
    fig_dir = base / "results" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    table_dir = base / "results" / "tables"

    rows = []
    for cond in args.conditions:
        p = table_dir / f"{cond}_pfba_fluxes.csv"
        if not p.exists():
            p = table_dir / f"{cond}_fba_fluxes.csv"
        if not p.exists():
            print("[WARN] missing flux:", cond)
            continue
        flux = load_flux(p)
        row = {"condition": cond}
        for ex in KEY_EX:
            row[ex] = flux.get(ex, 0)
        rows.append(row)

    df = pd.DataFrame(rows)
    if df.empty:
        print("[WARN] no data")
        return

    df.to_csv(table_dir / "summary_key_exchange_fluxes.csv", index=False, encoding="utf-8-sig")

    plot = df.set_index("condition")[KEY_EX].T
    fig, ax = plt.subplots(figsize=(11, 6))
    plot.plot(kind="bar", ax=ax)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Flux")
    ax.set_title("Key exchange flux comparison")
    ax.set_xticklabels(plot.index, rotation=45, ha="right")
    ax.grid(axis="y", linestyle=":", alpha=0.5)
    fig.tight_layout()
    out = fig_dir / "key_exchange_fluxes.png"
    fig.savefig(out, dpi=300)
    print("[OK] saved:", out)

if __name__ == "__main__":
    main()
