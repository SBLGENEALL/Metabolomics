import argparse
from pathlib import Path
import pandas as pd
import escher

def find_flux_file(base, condition):
    for p in [
        base / "results" / "tables" / f"{condition}_pfba_fluxes.csv",
        base / "results" / "tables" / f"{condition}_fba_fluxes.csv",
        base / "results" / "tables" / f"{condition}_relaxed_pfba_fluxes.csv",
        base / "results" / "tables" / f"{condition}_relaxed_fba_fluxes.csv",
    ]:
        if p.exists():
            return p
    raise FileNotFoundError(f"No flux file for {condition}")

def load_flux_data(path):
    df = pd.read_csv(path)
    id_col = "reaction" if "reaction" in df.columns else df.columns[0]
    flux_col = "flux" if "flux" in df.columns else df.columns[1]
    df[flux_col] = pd.to_numeric(df[flux_col], errors="coerce")
    df = df.dropna(subset=[flux_col])
    return dict(zip(df[id_col].astype(str), df[flux_col]))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--condition", required=True)
    ap.add_argument("--model-file", required=True)
    ap.add_argument("--map-json", default=None)
    args = ap.parse_args()

    base = Path(args.base)
    outdir = base / "results" / "escher"
    outdir.mkdir(parents=True, exist_ok=True)

    flux_file = find_flux_file(base, args.condition)
    reaction_data = load_flux_data(flux_file)

    kwargs = {
        "model_json": str(Path(args.model_file)),
        "reaction_data": reaction_data,
        "reaction_styles": ["color", "size", "text"],
    }
    if args.map_json:
        kwargs["map_json"] = str(Path(args.map_json))

    b = escher.Builder(**kwargs)
    out = outdir / f"escher_{args.condition}.html"
    b.save_html(str(out))
    print("[OK] saved:", out)

if __name__ == "__main__":
    main()
