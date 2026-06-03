import argparse
from pathlib import Path
import cobra
import pandas as pd

def load_model(model_file):
    p = Path(model_file)
    if p.suffix.lower() == ".json":
        return cobra.io.load_json_model(str(p))
    return cobra.io.read_sbml_model(str(p))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--model-file", required=True)
    args = ap.parse_args()

    base = Path(args.base)
    outdir = base / "results" / "tables"
    outdir.mkdir(parents=True, exist_ok=True)

    print("[INFO] Loading model:", args.model_file)
    m = load_model(args.model_file)
    print("[MODEL]")
    print("id          :", m.id)
    print("reactions   :", len(m.reactions))
    print("metabolites :", len(m.metabolites))
    print("genes       :", len(m.genes))
    print("objective   :", m.objective.expression)

    ex = []
    for r in m.exchanges:
        ex.append({
            "id": r.id,
            "name": r.name,
            "lower_bound": r.lower_bound,
            "upper_bound": r.upper_bound,
            "metabolites": "; ".join([met.id for met in r.metabolites]),
            "metabolite_names": "; ".join([met.name or "" for met in r.metabolites]),
        })
    pd.DataFrame(ex).to_csv(outdir / "exchange_reactions.csv", index=False, encoding="utf-8-sig")

    sol = m.optimize()
    print(f"[BASELINE FBA] status={sol.status}, objective={sol.objective_value}")
    print("[OK] exchange reactions saved:", outdir / "exchange_reactions.csv")

if __name__ == "__main__":
    main()
