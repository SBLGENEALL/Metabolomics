import argparse
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    args = ap.parse_args()

    base = Path(args.base)
    for d in [
        "data/metabolomics/raw",
        "data/metabolomics/processed",
        "model",
        "maps",
        "results/tables",
        "results/tables/batch",
        "results/figures",
        "results/figures_publication",
        "results/escher",
        "logs",
        "templates",
        "scripts",
    ]:
        p = base / d
        p.mkdir(parents=True, exist_ok=True)
        print("[DIR]", p)

    print("[OK] Project folders ready:", base)

if __name__ == "__main__":
    main()
