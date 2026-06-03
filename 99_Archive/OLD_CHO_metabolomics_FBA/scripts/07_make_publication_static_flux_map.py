import argparse
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle
from matplotlib import cm
from matplotlib.colors import Normalize

def find_flux_file(base, cond):
    for p in [
        base / "results" / "tables" / f"{cond}_pfba_fluxes.csv",
        base / "results" / "tables" / f"{cond}_fba_fluxes.csv",
        base / "results" / "tables" / f"{cond}_relaxed_pfba_fluxes.csv",
        base / "results" / "tables" / f"{cond}_relaxed_fba_fluxes.csv",
    ]:
        if p.exists():
            return p
    raise FileNotFoundError(cond)

def load_flux(base, cond):
    df = pd.read_csv(find_flux_file(base, cond))
    id_col = "reaction" if "reaction" in df.columns else df.columns[0]
    fl_col = "flux" if "flux" in df.columns else df.columns[1]
    return dict(zip(df[id_col].astype(str), pd.to_numeric(df[fl_col], errors="coerce").fillna(0)))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--condition", required=True)
    args = ap.parse_args()

    base = Path(args.base)
    cond = args.condition
    flux = load_flux(base, cond)

    outdir = base / "results" / "figures_publication"
    outdir.mkdir(parents=True, exist_ok=True)

    vals = [flux.get(k, 0.0) for k in ["EX_glc_e","EX_lac_L_e","EX_gln_L_e","EX_glu_L_e","EX_nh4_e","biomass_cho_prod"]]
    max_abs = max([abs(v) for v in vals] + [1.0])
    norm = Normalize(vmin=-max_abs, vmax=max_abs)
    cmap = cm.get_cmap("coolwarm")

    fig, ax = plt.subplots(figsize=(12, 7))
    ax.axis("off")
    ax.set_aspect("equal")

    pos = {
        "Glc": (-4, 1.5),
        "Pyr": (-1.5, 1.5),
        "Lac": (1.0, 2.5),
        "Gln": (-4, -1.0),
        "Glu": (-1.5, -1.0),
        "NH4": (0.8, -1.8),
        "AcCoA": (0.8, 0.8),
        "Cit": (2.0, 1.1),
        "aKG": (3.0, 0.0),
        "Succ": (2.0, -1.1),
        "Mal": (0.8, -0.8),
        "Biomass": (3.4, 1.8),
    }

    def metabolite(name, xy, color="#FAD7A0"):
        ax.add_patch(Circle(xy, 0.12, fc=color, ec="#7D6608", lw=1.2, zorder=4))
        ax.text(xy[0], xy[1] - 0.22, name, ha="center", va="top", fontsize=10, fontweight="bold")

    def arrow(a, b, val, label, rad=0):
        color = cmap(norm(val))
        lw = 1.3 + 4 * min(abs(val) / max_abs, 1)
        ax.add_patch(FancyArrowPatch(pos[a], pos[b], arrowstyle="-|>", mutation_scale=16, linewidth=lw, color=color, connectionstyle=f"arc3,rad={rad}", alpha=0.9))
        mx = (pos[a][0] + pos[b][0]) / 2
        my = (pos[a][1] + pos[b][1]) / 2
        ax.text(mx, my, f"{label}\n{val:.3g}", ha="center", va="center", fontsize=8, bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.2))

    for k, p in pos.items():
        metabolite(k, p, "#A9DFBF" if k == "Biomass" else "#FAD7A0")

    ax.text(-4.7, 3.15, f"CHO FBA central-flux schematic ({cond})", fontsize=16, fontweight="bold")
    ax.text(-4.7, 2.85, "Curated schematic: arrow color/width reflects model flux. Use Escher HTML for interactive QC.", fontsize=10)

    arrow("Glc", "Pyr", flux.get("EX_glc_e", 0), "Glc uptake")
    arrow("Pyr", "Lac", flux.get("EX_lac_L_e", 0), "Lac export", rad=0.1)
    arrow("Gln", "Glu", flux.get("EX_gln_L_e", 0), "Gln uptake")
    arrow("Glu", "NH4", flux.get("EX_nh4_e", 0), "NH4 export", rad=-0.1)
    arrow("Pyr", "AcCoA", flux.get("biomass_cho_prod", 0), "to TCA", rad=-0.05)

    for a, b, label in [("AcCoA","Cit","CS/TCA"),("Cit","aKG","IDH"),("aKG","Succ","AKGD"),("Succ","Mal","FUM/MDH"),("Mal","AcCoA","OAA cycle")]:
        arrow(a, b, flux.get("biomass_cho_prod", 0), label, rad=0.15)
    arrow("Cit", "Biomass", flux.get("biomass_cho_prod", 0), "objective", rad=0.1)

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Flux")

    ax.set_xlim(-5, 4.4)
    ax.set_ylim(-2.5, 3.4)

    png = outdir / f"central_flux_schematic_{cond}.png"
    pdf = outdir / f"central_flux_schematic_{cond}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    print("[OK] saved:", png)
    print("[OK] saved:", pdf)

if __name__ == "__main__":
    main()
