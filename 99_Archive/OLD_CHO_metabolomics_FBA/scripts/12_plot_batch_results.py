from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    args = ap.parse_args()

    base = Path(args.base)
    table_dir = base / "results" / "tables" / "batch"
    fig_dir = base / "results" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    summary_path = table_dir / "batch_objective_summary.csv"
    key_path = table_dir / "batch_key_exchange_summary.csv"
    if not summary_path.exists():
        raise FileNotFoundError(summary_path)

    summary = pd.read_csv(summary_path)
    if summary.empty:
        print("[WARN] empty summary")
        return

    numeric_cols = ["day_start","day_end","fba_objective","pfba_objective","EX_glc_e","EX_lac_L_e","EX_gln_L_e","EX_glu_L_e","EX_nh4_e","biomass_cho_prod","n_constraints","n_warnings"]
    for c in numeric_cols:
        if c in summary.columns:
            summary[c] = pd.to_numeric(summary[c], errors="coerce")

    clone_summary = summary.groupby("condition", as_index=False).agg(
        n_intervals=("interval", "count"),
        n_optimal=("fba_status", lambda x: (x.astype(str) == "optimal").sum()),
        mean_objective=("fba_objective", "mean"),
        max_objective=("fba_objective", "max"),
        mean_pfba_objective=("pfba_objective", "mean"),
        mean_glucose_flux=("EX_glc_e", "mean"),
        mean_lactate_flux=("EX_lac_L_e", "mean"),
        mean_glutamine_flux=("EX_gln_L_e", "mean"),
        mean_glutamate_flux=("EX_glu_L_e", "mean"),
        mean_ammonia_flux=("EX_nh4_e", "mean"),
        total_warnings=("n_warnings", "sum"),
    )
    clone_summary_path = table_dir / "batch_clone_level_summary.csv"
    clone_summary.to_csv(clone_summary_path, index=False, encoding="utf-8-sig")

    for metric, title, fname in [
        ("fba_objective", "FBA objective by clone and interval", "batch_objective_by_interval.png"),
        ("pfba_objective", "pFBA objective by clone and interval", "batch_pfba_objective_by_interval.png"),
    ]:
        fig, ax = plt.subplots(figsize=(12, 6))
        for cond, sub in summary.sort_values(["condition", "day_start"]).groupby("condition"):
            ax.plot(sub["day_end"], sub[metric], marker="o", label=str(cond))
        ax.set_xlabel("Culture day, interval end")
        ax.set_ylabel(metric)
        ax.set_title(title)
        ax.grid(True, linestyle=":", alpha=0.5)
        ax.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        out = fig_dir / fname
        fig.savefig(out, dpi=300)
        print("[OK] saved:", out)

    if key_path.exists():
        key = pd.read_csv(key_path)
        key["flux"] = pd.to_numeric(key["flux"], errors="coerce")
        key["day_end"] = pd.to_numeric(key["day_end"], errors="coerce")
        reactions = ["EX_glc_e", "EX_lac_L_e", "EX_gln_L_e", "EX_glu_L_e", "EX_nh4_e", "biomass_cho_prod"]
        for rxn in reactions:
            sub_rxn = key[key["reaction"] == rxn].copy()
            if sub_rxn.empty:
                continue
            fig, ax = plt.subplots(figsize=(12, 6))
            for cond, sub in sub_rxn.sort_values(["condition", "day_start"]).groupby("condition"):
                ax.plot(sub["day_end"], sub["flux"], marker="o", label=str(cond))
            ax.axhline(0, linewidth=0.8)
            ax.set_xlabel("Culture day, interval end")
            ax.set_ylabel("pFBA/FBA flux")
            ax.set_title(f"{rxn} flux by clone and interval")
            ax.grid(True, linestyle=":", alpha=0.5)
            ax.legend(fontsize=8, ncol=2)
            fig.tight_layout()
            out = fig_dir / f"batch_{rxn}_by_interval.png"
            fig.savefig(out, dpi=300)
            print("[OK] saved:", out)

    print("[OK] saved:", clone_summary_path)
    print(clone_summary.to_string(index=False))

if __name__ == "__main__":
    main()
