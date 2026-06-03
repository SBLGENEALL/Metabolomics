
"""
14_make_clone_decision_summary.py

Make a clone-level decision summary after running:
  09_convert_existing_excel_to_template.py
  13_recalculate_feed_corrected_rates.py
  11_batch_run_intervals.py
  12_plot_batch_results.py

It combines:
  - final titer / qP / viability from converted_process_data.csv
  - feed-corrected exchange rates from exchange_rates_feed_corrected.csv
  - batch FBA summary from results/tables/batch/batch_clone_level_summary.csv

Outputs:
  results/tables/clone_decision_summary.csv
  results/figures/clone_decision_score.png
  results/figures/final_titer_vs_lactate.png

Run:
  cd C:\CHO_POC_ChatGPT_260519
  python scripts\14_make_clone_decision_summary.py --base C:\CHO_POC_ChatGPT_260519
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def minmax_good_high(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    if s.notna().sum() == 0:
        return pd.Series([0.5] * len(s), index=s.index)
    lo, hi = s.min(), s.max()
    if hi == lo:
        return pd.Series([0.5] * len(s), index=s.index)
    return (s - lo) / (hi - lo)


def minmax_good_low(s: pd.Series) -> pd.Series:
    return 1 - minmax_good_high(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--process-csv", default="data/metabolomics/processed/converted_process_data.csv")
    ap.add_argument("--rates-csv", default="data/metabolomics/processed/exchange_rates_feed_corrected.csv")
    ap.add_argument("--batch-summary", default="results/tables/batch/batch_clone_level_summary.csv")
    args = ap.parse_args()

    base = Path(args.base)
    table_dir = base / "results" / "tables"
    fig_dir = base / "results" / "figures"
    table_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    process_path = Path(args.process_csv)
    if not process_path.is_absolute():
        process_path = base / process_path

    rates_path = Path(args.rates_csv)
    if not rates_path.is_absolute():
        rates_path = base / rates_path

    batch_path = Path(args.batch_summary)
    if not batch_path.is_absolute():
        batch_path = base / batch_path

    process = pd.read_csv(process_path)
    process["day"] = pd.to_numeric(process["day"], errors="coerce")
    for c in ["titer_g_L", "qP_pg_cell_day", "viability_pct", "vcd_10e6_cells_mL", "ivcd_10e6_cell_day_mL"]:
        if c in process.columns:
            process[c] = pd.to_numeric(process[c], errors="coerce")

    # Final day row per condition
    final_idx = process.sort_values(["condition", "day"]).groupby("condition")["day"].idxmax()
    final = process.loc[final_idx].copy()
    final = final.rename(columns={
        "day": "final_day",
        "titer_g_L": "final_titer_g_L",
        "qP_pg_cell_day": "last_interval_qP_pg_cell_day",
        "viability_pct": "final_viability_pct",
        "vcd_10e6_cells_mL": "final_vcd_10e6_cells_mL",
        "ivcd_10e6_cell_day_mL": "final_cumulative_ivcd",
    })
    final = final[[
        "condition", "final_day", "final_titer_g_L", "last_interval_qP_pg_cell_day",
        "final_viability_pct", "final_vcd_10e6_cells_mL", "final_cumulative_ivcd"
    ]]

    qp = process.groupby("condition", as_index=False).agg(
        mean_qP_pg_cell_day=("qP_pg_cell_day", "mean"),
        max_qP_pg_cell_day=("qP_pg_cell_day", "max"),
        mean_viability_pct=("viability_pct", "mean"),
        peak_vcd_10e6_cells_mL=("vcd_10e6_cells_mL", "max"),
    )

    summary = final.merge(qp, on="condition", how="left")

    # Rates summary
    if rates_path.exists():
        rates = pd.read_csv(rates_path)
        rates["rate"] = pd.to_numeric(rates["rate"], errors="coerce")
        rates["day_end"] = pd.to_numeric(rates["day_end"], errors="coerce")

        pivot_mean = rates.pivot_table(
            index="condition",
            columns="metabolite",
            values="rate",
            aggfunc="mean",
        ).reset_index()
        pivot_mean.columns = ["condition"] + [f"mean_rate_{c}" for c in pivot_mean.columns[1:]]
        summary = summary.merge(pivot_mean, on="condition", how="left")

        # late interval summary: use last interval per condition/metabolite
        late = (
            rates.sort_values(["condition", "metabolite", "day_end"])
            .groupby(["condition", "metabolite"], as_index=False)
            .tail(1)
        )
        late_pivot = late.pivot_table(
            index="condition",
            columns="metabolite",
            values="rate",
            aggfunc="mean",
        ).reset_index()
        late_pivot.columns = ["condition"] + [f"late_rate_{c}" for c in late_pivot.columns[1:]]
        summary = summary.merge(late_pivot, on="condition", how="left")
    else:
        print("[WARN] rates csv not found:", rates_path)

    # Batch FBA summary
    if batch_path.exists():
        batch = pd.read_csv(batch_path)
        summary = summary.merge(batch, on="condition", how="left", suffixes=("", "_batch"))
    else:
        print("[WARN] batch summary not found:", batch_path)

    # Decision score:
    # Higher final titer, qP, viability are good.
    # Lower lactate secretion, ammonia secretion, glucose uptake burden are generally good.
    # Here glucose uptake burden uses abs(mean_rate_Glucose), because rate is usually negative.
    for col in summary.columns:
        if col != "condition":
            summary[col] = pd.to_numeric(summary[col], errors="ignore")

    score = pd.Series(0.0, index=summary.index)
    weights = {}

    def add_score(name, series, weight, direction="high"):
        nonlocal score
        weights[name] = weight
        if direction == "high":
            score += weight * minmax_good_high(series)
        else:
            score += weight * minmax_good_low(series)

    add_score("final_titer", summary.get("final_titer_g_L"), 0.35, "high")
    add_score("mean_qP", summary.get("mean_qP_pg_cell_day"), 0.20, "high")
    add_score("final_viability", summary.get("final_viability_pct"), 0.15, "high")

    # lactate/ammonia lower is preferred
    if "mean_rate_Lactate" in summary.columns:
        add_score("mean_lactate_rate_low", summary["mean_rate_Lactate"], 0.10, "low")
    if "mean_rate_Ammonia" in summary.columns:
        add_score("mean_ammonia_rate_low", summary["mean_rate_Ammonia"], 0.10, "low")

    # Lower absolute glucose uptake burden can be interpreted as higher metabolic efficiency in this synthetic demo.
    if "mean_rate_Glucose" in summary.columns:
        add_score("glucose_burden_low", summary["mean_rate_Glucose"].abs(), 0.10, "low")

    summary["decision_score_0to1"] = score
    if sum(weights.values()) > 0:
        summary["decision_score_0to1"] = summary["decision_score_0to1"] / sum(weights.values())

    summary["rank"] = summary["decision_score_0to1"].rank(ascending=False, method="min").astype(int)

    # Simple recommendation labels
    def label(row):
        if row["rank"] <= 3:
            return "Top candidate"
        if row["rank"] <= 6:
            return "Middle candidate"
        return "Lower candidate"

    summary["recommendation"] = summary.apply(label, axis=1)

    summary = summary.sort_values(["rank", "condition"])

    out = table_dir / "clone_decision_summary.csv"
    summary.to_csv(out, index=False, encoding="utf-8-sig")
    print("[OK] saved:", out)

    # Score plot
    fig, ax = plt.subplots(figsize=(10, 5.5))
    plot_df = summary.sort_values("decision_score_0to1", ascending=True)
    ax.barh(plot_df["condition"], plot_df["decision_score_0to1"])
    ax.set_xlabel("Decision score, 0-1")
    ax.set_title("Clone decision score: titer/qP/viability/metabolic burden")
    ax.grid(axis="x", linestyle=":", alpha=0.5)
    fig.tight_layout()
    p = fig_dir / "clone_decision_score.png"
    fig.savefig(p, dpi=300)
    print("[OK] saved:", p)

    # Final titer vs lactate
    if "mean_rate_Lactate" in summary.columns:
        fig, ax = plt.subplots(figsize=(7, 5.5))
        ax.scatter(summary["mean_rate_Lactate"], summary["final_titer_g_L"])
        for _, r in summary.iterrows():
            ax.text(r["mean_rate_Lactate"], r["final_titer_g_L"], str(r["condition"]), fontsize=8, ha="left", va="bottom")
        ax.set_xlabel("Mean feed-corrected lactate exchange rate")
        ax.set_ylabel("Final titer, g/L")
        ax.set_title("Final titer vs lactate burden")
        ax.grid(True, linestyle=":", alpha=0.5)
        fig.tight_layout()
        p = fig_dir / "final_titer_vs_lactate.png"
        fig.savefig(p, dpi=300)
        print("[OK] saved:", p)

    cols_show = [
        "rank", "condition", "recommendation", "decision_score_0to1",
        "final_titer_g_L", "mean_qP_pg_cell_day", "final_viability_pct",
        "mean_rate_Glucose", "mean_rate_Lactate", "mean_rate_Glutamine",
        "mean_rate_Ammonia",
    ]
    cols_show = [c for c in cols_show if c in summary.columns]
    print("\n[DECISION SUMMARY]")
    print(summary[cols_show].to_string(index=False))


if __name__ == "__main__":
    main()
