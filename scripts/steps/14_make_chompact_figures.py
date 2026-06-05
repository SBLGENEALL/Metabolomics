"""
14_make_chompact_figures.py

Create publication-style CHOmpact interpretation figures from v1.1 tables.
"""
import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


def _repo_root() -> str:
    cur = os.path.dirname(os.path.abspath(__file__))
    while cur and cur != os.path.dirname(cur):
        if os.path.exists(os.path.join(cur, "src", "config.py")):
            return cur
        cur = os.path.dirname(cur)
    return os.getcwd()


ROOT = _repo_root()
sys.path.insert(0, ROOT)

from src.config import results_dir  # noqa: E402


COLORS = {
    "High": "#0072B2",
    "Mother": "#6A6A6A",
    "Moderate": "#009E73",
    "Mid": "#009E73",
    "Low": "#D55E00",
}


def read(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def setup_style() -> None:
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": "#333333",
        "axes.labelcolor": "#222222",
        "xtick.color": "#222222",
        "ytick.color": "#222222",
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "legend.frameon": False,
        "savefig.dpi": 300,
    })


def save_pathway_activity(fba: pd.DataFrame, fig_dir: str) -> str:
    if fba.empty:
        return ""
    use = fba.copy()
    if "mode" in use.columns and use["mode"].notna().any():
        preferred = "measured_demand" if "measured_demand" in set(use["mode"].astype(str)) else str(use["mode"].dropna().iloc[0])
        use = use[use["mode"].astype(str).eq(preferred)]
    top_paths = use.groupby("chompact_pathway")["median_abs_flux"].max().sort_values(ascending=False).head(12).index
    use = use[use["chompact_pathway"].isin(top_paths)]
    piv = use.pivot_table(index="chompact_pathway", columns="producer_group", values="median_abs_flux", aggfunc="median").fillna(0.0)
    piv = piv.loc[piv.max(axis=1).sort_values().index]
    fig, ax = plt.subplots(figsize=(8.5, max(4.2, 0.36 * len(piv) + 1.2)))
    y = np.arange(len(piv))
    width = 0.22
    groups = [g for g in ["High", "Mother", "Moderate", "Mid", "Low"] if g in piv.columns]
    for i, group in enumerate(groups):
        ax.barh(y + (i - (len(groups) - 1) / 2) * width, piv[group], height=width, color=COLORS.get(group, "#999999"), label=group)
    ax.set_yticks(y)
    ax.set_yticklabels(piv.index)
    ax.set_xlabel("Median absolute iCHO3K flux")
    ax.set_title("CHOmpact pathway activity from iCHO3K pFBA\nmodel-predicted, not measured qMet")
    ax.grid(True, axis="x", ls=":", alpha=0.45)
    ax.legend(loc="lower right")
    fig.tight_layout()
    path = os.path.join(fig_dir, "Fig12_chompact_pathway_activity.png")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def save_fva_overlap(sep: pd.DataFrame, fig_dir: str) -> str:
    if sep.empty:
        return ""
    use = sep[(sep.get("metric_family", "") == "FVA") & (sep.get("comparison", "") == "High_vs_Low")].copy()
    if use.empty:
        use = sep[sep.get("metric_family", "") == "FVA"].copy()
    if use.empty:
        return ""
    use["fva_non_overlap_score"] = pd.to_numeric(use["fva_non_overlap_score"], errors="coerce").fillna(0.0)
    use["label"] = use["chompact_pathway"].astype(str) + " | " + use["chompact_subpathway"].astype(str)
    top = use.sort_values("fva_non_overlap_score", ascending=False).head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.3, max(4.2, 0.34 * len(top) + 1.2)))
    ax.barh(top["label"], top["fva_non_overlap_score"], color="#CC79A7")
    ax.set_xlabel("FVA non-overlap score, High vs Low")
    ax.set_title("Robust feasible-range separation\nmodel-predicted FVA interval comparison")
    ax.set_xlim(0, max(1.0, float(top["fva_non_overlap_score"].max()) * 1.08))
    ax.grid(True, axis="x", ls=":", alpha=0.45)
    fig.tight_layout()
    path = os.path.join(fig_dir, "Fig13_chompact_fva_robustness.png")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def save_biomarker_ranking(ranked: pd.DataFrame, fig_dir: str) -> str:
    if ranked.empty:
        return ""
    use = ranked.head(15).copy().iloc[::-1]
    use["label"] = use["chompact_pathway"].astype(str) + " | " + use["chompact_subpathway"].astype(str)
    fig, ax = plt.subplots(figsize=(8.5, max(4.4, 0.36 * len(use) + 1.2)))
    ax.barh(use["label"], use["composite_biomarker_score"], color="#56B4E9")
    ax.set_xlabel("Composite pathway biomarker score")
    ax.set_title("CHOmpact pathway biomarker priority\nmeasured qMet + iCHO3K FBA/FVA evidence")
    ax.grid(True, axis="x", ls=":", alpha=0.45)
    fig.tight_layout()
    path = os.path.join(fig_dir, "Fig14_chompact_biomarker_ranking.png")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def write_summary(dataset: str, chompact_dir: str, fig_paths: list, ranked: pd.DataFrame) -> None:
    top_lines = []
    if not ranked.empty:
        for _, row in ranked.head(8).iterrows():
            top_lines.append(
                f"- {row['chompact_pathway']} / {row['chompact_subpathway']}: "
                f"score={row['composite_biomarker_score']:.3g}, evidence={row['primary_evidence_type']}"
            )
    else:
        top_lines.append("- No ranked pathway biomarkers were available.")
    rel_figs = [os.path.relpath(p, ROOT).replace("\\", "/") for p in fig_paths if p]
    text = f"""# CHOmpact v1.1 Executive Summary

Dataset: `{dataset}`

## Interpretation Boundary

- iCHO3K is the calculation engine.
- CHOmpact is used only as a pathway category and figure layer.
- Measured exchange-rate/qMet evidence, model-predicted pFBA flux, and FVA feasible intervals are reported separately.
- demand_scale is treated as sensitivity analysis, not a fixed biological truth.

## Top Pathway Biomarker Candidates

{chr(10).join(top_lines)}

## Figures

{chr(10).join(f'- `{p}`' for p in rel_figs)}
"""
    with open(os.path.join(chompact_dir, "CHOmpact_v1_1_executive_summary.md"), "w", encoding="utf-8") as handle:
        handle.write(text)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create CHOmpact interpretation figures")
    parser.add_argument("--dataset", default="practice_20aa", choices=["sowa2020", "practice_20aa", "own_experiment"])
    args = parser.parse_args()

    setup_style()
    tables = results_dir(args.dataset, "tables")
    chompact_dir = os.path.join(tables, "chompact")
    fig_dir = results_dir(args.dataset, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(chompact_dir, exist_ok=True)

    fba = read(os.path.join(chompact_dir, "chompact_pathway_fba_activity_scores.csv"))
    sep = read(os.path.join(chompact_dir, "chompact_pathway_high_low_separation.csv"))
    ranked = read(os.path.join(chompact_dir, "chompact_ranked_pathway_biomarkers.csv"))

    paths = [
        save_pathway_activity(fba, fig_dir),
        save_fva_overlap(sep, fig_dir),
        save_biomarker_ranking(ranked, fig_dir),
    ]
    write_summary(args.dataset, chompact_dir, paths, ranked)
    print("[saved] CHOmpact figures:")
    for path in paths:
        if path:
            print(f"  {os.path.relpath(path, ROOT)}")


if __name__ == "__main__":
    main()
