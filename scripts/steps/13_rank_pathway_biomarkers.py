"""
13_rank_pathway_biomarkers.py

Prioritize CHOmpact pathway biomarkers for clone-productivity interpretation.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd


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


def read(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def zscore(series: pd.Series) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce").fillna(0.0)
    std = x.std(ddof=0)
    if std <= 1e-12:
        return pd.Series(0.0, index=x.index)
    return (x - x.mean()) / std


def rank_biomarkers(sep: pd.DataFrame, demand: pd.DataFrame, qmet: pd.DataFrame) -> pd.DataFrame:
    keys = ["chompact_pathway", "chompact_subpathway"]
    rows = []
    if not sep.empty:
        use = sep.copy()
        if "comparison" in use.columns:
            use = use[use["comparison"].eq("High_vs_Low")]
        agg = use.groupby(keys, dropna=False).agg(
            fba_separation_score=("abs_delta", "max"),
            fva_non_overlap_score=("fva_non_overlap_score", "max"),
            direction_high_minus_low=("delta_a_minus_b", "median"),
            n_model_metrics=("metric_name", "count"),
        ).reset_index()
        rows.append(agg)
    if rows:
        base = rows[0]
    else:
        base = pd.DataFrame(columns=keys)

    if not demand.empty:
        d = demand.groupby(keys, dropna=False).agg(
            demand_top10_frequency=("top10_frequency", "max"),
            demand_rank_iqr=("rank_iqr", "median"),
            n_demand_settings=("n_demand_settings", "max"),
        ).reset_index()
        base = base.merge(d, on=keys, how="outer")
    if not qmet.empty and "median_abs_qmet" in qmet.columns:
        q = qmet.groupby(keys, dropna=False).agg(
            measured_qmet_score=("median_abs_qmet", "max"),
            n_measured_metabolites=("n_metabolites", "max"),
        ).reset_index()
        base = base.merge(q, on=keys, how="outer")

    for col in [
        "fba_separation_score",
        "fva_non_overlap_score",
        "demand_top10_frequency",
        "measured_qmet_score",
        "demand_rank_iqr",
        "n_model_metrics",
        "n_demand_settings",
        "n_measured_metabolites",
    ]:
        if col not in base.columns:
            base[col] = 0.0
        base[col] = pd.to_numeric(base[col], errors="coerce").fillna(0.0)

    base["z_fba_separation"] = zscore(base["fba_separation_score"]).clip(lower=0)
    base["z_fva_non_overlap"] = zscore(base["fva_non_overlap_score"]).clip(lower=0)
    base["z_measured_qmet"] = zscore(base["measured_qmet_score"]).clip(lower=0)
    base["demand_stability_score"] = np.where(
        base["n_demand_settings"] > 1,
        base["demand_top10_frequency"] / (1.0 + base["demand_rank_iqr"].clip(lower=0)),
        0.0,
    )
    base["composite_biomarker_score"] = (
        0.25 * base["z_measured_qmet"]
        + 0.30 * base["z_fba_separation"]
        + 0.30 * base["z_fva_non_overlap"]
        + 0.15 * base["demand_stability_score"]
    )
    base["primary_evidence_type"] = np.select(
        [
            (base["n_measured_metabolites"] > 0) & (base["fva_non_overlap_score"] > 0),
            base["fva_non_overlap_score"] > 0,
            base["fba_separation_score"] > 0,
            base["n_measured_metabolites"] > 0,
        ],
        [
            "measured_qMet_plus_FVA_non_overlap",
            "model_FVA_non_overlap",
            "model_FBA_flux_separation",
            "measured_qMet_only",
        ],
        default="weak_or_unavailable",
    )
    base["robustness_label"] = np.select(
        [
            (base["fva_non_overlap_score"] >= 0.5) & (base["composite_biomarker_score"] > 0),
            (base["fva_non_overlap_score"] > 0) | (base["fba_separation_score"] > 0),
        ],
        ["robust_candidate", "exploratory_candidate"],
        default="low_priority",
    )
    return base.sort_values("composite_biomarker_score", ascending=False)


def classification_potential(fba_scores: pd.DataFrame) -> pd.DataFrame:
    if fba_scores.empty or "producer_group" not in fba_scores.columns:
        return pd.DataFrame()
    keys = ["chompact_pathway", "chompact_subpathway"]
    rows = []
    for key, grp in fba_scores.groupby(keys, dropna=False):
        if not isinstance(key, tuple):
            key = (key,)
        piv = grp.pivot_table(index="producer_group", values="median_abs_flux", aggfunc="median")
        groups = piv.index.tolist()
        if "High" in groups and "Low" in groups:
            high = float(piv.loc["High", "median_abs_flux"])
            low = float(piv.loc["Low", "median_abs_flux"])
            separation = abs(high - low) / (abs(high) + abs(low) + 1e-12)
        else:
            separation = np.nan
        rows.append({
            "chompact_pathway": key[0],
            "chompact_subpathway": key[1],
            "n_groups_available": len(groups),
            "high_low_normalized_separation": separation,
            "classification_use": "exploratory_feature_only_needs_validation",
        })
    return pd.DataFrame(rows).sort_values("high_low_normalized_separation", ascending=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Rank CHOmpact pathway biomarkers")
    parser.add_argument("--dataset", default="practice_20aa", choices=["sowa2020", "practice_20aa", "own_experiment"])
    args = parser.parse_args()

    chompact_dir = os.path.join(results_dir(args.dataset, "tables"), "chompact")
    os.makedirs(chompact_dir, exist_ok=True)
    sep = read(os.path.join(chompact_dir, "chompact_pathway_high_low_separation.csv"))
    demand = read(os.path.join(chompact_dir, "chompact_robust_rank_across_demand_scales.csv"))
    qmet = read(os.path.join(chompact_dir, "chompact_pathway_measured_qmet_scores.csv"))
    fba = read(os.path.join(chompact_dir, "chompact_pathway_fba_activity_scores.csv"))

    ranked = rank_biomarkers(sep, demand, qmet)
    top = ranked.head(20).copy()
    evidence = ranked[[
        "chompact_pathway",
        "chompact_subpathway",
        "composite_biomarker_score",
        "primary_evidence_type",
        "robustness_label",
        "fba_separation_score",
        "fva_non_overlap_score",
        "measured_qmet_score",
        "demand_stability_score",
    ]] if not ranked.empty else pd.DataFrame()
    clf = classification_potential(fba)

    ranked.to_csv(os.path.join(chompact_dir, "chompact_ranked_pathway_biomarkers.csv"), index=False)
    top.to_csv(os.path.join(chompact_dir, "chompact_top_pathway_biomarkers_for_report.csv"), index=False)
    evidence.to_csv(os.path.join(chompact_dir, "chompact_biomarker_evidence_matrix.csv"), index=False)
    clf.to_csv(os.path.join(chompact_dir, "chompact_clone_classification_potential.csv"), index=False)
    print(f"[saved] {os.path.relpath(chompact_dir, ROOT)}")
    if not top.empty:
        print(top[["chompact_pathway", "chompact_subpathway", "composite_biomarker_score", "primary_evidence_type"]].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
