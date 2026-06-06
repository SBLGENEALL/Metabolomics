"""
13_rank_pathway_biomarkers.py

Create industrial decision-support candidate tables from separated evidence.

The output language is deliberately conservative:
- measured screening marker
- candidate pathway signature
- engineering target hypothesis
- demand-conditioned explanation

Product-demand-driven IgG reactions are excluded from predictive rankings.
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


PROVENANCE_COLUMNS = [
    "evidence_observability",
    "high_low_difference_source",
    "reconciliation_status",
    "decision_role",
    "interpretation_guardrail",
]


def read(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def percentile_score(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    result = pd.Series(np.nan, index=series.index, dtype=float)
    valid = values.notna()
    if valid.any():
        result.loc[valid] = 100.0 * values.loc[valid].rank(method="average", pct=True)
    return result


def mean_available(df: pd.DataFrame, columns: list) -> pd.Series:
    available = [col for col in columns if col in df.columns]
    if not available:
        return pd.Series(np.nan, index=df.index, dtype=float)
    return df[available].apply(pd.to_numeric, errors="coerce").mean(axis=1, skipna=True)


def ensure_columns(df: pd.DataFrame, columns: list) -> pd.DataFrame:
    out = df.copy()
    for col in columns:
        if col not in out.columns:
            out[col] = np.nan
    return out


def context_label(mode: object) -> str:
    text = str(mode).lower()
    if "no_igg" in text or "prediction" in text:
        return "no_igg_input_model_emergent"
    if "demand" in text or "explain" in text:
        return "measured_demand_conditioned"
    return "unspecified_model_context"


def select_predictive_context(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty or "mode" not in df.columns:
        return df
    predictive = df[df["mode"].astype(str).str.contains("no_igg|prediction", case=False, regex=True)]
    return predictive if not predictive.empty else df


def build_measured_markers(measured: pd.DataFrame) -> pd.DataFrame:
    if measured.empty:
        return pd.DataFrame()
    out = measured.copy()
    out["candidate_type"] = "measured_screening_marker"
    out["effect_magnitude"] = pd.to_numeric(out["standardized_high_low_effect"], errors="coerce")
    out["priority_score"] = percentile_score(out["effect_magnitude"])
    out["robustness_score"] = pd.to_numeric(out.get("robustness_score"), errors="coerce")
    out["evidence_coverage_score"] = np.where(out["robustness_score"].notna(), 100.0, 50.0)
    out["confidence_score"] = mean_available(
        out,
        ["evidence_coverage_score", "mapping_coverage_score", "robustness_score"],
    )
    out["hypothesis_ranking_score"] = out["priority_score"]
    out["confidence_basis"] = np.where(
        out["robustness_score"].notna(),
        "evidence coverage + mapping coverage + measured robustness",
        "evidence coverage + mapping coverage; reproducibility unavailable",
    )
    out["candidate_label"] = "candidate screening/feed-media marker"
    return out.sort_values(["priority_score", "confidence_score"], ascending=False)


def aggregate_model_pathways(separation: pd.DataFrame, demand: pd.DataFrame) -> pd.DataFrame:
    if separation.empty:
        return pd.DataFrame()
    use = separation.copy()
    use = use[use.get("comparison", "").eq("High_vs_Low")]
    use = use[use.get("high_low_difference_source", "").eq("model_emergent")]
    use = select_predictive_context(use)
    if use.empty:
        return pd.DataFrame()

    keys = ["chompact_pathway", "chompact_subpathway"]
    if "mode" in use.columns:
        keys.insert(0, "mode")
    fba = use[use["metric_family"].eq("FBA")].copy()
    fva = use[use["metric_family"].eq("FVA")].copy()

    fba_cols = keys + [
        "normalized_effect",
        "delta_a_minus_b",
        "abs_delta",
        "evidence_coverage_score",
        "mapping_coverage_score",
    ] + PROVENANCE_COLUMNS
    fva_cols = keys + [
        "fva_non_overlap_score",
        "fraction_reactions_non_overlap_gt_0_5",
        "robustness_score",
        "n_reactions_tested",
        "evidence_coverage_score",
        "mapping_coverage_score",
    ]
    fba = ensure_columns(fba, fba_cols)[fba_cols].rename(columns={
        "normalized_effect": "fba_normalized_effect",
        "delta_a_minus_b": "fba_signed_high_low_difference",
        "abs_delta": "fba_absolute_high_low_difference",
        "evidence_coverage_score": "fba_evidence_coverage_score",
        "mapping_coverage_score": "fba_mapping_coverage_score",
    })
    fva = ensure_columns(fva, fva_cols)[fva_cols].rename(columns={
        "evidence_coverage_score": "fva_evidence_coverage_score",
        "mapping_coverage_score": "fva_mapping_coverage_score",
        "robustness_score": "fva_robustness_score",
    })
    out = fba.merge(fva, on=keys, how="outer")

    if not demand.empty:
        demand_use = demand.copy()
        demand_use = demand_use[
            demand_use.get("comparison", "").eq("High_vs_Low")
            & demand_use.get("metric_family", "").isin(["FBA", "FVA"])
        ]
        if "high_low_difference_source" in demand_use.columns:
            demand_use = demand_use[demand_use["high_low_difference_source"].eq("model_emergent")]
        demand_summary = demand_use.groupby(
            ["chompact_pathway", "chompact_subpathway"],
            dropna=False,
        ).agg(
            demand_top10_frequency=("top10_frequency", "max"),
            demand_rank_iqr=("rank_iqr", "median"),
            n_demand_settings=("n_demand_settings", "max"),
        ).reset_index()
        out = out.merge(demand_summary, on=["chompact_pathway", "chompact_subpathway"], how="left")

    out["fba_priority_component"] = percentile_score(out.get("fba_normalized_effect", pd.Series(np.nan, index=out.index)))
    out["fva_priority_component"] = percentile_score(out.get("fva_non_overlap_score", pd.Series(np.nan, index=out.index)))
    out["priority_score"] = mean_available(out, ["fba_priority_component", "fva_priority_component"])
    out["hypothesis_ranking_score"] = out["priority_score"]
    out["evidence_coverage_score"] = 100.0 * out[
        ["fba_normalized_effect", "fva_non_overlap_score"]
    ].notna().mean(axis=1)
    out["mapping_coverage_score"] = mean_available(
        out,
        ["fba_mapping_coverage_score", "fva_mapping_coverage_score"],
    )
    out["robustness_score"] = pd.to_numeric(out.get("fva_robustness_score"), errors="coerce")
    if "n_demand_settings" in out.columns:
        stable = pd.to_numeric(out["n_demand_settings"], errors="coerce") > 1
        stability = np.where(
            stable,
            100.0
            * pd.to_numeric(out["demand_top10_frequency"], errors="coerce")
            / (1.0 + pd.to_numeric(out["demand_rank_iqr"], errors="coerce").clip(lower=0)),
            np.nan,
        )
        out["demand_stability_score"] = stability
    else:
        out["demand_stability_score"] = np.nan
    out["confidence_score"] = mean_available(
        out,
        [
            "evidence_coverage_score",
            "mapping_coverage_score",
            "robustness_score",
            "demand_stability_score",
        ],
    )
    out["candidate_type"] = "model_emergent_pathway_hypothesis"
    out["candidate_label"] = "candidate pathway signature"
    out["analysis_context"] = out["mode"].map(context_label) if "mode" in out.columns else "unspecified_model_context"
    for col, default in {
        "evidence_observability": "predicted_only",
        "high_low_difference_source": "model_emergent",
        "reconciliation_status": "untested",
        "decision_role": "engineering_hypothesis",
        "interpretation_guardrail": "Model-emergent iCHO3K hypothesis requiring independent validation.",
    }.items():
        if col not in out.columns:
            out[col] = default
        out[col] = out[col].fillna(default)
    out["confidence_basis"] = "available evidence coverage + mapping + FVA robustness + demand stability"
    return out.sort_values(["priority_score", "confidence_score"], ascending=False)


def build_demand_explanations(reaction_sep: pd.DataFrame) -> pd.DataFrame:
    if reaction_sep.empty:
        return pd.DataFrame()
    use = reaction_sep[
        reaction_sep.get("high_low_difference_source", "").eq("product_demand_driven")
    ].copy()
    if "mode" in use.columns:
        demand_mode = use[use["mode"].astype(str).str.contains("demand|explain", case=False, regex=True)]
        if not demand_mode.empty:
            use = demand_mode
    if use.empty:
        return use
    keys = ["reaction_id", "chompact_pathway", "chompact_subpathway"]
    if "mode" in use.columns:
        keys.insert(0, "mode")
    fba = use[use["metric_family"].eq("FBA")].copy()
    fva = use[use["metric_family"].eq("FVA")].copy()
    fba_keep = keys + ["normalized_effect", "delta_a_minus_b", "mapping_status"] + PROVENANCE_COLUMNS
    fva_keep = keys + ["fva_non_overlap_score", "normalized_effect"]
    fba = ensure_columns(fba, fba_keep)[fba_keep].rename(columns={
        "normalized_effect": "fba_normalized_effect",
        "delta_a_minus_b": "fba_signed_high_low_difference",
    })
    fva = ensure_columns(fva, fva_keep)[fva_keep].rename(columns={
        "normalized_effect": "fva_normalized_effect",
    })
    out = fba.merge(fva, on=keys, how="outer")
    out["effect_component"] = percentile_score(out["fba_normalized_effect"])
    out["fva_component"] = percentile_score(out["fva_non_overlap_score"])
    out["priority_score"] = mean_available(out, ["effect_component", "fva_component"])
    out["hypothesis_ranking_score"] = np.nan
    out["evidence_coverage_score"] = 100.0 * out[
        ["fba_normalized_effect", "fva_non_overlap_score"]
    ].notna().mean(axis=1)
    mapping_status = (
        out["mapping_status"]
        if "mapping_status" in out.columns
        else pd.Series("unknown", index=out.index)
    )
    out["mapping_coverage_score"] = np.where(mapping_status.eq("mapped"), 100.0, 0.0)
    out["robustness_score"] = 100.0 * pd.to_numeric(out["fva_non_overlap_score"], errors="coerce")
    out["confidence_score"] = mean_available(
        out,
        ["evidence_coverage_score", "mapping_coverage_score", "robustness_score"],
    )
    out["candidate_type"] = "demand_conditioned_explanation"
    out["candidate_label"] = "production-burden explanation"
    for col, default in {
        "evidence_observability": "predicted_only",
        "high_low_difference_source": "product_demand_driven",
        "reconciliation_status": "not_applicable",
        "decision_role": "do_not_rank_as_predictive",
        "interpretation_guardrail": (
            "Product-demand-driven reconstruction; explanatory only and excluded "
            "from independent predictive ranking."
        ),
    }.items():
        if col not in out.columns:
            out[col] = default
        out[col] = out[col].fillna(default)
    out["decision_role"] = "do_not_rank_as_predictive"
    out["predictive_ranking_eligible"] = False
    return out.sort_values(["priority_score", "confidence_score"], ascending=False)


def aggregate_model_reactions(reaction_sep: pd.DataFrame) -> pd.DataFrame:
    if reaction_sep.empty:
        return pd.DataFrame()
    use = reaction_sep[
        reaction_sep.get("high_low_difference_source", "").eq("model_emergent")
        & reaction_sep.get("comparison", "").eq("High_vs_Low")
    ].copy()
    use = select_predictive_context(use)
    if use.empty:
        return use
    keys = ["reaction_id", "chompact_pathway", "chompact_subpathway"]
    if "mode" in use.columns:
        keys.insert(0, "mode")
    fba = use[use["metric_family"].eq("FBA")].copy()
    fva = use[use["metric_family"].eq("FVA")].copy()
    fba_keep = keys + ["normalized_effect", "delta_a_minus_b", "abs_delta"] + PROVENANCE_COLUMNS
    fva_keep = keys + ["fva_non_overlap_score", "normalized_effect"]
    fba = ensure_columns(fba, fba_keep)[fba_keep].rename(columns={
        "normalized_effect": "fba_normalized_effect",
        "delta_a_minus_b": "fba_signed_high_low_difference",
        "abs_delta": "fba_absolute_high_low_difference",
    })
    fva = ensure_columns(fva, fva_keep)[fva_keep].rename(columns={
        "normalized_effect": "fva_normalized_effect",
    })
    out = fba.merge(fva, on=keys, how="outer")
    out["fba_priority_component"] = percentile_score(out["fba_normalized_effect"])
    out["fva_priority_component"] = percentile_score(out["fva_non_overlap_score"])
    out["priority_score"] = mean_available(out, ["fba_priority_component", "fva_priority_component"])
    out["hypothesis_ranking_score"] = out["priority_score"]
    out["evidence_coverage_score"] = 100.0 * out[["fba_normalized_effect", "fva_non_overlap_score"]].notna().mean(axis=1)
    out["mapping_coverage_score"] = 100.0
    out["robustness_score"] = 100.0 * pd.to_numeric(out["fva_non_overlap_score"], errors="coerce")
    out["confidence_score"] = mean_available(
        out,
        ["evidence_coverage_score", "mapping_coverage_score", "robustness_score"],
    )
    out["candidate_type"] = "engineering_target_hypothesis"
    out["candidate_label"] = "engineering target hypothesis"
    out["analysis_context"] = out["mode"].map(context_label) if "mode" in out.columns else "unspecified_model_context"
    for col, default in {
        "evidence_observability": "predicted_only",
        "high_low_difference_source": "model_emergent",
        "reconciliation_status": "untested",
        "decision_role": "engineering_hypothesis",
        "interpretation_guardrail": "Internal iCHO3K hypothesis requiring independent validation.",
    }.items():
        if col not in out.columns:
            out[col] = default
        out[col] = out[col].fillna(default)
    out["predictive_ranking_eligible"] = True
    return out.sort_values(["priority_score", "confidence_score"], ascending=False)


def pathway_measured_candidates(measured_markers: pd.DataFrame) -> pd.DataFrame:
    if measured_markers.empty:
        return pd.DataFrame()
    keys = ["chompact_pathway", "chompact_subpathway"]
    out = measured_markers.groupby(keys, dropna=False).agg(
        priority_score=("priority_score", "max"),
        hypothesis_ranking_score=("hypothesis_ranking_score", "max"),
        confidence_score=("confidence_score", "mean"),
        evidence_coverage_score=("evidence_coverage_score", "mean"),
        mapping_coverage_score=("mapping_coverage_score", "mean"),
        robustness_score=("robustness_score", "mean"),
        n_measured_reactions=("reaction_id", "nunique"),
        evidence_observability=("evidence_observability", "first"),
        high_low_difference_source=("high_low_difference_source", "first"),
        reconciliation_status=("reconciliation_status", "first"),
        decision_role=("decision_role", "first"),
        interpretation_guardrail=("interpretation_guardrail", "first"),
    ).reset_index()
    out["candidate_type"] = "measured_pathway_screening_signature"
    out["candidate_label"] = "candidate measured pathway signature"
    out["analysis_context"] = "measured_phenotype"
    return out


def classification_potential(fba_scores: pd.DataFrame) -> pd.DataFrame:
    if fba_scores.empty or "producer_group" not in fba_scores.columns:
        return pd.DataFrame()
    rows = []
    keys = ["chompact_pathway", "chompact_subpathway"]
    for key, grp in fba_scores.groupby(keys, dropna=False):
        key = key if isinstance(key, tuple) else (key,)
        piv = grp.pivot_table(index="producer_group", values="median_abs_flux", aggfunc="median")
        separation = np.nan
        if "High" in piv.index and "Low" in piv.index:
            high = float(piv.loc["High", "median_abs_flux"])
            low = float(piv.loc["Low", "median_abs_flux"])
            separation = abs(high - low) / (abs(high) + abs(low) + 1e-12)
        rows.append({
            "chompact_pathway": key[0],
            "chompact_subpathway": key[1],
            "n_groups_available": len(piv),
            "high_low_normalized_separation": separation,
            "classification_use": "exploratory_candidate_signature_not_validated_classifier",
        })
    return pd.DataFrame(rows).sort_values("high_low_normalized_separation", ascending=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Rank industrial CHOmpact candidate signatures and hypotheses")
    parser.add_argument("--dataset", default="practice_20aa", choices=["sowa2020", "practice_20aa", "own_experiment"])
    args = parser.parse_args()

    chompact_dir = os.path.join(results_dir(args.dataset, "tables"), "chompact")
    os.makedirs(chompact_dir, exist_ok=True)
    separation = read(os.path.join(chompact_dir, "chompact_pathway_high_low_separation.csv"))
    reaction_sep = read(os.path.join(chompact_dir, "chompact_reaction_high_low_separation.csv"))
    measured_effects = read(os.path.join(chompact_dir, "chompact_measured_high_low_effects.csv"))
    demand = read(os.path.join(chompact_dir, "chompact_robust_rank_across_demand_scales.csv"))
    fba_scores = read(os.path.join(chompact_dir, "chompact_pathway_fba_activity_scores.csv"))

    measured_markers = build_measured_markers(measured_effects)
    model_pathways = aggregate_model_pathways(separation, demand)
    demand_explanations = build_demand_explanations(reaction_sep)
    model_reactions = aggregate_model_reactions(reaction_sep)

    measured_reaction_candidates = measured_markers.copy()
    if not measured_reaction_candidates.empty:
        measured_reaction_candidates["predictive_ranking_eligible"] = True
    reaction_candidates = pd.concat(
        [measured_reaction_candidates, model_reactions],
        ignore_index=True,
        sort=False,
    )
    if not reaction_candidates.empty:
        reaction_candidates = reaction_candidates[
            ~reaction_candidates["high_low_difference_source"].eq("product_demand_driven")
        ].sort_values(["priority_score", "confidence_score"], ascending=False)

    measured_pathways = pathway_measured_candidates(measured_markers)
    pathway_candidates = pd.concat(
        [measured_pathways, model_pathways],
        ignore_index=True,
        sort=False,
    )
    if not pathway_candidates.empty:
        pathway_candidates = pathway_candidates[
            ~pathway_candidates["high_low_difference_source"].eq("product_demand_driven")
        ].sort_values(["priority_score", "confidence_score"], ascending=False)

    outputs = {
        "measured_screening_markers.csv": measured_markers,
        "model_emergent_pathway_hypotheses.csv": model_pathways,
        "demand_conditioned_explanations.csv": demand_explanations,
        "reaction_level_candidates.csv": reaction_candidates,
        "pathway_level_candidates.csv": pathway_candidates,
    }
    for name, df in outputs.items():
        df.to_csv(os.path.join(chompact_dir, name), index=False)

    # Backward-compatible filenames remain available, but use PR1 terminology
    # and the new hypothesis_ranking_score instead of composite biomarker score.
    ranked = pathway_candidates.copy()
    top = ranked.head(20).copy()
    evidence_cols = [
        "chompact_pathway",
        "chompact_subpathway",
        "candidate_type",
        "candidate_label",
        "hypothesis_ranking_score",
        "priority_score",
        "confidence_score",
        "evidence_coverage_score",
        "mapping_coverage_score",
        "robustness_score",
    ] + PROVENANCE_COLUMNS
    evidence = ensure_columns(ranked, evidence_cols)[evidence_cols] if not ranked.empty else pd.DataFrame(columns=evidence_cols)
    clf = classification_potential(fba_scores)

    ranked.to_csv(os.path.join(chompact_dir, "chompact_ranked_pathway_biomarkers.csv"), index=False)
    top.to_csv(os.path.join(chompact_dir, "chompact_top_pathway_biomarkers_for_report.csv"), index=False)
    evidence.to_csv(os.path.join(chompact_dir, "chompact_biomarker_evidence_matrix.csv"), index=False)
    clf.to_csv(os.path.join(chompact_dir, "chompact_clone_classification_potential.csv"), index=False)

    qc = pd.DataFrame([
        {
            "table": name,
            "n_rows": len(df),
            "n_product_demand_rows": int(
                df.get("high_low_difference_source", pd.Series(dtype=str))
                .eq("product_demand_driven")
                .sum()
            ),
        }
        for name, df in outputs.items()
    ])
    qc.to_csv(os.path.join(chompact_dir, "chompact_industrial_ranking_qc.csv"), index=False)

    print(f"[saved] {os.path.relpath(chompact_dir, ROOT)}")
    print(qc.to_string(index=False))
    if not model_pathways.empty:
        print("\nTop model-emergent candidate pathway signatures:")
        print(
            model_pathways[
                ["chompact_pathway", "chompact_subpathway", "priority_score", "confidence_score"]
            ].head(10).to_string(index=False)
        )


if __name__ == "__main__":
    main()
