"""
12_demand_sensitivity.py

Summarize pathway-ranking sensitivity across demand_scale strategies.

The step treats demand scaling as sensitivity analysis. If only one v1.0 run is
available, it exports a current-run baseline and marks multi-scale stability as
not assessed.
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


def load_runsheet(path: str) -> pd.DataFrame:
    if os.path.exists(path):
        return pd.read_csv(path)
    return pd.DataFrame(columns=["run_label", "dataset", "demand_scale", "chompact_dir"])


def summarize_single_run(dataset: str, chompact_dir: str) -> pd.DataFrame:
    sep = read(os.path.join(chompact_dir, "chompact_pathway_high_low_separation.csv"))
    if sep.empty:
        return pd.DataFrame()
    use = sep.copy()
    if "comparison" in use.columns:
        use = use[use["comparison"].isin(["High_vs_Low", "High_vs_Mother", "Mother_vs_Low"])]
    use["sensitivity_run"] = "current_run"
    use["dataset"] = dataset
    use["demand_scale"] = os.environ.get("CHO_DEMAND_SCALE", "current_or_auto")
    use["rank_metric"] = np.where(use.get("metric_family", "").eq("FVA"), "fva_non_overlap_score", "abs_delta")
    use["rank_value"] = np.where(
        use.get("metric_family", "").eq("FVA"),
        pd.to_numeric(use.get("fva_non_overlap_score", 0.0), errors="coerce").fillna(0.0),
        pd.to_numeric(use.get("abs_delta", 0.0), errors="coerce").fillna(0.0),
    )
    use["rank_within_run"] = use.groupby(["sensitivity_run", "comparison", "metric_family"])["rank_value"].rank(
        ascending=False,
        method="min",
    )
    return use


def summarize_runsheet(runsheet: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for _, row in runsheet.iterrows():
        cdir = str(row.get("chompact_dir", "")).strip()
        if not cdir:
            dataset = str(row.get("dataset", "practice_20aa"))
            cdir = os.path.join(results_dir(dataset, "tables"), "chompact")
        if not os.path.isabs(cdir):
            cdir = os.path.join(ROOT, cdir)
        df = summarize_single_run(str(row.get("dataset", "")), cdir)
        if df.empty:
            continue
        df["sensitivity_run"] = row.get("run_label", row.get("demand_scale", "run"))
        df["demand_scale"] = row.get("demand_scale", df["demand_scale"].iloc[0])
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def robust_rank(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    keys = ["chompact_pathway", "chompact_subpathway", "comparison", "metric_family"]
    aggregations = {
        "n_demand_settings": ("sensitivity_run", "nunique"),
        "median_rank": ("rank_within_run", "median"),
        "best_rank": ("rank_within_run", "min"),
        "worst_rank": ("rank_within_run", "max"),
        "rank_iqr": ("rank_within_run", lambda s: float(np.nanpercentile(s, 75) - np.nanpercentile(s, 25))),
        "median_rank_value": ("rank_value", "median"),
        "sign_consistency": (
            "delta_a_minus_b",
            lambda s: float((np.sign(s.dropna()) == np.sign(s.dropna()).mode().iloc[0]).mean())
            if len(s.dropna())
            else np.nan,
        ),
    }
    for col in PROVENANCE_COLUMNS:
        if col in df.columns:
            aggregations[col] = (
                col,
                lambda s: s.dropna().iloc[0]
                if s.dropna().nunique() == 1
                else ("mixed" if len(s.dropna()) else "unknown"),
            )
    out = df.groupby(keys, dropna=False).agg(**aggregations).reset_index()
    top10 = (
        df.assign(is_top10=df["rank_within_run"] <= 10)
        .groupby(keys, dropna=False)["is_top10"]
        .mean()
        .rename("top10_frequency")
        .reset_index()
    )
    out = out.merge(top10, on=keys, how="left")
    out["demand_stability_label"] = np.where(
        out["n_demand_settings"] <= 1,
        "single_setting_not_assessed",
        np.where((out["rank_iqr"] <= 5) & (out["sign_consistency"] >= 0.8), "stable", "sensitive"),
    )
    out = out.sort_values(["top10_frequency", "median_rank_value"], ascending=[False, False])
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Review demand_scale sensitivity for CHOmpact pathway rankings")
    parser.add_argument("--dataset", default="practice_20aa", choices=["sowa2020", "practice_20aa", "own_experiment"])
    parser.add_argument("--runsheet", default=os.path.join(ROOT, "data", "demand_sensitivity_runsheet.csv"))
    args = parser.parse_args()

    chompact_dir = os.path.join(results_dir(args.dataset, "tables"), "chompact")
    os.makedirs(chompact_dir, exist_ok=True)
    runsheet = load_runsheet(args.runsheet)
    sensitivity = summarize_runsheet(runsheet) if not runsheet.empty else summarize_single_run(args.dataset, chompact_dir)
    ranks = robust_rank(sensitivity)

    sensitivity.to_csv(os.path.join(chompact_dir, "chompact_demand_scale_sensitivity.csv"), index=False)
    ranks.to_csv(os.path.join(chompact_dir, "chompact_robust_rank_across_demand_scales.csv"), index=False)
    qc = pd.DataFrame([{
        "dataset": args.dataset,
        "runsheet": args.runsheet,
        "runsheet_found": os.path.exists(args.runsheet),
        "n_sensitivity_rows": len(sensitivity),
        "n_ranked_pathways": len(ranks),
        "interpretation": "demand_scale is sensitivity analysis, not fixed biological truth",
    }])
    qc.to_csv(os.path.join(chompact_dir, "chompact_demand_strategy_qc.csv"), index=False)
    print(f"[saved] {os.path.relpath(chompact_dir, ROOT)}")
    print(qc.to_string(index=False))


if __name__ == "__main__":
    main()
