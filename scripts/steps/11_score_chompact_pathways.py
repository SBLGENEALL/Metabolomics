"""
11_score_chompact_pathways.py

Create evidence-aware CHOmpact pathway scores from mapped v1.0 outputs.

Measured qMet, model-emergent FBA/FVA, exchange-constraint-driven outputs, and
product-demand-driven explanations remain explicitly separated. CHOmpact is an
interpretation schema only; no CHOmpact flux is reconstructed.
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


GROUP_ORDER = ["High", "Mother", "Moderate", "Mid", "Low"]
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


def infer_group(value: object) -> str:
    text = str(value)
    low = text.lower()
    if "high" in low:
        return "High"
    if "mother" in low or "parent" in low:
        return "Mother"
    if "moderate" in low:
        return "Moderate"
    if "low" in low:
        return "Low"
    if "mid" in low:
        return "Mid"
    return text


def add_group(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if out.empty:
        return out
    if "producer_group" not in out.columns:
        source = next(
            (col for col in ["group", "condition", "clone"] if col in out.columns),
            None,
        )
        out["producer_group"] = out[source].map(infer_group) if source else "Unknown"
    return out


def numeric(df: pd.DataFrame, col: str, default=np.nan) -> pd.Series:
    if col in df.columns:
        return pd.to_numeric(df[col], errors="coerce")
    return pd.Series(default, index=df.index, dtype=float)


def pathway_keys() -> list:
    return ["chompact_pathway", "chompact_subpathway"]


def unique_or_mixed(series: pd.Series, missing: str = "unknown") -> str:
    values = [str(x) for x in series.dropna().unique() if str(x).strip()]
    if not values:
        return missing
    if len(values) == 1:
        return values[0]
    return "mixed"


def provenance_summary(df: pd.DataFrame, keys: list) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=keys + PROVENANCE_COLUMNS)
    aggregations = {
        col: (col, lambda s, col=col: unique_or_mixed(s, "untested" if col == "reconciliation_status" else "unknown"))
        for col in PROVENANCE_COLUMNS
        if col in df.columns
    }
    if not aggregations:
        return pd.DataFrame(columns=keys + PROVENANCE_COLUMNS)
    summary = df.groupby(keys, dropna=False).agg(**aggregations).reset_index()
    for col in PROVENANCE_COLUMNS:
        if col not in summary.columns:
            summary[col] = "unknown"
    return summary


def mapping_coverage(df: pd.DataFrame, keys: list) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=keys + ["n_reactions", "n_mapped_reactions", "mapping_coverage_score"])
    work = df.copy()
    if "mapping_status" not in work.columns:
        work["mapping_status"] = "mapped"
    total = work.groupby(keys, dropna=False)["reaction_id"].nunique().rename("n_reactions")
    mapped = (
        work[work["mapping_status"].eq("mapped")]
        .groupby(keys, dropna=False)["reaction_id"]
        .nunique()
        .rename("n_mapped_reactions")
    )
    out = pd.concat([total, mapped], axis=1).fillna({"n_mapped_reactions": 0}).reset_index()
    out["mapping_coverage_score"] = np.where(
        out["n_reactions"] > 0,
        100.0 * out["n_mapped_reactions"] / out["n_reactions"],
        np.nan,
    )
    return out


def enrich_pathway_table(table: pd.DataFrame, source: pd.DataFrame, keys: list) -> pd.DataFrame:
    if table.empty:
        return table
    out = table.merge(provenance_summary(source, keys), on=keys, how="left")
    coverage = mapping_coverage(source, keys).drop(columns=["n_reactions"], errors="ignore")
    out = out.merge(coverage, on=keys, how="left")
    out["evidence_coverage_score"] = 100.0
    return out


def score_fba(flux: pd.DataFrame) -> pd.DataFrame:
    if flux.empty or "reaction_id" not in flux.columns:
        return pd.DataFrame()
    df = add_group(flux)
    if "flux" not in df.columns and "value" in df.columns:
        df = df.rename(columns={"value": "flux"})
    df["flux"] = numeric(df, "flux")
    df = df[df["flux"].notna()].copy()
    df["abs_flux"] = df["flux"].abs()
    group_cols = ["producer_group"] + pathway_keys()
    for opt in ["mode", "objective"]:
        if opt in df.columns:
            group_cols.insert(1, opt)
    agg = df.groupby(group_cols, dropna=False).agg(
        n_reactions=("reaction_id", "nunique"),
        median_flux=("flux", "median"),
        median_abs_flux=("abs_flux", "median"),
        mean_abs_flux=("abs_flux", "mean"),
        max_abs_flux=("abs_flux", "max"),
    ).reset_index()
    agg["source_type"] = "model_predicted_fba_pathway_score"
    agg = enrich_pathway_table(agg, df, pathway_keys())
    return agg


def score_fva(fva: pd.DataFrame) -> pd.DataFrame:
    if fva.empty or "reaction_id" not in fva.columns:
        return pd.DataFrame()
    df = add_group(fva)
    df["minimum"] = numeric(df, "minimum")
    df["maximum"] = numeric(df, "maximum")
    df = df[df["minimum"].notna() & df["maximum"].notna()].copy()
    df["fva_range"] = (df["maximum"] - df["minimum"]).abs()
    df["fva_midpoint"] = (df["maximum"] + df["minimum"]) / 2.0
    df["fva_tightness"] = 1.0 / (1.0 + df["fva_range"])
    group_cols = ["producer_group"] + pathway_keys()
    for opt in ["mode", "objective"]:
        if opt in df.columns:
            group_cols.insert(1, opt)
    agg = df.groupby(group_cols, dropna=False).agg(
        n_reactions=("reaction_id", "nunique"),
        median_fva_min=("minimum", "median"),
        median_fva_max=("maximum", "median"),
        median_fva_range=("fva_range", "median"),
        mean_fva_range=("fva_range", "mean"),
        median_fva_midpoint=("fva_midpoint", "median"),
        median_fva_tightness=("fva_tightness", "median"),
    ).reset_index()
    agg["source_type"] = "model_predicted_fva_pathway_score"
    agg = enrich_pathway_table(agg, df, pathway_keys())
    return agg


def interval_overlap(min_a: float, max_a: float, min_b: float, max_b: float) -> tuple:
    overlap = max(0.0, min(max_a, max_b) - max(min_a, min_b))
    union = max(max_a, max_b) - min(min_a, min_b)
    if union <= 1e-12:
        return overlap, union, 1.0, 0.0
    ratio = overlap / union
    return overlap, union, ratio, 1.0 - ratio


def reaction_fba_separation(flux: pd.DataFrame) -> pd.DataFrame:
    if flux.empty or "reaction_id" not in flux.columns:
        return pd.DataFrame()
    df = add_group(flux)
    if "flux" not in df.columns and "value" in df.columns:
        df = df.rename(columns={"value": "flux"})
    df["flux"] = numeric(df, "flux")
    df = df[df["flux"].notna()].copy()
    key_cols = ["reaction_id"] + pathway_keys()
    for opt in ["mode", "objective"]:
        if opt in df.columns:
            key_cols.insert(0, opt)
    rows = []
    for key, grp in df.groupby(key_cols, dropna=False):
        key = key if isinstance(key, tuple) else (key,)
        meta = dict(zip(key_cols, key))
        by_group = grp.groupby("producer_group")["flux"].median()
        if "High" not in by_group.index or "Low" not in by_group.index:
            continue
        high = float(by_group["High"])
        low = float(by_group["Low"])
        row = {
            **meta,
            "comparison": "High_vs_Low",
            "metric_family": "FBA",
            "metric_name": "reaction_flux",
            "group_a": "High",
            "group_b": "Low",
            "group_a_value": high,
            "group_b_value": low,
            "delta_a_minus_b": high - low,
            "abs_delta": abs(high - low),
            "normalized_effect": abs(high - low) / (0.5 * (abs(high) + abs(low)) + 1e-12),
            "fva_overlap_ratio": np.nan,
            "fva_non_overlap_score": np.nan,
            "source_type": "model_predicted_fba_separation",
            "mapping_status": unique_or_mixed(grp["mapping_status"]) if "mapping_status" in grp.columns else "unknown",
        }
        for col in PROVENANCE_COLUMNS:
            row[col] = unique_or_mixed(grp[col]) if col in grp.columns else "unknown"
        rows.append(row)
    return pd.DataFrame(rows)


def reaction_fva_separation(fva: pd.DataFrame) -> pd.DataFrame:
    if fva.empty or "reaction_id" not in fva.columns:
        return pd.DataFrame()
    df = add_group(fva)
    df["minimum"] = numeric(df, "minimum")
    df["maximum"] = numeric(df, "maximum")
    df = df[df["minimum"].notna() & df["maximum"].notna()].copy()
    key_cols = ["reaction_id"] + pathway_keys()
    for opt in ["mode", "objective"]:
        if opt in df.columns:
            key_cols.insert(0, opt)
    rows = []
    for key, grp in df.groupby(key_cols, dropna=False):
        key = key if isinstance(key, tuple) else (key,)
        meta = dict(zip(key_cols, key))
        by_group = grp.groupby("producer_group")[["minimum", "maximum"]].median()
        if "High" not in by_group.index or "Low" not in by_group.index:
            continue
        hmin, hmax = map(float, by_group.loc["High", ["minimum", "maximum"]])
        lmin, lmax = map(float, by_group.loc["Low", ["minimum", "maximum"]])
        overlap, union, ratio, non_overlap = interval_overlap(hmin, hmax, lmin, lmax)
        hmid = (hmin + hmax) / 2.0
        lmid = (lmin + lmax) / 2.0
        row = {
            **meta,
            "comparison": "High_vs_Low",
            "metric_family": "FVA",
            "metric_name": "reaction_feasible_interval",
            "group_a": "High",
            "group_b": "Low",
            "group_a_value": hmid,
            "group_b_value": lmid,
            "delta_a_minus_b": hmid - lmid,
            "abs_delta": abs(hmid - lmid),
            "normalized_effect": abs(hmid - lmid) / (0.5 * (abs(hmid) + abs(lmid)) + 1e-12),
            "overlap_width": overlap,
            "union_width": union,
            "fva_overlap_ratio": ratio,
            "fva_non_overlap_score": non_overlap,
            "source_type": "model_predicted_fva_overlap",
            "mapping_status": unique_or_mixed(grp["mapping_status"]) if "mapping_status" in grp.columns else "unknown",
        }
        for col in PROVENANCE_COLUMNS:
            row[col] = unique_or_mixed(grp[col]) if col in grp.columns else "unknown"
        rows.append(row)
    return pd.DataFrame(rows)


def pathway_separation(reaction_sep: pd.DataFrame) -> pd.DataFrame:
    if reaction_sep.empty:
        return pd.DataFrame()
    keys = pathway_keys() + ["comparison", "metric_family"]
    for opt in ["mode", "objective"]:
        if opt in reaction_sep.columns:
            keys.insert(0, opt)
    agg = reaction_sep.groupby(keys, dropna=False).agg(
        metric_name=("metric_name", lambda s: "reaction_level_pathway_summary"),
        group_a=("group_a", "first"),
        group_b=("group_b", "first"),
        group_a_value=("group_a_value", "median"),
        group_b_value=("group_b_value", "median"),
        delta_a_minus_b=("delta_a_minus_b", "median"),
        abs_delta=("abs_delta", "median"),
        normalized_effect=("normalized_effect", "median"),
        n_reactions_tested=("reaction_id", "nunique"),
        n_mapped_reactions=(
            "mapping_status",
            lambda s: int((s == "mapped").sum()) if len(s) else 0,
        ),
        fva_overlap_ratio=("fva_overlap_ratio", "median"),
        fva_non_overlap_score=("fva_non_overlap_score", "median"),
        fraction_reactions_non_overlap_gt_0_5=(
            "fva_non_overlap_score",
            lambda s: float((pd.to_numeric(s, errors="coerce").dropna() > 0.5).mean())
            if len(pd.to_numeric(s, errors="coerce").dropna())
            else np.nan,
        ),
    ).reset_index()
    provenance = provenance_summary(reaction_sep, keys)
    agg = agg.merge(provenance, on=keys, how="left")
    agg["source_type"] = np.where(
        agg["metric_family"].eq("FVA"),
        "model_predicted_fva_pathway_separation",
        "model_predicted_fba_pathway_separation",
    )
    agg["evidence_coverage_score"] = 100.0
    agg["mapping_coverage_score"] = np.where(
        agg["n_reactions_tested"] > 0,
        100.0 * agg["n_mapped_reactions"] / agg["n_reactions_tested"],
        np.nan,
    )
    agg["robustness_score"] = np.where(
        agg["metric_family"].eq("FVA"),
        100.0 * pd.to_numeric(agg["fva_non_overlap_score"], errors="coerce"),
        np.nan,
    )
    return agg


def measured_qmet_scores(rates: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "producer_group",
        "chompact_pathway",
        "chompact_subpathway",
        "n_metabolites",
        "median_qmet",
        "median_abs_qmet",
        "source_type",
    ] + PROVENANCE_COLUMNS
    if rates.empty:
        return pd.DataFrame(columns=columns)
    df = add_group(rates)
    rate_col = "rate_mmol_gDCWh" if "rate_mmol_gDCWh" in df.columns else None
    if not rate_col:
        return pd.DataFrame(columns=columns)
    df[rate_col] = numeric(df, rate_col)
    df = df[df[rate_col].notna()].copy()
    df["abs_qmet"] = df[rate_col].abs()
    keys = ["producer_group"] + pathway_keys()
    out = df.groupby(keys, dropna=False).agg(
        n_metabolites=("reaction_id", "nunique"),
        median_qmet=(rate_col, "median"),
        median_abs_qmet=("abs_qmet", "median"),
    ).reset_index()
    out["source_type"] = "measured_qmet_pathway_score"
    out = enrich_pathway_table(out, df, pathway_keys())
    return out


def measured_high_low_effects(rates: pd.DataFrame) -> pd.DataFrame:
    if rates.empty or "reaction_id" not in rates.columns:
        return pd.DataFrame()
    df = add_group(rates)
    rate_col = "rate_mmol_gDCWh" if "rate_mmol_gDCWh" in df.columns else None
    if not rate_col:
        return pd.DataFrame()
    df[rate_col] = numeric(df, rate_col)
    df = df[df[rate_col].notna()].copy()
    keys = ["reaction_id"] + pathway_keys()
    rows = []
    for key, grp in df.groupby(keys, dropna=False):
        key = key if isinstance(key, tuple) else (key,)
        meta = dict(zip(keys, key))
        by_group = grp.groupby("producer_group")[rate_col]
        medians = by_group.median()
        if "High" not in medians.index or "Low" not in medians.index:
            continue
        high = float(medians["High"])
        low = float(medians["Low"])
        scale = float(grp[rate_col].abs().median())
        row = {
            **meta,
            "comparison": "High_vs_Low",
            "high_median_qmet": high,
            "low_median_qmet": low,
            "signed_high_low_difference": high - low,
            "absolute_high_low_difference": abs(high - low),
            "standardized_high_low_effect": abs(high - low) / (scale + 1e-12),
            "n_high": int((grp["producer_group"] == "High").sum()),
            "n_low": int((grp["producer_group"] == "Low").sum()),
            "source_type": "measured_qmet_high_low_effect",
            "evidence_coverage_score": 100.0,
            "mapping_coverage_score": 100.0 if grp.get("mapping_status", pd.Series(["mapped"])).eq("mapped").all() else 0.0,
            "robustness_score": np.nan,
        }
        for col in PROVENANCE_COLUMNS:
            row[col] = unique_or_mixed(grp[col]) if col in grp.columns else "unknown"
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Score evidence-aware CHOmpact pathway categories")
    parser.add_argument("--dataset", default="practice_20aa", choices=["sowa2020", "practice_20aa", "own_experiment"])
    args = parser.parse_args()

    chompact_dir = os.path.join(results_dir(args.dataset, "tables"), "chompact")
    os.makedirs(chompact_dir, exist_ok=True)
    flux = read(os.path.join(chompact_dir, "chompact_mapped_flux.csv"))
    fva = read(os.path.join(chompact_dir, "chompact_mapped_fva.csv"))
    rates = read(os.path.join(chompact_dir, "chompact_mapped_measured_rates.csv"))

    fba_scores = score_fba(flux)
    fva_scores = score_fva(fva)
    qmet_scores = measured_qmet_scores(rates)
    measured_effects = measured_high_low_effects(rates)
    reaction_sep = pd.concat(
        [reaction_fba_separation(flux), reaction_fva_separation(fva)],
        ignore_index=True,
        sort=False,
    )
    pathway_sep = pathway_separation(reaction_sep)

    files = {
        "chompact_pathway_fba_activity_scores.csv": fba_scores,
        "chompact_pathway_fva_robustness_scores.csv": fva_scores,
        "chompact_pathway_measured_qmet_scores.csv": qmet_scores,
        "chompact_measured_high_low_effects.csv": measured_effects,
        "chompact_reaction_high_low_separation.csv": reaction_sep,
        "chompact_pathway_high_low_separation.csv": pathway_sep,
    }
    for name, df in files.items():
        df.to_csv(os.path.join(chompact_dir, name), index=False)

    qc = pd.DataFrame([
        {
            "table": name,
            "n_rows": len(df),
            "n_pathways": df["chompact_pathway"].nunique()
            if "chompact_pathway" in df.columns and len(df)
            else 0,
        }
        for name, df in files.items()
    ])
    qc.to_csv(os.path.join(chompact_dir, "chompact_pathway_score_qc.csv"), index=False)
    print(f"[saved] {os.path.relpath(chompact_dir, ROOT)}")
    print(qc.to_string(index=False))


if __name__ == "__main__":
    main()
