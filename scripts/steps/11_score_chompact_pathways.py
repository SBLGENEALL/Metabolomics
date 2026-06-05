"""
11_score_chompact_pathways.py

Score CHOmpact pathway categories from mapped v1.0 iCHO3K outputs.

Scores keep measured qMet, model-predicted flux, and FVA feasibility separate.
They are pathway summaries for interpretation, not reconstructed CHOmpact fluxes.
"""
import argparse
import os
import sys
from itertools import combinations

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


def read(path: str) -> pd.DataFrame:
    return pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()


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
        source = "condition" if "condition" in out.columns else ("clone" if "clone" in out.columns else None)
        out["producer_group"] = out[source].map(infer_group) if source else "Unknown"
    return out


def numeric(df: pd.DataFrame, col: str, default=np.nan) -> pd.Series:
    if col in df.columns:
        return pd.to_numeric(df[col], errors="coerce")
    return pd.Series(default, index=df.index)


def pathway_keys() -> list:
    return ["chompact_pathway", "chompact_subpathway"]


def score_fba(flux: pd.DataFrame) -> pd.DataFrame:
    if flux.empty or "reaction_id" not in flux.columns:
        return pd.DataFrame()
    df = add_group(flux)
    if "flux" not in df.columns and "value" in df.columns:
        df = df.rename(columns={"value": "flux"})
    df["flux"] = numeric(df, "flux", 0.0).fillna(0.0)
    df["abs_flux"] = df["flux"].abs()
    df["constraint_fraction"] = df.groupby(pathway_keys())["is_constraint_reaction"].transform(
        lambda s: pd.to_numeric(s, errors="coerce").fillna(False).astype(bool).mean()
    )
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
        constraint_fraction=("constraint_fraction", "median"),
    ).reset_index()
    agg["source_type"] = "model_predicted_fba_pathway_score"
    agg["interpretation_boundary"] = np.where(
        agg["constraint_fraction"] >= 0.5,
        "mostly_constraint_driven",
        "mostly_model_emergent",
    )
    return agg


def score_fva(fva: pd.DataFrame) -> pd.DataFrame:
    if fva.empty or "reaction_id" not in fva.columns:
        return pd.DataFrame()
    df = add_group(fva)
    df["minimum"] = numeric(df, "minimum", 0.0).fillna(0.0)
    df["maximum"] = numeric(df, "maximum", 0.0).fillna(0.0)
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
    return agg


def interval_overlap(min_a: float, max_a: float, min_b: float, max_b: float) -> tuple:
    overlap = max(0.0, min(max_a, max_b) - max(min_a, min_b))
    union = max(max_a, max_b) - min(min_a, min_b)
    if union <= 1e-12:
        return overlap, union, 1.0, 0.0
    ratio = overlap / union
    return overlap, union, ratio, 1.0 - ratio


def separation_scores(fba_scores: pd.DataFrame, fva_scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if not fba_scores.empty:
        base_cols = pathway_keys()
        mode_cols = [c for c in ["mode", "objective"] if c in fba_scores.columns]
        for keys, grp in fba_scores.groupby(mode_cols + base_cols, dropna=False) if mode_cols else fba_scores.groupby(base_cols, dropna=False):
            if not isinstance(keys, tuple):
                keys = (keys,)
            key_map = dict(zip(mode_cols + base_cols, keys))
            by_group = grp.set_index("producer_group")
            for a, b in combinations([g for g in GROUP_ORDER if g in by_group.index], 2):
                av = float(by_group.loc[a, "median_abs_flux"])
                bv = float(by_group.loc[b, "median_abs_flux"])
                rows.append({
                    **key_map,
                    "comparison": f"{a}_vs_{b}",
                    "metric_family": "FBA",
                    "metric_name": "median_abs_flux",
                    "group_a": a,
                    "group_b": b,
                    "group_a_value": av,
                    "group_b_value": bv,
                    "delta_a_minus_b": av - bv,
                    "abs_delta": abs(av - bv),
                    "log2_ratio_a_over_b": np.log2((av + 1e-12) / (bv + 1e-12)),
                    "fva_overlap_ratio": np.nan,
                    "fva_non_overlap_score": np.nan,
                })
    if not fva_scores.empty:
        base_cols = pathway_keys()
        mode_cols = [c for c in ["mode", "objective"] if c in fva_scores.columns]
        grouped = fva_scores.groupby(mode_cols + base_cols, dropna=False) if mode_cols else fva_scores.groupby(base_cols, dropna=False)
        for keys, grp in grouped:
            if not isinstance(keys, tuple):
                keys = (keys,)
            key_map = dict(zip(mode_cols + base_cols, keys))
            by_group = grp.set_index("producer_group")
            for a, b in combinations([g for g in GROUP_ORDER if g in by_group.index], 2):
                amin = float(by_group.loc[a, "median_fva_min"])
                amax = float(by_group.loc[a, "median_fva_max"])
                bmin = float(by_group.loc[b, "median_fva_min"])
                bmax = float(by_group.loc[b, "median_fva_max"])
                ov, union, ratio, non_overlap = interval_overlap(amin, amax, bmin, bmax)
                rows.append({
                    **key_map,
                    "comparison": f"{a}_vs_{b}",
                    "metric_family": "FVA",
                    "metric_name": "median_feasible_interval",
                    "group_a": a,
                    "group_b": b,
                    "group_a_value": (amin + amax) / 2.0,
                    "group_b_value": (bmin + bmax) / 2.0,
                    "delta_a_minus_b": ((amin + amax) / 2.0) - ((bmin + bmax) / 2.0),
                    "abs_delta": abs(((amin + amax) / 2.0) - ((bmin + bmax) / 2.0)),
                    "overlap_width": ov,
                    "union_width": union,
                    "fva_overlap_ratio": ratio,
                    "fva_non_overlap_score": non_overlap,
                })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["source_type"] = np.where(out["metric_family"].eq("FVA"), "model_predicted_fva_overlap", "model_predicted_fba_separation")
    return out


def measured_qmet_scores(rates: pd.DataFrame) -> pd.DataFrame:
    if rates.empty:
        return pd.DataFrame(columns=[
            "producer_group", "chompact_pathway", "chompact_subpathway", "n_metabolites",
            "median_qmet", "median_abs_qmet", "source_type",
        ])
    df = add_group(rates)
    rate_col = "rate_mmol_gDCWh" if "rate_mmol_gDCWh" in df.columns else None
    if not rate_col:
        return pd.DataFrame()
    df[rate_col] = numeric(df, rate_col, 0.0).fillna(0.0)
    df["abs_qmet"] = df[rate_col].abs()
    out = df.groupby(["producer_group"] + pathway_keys(), dropna=False).agg(
        n_metabolites=("reaction_id", "nunique"),
        median_qmet=(rate_col, "median"),
        median_abs_qmet=("abs_qmet", "median"),
    ).reset_index()
    out["source_type"] = "measured_qmet_pathway_score"
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Score CHOmpact pathway categories from mapped iCHO3K outputs")
    parser.add_argument("--dataset", default="practice_20aa", choices=["sowa2020", "practice_20aa", "own_experiment"])
    args = parser.parse_args()

    chompact_dir = os.path.join(results_dir(args.dataset, "tables"), "chompact")
    os.makedirs(chompact_dir, exist_ok=True)
    flux = read(os.path.join(chompact_dir, "chompact_mapped_flux.csv"))
    fva = read(os.path.join(chompact_dir, "chompact_mapped_fva.csv"))
    rates = read(os.path.join(chompact_dir, "chompact_mapped_measured_rates.csv"))

    fba_scores = score_fba(flux)
    fva_scores = score_fva(fva)
    qmet = measured_qmet_scores(rates)
    sep = separation_scores(fba_scores, fva_scores)

    files = {
        "chompact_pathway_fba_activity_scores.csv": fba_scores,
        "chompact_pathway_fva_robustness_scores.csv": fva_scores,
        "chompact_pathway_measured_qmet_scores.csv": qmet,
        "chompact_pathway_high_low_separation.csv": sep,
    }
    for name, df in files.items():
        df.to_csv(os.path.join(chompact_dir, name), index=False)

    qc = pd.DataFrame([
        {"table": name, "n_rows": len(df), "n_pathways": df["chompact_pathway"].nunique() if "chompact_pathway" in df.columns and len(df) else 0}
        for name, df in files.items()
    ])
    qc.to_csv(os.path.join(chompact_dir, "chompact_pathway_score_qc.csv"), index=False)
    print(f"[saved] {os.path.relpath(chompact_dir, ROOT)}")
    print(qc.to_string(index=False))


if __name__ == "__main__":
    main()
