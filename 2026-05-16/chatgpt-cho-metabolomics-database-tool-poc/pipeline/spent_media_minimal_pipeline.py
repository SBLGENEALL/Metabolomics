#!/usr/bin/env python
"""Minimal spent-media/extracellular metabolite pipeline for CHO studies.

Input can be either:
1) wide format: one row per sample, metabolite columns such as Glucose/lactate
2) long format: one row per sample-metabolite with columns metabolite,value

Required sample columns:
sample_id, passage_or_clone, day, replicate

Recommended optional columns:
producer_group, viable_cell_density_1e6_mL, viability_pct, titer_mg_L
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


SAMPLE_COLUMNS = [
    "sample_id",
    "passage_or_clone",
    "producer_group",
    "day",
    "replicate",
    "viable_cell_density_1e6_mL",
    "viability_pct",
    "titer_mg_L",
]

KNOWN_NON_METABOLITE = set(SAMPLE_COLUMNS + ["sample_type", "batch_id", "notes"])


def read_input(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".xlsx", ".xls"}:
        xl = pd.ExcelFile(path)
        sheet = "Experiment_Input" if "Experiment_Input" in xl.sheet_names else "Spent_Media_Input"
        df = pd.read_excel(path, sheet_name=sheet)
    else:
        df = pd.read_csv(path)

    if {"metabolite", "value"}.issubset(df.columns):
        long_df = df.copy()
    else:
        id_cols = [c for c in SAMPLE_COLUMNS if c in df.columns]
        required = {"sample_id", "passage_or_clone", "day", "replicate"}
        missing = required.difference(id_cols)
        if missing:
            raise ValueError(f"Missing required columns: {sorted(missing)}")
        metabolite_cols = [c for c in df.columns if c not in KNOWN_NON_METABOLITE]
        long_df = df.melt(
            id_vars=id_cols,
            value_vars=metabolite_cols,
            var_name="metabolite",
            value_name="value",
        )

    long_df["day"] = pd.to_numeric(long_df["day"], errors="coerce")
    long_df["replicate"] = pd.to_numeric(long_df["replicate"], errors="coerce")
    long_df["value"] = pd.to_numeric(long_df["value"], errors="coerce")
    for optional_numeric in ["viable_cell_density_1e6_mL", "viability_pct", "titer_mg_L"]:
        if optional_numeric in long_df.columns:
            long_df[optional_numeric] = pd.to_numeric(long_df[optional_numeric], errors="coerce")
    if "producer_group" not in long_df.columns:
        long_df["producer_group"] = long_df["passage_or_clone"].astype(str)
    return long_df.dropna(subset=["sample_id", "passage_or_clone", "day", "metabolite"])


def summarize_qc(long_df: pd.DataFrame) -> pd.DataFrame:
    qc = (
        long_df.groupby("metabolite")
        .agg(
            n=("value", "size"),
            n_missing=("value", lambda s: int(s.isna().sum())),
            mean_value=("value", "mean"),
            sd_value=("value", "std"),
        )
        .reset_index()
    )
    qc["missing_rate"] = qc["n_missing"] / qc["n"]
    qc["cv_overall"] = qc["sd_value"] / qc["mean_value"]
    qc["qc_flag"] = np.select(
        [qc["missing_rate"] > 0.2, qc["cv_overall"] > 1.5],
        ["High missing rate", "Very high overall CV"],
        default="Pass",
    )
    return qc


def baseline_delta(long_df: pd.DataFrame) -> pd.DataFrame:
    df = long_df.copy()
    keys = ["passage_or_clone", "replicate", "metabolite"]
    baseline = (
        df.sort_values("day")
        .groupby(keys, dropna=False)
        .first()["value"]
        .rename("day0_value")
        .reset_index()
    )
    df = df.merge(baseline, on=keys, how="left")
    df["delta_from_first_day"] = df["value"] - df["day0_value"]
    df["log2_ratio_from_first_day"] = np.log2((df["value"] + 1) / (df["day0_value"] + 1))
    return df


def interval_rates(delta_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_cols = ["passage_or_clone", "replicate", "metabolite"]
    for keys, sub in delta_df.sort_values("day").groupby(group_cols, dropna=False):
        sub = sub.sort_values("day")
        records = sub.to_dict("records")
        for prev, curr in zip(records[:-1], records[1:]):
            dt = curr["day"] - prev["day"]
            if not np.isfinite(dt) or dt <= 0:
                continue
            apparent_rate = (curr["value"] - prev["value"]) / dt
            vcd_prev = prev.get("viable_cell_density_1e6_mL", np.nan)
            vcd_curr = curr.get("viable_cell_density_1e6_mL", np.nan)
            vcd_values = [v for v in [vcd_prev, vcd_curr] if np.isfinite(v)]
            mean_vcd = float(np.mean(vcd_values)) if vcd_values else np.nan
            qmet = apparent_rate / mean_vcd if np.isfinite(mean_vcd) and mean_vcd > 0 else np.nan
            rows.append(
                {
                    "passage_or_clone": keys[0],
                    "replicate": keys[1],
                    "producer_group": curr.get("producer_group", prev.get("producer_group", "")),
                    "metabolite": keys[2],
                    "day_start": prev["day"],
                    "day_end": curr["day"],
                    "delta_value": curr["value"] - prev["value"],
                    "apparent_rate_per_day": apparent_rate,
                    "qmet_per_1e6_cells_day": qmet,
                }
            )
    return pd.DataFrame(rows)


def endpoint_contrast(
    delta_df: pd.DataFrame,
    reference_group: str | None = None,
    compare_group: str | None = None,
) -> pd.DataFrame:
    endpoint = (
        delta_df.sort_values("day")
        .groupby(["passage_or_clone", "replicate", "metabolite"], dropna=False)
        .tail(1)
    )
    summary = (
        endpoint.groupby(["producer_group", "metabolite"], dropna=False)
        .agg(
            n=("value", "size"),
            endpoint_mean=("value", "mean"),
            endpoint_sd=("value", "std"),
            delta_mean=("delta_from_first_day", "mean"),
            log2_ratio_mean=("log2_ratio_from_first_day", "mean"),
        )
        .reset_index()
    )

    groups = sorted(endpoint["producer_group"].dropna().astype(str).unique())
    if len(groups) >= 2:
        ref = reference_group or groups[0]
        comp = compare_group or groups[-1]
        missing_groups = [g for g in [ref, comp] if g not in groups]
        if missing_groups:
            raise ValueError(
                f"Requested group(s) not found: {missing_groups}. "
                f"Available producer_group values: {groups}"
            )
        pivot = summary.pivot(index="metabolite", columns="producer_group", values="log2_ratio_mean")
        contrast = pd.DataFrame(index=pivot.index)
        contrast[f"log2_ratio_mean_{ref}"] = pivot.get(ref)
        contrast[f"log2_ratio_mean_{comp}"] = pivot.get(comp)
        contrast[f"delta_log2_ratio_{comp}_vs_{ref}"] = pivot.get(comp) - pivot.get(ref)
        return summary, contrast.reset_index()
    return summary, pd.DataFrame()


def pathway_rollup(contrast: pd.DataFrame) -> pd.DataFrame:
    pathway_map = {
        "Glucose": "Carbon source/glycolysis",
        "glucose": "Carbon source/glycolysis",
        "lactate": "Overflow metabolism",
        "Lactate": "Overflow metabolism",
        "glutamine": "Glutamine metabolism",
        "Glutamine": "Glutamine metabolism",
        "glutamate": "Glutamine metabolism",
        "Glutamate": "Glutamine metabolism",
        "alanine": "Amino acid byproduct",
        "Alanine": "Amino acid byproduct",
        "pyruvate": "Glycolysis/TCA entry",
        "Pyruvate": "Glycolysis/TCA entry",
        "asparagine": "Amino acid consumption",
        "aspartate": "Amino acid/TCA anaplerosis",
        "citrate": "TCA-related extracellular signal",
        "acetate": "Short-chain acid/byproduct",
    }
    if contrast.empty:
        return contrast
    value_col = [c for c in contrast.columns if c.startswith("delta_log2_ratio")]
    if not value_col:
        return pd.DataFrame()
    value_col = value_col[0]
    tmp = contrast.copy()
    tmp["pathway_group"] = tmp["metabolite"].map(pathway_map).fillna("Other")
    return (
        tmp.groupby("pathway_group")
        .agg(
            n_metabolites=("metabolite", "count"),
            mean_delta_log2_ratio=(value_col, "mean"),
            metabolites=("metabolite", lambda s: ", ".join(s.astype(str))),
        )
        .reset_index()
        .sort_values("mean_delta_log2_ratio")
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--reference-group", help="Baseline group for endpoint contrast, e.g. Mother")
    parser.add_argument("--compare-group", help="Comparison group for endpoint contrast, e.g. High")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    long_df = read_input(args.input)
    qc = summarize_qc(long_df)
    delta = baseline_delta(long_df)
    rates = interval_rates(delta)
    summary, contrast = endpoint_contrast(delta, args.reference_group, args.compare_group)
    pathways = pathway_rollup(contrast)

    long_df.to_csv(args.outdir / "01_long_input.csv", index=False)
    qc.to_csv(args.outdir / "02_qc_summary.csv", index=False)
    delta.to_csv(args.outdir / "03_baseline_delta.csv", index=False)
    rates.to_csv(args.outdir / "04_interval_rates.csv", index=False)
    summary.to_csv(args.outdir / "05_endpoint_summary.csv", index=False)
    contrast.to_csv(args.outdir / "06_endpoint_contrast.csv", index=False)
    pathways.to_csv(args.outdir / "07_pathway_rollup.csv", index=False)

    print(f"Wrote minimal spent-media outputs to {args.outdir}")
    if not contrast.empty:
        print("\nEndpoint contrast preview:")
        print(contrast.head(12).to_string(index=False))


if __name__ == "__main__":
    main()
