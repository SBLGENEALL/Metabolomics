#!/usr/bin/env python
"""Prepare iCHO3K-ready multi-omics model input tables.

This script does not run FBA by itself. It converts experimental spent-media
rates and RNA-seq gene expression into tables that can be used by COBRApy,
RAVEN, COBRA Toolbox, or cameo once the iCHO3K SBML/MAT model is available.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import json

import numpy as np
import pandas as pd


def load_spent_media_rates(path: Path) -> pd.DataFrame:
    rates = pd.read_csv(path)
    required = {"metabolite", "day_start", "day_end", "apparent_rate_per_day"}
    missing = required.difference(rates.columns)
    if missing:
        raise ValueError(f"Missing rate columns: {sorted(missing)}")
    return rates


def make_exchange_constraints(rates: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    merged = rates.merge(mapping, left_on="metabolite", right_on="input_column", how="left", suffixes=("", "_map"))
    if merged["model_exchange_reaction_id"].isna().any():
        missing = sorted(merged.loc[merged["model_exchange_reaction_id"].isna(), "metabolite"].dropna().unique())
        print(f"Warning: no exchange mapping for {missing}")

    # Convention for most COBRA models: uptake is negative, secretion positive.
    flux = merged["qmet_per_1e6_cells_day"].where(
        merged["qmet_per_1e6_cells_day"].notna(), merged["apparent_rate_per_day"]
    )
    merged["measured_flux"] = flux
    merged["lower_bound"] = np.where(flux < 0, flux, 0.0)
    merged["upper_bound"] = np.where(flux > 0, flux, 0.0)
    merged.loc[merged["sign_convention"].eq("uptake_negative") & (flux > 0), "notes"] = (
        "Positive measured value for expected uptake; check concentration trend/sign"
    )
    cols = [
        "passage_or_clone",
        "producer_group",
        "replicate",
        "day_start",
        "day_end",
        "metabolite",
        "model_exchange_reaction_id",
        "measured_flux",
        "lower_bound",
        "upper_bound",
        "unit",
        "notes",
    ]
    return merged[[c for c in cols if c in merged.columns]]


def write_escher_reaction_data(constraints: pd.DataFrame, outdir: Path) -> None:
    if constraints.empty or "model_exchange_reaction_id" not in constraints.columns:
        return
    grouped = (
        constraints.dropna(subset=["model_exchange_reaction_id", "measured_flux"])
        .groupby("model_exchange_reaction_id", dropna=False)["measured_flux"]
        .mean()
        .reset_index()
        .rename(columns={"model_exchange_reaction_id": "reaction_id", "measured_flux": "mean_flux"})
    )
    grouped.to_csv(outdir / "escher_reaction_data_mean_flux.csv", index=False)
    reaction_data = dict(zip(grouped["reaction_id"], grouped["mean_flux"].astype(float)))
    (outdir / "escher_reaction_data_mean_flux.json").write_text(
        json.dumps(reaction_data, indent=2),
        encoding="utf-8",
    )
    html = """<!doctype html>
<html>
<head><meta charset="utf-8"><title>Escher overlay instructions</title></head>
<body style="font-family:Arial,sans-serif;max-width:880px;margin:40px auto;line-height:1.5">
<h1>iCHO3K Escher overlay data</h1>
<p>This folder contains <code>escher_reaction_data_mean_flux.json</code> and
<code>escher_reaction_data_mean_flux.csv</code>. Use the JSON as reaction data
in Escher or Escher Builder with an iCHO3K-compatible map.</p>
<ol>
<li>Open Escher Builder in the Python/conda environment that has <code>escher</code>.</li>
<li>Load or build an iCHO3K map for central carbon / amino acid exchange reactions.</li>
<li>Load <code>escher_reaction_data_mean_flux.json</code> as reaction data.</li>
</ol>
<p>Sign convention follows COBRA exchange flux: uptake is negative, secretion is positive.</p>
</body>
</html>
"""
    (outdir / "escher_overlay_instructions.html").write_text(html, encoding="utf-8")


def make_reaction_scores(transcriptomics: pd.DataFrame, gpr: pd.DataFrame) -> pd.DataFrame:
    if transcriptomics.empty:
        return pd.DataFrame()
    tx = transcriptomics.copy()
    tx["tpm"] = pd.to_numeric(tx["tpm"], errors="coerce")
    tx["log_tpm"] = np.log2(tx["tpm"].fillna(0) + 1)
    gene_summary = (
        tx.groupby(["condition", "timepoint", "gene_symbol"], dropna=False)
        .agg(mean_log_tpm=("log_tpm", "mean"), mean_tpm=("tpm", "mean"))
        .reset_index()
    )
    joined = gpr.merge(gene_summary, on="gene_symbol", how="left")
    # For exact GPR parsing, use the model's native GPR parser. This template gives a first pass.
    reaction_scores = (
        joined.groupby(["condition", "timepoint", "reaction_id", "reaction_name", "pathway"], dropna=False)
        .agg(
            n_mapped_genes=("gene_symbol", "nunique"),
            reaction_expression_score=("mean_log_tpm", "max"),
            mapped_genes=("gene_symbol", lambda s: ";".join(sorted(set(str(x) for x in s.dropna())))),
        )
        .reset_index()
    )
    return reaction_scores


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rates", required=True, type=Path, help="04_interval_rates.csv from spent_media_minimal_pipeline")
    parser.add_argument("--exchange-map", required=True, type=Path)
    parser.add_argument("--transcriptomics", type=Path)
    parser.add_argument("--gpr-map", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    rates = load_spent_media_rates(args.rates)
    exchange_map = pd.read_csv(args.exchange_map)
    constraints = make_exchange_constraints(rates, exchange_map)
    constraints.to_csv(args.outdir / "icho_exchange_constraints.csv", index=False)
    write_escher_reaction_data(constraints, args.outdir)

    if args.transcriptomics and args.gpr_map:
        transcriptomics = pd.read_csv(args.transcriptomics)
        gpr = pd.read_csv(args.gpr_map)
        scores = make_reaction_scores(transcriptomics, gpr)
        scores.to_csv(args.outdir / "reaction_expression_scores.csv", index=False)

    print(f"Wrote model-ready inputs to {args.outdir}")


if __name__ == "__main__":
    main()
