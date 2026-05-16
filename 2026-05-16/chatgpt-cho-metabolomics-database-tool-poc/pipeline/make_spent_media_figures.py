#!/usr/bin/env python
"""Create publication-style SVG figures from the minimal spent-media pipeline.

This intentionally uses only Python standard library + pandas/numpy because the
bundled runtime may not include matplotlib/seaborn.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


COLORS = {
    "Early_Passage": "#2563EB",
    "Late_Passage": "#DC2626",
    "High": "#2563EB",
    "Low": "#DC2626",
}
DEFAULT_COLORS = ["#2563EB", "#DC2626", "#059669", "#7C3AED", "#EA580C", "#0891B2"]


def svg_text(x, y, text, size=12, anchor="start", weight="normal", color="#111827", rotate=None):
    transform = f' transform="rotate({rotate} {x} {y})"' if rotate else ""
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial, sans-serif" font-size="{size}" '
        f'font-weight="{weight}" fill="{color}" text-anchor="{anchor}"{transform}>{text}</text>'
    )


def scale(values, out_min, out_max):
    finite = [v for v in values if np.isfinite(v)]
    if not finite:
        return lambda _: (out_min + out_max) / 2
    vmin, vmax = min(finite), max(finite)
    if math.isclose(vmin, vmax):
        return lambda _: (out_min + out_max) / 2
    return lambda v: out_min + (v - vmin) * (out_max - out_min) / (vmax - vmin)


def write_svg(path: Path, width: int, height: int, body: list[str]):
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#FFFFFF"/>',
        *body,
        "</svg>",
    ]
    path.write_text("\n".join(svg), encoding="utf-8")


def line_plot(long_df: pd.DataFrame, output: Path, metabolites: list[str]):
    metabolites = [m for m in metabolites if m in set(long_df["metabolite"])]
    if not metabolites:
        return
    width, height = 1120, 760
    left, right, top, bottom = 80, 30, 80, 80
    panel_w = (width - left - right - 30) / 2
    panel_h = (height - top - bottom - 40) / 2
    body = [
        svg_text(40, 38, "A. Spent-media time-course", 22, weight="bold"),
        svg_text(40, 60, "Group mean profiles from extracellular metabolite measurements", 12, color="#475569"),
    ]
    groups = list(long_df["producer_group"].dropna().unique())
    for idx, metabolite in enumerate(metabolites[:4]):
        px = left + (idx % 2) * (panel_w + 30)
        py = top + (idx // 2) * (panel_h + 40)
        sub = long_df[long_df["metabolite"] == metabolite]
        grouped = sub.groupby(["producer_group", "day"])["value"].mean().reset_index()
        x_map = scale(grouped["day"], px, px + panel_w)
        y_map = scale(grouped["value"], py + panel_h, py)
        body.append(f'<rect x="{px:.1f}" y="{py:.1f}" width="{panel_w:.1f}" height="{panel_h:.1f}" fill="#FFFFFF" stroke="#CBD5E1"/>')
        body.append(svg_text(px, py - 12, metabolite, 13, weight="bold"))
        for i in range(5):
            y = py + i * panel_h / 4
            body.append(f'<line x1="{px:.1f}" x2="{px+panel_w:.1f}" y1="{y:.1f}" y2="{y:.1f}" stroke="#E5E7EB"/>')
        for g_idx, group in enumerate(groups):
            g = grouped[grouped["producer_group"] == group].sort_values("day")
            if g.empty:
                continue
            color = COLORS.get(str(group), DEFAULT_COLORS[g_idx % len(DEFAULT_COLORS)])
            points = [(x_map(r.day), y_map(r.value)) for r in g.itertuples()]
            d = " ".join([f"{x:.1f},{y:.1f}" for x, y in points])
            body.append(f'<polyline points="{d}" fill="none" stroke="{color}" stroke-width="2.5"/>')
            for x, y in points:
                body.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3.5" fill="{color}"/>')
        days = sorted(grouped["day"].dropna().unique())
        for d in days[:: max(1, len(days) // 4)]:
            body.append(svg_text(x_map(d), py + panel_h + 18, str(int(d) if float(d).is_integer() else d), 10, anchor="middle", color="#475569"))
        body.append(svg_text(px + panel_w / 2, py + panel_h + 38, "Culture day", 11, anchor="middle", color="#475569"))
    lx = width - 230
    ly = 34
    for g_idx, group in enumerate(groups):
        color = COLORS.get(str(group), DEFAULT_COLORS[g_idx % len(DEFAULT_COLORS)])
        body.append(f'<line x1="{lx}" x2="{lx+26}" y1="{ly+g_idx*20}" y2="{ly+g_idx*20}" stroke="{color}" stroke-width="3"/>')
        body.append(svg_text(lx + 34, ly + 4 + g_idx * 20, str(group), 12, color="#334155"))
    write_svg(output, width, height, body)


def endpoint_bar(contrast: pd.DataFrame, output: Path):
    value_col = next((c for c in contrast.columns if c.startswith("delta_log2_ratio")), None)
    if value_col is None or contrast.empty:
        return
    df = contrast.sort_values(value_col)
    width, height = 920, max(420, 48 * len(df) + 120)
    left, right, top, bottom = 220, 70, 70, 50
    plot_w = width - left - right
    plot_h = height - top - bottom
    values = df[value_col].astype(float).to_numpy()
    max_abs = max(1e-9, np.nanmax(np.abs(values)))
    x_map = lambda v: left + plot_w / 2 + (v / max_abs) * (plot_w / 2)
    zero_x = x_map(0)
    row_h = plot_h / len(df)
    body = [
        svg_text(40, 34, "B. Endpoint contrast", 22, weight="bold"),
        svg_text(40, 56, value_col.replace("_", " "), 12, color="#475569"),
        f'<line x1="{zero_x:.1f}" x2="{zero_x:.1f}" y1="{top:.1f}" y2="{top+plot_h:.1f}" stroke="#111827" stroke-width="1"/>',
    ]
    for i, row in enumerate(df.itertuples()):
        y = top + i * row_h + row_h * 0.2
        val = getattr(row, value_col)
        color = "#DC2626" if val > 0 else "#2563EB"
        x0, x1 = zero_x, x_map(val)
        x = min(x0, x1)
        w = abs(x1 - x0)
        body.append(svg_text(left - 12, y + row_h * 0.45, str(row.metabolite), 12, anchor="end"))
        body.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{row_h*0.55:.1f}" fill="{color}" opacity="0.86"/>')
        body.append(svg_text(x1 + (6 if val >= 0 else -6), y + row_h * 0.42, f"{val:.2f}", 11, anchor=("start" if val >= 0 else "end"), color="#334155"))
    body.append(svg_text(left + plot_w / 2, height - 18, "Delta log2 endpoint ratio", 12, anchor="middle", color="#475569"))
    write_svg(output, width, height, body)


def heatmap(delta_df: pd.DataFrame, output: Path):
    top_metabolites = [
        "Glucose", "lactate", "glutamine", "glutamate", "alanine",
        "asparagine", "aspartate", "citrate", "acetate", "pyruvate",
    ]
    df = delta_df[delta_df["metabolite"].isin(top_metabolites)].copy()
    if df.empty:
        return
    summary = df.groupby(["producer_group", "metabolite"])["log2_ratio_from_first_day"].mean().reset_index()
    pivot = summary.pivot(index="metabolite", columns="producer_group", values="log2_ratio_from_first_day")
    pivot = pivot.reindex([m for m in top_metabolites if m in pivot.index])
    groups = list(pivot.columns)
    width, height = 620, 520
    left, top = 170, 80
    cell_w, cell_h = 150, 34
    vals = pivot.to_numpy(dtype=float).flatten()
    max_abs = max(1e-9, np.nanmax(np.abs(vals)))
    def color(v):
        if not np.isfinite(v):
            return "#F1F5F9"
        if v >= 0:
            strength = min(abs(v) / max_abs, 1)
            return f"rgb({255},{int(245-110*strength)},{int(235-180*strength)})"
        strength = min(abs(v) / max_abs, 1)
        return f"rgb({int(239-180*strength)},{int(246-120*strength)},{255})"
    body = [
        svg_text(40, 34, "C. Mean log2 change heatmap", 22, weight="bold"),
        svg_text(40, 56, "Relative to first sampling day within passage/replicate", 12, color="#475569"),
    ]
    for j, group in enumerate(groups):
        body.append(svg_text(left + j * cell_w + cell_w / 2, top - 14, str(group), 12, anchor="middle", weight="bold"))
    for i, metabolite in enumerate(pivot.index):
        y = top + i * cell_h
        body.append(svg_text(left - 12, y + 22, metabolite, 12, anchor="end"))
        for j, group in enumerate(groups):
            x = left + j * cell_w
            v = pivot.loc[metabolite, group]
            body.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{cell_w:.1f}" height="{cell_h:.1f}" fill="{color(v)}" stroke="#FFFFFF"/>')
            body.append(svg_text(x + cell_w / 2, y + 22, "" if not np.isfinite(v) else f"{v:.2f}", 11, anchor="middle", color="#111827"))
    write_svg(output, width, height, body)


def qmet_heatmap(rates: pd.DataFrame, output: Path):
    if "qmet_per_1e6_cells_day" not in rates.columns or rates.empty:
        return
    df = rates.copy()
    df["qmet_per_1e6_cells_day"] = pd.to_numeric(df["qmet_per_1e6_cells_day"], errors="coerce")
    df = df.dropna(subset=["qmet_per_1e6_cells_day"])
    if df.empty:
        return
    priority = [
        "Glucose", "lactate", "ammonia", "glutamine", "glutamate", "alanine",
        "asparagine", "aspartate", "arginine", "serine", "glycine", "cysteine",
        "leucine", "isoleucine", "valine", "pyruvate",
    ]
    present = [m for m in priority if m in set(df["metabolite"])]
    present += [m for m in sorted(set(df["metabolite"])) if m not in present]
    present = present[:18]
    pivot = (
        df[df["metabolite"].isin(present)]
        .groupby(["metabolite", "passage_or_clone"])["qmet_per_1e6_cells_day"]
        .mean()
        .reset_index()
        .pivot(index="metabolite", columns="passage_or_clone", values="qmet_per_1e6_cells_day")
        .reindex(present)
    )
    groups = list(pivot.columns)
    width = max(760, 190 + 92 * len(groups))
    height = max(560, 110 + 30 * len(pivot.index))
    left, top = 170, 86
    cell_w = min(92, (width - left - 40) / max(1, len(groups)))
    cell_h = 28
    vals = pivot.to_numpy(dtype=float).flatten()
    max_abs = max(1e-9, np.nanmax(np.abs(vals)))

    def color(v):
        if not np.isfinite(v):
            return "#F1F5F9"
        strength = min(abs(v) / max_abs, 1)
        if v >= 0:
            return f"rgb({255},{int(237-120*strength)},{int(213-160*strength)})"
        return f"rgb({int(219-160*strength)},{int(234-100*strength)},{254})"

    body = [
        svg_text(40, 34, "D. qMet uptake/secretion heatmap", 22, weight="bold"),
        svg_text(40, 56, "Mean interval rate normalized by VCD; uptake negative, secretion positive", 12, color="#475569"),
    ]
    for j, group in enumerate(groups):
        body.append(svg_text(left + j * cell_w + cell_w / 2, top - 12, str(group), 10, anchor="middle", weight="bold", rotate=-30))
    for i, metabolite in enumerate(pivot.index):
        y = top + i * cell_h
        body.append(svg_text(left - 12, y + 19, str(metabolite), 11, anchor="end"))
        for j, group in enumerate(groups):
            x = left + j * cell_w
            v = pivot.loc[metabolite, group]
            body.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{cell_w:.1f}" height="{cell_h:.1f}" fill="{color(v)}" stroke="#FFFFFF"/>')
            body.append(svg_text(x + cell_w / 2, y + 18, "" if not np.isfinite(v) else f"{v:.2g}", 9, anchor="middle"))
    write_svg(output, width, height, body)


def qmet_bar(rates: pd.DataFrame, output: Path):
    if "qmet_per_1e6_cells_day" not in rates.columns or rates.empty:
        return
    df = rates.copy()
    df["qmet_per_1e6_cells_day"] = pd.to_numeric(df["qmet_per_1e6_cells_day"], errors="coerce")
    df = df.dropna(subset=["qmet_per_1e6_cells_day"])
    if df.empty:
        return
    summary = (
        df.groupby("metabolite")["qmet_per_1e6_cells_day"]
        .mean()
        .reset_index()
        .assign(abs_flux=lambda x: x["qmet_per_1e6_cells_day"].abs())
        .sort_values("abs_flux", ascending=False)
        .head(16)
        .sort_values("qmet_per_1e6_cells_day")
    )
    width, height = 920, 660
    left, top, right, bottom = 210, 72, 70, 54
    plot_w, plot_h = width - left - right, height - top - bottom
    values = summary["qmet_per_1e6_cells_day"].to_numpy()
    max_abs = max(1e-9, np.nanmax(np.abs(values)))
    x_map = lambda v: left + plot_w / 2 + (v / max_abs) * (plot_w / 2)
    zero_x = x_map(0)
    row_h = plot_h / max(1, len(summary))
    body = [
        svg_text(40, 34, "E. Dominant exchange fluxes", 22, weight="bold"),
        svg_text(40, 56, "Largest mean qMet values across intervals and clones", 12, color="#475569"),
        f'<line x1="{zero_x:.1f}" x2="{zero_x:.1f}" y1="{top}" y2="{top+plot_h}" stroke="#111827"/>',
    ]
    for i, row in enumerate(summary.itertuples()):
        y = top + i * row_h + row_h * 0.18
        val = row.qmet_per_1e6_cells_day
        color = "#DC2626" if val > 0 else "#2563EB"
        x1 = x_map(val)
        x = min(zero_x, x1)
        body.append(svg_text(left - 12, y + row_h * 0.45, str(row.metabolite), 12, anchor="end"))
        body.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{abs(x1-zero_x):.1f}" height="{row_h*0.58:.1f}" fill="{color}" opacity="0.88"/>')
        body.append(svg_text(x1 + (6 if val >= 0 else -6), y + row_h * 0.43, f"{val:.3g}", 11, anchor=("start" if val >= 0 else "end")))
    body.append(svg_text(left + plot_w / 2, height - 18, "qMet per 1e6 cells per day", 12, anchor="middle", color="#475569"))
    write_svg(output, width, height, body)


def culture_profile(long_df: pd.DataFrame, output: Path):
    needed = ["viable_cell_density_1e6_mL", "viability_pct", "titer_mg_L"]
    present = [c for c in needed if c in long_df.columns]
    if not present:
        return
    sample_cols = ["sample_id", "producer_group", "day", *present]
    samples = long_df[sample_cols].drop_duplicates("sample_id")
    metric_labels = {
        "viable_cell_density_1e6_mL": "VCD (1e6 cells/mL)",
        "viability_pct": "Viability (%)",
        "titer_mg_L": "Titer (mg/L)",
    }
    width, height = 1120, 380
    left, top, right, bottom = 78, 76, 34, 62
    panel_w = (width - left - right - 36 * (len(present) - 1)) / len(present)
    panel_h = height - top - bottom
    body = [
        svg_text(40, 34, "F. Culture performance profile", 22, weight="bold"),
        svg_text(40, 56, "Mean process profiles by producer group", 12, color="#475569"),
    ]
    groups = list(samples["producer_group"].dropna().unique())
    for idx, metric in enumerate(present):
        px = left + idx * (panel_w + 36)
        py = top
        grouped = samples.groupby(["producer_group", "day"])[metric].mean().reset_index()
        grouped[metric] = pd.to_numeric(grouped[metric], errors="coerce")
        x_map = scale(grouped["day"], px, px + panel_w)
        y_map = scale(grouped[metric], py + panel_h, py)
        body.append(f'<rect x="{px:.1f}" y="{py:.1f}" width="{panel_w:.1f}" height="{panel_h:.1f}" fill="#FFFFFF" stroke="#CBD5E1"/>')
        body.append(svg_text(px, py - 12, metric_labels.get(metric, metric), 13, weight="bold"))
        for g_idx, group in enumerate(groups):
            g = grouped[grouped["producer_group"] == group].sort_values("day")
            color = COLORS.get(str(group), DEFAULT_COLORS[g_idx % len(DEFAULT_COLORS)])
            points = [(x_map(r.day), y_map(getattr(r, metric))) for r in g.itertuples() if np.isfinite(getattr(r, metric))]
            if len(points) < 1:
                continue
            body.append(f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in points)}" fill="none" stroke="{color}" stroke-width="2.5"/>')
            for x, y in points:
                body.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3.2" fill="{color}"/>')
        body.append(svg_text(px + panel_w / 2, py + panel_h + 38, "Culture day", 11, anchor="middle", color="#475569"))
    write_svg(output, width, height, body)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--figdir", type=Path)
    args = parser.parse_args()
    figdir = args.figdir or (args.outdir / "figures")
    figdir.mkdir(parents=True, exist_ok=True)

    long_df = pd.read_csv(args.outdir / "01_long_input.csv")
    delta_df = pd.read_csv(args.outdir / "03_baseline_delta.csv")
    rates = pd.read_csv(args.outdir / "04_interval_rates.csv")
    contrast = pd.read_csv(args.outdir / "06_endpoint_contrast.csv")

    line_plot(long_df, figdir / "figure_A_time_course.svg", ["Glucose", "lactate", "glutamine", "alanine"])
    endpoint_bar(contrast, figdir / "figure_B_endpoint_contrast.svg")
    heatmap(delta_df, figdir / "figure_C_log2_change_heatmap.svg")
    qmet_heatmap(rates, figdir / "figure_D_qmet_heatmap.svg")
    qmet_bar(rates, figdir / "figure_E_dominant_exchange_fluxes.svg")
    culture_profile(long_df, figdir / "figure_F_culture_profile.svg")
    print(f"Wrote figures to {figdir}")


if __name__ == "__main__":
    main()
