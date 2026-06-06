"""
14_make_chompact_figures.py

Create evidence-aware industrial decision-support figures and executive summary.
"""
import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402


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


EVIDENCE_COLORS = {
    "measured_pathway_screening_signature": "#0072B2",
    "model_emergent_pathway_hypothesis": "#D55E00",
    "demand_conditioned_explanation": "#8A8A8A",
}

PATHWAY_COLORS = {
    "TCA": "#0072B2",
    "PPP": "#009E73",
    "OXPHOS": "#CC79A7",
    "Glutamine": "#E69F00",
    "Exchange": "#56B4E9",
    "Other": "#777777",
}

EVIDENCE_MARKERS = {
    "model-emergent": "o",
    "measured": "^",
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


def compact_label(row: pd.Series) -> str:
    pathway = str(row.get("chompact_pathway", ""))
    subpathway = str(row.get("chompact_subpathway", ""))
    reaction = str(row.get("reaction_id", ""))
    if reaction and reaction != "nan":
        return f"{reaction} | {subpathway}"
    return f"{pathway} | {subpathway}"


def save_measured_markers(measured: pd.DataFrame, fig_dir: str) -> str:
    if measured.empty:
        return ""
    use = measured.copy()
    use["signed_high_low_difference"] = pd.to_numeric(
        use["signed_high_low_difference"],
        errors="coerce",
    )
    use = use[use["signed_high_low_difference"].notna()]
    if use.empty:
        return ""
    use["label"] = use.apply(compact_label, axis=1)
    use = use.sort_values("priority_score", ascending=False).head(15)
    use = use.sort_values("signed_high_low_difference")
    colors = np.where(use["signed_high_low_difference"] >= 0, "#0072B2", "#D55E00")
    fig, ax = plt.subplots(figsize=(8.8, max(4.5, 0.36 * len(use) + 1.4)))
    ax.barh(use["label"], use["signed_high_low_difference"], color=colors)
    ax.axvline(0, color="#222222", lw=1)
    ax.set_xlabel("Measured qMet difference: High - Low")
    ax.set_title("Figure 12. Measured Screening and Feed-Media Markers")
    ax.grid(True, axis="x", ls=":", alpha=0.45)
    fig.tight_layout()
    path = os.path.join(fig_dir, "Fig12_measured_screening_markers.png")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def save_model_robustness(model_pathways: pd.DataFrame, fig_dir: str) -> str:
    if model_pathways.empty:
        return ""
    use = model_pathways.copy()
    use = use[use["high_low_difference_source"].eq("model_emergent")]
    use["robustness_score"] = pd.to_numeric(use["robustness_score"], errors="coerce")
    use = use[use["robustness_score"].notna()]
    if use.empty:
        return ""
    use["label"] = use.apply(compact_label, axis=1)
    use = use.sort_values("robustness_score", ascending=False).head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.8, max(4.5, 0.36 * len(use) + 1.4)))
    ax.barh(use["label"], use["robustness_score"], color="#CC79A7")
    ax.set_xlabel("FVA robustness score (0-100)")
    ax.set_title("Figure 13. Model-Emergent Pathway Robustness")
    ax.set_xlim(0, 105)
    ax.grid(True, axis="x", ls=":", alpha=0.45)
    fig.tight_layout()
    path = os.path.join(fig_dir, "Fig13_model_emergent_fva_robustness.png")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def pathway_family(row: pd.Series) -> str:
    text = " ".join(
        str(row.get(col, ""))
        for col in ["chompact_pathway", "chompact_subpathway", "reaction_id"]
    ).lower()
    evidence_type = str(row.get("evidence_type", "")).lower()
    if evidence_type == "measured" or "exchange" in text or "uptake" in text:
        return "Exchange"
    if "pentose" in text or "ppp" in text:
        return "PPP"
    if "tca" in text or "citrate" in text or "anapler" in text:
        return "TCA"
    if any(token in text for token in ["oxidative phosphorylation", "complex i", "complex iii", "complex iv", "atp synthase"]):
        return "OXPHOS"
    if "glutamine" in text or "glutamate" in text:
        return "Glutamine"
    return "Other"


def evidence_type(row: pd.Series) -> str:
    candidate_type = str(row.get("candidate_type", "")).lower()
    observability = str(row.get("evidence_observability", "")).lower()
    if "measured" in candidate_type or observability == "measured":
        return "measured"
    return "model-emergent"


def marker_size(coverage: pd.Series) -> pd.Series:
    values = pd.to_numeric(coverage, errors="coerce").clip(lower=0, upper=100)
    return 45.0 + 1.55 * values


def select_annotation_rows(
    candidates: pd.DataFrame,
    top_overall: int = 4,
    top_model: int = 4,
    top_measured: int = 3,
) -> pd.DataFrame:
    selections = [
        candidates.nlargest(top_overall, "priority_score"),
        candidates[candidates["evidence_type"].eq("model-emergent")].nlargest(
            top_model,
            "priority_score",
        ),
        candidates[candidates["evidence_type"].eq("measured")].nlargest(
            top_measured,
            "priority_score",
        ),
    ]
    selected = pd.concat(selections, ignore_index=False)
    selected = selected.sort_values(
        ["priority_score", "confidence_score"],
        ascending=False,
        kind="stable",
    )
    return selected.drop_duplicates("label", keep="first")


def annotate_in_lanes(ax: plt.Axes, rows: pd.DataFrame) -> None:
    """Place selected labels in deterministic evidence-specific lanes."""
    lane_specs = {
        "model-emergent": {
            "x": 0.04,
            "start_y": 0.84,
            "step_y": -0.042,
            "ha": "left",
        },
        "measured": {
            "x": 0.96,
            "start_y": 0.61,
            "step_y": -0.042,
            "ha": "right",
        },
    }
    for evidence, group in rows.groupby("evidence_type", sort=False):
        spec = lane_specs.get(str(evidence), lane_specs["model-emergent"])
        ordered = group.sort_values(
            ["priority_score", "confidence_score"],
            ascending=False,
            kind="stable",
        )
        for lane_index, (_, row) in enumerate(ordered.iterrows()):
            color = PATHWAY_COLORS.get(
                str(row["pathway_family"]),
                PATHWAY_COLORS["Other"],
            )
            label = str(row.get("chompact_subpathway", row["label"]))
            ax.annotate(
                label,
                (float(row["priority_score"]), float(row["confidence_score"])),
                xytext=(
                    spec["x"],
                    spec["start_y"] + lane_index * spec["step_y"],
                ),
                textcoords="axes fraction",
                fontsize=7,
                ha=spec["ha"],
                va="center",
                bbox={
                    "boxstyle": "square,pad=0.12",
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.88,
                },
                arrowprops={
                    "arrowstyle": "-",
                    "color": color,
                    "lw": 0.65,
                    "alpha": 0.72,
                    "shrinkA": 3,
                    "shrinkB": 4,
                },
                annotation_clip=True,
            )


def metadata_caption(candidates: pd.DataFrame) -> str:
    fields = [
        ("fva_source", "FVA source"),
        ("fva_scope", "scope"),
        ("discovery_role", "role"),
    ]
    parts = []
    for column, label in fields:
        if column not in candidates.columns:
            value = "unavailable"
        else:
            values = sorted(
                {
                    str(value)
                    for value in candidates[column].dropna()
                    if str(value).strip()
                }
            )
            value = ", ".join(values) if values else "unavailable"
        parts.append(f"{label}: {value}")
    return " | ".join(parts)


def save_priority_confidence(
    candidates: pd.DataFrame,
    fig_dir: str,
    chompact_dir: str,
) -> str:
    if candidates.empty:
        return ""
    use = candidates.copy()
    use["priority_score"] = pd.to_numeric(use["priority_score"], errors="coerce")
    use["confidence_score"] = pd.to_numeric(use["confidence_score"], errors="coerce")
    use = use[use["priority_score"].notna() & use["confidence_score"].notna()]
    use = use[~use["high_low_difference_source"].eq("product_demand_driven")]
    if use.empty:
        return ""
    use = use.sort_values(
        ["priority_score", "confidence_score"],
        ascending=False,
        kind="stable",
    )
    use["label"] = use.apply(compact_label, axis=1)
    use["evidence_type"] = use.apply(evidence_type, axis=1)
    use["pathway_family"] = use.apply(pathway_family, axis=1)
    use["plot_size"] = marker_size(use["evidence_coverage_score"]).fillna(45.0)
    annotation_rows = select_annotation_rows(use)
    use["annotation_display"] = use.index.isin(annotation_rows.index)

    companion = use.rename(columns={
        "chompact_pathway": "pathway",
        "chompact_subpathway": "subpathway",
    })
    companion_columns = [
        "pathway",
        "subpathway",
        "priority_score",
        "confidence_score",
        "evidence_coverage_score",
        "mapping_coverage_score",
        "robustness_score",
        "evidence_type",
        "pathway_family",
        "fva_source",
        "fva_scope",
        "annotation_display",
    ]
    for col in companion_columns:
        if col not in companion.columns:
            companion[col] = np.nan
    companion[companion_columns].to_csv(
        os.path.join(chompact_dir, "Fig14_candidate_pathway_priority_confidence_data.csv"),
        index=False,
    )

    fig, ax = plt.subplots(figsize=(11.5, 7.2))
    for (evidence, family), group in use.groupby(
        ["evidence_type", "pathway_family"],
        dropna=False,
    ):
        color = PATHWAY_COLORS.get(str(family), PATHWAY_COLORS["Other"])
        ax.scatter(
            group["priority_score"],
            group["confidence_score"],
            s=group["plot_size"],
            color=color,
            marker=EVIDENCE_MARKERS.get(str(evidence), "o"),
            alpha=0.78,
            edgecolor="white",
            linewidth=0.8,
        )
    annotate_in_lanes(ax, annotation_rows)
    ax.axvline(60, color="#999999", lw=0.9, ls="--")
    ax.axhline(60, color="#999999", lw=0.9, ls="--")
    ax.set_xlim(0, 105)
    ax.set_ylim(0, 105)
    ax.set_xlabel("Priority score")
    ax.set_ylabel("Confidence score")
    ax.set_title(
        "Figure 14. Candidate Pathway Priority and Confidence",
        pad=28,
    )
    ax.text(
        0.5,
        1.012,
        metadata_caption(use),
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=8,
        color="#555555",
    )
    ax.grid(True, ls=":", alpha=0.35)
    evidence_handles = [
        Line2D(
            [0],
            [0],
            marker=marker,
            color="none",
            markerfacecolor="#555555",
            markeredgecolor="white",
            markersize=8,
            label=label,
        )
        for label, marker in EVIDENCE_MARKERS.items()
    ]
    pathway_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=color,
            markeredgecolor="white",
            markersize=8,
            label=family,
        )
        for family, color in PATHWAY_COLORS.items()
    ]
    coverage_handles = [
        ax.scatter([], [], s=float(marker_size(pd.Series([coverage])).iloc[0]), color="#999999", alpha=0.65, label=f"{coverage}%")
        for coverage in [25, 50, 100]
    ]
    evidence_legend = ax.legend(
        handles=evidence_handles,
        title="Evidence type",
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        fontsize=8,
    )
    ax.add_artist(evidence_legend)
    pathway_legend = ax.legend(
        handles=pathway_handles,
        title="Pathway family",
        loc="center left",
        bbox_to_anchor=(1.01, 0.53),
        fontsize=8,
    )
    ax.add_artist(pathway_legend)
    ax.legend(
        handles=coverage_handles,
        title="Evidence coverage",
        loc="lower left",
        bbox_to_anchor=(1.01, 0.0),
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0, 0.82, 1))
    path = os.path.join(fig_dir, "Fig14_candidate_pathway_priority_confidence.png")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def format_rows(df: pd.DataFrame, fields: list, limit: int = 6) -> list:
    if df.empty:
        return ["- No eligible result was available."]
    lines = []
    for _, row in df.head(limit).iterrows():
        label = compact_label(row)
        details = []
        for field, display in fields:
            value = row.get(field, np.nan)
            if pd.notna(value):
                if isinstance(value, (int, float, np.integer, np.floating)):
                    details.append(f"{display}={float(value):.1f}")
                else:
                    details.append(f"{display}={value}")
        lines.append(f"- {label}: " + ", ".join(details))
    return lines


def write_summary(
    dataset: str,
    chompact_dir: str,
    fig_paths: list,
    measured: pd.DataFrame,
    model: pd.DataFrame,
    demand: pd.DataFrame,
) -> None:
    measured_lines = format_rows(
        measured,
        [
            ("priority_score", "priority"),
            ("confidence_score", "confidence"),
            ("signed_high_low_difference", "High-Low qMet"),
        ],
    )
    model_lines = format_rows(
        model,
        [
            ("priority_score", "priority"),
            ("confidence_score", "confidence"),
            ("robustness_score", "FVA robustness"),
        ],
    )
    demand_lines = format_rows(
        demand,
        [
            ("priority_score", "explanation priority"),
            ("confidence_score", "confidence"),
        ],
    )
    rel_figs = [os.path.relpath(path, ROOT).replace("\\", "/") for path in fig_paths if path]
    text = f"""# CHOmpact v1.1 PR1 Industrial Decision Summary

Dataset: `{dataset}`

## A. Measured Phenotype

These are observed or measurement-derived screening/feed-media markers.

{chr(10).join(measured_lines)}

## B. Model-Emergent Oxidative and Metabolic Hypotheses

These are predicted-only iCHO3K hypotheses under measured constraints. They are
candidate pathway signatures and engineering target hypotheses, not measured
intracellular fluxes.

{chr(10).join(model_lines)}

## C. Demand-Conditioned Production-Burden Explanations

These results explain imposed measured IgG demand. They are excluded from
independent predictive ranking.

{chr(10).join(demand_lines)}

## Interpretation Warnings

- iCHO3K is the calculation engine; CHOmpact is an interpretation category layer.
- Exchange reactions are measured/exchange-constraint-driven markers.
- Internal ATP synthase, ETC, PDH, citrate synthase, TCA, and glutamine-catabolism signals are predicted-only model-emergent hypotheses.
- IgG demand, assembly, heavy-chain, and light-chain reactions are product-demand-driven explanatory outputs.
- Priority is not confidence. Missing evidence remains unavailable rather than being converted to zero.
- No output in this report is a validated biomarker or validated engineering target.

## Figures

{chr(10).join(f"- `{path}`" for path in rel_figs)}
"""
    with open(
        os.path.join(chompact_dir, "CHOmpact_v1_1_executive_summary.md"),
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(text)


def update_main_report(dataset: str, fig_paths: list) -> None:
    report_path = os.path.join(results_dir(dataset, "."), "REPORT_SUMMARY.md")
    if not os.path.exists(report_path):
        return
    start_marker = "<!-- PR1_FIGURES_START -->"
    end_marker = "<!-- PR1_FIGURES_END -->"
    with open(report_path, "r", encoding="utf-8") as handle:
        text = handle.read()
    if start_marker in text and end_marker in text:
        before = text.split(start_marker, 1)[0].rstrip()
        after = text.split(end_marker, 1)[1].lstrip()
        text = before + ("\n\n" + after if after else "")
    rel_figs = [
        os.path.relpath(path, ROOT).replace("\\", "/")
        for path in fig_paths
        if path
    ]
    section = [
        start_marker,
        "## PR1 evidence-aware decision figures",
        "",
        *(f"- `{path}`" for path in rel_figs),
        end_marker,
    ]
    text = text.rstrip() + "\n\n" + "\n".join(section) + "\n"
    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write(text)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create evidence-aware CHOmpact decision-support figures")
    parser.add_argument("--dataset", default="practice_20aa", choices=["sowa2020", "practice_20aa", "own_experiment"])
    args = parser.parse_args()

    setup_style()
    tables = results_dir(args.dataset, "tables")
    chompact_dir = os.path.join(tables, "chompact")
    fig_dir = results_dir(args.dataset, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(chompact_dir, exist_ok=True)

    measured = read(os.path.join(chompact_dir, "measured_screening_markers.csv"))
    model = read(os.path.join(chompact_dir, "model_emergent_pathway_hypotheses.csv"))
    demand = read(os.path.join(chompact_dir, "demand_conditioned_explanations.csv"))
    candidates = read(os.path.join(chompact_dir, "pathway_level_candidates.csv"))

    paths = [
        save_measured_markers(measured, fig_dir),
        save_model_robustness(model, fig_dir),
        save_priority_confidence(candidates, fig_dir, chompact_dir),
    ]
    write_summary(args.dataset, chompact_dir, paths, measured, model, demand)
    update_main_report(args.dataset, paths)
    print("[saved] CHOmpact PR1 decision-support figures:")
    for path in paths:
        if path:
            print(f"  {os.path.relpath(path, ROOT)}")


if __name__ == "__main__":
    main()
