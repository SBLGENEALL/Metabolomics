"""
10_map_to_chompact.py

Map validated v1.0 iCHO3K outputs onto CHOmpact-style pathway categories.

This step does not recompute flux and does not use CHOmpact as a model. It only
adds an interpretation schema to existing iCHO3K FBA/FVA outputs.
"""
import argparse
import os
import sys
from typing import Dict, List, Optional

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


REQUIRED_MAPPING_COLUMNS = [
    "reaction_id",
    "reaction_name",
    "subsystem",
    "chompact_pathway",
    "chompact_subpathway",
    "interpretation_note",
]

OPTIONAL_MAPPING_COLUMNS = [
    "reaction_class",
    "is_measured_exchange",
    "is_product_related",
    "is_constraint_reaction",
    "priority_for_figures",
]

REACTION_ID_ALIASES = ["rxn_id", "reaction", "Reaction ID", "reactionID", "id"]
MAPPING_VALUE_COLUMNS = REQUIRED_MAPPING_COLUMNS[1:] + OPTIONAL_MAPPING_COLUMNS
PROVENANCE_COLUMNS = [
    "evidence_observability",
    "high_low_difference_source",
    "reconciliation_status",
    "decision_role",
    "interpretation_guardrail",
]

PRODUCT_DEMAND_REACTIONS = {
    "DM_igg_g",
    "DM_for_igg",
    "igg_formation",
    "igg_hc",
    "igg_lc",
}

FEED_MEDIA_EXCHANGES = {
    "EX_glc_e",
    "EX_lac_L_e",
    "EX_gln_L_e",
    "EX_glu_L_e",
    "EX_nh4_e",
}


def read_csv_if_exists(path: Optional[str]) -> pd.DataFrame:
    if path and os.path.exists(path):
        try:
            return pd.read_csv(path)
        except pd.errors.EmptyDataError:
            return pd.DataFrame()
    return pd.DataFrame()


def collapse_duplicate_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse duplicate column labels by taking the first non-null value.

    Some upstream v1.0 tables may carry both `reaction_id` and aliases such as
    `rxn_id`, or duplicate labels after previous merges. Pandas refuses to merge
    when a merge key is not unique, so normalize tables defensively here.
    """
    if df.empty or not df.columns.has_duplicates:
        return df

    data = {}
    for col in dict.fromkeys(df.columns):
        same = df.loc[:, df.columns == col]
        if same.shape[1] == 1:
            data[col] = same.iloc[:, 0]
        else:
            data[col] = same.bfill(axis=1).iloc[:, 0]
    return pd.DataFrame(data)


def normalize_reaction_id(df: pd.DataFrame) -> pd.DataFrame:
    out = collapse_duplicate_columns(df.copy())

    if "reaction_id" not in out.columns:
        for alias in REACTION_ID_ALIASES:
            if alias in out.columns:
                out = out.rename(columns={alias: "reaction_id"})
                break
    else:
        # If aliases also exist, use them only to fill missing reaction_id values,
        # then drop them to keep a single unambiguous merge key.
        rid = out["reaction_id"].astype("string")
        missing = rid.isna() | rid.fillna("").str.strip().eq("")
        for alias in REACTION_ID_ALIASES:
            if alias in out.columns:
                alias_values = out[alias].astype("string")
                rid = rid.mask(missing, alias_values)
                missing = rid.isna() | rid.fillna("").str.strip().eq("")
        out["reaction_id"] = rid

    drop_aliases = [c for c in REACTION_ID_ALIASES if c in out.columns and c != "reaction_id"]
    if drop_aliases:
        out = out.drop(columns=drop_aliases)

    out = collapse_duplicate_columns(out)
    if "reaction_id" in out.columns:
        out["reaction_id"] = out["reaction_id"].astype(str)
    return out


def infer_reaction_class(reaction_id: str) -> str:
    rid = str(reaction_id)
    if rid.startswith("EX_"):
        return "exchange"
    if rid.startswith("DM_"):
        return "demand"
    if "igg" in rid.lower():
        return "product"
    return "internal"


def evidence_provenance(reaction_id: str, source_type: str, reaction_class: str) -> dict:
    """Return conservative industrial interpretation labels for one reaction.

    Observability describes the biological quantity represented by the row,
    while source_type still records which pipeline table produced the value.
    Measured exchange reactions remain screening/feed markers even when the row
    contains a model flux constrained by the measured exchange phenotype.
    """
    rid = str(reaction_id)
    rclass = str(reaction_class or infer_reaction_class(rid))
    is_product = rid in PRODUCT_DEMAND_REACTIONS or "igg" in rid.lower()
    is_exchange = rclass == "exchange" or rid.startswith("EX_")

    if is_product:
        return {
            "evidence_observability": "predicted_only",
            "high_low_difference_source": "product_demand_driven",
            "reconciliation_status": "not_applicable",
            "decision_role": "do_not_rank_as_predictive",
            "interpretation_guardrail": (
                "Product-demand-driven reconstruction; explanatory only and "
                "excluded from independent predictive ranking."
            ),
        }

    if is_exchange:
        role = "feed_media_marker" if rid in FEED_MEDIA_EXCHANGES else "clone_selection_marker"
        return {
            "evidence_observability": "measured",
            "high_low_difference_source": "exchange_constraint_driven",
            "reconciliation_status": "not_applicable",
            "decision_role": role,
            "interpretation_guardrail": (
                "Measured or measurement-derived exchange phenotype; use for "
                "screening and feed/media triage, not as a model-emergent mechanism."
            ),
        }

    return {
        "evidence_observability": "predicted_only",
        "high_low_difference_source": "model_emergent",
        "reconciliation_status": "untested",
        "decision_role": "engineering_hypothesis",
        "interpretation_guardrail": (
            "Internal iCHO3K prediction under measured constraints; treat as an "
            "engineering target hypothesis requiring independent validation."
        ),
    }


def load_mapping(path: str) -> pd.DataFrame:
    mapping = collapse_duplicate_columns(pd.read_csv(path))
    mapping = normalize_reaction_id(mapping)
    missing = [c for c in REQUIRED_MAPPING_COLUMNS if c not in mapping.columns]
    if missing:
        raise ValueError(f"Mapping file is missing required columns: {missing}")
    for col in OPTIONAL_MAPPING_COLUMNS:
        if col not in mapping.columns:
            mapping[col] = ""
    mapping["reaction_id"] = mapping["reaction_id"].astype(str)
    return mapping.drop_duplicates("reaction_id")


def attach_mapping(df: pd.DataFrame, mapping: pd.DataFrame, source_type: str) -> pd.DataFrame:
    if df.empty:
        return df
    out = normalize_reaction_id(df)
    if "reaction_id" not in out.columns:
        raise ValueError(f"{source_type} table does not contain rxn_id/reaction_id")

    # Avoid carrying stale mapping columns from prior runs into the new merge.
    stale_mapping_cols = [c for c in MAPPING_VALUE_COLUMNS if c in out.columns]
    if stale_mapping_cols:
        out = out.drop(columns=stale_mapping_cols)

    out = out.merge(mapping, on="reaction_id", how="left", suffixes=("", "_chompact"))
    out = collapse_duplicate_columns(out)
    out["source_type"] = source_type
    out["mapping_status"] = out["chompact_pathway"].notna().map({True: "mapped", False: "unmapped"})
    out["chompact_pathway"] = out["chompact_pathway"].fillna("Other / unmapped")
    out["chompact_subpathway"] = out["chompact_subpathway"].fillna("Unmapped")
    out["interpretation_note"] = out["interpretation_note"].fillna(
        "No CHOmpact pathway category assigned; keep as iCHO3K-specific result."
    )
    if "reaction_class" not in out.columns:
        out["reaction_class"] = ""
    out["reaction_class"] = out["reaction_class"].fillna("")
    empty_class = out["reaction_class"].eq("")
    out.loc[empty_class, "reaction_class"] = out.loc[empty_class, "reaction_id"].map(infer_reaction_class)
    for col in ["is_measured_exchange", "is_product_related", "is_constraint_reaction"]:
        if col not in out.columns:
            out[col] = False
        out[col] = out[col].fillna(False)
    if "priority_for_figures" not in out.columns:
        out["priority_for_figures"] = 3
    out["priority_for_figures"] = pd.to_numeric(out["priority_for_figures"], errors="coerce").fillna(3).astype(int)
    provenance = out.apply(
        lambda row: evidence_provenance(
            row["reaction_id"],
            source_type,
            row.get("reaction_class", ""),
        ),
        axis=1,
        result_type="expand",
    )
    for col in PROVENANCE_COLUMNS:
        out[col] = provenance[col]
    return out


def find_first(paths: List[str]) -> Optional[str]:
    for path in paths:
        if os.path.exists(path):
            return path
    return None


def make_qc(mapped: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for name, df in mapped.items():
        if df.empty:
            rows.append({
                "table": name,
                "n_rows": 0,
                "n_reactions": 0,
                "n_mapped_reactions": 0,
                "mapping_rate_reactions": 0.0,
            })
            continue
        df = normalize_reaction_id(df)
        unique = df.drop_duplicates("reaction_id")
        n_rxn = len(unique)
        n_mapped = int((unique["mapping_status"] == "mapped").sum())
        rows.append({
            "table": name,
            "n_rows": len(df),
            "n_reactions": n_rxn,
            "n_mapped_reactions": n_mapped,
            "mapping_rate_reactions": n_mapped / n_rxn if n_rxn else 0.0,
        })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Map v1.0 iCHO3K outputs to CHOmpact pathway categories")
    parser.add_argument("--dataset", default="practice_20aa", choices=["sowa2020", "practice_20aa", "own_experiment"])
    parser.add_argument("--mapping", default=os.path.join(ROOT, "data", "chompact_pathway_mapping.csv"))
    parser.add_argument("--flux_table", default=None)
    parser.add_argument("--fva_table", default=None)
    parser.add_argument("--full_fva_scope", default=os.environ.get("CHO_FVA_SCOPE", "all"))
    args = parser.parse_args()

    tables = results_dir(args.dataset, "tables")
    chompact_dir = os.path.join(tables, "chompact")
    os.makedirs(chompact_dir, exist_ok=True)

    mapping = load_mapping(args.mapping)
    flux_path = args.flux_table or find_first([
        os.path.join(tables, "central_mab_flux_state.csv"),
        os.path.join(tables, "fba_results.csv"),
    ])
    focused_fva_path = args.fva_table or find_first([
        os.path.join(tables, "central_mab_fva_report.csv"),
        os.path.join(tables, "full_fva", f"full_fva_{args.full_fva_scope}_combined_report.csv"),
        os.path.join(tables, "full_fva", "full_fva_all_combined_report.csv"),
    ])
    overlap_path = find_first([
        os.path.join(tables, "full_fva", f"full_fva_{args.full_fva_scope}_high_low_overlap.csv"),
        os.path.join(tables, "full_fva", "full_fva_all_high_low_overlap.csv"),
        os.path.join(tables, "central_mab_fva_high_low_overlap.csv"),
    ])
    rate_path = os.path.join(tables, "exchange_rates.csv")

    mapped_flux = attach_mapping(read_csv_if_exists(flux_path), mapping, "model_predicted_flux")
    mapped_fva = attach_mapping(read_csv_if_exists(focused_fva_path), mapping, "model_predicted_fva")
    mapped_overlap = attach_mapping(read_csv_if_exists(overlap_path), mapping, "model_predicted_fva_overlap")

    rates = read_csv_if_exists(rate_path)
    mapped_rates = pd.DataFrame()
    if not rates.empty:
        rates = rates.rename(columns={"exchange_id": "reaction_id"})
        mapped_rates = attach_mapping(rates, mapping, "measured_exchange_rate")

    outputs = {
        "mapped_flux": mapped_flux,
        "mapped_fva": mapped_fva,
        "mapped_fva_overlap": mapped_overlap,
        "mapped_measured_rates": mapped_rates,
    }
    filenames = {
        "mapped_flux": "chompact_mapped_flux.csv",
        "mapped_fva": "chompact_mapped_fva.csv",
        "mapped_fva_overlap": "chompact_mapped_fva_overlap.csv",
        "mapped_measured_rates": "chompact_mapped_measured_rates.csv",
    }
    for key, df in outputs.items():
        if not df.empty:
            df.to_csv(os.path.join(chompact_dir, filenames[key]), index=False)

    qc = make_qc(outputs)
    qc.to_csv(os.path.join(chompact_dir, "chompact_mapping_qc.csv"), index=False)
    non_empty = [normalize_reaction_id(df) for df in outputs.values() if not df.empty]
    unmapped = pd.concat(non_empty, ignore_index=True) if non_empty else pd.DataFrame()
    if not unmapped.empty:
        cols = [c for c in ["reaction_id", "source_type", "reaction_name", "subsystem", "mapping_status"] if c in unmapped.columns]
        unmapped = unmapped.loc[unmapped["mapping_status"] == "unmapped", cols].drop_duplicates()
        unmapped.to_csv(os.path.join(chompact_dir, "chompact_unmapped_reactions.csv"), index=False)

    print(f"[saved] {os.path.relpath(chompact_dir, ROOT)}")
    print(qc.to_string(index=False))


if __name__ == "__main__":
    main()
