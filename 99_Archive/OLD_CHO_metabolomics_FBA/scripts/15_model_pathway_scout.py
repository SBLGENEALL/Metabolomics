
"""
15_model_pathway_scout.py

Purpose
-------
Scan the iCHO3K / COBRA model and make candidate internal reaction sets
for pathway-level pFBA/FVA interpretation.

This is the missing step between:
  "I only plotted measured extracellular metabolite changes"
and:
  "I used exchange-constrained FBA/FVA to compare model-inferred pathway phenotypes."

Outputs
-------
results/tables/model_reaction_scout.csv
results/tables/objective_candidates.csv
results/tables/pathway_reaction_sets.csv

Run
---
cd C:\CHO_POC_ChatGPT_260519

python scripts\15_model_pathway_scout.py ^
  --base C:\CHO_POC_ChatGPT_260519 ^
  --model-file "C:\CHO_POC_ChatGPT_260519\model\iCHO3K-main\iCHO3K\Model\iCHO3K_cho_prod_generic_unblocked.json"

If your model is already in the default path, --model-file can be omitted.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import cobra
import pandas as pd


PATHWAY_KEYWORDS = {
    "glycolysis": [
        "hexokinase", "glucokinase", "phosphofructokinase", "fructose-bisphosphate",
        "glyceraldehyde", "phosphoglycerate", "enolase", "pyruvate kinase",
        "HK", "PFK", "FBA", "TPI", "GAPD", "PGK", "PGM", "ENO", "PYK",
        "glucose-6-phosphate", "fructose-6-phosphate", "pyruvate"
    ],
    "tca_cycle": [
        "citrate synthase", "aconitase", "isocitrate", "2-oxoglutarate",
        "alpha-ketoglutarate", "ketoglutarate", "succinate", "fumarate",
        "malate", "oxaloacetate", "CS", "ACONT", "ICDH", "AKGD", "SUCOAS",
        "FUM", "MDH", "citrate", "TCA"
    ],
    "lactate_pyruvate": [
        "lactate", "lactate dehydrogenase", "LDH", "pyruvate", "PYR",
        "EX_lac", "lac_L", "l-Lactate"
    ],
    "glutamine_glutamate": [
        "glutamine", "glutamate", "glutaminase", "glutamine synthetase",
        "glutamate dehydrogenase", "GLN", "GLU", "GLNS", "GLUD", "GLUN",
        "alpha-ketoglutarate", "2-oxoglutarate"
    ],
    "ammonia_nitrogen": [
        "ammonia", "ammonium", "NH4", "nitrogen", "glutaminase",
        "glutamate dehydrogenase", "urea", "EX_nh4"
    ],
    "oxidative_phosphorylation": [
        "NADH dehydrogenase", "cytochrome", "ubiquinone", "ATP synthase",
        "oxidative phosphorylation", "electron transport", "complex I",
        "complex II", "complex III", "complex IV", "ATPS", "CYT", "NADH"
    ],
    "pentose_phosphate": [
        "glucose-6-phosphate dehydrogenase", "6-phosphogluconate",
        "ribulose", "ribose", "transketolase", "transaldolase",
        "G6PD", "PGL", "GND", "TKT", "TALA", "pentose phosphate"
    ],
    "amino_acid_transport": [
        "exchange", "transport", "EX_", "alanine", "arginine", "asparagine",
        "aspartate", "histidine", "isoleucine", "leucine", "lysine",
        "methionine", "phenylalanine", "serine", "threonine", "tryptophan",
        "tyrosine", "valine"
    ],
    "biomass_product": [
        "biomass", "growth", "product", "igg", "antibody", "mab",
        "protein", "secreted", "demand"
    ],
}


OBJECTIVE_KEYWORDS = [
    "biomass", "growth", "product", "igg", "antibody", "mab",
    "protein", "demand", "DM_", "sink"
]


def pick_model(base: Path, model_file: str | None) -> Path:
    if model_file:
        p = Path(model_file)
        return p if p.is_absolute() else base / p

    candidates = [
        base / "model" / "iCHO3K-main" / "iCHO3K" / "Model" / "iCHO3K_cho_prod_generic_unblocked.json",
        base / "model" / "iCHO3K-main" / "iCHO3K" / "Model" / "iCHO3K_cho_prod_generic_unblocked.xml",
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError("Model file not found. Use --model-file")


def load_model(model_file: Path):
    if model_file.suffix.lower() == ".json":
        return cobra.io.load_json_model(str(model_file))
    return cobra.io.read_sbml_model(str(model_file))


def text_blob(rxn) -> str:
    mets = []
    met_names = []
    for m in rxn.metabolites:
        mets.append(m.id)
        met_names.append(m.name or "")
    parts = [
        rxn.id,
        rxn.name or "",
        rxn.reaction or "",
        " ".join(mets),
        " ".join(met_names),
        str(rxn.gene_reaction_rule or ""),
    ]
    return " | ".join(parts)


def score_keywords(blob: str, keywords: list[str]) -> tuple[int, list[str]]:
    b = blob.lower()
    hits = []
    for kw in keywords:
        if kw.lower() in b:
            hits.append(kw)
    return len(hits), hits


def is_exchange_like(rxn) -> bool:
    return rxn.boundary or rxn.id.startswith("EX_") or rxn.id.startswith("DM_") or rxn.id.startswith("SK_")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--model-file", default=None)
    ap.add_argument("--min-score", type=int, default=1)
    ap.add_argument("--include-boundary", action="store_true", help="include exchange/demand/sink reactions in pathway sets")
    args = ap.parse_args()

    base = Path(args.base)
    model_file = pick_model(base, args.model_file)
    model = load_model(model_file)

    outdir = base / "results" / "tables"
    outdir.mkdir(parents=True, exist_ok=True)

    print("[INFO] model:", model_file)
    print("[INFO] reactions:", len(model.reactions))

    scout_rows = []
    pathway_rows = []
    objective_rows = []

    for rxn in model.reactions:
        blob = text_blob(rxn)
        boundary = is_exchange_like(rxn)

        matched_pathways = []
        matched_keywords = []
        for pathway, kws in PATHWAY_KEYWORDS.items():
            score, hits = score_keywords(blob, kws)
            if score >= args.min_score:
                matched_pathways.append(pathway)
                matched_keywords += hits

                if args.include_boundary or not boundary or pathway in ["amino_acid_transport", "biomass_product"]:
                    pathway_rows.append({
                        "pathway": pathway,
                        "reaction": rxn.id,
                        "reaction_name": rxn.name,
                        "lower_bound": rxn.lower_bound,
                        "upper_bound": rxn.upper_bound,
                        "is_boundary": boundary,
                        "include": True,
                        "keyword_score": score,
                        "matched_keywords": "; ".join(hits),
                        "reaction_string": rxn.reaction,
                    })

        obj_score, obj_hits = score_keywords(blob, OBJECTIVE_KEYWORDS)
        if obj_score > 0:
            old_obj = model.objective
            try:
                model.objective = rxn.id
                sol = model.optimize()
                status = sol.status
                value = sol.objective_value
            except Exception as e:
                status = f"error: {e}"
                value = None
            finally:
                try:
                    model.objective = old_obj
                except Exception:
                    pass

            objective_rows.append({
                "reaction": rxn.id,
                "reaction_name": rxn.name,
                "objective_test_status": status,
                "objective_test_value": value,
                "is_boundary": boundary,
                "keyword_score": obj_score,
                "matched_keywords": "; ".join(obj_hits),
                "lower_bound": rxn.lower_bound,
                "upper_bound": rxn.upper_bound,
                "reaction_string": rxn.reaction,
            })

        scout_rows.append({
            "reaction": rxn.id,
            "reaction_name": rxn.name,
            "subsystem": getattr(rxn, "subsystem", ""),
            "lower_bound": rxn.lower_bound,
            "upper_bound": rxn.upper_bound,
            "is_boundary": boundary,
            "matched_pathways": "; ".join(sorted(set(matched_pathways))),
            "matched_keywords": "; ".join(sorted(set(matched_keywords))),
            "reaction_string": rxn.reaction,
            "gene_rule": rxn.gene_reaction_rule,
        })

    scout = pd.DataFrame(scout_rows)
    pathway = pd.DataFrame(pathway_rows).drop_duplicates(subset=["pathway", "reaction"])
    objective = pd.DataFrame(objective_rows)

    if not objective.empty:
        objective["objective_test_value_numeric"] = pd.to_numeric(objective["objective_test_value"], errors="coerce")
        objective = objective.sort_values(
            ["keyword_score", "objective_test_value_numeric"],
            ascending=[False, False],
            na_position="last",
        ).drop(columns=["objective_test_value_numeric"])

    scout.to_csv(outdir / "model_reaction_scout.csv", index=False, encoding="utf-8-sig")
    objective.to_csv(outdir / "objective_candidates.csv", index=False, encoding="utf-8-sig")
    pathway.to_csv(outdir / "pathway_reaction_sets.csv", index=False, encoding="utf-8-sig")

    print("[OK] saved:", outdir / "model_reaction_scout.csv")
    print("[OK] saved:", outdir / "objective_candidates.csv")
    print("[OK] saved:", outdir / "pathway_reaction_sets.csv")

    print("\n[OBJECTIVE CANDIDATES TOP 20]")
    if objective.empty:
        print("No objective candidates found.")
    else:
        cols = ["reaction", "reaction_name", "objective_test_status", "objective_test_value", "matched_keywords"]
        print(objective[cols].head(20).to_string(index=False))

    print("\n[PATHWAY REACTION COUNTS]")
    if pathway.empty:
        print("No pathway reactions found.")
    else:
        print(pathway.groupby("pathway")["reaction"].nunique().sort_values(ascending=False).to_string())


if __name__ == "__main__":
    main()
