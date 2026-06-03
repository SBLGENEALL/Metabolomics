import json
from pathlib import Path
import argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    args = ap.parse_args()

    base = Path(args.base)
    outdir = base / "maps"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "cho_exchange_flux_map.json"

    nodes, reactions, text_labels = {}, {}, {}
    nc = rc = tc = 0

    def text(x, y, t):
        nonlocal tc
        text_labels[str(tc)] = {"x": x, "y": y, "text": t}
        tc += 1

    def met(x, y, bigg, name):
        nonlocal nc
        nid = str(nc)
        nc += 1
        nodes[nid] = {
            "node_type": "metabolite",
            "x": x,
            "y": y,
            "bigg_id": bigg,
            "name": name,
            "label_x": x - 20,
            "label_y": y - 52,
            "node_is_primary": True,
        }
        return nid

    def marker(x, y, typ):
        nonlocal nc
        nid = str(nc)
        nc += 1
        nodes[nid] = {"node_type": typ, "x": x, "y": y}
        return nid

    def add_exchange(rxn, metid, name, x, y):
        nonlocal rc
        n1 = met(x, y, metid, name)
        n2 = marker(x + 150, y, "midmarker")
        n3 = marker(x + 300, y, "multimarker")
        key = str(rc)
        rc += 1
        reactions[key] = {
            "name": rxn,
            "bigg_id": rxn,
            "reversibility": True,
            "label_x": x + 55,
            "label_y": y + 58,
            "gene_reaction_rule": "",
            "genes": [],
            "metabolites": [{"bigg_id": metid, "coefficient": -1}],
            "segments": {
                f"{key}_0": {"from_node_id": n1, "to_node_id": n2, "b1": None, "b2": None},
                f"{key}_1": {"from_node_id": n2, "to_node_id": n3, "b1": None, "b2": None},
            },
        }

    text(-100, -170, "Carbon / Nitrogen")
    text(520, -170, "Amino acids I")
    text(1050, -170, "Amino acids II")
    text(1580, -170, "TCA-related / Objective")

    row = 165
    columns = [
        (0, [("EX_glc_e","glc_e","Glucose"),("EX_lac_L_e","lac_L_e","L-Lactate"),("EX_gln_L_e","gln_L_e","L-Glutamine"),("EX_glu_L_e","glu_L_e","L-Glutamate"),("EX_nh4_e","nh4_e","Ammonium"),("EX_ala_L_e","ala_L_e","L-Alanine")]),
        (560, [("EX_asn_L_e","asn_L_e","L-Asparagine"),("EX_asp_L_e","asp_L_e","L-Aspartate"),("EX_ser_L_e","ser_L_e","L-Serine"),("EX_gly_e","gly_e","Glycine"),("EX_pro_L_e","pro_L_e","L-Proline"),("EX_leu_L_e","leu_L_e","L-Leucine"),("EX_ile_L_e","ile_L_e","L-Isoleucine"),("EX_val_L_e","val_L_e","L-Valine")]),
        (1120, [("EX_lys_L_e","lys_L_e","L-Lysine"),("EX_arg_L_e","arg_L_e","L-Arginine"),("EX_his_L_e","his_L_e","L-Histidine"),("EX_thr_L_e","thr_L_e","L-Threonine"),("EX_phe_L_e","phe_L_e","L-Phenylalanine"),("EX_tyr_L_e","tyr_L_e","L-Tyrosine"),("EX_met_L_e","met_L_e","L-Methionine"),("EX_trp_L_e","trp_L_e","L-Tryptophan")]),
        (1680, [("EX_pyr_e","pyr_e","Pyruvate"),("EX_cit_e","cit_e","Citrate"),("EX_succ_e","succ_e","Succinate"),("EX_fum_e","fum_e","Fumarate"),("EX_mal_L_e","mal_L_e","L-Malate"),("biomass_cho_prod","biomass_c","Biomass/Productive demand")]),
    ]

    for x, items in columns:
        for i, (rxn, mid, name) in enumerate(items):
            add_exchange(rxn, mid, name, x, i * row)

    meta = {
        "map_name": "CHO exchange flux map",
        "map_id": "cho_exchange_flux_map",
        "map_description": "Exchange flux overview map",
        "homepage": "https://escher.github.io",
        "schema": "https://escher.github.io/escher/jsonschema/1-0-0#",
    }
    data = {
        "reactions": reactions,
        "nodes": nodes,
        "text_labels": text_labels,
        "canvas": {"x": -300, "y": -300, "width": 2400, "height": 1550},
    }

    with open(out, "w", encoding="utf-8") as f:
        json.dump([meta, data], f, indent=2)
    print("[OK] saved:", out)

if __name__ == "__main__":
    main()
