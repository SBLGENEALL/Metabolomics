"""
16_focused_escher_maps.py
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Create publication-friendly focused CHO/iCHO3K Escher maps and filtered flux/FVA
JSON files.

Why this exists
---------------
Large RECON/Escher maps are too broad for clone-comparison presentations. They
also mix many side branches and may show broken paths when reaction IDs do not
match iCHO3K exactly. This step builds small maps using iCHO3K reaction IDs, then
filters flux/FVA JSON data to only the reactions drawn on those maps.

Outputs
-------
results/<dataset>/escher_maps/focused/
  - CHO_focus_core_carbon_map.json
  - CHO_focus_glycolysis_lactate_map.json
  - CHO_focus_ppp_map.json
  - CHO_focus_tca_gln_map.json
  - focused_escher_guide.md

results/<dataset>/tables/focused_escher/
  - filtered flux JSON files for the focused maps
  - FVA range/mean JSON files when central_mab_fva_report.csv exists

results/<dataset>/figures/
  - Fig12_focused_core_fva_range.png when FVA data exists
"""
import argparse
import glob
import json
import os
import sys
from typing import Dict, Tuple, List

import numpy as np
import pandas as pd


def _find_root():
    cur = os.path.dirname(os.path.abspath(__file__))
    for _ in range(7):
        if os.path.exists(os.path.join(cur, "src", "config.py")):
            return cur
        cur = os.path.dirname(cur)
    return cur


ROOT = _find_root()
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from src.config import MODEL_PATH, results_dir  # noqa

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa")
parser.add_argument("--model", default=None, help="Optional COBRA JSON model path")
parser.add_argument("--analysis_mode", default=os.environ.get("CHO_ANALYSIS_MODE", "both"), choices=["predict", "explain", "both"])
args = parser.parse_args()

DATASET = args.dataset
MODEL_JSON = args.model or MODEL_PATH
if not MODEL_JSON or not os.path.exists(MODEL_JSON):
    raise FileNotFoundError(f"Model JSON not found: {MODEL_JSON}")

with open(MODEL_JSON, "r", encoding="utf-8") as f:
    model = json.load(f)
rxns: Dict[str, dict] = {r["id"]: r for r in model.get("reactions", [])}
mets: Dict[str, dict] = {m["id"]: m for m in model.get("metabolites", [])}

OUT_MAP_DIR = os.path.join(results_dir(DATASET, "escher_maps"), "focused")
OUT_TABLE_DIR = os.path.join(results_dir(DATASET, "tables"), "focused_escher")
OUT_FIG_DIR = results_dir(DATASET, "figures")
os.makedirs(OUT_MAP_DIR, exist_ok=True)
os.makedirs(OUT_TABLE_DIR, exist_ok=True)
os.makedirs(OUT_FIG_DIR, exist_ok=True)

# ---------------------------------------------------------------------
# Focused reaction templates. Reaction IDs are iCHO3K IDs.
# ---------------------------------------------------------------------
GLYCOLYSIS_LACTATE = [
    ("EX_glc_e", "glc_D_boundary", "glc_D_e"),
    ("GLCt1r", "glc_D_e", "glc_D_c"),
    ("HEX1", "glc_D_c", "g6p_c"),
    ("PGI", "g6p_c", "f6p_c"),
    ("PFK", "f6p_c", "fdp_c"),
    ("FBA", "fdp_c", "g3p_c"),
    ("TPI", "dhap_c", "g3p_c"),
    ("GAPD", "g3p_c", "13dpg_c"),
    ("PGK", "13dpg_c", "3pg_c"),
    ("PGM", "3pg_c", "2pg_c"),
    ("ENO", "2pg_c", "pep_c"),
    ("PYK", "pep_c", "pyr_c"),
    ("LDH_L", "pyr_c", "lac_L_c"),
    ("L_LACt2r", "lac_L_c", "lac_L_e"),
    ("EX_lac_L_e", "lac_L_e", "lac_L_boundary"),
]

PPP = [
    ("G6PDH2r", "g6p_c", "6pgl_c"),
    ("PGL", "6pgl_c", "6pgc_c"),
    ("GND", "6pgc_c", "ru5p_D_c"),
    ("RPI", "ru5p_D_c", "r5p_c"),
    ("RPE", "ru5p_D_c", "xu5p_D_c"),
    ("TKT1", "r5p_c", "sdhpt7p_c"),
    ("TALA", "sdhpt7p_c", "f6p_c"),
    ("TKT2", "xu5p_D_c", "g3p_c"),
]

TCA_GLN = [
    ("PYRt2m", "pyr_c", "pyr_m"),
    ("PCm", "pyr_m", "oaa_m"),
    ("CSm", "oaa_m", "cit_m"),
    ("ACONTam", "cit_m", "acon_C_m"),
    ("ACONTbm", "acon_C_m", "icit_m"),
    ("ICDHyrm", "icit_m", "akg_m"),
    ("AKGDm", "akg_m", "succoa_m"),
    ("SUCOASm", "succoa_m", "succ_m"),
    ("SUCD1m", "succ_m", "fum_m"),
    ("FUMm", "fum_m", "mal_L_m"),
    ("MDHm", "mal_L_m", "oaa_m"),
    ("EX_gln_L_e", "gln_L_boundary", "gln_L_e"),
    ("GLNt4rev", "gln_L_e", "gln_L_c"),
    ("GLNtm", "gln_L_c", "gln_L_m"),
    ("GLUNm", "gln_L_m", "glu_L_m"),
    ("GLUDym", "glu_L_m", "akg_m"),
    ("GLUDxm", "glu_L_m", "akg_m"),
    ("NH4t3r", "nh4_c", "nh4_e"),
    ("EX_nh4_e", "nh4_e", "nh4_boundary"),
]

# All-core map is intentionally smaller than official RECON maps: it keeps only
# the path a presenter can verbally follow.
CORE = GLYCOLYSIS_LACTATE + PPP + TCA_GLN

# Coordinate set shared by all maps. Smaller maps reuse subsets and auto-crop.
COORDS: Dict[str, Tuple[float, float]] = {
    # glycolysis vertical lane
    "glc_D_boundary": (0, 0),
    "glc_D_e": (0, 90),
    "glc_D_c": (0, 190),
    "g6p_c": (0, 310),
    "f6p_c": (0, 430),
    "fdp_c": (0, 550),
    "dhap_c": (-190, 655),
    "g3p_c": (0, 670),
    "13dpg_c": (0, 790),
    "3pg_c": (0, 910),
    "2pg_c": (0, 1030),
    "pep_c": (0, 1150),
    "pyr_c": (0, 1280),
    "lac_L_c": (-230, 1380),
    "lac_L_e": (-230, 1490),
    "lac_L_boundary": (-230, 1600),
    # PPP right branch
    "6pgl_c": (330, 310),
    "6pgc_c": (560, 310),
    "ru5p_D_c": (790, 310),
    "r5p_c": (1000, 210),
    "xu5p_D_c": (1000, 430),
    "sdhpt7p_c": (610, 560),
    "e4p_c": (850, 565),
    # pyruvate/TCA and glutamine
    "pyr_m": (330, 1280),
    "oaa_m": (700, 1180),
    "cit_m": (970, 1040),
    "acon_C_m": (1230, 1120),
    "icit_m": (1320, 1350),
    "akg_m": (1110, 1600),
    "succoa_m": (790, 1700),
    "succ_m": (540, 1580),
    "fum_m": (450, 1340),
    "mal_L_m": (540, 1120),
    "gln_L_boundary": (-590, 1030),
    "gln_L_e": (-590, 1140),
    "gln_L_c": (-590, 1260),
    "gln_L_m": (-240, 1560),
    "glu_L_m": (270, 1580),
    "nh4_c": (-590, 1440),
    "nh4_e": (-590, 1560),
    "nh4_boundary": (-590, 1680),
}

BOUNDARY_NAMES = {
    "glc_D_boundary": "Glucose boundary",
    "lac_L_boundary": "Lactate boundary",
    "gln_L_boundary": "Glutamine boundary",
    "nh4_boundary": "NH4 boundary",
}

MET_LABELS = {
    "glc_D_e": "Glucose[e]", "glc_D_c": "Glucose[c]", "g6p_c": "G6P", "f6p_c": "F6P",
    "fdp_c": "FBP", "dhap_c": "DHAP", "g3p_c": "G3P", "13dpg_c": "1,3-BPG",
    "3pg_c": "3PG", "2pg_c": "2PG", "pep_c": "PEP", "pyr_c": "Pyruvate[c]",
    "lac_L_c": "Lactate[c]", "lac_L_e": "Lactate[e]", "6pgl_c": "6PGL", "6pgc_c": "6PG",
    "ru5p_D_c": "Ru5P", "r5p_c": "R5P", "xu5p_D_c": "Xu5P", "sdhpt7p_c": "S7P",
    "e4p_c": "E4P", "pyr_m": "Pyruvate[m]", "oaa_m": "OAA[m]", "cit_m": "Citrate[m]",
    "acon_C_m": "Aconitate[m]", "icit_m": "Isocitrate[m]", "akg_m": "αKG[m]", "succoa_m": "Succinyl-CoA[m]",
    "succ_m": "Succinate[m]", "fum_m": "Fumarate[m]", "mal_L_m": "Malate[m]", "gln_L_e": "Gln[e]",
    "gln_L_c": "Gln[c]", "gln_L_m": "Gln[m]", "glu_L_m": "Glu[m]", "nh4_c": "NH4[c]", "nh4_e": "NH4[e]",
}


def met_name(met_id: str) -> str:
    return MET_LABELS.get(met_id) or BOUNDARY_NAMES.get(met_id) or mets.get(met_id, {}).get("name") or met_id


class MiniEscherMap:
    def __init__(self, name: str, reactions: List[Tuple[str, str, str]], labels: Dict[str, dict]):
        self.name = name
        self.template = reactions
        self.labels = labels
        self.next_id = 1
        self.nodes = {}
        self.reactions = {}
        self.met_node_id = {}
        self.drawn_rxns = []
        self.skipped = []

    def _id(self):
        out = str(self.next_id)
        self.next_id += 1
        return out

    def add_met(self, met_id):
        if met_id in self.met_node_id:
            return self.met_node_id[met_id]
        x, y = COORDS[met_id]
        nid = self._id()
        self.nodes[nid] = {
            "node_type": "metabolite",
            "x": x,
            "y": y,
            "bigg_id": met_id,
            "name": met_name(met_id),
            "label_x": x + 18,
            "label_y": y - 28,
            "node_is_primary": True,
        }
        self.met_node_id[met_id] = nid
        return nid

    def add_marker(self, node_type, x, y):
        nid = self._id()
        self.nodes[nid] = {"node_type": node_type, "x": x, "y": y}
        return nid

    def add_rxn(self, rid, start, end):
        if rid not in rxns:
            self.skipped.append(f"{rid} missing from model")
            return
        if start not in COORDS or end not in COORDS:
            self.skipped.append(f"{rid} missing coordinates")
            return
        x1, y1 = COORDS[start]
        x2, y2 = COORDS[end]
        start_id = self.add_met(start)
        end_id = self.add_met(end)
        m1 = self.add_marker("multimarker", x1 + 0.33 * (x2 - x1), y1 + 0.33 * (y2 - y1))
        mid = self.add_marker("midmarker", x1 + 0.50 * (x2 - x1), y1 + 0.50 * (y2 - y1))
        m2 = self.add_marker("multimarker", x1 + 0.67 * (x2 - x1), y1 + 0.67 * (y2 - y1))
        segments = {}
        for a, b in [(start_id, m1), (m1, mid), (mid, m2), (m2, end_id)]:
            sid = self._id()
            segments[sid] = {"from_node_id": a, "to_node_id": b, "b1": None, "b2": None}
        r = rxns[rid]
        lower = float(r.get("lower_bound", 0.0))
        upper = float(r.get("upper_bound", 0.0))
        escher_rid = self._id()
        self.reactions[escher_rid] = {
            "name": r.get("name") or rid,
            "bigg_id": rid,
            "reversibility": bool(lower < 0 and upper > 0),
            "label_x": (x1 + x2) / 2 + 12,
            "label_y": (y1 + y2) / 2 - 14,
            "gene_reaction_rule": r.get("gene_reaction_rule", ""),
            "metabolites": [
                {"bigg_id": m, "coefficient": float(c)}
                for m, c in r.get("metabolites", {}).items()
            ],
            "segments": segments,
        }
        self.drawn_rxns.append(rid)

    def build(self):
        for rid, start, end in self.template:
            self.add_rxn(rid, start, end)
        xs = [n.get("x") for n in self.nodes.values() if "x" in n]
        ys = [n.get("y") for n in self.nodes.values() if "y" in n]
        pad = 180
        canvas = {
            "x": min(xs) - pad if xs else -100,
            "y": min(ys) - pad if ys else -100,
            "width": (max(xs) - min(xs) + 2 * pad) if xs else 1000,
            "height": (max(ys) - min(ys) + 2 * pad) if ys else 1000,
        }
        info = {
            "map_name": self.name,
            "map_id": self.name,
            "map_description": "Focused CHO/iCHO3K map for clone flux/FVA comparison. Uses iCHO3K reaction IDs.",
            "homepage": "https://escher.github.io/",
            "schema": "https://escher.github.io/escher/jsonschema/1-0-0#",
        }
        body = {"reactions": self.reactions, "nodes": self.nodes, "text_labels": self.labels, "canvas": canvas}
        return [info, body]


def write_map(name, reactions, labels):
    m = MiniEscherMap(name, reactions, labels)
    out = m.build()
    path = os.path.join(OUT_MAP_DIR, f"{name}.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    rxn_path = os.path.join(OUT_MAP_DIR, f"{name}_reaction_ids.txt")
    with open(rxn_path, "w", encoding="utf-8") as f:
        for rid in m.drawn_rxns:
            f.write(rid + "\n")
        if m.skipped:
            f.write("\n# skipped\n")
            for s in m.skipped:
                f.write("# " + s + "\n")
    print(f"  [saved] {os.path.relpath(path, ROOT)} ({len(m.drawn_rxns)} reactions)")
    return path, set(m.drawn_rxns), m.skipped


maps = {}
map_specs = [
    ("CHO_focus_glycolysis_lactate_map", GLYCOLYSIS_LACTATE, {"1": {"x": -80, "y": -70, "text": "Glucose → Glycolysis → Pyruvate/Lactate"}}),
    ("CHO_focus_ppp_map", [("HEX1", "glc_D_c", "g6p_c"), ("PGI", "g6p_c", "f6p_c")] + PPP, {"1": {"x": 280, "y": 160, "text": "Pentose Phosphate Pathway + return to F6P/G3P"}}),
    ("CHO_focus_tca_gln_map", TCA_GLN, {"1": {"x": 620, "y": 930, "text": "TCA + Glutamine/Glutamate Anaplerosis"}}),
    ("CHO_focus_core_carbon_map", CORE, {
        "1": {"x": -100, "y": -70, "text": "Glucose/Glycolysis"},
        "2": {"x": 420, "y": 180, "text": "PPP / NADPH support"},
        "3": {"x": 650, "y": 930, "text": "TCA cycle"},
        "4": {"x": -680, "y": 950, "text": "Gln/Glu/NH4"},
        "5": {"x": -330, "y": 1330, "text": "Lactate"},
    }),
]

all_focus_rxns = set()
for name, reactions, labels in map_specs:
    path, rxn_set, skipped = write_map(name, reactions, labels)
    maps[name] = rxn_set
    all_focus_rxns.update(rxn_set)

# ---------------------------------------------------------------------
# Filter flux JSON files for focused maps.
# ---------------------------------------------------------------------
def filter_json_files():
    tables_dir = results_dir(DATASET, "tables")
    candidates = sorted(glob.glob(os.path.join(tables_dir, "escher_flux_*_cho*.json")))
    # keep CHO-ID files only; skip existing focused outputs
    n_written = 0
    for p in candidates:
        if "focused_escher" in p:
            continue
        base = os.path.basename(p)
        try:
            d = json.load(open(p, "r", encoding="utf-8"))
        except Exception:
            continue
        filt = {rid: float(v) for rid, v in d.items() if rid in all_focus_rxns and abs(float(v)) > 1e-12}
        out = os.path.join(OUT_TABLE_DIR, "focused_" + base)
        with open(out, "w", encoding="utf-8") as f:
            json.dump(filt, f, indent=2, ensure_ascii=False, sort_keys=True)
        n_written += 1
    print(f"  [saved] focused flux JSON files: {n_written}")


filter_json_files()

# ---------------------------------------------------------------------
# Export FVA range/mean JSON files and Fig12 when FVA exists.
# ---------------------------------------------------------------------
def export_fva():
    fva_path = os.path.join(results_dir(DATASET, "tables"), "central_mab_fva_report.csv")
    if not os.path.exists(fva_path):
        print("  [skip] central_mab_fva_report.csv not found — run step 14 to get FVA outputs")
        return
    fva = pd.read_csv(fva_path)
    if fva.empty or "rxn_id" not in fva.columns:
        print("  [skip] FVA table empty")
        return
    fva = fva[fva["rxn_id"].isin(all_focus_rxns)].copy()
    if fva.empty:
        print("  [skip] FVA table has no focused-map reactions")
        return
    # Normalize expected columns.
    fva["range"] = fva["maximum"].astype(float) - fva["minimum"].astype(float)
    fva["mean"] = (fva["maximum"].astype(float) + fva["minimum"].astype(float)) / 2.0

    written = 0
    for mode in sorted(fva["mode"].dropna().unique()):
        subm = fva[fva["mode"] == mode]
        for cond in sorted(subm["condition"].dropna().unique()):
            sub = subm[subm["condition"] == cond]
            for value_col in ["range", "mean", "minimum", "maximum"]:
                d = {r: round(float(v), 8) for r, v in zip(sub["rxn_id"], sub[value_col]) if np.isfinite(float(v))}
                out = os.path.join(OUT_TABLE_DIR, f"focused_fva_{mode}_{cond}_{value_col}_cho.json")
                with open(out, "w", encoding="utf-8") as f:
                    json.dump(d, f, indent=2, ensure_ascii=False, sort_keys=True)
                written += 1
        # High-Low FVA range delta.
        conds = set(subm["condition"])
        if "HighAvg" in conds and "LowAvg" in conds:
            h = subm[subm["condition"] == "HighAvg"].set_index("rxn_id")
            l = subm[subm["condition"] == "LowAvg"].set_index("rxn_id")
            common = sorted(set(h.index) & set(l.index))
            for value_col in ["range", "mean"]:
                d = {rid: round(float(h.loc[rid, value_col] - l.loc[rid, value_col]), 8) for rid in common}
                out = os.path.join(OUT_TABLE_DIR, f"focused_fva_{mode}_high_minus_low_{value_col}_cho.json")
                with open(out, "w", encoding="utf-8") as f:
                    json.dump(d, f, indent=2, ensure_ascii=False, sort_keys=True)
                written += 1
    print(f"  [saved] focused FVA JSON files: {written}")

    # Figure 12: FVA range heatmap on focused reactions.
    try:
        import matplotlib.pyplot as plt
        plot_mode = "measured_demand" if "measured_demand" in set(fva["mode"]) else sorted(fva["mode"].unique())[0]
        sub = fva[fva["mode"] == plot_mode].copy()
        order_rxns = [rid for rid, _, _ in CORE if rid in set(sub["rxn_id"])]
        # Use readable labels if present.
        if "label" in sub.columns:
            sub["row_label"] = sub["rxn_id"] + " | " + sub["label"].astype(str)
        else:
            sub["row_label"] = sub["rxn_id"]
        label_map = sub.drop_duplicates("rxn_id").set_index("rxn_id")["row_label"].to_dict()
        keep = [r for r in order_rxns if r in label_map]
        piv = sub.pivot_table(index="rxn_id", columns="condition", values="range", aggfunc="first").fillna(0.0)
        piv = piv.loc[[r for r in keep if r in piv.index]]
        # Prefer group-average columns first.
        pref = [c for c in ["HighAvg", "MotherAvg", "LowAvg"] if c in piv.columns]
        rest = [c for c in piv.columns if c not in pref]
        piv = piv[pref + rest]
        fig_h = max(7, 0.26 * len(piv) + 2)
        fig, ax = plt.subplots(figsize=(10.5, fig_h))
        im = ax.imshow(piv.values, aspect="auto")
        ax.set_yticks(np.arange(len(piv.index)))
        ax.set_yticklabels([label_map.get(r, r) for r in piv.index], fontsize=8)
        ax.set_xticks(np.arange(len(piv.columns)))
        ax.set_xticklabels(piv.columns, rotation=35, ha="right")
        ax.set_title(f"Figure 12. Focused core pathway FVA range ({plot_mode})", weight="bold")
        ax.set_xlabel("Condition")
        ax.set_ylabel("Reaction")
        cbar = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.02)
        cbar.set_label("FVA range = max - min")
        fig.tight_layout()
        outfig = os.path.join(OUT_FIG_DIR, "Fig12_focused_core_fva_range.png")
        fig.savefig(outfig, dpi=220, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {os.path.relpath(outfig, ROOT)}")

        # Figure 12B: row-wise z-score of FVA range. Raw FVA range is dominated
        # by reactions that are inherently flexible in all clones; z-scoring each
        # reaction emphasizes which clones/groups are relatively more constrained
        # or flexible for that reaction.
        try:
            row_mean = piv.mean(axis=1)
            row_std = piv.std(axis=1).replace(0, np.nan)
            zpiv = piv.sub(row_mean, axis=0).div(row_std, axis=0).fillna(0.0)
            fig, ax = plt.subplots(figsize=(10.5, fig_h))
            vmax = max(1.0, float(np.nanmax(np.abs(zpiv.values)))) if zpiv.size else 1.0
            im = ax.imshow(zpiv.values, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
            ax.set_yticks(np.arange(len(zpiv.index)))
            ax.set_yticklabels([label_map.get(r, r) for r in zpiv.index], fontsize=8)
            ax.set_xticks(np.arange(len(zpiv.columns)))
            ax.set_xticklabels(zpiv.columns, rotation=35, ha="right")
            ax.set_title(f"Figure 12B. Relative FVA range pattern across clones/groups ({plot_mode})", weight="bold")
            ax.set_xlabel("Condition")
            ax.set_ylabel("Reaction")
            cbar = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.02)
            cbar.set_label("Row-wise z-score of FVA range")
            fig.tight_layout()
            outfig2 = os.path.join(OUT_FIG_DIR, "Fig12B_focused_core_fva_range_zscore.png")
            fig.savefig(outfig2, dpi=220, bbox_inches="tight")
            plt.close(fig)
            print(f"  [saved] {os.path.relpath(outfig2, ROOT)}")
        except Exception as e:
            print(f"  [skip] Fig12B failed: {e}")

        # Figure 12C: HighAvg - LowAvg FVA range delta. This summarizes where
        # the high-producer state has more or less feasible-space flexibility.
        try:
            if "HighAvg" in piv.columns and "LowAvg" in piv.columns:
                delta = (piv["HighAvg"] - piv["LowAvg"]).sort_values(key=lambda x: x.abs(), ascending=False).head(25).iloc[::-1]
                fig, ax = plt.subplots(figsize=(9.5, max(5.5, 0.32 * len(delta) + 1.2)))
                ax.barh([label_map.get(r, r) for r in delta.index], delta.values)
                ax.axvline(0, color="black", lw=1)
                ax.set_xlabel("HighAvg - LowAvg FVA range")
                ax.set_title(f"Figure 12C. High-vs-Low FVA flexibility difference ({plot_mode})", weight="bold")
                ax.grid(True, axis="x", ls=":", alpha=0.5)
                fig.tight_layout()
                outfig3 = os.path.join(OUT_FIG_DIR, "Fig12C_focused_fva_high_low_delta.png")
                fig.savefig(outfig3, dpi=220, bbox_inches="tight")
                plt.close(fig)
                print(f"  [saved] {os.path.relpath(outfig3, ROOT)}")
        except Exception as e:
            print(f"  [skip] Fig12C failed: {e}")

    except Exception as e:
        print(f"  [skip] Fig12 failed: {e}")


export_fva()

# ---------------------------------------------------------------------
# Guide
# ---------------------------------------------------------------------
guide = f"""# Focused Escher map guide

These maps are the recommended alternative to large RECON/global maps when the
question is: **High vs Low CHO producer flux phenotype**.

## Load in Escher

1. Open <https://escher.github.io>
2. `Load map` → select one map from `results/{DATASET}/escher_maps/focused/`
3. `Load model` is optional. If loaded, use the iCHO3K JSON model.
4. `Load reaction data` → select one focused JSON from `results/{DATASET}/tables/focused_escher/`

## Most useful combinations

### High vs Low flux difference

Map:
`CHO_focus_core_carbon_map.json`

Data:
`focused_escher_flux_high_minus_low_cho.json`

Interpretation: positive = higher in HighAvg; negative = higher in LowAvg.

### Individual pathway maps

- Glycolysis/lactate: `CHO_focus_glycolysis_lactate_map.json`
- PPP: `CHO_focus_ppp_map.json`
- TCA + glutamine: `CHO_focus_tca_gln_map.json`

Use the same data file, for example `focused_escher_flux_high_minus_low_cho.json`.
The map only displays reactions drawn on that focused pathway.

## FVA outputs

If step 14 was run, this folder also contains files such as:

- `focused_fva_measured_demand_HighAvg_range_cho.json`
- `focused_fva_measured_demand_LowAvg_range_cho.json`
- `focused_fva_measured_demand_high_minus_low_range_cho.json`

These are not pFBA fluxes. They show feasible flux **range width** or range
center for the focused reactions. Escher can color them, but it cannot display
min-max bars directly; use `Fig12_focused_core_fva_range.png` and
`central_mab_fva_report.csv` for full FVA interpretation.

## Why this map may still show some apparent gaps

The iCHO3K prod model used here does not contain a simple canonical PDH reaction
ID in the form `PDHm`; TCA entry is represented through the reactions present in
this model, such as `PCm`, `CSm`, and anaplerotic glutamine/glutamate reactions.
Therefore, do not interpret absence of a drawn PDH arrow as absence of biology;
it is a model-reaction representation issue.
"""
with open(os.path.join(OUT_MAP_DIR, "focused_escher_guide.md"), "w", encoding="utf-8") as f:
    f.write(guide)
print(f"  [saved] {os.path.relpath(os.path.join(OUT_MAP_DIR, 'focused_escher_guide.md'), ROOT)}")
print("  OK 16_focused_escher_maps 완료")
