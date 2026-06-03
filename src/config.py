"""
config.py — CHO_METABOLOMICS global settings

2026-05 objective correction
---------------------------
The previous pipeline could select an arbitrary reaction such as BiGGRxn05 by
correlation screening. That is disabled for the main analysis. The production
objective is fixed to an IgG/mAb-producing reaction in the iCHO3K prod model.
"""
import os
import glob as _g


def _find_root():
    cur = os.path.dirname(os.path.abspath(__file__))
    for _ in range(6):
        if os.path.exists(os.path.join(cur, "src", "config.py")):
            return cur
        cur = os.path.dirname(cur)
    return cur


ROOT = _find_root()

# Paths
DATA_RAW = os.path.join(ROOT, "data", "raw")
DATA_PROCESSED = os.path.join(ROOT, "data", "processed")
MODEL_DIR = os.path.join(ROOT, "model", "iCHO3K")
LOGS_DIR = os.path.join(ROOT, "logs")

SOWA_MMC1 = os.path.join(DATA_RAW, "sowa2020", "mmc1.xlsx")
AA20_EXCEL = os.path.join(DATA_RAW, "practice_20aa", "CHO_raw_data_practice_20AA.xlsx")
AA20_TSV_DIR = os.path.join(DATA_RAW, "practice_20aa_tsv")
OWN_DIR = os.path.join(DATA_RAW, "own_experiment")

# TSV input filenames for offline/Linux workstations where Excel files may be
# inconvenient or DRM-wrapped. Put these files in either data/raw/practice_20aa_tsv
# or data/raw/own_experiment:
#   raw_timeseries.tsv, metabolite_map.tsv, feed_composition.tsv
RAW_TIMESERIES_TSV = "raw_timeseries.tsv"
METABOLITE_MAP_TSV = "metabolite_map.tsv"
FEED_COMPOSITION_TSV = "feed_composition.tsv"

# Model path: prefer prod model
_prod = _g.glob(os.path.join(MODEL_DIR, "*prod*.json"))
MODEL_PATH = _prod[0] if _prod else next(iter(_g.glob(os.path.join(MODEL_DIR, "*.json"))), None)


def results_dir(dataset: str, subdir: str = "tables") -> str:
    p = os.path.join(ROOT, "results", dataset, subdir)
    os.makedirs(p, exist_ok=True)
    return p


# FBA objective settings
# Main mAb production objective. DM_igg_g consumes the IgG produced by igg_formation.
PRODUCTION_OBJECTIVE = "DM_igg_g"
OBJ_RXN = PRODUCTION_OBJECTIVE  # backward-compatible name used by existing scripts

# Useful alternatives. The code resolves the first reaction that exists in the loaded model.
PRODUCTION_OBJECTIVE_CANDIDATES = [
    "DM_igg_g",
    "igg_formation",
    "igg_hc",
    "igg_lc",
]

# Growth/maintenance support. Default 0 means no minimum growth constraint is imposed.
# If you want growth-coupled production, try 0.05-0.20, but report it explicitly.
BIOMASS_RXN = "biomass_cho_prod"
BIOMASS_MIN_FRACTION = 0.0

FVA_FRACTION = 0.9
FVA_PROCESSES = int(os.environ.get("CHO_FVA_PROCESSES", "1"))  # Windows-safe default; set >1 on Linux/workstation.
FVA_SCOPE_DEFAULT = os.environ.get("CHO_FVA_SCOPE", "focused")  # exchange, focused, internal, all
FVA_FULL_TARGETS_DEFAULT = os.environ.get("CHO_FVA_TARGETS", "group_avg")  # group_avg, representative, all_clones

# Verified iCHO3K prod exchange IDs
EXCHANGE_IDS = {
    "Glucose": "EX_glc_e",
    "Lactate": "EX_lac_L_e",
    "Glutamine": "EX_gln_L_e",
    "Glutamate": "EX_glu_L_e",
    "Alanine": "EX_ala_L_e",
    "Arginine": "EX_arg_L_e",
    "Asparagine": "EX_asn_L_e",
    "Aspartate": "EX_asp_L_e",
    "Cysteine": "EX_cys_L_e",
    "Glycine": "EX_gly_e",
    "Histidine": "EX_his_L_e",
    "Isoleucine": "EX_ile_L_e",
    "Leucine": "EX_leu_L_e",
    "Lysine": "EX_lys_L_e",
    "Methionine": "EX_met_L_e",
    "Phenylalanine": "EX_phe_L_e",
    "Proline": "EX_pro_L_e",
    "Serine": "EX_ser_L_e",
    "Threonine": "EX_thr_L_e",
    "Tryptophan": "EX_trp_L_e",
    "Tyrosine": "EX_tyr_L_e",
    "Valine": "EX_val_L_e",
    "NH4+": "EX_nh4_e",
    "Pyruvate": "EX_pyr_e",
    "Citrate": "EX_cit_e",
    "Fumarate": "EX_fum_e",
    "Succinate": "EX_succ_e",
    "Malate": "EX_mal_L_e",
}

# Literature fallback constraints, used only when explicitly requested as fallback.
LITERATURE_HIGH = {
    "EX_glc_e": (-0.048, 0.0),
    "EX_lac_L_e": (0.0, 0.040),
    "EX_gln_L_e": (-0.012, 0.0),
    "EX_nh4_e": (0.0, 0.010),
    "EX_ala_L_e": (0.0, 0.006),
    "EX_glu_L_e": (-0.002, 0.002),
}
LITERATURE_LOW = {
    "EX_glc_e": (-0.055, 0.0),
    "EX_lac_L_e": (0.0, 0.090),
    "EX_gln_L_e": (-0.015, 0.0),
    "EX_nh4_e": (0.0, 0.020),
    "EX_ala_L_e": (0.0, 0.012),
    "EX_glu_L_e": (0.0, 0.006),
}

PALETTE = {
    "high": "#D62728",
    "low": "#1F77B4",
    "mid": ["#FF7F0E", "#2CA02C", "#9467BD", "#8C564B", "#E377C2", "#7F8C8D", "#BCBD22", "#17BECF"],
    "impr": "#2CA02C",
    "neutral": "#AAAAAA",
    "danger": "#C0392B",
}

for _d in [DATA_PROCESSED, LOGS_DIR]:
    os.makedirs(_d, exist_ok=True)

if __name__ == "__main__":
    print(f"ROOT       : {ROOT}")
    print(f"MODEL      : {MODEL_PATH}  {'OK' if MODEL_PATH and os.path.exists(MODEL_PATH) else 'MISSING'}")
    print(f"OBJECTIVE  : {OBJ_RXN}")
    print(f"BIOMASS    : {BIOMASS_RXN}  min_fraction={BIOMASS_MIN_FRACTION}")
    print(f"20AA DATA  : {AA20_EXCEL}  {'OK' if os.path.exists(AA20_EXCEL) else 'MISSING'}")

# IgG/mAb production zero-aware settings
# If DM_igg_g / igg_formation maximization is zero under strict constraints,
# this means the model cannot predict production under those constraints.
# Do NOT replace it with an unrelated objective. Instead, either report zero
# capacity or enforce measured IgG production as a demand lower bound.
IGG_MW_G_PER_MOL = 150000.0
PRODUCTION_ZERO_EPS = 1e-12
MEASURED_IGG_DEMAND_FRACTION = 0.80  # enforce 80% of observed interval qIgG for feasibility


# Constraint policy settings for production-objective analyses
CONSTRAINT_POLICY_DEFAULT = "production_relaxed"
CONSTRAINT_BUFFER = 0.20
CORE_HARD_METABOLITES = {"Glucose", "Lactate", "NH4+", "Ammonia"}
FLEXIBLE_NUTRIENT_METABOLITES = {
    "Glutamine", "Glutamate", "Alanine", "Arginine", "Asparagine", "Aspartate",
    "Cysteine", "Glycine", "Histidine", "Isoleucine", "Leucine", "Lysine",
    "Methionine", "Phenylalanine", "Proline", "Serine", "Threonine",
    "Tryptophan", "Tyrosine", "Valine",
}
MIN_NUTRIENT_UPTAKE_CAPACITY = 0.02


# Measured IgG-demand pFBA settings
# The synthetic practice file can have observed qIgG units larger than the model's
# DM_igg_g capacity. In that case, do NOT maximize the capped objective. Instead
# scale all observed qIgG demands by one global factor that keeps every clone
# feasible while preserving relative production differences.
MEASURED_IGG_DEMAND_SCALE = "auto"  # auto or numeric string/float, e.g. 0.20
MEASURED_IGG_DEMAND_SCALE_SAFETY = 0.90
FIXED_DEMAND_TOL = 1e-12


# ─────────────────────────────────────────────────────────────────────────────
# ML-ready prediction settings
# ─────────────────────────────────────────────────────────────────────────────
# FBA is used as a mechanistic feature generator. ML prediction is only trained
# when enough independent samples/runs/clones are available. With fewer samples,
# the pipeline exports feature tables only and reports that ML is not statistically
# meaningful yet.
ML_MIN_TRAIN_SAMPLES = 12
ML_TARGET_DEFAULT = "IgG_day14"  # alternatives: qIgG_interval, productivity_class
ML_FEATURE_EXCLUDE_PATTERNS = ["IgG", "qIgG", "target", "measured_igg", "scaled_igg", "production_flux"]
