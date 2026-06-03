"""
config.py — CHO_METABOLOMICS 전역 설정
모든 스크립트에서 from src.config import * 로 사용
"""
import os, glob as _g

# ── 프로젝트 루트 ────────────────────────────────────
def _find_root():
    cur = os.path.dirname(os.path.abspath(__file__))
    for _ in range(6):
        if os.path.exists(os.path.join(cur, "src", "config.py")):
            return cur
        cur = os.path.dirname(cur)
    return cur

ROOT = _find_root()

# ── 경로 ────────────────────────────────────────────
DATA_RAW       = os.path.join(ROOT, "data", "raw")
DATA_PROCESSED = os.path.join(ROOT, "data", "processed")
MODEL_DIR      = os.path.join(ROOT, "model", "iCHO3K")
LOGS_DIR       = os.path.join(ROOT, "logs")

SOWA_MMC1  = os.path.join(DATA_RAW, "sowa2020", "mmc1.xlsx")
AA20_EXCEL = os.path.join(DATA_RAW, "practice_20aa",
                          "CHO_raw_data_practice_20AA.xlsx")
OWN_DIR    = os.path.join(DATA_RAW, "own_experiment")

# ── 모델 경로 (prod 우선) ────────────────────────────
_prod = _g.glob(os.path.join(MODEL_DIR, "*prod*.json"))
MODEL_PATH = _prod[0] if _prod else next(
    iter(_g.glob(os.path.join(MODEL_DIR, "*.json"))), None)

# ── 결과 폴더 ───────────────────────────────────────
def results_dir(dataset: str, subdir: str = "tables") -> str:
    p = os.path.join(ROOT, "results", dataset, subdir)
    os.makedirs(p, exist_ok=True)
    return p

# ── FBA 파라미터 ─────────────────────────────────────
# Ref: Orth et al. (2010) Nat Biotechnol 28:245
# biomass_cho_prod: igg_formation이 모델에서 연결 끊김 → prod biomass 사용
# Ref: Gopalakrishnan et al. (2024) Metab Eng 82:110
OBJ_RXN       = "biomass_cho_prod"
FVA_FRACTION  = 0.9   # Ref: Mahadevan & Schilling (2003) Metab Eng 5:264
FVA_PROCESSES = 1     # Windows=1, Linux=os.cpu_count()

# ── 검증된 iCHO3K prod Exchange ID ──────────────────
# Ref: 실제 모델 파일 검증 (find_exchange_ids.py)
EXCHANGE_IDS = {
    "Glucose":       "EX_glc_e",     # D-glucose (기본 lb=0 → apply_bounds 필수)
    "Lactate":       "EX_lac_L_e",   # (S)-lactate
    "Glutamine":     "EX_gln_L_e",
    "Glutamate":     "EX_glu_L_e",
    "Alanine":       "EX_ala_L_e",
    "Arginine":      "EX_arg_L_e",
    "Asparagine":    "EX_asn_L_e",
    "Aspartate":     "EX_asp_L_e",
    "Cysteine":      "EX_cys_L_e",
    "Glycine":       "EX_gly_e",
    "Histidine":     "EX_his_L_e",
    "Isoleucine":    "EX_ile_L_e",
    "Leucine":       "EX_leu_L_e",
    "Lysine":        "EX_lys_L_e",
    "Methionine":    "EX_met_L_e",
    "Phenylalanine": "EX_phe_L_e",
    "Proline":       "EX_pro_L_e",
    "Serine":        "EX_ser_L_e",
    "Threonine":     "EX_thr_L_e",
    "Tryptophan":    "EX_trp_L_e",
    "Tyrosine":      "EX_tyr_L_e",
    "Valine":        "EX_val_L_e",
    "NH4+":          "EX_nh4_e",
    "Pyruvate":      "EX_pyr_e",
    "Citrate":       "EX_cit_e",
    "Fumarate":      "EX_fum_e",
    "Succinate":     "EX_succ_e",
    "Malate":        "EX_mal_L_e",
}

# ── 문헌 Constraints ─────────────────────────────────
# Ref: Fouladiha et al. (2020) Bioprocess Biosyst Eng 44:1
LITERATURE_HIGH = {
    "EX_glc_e":   (-0.048, 0.0),   # uptake only → ub=0
    "EX_lac_L_e": (0.0,    0.040),  # secretion only → lb=0
    "EX_gln_L_e": (-0.012, 0.0),
    "EX_nh4_e":   (0.0,    0.010),
    "EX_ala_L_e": (0.0,    0.006),
    "EX_glu_L_e": (-0.002, 0.002),
}
LITERATURE_LOW = {
    "EX_glc_e":   (-0.055, 0.0),
    "EX_lac_L_e": (0.0,    0.090),
    "EX_gln_L_e": (-0.015, 0.0),
    "EX_nh4_e":   (0.0,    0.020),
    "EX_ala_L_e": (0.0,    0.012),
    "EX_glu_L_e": (0.0,    0.006),
}

# ── 색상 팔레트 (Nature/Science 스타일) ──────────────
PALETTE = {
    "high":    "#D62728",
    "low":     "#1F77B4",
    "mid":     ["#FF7F0E","#2CA02C","#9467BD","#8C564B",
                "#E377C2","#7F8C8D","#BCBD22","#17BECF"],
    "impr":    "#2CA02C",
    "neutral": "#AAAAAA",
    "danger":  "#C0392B",
}

# ── 폴더 생성 ────────────────────────────────────────
for _d in [DATA_PROCESSED, LOGS_DIR]:
    os.makedirs(_d, exist_ok=True)

if __name__ == "__main__":
    print(f"ROOT  : {ROOT}")
    print(f"MODEL : {MODEL_PATH}  {'OK' if MODEL_PATH and os.path.exists(MODEL_PATH) else 'MISSING'}")
    print(f"SOWA  : {SOWA_MMC1}  {'OK' if os.path.exists(SOWA_MMC1) else 'MISSING'}")
    print(f"20AA  : {AA20_EXCEL}  {'OK' if os.path.exists(AA20_EXCEL) else 'MISSING'}")
