"""
create_project.py
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
CHO_METABOLOMICS 팀 공유 프로젝트 완전 새로 생성

사용:
  python create_project.py
  python create_project.py --root C:\\Projects\\CHO_METABOLOMICS

- 팀원 2~5명 공유 구조
- 단계별 스크립트 (01_, 02_, 03_ ...)
- 데이터셋별 결과 분리
- 자체 실험 데이터 슬롯 포함
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
import os, sys, shutil, argparse, platform, textwrap

parser = argparse.ArgumentParser()
parser.add_argument("--root", default=None)
parser.add_argument("--old",  default=None, help="기존 CHO_FBA 경로 (기본: 자동 탐색)")
args = parser.parse_args()

IS_WIN = platform.system() == "Windows"
ROOT   = args.root or (r"C:\CHO_METABOLOMICS" if IS_WIN else os.path.expanduser("~/CHO_METABOLOMICS"))

OLD_CANDIDATES = [r"C:\CHO_FBA", r"C:\CHO_POC",
                  os.path.expanduser("~/CHO_FBA"), os.path.expanduser("~/CHO_POC")]
OLD = args.old or next((p for p in OLD_CANDIDATES if os.path.exists(p)), None)

SEP = "─" * 60

def section(title):
    print(f"\n{SEP}")
    print(f"  {title}")
    print(SEP)

def ok(msg):  print(f"  ✔  {msg}")
def warn(msg):print(f"  ⚠  {msg}")
def skip(msg):print(f"  –  {msg} (건너뜀)")

def make(rel, comment=""):
    full = os.path.join(ROOT, rel)
    os.makedirs(full, exist_ok=True)
    open(os.path.join(full, ".gitkeep"), "w").close()
    pad = " " * max(1, 45 - len(rel))
    print(f"  📁 {rel}{pad}{comment}")
    return full

def write(rel, content, note=""):
    full = os.path.join(ROOT, rel)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, "w", encoding="utf-8") as f:
        f.write(textwrap.dedent(content).lstrip())
    ok(f"{rel}  {note}")

def copy_file(src, dst_dir, rename=None):
    if not os.path.exists(src): return False
    os.makedirs(dst_dir, exist_ok=True)
    dst = os.path.join(dst_dir, rename or os.path.basename(src))
    if os.path.exists(dst):
        skip(os.path.basename(dst))
        return True
    shutil.copy2(src, dst)
    ok(f"복사: {os.path.basename(src)}")
    return True

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
print(f"\n{'━'*60}")
print("  CHO_METABOLOMICS 팀 프로젝트 생성")
print(f"{'━'*60}")
print(f"  루트   : {ROOT}")
print(f"  기존   : {OLD or '없음'}")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 1: 폴더 구조
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
section("STEP 1. 폴더 구조 생성")

# 데이터 (원본 절대 수정 금지)
make("data/raw/sowa2020",          "Sowa et al. 2020  ← 원본, 수정 금지")
make("data/raw/practice_20aa",     "CHO_raw_data_practice_20AA.xlsx")
make("data/raw/own_experiment",    "자체 fed-batch 실험 원본 (추후)")
make("data/raw/literature_rates",  "Fouladiha 2020, Gopalakrishnan 2024 CSV")
make("data/processed",             "rate 계산 후 저장 (pkl / csv)")

# 모델
make("model/iCHO3K",               "iCHO3K JSON 파일만 (repo 전체 X)")

# 공통 모듈
make("src",                        "패키지처럼 import 가능한 공통 코드")

# 단계별 스크립트 (데이터셋 독립)
make("scripts",                    "진입점 스크립트 (run_*.py)")
make("scripts/steps",              "01_ ~ 06_ 단계별 분석 모듈")
make("scripts/diagnostics",        "모델 탐색·디버그 유틸")

# 결과 (데이터셋 × 분석 유형)
for ds in ["sowa2020", "practice_20aa", "own_experiment"]:
    make(f"results/{ds}/tables",   "CSV / Excel 결과")
    make(f"results/{ds}/figures",  "PNG 300 DPI 논문 Figure")
    make(f"results/{ds}/reports",  "요약 리포트 (md / html)")

# 노트북, 문서, 로그
make("notebooks",                  "탐색 EDA Jupyter (공유 전 검토)")
make("docs",                       "팀 메모·분석 메서드 노트")
make("logs",                       "파이프라인 실행 로그")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 2: 기존 데이터 복사
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
section("STEP 2. 기존 파일 복사")

if OLD:
    # mmc1.xlsx
    for c in [os.path.join(OLD, "data", "sowa2020", "mmc1.xlsx"),
              os.path.join(OLD, "data", "mmc1.xlsx")]:
        if copy_file(c, os.path.join(ROOT, "data", "raw", "sowa2020")): break

    # 20AA Excel
    for c in [os.path.join(OLD, "data", "CHO_raw_data_practice_20AA.xlsx"),
              os.path.join(OLD, "data", "metabolomics", "raw", "CHO_raw_data_practice_20AA.xlsx")]:
        if copy_file(c, os.path.join(ROOT, "data", "raw", "practice_20aa")): break

    # iCHO3K JSON 파일만
    for sub in [os.path.join("model","iCHO3K-main","iCHO3K","Model"),
                os.path.join("model","iCHO3K-main","Model")]:
        d = os.path.join(OLD, sub)
        if os.path.exists(d):
            for f in os.listdir(d):
                if f.endswith(".json"):
                    copy_file(os.path.join(d, f), os.path.join(ROOT, "model", "iCHO3K"))
            break

    # 기존 분석 스크립트 → scripts/steps/ 로 이동
    old_scripts = os.path.join(OLD, "scripts")
    if os.path.exists(old_scripts):
        for f in sorted(os.listdir(old_scripts)):
            if f.endswith(".py") and (f[0].isdigit() or f.startswith("run_")):
                copy_file(os.path.join(old_scripts, f),
                          os.path.join(ROOT, "scripts", "steps"))

    # 진단 스크립트
    for f in ["find_exchange_ids.py","fix_constraints.py","debug_model.py","debug_excel.py"]:
        copy_file(os.path.join(OLD, f), os.path.join(ROOT, "scripts", "diagnostics"))
else:
    warn("기존 CHO_FBA 폴더 없음 — 파일 복사 건너뜀")
    warn("모델은 수동 배치 필요: model/iCHO3K/")
    warn("데이터는 수동 배치 필요: data/raw/")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 3: 핵심 소스 파일 생성
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
section("STEP 3. 공통 모듈 생성 (src/)")

# src/__init__.py
write("src/__init__.py", "")

# src/config.py
write("src/config.py", f'''
    """
    config.py — CHO_METABOLOMICS 프로젝트 전역 설정
    팀원 모두 여기에서 경로·파라미터를 참조합니다.
    """
    import os, glob as _g

    ROOT = r"{ROOT}"

    # ── 경로 ───────────────────────────────────────────
    DATA_RAW        = os.path.join(ROOT, "data", "raw")
    DATA_PROCESSED  = os.path.join(ROOT, "data", "processed")
    MODEL_DIR       = os.path.join(ROOT, "model", "iCHO3K")
    LOGS_DIR        = os.path.join(ROOT, "logs")

    # 원본 데이터 경로
    SOWA_MMC1  = os.path.join(DATA_RAW, "sowa2020", "mmc1.xlsx")
    AA20_EXCEL = os.path.join(DATA_RAW, "practice_20aa",
                              "CHO_raw_data_practice_20AA.xlsx")
    OWN_DIR    = os.path.join(DATA_RAW, "own_experiment")   # 추후

    # 모델 (prod 우선)
    _prod = _g.glob(os.path.join(MODEL_DIR, "*prod*.json"))
    MODEL_PATH = _prod[0] if _prod else next(
        iter(_g.glob(os.path.join(MODEL_DIR, "*.json"))), None)

    # 결과 폴더 함수
    def results_dir(dataset: str, subdir: str = "tables") -> str:
        """
        dataset: "sowa2020" | "practice_20aa" | "own_experiment"
        subdir:  "tables" | "figures" | "reports"
        """
        p = os.path.join(ROOT, "results", dataset, subdir)
        os.makedirs(p, exist_ok=True)
        return p

    # ── FBA 파라미터 ────────────────────────────────────
    OBJ_RXN       = "biomass_cho_prod"  # igg_formation은 연결 끊김
    FVA_FRACTION  = 0.9
    FVA_PROCESSES = 1                   # Windows=1, Linux=코어수

    # ── 검증된 iCHO3K prod Exchange ID ─────────────────
    EXCHANGE_IDS = {{
        "Glucose":       "EX_glc_e",
        "Lactate":       "EX_lac_L_e",
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
    }}

    # ── 문헌 Constraints (Fouladiha 2020) ───────────────
    LITERATURE_HIGH = {{
        "EX_glc_e":   (-0.048, -0.032), "EX_lac_L_e": ( 0.024,  0.040),
        "EX_gln_L_e": (-0.012, -0.007), "EX_nh4_e":   ( 0.005,  0.010),
        "EX_ala_L_e": ( 0.002,  0.006), "EX_glu_L_e": (-0.002,  0.002),
    }}
    LITERATURE_LOW = {{
        "EX_glc_e":   (-0.055, -0.040), "EX_lac_L_e": ( 0.060,  0.090),
        "EX_gln_L_e": (-0.015, -0.010), "EX_nh4_e":   ( 0.012,  0.020),
        "EX_ala_L_e": ( 0.005,  0.012), "EX_glu_L_e": ( 0.001,  0.006),
    }}

    # ── 색상 팔레트 (Nature/Science 스타일) ─────────────
    PALETTE = {{
        "high":    "#D62728", "low":     "#1F77B4",
        "mid":     ["#FF7F0E","#2CA02C","#9467BD","#8C564B",
                   "#E377C2","#7F7F7F","#BCBD22","#17BECF"],
        "impr":    "#2CA02C", "neutral": "#AAAAAA", "danger":  "#C0392B",
    }}

    if __name__ == "__main__":
        for k, v in {{"MODEL": MODEL_PATH, "SOWA": SOWA_MMC1, "20AA": AA20_EXCEL}}.items():
            status = "✔" if v and os.path.exists(v) else "✘"
            print(f"  {{status}}  {{k}}: {{v}}")
''')

# src/fba_utils.py
write("src/fba_utils.py", '''
    """
    fba_utils.py — FBA/FVA 공통 함수
    모든 스크립트에서 from src.fba_utils import * 로 사용
    """
    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    from src.config import *

    import cobra, pandas as pd, numpy as np
    from cobra.flux_analysis import flux_variability_analysis


    def load_model(verbose=True):
        """iCHO3K prod 모델 로딩 + 목적함수 설정"""
        if not MODEL_PATH or not os.path.exists(MODEL_PATH):
            raise FileNotFoundError(
                f"모델 없음: {MODEL_PATH}\\n"
                "model/iCHO3K/ 폴더에 iCHO3K_cho_prod_generic_unblocked.json 배치 필요"
            )
        m = cobra.io.load_json_model(MODEL_PATH) if MODEL_PATH.endswith(".json") \\
            else cobra.io.read_sbml_model(MODEL_PATH)
        m.objective = OBJ_RXN
        if verbose:
            print(f"  Model: {os.path.basename(MODEL_PATH)}")
            print(f"  Reactions: {len(m.reactions)} | Genes: {len(m.genes)}")
            print(f"  Objective: {OBJ_RXN}")
        return m


    def apply_bounds(model, constraints: dict) -> int:
        """
        constraints를 모델에 적용
        clamping 없이 직접 설정 — EX_glc_e(lb=0 기본)도 음수 설정 가능

        Args:
            constraints: {rxn_id: (lb, ub)} or {met_name: (lb, ub)}
                         met_name은 EXCHANGE_IDS에서 자동 변환
        Returns:
            적용된 constraint 수
        """
        rxn_ids = {r.id for r in model.reactions}
        n = 0
        for key, (lb, ub) in constraints.items():
            rxn_id = EXCHANGE_IDS.get(key, key)
            if rxn_id not in rxn_ids:
                continue
            rxn = model.reactions.get_by_id(rxn_id)
            rxn.lower_bound = float(lb)
            rxn.upper_bound = float(ub)
            n += 1
        return n


    def rates_to_constraints(rates: dict, buffer: float = 0.20) -> dict:
        """
        측정 exchange rate → FBA bounds 변환
        Args:
            rates:  {met_name: rate_mmol_gDCWh}  (neg=uptake, pos=secretion)
            buffer: 측정 불확실성 ±% (기본 ±20%)
        Returns:
            {exchange_id: (lb, ub)}
        """
        cst = {}
        for key, q in rates.items():
            rxn_id = EXCHANGE_IDS.get(key, key)
            if q == 0:
                continue
            if q < 0:   # uptake
                lb, ub = q * (1 + buffer), q * (1 - buffer)
            else:       # secretion
                lb, ub = q * (1 - buffer), q * (1 + buffer)
            cst[rxn_id] = (min(lb, ub), max(lb, ub))
        return cst


    def run_fba(model, constraints: dict, label: str = "") -> dict:
        """
        FBA 실행 (context manager 사용 → 모델 원본 보존)
        Returns:
            {"label", "obj", "status", "ex_flux", "n_constraints"}
        """
        with model:
            n = apply_bounds(model, constraints)
            model.objective = OBJ_RXN
            sol = model.optimize()
            obj = sol.objective_value if sol.status == "optimal" else 0.0
            ex  = {r.id: sol.fluxes[r.id] for r in model.exchanges
                   if abs(sol.fluxes.get(r.id, 0)) > 1e-9} if sol.status == "optimal" else {}
        if label:
            print(f"  {label:42s} n={n:2d}  {sol.status:10s}  obj={obj:.5f}")
        return {"label": label, "obj": obj, "status": sol.status,
                "ex_flux": ex, "n_constraints": n}


    def run_fva(model, constraints: dict, label: str = "",
                exchange_ids: list = None) -> pd.DataFrame:
        """
        FVA 실행 (fraction=FVA_FRACTION)
        Returns:
            DataFrame with columns: minimum, maximum, mean, range
        """
        with model:
            apply_bounds(model, constraints)
            model.objective = OBJ_RXN
            rxn_list = None
            if exchange_ids:
                rxn_list = [model.reactions.get_by_id(r)
                            for r in exchange_ids
                            if r in {x.id for x in model.reactions}]
            try:
                fva = flux_variability_analysis(
                    model, fraction_of_optimum=FVA_FRACTION,
                    processes=FVA_PROCESSES,
                    reaction_list=rxn_list or list(model.exchanges)
                )
                fva["mean"]  = (fva["minimum"] + fva["maximum"]) / 2
                fva["range"] = fva["maximum"] - fva["minimum"]
                if label:
                    print(f"  FVA {label}: OK ({len(fva)} reactions)")
                return fva
            except Exception as e:
                print(f"  FVA {label}: {e}")
                return pd.DataFrame()


    def compute_ivcd(vcd1, vcd2, day_start, day_end,
                     volume_mL, cell_dw_g=8e-12) -> float:
        """
        IVCD [gDCW·h] = avg_VCD × vol × Δt × cell_weight
        Args:
            vcd1, vcd2: viable cell density at day_start, day_end [10^6 cells/mL]
        """
        return ((vcd1 + vcd2) / 2) * 1e6 * volume_mL * (day_end - day_start) * 24 * cell_dw_g


    def save_figure(fig, name: str, dataset: str, subdir: str = "figures") -> str:
        """
        Figure를 results/{dataset}/{subdir}/{name} 에 300 DPI 저장
        사용: save_figure(fig, "Fig1_timecourse.png", dataset="practice_20aa")
        """
        import matplotlib.pyplot as plt
        out = results_dir(dataset, subdir)
        path = os.path.join(out, name)
        fig.savefig(path, dpi=300, bbox_inches="tight",
                    facecolor="white", edgecolor="none")
        plt.close(fig)
        print(f"  [saved] results/{dataset}/{subdir}/{name}")
        return path


    def save_table(df: pd.DataFrame, name: str, dataset: str) -> str:
        """
        DataFrame을 results/{dataset}/tables/{name} 에 저장
        """
        out = results_dir(dataset, "tables")
        path = os.path.join(out, name)
        if name.endswith(".csv"):
            df.to_csv(path, index=False)
        elif name.endswith(".xlsx"):
            df.to_excel(path, index=False)
        print(f"  [saved] results/{dataset}/tables/{name}")
        return path
''')

print()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 4: 단계별 스크립트 뼈대 생성
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
section("STEP 4. 단계별 스크립트 생성 (scripts/steps/)")

STEPS = [
    ("01_load_data.py",      "데이터 로딩 + Feed-Corrected Rate 계산"),
    ("02_map_metabolites.py","대사물질명 → iCHO3K Exchange ID 매핑"),
    ("03_run_fba.py",        "FBA 실행 (클론별 / 데이터셋별)"),
    ("04_run_fva.py",        "FVA 실행 (High vs Low producer)"),
    ("05_ko_screening.py",   "Reaction KO In Silico 스크리닝"),
    ("06_figures.py",        "논문 Figure 생성 (8종)"),
]

for fname, desc in STEPS:
    step_num = fname[:2]
    path = os.path.join(ROOT, "scripts", "steps", fname)
    content = f'''"""
{fname} — {desc}

사용:
  cd {ROOT}
  python scripts/steps/{fname} --dataset practice_20aa

  dataset 옵션: sowa2020 | practice_20aa | own_experiment
"""
import sys, os, argparse
sys.path.insert(0, r"{ROOT}")
from src.config import *
from src.fba_utils import *

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa",
                    choices=["sowa2020", "practice_20aa", "own_experiment"])
args = parser.parse_args()
DATASET = args.dataset

print("=" * 60)
print(f"  {fname}")
print(f"  Dataset: {{DATASET}}")
print("=" * 60)

# ── TODO: 이 아래에 분석 코드 작성 ─────────────────────
# 이전 스크립트 결과 로딩 예시:
# import pickle
# with open(os.path.join(DATA_PROCESSED, f"rates_{{DATASET}}.pkl"), "rb") as f:
#     data = pickle.load(f)

# 결과 저장 예시:
# save_table(df, "exchange_rates.csv", dataset=DATASET)
# save_figure(fig, "Fig1_timecourse.png", dataset=DATASET)

raise NotImplementedError("{desc} — scripts/steps/{fname} 에 코드를 작성하세요")
'''
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"  📄 scripts/steps/{fname:<30s} # {desc}")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 5: 실행 진입점 스크립트
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
section("STEP 5. 실행 진입점 생성")

# run_pipeline.py (루트 레벨)
write("run_pipeline.py", f'''
    """
    run_pipeline.py — 전체 파이프라인 한 번에 실행
    사용: python run_pipeline.py --dataset practice_20aa
    """
    import subprocess, sys, os, time, argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", default="practice_20aa",
                        choices=["sowa2020", "practice_20aa", "own_experiment"])
    parser.add_argument("--steps", default="1,2,3,4,5,6",
                        help="실행할 스텝 번호 (예: 1,2,3)")
    args = parser.parse_args()

    ROOT    = r"{ROOT}"
    STEPS_DIR = os.path.join(ROOT, "scripts", "steps")
    LOGS_DIR  = os.path.join(ROOT, "logs")
    os.makedirs(LOGS_DIR, exist_ok=True)

    STEPS = {{
        1: "01_load_data.py",
        2: "02_map_metabolites.py",
        3: "03_run_fba.py",
        4: "04_run_fva.py",
        5: "05_ko_screening.py",
        6: "06_figures.py",
    }}
    DESCS = {{
        1: "데이터 로딩 + Rate 계산",
        2: "대사물질 매핑",
        3: "FBA 실행",
        4: "FVA 실행",
        5: "KO 스크리닝",
        6: "Figure 생성",
    }}

    to_run = [int(x) for x in args.steps.split(",")]
    env    = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"

    print("=" * 60)
    print(f"  CHO_METABOLOMICS Pipeline — {{args.dataset}}")
    print("=" * 60)
    t_total = time.time()

    for num in to_run:
        script = STEPS.get(num)
        if not script:
            continue
        path = os.path.join(STEPS_DIR, script)
        print(f"\\n[{{num}}/6] {{DESCS[num]}}")
        t0 = time.time()
        log_path = os.path.join(LOGS_DIR, script.replace(".py", f"_{{args.dataset}}.log"))
        result = subprocess.run(
            [sys.executable, path, "--dataset", args.dataset],
            cwd=ROOT, capture_output=True, text=True,
            encoding="utf-8", errors="replace", env=env
        )
        with open(log_path, "w", encoding="utf-8") as lf:
            lf.write(result.stdout + "\\n" + result.stderr)
        for line in result.stdout.split("\\n"):
            if any(k in line for k in ["✔","✘","saved","obj=","status=","Error","완료"]):
                print(f"    {{line.rstrip()}}")
        elapsed = time.time() - t0
        if result.returncode == 0:
            print(f"  ✔ 완료 ({{elapsed:.1f}}s)  로그: {{os.path.basename(log_path)}}")
        else:
            print(f"  ✘ 오류 (code={{result.returncode}}, {{elapsed:.1f}}s)")
            print(f"     로그: {{log_path}}")
            if input("  계속? (y/n): ").strip().lower() != "y":
                break

    print(f"\\n{'='*60}")
    print(f"  완료! 소요: {{time.time()-t_total:.1f}}초")
    print(f"  결과: {{ROOT}}\\\\results\\\\{{args.dataset}}")
    print(f"{'='*60}")
''')

# setup_env.py (환경 자동 설정)
write("setup_env.py", f'''
    """
    setup_env.py — 팀원 환경 초기 설정
    처음 한 번만 실행: python setup_env.py
    """
    import subprocess, sys, os

    print("=" * 60)
    print("  CHO_METABOLOMICS 환경 설정")
    print("=" * 60)

    PKGS = ["cobra","pandas","numpy","scipy","matplotlib",
            "seaborn","openpyxl","requests","tqdm","scikit-learn"]

    missing = []
    for pkg in PKGS:
        imp = {{"scikit-learn": "sklearn"}}.get(pkg, pkg)
        try:
            __import__(imp); print(f"  ✔  {{pkg}}")
        except ImportError:
            print(f"  ✘  {{pkg}}"); missing.append(pkg)

    if missing:
        print(f"\\n  설치 중: {{missing}}")
        subprocess.check_call([sys.executable, "-m", "pip", "install"] + missing)
        print("  ✔ 설치 완료")

    # iCHO3K 모델 다운로드
    import zipfile, io, requests
    MODEL_DIR = r"{os.path.join(ROOT, 'model', 'iCHO3K')}"
    os.makedirs(MODEL_DIR, exist_ok=True)
    prod_json = next((f for f in os.listdir(MODEL_DIR) if "prod" in f and f.endswith(".json")), None)
    if prod_json:
        print(f"\\n  ✔ 모델 이미 있음: {{prod_json}}")
    else:
        print("\\n  iCHO3K 다운로드 중...")
        url = "https://github.com/LewisLabUCSD/iCHO3K/archive/refs/heads/main.zip"
        try:
            r = requests.get(url, timeout=60)
            r.raise_for_status()
            with zipfile.ZipFile(io.BytesIO(r.content)) as z:
                for member in z.namelist():
                    if "iCHO3K/Model/" in member and member.endswith(".json"):
                        fname = os.path.basename(member)
                        with z.open(member) as src, open(os.path.join(MODEL_DIR, fname), "wb") as dst:
                            dst.write(src.read())
                        print(f"  ✔ {{fname}}")
        except Exception as e:
            print(f"  ✘ 자동 다운로드 실패: {{e}}")
            print("  수동: github.com/LewisLabUCSD/iCHO3K → JSON 파일 → model/iCHO3K/")

    # config 검증
    sys.path.insert(0, r"{ROOT}")
    from src.config import MODEL_PATH, SOWA_MMC1, AA20_EXCEL
    print("\\n  [파일 확인]")
    for k, v in [("모델", MODEL_PATH), ("Sowa MMC1", SOWA_MMC1), ("20AA Excel", AA20_EXCEL)]:
        status = "✔" if v and os.path.exists(v) else "✘"
        print(f"  {{status}}  {{k}}: {{v or '경로 없음'}}")

    print("\\n  설정 완료! 다음 실행:")
    print(f"  python run_pipeline.py --dataset practice_20aa")
''')

print()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 6: 팀 공유 파일
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
section("STEP 6. 팀 공유 파일 생성")

# requirements.txt
write("requirements.txt", '''
    cobra>=0.26
    pandas>=2.0
    numpy>=1.24
    scipy>=1.10
    matplotlib>=3.7
    seaborn>=0.12
    openpyxl>=3.1
    scikit-learn>=1.3
    requests>=2.28
    tqdm>=4.65
''')

# .gitignore
write(".gitignore", f'''
    # Python
    __pycache__/
    *.pyc
    *.pyo
    .env
    *.egg-info/

    # 대용량 원본 데이터 (Git LFS 또는 별도 공유 드라이브 사용)
    data/raw/
    model/

    # 중간 결과 (pkl은 용량이 클 수 있음)
    data/processed/*.pkl

    # 결과 PNG (선택 — 공유하려면 주석 처리)
    # results/**/*.png

    # Jupyter 체크포인트
    .ipynb_checkpoints/
    notebooks/.ipynb_checkpoints/

    # OS
    .DS_Store
    Thumbs.db
    logs/
''')

# README.md
write("README.md", f'''
    # CHO_METABOLOMICS

    CHO 항체 생산 세포주 Metabolomics + FBA/FVA 분석 팀 프로젝트

    ## 빠른 시작

    ```bash
    # 1. 환경 설정 (처음 한 번만)
    conda create -n cho python=3.10 -y
    conda activate cho
    pip install -r requirements.txt
    python setup_env.py        # 패키지 확인 + iCHO3K 자동 다운로드

    # 2. 데이터 배치 (공유 드라이브 또는 수동)
    #    data/raw/sowa2020/mmc1.xlsx
    #    data/raw/practice_20aa/CHO_raw_data_practice_20AA.xlsx

    # 3. 파이프라인 실행
    python run_pipeline.py --dataset practice_20aa
    python run_pipeline.py --dataset sowa2020
    ```

    ## 폴더 구조

    ```
    CHO_METABOLOMICS/
    ├── data/
    │   ├── raw/
    │   │   ├── sowa2020/           Sowa et al. 2020 (원본 수정 금지)
    │   │   ├── practice_20aa/      20AA 연습 데이터
    │   │   ├── own_experiment/     자체 실험 원본
    │   │   └── literature_rates/   Fouladiha/Gopalakrishnan CSV
    │   └── processed/              rate 계산 중간 결과 (pkl)
    ├── model/
    │   └── iCHO3K/                 iCHO3K JSON 파일 (setup_env.py로 자동 다운로드)
    ├── src/
    │   ├── config.py               경로·파라미터 전역 설정
    │   └── fba_utils.py            FBA 공통 함수 (load_model, run_fba, ...)
    ├── scripts/
    │   ├── steps/
    │   │   ├── 01_load_data.py     데이터 로딩 + rate 계산
    │   │   ├── 02_map_metabolites.py
    │   │   ├── 03_run_fba.py
    │   │   ├── 04_run_fva.py
    │   │   ├── 05_ko_screening.py
    │   │   └── 06_figures.py
    │   └── diagnostics/            모델 탐색·디버그
    ├── results/
    │   ├── sowa2020/tables/        exchange_rates.csv, fba_results.csv ...
    │   ├── sowa2020/figures/       Fig1~8 PNG (300 DPI)
    │   ├── practice_20aa/
    │   └── own_experiment/         자체 실험 결과 (추후)
    ├── notebooks/                  EDA Jupyter (팀 검토용)
    ├── docs/                       분석 메서드 노트
    ├── run_pipeline.py             전체 파이프라인 실행
    ├── setup_env.py                환경 초기 설정
    └── requirements.txt
    ```

    ## 스크립트 단계

    | 스텝 | 파일 | 역할 |
    |------|------|------|
    | 01 | `01_load_data.py` | Excel 파싱, feed 보정, IVCD 기반 specific rate 계산 |
    | 02 | `02_map_metabolites.py` | 대사물질명 → iCHO3K Exchange ID 매핑 |
    | 03 | `03_run_fba.py` | 클론별 FBA, 목적함수=biomass_cho_prod |
    | 04 | `04_run_fva.py` | High vs Low FVA, 90% optimality |
    | 05 | `05_ko_screening.py` | 반응 KO 스크리닝 |
    | 06 | `06_figures.py` | 논문 수준 Figure 8종 생성 |

    ## 중요: Exchange Rate 계산

    ```
    q [mmol/gDCW/h] = (ΔC_obs + ΔC_feed) × Volume / IVCD
    IVCD = avg_VCD × Volume × Δt × 8e-12 g/cell
    음수(-)  = uptake (세포가 흡수)
    양수(+)  = secretion (세포가 분비)
    ```

    ## 모델 핵심 파라미터

    - 파일: `iCHO3K_cho_prod_generic_unblocked.json`
    - 목적함수: `biomass_cho_prod` (`igg_formation`은 연결 끊김)
    - Glucose: `EX_glc_e` (기본 lb=0 → apply_bounds로 음수 설정 필수)
    - Lactate: `EX_lac_L_e`

    ## 팀 컨벤션

    - `data/raw/` — **절대 수정 금지** (원본 보존)
    - 결과는 반드시 `save_table()` / `save_figure()` 함수로 저장
    - 새 데이터셋 추가 시: `data/raw/<name>/` + `results/<name>/` 폴더 생성
    - 스크립트는 `--dataset` 인수로 데이터셋 분기
''')

# docs/method_notes.md
write("docs/method_notes.md", '''
    # 분석 메서드 노트

    ## FBA 개요

    iCHO3K prod 모델을 사용한 Flux Balance Analysis.
    - 목적함수: `biomass_cho_prod` (항체 생산 특화 biomass)
    - Exchange constraint: feed-corrected specific rate ± 20%
    - Infeasibility 원인 대부분: `EX_glc_e`의 기본 lb=0 → `apply_bounds()` 필수

    ## Exchange Rate 계산

    ```
    q = (ΔC_obs - ΔC_feed) / IVCD
    IVCD = avg_VCD × volume × Δt × 8e-12 g/cell
    ```

    ## 검증된 Exchange ID (iCHO3K prod)

    | 대사물질 | Exchange ID | 주의 |
    |---------|-------------|------|
    | Glucose | `EX_glc_e` | lb=0 기본, 반드시 apply_bounds 사용 |
    | Lactate | `EX_lac_L_e` | |
    | Glutamine | `EX_gln_L_e` | |
    | NH4+ | `EX_nh4_e` | |

    ## 알려진 이슈

    - `igg_formation`: prod 모델에서 대사망 연결 끊김 → 사용 불가
    - fuzzy 매핑 오류: glucose→dopa3glcur_e 등 엉뚱한 ID 반환 가능 → EXCHANGE_IDS 직접 매핑 사용
    - FBA infeasible 원인 1순위: glucose uptake constraint 미적용
''')

print()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 최종 리포트
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
section("완료!")
print(f"  루트: {ROOT}\n")

CHECKS = [
    ("data/raw/sowa2020",                   "Sowa 2020 폴더"),
    ("data/raw/practice_20aa",              "20AA 폴더"),
    ("model/iCHO3K",                        "iCHO3K 모델 폴더"),
    ("src/config.py",                       "config.py"),
    ("src/fba_utils.py",                    "fba_utils.py"),
    ("scripts/steps/01_load_data.py",       "01_load_data.py"),
    ("scripts/steps/06_figures.py",         "06_figures.py"),
    ("run_pipeline.py",                     "run_pipeline.py"),
    ("setup_env.py",                        "setup_env.py"),
    ("requirements.txt",                    "requirements.txt"),
    ("README.md",                           "README.md"),
]

all_ok = True
for rel, desc in CHECKS:
    full = os.path.join(ROOT, rel)
    if os.path.isdir(full):
        n = len([f for f in os.listdir(full) if not f.startswith(".")])
        s = f"✔  ({n}개)" if n > 0 else "⚠  비어 있음"
        if n == 0: all_ok = False
    else:
        s = "✔ " if os.path.exists(full) else "✘ 없음"
        if not os.path.exists(full): all_ok = False
    print(f"  {s:12s} {rel}")

print()
print("  ─── 다음 단계 ───────────────────────────────")
print(f"  cd {ROOT}")
print( "  python setup_env.py              # 환경 설정 + 모델 다운로드")
print( "  python run_pipeline.py --dataset practice_20aa")
print()
if not all_ok:
    warn("일부 파일 없음 — setup_env.py 실행으로 자동 보완")
