"""
build_final_pipeline.py
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
CHO_METABOLOMICS 최종 파이프라인 빌드 스크립트
이 스크립트 하나를 실행하면 전체 프로젝트 파일이 생성됩니다.

사용:
  python build_final_pipeline.py
  python build_final_pipeline.py --root C:\\MyProject\\CHO_METABOLOMICS
  python build_final_pipeline.py --old C:\\CHO_FBA   (기존 프로젝트 경로)
"""
import os, sys, shutil, argparse, platform, textwrap

parser = argparse.ArgumentParser()
parser.add_argument("--root", default=None,
                    help="새 프로젝트 경로 (기본: C:\\CHO_METABOLOMICS)")
parser.add_argument("--old",  default=None,
                    help="기존 CHO_FBA 경로 (자동 탐색)")
args = parser.parse_args()

IS_WIN = platform.system() == "Windows"
ROOT   = args.root or (r"C:\CHO_METABOLOMICS" if IS_WIN
                       else os.path.expanduser("~/CHO_METABOLOMICS"))
OLD_CANDIDATES = [r"C:\CHO_FBA", r"C:\CHO_POC",
                  os.path.expanduser("~/CHO_FBA")]
OLD = args.old or next((p for p in OLD_CANDIDATES if os.path.exists(p)), None)

print("="*60)
print("  CHO_METABOLOMICS Final Pipeline Build")
print("="*60)
print(f"  ROOT : {ROOT}")
print(f"  OLD  : {OLD or '없음'}")

def mkd(rel):
    p = os.path.join(ROOT, rel)
    os.makedirs(p, exist_ok=True)
    keep = os.path.join(p, ".gitkeep")
    if not os.listdir(p):
        open(keep, "w").close()

def write(rel, content):
    p = os.path.join(ROOT, rel)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        f.write(textwrap.dedent(content).lstrip())
    print(f"  [OK] {rel}")

def copy_f(src, dst_dir, rename=None):
    if not os.path.exists(src): return False
    os.makedirs(dst_dir, exist_ok=True)
    dst = os.path.join(dst_dir, rename or os.path.basename(src))
    if not os.path.exists(dst):
        shutil.copy2(src, dst)
        print(f"  [CP] {os.path.basename(src)}")
    return True

# ════════════════════════════════════════════════════
# 1. 폴더 구조
# ════════════════════════════════════════════════════
print("\n[1/6] 폴더 구조 생성")
for d in [
    "data/raw/sowa2020", "data/raw/practice_20aa",
    "data/raw/own_experiment", "data/raw/literature_rates",
    "data/processed", "model/iCHO3K", "src",
    "scripts/steps", "scripts/diagnostics",
    "results/sowa2020/tables",      "results/sowa2020/figures",
    "results/practice_20aa/tables", "results/practice_20aa/figures",
    "results/own_experiment/tables","results/own_experiment/figures",
    "notebooks", "docs", "logs",
]:
    mkd(d)

# ════════════════════════════════════════════════════
# 2. 기존 데이터 복사
# ════════════════════════════════════════════════════
print("\n[2/6] 기존 파일 복사")
if OLD:
    copy_f(os.path.join(OLD, "data", "sowa2020", "mmc1.xlsx"),
           os.path.join(ROOT, "data", "raw", "sowa2020"))
    copy_f(os.path.join(OLD, "data", "CHO_raw_data_practice_20AA.xlsx"),
           os.path.join(ROOT, "data", "raw", "practice_20aa"))
    for sub in ["model/iCHO3K-main/iCHO3K/Model",
                "model/iCHO3K-main/Model"]:
        d = os.path.join(OLD, sub)
        if os.path.exists(d):
            for f in os.listdir(d):
                if f.endswith(".json"):
                    copy_f(os.path.join(d, f),
                           os.path.join(ROOT, "model", "iCHO3K"))
            break
else:
    print("  -- 기존 프로젝트 없음, 데이터 수동 배치 필요")

# ════════════════════════════════════════════════════
# 3. src/ 공통 모듈
# ════════════════════════════════════════════════════
print("\n[3/6] src/ 공통 모듈 생성")

write("src/__init__.py", "")

write("src/config.py", f'''
    """
    config.py — CHO_METABOLOMICS 전역 설정
    Ref: Orth et al. (2010) Nat Biotechnol 28:245
    """
    import os, glob as _g

    def _find_root():
        cur = os.path.dirname(os.path.abspath(__file__))
        for _ in range(6):
            if os.path.exists(os.path.join(cur, "src", "config.py")):
                return cur
            cur = os.path.dirname(cur)
        return cur

    ROOT = _find_root()

    DATA_RAW       = os.path.join(ROOT, "data", "raw")
    DATA_PROCESSED = os.path.join(ROOT, "data", "processed")
    MODEL_DIR      = os.path.join(ROOT, "model", "iCHO3K")
    LOGS_DIR       = os.path.join(ROOT, "logs")

    SOWA_MMC1  = os.path.join(DATA_RAW, "sowa2020", "mmc1.xlsx")
    AA20_EXCEL = os.path.join(DATA_RAW, "practice_20aa",
                              "CHO_raw_data_practice_20AA.xlsx")
    OWN_DIR    = os.path.join(DATA_RAW, "own_experiment")

    _prod = _g.glob(os.path.join(MODEL_DIR, "*prod*.json"))
    MODEL_PATH = _prod[0] if _prod else next(
        iter(_g.glob(os.path.join(MODEL_DIR, "*.json"))), None)

    def results_dir(dataset: str, subdir: str = "tables") -> str:
        p = os.path.join(ROOT, "results", dataset, subdir)
        os.makedirs(p, exist_ok=True)
        return p

    # Ref: Gopalakrishnan et al. (2024) Metab Eng 82:110
    # biomass_cho_prod: igg_formation 모델에서 연결 끊김 → prod biomass 사용
    OBJ_RXN       = "biomass_cho_prod"
    FVA_FRACTION  = 0.9   # Ref: Mahadevan & Schilling (2003) Metab Eng 5:264
    FVA_PROCESSES = 1

    # Ref: 실제 iCHO3K prod 모델 검증 (EX_glc_e lb=0 기본 → apply_bounds 필수)
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
        "Malate":        "EX_mal_L_e",
    }}

    # Ref: Fouladiha et al. (2020) Bioprocess Biosyst Eng 44:1
    # ub=0 / lb=0 단방향: Varma & Palsson (1994) Appl Environ Microbiol 60:3724
    LITERATURE_HIGH = {{
        "EX_glc_e":   (-0.048,  0.0),
        "EX_lac_L_e": ( 0.0,    0.040),
        "EX_gln_L_e": (-0.012,  0.0),
        "EX_nh4_e":   ( 0.0,    0.010),
        "EX_ala_L_e": ( 0.0,    0.006),
        "EX_glu_L_e": (-0.002,  0.002),
    }}
    LITERATURE_LOW = {{
        "EX_glc_e":   (-0.055,  0.0),
        "EX_lac_L_e": ( 0.0,    0.090),
        "EX_gln_L_e": (-0.015,  0.0),
        "EX_nh4_e":   ( 0.0,    0.020),
        "EX_ala_L_e": ( 0.0,    0.012),
        "EX_glu_L_e": ( 0.0,    0.006),
    }}

    PALETTE = {{
        "high":    "#D62728", "low":     "#1F77B4",
        "mid":     ["#FF7F0E","#2CA02C","#9467BD","#8C564B",
                    "#E377C2","#7F8C8D","#BCBD22","#17BECF"],
        "impr":    "#2CA02C", "neutral": "#AAAAAA", "danger": "#C0392B",
    }}

    for _d in [DATA_PROCESSED, LOGS_DIR]:
        os.makedirs(_d, exist_ok=True)

    if __name__ == "__main__":
        for k, v in {{"MODEL": MODEL_PATH, "SOWA": SOWA_MMC1,
                      "20AA": AA20_EXCEL}}.items():
            print(f"  {{\'OK\' if v and os.path.exists(v) else \'--\':2s}} {{k}}: {{v}}")
''')

write("src/fba_utils.py", '''
    """
    fba_utils.py — FBA/FVA 공통 함수 (최종본)
    핵심 수정:
      apply_bounds  : cobra bound 순서 충돌 방지
      fix_one_sided : lb/ub 단방향 보정
    Ref: Varma & Palsson (1994) Appl Environ Microbiol 60:3724
    Ref: Orth et al. (2010) Nat Biotechnol 28:245
    """
    import os, sys
    def _find_root():
        cur = os.path.dirname(os.path.abspath(__file__))
        for _ in range(6):
            if os.path.exists(os.path.join(cur, "src", "config.py")): return cur
            cur = os.path.dirname(cur)
        return cur
    ROOT = _find_root()
    if ROOT not in sys.path: sys.path.insert(0, ROOT)

    from src.config import *
    import cobra, pandas as pd, numpy as np
    from cobra.flux_analysis import flux_variability_analysis

    def load_model(verbose=True):
        if not MODEL_PATH or not os.path.exists(MODEL_PATH):
            raise FileNotFoundError(f"모델 없음: {MODEL_PATH}\\npython setup_env.py 실행")
        m = cobra.io.load_json_model(MODEL_PATH) if MODEL_PATH.endswith(".json") \\
            else cobra.io.read_sbml_model(MODEL_PATH)
        m.objective = OBJ_RXN
        if verbose:
            print(f"  Model : {os.path.basename(MODEL_PATH)}")
            print(f"  Rxns  : {len(m.reactions)} | Genes: {len(m.genes)}")
            print(f"  Obj   : {OBJ_RXN}")
        return m

    def apply_bounds(model, constraints: dict) -> int:
        """cobra bound 순서 충돌 방지 — new lb > current ub 시 ub 먼저"""
        rxn_ids = {r.id for r in model.reactions}
        n = 0
        for key, (lb, ub) in constraints.items():
            rxn_id = EXCHANGE_IDS.get(key, key)
            if rxn_id not in rxn_ids: continue
            rxn = model.reactions.get_by_id(rxn_id)
            lb_f, ub_f = float(lb), float(ub)
            if lb_f > ub_f: lb_f, ub_f = ub_f, lb_f
            if lb_f > rxn.upper_bound:
                rxn.upper_bound = ub_f; rxn.lower_bound = lb_f
            elif ub_f < rxn.lower_bound:
                rxn.lower_bound = lb_f; rxn.upper_bound = ub_f
            else:
                rxn.lower_bound = lb_f; rxn.upper_bound = ub_f
            n += 1
        return n

    def fix_one_sided(constraints: dict) -> dict:
        """lb/ub 둘 다 음수 → (lb, 0), 둘 다 양수 → (0, ub)"""
        fixed = {}
        for rxn_id, (lb, ub) in constraints.items():
            if lb < 0 and ub < 0: fixed[rxn_id] = (lb, 0.0)
            elif lb > 0 and ub > 0: fixed[rxn_id] = (0.0, ub)
            else: fixed[rxn_id] = (lb, ub)
        return fixed

    def rates_to_constraints(rates: dict, buffer: float = 0.20) -> dict:
        """rate → FBA bounds 변환 (Goudar et al. 2005 Biotechnol Prog 21:1193)"""
        cst = {}
        for key, q in rates.items():
            rxn_id = EXCHANGE_IDS.get(key, key)
            if abs(q) < 1e-9: continue
            if q < 0: lb, ub = q*(1+buffer), q*(1-buffer)
            else:     lb, ub = q*(1-buffer), q*(1+buffer)
            cst[rxn_id] = (min(lb,ub), max(lb,ub))
        return fix_one_sided(cst)

    def avg_constraints(clone_list, all_constraints):
        cst_list = [all_constraints[c] for c in clone_list if c in all_constraints]
        if not cst_list: return {}
        ids = set()
        for c in cst_list: ids.update(c.keys())
        return {rid: (np.mean([c[rid][0] for c in cst_list if rid in c]),
                       np.mean([c[rid][1] for c in cst_list if rid in c]))
                for rid in ids}

    def run_fba(model, constraints, label=""):
        with model:
            n = apply_bounds(model, constraints)
            model.objective = OBJ_RXN
            sol = model.optimize()
            obj = sol.objective_value if sol.status=="optimal" else 0.0
            ex  = {r.id: sol.fluxes[r.id] for r in model.exchanges
                   if abs(sol.fluxes.get(r.id,0)) > 1e-9} if sol.status=="optimal" else {}
        if label:
            print(f"  {label:42s} n={n:2d}  {sol.status:10s}  obj={obj:.5f}")
        return {"label":label,"obj":obj,"status":sol.status,
                "ex_flux":ex,"n_constraints":n}

    def run_fva(model, constraints, label="", exchange_ids=None):
        """FVA (Mahadevan & Schilling 2003 Metab Eng 5:264)"""
        with model:
            apply_bounds(model, constraints)
            model.objective = OBJ_RXN
            rxn_ids_set = {r.id for r in model.reactions}
            rxn_list = [model.reactions.get_by_id(r) for r in (exchange_ids or [])
                        if r in rxn_ids_set] or None
            try:
                fva = flux_variability_analysis(
                    model, fraction_of_optimum=FVA_FRACTION,
                    processes=FVA_PROCESSES,
                    reaction_list=rxn_list or list(model.exchanges))
                fva["mean"]  = (fva["minimum"]+fva["maximum"])/2
                fva["range"] = fva["maximum"]-fva["minimum"]
                if label: print(f"  FVA {label}: OK ({len(fva)} reactions)")
                return fva
            except Exception as e:
                print(f"  FVA {label}: {e}")
                return pd.DataFrame()

    def save_figure(fig, name, dataset, subdir="figures"):
        import matplotlib.pyplot as plt
        out  = results_dir(dataset, subdir)
        path = os.path.join(out, name)
        fig.savefig(path, dpi=300, bbox_inches="tight",
                    facecolor="white", edgecolor="none")
        plt.close(fig)
        print(f"  [saved] results/{dataset}/{subdir}/{name}")
        return path

    def save_table(df, name, dataset):
        out  = results_dir(dataset, "tables")
        path = os.path.join(out, name)
        df.to_csv(path, index=False) if name.endswith(".csv") else df.to_excel(path, index=False)
        print(f"  [saved] results/{dataset}/tables/{name}")
        return path

    def ax_style(ax):
        ax.xaxis.grid(True,ls=":",color="0.88",zorder=0)
        ax.yaxis.grid(True,ls=":",color="0.88",zorder=0)
        ax.set_axisbelow(True)
        ax.spines[["top","right"]].set_visible(False)
''')

# ════════════════════════════════════════════════════
# 4. 단계별 스크립트 복사
# ════════════════════════════════════════════════════
print("\n[4/6] 단계별 스크립트 복사")

SCRIPT_MAP = {
    "scripts/steps/00_auto_objective.py": None,    # 별도 파일에서
    "scripts/steps/01_load_data.py":      None,
    "scripts/steps/02_map_metabolites.py":None,
    "scripts/steps/03_run_fba.py":        None,
    "scripts/steps/04_run_fva.py":        None,
    "scripts/steps/05_ko_screening.py":   None,
    "scripts/steps/06_figures.py":        None,
    "scripts/steps/fig_qp_timeseries.py": None,
    "scripts/steps/fig_pathway_flux_map.py":None,
    "scripts/steps/multi_interval_fba.py":  None,
}

# 스크립트들의 공통 헤더 수정 함수
def normalize_header(content, script_name):
    """모든 스크립트의 import 헤더를 통일된 형식으로 교체"""
    # 기존 헤더 패턴들 제거 후 표준 헤더로 교체
    import re
    std_header = '''import sys, os, argparse, warnings, pickle
warnings.filterwarnings("ignore")

def _find_root():
    cur = os.path.dirname(os.path.abspath(__file__))
    for _ in range(6):
        if os.path.exists(os.path.join(cur, "src", "config.py")): return cur
        cur = os.path.dirname(cur)
    return cur
ROOT = _find_root()
if ROOT not in sys.path: sys.path.insert(0, ROOT)
from src.config import *
from src.fba_utils import (load_model, apply_bounds, fix_one_sided,
                            rates_to_constraints, avg_constraints,
                            run_fba, run_fva, save_figure, save_table, ax_style)
import pandas as pd
import numpy as np
'''
    # find first import after docstring
    lines = content.split("\n")
    doc_end = 0
    if lines[0].startswith('"""'):
        for i,l in enumerate(lines[1:],1):
            if '"""' in l:
                doc_end = i+1; break
    # replace header up to first real code
    header_end = doc_end
    for i in range(doc_end, min(doc_end+30, len(lines))):
        l = lines[i]
        if (l.startswith("import sys") or l.startswith("def _find_root") or
            l.startswith("ROOT =") or l.startswith("from src")):
            header_end = i; break
    # find where actual code starts (after all imports/setup)
    code_start = header_end
    for i in range(header_end, min(header_end+40, len(lines))):
        l = lines[i]
        if (l.startswith("parser") or l.startswith("print(") or
            l.startswith("# ──") or l.startswith("DATASET")):
            code_start = i; break
    
    docstring = "\n".join(lines[:doc_end])
    code_body  = "\n".join(lines[code_start:])
    return docstring + "\n" + std_header + "\n" + code_body

# 각 스크립트 소스 파일 경로
SOURCES = {
    "00_auto_objective.py":    "/home/claude/cho_figures/00_auto_objective_v2.py",
    "01_load_data.py":         "/mnt/user-data/outputs/01_load_data.py",
    "02_map_metabolites.py":   "/home/claude/cho_steps/02_map_metabolites.py",
    "03_run_fba.py":           "/home/claude/cho_steps/03_run_fba.py",
    "04_run_fva.py":           "/home/claude/cho_steps/04_run_fva.py",
    "05_ko_screening.py":      "/home/claude/cho_steps/05_ko_screening.py",
    "06_figures.py":           "/home/claude/cho_steps/06_figures.py",
    "fig_qp_timeseries.py":    "/mnt/user-data/outputs/fig_qp_timeseries.py",
    "fig_pathway_flux_map.py": "/mnt/user-data/outputs/fix_fig7_and_detailed_map.py",
    "multi_interval_fba.py":   "/mnt/user-data/outputs/multi_interval_fba.py",
}

for name, src in SOURCES.items():
    if os.path.exists(src):
        content = open(src, encoding="utf-8").read()
        # ROOT 경로를 새 프로젝트 루트로 업데이트
        content = content.replace(
            r'r"C:\\CHO_METABOLOMICS"', f'r"{ROOT}"')
        dst_path = os.path.join(ROOT, "scripts", "steps", name)
        with open(dst_path, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"  [OK] scripts/steps/{name}")
    else:
        print(f"  [--] {name} 소스 없음: {src}")

# ════════════════════════════════════════════════════
# 5. 루트 레벨 스크립트
# ════════════════════════════════════════════════════
print("\n[5/6] 루트 스크립트 생성")

write("run_pipeline.py", f'''
    """
    run_pipeline.py — 전체 파이프라인 실행
    사용: python run_pipeline.py --dataset practice_20aa
    """
    import subprocess, sys, os, time, argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", default="practice_20aa",
                        choices=["sowa2020","practice_20aa","own_experiment"])
    parser.add_argument("--steps",   default="0,1,2,3,4,5,6",
                        help="실행할 스텝 (예: 1,2,3)")
    parser.add_argument("--rate_days", default="7,10",
                        help="Rate 계산 구간 (기본: 7,10 = 전환기)")
    args = parser.parse_args()

    ROOT      = r"{ROOT}"
    STEPS_DIR = os.path.join(ROOT, "scripts", "steps")
    LOGS_DIR  = os.path.join(ROOT, "logs")
    os.makedirs(LOGS_DIR, exist_ok=True)

    STEP_MAP = {{
        0: ("00_auto_objective.py",    "Objective Auto-Screening"),
        1: ("01_load_data.py",         "Data Load + Rate Calculation"),
        2: ("02_map_metabolites.py",   "Metabolite Mapping"),
        3: ("03_run_fba.py",           "FBA Run"),
        4: ("04_run_fva.py",           "FVA Run"),
        5: ("05_ko_screening.py",      "Reaction KO Screening"),
        6: ("06_figures.py",           "Publication Figures"),
    }}
    EXTRA_STEPS = {{
        "qp":    ("fig_qp_timeseries.py",    "qP Timeseries + Trade-off"),
        "flux":  ("fig_pathway_flux_map.py", "Detailed Flux Map"),
        "multi": ("multi_interval_fba.py",   "Multi-interval FBA"),
    }}

    to_run = [int(x) for x in args.steps.split(",") if x.strip().isdigit()]
    env    = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"

    print("="*60)
    print(f"  CHO_METABOLOMICS Pipeline — {{args.dataset}}")
    print(f"  Rate interval: Day {{args.rate_days.replace(\",\",\"→\")}}")
    print("="*60)
    t0 = time.time()

    for num in to_run:
        script, desc = STEP_MAP.get(num, (None, None))
        if not script: continue
        path = os.path.join(STEPS_DIR, script)
        print(f"\\n[{{num}}/6] {{desc}}")
        cmd = [sys.executable, path, "--dataset", args.dataset]
        if num == 1: cmd += ["--rate_days", args.rate_days]
        log = os.path.join(LOGS_DIR, f"{{script.replace(\'.py\',\'\')}}_{{args.dataset}}.log")
        r = subprocess.run(cmd, cwd=ROOT, capture_output=True,
                           text=True, encoding="utf-8", errors="replace", env=env)
        with open(log, "w", encoding="utf-8") as lf:
            lf.write(r.stdout + "\\n" + r.stderr)
        kw = ["saved","obj=","status=","optimal","OK","ERROR","오류","완료","★","✔","✘"]
        for line in r.stdout.split("\\n"):
            if any(k in line for k in kw) and line.strip():
                print(f"  {{line.rstrip()}}")
        if r.returncode == 0:
            print(f"  ✔ 완료")
        else:
            print(f"  ✘ 오류 (code={{r.returncode}}) → {{log}}")
            tail = (r.stdout+r.stderr).strip().split("\\n")
            for l in tail[-5:]:
                if l.strip(): print(f"    {{l}}")
            if input("  계속? (y/n): ").lower() != "y": break

    # 추가 Figure
    extras = [k for k in EXTRA_STEPS if k in args.steps]
    for key in extras:
        script, desc = EXTRA_STEPS[key]
        print(f"\\n[extra] {{desc}}")
        path = os.path.join(STEPS_DIR, script)
        r = subprocess.run([sys.executable, path, "--dataset", args.dataset],
                           cwd=ROOT, capture_output=True, text=True,
                           encoding="utf-8", errors="replace", env=env)
        for line in r.stdout.split("\\n"):
            if "saved" in line and line.strip():
                print(f"  {{line.rstrip()}}")

    elapsed = time.time() - t0
    print(f"\\n{{'='*60}}")
    print(f"  완료! 소요: {{elapsed:.1f}}초")
    print(f"  결과: {{ROOT}}/results/{{args.dataset}}")

    # 생성 Figure 목록
    import glob
    pngs = sorted(glob.glob(os.path.join(ROOT,"results",args.dataset,"figures","*.png")))
    if pngs:
        print(f"\\n  생성 Figure ({{len(pngs)}}개):")
        for p in pngs:
            print(f"    {{os.path.basename(p)}}  ({{os.path.getsize(p)//1024}} KB)")
''')

write("setup_env.py", f'''
    """
    setup_env.py — 환경 초기 설정 (처음 한 번만 실행)
    사용: python setup_env.py
    """
    import subprocess, sys, os, zipfile, io

    ROOT = r"{ROOT}"
    print("="*55)
    print("  CHO_METABOLOMICS 환경 설정")
    print("="*55)

    # 패키지 설치
    PKGS = ["cobra","pandas","numpy","scipy","matplotlib",
            "seaborn","openpyxl","requests","scikit-learn"]
    missing = []
    for pkg in PKGS:
        imp = {{"scikit-learn":"sklearn"}}.get(pkg, pkg)
        try: __import__(imp); print(f"  ✔ {{pkg}}")
        except ImportError: missing.append(pkg); print(f"  ✘ {{pkg}}")
    if missing:
        subprocess.check_call([sys.executable,"-m","pip","install"]+missing)

    # iCHO3K 모델 다운로드
    MODEL_DIR = os.path.join(ROOT, "model", "iCHO3K")
    os.makedirs(MODEL_DIR, exist_ok=True)
    prod = next((f for f in os.listdir(MODEL_DIR) if "prod" in f and f.endswith(".json")), None)
    if prod:
        print(f"\\n  ✔ 모델 있음: {{prod}}")
    else:
        print("\\n  iCHO3K 다운로드 중...")
        try:
            import requests
            url = "https://github.com/LewisLabUCSD/iCHO3K/archive/refs/heads/main.zip"
            r = requests.get(url, timeout=120)
            r.raise_for_status()
            with zipfile.ZipFile(io.BytesIO(r.content)) as z:
                for member in z.namelist():
                    if "iCHO3K/Model/" in member and member.endswith(".json"):
                        fname = os.path.basename(member)
                        with z.open(member) as src, open(os.path.join(MODEL_DIR,fname),"wb") as dst:
                            dst.write(src.read())
                        print(f"  ✔ {{fname}}")
        except Exception as e:
            print(f"  ✘ {{e}}")
            print("  수동: github.com/LewisLabUCSD/iCHO3K → JSON → model/iCHO3K/")

    sys.path.insert(0, ROOT)
    from src.config import MODEL_PATH, SOWA_MMC1, AA20_EXCEL
    print("\\n  [파일 확인]")
    for k,v in [("모델",MODEL_PATH),("Sowa",SOWA_MMC1),("20AA",AA20_EXCEL)]:
        ok = "✔" if v and os.path.exists(v) else "✘"
        print(f"  {{ok}} {{k}}: {{v or '없음'}}")

    print("\\n  설정 완료!")
    print("  python run_pipeline.py --dataset practice_20aa")
''')

# .gitignore
write(".gitignore", '''
    # Python
    __pycache__/
    *.pyc
    .env

    # 대용량 (별도 공유)
    data/raw/
    model/
    data/processed/*.pkl

    # 로그
    logs/

    # OS
    .DS_Store
    Thumbs.db
    .ipynb_checkpoints/
''')

# README
write("README.md", f'''
    # CHO_METABOLOMICS — FBA/FVA Pipeline

    항체 고생산 CHO 세포주 대사 분석 + 개선 타겟 탐색

    ## 빠른 시작

    ```bash
    conda create -n cho python=3.10 -y && conda activate cho
    pip install -r requirements.txt
    python setup_env.py

    # 전체 파이프라인 (Day 7→10 전환기 기준)
    python run_pipeline.py --dataset practice_20aa --rate_days 7,10

    # 추가 Figure
    python scripts/steps/fig_qp_timeseries.py   --dataset practice_20aa
    python scripts/steps/fig_pathway_flux_map.py --dataset practice_20aa
    python scripts/steps/multi_interval_fba.py   --dataset practice_20aa
    ```

    ## 파이프라인 단계

    | Step | 파일 | 내용 | 근거 논문 |
    |------|------|------|-----------|
    | 0 | 00_auto_objective.py | Objective 자동 선택 | Feist & Palsson 2010 |
    | 1 | 01_load_data.py | 데이터 로딩 + Feed-corrected rate | Goudar et al. 2005 |
    | 2 | 02_map_metabolites.py | 대사물질 → Exchange ID 매핑 | — |
    | 3 | 03_run_fba.py | FBA (클론별) | Orth et al. 2010 |
    | 4 | 04_run_fva.py | FVA (High vs Low) | Mahadevan & Schilling 2003 |
    | 5 | 05_ko_screening.py | Reaction KO 스크리닝 | Lim et al. 2010 |
    | 6 | 06_figures.py | 논문 Figure 1~8 | — |
    | + | fig_qp_timeseries.py | qP + Growth-Production trade-off | Templeton 2013 |
    | + | fig_pathway_flux_map.py | 상세 Flux pathway map | King et al. 2015 |
    | + | multi_interval_fba.py | 구간별 FBA → 최적 구간 탐색 | Mulukutla 2012 |

    ## 핵심 수정사항 (이전 버전 대비)

    - **apply_bounds**: cobra bound 순서 충돌 방지 (lb>current_ub → ub 먼저)
    - **fix_one_sided**: lb/ub 단방향 보정 (Day 7→10 lactate 재소비 처리)
    - **Fig3**: all_rates 실측값 기반 Lac/Glc ratio (constraint 아님)
    - **Fig7 KO**: 3단계 반응 탐색 + 자동 크기 조정
    - **Objective**: IgG titer와 Pearson r 기반 자동 선택

    ## 새 데이터 추가

    ```bash
    data/raw/own_experiment/내실험.xlsx 배치 후:
    python run_pipeline.py --dataset own_experiment --rate_days 7,10
    ```

    ## 참고 논문

    | 분석 | 논문 |
    |------|------|
    | FBA | Orth et al. 2010 Nat Biotechnol 28:245 |
    | iCHO3K | Hernandez Bort et al. 2023 NPJ Syst Biol |
    | Rate | Goudar et al. 2005 Biotechnol Prog 21:1193 |
    | Lac/Glc | Zagari et al. 2013 New Biotechnology 30:238 |
    | 구간 | Mulukutla et al. 2012 Trends Biotechnol 30:616 |
    | LDH KO | Lim et al. 2010 Metab Eng 12:614 |
    | FVA | Mahadevan & Schilling 2003 Metab Eng 5:264 |
    | qP | Templeton et al. 2013 Biotechnol Bioeng 110:2508 |
    | Escher | King et al. 2015 PLOS Comput Biol 11:e1004321 |
    | Fouladiha | Fouladiha et al. 2020 Bioprocess Biosyst Eng 44:1 |
''')

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
    escher>=1.7
''')

# ════════════════════════════════════════════════════
# 6. 완료 리포트
# ════════════════════════════════════════════════════
print("\n[6/6] 완료 확인")

CHECKS = [
    ("src/config.py",                             "설정 파일"),
    ("src/fba_utils.py",                          "FBA 유틸 (apply_bounds 수정)"),
    ("scripts/steps/00_auto_objective.py",        "Objective 자동 선택"),
    ("scripts/steps/01_load_data.py",             "데이터 로딩"),
    ("scripts/steps/02_map_metabolites.py",       "대사물질 매핑"),
    ("scripts/steps/03_run_fba.py",               "FBA"),
    ("scripts/steps/04_run_fva.py",               "FVA"),
    ("scripts/steps/05_ko_screening.py",          "KO 스크리닝"),
    ("scripts/steps/06_figures.py",               "Figure 1~8"),
    ("scripts/steps/fig_qp_timeseries.py",        "qP + Trade-off"),
    ("scripts/steps/fig_pathway_flux_map.py",     "상세 Flux Map"),
    ("scripts/steps/multi_interval_fba.py",       "Multi-interval FBA"),
    ("run_pipeline.py",                           "전체 파이프라인 실행"),
    ("setup_env.py",                              "환경 설정"),
    ("requirements.txt",                          "패키지 목록"),
    ("README.md",                                 "문서"),
]

all_ok = True
for rel, desc in CHECKS:
    exists = os.path.exists(os.path.join(ROOT, rel))
    icon   = "✔" if exists else "✘"
    if not exists: all_ok = False
    print(f"  {icon}  {rel:<48s} {desc}")

print(f"\n{'='*60}")
print(f"  빌드 {'완료!' if all_ok else '일부 실패 — 위 ✘ 확인'}")
print(f"  ROOT: {ROOT}")
print(f"\n  실행 순서:")
print(f"  1. python setup_env.py              # 환경 설정 + 모델 다운로드")
print(f"  2. python run_pipeline.py --dataset practice_20aa --rate_days 7,10")
print(f"     (rate_days: 3,7=성장기 | 7,10=전환기★ | 10,14=생산기)")
print(f"  3. python scripts/steps/multi_interval_fba.py --dataset practice_20aa")
print(f"     → 최적 구간 자동 탐색")
print(f"{'='*60}")
