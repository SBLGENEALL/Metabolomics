"""
fba_utils.py — FBA/FVA 공통 함수 (최종본)

핵심 수정사항:
  apply_bounds: cobra bound 순서 충돌 방지
    Ref: Varma & Palsson (1994) Appl Environ Microbiol 60:3724
  fix_one_sided: lb/ub 둘 다 음수/양수 → 단방향 constraint 변환
    Ref: Orth et al. (2010) Nat Biotechnol 28:245
"""
import os, sys

def _find_root():
    cur = os.path.dirname(os.path.abspath(__file__))
    for _ in range(6):
        if os.path.exists(os.path.join(cur, "src", "config.py")):
            return cur
        cur = os.path.dirname(cur)
    return cur

ROOT = _find_root()
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from src.config import *
import cobra
import pandas as pd
import numpy as np
from cobra.flux_analysis import flux_variability_analysis


def load_model(verbose=True):
    """iCHO3K prod 모델 로딩"""
    if not MODEL_PATH or not os.path.exists(MODEL_PATH):
        raise FileNotFoundError(
            f"모델 없음: {MODEL_PATH}\n"
            "model/iCHO3K/ 에 iCHO3K_cho_prod_generic_unblocked.json 필요\n"
            "python setup_env.py 실행 시 자동 다운로드"
        )
    m = cobra.io.load_json_model(MODEL_PATH) if MODEL_PATH.endswith(".json") \
        else cobra.io.read_sbml_model(MODEL_PATH)
    m.objective = OBJ_RXN
    if verbose:
        print(f"  Model : {os.path.basename(MODEL_PATH)}")
        print(f"  Rxns  : {len(m.reactions)} | Genes: {len(m.genes)}")
        print(f"  Obj   : {OBJ_RXN}")
    return m


def apply_bounds(model, constraints: dict) -> int:
    """
    cobra bound 순서 충돌 방지 버전
    Ref: Varma & Palsson (1994) Appl Environ Microbiol 60:3724
      - new lb > current ub → ub 먼저 설정
      - new ub < current lb → lb 먼저 설정
    """
    rxn_ids = {r.id for r in model.reactions}
    n = 0
    for key, (lb, ub) in constraints.items():
        rxn_id = EXCHANGE_IDS.get(key, key)
        if rxn_id not in rxn_ids:
            continue
        rxn = model.reactions.get_by_id(rxn_id)
        lb_f, ub_f = float(lb), float(ub)
        if lb_f > ub_f:
            lb_f, ub_f = ub_f, lb_f
        if lb_f > rxn.upper_bound:
            rxn.upper_bound = ub_f
            rxn.lower_bound = lb_f
        elif ub_f < rxn.lower_bound:
            rxn.lower_bound = lb_f
            rxn.upper_bound = ub_f
        else:
            rxn.lower_bound = lb_f
            rxn.upper_bound = ub_f
        n += 1
    return n


def fix_one_sided(constraints: dict) -> dict:
    """
    lb/ub 둘 다 음수(uptake) → lb 유지, ub=0
    lb/ub 둘 다 양수(secretion) → lb=0, ub 유지
    Ref: Orth et al. (2010) Nat Biotechnol 28:245
      "Exchange constraints should be directional"
    """
    fixed = {}
    for rxn_id, (lb, ub) in constraints.items():
        if lb < 0 and ub < 0:
            fixed[rxn_id] = (lb, 0.0)
        elif lb > 0 and ub > 0:
            fixed[rxn_id] = (0.0, ub)
        else:
            fixed[rxn_id] = (lb, ub)
    return fixed


def rates_to_constraints(rates: dict, buffer: float = 0.20) -> dict:
    """
    exchange rate → FBA bounds 변환 + 단방향 보정
    Ref: Goudar et al. (2005) Biotechnol Prog 21:1193
    """
    cst = {}
    for key, q in rates.items():
        rxn_id = EXCHANGE_IDS.get(key, key)
        if abs(q) < 1e-9:
            continue
        if q < 0:
            lb, ub = q * (1 + buffer), q * (1 - buffer)
        else:
            lb, ub = q * (1 - buffer), q * (1 + buffer)
        cst[rxn_id] = (min(lb, ub), max(lb, ub))
    return fix_one_sided(cst)


def avg_constraints(clone_list: list, all_constraints: dict) -> dict:
    """여러 클론 constraints 평균"""
    cst_list = [all_constraints[c] for c in clone_list if c in all_constraints]
    if not cst_list:
        return {}
    all_ids = set()
    for c in cst_list:
        all_ids.update(c.keys())
    return {rid: (np.mean([c[rid][0] for c in cst_list if rid in c]),
                  np.mean([c[rid][1] for c in cst_list if rid in c]))
            for rid in all_ids}


def run_fba(model, constraints: dict, label: str = "") -> dict:
    """FBA 실행 (context manager → 모델 원본 보존)"""
    with model:
        n = apply_bounds(model, constraints)
        model.objective = OBJ_RXN
        sol = model.optimize()
        obj = sol.objective_value if sol.status == "optimal" else 0.0
        ex = {r.id: sol.fluxes[r.id] for r in model.exchanges
              if abs(sol.fluxes.get(r.id, 0)) > 1e-9} \
             if sol.status == "optimal" else {}
    if label:
        print(f"  {label:42s} n={n:2d}  {sol.status:10s}  obj={obj:.5f}")
    return {"label": label, "obj": obj, "status": sol.status,
            "ex_flux": ex, "n_constraints": n}


def run_fva(model, constraints: dict, label: str = "",
            exchange_ids: list = None) -> pd.DataFrame:
    """
    FVA 실행
    Ref: Mahadevan & Schilling (2003) Metab Eng 5:264
    """
    with model:
        apply_bounds(model, constraints)
        model.objective = OBJ_RXN
        rxn_ids_set = {r.id for r in model.reactions}
        rxn_list = None
        if exchange_ids:
            rxn_list = [model.reactions.get_by_id(r)
                        for r in exchange_ids if r in rxn_ids_set]
        try:
            fva = flux_variability_analysis(
                model, fraction_of_optimum=FVA_FRACTION,
                processes=FVA_PROCESSES,
                reaction_list=rxn_list or list(model.exchanges),
            )
            fva["mean"]  = (fva["minimum"] + fva["maximum"]) / 2
            fva["range"] = fva["maximum"] - fva["minimum"]
            if label:
                print(f"  FVA {label}: OK ({len(fva)} reactions)")
            return fva
        except Exception as e:
            print(f"  FVA {label}: {e}")
            return pd.DataFrame()


def save_figure(fig, name: str, dataset: str, subdir: str = "figures") -> str:
    """Figure 저장 300 DPI"""
    import matplotlib.pyplot as plt
    out = results_dir(dataset, subdir)
    path = os.path.join(out, name)
    fig.savefig(path, dpi=300, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"  [saved] results/{dataset}/{subdir}/{name}")
    return path


def save_table(df: pd.DataFrame, name: str, dataset: str) -> str:
    """DataFrame 저장"""
    out = results_dir(dataset, "tables")
    path = os.path.join(out, name)
    if name.endswith(".csv"):
        df.to_csv(path, index=False)
    elif name.endswith(".xlsx"):
        df.to_excel(path, index=False)
    print(f"  [saved] results/{dataset}/tables/{name}")
    return path


def ax_style(ax):
    """논문 스타일 axes"""
    ax.xaxis.grid(True, ls=":", color="0.88", zorder=0)
    ax.yaxis.grid(True, ls=":", color="0.88", zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
