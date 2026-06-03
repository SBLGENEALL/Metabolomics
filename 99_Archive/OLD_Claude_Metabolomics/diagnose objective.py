"""
diagnose_objective.py
왜 feasible objective가 없는지 단계별 진단
python scripts/diagnostics/diagnose_objective.py --dataset practice_20aa
"""
import sys, os, argparse, warnings, pickle
warnings.filterwarnings("ignore")

def _find_root():
    cur = os.path.dirname(os.path.abspath(__file__))
    for _ in range(5):
        if os.path.exists(os.path.join(cur, "src", "config.py")):
            return cur
        cur = os.path.dirname(cur)
    return cur

ROOT = _find_root()
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from src.config import *
from src.fba_utils import load_model

import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa")
args = parser.parse_args()
DATASET = args.dataset

print("=" * 65)
print(f"  Objective Diagnosis — {DATASET}")
print("=" * 65)

# ── 데이터 로딩 ────────────────────────────────────────
pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
with open(pkl_path, "rb") as f:
    data = pickle.load(f)

all_constraints = data.get("all_constraints", {})
high_clones     = data["high_clones"]
clones          = data["clones"]

model = load_model(verbose=False)
all_rxn_ids = {r.id for r in model.reactions}

def apply_bounds(model, constraints):
    n = 0
    for key, (lb, ub) in constraints.items():
        rxn_id = EXCHANGE_IDS.get(key, key)
        if rxn_id not in all_rxn_ids:
            continue
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

def avg_cst(clone_list):
    cst_list = [all_constraints[c] for c in clone_list if c in all_constraints]
    if not cst_list: return {}
    all_ids = set()
    for c in cst_list: all_ids.update(c.keys())
    return {rid: (np.mean([c[rid][0] for c in cst_list if rid in c]),
                  np.mean([c[rid][1] for c in cst_list if rid in c]))
            for rid in all_ids}

HIGH_CST = avg_cst(high_clones)

# ── 진단 1: constraints 자체가 feasible한지 ────────────
print("\n[진단 1] Constraint 적용 후 기본 FBA feasibility")
for obj_id in ["biomass_cho_prod", "biomass_cho", "biomass_cho_s"]:
    if obj_id not in all_rxn_ids:
        print(f"  {obj_id:35s} — 모델에 없음")
        continue
    with model:
        n = apply_bounds(model, HIGH_CST)
        model.objective = obj_id
        sol = model.optimize()
        print(f"  {obj_id:35s}  n={n}  {sol.status:12s}  "
              f"obj={sol.objective_value:.6f}" if sol.status=="optimal"
              else f"  {obj_id:35s}  n={n}  {sol.status}")

# ── 진단 2: constraint 없이 FBA ────────────────────────
print("\n[진단 2] Constraint 없이 FBA (모델 자체 feasibility)")
for obj_id in ["biomass_cho_prod", "igg_formation", "DM_igg_g",
               "igg_hc", "igg_lc"]:
    if obj_id not in all_rxn_ids:
        print(f"  {obj_id:35s} — 모델에 없음")
        continue
    with model:
        model.objective = obj_id
        sol = model.optimize()
        obj = sol.objective_value if sol.status == "optimal" else 0
        print(f"  {obj_id:35s}  {sol.status:12s}  obj={obj:.6f}")

# ── 진단 3: 어떤 constraint가 문제인지 ─────────────────
print("\n[진단 3] 문제 constraint 찾기 (하나씩 제거)")
if HIGH_CST:
    with model:
        apply_bounds(model, HIGH_CST)
        model.objective = "biomass_cho_prod"
        sol_full = model.optimize()
        print(f"  전체 {len(HIGH_CST)}개 적용: {sol_full.status}")

    # constraint 하나씩 제거하며 infeasible 원인 탐색
    if sol_full.status != "optimal":
        print("\n  문제 constraint 탐색 중...")
        for skip_id in list(HIGH_CST.keys()):
            cst_without = {k: v for k, v in HIGH_CST.items() if k != skip_id}
            with model:
                apply_bounds(model, cst_without)
                model.objective = "biomass_cho_prod"
                sol_t = model.optimize()
                if sol_t.status == "optimal":
                    rxn_name = model.reactions.get_by_id(skip_id).name \
                               if skip_id in all_rxn_ids else skip_id
                    print(f"  !! 이 constraint 제거 시 optimal: {skip_id} ({rxn_name})")
                    lb, ub = HIGH_CST[skip_id]
                    print(f"     lb={lb:.5f}, ub={ub:.5f}")

# ── 진단 4: constraints 적용된 exchange rate 확인 ──────
print("\n[진단 4] 실제 적용된 High constraints 값")
print(f"  {'Exchange ID':30s} {'lb':>10s} {'ub':>10s}  {'방향'}")
print("  " + "─" * 58)
for rxn_id, (lb, ub) in sorted(HIGH_CST.items()):
    direction = "uptake" if lb < 0 and ub <= 0 else \
                "secretion" if lb >= 0 else "both"
    met_name = next((name for name, eid in EXCHANGE_IDS.items()
                     if eid == rxn_id), rxn_id)
    print(f"  {rxn_id:30s} {lb:10.5f} {ub:10.5f}  {direction} ({met_name})")

# ── 진단 5: rate interval 권고 ─────────────────────────
print("\n[진단 5] 현재 Rate interval")
interval = data.get("day_interval", "unknown")
print(f"  현재: Day {interval[0]} → Day {interval[1]}")
if interval[0] <= 5:
    print("  ⚠ 성장기 구간 — FBA infeasible 가능성 높음")
    print("  권고: --rate_days 7,10 또는 --rate_days 7,12 로 재실행")
    print()
    print("  재실행 명령:")
    print(f"  python scripts/steps/01_load_data.py "
          f"--dataset {DATASET} --rate_days 7,10")
    print(f"  python scripts/steps/02_map_metabolites.py --dataset {DATASET}")
    print(f"  python scripts/steps/00_auto_objective.py --dataset {DATASET}")
