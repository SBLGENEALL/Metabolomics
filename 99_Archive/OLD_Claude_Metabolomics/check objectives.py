"""
check_objectives.py
iCHO3K prod 모델에서 항체 생산 관련 objective 확인

실행: python scripts/diagnostics/check_objectives.py
"""
import sys, os
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.config import *
from src.fba_utils import load_model, apply_bounds

import cobra

model = load_model(verbose=False)
all_rxn_ids = {r.id for r in model.reactions}

print("=" * 65)
print("  iCHO3K prod — 항체 관련 반응 목록")
print("=" * 65)

# IgG/antibody 관련 반응 탐색
keywords = ["igg", "antibody", "mab", "igG", "IgG",
            "heavy chain", "light chain", "igg_hc", "igg_lc",
            "DM_igg", "biomass_cho_prod", "biomass_cho"]

print("\n[항체/IgG 관련 반응]")
candidates = []
for r in model.reactions:
    r_str = f"{r.id} {r.name or ''}".lower()
    if any(k.lower() in r_str for k in keywords):
        candidates.append(r)
        print(f"  {r.id:35s}  {r.name or '':40s}")
        print(f"  {'':35s}  lb={r.lower_bound:.1f}, ub={r.upper_bound:.1f}")

# 각 후보를 objective로 설정해서 feasibility 테스트
print("\n[Objective 후보별 FBA 결과 — constraint 없이]")
print(f"  {'Reaction':35s} {'status':12s} {'obj':>10s}  {'의미'}")
print("  " + "─" * 75)

for rxn in candidates:
    try:
        with model:
            model.objective = rxn.id
            sol = model.optimize()
            obj = sol.objective_value if sol.status == "optimal" else 0
            print(f"  {rxn.id:35s} {sol.status:12s} {obj:10.5f}  {rxn.name or ''}")
    except Exception as e:
        print(f"  {rxn.id:35s} ERROR: {e}")

# Fouladiha 문헌 constraints 적용 후 재테스트
print("\n[Fouladiha High constraints 적용 후 FBA]")
print(f"  {'Reaction':35s} {'status':12s} {'obj':>10s}")
print("  " + "─" * 60)

for rxn in candidates:
    try:
        with model:
            apply_bounds(model, LITERATURE_HIGH)
            model.objective = rxn.id
            sol = model.optimize()
            obj = sol.objective_value if sol.status == "optimal" else 0
            star = " ★" if (sol.status == "optimal" and abs(obj) > 1e-6) else ""
            print(f"  {rxn.id:35s} {sol.status:12s} {obj:10.5f}{star}")
    except Exception as e:
        print(f"  {rxn.id:35s} ERROR: {e}")

print("""
결론:
  obj > 0  → 사용 가능 (항체 생산 proxy)
  obj = 0  → 연결 끊김 (사용 불가)

  권장 순서:
    1. igg_formation     → 직접적 IgG 합성 flux
    2. DM_igg_g          → IgG demand reaction
    3. igg_hc / igg_lc   → heavy/light chain 생산
    4. biomass_cho_prod  → 간접 proxy (지금 사용 중)
""")
