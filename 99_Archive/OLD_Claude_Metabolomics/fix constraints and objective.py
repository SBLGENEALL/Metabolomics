"""
fix_constraints_and_objective.py
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
constraint 값 실제로 출력 + 모델 허용 범위와 비교 + 자동 수정

python scripts/diagnostics/fix_constraints_and_objective.py --dataset practice_20aa
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
from src.fba_utils import load_model, save_table

import pandas as pd
import numpy as np
from scipy import stats

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa")
args = parser.parse_args()
DATASET = args.dataset

print("="*65)
print(f"  fix_constraints_and_objective.py — {DATASET}")
print("="*65)

# ── 로딩 ──────────────────────────────────────────────
pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
with open(pkl_path, "rb") as f:
    data = pickle.load(f)

all_constraints = data.get("all_constraints", {})
all_rates       = data.get("all_rates", {})
igG_day14       = data.get("igG_day14", {})
high_clones     = data["high_clones"]
low_clones      = data["low_clones"]
clones          = data["clones"]
interval        = data.get("day_interval", (7,10))

model = load_model(verbose=False)
all_rxn_ids = {r.id for r in model.reactions}

def apply_bounds(model, constraints):
    n = 0
    for key, (lb, ub) in constraints.items():
        rxn_id = EXCHANGE_IDS.get(key, key)
        if rxn_id not in all_rxn_ids: continue
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

# ── 1. 실제 constraint 값 출력 ────────────────────────
print(f"\n[1] 실제 constraint 값 (High clones 평균)")
print(f"    Rate interval: Day {interval[0]} → Day {interval[1]}")
print()

# High clone 평균 constraints
high_csts = [all_constraints[c] for c in high_clones if c in all_constraints]
all_ids = set()
for c in high_csts: all_ids.update(c.keys())
avg_cst = {rid: (np.mean([c[rid][0] for c in high_csts if rid in c]),
                  np.mean([c[rid][1] for c in high_csts if rid in c]))
           for rid in all_ids}

print(f"  {'Exchange ID':30s} {'lb':>10s} {'ub':>10s}  "
      f"{'Model lb':>9s} {'Model ub':>9s}  {'상태'}")
print("  " + "─"*80)

problems = []
good = []
for rxn_id, (lb, ub) in sorted(avg_cst.items()):
    if rxn_id not in all_rxn_ids:
        continue
    rxn = model.reactions.get_by_id(rxn_id)
    m_lb, m_ub = rxn.lower_bound, rxn.upper_bound

    # 문제 유형 판단
    if lb > m_ub or ub < m_lb:
        status = "!! 범위 밖"
        problems.append(rxn_id)
    elif abs(lb) > 100 or abs(ub) > 100:
        status = "!! 비정상 크기"
        problems.append(rxn_id)
    else:
        status = "OK"
        good.append(rxn_id)

    icon = "!!" if rxn_id in problems else "  "
    print(f"  {icon} {rxn_id:30s} {lb:10.4f} {ub:10.4f}  "
          f"{m_lb:9.1f} {m_ub:9.1f}  {status}")

print(f"\n  문제 constraints: {len(problems)}개")
print(f"  정상 constraints: {len(good)}개")

# ── 2. 문제 원인 진단 ─────────────────────────────────
print(f"\n[2] 문제 원인 분석")

# 원본 rate 값 확인
if all_rates:
    sample_clone = high_clones[0] if high_clones else clones[0]
    if sample_clone in all_rates:
        print(f"\n  {sample_clone} 원본 rate 값 (Day {interval[0]}→{interval[1]}):")
        print(f"  {'대사물질':20s} {'rate(mmol/gDCW/h)':>18s}  {'Exchange ID'}")
        print("  " + "─"*60)
        rates = all_rates[sample_clone]
        for col, info in rates.items():
            q = info.get("rate_mmol_gDCWh", 0)
            ex_id = info.get("exchange_id", "")
            if abs(q) > 0.001:
                flag = " !!" if abs(q) > 10 else ""
                print(f"  {col:20s} {q:18.6f}  {ex_id}{flag}")

# ── 3. 올바른 생리적 범위로 재정규화 ─────────────────
print(f"\n[3] 생리적 범위 기반 Constraints 재생성")
print("""
  근거: Fouladiha et al. (2020) Bioprocess Biosyst Eng
        CHO fed-batch 정상 범위:
          Glucose uptake:  0.02 ~ 0.08 mmol/gDCW/h
          Lactate secr:    0.01 ~ 0.10 mmol/gDCW/h
          Glutamine uptake:0.005~ 0.025 mmol/gDCW/h
          NH4+ secretion:  0.003~ 0.020 mmol/gDCW/h
""")

# 생리적 허용 범위
PHYSIO_BOUNDS = {
    "EX_glc_e":   (-0.15, 0.0),    # glucose uptake only
    "EX_lac_L_e": (-0.05, 0.15),   # lactate mostly secreted
    "EX_gln_L_e": (-0.05, 0.01),   # glutamine mostly uptake
    "EX_glu_L_e": (-0.03, 0.03),
    "EX_nh4_e":   (-0.01, 0.05),   # NH4+ mostly secreted
    "EX_ala_L_e": (-0.02, 0.05),
    "EX_asp_L_e": (-0.02, 0.02),
    "EX_asn_L_e": (-0.02, 0.02),
    "EX_leu_L_e": (-0.03, 0.01),
    "EX_ile_L_e": (-0.03, 0.01),
    "EX_val_L_e": (-0.03, 0.01),
    "EX_ser_L_e": (-0.03, 0.01),
    "EX_pyr_e":   (-0.02, 0.05),
    "EX_cit_e":   (-0.02, 0.02),
}

def clamp_to_physio(lb, ub, rxn_id):
    """측정값을 생리적 허용 범위 내로 클램핑"""
    if rxn_id not in PHYSIO_BOUNDS:
        # 일반적인 exchange reaction 범위
        lb_c = max(lb, -0.2)
        ub_c = min(ub,  0.2)
    else:
        p_lb, p_ub = PHYSIO_BOUNDS[rxn_id]
        lb_c = max(lb, p_lb)
        ub_c = min(ub, p_ub)
    if lb_c > ub_c:
        # 클램핑 후에도 역전이면 중앙값 ±10% 사용
        mid = (lb + ub) / 2
        mid_c = max(PHYSIO_BOUNDS.get(rxn_id, (-0.05, 0.05))[0],
                    min(PHYSIO_BOUNDS.get(rxn_id, (-0.05, 0.05))[1], mid))
        lb_c = mid_c * 1.1 if mid_c < 0 else mid_c * 0.9
        ub_c = mid_c * 0.9 if mid_c < 0 else mid_c * 1.1
        if lb_c > ub_c: lb_c, ub_c = ub_c, lb_c
    return lb_c, ub_c

# 정규화 방법: Glucose rate 기준 스케일링
# Ref: Goudar et al. (2005) Biotechnol Prog 21:1193
def normalize_by_glucose(constraints, target_glc=-0.040):
    """Glucose uptake rate를 기준으로 전체 스케일 조정"""
    glc_id = "EX_glc_e"
    if glc_id not in constraints:
        return constraints
    glc_lb, glc_ub = constraints[glc_id]
    glc_val = (glc_lb + glc_ub) / 2  # 음수
    if abs(glc_val) < 1e-9:
        return constraints
    scale = target_glc / glc_val  # 두 음수의 비 → 양수
    if abs(scale - 1.0) < 0.1:
        return constraints  # 이미 정상 스케일
    print(f"  Glucose 기준 스케일 조정: scale={scale:.4f} "
          f"(현재 Glc rate={glc_val:.6f})")
    return {k: (lb*scale, ub*scale) for k,(lb,ub) in constraints.items()}

# 클론별 재생성
new_all_constraints = {}
print(f"\n  {'Clone':12s} {'원본 n':>7s} {'수정 n':>7s}  feasible")
print("  " + "─"*45)

for clone in clones:
    if clone not in all_constraints:
        continue
    raw_cst = all_constraints[clone]

    # Step A: Glucose 기준 정규화
    norm_cst = normalize_by_glucose(raw_cst)

    # Step B: 생리적 범위 클램핑
    clamped = {}
    for rxn_id, (lb, ub) in norm_cst.items():
        if rxn_id not in all_rxn_ids: continue
        lb_c, ub_c = clamp_to_physio(lb, ub, rxn_id)
        clamped[rxn_id] = (lb_c, ub_c)

    # feasibility 테스트
    with model:
        apply_bounds(model, clamped)
        model.objective = "biomass_cho_prod"
        sol = model.optimize()
        feasible = sol.status == "optimal"
        obj = sol.objective_value if feasible else 0.0

    new_all_constraints[clone] = clamped
    icon = "✔" if feasible else "✘"
    print(f"  {icon} {clone:12s} {len(raw_cst):7d} {len(clamped):7d}  "
          f"{'optimal' if feasible else 'infeasible'}  "
          f"{'obj='+str(round(obj,5)) if feasible else ''}")

# 전체 feasible 개수
n_feasible = sum(
    1 for clone, cst in new_all_constraints.items()
    if (lambda c: model.optimize().status=="optimal" if
        apply_bounds(model, c) >= 0 and
        setattr(model, 'objective', 'biomass_cho_prod') is None
        else False)(cst)
)

# ── 4. Objective 탐색 (수정된 constraints로) ──────────
print(f"\n[4] 수정된 constraints로 Objective 탐색")

# High 평균 (수정된 것)
high_new = [new_all_constraints[c] for c in high_clones
            if c in new_all_constraints]
avg_new = {}
if high_new:
    all_ids2 = set()
    for c in high_new: all_ids2.update(c.keys())
    avg_new = {rid: (np.mean([c[rid][0] for c in high_new if rid in c]),
                      np.mean([c[rid][1] for c in high_new if rid in c]))
               for rid in all_ids2}

CANDS = [r for r in model.reactions
         if any(k in f"{r.id} {r.name or ''}".lower()
                for k in ["igg","biomass_cho_prod","biomass_cho",
                           "antibody","dm_igg"])]

print(f"\n  {'Reaction':35s} {'obj':>10s}  feasible")
print("  " + "─"*55)

screening = []
for rxn in CANDS:
    with model:
        apply_bounds(model, avg_new)
        model.objective = rxn.id
        sol = model.optimize()
        obj = sol.objective_value if sol.status == "optimal" else 0.0
    feasible = obj > 1e-9
    icon = "★" if feasible else " "
    print(f"  {icon} {rxn.id:35s} {obj:10.5f}  "
          f"{'YES' if feasible else 'no'}")
    screening.append({"id": rxn.id, "obj": obj, "feasible": feasible})

screen_df = pd.DataFrame(screening)
feasible_cands = [s["id"] for s in screening if s["feasible"]]

# ── 5. IgG 상관관계 + 최적 objective 선택 ─────────────
corr_rows = []
if feasible_cands and igG_day14:
    print(f"\n[5] IgG 상관관계")
    print(f"  {'Reaction':35s} {'r':>7s} {'p':>7s}  sep")
    print("  " + "─"*60)
    for rxn_id in feasible_cands:
        clone_objs = {}
        for clone in clones:
            if clone not in new_all_constraints: continue
            with model:
                apply_bounds(model, new_all_constraints[clone])
                model.objective = rxn_id
                sol = model.optimize()
                clone_objs[clone] = sol.objective_value \
                    if sol.status == "optimal" else 0.0

        pairs = [(igG_day14[c], clone_objs[c])
                 for c in clones if c in igG_day14 and c in clone_objs]
        if len(pairs) >= 3:
            r, p = stats.pearsonr([x[0] for x in pairs],
                                   [x[1] for x in pairs])
        else:
            r, p = 0.0, 1.0

        h_mean = np.mean([clone_objs.get(c,0) for c in high_clones])
        l_mean = np.mean([clone_objs.get(c,0) for c in low_clones])
        sep    = h_mean - l_mean
        sig = "★★" if p<0.01 else ("★" if p<0.05 else "  ")
        print(f"  {rxn_id:35s} {r:+7.3f} {p:7.3f} {sig}  {sep:+.5f}")
        corr_rows.append({"id":rxn_id,"pearson_r":r,"p_value":p,
                           "separation":sep,"h_mean":h_mean,"l_mean":l_mean})

# 최적 선택
if corr_rows:
    cdf = pd.DataFrame(corr_rows)
    pos = cdf[(cdf["pearson_r"]>0)&(cdf["separation"]>0)].sort_values(
          "pearson_r", ascending=False)
    chosen = pos.iloc[0]["id"] if len(pos)>0 else \
             cdf.reindex(cdf["pearson_r"].abs().sort_values(
                         ascending=False).index).iloc[0]["id"]
elif feasible_cands:
    chosen = feasible_cands[0]
else:
    chosen = "biomass_cho_prod"

print(f"\n  ★ 선택된 Objective: {chosen}")

# ── 6. pkl + config.py 업데이트 ───────────────────────
import re as _re
data["all_constraints"]         = new_all_constraints  # 수정된 버전으로 교체
data["all_constraints_original"]= all_constraints      # 원본 보존
data["selected_objective"]      = chosen
with open(pkl_path, "wb") as f:
    pickle.dump(data, f)

cfg_path = os.path.join(ROOT, "src", "config.py")
with open(cfg_path, "r", encoding="utf-8") as f:
    cfg = f.read()
pat = _re.compile(r'OBJ_RXN\s*=\s*"[^"]*"')
if pat.search(cfg):
    cfg = pat.sub(f'OBJ_RXN = "{chosen}"', cfg)
else:
    cfg += f'\nOBJ_RXN = "{chosen}"\n'
with open(cfg_path, "w", encoding="utf-8") as f:
    f.write(cfg)

if corr_rows:
    save_table(pd.DataFrame(corr_rows), "objective_screening.csv", DATASET)
save_table(pd.DataFrame(screening),     "objective_candidates.csv", DATASET)

print(f"""
{'='*65}
  완료!
  - Constraints: 생리적 범위로 재정규화 + 클램핑
  - Objective  : {chosen}
  - pkl 업데이트: all_constraints (수정됨), 원본은 all_constraints_original 보존

  다음 실행:
  python scripts/steps/03_run_fba.py --dataset {DATASET}
  python scripts/steps/04_run_fva.py --dataset {DATASET}
{'='*65}
""")
