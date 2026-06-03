"""
fix_one_sided_constraints.py
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
lb/ub 둘 다 음수 (또는 양수) 인 constraint를 단방향으로 수정

문제:
  EX_glc_e:   lb=-0.4233, ub=-0.2822  → 양쪽 다 음수
  EX_lac_L_e: lb=-0.0826, ub=-0.0551  → 양쪽 다 음수
  → 모델이 정확히 이 범위 내에서만 flux 허용
  → 다른 반응 stoichiometry와 충돌 → infeasible

해결 (Ref: Varma & Palsson 1994 Appl Environ Microbiol):
  uptake (둘 다 음수):
    lb = 측정값 lb (최대 uptake 제한)
    ub = 0         (uptake 방향만 강제, 크기는 유연)

  secretion (둘 다 양수):
    lb = 0         (secretion 방향만 강제)
    ub = 측정값 ub (최대 secretion 제한)

  이렇게 하면 모델이 나머지 반응들과 균형을 맞출 수 있음

python scripts/diagnostics/fix_one_sided_constraints.py --dataset practice_20aa
"""
import sys, os, argparse, warnings, pickle, re
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
args   = parser.parse_args()
DATASET = args.dataset

print("="*65)
print(f"  fix_one_sided_constraints.py — {DATASET}")
print("="*65)
print()
print("  Ref: Varma & Palsson (1994) Appl Environ Microbiol 60:3724")
print("  Ref: Orth et al. (2010) Nat Biotechnol 28:245")
print("  Ref: Fouladiha et al. (2020) Bioprocess Biosyst Eng 44:1")

# ── 로딩 ──────────────────────────────────────────────
pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
with open(pkl_path, "rb") as f:
    data = pickle.load(f)

all_constraints = data.get("all_constraints", {})
igG_day14       = data.get("igG_day14", {})
high_clones     = data["high_clones"]
low_clones      = data["low_clones"]
clones          = data["clones"]

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

# ── 핵심 수정 함수 ─────────────────────────────────────
def fix_constraint(lb, ub):
    """
    lb/ub가 같은 부호일 때 단방향으로 수정
    Ref: Varma & Palsson 1994 — directional constraints
    """
    if lb < 0 and ub < 0:
        # 둘 다 음수 = uptake 방향
        # lb만 유지 (최대 uptake 제한), ub=0 (방향만 강제)
        return lb, 0.0

    elif lb > 0 and ub > 0:
        # 둘 다 양수 = secretion 방향
        # lb=0, ub만 유지 (최대 secretion 제한)
        return 0.0, ub

    else:
        # 이미 양방향 → 그대로 유지
        return lb, ub


# ── 클론별 constraints 수정 ──────────────────────────
print(f"\n[1] Constraint 수정 전/후 비교 (High clones)")
print(f"  {'Exchange ID':30s} {'원본 lb':>10s} {'원본 ub':>10s}  "
      f"{'수정 lb':>10s} {'수정 ub':>10s}")
print("  " + "─"*75)

# 첫 High clone으로 예시 출력
sample = high_clones[0] if high_clones else clones[0]
if sample in all_constraints:
    for rxn_id, (lb, ub) in sorted(all_constraints[sample].items()):
        new_lb, new_ub = fix_constraint(lb, ub)
        changed = "← 수정" if (new_lb != lb or new_ub != ub) else ""
        print(f"  {rxn_id:30s} {lb:10.4f} {ub:10.4f}  "
              f"{new_lb:10.4f} {new_ub:10.4f}  {changed}")

# 전체 클론에 적용
new_all_constraints = {}
for clone, cst in all_constraints.items():
    new_cst = {}
    for rxn_id, (lb, ub) in cst.items():
        new_lb, new_ub = fix_constraint(lb, ub)
        new_cst[rxn_id] = (new_lb, new_ub)
    new_all_constraints[clone] = new_cst

# ── Feasibility 확인 ──────────────────────────────────
print(f"\n[2] 수정 후 Feasibility 확인")
print(f"  {'Clone':12s} {'status':12s} {'obj':>10s}")
print("  " + "─"*40)

n_ok = 0
for clone in clones:
    if clone not in new_all_constraints:
        continue
    with model:
        apply_bounds(model, new_all_constraints[clone])
        model.objective = "biomass_cho_prod"
        sol = model.optimize()
        obj = sol.objective_value if sol.status == "optimal" else 0.0
        tag = "★" if clone in high_clones else \
              ("▼" if clone in low_clones else " ")
        print(f"  {tag} {clone:12s} {sol.status:12s} {obj:10.5f}")
        if sol.status == "optimal":
            n_ok += 1

print(f"\n  Feasible: {n_ok}/{len(clones)}")

# ── Objective 탐색 ─────────────────────────────────────
print(f"\n[3] Objective 후보 탐색 + IgG 상관관계")
print(f"  Ref: Feist & Palsson (2010) Nat Methods 7:482")
print()

# High 평균 constraints (수정된 것)
high_new = [new_all_constraints[c] for c in high_clones
            if c in new_all_constraints]
avg_new = {}
if high_new:
    all_ids = set()
    for c in high_new: all_ids.update(c.keys())
    avg_new = {rid: (np.mean([c[rid][0] for c in high_new if rid in c]),
                      np.mean([c[rid][1] for c in high_new if rid in c]))
               for rid in all_ids}

# IgG/biomass 관련 반응 탐색
CANDS = [r for r in model.reactions
         if any(k in f"{r.id} {r.name or ''}".lower()
                for k in ["igg","biomass_cho_prod","biomass_cho",
                           "antibody","dm_igg","igg_hc","igg_lc"])]

print(f"  {'Reaction':35s} {'obj':>9s} {'r':>7s} {'p':>7s}  sep")
print("  " + "─"*70)

corr_rows = []
for rxn in CANDS:
    # 작동 확인
    with model:
        apply_bounds(model, avg_new)
        model.objective = rxn.id
        sol = model.optimize()
        obj = sol.objective_value if sol.status == "optimal" else 0.0

    if obj < 1e-9:
        print(f"    {rxn.id:35s} {obj:9.5f}  (no flux)")
        continue

    # 클론별 FBA → IgG 상관관계
    clone_objs = {}
    for clone in clones:
        if clone not in new_all_constraints: continue
        with model:
            apply_bounds(model, new_all_constraints[clone])
            model.objective = rxn.id
            sol2 = model.optimize()
            clone_objs[clone] = sol2.objective_value \
                if sol2.status == "optimal" else 0.0

    pairs = [(igG_day14[c], clone_objs[c])
             for c in clones if c in igG_day14 and c in clone_objs]

    if len(pairs) >= 3:
        r, p = stats.pearsonr([x[0] for x in pairs], [x[1] for x in pairs])
    else:
        r, p = 0.0, 1.0

    h_mean = np.mean([clone_objs.get(c, 0) for c in high_clones])
    l_mean = np.mean([clone_objs.get(c, 0) for c in low_clones])
    sep    = h_mean - l_mean

    sig = "★★" if p < 0.01 else ("★" if p < 0.05 else "  ")
    print(f"  ★ {rxn.id:35s} {obj:9.5f} {r:+7.3f} {p:7.3f}{sig}  {sep:+.5f}")

    corr_rows.append({
        "id": rxn.id, "name": rxn.name or "",
        "obj_high_avg": obj,
        "pearson_r": r, "p_value": p, "separation": sep,
        "h_mean": h_mean, "l_mean": l_mean,
    })

# 최적 선택
if corr_rows:
    cdf = pd.DataFrame(corr_rows)
    pos = cdf[(cdf["pearson_r"] > 0) & (cdf["separation"] > 0)].sort_values(
          "pearson_r", ascending=False)
    chosen = pos.iloc[0]["id"] if len(pos) > 0 else \
             cdf.reindex(cdf["pearson_r"].abs().sort_values(
                         ascending=False).index).iloc[0]["id"]
    save_table(cdf, "objective_screening.csv", DATASET)
else:
    chosen = "biomass_cho_prod"

print(f"\n  ★ 선택된 Objective: {chosen}")

# ── pkl + config.py 업데이트 ──────────────────────────
data["all_constraints"]          = new_all_constraints
data["all_constraints_original"] = all_constraints
data["selected_objective"]       = chosen
with open(pkl_path, "wb") as f:
    pickle.dump(data, f)

cfg_path = os.path.join(ROOT, "src", "config.py")
with open(cfg_path, "r", encoding="utf-8") as f:
    cfg = f.read()
pat = re.compile(r'OBJ_RXN\s*=\s*"[^"]*"')
cfg = pat.sub(f'OBJ_RXN = "{chosen}"', cfg) if pat.search(cfg) \
      else cfg + f'\nOBJ_RXN = "{chosen}"\n'
with open(cfg_path, "w", encoding="utf-8") as f:
    f.write(cfg)

print(f"""
{'='*65}
  완료!

  수정 내용:
    lb/ub 둘 다 음수 → lb 유지, ub=0  (uptake 방향 강제)
    lb/ub 둘 다 양수 → lb=0, ub 유지  (secretion 방향 강제)

  Objective: {chosen}
  저장: data/processed/rates_{DATASET}.pkl (원본은 all_constraints_original)

  다음 실행:
  python scripts/steps/03_run_fba.py --dataset {DATASET}
  python scripts/steps/04_run_fva.py --dataset {DATASET}
{'='*65}
""")
