"""
00_auto_objective.py  v2
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Objective Auto-Screening + Feasibility Auto-Recovery

infeasible 시 자동 복구 전략:
  1단계: buffer 확대 (20% → 50% → 100%)
  2단계: 핵심 constraints만 유지 (Glc, Lac, Gln, NH4+)
  3단계: 모든 constraints 제거 후 모델만으로 테스트
  → 작동하는 최소 feasible set 탐색

근거:
  Orth et al. (2010) Nat Biotechnol 28:245
    "FBA requires consistent constraints"
  Feist & Palsson (2010) Nat Methods 7:482
    "Objective selection must account for model feasibility"
  Reed & Palsson (2004) Genome Res 14:1797
    "Constraint relaxation for feasibility recovery"

사용:
  python scripts/steps/00_auto_objective.py --dataset practice_20aa
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
parser.add_argument("--dataset", default="practice_20aa",
                    choices=["sowa2020", "practice_20aa", "own_experiment"])
parser.add_argument("--force", default=None)
args = parser.parse_args()
DATASET = args.dataset

print("=" * 65)
print(f"  00_auto_objective.py v2 — {DATASET}")
print("=" * 65)
print("  Ref: Orth et al. (2010) Nat Biotechnol 28:245")
print("  Ref: Feist & Palsson (2010) Nat Methods 7:482")
print("  Ref: Reed & Palsson (2004) Genome Res 14:1797")

# ── 로딩 ──────────────────────────────────────────────
pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
if not os.path.exists(pkl_path):
    print(f"\n  !! {pkl_path} 없음 — 01 & 02 먼저 실행")
    sys.exit(1)
with open(pkl_path, "rb") as f:
    data = pickle.load(f)

all_constraints = data.get("all_constraints", {})
igG_day14       = data.get("igG_day14", {})
high_clones     = data["high_clones"]
low_clones      = data["low_clones"]
clones          = data["clones"]

model = load_model(verbose=True)
all_rxn_ids = {r.id for r in model.reactions}

# ── Safe apply_bounds ──────────────────────────────────
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

def is_feasible(cst, obj="biomass_cho_prod"):
    """주어진 constraint + objective가 feasible한지 확인"""
    if obj not in all_rxn_ids:
        return False, 0.0
    with model:
        apply_bounds(model, cst)
        model.objective = obj
        sol = model.optimize()
        obj_val = sol.objective_value if sol.status == "optimal" else 0.0
        return sol.status == "optimal", obj_val

def widen_bounds(cst, buffer):
    """
    constraint bounds를 buffer 비율로 확장
    Ref: Reed & Palsson 2004 — constraint relaxation
    """
    widened = {}
    for rxn_id, (lb, ub) in cst.items():
        if lb < 0:   # uptake: lb를 더 음수로
            new_lb = lb * (1 + buffer)
            new_ub = ub * max(0, (1 - buffer))
        else:        # secretion: ub를 더 크게
            new_lb = lb * max(0, (1 - buffer))
            new_ub = ub * (1 + buffer)
        widened[rxn_id] = (min(new_lb, new_ub), max(new_lb, new_ub))
    return widened

# ── STEP 1: Feasibility 자동 복구 ─────────────────────
print(f"\n{'─'*65}")
print("  [STEP 1] Feasibility 자동 복구")
print(f"{'─'*65}")

CORE_IDS = [
    EXCHANGE_IDS.get("Glucose",   "EX_glc_e"),
    EXCHANGE_IDS.get("Lactate",   "EX_lac_L_e"),
    EXCHANGE_IDS.get("Glutamine", "EX_gln_L_e"),
    EXCHANGE_IDS.get("NH4+",      "EX_nh4_e"),
]

RECOVERY_STRATEGIES = [
    ("원본 (±20%)",           HIGH_CST,                                  False),
    ("완화 ±50%",             widen_bounds(HIGH_CST, 0.50),              False),
    ("완화 ±100%",            widen_bounds(HIGH_CST, 1.00),              False),
    ("완화 ±200%",            widen_bounds(HIGH_CST, 2.00),              False),
    ("핵심 4개만 (±50%)",
     {k: v for k, v in widen_bounds(HIGH_CST, 0.50).items()
      if k in CORE_IDS},                                                  False),
    ("핵심 4개만 (±100%)",
     {k: v for k, v in widen_bounds(HIGH_CST, 1.00).items()
      if k in CORE_IDS},                                                  False),
    ("Glucose + Lactate만",
     {k: v for k, v in widen_bounds(HIGH_CST, 0.50).items()
      if k in CORE_IDS[:2]},                                              False),
    ("Constraints 없음",      {},                                          True),
]

working_cst   = None
working_label = None
fallback_no_cst = False

for label, cst, no_cst_flag in RECOVERY_STRATEGIES:
    feasible, obj_val = is_feasible(cst)
    status_str = f"optimal  obj={obj_val:.5f}" if feasible else "infeasible"
    icon = "✔" if feasible else "✘"
    print(f"  {icon}  {label:30s}  {status_str}")

    if feasible and working_cst is None:
        working_cst   = cst
        working_label = label
        fallback_no_cst = no_cst_flag

if working_cst is None:
    print("\n  !! 모든 전략 실패 — 모델 자체 문제")
    print("  → iCHO3K_cho_generic_unblocked.json 사용 권고")
    sys.exit(1)

print(f"\n  ★ 사용할 Constraints: '{working_label}'")
if fallback_no_cst:
    print("  ⚠ Constraint 없이 FBA — objective 탐색만 가능")
    print("    Rate interval 재검토 권장 (multi_interval_fba.py 실행)")

# ── STEP 2: 문제 constraint 진단 리포트 ───────────────
if working_label != "원본 (±20%)":
    print(f"\n{'─'*65}")
    print("  [STEP 2] 문제 Constraint 진단")
    print(f"{'─'*65}")
    print("  원본 infeasible 원인 탐색:")
    found_culprit = False
    for skip_id in list(HIGH_CST.keys()):
        cst_without = {k: v for k, v in HIGH_CST.items() if k != skip_id}
        feasible_t, _ = is_feasible(cst_without)
        if feasible_t:
            rxn_name = model.reactions.get_by_id(skip_id).name \
                       if skip_id in all_rxn_ids else "unknown"
            lb, ub = HIGH_CST[skip_id]
            # 모델 기본값과 비교
            if skip_id in all_rxn_ids:
                rxn = model.reactions.get_by_id(skip_id)
                print(f"\n  !! 원인: {skip_id} ({rxn_name})")
                print(f"     측정값:  lb={lb:.5f}, ub={ub:.5f}")
                print(f"     모델기본: lb={rxn.lower_bound:.1f}, ub={rxn.upper_bound:.1f}")
                if ub < rxn.lower_bound:
                    print(f"     → 측정 ub({ub:.5f}) < 모델 lb({rxn.lower_bound:.1f})")
                    print(f"        측정값이 모델 허용 범위 밖")
            found_culprit = True
            break
    if not found_culprit:
        print("  복수 constraint의 조합이 원인 (단일 제거로 해결 안 됨)")
        print("  → buffer 완화가 올바른 해결책")

# ── STEP 3: Objective 후보 탐색 ───────────────────────
print(f"\n{'─'*65}")
print("  [STEP 3] Objective 후보 탐색")
print(f"{'─'*65}")

ANTIBODY_KW = ["igg","antibody","mab","heavy chain","light chain",
               "DM_igg","BiGGEx","immunoglobulin"]
BIOMASS_KW  = ["biomass_cho_prod","biomass_cho"]

candidates = []
for r in model.reactions:
    r_str = f"{r.id} {r.name or ''}".lower()
    if any(k.lower() in r_str for k in ANTIBODY_KW + BIOMASS_KW):
        candidates.append(r)

print(f"  {'Reaction':35s} {'no_cst':>9s} {'w_cst':>9s}  feasible")
print("  " + "─" * 60)

screening = []
for rxn in candidates:
    # constraint 없이
    with model:
        model.objective = rxn.id
        sol_nc = model.optimize()
        obj_nc = sol_nc.objective_value if sol_nc.status == "optimal" else 0.0

    # 작동하는 constraint로
    with model:
        apply_bounds(model, working_cst)
        model.objective = rxn.id
        sol_wc = model.optimize()
        obj_wc = sol_wc.objective_value if sol_wc.status == "optimal" else 0.0

    feasible = obj_wc > 1e-9 or obj_nc > 1e-9
    icon = "★" if feasible else " "
    print(f"  {icon} {rxn.id:35s} {obj_nc:9.5f} {obj_wc:9.5f}  "
          f"{'YES' if feasible else 'no'}")
    screening.append({
        "id":          rxn.id,
        "name":        rxn.name or "",
        "obj_no_cst":  obj_nc,
        "obj_w_cst":   obj_wc,
        "feasible":    feasible,
    })

# ── STEP 4: IgG 상관관계 계산 ─────────────────────────
print(f"\n{'─'*65}")
print("  [STEP 4] IgG titer 상관관계 계산")
print(f"{'─'*65}")

feasible_cands = [s for s in screening if s["feasible"]]

if not feasible_cands:
    print("  !! 작동하는 objective 없음")
    chosen = "biomass_cho_prod"
    print(f"  → 기본값 사용: {chosen}")
else:
    corr_rows = []
    for cand in feasible_cands:
        rxn_id = cand["id"]
        clone_objs = {}
        for clone in clones:
            if clone not in all_constraints:
                continue
            # 클론별: 작동하는 전략과 동일한 buffer 적용
            if working_label == "원본 (±20%)":
                clone_cst = all_constraints[clone]
            elif "핵심" in working_label:
                buf = 0.50 if "50%" in working_label else 1.00
                clone_cst = {k: v for k, v in
                             widen_bounds(all_constraints[clone], buf).items()
                             if k in CORE_IDS}
            elif "Constraints 없음" in working_label:
                clone_cst = {}
            else:
                buf_map = {"±50%": 0.50, "±100%": 1.00, "±200%": 2.00}
                buf = next((v for k, v in buf_map.items()
                            if k in working_label), 0.50)
                clone_cst = widen_bounds(all_constraints[clone], buf)

            with model:
                apply_bounds(model, clone_cst)
                model.objective = rxn_id
                sol = model.optimize()
                clone_objs[clone] = sol.objective_value \
                    if sol.status == "optimal" else 0.0

        pairs = [(igG_day14[c], clone_objs[c])
                 for c in clones if c in igG_day14 and c in clone_objs]
        if len(pairs) >= 3:
            r, p_val = stats.pearsonr([p[0] for p in pairs],
                                       [p[1] for p in pairs])
        else:
            r, p_val = 0.0, 1.0

        high_mean = np.mean([clone_objs.get(c, 0) for c in high_clones])
        low_mean  = np.mean([clone_objs.get(c, 0) for c in low_clones])
        sep = high_mean - low_mean

        sig = "★★" if p_val < 0.01 else ("★" if p_val < 0.05 else "  ")
        direction = "↑High" if sep > 0 else "↓High"
        print(f"  {rxn_id:35s}  r={r:+.3f} {sig}  "
              f"sep={sep:+.5f}({direction})")

        corr_rows.append({
            "id": rxn_id, "name": cand["name"],
            "pearson_r": r, "p_value": p_val,
            "high_mean": high_mean, "low_mean": low_mean,
            "separation": sep,
        })

    corr_df = pd.DataFrame(corr_rows)

    # 선택: 정상관 + High > Low 우선
    positive = corr_df[(corr_df["pearson_r"] > 0) &
                       (corr_df["separation"] > 0)].sort_values(
                       "pearson_r", ascending=False)

    if len(positive) > 0:
        best = positive.iloc[0]
        chosen = best["id"]
        print(f"\n  ★ 선택: {chosen}  r={best['pearson_r']:.3f}  "
              f"p={best['p_value']:.3f}")
    else:
        # 절대값 r 최대
        best = corr_df.reindex(
            corr_df["pearson_r"].abs().sort_values(ascending=False).index
        ).iloc[0]
        chosen = best["id"]
        print(f"\n  ~ 정상관 없음 → 절대값 r 최대: {chosen}  "
              f"r={best['pearson_r']:.3f}")
        if best["pearson_r"] < 0:
            print("  ⚠ 역상관 — Rate interval 변경 권고")
            print("    python scripts/steps/multi_interval_fba.py "
                  f"--dataset {DATASET}")

    save_table(corr_df, "objective_screening.csv", DATASET)

# ── STEP 5: config.py 업데이트 ────────────────────────
print(f"\n{'─'*65}")
print("  [STEP 5] config.py 업데이트")
print(f"{'─'*65}")

cfg_path = os.path.join(ROOT, "src", "config.py")
with open(cfg_path, "r", encoding="utf-8") as f:
    cfg = f.read()

old_pat = re.compile(r'OBJ_RXN\s*=\s*"[^"]*"')
new_line = f'OBJ_RXN = "{chosen}"'
if old_pat.search(cfg):
    old = old_pat.search(cfg).group()
    cfg = old_pat.sub(new_line, cfg)
    print(f"  {old}  →  {new_line}")
else:
    cfg += f"\n{new_line}\n"
    print(f"  추가: {new_line}")

with open(cfg_path, "w", encoding="utf-8") as f:
    f.write(cfg)

# pkl 저장
data["selected_objective"]      = chosen
data["working_constraint_label"] = working_label
data["working_cst_high"]        = working_cst
# 클론별 완화된 constraints도 저장
if working_label != "원본 (±20%)":
    buf_map = {"±50%": 0.50, "±100%": 1.00, "±200%": 2.00}
    buf = next((v for k, v in buf_map.items() if k in working_label), 0.50)
    if "핵심" in working_label:
        data["all_constraints_relaxed"] = {
            clone: {k: v for k, v in widen_bounds(cst, buf).items()
                    if k in CORE_IDS}
            for clone, cst in all_constraints.items()
        }
    elif "없음" not in working_label:
        data["all_constraints_relaxed"] = {
            clone: widen_bounds(cst, buf)
            for clone, cst in all_constraints.items()
        }
    print(f"\n  완화된 constraints도 pkl에 저장 (키: all_constraints_relaxed)")

with open(pkl_path, "wb") as f:
    pickle.dump(data, f)

print(f"\n{'='*65}")
print(f"  완료!")
print(f"{'='*65}")
print(f"""
  선택된 Objective : {chosen}
  사용된 Constraints: {working_label}

  다음 실행:
  python scripts/steps/03_run_fba.py --dataset {DATASET}
  python scripts/steps/04_run_fva.py --dataset {DATASET}
""")
