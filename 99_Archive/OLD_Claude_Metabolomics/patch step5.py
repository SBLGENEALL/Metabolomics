"""
patch_step5.py
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Step 5 baseline infeasible 자동 수정
  1. Step 3 결과에서 실제로 optimal인 클론만 추출
  2. 그 클론들의 constraints로 baseline 재설정
  3. 수정된 05_ko_screening.py 저장

python patch_step5.py --dataset practice_20aa
"""
import sys, os, argparse, pickle, warnings, re
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
from src.fba_utils import load_model, apply_bounds, avg_constraints

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa")
args = parser.parse_args()
DATASET = args.dataset

print("="*60)
print(f"  patch_step5.py — {DATASET}")
print("="*60)

# ── pkl 로딩 ──────────────────────────────────────────
pkl = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
with open(pkl, "rb") as f:
    data = pickle.load(f)

all_constraints = data.get("all_constraints", {})
fba_results     = data.get("fba_results", {})
high_clones     = data["high_clones"]
low_clones      = data["low_clones"]
clones          = data["clones"]
OBJ = data.get("selected_objective", OBJ_RXN)

print(f"\n  Step 3 결과:")
print(f"  {'Clone':12s} {'status':12s} {'obj':>10s}  group")
print("  " + "─"*45)
for clone in clones:
    r = fba_results.get(clone, {})
    tag = "★High" if clone in high_clones else \
          ("▼Low"  if clone in low_clones  else "  mid")
    print(f"  {clone:12s} {r.get('status','N/A'):12s} "
          f"{r.get('obj',0):10.5f}  {tag}")

# optimal인 클론 선별
optimal_clones = [c for c in clones
                  if fba_results.get(c,{}).get("status")=="optimal"
                  and abs(fba_results.get(c,{}).get("obj",0)) > 1e-9]
print(f"\n  Optimal 클론: {optimal_clones} ({len(optimal_clones)}/{len(clones)})")

# ── 모델 로딩 ─────────────────────────────────────────
model = load_model(verbose=True)
all_rxn_ids = {r.id for r in model.reactions}

def widen(cst, buf):
    out = {}
    for k,(lb,ub) in cst.items():
        if lb < 0: out[k] = (lb*(1+buf), 0.0)
        elif ub > 0: out[k] = (0.0, ub*(1+buf))
        else: out[k] = (lb,ub)
    return out

def test_feasible(cst, obj=OBJ):
    if obj not in all_rxn_ids: return False, 0.0
    with model:
        apply_bounds(model, cst)
        model.objective = obj
        sol = model.optimize()
        return sol.status=="optimal", (sol.objective_value if sol.status=="optimal" else 0.0)

# ── 여러 전략으로 feasible baseline 탐색 ──────────────
print(f"\n  Baseline 탐색:")

# 전략 1: optimal 클론만으로 High avg
if optimal_clones:
    opt_high = [c for c in optimal_clones if c in high_clones]
    opt_any  = optimal_clones
    for label, cl in [
        ("High optimal 클론만", opt_high or opt_any[:len(opt_any)//2+1]),
        ("Optimal 전체 avg",    opt_any),
        ("모든 클론 avg",       clones),
    ]:
        cl2 = [c for c in cl if c in all_constraints]
        if not cl2: continue
        cst = avg_constraints(cl2, all_constraints)
        ok, obj_val = test_feasible(cst)
        print(f"  {'✔' if ok else '✘'} {label:30s} obj={obj_val:.5f}  ({ok})")
        if ok:
            best_cst   = cst
            best_label = label
            best_obj   = obj_val
            break
    else:
        # 전략 2: buffer 확대
        for buf, label in [(0.5,"±50%"),(1.0,"±100%"),(2.0,"±200%")]:
            cl2 = [c for c in opt_any if c in all_constraints]
            cst = widen(avg_constraints(cl2, all_constraints), buf)
            ok, obj_val = test_feasible(cst)
            print(f"  {'✔' if ok else '✘'} buffer {label:8s}  obj={obj_val:.5f}")
            if ok:
                best_cst   = cst
                best_label = f"buffer {label}"
                best_obj   = obj_val
                break
        else:
            # 전략 3: 문헌값 사용
            print("  ~ 문헌값(Fouladiha) 사용")
            best_cst   = LITERATURE_HIGH
            best_label = "Fouladiha 2020 문헌값"
            ok, best_obj = test_feasible(LITERATURE_HIGH)
            print(f"  {'✔' if ok else '✘'} Fouladiha High  obj={best_obj:.5f}")
else:
    print("  !! Optimal 클론 없음 → 문헌값 사용")
    best_cst   = LITERATURE_HIGH
    best_label = "Fouladiha 2020 문헌값"
    _, best_obj = test_feasible(LITERATURE_HIGH)

print(f"\n  ★ 선택된 baseline: {best_label}  obj={best_obj:.5f}")

# pkl에 저장
data["ko_baseline_cst"]   = best_cst
data["ko_baseline_label"] = best_label
data["ko_baseline_obj"]   = best_obj
with open(pkl, "wb") as f:
    pickle.dump(data, f)
print(f"  pkl 업데이트 완료 (ko_baseline_cst)")

# ── 05_ko_screening.py 패치 ───────────────────────────
f05 = os.path.join(ROOT, "scripts", "steps", "05_ko_screening.py")
with open(f05, "r", encoding="utf-8") as f:
    content = f.read()

# HIGH_CST 계산 부분을 pkl의 ko_baseline_cst 사용으로 교체
old_block = """HIGH_CST = avg_cst(high_clones)
print(f"  High avg constraints: {len(HIGH_CST)}개")"""

new_block = """# ko_baseline_cst: patch_step5.py가 탐색한 feasible constraint
if "ko_baseline_cst" in data and data["ko_baseline_cst"]:
    HIGH_CST = data["ko_baseline_cst"]
    print(f"  Baseline: {data.get('ko_baseline_label','custom')}  "
          f"({len(HIGH_CST)}개 constraints)")
else:
    HIGH_CST = avg_cst(high_clones)
    print(f"  High avg constraints: {len(HIGH_CST)}개")"""

if old_block in content:
    content = content.replace(old_block, new_block)
    with open(f05, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"\n  05_ko_screening.py 패치 완료")
else:
    print(f"\n  ⚠ 패턴 못 찾음 — 수동 확인: {f05}")
    # fallback: 앞에 삽입
    insert = """
# ── patch_step5 적용 ──────────────────────────────────
if "ko_baseline_cst" in data and data["ko_baseline_cst"]:
    HIGH_CST = data["ko_baseline_cst"]
    print(f"  [patch] Baseline: {data.get('ko_baseline_label','custom')}")
else:
    HIGH_CST = avg_cst(high_clones)
"""
    content = content.replace(
        "# ── Baseline ──",
        insert + "# ── Baseline ──"
    )
    with open(f05, "w", encoding="utf-8") as f:
        f.write(content)

print(f"""
{'='*60}
  완료!

  1. pkl에 feasible baseline constraint 저장
  2. 05_ko_screening.py 자동 패치

  다음 실행:
  python scripts/steps/05_ko_screening.py --dataset {DATASET}

  또는 전체:
  python run_pipeline.py --dataset {DATASET} --steps 5,6
{'='*60}
""")
