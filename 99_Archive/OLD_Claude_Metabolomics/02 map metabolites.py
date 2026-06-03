"""
02_map_metabolites.py — 대사물질명 → iCHO3K Exchange ID 매핑 확인
01에서 생성된 pkl을 읽어 각 대사물질이 모델에 존재하는지 검증하고
FBA constraints dict를 생성한다.

사용:
  python scripts/steps/02_map_metabolites.py --dataset practice_20aa
"""
import sys, os, argparse, warnings, pickle
warnings.filterwarnings("ignore")

def _find_root():
    cur = os.path.dirname(os.path.abspath(__file__))
    for _ in range(5):
        if os.path.exists(os.path.join(cur, 'src', 'config.py')):
            return cur
        cur = os.path.dirname(cur)
    return cur
ROOT = _find_root()
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from src.config import *
from src.fba_utils import load_model

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa",
                    choices=["sowa2020", "practice_20aa", "own_experiment"])
args = parser.parse_args()
DATASET = args.dataset

print("=" * 60)
print(f"  02_map_metabolites.py — {DATASET}")
print("=" * 60)

# ── 01 결과 로딩 ──────────────────────────────────
pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
if not os.path.exists(pkl_path):
    print(f"  !! {pkl_path} 없음 — 01_load_data.py 먼저 실행")
    sys.exit(1)
with open(pkl_path, "rb") as f:
    data = pickle.load(f)

all_rates = data["all_rates"]
COL2EX    = data.get("COL2EX", {})
COL2NAME  = data.get("COL2NAME", {})
MET_COLS  = data.get("MET_COLS", [])

# ── 모델 로딩 ─────────────────────────────────────
model = load_model(verbose=True)
model_rxn_ids = {r.id for r in model.reactions}

# ── 매핑 검증 ─────────────────────────────────────
print(f"\n  [Exchange ID 검증]")
print(f"  {'대사물질':20s} {'Exchange ID':30s} {'상태':10s} {'lb/ub':20s}")
print("  " + "─" * 82)

mapping_records = []
confirmed = {}   # met_name → exchange_id (모델에 존재하는 것만)

for col in MET_COLS:
    met_name = COL2NAME.get(col, col)
    # 우선순위: EXCHANGE_IDS(검증) → COL2EX(파일 내 매핑) → 없음
    ex_id = EXCHANGE_IDS.get(met_name) \
         or EXCHANGE_IDS.get(col) \
         or COL2EX.get(col, "")

    if ex_id and ex_id in model_rxn_ids:
        rxn = model.reactions.get_by_id(ex_id)
        status = "OK"
        bounds = f"lb={rxn.lower_bound:.0f}, ub={rxn.upper_bound:.0f}"
        confirmed[col] = ex_id
    elif ex_id:
        status = "NOT IN MODEL"
        bounds = ""
    else:
        status = "NO MAPPING"
        bounds = ""

    icon = "✔" if status == "OK" else "✘"
    print(f"  {icon} {col:20s} {ex_id:30s} {status:10s} {bounds}")
    mapping_records.append({
        "column": col, "metabolite_name": met_name,
        "exchange_id": ex_id, "status": status, "in_model": status == "OK",
    })

mapping_df = pd.DataFrame(mapping_records)
n_ok = (mapping_df["status"] == "OK").sum()
print(f"\n  매핑 성공: {n_ok}/{len(mapping_df)}개")

# ── FBA Constraints 생성 ──────────────────────────
print(f"\n  [FBA Constraints 생성]")
all_constraints = {}
BUFFER = 0.20

for clone, rates in all_rates.items():
    cst = {}
    for col, info in rates.items():
        if col not in confirmed:
            continue
        ex_id = confirmed[col]
        q = info["rate_mmol_gDCWh"]
        if q == 0 or abs(q) < 1e-9:
            continue
        if q < 0:   # uptake
            lb, ub = q * (1 + BUFFER), q * (1 - BUFFER)
        else:       # secretion
            lb, ub = q * (1 - BUFFER), q * (1 + BUFFER)
        cst[ex_id] = (min(lb, ub), max(lb, ub))

    all_constraints[clone] = cst
    print(f"  {clone:12s}: {len(cst)}개 constraint")

# ── 저장 ──────────────────────────────────────────
data["confirmed_mapping"] = confirmed
data["all_constraints"]   = all_constraints
data["mapping_df"]        = mapping_df

with open(pkl_path, "wb") as f:
    pickle.dump(data, f)

out = results_dir(DATASET, "tables")
mapping_df.to_csv(os.path.join(out, "metabolite_mapping.csv"), index=False)
print(f"\n  [saved] results/{DATASET}/tables/metabolite_mapping.csv")
print(f"  [saved] data/processed/rates_{DATASET}.pkl (updated)")
print(f"\n  OK  02_map_metabolites 완료")
print(f"  다음: python scripts/steps/03_run_fba.py --dataset {DATASET}")
