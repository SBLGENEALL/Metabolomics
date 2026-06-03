
"""
02_map_metabolites.py — metabolite-to-iCHO3K exchange mapping + constraint policy

Default policy: production_relaxed. This prevents apparent positive amino-acid
net rates from being treated as secretion-only constraints that block IgG/mAb
synthesis.
"""
import sys, os, argparse, warnings, pickle
warnings.filterwarnings("ignore")

def _find_root():
    cur = os.path.dirname(os.path.abspath(__file__))
    for _ in range(6):
        if os.path.exists(os.path.join(cur, 'src', 'config.py')):
            return cur
        cur = os.path.dirname(cur)
    return cur
ROOT = _find_root()
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from src.config import *  # noqa
from src.fba_utils import load_model
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa", choices=["sowa2020", "practice_20aa", "own_experiment"])
parser.add_argument("--constraint_policy", default=os.environ.get("CHO_CONSTRAINT_POLICY", CONSTRAINT_POLICY_DEFAULT),
                    choices=["strict", "production_relaxed"],
                    help="strict=all rates directional; production_relaxed=AA nutrients flexible if apparent positive")
args = parser.parse_args()
DATASET = args.dataset
POLICY = args.constraint_policy

print("=" * 65)
print(f"  02_map_metabolites.py — {DATASET}")
print("=" * 65)
print(f"  Constraint policy: {POLICY}")

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

model = load_model(verbose=True)
model_rxn_ids = {r.id for r in model.reactions}

print(f"\n  [Exchange ID 검증]")
print(f"  {'대사물질':20s} {'Exchange ID':30s} {'상태':10s} {'lb/ub':20s}")
print("  " + "─" * 82)

mapping_records = []
confirmed = {}
for col in MET_COLS:
    met_name = COL2NAME.get(col, col)
    ex_id = EXCHANGE_IDS.get(met_name) or EXCHANGE_IDS.get(col) or COL2EX.get(col, "")
    if ex_id and ex_id in model_rxn_ids:
        rxn = model.reactions.get_by_id(ex_id)
        status = "OK"
        bounds = f"lb={rxn.lower_bound:.0f}, ub={rxn.upper_bound:.0f}"
        confirmed[col] = ex_id
    elif ex_id:
        status = "NOT IN MODEL"; bounds = ""
    else:
        status = "NO MAPPING"; bounds = ""
    icon = "✔" if status == "OK" else "✘"
    print(f"  {icon} {col:20s} {ex_id:30s} {status:10s} {bounds}")
    mapping_records.append({"column": col, "metabolite_name": met_name, "exchange_id": ex_id, "status": status, "in_model": status == "OK"})

mapping_df = pd.DataFrame(mapping_records)
print(f"\n  매핑 성공: {(mapping_df['status'] == 'OK').sum()}/{len(mapping_df)}개")

BUFFER = float(globals().get("CONSTRAINT_BUFFER", 0.20))
MIN_NUT = float(globals().get("MIN_NUTRIENT_UPTAKE_CAPACITY", 0.02))
CORE_HARD = set(globals().get("CORE_HARD_METABOLITES", {"Glucose", "Lactate", "NH4+", "Ammonia"}))
FLEX_NUT = set(globals().get("FLEXIBLE_NUTRIENT_METABOLITES", set()))

def strict_bound(q):
    if q < 0:
        return (q * (1 + BUFFER), 0.0)
    return (0.0, q * (1 + BUFFER))

def production_relaxed_bound(met_name, q):
    if met_name in CORE_HARD:
        return strict_bound(q)
    if met_name in FLEX_NUT:
        cap = max(abs(q) * (1 + BUFFER), MIN_NUT)
        if q < 0:
            return (q * (1 + BUFFER), 0.0)
        return (-cap, q * (1 + BUFFER))
    return strict_bound(q)

print(f"\n  [FBA Constraints 생성]")
print("  strict constraints are saved as all_constraints_strict")
print(f"  active constraints use policy={POLICY}")
all_constraints_strict = {}
all_constraints_active = {}
constraint_records = []
for clone, rates in all_rates.items():
    cst_strict = {}
    cst_active = {}
    clone_relaxed = 0
    for col, info in rates.items():
        if col not in confirmed:
            continue
        ex_id = confirmed[col]
        met_name = COL2NAME.get(col, col)
        q = float(info["rate_mmol_gDCWh"])
        if abs(q) < 1e-9:
            continue
        sb = strict_bound(q)
        ab = sb if POLICY == "strict" else production_relaxed_bound(met_name, q)
        cst_strict[ex_id] = sb
        cst_active[ex_id] = ab
        note = "directional_measured"
        if ab != sb:
            note = "relaxed_positive_nutrient" if q > 0 else "relaxed_nutrient"
            clone_relaxed += 1
        constraint_records.append({
            "clone": clone, "column": col, "metabolite_name": met_name,
            "exchange_id": ex_id, "rate": q,
            "strict_lb": sb[0], "strict_ub": sb[1],
            "active_lb": ab[0], "active_ub": ab[1],
            "policy": POLICY, "note": note,
        })
    all_constraints_strict[clone] = cst_strict
    all_constraints_active[clone] = cst_active
    print(f"  {clone:12s}: {len(cst_active)} active constraints ({clone_relaxed} relaxed nutrient bounds)")

constraint_policy_df = pd.DataFrame(constraint_records)
data["confirmed_mapping"] = confirmed
data["all_constraints_strict"] = all_constraints_strict
data["all_constraints"] = all_constraints_active
data["constraint_policy"] = POLICY
data["constraint_policy_df"] = constraint_policy_df
data["mapping_df"] = mapping_df
with open(pkl_path, "wb") as f:
    pickle.dump(data, f)

out = results_dir(DATASET, "tables")
mapping_df.to_csv(os.path.join(out, "metabolite_mapping.csv"), index=False)
constraint_policy_df.to_csv(os.path.join(out, "constraint_policy_bounds.csv"), index=False)
print(f"\n  [saved] results/{DATASET}/tables/metabolite_mapping.csv")
print(f"  [saved] results/{DATASET}/tables/constraint_policy_bounds.csv")
print(f"  [saved] data/processed/rates_{DATASET}.pkl (updated)")
print(f"\n  OK  02_map_metabolites 완료")
print(f"  다음: python scripts/steps/03_run_fba.py --dataset {DATASET}")
