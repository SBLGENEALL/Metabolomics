"""
03_run_fba.py — 클론별 FBA 실행
02에서 생성된 constraints로 iCHO3K prod 모델 FBA 실행.
목적함수: biomass_cho_prod

사용:
  python scripts/steps/03_run_fba.py --dataset practice_20aa
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
from src.fba_utils import load_model, apply_bounds, save_table

import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa",
                    choices=["sowa2020", "practice_20aa", "own_experiment"])
args = parser.parse_args()
DATASET = args.dataset

print("=" * 60)
print(f"  03_run_fba.py — {DATASET}")
print("=" * 60)

# ── 데이터 로딩 ───────────────────────────────────
pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
if not os.path.exists(pkl_path):
    print(f"  !! {pkl_path} 없음 — 02_map_metabolites.py 먼저 실행")
    sys.exit(1)
with open(pkl_path, "rb") as f:
    data = pickle.load(f)

all_constraints = data.get("all_constraints", {})
all_rates       = data["all_rates"]
clones          = data["clones"]
igG_day14       = data.get("igG_day14", {})
high_clones     = data["high_clones"]
low_clones      = data["low_clones"]

if not all_constraints:
    print("  !! constraints 없음 — 02_map_metabolites.py 먼저 실행")
    sys.exit(1)

# ── 모델 로딩 ─────────────────────────────────────
model = load_model(verbose=True)
model_rxn_ids = {r.id for r in model.reactions}

# ── FBA 실행 ──────────────────────────────────────
print(f"\n  [FBA 실행 — {OBJ_RXN}]")
print(f"  {'':2s} {'Clone':12s} {'n_cst':5s} {'status':10s} {'obj':>10s} {'IgG(D14)':>10s}")
print("  " + "─" * 55)

fba_results = {}

for clone in clones:
    if clone not in all_constraints:
        continue
    cst = all_constraints[clone]
    with model:
        n = apply_bounds(model, cst)
        model.objective = OBJ_RXN
        sol = model.optimize()
        obj = sol.objective_value if sol.status == "optimal" else 0.0
        ex_flux = {}
        if sol.status == "optimal":
            ex_flux = {r.id: sol.fluxes[r.id]
                       for r in model.exchanges
                       if abs(sol.fluxes.get(r.id, 0)) > 1e-9}

    igG  = igG_day14.get(clone, 0)
    tag  = "★" if clone in high_clones else ("▼" if clone in low_clones else " ")
    print(f"  {tag} {clone:12s} {n:5d} {sol.status:10s} {obj:10.5f} {igG:10.1f}")

    fba_results[clone] = {
        "clone": clone, "obj": obj, "status": sol.status,
        "n_constraints": n, "igG_day14": igG,
        "ex_flux": ex_flux, "constraints": cst,
    }

# ── 통계 요약 ─────────────────────────────────────
fba_df = pd.DataFrame([{
    "clone":         v["clone"],
    "obj":           v["obj"],
    "status":        v["status"],
    "n_constraints": v["n_constraints"],
    "igG_day14":     v["igG_day14"],
    "group":         "High" if v["clone"] in high_clones
                     else ("Low" if v["clone"] in low_clones else "Mid"),
} for v in fba_results.values()])
fba_df = fba_df.sort_values("igG_day14", ascending=False)

optimal_n = (fba_df["status"] == "optimal").sum()
print(f"\n  Optimal: {optimal_n}/{len(fba_df)}개")

if optimal_n > 0:
    h_mean = fba_df[fba_df["group"]=="High"]["obj"].mean()
    l_mean = fba_df[fba_df["group"]=="Low"]["obj"].mean()
    print(f"  High producer 평균 obj: {h_mean:.5f}")
    print(f"  Low  producer 평균 obj: {l_mean:.5f}")
    if l_mean > 0:
        print(f"  High/Low 비율:          {h_mean/l_mean:.2f}×")

# Lac/Glc ratio 계산
G, L = "EX_glc_e", "EX_lac_L_e"
lac_glc_rows = []
for clone, cst in all_constraints.items():
    glc = abs(cst.get(G, (0, 0))[0])
    lac = cst.get(L, (0, 0))[0]
    ratio = lac / glc if (glc > 0 and lac >= 0) else 0
    lac_glc_rows.append({"clone": clone, "lac_glc_ratio": ratio,
                          "igG_day14": igG_day14.get(clone, 0)})
lac_glc_df = pd.DataFrame(lac_glc_rows).sort_values("igG_day14", ascending=False)

print(f"\n  [Lac/Glc ratio]  (<1.0=TCA 효율, >1.0=Warburg)")
for _, row in lac_glc_df.iterrows():
    bar = "★" if row["clone"] in high_clones else ("▼" if row["clone"] in low_clones else " ")
    print(f"  {bar} {row['clone']:12s} ratio={row['lac_glc_ratio']:.3f}  IgG={row['igG_day14']:.0f}")

# ── Exchange flux 비교표 ──────────────────────────
key_ex = [G, L, "EX_gln_L_e", "EX_nh4_e", "EX_ala_L_e", "EX_glu_L_e"]
ex_compare = pd.DataFrame({
    clone: {ex: res["ex_flux"].get(ex, 0) for ex in key_ex}
    for clone, res in fba_results.items()
}).T
ex_compare.index.name = "clone"

# ── 저장 ──────────────────────────────────────────
data["fba_results"]  = fba_results
data["fba_df"]       = fba_df
data["lac_glc_df"]   = lac_glc_df
data["ex_compare"]   = ex_compare
with open(pkl_path, "wb") as f:
    pickle.dump(data, f)

save_table(fba_df,       "fba_results.csv",       DATASET)
save_table(lac_glc_df,   "lac_glc_ratio.csv",     DATASET)
save_table(ex_compare.reset_index(), "fba_exchange_compare.csv", DATASET)

print(f"\n  OK  03_run_fba 완료")
print(f"  다음: python scripts/steps/04_run_fva.py --dataset {DATASET}")
