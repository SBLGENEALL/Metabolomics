"""
04_run_fva.py — FVA (Flux Variability Analysis)
High vs Low producer 평균 constraints로 FVA 실행.
각 Exchange reaction의 feasible flux 범위 계산.

사용:
  python scripts/steps/04_run_fva.py --dataset practice_20aa
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
from cobra.flux_analysis import flux_variability_analysis

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa",
                    choices=["sowa2020", "practice_20aa", "own_experiment"])
parser.add_argument("--fraction", type=float, default=0.9,
                    help="FVA optimality fraction (기본 0.9)")
args = parser.parse_args()
DATASET  = args.dataset
FRACTION = args.fraction

print("=" * 60)
print(f"  04_run_fva.py — {DATASET}  (fraction={FRACTION})")
print("=" * 60)

# ── 데이터 로딩 ───────────────────────────────────
pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
if not os.path.exists(pkl_path):
    print(f"  !! {pkl_path} 없음 — 03_run_fba.py 먼저 실행")
    sys.exit(1)
with open(pkl_path, "rb") as f:
    data = pickle.load(f)

all_constraints = data["all_constraints"]
high_clones     = data["high_clones"]
low_clones      = data["low_clones"]
COL2EX          = data.get("COL2EX", {})
COL2NAME        = data.get("COL2NAME", {})

# ── 모델 로딩 ─────────────────────────────────────
model = load_model(verbose=True)
model_rxn_ids = {r.id for r in model.reactions}

# ── High/Low 평균 constraints 계산 ────────────────
def avg_constraints(clone_list):
    cst_list = [all_constraints[c] for c in clone_list if c in all_constraints]
    if not cst_list:
        return {}
    all_ids = set()
    for c in cst_list:
        all_ids.update(c.keys())
    avg = {}
    for rid in all_ids:
        lbs = [c[rid][0] for c in cst_list if rid in c]
        ubs = [c[rid][1] for c in cst_list if rid in c]
        if lbs:
            avg[rid] = (np.mean(lbs), np.mean(ubs))
    return avg

high_cst = avg_constraints(high_clones)
low_cst  = avg_constraints(low_clones)
print(f"\n  High avg constraints: {len(high_cst)}개")
print(f"  Low  avg constraints: {len(low_cst)}개")

# Exchange ID 목록 (모델에 있는 것만)
ex_ids_in_model = [r.id for r in model.exchanges]

# ── FVA 실행 ──────────────────────────────────────
print(f"\n  [FVA 실행 — Exchange reactions only]")

fva_results = {}
for label, cst in [("High", high_cst), ("Low", low_cst)]:
    print(f"  → {label} producer ({len(cst)} constraints)...")
    with model:
        apply_bounds(model, cst)
        model.objective = OBJ_RXN
        # feasibility 먼저 확인
        pre = model.optimize()
        if pre.status != "optimal":
            print(f"    !! {label}: infeasible — FVA 건너뜀")
            continue
        print(f"    baseline obj = {pre.objective_value:.5f}")
        try:
            fva = flux_variability_analysis(
                model,
                fraction_of_optimum=FRACTION,
                processes=FVA_PROCESSES,
                reaction_list=list(model.exchanges),
            )
            fva["mean"]  = (fva["minimum"] + fva["maximum"]) / 2
            fva["range"] = fva["maximum"]  - fva["minimum"]
            fva["label"] = label
            fva_results[label] = fva
            print(f"    OK: {len(fva)} exchange reactions")
        except Exception as e:
            print(f"    !! FVA 오류: {e}")

# ── 비교표 생성 ───────────────────────────────────
if len(fva_results) == 2:
    fva_h = fva_results["High"][["minimum","maximum","mean","range"]].rename(
        columns={c: f"{c}_High" for c in ["minimum","maximum","mean","range"]})
    fva_l = fva_results["Low"][["minimum","maximum","mean","range"]].rename(
        columns={c: f"{c}_Low" for c in ["minimum","maximum","mean","range"]})
    fva_compare = fva_h.join(fva_l, how="outer")
    fva_compare["delta_mean"] = (fva_compare["mean_High"] - fva_compare["mean_Low"])
    fva_compare["abs_delta"]  = fva_compare["delta_mean"].abs()
    fva_compare = fva_compare.sort_values("abs_delta", ascending=False)

    # 주요 exchange 상위 15개 출력
    key_ex_names = {
        "EX_glc_e":"Glucose", "EX_lac_L_e":"Lactate",
        "EX_gln_L_e":"Glutamine", "EX_nh4_e":"NH4+",
        "EX_ala_L_e":"Alanine",   "EX_glu_L_e":"Glutamate",
    }
    print(f"\n  [주요 Exchange FVA 결과]")
    print(f"  {'반응':25s} {'High mean':>10s} {'Low mean':>10s} {'Delta':>10s}")
    print("  " + "─" * 58)
    for eid, name in key_ex_names.items():
        if eid in fva_compare.index:
            row = fva_compare.loc[eid]
            h_m = row.get("mean_High", 0)
            l_m = row.get("mean_Low",  0)
            d   = h_m - l_m
            print(f"  {name:25s} {h_m:10.4f} {l_m:10.4f} {d:+10.4f}")

    data["fva_compare"] = fva_compare
    save_table(fva_compare.reset_index(), "fva_comparison.csv", DATASET)

# 개별 저장
for label, fva in fva_results.items():
    save_table(fva.reset_index(), f"fva_{label.lower()}.csv", DATASET)

data["fva_results"] = fva_results
with open(pkl_path, "wb") as f:
    pickle.dump(data, f)

print(f"\n  OK  04_run_fva 완료")
print(f"  다음: python scripts/steps/05_ko_screening.py --dataset {DATASET}")
