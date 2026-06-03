"""
05_ko_screening.py — Reaction KO In Silico 스크리닝
High producer 평균 constraints 기준으로
핵심 대사 반응을 knock-out했을 때 objective 변화율 계산.

사용:
  python scripts/steps/05_ko_screening.py --dataset practice_20aa
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
print(f"  05_ko_screening.py — {DATASET}")
print("=" * 60)

# ── 데이터 로딩 ───────────────────────────────────
pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
if not os.path.exists(pkl_path):
    print(f"  !! {pkl_path} 없음 — 04_run_fva.py 먼저 실행")
    sys.exit(1)
with open(pkl_path, "rb") as f:
    data = pickle.load(f)

all_constraints = data["all_constraints"]
high_clones     = data["high_clones"]

# ── 모델 로딩 ─────────────────────────────────────
model = load_model(verbose=True)

# ── High 평균 constraints ─────────────────────────
def avg_constraints(clone_list):
    cst_list = [all_constraints[c] for c in clone_list if c in all_constraints]
    if not cst_list: return {}
    all_ids = set()
    for c in cst_list: all_ids.update(c.keys())
    return {rid: (np.mean([c[rid][0] for c in cst_list if rid in c]),
                  np.mean([c[rid][1] for c in cst_list if rid in c]))
            for rid in all_ids}

high_cst = avg_constraints(high_clones)
print(f"\n  High avg constraints: {len(high_cst)}개")

# ── Baseline ──────────────────────────────────────
with model:
    apply_bounds(model, high_cst)
    model.objective = OBJ_RXN
    base_sol = model.optimize()
    base_obj = base_sol.objective_value if base_sol.status == "optimal" else 0.0

print(f"  Baseline (High avg): {base_sol.status}  obj={base_obj:.5f}")
if base_sol.status != "optimal":
    print("  !! baseline infeasible — constraints 확인 필요")
    print("  !! 03_run_fba.py 결과를 먼저 확인하세요")
    sys.exit(1)

# ── KO 타겟 정의 ──────────────────────────────────
RXN_KO_TARGETS = {
    "LDH (Lactate DH)":      ["lactate dehydrogenase", "LDH_"],
    "PDH kinase":            ["pyruvate dehydrogenase kinase", "PDHK", "PDC"],
    "Glutaminase":           ["glutaminase", "GLS_"],
    "Pyruvate CX":           ["pyruvate carboxylase", "PC_"],
    "Isocitrate DH":         ["isocitrate dehydrogenase", "IDH"],
    "Malate DH":             ["malate dehydrogenase", "MDH"],
    "Citrate synthase":      ["citrate synthase", "CS_"],
    "Malic enzyme":          ["malic enzyme", "ME1_", "ME2_"],
    "Fatty acid syn (FASN)": ["fatty acid synthase", "FASN"],
    "ATP-citrate lyase":     ["ATP-citrate", "ACLY"],
    "Glutamate DH":          ["glutamate dehydrogenase", "GLUD"],
    "Pyruvate kinase":       ["pyruvate kinase", "PKM"],
    "G6P DH":                ["glucose-6-phosphate dehydrogenase", "G6PD"],
}

def find_rxns(kws):
    return [r for r in model.reactions
            if any(k.lower() in f"{r.id} {r.name or ''}".lower() for k in kws)]

# ── KO 스크리닝 ───────────────────────────────────
print(f"\n  [Reaction KO 스크리닝]")
print(f"  {'':2s} {'타겟':25s} {'n_rxns':6s} {'ko_obj':>10s} {'delta(%)':>10s}  status")
print("  " + "─" * 68)

ko_rows = []
for label, kws in RXN_KO_TARGETS.items():
    rxns = find_rxns(kws)
    if not rxns:
        print(f"  ? {label:25s} — 모델에서 찾을 수 없음")
        ko_rows.append({"target": label, "n_rxns": 0,
                        "ko_obj": None, "delta_pct": None, "status": "not_found"})
        continue
    try:
        with model:
            apply_bounds(model, high_cst)
            model.objective = OBJ_RXN
            for r in rxns:
                r.knock_out()
            sol = model.optimize()
            ko_obj = sol.objective_value if sol.status == "optimal" else 0.0
            dp = (ko_obj - base_obj) / base_obj * 100 if base_obj else 0.0
            icon = "+" if dp > 0.5 else ("-" if dp < -2 else "=")
            print(f"  {icon} {label:25s} {len(rxns):6d} {ko_obj:10.5f} {dp:+10.1f}%  {sol.status}")
            ko_rows.append({"target": label, "n_rxns": len(rxns),
                            "ko_obj": ko_obj, "baseline_obj": base_obj,
                            "delta_pct": dp, "status": sol.status})
    except Exception as e:
        print(f"  ! {label}: {e}")
        ko_rows.append({"target": label, "n_rxns": len(rxns),
                        "ko_obj": None, "delta_pct": None, "status": "error"})

ko_df = pd.DataFrame(ko_rows)

# ── 결과 요약 ─────────────────────────────────────
ko_valid = ko_df[ko_df["delta_pct"].notna()]
improved = ko_valid[ko_valid["delta_pct"] > 0.5].sort_values("delta_pct", ascending=False)
decreased = ko_valid[ko_valid["delta_pct"] < -2].sort_values("delta_pct")

if len(improved) > 0:
    print(f"\n  ★ 생산성 향상 타겟 ({len(improved)}개):")
    for _, r in improved.iterrows():
        print(f"    {r['target']:25s} +{r['delta_pct']:.1f}%")
if len(decreased) > 0:
    print(f"\n  ▼ 필수 반응 (건드리면 안 됨, {len(decreased)}개):")
    for _, r in decreased.iterrows():
        print(f"    {r['target']:25s} {r['delta_pct']:.1f}%")

# ── 개선 시나리오 시뮬레이션 ──────────────────────
print(f"\n  [개선 시나리오 시뮬레이션]")
top_ko = improved.head(3)["target"].tolist()

media_extra = {
    "EX_leu_L_e": (-0.015, 0), "EX_ile_L_e": (-0.012, 0),
    "EX_val_L_e": (-0.012, 0), "EX_glu_L_e": (-0.010, 0),
    "EX_asn_L_e": (-0.010, 0),
}

# Low producer constraints
low_clones = data["low_clones"]
low_cst = avg_constraints(low_clones)

impr_scenarios = [
    ("Baseline (High)",       high_cst, []),
    ("Reference (Low)",       low_cst,  []),
]
for t in top_ko[:2]:
    impr_scenarios.append((f"KO: {t[:18]}", high_cst, [t]))
if len(top_ko) >= 2:
    impr_scenarios.append((f"KO: {top_ko[0][:10]}+{top_ko[1][:10]}",
                           high_cst, top_ko[:2]))
impr_scenarios.append(("Media (BCAA+Glu+Asn)", high_cst, []))
if top_ko:
    impr_scenarios.append(("KO + Media", high_cst, top_ko[:1]))

impr_rows = []
print(f"  {'시나리오':30s} {'obj':>10s} {'vs_High(%)':>12s}  status")
print("  " + "─" * 60)
for label, base_c, ko_labels in impr_scenarios:
    ko_rxns = []
    for kl in ko_labels:
        ko_rxns += find_rxns(RXN_KO_TARGETS.get(kl, []))
    with model:
        apply_bounds(model, base_c)
        if "Media" in label:
            apply_bounds(model, media_extra)
        model.objective = OBJ_RXN
        for r in ko_rxns:
            r.knock_out()
        sol = model.optimize()
        obj = sol.objective_value if sol.status == "optimal" else 0.0
    vs = (obj / base_obj - 1) * 100 if base_obj else 0
    print(f"  {label:30s} {obj:10.5f} {vs:+12.1f}%  {sol.status}")
    impr_rows.append({"scenario": label, "obj": obj,
                      "vs_high_pct": vs, "status": sol.status})

impr_df = pd.DataFrame(impr_rows)

# ── 저장 ──────────────────────────────────────────
data["ko_df"]        = ko_df
data["impr_df"]      = impr_df
data["base_obj"]     = base_obj
data["top_ko"]       = top_ko
data["high_cst_avg"] = high_cst
data["low_cst_avg"]  = low_cst
with open(pkl_path, "wb") as f:
    pickle.dump(data, f)

save_table(ko_df,     "ko_screening.csv",   DATASET)
save_table(impr_df,   "improvement_scenarios.csv", DATASET)

print(f"\n  OK  05_ko_screening 완료")
print(f"  다음: python scripts/steps/06_figures.py --dataset {DATASET}")
