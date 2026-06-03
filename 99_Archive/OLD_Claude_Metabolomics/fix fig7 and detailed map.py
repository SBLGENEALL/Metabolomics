"""
fix_fig7_and_detailed_map.py
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
1. Fig 7 KO Screening 수정
   - 반응 탐색 로직 개선 (iCHO3K 실제 ID 기반)
   - 모든 타겟 표시 (not found 포함)
   - 그래프 크기 자동 조정

2. 상세 Flux Pathway Map
   - 40+ 노드: Glycolysis, TCA, PPP, Gln, AA, FA, IgG
   - 실제 flux 값 + 배율 표시
   - High vs Low 나란히 비교

python scripts/steps/fix_fig7_and_detailed_map.py --dataset practice_20aa
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
from src.fba_utils import load_model, apply_bounds, save_figure, save_table

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.cm as cm
from scipy import stats

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa")
args   = parser.parse_args()
DATASET = args.dataset

print("="*65)
print(f"  fix_fig7_and_detailed_map.py — {DATASET}")
print("="*65)

# ── 로딩 ──────────────────────────────────────────────
pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
with open(pkl_path, "rb") as f:
    data = pickle.load(f)

all_constraints = data.get("all_constraints", {})
igG_day14       = data.get("igG_day14", {})
high_clones     = data["high_clones"]
low_clones      = data["low_clones"]
clones          = data["clones"]
interval        = data.get("day_interval", (7,10))
OBJ = data.get("selected_objective", OBJ_RXN)

P = PALETTE

model = load_model(verbose=False)
all_rxn_ids = {r.id for r in model.reactions}

def safe_apply(model, constraints):
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

def avg_cst(clone_list):
    cst_list = [all_constraints[c] for c in clone_list if c in all_constraints]
    if not cst_list: return {}
    ids = set()
    for c in cst_list: ids.update(c.keys())
    return {rid: (np.mean([c[rid][0] for c in cst_list if rid in c]),
                  np.mean([c[rid][1] for c in cst_list if rid in c]))
            for rid in ids}

HIGH_CST = avg_cst(high_clones)

# ── FBA flux 추출 ──────────────────────────────────────
flux_high, flux_low = {}, {}
for label, cl, fd in [("High", high_clones, flux_high),
                       ("Low",  low_clones,  flux_low)]:
    with model:
        safe_apply(model, avg_cst(cl))
        model.objective = OBJ
        sol = model.optimize()
        if sol.status == "optimal":
            for r in model.reactions:
                fd[r.id] = sol.fluxes[r.id]
            print(f"  {label}: obj={sol.objective_value:.5f}")
        else:
            print(f"  {label}: {sol.status}")

# ════════════════════════════════════════════════════════
# FIG 7 수정: KO Screening
# ════════════════════════════════════════════════════════
print("\n[1] Fig 7 KO Screening 수정")
print("    Ref: Lim et al. (2010) Metab Eng 12:614")
print("    Ref: Wilkens et al. (2011) J Biotechnol 153:51")

# iCHO3K 실제 반응명 기반 KO 타겟
# 모델에서 직접 검색한 정확한 키워드
RXN_KO_TARGETS = {
    "LDH\n(Lactate DH)":          ["LDH_","LDHA","LDHB","lactate dehydro"],
    "PDH kinase\n(PDHK1/2)":      ["PDC","PDHK","pyruvate dehydrogenase kinase"],
    "Glutaminase\n(GLS)":          ["GLS","glutaminase","GLNase"],
    "Pyruvate CX\n(PC)":           ["PC_","pyruvate carboxylase","PYC"],
    "Isocitrate DH\n(IDH1/2)":     ["IDH1","IDH2","IDH_","isocitrate dehydro"],
    "Malate DH\n(MDH1/2)":         ["MDH1","MDH2","MDH_","malate dehydro"],
    "Citrate syn\n(CS)":           ["CS_","citrate syn","CITS"],
    "Malic enzyme\n(ME1/2)":       ["ME1","ME2","malic enzyme","MALA"],
    "Fatty acid syn\n(FASN)":      ["FASN","FAS","fatty acid syn"],
    "ATP-citrate ly\n(ACLY)":      ["ACLY","ATP cit","citrate lyase","ATPCL"],
    "Glutamate DH\n(GLUD)":        ["GLUD","glutamate dehydro","GLDH"],
    "Pyruvate kin\n(PKM)":         ["PKM","pyruvate kin","PYK"],
    "G6P DH\n(G6PD)":              ["G6PD","G6PDH","glucose-6-phos dehydro"],
    "Phosphogluco\n(PGI)":         ["PGI","phosphoglucose isomerase","GPI"],
    "Aspartate AT\n(GOT)":         ["GOT","ASPTA","aspartate amino"],
    "PSAT1\n(Serine syn)":         ["PSAT","serine","phosphoserine amino"],
    "Alanine AT\n(ALAT)":          ["ALAT","alanine amino","ALATC"],
}

def find_rxns_improved(kws):
    """
    iCHO3K 모델에서 반응 탐색 - 다단계 매칭
    1) 정확한 ID 매칭
    2) 이름 포함 매칭  
    3) subsystem 매칭
    """
    found = set()
    for r in model.reactions:
        r_str = f"{r.id} {r.name or ''} {getattr(r,'subsystem','') or ''}".lower()
        for kw in kws:
            if kw.lower() in r_str:
                found.add(r.id)
                break
    return [model.reactions.get_by_id(rid) for rid in found]

# baseline
with model:
    safe_apply(model, HIGH_CST)
    model.objective = OBJ
    base_sol = model.optimize()
    base_obj = base_sol.objective_value if base_sol.status=="optimal" else 0

print(f"  Baseline: {base_sol.status}  obj={base_obj:.5f}")

ko_rows = []
print(f"\n  {'Target':25s} {'n':>4s} {'delta%':>8s}")
print("  " + "─"*42)

for label, kws in RXN_KO_TARGETS.items():
    rxns = find_rxns_improved(kws)
    display = label.replace("\n", " ")

    if not rxns:
        print(f"  {'not found':4s} {display:30s}  n=0")
        ko_rows.append({"target": label, "n_rxns": 0,
                        "delta_pct": None, "status": "not_found",
                        "ko_obj": None})
        continue

    try:
        with model:
            safe_apply(model, HIGH_CST)
            model.objective = OBJ
            for r in rxns: r.knock_out()
            sol = model.optimize()
            ko_obj = sol.objective_value if sol.status=="optimal" else 0
            dp = (ko_obj - base_obj)/base_obj*100 if base_obj else 0
            icon = "+" if dp > 0.5 else ("-" if dp < -2 else "=")
            print(f"  {icon}    {display:30s}  n={len(rxns)}  {dp:+.1f}%")
            ko_rows.append({"target": label, "n_rxns": len(rxns),
                            "delta_pct": dp, "ko_obj": ko_obj,
                            "status": sol.status})
    except Exception as e:
        ko_rows.append({"target": label, "n_rxns": len(rxns),
                        "delta_pct": None, "status": "error", "ko_obj": None})

ko_df = pd.DataFrame(ko_rows)
save_table(ko_df, "ko_screening.csv", DATASET)

# ── Fig 7 그리기 ───────────────────────────────────────
all_targets = ko_df["target"].tolist()
n_total     = len(all_targets)

# 높이 자동 조정: 타겟 수에 맞게
fig_h = max(8, n_total * 0.55)
fig, ax = plt.subplots(figsize=(11, fig_h))

y_pos = np.arange(n_total)

# 색상 및 bar 값
bar_vals = []
bar_cols = []
bar_patterns = []

for _, row in ko_df.iterrows():
    dp = row["delta_pct"]
    status = row["status"]

    if status == "not_found":
        bar_vals.append(0)
        bar_cols.append("#E8E8E8")
        bar_patterns.append("////")    # 빗금 = 모델에 없음
    elif dp is None:
        bar_vals.append(0)
        bar_cols.append("#E8E8E8")
        bar_patterns.append("....")
    elif dp > 0.5:
        bar_vals.append(dp)
        bar_cols.append(P["impr"])
        bar_patterns.append("")
    elif dp < -2:
        bar_vals.append(dp)
        bar_cols.append(P["danger"])
        bar_patterns.append("")
    else:
        bar_vals.append(dp)
        bar_cols.append(P["neutral"])
        bar_patterns.append("")

bars = ax.barh(y_pos, bar_vals, color=bar_cols,
               edgecolor="k", lw=0.5, zorder=3, height=0.7)

# 빗금 패턴 (모델에 없는 타겟)
for bar, pat in zip(bars, bar_patterns):
    if pat:
        bar.set_hatch(pat)
        bar.set_alpha(0.4)

ax.axvline(0, color="k", lw=1)
ax.set_yticks(y_pos)
ax.set_yticklabels(all_targets, fontsize=8)
ax.set_xlabel(f"Δ Objective ({OBJ}) [%]\n"
              "(positive = improved production after KO)",
              fontsize=9)
ax.set_title(
    "Figure 7. Reaction KO In Silico Screening\n"
    f"(High producer avg constraints | obj={OBJ})\n"
    "Ref: Lim et al. 2010 Metab Eng | Wilkens et al. 2011 J Biotechnol",
    fontsize=10, fontweight="bold"
)
ax.xaxis.grid(True, ls=":", color="0.88"); ax.set_axisbelow(True)
ax.spines[["top","right"]].set_visible(False)

# 값 라벨
max_abs = max(abs(v) for v in bar_vals if v) if any(bar_vals) else 1
for i, (row_data, v) in enumerate(zip(ko_df.itertuples(), bar_vals)):
    status = row_data.status
    if status == "not_found":
        ax.text(0.01 * max_abs, i, "not in model",
                va="center", ha="left", fontsize=7, color="#999999",
                style="italic")
    elif v != 0:
        x = v + max_abs*0.02 if v >= 0 else v - max_abs*0.02
        ax.text(x, i, f"{v:+.1f}%", va="center",
                ha="left" if v >= 0 else "right",
                fontsize=8, fontweight="bold" if abs(v) > 1 else "normal",
                color=P["impr"] if v > 0.5 else
                      P["danger"] if v < -2 else "0.4")
    else:
        ax.text(0.01 * max_abs, i, "neutral",
                va="center", ha="left", fontsize=7, color="0.5")

# 범례
patches = [
    mpatches.Patch(color=P["impr"],    label="Improved (>0.5%)"),
    mpatches.Patch(color=P["neutral"], label="Neutral (−2~+0.5%)"),
    mpatches.Patch(color=P["danger"],  label="Reduced (<−2%)"),
    mpatches.Patch(fc="#E8E8E8", ec="k", hatch="////",
                   label="Not in model"),
]
ax.legend(handles=patches, loc="lower right", fontsize=8)
plt.tight_layout()
save_figure(fig, "Fig7_ko_screen.png", DATASET)

# ════════════════════════════════════════════════════════
# 상세 FLUX MAP
# ════════════════════════════════════════════════════════
print("\n[2] 상세 Flux Pathway Map 생성")
print("    Ref: Orth et al. (2010) Nat Biotechnol 28:245")
print("    Ref: Gopalakrishnan et al. (2024) Metab Eng 82:110")
print("    Ref: Mulukutla et al. (2012) Trends Biotechnol 30:616")

def get_f(fd, kws, default=0.0):
    vals = []
    for rxn_id, v in fd.items():
        r = model.reactions.get_by_id(rxn_id)
        s = f"{rxn_id} {r.name or ''} {getattr(r,'subsystem','') or ''}".lower()
        if any(k.lower() in s for k in kws):
            vals.append(v)
    return np.mean(vals) if vals else default

# ── 노드 정의 (x, y, label, color) ───────────────────
# 좌표계: x=0~16, y=0~18
NODES = {
    # 세포외
    "Glc_ex":   (1.0, 17.0, "Glucose\n(ext)",   "#AED6F1"),
    "Lac_ex":   (5.5, 17.0, "Lactate\n(ext)",   "#F1948A"),
    "Gln_ex":   (13.5, 11.0, "Glutamine\n(ext)", "#A9DFBF"),
    "NH4_ex":   (13.5, 7.5,  "NH4+\n(ext)",      "#F9E79F"),
    "Ala_ex":   (13.5, 5.0,  "Alanine\n(ext)",   "#D7BDE2"),
    "Glu_ex":   (13.5, 9.0,  "Glutamate\n(ext)", "#FDEBD0"),
    # Glycolysis
    "Glc":      (3.0, 17.0,  "Glucose",          "#85C1E9"),
    "G6P":      (3.0, 15.5,  "G-6-P",            "#85C1E9"),
    "F6P":      (3.0, 14.0,  "F-6-P",            "#85C1E9"),
    "FBP":      (3.0, 12.5,  "F-1,6-BP",         "#85C1E9"),
    "DHAP":     (1.5, 11.0,  "DHAP",             "#85C1E9"),
    "GAP":      (3.0, 11.0,  "GAP",              "#85C1E9"),
    "3PG":      (3.0, 9.5,   "3-PG",             "#85C1E9"),
    "2PG":      (3.0, 8.5,   "2-PG",             "#85C1E9"),
    "PEP":      (3.0, 7.5,   "PEP",              "#85C1E9"),
    "Pyr":      (3.0, 6.0,   "Pyruvate",         "#85C1E9"),
    # LDH branch
    "Lac":      (5.0, 6.0,   "Lactate",          "#F1948A"),
    # PPP
    "6PGL":     (1.0, 14.5,  "6-PGL\n(PPP)",     "#FAD7A0"),
    "R5P":      (1.0, 13.0,  "R-5-P\n(PPP)",     "#FAD7A0"),
    # TCA
    "AcCoA":    (3.0, 4.5,   "Acetyl-CoA",       "#A9DFBF"),
    "OAA":      (5.0, 4.5,   "OAA",              "#A9DFBF"),
    "Cit":      (6.5, 5.5,   "Citrate",          "#A9DFBF"),
    "IsoCit":   (8.0, 5.5,   "Isocitrate",       "#A9DFBF"),
    "aKG":      (9.0, 4.5,   "α-KG",             "#A9DFBF"),
    "SucCoA":   (9.0, 3.0,   "Succinyl-CoA",     "#A9DFBF"),
    "Suc":      (7.5, 2.0,   "Succinate",        "#A9DFBF"),
    "Fum":      (6.0, 2.0,   "Fumarate",         "#A9DFBF"),
    "Mal":      (5.0, 2.0,   "Malate",           "#A9DFBF"),
    # Pyruvate anaplerosis
    "OAA2":     (3.0, 3.0,   "OAA\n(PC)",        "#A9DFBF"),
    # Gln/Glu
    "Gln":      (11.5, 11.0, "Glutamine",        "#A9DFBF"),
    "Glu":      (11.5, 9.0,  "Glutamate",        "#A9DFBF"),
    "Pro":      (11.5, 7.0,  "Proline",          "#D7BDE2"),
    "Ala":      (11.5, 5.0,  "Alanine",          "#D7BDE2"),
    # Asp/Asn
    "Asp":      (7.0, 7.5,   "Aspartate",        "#D7BDE2"),
    "Asn":      (7.0, 9.0,   "Asparagine",       "#D7BDE2"),
    # Ser/Gly (one-carbon)
    "3PG_ser":  (5.0, 9.5,   "3-PG→Ser",         "#D7BDE2"),
    "Ser":      (5.0, 10.5,  "Serine",           "#D7BDE2"),
    # FA synthesis
    "MalCoA":   (1.0, 4.0,   "Malonyl-CoA",      "#F9E79F"),
    "FA":       (1.0, 2.5,   "Fatty acids",      "#F9E79F"),
    # Biomass/IgG
    "IgG":      (6.0, 0.5,   f"IgG/Biomass\n({OBJ[:15]})", "#2C3E50"),
}

# ── 엣지 정의 (from, to, keywords, pathway, label) ──
EDGES = [
    # Uptake/Exchange
    ("Glc_ex", "Glc",     ["EX_glc"],                  "exchange", "Glc uptake"),
    ("Lac",    "Lac_ex",  ["EX_lac"],                   "exchange", "Lac secretion"),
    ("Gln_ex", "Gln",     ["EX_gln"],                   "exchange", "Gln uptake"),
    ("NH4_ex", "NH4_ex",  ["EX_nh4"],                   "exchange", "NH4 secr"),
    ("Glu",    "Glu_ex",  ["EX_glu"],                   "exchange", "Glu"),
    ("Ala",    "Ala_ex",  ["EX_ala"],                   "exchange", "Ala secr"),
    # Glycolysis
    ("Glc",    "G6P",     ["HEX","hexokinase"],          "glycolysis", "HK"),
    ("G6P",    "F6P",     ["PGI","phosphoglucose iso"],  "glycolysis", "PGI"),
    ("F6P",    "FBP",     ["PFK","phosphofructo"],       "glycolysis", "PFK"),
    ("FBP",    "GAP",     ["FBA","aldolase"],            "glycolysis", "FBA"),
    ("FBP",    "DHAP",    ["FBA","aldolase"],            "glycolysis", ""),
    ("DHAP",   "GAP",     ["TPI","triose"],              "glycolysis", "TPI"),
    ("GAP",    "3PG",     ["GAPD","PGK","glyceral"],     "glycolysis", "GAPD"),
    ("3PG",    "2PG",     ["PGM","phosphoglycerate m"],  "glycolysis", "PGM"),
    ("2PG",    "PEP",     ["ENO","enolase"],             "glycolysis", "ENO"),
    ("PEP",    "Pyr",     ["PYK","pyruvate kin"],        "glycolysis", "PK"),
    # LDH
    ("Pyr",    "Lac",     ["LDH","lactate dehydro"],     "ldh",       "LDH"),
    # PPP
    ("G6P",    "6PGL",    ["G6PD","glucose-6-phos d"],   "ppp",       "G6PDH"),
    ("6PGL",   "R5P",     ["GND","6PGL","phosphoglucon"],"ppp",       "→R5P"),
    # TCA entry
    ("Pyr",    "AcCoA",   ["PDH","pyruvate dehydro"],    "tca",       "PDH"),
    ("Pyr",    "OAA2",    ["PC","pyruvate carboxyl"],    "tca",       "PC"),
    ("OAA2",   "OAA",     [],                             "tca",       ""),
    # TCA cycle
    ("AcCoA",  "Cit",     ["CS","citrate syn"],          "tca",       "CS"),
    ("OAA",    "Cit",     ["CS","citrate syn"],          "tca",       ""),
    ("Cit",    "IsoCit",  ["ACO","aconit"],              "tca",       "ACO"),
    ("IsoCit", "aKG",     ["IDH","isocitrate dehydro"],  "tca",       "IDH"),
    ("aKG",    "SucCoA",  ["AKGDH","oxoglutarate"],      "tca",       "AKGDH"),
    ("SucCoA", "Suc",     ["SUCOAS","succinyl"],         "tca",       "SUCS"),
    ("Suc",    "Fum",     ["SUCD","succinate dehydro"],  "tca",       "SDH"),
    ("Fum",    "Mal",     ["FUM","fumarase"],            "tca",       "FUM"),
    ("Mal",    "OAA",     ["MDH","malate dehydro"],      "tca",       "MDH"),
    # Gln/Glu
    ("Gln",    "Glu",     ["GLS","glutaminase"],         "gln",       "GLS"),
    ("Glu",    "aKG",     ["GLUD","glutamate dehydro"],  "gln",       "GLUD"),
    ("Glu",    "Ala",     ["ALAT","alanine amino"],      "gln",       "ALAT"),
    ("Glu",    "Asp",     ["GOT","ASPTA","aspartate am"],"gln",       "GOT"),
    ("Glu",    "Pro",     ["P5CS","proline"],            "gln",       "→Pro"),
    # Asp/Asn
    ("OAA",    "Asp",     ["GOT","aspartate amino"],     "aa",        ""),
    ("Asp",    "Asn",     ["ASNS","asparagine syn"],     "aa",        "ASNS"),
    # Ser
    ("3PG",    "3PG_ser", ["PGCD","PSAT","phosphoser"],  "aa",        "PSAT"),
    ("3PG_ser","Ser",     ["PSP","PSPH","ser syn"],      "aa",        ""),
    # Malic enzyme
    ("Mal",    "Pyr",     ["ME1","ME2","malic enzyme"],  "me",        "ME"),
    # ACLY
    ("Cit",    "AcCoA",   ["ACLY","ATP cit"],            "fa",        "ACLY"),
    # FA
    ("AcCoA",  "MalCoA",  ["ACACA","ACC","acetyl-CoA c"],"fa",        "ACC"),
    ("MalCoA", "FA",      ["FASN","fatty acid syn"],     "fa",        "FASN"),
    # Biomass
    ("IgG",    "IgG",     [OBJ],                         "biomass",   OBJ[:12]),
    ("AcCoA",  "IgG",     [OBJ],                         "biomass",   ""),
    ("OAA",    "IgG",     [OBJ],                         "biomass",   ""),
    ("Glu",    "IgG",     [OBJ],                         "biomass",   ""),
]

PATH_COLORS = {
    "exchange":  "#7F8C8D",
    "glycolysis":"#2980B9",
    "ldh":       "#E74C3C",
    "ppp":       "#F39C12",
    "tca":       "#27AE60",
    "gln":       "#8E44AD",
    "aa":        "#9B59B6",
    "me":        "#16A085",
    "fa":        "#D35400",
    "biomass":   "#2C3E50",
}

def draw_detailed_map(ax, flux_dict, title, base_color):
    ax.set_xlim(-0.5, 15); ax.set_ylim(-0.5, 18.5)
    ax.axis("off"); ax.set_facecolor("white")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=8)

    # 구획 배경
    # Glycolysis zone
    ax.add_patch(FancyBboxPatch(xy=(2.2, 5.5), width=1.5, height=12,
        boxstyle="round,pad=0.2", fc="#EBF5FB", ec="#AED6F1",
        lw=1.5, ls="-", zorder=0, alpha=0.5))
    ax.text(3.0, 18.0, "Glycolysis", ha="center", fontsize=8,
            color="#2980B9", fontweight="bold")

    # TCA zone
    ax.add_patch(Circle((6.5, 3.5), 3.8,
        fc="#EAFAF1", ec="#27AE60", lw=1.5, ls="-", zorder=0, alpha=0.5))
    ax.text(6.5, 0.8, "TCA Cycle", ha="center", fontsize=8,
            color="#27AE60", fontweight="bold")

    # FA zone
    ax.add_patch(FancyBboxPatch(xy=(0.3, 1.5), width=2.0, height=3.5,
        boxstyle="round,pad=0.2", fc="#FEF9E7", ec="#F39C12",
        lw=1, ls="--", zorder=0, alpha=0.5))
    ax.text(1.3, 5.3, "Fatty Acid\nSyn", ha="center", fontsize=7,
            color="#D35400", fontweight="bold")

    # PPP zone
    ax.add_patch(FancyBboxPatch(xy=(0.3, 12.5), width=1.5, height=3.5,
        boxstyle="round,pad=0.2", fc="#FEF5E7", ec="#F39C12",
        lw=1, ls="--", zorder=0, alpha=0.4))
    ax.text(1.0, 16.3, "PPP", ha="center", fontsize=7,
            color="#F39C12", fontweight="bold")

    # AA zone
    ax.add_patch(FancyBboxPatch(xy=(10.5, 4.0), width=3.5, height=8,
        boxstyle="round,pad=0.2", fc="#F5EEF8", ec="#9B59B6",
        lw=1, ls="--", zorder=0, alpha=0.4))
    ax.text(12.5, 12.3, "Amino Acid\nMetabolism", ha="center", fontsize=7,
            color="#8E44AD", fontweight="bold")

    # max flux 계산
    fv_all = []
    for _, _, kws, _, _ in EDGES:
        if kws:
            fv_all.append(abs(get_f(flux_dict, kws)))
    max_fv = np.percentile([v for v in fv_all if v > 1e-6], 95) if any(v > 1e-6 for v in fv_all) else 1

    # 엣지 그리기
    drawn_edges = set()
    for src, dst, kws, pathway, edge_label in EDGES:
        if src not in NODES or dst not in NODES: continue
        if src == dst: continue
        edge_key = (src, dst)
        if edge_key in drawn_edges: continue
        drawn_edges.add(edge_key)

        x1, y1 = NODES[src][:2]
        x2, y2 = NODES[dst][:2]
        fv = abs(get_f(flux_dict, kws)) if kws else 0
        lw = max(0.8, min(10, (fv/max_fv)*9)) if max_fv > 0 else 1
        color = PATH_COLORS.get(pathway, "#888888")
        alpha = min(0.9, max(0.2, fv/max_fv*1.5)) if max_fv > 0 else 0.4

        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
            arrowprops=dict(
                arrowstyle=f"-|>,head_width={max(0.04, lw*0.025)},head_length=0.05",
                color=color, lw=lw, alpha=alpha,
                connectionstyle="arc3,rad=0.05",
            ), zorder=3)

        # 엣지 라벨 + flux 값
        if edge_label and fv > 1e-4:
            mx = (x1+x2)/2; my = (y1+y2)/2
            offset_x = 0.15 if x2 > x1 else -0.15
            ax.text(mx + offset_x, my,
                    f"{edge_label}\n{fv:.3f}",
                    ha="center", va="center", fontsize=5.5, color=color,
                    bbox=dict(boxstyle="round,pad=0.1", fc="white",
                              ec="none", alpha=0.75),
                    zorder=4)
        elif edge_label:
            mx = (x1+x2)/2; my = (y1+y2)/2
            ax.text(mx, my, edge_label, ha="center", va="center",
                    fontsize=5, color=color, alpha=0.6, zorder=4)

    # 노드 그리기
    for name, (nx, ny, label, node_color) in NODES.items():
        is_ext  = "_ex" in name
        size    = 800 if "IgG" in name else (400 if is_ext else 550)
        shape   = "*" if "IgG" in name else ("D" if is_ext else "o")
        ec      = base_color if "IgG" in name else "#555555"
        lw_node = 2.0 if "IgG" in name else 0.8

        ax.scatter(nx, ny, s=size, c=node_color, marker=shape,
                   edgecolors=ec, lw=lw_node, zorder=5)

        # 라벨 위치 자동 조정
        offset_y = 0.5 if ny > 9 else -0.5
        if "IgG" in name: offset_y = -0.7
        ax.text(nx, ny + offset_y, label,
                ha="center", va="center",
                fontsize=6 if "\n" in label else 6.5,
                fontweight="bold" if "IgG" in name else "normal",
                zorder=6,
                bbox=dict(boxstyle="round,pad=0.1", fc="white",
                          ec="none", alpha=0.6) if not is_ext else None)

# ── 3패널 그리기 ──────────────────────────────────────
fig2, axes = plt.subplots(1, 3, figsize=(24, 20))
fig2.patch.set_facecolor("white")

draw_detailed_map(axes[0], flux_high,
                  f"High Producer\n(avg: {', '.join(high_clones)})",
                  P["high"])
draw_detailed_map(axes[1], flux_low,
                  f"Low Producer\n(avg: {', '.join(low_clones)})",
                  P["low"])

# Panel 3: Fold Change Map
axes[2].set_xlim(-0.5, 15); axes[2].set_ylim(-0.5, 18.5)
axes[2].axis("off"); axes[2].set_facecolor("white")
axes[2].set_title("Flux Fold Change\nHigh / Low\n(Color = direction, Thickness = magnitude)",
                   fontsize=12, fontweight="bold")

# 배경 구획 (동일)
axes[2].add_patch(FancyBboxPatch(xy=(2.2, 5.5), width=1.5, height=12,
    boxstyle="round,pad=0.2", fc="#EBF5FB", ec="#AED6F1", lw=1.5, zorder=0, alpha=0.4))
axes[2].add_patch(Circle((6.5, 3.5), 3.8,
    fc="#EAFAF1", ec="#27AE60", lw=1.5, zorder=0, alpha=0.4))

cmap = cm.RdYlGn  # 빨강(Low↑) → 초록(High↑)
norm = Normalize(vmin=0.3, vmax=2.0)

drawn2 = set()
for src, dst, kws, pathway, edge_label in EDGES:
    if src not in NODES or dst not in NODES or src == dst: continue
    if (src,dst) in drawn2: continue
    drawn2.add((src,dst))
    x1,y1 = NODES[src][:2]; x2,y2 = NODES[dst][:2]
    fh = abs(get_f(flux_high, kws)) if kws else 0
    fl = abs(get_f(flux_low,  kws)) if kws else 0
    fc = fh/fl if fl > 1e-6 else (2.0 if fh > 1e-6 else 1.0)
    color = cmap(norm(fc))
    lw = max(0.5, min(10, abs(fc-1)*7 + 1))
    axes[2].annotate("", xy=(x2,y2), xytext=(x1,y1),
        arrowprops=dict(
            arrowstyle=f"-|>,head_width={max(0.04,lw*0.025)},head_length=0.05",
            color=color, lw=lw, alpha=0.8,
            connectionstyle="arc3,rad=0.05",
        ), zorder=3)
    if edge_label and abs(fc-1) > 0.15:
        mx,my = (x1+x2)/2,(y1+y2)/2
        axes[2].text(mx, my, f"{fc:.1f}×", ha="center", va="center",
                     fontsize=5.5, fontweight="bold",
                     color="darkgreen" if fc>1.2 else "darkred",
                     bbox=dict(boxstyle="round,pad=0.1", fc="white",
                               ec="none", alpha=0.7), zorder=4)

for name,(nx,ny,label,nc) in NODES.items():
    axes[2].scatter(nx, ny, s=450, c=nc, marker="o" if "_ex" not in name else "D",
                    edgecolors="#555555", lw=0.8, zorder=5)
    axes[2].text(nx, ny+0.45, label, ha="center", va="center",
                 fontsize=6, zorder=6,
                 bbox=dict(fc="white", ec="none", alpha=0.6, pad=0.1))

# Colorbar
sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig2.colorbar(sm, ax=axes[2], shrink=0.4, pad=0.02, location="bottom")
cbar.set_label("Flux Fold Change (High/Low)\n"
               "Green > 1: Higher in High | Red < 1: Higher in Low",
               fontsize=8)

# 경로 범례
leg_handles = [mpatches.Patch(color=c, label=p)
               for p, c in PATH_COLORS.items()]
fig2.legend(handles=leg_handles, loc="lower center", ncol=5,
            fontsize=8, bbox_to_anchor=(0.5, -0.01),
            title="Metabolic Pathways", title_fontsize=9)

fig2.suptitle(
    "Figure E2. Detailed Central Carbon Metabolism Flux Map\n"
    f"Objective: {OBJ} | "
    f"Day {interval[0]}→{interval[1]} constraints\n"
    "Ref: Orth et al. 2010 Nat Biotechnol | "
    "King et al. 2015 PLOS Comput Biol (Escher) | "
    "Mulukutla et al. 2012 Trends Biotechnol",
    fontsize=12, fontweight="bold", y=1.01
)
plt.tight_layout(rect=[0, 0.04, 1, 1])
save_figure(fig2, "FigE2_detailed_flux_map.png", DATASET)

print(f"\n{'='*65}")
print("  완료!")
print(f"  Fig7: results/{DATASET}/figures/Fig7_ko_screen.png")
print(f"  FigE2: results/{DATASET}/figures/FigE2_detailed_flux_map.png")
print(f"{'='*65}")
