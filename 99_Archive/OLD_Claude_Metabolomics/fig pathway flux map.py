"""
fig_pathway_flux_map.py
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Central Carbon Metabolism Flux Map
FBA로 구한 flux를 TCA/Glycolysis 경로 그림 위에
화살표 두께 + 색상으로 시각화

근거 논문:
  Orth et al. (2010) Nature Biotechnology 28:245
    "FBA flux maps visualize metabolic state"
  Goudar et al. (2005) Biotechnol Prog 21:1193
    "Central carbon metabolism flux visualization"
  Mulukutla et al. (2012) Trends Biotechnol 30:616
    "Glycolysis → TCA shift is hallmark of high producers"
  Gopalakrishnan et al. (2024) Metab Eng 82:110
    "Flux maps reveal TCA upregulation in high producers"

사용:
  cd C:\\CHO_METABOLOMICS
  python scripts/steps/fig_pathway_flux_map.py --dataset practice_20aa
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
from src.fba_utils import load_model, save_figure

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch
import pickle

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa",
                    choices=["sowa2020", "practice_20aa", "own_experiment"])
args = parser.parse_args()
DATASET = args.dataset

print("=" * 65)
print(f"  fig_pathway_flux_map.py — {DATASET}")
print("  Central Carbon Metabolism Flux Map")
print("=" * 65)
print()
print("  Ref: Orth et al. (2010) Nat Biotechnol 28:245")
print("  Ref: Mulukutla et al. (2012) Trends Biotechnol 30:616")
print("  Ref: Gopalakrishnan et al. (2024) Metab Eng 82:110")

# ── 데이터 로딩 ────────────────────────────────────────
pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
if not os.path.exists(pkl_path):
    print(f"\n  !! {pkl_path} 없음 — 03_run_fba.py 먼저 실행")
    sys.exit(1)
with open(pkl_path, "rb") as f:
    data = pickle.load(f)

fba_results  = data.get("fba_results", {})
high_clones  = data["high_clones"]
low_clones   = data["low_clones"]
all_constraints = data.get("all_constraints", {})

if not fba_results:
    print("  !! FBA 결과 없음 — 03_run_fba.py 먼저 실행")
    sys.exit(1)

# ── 모델 로딩 + FBA 재실행 ─────────────────────────────
# High / Low 평균 constraints로 FBA → 전체 flux 추출
model = load_model(verbose=True)
all_rxn_ids = {r.id for r in model.reactions}

def apply_bounds_safe(model, constraints):
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
    all_ids = set()
    for c in cst_list: all_ids.update(c.keys())
    return {rid: (np.mean([c[rid][0] for c in cst_list if rid in c]),
                  np.mean([c[rid][1] for c in cst_list if rid in c]))
            for rid in all_ids}

# High/Low FBA 실행 → 전체 flux solution
flux_high = {}
flux_low  = {}
for label, clone_list, flux_dict in [
    ("High", high_clones, flux_high),
    ("Low",  low_clones,  flux_low)
]:
    cst = avg_cst(clone_list)
    with model:
        apply_bounds_safe(model, cst)
        model.objective = OBJ_RXN
        sol = model.optimize()
        if sol.status == "optimal":
            for r in model.reactions:
                flux_dict[r.id] = sol.fluxes[r.id]
            print(f"  {label}: optimal, obj={sol.objective_value:.5f}")
        else:
            print(f"  {label}: {sol.status}")

# ── 핵심 경로 flux 추출 함수 ──────────────────────────
def get_flux(flux_dict, keywords, default=0.0):
    """키워드로 반응을 찾아 flux 반환 (절대값 평균)"""
    matches = []
    for rxn_id, v in flux_dict.items():
        rxn_str = f"{rxn_id} {model.reactions.get_by_id(rxn_id).name or ''}".lower()
        if any(k.lower() in rxn_str for k in keywords):
            matches.append(abs(v))
    return np.mean(matches) if matches else default

# 핵심 경로별 flux 추출
# Ref: Mulukutla 2012, Goudar 2005 — 중요 경로 선정 근거
PATHWAYS = {
    # Glycolysis
    "Glc uptake":      (["EX_glc_e"],           "exchange"),
    "Hexokinase (HK)": (["HEX1","hexokinase"],   "reaction"),
    "PFK":             (["PFK","phosphofructo"],  "reaction"),
    "Pyruvate kinase": (["pyruvate kinase","PKM"],"reaction"),
    # Lactate branch
    "LDH (→Lac)":      (["lactate dehydrogenase","LDH_"], "reaction"),
    "Lac secretion":   (["EX_lac_L_e"],           "exchange"),
    # TCA entry
    "PDH":             (["pyruvate dehydrogenase","PDHm"], "reaction"),
    "Citrate syn":     (["citrate synthase","CS_"],        "reaction"),
    # TCA cycle
    "Isocitrate DH":   (["isocitrate dehydrogenase","IDH"],"reaction"),
    "α-KG DH":         (["oxoglutarate dehydrogenase","AKGDH"],"reaction"),
    "Malate DH":       (["malate dehydrogenase","MDH"],    "reaction"),
    # Glutamine anaplerosis
    "Gln uptake":      (["EX_gln_L_e"],           "exchange"),
    "Glutaminase":     (["glutaminase","GLS_"],    "reaction"),
    "Glutamate DH":    (["glutamate dehydrogenase","GLUD"],"reaction"),
    # NH4+ / byproduct
    "NH4+ secretion":  (["EX_nh4_e"],             "exchange"),
    # PPP
    "G6P DH (PPP)":    (["glucose-6-phosphate dehydrogenase","G6PD"],"reaction"),
    # Fatty acid
    "Fatty acid syn":  (["fatty acid synthase","FASN"],    "reaction"),
    # Biomass/IgG
    "Biomass":         ([OBJ_RXN],                "reaction"),
}

pathway_data = {}
for name, (kws, ptype) in PATHWAYS.items():
    h_flux = get_flux(flux_high, kws) if flux_high else 0
    l_flux = get_flux(flux_low,  kws) if flux_low  else 0
    pathway_data[name] = {"high": h_flux, "low": l_flux,
                          "ratio": h_flux/l_flux if l_flux > 1e-8 else 1.0}

# ════════════════════════════════════════════════════════
# FIGURE C: Flux Map (경로 다이어그램)
# ════════════════════════════════════════════════════════
print("\n  Fig C: Pathway flux map...")

P = PALETTE

def draw_arrow(ax, x1, y1, x2, y2, flux_val, max_flux,
               color, label=None, lw_base=2):
    """
    화살표 두께 = flux 크기에 비례
    Ref: Orth et al. 2010 — flux map convention
    """
    lw = max(0.5, (flux_val / max_flux) * 12) if max_flux > 0 else 1
    ax.annotate("",
                xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(
                    arrowstyle=f"-|>,head_width={lw*0.04},head_length=0.03",
                    color=color, lw=lw,
                    connectionstyle="arc3,rad=0",
                ))
    if label:
        mx, my = (x1+x2)/2, (y1+y2)/2
        ax.text(mx, my, f"{label}\n{flux_val:.3f}",
                ha="center", va="center", fontsize=6.5,
                color=color, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))

# ── 두 패널 (High / Low) ──────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 10))
fig.patch.set_facecolor("white")

def draw_flux_map(ax, flux_dict, title, color_scheme):
    """중심 탄소 대사 flux map 그리기"""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 14)
    ax.axis("off")
    ax.set_facecolor("#FAFAFA")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=12)

    # ── 노드 위치 정의 ──────────────────────────────────
    nodes = {
        # 세포 외부
        "Glucose (ext)":  (1.0, 13.0),
        "Lactate (ext)":  (9.0, 10.5),
        "Gln (ext)":      (9.0, 6.5),
        "NH4+ (ext)":     (9.0, 4.5),
        # 해당과정
        "Glucose (int)":  (2.5, 13.0),
        "G6P":            (2.5, 11.5),
        "F6P/FBP":        (2.5, 10.0),
        "Pyruvate":       (2.5, 7.5),
        # 젖산
        "Lactate (int)":  (4.5, 7.5),
        # TCA cycle
        "Acetyl-CoA":     (2.5, 6.0),
        "Citrate":        (4.5, 5.0),
        "α-KG":           (4.5, 3.5),
        "Succinate":      (4.5, 2.0),
        "Malate":         (2.5, 2.0),
        "OAA":            (2.5, 3.5),
        # Gln 경로
        "Glutamine":      (7.5, 6.5),
        "Glutamate":      (6.0, 5.0),
        # PPP
        "PPP":            (0.8, 11.0),
        # 산물
        "Biomass/IgG":    (2.5, 0.5),
        # FA
        "Fatty acid":     (0.8, 5.5),
    }

    # 노드 그리기
    node_styles = {
        "Glucose (ext)":  {"color":"#AED6F1","size":700,"shape":"s"},
        "Lactate (ext)":  {"color":"#F1948A","size":700,"shape":"s"},
        "Gln (ext)":      {"color":"#A9DFBF","size":700,"shape":"s"},
        "NH4+ (ext)":     {"color":"#F9E79F","size":600,"shape":"s"},
        "Biomass/IgG":    {"color":color_scheme,"size":900,"shape":"*"},
    }
    default_style = {"color":"#ECF0F1","size":500,"shape":"o"}

    for name, (nx, ny) in nodes.items():
        style = node_styles.get(name, default_style)
        ax.scatter(nx, ny, s=style["size"], c=style["color"],
                   marker=style["shape"], edgecolors="k", lw=0.8,
                   zorder=5)
        offset_y = 0.5 if ny < 13 else -0.5
        ax.text(nx, ny + (0.4 if name not in ["Biomass/IgG"] else -0.5),
                name, ha="center", va="center",
                fontsize=6.5, fontweight="bold", zorder=6)

    # ── 화살표 연결 ─────────────────────────────────────
    # max flux 정규화 기준값
    all_fluxes = [abs(v) for v in flux_dict.values() if abs(v) > 1e-6]
    max_f = np.percentile(all_fluxes, 90) if all_fluxes else 1

    def arrow(from_node, to_node, kws, label_text=""):
        x1, y1 = nodes[from_node]
        x2, y2 = nodes[to_node]
        fv = get_flux(flux_dict, kws)
        # 방향 색: 높으면 진한 색, 낮으면 회색
        ratio = fv / max_f if max_f > 0 else 0
        c = color_scheme if ratio > 0.1 else "#AAAAAA"
        draw_arrow(ax, x1, y1, x2, y2, fv, max_f, c,
                   label=f"{fv:.3f}" if fv > 1e-4 else None)

    # Glycolysis
    arrow("Glucose (ext)",  "Glucose (int)", ["EX_glc_e"])
    arrow("Glucose (int)",  "G6P",           ["HEX1","hexokinase"])
    arrow("G6P",            "F6P/FBP",       ["PFK","phosphofructo"])
    arrow("F6P/FBP",        "Pyruvate",      ["pyruvate kinase","PKM",
                                               "enolase","phosphoglycerate"])
    # LDH branch
    arrow("Pyruvate",       "Lactate (int)", ["lactate dehydrogenase","LDH_"])
    arrow("Lactate (int)",  "Lactate (ext)", ["EX_lac_L_e"])
    # TCA entry
    arrow("Pyruvate",       "Acetyl-CoA",    ["pyruvate dehydrogenase","PDHm"])
    arrow("Acetyl-CoA",     "Citrate",       ["citrate synthase","CS_"])
    # TCA cycle
    arrow("Citrate",        "α-KG",          ["isocitrate dehydrogenase","IDH"])
    arrow("α-KG",           "Succinate",     ["oxoglutarate dehydrogenase","AKGDH"])
    arrow("Succinate",      "Malate",        ["succinate dehydrogenase","SUCD"])
    arrow("Malate",         "OAA",           ["malate dehydrogenase","MDH"])
    arrow("OAA",            "Citrate",       ["citrate synthase","CS_"])
    # Gln anaplerosis
    arrow("Gln (ext)",      "Glutamine",     ["EX_gln_L_e"])
    arrow("Glutamine",      "Glutamate",     ["glutaminase","GLS_"])
    arrow("Glutamate",      "α-KG",          ["glutamate dehydrogenase","GLUD"])
    # NH4+
    arrow("Glutamate",      "NH4+ (ext)",    ["EX_nh4_e"])
    # PPP
    arrow("G6P",            "PPP",           ["glucose-6-phosphate dehydrogenase","G6PD"])
    # Biomass
    arrow("OAA",            "Biomass/IgG",   [OBJ_RXN])
    arrow("Malate",         "Biomass/IgG",   [OBJ_RXN])
    # Fatty acid
    arrow("Acetyl-CoA",     "Fatty acid",    ["fatty acid synthase","FASN",
                                               "fatty acid","FA synthesis"])

    # 구획 표시 (세포막)
    cell_border = plt.Polygon(
        [(0.3, 0.1), (9.7, 0.1), (9.7, 12.0), (0.3, 12.0)],
        fill=True, facecolor="#F8F9FA", edgecolor="#BDC3C7",
        lw=1.5, ls="--", zorder=0, alpha=0.5
    )
    ax.add_patch(cell_border)
    ax.text(0.5, 11.7, "Cell", fontsize=9, color="#7F8C8D",
            fontstyle="italic", zorder=4)

    # 경로 영역 배경
    # Glycolysis
    gl_bg = plt.Polygon([(1.8,6.5),(3.2,6.5),(3.2,13.5),(1.8,13.5)],
                         fill=True, facecolor="#EBF5FB", edgecolor="none",
                         zorder=1, alpha=0.6)
    ax.add_patch(gl_bg)
    ax.text(2.5, 12.5, "Glycolysis", ha="center", fontsize=7.5,
            color="#2980B9", fontweight="bold", zorder=4)

    # TCA
    tca_bg = plt.Circle((3.5, 3.5), 2.2,
                          fill=True, facecolor="#EAFAF1",
                          edgecolor="#27AE60", lw=1, ls="--",
                          zorder=1, alpha=0.5)
    ax.add_patch(tca_bg)
    ax.text(3.5, 1.2, "TCA Cycle", ha="center", fontsize=7.5,
            color="#27AE60", fontweight="bold", zorder=4)

for ax_idx, (ax, flux_dict, title, cs) in enumerate([
    (axes[0], flux_high, f"High Producer\n(avg of {high_clones})", P["high"]),
    (axes[1], flux_low,  f"Low Producer\n(avg of {low_clones})", P["low"]),
    (axes[2], {}, "Flux Ratio\n(High / Low)", "#444444"),
]):
    if ax_idx < 2:
        draw_flux_map(ax, flux_dict, title, cs)
    else:
        # 3번째 패널: flux ratio 비교 bar chart
        ax.set_title("Flux Ratio (High / Low)\nby Pathway",
                     fontsize=11, fontweight="bold")
        pnames = list(pathway_data.keys())
        ratios = [pathway_data[n]["ratio"] for n in pnames]
        # 1 기준으로 색상
        bar_cols = [P["high"] if r > 1.1 else
                    (P["low"] if r < 0.9 else P["neutral"])
                    for r in ratios]
        y = np.arange(len(pnames))
        ax.barh(y, [r-1 for r in ratios], color=bar_cols,
                edgecolor="k", lw=0.4, zorder=3)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(pnames, fontsize=8)
        ax.set_xlabel("Flux Ratio − 1\n(positive = higher in High producer)")
        ax.xaxis.grid(True, ls=":", color="0.88"); ax.set_axisbelow(True)
        ax.spines[["top","right"]].set_visible(False)

        # 값 라벨
        for i, (r, b) in enumerate(zip(ratios, bar_cols)):
            val = r - 1
            xpos = val + 0.02 if val >= 0 else val - 0.02
            ax.text(xpos, i, f"{r:.2f}×",
                    va="center", ha="left" if val >= 0 else "right",
                    fontsize=7.5, color=b if b != P["neutral"] else "0.4")

        patches = [
            mpatches.Patch(color=P["high"],    label="Higher in High (>1.1×)"),
            mpatches.Patch(color=P["low"],     label="Higher in Low (<0.9×)"),
            mpatches.Patch(color=P["neutral"], label="Similar (0.9–1.1×)"),
        ]
        ax.legend(handles=patches, loc="lower right", fontsize=7.5)

fig.suptitle("Figure C. Central Carbon Metabolism Flux Map\n"
             "Arrow thickness ∝ flux magnitude | "
             "Ref: Orth et al. 2010 Nat Biotechnol | "
             "Mulukutla et al. 2012 Trends Biotechnol",
             fontsize=11, fontweight="bold", y=1.01)
plt.tight_layout()
save_figure(fig, "FigC_pathway_flux_map.png", DATASET)

# ════════════════════════════════════════════════════════
# FIGURE D: Subsystem별 flux 비교 (논문에서 자주 쓰는 형식)
# Ref: Gopalakrishnan 2024 — subsystem flux comparison
# ════════════════════════════════════════════════════════
print("  Fig D: Subsystem flux comparison...")

subsystems_of_interest = {
    "Glycolysis/Gluconeogenesis": ["glycolysis","gluconeogenesis","HEX","PFK","PKM","GAPD","PGK","ENO"],
    "TCA cycle":                  ["TCA","citrate synthase","isocitrate","oxoglutarate","succinate","malate","fumarate"],
    "Oxidative phosphorylation":  ["NADH","ATP synthase","complex","OXPHOS","cytochrome"],
    "Amino acid metabolism":      ["amino acid","transaminase","aminotransfer","amine"],
    "Fatty acid synthesis":       ["fatty acid","FASN","acetyl-CoA carboxylase","malonyl"],
    "Pentose phosphate":          ["pentose","G6PD","6PGL","GND","RPE","TKT"],
    "Nucleotide synthesis":       ["nucleotide","purine","pyrimidine","IMP","AMP","GMP"],
    "Glutamine metabolism":       ["glutamine","glutamate","glutaminase","GLNS","GLS"],
}

# 서브시스템별 총 flux 계산
subflux_high = {}
subflux_low  = {}
for subname, kws in subsystems_of_interest.items():
    fv_h = sum(abs(v) for rxn_id, v in flux_high.items()
               if any(k.lower() in f"{rxn_id} {model.reactions.get_by_id(rxn_id).name or ''}".lower()
                      for k in kws))
    fv_l = sum(abs(v) for rxn_id, v in flux_low.items()
               if any(k.lower() in f"{rxn_id} {model.reactions.get_by_id(rxn_id).name or ''}".lower()
                      for k in kws))
    subflux_high[subname] = fv_h
    subflux_low[subname]  = fv_l

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

# D1: Grouped bar — subsystem flux 비교
ax = axes[0]
subs = list(subsystems_of_interest.keys())
short_labels = [s.split("/")[0][:18] for s in subs]
x = np.arange(len(subs)); w = 0.35
fh = [subflux_high[s] for s in subs]
fl = [subflux_low[s]  for s in subs]

bars1 = ax.bar(x - w/2, fh, w, color=P["high"], label=f"High (avg {high_clones})",
               edgecolor="k", lw=0.5, zorder=3, alpha=0.85)
bars2 = ax.bar(x + w/2, fl, w, color=P["low"],  label=f"Low (avg {low_clones})",
               edgecolor="k", lw=0.5, zorder=3, alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels(short_labels, rotation=40, ha="right", fontsize=8.5)
ax.set_ylabel("Total Subsystem Flux\n(sum of |flux| in pathway)")
ax.set_title("(A) Subsystem-level Flux Comparison\n"
             "High vs Low Producer")
ax.legend(fontsize=8)
ax.yaxis.grid(True, ls=":", color="0.88"); ax.set_axisbelow(True)
ax.spines[["top","right"]].set_visible(False)

# D2: Fold-change bar (High/Low ratio)
ax = axes[1]
fc = [fh[i]/fl[i] if fl[i] > 1e-6 else 1.0 for i in range(len(subs))]
fc_cols = [P["high"] if f > 1.1 else (P["low"] if f < 0.9 else P["neutral"])
           for f in fc]
bars = ax.bar(x, [f-1 for f in fc], color=fc_cols,
              edgecolor="k", lw=0.5, zorder=3, alpha=0.85)
ax.axhline(0, color="k", lw=0.8)
ax.set_xticks(x)
ax.set_xticklabels(short_labels, rotation=40, ha="right", fontsize=8.5)
ax.set_ylabel("Fold Change − 1  (High / Low)\n(positive = upregulated in High)")
ax.set_title("(B) Metabolic Reprogramming\n"
             "Upregulated vs Downregulated Pathways in High Producer\n"
             "[Ref: Gopalakrishnan et al. 2024 Metab Eng]")
ax.yaxis.grid(True, ls=":", color="0.88"); ax.set_axisbelow(True)
ax.spines[["top","right"]].set_visible(False)

# 값 라벨
for bar, f_val in zip(bars, fc):
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2,
            h + 0.01 if h >= 0 else h - 0.01,
            f"{f_val:.2f}×", ha="center",
            va="bottom" if h >= 0 else "top", fontsize=7.5)

patches = [
    mpatches.Patch(color=P["high"],    label="Upregulated in High (>1.1×)"),
    mpatches.Patch(color=P["low"],     label="Downregulated in High (<0.9×)"),
    mpatches.Patch(color=P["neutral"], label="Similar"),
]
ax.legend(handles=patches, loc="upper right", fontsize=7.5)

fig.suptitle("Figure D. Metabolic Subsystem Flux Comparison — High vs Low Producer\n"
             "Ref: Gopalakrishnan et al. 2024 Metab Eng | "
             "Mulukutla et al. 2012 Trends Biotechnol",
             fontsize=11, fontweight="bold")
plt.tight_layout()
save_figure(fig, "FigD_subsystem_flux_comparison.png", DATASET)

print(f"\n  OK  Flux map figures 완료")
print(f"\n  생성 Figure:")
for fname in ["FigC_pathway_flux_map.png", "FigD_subsystem_flux_comparison.png"]:
    p = os.path.join(ROOT, "results", DATASET, "figures", fname)
    if os.path.exists(p):
        print(f"    {fname}  ({os.path.getsize(p)//1024} KB)")
