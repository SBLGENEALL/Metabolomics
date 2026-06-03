"""
fix_fig3_and_escher.py
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
1. Fig 3 Lac/Glc ratio 수정
   - constraint bounds가 아닌 all_rates 실측값 사용
   - Ref: Zagari et al. (2013) New Biotechnology 30:238

2. Escher 스타일 Flux Pathway Map (HTML + PNG)
   - Escher (King et al. 2015 PLOS Comput Biol) 형식
   - 설치: conda activate cho_fba && pip install escher
   - Central carbon metabolism map on CHO FBA flux

사용:
  pip install escher          # 처음 한 번만
  cd C:\\CHO_METABOLOMICS
  python scripts/steps/fix_fig3_and_escher.py --dataset practice_20aa
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
from scipy import stats

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa",
                    choices=["sowa2020", "practice_20aa", "own_experiment"])
args   = parser.parse_args()
DATASET = args.dataset

print("="*65)
print(f"  fix_fig3_and_escher.py — {DATASET}")
print("="*65)

# ── 데이터 로딩 ────────────────────────────────────────
pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
with open(pkl_path, "rb") as f:
    data = pickle.load(f)

all_rates       = data.get("all_rates", {})
all_constraints = data.get("all_constraints", {})
fba_results     = data.get("fba_results", {})
igG_day14       = data.get("igG_day14", {})
high_clones     = data["high_clones"]
low_clones      = data["low_clones"]
clones          = data["clones"]
sorted_clones   = data.get("sorted_clones", clones)
interval        = data.get("day_interval", (7, 10))
COL2NAME        = data.get("COL2NAME", {})

P = PALETTE
clone_colors = {}
for c in clones:
    if c in high_clones:   clone_colors[c] = P["high"]
    elif c in low_clones:  clone_colors[c] = P["low"]
    else:                  clone_colors[c] = P["mid"][clones.index(c) % len(P["mid"])]

def ax_style(ax):
    ax.xaxis.grid(True, ls=":", color="0.88", zorder=0)
    ax.yaxis.grid(True, ls=":", color="0.88", zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top","right"]].set_visible(False)

# ════════════════════════════════════════════════════════
# FIG 3 수정: all_rates 기반 Lac/Glc ratio
# ════════════════════════════════════════════════════════
# Ref: Zagari et al. (2013) New Biotechnology 30:238
#   "Lac/Glc molar ratio < 1 indicates efficient TCA utilization"
print("\n[1] Fig 3 수정 — all_rates 기반 Lac/Glc ratio")

def get_lac_glc_from_rates(clone):
    """
    constraint bounds가 아닌 실측 rate에서 Lac/Glc 계산
    Ref: Zagari et al. (2013)
    """
    if clone not in all_rates:
        return None
    rates = all_rates[clone]
    # Glucose: 'Gluc' 또는 'Glucose' 컬럼
    glc_col = next((c for c in rates if "gluc" in c.lower()), None)
    lac_col = next((c for c in rates if "lac" in c.lower()), None)
    if not glc_col or not lac_col:
        return None
    glc_rate = rates[glc_col].get("rate_mmol_gDCWh", 0)
    lac_rate  = rates[lac_col].get("rate_mmol_gDCWh", 0)
    if abs(glc_rate) < 1e-9:
        return None
    # Lac/Glc: 양수면 lactate 분비(비효율), 음수면 lactate 소비(효율)
    return lac_rate / abs(glc_rate)

lgc_rows = []
for clone in clones:
    ratio = get_lac_glc_from_rates(clone)
    lgc_rows.append({
        "clone":    clone,
        "lac_glc":  ratio if ratio is not None else 0.0,
        "igG":      igG_day14.get(clone, 0),
        "group":    "High" if clone in high_clones else
                    ("Low" if clone in low_clones else "Mid"),
        "valid":    ratio is not None,
    })
lgc_df = pd.DataFrame(lgc_rows)
save_table(lgc_df, "lac_glc_ratio.csv", DATASET)

print(f"  {'Clone':12s} {'Lac/Glc':>10s} {'IgG':>10s}  해석")
print("  " + "─"*48)
for _, row in lgc_df.sort_values("igG", ascending=False).iterrows():
    interp = "TCA 활성(효율)" if row["lac_glc"] < 0 else \
             ("효율" if row["lac_glc"] < 0.5 else "비효율")
    tag = "★" if row["clone"] in high_clones else \
          ("▼" if row["clone"] in low_clones else " ")
    print(f"  {tag} {row['clone']:12s} {row['lac_glc']:10.3f} "
          f"{row['igG']:10.0f}  {interp}")

# Fig 3 그리기
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# 3A: Bar chart
ax = axes[0]
lgc_s = lgc_df.sort_values("igG", ascending=False)
cols  = [clone_colors[c] for c in lgc_s["clone"]]
bars  = ax.bar(lgc_s["clone"], lgc_s["lac_glc"],
               color=cols, edgecolor="k", lw=0.5, zorder=3)
ax.axhline(0, color="k", lw=0.8, ls="--")
ax.axhline(1.0, color="gray", lw=0.8, ls=":", alpha=0.6, label="Ratio=1.0")
# 음영 영역
ax.axhspan(-0.5, 0, alpha=0.06, color="green")
ax.axhspan(0, 0.5, alpha=0.04, color="orange")
ax.axhspan(0.5, 2.0, alpha=0.04, color="red")
ax.text(0.02, -0.25, "Lactate consumed\n(TCA active)",
        transform=ax.get_yaxis_transform(), fontsize=7.5,
        color="green", ha="left")
ax.text(0.02, 0.25, "Efficient",
        transform=ax.get_yaxis_transform(), fontsize=7.5,
        color="orange", ha="left")
ax.set_ylabel("Lac/Glc Rate Ratio\n(negative = lactate consumed)")
ax.set_title("(A) Lac/Glc Ratio per Clone\n"
             "[Ref: Zagari et al. 2013 New Biotechnol]")
ax.tick_params(axis="x", rotation=35)
ax_style(ax)
ax.legend(fontsize=8)
for bar, v in zip(bars, lgc_s["lac_glc"]):
    ypos = v + 0.02 if v >= 0 else v - 0.05
    ax.text(bar.get_x() + bar.get_width()/2, ypos,
            f"{v:.2f}", ha="center", fontsize=7.5)

# 3B: Scatter + regression
ax = axes[1]
valid = lgc_df[lgc_df["valid"]]
ax.scatter(valid["lac_glc"], valid["igG"],
           c=[clone_colors[c] for c in valid["clone"]],
           s=90, edgecolors="k", lw=0.6, zorder=3)

# x값이 모두 동일한지 확인 후 regression
x_vals = valid["lac_glc"].values
if len(np.unique(x_vals)) > 1:
    sl, ic, r, p, _ = stats.linregress(x_vals, valid["igG"].values)
    xl = np.linspace(x_vals.min(), x_vals.max(), 50)
    ax.plot(xl, sl*xl + ic, "k--", lw=1.2,
            label=f"r={r:.2f}, p={p:.3f}")
    ax.legend(fontsize=8)
else:
    ax.text(0.5, 0.5, "All x values identical\n(rate interval 재검토 필요)",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=9, color="gray")

for _, row in valid.iterrows():
    ax.annotate(row["clone"], (row["lac_glc"], row["igG"]),
                xytext=(4, 2), textcoords="offset points", fontsize=7.5)

ax.axvline(0, color="k", lw=0.8, ls="--", alpha=0.4)
ax.set_xlabel("Lac/Glc Ratio\n(negative = lactate consumed = TCA active)")
ax.set_ylabel("IgG Day 14 (mg/L)")
ax.set_title("(B) Lac/Glc Ratio vs IgG Titer\n"
             "[Ref: Zagari 2013 | Mulukutla 2012]")
ax_style(ax)

fig.suptitle("Figure 3. Lactate/Glucose Ratio — Metabolic Efficiency Indicator\n"
             f"Day {interval[0]}→{interval[1]} rates | "
             "negative ratio = lactate re-consumption (TCA active phase)",
             fontsize=11, fontweight="bold")
plt.tight_layout()
save_figure(fig, "Fig3_lac_glc_ratio.png", DATASET)
print("  Fig 3 저장 완료")

# ════════════════════════════════════════════════════════
# Escher Flux Map
# Ref: King et al. (2015) PLOS Comput Biol 11:e1004321
#   "Escher: A web application for building, sharing,
#    and embedding data-rich visualizations of metabolic pathways"
# ════════════════════════════════════════════════════════
print("\n[2] Escher Flux Map")
print("    Ref: King et al. (2015) PLOS Comput Biol 11:e1004321")

# FBA flux 가져오기
model = load_model(verbose=False)
all_rxn_ids = {r.id for r in model.reactions}

def avg_cst(clone_list):
    cst_list = [all_constraints[c] for c in clone_list if c in all_constraints]
    if not cst_list: return {}
    all_ids = set()
    for c in cst_list: all_ids.update(c.keys())
    return {rid: (np.mean([c[rid][0] for c in cst_list if rid in c]),
                  np.mean([c[rid][1] for c in cst_list if rid in c]))
            for rid in all_ids}

OBJ = data.get("selected_objective", OBJ_RXN)
flux_high, flux_low = {}, {}

for label, clone_list, flux_dict in [("High", high_clones, flux_high),
                                      ("Low",  low_clones,  flux_low)]:
    cst = avg_cst(clone_list)
    with model:
        apply_bounds(model, cst)
        model.objective = OBJ
        sol = model.optimize()
        if sol.status == "optimal":
            for r in model.reactions:
                flux_dict[r.id] = sol.fluxes[r.id]
            print(f"  {label}: optimal obj={sol.objective_value:.5f}")
        else:
            print(f"  {label}: {sol.status}")

# ── Escher HTML 생성 ──────────────────────────────────
try:
    import escher
    ESCHER_AVAILABLE = True
    print("  Escher 패키지 감지됨 → HTML flux map 생성")
except ImportError:
    ESCHER_AVAILABLE = False
    print("  Escher 미설치 → matplotlib flux map으로 대체")
    print("  설치: pip install escher")

if ESCHER_AVAILABLE:
    # Escher HTML 생성
    # Ref: King et al. 2015 — BiGG central metabolism map 사용
    for label, flux_dict in [("High", flux_high), ("Low", flux_low)]:
        if not flux_dict:
            continue
        try:
            b = escher.Builder(
                map_name="e_coli_core.Core metabolism",   # 범용 중심대사 맵
                reaction_data=flux_dict,
                # 색상: 양수=파랑(분비), 음수=빨강(흡수)
                reaction_scale=[
                    {"type": "min",    "color": "#D62728", "size": 20},
                    {"type": "value",  "value": 0, "color": "#AAAAAA", "size": 5},
                    {"type": "max",    "color": "#1F77B4", "size": 20},
                ],
            )
            out_dir  = results_dir(DATASET, "figures")
            html_path = os.path.join(out_dir, f"FigE_escher_flux_{label}.html")
            b.save_html(html_path)
            print(f"  [saved] FigE_escher_flux_{label}.html")
            print(f"    → 브라우저로 열면 인터랙티브 flux map 확인 가능")
        except Exception as e:
            print(f"  Escher 오류: {e}")
            ESCHER_AVAILABLE = False

# Escher 없으면 matplotlib 버전 + 설치 안내
if not ESCHER_AVAILABLE:
    # ── Matplotlib 버전 Escher 스타일 flux map ──────────
    print("\n  Matplotlib Escher 스타일 flux map 생성...")
    print("  (Escher 설치 후 인터랙티브 버전 사용 권장)")

    def get_flux(fd, kws):
        for rxn_id, v in fd.items():
            if any(k.lower() in rxn_id.lower() for k in kws):
                return v
        return 0.0

    # 경로 노드/엣지 정의
    # Ref: Mulukutla et al. (2012) Trends Biotechnol
    NODES = {
        "Glc(ext)":  (1, 9), "Glc-6P":   (2.5, 9),
        "F-1,6BP":   (2.5, 7.5), "GAP":  (2.5, 6),
        "PEP":       (2.5, 4.5), "Pyr":  (2.5, 3),
        "Lac(ext)":  (4.5, 3),  "Lac":   (3.8, 3),
        "AcCoA":     (2.5, 1.5),
        "Cit":       (4.5, 1), "IsoCit": (5.5, 2),
        "aKG":       (5.5, 3.5), "SucCoA":(5.5, 5),
        "Suc":       (4.5, 6), "Fum":    (3.8, 6.5),
        "Mal":       (3.2, 7), "OAA":    (2.5, 7.5),
        "Gln(ext)":  (7, 4), "Gln":      (6.2, 4),
        "Glu":       (6.2, 3.5), "NH4":  (7, 2.5),
        "PPP":       (1.2, 7.5),
        "Biomass":   (2.5, 0.3),
    }

    EDGES = [
        # Glycolysis
        ("Glc(ext)", "Glc-6P",  ["HEX1","hexokinase"],         "Glycolysis"),
        ("Glc-6P",   "F-1,6BP", ["PFK","phosphofructo"],        "Glycolysis"),
        ("F-1,6BP",  "GAP",     ["FBA","aldolase"],             "Glycolysis"),
        ("GAP",      "PEP",     ["GAPD","PGK","ENO"],           "Glycolysis"),
        ("PEP",      "Pyr",     ["PYK","pyruvate kinase"],      "Glycolysis"),
        # LDH
        ("Pyr",      "Lac",     ["LDH","lactate dehydrogen"],   "LDH"),
        ("Lac",      "Lac(ext)",["EX_lac"],                     "Exchange"),
        # PDH→TCA
        ("Pyr",      "AcCoA",   ["PDH","pyruvate dehydrogen"],  "TCA"),
        ("AcCoA",    "Cit",     ["CS","citrate syn"],           "TCA"),
        ("Cit",      "IsoCit",  ["ACO","aconit"],               "TCA"),
        ("IsoCit",   "aKG",     ["IDH","isocitrate"],           "TCA"),
        ("aKG",      "SucCoA",  ["AKGDH","oxoglutarate"],       "TCA"),
        ("SucCoA",   "Suc",     ["SUCOAS","succinyl"],          "TCA"),
        ("Suc",      "Fum",     ["SUCD","succinate dehydro"],   "TCA"),
        ("Fum",      "Mal",     ["FUM","fumarase"],             "TCA"),
        ("Mal",      "OAA",     ["MDH","malate dehydro"],       "TCA"),
        ("OAA",      "Cit",     ["CS","citrate syn"],           "TCA"),
        # Anaplerosis
        ("Gln(ext)", "Gln",     ["EX_gln"],                     "Gln"),
        ("Gln",      "Glu",     ["GLS","glutaminase"],          "Gln"),
        ("Glu",      "aKG",     ["GLUD","glutamate dehydro"],   "Gln"),
        ("Glu",      "NH4",     ["EX_nh4"],                     "Exchange"),
        # Exchange
        ("Glc(ext)", "Glc-6P",  ["EX_glc"],                     "Exchange"),
        # PPP
        ("Glc-6P",   "PPP",     ["G6PD","glucose-6-phosphate"], "PPP"),
        # Biomass
        ("OAA",      "Biomass", [OBJ],                          "Biomass"),
    ]

    PATH_COLORS = {
        "Glycolysis": "#2980B9",
        "TCA":        "#27AE60",
        "LDH":        "#E74C3C",
        "Gln":        "#8E44AD",
        "Exchange":   "#7F8C8D",
        "PPP":        "#F39C12",
        "Biomass":    "#2C3E50",
    }

    fig, axes = plt.subplots(1, 3, figsize=(20, 11))
    fig.patch.set_facecolor("white")

    def draw_map(ax, flux_dict, title, base_color):
        ax.set_xlim(0.5, 8); ax.set_ylim(-0.2, 10.5)
        ax.axis("off"); ax.set_facecolor("#F8F9FA")
        ax.set_title(title, fontsize=12, fontweight="bold", pad=10)

        # 모든 flux 중 최대값 (화살표 두께 정규화)
        all_fv = [abs(get_flux(flux_dict, e[2])) for e in EDGES]
        max_fv = np.percentile([v for v in all_fv if v > 0], 90) if any(v>0 for v in all_fv) else 1

        # 배경 영역
        from matplotlib.patches import FancyBboxPatch
        regions = [
            ((1.8,2.7),(1.2,7.0), "#EBF5FB", "Glycolysis"),
            ((3.5,3.5), 2.8,      "#EAFAF1", "TCA Cycle"),
        ]
        ax.add_patch(FancyBboxPatch((1.8, 2.5), 1.2, 7.2,
            boxstyle="round,pad=0.1", fc="#EBF5FB", ec="none", zorder=0, alpha=0.6))
        from matplotlib.patches import Circle
        ax.add_patch(Circle((4.5, 4), 2.4,
            fc="#EAFAF1", ec="#27AE60", lw=1, ls="--", zorder=0, alpha=0.5))

        ax.text(2.4, 9.8, "Glycolysis", fontsize=8, color="#2980B9",
                fontweight="bold", ha="center")
        ax.text(4.5, 1.3, "TCA Cycle", fontsize=8, color="#27AE60",
                fontweight="bold", ha="center")

        # 엣지 (화살표)
        for src, dst, kws, pathway in EDGES:
            if src not in NODES or dst not in NODES: continue
            x1, y1 = NODES[src]; x2, y2 = NODES[dst]
            fv = abs(get_flux(flux_dict, kws))
            lw = max(0.8, min(12, (fv/max_fv)*10)) if max_fv > 0 else 1
            color = PATH_COLORS.get(pathway, "#888888")
            alpha = 0.7 if fv > 1e-4 else 0.2

            ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(
                    arrowstyle=f"-|>,head_width={max(0.05,lw*0.03)},head_length=0.06",
                    color=color, lw=lw, alpha=alpha,
                    connectionstyle="arc3,rad=0.05",
                ))
            # flux 값 표시 (주요 반응만)
            if fv > 1e-4 and pathway in ["Glycolysis","TCA","LDH"]:
                mx, my = (x1+x2)/2 + 0.08, (y1+y2)/2
                ax.text(mx, my, f"{fv:.3f}", fontsize=5.5,
                        color=color, ha="center", va="center",
                        bbox=dict(fc="white", ec="none", alpha=0.6, pad=0.5))

        # 노드
        for name, (nx, ny) in NODES.items():
            is_ext = "(ext)" in name
            size   = 500 if not is_ext else 350
            color  = "#AED6F1" if "Glc" in name else \
                     "#F1948A" if "Lac" in name else \
                     "#A9DFBF" if "Gln" in name else \
                     "#F9E79F" if "NH4" in name else \
                     base_color if "Biomass" in name else "#ECF0F1"
            shape  = "*" if "Biomass" in name else ("s" if is_ext else "o")
            ax.scatter(nx, ny, s=size, c=color, marker=shape,
                       edgecolors="k", lw=0.8, zorder=5)
            offset = 0.35 if ny < 5 else -0.35
            ax.text(nx, ny + offset, name.replace("(ext)", "\n(ext)"),
                    ha="center", va="center", fontsize=6.5,
                    fontweight="bold" if "Biomass" in name else "normal",
                    zorder=6)

        # 범례 (경로별 색상)
        handles = [mpatches.Patch(color=c, label=p)
                   for p, c in PATH_COLORS.items()]
        ax.legend(handles=handles, loc="upper right", fontsize=6.5,
                  framealpha=0.9, ncol=1)

    # Panel 1: High producer
    draw_map(axes[0], flux_high,
             f"High Producer\n(avg: {', '.join(high_clones)})",
             P["high"])

    # Panel 2: Low producer
    draw_map(axes[1], flux_low,
             f"Low Producer\n(avg: {', '.join(low_clones)})",
             P["low"])

    # Panel 3: Fold change (High/Low)
    axes[2].set_title("Flux Fold Change\n(High / Low)",
                      fontsize=12, fontweight="bold")
    axes[2].set_xlim(0.5, 8); axes[2].set_ylim(-0.2, 10.5)
    axes[2].axis("off"); axes[2].set_facecolor("#F8F9FA")

    for src, dst, kws, pathway in EDGES:
        if src not in NODES or dst not in NODES: continue
        x1, y1 = NODES[src]; x2, y2 = NODES[dst]
        fh = abs(get_flux(flux_high, kws))
        fl = abs(get_flux(flux_low,  kws))
        fc = fh/fl if fl > 1e-6 else (2.0 if fh > 1e-6 else 1.0)
        color = P["high"] if fc > 1.2 else (P["low"] if fc < 0.8 else "#AAAAAA")
        lw    = max(0.5, min(10, abs(fc-1)*8))
        axes[2].annotate("", xy=(x2, y2), xytext=(x1, y1),
            arrowprops=dict(
                arrowstyle=f"-|>,head_width={max(0.04,lw*0.03)},head_length=0.06",
                color=color, lw=max(0.5, lw),
                connectionstyle="arc3,rad=0.05",
            ))
        if abs(fc - 1) > 0.15:
            mx, my = (x1+x2)/2 + 0.08, (y1+y2)/2
            axes[2].text(mx, my, f"{fc:.1f}×", fontsize=6,
                         color=color, ha="center", va="center",
                         fontweight="bold",
                         bbox=dict(fc="white", ec="none", alpha=0.7, pad=0.5))

    for name, (nx, ny) in NODES.items():
        axes[2].scatter(nx, ny, s=400, c="#ECF0F1", marker="o",
                        edgecolors="k", lw=0.6, zorder=5)
        axes[2].text(nx, ny+0.3, name, ha="center", va="center",
                     fontsize=6, zorder=6)

    fold_patches = [
        mpatches.Patch(color=P["high"], label="Higher in High (>1.2×)"),
        mpatches.Patch(color=P["low"],  label="Higher in Low (<0.8×)"),
        mpatches.Patch(color="#AAAAAA", label="Similar (0.8–1.2×)"),
    ]
    axes[2].legend(handles=fold_patches, loc="upper right",
                   fontsize=7, framealpha=0.9)

    fig.suptitle(
        "Figure E. Central Carbon Metabolism Flux Map\n"
        "Arrow thickness ∝ flux magnitude | "
        f"Objective: {OBJ} | Day {interval[0]}→{interval[1]}\n"
        "Ref: King et al. 2015 PLOS Comput Biol (Escher) | "
        "Orth et al. 2010 Nat Biotechnol",
        fontsize=11, fontweight="bold", y=1.01
    )
    plt.tight_layout()
    save_figure(fig, "FigE_flux_map_escher_style.png", DATASET)

# ── Escher 설치 안내 (HTML 버전) ──────────────────────
if not ESCHER_AVAILABLE:
    escher_guide = os.path.join(results_dir(DATASET, "figures"),
                                "escher_setup_guide.txt")
    guide = f"""
Escher Interactive Flux Map 설치 방법
══════════════════════════════════════
Ref: King et al. (2015) PLOS Comput Biol 11:e1004321

1. 설치:
   conda activate cho_fba
   pip install escher

2. 이 스크립트 재실행:
   python scripts/steps/fix_fig3_and_escher.py --dataset {DATASET}
   → FigE_escher_flux_High.html, FigE_escher_flux_Low.html 생성

3. HTML 파일을 브라우저로 열면:
   - 반응을 클릭해서 flux 값 확인
   - 줌/패닝 가능
   - PNG로 내보내기 가능

4. 또는 온라인 Escher:
   https://escher.github.io
   → Load Map → Central Metabolism
   → Load Data → reaction_fluxes.json 업로드

reaction_fluxes.json 경로:
   results/{DATASET}/tables/reaction_fluxes_high.json
   results/{DATASET}/tables/reaction_fluxes_low.json
"""
    with open(escher_guide, "w", encoding="utf-8") as f:
        f.write(guide)
    print(f"\n  Escher 설치 안내: {escher_guide}")

    # flux를 JSON으로 저장 (온라인 Escher에서 사용 가능)
    import json
    for label, fd in [("high", flux_high), ("low", flux_low)]:
        if fd:
            json_path = os.path.join(
                results_dir(DATASET, "tables"),
                f"reaction_fluxes_{label}.json"
            )
            # Escher 형식: {"reaction_id": flux_value}
            escher_data = {k: round(v, 6) for k, v in fd.items()
                           if abs(v) > 1e-6}
            with open(json_path, "w") as f:
                json.dump(escher_data, f, indent=2)
            print(f"  [saved] reaction_fluxes_{label}.json "
                  f"({len(escher_data)}개 반응)")
    print("  → https://escher.github.io 에서 이 JSON을 업로드하면")
    print("    인터랙티브 flux map 사용 가능")

print(f"\n{'='*65}")
print("  완료!")
print(f"  Fig 3: results/{DATASET}/figures/Fig3_lac_glc_ratio.png")
print(f"  Fig E: results/{DATASET}/figures/FigE_flux_map_escher_style.png")
print(f"  JSON : results/{DATASET}/tables/reaction_fluxes_high/low.json")
print(f"{'='*65}")
