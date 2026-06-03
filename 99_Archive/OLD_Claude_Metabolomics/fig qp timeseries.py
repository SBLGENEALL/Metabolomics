"""
fig_qp_timeseries.py
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Figure A: Specific Productivity (qP) 시계열
Figure B: Growth-Production Trade-off scatter
Figure C: qP vs VCD 구간별 분포

근거 논문:
  Templeton et al. (2013) Biotechnol Bioeng 110:2508
    "High-producing CHO clones show reduced growth rate
     but elevated specific productivity (qP)"
  Gopalakrishnan et al. (2024) Metab Eng 82:110
    "High-producing clones do not necessarily need to grow to
     high cell densities; low cell density clones can achieve
     high titers through high specific antibody productivity"
  Doolan et al. (2010) Biotechnol J 5:1085
    qP [pg/cell/day] = ΔProduct / (avg_VCD × Δt)

사용:
  cd C:\\CHO_METABOLOMICS
  python scripts/steps/fig_qp_timeseries.py --dataset practice_20aa
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
from src.fba_utils import save_figure, save_table

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from scipy import stats

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa",
                    choices=["sowa2020", "practice_20aa", "own_experiment"])
args = parser.parse_args()
DATASET = args.dataset

print("=" * 65)
print(f"  fig_qp_timeseries.py — {DATASET}")
print("  Specific Productivity (qP) + Growth-Production Trade-off")
print("=" * 65)
print()
print("  Ref: Templeton et al. (2013) Biotechnol Bioeng 110:2508")
print("  Ref: Gopalakrishnan et al. (2024) Metab Eng 82:110")
print("  Ref: Doolan et al. (2010) Biotechnol J 5:1085")

# ── 데이터 로딩 ────────────────────────────────────────
pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
if not os.path.exists(pkl_path):
    print(f"\n  !! {pkl_path} 없음 — 01_load_data.py 먼저 실행")
    sys.exit(1)
with open(pkl_path, "rb") as f:
    data = pickle.load(f)

df_raw      = data["df_raw"]
high_clones = data["high_clones"]
low_clones  = data["low_clones"]
sorted_clones = data["sorted_clones"]
clones      = data["clones"]
igG_day14   = data.get("igG_day14", {})

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

# ── qP 계산 ────────────────────────────────────────────
# Ref: Doolan et al. (2010)
# qP [pg/cell/day] = Δ[IgG](mg/L) × Volume(L) × 10^9 pg/mg
#                    / (avg_VCD [cells/mL] × Volume(mL) × Δt [day])
# = Δ[IgG] × 10^6 / (avg_VCD × Δt)   [pg/cell/day]

qP_records = []
days = sorted(df_raw["DAY"].unique())

for clone in clones:
    sub = df_raw[df_raw["Sample ID"] == clone].sort_values("DAY")
    rows = {int(r["DAY"]): r for _, r in sub.iterrows()}

    for i in range(len(days) - 1):
        d1, d2 = days[i], days[i+1]
        if d1 not in rows or d2 not in rows:
            continue
        r1, r2 = rows[d1], rows[d2]

        igG1 = float(r1["IgG"])   if pd.notna(r1.get("IgG"))  else 0
        igG2 = float(r2["IgG"])   if pd.notna(r2.get("IgG"))  else 0
        vcd1 = float(r1["Viable Density"]) if pd.notna(r1.get("Viable Density")) else 0
        vcd2 = float(r2["Viable Density"]) if pd.notna(r2.get("Viable Density")) else 0
        dt   = d2 - d1  # days

        avg_vcd = (vcd1 + vcd2) / 2  # 10^6 cells/mL
        delta_igG = igG2 - igG1       # mg/L

        # qP [pg/cell/day] = ΔIgG[mg/L] × 10^6[pg/mg] / (avg_VCD[10^6 cell/mL] × Δt)
        qP = (delta_igG * 1e6) / (avg_vcd * 1e6 * dt) if avg_vcd > 0 else 0
        # 단위: mg/L / (10^6 cells/mL × day) × 10^6 = pg/cell/day

        qP_records.append({
            "clone":    clone,
            "day_mid":  (d1 + d2) / 2,   # 구간 중간 day
            "day_start": d1,
            "day_end":   d2,
            "qP_pg_cell_day": qP,
            "avg_vcd_M_cells_mL": avg_vcd,
            "delta_igG_mg_L": delta_igG,
            "igG_start": igG1,
            "igG_end":   igG2,
            "group": "High" if clone in high_clones else
                     ("Low" if clone in low_clones else "Mid"),
        })

qP_df = pd.DataFrame(qP_records)
save_table(qP_df, "qP_timeseries.csv", DATASET)
print(f"\n  qP 계산 완료: {len(qP_df)} 레코드")

# ── SGR (specific growth rate) 계산 ───────────────────
# Ref: Mulukutla et al. (2012) Trends Biotechnol 30:616
# μ [day^-1] = ln(VCD2/VCD1) / Δt

sgr_records = []
for clone in clones:
    sub = df_raw[df_raw["Sample ID"] == clone].sort_values("DAY")
    rows = {int(r["DAY"]): r for _, r in sub.iterrows()}
    for i in range(len(days) - 1):
        d1, d2 = days[i], days[i+1]
        if d1 not in rows or d2 not in rows: continue
        vcd1 = float(rows[d1]["Viable Density"]) if pd.notna(rows[d1].get("Viable Density")) else 0
        vcd2 = float(rows[d2]["Viable Density"]) if pd.notna(rows[d2].get("Viable Density")) else 0
        dt   = d2 - d1
        mu   = np.log(vcd2/vcd1)/dt if (vcd1 > 0 and vcd2 > 0) else 0
        sgr_records.append({
            "clone": clone, "day_mid": (d1+d2)/2,
            "day_start": d1, "day_end": d2,
            "mu_per_day": mu,
            "avg_vcd": (vcd1+vcd2)/2,
            "group": "High" if clone in high_clones else
                     ("Low" if clone in low_clones else "Mid"),
        })

sgr_df = pd.DataFrame(sgr_records)

# ════════════════════════════════════════════════════════
# FIGURE A: qP 시계열 (논문 표준 형식)
# ════════════════════════════════════════════════════════
# Ref: Gopalakrishnan et al. (2024) Fig 1 형식
print("\n  Fig A: qP timeseries...")

fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# A1: IgG titer 시계열
ax = axes[0, 0]
for clone in sorted_clones:
    sub = df_raw[df_raw["Sample ID"] == clone].sort_values("DAY")
    lw = 2.5 if clone in high_clones else (1.0 if clone in low_clones else 1.3)
    ls = "-"  if clone in high_clones else (":" if clone in low_clones else "--")
    ax.plot(sub["DAY"], sub["IgG"], color=clone_colors[clone],
            lw=lw, ls=ls, marker="o", ms=5, label=clone, zorder=3)
ax.set_xlabel("Culture Day")
ax.set_ylabel("IgG Titer (mg/L)")
ax.set_title("(A) IgG Titer Time Course")
ax_style(ax)
ax.legend(fontsize=7, ncol=2, loc="upper left")
# 최종 titer 주석
for clone in high_clones + low_clones:
    sub = df_raw[df_raw["Sample ID"] == clone]
    last = sub.loc[sub["DAY"].idxmax()]
    ax.annotate(f"{last['IgG']:.0f}",
                xy=(last["DAY"], last["IgG"]),
                xytext=(4, 0), textcoords="offset points",
                fontsize=7, color=clone_colors[clone], va="center")

# A2: VCD 시계열
ax = axes[0, 1]
for clone in sorted_clones:
    sub = df_raw[df_raw["Sample ID"] == clone].sort_values("DAY")
    lw = 2.5 if clone in high_clones else (1.0 if clone in low_clones else 1.3)
    ls = "-"  if clone in high_clones else (":" if clone in low_clones else "--")
    ax.plot(sub["DAY"], sub["Viable Density"], color=clone_colors[clone],
            lw=lw, ls=ls, marker="o", ms=5, zorder=3)
ax.set_xlabel("Culture Day")
ax.set_ylabel("VCD (×10⁶ cells/mL)")
ax.set_title("(B) Viable Cell Density Time Course")
ax_style(ax)

# A3: qP 시계열 — 핵심 Figure
# Ref: Templeton et al. (2013): qP 높은 클론이 반드시 VCD 높지 않음
ax = axes[1, 0]
for clone in sorted_clones:
    cdf = qP_df[qP_df["clone"] == clone].sort_values("day_mid")
    if len(cdf) == 0: continue
    lw = 2.5 if clone in high_clones else (1.0 if clone in low_clones else 1.3)
    ls = "-"  if clone in high_clones else (":" if clone in low_clones else "--")
    ax.plot(cdf["day_mid"], cdf["qP_pg_cell_day"],
            color=clone_colors[clone], lw=lw, ls=ls,
            marker="s", ms=5, label=clone, zorder=3)
ax.axhline(0, color="k", lw=0.8, ls="--", alpha=0.4)
ax.set_xlabel("Culture Day (interval midpoint)")
ax.set_ylabel("Specific Productivity qP\n(pg IgG/cell/day)")
ax.set_title("(C) Specific Productivity (qP) Time Course\n"
             "[Ref: Doolan et al. 2010; Templeton et al. 2013]")
ax_style(ax)
# 음영: 전환기 표시
ax.axvspan(7, 10, alpha=0.08, color="gold", label="Metabolic shift")
ax.legend(fontsize=7, ncol=2, loc="upper right")

# A4: Specific growth rate (μ) 시계열
ax = axes[1, 1]
for clone in sorted_clones:
    sdf = sgr_df[sgr_df["clone"] == clone].sort_values("day_mid")
    if len(sdf) == 0: continue
    lw = 2.5 if clone in high_clones else (1.0 if clone in low_clones else 1.3)
    ls = "-"  if clone in high_clones else (":" if clone in low_clones else "--")
    ax.plot(sdf["day_mid"], sdf["mu_per_day"],
            color=clone_colors[clone], lw=lw, ls=ls,
            marker="^", ms=5, zorder=3)
ax.axhline(0, color="k", lw=0.8, ls="--", alpha=0.4)
ax.set_xlabel("Culture Day (interval midpoint)")
ax.set_ylabel("Specific Growth Rate μ (day⁻¹)")
ax.set_title("(D) Specific Growth Rate Time Course\n"
             "[Ref: Mulukutla et al. 2012 Trends Biotechnol]")
ax_style(ax)
ax.axvspan(7, 10, alpha=0.08, color="gold")

# 범례
high_p = mpatches.Patch(color=P["high"], label=f"High producer (n={len(high_clones)})")
low_p  = mpatches.Patch(color=P["low"],  label=f"Low producer (n={len(low_clones)})")
mid_p  = mpatches.Patch(color=P["mid"][0], label="Mid")
shade_p = mpatches.Patch(color="gold", alpha=0.3, label="Metabolic shift zone (Day 7-10)")
fig.legend(handles=[high_p, low_p, mid_p, shade_p],
           loc="lower center", ncol=4, fontsize=8.5,
           bbox_to_anchor=(0.5, -0.02), frameon=True)

fig.suptitle("Figure A. Culture Performance Profiles\n"
             "IgG Titer / VCD / Specific Productivity (qP) / Specific Growth Rate (μ)\n"
             "Ref: Templeton et al. 2013 Biotechnol Bioeng | "
             "Gopalakrishnan et al. 2024 Metab Eng",
             fontsize=11, fontweight="bold")
plt.tight_layout(rect=[0, 0.05, 1, 1])
save_figure(fig, "FigA_culture_performance_qP.png", DATASET)

# ════════════════════════════════════════════════════════
# FIGURE B: Growth-Production Trade-off
# ════════════════════════════════════════════════════════
# Ref: Templeton et al. (2013): High producers have LOW growth, HIGH qP
# Ref: Gopalakrishnan et al. (2024): Fig 1 — VCD vs qP scatter
print("  Fig B: Growth-Production trade-off...")

# Day 14 기준 최종 VCD, qP 집계
final_data = []
for clone in clones:
    sub  = df_raw[df_raw["Sample ID"] == clone].sort_values("DAY")
    last = sub.iloc[-1]
    # 전체 배양 평균 qP
    clone_qP = qP_df[qP_df["clone"] == clone]["qP_pg_cell_day"]
    # 후반부 (Day 7 이후) 평균 qP가 더 의미있음
    late_qP = qP_df[(qP_df["clone"] == clone) &
                    (qP_df["day_start"] >= 7)]["qP_pg_cell_day"]
    # 최대 VCD
    max_vcd = sub["Viable Density"].max()
    final_data.append({
        "clone":          clone,
        "igG_day14":      igG_day14.get(clone, 0),
        "final_vcd":      float(last["Viable Density"]),
        "max_vcd":        max_vcd,
        "mean_qP_all":    clone_qP.mean()  if len(clone_qP) > 0 else 0,
        "mean_qP_late":   late_qP.mean()   if len(late_qP) > 0 else 0,
        "group":          "High" if clone in high_clones else
                          ("Low" if clone in low_clones else "Mid"),
    })

final_df = pd.DataFrame(final_data)
save_table(final_df, "growth_production_summary.csv", DATASET)

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# B1: Max VCD vs IgG titer
# Ref: Gopalakrishnan 2024: "high-producing clones do NOT necessarily
#      need to grow to high cell densities"
ax = axes[0]
for _, row in final_df.iterrows():
    ax.scatter(row["max_vcd"], row["igG_day14"],
               c=clone_colors[row["clone"]], s=120,
               edgecolors="k", lw=0.8, zorder=3)
    ax.annotate(row["clone"], (row["max_vcd"], row["igG_day14"]),
                xytext=(4, 3), textcoords="offset points", fontsize=7.5)
if len(final_df) > 2:
    sl, ic, r, p, _ = stats.linregress(final_df["max_vcd"], final_df["igG_day14"])
    xl = np.linspace(final_df["max_vcd"].min(), final_df["max_vcd"].max(), 50)
    ax.plot(xl, sl*xl+ic, "k--", lw=1.2, alpha=0.6,
            label=f"r={r:.2f}, p={p:.3f}")
    ax.legend(fontsize=8)
ax.set_xlabel("Max VCD (×10⁶ cells/mL)")
ax.set_ylabel("IgG Titer Day 14 (mg/L)")
ax.set_title("(A) Max VCD vs IgG Titer\n"
             "★ High producers can achieve high titers\nat lower cell densities")
ax_style(ax)

# B2: Late-phase qP vs IgG titer — 핵심
ax = axes[1]
for _, row in final_df.iterrows():
    ax.scatter(row["mean_qP_late"], row["igG_day14"],
               c=clone_colors[row["clone"]], s=120,
               edgecolors="k", lw=0.8, zorder=3)
    ax.annotate(row["clone"], (row["mean_qP_late"], row["igG_day14"]),
                xytext=(4, 3), textcoords="offset points", fontsize=7.5)
if len(final_df) > 2:
    valid = final_df[final_df["mean_qP_late"].abs() > 0]
    if len(valid) > 2:
        sl, ic, r, p, _ = stats.linregress(valid["mean_qP_late"], valid["igG_day14"])
        xl = np.linspace(valid["mean_qP_late"].min(), valid["mean_qP_late"].max(), 50)
        ax.plot(xl, sl*xl+ic, "k--", lw=1.2, alpha=0.6,
                label=f"r={r:.2f}, p={p:.3f}")
        ax.legend(fontsize=8)
ax.set_xlabel("Mean qP Day 7-14 (pg/cell/day)")
ax.set_ylabel("IgG Titer Day 14 (mg/L)")
ax.set_title("(B) Late-phase qP vs IgG Titer\n"
             "[Ref: Doolan et al. 2010 Biotechnol J]")
ax_style(ax)

# B3: Growth-Production Trade-off 사분면 그래프
# Ref: Templeton et al. (2013): 고생산 클론은 qP 높고 μ 낮음
ax = axes[2]
# Day 7~14 평균 SGR
late_sgr = sgr_df[sgr_df["day_start"] >= 7].groupby("clone")["mu_per_day"].mean()
final_df["late_mu"] = final_df["clone"].map(late_sgr)

for _, row in final_df.iterrows():
    if pd.isna(row.get("late_mu")): continue
    ax.scatter(row["late_mu"], row["mean_qP_late"],
               c=clone_colors[row["clone"]], s=120,
               edgecolors="k", lw=0.8, zorder=3)
    ax.annotate(row["clone"], (row["late_mu"], row["mean_qP_late"]),
                xytext=(4, 3), textcoords="offset points", fontsize=7.5)

# 사분면 구분선
mu_med  = final_df["late_mu"].median()
qP_med  = final_df["mean_qP_late"].median()
ax.axvline(mu_med, color="gray", lw=1, ls="--", alpha=0.5)
ax.axhline(qP_med, color="gray", lw=1, ls="--", alpha=0.5)

# 사분면 라벨
xlim = ax.get_xlim(); ylim = ax.get_ylim()
ax.text(mu_med * 0.5 if mu_med > 0 else mu_med * 1.5,
        qP_med * 1.5 if qP_med > 0 else qP_med * 0.5,
        "High qP\nLow growth\n(IDEAL)", ha="center", fontsize=8,
        color=P["high"], fontweight="bold", alpha=0.7)

ax.set_xlabel("Late-phase μ (day⁻¹, Day 7-14)")
ax.set_ylabel("Late-phase qP (pg/cell/day)")
ax.set_title("(C) Growth-Production Trade-off\n"
             "[Ref: Templeton et al. 2013 Biotechnol Bioeng]")
ax_style(ax)

# 범례
for label, color in [(f"High (n={len(high_clones)})", P["high"]),
                     (f"Low (n={len(low_clones)})", P["low"]),
                     ("Mid", P["mid"][0])]:
    ax.scatter([], [], c=color, s=80, label=label)
ax.legend(fontsize=8, loc="best")

fig.suptitle("Figure B. Growth-Production Trade-off Analysis\n"
             "High producers: LOW growth rate (μ), HIGH specific productivity (qP)\n"
             "Ref: Templeton et al. 2013 | Gopalakrishnan et al. 2024",
             fontsize=11, fontweight="bold")
plt.tight_layout()
save_figure(fig, "FigB_growth_production_tradeoff.png", DATASET)

print(f"\n  OK  qP + trade-off figures 완료")
print(f"  results/{DATASET}/figures/ 에 저장됨")
