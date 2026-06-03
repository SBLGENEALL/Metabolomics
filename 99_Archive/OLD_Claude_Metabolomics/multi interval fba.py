"""
multi_interval_fba.py
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
전체 배양 구간을 나눠서 각 구간별 FBA 실행
→ IgG titer와 상관관계가 가장 높은 구간 탐색

CHO fed-batch 14일 구간:
  Day 0→3   적응기
  Day 3→5   성장 초기
  Day 5→7   성장 후기
  Day 7→10  전환기 ← 보통 가장 중요
  Day 10→12 생산기 초기
  Day 12→14 생산기 후기
  Day 3→7   성장기 전체 (현재 사용)
  Day 7→14  생산기 전체
  Day 0→14  전체 ISR

사용:
  cd C:\\CHO_METABOLOMICS
  python scripts\\steps\\multi_interval_fba.py --dataset practice_20aa
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
from src.fba_utils import load_model, save_figure, save_table

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa",
                    choices=["sowa2020", "practice_20aa", "own_experiment"])
args   = parser.parse_args()
DATASET = args.dataset

print("=" * 65)
print(f"  Multi-Interval FBA Analysis — {DATASET}")
print("=" * 65)

# ── 데이터 로딩 ────────────────────────────────────────
pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
if not os.path.exists(pkl_path):
    print(f"  !! {pkl_path} 없음 — 01_load_data.py 먼저 실행")
    sys.exit(1)

with open(pkl_path, "rb") as f:
    base_data = pickle.load(f)

df_raw      = base_data["df_raw"]
igG_day14   = base_data["igG_day14"]
high_clones = base_data["high_clones"]
low_clones  = base_data["low_clones"]
clones      = base_data["clones"]
MET_COLS    = base_data["MET_COLS"]
COL2EX      = base_data["COL2EX"]
COL2NAME    = base_data["COL2NAME"]

# ── 모델 로딩 ──────────────────────────────────────────
model = load_model(verbose=True)
all_rxn_ids = {r.id for r in model.reactions}

# ── Safe apply_bounds ──────────────────────────────────
def apply_bounds(model, constraints):
    n = 0
    for key, (lb, ub) in constraints.items():
        rxn_id = EXCHANGE_IDS.get(key, key)
        if rxn_id not in all_rxn_ids:
            continue
        rxn = model.reactions.get_by_id(rxn_id)
        lb_f, ub_f = float(lb), float(ub)
        if lb_f > ub_f:
            lb_f, ub_f = ub_f, lb_f
        if lb_f > rxn.upper_bound:
            rxn.upper_bound = ub_f
            rxn.lower_bound = lb_f
        elif ub_f < rxn.lower_bound:
            rxn.lower_bound = lb_f
            rxn.upper_bound = ub_f
        else:
            rxn.lower_bound = lb_f
            rxn.upper_bound = ub_f
        n += 1
    return n

# ── 구간 정의 ──────────────────────────────────────────
days_in_data = sorted(df_raw["DAY"].unique())
print(f"\n  데이터 내 Days: {days_in_data}")

# 연속 구간 자동 생성
INTERVALS = []
for i in range(len(days_in_data) - 1):
    INTERVALS.append((days_in_data[i], days_in_data[i+1]))

# 주요 복합 구간 추가
d = days_in_data
if len(d) >= 3:
    INTERVALS.append((d[0],  d[2]))   # Day 0→6  또는 비슷한 구간
if len(d) >= 4:
    INTERVALS.append((d[1],  d[3]))   # Day 3→8
    INTERVALS.append((d[2],  d[-1]))  # Day 6→14 (생산기 전체)
if len(d) >= 5:
    INTERVALS.append((d[1],  d[-2]))  # Day 3→12
INTERVALS.append((d[0], d[-1]))       # 전체 구간

# 중복 제거
INTERVALS = sorted(list(set(INTERVALS)))
print(f"  분석 구간: {INTERVALS}")

# ── Feed 조성 로딩 ────────────────────────────────────
try:
    xl = pd.ExcelFile(AA20_EXCEL if DATASET == "practice_20aa" else SOWA_MMC1)
    df_feed = pd.read_excel(AA20_EXCEL if DATASET == "practice_20aa" else SOWA_MMC1,
                             sheet_name="Feed_Composition").set_index("Feed")
    has_feed = True
except Exception:
    df_feed  = pd.DataFrame()
    has_feed = False

FEED_COLS = {
    "FunctionMax":  "FunctionMax mL",
    "CellBoost 7A": "CellBoost 7A mL",
    "CellBoost 7B": "CellBoost 7B mL",
}

# ── Rate 계산 함수 ────────────────────────────────────
def compute_rate_interval(clone, d1, d2):
    sub  = df_raw[df_raw["Sample ID"] == clone].sort_values("DAY")
    rows = {int(r["DAY"]): r for _, r in sub.iterrows()}
    if d1 not in rows or d2 not in rows:
        return {}
    r1, r2  = rows[d1], rows[d2]
    vcd1    = float(r1["Viable Density"])
    vcd2    = float(r2["Viable Density"])
    v_end   = float(r2["Culture Volume mL"])
    ivcd    = ((vcd1 + vcd2) / 2) * 1e6 * v_end * (d2 - d1) * 24 * 8e-12
    if ivcd <= 0:
        return {}

    rates = {}
    for col in MET_COLS:
        if col not in df_raw.columns:
            continue
        c1 = float(r1[col]) if pd.notna(r1.get(col)) else 0.0
        c2 = float(r2[col]) if pd.notna(r2.get(col)) else 0.0
        v_avg   = (float(r1.get("Culture Volume mL", v_end)) + v_end) / 2
        delta   = (c2 - c1) * v_avg
        feed_add = 0.0
        if has_feed:
            for fn, fc_col in FEED_COLS.items():
                if fc_col not in df_raw.columns: continue
                fv1 = float(r1.get(fc_col, 0)) if pd.notna(r1.get(fc_col, 0)) else 0.0
                fv2 = float(r2.get(fc_col, 0)) if pd.notna(r2.get(fc_col, 0)) else 0.0
                vol = max(0, fv2 - fv1)
                if fn in df_feed.index and col in df_feed.columns:
                    feed_add += vol * float(df_feed.loc[fn, col])
        q = (delta - feed_add) / 1000 / ivcd
        rates[col] = q
    return rates

def rates_to_cst(rates, buffer=0.20):
    cst = {}
    for col, q in rates.items():
        ex_id = EXCHANGE_IDS.get(COL2NAME.get(col, col),
                  EXCHANGE_IDS.get(col, COL2EX.get(col, "")))
        if not ex_id or ex_id not in all_rxn_ids or q == 0:
            continue
        if q < 0:
            lb, ub = q*(1+buffer), q*(1-buffer)
        else:
            lb, ub = q*(1-buffer), q*(1+buffer)
        cst[ex_id] = (min(lb, ub), max(lb, ub))
    return cst

# ── 구간별 FBA 실행 ───────────────────────────────────
print(f"\n  {'구간':12s} {'클론':12s} {'status':10s} {'obj':>10s}")
print("  " + "─" * 48)

interval_results = {}   # {interval: {clone: obj}}

for (d1, d2) in INTERVALS:
    label = f"D{d1}→D{d2}"
    clone_obj = {}
    n_optimal = 0

    for clone in clones:
        rates = compute_rate_interval(clone, d1, d2)
        if not rates:
            clone_obj[clone] = None
            continue
        cst = rates_to_cst(rates)
        with model:
            apply_bounds(model, cst)
            model.objective = OBJ_RXN
            sol = model.optimize()
            obj = sol.objective_value if sol.status == "optimal" else None
            clone_obj[clone] = obj
            if sol.status == "optimal":
                n_optimal += 1

    interval_results[label] = clone_obj
    # 상관관계 계산
    pairs = [(igG_day14[c], clone_obj[c])
             for c in clones if clone_obj.get(c) is not None]
    if len(pairs) >= 3:
        igG_vals = [p[0] for p in pairs]
        obj_vals = [p[1] for p in pairs]
        r, p = stats.pearsonr(igG_vals, obj_vals)
        sig = "★★" if p < 0.01 else ("★" if p < 0.05 else "  ")
        print(f"  {label:12s} optimal={n_optimal}/{len(clones)}  r={r:+.3f}  p={p:.3f} {sig}")
    else:
        print(f"  {label:12s} optimal={n_optimal}/{len(clones)}  (데이터 부족)")

# ── 결과 정리 ─────────────────────────────────────────
# 구간별 r값 계산
corr_rows = []
for label, clone_obj in interval_results.items():
    pairs = [(igG_day14[c], clone_obj[c])
             for c in clones if clone_obj.get(c) is not None]
    if len(pairs) < 3:
        continue
    igG_v = [p[0] for p in pairs]
    obj_v = [p[1] for p in pairs]
    r, p  = stats.pearsonr(igG_v, obj_v)
    corr_rows.append({
        "interval": label,
        "n_optimal": sum(1 for v in clone_obj.values() if v is not None),
        "pearson_r": r,
        "p_value":   p,
        "significant": p < 0.05,
    })

corr_df = pd.DataFrame(corr_rows).sort_values("pearson_r", ascending=False)
print(f"\n  [상관관계 순위 — IgG vs FBA objective]")
print(corr_df.to_string(index=False))

best_interval = corr_df.iloc[0]["interval"] if len(corr_df) > 0 else "D7→D10"
print(f"\n  ★ 최적 구간: {best_interval}  (r={corr_df.iloc[0]['pearson_r']:+.3f})")

# ── Figure 생성 ────────────────────────────────────────
print(f"\n  Figures 생성 중...")

P = PALETTE
clone_colors = {}
for c in clones:
    if c in high_clones:   clone_colors[c] = P["high"]
    elif c in low_clones:  clone_colors[c] = P["low"]
    else:                  clone_colors[c] = P["mid"][clones.index(c) % len(P["mid"])]

# ── FIG A: 구간별 Pearson r 비교 ──────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

ax = axes[0]
colors_bar = [P["high"] if r > 0 else P["low"] for r in corr_df["pearson_r"]]
bars = ax.barh(corr_df["interval"], corr_df["pearson_r"],
               color=colors_bar, edgecolor="k", lw=0.5, zorder=3)
ax.axvline(0, color="k", lw=0.8)
# 유의한 구간 표시
for i, (_, row) in enumerate(corr_df.iterrows()):
    if row["p_value"] < 0.05:
        ax.text(row["pearson_r"] + 0.01 if row["pearson_r"] >= 0
                else row["pearson_r"] - 0.01,
                i, "★ p<0.05", va="center",
                ha="left" if row["pearson_r"] >= 0 else "right",
                fontsize=8, color=P["high"] if row["pearson_r"] > 0 else P["low"])
ax.set_xlabel("Pearson r (IgG Day14 vs FBA Objective)")
ax.set_title("(A) 구간별 IgG–FBA 상관계수\n(양수=정상관: FBA↑→IgG↑, 음수=역상관)")
ax.xaxis.grid(True, ls=":", color="0.88"); ax.set_axisbelow(True)
ax.set_xlim(-1.1, 1.1)

# ── FIG B: 최적 구간 scatter ──────────────────────────
ax = axes[1]
best_obj = interval_results.get(best_interval, {})
plot_pts  = [(c, igG_day14[c], best_obj[c])
             for c in clones if best_obj.get(c) is not None]
if plot_pts:
    xs = [p[1] for p in plot_pts]
    ys = [p[2] for p in plot_pts]
    ax.scatter(xs, ys,
               c=[clone_colors[p[0]] for p in plot_pts],
               s=90, edgecolors="k", lw=0.5, zorder=3)
    if len(plot_pts) > 2:
        sl, ic, r, p_val, _ = stats.linregress(xs, ys)
        xl = np.linspace(min(xs), max(xs), 50)
        ax.plot(xl, sl*xl + ic, "k--", lw=1.2,
                label=f"r={r:.2f}, p={p_val:.3f}")
        ax.legend(fontsize=9)
    for c, x, y in plot_pts:
        ax.annotate(c, (x, y), xytext=(4, 2),
                    textcoords="offset points", fontsize=7.5)
ax.set_xlabel("IgG Day 14 (mg/L)")
ax.set_ylabel(f"FBA Objective ({OBJ_RXN})")
ax.set_title(f"(B) 최적 구간 {best_interval}\nIgG vs FBA Objective")
ax.xaxis.grid(True, ls=":", color="0.88")
ax.yaxis.grid(True, ls=":", color="0.88")
ax.set_axisbelow(True)

fig.suptitle(f"Multi-Interval FBA Analysis — {DATASET}\n"
             f"각 구간 rate 기반 FBA → IgG 상관관계 탐색",
             fontsize=12, fontweight="bold")
plt.tight_layout()
save_figure(fig, "FigX_multi_interval_correlation.png", DATASET)

# ── FIG C: 구간별 × 클론별 heatmap ───────────────────
obj_matrix = pd.DataFrame(
    {label: {c: interval_results[label].get(c) for c in clones}
     for label in interval_results}
)
obj_matrix = obj_matrix.loc[
    sorted(obj_matrix.index, key=lambda c: igG_day14.get(c, 0), reverse=True)
]

import seaborn as sns
fig2, ax2 = plt.subplots(figsize=(max(10, len(INTERVALS)*0.9), 5))
sns.heatmap(obj_matrix.astype(float), ax=ax2,
            cmap="RdYlGn", annot=True, fmt=".4f",
            linewidths=0.4, linecolor="white",
            cbar_kws={"label": f"FBA Objective ({OBJ_RXN})", "shrink": 0.7},
            annot_kws={"size": 7.5})
ax2.set_title(f"Multi-Interval FBA Objective Heatmap\n"
              f"(Rows=clones sorted by IgG, Cols=time intervals)\n"
              f"Green=high objective, Red=low objective")
ax2.set_xlabel("Time Interval")
ax2.set_ylabel("Clone (sorted by IgG titer)")
ax2.tick_params(axis="x", rotation=35)
ax2.tick_params(axis="y", rotation=0)
plt.tight_layout()
save_figure(fig2, "FigX_multi_interval_heatmap.png", DATASET)

# ── FIG D: 구간별 클론 obj 선 그래프 ─────────────────
# 연속 구간만 (복합 구간 제외)
simple_intervals = [(d1, d2) for (d1, d2) in INTERVALS
                    if d2 - d1 <= 3]   # 짧은 구간만
if len(simple_intervals) >= 3:
    fig3, ax3 = plt.subplots(figsize=(10, 5))
    for clone in clones:
        xs, ys = [], []
        for (d1, d2) in simple_intervals:
            label = f"D{d1}→D{d2}"
            obj   = interval_results.get(label, {}).get(clone)
            if obj is not None:
                xs.append(f"D{d1}→{d2}")
                ys.append(obj)
        if xs:
            lw  = 2.5 if clone in high_clones else (1.0 if clone in low_clones else 1.5)
            ls  = "-"  if clone in high_clones else (":" if clone in low_clones else "--")
            ax3.plot(xs, ys, color=clone_colors[clone],
                     lw=lw, ls=ls, marker="o", ms=5, label=clone, zorder=3)
    ax3.set_xlabel("Time Interval")
    ax3.set_ylabel(f"FBA Objective ({OBJ_RXN})")
    ax3.set_title("FBA Objective Across Time Intervals\n"
                  "(Red=High producer, Blue=Low producer)")
    ax3.xaxis.grid(True, ls=":", color="0.88")
    ax3.yaxis.grid(True, ls=":", color="0.88")
    ax3.set_axisbelow(True)
    ax3.legend(fontsize=7.5, ncol=2, loc="upper right")
    plt.tight_layout()
    save_figure(fig3, "FigX_multi_interval_trend.png", DATASET)

# ── 저장 ──────────────────────────────────────────────
save_table(corr_df, "multi_interval_correlation.csv", DATASET)

obj_matrix_save = obj_matrix.copy()
obj_matrix_save.index.name = "clone"
save_table(obj_matrix_save.reset_index(), "multi_interval_obj_matrix.csv", DATASET)

print(f"\n{'='*65}")
print(f"  완료!")
print(f"\n  권장 구간: {best_interval}  (r={corr_df.iloc[0]['pearson_r']:+.3f})")
print(f"\n  다음 액션:")
print(f"  → 01_load_data.py --rate_days {best_interval.replace('D','').replace('→',',')}")
print(f"     위 구간으로 전체 파이프라인 재실행")
print(f"\n  생성 Figure:")
for f in ["FigX_multi_interval_correlation.png",
          "FigX_multi_interval_heatmap.png",
          "FigX_multi_interval_trend.png"]:
    p = os.path.join(ROOT, "results", DATASET, "figures", f)
    if os.path.exists(p):
        print(f"    {f}  ({os.path.getsize(p)//1024} KB)")
