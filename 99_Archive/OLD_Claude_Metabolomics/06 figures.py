"""
06_figures.py — 논문 수준 Figure 생성 (8종 300 DPI)

Fig 1. IgG titer 시계열
Fig 2. Exchange rate 히트맵 (클론 × 대사물질)
Fig 3. Lac/Glc ratio 비교 + IgG 상관관계
Fig 4. FBA objective 클론별 + IgG 상관관계
Fig 5. High vs Low Exchange rate 비교
Fig 6. FVA Forest plot
Fig 7. Reaction KO 스크리닝
Fig 8. 종합 Summary Panel

사용:
  python scripts/steps/06_figures.py --dataset practice_20aa
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
from src.fba_utils import save_figure

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats

parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="practice_20aa",
                    choices=["sowa2020", "practice_20aa", "own_experiment"])
args = parser.parse_args()
DATASET = args.dataset

print("=" * 60)
print(f"  06_figures.py — {DATASET}")
print("=" * 60)

# ── 스타일 설정 ───────────────────────────────────
plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "font.size":         9,
    "axes.labelsize":    10,
    "axes.titlesize":    10,
    "axes.titleweight":  "bold",
    "axes.linewidth":    0.8,
    "xtick.labelsize":   8.5,
    "ytick.labelsize":   8.5,
    "legend.fontsize":   8,
    "legend.framealpha": 0.85,
    "figure.dpi":        150,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
    "lines.linewidth":   1.5,
    "patch.linewidth":   0.6,
})

P = PALETTE

# ── 데이터 로딩 ───────────────────────────────────
pkl_path = os.path.join(DATA_PROCESSED, f"rates_{DATASET}.pkl")
if not os.path.exists(pkl_path):
    print(f"  !! {pkl_path} 없음 — 05_ko_screening.py 먼저 실행")
    sys.exit(1)
with open(pkl_path, "rb") as f:
    data = pickle.load(f)

df_raw          = data["df_raw"]
all_rates       = data["all_rates"]
all_constraints = data["all_constraints"]
igG_day14       = data.get("igG_day14", {})
high_clones     = data["high_clones"]
low_clones      = data["low_clones"]
sorted_clones   = data["sorted_clones"]
clones          = data["clones"]
COL2NAME        = data.get("COL2NAME", {})
MET_COLS        = data.get("MET_COLS", [])
fba_df          = data.get("fba_df", pd.DataFrame())
lac_glc_df      = data.get("lac_glc_df", pd.DataFrame())
fva_results     = data.get("fva_results", {})
ko_df           = data.get("ko_df", pd.DataFrame())
impr_df         = data.get("impr_df", pd.DataFrame())
base_obj        = data.get("base_obj", 0)

# 클론별 색상
clone_colors = {}
for i, c in enumerate(clones):
    if c in high_clones:
        clone_colors[c] = P["high"]
    elif c in low_clones:
        clone_colors[c] = P["low"]
    else:
        clone_colors[c] = P["mid"][i % len(P["mid"])]

def ax_grid(ax):
    ax.xaxis.grid(True, ls=":", color="0.88")
    ax.yaxis.grid(True, ls=":", color="0.88")
    ax.set_axisbelow(True)

# ════════════════════════════════════════════════════
# FIG 1: IgG titer 시계열
# ════════════════════════════════════════════════════
print("  Fig 1: IgG timecourse...")
try:
    days = sorted(df_raw["DAY"].unique())
    fig, ax = plt.subplots(figsize=(9, 5))
    for clone in sorted_clones:
        sub = df_raw[df_raw["Sample ID"] == clone].sort_values("DAY")
        lw  = 2.5 if clone in high_clones else (1.0 if clone in low_clones else 1.5)
        ls  = "-"  if clone in high_clones else (":" if clone in low_clones else "--")
        ms  = 7    if clone in (high_clones + low_clones) else 5
        ax.plot(sub["DAY"], sub["IgG"], color=clone_colors[clone],
                lw=lw, ls=ls, marker="o", ms=ms, label=clone, zorder=3)
    ax.set_xlabel("Culture Day")
    ax.set_ylabel("IgG Titer (mg/L)")
    ax.set_title(f"Figure 1. IgG Production Time Course\n"
                 f"(Red=High producer top {len(high_clones)}, "
                 f"Blue=Low producer bottom {len(low_clones)})")
    ax_grid(ax)
    ax.legend(loc="upper left", ncol=2, fontsize=8)
    # 최종 titer 주석
    for clone in high_clones + low_clones:
        sub = df_raw[df_raw["Sample ID"] == clone]
        if len(sub) == 0: continue
        last = sub.loc[sub["DAY"].idxmax()]
        ax.annotate(f"{last['IgG']:.0f}",
                    xy=(last["DAY"], last["IgG"]),
                    xytext=(4, 0), textcoords="offset points",
                    fontsize=7.5, color=clone_colors[clone], va="center")
    plt.tight_layout()
    save_figure(fig, "Fig1_IgG_timecourse.png", DATASET)
except Exception as e:
    print(f"  !! Fig 1 오류: {e}")

# ════════════════════════════════════════════════════
# FIG 2: Exchange Rate 히트맵
# ════════════════════════════════════════════════════
print("  Fig 2: Exchange rate heatmap...")
try:
    rate_matrix = pd.DataFrame({
        clone: {COL2NAME.get(col, col): info["rate_mmol_gDCWh"]
                for col, info in rates.items()}
        for clone, rates in all_rates.items()
    }).T
    # IgG 기준 정렬
    rate_matrix = rate_matrix.loc[
        sorted(rate_matrix.index, key=lambda c: igG_day14.get(c, 0), reverse=True)
    ]
    fig, ax = plt.subplots(figsize=(max(12, len(MET_COLS) * 0.7), 5))
    sns.heatmap(rate_matrix, ax=ax, cmap="RdBu_r", center=0,
                annot=True, fmt=".3f", linewidths=0.4, linecolor="white",
                cbar_kws={"label": "Rate (mmol/gDCW/h)\nneg=uptake, pos=secretion",
                           "shrink": 0.7},
                annot_kws={"size": 7.5})
    ax.set_title(f"Figure 2. Exchange Rate Heatmap\n"
                 f"(Day {data['day_interval'][0]}→{data['day_interval'][1]}, "
                 f"feed-corrected, sorted by IgG titer)")
    ax.set_xlabel("Metabolite")
    ax.set_ylabel("Clone (sorted by IgG)")
    ax.tick_params(axis="x", rotation=45)
    ax.tick_params(axis="y", rotation=0)
    plt.tight_layout()
    save_figure(fig, "Fig2_rate_heatmap.png", DATASET)
except Exception as e:
    print(f"  !! Fig 2 오류: {e}")

# ════════════════════════════════════════════════════
# FIG 3: Lac/Glc ratio + IgG 상관관계
# ════════════════════════════════════════════════════
print("  Fig 3: Lac/Glc ratio...")
try:
    G = "EX_glc_e"; L = "EX_lac_L_e"
    lgc_rows = []
    for clone, cst in all_constraints.items():
        glc = abs(cst.get(G, (0, 0))[0])
        lac = cst.get(L, (0, 0))[0]
        ratio = lac / glc if (glc > 0 and lac >= 0) else 0
        lgc_rows.append({"clone": clone, "lac_glc": ratio,
                          "igG": igG_day14.get(clone, 0)})
    lgdf = pd.DataFrame(lgc_rows).sort_values("igG", ascending=False)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    ax = axes[0]
    cols_bar = [clone_colors[c] for c in lgdf["clone"]]
    bars = ax.bar(lgdf["clone"], lgdf["lac_glc"], color=cols_bar,
                  edgecolor="k", lw=0.5, zorder=3)
    ax.axhline(1.0, color="k", lw=1, ls="--", label="Ratio = 1.0")
    ax.set_ylabel("Lac/Glc Ratio (mmol/mmol)")
    ax.set_title("(A) Lac/Glc Ratio per Clone\n(<1.0 = efficient TCA)")
    ax.tick_params(axis="x", rotation=35)
    ax_grid(ax); ax.legend(fontsize=8)
    for bar, v in zip(bars, lgdf["lac_glc"]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f"{v:.2f}", ha="center", fontsize=7.5)

    ax = axes[1]
    sc_cols = [clone_colors[c] for c in lgdf["clone"]]
    ax.scatter(lgdf["lac_glc"], lgdf["igG"], c=sc_cols,
               s=80, edgecolors="k", lw=0.5, zorder=3)
    if len(lgdf) > 2:
        sl, ic, r, p, _ = stats.linregress(lgdf["lac_glc"], lgdf["igG"])
        x_line = np.linspace(lgdf["lac_glc"].min(), lgdf["lac_glc"].max(), 50)
        ax.plot(x_line, sl*x_line + ic, "k--", lw=1.2, alpha=0.7,
                label=f"r={r:.2f}, p={p:.3f}")
        ax.legend(fontsize=8)
    for _, row in lgdf.iterrows():
        ax.annotate(row["clone"], (row["lac_glc"], row["igG"]),
                    xytext=(4, 2), textcoords="offset points", fontsize=7.5)
    ax.set_xlabel("Lac/Glc Ratio")
    ax.set_ylabel("IgG Day 14 (mg/L)")
    ax.set_title("(B) Lac/Glc Ratio vs IgG Titer")
    ax_grid(ax)

    fig.suptitle("Figure 3. Lactate/Glucose Ratio — Metabolic Efficiency",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    save_figure(fig, "Fig3_lac_glc_ratio.png", DATASET)
except Exception as e:
    print(f"  !! Fig 3 오류: {e}")

# ════════════════════════════════════════════════════
# FIG 4: FBA objective + IgG 상관관계
# ════════════════════════════════════════════════════
print("  Fig 4: FBA results...")
try:
    if len(fba_df) > 0:
        fba_s = fba_df.sort_values("igG_day14", ascending=False)
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

        ax = axes[0]
        cols4 = [clone_colors.get(c, P["neutral"]) for c in fba_s["clone"]]
        bars = ax.bar(fba_s["clone"], fba_s["obj"], color=cols4,
                      edgecolor="k", lw=0.5, zorder=3)
        ax.set_ylabel(f"FBA Objective ({OBJ_RXN})")
        ax.set_title("(A) FBA Objective per Clone")
        ax.tick_params(axis="x", rotation=35)
        ax_grid(ax)
        for bar, v in zip(bars, fba_s["obj"]):
            if abs(v) > 1e-8:
                ax.text(bar.get_x() + bar.get_width()/2,
                        bar.get_height() * 1.02,
                        f"{v:.4f}", ha="center", fontsize=7)

        ax = axes[1]
        opt = fba_s[(fba_s["status"] == "optimal") & (fba_s["obj"].abs() > 1e-8)]
        if len(opt) > 2:
            ax.scatter(opt["igG_day14"], opt["obj"],
                       c=[clone_colors.get(c, P["neutral"]) for c in opt["clone"]],
                       s=80, edgecolors="k", lw=0.5, zorder=3)
            sl, ic, r, p, _ = stats.linregress(opt["igG_day14"], opt["obj"])
            xl = np.linspace(opt["igG_day14"].min(), opt["igG_day14"].max(), 50)
            ax.plot(xl, sl*xl + ic, "k--", lw=1.2, label=f"r={r:.2f}")
            for _, row in opt.iterrows():
                ax.annotate(row["clone"], (row["igG_day14"], row["obj"]),
                            xytext=(4, 2), textcoords="offset points", fontsize=7.5)
            ax.legend(fontsize=8)
        ax.set_xlabel("IgG Day 14 (mg/L)")
        ax.set_ylabel("FBA Objective")
        ax.set_title("(B) FBA Objective vs IgG Titer")
        ax_grid(ax)

        fig.suptitle(f"Figure 4. FBA Results — {OBJ_RXN}\n"
                     "iCHO3K prod model + feed-corrected constraints",
                     fontsize=11, fontweight="bold")
        plt.tight_layout()
        save_figure(fig, "Fig4_fba_results.png", DATASET)
except Exception as e:
    print(f"  !! Fig 4 오류: {e}")

# ════════════════════════════════════════════════════
# FIG 5: High vs Low Exchange Rate 비교
# ════════════════════════════════════════════════════
print("  Fig 5: High vs Low exchange rates...")
try:
    high_avg = {col: np.mean([all_rates[c][col]["rate_mmol_gDCWh"]
                               for c in high_clones if c in all_rates and col in all_rates[c]])
                for col in MET_COLS}
    low_avg  = {col: np.mean([all_rates[c][col]["rate_mmol_gDCWh"]
                               for c in low_clones  if c in all_rates and col in all_rates[c]])
                for col in MET_COLS}
    diffs = {c: abs(high_avg.get(c, 0) - low_avg.get(c, 0)) for c in MET_COLS}
    sorted_mets = sorted(diffs, key=diffs.get, reverse=True)
    met_labels  = [COL2NAME.get(m, m) for m in sorted_mets]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    y = np.arange(len(sorted_mets)); w = 0.36

    ax = axes[0]
    vh = [high_avg.get(m, 0) for m in sorted_mets]
    vl = [low_avg.get(m, 0)  for m in sorted_mets]
    ax.barh(y + w/2, vh, w, color=P["high"], label=f"High (n={len(high_clones)})",
            edgecolor="k", lw=0.4, zorder=3)
    ax.barh(y - w/2, vl, w, color=P["low"],  label=f"Low (n={len(low_clones)})",
            edgecolor="k", lw=0.4, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(met_labels, fontsize=8.5)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("Exchange Rate (mmol/gDCW/h)\n(neg=uptake, pos=secretion)")
    ax.set_title("(A) High vs Low Producer")
    ax.legend(); ax_grid(ax)

    ax = axes[1]
    delta = [high_avg.get(m, 0) - low_avg.get(m, 0) for m in sorted_mets]
    d_cols = [P["high"] if d > 0 else P["low"] for d in delta]
    ax.barh(y, delta, color=d_cols, edgecolor="k", lw=0.4, zorder=3)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels(met_labels, fontsize=8.5)
    ax.set_xlabel("Delta Rate = High - Low (mmol/gDCW/h)")
    ax.set_title("(B) Difference (High - Low)")
    patches = [mpatches.Patch(color=P["high"], label="Higher in High"),
               mpatches.Patch(color=P["low"],  label="Lower in High")]
    ax.legend(handles=patches, loc="lower right")
    ax_grid(ax)

    fig.suptitle(f"Figure 5. Exchange Rate Comparison: High vs Low Producer\n"
                 f"(High={high_clones}, Low={low_clones})",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    save_figure(fig, "Fig5_high_vs_low_rates.png", DATASET)
except Exception as e:
    print(f"  !! Fig 5 오류: {e}")

# ════════════════════════════════════════════════════
# FIG 6: FVA Forest Plot
# ════════════════════════════════════════════════════
print("  Fig 6: FVA forest plot...")
try:
    if len(fva_results) >= 1:
        KEY_EX = {
            "EX_glc_e": "Glucose",     "EX_lac_L_e": "Lactate",
            "EX_gln_L_e": "Glutamine", "EX_nh4_e":   "NH4+",
            "EX_ala_L_e": "Alanine",   "EX_glu_L_e": "Glutamate",
            "EX_leu_L_e": "Leucine",   "EX_ile_L_e": "Isoleucine",
            "EX_val_L_e": "Valine",    "EX_asp_L_e": "Aspartate",
        }
        fva_cols = {"High": P["high"], "Low": P["low"]}
        fig, ax = plt.subplots(figsize=(10, 6))
        y = np.arange(len(KEY_EX))

        for i, (label, fdf) in enumerate(fva_results.items()):
            color  = fva_cols.get(label, P["mid"][i])
            offset = (i - (len(fva_results) - 1) / 2) * 0.25
            for j, (eid, name) in enumerate(KEY_EX.items()):
                if eid not in fdf.index: continue
                mn  = fdf.loc[eid, "minimum"]
                mx  = fdf.loc[eid, "maximum"]
                mid = fdf.loc[eid, "mean"]
                ax.plot([mn, mx], [y[j] + offset]*2, color=color,
                        lw=4, alpha=0.55, solid_capstyle="round",
                        label=label if j == 0 else "")
                ax.plot(mid, y[j] + offset, "D", color=color, ms=5, zorder=5)
                ax.text(mid, y[j] + offset + 0.18,
                        f"{mid:.3f}", ha="center", va="bottom",
                        fontsize=6.5, color=color)

        ax.set_yticks(y)
        ax.set_yticklabels(list(KEY_EX.values()))
        ax.axvline(0, color="k", lw=0.8, ls="--")
        ax.set_xlabel("Flux (mmol/gDCW/h)")
        ax.set_title("Figure 6. FVA Flux Ranges — High vs Low Producer\n"
                     "(diamond=mean, bar=feasible range at 90% optimality)")
        handles = [mpatches.Patch(color=fva_cols.get(l, P["mid"][i]),
                                   label=f"{l} (n={len(high_clones if l=='High' else low_clones)})")
                   for i, l in enumerate(fva_results)]
        ax.legend(handles=handles, loc="lower right")
        ax_grid(ax)
        plt.tight_layout()
        save_figure(fig, "Fig6_fva_forest.png", DATASET)
    else:
        print("  -- FVA 결과 없음 — 04_run_fva.py 실행 후 다시 시도")
except Exception as e:
    print(f"  !! Fig 6 오류: {e}")

# ════════════════════════════════════════════════════
# FIG 7: KO Screening
# ════════════════════════════════════════════════════
print("  Fig 7: KO screening...")
try:
    if len(ko_df) > 0 and "delta_pct" in ko_df.columns:
        ko_valid = ko_df[ko_df["delta_pct"].notna()].sort_values("delta_pct")
        fig, ax = plt.subplots(figsize=(9, max(4.5, len(ko_valid) * 0.5)))
        cols_k = [P["impr"]   if v > 0.5 else
                  P["neutral"] if v > -2  else P["danger"]
                  for v in ko_valid["delta_pct"]]
        ax.barh(ko_valid["target"], ko_valid["delta_pct"],
                color=cols_k, edgecolor="k", lw=0.4, zorder=3)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel(f"Delta {OBJ_RXN} (%)\n(positive = improved production)")
        ax.set_title(f"Figure 7. Reaction KO In Silico Screening\n"
                     f"(High producer constraints + iCHO3K prod, obj={OBJ_RXN})")
        ax_grid(ax)
        rng = max(abs(ko_valid["delta_pct"].max()),
                  abs(ko_valid["delta_pct"].min()), 0.5)
        for _, row in ko_valid.iterrows():
            v = row["delta_pct"]
            x = v + rng * 0.03 if v >= 0 else v - rng * 0.03
            ax.text(x, row["target"], f"{v:+.1f}%", va="center",
                    ha="left" if v >= 0 else "right", fontsize=8.5,
                    color=P["impr"] if v > 0.5 else
                          P["danger"] if v < -2 else "0.4")
        patches = [mpatches.Patch(color=P["impr"],    label="Improved (>0.5%)"),
                   mpatches.Patch(color=P["neutral"],  label="Neutral"),
                   mpatches.Patch(color=P["danger"],   label="Reduced (<-2%)")]
        ax.legend(handles=patches, loc="lower right")
        plt.tight_layout()
        save_figure(fig, "Fig7_ko_screen.png", DATASET)
except Exception as e:
    print(f"  !! Fig 7 오류: {e}")

# ════════════════════════════════════════════════════
# FIG 8: Summary Panel (논문 main figure)
# ════════════════════════════════════════════════════
print("  Fig 8: Summary panel...")
try:
    fig = plt.figure(figsize=(16, 10))
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.55, wspace=0.42)

    # 8A: IgG timecourse
    ax8a = fig.add_subplot(gs[0, 0])
    for clone in sorted_clones:
        sub = df_raw[df_raw["Sample ID"] == clone].sort_values("DAY")
        lw = 2.2 if clone in high_clones else (1.0 if clone in low_clones else 1.0)
        ls = "-"  if clone in high_clones else (":" if clone in low_clones else "--")
        ax8a.plot(sub["DAY"], sub["IgG"], color=clone_colors[clone],
                  lw=lw, ls=ls, marker="o", ms=4, label=clone, zorder=3)
    ax8a.set_xlabel("Day"); ax8a.set_ylabel("IgG (mg/L)")
    ax8a.set_title("(A) IgG Time Course")
    ax_grid(ax8a); ax8a.legend(fontsize=6, ncol=2)

    # 8B: Lac/Glc ratio
    ax8b = fig.add_subplot(gs[0, 1])
    lgdf2 = pd.DataFrame([
        {"clone": c, "lac_glc": all_constraints[c].get(L, (0,0))[0] /
                                  abs(all_constraints[c].get(G, (-1,))[0])
                                  if abs(all_constraints[c].get(G, (-1,))[0]) > 0 else 0,
         "igG": igG_day14.get(c, 0)}
        for c in sorted_clones if c in all_constraints
    ]).sort_values("igG", ascending=False)
    b8 = ax8b.bar(lgdf2["clone"], lgdf2["lac_glc"],
                   color=[clone_colors.get(c, P["neutral"]) for c in lgdf2["clone"]],
                   edgecolor="k", lw=0.5, zorder=3)
    ax8b.axhline(1.0, color="k", lw=1, ls="--")
    ax8b.set_ylabel("Lac/Glc Ratio")
    ax8b.set_title("(B) Lac/Glc Ratio")
    ax8b.tick_params(axis="x", rotation=35)
    ax_grid(ax8b)
    for bar, v in zip(b8, lgdf2["lac_glc"]):
        ax8b.text(bar.get_x() + bar.get_width()/2, v + 0.01,
                  f"{v:.2f}", ha="center", fontsize=7)

    # 8C: FBA objective
    ax8c = fig.add_subplot(gs[0, 2])
    if len(fba_df) > 0:
        fba_s2 = fba_df.sort_values("igG_day14", ascending=False)
        ax8c.bar(fba_s2["clone"], fba_s2["obj"],
                  color=[clone_colors.get(c, P["neutral"]) for c in fba_s2["clone"]],
                  edgecolor="k", lw=0.5, zorder=3)
    ax8c.set_ylabel("FBA Obj"); ax8c.set_title(f"(C) FBA Obj ({OBJ_RXN[:15]}...)")
    ax8c.tick_params(axis="x", rotation=35); ax_grid(ax8c)

    # 8D: Exchange rate High vs Low (상위 8종)
    ax8d = fig.add_subplot(gs[1, :2])
    top8 = sorted_mets[:8]; top8_labs = [COL2NAME.get(m, m) for m in top8]
    yd = np.arange(len(top8)); wd = 0.32
    vh2 = [high_avg.get(m, 0) for m in top8]
    vl2 = [low_avg.get(m, 0)  for m in top8]
    ax8d.barh(yd + wd/2, vh2, wd, color=P["high"],
               label=f"High (n={len(high_clones)})", edgecolor="k", lw=0.4, zorder=3)
    ax8d.barh(yd - wd/2, vl2, wd, color=P["low"],
               label=f"Low (n={len(low_clones)})", edgecolor="k", lw=0.4, zorder=3)
    ax8d.set_yticks(yd); ax8d.set_yticklabels(top8_labs, fontsize=8)
    ax8d.axvline(0, color="k", lw=0.8)
    ax8d.set_xlabel("Rate (mmol/gDCW/h)")
    ax8d.set_title("(D) Exchange Rates — High vs Low (Top 8)")
    ax8d.legend(fontsize=8); ax_grid(ax8d)

    # 8E: KO screen
    ax8e = fig.add_subplot(gs[1, 2])
    if len(ko_df) > 0 and "delta_pct" in ko_df.columns:
        kv3 = ko_df[ko_df["delta_pct"].notna()].sort_values("delta_pct", ascending=True)
        ck3 = [P["impr"] if v > 0.5 else P["neutral"] if v > -2 else P["danger"]
               for v in kv3["delta_pct"]]
        ax8e.barh(kv3["target"], kv3["delta_pct"], color=ck3, edgecolor="k", lw=0.4, zorder=3)
        ax8e.axvline(0, color="k", lw=0.8)
    ax8e.set_xlabel("Delta Obj (%)"); ax8e.set_title("(E) Reaction KO")
    ax_grid(ax8e)

    day_s, day_e = data.get("day_interval", (3, 7))
    fig.suptitle(f"CHO Fed-Batch FBA/FVA Summary — {DATASET}\n"
                 f"iCHO3K prod model | Day {day_s}→{day_e} rates | "
                 f"High={high_clones} | Low={low_clones}",
                 fontsize=12, fontweight="bold", y=1.01)
    save_figure(fig, "Fig8_summary_panel.png", DATASET)
except Exception as e:
    print(f"  !! Fig 8 오류: {e}")

# ── 결과 목록 ────────────────────────────────────
fig_dir = results_dir(DATASET, "figures")
pngs = sorted(f for f in os.listdir(fig_dir) if f.endswith(".png"))
print(f"\n  [생성된 Figure — results/{DATASET}/figures/]")
for f in pngs:
    sz = os.path.getsize(os.path.join(fig_dir, f)) // 1024
    print(f"    {f:<45s} ({sz:>4d} KB)")

print(f"\n  OK  06_figures 완료 — {len(pngs)}개 Figure 생성")
