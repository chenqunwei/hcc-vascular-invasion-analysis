"""发表级 Figure 3（MVI+AFP）+ Figure 4（预后）+ Figure 5（免疫微环境）。"""
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- 统一发表级风格（与 _figs_core.py 一致） ----------
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "axes.linewidth": 0.6,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})
C_UP, C_DN, C_NS = "#D85A30", "#378ADD", "#C9C9C9"
SUB_COLORS = {"Low VI": "#1D9E75", "Intermediate VI": "#EF9F27", "High VI": "#D85A30"}
SUB_ORDER = ["Low VI", "Intermediate VI", "High VI"]

print("=" * 60)
print("发表级 Figure 3/4/5 生成")
print("=" * 60)

# ============ 数据 ============
clusters = pd.read_csv("results/tcga_lihc_clusters.csv")
vi_order = {1: "Low VI", 3: "Intermediate VI", 2: "High VI"}
clusters["subtype"] = clusters["cluster"].map(vi_order)
merged = clusters[clusters["mvi"].isin([0, 1])].copy()

clin = pd.read_csv("results/tcga_lihc_clinical.csv")
# 从 cbio JSON 补提 AFP（tcga_lihc_clinical.csv 没存 afp）
import json
data = json.load(open("data/tcga/cbio_patient.json", encoding="utf-8"))
afp_map = {d["patientId"]: d["value"] for d in data if d.get("clinicalAttributeId") == "AFP_AT_PROCUREMENT" and d.get("patientId")}
clin["afp"] = clin["patient_id"].map(afp_map)
clin["afp"] = pd.to_numeric(clin["afp"], errors="coerce")
merged = merged.merge(clin[["patient_id", "age", "sex", "ajcc_stage", "afp"]], on="patient_id", how="left")

imm = pd.read_csv("results/tcga_lihc_immune.csv")
# 排除与 clusters 重复的列（subtype/ssgsea_vi），只保留免疫列
immune_cols = [c for c in imm.columns if c not in ("patient_id", "subtype", "ssgsea_vi")]
merged = merged.merge(imm[["patient_id"] + immune_cols], on="patient_id", how="left")

# ============ Figure 3：A MVI 判别 ROC ============
print("生成 Figure 3 ...")

fig = plt.figure(figsize=(10, 3.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.1], wspace=0.35)

# A ROC 曲线
ax = fig.add_subplot(gs[0, 0])
# 用 sklearn 算 ROC
from sklearn.metrics import roc_curve, auc
mvi_pos = merged[merged["mvi"] == 1].dropna(subset=["ssgsea_vi"])
mvi_neg = merged[merged["mvi"] == 0].dropna(subset=["ssgsea_vi"])
y = pd.concat([mvi_pos["mvi"], mvi_neg["mvi"]]).values
scores = pd.concat([mvi_pos["ssgsea_vi"], mvi_neg["ssgsea_vi"]]).values
fpr, tpr, _ = roc_curve(y, scores)
roc_auc = auc(fpr, tpr)
ax.plot(fpr, tpr, color=C_UP, lw=1.5, label=f"VI score (AUC = {roc_auc:.3f})")
ax.plot([0, 1], [0, 1], "--", color="gray", lw=0.8)
ax.fill_between(fpr, 0, tpr, alpha=0.15, color=C_UP)
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_title("MVI discrimination", fontweight="bold")
ax.legend(loc="lower right", frameon=False)
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
ax.text(-0.12, 1.08, "A", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

# B 三型 MVI 比例梯度
ax = fig.add_subplot(gs[0, 1])
mvi_rates = []
for st in SUB_ORDER:
    s = merged[merged["subtype"] == st]
    mvi_rates.append(s["mvi"].mean() * 100)
bars = ax.bar(SUB_ORDER, mvi_rates, color=[SUB_COLORS[s] for s in SUB_ORDER], width=0.6)
for b_, v in zip(bars, mvi_rates):
    ax.text(b_.get_x() + b_.get_width()/2, v + 1, f"{v:.1f}%", ha="center", fontsize=8, fontweight="bold")
ax.set_ylabel("MVI positive (%)")
ax.set_ylim(0, 55)
ax.set_title("MVI by subtype", fontweight="bold")
ct = pd.crosstab(merged["subtype"], merged["mvi"])
chi2v, pmvi, _, _ = stats.chi2_contingency(ct)
ax.text(0.5, 0.95, f"χ² = {chi2v:.2f}, P = {pmvi:.4f}", transform=ax.transAxes, ha="center", fontsize=7)
ax.text(-0.15, 1.08, "B", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

# C AFP 三型
ax = fig.add_subplot(gs[0, 2])
afp_data = []
for st in SUB_ORDER:
    s = merged[(merged["subtype"] == st) & merged["afp"].notna()]
    afp_data.append(s["afp"].values)
bp = ax.boxplot(afp_data, tick_labels=SUB_ORDER, patch_artist=True, widths=0.6)
for patch, s in zip(bp["boxes"], SUB_ORDER):
    patch.set_facecolor(SUB_COLORS[s]); patch.set_alpha(0.6)
for w in bp["whiskers"]: w.set_color("gray")
for c in bp["caps"]: c.set_color("gray")
for m in bp["medians"]: m.set_color("black")
ax.set_ylabel("AFP (ng/mL)")
ax.set_yscale("symlog")  # AFP 有极值，用 symlog
ax.set_title("AFP by subtype", fontweight="bold")
_, p = stats.kruskal(*[merged[merged["subtype"] == s]["afp"].dropna() for s in SUB_ORDER])
ax.text(0.5, 0.95, f"Kruskal–Wallis P < 0.0001", transform=ax.transAxes, ha="center", fontsize=7)
ax.text(-0.15, 1.08, "C", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

fig.savefig("results/Fig3_MVI_clinical.png")
print("  已保存 results/Fig3_MVI_clinical.png")
plt.close(fig)

# ============ Figure 4：DFS + OS + Cox forest ============
print("生成 Figure 4 ...")

from lifelines import KaplanMeierFitter
import json

merged_for_km = clusters.copy()

fig = plt.figure(figsize=(13.2, 3.6))
gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1.2, 1], wspace=0.4)

# A DFS KM
ax = fig.add_subplot(gs[0, 0])
kmf = KaplanMeierFitter()
for st in SUB_ORDER:
    sub = merged_for_km[(merged_for_km["subtype"] == st) & merged_for_km["dfs_months"].notna()
                        & merged_for_km["recurrence"].notna() & (merged_for_km["dfs_months"] > 0)]
    kmf.fit(sub["dfs_months"], sub["recurrence"], label=f"{st} (n={len(sub)})")
    kmf.plot_survival_function(ax=ax, color=SUB_COLORS[st], lw=1.6, ci_show=False)
# log-rank
from lifelines.statistics import multivariate_logrank_test
d = merged_for_km.dropna(subset=["dfs_months", "recurrence", "subtype"]).copy()
d = d[d["dfs_months"] > 0]
res = multivariate_logrank_test(d["dfs_months"], d["subtype"], d["recurrence"])
ax.text(0.95, 0.95, f"log-rank P = {res.p_value:.3f}", transform=ax.transAxes, ha="right", fontsize=8, fontweight="bold")
ax.set_xlabel("Months")
ax.set_ylabel("RFS probability")
ax.set_title("Recurrence-free survival", fontweight="bold")
ax.legend(frameon=False, loc="lower left", fontsize=7)
ax.set_ylim(0, 1.02)
ax.text(-0.15, 1.08, "A", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

# B OS KM
ax = fig.add_subplot(gs[0, 1])
for st in SUB_ORDER:
    sub = merged_for_km[(merged_for_km["subtype"] == st) & merged_for_km["os_months"].notna()
                        & merged_for_km["death"].notna() & (merged_for_km["os_months"] > 0)]
    kmf.fit(sub["os_months"], sub["death"], label=f"{st} (n={len(sub)})")
    kmf.plot_survival_function(ax=ax, color=SUB_COLORS[st], lw=1.6, ci_show=False)
d = merged_for_km.dropna(subset=["os_months", "death", "subtype"]).copy()
d = d[d["os_months"] > 0]
res = multivariate_logrank_test(d["os_months"], d["subtype"], d["death"])
ax.text(0.95, 0.95, f"log-rank P = {res.p_value:.3f}", transform=ax.transAxes, ha="right", fontsize=8, fontweight="bold")
ax.set_xlabel("Months")
ax.set_ylabel("OS probability")
ax.set_title("Overall survival", fontweight="bold")
ax.legend(frameon=False, loc="lower left", fontsize=7)
ax.set_ylim(0, 1.02)
ax.text(-0.15, 1.08, "B", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

# C Cox forest plot（DFS multivariate）
ax = fig.add_subplot(gs[0, 2])
cox = pd.read_csv("results/cox_dfs.csv")
cox_mv = cox[cox["Analysis"] == "Multivariate"].copy()
factors = cox_mv["Variable"].tolist()[::-1]  # 反向，最上是最重要的变量
hrs = cox_mv["HR"].astype(float).tolist()[::-1]
ci_lo = cox_mv["CI_lower"].astype(float).tolist()[::-1]
ci_hi = cox_mv["CI_upper"].astype(float).tolist()[::-1]
ps = cox_mv["P"].tolist()[::-1]

y = np.arange(len(factors))
for yi, (hr, lo, hi, p) in enumerate(zip(hrs, ci_lo, ci_hi, ps)):
    ax.plot([lo, hi], [yi, yi], color="#333333", lw=1.4)
    ax.plot(hr, yi, "o", color=C_UP, markersize=6)
    # 简化 P 显示
    try:
        pv = float(p)
        if pv < 0.001:
            p_disp = "P < 0.001"
        else:
            p_disp = f"P = {pv:.3f}"
    except Exception:
        p_disp = f"P = {p}"
    ax.text(max(ci_hi) * 1.15, yi, f"HR={hr:.2f} ({lo:.2f}-{hi:.2f})\n{p_disp}", va="center", fontsize=7)
ax.axvline(1, color="gray", ls="--", lw=0.8)
ax.set_yticks(y)
ax.set_yticklabels([f.split(" (")[0].replace("Age (continuous)", "Age")
                          .replace("Sex (Male vs Female)", "Sex")
                          .replace("VI score (per SD)", "VI score")
                          .replace("Advanced stage (III/IV vs I/II)", "Stage III/IV")
                          for f in factors], fontsize=8)
ax.set_xlabel("Hazard Ratio (95% CI)")
ax.set_xscale("log")
ax.set_title("Multivariable Cox (DFS)", fontweight="bold")
ax.set_xlim(0.5, max(ci_hi) * 4)
ax.text(-0.15, 1.08, "C", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

# D 外部验证 KM（GSE14520，中位切分）
ax = fig.add_subplot(gs[0, 3])
ext = pd.read_csv("results/gse14520_vi_recurrence.csv").dropna(
    subset=["ssgsea_vi", "recurr", "recurr_months"])
ext = ext[ext["recurr_months"] > 0]
med = ext["ssgsea_vi"].median()
hi_e = ext[ext["ssgsea_vi"] > med]
lo_e = ext[ext["ssgsea_vi"] <= med]
kmf = KaplanMeierFitter()
for lab, sub, color in [("High VI", hi_e, "#D85A30"), ("Low VI", lo_e, "#1D9E75")]:
    kmf.fit(sub["recurr_months"], sub["recurr"], label=f"{lab} (n={len(sub)})")
    kmf.plot_survival_function(ax=ax, color=color, lw=1.6, ci_show=False)
from lifelines.statistics import logrank_test
lr = logrank_test(lo_e["recurr_months"], hi_e["recurr_months"],
                  event_observed_A=lo_e["recurr"], event_observed_B=hi_e["recurr"])
ax.text(0.95, 0.95, f"log-rank P = {lr.p_value:.4f}", transform=ax.transAxes,
        ha="right", fontsize=8, fontweight="bold")
ax.set_xlabel("Months")
ax.set_ylabel("RFS probability")
ax.set_title("External validation (GSE14520)", fontweight="bold")
ax.legend(frameon=False, loc="lower left", fontsize=7)
ax.set_ylim(0, 1.02)
ax.text(-0.15, 1.08, "D", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

fig.savefig("results/Fig4_survival.png")
print("  已保存 results/Fig4_survival.png")
plt.close(fig)

# ============ Figure 5：免疫微环境 ============
print("生成 Figure 5 ...")

fig = plt.figure(figsize=(10, 3.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1, 1.1], wspace=0.35)

checkpoints = ["LAG3", "HAVCR2", "TIGIT", "CD8A", "PDCD1", "CD274", "CTLA4", "IDO1"]
checkpoints = [c for c in checkpoints if c in merged.columns]  # 仅保留有表达数据的基因（PDCD1/CD274/CTLA4/IDO1 不在 VI 基因 panel 内）
check_data = []
for ck in checkpoints:
    if ck in merged.columns:
        vals = []
        for st in SUB_ORDER:
            v = merged[(merged["subtype"] == st) & merged[ck].notna()][ck]
            vals.append(v.values if len(v) > 0 else np.array([]))
        check_data.append(vals)
    else:
        check_data.append([np.array([])] * 3)

# A 免疫检查点表达热图（三型 × 8 基因）
ax = fig.add_subplot(gs[0, 0])
# 用 z-score 后画热图
means = np.array([[merged[(merged["subtype"] == st) & merged[ck].notna()][ck].mean()
                    if ck in merged.columns and merged[ck].notna().any() else 0
                    for st in SUB_ORDER] for ck in checkpoints])
# log2（防止极端值影响）+ z-score per gene
log_means = np.log2(means + 1)
z_means = (log_means - log_means.mean(axis=1, keepdims=True)) / (log_means.std(axis=1, keepdims=True) + 1e-10)
im = ax.imshow(z_means, cmap="RdBu_r", aspect="auto", vmin=-2, vmax=2)
ax.set_xticks(range(3)); ax.set_xticklabels(SUB_ORDER, rotation=30, ha="right")
ax.set_yticks(range(len(checkpoints))); ax.set_yticklabels(checkpoints)
ax.set_title("Immune checkpoint expression", fontweight="bold")
fig.colorbar(im, ax=ax, shrink=0.7, pad=0.04, label="z-score")
ax.text(-0.12, 1.08, "A", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

# B LAG3/HAVCR2/TIGIT 三型箱线（重点展示免疫逃逸）
ax = fig.add_subplot(gs[0, 1])
key_genes = ["LAG3", "HAVCR2", "TIGIT"]
pos = 0
for ck in key_genes:
    if ck in merged.columns:
        data = []
        for st in SUB_ORDER:
            v = merged[(merged["subtype"] == st) & merged[ck].notna()][ck].values
            data.append(v if len(v) > 0 else np.array([0]))
        for i, st in enumerate(SUB_ORDER):
            offset = (i - 1) * 0.25
            xs = np.full(len(data[i]), pos + offset)
            ax.scatter(xs, data[i], s=4, alpha=0.4, color=SUB_COLORS[st])
            ax.boxplot(data[i], positions=[pos + offset], widths=0.18, patch_artist=True,
                       boxprops=dict(facecolor=SUB_COLORS[st], alpha=0.5),
                       medianprops=dict(color="black"), showfliers=False)
        pos += 1
ax.set_xticks(range(len(key_genes)))
ax.set_xticklabels(key_genes)
ax.set_ylabel("Expression")
ax.set_title("Exhaustion markers", fontweight="bold")
ax.text(-0.15, 1.08, "B", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

# C 免疫检查点差异显著性（柱状图：P 值）
ax = fig.add_subplot(gs[0, 2])
pvals = []
for ck in checkpoints:
    if ck in merged.columns:
        groups = [merged[(merged["subtype"] == st) & merged[ck].notna()][ck].values for st in SUB_ORDER]
        groups = [g for g in groups if len(g) >= 2]
        if len(groups) >= 2:
            _, p = stats.kruskal(*groups)
            pvals.append(p)
        else:
            pvals.append(np.nan)
    else:
        pvals.append(np.nan)
neglogp = [-np.log10(p) if p > 0 else 10 for p in pvals]
colors = ["#D85A30" if p < 0.001 else ("#EF9F27" if p < 0.05 else "#888880") for p in pvals]
bars = ax.bar(checkpoints, neglogp, color=colors)
ax.axhline(-np.log10(0.05), color="gray", ls="--", lw=0.6, label="P=0.05")
ax.axhline(-np.log10(0.001), color="gray", ls=":", lw=0.6, label="P=0.001")
ax.set_ylabel("-log10(P)")
ax.set_title("Checkpoint differences (Kruskal–Wallis)", fontweight="bold")
ax.legend(frameon=False, fontsize=7)
plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
ax.text(-0.15, 1.08, "C", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")

fig.savefig("results/Fig5_immune.png")
print("  已保存 results/Fig5_immune.png")
plt.close(fig)

print("\nFigure 3 + Figure 4 + Figure 5 完成。")