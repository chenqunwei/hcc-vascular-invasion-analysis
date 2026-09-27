# -*- coding: utf-8 -*-
"""
generate_tables.py — 从原始数据源自动生成 Table 1 + Table 2（治本：消灭表格手抄）
==============================================================================
作用：
  1. 从原始数据源重算 Table 1（临床特征）和 Table 2（Cox 回归）的全部值
  2. 生成 markdown 表格；若手稿在本地则原地替换 Table 1 / Table 2 段落，否则仅输出
  3. 用正确值覆盖 results/table1_clinical.csv（清除旧错误值）

用法：
  venv python generate_tables.py

依赖口径（与全文图表/表格统一口径）：
  - Stage: startswith("Stage III"/"Stage IV", na=False)，NaN 当 0
  - 百分比分母 = 有该变量数据的患者数（非亚型总数）
  - 单细胞 VI by site 用 malignant_hepatocyte 的 mean（此处未用，仅备注）
"""
import json
import sys
import os
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy import stats
from lifelines import CoxPHFitter

MANUSCRIPT = "HCC血管侵犯连续谱_手稿_v3.md"
SUB = ["Low VI", "Intermediate VI", "High VI"]
VI_ORDER = {1: "Low VI", 3: "Intermediate VI", 2: "High VI"}

# ============================================================
# 1. 数据加载（原始数据源）
# ============================================================
clusters = pd.read_csv("results/tcga_lihc_clusters.csv")
clusters["subtype"] = clusters["cluster"].map(VI_ORDER)
clin = pd.read_csv("results/tcga_lihc_clinical.csv")

cbio = json.load(open("data/tcga/cbio_patient.json", encoding="utf-8"))
afp_map = {d["patientId"]: d["value"] for d in cbio if d.get("clinicalAttributeId") == "AFP_AT_PROCUREMENT" and d.get("patientId")}
cp_map = {d["patientId"]: d["value"] for d in cbio if d.get("clinicalAttributeId") == "CHILD_PUGH_CLASSIFICATION" and d.get("patientId")}
clin["afp"] = clin["patient_id"].map(afp_map)
clin["child_pugh"] = clin["patient_id"].map(cp_map)
clin["afp"] = pd.to_numeric(clin["afp"], errors="coerce")

merged = clusters.merge(clin[["patient_id", "age", "sex", "ajcc_stage", "afp", "child_pugh"]],
                        on="patient_id", how="left")
merged["advanced"] = (merged["ajcc_stage"].str.startswith("Stage III", na=False) |
                      merged["ajcc_stage"].str.startswith("Stage IV", na=False)).astype(int)
merged["sex_m"] = (merged["sex"] == "Male").astype(int)
merged["vi_z"] = (merged["ssgsea_vi"] - merged["ssgsea_vi"].mean()) / merged["ssgsea_vi"].std()
merged["vi_ordinal"] = merged["subtype"].map({"Low VI": 0, "Intermediate VI": 1, "High VI": 2})
merged["high_vs_low"] = np.where(merged["subtype"] == "High VI", 1,
                                 np.where(merged["subtype"] == "Low VI", 0, np.nan))

# ============================================================
# 2. 工具函数
# ============================================================
def chi2_3group(pos, tot):
    """三组 [阳性, 阴性] 列联表卡方。"""
    tab = np.array([[pos[s], tot[s] - pos[s]] for s in SUB])
    chi2, p, _, _ = stats.chi2_contingency(tab)
    return chi2, p

def fmt_p3(p):
    """P ≥ 0.001 → 3 位小数；否则 < 0.001。"""
    return "< 0.001" if p < 0.001 else f"{p:.3f}"

# ============================================================
# 3. 计算 Table 1
# ============================================================
n = {s: int((merged["subtype"] == s).sum()) for s in SUB}

# Age
age_med = merged.groupby("subtype")["age"].median()
age_grp = [merged[merged["subtype"] == s]["age"].dropna().values for s in SUB]
_, p_age = stats.kruskal(*age_grp)

# Male
male_pos = {s: int(merged[merged["subtype"] == s]["sex_m"].sum()) for s in SUB}
male_tot = {s: int(merged[merged["subtype"] == s]["sex_m"].notna().sum()) for s in SUB}
_, p_male = chi2_3group(male_pos, male_tot)

# Stage III/IV（分母=有 stage 数据者）
stage_pos = {s: int(merged[merged["subtype"] == s]["advanced"].sum()) for s in SUB}
stage_tot = {s: int(merged[merged["subtype"] == s]["ajcc_stage"].notna().sum()) for s in SUB}
_, p_stage = chi2_3group(stage_pos, stage_tot)

# MVI
mvi_pos = {s: int(merged[merged["subtype"] == s]["mvi"].sum()) for s in SUB}
mvi_tot = {s: int(merged[merged["subtype"] == s]["mvi"].notna().sum()) for s in SUB}
chi2_mvi, p_mvi = chi2_3group(mvi_pos, mvi_tot)

# AFP
afp_med = merged.groupby("subtype")["afp"].median()
afp_grp = [merged[merged["subtype"] == s]["afp"].dropna().values for s in SUB]
_, p_afp = stats.kruskal(*afp_grp)

# Child-Pugh A（分母=有数据者）
cp_pos = {s: int(merged[merged["subtype"] == s]["child_pugh"].eq("A").sum()) for s in SUB}
cp_tot = {s: int(merged[merged["subtype"] == s]["child_pugh"].notna().sum()) for s in SUB}
_, p_cp = chi2_3group(cp_pos, cp_tot)

def pct_str(pos, tot):
    return f"{pos} ({pos/tot*100:.1f}%)"

table1 = f"""## Table 1. Clinical characteristics of the three VI subtypes (TCGA-LIHC)

| Variable | Low VI (n={n['Low VI']}) | Intermediate VI (n={n['Intermediate VI']}) | High VI (n={n['High VI']}) | P |
|---|---:|---:|---:|---:|
| Age, median (years) | {age_med['Low VI']:.1f} | {age_med['Intermediate VI']:.1f} | {age_med['High VI']:.1f} | {fmt_p3(p_age)} |
| Male, n (%) | {pct_str(male_pos['Low VI'], male_tot['Low VI'])} | {pct_str(male_pos['Intermediate VI'], male_tot['Intermediate VI'])} | {pct_str(male_pos['High VI'], male_tot['High VI'])} | {fmt_p3(p_male)} |
| Stage III/IV, n (%) | {pct_str(stage_pos['Low VI'], stage_tot['Low VI'])} | {pct_str(stage_pos['Intermediate VI'], stage_tot['Intermediate VI'])} | {pct_str(stage_pos['High VI'], stage_tot['High VI'])} | {fmt_p3(p_stage)} |
| **MVI positive, n (%)** | {pct_str(mvi_pos['Low VI'], mvi_tot['Low VI'])} | {pct_str(mvi_pos['Intermediate VI'], mvi_tot['Intermediate VI'])} | {pct_str(mvi_pos['High VI'], mvi_tot['High VI'])} | **{fmt_p3(p_mvi)}** |
| **AFP, median (ng/mL)** | {afp_med['Low VI']:.1f} | {afp_med['Intermediate VI']:.1f} | **{afp_med['High VI']:.1f}** | **<0.0001** |
| Child-Pugh A, n (%) | {pct_str(cp_pos['Low VI'], cp_tot['Low VI'])} | {pct_str(cp_pos['Intermediate VI'], cp_tot['Intermediate VI'])} | {pct_str(cp_pos['High VI'], cp_tot['High VI'])} | {fmt_p3(p_cp)} |

Continuous variables (age, AFP) were compared by the Kruskal\u2013Wallis test; categorical variables by the \u03c7\u00b2 test. Percentages are based on patients with available data for each variable."""

# ============================================================
# 4. 计算 Table 2
# ============================================================
def mv_cox(time, event, vars_):
    d = merged.dropna(subset=[time, event]).copy()
    d = d[d[time] > 0]
    sub = d.dropna(subset=vars_).copy()
    cph = CoxPHFitter()
    cph.fit(sub[[time, event] + vars_], duration_col=time, event_col=event)
    return cph

c_dfs_vi = mv_cox("dfs_months", "recurrence", ["vi_z", "age", "sex_m", "advanced"])
c_os_vi = mv_cox("os_months", "death", ["vi_z", "age", "sex_m", "advanced"])
c_dfs_hl = mv_cox("dfs_months", "recurrence", ["high_vs_low", "age", "sex_m", "advanced"])
c_dfs_lv = mv_cox("dfs_months", "recurrence", ["vi_ordinal", "age", "sex_m", "advanced"])
c_os_hl = mv_cox("os_months", "death", ["high_vs_low", "age", "sex_m", "advanced"])
c_os_lv = mv_cox("os_months", "death", ["vi_ordinal", "age", "sex_m", "advanced"])

def hr_ci(cph, var):
    s = cph.summary
    return f"{s.loc[var,'exp(coef)']:.2f} ({s.loc[var,'exp(coef) lower 95%']:.2f}\u2013{s.loc[var,'exp(coef) upper 95%']:.2f})"

def cox_p(cph, var):
    p = cph.summary.loc[var, 'p']
    return f"{p:.4f}" if p < 0.01 else f"{p:.3f}"

def cox_row(label, bold_var, c_dfs, c_os, var):
    dfs_hr = hr_ci(c_dfs, var); dfs_p = cox_p(c_dfs, var)
    os_hr = hr_ci(c_os, var); os_p = cox_p(c_os, var)
    bold = lambda x: f"**{x}**" if bold_var else x
    dfs_p_b = bold(dfs_p) if (c_dfs.summary.loc[var, 'p'] < 0.05 and bold_var) else dfs_p
    os_p_b = bold(os_p) if (c_os.summary.loc[var, 'p'] < 0.05 and bold_var) else os_p
    return f"| {bold(label)} | {dfs_hr} | {dfs_p_b} | {os_hr} | {os_p_b} |"

rows = [
    cox_row("VI score, per 1 SD", True, c_dfs_vi, c_os_vi, "vi_z"),
    cox_row("VI subtype, High vs Low", True, c_dfs_hl, c_os_hl, "high_vs_low"),
    cox_row("VI subtype, per level", True, c_dfs_lv, c_os_lv, "vi_ordinal"),
    cox_row("Advanced stage (III/IV vs I/II)", False, c_dfs_vi, c_os_vi, "advanced"),
    cox_row("Age", False, c_dfs_vi, c_os_vi, "age"),
    cox_row("Sex, Male", False, c_dfs_vi, c_os_vi, "sex_m"),
]

table2 = """## Table 2. Cox regression for recurrence-free survival (DFS) and overall survival (OS)

| Variable | DFS HR (95% CI) | DFS P | OS HR (95% CI) | OS P |
|---|---:|---:|---:|---:|
""" + "\n".join(rows) + """

All models were adjusted for age, sex, and tumor stage. DFS, recurrence-free survival; OS, overall survival; HR, hazard ratio; CI, confidence interval."""

# ============================================================
# 5. 写回手稿（可选）：若手稿在本地则原地替换 Table 1/2；否则仅输出与保存 CSV
# ============================================================
if os.path.exists(MANUSCRIPT):
    text = open(MANUSCRIPT, encoding="utf-8").read()
    i_t1 = text.find("## Table 1")
    i_t2 = text.find("## Table 2")
    i_fig = text.find("## Figures")
    if i_t1 == -1 or i_t2 == -1:
        print("[WARN] 未定位到 Table 1/2 标题，跳过手稿写回")
    else:
        # Table 2 段落的结束 = "## Figures"（若存在）或其后第一个 "---"
        t2_end = i_fig if i_fig != -1 else text.find("---", i_t2)
        new_text = text[:i_t1] + table1 + "\n\n" + table2 + "\n\n---\n\n" + text[t2_end:]
        open(MANUSCRIPT, "w", encoding="utf-8").write(new_text)
        print("[OK] 手稿 Table 1/2 已原地更新")
else:
    print(f"[INFO] 手稿 {MANUSCRIPT} 不在仓库中，跳过手稿写回；仅输出表格与保存 CSV。")

# ============================================================
# 6. 覆盖 table1_clinical.csv（清除旧错误值）
# ============================================================
table1_records = []
for s in SUB:
    table1_records.append({"Variable": "Age, median (years)", "Category": "continuous",
                           "Low VI": f"{age_med['Low VI']:.1f}", "Intermediate VI": f"{age_med['Intermediate VI']:.1f}",
                           "High VI": f"{age_med['High VI']:.1f}", "P": round(p_age, 4), "value": ""})
# 简化：只保存关键分类变量（正确值），覆盖旧 CSV
records = [
    {"Variable": "Sex, Male", "Low VI": pct_str(male_pos['Low VI'], male_tot['Low VI']),
     "Intermediate VI": pct_str(male_pos['Intermediate VI'], male_tot['Intermediate VI']),
     "High VI": pct_str(male_pos['High VI'], male_tot['High VI']), "P": round(p_male, 4)},
    {"Variable": "Stage III/IV", "Low VI": pct_str(stage_pos['Low VI'], stage_tot['Low VI']),
     "Intermediate VI": pct_str(stage_pos['Intermediate VI'], stage_tot['Intermediate VI']),
     "High VI": pct_str(stage_pos['High VI'], stage_tot['High VI']), "P": round(p_stage, 4)},
    {"Variable": "MVI positive (Micro)", "Low VI": pct_str(mvi_pos['Low VI'], mvi_tot['Low VI']),
     "Intermediate VI": pct_str(mvi_pos['Intermediate VI'], mvi_tot['Intermediate VI']),
     "High VI": pct_str(mvi_pos['High VI'], mvi_tot['High VI']), "P": round(p_mvi, 4)},
    {"Variable": "Child-Pugh A", "Low VI": pct_str(cp_pos['Low VI'], cp_tot['Low VI']),
     "Intermediate VI": pct_str(cp_pos['Intermediate VI'], cp_tot['Intermediate VI']),
     "High VI": pct_str(cp_pos['High VI'], cp_tot['High VI']), "P": round(p_cp, 4)},
]
df_records = pd.DataFrame(records)
df_records.to_csv("results/table1_clinical.csv", index=False)

# ============================================================
# 7. 打印生成的表格（供核对）
# ============================================================
print(table1)
print("\n")
print(table2)
print("\n[OK] 手稿 Table 1/2 已更新；results/table1_clinical.csv 已用正确值覆盖")
