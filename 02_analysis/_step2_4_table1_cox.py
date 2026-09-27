import json
import numpy as np
import pandas as pd
from scipy import stats
from lifelines import CoxPHFitter

print("=" * 60)
print("Table 1（三型临床基线）+ 单/多因素 Cox 回归")
print("=" * 60)

# ---------- 1. 提取扩展临床字段 ----------
print("[1] 提取扩展临床字段 ...")
data = json.load(open("data/tcga/cbio_patient.json", encoding="utf-8"))
patients = {}
for d in data:
    pid = d.get("patientId"); attr = d.get("clinicalAttributeId"); val = d.get("value")
    if pid and attr:
        patients.setdefault(pid, {})[attr] = val

field_map = {
    "AGE": "age", "SEX": "sex", "VASCULAR_INVASION": "vascular_invasion",
    "AJCC_PATHOLOGIC_TUMOR_STAGE": "stage", "AFP_AT_PROCUREMENT": "afp",
    "CHILD_PUGH_CLASSIFICATION": "child_pugh", "HISTORY_HEPATO_CARCINOMA_RISK_FACTORS": "risk_factors",
    "ISHAK_FIBROSIS_SCORE": "fibrosis", "GRADE": "grade",
    "OS_STATUS": "os_status", "OS_MONTHS": "os_months",
    "DFS_STATUS": "dfs_status", "DFS_MONTHS": "dfs_months",
}
rows = []
for pid, p in patients.items():
    row = {"patient_id": pid}
    for src, dst in field_map.items():
        row[dst] = p.get(src)
    rows.append(row)
clin = pd.DataFrame(rows)
clin["age"] = pd.to_numeric(clin["age"], errors="coerce")
clin["afp"] = pd.to_numeric(clin["afp"], errors="coerce")
clin["os_months"] = pd.to_numeric(clin["os_months"], errors="coerce")
clin["dfs_months"] = pd.to_numeric(clin["dfs_months"], errors="coerce")
# 二分类
clin["death"] = clin["os_status"].map(lambda x: 1 if x and x.startswith("1:") else (0 if x and x.startswith("0:") else None))
clin["recurrence"] = clin["dfs_status"].map(lambda x: 1 if x and x.startswith("1:") else (0 if x and x.startswith("0:") else None))
clin["mvi"] = clin["vascular_invasion"].map({"Micro": 1, "None": 0, "Macro": 2})
# 简化 stage
def simplify_stage(s):
    if not isinstance(s, str):
        return None
    s = s.strip()
    if s.startswith("Stage I "): return "Stage I"
    if s.startswith("Stage II"): return "Stage II"
    if s.startswith("Stage III"): return "Stage III"
    if s.startswith("Stage IV"): return "Stage IV"
    return s
clin["stage_simple"] = clin["stage"].apply(simplify_stage)
print(f"  临床表: {len(clin)} 例")

# ---------- 2. 合并分型 ----------
print("[2] 合并分型 ...")
clusters = pd.read_csv("results/tcga_lihc_clusters.csv")
vi_order = {1: "Low VI", 3: "Intermediate VI", 2: "High VI"}
clusters["subtype"] = clusters["cluster"].map(vi_order)
order = ["Low VI", "Intermediate VI", "High VI"]

merged = clusters[["patient_id", "subtype", "ssgsea_vi"]].merge(clin, on="patient_id", how="left")
merged = merged[merged["subtype"].notna()]
print(f"  合并后: {len(merged)} 例")

# ---------- 3. Table 1 ----------
print("\n[3] 生成 Table 1 ...")
table1_rows = []

def add_row(name, category, value):
    table1_rows.append({"Variable": name, "Category": category, **value})

for st in order:
    pass  # 占位

# 连续变量：age, afp（中位数 + IQR，Kruskal 检验）
for var, name in [("age", "Age (median)"), ("afp", "AFP (median, ng/mL)")]:
    vals = {}
    medians = {}
    for st in order:
        s = merged[merged["subtype"] == st][var].dropna()
        medians[st] = f"{s.median():.1f}"
    groups = [merged[merged['subtype']==st][var].dropna() for st in order]
    groups = [g for g in groups if len(g) >= 2]
    p = stats.kruskal(*groups).pvalue if len(groups) >= 2 else np.nan
    add_row(name, "continuous", {st: medians[st] for st in order} | {"P": f"{p:.4f}"})

# 分类变量：sex, stage, grade, mvi, child_pugh, fibrosis
def cat_table(var, name, mapping=None):
    vals = {}
    p = np.nan
    # 卡方
    ct = pd.crosstab(merged["subtype"], merged[var])
    if ct.shape[0] >= 2 and ct.shape[1] >= 2:
        try:
            from scipy.stats import chi2_contingency
            p = chi2_contingency(ct)[1]
        except Exception:
            p = np.nan
    for st in order:
        sub = merged[merged["subtype"] == st]
        cnt = sub[var].value_counts()
        total = len(sub)
        # 取该变量最常见的类别做展示
    return p, ct

# sex
p, ct = cat_table("sex", "Sex")
for st in order:
    sub = merged[merged["subtype"] == st]
    n_male = (sub["sex"] == "Male").sum()
    table1_rows.append({"Variable": "Sex (Male, %)", "Category": st, "value": f"{n_male}/{len(sub)} ({n_male/len(sub)*100:.0f}%)"})

# 简化：直接构建 Table 1 DataFrame
t1 = pd.DataFrame(merged)
t1_p = {}

# 连续变量统计
cont_vars = [("age", "Age (years)"), ("afp", "AFP (ng/mL)")]
for var, label in cont_vars:
    row = {"Variable": label}
    for st in order:
        s = merged[merged["subtype"] == st][var].dropna()
        row[st] = f"{s.median():.1f}"
    g = [merged[merged['subtype']==st][var].dropna() for st in order]
    g = [x for x in g if len(x) >= 2]
    row["P"] = f"{stats.kruskal(*g).pvalue:.4f}" if len(g) >= 2 else "-"
    table1_rows.append(row)

# 分类变量统计（计数+比例）
cat_vars = [
    ("sex", "Male", "Sex, Male"),
    ("stage_simple", "Stage III/IV", "Advanced stage (III/IV)"),
    ("mvi", 1.0, "MVI positive"),
    ("grade", None, None),  # 特殊处理
]
# sex
row = {"Variable": "Sex, Male"}
for st in order:
    s = merged[merged["subtype"] == st]
    n = (s["sex"] == "Male").sum()
    row[st] = f"{n}/{len(s)} ({n/len(s)*100:.0f}%)"
ct = pd.crosstab(merged["subtype"], merged["sex"])
row["P"] = f"{stats.chi2_contingency(ct)[1]:.4f}"
table1_rows.append(row)

# stage III/IV
row = {"Variable": "Stage III/IV"}
merged["advanced"] = merged["stage_simple"].isin(["Stage III", "Stage IV"])
for st in order:
    s = merged[merged["subtype"] == st]
    n = s["advanced"].sum()
    row[st] = f"{n}/{len(s)} ({n/len(s)*100:.0f}%)"
ct = pd.crosstab(merged["subtype"], merged["advanced"])
row["P"] = f"{stats.chi2_contingency(ct)[1]:.4f}"
table1_rows.append(row)

# MVI positive
row = {"Variable": "MVI positive (Micro)"}
for st in order:
    s = merged[merged["subtype"] == st]
    n = (s["mvi"] == 1).sum()
    row[st] = f"{n}/{len(s)} ({n/len(s)*100:.0f}%)"
ct = pd.crosstab(merged["subtype"], merged["mvi"] == 1)
row["P"] = f"{stats.chi2_contingency(ct)[1]:.4f}"
table1_rows.append(row)

# Child-Pugh A
row = {"Variable": "Child-Pugh A"}
for st in order:
    s = merged[merged["subtype"] == st]
    n = (s["child_pugh"] == "A").sum()
    row[st] = f"{n}/{len(s)} ({n/len(s)*100:.0f}%)"
ct = pd.crosstab(merged["subtype"], merged["child_pugh"] == "A")
row["P"] = f"{stats.chi2_contingency(ct)[1]:.4f}"
table1_rows.append(row)

table1 = pd.DataFrame(table1_rows)
table1.to_csv("results/table1_clinical.csv", index=False)
print("  Table 1:")
print(table1.to_string(index=False))

# ---------- 4. Cox 回归 ----------
print("\n[4] Cox 回归（OS + DFS）...")

def cox_analysis(df, time_col, event_col, label):
    d = df.dropna(subset=[time_col, event_col]).copy()
    d = d[d[time_col] > 0]
    results = []
    # 单因素
    factors = {
        "VI score (continuous)": "ssgsea_vi",
        "Age (continuous)": "age",
        "Sex (Male vs Female)": "sex",
        "Advanced stage (III/IV vs I/II)": "advanced",
    }
    for name, var in factors.items():
        sub = d.dropna(subset=[var])
        if var == "sex":
            sub = sub.copy(); sub["sex_m"] = (sub["sex"] == "Male").astype(int)
            var = "sex_m"
        elif var == "advanced":
            sub = sub.copy(); sub["advanced"] = sub["advanced"].astype(int)
        if len(sub) < 20:
            results.append({"Factor": name, "HR": np.nan, "CI": "-", "P": np.nan})
            continue
        cph = CoxPHFitter()
        try:
            cph.fit(sub[[time_col, event_col, var]], duration_col=time_col, event_col=event_col)
            hr = np.exp(cph.params_[var])
            ci = np.exp(cph.confidence_intervals_.loc[var])
            p = cph.summary.loc[var, "p"]
            results.append({"Factor": name, "HR": f"{hr:.2f}", "CI": f"{ci[0]:.2f}-{ci[1]:.2f}", "P": f"{p:.4f}"})
        except Exception as e:
            results.append({"Factor": name, "HR": "err", "CI": "-", "P": "-"})
    # 多因素
    sub = d.dropna(subset=["ssgsea_vi", "age", "sex", "advanced"]).copy()
    sub["sex_m"] = (sub["sex"] == "Male").astype(int)
    sub["advanced"] = sub["advanced"].astype(int)
    cph = CoxPHFitter()
    try:
        cph.fit(sub[[time_col, event_col, "ssgsea_vi", "age", "sex_m", "advanced"]],
                duration_col=time_col, event_col=event_col)
        for var, name in [("ssgsea_vi", "VI score (multivariate)"), ("age", "Age (multivariate)"),
                          ("sex_m", "Sex (multivariate)"), ("advanced", "Stage (multivariate)")]:
            hr = np.exp(cph.params_[var])
            ci = np.exp(cph.confidence_intervals_.loc[var])
            p = cph.summary.loc[var, "p"]
            results.append({"Factor": name, "HR": f"{hr:.2f}", "CI": f"{ci[0]:.2f}-{ci[1]:.2f}", "P": f"{p:.4f}"})
    except Exception as e:
        print(f"  多因素 Cox 失败: {e}")
    res_df = pd.DataFrame(results)
    print(f"\n  === {label} Cox 回归 ===")
    print(res_df.to_string(index=False))
    return res_df

os_cox = cox_analysis(merged, "os_months", "death", "OS")
dfs_cox = cox_analysis(merged, "dfs_months", "recurrence", "DFS")
os_cox.to_csv("results/cox_os.csv", index=False)
dfs_cox.to_csv("results/cox_dfs.csv", index=False)

print("\n[5] 已保存: results/table1_clinical.csv, cox_os.csv, cox_dfs.csv")
print("\nTable 1 + Cox 完成。")
