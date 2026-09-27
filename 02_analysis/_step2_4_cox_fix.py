import warnings
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

warnings.filterwarnings("ignore")

print("重新生成 Cox 回归表（修正 VI score scale 问题）")

# ---------- 数据 ----------
clusters = pd.read_csv("results/tcga_lihc_clusters.csv")
clin = pd.read_csv("results/tcga_lihc_clinical.csv")
merged = clusters.merge(clin[["patient_id", "age", "sex", "ajcc_stage"]], on="patient_id", how="left")
vi_order = {1: "Low VI", 3: "Intermediate VI", 2: "High VI"}
merged["subtype"] = merged["cluster"].map(vi_order)
merged["advanced"] = (merged["ajcc_stage"].str.startswith("Stage III", na=False) |
                      merged["ajcc_stage"].str.startswith("Stage IV", na=False)).astype(int)
merged["sex_m"] = (merged["sex"] == "Male").astype(int)
# 标准化 VI score（每 1 SD）
merged["vi_z"] = (merged["ssgsea_vi"] - merged["ssgsea_vi"].mean()) / merged["ssgsea_vi"].std()
# 分型有序（Low=0, Int=1, High=2）
merged["vi_ordinal"] = merged["subtype"].map({"Low VI": 0, "Intermediate VI": 1, "High VI": 2})
# High vs Low 二分类
merged["high_vs_low"] = np.where(merged["subtype"] == "High VI", 1,
                                 np.where(merged["subtype"] == "Low VI", 0, np.nan))

def run_cox(df, time_col, event_col):
    d = df.dropna(subset=[time_col, event_col]).copy()
    d = d[d[time_col] > 0]
    rows = []

    # 单因素
    uni_vars = [
        ("vi_z", "VI score (per SD)"),
        ("vi_ordinal", "VI subtype (per level)"),
        ("age", "Age (continuous)"),
        ("sex_m", "Sex (Male vs Female)"),
        ("advanced", "Advanced stage (III/IV vs I/II)"),
    ]
    for var, name in uni_vars:
        sub = d.dropna(subset=[var])
        cph = CoxPHFitter()
        try:
            cph.fit(sub[[time_col, event_col, var]], duration_col=time_col, event_col=event_col)
            hr = float(np.exp(cph.params_[var]))
            ci = np.exp(cph.confidence_intervals_.loc[var].values)
            p = float(cph.summary.loc[var, "p"])
            rows.append({"Analysis": "Univariate", "Variable": name, "HR": round(hr, 3),
                         "CI_lower": round(float(ci[0]), 3), "CI_upper": round(float(ci[1]), 3),
                         "P": f"{p:.4f}"})
        except Exception as e:
            rows.append({"Analysis": "Univariate", "Variable": name, "HR": np.nan,
                         "CI_lower": np.nan, "CI_upper": np.nan, "P": "NA"})

    # High vs Low 单独（分类 HR）
    sub = d.dropna(subset=["high_vs_low"])
    sub = sub[sub["high_vs_low"].isin([0, 1])]
    if len(sub) > 20:
        cph = CoxPHFitter()
        cph.fit(sub[[time_col, event_col, "high_vs_low"]], duration_col=time_col, event_col=event_col)
        hr = float(np.exp(cph.params_["high_vs_low"]))
        ci = np.exp(cph.confidence_intervals_.loc["high_vs_low"].values)
        p = float(cph.summary.loc["high_vs_low", "p"])
        rows.append({"Analysis": "Univariate", "Variable": "VI subtype (High vs Low)", "HR": round(hr, 3),
                     "CI_lower": round(float(ci[0]), 3), "CI_upper": round(float(ci[1]), 3), "P": f"{p:.4f}"})

    # 多因素（VI score + age + sex + stage）
    sub = d.dropna(subset=["vi_z", "age", "sex_m", "advanced"]).copy()
    cph = CoxPHFitter()
    try:
        cph.fit(sub[[time_col, event_col, "vi_z", "age", "sex_m", "advanced"]],
                duration_col=time_col, event_col=event_col)
        mv_vars = [
            ("vi_z", "VI score (per SD)"),
            ("age", "Age (continuous)"),
            ("sex_m", "Sex (Male vs Female)"),
            ("advanced", "Advanced stage (III/IV vs I/II)"),
        ]
        for var, name in mv_vars:
            hr = float(np.exp(cph.params_[var]))
            ci = np.exp(cph.confidence_intervals_.loc[var].values)
            p = float(cph.summary.loc[var, "p"])
            rows.append({"Analysis": "Multivariate", "Variable": name, "HR": round(hr, 3),
                         "CI_lower": round(float(ci[0]), 3), "CI_upper": round(float(ci[1]), 3), "P": f"{p:.4f}"})
    except Exception as e:
        print(f"  多因素 Cox 失败: {e}")

    return pd.DataFrame(rows)

print("\n=== OS Cox ===")
os_df = run_cox(merged, "os_months", "death")
print(os_df.to_string(index=False))
os_df.to_csv("results/cox_os.csv", index=False)

print("\n=== DFS Cox ===")
dfs_df = run_cox(merged, "dfs_months", "recurrence")
print(dfs_df.to_string(index=False))
dfs_df.to_csv("results/cox_dfs.csv", index=False)

print("\n已重存 results/cox_os.csv 和 results/cox_dfs.csv（修正版，HR 有意义）。")
