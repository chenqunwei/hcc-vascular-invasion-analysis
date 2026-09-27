import json
import pandas as pd

data = json.load(open("data/tcga/cbio_patient.json", encoding="utf-8"))
patients = {}
for d in data:
    pid = d.get("patientId")
    attr = d.get("clinicalAttributeId")
    val = d.get("value")
    if pid and attr:
        patients.setdefault(pid, {})[attr] = val

rows = []
for pid, p in patients.items():
    rows.append({
        "patient_id": pid,
        "vascular_invasion": p.get("VASCULAR_INVASION"),
        "dfs_status": p.get("DFS_STATUS"),
        "dfs_months": p.get("DFS_MONTHS"),
        "os_status": p.get("OS_STATUS"),
        "os_months": p.get("OS_MONTHS"),
        "ajcc_stage": p.get("AJCC_PATHOLOGIC_TUMOR_STAGE"),
        "age": p.get("AGE"),
        "sex": p.get("SEX"),
    })

df = pd.DataFrame(rows)
df["dfs_months"] = pd.to_numeric(df["dfs_months"], errors="coerce")
df["os_months"] = pd.to_numeric(df["os_months"], errors="coerce")
df["age"] = pd.to_numeric(df["age"], errors="coerce")

# MVI 三分类 → 二分类（Micro=1 为 MVI+，None=0 为 MVI-，Macro 单独标记）
df["mvi"] = df["vascular_invasion"].map({"Micro": 1, "None": 0, "Macro": 2})
# 复发二分类（值格式 "1:Recurred/Progressed" / "0:DiseaseFree"）
df["recurrence"] = df["dfs_status"].map(lambda x: 1 if x and x.startswith("1:") else (0 if x and x.startswith("0:") else None))
# 死亡二分类（值格式 "1:DECEASED" / "0:LIVING"）
df["death"] = df["os_status"].map(lambda x: 1 if x and x.startswith("1:") else (0 if x and x.startswith("0:") else None))

df.to_csv("results/tcga_lihc_clinical.csv", index=False)

print("临床表已保存 results/tcga_lihc_clinical.csv")
print("总样本:", len(df))
print("\n血管侵犯分布:")
print(df["vascular_invasion"].value_counts(dropna=False).to_string())
print("\nMVI 二分类 (Micro=1/None=0/Macro=2):")
print(df["mvi"].value_counts(dropna=False).to_string())
print("\n复发分布:")
print(df["dfs_status"].value_counts(dropna=False).to_string())
print("\nMVI 有值 + 复发有值 的样本数:",
      df[df["mvi"].notna() & df["recurrence"].notna()].shape[0])
print("\nMVI 有值 + OS 有值 的样本数:",
      df[df["mvi"].notna() & df["death"].notna()].shape[0])
