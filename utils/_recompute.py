# -*- coding: utf-8 -*-
"""Recompute ALL manuscript numbers from the NEW deterministic clusters.csv (seed=42).
Outputs pure-ASCII report so it can be read reliably.
"""
import json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy import stats
from lifelines import CoxPHFitter

SUB = ["Low VI", "Intermediate VI", "High VI"]
VI_ORDER = {1: "Low VI", 3: "Intermediate VI", 2: "High VI"}

clusters = pd.read_csv("results/tcga_lihc_clusters.csv")
clusters["subtype"] = clusters["cluster"].map(VI_ORDER)
clin = pd.read_csv("results/tcga_lihc_clinical.csv")

cbio = json.load(open("data/tcga/cbio_patient.json", encoding="utf-8"))
afp_map = {d["patientId"]: d["value"] for d in cbio if d.get("clinicalAttributeId") == "AFP_AT_PROCUREMENT" and d.get("patientId")}
cp_map = {d["patientId"]: d["value"] for d in cbio if d.get("clinicalAttributeId") == "CHILD_PUGH_CLASSIFICATION" and d.get("patientId")}
clin["afp"] = clin["patient_id"].map(afp_map)
clin["child_pugh"] = clin["patient_id"].map(cp_map)
clin["afp"] = pd.to_numeric(clin["afp"], errors="coerce")

merged = clusters.merge(clin[["patient_id", "age", "sex", "ajcc_stage", "afp", "child_pugh"]], on="patient_id", how="left")
merged["advanced"] = (merged["ajcc_stage"].str.startswith("Stage III", na=False) |
                      merged["ajcc_stage"].str.startswith("Stage IV", na=False)).astype(int)
merged["sex_m"] = (merged["sex"] == "Male").astype(int)
merged["vi_z"] = (merged["ssgsea_vi"] - merged["ssgsea_vi"].mean()) / merged["ssgsea_vi"].std()
merged["vi_ordinal"] = merged["subtype"].map({"Low VI": 0, "Intermediate VI": 1, "High VI": 2})
merged["high_vs_low"] = np.where(merged["subtype"] == "High VI", 1, np.where(merged["subtype"] == "Low VI", 0, np.nan))

def chi2_3(pos, tot):
    tab = np.array([[pos[s], tot[s]-pos[s]] for s in SUB])
    c, p, _, _ = stats.chi2_contingency(tab)
    return c, p

out = {}
sub_n = merged["subtype"].value_counts()
out["n_low"] = int(sub_n["Low VI"]); out["n_int"] = int(sub_n["Intermediate VI"]); out["n_high"] = int(sub_n["High VI"])
med_vi = merged.groupby("subtype")["ssgsea_vi"].median()
out["med_vi_low"] = round(float(med_vi["Low VI"]),1); out["med_vi_int"]=round(float(med_vi["Intermediate VI"]),1); out["med_vi_high"]=round(float(med_vi["High VI"]),1)

mvi_pos = {s: int(merged[merged["subtype"]==s]["mvi"].sum()) for s in SUB}
mvi_cnt = {s: int(merged[merged["subtype"]==s]["mvi"].notna().sum()) for s in SUB}
rate = {s: mvi_pos[s]/mvi_cnt[s]*100 for s in SUB}
c_mvi, p_mvi = chi2_3(mvi_pos, mvi_cnt)
out["mvi_low"]=f"{mvi_pos['Low VI']} ({rate['Low VI']:.1f}%)"; out["mvi_int"]=f"{mvi_pos['Intermediate VI']} ({rate['Intermediate VI']:.1f}%)"; out["mvi_high"]=f"{mvi_pos['High VI']} ({rate['High VI']:.1f}%)"
out["mvi_chi2"]=round(float(c_mvi),2); out["mvi_p"]=round(float(p_mvi),4)

afp_med = merged.groupby("subtype")["afp"].median()
ag = [merged[merged["subtype"]==s]["afp"].dropna().values for s in SUB]
_, p_afp = stats.kruskal(*ag)
out["afp_low"]=round(float(afp_med["Low VI"]),1); out["afp_int"]=round(float(afp_med["Intermediate VI"]),1); out["afp_high"]=round(float(afp_med["High VI"]),1); out["afp_p"]=f"{p_afp:.2e}"

def adv_by(s):
    g=merged[merged["subtype"]==s]; return int(g["advanced"].sum()), int(g["ajcc_stage"].notna().sum())
adv={s:adv_by(s) for s in SUB}
c_a,p_a,_ ,_=stats.chi2_contingency(np.array([[adv[s][0],adv[s][1]-adv[s][0]] for s in SUB]))
out["stage_low"]=f"{adv['Low VI'][0]} ({adv['Low VI'][0]/adv['Low VI'][1]*100:.1f}%)"
out["stage_int"]=f"{adv['Intermediate VI'][0]} ({adv['Intermediate VI'][0]/adv['Intermediate VI'][1]*100:.1f}%)"
out["stage_high"]=f"{adv['High VI'][0]} ({adv['High VI'][0]/adv['High VI'][1]*100:.1f}%)"
out["stage_p"]=round(float(p_a),4)

def cp_by(s):
    g=merged[merged["subtype"]==s]; return int(g["child_pugh"].eq("A").sum()), int(g["child_pugh"].notna().sum())
cpa={s:cp_by(s) for s in SUB}
c_c,p_c,_,_=stats.chi2_contingency(np.array([[cpa[s][0],cpa[s][1]-cpa[s][0]] for s in SUB]))
out["cp_low"]=f"{cpa['Low VI'][0]} ({cpa['Low VI'][0]/cpa['Low VI'][1]*100:.1f}%)"
out["cp_int"]=f"{cpa['Intermediate VI'][0]} ({cpa['Intermediate VI'][0]/cpa['Intermediate VI'][1]*100:.1f}%)"
out["cp_high"]=f"{cpa['High VI'][0]} ({cpa['High VI'][0]/cpa['High VI'][1]*100:.1f}%)"
out["cp_p"]=round(float(p_c),4)

male_by={s:(int(merged[merged["subtype"]==s]["sex_m"].sum()),int(merged[merged["subtype"]==s]["sex_m"].notna().sum())) for s in SUB}
c_m,p_m,_,_=stats.chi2_contingency(np.array([[male_by[s][0],male_by[s][1]-male_by[s][0]] for s in SUB]))
out["male_low"]=f"{male_by['Low VI'][0]} ({male_by['Low VI'][0]/male_by['Low VI'][1]*100:.1f}%)"
out["male_int"]=f"{male_by['Intermediate VI'][0]} ({male_by['Intermediate VI'][0]/male_by['Intermediate VI'][1]*100:.1f}%)"
out["male_high"]=f"{male_by['High VI'][0]} ({male_by['High VI'][0]/male_by['High VI'][1]*100:.1f}%)"
out["male_p"]=round(float(p_m),4)

age_med=merged.groupby("subtype")["age"].median()
ag2=[merged[merged["subtype"]==s]["age"].dropna().values for s in SUB]
_,p_age=stats.kruskal(*ag2)
out["age_low"]=round(float(age_med["Low VI"]),1); out["age_int"]=round(float(age_med["Intermediate VI"]),1); out["age_high"]=round(float(age_med["High VI"]),1); out["age_p"]=round(float(p_age),4)

def mv_cox(time,event,vars_):
    d=merged.dropna(subset=[time,event]).copy(); d=d[d[time]>0]; s=d.dropna(subset=vars_).copy()
    cph=CoxPHFitter(); cph.fit(s[[time,event]+vars_],duration_col=time,event_col=event); return cph
c_dfs_vi=mv_cox("dfs_months","recurrence",["vi_z","age","sex_m","advanced"])
c_os_vi=mv_cox("os_months","death",["vi_z","age","sex_m","advanced"])
c_dfs_hl=mv_cox("dfs_months","recurrence",["high_vs_low","age","sex_m","advanced"])
c_dfs_lv=mv_cox("dfs_months","recurrence",["vi_ordinal","age","sex_m","advanced"])
c_os_hl=mv_cox("os_months","death",["high_vs_low","age","sex_m","advanced"])
c_os_lv=mv_cox("os_months","death",["vi_ordinal","age","sex_m","advanced"])
def hr(cph,v):
    s=cph.summary; return f"{s.loc[v,'exp(coef)']:.2f} ({s.loc[v,'exp(coef) lower 95%']:.2f}-{s.loc[v,'exp(coef) upper 95%']:.2f})"
def pp(cph,v):
    p=cph.summary.loc[v,'p']; return f"{p:.4f}" if p<0.01 else f"{p:.3f}"
out["cox_vi_dfs"]=hr(c_dfs_vi,"vi_z"); out["cox_vi_dfs_p"]=pp(c_dfs_vi,"vi_z")
out["cox_vi_os"]=hr(c_os_vi,"vi_z"); out["cox_vi_os_p"]=pp(c_os_vi,"vi_z")
out["cox_hl_dfs"]=hr(c_dfs_hl,"high_vs_low"); out["cox_hl_dfs_p"]=pp(c_dfs_hl,"high_vs_low")
out["cox_lv_dfs"]=hr(c_dfs_lv,"vi_ordinal"); out["cox_lv_dfs_p"]=pp(c_dfs_lv,"vi_ordinal")
out["cox_hl_os"]=hr(c_os_hl,"high_vs_low"); out["cox_hl_os_p"]=pp(c_os_hl,"high_vs_low")
out["cox_lv_os"]=hr(c_os_lv,"vi_ordinal"); out["cox_lv_os_p"]=pp(c_os_lv,"vi_ordinal")

# survival log-rank
surv=pd.read_csv("results/survival_summary.csv")
out["dfs_lr"]=round(float(surv["DFS_logrank_p_3group"][0]),4)
out["os_lr"]=round(float(surv["OS_logrank_p_3group"][0]),4)

# external
ext=pd.read_csv("results/gse14520_vi_recurrence.csv").dropna(subset=["ssgsea_vi","recurr","recurr_months"])
ext=ext[ext["recurr_months"]>0]
out["gse14520_n"]=int(len(ext))

with open("_recompute_out.txt","w",encoding="ascii") as f:
    for k,v in out.items():
        f.write(f"{k} = {v}\n")
print("OK wrote _recompute_out.txt")
