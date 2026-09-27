# -*- coding: utf-8 -*-
"""Compute authoritative values for every manuscript number directly from data."""
import json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from lifelines import CoxPHFitter
from lifelines.statistics import multivariate_logrank_test
from collections import Counter

SUB = ["Low VI", "Intermediate VI", "High VI"]
VI_ORDER = {1: "Low VI", 3: "Intermediate VI", 2: "High VI"}

def hdr(t): print("\n" + "="*66 + f"\n{t}\n" + "="*66)

# ---------- load ----------
clusters = pd.read_csv("results/tcga_lihc_clusters.csv")
clusters["subtype"] = clusters["cluster"].map(VI_ORDER)
clin = pd.read_csv("results/tcga_lihc_clinical.csv")
cbio = json.load(open("data/tcga/cbio_patient.json", encoding="utf-8"))
afp_map = {d["patientId"]: d["value"] for d in cbio if d.get("clinicalAttributeId")=="AFP_AT_PROCUREMENT" and d.get("patientId")}
cp_map = {d["patientId"]: d["value"] for d in cbio if d.get("clinicalAttributeId")=="CHILD_PUGH_CLASSIFICATION" and d.get("patientId")}
clin["afp"] = pd.to_numeric(clin["patient_id"].map(afp_map), errors="coerce")
clin["child_pugh"] = clin["patient_id"].map(cp_map)
merged = clusters.merge(clin[["patient_id","age","sex","ajcc_stage","afp","child_pugh"]], on="patient_id", how="left")
merged["advanced"] = (merged["ajcc_stage"].str.startswith("Stage III", na=False) |
                      merged["ajcc_stage"].str.startswith("Stage IV", na=False)).astype(int)
merged["sex_m"] = (merged["sex"]=="Male").astype(int)
merged["vi_z"] = (merged["ssgsea_vi"]-merged["ssgsea_vi"].mean())/merged["ssgsea_vi"].std()
merged["vi_ordinal"] = merged["subtype"].map({"Low VI":0,"Intermediate VI":1,"High VI":2})
merged["high_vs_low"] = np.where(merged["subtype"]=="High VI",1,np.where(merged["subtype"]=="Low VI",0,np.nan))

hdr("0. cohort shape")
print("clusters rows:", len(clusters), " cols:", list(clusters.columns))
print("cluster value counts:", clusters["cluster"].value_counts().to_dict())
print("subtype counts (all):", merged["subtype"].value_counts().to_dict())
print("mvi values:", merged["mvi"].value_counts(dropna=False).to_dict())
mvi_sub = merged[merged["mvi"].isin([0,1])].copy()
print("mvi-subset rows:", len(mvi_sub), " subtype counts:", mvi_sub["subtype"].value_counts().to_dict())

hdr("1. core gene counts")
up_genes = pd.read_csv("results/vascular_invasion_up_genes.csv"); dn_genes = pd.read_csv("results/vascular_invasion_dn_genes.csv")
print("up:", len(up_genes), " dn:", len(dn_genes))

hdr("2. subtype n + median VI")
for s in SUB:
    g = merged[merged["subtype"]==s]
    print(f"  {s:18} n={len(g):3d}  median VI={g['ssgsea_vi'].median():.1f}  mean={g['ssgsea_vi'].mean():.1f}")
print("  median rounded:", [f"{merged[merged['subtype']==s]['ssgsea_vi'].median():.0f}" for s in SUB])

hdr("3. MVI rate + chi2 (mvi subset)")
rates = mvi_sub.groupby("subtype")["mvi"].mean()*100
pos = {s:int(mvi_sub[mvi_sub["subtype"]==s]["mvi"].sum()) for s in SUB}
cnt = {s:int(mvi_sub[mvi_sub["subtype"]==s]["mvi"].notna().sum()) for s in SUB}
tab = np.array([[pos[s], cnt[s]-pos[s]] for s in SUB])
chi2,p,_,_ = stats.chi2_contingency(tab)
print("  counts pos:", pos, " totals:", cnt)
print("  rates:", {s:f"{rates[s]:.1f}%" for s in SUB})
print(f"  chi2={chi2:.3f}  P={p:.4f}")

hdr("4. AUC + Mann-Whitney (mvi subset)")
auc = roc_auc_score(mvi_sub["mvi"], mvi_sub["ssgsea_vi"])
u,pmw = stats.mannwhitneyu(mvi_sub[mvi_sub["mvi"]==1]["ssgsea_vi"], mvi_sub[mvi_sub["mvi"]==0]["ssgsea_vi"])
print(f"  AUC={auc:.4f}  MW P={pmw:.3e}")

hdr("5. AFP median + KW")
afp_med = merged.groupby("subtype")["afp"].median()
grp=[merged[merged["subtype"]==s]["afp"].dropna().values for s in SUB]
_,pafp = stats.kruskal(*grp)
print("  medians:", {s:round(float(afp_med[s]),1) for s in SUB}, f" KW P={pafp:.3e}")

hdr("6. Stage III/IV + Child-Pugh A + Male + Age")
def pct_by(col, pred):
    out={}
    for s in SUB:
        g=merged[merged["subtype"]==s]
        out[s]=(int(pred(g).sum()), int(g[col].notna().sum() if col else pred(g).notna().sum()))
    return out
# stage
adv={s:(int(merged[merged["subtype"]==s]["advanced"].sum()), int(merged[merged["subtype"]==s]["ajcc_stage"].notna().sum())) for s in SUB}
tab=np.array([[adv[s][0], adv[s][1]-adv[s][0]] for s in SUB]); c2,padv,_,_=stats.chi2_contingency(tab)
print("  Stage III/IV:", {s:f"{adv[s][0]} ({adv[s][0]/adv[s][1]*100:.1f}%)" for s in SUB}, f" P={padv:.4f}")
cp={s:(int(merged[merged["subtype"]==s]["child_pugh"].eq('A').sum()), int(merged[merged["subtype"]==s]["child_pugh"].notna().sum())) for s in SUB}
tab=np.array([[cp[s][0], cp[s][1]-cp[s][0]] for s in SUB]); c2,pcp,_,_=stats.chi2_contingency(tab)
print("  Child-Pugh A:", {s:f"{cp[s][0]} ({cp[s][0]/cp[s][1]*100:.1f}%)" for s in SUB}, f" P={pcp:.4f}")
mal={s:(int(merged[merged["subtype"]==s]["sex_m"].sum()), int(merged[merged["subtype"]==s]["sex_m"].notna().sum())) for s in SUB}
tab=np.array([[mal[s][0], mal[s][1]-mal[s][0]] for s in SUB]); c2,pmal,_,_=stats.chi2_contingency(tab)
print("  Male:", {s:f"{mal[s][0]} ({mal[s][0]/mal[s][1]*100:.1f}%)" for s in SUB}, f" P={pmal:.4f}")
age_med=merged.groupby("subtype")["age"].median()
grp=[merged[merged["subtype"]==s]["age"].dropna().values for s in SUB]
_,page=stats.kruskal(*grp)
print("  Age medians:", {s:round(float(age_med[s]),1) for s in SUB}, f" P={page:.4f}")

hdr("7. Survival log-rank (3 group)")
surv=pd.read_csv("results/survival_summary.csv")
print("  DFS logrank P=", surv["DFS_logrank_p_3group"][0], " OS logrank P=", surv["OS_logrank_p_3group"][0])
print("  DFS cohort n per subtype:")
d=merged.dropna(subset=["dfs_months","recurrence"]); d=d[d["dfs_months"]>0]
print("   DFS n:", len(d), d["subtype"].value_counts().to_dict())
d2=merged.dropna(subset=["os_months","death"]); d2=d2[d2["os_months"]>0]
print("   OS n:", len(d2), d2["subtype"].value_counts().to_dict())
ext=pd.read_csv("results/gse14520_vi_recurrence.csv").dropna(subset=["ssgsea_vi","recurr","recurr_months"]); ext=ext[ext["recurr_months"]>0]
print("  GSE14520 n=", len(ext))

hdr("8. Cox multivariate (lifelines)")
def mv(time, event, vars_):
    d=merged.dropna(subset=[time,event]).copy(); d=d[d[time]>0]; d=d.dropna(subset=vars_).copy()
    cph=CoxPHFitter(); cph.fit(d[[time,event]+vars_], duration_col=time, event_col=event)
    return cph
def show(tag, cph, var):
    s=cph.summary.loc[var]
    print(f"  {tag:22} {var:14} HR={s['exp(coef)']:.3f} ({s['exp(coef) lower 95%']:.3f}-{s['exp(coef) upper 95%']:.3f}) P={s['p']:.4f}")
for lab,time,ev in [("DFS","dfs_months","recurrence"),("OS","os_months","death")]:
    c=mv(time,ev,["vi_z","age","sex_m","advanced"]); show(lab+" VIscore/SD", c,"vi_z"); show(lab+" Age",c,"age"); show(lab+" Sex",c,"sex_m"); show(lab+" AdvStage",c,"advanced")
    ch=mv(time,ev,["high_vs_low","age","sex_m","advanced"]); show(lab+" HighVsLow", ch,"high_vs_low")
    cl=mv(time,ev,["vi_ordinal","age","sex_m","advanced"]); show(lab+" PerLevel", cl,"vi_ordinal")

hdr("9. GSE14520 Cox")
extc=ext.copy(); extc["vi_z"]=(extc["ssgsea_vi"]-extc["ssgsea_vi"].mean())/extc["ssgsea_vi"].std()
cph=CoxPHFitter(); cph.fit(extc[["recurr_months","recurr","vi_z"]],duration_col="recurr_months",event_col="recurr")
s=cph.summary.loc["vi_z"]; print(f"  HR/SD={s['exp(coef)']:.3f} ({s['exp(coef) lower 95%']:.3f}-{s['exp(coef) upper 95%']:.3f}) P={s['p']:.4f}")

hdr("10. Single-cell VI by site (mean) + cluster2")
mh=pd.read_csv("results/malignant_hepatocyte_vi_score.csv")
site=mh.groupby("site")["VI_score"].mean()
print("  mean VI by site:", {k:round(float(v),3) for k,v in site.items()})
mc=pd.read_csv("results/malig_cluster_vi_score.csv")
c2=mc[mc["leiden"]==2]
print(f"  cluster2 VI mean={c2['VI_score'].mean():.3f}  PVTT%={ (c2['site']=='PVTT').mean()*100:.1f}%  overallPVTT%={(mc['site']=='PVTT').mean()*100:.1f}%  fold={(c2['site']=='PVTT').mean()/(mc['site']=='PVTT').mean():.2f}")

hdr("11. Trajectory pattern counts")
for name in ["up","dn"]:
    pat=pd.read_csv(f"results/trajectory_{name}_patterns.csv", index_col=0)["pattern"]
    print(f"  {name}:", dict(sorted(Counter(pat).items())))
