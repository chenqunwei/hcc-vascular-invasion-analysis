import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

print("=" * 60)
print("功能富集分析（GO + KEGG）：血管侵犯上调/下调基因集")
print("=" * 60)

up = pd.read_csv("results/vascular_invasion_up_genes_symbol.csv")["symbol"].dropna().tolist()
dn = pd.read_csv("results/vascular_invasion_dn_genes_symbol.csv")["symbol"].dropna().tolist()
print(f"上调基因: {len(up)}, 下调基因: {len(dn)}")

LIBS = ["GO_Biological_Process_2023", "GO_Cellular_Component_2023",
        "GO_Molecular_Function_2023", "KEGG_2021_Human"]

def run_enrichr(gene_list, name):
    """对基因列表做多库富集，返回合并结果。"""
    all_res = []
    for lib in LIBS:
        try:
            enr = gp.enrichr(gene_list=gene_list, gene_sets=[lib], organism="human", outdir=None)
            res = enr.results.copy()
            res["Library"] = lib
            all_res.append(res)
        except Exception as e:
            print(f"  {lib} 失败: {e}")
    if all_res:
        return pd.concat(all_res, ignore_index=True)
    return pd.DataFrame()

# ---------- 上调基因富集 ----------
print("\n[1] 上调基因富集 ...")
up_res = run_enrichr(up, "up")
up_res.to_csv("results/enrichment_up.csv", index=False)
print(f"  上调富集结果: {len(up_res)} 条")
# top 显著
top_up = up_res[up_res["Adjusted P-value"] < 0.05].sort_values("Adjusted P-value")
print(f"  上调显著（FDR<0.05）: {len(top_up)} 条")
print("  Top 10 GO BP + KEGG:")
bp_up = top_up[top_up["Library"] == "GO_Biological_Process_2023"].head(10)
print(bp_up[["Term", "Overlap", "Adjusted P-value"]].to_string(index=False))

# ---------- 下调基因富集 ----------
print("\n[2] 下调基因富集 ...")
dn_res = run_enrichr(dn, "dn")
dn_res.to_csv("results/enrichment_dn.csv", index=False)
print(f"  下调富集结果: {len(dn_res)} 条")
top_dn = dn_res[dn_res["Adjusted P-value"] < 0.05].sort_values("Adjusted P-value")
print(f"  下调显著（FDR<0.05）: {len(top_dn)} 条")
print("  Top 10 GO BP + KEGG:")
bp_dn = top_dn[top_dn["Library"] == "GO_Biological_Process_2023"].head(10)
print(bp_dn[["Term", "Overlap", "Adjusted P-value"]].to_string(index=False))

# ---------- 画富集气泡图 ----------
def plot_bubble(res, title, out_png, n=15):
    res = res[res["Adjusted P-value"] < 0.05].copy()
    res = res.sort_values("Adjusted P-value").head(n)
    if len(res) == 0:
        print(f"  {title}: 无显著富集，跳过绘图")
        return
    # 解析 Overlap "x/y" -> gene ratio
    res["overlap_n"] = res["Overlap"].apply(lambda s: int(s.split("/")[0]))
    res["total_n"] = res["Overlap"].apply(lambda s: int(s.split("/")[1]))
    res["gene_ratio"] = res["overlap_n"] / res["total_n"]
    res = res.sort_values("gene_ratio", ascending=True)
    res["-log10padj"] = -np.log10(res["Adjusted P-value"].clip(lower=1e-30))

    fig, ax = plt.subplots(figsize=(9, max(4, n * 0.38)))
    sc = ax.scatter(res["gene_ratio"], range(len(res)), s=res["overlap_n"] * 15,
                    c=res["-log10padj"], cmap="Reds", alpha=0.8, edgecolors="black", linewidths=0.3)
    ax.set_yticks(range(len(res)))
    ax.set_yticklabels(res["Term"], fontsize=9)
    ax.set_xlabel("Gene ratio", fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    cbar = fig.colorbar(sc, ax=ax, pad=0.01)
    cbar.set_label("-log10(adjusted P)", fontsize=9)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  已保存 {out_png}")

print("\n[3] 绘制富集气泡图 ...")
plot_bubble(up_res[up_res["Library"] == "GO_Biological_Process_2023"],
            "Up-regulated in vascular invasion: GO BP", "results/enrich_up_GOBP.png")
plot_bubble(dn_res[dn_res["Library"] == "GO_Biological_Process_2023"],
            "Down-regulated in vascular invasion: GO BP", "results/enrich_dn_GOBP.png")
plot_bubble(up_res[up_res["Library"] == "KEGG_2021_Human"],
            "Up-regulated: KEGG", "results/enrich_up_KEGG.png")
plot_bubble(dn_res[dn_res["Library"] == "KEGG_2021_Human"],
            "Down-regulated: KEGG", "results/enrich_dn_KEGG.png")

print("\n功能富集完成。")
