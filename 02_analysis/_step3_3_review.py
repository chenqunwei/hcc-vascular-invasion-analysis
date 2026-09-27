# -*- coding: utf-8 -*-
"""Step 3.3 方法学重新审核 —— 核查重做方案的可行性
核心问题：当前 core 程序 = Normal->Tumor 显著，混杂了"正常 vs 恶性"差异。
重做方向：仅在恶性肝细胞(Tumor/PVTT/Lymph)内部，用"Tumor 中已表达比例"定义 core，
         用"PVTT vs Tumor 差异"定义 stage。
本脚本核查：血管侵犯基因在 Tumor 中的表达比例分布，能否支撑"早期已启动 vs 晚期才启动"的区分。
"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, scanpy as sc

adata = sc.read_h5ad('results/gse149614_processed_qc.h5ad')
heps = adata[adata.obs['celltype'] == 'Hepatocyte'].copy()
malig = heps[heps.obs['site'].isin(['Tumor', 'PVTT', 'Lymph'])].copy()
print('恶性肝细胞总数:', malig.shape[0])
print('site 分布:', malig.obs['site'].value_counts().to_dict())
print('Normal 肝细胞(应排除):', (heps.obs['site']=='Normal').sum())

up_sym = pd.read_csv('results/vascular_invasion_up_genes_symbol.csv')['symbol'].dropna().tolist()
dn_sym = pd.read_csv('results/vascular_invasion_dn_genes_symbol.csv')['symbol'].dropna().tolist()
up_hit = [g for g in up_sym if g in malig.var_names]
dn_hit = [g for g in dn_sym if g in malig.var_names]
print(f'\nup 命中 {len(up_hit)}/{len(up_sym)}, dn 命中 {len(dn_hit)}/{len(dn_sym)}')

tumor_cells = malig[malig.obs['site'] == 'Tumor']
pvtt_cells = malig[malig.obs['site'] == 'PVTT']

def detect_frac(X):
    return float((np.asarray(X) > 0).mean())

def gene_fracs(genes, cells):
    """返回每个基因在给定细胞中的表达比例(>0 的细胞占比)。"""
    fracs = {}
    for g in genes:
        gi = list(cells.var_names).index(g)
        fracs[g] = detect_frac(cells[:, gi].X.toarray().ravel())
    return fracs

print('\n计算 up 基因在 Tumor / PVTT 中的表达比例 ...')
frac_t = gene_fracs(up_hit, tumor_cells)
frac_p = gene_fracs(up_hit, pvtt_cells)

df = pd.DataFrame({'gene': up_hit,
                   'frac_tumor': [frac_t[g] for g in up_hit],
                   'frac_pvtt': [frac_p[g] for g in up_hit]})
print('\nup 基因表达比例分布:')
print(df[['frac_tumor', 'frac_pvtt']].describe().round(3).to_string())

# 尝试区分
core_like = df[df.frac_tumor > 0.5]
stage_like = df[(df.frac_tumor < 0.3) & (df.frac_pvtt > 0.5)]
mid = df[~df.index.isin(core_like.index) & ~df.index.isin(stage_like.index)]
print(f'\n"早期已启动"(Tumor 表达比例>50%): {len(core_like)} 个')
print(f'"晚期才启动"(Tumor<30% 且 PVTT>50%): {len(stage_like)} 个')
print(f'"中间"(未分类): {len(mid)} 个')

# 展示几个代表基因
print('\n代表基因示例:')
print('  早期已启动(前5):', core_like.sort_values('frac_tumor', ascending=False)['gene'].head(5).tolist())
print('  晚期才启动(前5):', stage_like.sort_values('frac_pvtt', ascending=False)['gene'].head(5).tolist())
