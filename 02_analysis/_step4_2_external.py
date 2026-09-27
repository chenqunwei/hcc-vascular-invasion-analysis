# -*- coding: utf-8 -*-
"""Step 4.2 复发预测外部验证（GSE14520）
定位：证明「血管侵犯分子程序(VI score)预测复发」在独立队列(GSE14520)中可跨数据集复现。
方法：每个平台(GPL3921/GPL571)内部算 VI score 后合并，对复发做 KM + Cox。
VI score 定义与 Part A 下游一致：**ssGSEA 为主**（es_up − es_dn），mean z-score 为验证版。
数据：series matrix 已含 RMA 归一化表达值；临床复发标注来自 Extra_Supplement。
修复记录：
  2026-09-16a) GPL annot 列名精确等值匹配（'gene symbol' 不被 'unigene symbol' 覆盖）；
  2026-09-16b) series matrix GSM/probe ID strip 双引号；
  2026-09-16c) /// 多基因拆分；
  2026-09-16d) **VI score 主方法改用 ssGSEA**（与 Step 2.3/2.4 下游一致），mean z-score 降为验证版。
"""
import warnings; warnings.filterwarnings('ignore')
import gzip, os, re
import numpy as np, pandas as pd
import gseapy as gp
from scipy import stats
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter


def load_series_matrix(path):
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    sample_ids = None
    table_start = None
    for i, line in enumerate(lines):
        line = line.rstrip('\n')
        if line.startswith('!Sample_geo_accession'):
            sample_ids = [x.strip().strip('"') for x in line.split('\t')[1:]]
        if line.startswith('!series_matrix_table_begin'):
            table_start = i
            break
    data_rows = []
    for line in lines[table_start + 2:]:
        line = line.rstrip('\n')
        if line.startswith('!series_matrix_table_end'):
            break
        parts = line.split('\t')
        probe = parts[0].strip().strip('"')
        try:
            vals = [float(x.strip().strip('"')) for x in parts[1:]]
        except ValueError:
            continue
        if len(vals) != len(sample_ids):
            continue
        data_rows.append((probe, vals))
    expr = pd.DataFrame({p: v for p, v in data_rows}, index=sample_ids).T
    return expr, sample_ids


def load_gpl_annot(path):
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        txt = f.read()
    m = re.search(r'!platform_table_begin\n(.*?)!platform_table_end', txt, re.S)
    if not m:
        return {}
    table = m.group(1)
    lines = table.strip().split('\n')
    header = lines[0].split('\t')
    id_col = sym_col = None
    for i, h in enumerate(header):
        hl = h.strip().lower()
        if hl == 'id':
            id_col = i
        elif hl == 'gene symbol':
            sym_col = i
    if sym_col is None:
        for i, h in enumerate(header):
            if h.strip().lower() == 'symbol':
                sym_col = i
                break
    if id_col is None or sym_col is None:
        return {}
    mapping = {}
    for line in lines[1:]:
        parts = line.split('\t')
        if len(parts) <= max(id_col, sym_col):
            continue
        pid = parts[id_col].strip().strip('"')
        sym_raw = parts[sym_col].strip().strip('"')
        if not pid or not sym_raw or sym_raw in ('---', 'NA', ''):
            continue
        for s in sym_raw.split('///'):
            s = s.strip()
            if s:
                mapping.setdefault(pid, []).append(s)
    return mapping


# ---------- 1. 血管侵犯基因集（symbol）----------
print('[1] 加载血管侵犯基因集 ...')
up_sym = pd.read_csv('results/vascular_invasion_up_genes_symbol.csv')['symbol'].dropna().astype(str).tolist()
dn_sym = pd.read_csv('results/vascular_invasion_dn_genes_symbol.csv')['symbol'].dropna().astype(str).tolist()
print(f'    up={len(up_sym)}, dn={len(dn_sym)}')

# ---------- 2. 逐平台算 VI score（ssGSEA 主 + mean z-score 验证）----------
all_vi = []
for gpl in ['GPL3921', 'GPL571']:
    mat_path = f'data/geo/GSE14520-{gpl}_series_matrix.txt.gz'
    annot_path = f'data/geo/{gpl}.annot.gz'
    if not os.path.exists(mat_path):
        continue
    print(f'[2] 处理 {gpl} ...')
    expr, samples = load_series_matrix(mat_path)
    print(f'    表达矩阵: {expr.shape[0]} 探针 × {expr.shape[1]} 样本')

    mapping = load_gpl_annot(annot_path) if os.path.exists(annot_path) else {}
    print(f'    GPL 注释探针数: {len(mapping)}')
    rows = [(p, s) for p in expr.index for s in mapping.get(p, [])]
    probe_symbol = pd.DataFrame(rows, columns=['probe', 'symbol'])
    expr = expr.merge(probe_symbol, left_index=True, right_on='probe', how='inner')
    expr = expr.drop(columns='probe').groupby('symbol').mean()
    print(f'    映射后基因数: {expr.shape[0]}')

    up_hit = [g for g in up_sym if g in expr.index]
    dn_hit = [g for g in dn_sym if g in expr.index]
    print(f'    up 命中 {len(up_hit)}/{len(up_sym)}, dn 命中 {len(dn_hit)}/{len(dn_sym)}')

    # --- (A) ssGSEA 主方法 ---
    gene_sets = {"VI_up": up_hit, "VI_dn": dn_hit}
    ss = gp.ssgsea(data=expr, gene_sets=gene_sets, outdir=None,
                   min_size=5, max_size=10000, sample_norm_method="rank",
                   no_plot=True, threads=1)
    res = ss.res2d.copy()
    if "Term" in res.columns:
        es_mat = res.pivot(index="Name", columns="Term", values="ES")
    else:
        es_mat = res.set_index("Name")
    ssgsea_vi = es_mat["VI_up"] - es_mat["VI_dn"]

    # --- (B) mean z-score 验证版 ---
    z = (expr - expr.mean(axis=1).values.reshape(-1, 1)) / expr.std(axis=1).values.reshape(-1, 1)
    z = z.replace([np.inf, -np.inf], 0).fillna(0)
    meanz_vi = z.loc[up_hit].mean(axis=0) - z.loc[dn_hit].mean(axis=0)

    for gsm in samples:
        all_vi.append((gsm, float(ssgsea_vi[gsm]), float(meanz_vi[gsm])))
    print(f'    {gpl} VI score 完成（ssGSEA 主 + mean-z 验证），样本 {len(samples)}')

vi_df = pd.DataFrame(all_vi, columns=['gsm', 'ssgsea_vi', 'meanz_vi'])
print(f'\n[3] 合并两平台: {len(vi_df)} 样本, ssGSEA 非 NaN {vi_df.ssgsea_vi.notna().sum()}')

# ssGSEA vs mean-z 一致性（GSE14520 上）
r, rp = stats.spearmanr(vi_df['ssgsea_vi'], vi_df['meanz_vi'])
print(f'    ssGSEA vs mean-z Spearman ρ = {r:.3f} (p={rp:.2e})')

# ---------- 3. 合并临床复发标注 ----------
print('[4] 合并临床复发标注 ...')
clin = pd.read_csv('data/geo/GSE14520_Extra_Supplement.txt.gz', sep='\t', compression='gzip')
clin = clin.rename(columns={'Affy_GSM': 'gsm'})
clin_tumor = clin[clin['Tissue Type'] == 'Tumor'].copy()
clin_tumor['recurr'] = pd.to_numeric(clin_tumor['Recurr status'], errors='coerce')
clin_tumor['recurr_months'] = pd.to_numeric(clin_tumor['Recurr months'], errors='coerce')
clin_tumor = clin_tumor.dropna(subset=['recurr', 'recurr_months'])
print(f'    肿瘤样本有复发标注: {len(clin_tumor)} 例')

merged = vi_df.merge(clin_tumor[['gsm', 'recurr', 'recurr_months']], on='gsm', how='inner')
merged = merged.dropna(subset=['ssgsea_vi'])
print(f'    VI score 与复发标注匹配: {len(merged)} 例')
print(f'    复发分布: {merged.recurr.value_counts().to_dict()}')

# ---------- 4. KM + Cox（主 = ssGSEA）----------
def cox_km(df, vi_col, tag):
    df = df.copy()
    df['vi_group'] = (df[vi_col] > df[vi_col].median()).astype(int)
    kmf_hi = KaplanMeierFitter(); kmf_lo = KaplanMeierFitter()
    hi = df[df['vi_group'] == 1]; lo = df[df['vi_group'] == 0]
    kmf_hi.fit(hi['recurr_months'], event_observed=hi['recurr'])
    kmf_lo.fit(lo['recurr_months'], event_observed=lo['recurr'])
    lr = logrank_test(lo['recurr_months'], hi['recurr_months'],
                      event_observed_A=lo['recurr'], event_observed_B=hi['recurr'])
    df['vi_z'] = (df[vi_col] - df[vi_col].mean()) / df[vi_col].std()
    cph = CoxPHFitter()
    cph.fit(df[['recurr_months', 'recurr', 'vi_z']], duration_col='recurr_months', event_col='recurr')
    hr = np.exp(cph.params_['vi_z'])
    ci = np.exp(cph.confidence_intervals_.loc['vi_z'].values)
    p = cph.summary.loc['vi_z', 'p']
    print(f'    [{tag}] KM log-rank p={lr.p_value:.4f} | Cox HR={hr:.2f} ({ci[0]:.2f}-{ci[1]:.2f}), p={p:.4f}')
    return kmf_hi, kmf_lo, hi, lo, lr, hr, ci, p

print('\n[5] 复发预测 KM + Cox ...')
kmf_hi, kmf_lo, hi, lo, lr, hr, ci, p = cox_km(merged, 'ssgsea_vi', 'ssGSEA 主')
cox_km(merged, 'meanz_vi', 'mean-z 验证')

# 复发组 vs 未复发组 ssGSEA 均值
print(f'    复发组 ssgsea_vi 均值 {merged[merged.recurr==1].ssgsea_vi.mean():.3f} '
      f'vs 未复发组 {merged[merged.recurr==0].ssgsea_vi.mean():.3f}')

# ---------- 5. 保存 + 画图 ----------
merged.to_csv('results/gse14520_vi_recurrence.csv', index=False)
print(f'\n[6] 已保存 results/gse14520_vi_recurrence.csv')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(4.5, 3.5))
kmf_hi.plot_survival_function(ax=ax, label=f'High VI (n={len(hi)})', color='#D85A30')
kmf_lo.plot_survival_function(ax=ax, label=f'Low VI (n={len(lo)})', color='#1D9E75')
ax.set_xlabel('Recurrence-free months', fontsize=9)
ax.set_ylabel('DFS probability', fontsize=9)
ax.set_title(f'GSE14520 recurrence (ssGSEA VI, log-rank P = {lr.p_value:.3f})', fontsize=9)
ax.legend(fontsize=7, frameon=False)
plt.tight_layout()
fig.savefig('results/FigS_gse14520_recurrence_km.png', dpi=300, bbox_inches='tight')
plt.close(fig)
print('    已保存 results/FigS_gse14520_recurrence_km.png')
print('\n完成 Step 4.2 复发预测外部验证（ssGSEA 主方法）。')
