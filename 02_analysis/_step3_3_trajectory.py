# -*- coding: utf-8 -*-
"""Step 3.3（重做版 v2）：血管侵犯基因的四层次轨迹模式聚类
保持四层次框架（Normal->Tumor->PVTT->Lymph），对血管侵犯 up/dn 基因分别
做"表达轨迹"模式聚类，展示基因激活时序的连续分布（不再做 core/stage 二分）。
"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

SITES = ['Normal', 'Tumor', 'PVTT', 'Lymph']

print('[1] 加载 h5ad ...')
adata = sc.read_h5ad('results/gse149614_processed_qc.h5ad')
heps = adata[adata.obs['celltype'] == 'Hepatocyte'].copy()
print('    肝细胞:', heps.shape[0], '| site:', heps.obs['site'].value_counts().to_dict())

up = pd.read_csv('results/vascular_invasion_up_genes_symbol.csv')['symbol'].dropna().tolist()
dn = pd.read_csv('results/vascular_invasion_dn_genes_symbol.csv')['symbol'].dropna().tolist()
up_hit = [g for g in up if g in heps.var_names]
dn_hit = [g for g in dn if g in heps.var_names]
print(f'    up 命中 {len(up_hit)}/{len(up)}, dn 命中 {len(dn_hit)}/{len(dn)}')

print('[2] 计算每基因四 site 平均表达 ...')
def site_mean_matrix(genes):
    var_idx = {g: i for i, g in enumerate(heps.var_names)}
    mat = np.zeros((len(genes), 4))
    site_masks = {s: (heps.obs['site'] == s).values for s in SITES}
    for gi, g in enumerate(genes):
        col = heps[:, var_idx[g]].X
        if hasattr(col, 'toarray'):
            col = col.toarray().ravel()
        else:
            col = np.asarray(col).ravel()
        for si, s in enumerate(SITES):
            mat[gi, si] = col[site_masks[s]].mean()
    return mat

mat_up = site_mean_matrix(up_hit)
mat_dn = site_mean_matrix(dn_hit)

def zscore_rows(mat):
    mu = mat.mean(axis=1, keepdims=True)
    sd = mat.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    return (mat - mu) / sd

z_up = zscore_rows(mat_up)
z_dn = zscore_rows(mat_dn)

print('[3] 选择最优模式数 k（silhouette）...')
def choose_k(zmat, name, k_range=(3, 7)):
    best_k, best_s = None, -1
    for k in range(*k_range):
        km = KMeans(n_clusters=k, n_init=20, random_state=0).fit(zmat)
        s = silhouette_score(zmat, km.labels_)
        print(f'    {name} k={k}: silhouette={s:.3f}')
        if s > best_s:
            best_k, best_s = k, s
    return best_k

k_up = choose_k(z_up, 'up')
k_dn = choose_k(z_dn, 'dn')

print('[4] 聚类 + 按轨迹形状排序模式 ...')
def cluster_and_order(zmat, k, name):
    km = KMeans(n_clusters=k, n_init=50, random_state=0).fit(zmat)
    labels = km.labels_
    # 每个 cluster 的平均轨迹
    mean_traj = np.array([zmat[labels == c].mean(axis=0) for c in range(k)])
    # 按"Tumor 位置的 z-score"排序（up: 越大越早激活；dn: 越小越早下降）
    # 统一用 Tumor 列值排序
    order = np.argsort(mean_traj[:, 1])[::-1]  # Tumor 列从高到低
    # 重映射 label -> 排序后的序号(0=最早)
    remap = {old: new for new, old in enumerate(order)}
    new_labels = np.array([remap[l] for l in labels])
    return new_labels, mean_traj[order]

labels_up, traj_up = cluster_and_order(z_up, k_up, 'up')
labels_dn, traj_dn = cluster_and_order(z_dn, k_dn, 'dn')

print('\n    up 模式平均轨迹（行=模式, 列=Normal/Tumor/PVTT/Lymph z-score）:')
print(np.round(traj_up, 2))
print('    dn 模式平均轨迹:')
print(np.round(traj_dn, 2))

print('[5] 保存分类结果 ...')
df_up = pd.DataFrame({'gene': up_hit, 'pattern': labels_up})
df_dn = pd.DataFrame({'gene': dn_hit, 'pattern': labels_dn})
df_up.to_csv('results/trajectory_up_patterns.csv', index=False)
df_dn.to_csv('results/trajectory_dn_patterns.csv', index=False)

print('[6] 出图 ...')
# 配色：模式用连续色
import matplotlib.cm as cm
def plot_trajectory_heatmap(zmat, labels, genes, name, cmap_name):
    # 按 label 分组，组内按第一主方向排序
    order_idx = []
    for c in sorted(set(labels)):
        idx = np.where(labels == c)[0]
        # 组内按 Tumor->PVTT 变化量排序
        key = zmat[idx, 2] - zmat[idx, 1]
        idx = idx[np.argsort(key)[::-1]]
        order_idx.extend(idx)
    zmat_ord = zmat[order_idx]
    labels_ord = labels[order_idx]
    k = len(set(labels))
    colors = cm.viridis(np.linspace(0, 0.85, k))
    row_colors = [colors[c] for c in labels_ord]

    fig, ax = plt.subplots(figsize=(2.6, 7))
    im = ax.imshow(zmat_ord, aspect='auto', cmap=cmap_name, vmin=-2, vmax=2)
    ax.set_xticks(range(4)); ax.set_xticklabels(SITES, fontsize=8, rotation=30, ha='right')
    ax.set_yticks([])
    ax.set_title(f'{name} ({len(genes)} genes)', fontsize=9)
    # 右侧模式颜色条
    ax2 = fig.add_axes([0.93, 0.15, 0.03, 0.7])
    ax2.imshow(np.array(row_colors).reshape(-1, 1, 3) if False else [[c] for c in row_colors],
               aspect='auto', interpolation='nearest')
    ax2.set_yticks([]); ax2.set_xticks([])
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.04)
    cbar.ax.tick_params(labelsize=7)
    plt.tight_layout()
    fig.savefig(f'results/FigS_trajectory_heatmap_{name.lower()}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

plot_trajectory_heatmap(z_up, labels_up, up_hit, 'Up', 'Reds')
plot_trajectory_heatmap(z_dn, labels_dn, dn_hit, 'Dn', 'Blues')

# 平均轨迹线图
def plot_mean_trajectory(traj, name, title):
    fig, ax = plt.subplots(figsize=(3.4, 2.8))
    k = traj.shape[0]
    colors = cm.viridis(np.linspace(0, 0.85, k))
    x = np.arange(4)
    for c in range(k):
        ax.plot(x, traj[c], marker='o', markersize=4, color=colors[c],
                label=f'P{c+1} (n={int((labels_up==c).sum()) if name=="Up" else int((labels_dn==c).sum())})', linewidth=1.5)
    ax.axhline(0, color='grey', lw=0.6, ls='--')
    ax.set_xticks(x); ax.set_xticklabels(SITES, fontsize=8)
    ax.set_ylabel('Mean z-score', fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=6, loc='best', frameon=False)
    plt.tight_layout()
    fig.savefig(f'results/FigS_trajectory_mean_{name.lower()}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

plot_mean_trajectory(traj_up, 'Up', 'Up-regulated program trajectories')
plot_mean_trajectory(traj_dn, 'Dn', 'Down-regulated program trajectories')

print('\n完成。输出:')
print('  results/trajectory_up_patterns.csv / trajectory_dn_patterns.csv')
print('  results/FigS_trajectory_heatmap_up.png / dn.png')
print('  results/FigS_trajectory_mean_up.png / dn.png')
