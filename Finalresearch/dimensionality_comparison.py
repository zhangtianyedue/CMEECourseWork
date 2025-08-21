import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import param
import pandas as pd
import seaborn as sns

# ========== 无量纲化函数 ==========
def dimensionless_parameters():
    """
    返回无量纲化的参考尺度
    """
    R0 = 2.0      # 资源参考浓度
    C0 = 1.0      # 物种参考浓度  
    T0 = 1.0      # 时间参考尺度（小时）
    L0 = 1.0      # 空间参考尺度（patch间距）
    
    return R0, C0, T0, L0

def normalize_parameters(community_params, R0, C0, T0, L0):
    """
    将所有参数无量纲化
    """
    normalized_params = []
    
    for cp in community_params:
        norm_cp = cp.copy()
        
        # 无量纲化资源输入率: ρ̃ = ρ / (R0/T0)
        norm_cp['rho'] = cp['rho'] / (R0 / T0)
        
        # 无量纲化资源消耗率: ω̃ = ω * T0
        norm_cp['omega'] = cp['omega'] * T0
        
        # 无量纲化物种死亡率: m̃ = m * T0
        norm_cp['m'] = cp['m'] * T0
        
        # 无量纲化消耗矩阵: ũ = u * C0 * R0 * T0
        norm_cp['u'] = cp['u'] * C0 * R0 * T0
        
        # 泄漏张量保持比例形式
        norm_cp['l'] = cp['l']
        
        normalized_params.append(norm_cp)
    
    return normalized_params

def normalize_diffusion_coefficients(D0, D_species, T0, L0):
    """
    无量纲化扩散系数: D̃ = D * T0 / L0^2
    """
    D0_norm = D0 * T0 / (L0**2)
    D_species_norm = D_species * T0 / (L0**2)
    return D0_norm, D_species_norm

def normalize_initial_conditions(Y0, R0, C0, K, N_pool, M):
    """
    无量纲化初始条件
    """
    Y0_norm = Y0.copy()
    
    # 物种初始条件归一化: C* = C / C0
    for k in range(K):
        Y0_norm[k*N_pool:(k+1)*N_pool] /= C0
    
    # 资源初始条件归一化: R* = R / R0
    for k in range(K):
        Y0_norm[K*N_pool + k*M:K*N_pool + (k+1)*M] /= R0
    
    return Y0_norm

# 1. 支持不同空间维度的patch坐标和距离矩阵生成
def generate_coords_and_dist(K, dim):
    if dim == 1:
        coords = np.arange(K).reshape(-1, 1)
    elif dim == 2:
        side = int(np.round(K ** (1/2)))
        coords = np.array([[i, j] for i in range(side) for j in range(side)])
    elif dim == 3:
        # 对于3D，使用5x5x1的结构，保持25个patches但具有3D坐标
        side_2d = int(np.round(K ** (1/2)))
        coords = np.array([[i, j, 0] for i in range(side_2d) for j in range(side_2d)])
    else:
        raise ValueError('dim must be 1, 2, or 3')
    D_ij = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    return coords, D_ij

# 2. 无量纲化参数设置
R0, C0, T0, L0 = dimensionless_parameters()

K_1d = 9
K_2d = 9  # 5x5
K_3d = 9  # 5x5x1 (为了保持相同patch数量，3D使用5x5x1的扁平结构，但仍然是3D坐标)

# 无量纲化后的扩散系数（数值更稳定）
D0 = 1.0
D_species = 1.0
lambda_decay = 1.0  # 无量纲距离衰减参数
beta = 1.0
N_pool = 50  # 全局物种池
M = 10     # 每个patch资源数

def run_sim(K, dim, D0, D_species, with_diffusion, seed):
    np.random.seed(seed)
    
    # patch参数和初始条件（有物种池结构）
    patch_species_list = []
    N_k_list = []
    for k in range(K):
        N_k = np.random.randint(3, N_pool+1)
        N_k_list.append(N_k)
        patch_species = np.sort(np.random.choice(np.arange(N_pool), N_k, replace=False))
        patch_species_list.append(patch_species)
    
    # 原始参数设置
    community_params = []
    for k in range(K):
        N_k = N_k_list[k]
        patch_species = patch_species_list[k]
        u_k = param.modular_uptake(N_k, M, N_modules=1, s_ratio=10.0)
        l_k = param.generate_l_tensor(N_k, M, N_modules=2, s_ratio=10.0, λ=0.3)
        if k == 0:
            rho_k = np.ones(M) * 10.0
        elif k == K-1:
            rho_k = np.ones(M) * 5.0
        else:
            rho_k = np.ones(M) * (5.0 - (4.0 * k / (K-1)))
        omega_k = np.ones(M) * 0.05
        m_k = np.ones(N_k) * 0.09
        community_params.append({'N': N_k, 'M': M, 'u': u_k, 'l': l_k, 'rho': rho_k, 'omega': omega_k, 'm': m_k, 'species_idx': patch_species})
    
    # 无量纲化所有参数
    community_params = normalize_parameters(community_params, R0, C0, T0, L0)
    
    # 变量索引
    C_idx = []
    R_idx = []
    for k in range(K):
        C_idx.append((k*N_pool, (k+1)*N_pool))
    for k in range(K):
        R_idx.append((K*N_pool + k*M, K*N_pool + (k+1)*M))
    n_total = K*N_pool + K*M
    
    # 初始条件
    Y0 = np.zeros(n_total)
    for k, cp in enumerate(community_params):
        N_k = cp['N']
        patch_species = cp['species_idx']
        Y0[C_idx[k][0]:C_idx[k][1]][patch_species] = np.random.uniform(0.5, 1.5, N_k)
        Y0[R_idx[k][0]:R_idx[k][1]] = np.ones(M) * 2.0
    
    # 无量纲化初始条件
    Y0 = normalize_initial_conditions(Y0, R0, C0, K, N_pool, M)
    
    # 距离矩阵（无量纲化）
    coords, D_ij = generate_coords_and_dist(K, dim)
    # 无量纲化距离矩阵（除以参考长度）
    D_ij = D_ij / L0
    
    k_j = np.ones(K)
    A_kj = np.zeros((K, K))
    for k1 in range(K):
        for j in range(K):
            if k1 != j:
                A_kj[k1, j] = k_j[j]**beta * np.exp(-lambda_decay * D_ij[k1, j])
            else:
                A_kj[k1, j] = 0
    
    # 预处理：统计每个物种在哪些patch中出现
    species_in_patches = [[] for _ in range(N_pool)]
    for k, cp in enumerate(community_params):
        for sp in cp['species_idx']:
            species_in_patches[sp].append(k)
    
    def get_abund(sol):
        abund = np.zeros(N_pool)
        for k, cp in enumerate(community_params):
            patch_species = cp['species_idx']
            C_patch = sol.y[C_idx[k][0]:C_idx[k][1], -1]
            for idx, sp in enumerate(patch_species):
                abund[sp] += C_patch[sp]
        return abund
    
    def global_ode(t, y, D0, D_species, with_diffusion=True):
        # 无量纲化扩散系数
        D0_norm, D_species_norm = normalize_diffusion_coefficients(D0, D_species, T0, L0)
        
        dydt = np.zeros_like(y)
        for k, cp in enumerate(community_params):
            N_k = cp['N']
            patch_species = cp['species_idx']
            u, l, rho, omega, m = cp['u'], cp['l'], cp['rho'], cp['omega'], cp['m']
            C_patch = y[C_idx[k][0]:C_idx[k][1]]
            R = y[R_idx[k][0]:R_idx[k][1]]
            dCdt = np.zeros(N_k)
            dRdt = np.zeros(M)
            for i in range(N_k):
                dCdt[i] = sum(C_patch[patch_species[i]] * R[alpha] * u[i, alpha] * (1 - 0.3) for alpha in range(M)) - C_patch[patch_species[i]] * m[i]
            for alpha in range(M):
                dRdt[alpha] = rho[alpha] - R[alpha] * omega[alpha]
                dRdt[alpha] -= sum(C_patch[patch_species[i]] * R[alpha] * u[i, alpha] for i in range(N_k))
                dRdt[alpha] += sum(sum(C_patch[patch_species[i]] * R[beta] * u[i, beta] * l[i, beta, alpha] for beta in range(M)) for i in range(N_k))
            for idx, sp in enumerate(patch_species):
                dydt[C_idx[k][0] + sp] = dCdt[idx]
            dydt[R_idx[k][0]:R_idx[k][1]] = dRdt
        # 物种扩散项（A_kj全局扩散）
        if with_diffusion and D_species_norm > 0:
            for i in range(N_pool):
                patch_list = species_in_patches[i]
                if len(patch_list) < 2:
                    continue
                for k1 in patch_list:
                    for k2 in patch_list:
                        if k1 == k2:
                            continue
                        dydt[C_idx[k1][0] + i] += D_species_norm * A_kj[k1, k2] * (y[C_idx[k2][0] + i] - y[C_idx[k1][0] + i])
        # 资源扩散项（A_kj全局扩散）
        if with_diffusion and D0_norm > 0:
            for alpha in range(M):
                for i in range(K):
                    for j in range(K):
                        if i == j:
                            continue
                        dydt[R_idx[i][0] + alpha] += D0_norm * A_kj[i, j] * (y[R_idx[j][0] + alpha] - y[R_idx[i][0] + alpha])
        return dydt
    
    t_span = (0, 500)
    t_eval = np.linspace(*t_span, 500)
    sol = solve_ivp(lambda t, y: global_ode(t, y, D0, D_species, with_diffusion=with_diffusion), t_span, Y0, t_eval=t_eval, method='LSODA', rtol=1e-5, atol=1e-7)
    abund = get_abund(sol)
    return abund

# 3. 进行50次模拟
n_simulations = 50
print(f"开始进行{n_simulations}次模拟...")

# 存储所有结果
results = {
    '1D_no_diff': [],
    '1D_with_diff': [],
    '2D_no_diff': [],
    '2D_with_diff': [],
    '3D_no_diff': [],
    '3D_with_diff': []
}

# 进行模拟
for sim in range(n_simulations):
    print(f"模拟进度: {sim+1}/{n_simulations}")
    
    # 1D空间
    abund_1d_no = run_sim(K_1d, 1, D0=0, D_species=0, with_diffusion=False, seed=42+sim)
    abund_1d_with = run_sim(K_1d, 1, D0=1, D_species=1, with_diffusion=True, seed=42+sim)
    
    # 2D空间
    abund_2d_no = run_sim(K_2d, 2, D0=0, D_species=0, with_diffusion=False, seed=42+sim)
    abund_2d_with = run_sim(K_2d, 2, D0=1, D_species=1, with_diffusion=True, seed=42+sim)
    
    # 3D空间
    abund_3d_no = run_sim(K_3d, 3, D0=0, D_species=0, with_diffusion=False, seed=42+sim)
    abund_3d_with = run_sim(K_3d, 3, D0=1, D_species=1, with_diffusion=True, seed=42+sim)
    
    # 存储结果
    results['1D_no_diff'].append(abund_1d_no)
    results['1D_with_diff'].append(abund_1d_with)
    results['2D_no_diff'].append(abund_2d_no)
    results['2D_with_diff'].append(abund_2d_with)
    results['3D_no_diff'].append(abund_3d_no)
    results['3D_with_diff'].append(abund_3d_with)

print("所有模拟完成！")

# 4. 创建可视化函数
def create_rank_abundance_plot(results, title, save_name=None):
    """创建rank-abundance图，包含箱线图和散点"""
    
    # 准备数据
    plot_data = []
    
    for condition, abund_list in results.items():
        for sim_idx, abund in enumerate(abund_list):
            # 计算rank-abundance
            sorted_abund = np.sort(abund)[::-1]  # 从大到小排序
            ranks = np.arange(1, len(sorted_abund)+1)
            
            # 只取前10个rank（避免图太复杂）
            for rank, abundance in zip(ranks[:10], sorted_abund[:10]):
                plot_data.append({
                    'Condition': condition,
                    'Rank': rank,
                    'Abundance': abundance,
                    'Simulation': sim_idx
                })
    
    df = pd.DataFrame(plot_data)
    
    # 创建图
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # 颜色设置
    colors = ['#E74C3C', '#3498DB', '#2ECC71', '#F39C12', '#9B59B6', '#1ABC9C']
    
    # 为每个rank创建子图
    for rank in range(1, 11):
        row = (rank - 1) // 3
        col = (rank - 1) % 3
        ax = axes[row, col]
        
        # 筛选当前rank的数据
        rank_data = df[df['Rank'] == rank]
        
        # 创建箱线图
        sns.boxplot(x='Condition', y='Abundance', data=rank_data, ax=ax, palette=colors)
        
        # 添加散点
        sns.stripplot(x='Condition', y='Abundance', data=rank_data, ax=ax, 
                     color='black', size=3, alpha=0.6, jitter=True)
        
        ax.set_title(f'Rank {rank}')
        ax.set_ylabel('Abundance')
        ax.set_xlabel('')
        ax.tick_params(axis='x', rotation=45)
        ax.set_yscale('log')
    
    plt.suptitle(title, fontsize=16)
    plt.tight_layout()
    
    if save_name:
        plt.savefig(save_name, dpi=300, bbox_inches='tight')
    
    plt.show()

# 5. 创建对比图
def create_comparison_plots():
    """创建1D/2D/3D空间对比图（分开的两张图）"""
    
    # 设置字体大小
    plt.rcParams.update({'font.size': 11})
    plt.rcParams.update({'axes.titlesize': 12})
    plt.rcParams.update({'axes.labelsize': 24})  # 增大坐标轴标签字体
    plt.rcParams.update({'xtick.labelsize': 16})  # 增大X轴刻度数字字体
    plt.rcParams.update({'ytick.labelsize': 16})  # 增大Y轴刻度数字字体
    plt.rcParams.update({'legend.fontsize': 14})  # 增大图例字体
    
    # 准备数据
    plot_data = []
    
    # 有扩散的情况
    for condition in ['1D_with_diff', '2D_with_diff', '3D_with_diff']:
        abund_list = results[condition]
        for sim_idx, abund in enumerate(abund_list):
            sorted_abund = np.sort(abund)[::-1]
            ranks = np.arange(1, len(sorted_abund)+1)
            
            # 只取前10个rank
            for rank, abundance in zip(ranks[:10], sorted_abund[:10]):
                dim = condition[0]
                plot_data.append({
                    'Condition': f'{dim}D with diffusion',
                    'Rank': rank,
                    'Abundance': abundance,
                    'Simulation': sim_idx
                })
    
    # 无扩散的情况（所有维度都一样）
    abund_list = results['1D_no_diff']  # 用1D的结果代表所有无扩散情况
    for sim_idx, abund in enumerate(abund_list):
        sorted_abund = np.sort(abund)[::-1]
        ranks = np.arange(1, len(sorted_abund)+1)
        
        # 只取前10个rank
        for rank, abundance in zip(ranks[:10], sorted_abund[:10]):
            plot_data.append({
                'Condition': 'No diffusion (baseline)',
                'Rank': rank,
                'Abundance': abundance,
                'Simulation': sim_idx
            })
    
    df = pd.DataFrame(plot_data)
    
    # 色盲友好的配色
    colorblind_palette = ['#0173B2', '#DE8F05', '#029E73', '#D55E00']

    # 第一张图：箱线图
    plt.figure(figsize=(14, 8))
    sns.boxplot(x='Rank', y='Abundance', hue='Condition', data=df, 
            palette=colorblind_palette)
    plt.yscale('log')
    plt.title('Species Abundance Distribution by Rank', fontsize=16, fontweight='bold')
    plt.xlabel('Species Rank', fontsize=24)
    plt.ylabel('Abundance', fontsize=24)
    plt.legend(title='Spatial Configuration', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=14)
    plt.tight_layout()
    plt.savefig('species_abundance_boxplots.png', dpi=300, bbox_inches='tight')
    plt.show()

    # 第二张图：趋势线图（带error bar）
    plt.figure(figsize=(14, 8))
    conditions = ['1D with diffusion', '2D with diffusion', '3D with diffusion', 'No diffusion (baseline)']
    for i, condition in enumerate(conditions):
        condition_data = df[df['Condition'] == condition]
        medians = condition_data.groupby('Rank')['Abundance'].median()
        stds = condition_data.groupby('Rank')['Abundance'].std()
        
        # 使用errorbar绘制带误差条的趋势线
        plt.errorbar(medians.index, medians.values, yerr=stds.values,
                     color=colorblind_palette[i], linewidth=3, marker='o', 
                     markersize=8, label=condition, capsize=5, capthick=2)

    plt.yscale('log')
    plt.xlabel('Species Rank', fontsize=24)
    plt.ylabel('Median Species Abundance', fontsize=24)
    plt.title('Median Abundance Trends', fontsize=16, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('median_abundance_trends.png', dpi=300, bbox_inches='tight')
    plt.show()

# 6. 执行可视化
print("Generating visualization plots...")

# 生成对比图（箱线图版本）
create_comparison_plots()

# 打印统计摘要
print("\n" + "="*80)
print("Statistical Summary (50 simulations)")
print("="*80)

# 检查无扩散情况下不同维度是否真的相同
print("Verifying consistency across dimensions (no diffusion):")
no_diff_1d = np.mean(results['1D_no_diff'], axis=0)
no_diff_2d = np.mean(results['2D_no_diff'], axis=0)
no_diff_3d = np.mean(results['3D_no_diff'], axis=0)

diff_1d_2d = np.mean(np.abs(no_diff_1d - no_diff_2d))
diff_1d_3d = np.mean(np.abs(no_diff_1d - no_diff_3d))
diff_2d_3d = np.mean(np.abs(no_diff_2d - no_diff_3d))

print(f"1D vs 2D no diffusion mean difference: {diff_1d_2d:.6f}")
print(f"1D vs 3D no diffusion mean difference: {diff_1d_3d:.6f}")
print(f"2D vs 3D no diffusion mean difference: {diff_2d_3d:.6f}")

if diff_1d_2d < 1e-10 and diff_1d_3d < 1e-10 and diff_2d_3d < 1e-10:
    print("✓ No diffusion results are identical across all dimensions (within numerical precision)")
else:
    print("⚠ Small differences detected in no diffusion results across dimensions")

print("\nStatistical results by condition:")

# 只打印有扩散的情况和无扩散的代表情况
conditions_to_print = ['1D_with_diff', '2D_with_diff', '3D_with_diff', '1D_no_diff']

for condition in conditions_to_print:
    abund_list = results[condition]
    mean_abund = np.mean(abund_list, axis=0)
    std_abund = np.std(abund_list, axis=0)
    
    # 计算多样性指标
    total_abund = np.sum(mean_abund)
    shannon = -np.sum((mean_abund/total_abund) * np.log(mean_abund/total_abund + 1e-10))
    
    condition_name = condition.replace('_', ' ').title()
    if condition == '1D_no_diff':
        condition_name = 'No Diffusion (all dimensions)'
    elif condition == '1D_with_diff':
        condition_name = '1D With Diffusion'
    elif condition == '2D_with_diff':
        condition_name = '2D With Diffusion'
    elif condition == '3D_with_diff':
        condition_name = '3D With Diffusion'
    
    print(f"\n{condition_name}:")
    print(f"  Mean total abundance: {total_abund:.3f}")
    print(f"  Shannon diversity: {shannon:.3f}")
    print(f"  Maximum abundance: {np.max(mean_abund):.3f}")
    print(f"  Minimum abundance: {np.min(mean_abund):.3f}")
    print(f"  Abundance standard deviation: {np.std(mean_abund):.3f}")

print("\nAll visualization plots have been generated!") 