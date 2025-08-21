import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy import stats
from scipy.optimize import minimize
import param


# ========== 无量纲化函数 ==========
def dimensionless_parameters():
    """
    返回无量纲化的参考尺度
    所有参数都将基于这些尺度进行归一化
    """
    # 参考尺度设定
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
        # 复制原始参数
        norm_cp = cp.copy()
        
        # 无量纲化资源输入率: ρ̃ = ρ / (R0/T0)
        norm_cp['rho'] = cp['rho'] / (R0 / T0)
        
        # 无量纲化资源消耗率: ω̃ = ω * T0
        norm_cp['omega'] = cp['omega'] * T0
        
        # 无量纲化物种死亡率: m̃ = m * T0
        norm_cp['m'] = cp['m'] * T0
        
        # 无量纲化消耗矩阵: ũ = u * C0 * R0 * T0
        norm_cp['u'] = cp['u'] * C0 * R0 * T0
        
        # 泄漏张量保持比例形式（如果本身就是无量纲的）
        norm_cp['l'] = cp['l']  # 假设l已经是无量纲的
        
        normalized_params.append(norm_cp)
    
    return normalized_params

def normalize_diffusion_coefficients(D0, D_species, T0, L0):
    """
    无量纲化扩散系数: D̃ = D * T0 / L0^2
    """
    D0_norm = D0 * T0 / (L0**2)
    D_species_norm = D_species * T0 / (L0**2)
    return D0_norm, D_species_norm

def normalize_initial_conditions(Y0, R0, C0):
    """
    无量纲化初始条件
    """
    Y0_norm = Y0.copy()
    
    # 物种初始条件归一化: C* = C / C0
    for k in range(K):
        Y0_norm[k*N_pool:(k+1)*N_pool] /= C0
    
    # 资源初始条件归一化: R* = R / R0
    M = 15  # 使用实际的M值
    for k in range(K):
        Y0_norm[K*N_pool + k*M:K*N_pool + (k+1)*M] /= R0
    
    return Y0_norm

def create_distance_decay_matrix(K, lambda_spatial, beta):
    # 1D线性布局
    coords = np.arange(K)
    D_kj = np.abs(coords[:, None] - coords[None, :])  # 线性距离
    
    # 连接度：每个patch都有K-1个邻居
    k_j = np.ones(K) * (K - 1)
    
    A = np.zeros((K, K))
    for k in range(K):
        for j in range(K):
            if k != j:
                A[k, j] = (k_j[j] ** beta) * np.exp(-lambda_spatial * D_kj[k, j])
    
    return A, coords

# ========== 无量纲化参数设置 ==========
R0, C0, T0, L0 = dimensionless_parameters()

# ========== 分析函数 ==========
def gini_coefficient(x):
    x = np.sort(np.array(x))
    n = len(x)
    if np.sum(x) == 0:
        return 0
    return (2 * np.sum((np.arange(1, n+1) * x)) / (n * np.sum(x))) - (n + 1) / n

def shannon_diversity(x):
    x = np.array(x)
    x = x[x > 0]
    if len(x) == 0:
        return 0
    p = x / np.sum(x)
    return -np.sum(p * np.log(p))

def simpson_index(x):
    """Simpson多样性指数 (1-D)"""
    x = np.array(x)
    x = x[x > 0]
    if len(x) == 0:
        return 0
    p = x / np.sum(x)
    return 1 - np.sum(p**2)

def berger_parker_index(x):
    """Berger-Parker优势度指数"""
    x = np.array(x)
    x = x[x > 0]
    if len(x) == 0:
        return 0
    return np.max(x) / np.sum(x)

def shannon_evenness(x):
    """Shannon均匀度指数"""
    x = np.array(x)
    x = x[x > 0]
    if len(x) == 0:
        return 0
    H = shannon_diversity(x)
    S = len(x)  # 物种数
    if S <= 1:
        return 1.0
    return H / np.log(S)

def calculate_all_diversity_metrics(abund):
    """计算所有多样性指标"""
    # 过滤掉零丰度物种
    non_zero_abund = abund[abund > 0]
    
    if len(non_zero_abund) == 0:
        return {
            'shannon': 0,
            'simpson': 0,
            'berger_parker': 0,
            'shannon_evenness': 0,
            'species_richness': 0
        }
    
    # 计算各指标
    shannon = shannon_diversity(non_zero_abund)
    simpson = simpson_index(non_zero_abund)
    berger_parker = berger_parker_index(non_zero_abund)
    evenness = shannon_evenness(non_zero_abund)
    richness = len(non_zero_abund)
    
    return {
        'shannon': shannon,
        'simpson': simpson,
        'berger_parker': berger_parker,
        'shannon_evenness': evenness,
        'species_richness': richness
    }

def get_abund(sol):
    # 统计全局物种池每个物种的总丰度（所有patch加和）
    abund = np.zeros(N_pool)
    for k, cp in enumerate(community_params):
        patch_species = cp['species_idx']
        C_patch = sol.y[C_idx[k][0]:C_idx[k][1], -1]  # 取最后一个时间点
        for sp in patch_species:
            abund[sp] += C_patch[sp]
    return abund

np.random.seed(42)
K = 50  # patch数
N_pool = 50# 全局物种池总数
M = 20 # 每个patch的资源数（可固定）

# 每个patch随机抽取物种子集
patch_species_list = []  # 每个patch的物种编号列表
N_k_list = []            # 每个patch实际物种数
for k in range(K):
    N_k = np.random.randint(3, N_pool+1)  # 每个patch物种数3~N_pool
    N_k_list.append(N_k)
    patch_species = np.sort(np.random.choice(np.arange(N_pool), N_k, replace=False))
    patch_species_list.append(patch_species)
M_k_list = [M for _ in range(K)]

# 原始参数设置（保持原有逻辑）
community_params = []
for k in range(K):
    N_k = N_k_list[k]
    patch_species = patch_species_list[k]
    u_k = param.modular_uptake(N_k, M, N_modules=1, s_ratio=10.0)
    l_k = param.generate_l_tensor(N_k, M, N_modules=2, s_ratio=10.0, λ=0.3)
    rho_k = np.ones(M) * (5.0 + 5.0 * np.sin(np.pi * k / (K-1)))
    omega_k = np.ones(M) * (0.05 + 0.01 * k)
    m_k = np.ones(N_k) * (0.09 + 0.02 * k)
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

# 无量纲化初始条件
Y0 = np.zeros(n_total)
for k, cp in enumerate(community_params):
    N_k = cp['N']
    patch_species = cp['species_idx']
    # 只给本patch有的物种赋初值，其余为0
    Y0[C_idx[k][0]:C_idx[k][1]][patch_species] = np.random.uniform(0.5, 1.5, N_k)
    Y0[R_idx[k][0]:R_idx[k][1]] = np.ones(M) * 2.0

# 无量纲化初始条件
Y0 = normalize_initial_conditions(Y0, R0, C0)

# ========== 预处理：统计每个物种在哪些patch中出现 ==========
species_in_patches = [[] for _ in range(N_pool)]
for k, cp in enumerate(community_params):
    for sp in cp['species_idx']:
        species_in_patches[sp].append(k)

# ========== ODE系统（使用距离衰减矩阵） ==========
def global_ode_distance_decay(t, y, with_diffusion=True, D0=0.1, D_species=0.01, lambda_spatial=1.0, beta=1.0):
    # 无量纲化扩散系数
    D0_norm, D_species_norm = normalize_diffusion_coefficients(D0, D_species, T0, L0)
    
    # 创建距离衰减矩阵
    A, positions = create_distance_decay_matrix(K, lambda_spatial, beta)
    
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
        # 只更新本patch有的物种
        for idx, sp in enumerate(patch_species):
            dydt[C_idx[k][0] + sp] = dCdt[idx]
        dydt[R_idx[k][0]:R_idx[k][1]] = dRdt
    
    # 物种扩散项（使用距离衰减矩阵）
    if with_diffusion and D_species_norm > 0:
        for i in range(N_pool):
            patch_list = species_in_patches[i]
            if len(patch_list) < 2:
                continue  # 只在多个patch出现才扩散
            
            # 计算这个物种在所有patch中的扩散
            for k in patch_list:
                for j in range(K):
                    if j != k and A[k, j] > 0:  # 使用距离衰减矩阵
                        # 检查物种i是否在patch j中存在
                        if i in community_params[j]['species_idx']:
                            C_k = y[C_idx[k][0] + i]
                            C_j = y[C_idx[j][0] + i]
                            dydt[C_idx[k][0] + i] += D_species_norm * A[k, j] * (C_j - C_k)
    
    # 资源扩散项（使用距离衰减矩阵）
    if with_diffusion:
        for alpha in range(M):
            for k in range(K):
                for j in range(K):
                    if j != k and A[k, j] > 0:  # 使用距离衰减矩阵
                        R_k_alpha = y[R_idx[k][0] + alpha]
                        R_j_alpha = y[R_idx[j][0] + alpha]
                        dydt[R_idx[k][0] + alpha] += D0_norm * A[k, j] * (R_j_alpha - R_k_alpha)
    
    return dydt

# ========== 数值积分 ==========
t_span = (0, 1500)
t_eval = np.linspace(*t_span, 1500)

# 距离衰减参数设置（固定为1）
lambda_spatial = 1.0  # 空间衰减参数
beta = 0.5            # 网络拓扑参数

# 扩散系数设置
D0_list = [1.0]  # 固定资源扩散系数
D_species_list = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0] # 更多物种扩散系数取值

# 进行10次模拟
n_simulations = 50

# 存储所有模拟结果
all_rank_abundances = []  # 存储所有rank-abundance曲线数据
all_gini = np.zeros((n_simulations, len(D0_list), len(D_species_list)))
all_shannon = np.zeros((n_simulations, len(D0_list), len(D_species_list)))
all_simpson = np.zeros((n_simulations, len(D0_list), len(D_species_list)))
all_berger_parker = np.zeros((n_simulations, len(D0_list), len(D_species_list)))
all_evenness = np.zeros((n_simulations, len(D0_list), len(D_species_list)))
all_richness = np.zeros((n_simulations, len(D0_list), len(D_species_list)))

# 存储分布拟合结果
all_distribution_fits = []  # 存储所有分布拟合结果

# 分布拟合函数定义
def fit_log_normal(abundances):
    """
    拟合log-normal分布
    返回拟合参数和AIC
    """
    # 过滤掉零值
    valid_abundances = abundances[abundances > 0]
    if len(valid_abundances) < 3:
        return None, None, None
    
    # 对丰度取对数
    log_abundances = np.log(valid_abundances)
    
    # 拟合正态分布（对数的丰度）
    mu, sigma = stats.norm.fit(log_abundances)
    
    # 计算AIC
    log_likelihood = np.sum(stats.norm.logpdf(log_abundances, mu, sigma))
    aic = 2 * 2 - 2 * log_likelihood  # 2个参数
    
    return mu, sigma, aic

def fit_log_series(abundances):
    """
    拟合log-series分布
    返回α参数和AIC
    """
    # 过滤掉零值
    valid_abundances = abundances[abundances > 0]
    if len(valid_abundances) < 3:
        return None, None
    
    # 计算物种总数和个体总数
    S = len(valid_abundances)  # 物种数
    N = np.sum(valid_abundances)  # 个体总数
    
    if N <= S:
        return None, None
    
    # 使用Fisher's α估计
    # 通过数值方法求解: S = α * ln(1 + N/α)
    def equation(alpha):
        return S - alpha * np.log(1 + N / alpha)
    
    # 初始猜测
    alpha_guess = N / np.log(N)
    
    try:
        # 使用scipy.optimize求解
        result = minimize(lambda x: abs(equation(x)), alpha_guess, 
                         bounds=[(0.1, N)], method='L-BFGS-B')
        alpha = result.x[0]
        
        if not result.success:
            return None, None
        
        # 计算AIC（简化版本）
        # log-series的似然函数比较复杂，这里用简化方法
        expected_abundances = N * (1/alpha) * (1 - 1/alpha)**(np.arange(1, S+1) - 1)
        expected_abundances = expected_abundances[:len(valid_abundances)]
        
        # 计算卡方统计量作为拟合优度
        chi_square = np.sum((valid_abundances - expected_abundances)**2 / expected_abundances)
        aic = chi_square + 2  # 1个参数
        
        return alpha, aic
        
    except:
        return None, None

def fit_distributions(abundances):
    """
    对物种丰度分布拟合log-normal和log-series
    返回拟合结果字典
    """
    results = {}
    
    # 拟合log-normal
    mu, sigma, aic_ln = fit_log_normal(abundances)
    results['log_normal'] = {
        'mu': mu,
        'sigma': sigma,
        'aic': aic_ln,
        'success': aic_ln is not None
    }
    
    # 拟合log-series
    alpha, aic_ls = fit_log_series(abundances)
    results['log_series'] = {
        'alpha': alpha,
        'aic': aic_ls,
        'success': aic_ls is not None
    }
    
    # 判断哪个模型更好
    if aic_ln is not None and aic_ls is not None:
        if aic_ln < aic_ls:
            results['better_model'] = 'log_normal'
            results['aic_difference'] = aic_ls - aic_ln
        else:
            results['better_model'] = 'log_series'
            results['aic_difference'] = aic_ln - aic_ls
    elif aic_ln is not None:
        results['better_model'] = 'log_normal'
        results['aic_difference'] = None
    elif aic_ls is not None:
        results['better_model'] = 'log_series'
        results['aic_difference'] = None
    else:
        results['better_model'] = 'none'
        results['aic_difference'] = None
    
    return results

print(f"开始进行{n_simulations}次模拟...")

for sim in range(n_simulations):
    print(f"模拟进度: {sim+1}/{n_simulations}")
    
    # 重新生成随机参数（保持相同的随机种子结构）
    np.random.seed(42 + sim)
    
    # 重新生成patch物种分配
    patch_species_list = []
    N_k_list = []
    for k in range(K):
        N_k = np.random.randint(3, N_pool+1)
        N_k_list.append(N_k)
        patch_species = np.sort(np.random.choice(np.arange(N_pool), N_k, replace=False))
        patch_species_list.append(patch_species)
    
    # 重新生成社区参数
    community_params = []
    for k in range(K):
        N_k = N_k_list[k]
        patch_species = patch_species_list[k]
        u_k = param.modular_uptake(N_k, M, N_modules=1, s_ratio=10.0)
        l_k = param.generate_l_tensor(N_k, M, N_modules=2, s_ratio=10.0, λ=0.3)
        rho_k = np.ones(M) * (5.0 + 5.0 * np.sin(np.pi * k / (K-1)))
        omega_k = np.ones(M) * (0.05 + 0.01 * k)
        m_k = np.ones(N_k) * (0.09 + 0.02 * k)
        community_params.append({'N': N_k, 'M': M, 'u': u_k, 'l': l_k, 'rho': rho_k, 'omega': omega_k, 'm': m_k, 'species_idx': patch_species})
    
    # 无量纲化参数
    community_params = normalize_parameters(community_params, R0, C0, T0, L0)
    
    # 重新生成初始条件
    Y0 = np.zeros(n_total)
    for k, cp in enumerate(community_params):
        N_k = cp['N']
        patch_species = cp['species_idx']
        Y0[C_idx[k][0]:C_idx[k][1]][patch_species] = np.random.uniform(0.5, 1.5, N_k)
        Y0[R_idx[k][0]:R_idx[k][1]] = np.ones(M) * 2.0
    
    Y0 = normalize_initial_conditions(Y0, R0, C0)
    
    # 重新统计物种在patch中的分布
    species_in_patches = [[] for _ in range(N_pool)]
    for k, cp in enumerate(community_params):
        for sp in cp['species_idx']:
            species_in_patches[sp].append(k)
    
    # 对每个扩散系数组合进行模拟
    for i, D0 in enumerate(D0_list):
        for j, D_species in enumerate(D_species_list):
            try:
                sol = solve_ivp(lambda t, y: global_ode_distance_decay(t, y, with_diffusion=True, D0=D0, D_species=D_species, lambda_spatial=lambda_spatial, beta=beta),
                                t_span, Y0, t_eval=t_eval, method='LSODA', rtol=1e-5, atol=1e-7)
                
                abund = get_abund(sol)
                
                # 存储完整的rank-abundance数据
                sorted_abund = np.sort(abund)[::-1]  # 从大到小排序
                rank_data = {
                    'sim': sim,
                    'D0': D0,
                    'D_species': D_species,
                    'ranks': np.arange(1, len(sorted_abund)+1),
                    'abundances': sorted_abund
                }
                all_rank_abundances.append(rank_data)
                
                # 计算全局多样性指标
                diversity_metrics = calculate_all_diversity_metrics(abund)
                all_shannon[sim, i, j] = diversity_metrics['shannon']
                all_simpson[sim, i, j] = diversity_metrics['simpson']
                all_berger_parker[sim, i, j] = diversity_metrics['berger_parker']
                all_evenness[sim, i, j] = diversity_metrics['shannon_evenness']
                all_richness[sim, i, j] = diversity_metrics['species_richness']
                
                # 计算空间分布的不平等性（基于所有物种）
                patch_totals = []
                for k, cp in enumerate(community_params):
                    patch_species = cp['species_idx']
                    C_patch = sol.y[C_idx[k][0]:C_idx[k][1], -1]
                    patch_total = np.sum([C_patch[sp] for sp in patch_species])
                    patch_totals.append(patch_total)
                
                all_gini[sim, i, j] = gini_coefficient(patch_totals)
                
                # 进行分布拟合
                fit_results = fit_distributions(abund)
                fit_data = {
                    'sim': sim,
                    'D0': D0,
                    'D_species': D_species,
                    'fit_results': fit_results
                }
                all_distribution_fits.append(fit_data)
                
            except Exception as e:
                print(f"模拟 {sim+1}, D0={D0}, D_species={D_species} 失败: {e}")
                all_gini[sim, i, j] = np.nan
                all_shannon[sim, i, j] = np.nan
                all_simpson[sim, i, j] = np.nan
                all_berger_parker[sim, i, j] = np.nan
                all_evenness[sim, i, j] = np.nan
                all_richness[sim, i, j] = np.nan

print("所有模拟完成！")

# 计算平均值和标准差
mean_gini = np.nanmean(all_gini, axis=0)
std_gini = np.nanstd(all_gini, axis=0)
mean_shannon = np.nanmean(all_shannon, axis=0)
std_shannon = np.nanstd(all_shannon, axis=0)
mean_simpson = np.nanmean(all_simpson, axis=0)
std_simpson = np.nanstd(all_simpson, axis=0)
mean_berger_parker = np.nanmean(all_berger_parker, axis=0)
std_berger_parker = np.nanstd(all_berger_parker, axis=0)
mean_evenness = np.nanmean(all_evenness, axis=0)
std_evenness = np.nanstd(all_evenness, axis=0)
mean_richness = np.nanmean(all_richness, axis=0)
std_richness = np.nanstd(all_richness, axis=0)

# 计算斜率（使用对数坐标）
def calculate_slope(x, y):
    """计算对数坐标下的斜率"""
    log_x = np.log10(x)
    log_y = np.log10(y)
    
    # 使用线性回归计算斜率
    n = len(log_x)
    if n < 2:
        return np.nan
    
    # 计算线性回归
    x_mean = np.mean(log_x)
    y_mean = np.mean(log_y)
    
    numerator = np.sum((log_x - x_mean) * (log_y - y_mean))
    denominator = np.sum((log_x - x_mean) ** 2)
    
    if denominator == 0:
        return np.nan
    
    slope = numerator / denominator
    return slope

def fit_log_normal(abundances):
    """
    拟合log-normal分布
    返回拟合参数和AIC
    """
    # 过滤掉零值
    valid_abundances = abundances[abundances > 0]
    if len(valid_abundances) < 3:
        return None, None, None
    
    # 对丰度取对数
    log_abundances = np.log(valid_abundances)
    
    # 拟合正态分布（对数的丰度）
    mu, sigma = stats.norm.fit(log_abundances)
    
    # 计算AIC
    log_likelihood = np.sum(stats.norm.logpdf(log_abundances, mu, sigma))
    aic = 2 * 2 - 2 * log_likelihood  # 2个参数
    
    return mu, sigma, aic

def fit_log_series(abundances):
    """
    拟合log-series分布
    返回α参数和AIC
    """
    # 过滤掉零值
    valid_abundances = abundances[abundances > 0]
    if len(valid_abundances) < 3:
        return None, None
    
    # 计算物种总数和个体总数
    S = len(valid_abundances)  # 物种数
    N = np.sum(valid_abundances)  # 个体总数
    
    if N <= S:
        return None, None
    
    # 使用Fisher's α估计
    # 通过数值方法求解: S = α * ln(1 + N/α)
    def equation(alpha):
        return S - alpha * np.log(1 + N / alpha)
    
    # 初始猜测
    alpha_guess = N / np.log(N)
    
    try:
        # 使用scipy.optimize求解
        result = minimize(lambda x: abs(equation(x)), alpha_guess, 
                         bounds=[(0.1, N)], method='L-BFGS-B')
        alpha = result.x[0]
        
        if not result.success:
            return None, None
        
        # 计算AIC（简化版本）
        # log-series的似然函数比较复杂，这里用简化方法
        expected_abundances = N * (1/alpha) * (1 - 1/alpha)**(np.arange(1, S+1) - 1)
        expected_abundances = expected_abundances[:len(valid_abundances)]
        
        # 计算卡方统计量作为拟合优度
        chi_square = np.sum((valid_abundances - expected_abundances)**2 / expected_abundances)
        aic = chi_square + 2  # 1个参数
        
        return alpha, aic
        
    except:
        return None, None

def fit_distributions(abundances):
    """
    对物种丰度分布拟合log-normal和log-series
    返回拟合结果字典
    """
    results = {}
    
    # 拟合log-normal
    mu, sigma, aic_ln = fit_log_normal(abundances)
    results['log_normal'] = {
        'mu': mu,
        'sigma': sigma,
        'aic': aic_ln,
        'success': aic_ln is not None
    }
    
    # 拟合log-series
    alpha, aic_ls = fit_log_series(abundances)
    results['log_series'] = {
        'alpha': alpha,
        'aic': aic_ls,
        'success': aic_ls is not None
    }
    
    # 判断哪个模型更好
    if aic_ln is not None and aic_ls is not None:
        if aic_ln < aic_ls:
            results['better_model'] = 'log_normal'
            results['aic_difference'] = aic_ls - aic_ln
        else:
            results['better_model'] = 'log_series'
            results['aic_difference'] = aic_ln - aic_ls
    elif aic_ln is not None:
        results['better_model'] = 'log_normal'
        results['aic_difference'] = None
    elif aic_ls is not None:
        results['better_model'] = 'log_series'
        results['aic_difference'] = None
    else:
        results['better_model'] = 'none'
        results['aic_difference'] = None
    
    return results

# 计算每个扩散系数组合的斜率
slopes_gini = np.zeros((len(D0_list), len(D_species_list)))
slopes_shannon = np.zeros((len(D0_list), len(D_species_list)))
slopes_simpson = np.zeros((len(D0_list), len(D_species_list)))

for i, D0 in enumerate(D0_list):
    for j, D_species in enumerate(D_species_list):
        # 计算Gini指数的斜率
        x_data = D_species_list
        y_data_gini = mean_gini[i, :]
        slopes_gini[i, j] = calculate_slope(x_data, y_data_gini)
        
        # 计算Shannon指数的斜率
        y_data_shannon = mean_shannon[i, :]
        slopes_shannon[i, j] = calculate_slope(x_data, y_data_shannon)
        
        # 计算Simpson指数的斜率
        y_data_simpson = mean_simpson[i, :]
        slopes_simpson[i, j] = calculate_slope(x_data, y_data_simpson)

# 计算每个D0对应的平均斜率
mean_slopes_gini = np.nanmean(slopes_gini, axis=1)
mean_slopes_shannon = np.nanmean(slopes_shannon, axis=1)
mean_slopes_simpson = np.nanmean(slopes_simpson, axis=1)

# 计算分段斜率函数
def fit_segment_slopes(ranks, abundances, n_segments=3):
    """分段拟合斜率"""
    n_species = len(ranks)
    segment_size = n_species // n_segments
    slopes = []
    
    for seg in range(n_segments):
        start_idx = seg * segment_size
        end_idx = (seg + 1) * segment_size if seg < n_segments - 1 else n_species
        
        seg_ranks = ranks[start_idx:end_idx]
        seg_abund = abundances[start_idx:end_idx]
        
        if len(seg_ranks) > 1:
            slope = calculate_slope(seg_ranks, seg_abund)
            slopes.append(slope)
        else:
            slopes.append(np.nan)
    
    return slopes

# 计算分段斜率
print("计算分段斜率...")
segment_slopes = {}
for i, D0 in enumerate(D0_list):
    for j, D_species in enumerate(D_species_list):
        param_data = [data for data in all_rank_abundances 
                     if data['D0'] == D0 and data['D_species'] == D_species]
        
        if param_data:
            all_abundances = np.array([data['abundances'] for data in param_data])
            mean_abundances = np.nanmean(all_abundances, axis=0)
            ranks = param_data[0]['ranks']
            
            slopes = fit_segment_slopes(ranks, mean_abundances, n_segments=3)
            segment_slopes[f'D0={D0}, D_s={D_species}'] = slopes

# 绘制完整的Rank-Abundance分布（对数坐标）
plt.figure(figsize=(12, 8))

# 设置字体大小（至少11号）
plt.rcParams.update({'font.size': 11})
plt.rcParams.update({'axes.titlesize': 20})
plt.rcParams.update({'axes.labelsize': 24})  # 从14改为18，让坐标轴标签更大
plt.rcParams.update({'xtick.labelsize': 20})
plt.rcParams.update({'ytick.labelsize': 20})
plt.rcParams.update({'legend.fontsize': 10})

# 首先添加无扩散基线（黑色实线）
print("计算无扩散基线...")
no_diff_abundances = []
for sim in range(n_simulations):
    np.random.seed(42 + sim)
    # 重新生成参数（与主循环相同）
    patch_species_list = []
    N_k_list = []
    for k in range(K):
        N_k = np.random.randint(3, N_pool+1)
        N_k_list.append(N_k)
        patch_species = np.sort(np.random.choice(np.arange(N_pool), N_k, replace=False))
        patch_species_list.append(patch_species)
    
    community_params = []
    for k in range(K):
        N_k = N_k_list[k]
        patch_species = patch_species_list[k]
        u_k = param.modular_uptake(N_k, M, N_modules=1, s_ratio=10.0)
        l_k = param.generate_l_tensor(N_k, M, N_modules=2, s_ratio=10.0, λ=0.3)
        rho_k = np.ones(M) * (5.0 + 5.0 * np.sin(np.pi * k / (K-1)))
        omega_k = np.ones(M) * (0.05 + 0.01 * k)
        m_k = np.ones(N_k) * (0.09 + 0.02 * k)
        community_params.append({'N': N_k, 'M': M, 'u': u_k, 'l': l_k, 'rho': rho_k, 'omega': omega_k, 'm': m_k, 'species_idx': patch_species})
    
    community_params = normalize_parameters(community_params, R0, C0, T0, L0)
    
    Y0 = np.zeros(n_total)
    for k, cp in enumerate(community_params):
        N_k = cp['N']
        patch_species = cp['species_idx']
        Y0[C_idx[k][0]:C_idx[k][1]][patch_species] = np.random.uniform(0.5, 1.5, N_k)
        Y0[R_idx[k][0]:R_idx[k][1]] = np.ones(M) * 2.0
    
    Y0 = normalize_initial_conditions(Y0, R0, C0)
    
    # 无扩散模拟
    sol_no_diff = solve_ivp(lambda t, y: global_ode_distance_decay(t, y, with_diffusion=False),
                            t_span, Y0, t_eval=t_eval, method='LSODA', rtol=1e-5, atol=1e-7)
    abund_no_diff = get_abund(sol_no_diff)
    sorted_abund_no_diff = np.sort(abund_no_diff)[::-1]
    no_diff_abundances.append(sorted_abund_no_diff)

# 计算无扩散的平均rank-abundance曲线和误差条
no_diff_mean = np.nanmean(no_diff_abundances, axis=0)
no_diff_std = np.nanstd(no_diff_abundances, axis=0)
ranks_no_diff = np.arange(1, len(no_diff_mean)+1)

# 绘制无扩散基线（黑色实线）带误差条
plt.plot(ranks_no_diff, no_diff_mean, 'k-', linewidth=3, label='No Diffusion (baseline)', alpha=0.9)
plt.fill_between(ranks_no_diff, no_diff_mean - no_diff_std, no_diff_mean + no_diff_std, 
                 color='black', alpha=0.2)

# 注释掉完全扩散的计算和绘制
# print("计算完全扩散状态...")
# ... (完全扩散相关代码已注释)

# 绘制其他扩散参数组合
print(f"all_rank_abundances中有{len(all_rank_abundances)}条数据")
print("参数组合:")
for data in all_rank_abundances:
    print(f"  D0={data['D0']}, D_s={data['D_species']}")

# 使用不同的颜色和线型来区分所有10条线
colors = ['blue', 'green', 'red', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
linestyles = ['--', '-.', ':', '--', '-.', ':', '--', '-.', ':', '--']

for i, D0 in enumerate(D0_list):
    for j, D_species in enumerate(D_species_list):
        # 收集该参数组合的所有rank-abundance数据
        param_data = [data for data in all_rank_abundances 
                     if data['D0'] == D0 and data['D_species'] == D_species]
        
        print(f"D0={D0}, D_s={D_species}: 找到{len(param_data)}条数据")
        
        if param_data:
            # 计算平均rank-abundance曲线和误差条
            all_abundances = np.array([data['abundances'] for data in param_data])
            mean_abundances = np.nanmean(all_abundances, axis=0)
            std_abundances = np.nanstd(all_abundances, axis=0)
            ranks = param_data[0]['ranks']
            
            # 使用不同的颜色和线型
            line_idx = j  # 因为只有一个D0值，所以直接用j作为索引
            
            # 检查数据是否适合对数坐标
            if np.any(mean_abundances <= 0):
                print(f"警告: D_species={D_species} 有非正值，不适合对数坐标")
                # 将0值替换为很小的正数
                mean_abundances = np.where(mean_abundances <= 0, 1e-10, mean_abundances)
            
            plt.plot(ranks, mean_abundances, 
                    color=colors[line_idx], alpha=0.8, linewidth=2, 
                    linestyle=linestyles[line_idx],
                    label=f'D_species={D_species}')
            
            # 添加误差条
            plt.fill_between(ranks, mean_abundances - std_abundances, mean_abundances + std_abundances, 
                           color=colors[line_idx], alpha=0.2)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Species Rank', fontsize=18)  # 从14改为18
plt.ylabel('Species Abundance', fontsize=18)  # 从14改为18
plt.title('Complete Rank-Abundance Distributions', fontsize=12, fontweight='bold')
plt.legend(loc='lower left', fontsize=10)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('rank_abundance_distributions.png', dpi=300, bbox_inches='tight')
plt.show()

# 绘制多样性指标随D_species变化的图 - 分别绘制6个独立图
# 设置全局字体大小
plt.rcParams.update({'font.size': 11})
plt.rcParams.update({'axes.titlesize': 12})
plt.rcParams.update({'axes.labelsize': 24})  # 从20改为24，让坐标轴标签更大
plt.rcParams.update({'xtick.labelsize': 16})  # 从11改为16，让X轴刻度数字更大
plt.rcParams.update({'ytick.labelsize': 16})  # 从11改为16，让Y轴刻度数字更大
plt.rcParams.update({'legend.fontsize': 10})

# 提取数据
D_species_values = []
shannon_values = []
simpson_values = []
berger_parker_values = []
gini_values = []
evenness_values = []

# 计算误差条数据
shannon_stds = []
simpson_stds = []
berger_parker_stds = []
gini_stds = []
evenness_stds = []

for i, D0 in enumerate(D0_list):
    for j, D_species in enumerate(D_species_list):
        D_species_values.append(D_species)
        shannon_values.append(mean_shannon[i, j])
        simpson_values.append(mean_simpson[i, j])
        berger_parker_values.append(mean_berger_parker[i, j])
        gini_values.append(mean_gini[i, j])
        evenness_values.append(mean_evenness[i, j])
        
        # 计算标准差作为误差条
        shannon_stds.append(std_shannon[i, j])
        simpson_stds.append(std_simpson[i, j])
        berger_parker_stds.append(std_berger_parker[i, j])
        gini_stds.append(std_gini[i, j])
        evenness_stds.append(std_evenness[i, j])


# 1. Shannon指数
plt.figure(figsize=(10, 6))
plt.errorbar(D_species_values, shannon_values, yerr=shannon_stds, 
             fmt='bo-', linewidth=2, markersize=6, capsize=5, capthick=2)
plt.xlabel('Species Diffusion Coefficient (D_species)', fontsize=24)
plt.ylabel('Shannon Diversity Index', fontsize=24)
plt.title('Shannon Diversity vs Species Diffusion', fontsize=12, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.xscale('log')
plt.tight_layout()
plt.savefig('shannon_diversity.png', dpi=300, bbox_inches='tight')
plt.show()

# 2. Simpson指数
plt.figure(figsize=(10, 6))
plt.errorbar(D_species_values, simpson_values, yerr=simpson_stds, 
             fmt='ro-', linewidth=2, markersize=6, capsize=5, capthick=2)
plt.xlabel('Species Diffusion Coefficient (D_species)', fontsize=24)
plt.ylabel('Simpson Diversity Index', fontsize=24)
plt.title('Simpson Diversity vs Species Diffusion', fontsize=12, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.xscale('log')
plt.tight_layout()
plt.savefig('simpson_diversity.png', dpi=300, bbox_inches='tight')
plt.show()

# 3. Berger-Parker指数
plt.figure(figsize=(10, 6))
plt.errorbar(D_species_values, berger_parker_values, yerr=berger_parker_stds, 
             fmt='go-', linewidth=2, markersize=6, capsize=5, capthick=2)
plt.xlabel('Species Diffusion Coefficient (D_species)', fontsize=24)
plt.ylabel('Berger-Parker Dominance Index', fontsize=24)
plt.title('Dominance vs Species Diffusion', fontsize=12, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.xscale('log')
plt.tight_layout()
plt.savefig('berger_parker_dominance.png', dpi=300, bbox_inches='tight')
plt.show()

# 4. Gini指数
plt.figure(figsize=(10, 6))
plt.errorbar(D_species_values, gini_values, yerr=gini_stds, 
             fmt='mo-', linewidth=2, markersize=6, capsize=5, capthick=2)
plt.xlabel('Species Diffusion Coefficient (D_species)', fontsize=24)
plt.ylabel('Gini Coefficient', fontsize=24)
plt.title('Spatial Inequality vs Species Diffusion', fontsize=12, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.xscale('log')
plt.tight_layout()
plt.savefig('gini_coefficient.png', dpi=300, bbox_inches='tight')
plt.show()

# 5. Evenness
plt.figure(figsize=(10, 6))
plt.errorbar(D_species_values, evenness_values, yerr=evenness_stds, 
             fmt='co-', linewidth=2, markersize=6, capsize=5, capthick=2)
plt.xlabel('Species Diffusion Coefficient (D_species)', fontsize=24)
plt.ylabel('Shannon Evenness', fontsize=24)
plt.title('Evenness vs Species Diffusion', fontsize=12, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.xscale('log')
plt.tight_layout()
plt.savefig('shannon_evenness.png', dpi=300, bbox_inches='tight')
plt.show()

# 最后一张图：三个多样性指数的组合图
plt.figure(figsize=(12, 8))
plt.errorbar(D_species_values, shannon_values, yerr=shannon_stds, 
             marker='o', linewidth=2, markersize=8, capsize=5, capthick=2, 
             color='blue', label='Shannon')
plt.errorbar(D_species_values, simpson_values, yerr=simpson_stds, 
             marker='s', linewidth=2, markersize=8, capsize=5, capthick=2, 
             color='red', label='Simpson')
plt.errorbar(D_species_values, berger_parker_values, yerr=berger_parker_stds, 
             marker='^', linewidth=2, markersize=8, capsize=5, capthick=2, 
             color='green', label='Berger-Parker')

plt.xlabel('Species Diffusion Coefficient (D_species)', fontsize=24)
plt.ylabel('Index Value', fontsize=24)
plt.title('All Diversity Indices vs Species Diffusion', fontsize=20, fontweight='bold')
plt.xscale('log')  # X轴对数刻度
plt.legend(fontsize=14)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('all_diversity_indices.png', dpi=300, bbox_inches='tight')
plt.show()

# 打印所有多样性指标的详细结果
print("\n" + "="*80)
print("多样性指标矩阵结果 (50次模拟的平均值 ± 标准差)")
print("="*80)

# 打印表头
print(f"{'D0':<6} {'D_s':<6} {'Gini':<10} {'Shannon':<10} {'Simpson':<10} {'B-Parker':<10} {'Evenness':<10} {'Richness':<10}")
print("-"*90)

# 打印每个参数组合的结果
for i, D0 in enumerate(D0_list):
    for j, D_species in enumerate(D_species_list):
        print(f"{D0:<6.1f} {D_species:<6.2f} "
              f"{mean_gini[i,j]:<10.3f}±{std_gini[i,j]:<6.3f} "
              f"{mean_shannon[i,j]:<10.3f}±{std_shannon[i,j]:<6.3f} "
              f"{mean_simpson[i,j]:<10.3f}±{std_simpson[i,j]:<6.3f} "
              f"{mean_berger_parker[i,j]:<10.3f}±{std_berger_parker[i,j]:<6.3f} "
              f"{mean_evenness[i,j]:<10.3f}±{std_evenness[i,j]:<6.3f} "
              f"{mean_richness[i,j]:<10.1f}±{std_richness[i,j]:<6.1f}")

# 打印分段斜率分析结果
print("\n" + "="*80)
print("分段斜率分析结果 (Dominant/Intermediate/Rare segments)")
print("="*80)
for key, slopes in segment_slopes.items():
    print(f"{key}: Dominant={slopes[0]:.4f}, Intermediate={slopes[1]:.4f}, Rare={slopes[2]:.4f}")

print("\n" + "="*80)
print("斜率分析结果 (对数坐标下的变化率)")
print("="*80)

# 打印斜率统计信息
print("\n=== 斜率统计信息 ===")
print("Gini指数斜率:")
for i, D0 in enumerate(D0_list):
    print(f"D0={D0}: {mean_slopes_gini[i]:.4f}")

print("\nShannon指数斜率:")
for i, D0 in enumerate(D0_list):
    print(f"D0={D0}: {mean_slopes_shannon[i]:.4f}")

print("\nSimpson指数斜率:")
for i, D0 in enumerate(D0_list):
    print(f"D0={D0}: {mean_slopes_simpson[i]:.4f}")

# 分析分布拟合结果
print("\n" + "="*80)
print("Log-Normal vs Log-Series 分布拟合结果")
print("="*80)

# 按参数组合统计拟合结果
fit_summary = {}
for fit_data in all_distribution_fits:
    D0 = fit_data['D0']
    D_species = fit_data['D_species']
    key = f"D0={D0}, D_s={D_species}"
    
    if key not in fit_summary:
        fit_summary[key] = {
            'log_normal_success': 0,
            'log_series_success': 0,
            'log_normal_aic': [],
            'log_series_aic': [],
            'log_normal_mu': [],
            'log_normal_sigma': [],
            'log_series_alpha': [],
            'better_model_counts': {'log_normal': 0, 'log_series': 0, 'none': 0}
        }
    
    fit_results = fit_data['fit_results']
    
    # 统计log-normal拟合
    if fit_results['log_normal']['success']:
        fit_summary[key]['log_normal_success'] += 1
        fit_summary[key]['log_normal_aic'].append(fit_results['log_normal']['aic'])
        fit_summary[key]['log_normal_mu'].append(fit_results['log_normal']['mu'])
        fit_summary[key]['log_normal_sigma'].append(fit_results['log_normal']['sigma'])
    
    # 统计log-series拟合
    if fit_results['log_series']['success']:
        fit_summary[key]['log_series_success'] += 1
        fit_summary[key]['log_series_aic'].append(fit_results['log_series']['aic'])
        fit_summary[key]['log_series_alpha'].append(fit_results['log_series']['alpha'])
    
    # 统计更好的模型
    fit_summary[key]['better_model_counts'][fit_results['better_model']] += 1

# 打印拟合结果摘要
print("\n=== 分布拟合成功率 ===")
for key, summary in fit_summary.items():
    total_sims = n_simulations
    ln_success_rate = summary['log_normal_success'] / total_sims * 100
    ls_success_rate = summary['log_series_success'] / total_sims * 100
    
    print(f"\n{key}:")
    print(f"  Log-Normal 拟合成功率: {ln_success_rate:.1f}% ({summary['log_normal_success']}/{total_sims})")
    print(f"  Log-Series 拟合成功率: {ls_success_rate:.1f}% ({summary['log_series_success']}/{total_sims})")

print("\n=== 更好的模型统计 ===")
for key, summary in fit_summary.items():
    total = sum(summary['better_model_counts'].values())
    ln_better = summary['better_model_counts']['log_normal'] / total * 100
    ls_better = summary['better_model_counts']['log_series'] / total * 100
    none_better = summary['better_model_counts']['none'] / total * 100
    
    print(f"\n{key}:")
    print(f"  Log-Normal 更好: {ln_better:.1f}%")
    print(f"  Log-Series 更好: {ls_better:.1f}%")
    print(f"  无有效拟合: {none_better:.1f}%")

print("\n=== 拟合参数统计 ===")
for key, summary in fit_summary.items():
    if summary['log_normal_aic']:
        mean_mu = np.mean(summary['log_normal_mu'])
        mean_sigma = np.mean(summary['log_normal_sigma'])
        mean_aic_ln = np.mean(summary['log_normal_aic'])
        print(f"\n{key} - Log-Normal:")
        print(f"  μ = {mean_mu:.3f}, σ = {mean_sigma:.3f}, AIC = {mean_aic_ln:.3f}")
    
    if summary['log_series_aic']:
        mean_alpha = np.mean(summary['log_series_alpha'])
        mean_aic_ls = np.mean(summary['log_series_aic'])
        print(f"{key} - Log-Series:")
        print(f"  α = {mean_alpha:.3f}, AIC = {mean_aic_ls:.3f}")

print("\n=== 群落稳定性分析 ===")
print("根据拟合结果判断群落状态:")
print("- Log-Normal 拟合更好: 群落相对稳定，物种丰度分布相对均匀")
print("- Log-Series 拟合更好: 群落可能受到扰动，稀有物种较多")
print("- 拟合失败率高: 群落结构复杂，可能需要其他模型")

# 分析扩散对群落稳定性的影响
print("\n=== 扩散对群落稳定性的影响 ===")
for key, summary in fit_summary.items():
    if 'better_model_counts' in summary:
        ln_better = summary['better_model_counts']['log_normal']
        ls_better = summary['better_model_counts']['log_series']
        total = ln_better + ls_better
        
        if total > 0:
            stability_ratio = ln_better / total
            if stability_ratio > 0.6:
                status = "稳定"
            elif stability_ratio < 0.4:
                status = "受扰动"
            else:
                status = "中等"
            
            print(f"{key}: 稳定性比例 = {stability_ratio:.2f} ({status})")

# 绘制分布拟合图（类似您提供的图片）
print("\n" + "="*80)
print("生成分布拟合可视化图...")
print("="*80)

# 为每个参数组合生成一个图
for key, summary in fit_summary.items():
    if summary['log_normal_aic'] and summary['log_series_aic']:
        # 获取该参数组合的所有拟合数据
        param_fits = [fit for fit in all_distribution_fits 
                     if fit['D0'] == float(key.split(',')[0].split('=')[1]) and 
                        fit['D_species'] == float(key.split(',')[1].split('=')[1])]
        
        if len(param_fits) > 0:
            # 选择一个成功的拟合作为示例
            example_fit = None
            for fit in param_fits:
                if (fit['fit_results']['log_normal']['success'] and 
                    fit['fit_results']['log_series']['success']):
                    example_fit = fit
                    break
            
            if example_fit:
                # 获取对应的丰度数据
                abund_data = None
                for rank_data in all_rank_abundances:
                    if (rank_data['D0'] == example_fit['D0'] and 
                        rank_data['D_species'] == example_fit['D_species'] and
                        rank_data['sim'] == example_fit['sim']):
                        abund_data = rank_data
                        break
                
                if abund_data:
                    # 创建丰度八度图
                    abundances = abund_data['abundances']
                    abundances = abundances[abundances > 0]  # 过滤零值
                    
                    if len(abundances) > 0:
                        # 计算丰度八度
                        log_abundances = np.log2(abundances)
                        octaves = np.floor(log_abundances).astype(int)
                        unique_octaves, counts = np.unique(octaves, return_counts=True)
                        
                        # 创建图
                        plt.figure(figsize=(10, 6))
                        
                        # 设置字体大小
                        plt.rcParams.update({'font.size': 11})
                        plt.rcParams.update({'axes.titlesize': 12})
                        plt.rcParams.update({'axes.labelsize': 14})  # 增大坐标轴标签字体
                        plt.rcParams.update({'xtick.labelsize': 11})
                        plt.rcParams.update({'ytick.labelsize': 11})
                        plt.rcParams.update({'legend.fontsize': 10})
                        
                        # 绘制观察数据（绿色柱状图）
                        plt.bar(unique_octaves, counts, color='green', alpha=0.7, 
                               label='Observed', width=0.8)
                        
                        # 绘制拟合曲线
                        fit_results = example_fit['fit_results']
                        
                        # Log-Normal拟合（红色线）
                        if fit_results['log_normal']['success']:
                            mu = fit_results['log_normal']['mu']
                            sigma = fit_results['log_normal']['sigma']
                            aic_ln = fit_results['log_normal']['aic']
                            
                            # 生成拟合曲线
                            x_fit = np.linspace(min(unique_octaves), max(unique_octaves), 100)
                            log_x_fit = x_fit * np.log(2)  # 转换回对数尺度
                            y_fit_ln = len(abundances) * stats.norm.pdf(log_x_fit, mu, sigma) * np.log(2)
                            
                            plt.plot(x_fit, y_fit_ln, 'r-', linewidth=2, 
                                   label=f'Lognormal (AIC={aic_ln:.1f})')
                        
                        # Log-Series拟合（蓝色线）
                        if fit_results['log_series']['success']:
                            alpha = fit_results['log_series']['alpha']
                            aic_ls = fit_results['log_series']['aic']
                            
                            # 生成拟合曲线（简化版本）
                            x_fit = np.linspace(min(unique_octaves), max(unique_octaves), 100)
                            y_fit_ls = len(abundances) * (1/alpha) * (1 - 1/alpha)**(x_fit - min(unique_octaves))
                            
                            plt.plot(x_fit, y_fit_ls, 'b-', linewidth=2, 
                                   label=f'Logseries (AIC={aic_ls:.1f})')
                        
                        plt.xlabel('Abundance Octave', fontsize=14)
                        plt.ylabel('Number of Species', fontsize=14)
                        plt.title(f'Species Abundance Distribution - {key}', fontsize=12, fontweight='bold')
                        plt.legend(fontsize=10)
                        plt.grid(True, alpha=0.3)
                        plt.tight_layout()
                        plt.savefig(f'distribution_fit_{key.replace("=", "_").replace(", ", "_")}.png', dpi=300, bbox_inches='tight')
                        plt.show()
                        
                        print(f"已生成 {key} 的分布拟合图")
                        break  # 只生成一个示例图
