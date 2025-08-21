import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats
from scipy.stats import shapiro, normaltest, anderson
import warnings
warnings.filterwarnings('ignore')

# 运行原始代码来获取数据
exec(open('compare_space_structure_shannon.py').read())

# 分析分布特征
def analyze_distribution(data, name, metric):
    """分析数据分布特征"""
    print(f"\n=== {name} - {metric} 分布分析 ===")
    
    # 基本统计量
    mean_val = np.mean(data)
    std_val = np.std(data)
    skewness = stats.skew(data)
    kurtosis = stats.kurtosis(data)
    
    print(f"样本数: {len(data)}")
    print(f"均值: {mean_val:.4f}")
    print(f"标准差: {std_val:.4f}")
    print(f"偏度: {skewness:.4f}")
    print(f"峰度: {kurtosis:.4f}")
    
    # 正态性检验
    print("\n正态性检验:")
    
    # Shapiro-Wilk检验
    shapiro_stat, shapiro_p = shapiro(data)
    print(f"Shapiro-Wilk检验: 统计量={shapiro_stat:.4f}, p值={shapiro_p:.4f}")
    
    # D'Agostino K^2检验
    dagostino_stat, dagostino_p = normaltest(data)
    print(f"D'Agostino K²检验: 统计量={dagostino_stat:.4f}, p值={dagostino_p:.4f}")
    
    # Anderson-Darling检验
    anderson_result = anderson(data, dist='norm')
    print(f"Anderson-Darling检验: 统计量={anderson_result.statistic:.4f}")
    print(f"临界值: {anderson_result.critical_values}")
    print(f"显著性水平: {anderson_result.significance_level}")
    
    # 判断是否为正态分布
    is_normal = shapiro_p > 0.05 and dagostino_p > 0.05
    print(f"\n是否为正态分布: {'是' if is_normal else '否'}")
    
    return {
        'mean': mean_val,
        'std': std_val,
        'skewness': skewness,
        'kurtosis': kurtosis,
        'shapiro_p': shapiro_p,
        'dagostino_p': dagostino_p,
        'anderson_stat': anderson_result.statistic,
        'is_normal': is_normal
    }

# 创建分布可视化
def plot_distributions():
    """绘制分布图"""
    structure_names = list(space_structures.keys())
    metrics = ['alpha_shannon', 'gamma_shannon', 'beta_shannon']
    
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))
    
    for i, metric in enumerate(metrics):
        for j, name in enumerate(structure_names):
            ax = axes[i, j]
            data = results[name][metric]
            
            # 直方图
            ax.hist(data, bins=15, density=True, alpha=0.7, color='skyblue', edgecolor='black')
            
            # 核密度估计
            sns.kdeplot(data=data, ax=ax, color='red', linewidth=2)
            
            # 理论正态分布
            x = np.linspace(min(data), max(data), 100)
            normal_dist = stats.norm.pdf(x, np.mean(data), np.std(data))
            ax.plot(x, normal_dist, 'g--', linewidth=2, label='理论正态分布')
            
            ax.set_title(f'{name} - {metric}')
            ax.set_xlabel('值')
            ax.set_ylabel('密度')
            ax.legend()
            ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('distribution_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

# 执行分析
print("开始分析Shannon多样性指数的分布特征...")

# 分析每个空间结构和指标的分布
analysis_results = {}
for name in space_structures.keys():
    analysis_results[name] = {}
    for metric in ['alpha_shannon', 'gamma_shannon', 'beta_shannon']:
        data = results[name][metric]
        analysis_results[name][metric] = analyze_distribution(data, name, metric)

# 绘制分布图
plot_distributions()

# 总结
print("\n" + "="*60)
print("分布分析总结:")
print("="*60)

for name in space_structures.keys():
    print(f"\n{name}:")
    for metric in ['alpha_shannon', 'gamma_shannon', 'beta_shannon']:
        result = analysis_results[name][metric]
        print(f"  {metric}: {'正态分布' if result['is_normal'] else '非正态分布'} "
              f"(偏度: {result['skewness']:.3f}, 峰度: {result['kurtosis']:.3f})")

# 创建QQ图
fig, axes = plt.subplots(3, 3, figsize=(18, 15))
structure_names = list(space_structures.keys())
metrics = ['alpha_shannon', 'gamma_shannon', 'beta_shannon']

for i, metric in enumerate(metrics):
    for j, name in enumerate(structure_names):
        ax = axes[i, j]
        data = results[name][metric]
        
        # QQ图
        stats.probplot(data, dist="norm", plot=ax)
        ax.set_title(f'{name} - {metric} QQ图')
        ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('qq_plots.png', dpi=300, bbox_inches='tight')
plt.show() 