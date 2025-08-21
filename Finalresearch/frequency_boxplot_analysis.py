import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats

# 运行原始代码来获取数据
exec(open('compare_space_structure_shannon.py').read())

# 设置字体
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 创建组合图
def create_combined_analysis():
    """创建频率分布和箱线图的组合分析"""
    
    structure_names = list(space_structures.keys())
    metrics = ['alpha_shannon', 'gamma_shannon', 'beta_shannon']
    # 使用更容易区分的颜色：红色、绿色、蓝色
    colors = ['#E74C3C', '#27AE60', '#3498DB']  # 红色、绿色、蓝色
    
    # 创建3x2的子图布局
    fig, axes = plt.subplots(3, 2, figsize=(16, 18))
    
    for i, metric in enumerate(metrics):
        # 左侧：频率分布直方图 + 核密度估计线
        ax_left = axes[i, 0]
        
        for j, name in enumerate(structure_names):
            data = results[name][metric]
            
            # 绘制直方图
            ax_left.hist(data, bins=15, density=True, alpha=0.3, 
                        color=colors[j], edgecolor='black', linewidth=0.5, label=f'{name} (Histogram)')
            
            # 绘制核密度估计线（拟合线）
            sns.kdeplot(data=data, ax=ax_left, color=colors[j], linewidth=2, label=f'{name} (Density)')
        
        ax_left.set_xlabel('Shannon Diversity Index Value')
        ax_left.set_ylabel('Density')
        ax_left.set_title(f'{metric} Frequency Distribution')
        ax_left.legend()
        ax_left.grid(True, alpha=0.3)
        
        # 右侧：箱线图 + 散点
        ax_right = axes[i, 1]
        
        # 准备箱线图数据
        box_data = []
        box_labels = []
        for name in structure_names:
            box_data.append(results[name][metric])
            box_labels.append(name)
        
        # 绘制箱线图
        bp = ax_right.boxplot(box_data, labels=box_labels, patch_artist=True)
        
        # 设置箱线图颜色
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        # 添加散点
        for j, name in enumerate(structure_names):
            data = results[name][metric]
            x_pos = j + 1
            ax_right.scatter([x_pos] * len(data), data, 
                           color=colors[j], alpha=0.6, s=30, edgecolors='black', linewidth=0.5)
        
        ax_right.set_ylabel('Shannon Diversity Index Value')
        ax_right.set_title(f'{metric} Box Plot')
        ax_right.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('combined_frequency_boxplot.png', dpi=300, bbox_inches='tight')
    plt.show()

# 创建统计摘要
def print_statistical_summary():
    """打印统计摘要"""
    print("="*80)
    print("Shannon Diversity Index Statistical Summary")
    print("="*80)
    
    structure_names = list(space_structures.keys())
    metrics = ['alpha_shannon', 'gamma_shannon', 'beta_shannon']
    
    for metric in metrics:
        print(f"\n{metric}:")
        print("-" * 50)
        for name in structure_names:
            data = np.array(results[name][metric])
            print(f"{name}:")
            print(f"  Sample size: {len(data)}")
            print(f"  Mean: {data.mean():.4f}")
            print(f"  Std: {data.std():.4f}")
            print(f"  Median: {np.median(data):.4f}")
            print(f"  Min: {data.min():.4f}")
            print(f"  Max: {data.max():.4f}")
            print(f"  Skewness: {stats.skew(data):.4f}")
            print(f"  Kurtosis: {stats.kurtosis(data):.4f}")
            print()

# 执行分析
if __name__ == "__main__":
    print("Creating combined frequency distribution and box plot analysis...")
    
    # 创建组合图
    create_combined_analysis()
    
    # 打印统计摘要
    print_statistical_summary() 