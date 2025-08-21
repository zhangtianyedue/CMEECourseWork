import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import param
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# ========== Reference Scales for Nondimensionalization ==========
R0 = 2.0      # Resource reference concentration
C0 = 1.0      # Species reference concentration
T0 = 1.0      # Time reference scale
L0 = 1.0      # Space reference scale (patch distance)

# ========== Parameter Settings (dimensionless) ==========
D0_raw = 1       # Resource diffusion coefficient (raw)
D_species_raw = 0.5  # Species diffusion coefficient (raw)
K = 10      # Number of patches
N_pool = 50     # Species pool size
M = 5         # Number of resources per patch
n_repeat = 100# Number of simulation repeats

# Nondimensionalize diffusion coefficients
D0 = D0_raw * T0 / (L0**2)
D_species = D_species_raw * T0 / (L0**2)

# ========== Spatial Structure Generation Functions ==========
def generate_distance_decay_A_kj(K, beta=1, lambda_decay=1, L0=1.0):
    coords = np.arange(K)
    D_ij = np.abs(coords[:, None] - coords[None, :]) / L0
    k_j = np.ones(K) * (K - 1)
    A = np.zeros((K, K))
    for i in range(K):
        for j in range(K):
            if i != j:
                A[i, j] = k_j[j]**beta * np.exp(-lambda_decay * D_ij[i, j])
    row_sums = A.sum(axis=1)
    row_sums[row_sums == 0] = 1
    A = A / row_sums[:, np.newaxis]
    return A

def generate_wellmixed_A_kj(K):
    A = np.ones((K, K))
    np.fill_diagonal(A, 0)
    A = A / (K - 1)
    return A

def generate_no_structure_A_kj(K):
    return np.eye(K)

space_structures = {
    'No Structure': generate_no_structure_A_kj,
    'Distance Decay': lambda K: generate_distance_decay_A_kj(K, beta=1, lambda_decay=1),
    'Well-mixed': generate_wellmixed_A_kj
}

# ========== Diversity Index Functions ==========
def shannon_index(abundances):
    abundances = np.array(abundances)
    total = np.sum(abundances)
    if total == 0:
        return 0
    p = abundances[abundances > 0] / total
    return -np.sum(p * np.log(p))

def evenness(H, S):
    return H / np.log(S) if S > 1 else 0

# ========== Main Loop ==========
results = {name: {'alpha_shannon': [], 'gamma_shannon': [], 'beta_shannon': [],
                  'alpha_S': [], 'gamma_S': [], 'beta_S': [],
                  'alpha_E': [], 'gamma_E': [], 'beta_E': [],
                  'slope': 0} for name in space_structures}

for name, gen_A in space_structures.items():
    print(f"Simulating spatial structure: {name}")
    S_list, H_list = [], []
    for rep in range(n_repeat):
        np.random.seed(rep)
        # Species assignment
        patch_species_list = []
        N_k_list = []
        for k in range(K):
            prob = np.random.dirichlet(np.ones(N_pool) * np.random.uniform(0.5, 2.0))
            N_k = np.random.randint(3, 8)
            patch_species = np.sort(np.random.choice(np.arange(N_pool), N_k, replace=False, p=prob))
            patch_species_list.append(patch_species)
            N_k_list.append(N_k)
        # Community parameters (nondimensionalized)
        community_params = []
        for k in range(K):
            N_k = N_k_list[k]
            patch_species = patch_species_list[k]
            u_k = param.modular_uptake(N_k, M, N_modules=1, s_ratio=10.0) * C0 * R0 * T0
            l_k = param.generate_l_tensor(N_k, M, N_modules=2, s_ratio=10.0, λ=0.3)  # 假设l已无量纲
            rho_k = (np.ones(M) * (5.0 + 5.0 * np.sin(np.pi * k / (K-1)))) / (R0 / T0)
            omega_k = np.ones(M) * (0.05 + 0.01 * k) * T0
            m_k = np.ones(N_k) * (0.09 + 0.02 * k) * T0
            community_params.append({'N': N_k, 'M': M, 'u': u_k, 'l': l_k, 'rho': rho_k, 'omega': omega_k, 'm': m_k, 'species_idx': patch_species})
        # Initial conditions (nondimensionalized)
        n_total = K*N_pool + K*M
        Y0 = np.zeros(n_total)
        for k, cp in enumerate(community_params):
            N_k = cp['N']
            patch_species = cp['species_idx']
            Y0[k*N_pool:(k+1)*N_pool][patch_species] = np.random.uniform(0.5, 1.5, N_k) / C0
            Y0[K*N_pool + k*M:K*N_pool + (k+1)*M] = np.ones(M) * 2.0 / R0
        # Record which patches each species appears in
        species_in_patches = [[] for _ in range(N_pool)]
        for k, cp in enumerate(community_params):
            for sp in cp['species_idx']:
                species_in_patches[sp].append(k)
        # Spatial diffusion matrix
        A = gen_A(K)
        # ODE system (nondimensionalized)
        def global_ode(t, y):
            dydt = np.zeros_like(y)
            for k, cp in enumerate(community_params):
                N_k = cp['N']
                patch_species = cp['species_idx']
                u, l, rho, omega, m = cp['u'], cp['l'], cp['rho'], cp['omega'], cp['m']
                C_patch = y[k*N_pool:(k+1)*N_pool]
                R = y[K*N_pool + k*M:K*N_pool + (k+1)*M]
                for i, sp in enumerate(patch_species):
                    consumption = sum(C_patch[sp] * R[alpha] * u[i, alpha] * (1 - 0.3) for alpha in range(M))
                    mortality = C_patch[sp] * m[i]
                    dydt[k*N_pool + sp] = consumption - mortality
                for alpha in range(M):
                    input_rate = rho[alpha]
                    output_rate = R[alpha] * omega[alpha]
                    consumption = sum(C_patch[patch_species[i]] * R[alpha] * u[i, alpha] for i in range(N_k))
                    leakage = sum(sum(C_patch[patch_species[i]] * R[beta] * u[i, beta] * l[i, beta, alpha] for beta in range(M)) for i in range(N_k))
                    dydt[K*N_pool + k*M + alpha] = input_rate - output_rate - consumption + leakage
            # Species diffusion
            for i in range(N_pool):
                patch_list = species_in_patches[i]
                if len(patch_list) >= 2:
                    for k in patch_list:
                        for j in range(K):
                            if A[k, j] > 0 and i in community_params[j]['species_idx']:
                                C_k = y[k*N_pool + i]
                                C_j = y[j*N_pool + i]
                                dydt[k*N_pool + i] += D_species * A[k, j] * (C_j - C_k)
            # Resource diffusion
            for alpha in range(M):
                for k in range(K):
                    for j in range(K):
                        if A[k, j] > 0:
                            R_k_alpha = y[K*N_pool + k*M + alpha]
                            R_j_alpha = y[K*N_pool + j*M + alpha]
                            dydt[K*N_pool + k*M + alpha] += D0 * A[k, j] * (R_j_alpha - R_k_alpha)
            return dydt
        # Numerical integration
        from scipy.integrate import solve_ivp
        t_span = (0, 500)
        t_eval = np.linspace(*t_span, 500)
        sol = solve_ivp(global_ode, t_span, Y0, t_eval=t_eval, method='BDF', rtol=1e-5, atol=1e-7)
        # Calculate diversity indices
        # α-Shannon
        alpha_shannons = []
        alpha_Ss = []
        alpha_Es = []
        for k, cp in enumerate(community_params):
            patch_species = cp['species_idx']
            C_patch = sol.y[k*N_pool:(k+1)*N_pool, -1]
            abund = C_patch[patch_species]
            S = np.sum(abund > 0)
            H = shannon_index(abund)
            E = evenness(H, S)
            alpha_shannons.append(H)
            alpha_Ss.append(S)
            alpha_Es.append(E)
        alpha_shannon = np.mean(alpha_shannons)
        alpha_S = np.mean(alpha_Ss)
        alpha_E = np.mean(alpha_Es)
        # γ-Shannon
        all_abund = np.zeros(N_pool)
        for k, cp in enumerate(community_params):
            patch_species = cp['species_idx']
            C_patch = sol.y[k*N_pool:(k+1)*N_pool, -1]
            all_abund[patch_species] += C_patch[patch_species]
        gamma_S = np.sum(all_abund > 0)
        gamma_H = shannon_index(all_abund)
        gamma_E = evenness(gamma_H, gamma_S)
        # β-Shannon
        beta_shannon = gamma_H - alpha_shannon
        beta_S = gamma_S / alpha_S if alpha_S > 0 else 0
        beta_E = gamma_E - alpha_E
        # Save
        results[name]['alpha_shannon'].append(alpha_shannon)
        results[name]['gamma_shannon'].append(gamma_H)
        results[name]['beta_shannon'].append(beta_shannon)
        results[name]['alpha_S'].append(alpha_S)
        results[name]['gamma_S'].append(gamma_S)
        results[name]['beta_S'].append(beta_S)
        results[name]['alpha_E'].append(alpha_E)
        results[name]['gamma_E'].append(gamma_E)
        results[name]['beta_E'].append(beta_E)
        S_list.append(gamma_S)
        H_list.append(gamma_H)
    # Linear regression: Shannon vs S
    try:
        if len(set(S_list)) == 1:
            slope = np.nan
        else:
            slope, _, _, _, _ = linregress(S_list, H_list)
    except Exception as e:
        print(f'Linear regression failed: {e}')
        slope = np.nan
    results[name]['slope'] = slope
    results[name]['gamma_S_list'] = S_list
    results[name]['gamma_H_list'] = H_list

# ========== Results Output ==========
print("\n===== Mean ± Std of Diversity Indices under Different Spatial Structures =====")
for name in space_structures:
    print(f"\n{name}:")
    for key in ['alpha_shannon', 'gamma_shannon', 'beta_shannon', 'alpha_S', 'gamma_S', 'beta_S', 'alpha_E', 'gamma_E', 'beta_E']:
        arr = np.array(results[name][key])
        print(f"  {key}: {arr.mean():.3f} ± {arr.std():.3f}")
    print(f"  Shannon vs S slope: {results[name]['slope']:.3f}")

# ========== Visualization ==========
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

structure_names = list(space_structures.keys())
metrics = ['alpha_shannon', 'gamma_shannon', 'beta_shannon']

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
if not isinstance(axes, np.ndarray):
    axes = [axes]
for idx, metric in enumerate(metrics):
    ax = axes[idx]
    data = []
    for name in structure_names:
        for value in results[name][metric]:
            data.append({'structure': name, 'value': value})
    df = pd.DataFrame(data)
    sns.boxplot(x='structure', y='value', data=df, ax=ax, showfliers=True, width=0.5)
    sns.stripplot(x='structure', y='value', data=df, ax=ax, color='black', size=4, jitter=True, alpha=0.7)
    ax.set_title(metric)
    ax.set_ylabel('Value')
    ax.set_xlabel('')
    ax.set_xticklabels(structure_names, rotation=30, ha='right')
plt.tight_layout()
plt.show()