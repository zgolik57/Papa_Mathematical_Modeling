import numpy as np
import pandas as pd
import itertools
from scipy.special import comb
from collections import defaultdict

#gamma parameters
gamma1_values = [0.25, 0.5, 1, 2, 4]
gamma2_values = [0.125, 0.5, 1, 2]

#valid pentamer compositions
compositions = []
for wt in range(6):
    for h in range(6 - wt):
        for f in range(6 - wt - h):
            if wt + h + f == 5:
                compositions.append((wt, h, f))

#compute circular distance
def circular_distance(i, j):
    return min(abs(i - j), 5 - abs(i - j))

#compute total gamma
def compute_gamma_total(config, gamma1, gamma2):
    total = 0
    for i, p1 in enumerate(config):
        if p1 != 'HALO':
            continue
        for j, p2 in enumerate(config):
            if p2 != 'FAST' or i == j:
                continue
            d = circular_distance(i, j)
            total += gamma1 if d == 1 else gamma2 if d == 2 else 0
    return total

#unique permutations for each composition
composition_perms = {}
for wt, h, f in compositions:
    items = ['WT'] * wt + ['HALO'] * h + ['FAST'] * f
    if len(items) == 5:
        perms = set(itertools.permutations(items))
        composition_perms[(wt, h, f)] = list(perms)

#sampling
n_samples = 100
std_dev = 0.05
conc_range = np.arange(0.0, 1.01, 0.05)
results = []

for wt_med in conc_range:
    for halo_med in conc_range:
        for fast_med in conc_range:
            total = wt_med + halo_med + fast_med
            if not np.isclose(total, 1.0):
                continue

            # Sample concentrations with variation
            wt_samples = np.clip(np.random.normal(wt_med, std_dev, n_samples), 0, 1)
            halo_samples = np.clip(np.random.normal(halo_med, std_dev, n_samples), 0, 1)
            fast_samples = np.clip(np.random.normal(fast_med, std_dev, n_samples), 0, 1)

            for gamma1 in gamma1_values:
                for gamma2 in gamma2_values:
                    if gamma1 <= gamma2:
                        continue

                    gvr_values = []
                    for i in range(n_samples):
                        wt = wt_samples[i]
                        halo = halo_samples[i]
                        fast = fast_samples[i]
                        norm_total = wt + halo + fast
                        if norm_total == 0:
                            continue
                        wt, halo, fast = wt / norm_total, halo / norm_total, fast / norm_total

                        weighted_signal = 0
                        for (wt_n, h_n, f_n), perms in composition_perms.items():
                            prob = (
                                comb(5, wt_n)
                                * comb(5 - wt_n, h_n)
                                * comb(5 - wt_n - h_n, f_n)
                                * (wt ** wt_n)
                                * (halo ** h_n)
                                * (fast ** f_n)
                            )

                            gamma_sum = sum(compute_gamma_total(p, gamma1, gamma2) for p in perms)
                            avg_gamma = gamma_sum / len(perms)
                            weighted_signal += prob * avg_gamma

                        gvr = weighted_signal / halo if halo > 0 else np.nan
                        gvr_values.append(gvr)

                    gvr_array = np.array(gvr_values)
                    results.append({
                        '[NPM1]_WT_median': wt_med,
                        '[NPM1]_Halo_median': halo_med,
                        '[NPM1]_FAST_median': fast_med,
                        'Gamma1': gamma1,
                        'Gamma2': gamma2,
                        'GVR_median': round(np.nanmedian(gvr_array), 3),
                        'GVR_std': round(np.nanstd(gvr_array), 3),
                    })

#output as df
df_gvr_var = pd.DataFrame(results)
df_gvr_var=df_gvr_var.round(3)
df_gvr_var.to_csv("npm1_gvr_variation_output.csv", index=False)