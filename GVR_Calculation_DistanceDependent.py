import numpy as np
import pandas as pd
import itertools
from scipy.special import comb

#define gamma values
gamma1_values = [0.25, 0.5, 1, 2, 4]
gamma2_values = [0.125, 0.25, 0.5, 1, 2]

#generate pentamer compositions of WT, HALO, FAST
compositions = []
for wt in range(6):
    for h in range(6 - wt):
        for f in range(6 - wt - h):
            if wt + h + f == 5:
                compositions.append((wt, h, f))

#compute circular distance in a pentamer
def circular_distance(i, j):
    return min(abs(i - j), 5 - abs(i - j))

#compute gamma_total for a given configuration and gamma1/gamma2
def compute_gamma_total(config, gamma1, gamma2):
    gamma_sum = 0
    for i, pos1 in enumerate(config):
        if pos1 != 'HALO':
            continue
        for j, pos2 in enumerate(config):
            if pos2 != 'FAST' or i == j:
                continue
            dist = circular_distance(i, j)
            if dist == 1:
                gamma_sum += gamma1
            elif dist == 2:
                gamma_sum += gamma2
    return gamma_sum

#sample concentrations
results = []
conc_range = np.arange(0.0, 1.01, 0.05)

for wt_conc in conc_range:
    for halo_conc in conc_range:
        for fast_conc in conc_range:
            total = wt_conc + halo_conc + fast_conc
            if not np.isclose(total, 1.0):
                continue

            p_wt = wt_conc
            p_halo = halo_conc
            p_fast = fast_conc

            for gamma1 in gamma1_values:
                for gamma2 in gamma2_values:
                    if gamma1 <= gamma2:
                        continue

                    weighted_signal = 0
                    for wt_count, halo_count, fast_count in compositions:
                        p_config = (
                            comb(5, wt_count)
                            * comb(5 - wt_count, halo_count)
                            * comb(5 - wt_count - halo_count, fast_count)
                            * (p_wt ** wt_count)
                            * (p_halo ** halo_count)
                            * (p_fast ** fast_count)
                        )

                        items = ['WT'] * wt_count + ['HALO'] * halo_count + ['FAST'] * fast_count
                        if len(items) != 5:
                            continue

                        perms = set(itertools.permutations(items))
                        gamma_total_sum = sum(compute_gamma_total(p, gamma1, gamma2) for p in perms)
                        avg_gamma_total = gamma_total_sum / len(perms) if perms else 0
                        weighted_signal += avg_gamma_total * p_config

                    violet_signal = halo_conc
                    gvr = weighted_signal / violet_signal if violet_signal > 0 else np.nan

                    results.append({
                        '[NPM1]_WT': wt_conc,
                        '[NPM1]_Halo': halo_conc,
                        '[NPM1]_FAST': fast_conc,
                        'Gamma1': gamma1,
                        'Gamma2': gamma2,
                        'GVR': round(gvr, 3)
                    })

#convert to df
df_distance_weighted = pd.DataFrame(results)
df_distance_weighted = df_distance_weighted.round(3)

df_distance_weighted.to_csv("npm1_gvr_distance_output.csv", index=False)