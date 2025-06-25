import numpy as np
import pandas as pd
from scipy.special import comb

#compute pentamer compositions
compositions = []
for wt in range(6):
    for h in range(6 - wt):
        for nf in range(6 - wt - h):
            cf = 5 - wt - h - nf
            if 0 <= cf <= 5:
                compositions.append((wt, h, nf, cf))

compositions = np.array(compositions)
wt_vec = compositions[:, 0]
h_vec = compositions[:, 1]
nf_vec = compositions[:, 2]
cf_vec = compositions[:, 3]
multinomial_coeffs = comb(5, wt_vec) * comb(5 - wt_vec, h_vec) * comb(5 - wt_vec - h_vec, nf_vec)

#Split-FAST GVR calculation with gamma
def compute_split_fast_gvr(wt_conc, halo_conc, nfast_conc, cfast_conc, gamma=1.0, alpha=1.0, beta=1.0):
    total_conc = wt_conc + halo_conc + nfast_conc + cfast_conc
    if total_conc == 0:
        return np.nan

    p_wt = wt_conc / total_conc
    p_h = halo_conc / total_conc
    p_nf = nfast_conc / total_conc
    p_cf = cfast_conc / total_conc

    probs = multinomial_coeffs * (p_wt ** wt_vec) * (p_h ** h_vec) * (p_nf ** nf_vec) * (p_cf ** cf_vec)
    paired_fast = np.minimum(nf_vec, cf_vec)
    papa_weight = h_vec * paired_fast

    green_signal = beta * np.sum(papa_weight * probs)
    violet_signal = alpha * halo_conc

    gvr = gamma * (green_signal / violet_signal) if violet_signal > 0 else np.nan
    return gvr

# Sweep over concentrations and gamma values
step = 0.05
concentration_range = np.arange(0.0, 1.01, step)
gamma_values = [0.25, 0.5, 1.0, 2.0, 4.0]

results = []

for wt in concentration_range:
    for h in concentration_range:
        for nf in concentration_range:
            for cf in concentration_range:
                total = wt + h + nf + cf
                if np.isclose(total, 1.0):
                    for gamma in gamma_values:
                        gvr = compute_split_fast_gvr(wt, h, nf, cf, gamma=gamma)
                        results.append({
                            '[NPM1]_WT': wt,
                            '[NPM1]_Halo': h,
                            '[NPM1]_NFAST': nf,
                            '[NPM1]_CFAST': cf,
                            'Gamma': gamma,
                            'GVR': gvr
                        })

# Create and round DataFrame
df = pd.DataFrame(results)
df = df.round(3)

# Show or export
print(df.head())
df.to_csv("split_fast_gvr_with_gamma.csv", index=False)