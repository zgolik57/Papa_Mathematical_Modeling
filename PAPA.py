import numpy as np
import pandas as pd
from scipy.special import comb

#compute valid pentamer combos
compositions = []
for w in range(6):
    for h in range(6 - w):
        f = 5 - w - h
        compositions.append((w, h, f))

compositions = np.array(compositions)
w_vec = compositions[:, 0]
h_vec = compositions[:, 1]
f_vec = compositions[:, 2]
multinomial_coeffs = comb(5, w_vec) * comb(5 - w_vec, h_vec)

#compute GVR (before gamma scaling)
def compute_base_gvr(wt_conc, halo_conc, fast_conc, alpha=1.0, beta=1.0):
    total_conc = wt_conc + halo_conc + fast_conc
    if total_conc == 0:
        return np.nan

    #probabilities
    p_w = wt_conc / total_conc
    p_h = halo_conc / total_conc
    p_f = fast_conc / total_conc

    #multinomial probabilities and PAPA weights
    probs = multinomial_coeffs * (p_w ** w_vec) * (p_h ** h_vec) * (p_f ** f_vec)
    w_papa = h_vec * f_vec

    #G and V signal
    V = alpha * halo_conc
    G = beta * np.sum(w_papa * probs)

    return G / V if V > 0 else np.nan

#sweep over concentrations and gamma values
step = 0.05
concentration_range = np.arange(0.0, 1.01, step)
gamma_values = [0.25, 0.5, 1.0, 2.0, 4.0]

rows = []

for wt in concentration_range:
    for halo in concentration_range:
        for fast in concentration_range:
            if np.isclose(wt + halo + fast, 1.0):
                base_gvr = compute_base_gvr(wt, halo, fast)
                for gamma in gamma_values:
                    scaled_gvr = gamma * base_gvr if not np.isnan(base_gvr) else np.nan
                    rows.append({
                        '[NPM1]_WT': wt,
                        '[NPM1]_Halo': halo,
                        '[NPM1]_FAST': fast,
                        'Gamma': gamma,
                        'GVR': scaled_gvr
                    })

#convert to DataFrame
gvr_df = pd.DataFrame(rows)
gvr_df = gvr_df.round(3)


#preview
#print(gvr_df.head())

#save to CSV
gvr_df.to_csv("npm1_gvr_gamma_output.csv", index=False)