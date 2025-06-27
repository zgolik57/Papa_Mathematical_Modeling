import numpy as np
import pandas as pd #^numerical tools and dataframes
import itertools #generates unique permutations
from scipy.special import comb #combinatorics choose function
from multiprocessing import Pool, cpu_count #allows parallelization
from tqdm import tqdm #progress bar

#define 2D positions on unit circle for pentamer
angles = [2 * np.pi * i / 5 for i in range(5)]
circle_coords = [(np.cos(a), np.sin(a)) for a in angles]

#assign gamma based on distance from midpoint of fast pair to halo
def assign_gamma(distance, gamma_map):
    if distance < 1.0:
        return gamma_map['gamma1']
    elif distance < 1.3:
        return gamma_map['gamma2']
    elif distance < 1.5:
        return gamma_map['gamma3']
    else:
        return gamma_map['gamma4']

# Compute total gamma for a given permutation
def compute_total_gamma(pentamer, gamma_map):
    halo_idx = [i for i, x in enumerate(pentamer) if x == 'HALO']
    nfast_idx = [i for i, x in enumerate(pentamer) if x == 'NFAST']
    cfast_idx = [i for i, x in enumerate(pentamer) if x == 'CFAST']
#^collects indicies of halo, nfast , and cfast, if not all are present, GVR=0
    if not halo_idx or not nfast_idx or not cfast_idx:
        return 0

    min_len = min(len(nfast_idx), len(cfast_idx))
    cfast_perms = list(itertools.permutations(cfast_idx, min_len))

    total_gamma_across_perms = 0
    for c_perm in cfast_perms:
        pairings = list(zip(nfast_idx[:min_len], c_perm))
        gamma_sum = 0
        for nf, cf in pairings:
            midpoint_idx = (nf + cf) / 2 % 5
            mid_angle = 2 * np.pi * midpoint_idx / 5
            midpoint_pos = (np.cos(mid_angle), np.sin(mid_angle))
            for h in halo_idx:
                halo_pos = circle_coords[h]
                dist = np.linalg.norm(np.array(halo_pos) - np.array(midpoint_pos))
                gamma_sum += assign_gamma(dist, gamma_map)
        total_gamma_across_perms += gamma_sum

    avg_gamma = total_gamma_across_perms / len(cfast_perms)
    return avg_gamma

#generate all valid compositions of a pentamer
def generate_compositions():
    comps = []
    for wt in range(6):
        for h in range(6 - wt):
            for nf in range(6 - wt - h):
                for cf in range(6 - wt - h - nf):
                    if wt + h + nf + cf == 5:
                        comps.append((wt, h, nf, cf))
    return comps
#^loops thru and stores combinations of wt, h, nf and cf that add to 5

#cache all unique permutations for each composition
def generate_permutations(compositions):
    composition_perms = {}
    for wt, h, nf, cf in compositions:
        tags = ['WT'] * wt + ['HALO'] * h + ['NFAST'] * nf + ['CFAST'] * cf
        if len(tags) == 5:
            composition_perms[(wt, h, nf, cf)] = list(set(itertools.permutations(tags)))
    return composition_perms

#function for a single gamma/concentration combination
def compute_gvr(args):
    wt, h, nf, cf, gamma1, gamma2, gamma3, gamma4, composition_perms = args
    gamma_map = {'gamma1': gamma1, 'gamma2': gamma2, 'gamma3': gamma3, 'gamma4': gamma4}
    total_signal = 0

    for (wt_n, h_n, nf_n, cf_n), perms in composition_perms.items():
        if not perms:
            continue
        prob = (
            comb(5, wt_n)
            * comb(5 - wt_n, h_n)
            * comb(5 - wt_n - h_n, nf_n)
            * (wt ** wt_n)
            * (h ** h_n)
            * (nf ** nf_n)
            * (cf ** cf_n)
        )

        gamma_sum = sum(compute_total_gamma(p, gamma_map) for p in perms)
        avg_gamma = gamma_sum / len(perms)
        total_signal += prob * avg_gamma

    gvr = total_signal / h if h > 0 else np.nan
    return {
        '[WT]': wt, '[HALO]': h, '[NFAST]': nf, '[CFAST]': cf,
        'Gamma1': gamma1, 'Gamma2': gamma2, 'Gamma3': gamma3, 'Gamma4': gamma4,
        'GVR': round(gvr, 3)
    }

#main execution
if __name__ == '__main__':
    concentrations = np.round(np.arange(0.0, 1.01, 0.1), 2)
    gamma_vals = np.round(np.arange(0.2, 2.01, 0.2), 2)
    compositions = generate_compositions()
    composition_perms = generate_permutations(compositions)

    tasks = []
    for wt in concentrations:
        for h in concentrations:
            for nf in concentrations:
                for cf in concentrations:
                    if np.isclose(wt + h + nf + cf, 1.0):
                        for g1 in gamma_vals:
                            for g2 in gamma_vals:
                                for g3 in gamma_vals:
                                    for g4 in gamma_vals:
                                        tasks.append((wt, h, nf, cf, g1, g2, g3, g4, composition_perms))

    print(f"Total tasks: {len(tasks)} — running on {cpu_count()} cores...")

    with Pool(cpu_count()) as pool:
        results = list(tqdm(pool.imap(compute_gvr, tasks), total=len(tasks)))

    df = pd.DataFrame(results).round(3)
    df.to_csv("papa_gvr_output.csv", index=False)
    print("Done. Results saved to papa_gvr_output.csv")