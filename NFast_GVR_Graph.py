import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("split_fast_gvr_with_gamma.csv")

#delete null GVRs
df = df.dropna(subset=["GVR"])

#group by fast concentration
grouped = df.groupby("[NPM1]_NFAST")

#compute mean and std of GVR at each fast conc
summary = grouped["GVR"].agg(['mean', 'std']).reset_index()

#plot with error bars
plt.figure(figsize=(8, 6))
plt.errorbar(summary["[NPM1]_NFAST"], summary["mean"], yerr=summary["std"],
             fmt='o-', ecolor='gray', capsize=4, linewidth=2)

plt.xlabel("[NPM1]_NFAST")
plt.ylabel("Green-to-Violet Ratio (GVR)")
plt.title("GVR vs [NPM1]_NFAST (error bars from concentration variation and gamma)")
plt.grid(True)
plt.tight_layout()
plt.show()