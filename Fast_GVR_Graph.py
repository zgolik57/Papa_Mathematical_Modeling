import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("npm1_gvr_gamma_output.csv")

#delete null GVRs
df = df.dropna(subset=["GVR"])

#group by fast concentration
grouped = df.groupby("Gamma")

#compute mean and std of GVR at each fast conc
summary = grouped["GVR"].agg(['mean', 'std']).reset_index()

#plot with error bars
plt.figure(figsize=(8, 6))
plt.errorbar(summary["Gamma"], summary["mean"], yerr=summary["std"],
             fmt='o-', ecolor='gray', capsize=4, linewidth=2)

plt.xlabel("Gamma constant value")
plt.ylabel("Green-to-Violet Ratio (GVR)")
plt.title("GVR vs Gamma value (error bars from concentration variation)")
plt.grid(True)
plt.tight_layout()
plt.show()