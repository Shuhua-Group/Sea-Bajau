"""Analyze digestive-enzyme copy numbers and haplotypes."""

import csv
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import itertools
import numpy as np


# === Read the mapping from sample to Superpop assignments ===
sample_to_superpop = {}
with open("/path/to/HGDP.info.list") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sample_to_superpop[row["Sample_ID"]] = row["Superpop"]

# === Read hap_pga_paths.txt, Collect each hap Of PGA Counts and assignments: Superpop ===
data = []
with open("/path/to/HGDP_statistic.merged_output.pattern.hap_PRSS_counts.txt") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) != 3:
            continue
        hap_name, copy3, copy5 = parts
        match = re.match(r"([A-Z0-9]+)_hap[12]", hap_name)
        print(match)
        if not match:
            continue
        sample_id = match.group(1)
        superpop = sample_to_superpop.get(sample_id)
        data.append({
            "Haplotype": hap_name,
            "Sample": sample_id,
            "Superpop": superpop,
            "3Copy": int(copy3) ,
            "5Copy": int(copy5),
        })

# === Convert to DataFrame Simplify processing ===
df = pd.DataFrame(data)

# === Calculate means ===
mean_df = df.groupby("Superpop")[["3Copy","5Copy"]].mean()
print("===  Each Superpop Mean  ===")
print(mean_df)

# === Statistical tests ===
print("\n=== Kruskal-Wallis Test results  ===")
for copy_type in ["3Copy","5Copy"]:
    groups = [group[copy_type].values for name, group in df.groupby("Superpop")]
    stat, pvalue = stats.kruskal(*groups)
    print(f"{copy_type}: H={stat:.3f}, p={pvalue:.3e}")

# === Save results ===
mean_df.to_csv("HGDP.superpop_prss_means.csv")
print("\n Save mean values to  HGDP.superpop_prss_means.csv")


# sns.set(style="whitegrid")
# fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# copy_types = ["3Copy", "5Copy"]

# # Calculate test statistics: p Value (Kruskal-Wallis)
# p_values = {}
# for copy in copy_types:
# groups = [group[copy].values for name, group in df.groupby("Superpop")]
# stat, p = stats.kruskal(*groups)
# p_values[copy] = p

# # Overlay box plots and individual points
# for i, pga in enumerate(pga_types):
# ax = axes[i]
# sns.boxplot(x="Superpop", y=pga, data=df, ax=ax, showfliers=False, palette="Set3")
# sns.stripplot(x="Superpop", y=pga, data=df, ax=ax, color='black', size=3, jitter=True, alpha=0.6)
# ax.set_title(f"{pga} by Superpop\nKruskal-Wallis p = {p_values[pga]:.3e}")
# ax.set_xlabel("Superpop")
# ax.set_ylabel("Count")

# plt.tight_layout()
# plt.savefig("PGA3_4_5_by_Superpop_boxplot.pdf", dpi=300)
