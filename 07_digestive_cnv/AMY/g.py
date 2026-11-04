"""Analyze digestive-enzyme copy numbers and haplotypes."""

import csv
import re
from collections import defaultdict

# === Step 1: Load samples -> Superpop Mapping ===
sample_to_superpop = {}
with open("sampleinfo.Malaylike.txt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sample_to_superpop[row["Sample_ID"]] = row["Superpop"]

# === Step 2: Superpop -> List of (hap_name, PGA3, PGA4, PGA5) ===
superpop_data = defaultdict(list)

# === Step 3: Process hap-AMY File ===
with open("AMY.GT.out") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) != 2:
            continue
        hap_name, path_str = parts
        # Extract sample_id(Supports . Or _ Delimited)
        match = re.match(r"([A-Z0-9]+)[._]hap[12]", hap_name)
        if not match:
            continue
        sample_id = match.group(1)
        superpop = sample_to_superpop.get(sample_id)
        if not superpop:
            continue
        # Summarize PGA3/4/5 occurrences
        tags = path_str.split(",")
        counts = {"AMY1": 0, "AMY2A": 0, "AMY2Ap": 0, "AMY2B": 0}
        for tag in tags:
            if tag in counts:
                counts[tag] += 1
        superpop_data[superpop].append(
            (hap_name, counts["AMY1"], counts["AMY2A"], counts["AMY2Ap"], counts["AMY2B"])
        )

# === Step 4: Write each Superpop file ===
for superpop, records in superpop_data.items():
    output_file = f"{superpop}.tsv"
    with open(output_file, "w") as out:
        out.write("Haplotype\tAMY1\tAMY2A\tAMY2Ap\tAMY2B\n")
        for hap_name, a1, a2a, a2ap, a2b in records:
            out.write(f"{hap_name}\t{a1}\t{a2a}\t{a2ap}\t{a2b}\n")

print("✅  Each  Superpop  files generated successfully .")


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
with open("sampleinfo.txt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sample_to_superpop[row["Sample_ID"]] = row["Superpop"]

# === Read hap_pga_paths.txt, Collect each hap Of PGA Counts and assignments: Superpop ===
data = []
with open("AMY.GT.out") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) != 2:
            continue
        hap_name, path_str = parts
        match = re.match(r"([A-Z0-9]+)[._]hap[12]", hap_name)
        if not match:
            continue
        sample_id = match.group(1)
        superpop = sample_to_superpop.get(sample_id)
        if not superpop:
            continue
        tags = path_str.split(",")
        counts = {"AMY1": 0, "AMY2A": 0, "AMY2Ap": 0, "AMY2B": 0}
        for tag in tags:
            if tag in counts:
                counts[tag] += 1
        data.append({
            "Haplotype": hap_name,
            "Sample": sample_id,
            "Superpop": superpop,
            "AMY1": counts["AMY1"],
            "AMY2A": counts["AMY2A"],
            "AMY2Ap": counts["AMY2Ap"],
            "AMY2B": counts["AMY2B"],
        })

# === Convert to DataFrame Simplify processing ===
df = pd.DataFrame(data)

# === Calculate means ===
mean_df = df.groupby("Superpop")[["AMY1","AMY2A","AMY2Ap","AMY2B"]].mean()
print("===  Each Superpop Mean  ===")
print(mean_df)

# === Statistical tests ===
print("\n=== Kruskal-Wallis Test results  ===")
for pga_type in ["PGA3","PGA4","PGA5"]:
    groups = [group[pga_type].values for name, group in df.groupby("Superpop")]
    stat, pvalue = stats.kruskal(*groups)
    print(f"{pga_type}: H={stat:.3f}, p={pvalue:.3e}")

# === Save results ===
mean_df.to_csv("superpop_pga_means.csv")
print("\n Save mean values to  superpop_pga_means.csv")


sns.set(style="whitegrid")
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

pga_types = ["PGA3", "PGA4", "PGA5"]

# Calculate test statistics: p Value (Kruskal-Wallis)
p_values = {}
for pga in pga_types:
    groups = [group[pga].values for name, group in df.groupby("Superpop")]
    stat, p = stats.kruskal(*groups)
    p_values[pga] = p

# Overlay box plots and individual points
for i, pga in enumerate(pga_types):
    ax = axes[i]
    sns.boxplot(x="Superpop", y=pga, data=df, ax=ax, showfliers=False, palette="Set3")
    sns.stripplot(x="Superpop", y=pga, data=df, ax=ax, color='black', size=3, jitter=True, alpha=0.6)
    ax.set_title(f"{pga} by Superpop\nKruskal-Wallis p = {p_values[pga]:.3e}")
    ax.set_xlabel("Superpop")
    ax.set_ylabel("Count")

plt.tight_layout()
plt.savefig("PGA3_4_5_by_Superpop_boxplot.pdf", dpi=300)
