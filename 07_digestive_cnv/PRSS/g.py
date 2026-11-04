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

# === Step 2: Superpop -> List of (hap_name, 5copy, 3copy) ===
superpop_data = defaultdict(list)

# === Step 3: Process hap-PGA File ===
with open("PRSS.GT.out3") as f:
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
        # Summarize 3copy/5copy occurrences path_str Format: 0,5copy_format1,2
        # Determine whether this row is 5copy Still 3copy
        if "5copy" in path_str:
            p5 = 1
        else:
            p5 = 0
        if "3copy" in path_str:
            p3 = 1
        else:
            p3 = 0
        superpop_data[superpop].append(
            (hap_name, p3, p5))



# === Step 4: Write each Superpop file ===
for superpop, records in superpop_data.items():
    output_file = f"{superpop}.tsv"
    with open(output_file, "w") as out:
        out.write("Haplotype\t3copy\t5copy\n")
        for hap_name, p3, p5 in records:
            out.write(f"{hap_name}\t{p3}\t{p5}\n")

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
with open("sampleinfo.Malaylike.txt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sample_to_superpop[row["Sample_ID"]] = row["Superpop"]

# === Read hap_pga_paths.txt, Collect each hap Of PGA Counts and assignments: Superpop ===
data = []
with open("PRSS.GT.out3") as f:
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
        if "5copy" in path_str:
            p5 = 1
        else:
            p5 = 0
        if "3copy" in path_str:
            p3 = 1
        else:
            p3 = 0
        data.append({
            "Haplotype": hap_name,
            "Sample": sample_id,
            "Superpop": superpop,
            "3copy": p3,
            "5copy": p5,
        })

# === Convert to DataFrame Simplify processing ===
df = pd.DataFrame(data)
df.to_csv("PRSS.GT.out3.tsv", index=False, sep="\t")
# === Calculate means ===
mean_df = df.groupby("Superpop")[["3copy","5copy"]].mean()
print("===  Each Superpop Mean  ===")
print(mean_df)

# === Statistical tests ===
print("\n=== Kruskal-Wallis Test results  ===")
for pga_type in ["3copy","5copy"]:
    groups = [group[pga_type].values for name, group in df.groupby("Superpop")]
    stat, pvalue = stats.kruskal(*groups)
    print(f"{pga_type}: H={stat:.3f}, p={pvalue:.3e}")

# === Save results ===
mean_df.to_csv("superpop_pga_means.csv")
print("\n Save mean values to  superpop_pga_means.csv")


sns.set(style="whitegrid")
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

pga_types = ["3copy", "5copy"]

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
plt.savefig("3_5_by_Superpop_boxplot.pdf", dpi=300)
