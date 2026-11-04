"""Analyze mitochondrial and Y-chromosome lineage variation."""

import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import rcParams

# Set embedding mode to TrueType Font (Type 42)
rcParams['pdf.fonttype'] = 42

# Input file
hap_file = "/path/to/25malay_8bornean_29purebajau_7mixedbajau.GRCh38.dv161.chrY.YLT.hapresult.hg"

# Read
df_hap = pd.read_csv(hap_file, sep="\t")
df_info = pd.read_csv(info_file)
df = df_hap.merge(df_info, left_on="SampleID", right_on="IID")

# Extract hierarchy levels
df["Level1"] = df["Haplogroup"].str[0:1]
df["Level2"] = df["Haplogroup"].str[0:2]
df["Level3"] = df["Haplogroup"].str[0:3]
df["Level4"] = df["cluster2"]

df_counts = df.groupby(["Level1", "Level2", "Level3", "Level4"]).size().reset_index(name="Count")
total = df_counts["Count"].sum()

cluster_colors = {"Malay":"#ED9364", "Sea-Bajau":"#4A68D9", "Bornean":"#C54A53"}
default_colors = plt.cm.tab20.colors

# Level4 Outer ring
outer_labels, outer_sizes, outer_colors = [], [], []
for _, row in df_counts.iterrows():
    outer_labels.append(f"{row['Level4']}\n{row['Count']/total*100:.1f}%")
    outer_sizes.append(row["Count"])
    outer_colors.append(cluster_colors.get(row["Level4"], "grey"))

# Level3
mid3_counts = df_counts.groupby("Level3")["Count"].sum().reset_index()
mid3_labels = [f"{row['Level3']}\n{row['Count']/total*100:.1f}%" for _, row in mid3_counts.iterrows()]
mid3_sizes = mid3_counts["Count"].tolist()
mid3_colors = default_colors[:len(mid3_sizes)]

# Level2
mid2_counts = df_counts.groupby("Level2")["Count"].sum().reset_index()
mid2_labels = [f"{row['Level2']}\n{row['Count']/total*100:.1f}%" for _, row in mid2_counts.iterrows()]
mid2_sizes = mid2_counts["Count"].tolist()
mid2_colors = default_colors[len(mid3_sizes):len(mid3_sizes)+len(mid2_counts)]

# Level1
inner_counts = df_counts.groupby("Level1")["Count"].sum().reset_index()
inner_labels = [f"{row['Level1']}\n{row['Count']/total*100:.1f}%" for _, row in inner_counts.iterrows()]
inner_sizes = inner_counts["Count"].tolist()
inner_colors = default_colors[len(mid3_sizes)+len(mid2_counts):len(mid3_sizes)+len(mid2_counts)+len(inner_counts)]

fig, ax = plt.subplots(figsize=(9,9))

# Inner ring Level1
ax.pie(inner_sizes, labels=inner_labels, colors=inner_colors,
       radius=0.3, wedgeprops=dict(width=0.2, edgecolor='white'))

# Second ring Level2
ax.pie(mid2_sizes, labels=mid2_labels, colors=mid2_colors,
       radius=0.5, wedgeprops=dict(width=0.2, edgecolor='white'))

# Third ring Level3
ax.pie(mid3_sizes, labels=mid3_labels, colors=mid3_colors,
       radius=0.7, wedgeprops=dict(width=0.2, edgecolor='white'))

# Outer ring Level4
ax.pie(outer_sizes, labels=outer_labels, colors=outer_colors,
       radius=0.9, wedgeprops=dict(width=0.15, edgecolor='white'))

ax.set(aspect="equal")
#plt.title("Y-haplogroup with 4 levels", fontsize=16)

plt.savefig("/path/to/Y_haplogroup_clusters.pdf", bbox_inches="tight")
plt.close()
