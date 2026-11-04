"""Summarize population variation and demographic history."""

import pandas as pd
import glob

# Population labels to process (Ensure .ibd Filenames contain these prefixes)
groups = ["7unpure.Bajau", "29pure.Bajau", "20Negrito.ngs", "33Malay"]

all_stats = []

for group in groups:
    ibd_file = f"{group}.chrA.QC.shapeited.map.hapIBD.ibd"
    try:
        # Load IBD Data
        colnames = ["sample1", "hap1", "sample2", "hap2", "chrom", "start", "end", "length_cM"]
        df = pd.read_csv(ibd_file, sep=r"\s+", header=None, names=colnames)
        df = df.dropna()
        # Add physical lengths
        df["length_bp"] = df["end"] - df["start"]
        # Pairs of individuals (Ignore haplotype Order)
        df["pair"] = df.apply(lambda row: tuple(sorted([row["sample1"], row["sample2"]])), axis=1)
        # Summarize
        pair_stats = df.groupby("pair").agg(
            ibd_segment_count = ("length_cM", "count"),
            total_ibd_length_cM = ("length_cM", "sum"),
            average_ibd_length_cM = ("length_cM", "mean"),
            total_ibd_length_bp = ("length_bp", "sum"),
            total_ibd_length_Mb = ("length_bp", lambda x: x.sum() / 1e6)
        ).reset_index()
        # Split pair
        pair_stats[["sample1", "sample2"]] = pd.DataFrame(pair_stats["pair"].tolist(), index=pair_stats.index)
        pair_stats["group"] = group
        pair_stats = pair_stats.drop(columns=["pair"])
        all_stats.append(pair_stats)
        print(f"[✓] Processed: {group}")
    except Exception as e:
        print(f"[✗] Failed for {group}: {e}")

# Merge results across populations
df_all = pd.concat(all_stats, ignore_index=True)
df_all.to_csv("all_groups.ibd_pairwise_stats.tsv", sep="\t", index=False)
print("[✓]  Merged output saved to : all_groups.ibd_pairwise_stats.tsv")
