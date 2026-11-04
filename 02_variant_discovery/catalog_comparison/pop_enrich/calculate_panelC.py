#!/usr/bin/env python3
"""Process variant calls and summarize variant catalogs."""

import argparse
import pandas as pd
from cyvcf2 import VCF
from tqdm import tqdm
from collections import defaultdict

parser = argparse.ArgumentParser(description="Calculate individual variant burden per class")
parser.add_argument("--vcf", required=True, help="Input APC-P1-specific VCF.gz")
parser.add_argument("--meta", required=True, help="Sample info TSV (IID, Location)")
parser.add_argument("--vclass", required=True, help="Variant class TSV (CHROM, POS, REF, ALT, CLASS)")
parser.add_argument("--out", required=True, help="Output burden TSV")
args = parser.parse_args()

print("1. Loading sample metadata...")
meta_df = pd.read_csv(args.meta, sep="\t")
sample_to_loc = dict(zip(meta_df["IID"], meta_df["Location"]))

print("2. Loading variant class annotations...")
# Read Class Read annotations and build a hash lookup, Implement O(1) Lookup
class_df = pd.read_csv(args.vclass, sep="\t", usecols=["CHROM", "POS", "REF", "ALT", "CLASS"])
class_dict = {}
for _, row in class_df.iterrows():
    class_dict[(str(row["CHROM"]), int(row["POS"]), str(row["REF"]), str(row["ALT"]))] = row["CLASS"]

print("3. Processing VCF and counting burdens...")
vcf = VCF(args.vcf)
samples = list(vcf.header.samples)

# Initialize Burden Dictionary: burden[sample][class] = count
burden = {s: defaultdict(int) for s in samples}

for variant in tqdm(vcf):
    chrom = variant.CHROM
    pos = variant.POS
    ref = variant.REF

    # Handle multiallelic or decomposed sites
    for alt_idx, alt in enumerate(variant.ALT):
        var_key = (chrom, pos, ref, alt)
        var_class = class_dict.get(var_key, "UNKNOWN")
        if var_class == "UNKNOWN":
            continue

        # Extract the current ALT allele (In cyvcf2 In, ALT1 Is 1, ALT2 Is 2)
        target_allele = alt_idx + 1

        # Iterate over individual genotypes
        for i, gt in enumerate(variant.genotypes):
            # gt For example [0, 1, False] list
            # If this individual carries target_allele, Burden Then +1
            if target_allele in gt[:2]:
                burden[samples[i]][var_class] += 1

print("4. Formatting output...")
out_data = []
for s in samples:
    loc = sample_to_loc.get(s, "UNKNOWN")
    if loc == "UNKNOWN":
        continue
    for var_class, count in burden[s].items():
        out_data.append([s, loc, var_class, count])

out_df = pd.DataFrame(out_data, columns=["IID", "Location", "CLASS", "Burden"])
out_df.to_csv(args.out, sep="\t", index=False)
print(f"Done! Saved to {args.out}")
