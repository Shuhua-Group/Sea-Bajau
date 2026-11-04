#!/bin/bash

# Keep only GujaratiA
# Process variant calls and summarize variant catalogs.
# Set input, reference and tool paths for your environment before execution.

cd /path/to/admixture || exit 1
plink --vcf admixture_20250108.102Malay.vcf.gz --vcf-half-call missing --geno 0.05 --bp-space 50000 --maf 0.05 --make-bed --out admixture_20250108.102Malay


for K in {1..10}; do ( /path/to/admixture --cv admixture_20250108.102Malay.bed $K | tee admixture_20250108.102Malay.admixture_log${K}.out &); done

