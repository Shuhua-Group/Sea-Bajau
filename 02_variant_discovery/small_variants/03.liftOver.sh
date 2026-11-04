#!/bin/bash

# Process variant calls and summarize variant catalogs.
# Set input, reference and tool paths for your environment before execution.

cd /path/to/67West_51East || exit 1

bcftools query -f '%CHROM\t%POS0\t%POS\t%ID\n' 67West_51East.chm13v2.dv161.chr1_22.snp.vcf.gz \
> 67West_51East.chm13v2.dv161.chr1_22.snp.bed

# Use CHM13 Convert coordinates to GRCh38, Compare sites with the reference populations
liftOver 67West_51East.chm13v2.dv161.chr1_22.snp.bed \
/path/to/chm13v2-grch38.chain \
67West_51East.grch38.dv161.chr1_22.snp.bed \
67West_51East.chm13v2togrch38.unmap.dv161.chr1_22.snp.bed

## Number of sites
# CHM13: 16915773
# GRCh38: 12750893
awk 'NR==FNR {coords[$4]=$2; next} $1 ~ /^#/ || $3 in coords {if ($1 !~ /^#/) $2=coords[$3]; print}' \
    67West_51East.grch38.dv161.chr1_22.snp.bed <(zcat 67West_51East.chm13v2.dv161.chr1_22.snp.vcf.gz) > 67West_51East.grch38.dv161.chr1_22.snp.vcf


# This step requires reheader, Intermediate steps omitted
bgzip 67West_51East.grch38.dv161.chr1_22.snp.vcf
tabix -f 67West_51East.grch38.dv161.chr1_22.snp.vcf.gz


# Alternatively use crossmap Convert directly


