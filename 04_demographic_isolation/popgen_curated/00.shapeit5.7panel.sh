#!/bin/bash

# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

conda activate shapeit5_env
cd /path/to/re_merge || exit 1

for i in {1..22};do
SHAPEIT5_phase_common --input 7panel.chr${i}.vcf.gz \
--region chr${i} \
--map /path/to/chr${i}.genetic_map_hg38.txt2 \
--output 7panel.chr${i}.shapeit5.bcf \
--thread 32
bcftools view "7panel.chr${i}.shapeit5.bcf" -Oz -o "7panel.chr${i}.shapeit5.vcf.gz"
tabix -p vcf "7panel.chr${i}.shapeit5.vcf.gz"
done

bcftools concat 7panel.chr*.shapeit5.vcf.gz -Oz -o 7panel.chrA.shapeit5.vcf.gz
tabix -p vcf 7panel.chrA.shapeit5.vcf.gz
