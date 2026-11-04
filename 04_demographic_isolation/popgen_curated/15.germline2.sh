#!/bin/bash

# Use germline2 Identify IBD

# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

wd=/path/to/GERMLINE2
cd "$wd" || exit 1

for list in `ls pop/*.list|grep -v 2018Cell|grep -v ho1926`;do
ind=$(basename "${list}" .list)
mkdir -p input
for i in {1..22};do
bcftools view -S "${list}" "4panel.chr${i}.shapeit5.vcf.gz" \
-Oz -o "input/4panel.${ind}.chr${i}.vcf.gz"
/path/to/RelateFileFormats \
                --mode ConvertFromVcf \
                --haps "input/4panel.${ind}.chr${i}.haps" \
                --sample "input/4panel.${ind}.chr${i}.sample" \
                -i "input/4panel.${ind}.chr${i}"
# Run g2, This step is fast
mkdir -p 4panel_g2out
## -h is strict mode, hap Do not mix between trans
/path/to/g2 \
-h \
input/4panel.${ind}.chr${i}.haps \
input/4panel.${ind}.chr${i}.sample \
/path/to/chr${i}.genetic_map_hg38.txt2 \
4panel_g2out/4panel.${ind}.chr${i}

## Run an additional comparison without -h Of
/path/to/g2 \
input/4panel.${ind}.chr${i}.haps \
input/4panel.${ind}.chr${i}.sample \
/path/to/chr${i}.genetic_map_hg38.txt2 \
4panel_g2out/4panel.${ind}.chr${i}.nohap
done &
done

# g2 results compared with HapIBD differ
