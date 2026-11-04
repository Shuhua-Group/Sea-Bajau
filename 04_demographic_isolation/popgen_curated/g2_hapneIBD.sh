#!/bin/bash

# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/hapNe"

# Start with hapIBD Test

for list in `ls pop/*.list|grep -v 2018Cell|grep -v ho1926`;do
ind=$(basename "${list}" .list)
mkdir -p "4panel_splitvcf/${ind}"
# Step1: Split VCF
# python step1.splitvcf.py \
# --vcf_file 4panel_vcf/4panel.${ind}.chrA \
# --save_in 4panel_splitvcf/${ind} \
# --genome_build grch38
# Step2: Run germline2
mkdir -p "4panel_g2_m1/${ind}"
for file in 4panel_splitvcf/${ind}/*.vcf.gz
do
    CHR=$(echo "$file" | awk -F'chr' '{print $2}' | awk -F'.' '{print $1}')
    PREFIX=$(echo "$file" | awk -F'.vcf.gz' '{print $1}' | awk -F'/' '{print $NF}')
    /path/to/RelateFileFormats \
                --mode ConvertFromVcf \
                --haps "4panel_g2_m1/${ind}/$PREFIX.haps" \
                --sample "4panel_g2_m1/${ind}/$PREFIX.sample" \
                -i "4panel_splitvcf/${ind}/$PREFIX"
    # Run g2, This step is fast
    /path/to/g2 \
    -h \
    -m 2 \
    "4panel_g2_m1/${ind}/$PREFIX.haps" \
    "4panel_g2_m1/${ind}/$PREFIX.sample" \
    "/path/to/chr${i}.genetic_map_hg38.txt2" \
    "4panel_g2_m1/${ind}/$PREFIX"
    echo "GERMLINE2 done for $file"
    python germline2ibd.py "4panel_g2_m1/${ind}/$PREFIX"
done

# Step3: Run HapNe-IBD
mkdir -p "4panel_g2_m1_HapNeIBD/${ind}"
n_sample=$(bcftools query -l 4panel_vcf/4panel.${ind}.chrA.vcf.gz | wc -l)
configfile="[CONFIG]
output_folder=4panel_g2_m1_HapNeIBD/${ind}
nb_samples=${n_sample}
population_name=$(basename 4panel_vcf/4panel.${ind}.chrA.vcf.gz .chrA.vcf.gz)
ibd_files=4panel_g2_m1/${ind}
column_cm_length=8
genome_build=grch38"
echo -e "$configfile" > "4panel_g2_m1_HapNeIBD/${ind}/config.txt"
python hapne_ibd.py --config_file "4panel_g2_m1_HapNeIBD/${ind}/config.txt"
done
