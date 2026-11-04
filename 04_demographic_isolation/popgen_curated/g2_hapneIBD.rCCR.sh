#!/bin/bash

# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/hapNe"

# Start with hapIBD Test

for list in `ls pop/*.list|grep -v 2018Cell|grep -v ho1926|grep -v 36Bajau`
do
ind=$(basename "${list}" .list)
mkdir -p "4panel_splitvcf/${ind}_36Bajau.combine"
# Step1: Split VCF
python step1.splitvcf.py \
--vcf_file 4panel_vcf/4panel.${ind}_36Bajau.combine.chrA \
--save_in 4panel_splitvcf/${ind}_36Bajau.combine \
--genome_build grch38
# Step2: Run germline2
mkdir -p "4panel_g2_m1/${ind}_36Bajau.combine"
for file in 4panel_splitvcf/${ind}_36Bajau.combine/*.vcf.gz
do
    CHR=$(echo "$file" | awk -F'chr' '{print $2}' | awk -F'.' '{print $1}')
    PREFIX=$(echo "$file" | awk -F'.vcf.gz' '{print $1}' | awk -F'/' '{print $NF}')
    /path/to/RelateFileFormats \
                --mode ConvertFromVcf \
                --haps "4panel_g2_m1/${ind}_36Bajau.combine/$PREFIX.haps" \
                --sample "4panel_g2_m1/${ind}_36Bajau.combine/$PREFIX.sample" \
                -i "4panel_splitvcf/${ind}_36Bajau.combine/$PREFIX"
    # Run g2, This step is fast
    /path/to/g2 \
    -h \
    -m 1 \
    "4panel_g2_m1/${ind}_36Bajau.combine/$PREFIX.haps" \
    "4panel_g2_m1/${ind}_36Bajau.combine/$PREFIX.sample" \
    "/path/to/chr${i}.genetic_map_hg38.txt2" \
    "4panel_g2_m1/${ind}_36Bajau.combine/$PREFIX"
    echo "GERMLINE2 done for $file"
    python germline2ibd.py "4panel_g2_m1/${ind}_36Bajau.combine/$PREFIX"
done

# Step3: Add a preprocessing step, Select ibd Segments shared across populations
mkdir -p "4panel_g2_m1/${ind}_36Bajau.combine_cross"
popA="pop/36Bajau.list"
popB="$list"
for file in 4panel_splitvcf/${ind}_36Bajau.combine/*.vcf.gz
do
    CHR=$(echo "$file" | awk -F'chr' '{print $2}' | awk -F'.' '{print $1}')
    PREFIX=$(echo "$file" | awk -F'.vcf.gz' '{print $1}' | awk -F'/' '{print $NF}')
    ibd_file=4panel_g2_m1/${ind}_36Bajau.combine/$PREFIX.ibd.gz
    cross_file=4panel_g2_m1/${ind}_36Bajau.combine_cross/$PREFIX.ibd.gz
    # Truncate the output file
    > "$cross_file"
    zcat "$ibd_file" | awk -v popA="$popA" -v popB="$popB" '
    BEGIN{
        while ((getline lineA < popA) > 0) a[lineA]=1
        while ((getline lineB < popB) > 0) b[lineB]=1
    }
    {
        if (($1 in a && $3 in b) || ($1 in b && $3 in a)) print $0
    }' | gzip -c > "$cross_file"
done

# Step4: Run HapNe-IBD
mkdir -p "4panel_g2_m1_HapNeIBD/${ind}_36Bajau.combine"
n_sample=$(bcftools query -l 4panel_vcf/4panel.${ind}_36Bajau.combine.chrA.vcf.gz | wc -l)
configfile="[CONFIG]
output_folder=4panel_g2_m1_HapNeIBD/${ind}_36Bajau.combine
nb_samples=${n_sample}
population_name=$(basename 4panel_vcf/4panel.${ind}_36Bajau.combine.chrA.vcf.gz .chrA.vcf.gz)
ibd_files=4panel_g2_m1/${ind}_36Bajau.combine_cross
column_cm_length=8
genome_build=grch38"
echo -e "$configfile" > "4panel_g2_m1_HapNeIBD/${ind}_36Bajau.combine/config.txt"
python hapne_ibd.py --config_file "4panel_g2_m1_HapNeIBD/${ind}_36Bajau.combine/config.txt"
done



for list in `ls pop/*.list|grep -v 2018Cell|grep -v ho1926|grep -v 36Bajau`;do
    ind=$(basename "${list}" .list)
    outdir="4panel_rCCR/${ind}_36Bajau.g2_m1.combine_cross"
    rm -rf "$outdir"
    mkdir -p "$outdir"
    hapneA="4panel_g2_m1_HapNeIBD/36Bajau/HapNe/hapne.csv"
    hapneB="4panel_g2_m1_HapNeIBD/${ind}/HapNe/hapne.csv"
    hapneABcross="4panel_g2_m1_HapNeIBD/${ind}_36Bajau.combine/HapNe/hapne.csv"
    # Write the header
    echo -e "Generation\t36Bajau_Ne\t${ind}_Ne\t${ind}_36Bajau_cross_Ne\t${ind}_36Bajau_rCCR" > "${outdir}/${ind}_36Bajau_rCCR.g2_m1_hapne.txt"
    # Calculate rCCR
    paste <(cut -d, -f1,2 "$hapneA" | sed 's/,/\t/') \
          <(cut -d, -f2 "$hapneB") \
          <(cut -d, -f2 "$hapneABcross") | \
          tail -100 | \
          awk '{rCCR = 2 / ($4 * (1/$2 + 1/$3)); print $1,$2,$3,$4,rCCR}' |tr " " "\t" >> "${outdir}/${ind}_36Bajau_rCCR.g2_m1_hapne.txt"
    cut -f5 "${outdir}/${ind}_36Bajau_rCCR.g2_m1_hapne.txt" > "${outdir}/${ind}_36Bajau.g2_m1_hapne.rCCR_only.txt"
done

paste 4panel_rCCR/*/*_36Bajau.g2_m1_hapne.rCCR_only.txt > 4panel_rCCR/4panel.g2_m1_hapne.rCCR_only.sum.txt

