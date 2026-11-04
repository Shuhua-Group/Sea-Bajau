#!/bin/bash

# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/hapNe"
cd "$wd" || exit 1

# Start with hapIBD Test
for min in 2 2.5
do
for max in 6 8 10
do
for list in `ls pop/*.list|grep -v 2018Cell|grep -v ho1926|tail -13`;do
ind=$(basename "${list}" .list)
mkdir -p 4panel_splitvcf/${ind}
# Step1: Split VCF
python step1.splitvcf.py \
--vcf_file 4panel_vcf/4panel.${ind}.chrA \
--save_in 4panel_splitvcf/${ind} \
--genome_build grch38
# Step2: Run hap-IBD
mkdir -p 4panel_hapIBD/${ind}
for file in 4panel_splitvcf/${ind}/*.vcf.gz
do
CHR=$(echo $file | awk -F'chr' '{print $2}' | awk -F'.' '{print $1}')
PREFIX=$(echo $file | awk -F'.vcf.gz' '{print $1}' | awk -F'/' '{print $NF}')
java -jar hap-ibd.jar gt=$file map=plink.chr$CHR.GRCh38.map out=IBD/$PREFIX
echo "Hap-IBD done for $file"
/path/to/java -jar \
/path/to/hap-ibd.jar \
gt=$file \
map=/path/to/plink.chr${CHR}.GRCh38.map \
out=4panel_hapIBD/${ind}/$PREFIX
done

# Step3: Filter IBD segments by length
python /path/to/filter_ibd_by_length.py \
    -i "4panel_hapIBD/${ind}/" \
    -o "4panel_hapIBD/${ind}_min${min}_max${max}cM/" \
    --min_cm "$min" \
    --max_cm "$max"

# Step4: Run HapNe-IBD
mkdir -p "4panel_HapNeIBD/${ind}_min${min}_max${max}cM"
n_sample=$(bcftools query -l 4panel_vcf/4panel.${ind}.chrA.vcf.gz | wc -l)
configfile="[CONFIG]
output_folder=4panel_HapNeIBD/${ind}_min${min}_max${max}cM
nb_samples=${n_sample}
population_name=$(basename 4panel_vcf/4panel.${ind}.chrA.vcf.gz .chrA.vcf.gz)
ibd_files=4panel_hapIBD/${ind}_min${min}_max${max}cM
column_cm_length=8
genome_build=grch38"
echo -e "$configfile" > "4panel_HapNeIBD/${ind}_min${min}_max${max}cM/config.txt"
python hapne_ibd.py --config_file "4panel_HapNeIBD/${ind}_min${min}_max${max}cM/config.txt"
done
done
done
