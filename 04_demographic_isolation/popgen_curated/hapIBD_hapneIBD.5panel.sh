#!/bin/bash

# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/hapNe"
cd "${wd}" || exit 1
# Start with hapIBD Test
for list in `ls pop/*.list|grep -v ho1926`;do
ind=$(basename "${list}" .list)
bcftools concat \
/path/to/5panel.${ind}.chr*.vcf.gz \
-Oz -o "5panel_vcf/5panel.${ind}.chrA.vcf.gz"
done

for list in `ls pop/*.list|grep -v ho1926|head -16|tail -3`;do
ind=$(basename "${list}" .list)
mkdir -p "5panel_splitvcf/${ind}"
# Step1: Split VCF
python step1.splitvcf.py \
--vcf_file 5panel_vcf/5panel.${ind}.chrA \
--save_in 5panel_splitvcf/${ind} \
--genome_build grch38

# Step2: Run hap-IBD
mkdir -p 5panel_hapIBD/${ind}
for file in 5panel_splitvcf/${ind}/*.vcf.gz
do
CHR=$(echo $file | awk -F'chr' '{print $2}' | awk -F'.' '{print $1}')
PREFIX=$(echo $file | awk -F'.vcf.gz' '{print $1}' | awk -F'/' '{print $NF}')
# java -jar hap-ibd.jar gt=$file map=plink.chr$CHR.GRCh38.map out=IBD/$PREFIX
echo "Hap-IBD done for $file"
/path/to/java -jar \
/path/to/hap-ibd.jar \
gt=$file \
map=/path/to/plink.chr${CHR}.GRCh38.map \
out=5panel_hapIBD/${ind}/$PREFIX
done


# Step3: Run HapNe-IBD
mkdir -p "5panel_HapNeIBD/${ind}"
n_sample=$(bcftools query -l 5panel_vcf/5panel.${ind}.chrA.vcf.gz | wc -l)
configfile="[CONFIG]
output_folder=5panel_HapNeIBD/${ind}
nb_samples=${n_sample}
population_name=$(basename 5panel_vcf/5panel.${ind}.chrA.vcf.gz .chrA.vcf.gz)
ibd_files=5panel_hapIBD/${ind}
column_cm_length=8
genome_build=grch38"
echo -e "$configfile" > "5panel_HapNeIBD/${ind}/config.txt"

python hapne_ibd.py --config_file "5panel_HapNeIBD/${ind}/config.txt"
done



python hapNeIBD.plot.multi.py -i \
5panel_HapNeIBD/2018Cell.Bajau/HapNe/hapne.csv \
5panel_HapNeIBD/2018Cell.Saluan/HapNe/hapne.csv \
5panel_HapNeIBD/25Malay/HapNe/hapne.csv \
5panel_HapNeIBD/36Bajau/HapNe/hapne.csv \
5panel_HapNeIBD/EGDP.Borneo/HapNe/hapne.csv \
5panel_HapNeIBD/EGDP.Luzon/HapNe/hapne.csv \
5panel_HapNeIBD/EGDP.Myanmar/HapNe/hapne.csv \
5panel_HapNeIBD/EGDP.Vietnam/HapNe/hapne.csv \
5panel_HapNeIBD/Pacific.Barito/HapNe/hapne.csv \
5panel_HapNeIBD/Pacific.Besemah/HapNe/hapne.csv \
5panel_HapNeIBD/Pacific.Kalimantan/HapNe/hapne.csv \
5panel_HapNeIBD/Pacific.Mamanwa/HapNe/hapne.csv \
5panel_HapNeIBD/Pacific.Manobo/HapNe/hapne.csv \
5panel_HapNeIBD/Pacific.Semende/HapNe/hapne.csv \
5panel_HapNeIBD/PanMLS.Jakun/HapNe/hapne.csv \
5panel_HapNeIBD/PanMLS.Temiar/HapNe/hapne.csv \
-o 5panel_HapNeIBD.ne.pdf --min_time 1 --max_time 100 --figsize 10 7


for list in `ls pop/*.list|grep -v ho1926|grep -v 36Bajau`;do
ind=$(basename "${list}" .list)
python hapNeIBD.plot.multi.py -i \
5panel_HapNeIBD/36Bajau/HapNe/hapne.csv \
5panel_HapNeIBD/${ind}/HapNe/hapne.csv \
-o "5panel_HapNeIBD.${ind}_36Bajau.ne.pdf" --min_time 1 --max_time 100 --figsize 10 7
done
