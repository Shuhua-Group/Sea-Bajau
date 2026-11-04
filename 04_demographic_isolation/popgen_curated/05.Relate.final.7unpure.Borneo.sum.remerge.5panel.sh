#!/bin/bash

# Run RELATE Estimate effective population size
# Select the target populations, Merge again

# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/re_merge_5panel"
mkdir -p "$wd"
cd "$wd" || exit 1
for i in {1..22};do
# Note:, All samples must index
# Set1. This study
set1="/path/to/102Malay.GRCh38.dv161.chr${i}.norm.QC.SNP.shapeited.reID2.vcf.gz"
# Set2. Ho1926
# set2="/path/to/Ho.1926.merge.chr${i}.vcf.gz"
# Set3. PanMLS
set3="/path/to/PanMLS461.VQSR_PASS.SNP.biallelic.universal_mask_low_complexity_heng75.chr1-22.hg38.chr${i}.sorted.vcf.gz"
# Set4. EGDP
set4="/path/to/EGDP_PaganiEtAl2016_Release_Build37_Chr${i}.94sample.hg38.sort.vcf.gz"
# # Set5. 2018Cell, Exclude for now
set5="/path/to/90combined.raw_variants.snps.PASS.hg38.chr${i}.vcf.gz"
# Set6. Pacific
set6="/path/to/Pacific.uni.hg38.chr${i}.sorted.selected.vcf.gz"
tmp_dir="${wd}/tmp_dir/chr${i}"

mkdir -p "${tmp_dir}"
# Extract the samples
for set in set1 set3 set4 set5 set6;do
bcftools view -S ${set}.list \
${!set} \
-Oz -o ${tmp_dir}/${set}.chr${i}.selected.vcf.gz \
&& tabix -p vcf ${tmp_dir}/${set}.chr${i}.selected.vcf.gz &
done

/path/to/bcftools isec --threads 4 \
${set1} ${set3} ${set4} ${set5} ${set6} \
-n=5 -p ${tmp_dir} &
bgzip ${tmp_dir}/0000.vcf -c > ${tmp_dir}/0000.vcf.gz && tabix -f ${tmp_dir}/0000.vcf.gz &
bgzip ${tmp_dir}/0001.vcf -c > ${tmp_dir}/0001.vcf.gz && tabix -f ${tmp_dir}/0001.vcf.gz &
bgzip ${tmp_dir}/0002.vcf -c > ${tmp_dir}/0002.vcf.gz && tabix -f ${tmp_dir}/0002.vcf.gz &
bgzip ${tmp_dir}/0003.vcf -c > ${tmp_dir}/0003.vcf.gz && tabix -f ${tmp_dir}/0003.vcf.gz &
bgzip ${tmp_dir}/0004.vcf -c > ${tmp_dir}/0004.vcf.gz && tabix -f ${tmp_dir}/0004.vcf.gz &

/path/to/bcftools merge "${tmp_dir}/0000.vcf.gz" "${tmp_dir}/0001.vcf.gz" \
${tmp_dir}/0002.vcf.gz "${tmp_dir}/0003.vcf.gz" "${tmp_dir}/0004.vcf.gz" \
--force-samples \
-Oz -o "${wd}/5panel.chr${i}.vcf.gz"  \
&& tabix -p vcf "${wd}/5panel.chr${i}.vcf.gz" &
done

# Merge all chromosomes
bcftools concat 5panel.chr*.vcf.gz \
-Oz -o 5panel.chrA.vcf.gz
tabix -f 5panel.chrA.vcf.gz
# 86885 sites, Approximate ..

conda activate shapeit5_env
for i in {1..22};do
SHAPEIT5_phase_common --input "5panel.chr${i}.vcf.gz" \
--region "chr${i}" \
--map "/path/to/chr${i}.genetic_map_hg38.txt2" \
--output "5panel.chr${i}.shapeit5.bcf" \
--thread 32
bcftools view "5panel.chr${i}.shapeit5.bcf" -Oz -o "5panel.chr${i}.shapeit5.vcf.gz"
tabix -p vcf "5panel.chr${i}.shapeit5.vcf.gz"
done &

bcftools concat 5panel.chr*.shapeit5.vcf.gz -Oz -o 5panel.chrA.shapeit5.vcf.gz
tabix -p vcf 5panel.chrA.shapeit5.vcf.gz

# cd ../
# # Calculate first IBD
mkdir -p ibd_5panel


for list in `ls re_merge_5panel/*list`;do
ind=$(basename ${list} .list)
for i in {1..22};do
bcftools view -S ${list} \
re_merge_5panel/5panel.chr${i}.shapeit5.vcf.gz \
-Oz -o ibd_5panel/5panel.chr${i}.shapeit5.${ind}.vcf.gz
tabix -p vcf ibd_5panel/5panel.chr${i}.shapeit5.${ind}.vcf.gz
plink --vcf ibd_5panel/5panel.chr${i}.shapeit5.${ind}.vcf.gz \
--cm-map /path/to/chr${i}.genetic_map_hg38.txt2 chr${i} \
--double-id \
--out ibd_5panel/5panel.chr${i}.shapeit5.${ind}.map \
--recode \
--set-missing-var-ids @:#

sed -i "s/^/chr/g" ibd_5panel/5panel.chr${i}.shapeit5.${ind}.map.map
# /path/to/java -jar \
# /path/to/hap-ibd.jar \
# gt=ibd/5panel.chr${i}.shapeit5.${ind}.vcf.gz \
# map=ibd/5panel.chr${i}.shapeit5.${ind}.map.map \
# out=ibd/5panel.chr${i}.shapeit5.${ind}.map.hapIBD
done
done
