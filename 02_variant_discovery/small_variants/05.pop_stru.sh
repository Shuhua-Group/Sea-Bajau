#!/bin/bash

# Merge the updated panel [ All are GRCh38 Coordinates ]
# Input data
# 1. 102 Malaysia Long-read data
# 2. ho_1926 Short-read data
# 3. PanMLS Short-read data
# 4. EGDP Short-read data
# 5. 2018cell Low-coverage short-read data
# 6. Pacific Array data
# 7. SG10K Singapore data

# Process variant calls and summarize variant catalogs.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/re_merge"
for i in {1..22};do
# Note:, All samples must index
set1="/path/to/102Malay.GRCh38.dv161.chr${i}.norm.QC.SNP.shapeited.reID2.vcf.gz"
set2="/path/to/Ho.1926.merge.chr${i}.vcf.gz"
set3="/path/to/PanMLS461.VQSR_PASS.SNP.biallelic.universal_mask_low_complexity_heng75.chr1-22.hg38.chr${i}.sorted.vcf.gz"
set4="/path/to/EGDP_PaganiEtAl2016_Release_Build37_Chr${i}.94sample.hg38.sort.vcf.gz"
set5="/path/to/90combined.raw_variants.snps.PASS.hg38.chr${i}.vcf.gz"
set6="/path/to/Pacific.uni.hg38.chr${i}.sorted.selected.vcf.gz"
set7="/path/to/SG10K.chr${i}.SNV.b38.correct_ref.biallelic.rm_multi.Peranakan_79.Indian_50.vcf.gz"
tmp_dir="${wd}/tmp_dir/chr${i}"
mkdir -p "${tmp_dir}"
/path/to/bcftools isec --threads 4 \
${set1} "${set2}" "${set3}" "${set4}" \
${set5} "${set6}" "${set7}" \
-n=7 -p "${tmp_dir}"
bgzip "${tmp_dir}/0000.vcf" -c > "${tmp_dir}/0000.vcf.gz" && tabix -f "${tmp_dir}/0000.vcf.gz"
bgzip "${tmp_dir}/0001.vcf" -c > "${tmp_dir}/0001.vcf.gz" && tabix -f "${tmp_dir}/0001.vcf.gz"
bgzip "${tmp_dir}/0002.vcf" -c > "${tmp_dir}/0002.vcf.gz" && tabix -f "${tmp_dir}/0002.vcf.gz"
bgzip "${tmp_dir}/0003.vcf" -c > "${tmp_dir}/0003.vcf.gz" && tabix -f "${tmp_dir}/0003.vcf.gz"
bgzip "${tmp_dir}/0004.vcf" -c > "${tmp_dir}/0004.vcf.gz" && tabix -f "${tmp_dir}/0004.vcf.gz"
bgzip "${tmp_dir}/0005.vcf" -c > "${tmp_dir}/0005.vcf.gz" && tabix -f "${tmp_dir}/0005.vcf.gz"
bgzip "${tmp_dir}/0006.vcf" -c > "${tmp_dir}/0006.vcf.gz" && tabix -f "${tmp_dir}/0006.vcf.gz"
/path/to/bcftools merge "${tmp_dir}/0000.vcf.gz" "${tmp_dir}/0001.vcf.gz" \
${tmp_dir}/0002.vcf.gz "${tmp_dir}/0003.vcf.gz" \
${tmp_dir}/0004.vcf.gz "${tmp_dir}/0005.vcf.gz" \
${tmp_dir}/0006.vcf.gz \
--force-samples \
-Oz -o "${wd}/7panel.chr${i}.vcf.gz"
done &

# Merge all chromosomes
bcftools concat 7panel.chr*.vcf.gz \
-Oz -o 7panel.chrA.vcf.gz
tabix -f 7panel.chrA.vcf.gz

# Not all sites are used
# Repeat sample selection

# 1. 102 Malaysia Long-read data [ Keep all records ]
# 2. ho_1926 Short-read data [ Keep all, Keep East and South Asian populations, Keep only CHS]
# 3. PanMLS Short-read data [../../reference/PanMLS461/PanMLS461.select.list]
# 4. EGDP Short-read data [ A single level ../../reference/EGDP/sample_info.txt]
# 5. 2018cell Low-coverage short-read data [ Keep all records ]
# 6. Pacific Array data [ Keep all records initially ]
# 7. SG10K Singapore data [ A single level Peranakan_79.Indian_50]

#

# Define several level, Use successively smaller sample sets
# 1. Keep unfiltered data, All samples
cd plink || exit 1

plink --vcf ../7panel.chrA.vcf.gz \
      --vcf-half-call missing \
      --make-bed \
      --const-fid 0 \
      --out 7panel.chrA
# 79767 variants and 4110 people pass filters and QC.

plink --bfile 7panel.chrA --missing --out 7panel.chrA.qc_results

king -b 7panel.chrA.bed --kinship --prefix 7panel.chrA.king
awk '{if ($9 == 0.5) print $3}' 7panel.chrA.king.kin > sample_duplicated.list
echo "IID FID" > sample_duplicated.txt
awk '{print "0", $0}' sample_duplicated.list >> sample_duplicated.txt
# sed -i "s/\t/ /g" sample_duplicated.txt

plink --bfile 7panel.chrA \
      --remove sample_duplicated.txt \
      --make-bed \
      --out 7panel.chrA.rm_duplicated
# 3985 people remaining.


plink --bfile 7panel.chrA \
--geno 0.01 --maf 0.05 --mind 0.1 \
--indep-pairwise 200kb 50 0.4 --make-bed \
--out 7panel.chrA.geno001.maf005.ld
# 77294 variants and 4065 people pass filters and QC.
# 200kb 50 0.4: 3390 of 77294 variants removed

plink --bfile 7panel.chrA.rm_duplicated \
--geno 0.01 --maf 0.05 --mind 0.1 \
--indep-pairwise 200kb 50 0.4 --make-bed \
--out 7panel.chrA.rm_duplicated.geno001.maf005.ld
# 77276 variants and 3940 people pass filters and QC.
# 200kb 50 0.4: 3404 of 77276 variants removed.


plink --bfile 7panel.chrA.geno001.maf005.ld \
--pca --out 7panel.chrA.geno001.maf005.ld.PCA_out
flashpca --bfile 7panel.chrA.geno001.maf005.ld --suffix 7panel.chrA.geno001.maf005.ld

plink --bfile 7panel.chrA.rm_duplicated.geno001.maf005.ld \
--pca --out 7panel.chrA.rm_duplicated.geno001.maf005.ld.PCA_out

plink --bfile 7panel.chrA.rm_duplicated.geno001.maf005.ld \
 --recode vcf-iid --out 7panel.chrA.rm_duplicated.geno001.maf005.ld
# Number of samples in VCF: 3940
python vcf2phy.py -i 7panel.chrA.rm_duplicated.geno001.maf005.ld.vcf \
-o 7panel.chrA.rm_duplicated.geno001.maf005.ld


# 2. level1 Use ho1926 Exclude distant populations 7panel.chrA.level1.vcf.gz
# 7panel.chrA.level1.vcf.gz 2040 samples
bcftools view ../7panel.chrA.vcf.gz \
-S sample_info/102Malay.ho_1926.PanMLS.SG10K.EGDP.Pacific.2018Cell.list \
-Oz -o 7panel.chrA.level1.vcf.gz
plink --vcf 7panel.chrA.level1.vcf.gz \
      --vcf-half-call missing \
      --make-bed \
      --const-fid 0 \
      --out 7panel.chrA.level1
# 79767 variants and 2040 people pass filters and QC.


plink --vcf 7panel.chrA.level1.vcf.gz \
      --vcf-half-call missing \
      --make-bed \
      --remove sample_duplicated.txt \
      --const-fid 0 \
      --out 7panel.chrA.rm_duplicated.level1
# 79767 variants and 2020 people pass filters and QC.



plink --bfile 7panel.chrA.level1 \
--geno 0.01 --maf 0.05 --mind 0.1 \
--indep-pairwise 200kb 50 0.4 --make-bed \
--out 7panel.chrA.level1.geno001.maf005.ld
# 74278 variants and 2037 people pass filters and QC.
# 200kb 50 0.4: 3728 of 74278 variants removed.

plink --bfile 7panel.chrA.rm_duplicated.level1 \
--geno 0.01 --maf 0.05 --mind 0.1 \
--indep-pairwise 200kb 50 0.4 --make-bed \
--out 7panel.chrA.rm_duplicated.level1.geno001.maf005.ld
# 74302 variants and 2017 people pass filters and QC.
# 200kb 50 0.4: 3686 of 74302 variants removed.


plink --bfile 7panel.chrA.level1.geno001.maf005.ld \
--pca --out 7panel.chrA.level1.geno001.maf005.ld.PCA_out
flashpca --bfile 7panel.chrA.level1.geno001.maf005.ld --suffix 7panel.chrA.level1.geno001.maf005.ld

plink --bfile 7panel.chrA.rm_duplicated.level1.geno001.maf005.ld \
--pca --out 7panel.chrA.rm_duplicated.level1.geno001.maf005.ld.PCA_out
flashpca --bfile 7panel.chrA.rm_duplicated.level1.geno001.maf005.ld --suffix 7panel.chrA.rm_duplicated.level1.geno001.maf005.ld

# $ which flashpca
# /path/to/flashpca

bcftools view ../7panel.chrA.vcf.gz \
-S sample_info/Bajau_Malay_Africa.for_tree.list \
-Oz -o 7panel.chrA.rm_duplicated.for_tree.vcf.gz
plink --vcf 7panel.chrA.rm_duplicated.for_tree.vcf.gz \
      --vcf-half-call missing \
      --make-bed \
      --const-fid 0 \
      --geno 0.01 --maf 0.05 --mind 0.1 \
      --indep-pairwise 200kb 50 0.4 \
      --out 7panel.chrA.rm_duplicated.for_tree.geno001.maf005.ld

plink --bfile 7panel.chrA.rm_duplicated.for_tree.geno001.maf005.ld \
 --recode vcf-iid --out 7panel.chrA.rm_duplicated.for_tree.geno001.maf005.ld
# Number of samples in VCF: 3940
python vcf2phy.py -i 7panel.chrA.rm_duplicated.for_tree.geno001.maf005.ld.vcf \
-o 7panel.chrA.rm_duplicated.for_tree.geno001.maf005.ld

# 3. level2 Keep only ho1926 selected southern Chinese minority populations +CHS(HGDP Of Han)
# Additionally, Here use SG10K Include only Peranakan_Chinese
# 7panel.chrA.level2.vcf.gz 599 samples
bcftools view ../7panel.chrA.vcf.gz \
-S sample_info/102Malay.ho_1926.PanMLS.SG10K.EGDP.Pacific.2018Cell.EAS.list \
-Oz -o 7panel.chrA.level2.vcf.gz

plink --vcf 7panel.chrA.level2.vcf.gz \
      --vcf-half-call missing \
      --make-bed \
      --const-fid 0 \
      --out 7panel.chrA.level2
# 79767 variants and 599 people pass filters and QC.


plink --vcf 7panel.chrA.level2.vcf.gz \
      --vcf-half-call missing \
      --remove sample_duplicated.txt \
      --make-bed \
      --const-fid 0 \
      --out 7panel.chrA.rm_duplicated.level2
# 79767 variants and 589 people pass filters and QC.

plink --bfile 7panel.chrA.level2 \
--geno 0.01 --maf 0.05 --mind 0.1 \
--indep-pairwise 200kb 50 0.4 --make-bed \
--out 7panel.chrA.level2.geno001.maf005.ld
# 62453 variants and 599 people pass filters and QC.
# 200kb 50 0.4: 2591 of 62453 variants removed.

plink --bfile 7panel.chrA.rm_duplicated.level2 \
--geno 0.01 --maf 0.05 --mind 0.1 \
--indep-pairwise 200kb 50 0.4 --make-bed \
--out 7panel.chrA.rm_duplicated.level2.geno001.maf005.ld
# 62474 variants and 589 people pass filters and QC.
# 200kb 50 0.4: 2576 of 62474 variants removed.


plink --bfile 7panel.chrA.level2.geno001.maf005.ld \
--pca --out 7panel.chrA.level2.geno001.maf005.ld.PCA_out
flashpca --bfile 7panel.chrA.level2.geno001.maf005.ld --suffix 7panel.chrA.level2.geno001.maf005.ld

plink --bfile 7panel.chrA.rm_duplicated.level2.geno001.maf005.ld \
--pca --out 7panel.chrA.rm_duplicated.level2.geno001.maf005.ld.PCA_out
flashpca --bfile 7panel.chrA.rm_duplicated.level2.geno001.maf005.ld --suffix 7panel.chrA.rm_duplicated.level2.geno001.maf005.ld

# 4. level3 Keep only Bajau,2018Cell, Philippines
bcftools view ../7panel.chrA.vcf.gz \
-S sample_info/102Malay.ho_1926.PanMLS.SG10K.EGDP.Pacific.2018Cell.ESA.list \
-Oz -o 7panel.chrA.level3.vcf.gz

plink --vcf 7panel.chrA.level3.vcf.gz \
      --vcf-half-call missing \
      --make-bed \
      --const-fid 0 \
      --out 7panel.chrA.level3
# 79767 variants and 130 people pass filters and QC.

plink --vcf 7panel.chrA.level3.vcf.gz \
      --vcf-half-call missing \
      --remove sample_duplicated.txt \
      --make-bed \
      --const-fid 0 \
      --out 7panel.chrA.rm_duplicated.level3
# 79767 variants and 130 people pass filters and QC.
# Identical, Skip the remaining steps for this branch

plink --bfile 7panel.chrA.level3 \
--geno 0.01 --maf 0.05 --mind 0.1 \
--indep-pairwise 200kb 50 0.4 --make-bed \
--out 7panel.chrA.level3.geno001.maf005.ld
# 29938 variants and 130 people pass filters and QC.
# 200kb 50 0.4: 394 of 29938 variants removed.
plink --bfile 7panel.chrA.level3.geno001.maf005.ld \
--pca --out 7panel.chrA.level3.geno001.maf005.ld.PCA_out
flashpca --bfile 7panel.chrA.level3.geno001.maf005.ld --suffix 7panel.chrA.level3.geno001.maf005.ld
