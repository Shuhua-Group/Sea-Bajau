#!/bin/bash

# Process variant calls and summarize variant catalogs.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/EGDP_PaganiEtAl2016_Release_Build37_PLINK"

cd /path/to/EGDP || exit 1
for i in {1..22};
do
( bcftools view -S sample.list \
$wd/EGDP_PaganiEtAl2016_Release_Build37_Chr${i}.vcf.gz \
-Oz -o "./EGDP_PaganiEtAl2016_Release_Build37_chr${i}.94sample.vcf.gz" &)
done



for i in {6..7};
do
CrossMap vcf /path/to/hg19ToHg38.over.chain.gz \
EGDP_PaganiEtAl2016_Release_Build37_Chr${i}.94sample.vcf.gz \
/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
EGDP_PaganiEtAl2016_Release_Build37_Chr${i}.94sample.hg38.vcf
bgzip "EGDP_PaganiEtAl2016_Release_Build37_Chr${i}.94sample.hg38.vcf"
bcftools sort "EGDP_PaganiEtAl2016_Release_Build37_Chr${i}.94sample.hg38.vcf.gz" \
-Oz -o "EGDP_PaganiEtAl2016_Release_Build37_Chr${i}.94sample.hg38.sorted.vcf.gz"
# Add chromosome prefixes
bcftools annotate \
  --rename-chrs /path/to/chr_name_mapping.txt \
  --threads 16 \
  -Oz -o "EGDP_PaganiEtAl2016_Release_Build37_Chr${i}.94sample.hg38.sort.vcf.gz" \
  "EGDP_PaganiEtAl2016_Release_Build37_Chr${i}.94sample.hg38.sorted.vcf.gz"
tabix -f "EGDP_PaganiEtAl2016_Release_Build37_Chr${i}.94sample.hg38.sort.vcf.gz"
done &

# Final processed data are stored in:
# /path/to/EGDP_PaganiEtAl2016_Release_Build37_Chr${i}.94sample.hg38.sort.vcf.gz

# Previous processing stopped here
# Merge with the previously combined dataset
wd="/path/to/pop_structure"
for i in {1..22};do
tmp_dir="${wd}/tmp_dir/chr${i}"
mkdir -p "$tmp_dir"
malay_ho1926_SG10K_PanMLS_2018Cell="${wd}/malay_ho1926_SG10K_PanMLS_BajauSaluan.chr${i}.vcf.gz"
EGDP="${wd}/../reference/EGDP/EGDP_PaganiEtAl2016_Release_Build37_Chr${i}.94sample.hg38.sort.vcf.gz"
bcftools isec --threads 8 "${malay_ho1926_SG10K_PanMLS_2018Cell}" "${EGDP}" -n=2 -p "${tmp_dir}"
bgzip "${tmp_dir}/0000.vcf" -c > "${tmp_dir}/0000.vcf.gz" && tabix -f "${tmp_dir}/0000.vcf.gz"
bgzip "${tmp_dir}/0001.vcf" -c > "${tmp_dir}/0001.vcf.gz" && tabix -f "${tmp_dir}/0001.vcf.gz"
done &

for i in {1..22};do
tmp_dir="${wd}/tmp_dir/chr${i}"
(bcftools merge "${tmp_dir}/0000.vcf.gz" "${tmp_dir}/0001.vcf.gz" -Oz -o "${wd}/malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr${i}.vcf.gz" && \
tabix -f "${wd}/malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr${i}.vcf.gz" &)
done

bcftools concat malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr*.vcf.gz \
--threads 32 \
-Oz -o malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.vcf.gz

# Use Malaysia And 2018Cell And EGDP Select, Add selected Southeast Asian populations, Add Oceanian minority populations
bcftools view malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.vcf.gz --threads 48 \
-S /path/to/merge.Malaysia_sample.90_2018Cell.EGDP.EastSouth.info.list \
-Oz -o malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA.vcf.gz

plink --vcf malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA.vcf.gz \
--vcf-half-call missing \
--make-bed \
--const-fid 0 \
--out malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA

plink --bfile malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA \
--geno 0.05 --maf 0.05 \
--bp-space 50000 --make-bed \
--out malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA.snp.geno005.maf005.ld50k

plink --bfile malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA.snp.geno005.maf005.ld50k \
--pca --out malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA.snp.geno005.maf005.ld50k.PCA_out

# Inspect the results, Use Bougainville NAN_Melanesian Sahul Exclude and inspect again
bcftools view malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA.vcf.gz --threads 32 \
-S /path/to/merge.Malaysia_sample.90_2018Cell.EGDP.EastSouth.info.v2.list \
-Oz -o malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA.v2.vcf.gz

plink --vcf malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA.v2.vcf.gz \
--vcf-half-call missing \
--make-bed \
--const-fid 0 \
--out malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA.v2

plink --bfile malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA.v2 \
--geno 0.05 --maf 0.05 \
--bp-space 50000 --make-bed \
--out malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA.v2.snp.geno005.maf005.ld50k

plink --bfile malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA.v2.snp.geno005.maf005.ld50k \
--pca --out malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA.v2.snp.geno005.maf005.ld50k.PCA_out
