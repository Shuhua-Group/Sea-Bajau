#!/bin/bash
# Process variant calls and summarize variant catalogs.
# Set input, reference and tool paths for your environment before execution.

cd /path/to/67West_51East || exit 1
bcftools filter -e 'abs(ILEN) >= 1' 67West_51East.chm13v2.dv161.chr1_22.vcf.gz |bgzip -@8 > 67West_51East.chm13v2.dv161.chr1_22.snp.vcf.gz
bcftools filter -i 'abs(ILEN) >= 1' 67West_51East.chm13v2.dv161.chr1_22.vcf.gz |bgzip -@8 > 67West_51East.chm13v2.dv161.chr1_22.indel.vcf.gz
cd plink || exit 1
plink --vcf ../67West_51East.chm13v2.dv161.chr1_22.snp.vcf.gz \
      --vcf-half-call missing \
      --make-bed \
      --out 67West_51East.chm13v2.dv161.chr1_22.snp

# Update this setting fam, Add sex metadata (Not yet included)

plink --bfile 67West_51East.chm13v2.dv161.chr1_22.snp --distance square --out 67West_51East.chm13v2.dv161.chr1_22.snp.similarity

# Use KING Estimate relatedness
# (, 0.0442] Unrelated
# (0.0442, 0.0884] 3rd Degree Relationship
# (0.0884, 0.1774] 2nd Degree Relationship
# (0.1774, 0.354] 1st Degree Relationship
# (0.354,) Duplicated
king -b 67West_51East.chm13v2.dv161.chr1_22.snp.bed --kinship --prefix 67West_51East.chm13v2.dv161.chr1_22.snp.king

####### 1st Degree Relationship #######
awk '{if ($8 > 0.1774 && $8 < 0.354) print $1,"\t",$3}' 67West_51East.chm13v2.dv161.chr1_22.snp.king.kin0|less
# JR008 KK024
# JR008 JR018
# KK024 KK030
# NB2340 NB2346
# NB2313 NB2314
# NB2336 NB2381
# KK037 KK054
# JR019 JR022
# JR006 JR00
# KK010 KK054
# KK014 KK037
# JR008 JR021
# KK004 KK021
# NB2336 NB2338
# KK010 KK037
# NB2314 NB2323
# NB2323 NB2324
# NB2313 NB2323
# JR018 JR021
# NB2322 NB2325
# KK014 KK021
# KK010 KK014
# NB2322 NB2328
# NB2336 NB2337
# NB2381 NB2382
# NB2315 NB2326
# KK014 KK054
####### 2nd Degree Relationship #######
awk '{if ($8 > 0.0884 && $8 < 0.1774) print $1,"\t",$3}' 67West_51East.chm13v2.dv161.chr1_22.snp.king.kin0
# JR002 JR006
# JR007 JR008
# JR007 KK003
# JR010 JR014
# JR010 JR015
# JR011 JR014
# JR016 KK003
# JR007 KK024
# JR008 KK030
# JR018 KK024
# JR021 KK024
# KK012 KK056
# KK021 KK037
# NB2313 NB2324
# NB2313 NB2326
# NB2314 NB2315
# NB2314 NB2318
# NB2314 NB2324
# NB2314 NB2326
# NB2320 NB2327
# NB2338 NB2381
####### 3rd Degree Relationship #######
awk '{if ($8 > 0.0442 && $8 < 0.0884) print $1,"\t",$3}' 67West_51East.chm13v2.dv161.chr1_22.snp.king.kin0|less
# JR002 JR008
# JR002 JR009
# JR002 JR017
# JR002 JR021
# JR006 JR013
# JR008 JR014
# JR009 JR013
# JR010 JR017
# JR010 JR018
# JR014 JR015
# JR014 JR016
# JR015 JR017
# JR015 JR021
# JR017 JR018
# JR017 JR020
# JR017 JR021
# JR017 KK009
# JR020 KK010
# JR021 KK009
# KK001 KK003
# JR005 KK055
# JR005 KK056
# JR007 KK030
# JR014 KK024
# JR021 KK030
# KK012 KK037
# KK012 KK054
# NB2313 NB2315
# NB2315 NB2320
# NB2316 NB2326
# NB2318 NB2323
# NB2320 NB2326
# NB2337 NB2338
# NB2336 NB2382
# NB2337 NB2381
# NB2339 NB2379
# NB2340 NB2341
# NB2341 NB2346

### 1.1 Standard filtering (Keep related individuals)
plink --bfile 67West_51East.chm13v2.dv161.chr1_22.snp \
--geno 0.05 --maf 0.05 \
--bp-space 50000 --make-bed \
--out 67West_51East.chm13v2.dv161.chr1_22.snp.geno005.maf005.ld50k
### 1.2 PCA
plink --bfile 67West_51East.chm13v2.dv161.chr1_22.snp.geno005.maf005.ld50k \
--pca --out 67West_51East.chm13v2.dv161.chr1_22.snp.geno005.maf005.ld50k.PCA_out

## 1.3 Output GT
plink --bfile 67West_51East.chm13v2.dv161.chr1_22.snp.geno005.maf005.ld50k --recode A --out 67West_51East.chm13v2.dv161.chr1_22.snp.geno005.maf005.ld50k.GT

## 1.4 Calculate the distance matrix
plink --bfile 67West_51East.chm13v2.dv161.chr1_22.snp.geno005.maf005.ld50k --distance-matrix --recode --out 67West_51East.chm13v2.dv161.chr1_22.snp.geno005.maf005.ld50k.Distance.out

## 1.5 Convert phylip
## Write another vcf
plink --bfile 67West_51East.chm13v2.dv161.chr1_22.snp.geno005.maf005.ld50k --recode vcf-iid --out 67West_51East.chm13v2.dv161.chr1_22.snp.geno005.maf005.ld50k
## Convert phy
python vcf2phy.py -i 67West_51East.chm13v2.dv161.chr1_22.snp.geno005.maf005.ld50k.vcf -o 67West_51East.chm13v2.dv161.chr1_22.snp.geno005.maf005.ld50k.phy
# Then Mega Build a phylogeny, Reference https://zhuanlan.zhihu.com/p/142570942

### 2.1 Exclude related individuals
cd remove_1stking || exit 1
plink --bfile ../67West_51East.chm13v2.dv161.chr1_22.snp \
      --remove 16samples_1stkinship_torm.txt \
      --make-bed \
      --out 67West_51East.rm16
plink --bfile 67West_51East.rm16 \
--geno 0.05 --maf 0.05 \
--bp-space 50000 --make-bed \
--out 67West_51East.rm16.geno005.maf005.ld50k

## 2.2 PCA
plink --bfile 67West_51East.rm16.geno005.maf005.ld50k \
--pca --out 67West_51East.rm16.geno005.maf005.ld50k.PCA_out

## 2.3 Output GT
plink --bfile 67West_51East.rm16.geno005.maf005.ld50k \
--recode A \
--out 67West_51East.rm16.geno005.maf005.ld50k.GT

## 2.4 Calculate the distance matrix
plink --bfile 67West_51East.rm16.geno005.maf005.ld50k  \
--distance-matrix --recode --out \
67West_51East.rm16.geno005.maf005.ld50k.Distance.out

## 2.5 Convert phylip
## Write another vcf
plink --bfile 67West_51East.rm16.geno005.maf005.ld50k \
--recode vcf-iid \
--out 67West_51East.rm16.geno005.maf005.ld50k
## Convert phy
python /path/to/vcf2phy.py vcf2phy.py 67West_51East.rm16.geno005.maf005.ld50k.vcf \
-o 67West_51East.rm16.geno005.maf005.ld50k.phy
# Then Mega Build a phylogeny, Reference https://zhuanlan.zhihu.com/p/142570942

## 3. Add reference populations to inspect genetic structure
#3.1 Merge first ho_1926 Inspect
wd="/path/to/pop_structure"
tmp_dir="${wd}/tmp_dir"
malay="/path/to/67West_51East.GRCh38.dv161.chr1-22.snp.vcf.gz"
ho_1926="/path/to/Ho.1926.merge.vcf.gz"
bcftools isec --threads 8 "${malay}" "${ho_1926}" -n=2 -p "${tmp_dir}"
bgzip "${tmp_dir}/0000.vcf" -c > "${tmp_dir}/0000.vcf.gz" && tabix -f "${tmp_dir}/0000.vcf.gz"
bgzip "${tmp_dir}/0001.vcf" -c > "${tmp_dir}/0001.vcf.gz" && tabix -f "${tmp_dir}/0001.vcf.gz"
bcftools merge "${tmp_dir}/0000.vcf.gz" "${tmp_dir}/0001.vcf.gz" -Oz -o "${wd}/malay_ho_1926.vcf.gz"

#3.1.1 Then plink, Run PCA
plink --vcf "${wd}/malay_ho_1926.vcf.gz" \
      --vcf-half-call missing \
      --make-bed \
      --const-fid 0 \
      --out malay_ho_1926

plink --bfile malay_ho_1926 \
--geno 0.05 --maf 0.05 \
--bp-space 50000 --make-bed \
--out malay_ho_1926.snp.geno005.maf005.ld50k

plink --bfile malay_ho_1926.snp.geno005.maf005.ld50k \
--pca --out malay_ho_1926.snp.geno005.maf005.ld50k.PCA_out

# Then continue with plot script
# /path/to/02.PCA.pop_stru.R
# /path/to/02.PCA.pop_stru2.R

# 3.2 Merge ho_1926 And PanMLS
wd="/path/to/pop_structure"
tmp_dir="${wd}/tmp_dir"
malay="/path/to/67West_51East.GRCh38.dv161.chr1-22.snp.vcf.gz"
ho_1926="/path/to/Ho.1926.merge.vcf.gz"
panmls="/path/to/PanMLS461.VQSR_PASS.SNP.biallelic.universal_mask_low_complexity_heng75.chr1-22.hg38.chr.sorted.vcf.gz"
bcftools isec --threads 8 "${malay}" "${ho_1926}" "${panmls}" -n=3 -p "${tmp_dir}"
bgzip "${tmp_dir}/0000.vcf" -c > "${tmp_dir}/0000.vcf.gz" && tabix -f "${tmp_dir}/0000.vcf.gz"
bgzip "${tmp_dir}/0001.vcf" -c > "${tmp_dir}/0001.vcf.gz" && tabix -f "${tmp_dir}/0001.vcf.gz"
bgzip "${tmp_dir}/0002.vcf" -c > "${tmp_dir}/0002.vcf.gz" && tabix -f "${tmp_dir}/0002.vcf.gz"
# ho_1926 And PanMLS Contains duplicate samples --force-samples
bcftools merge "${tmp_dir}/0000.vcf.gz" \
${tmp_dir}/0001.vcf.gz "${tmp_dir}/0002.vcf.gz" \
--force-samples \
-Oz -o "${wd}/malay_ho_1926_panmls.vcf.gz"

# 3.2.1 Then plink, Run PCA
plink --vcf "${wd}/malay_ho_1926_panmls.vcf.gz" \
      --vcf-half-call missing \
      --make-bed \
      --const-fid 0 \
      --out malay_ho_1926_panmls

plink --bfile malay_ho_1926_panmls \
--geno 0.05 --maf 0.05 \
--bp-space 50000 --make-bed \
--out malay_ho_1926_panmls.snp.geno005.maf005.ld50k

plink --bfile malay_ho_1926_panmls.snp.geno005.maf005.ld50k \
--pca --out malay_ho_1926_panmls.snp.geno005.maf005.ld50k.PCA_out

# 3.3 Merge PanMLS
wd="/path/to/pop_structure"
tmp_dir="${wd}/tmp_dir"
malay="/path/to/67West_51East.GRCh38.dv161.chr1-22.snp.vcf.gz"
panmls="/path/to/PanMLS461.VQSR_PASS.SNP.biallelic.universal_mask_low_complexity_heng75.chr1-22.hg38.chr.sorted.vcf.gz"
bcftools isec --threads 16 "${malay}" "${panmls}" -n=2 -p "${tmp_dir}"
bgzip "${tmp_dir}/0000.vcf" -c > "${tmp_dir}/0000.vcf.gz" && tabix -f "${tmp_dir}/0000.vcf.gz"
bgzip "${tmp_dir}/0001.vcf" -c > "${tmp_dir}/0001.vcf.gz" && tabix -f "${tmp_dir}/0001.vcf.gz"

bcftools merge "${tmp_dir}/0000.vcf.gz" \
${tmp_dir}/0001.vcf.gz \
-Oz -o "${wd}/malay_panmls.vcf.gz"
tabix -f "${wd}/malay_panmls.vcf.gz"

# 3.3.1 Then plink, Run PCA
plink --vcf "${wd}/malay_panmls.vcf.gz" \
      --vcf-half-call missing \
      --make-bed \
      --const-fid 0 \
      --out malay_panmls

plink --bfile malay_panmls \
--geno 0.05 --maf 0.05 \
--bp-space 50000 --make-bed \
--out malay_panmls.snp.geno005.maf005.ld50k

plink --bfile malay_panmls.snp.geno005.maf005.ld50k \
--pca --out malay_panmls.snp.geno005.maf005.ld50k.PCA_out

# 4. Refine ho_1926 scale
# 4.1 First scale: Remove Africa & WestEurasia
cat /path/to/118Allsamples.list >> ../reference/ho_1926/ho_1926.kp.info.rm_AFR_WestEura.list
bcftools view -S ../reference/ho_1926/ho_1926.kp.info.rm_AFR_WestEura.list malay_ho_1926.vcf.gz -Oz -o malay_ho_1926.rm_AFR_WestEura.vcf.gz
plink --vcf malay_ho_1926.rm_AFR_WestEura.vcf.gz \
      --vcf-half-call missing \
      --make-bed \
      --const-fid 0 \
      --out malay_ho_1926.rm_AFR_WestEura

plink --bfile malay_ho_1926.rm_AFR_WestEura \
--geno 0.05 --maf 0.05 \
--bp-space 50000 --make-bed \
--out malay_ho_1926.rm_AFR_WestEura.snp.geno005.maf005.ld50k

plink --bfile malay_ho_1926.rm_AFR_WestEura.snp.geno005.maf005.ld50k \
--pca --out malay_ho_1926.rm_AFR_WestEura.snp.geno005.maf005.ld50k.PCA_out


# 4.2 Second scale: Then remove America,CentralAsiaSiberia
cat /path/to/118Allsamples.list >> ../reference/ho_1926/ho_1926.kp.info.rm_AFR_WestEura_AMR_Siberia.list
bcftools view -S ../reference/ho_1926/ho_1926.kp.info.rm_AFR_WestEura_AMR_Siberia.list malay_ho_1926.vcf.gz -Oz -o malay_ho_1926.rm_AFR_WestEura_AMR_Siberia.vcf.gz
plink --vcf malay_ho_1926.rm_AFR_WestEura_AMR_Siberia.vcf.gz \
      --vcf-half-call missing \
      --make-bed \
      --const-fid 0 \
      --out malay_ho_1926.rm_AFR_WestEura_AMR_Siberia

plink --bfile malay_ho_1926.rm_AFR_WestEura_AMR_Siberia \
--geno 0.05 --maf 0.05 \
--bp-space 50000 --make-bed \
--out malay_ho_1926.rm_AFR_WestEura_AMR_Siberia.snp.geno005.maf005.ld50k

plink --bfile malay_ho_1926.rm_AFR_WestEura_AMR_Siberia.snp.geno005.maf005.ld50k \
--pca --out malay_ho_1926.rm_AFR_WestEura_AMR_Siberia.snp.geno005.maf005.ld50k.PCA_out


# Restart this step, Merge Singapore data with PanMLS selected data from
wd="/path/to/pop_structure"
for i in {1..22};do
tmp_dir="${wd}/tmp_dir/chr${i}"
mkdir -p "${tmp_dir}"
malay="/path/to/67West_51East.GRCh38.dv161.chr${i}.snp.vcf.gz"
ho_1926="/path/to/Ho.1926.merge.chr${i}.vcf.gz"
panmls="/path/to/PanMLS461.VQSR_PASS.SNP.biallelic.universal_mask_low_complexity_heng75.chr1-22.hg38.chr${i}.sorted.vcf.gz"
sg10k="/path/to/SG10K.chr${i}.SNV.b38.correct_ref.biallelic.rm_multi.chr.vcf.gz"
bcftools isec --threads 8 "${malay}" "${ho_1926}" "${panmls}" "${sg10k}" -n=4 -p "${tmp_dir}"
bgzip "${tmp_dir}/0000.vcf" -c > "${tmp_dir}/0000.vcf.gz" && tabix -f "${tmp_dir}/0000.vcf.gz"
bgzip "${tmp_dir}/0001.vcf" -c > "${tmp_dir}/0001.vcf.gz" && tabix -f "${tmp_dir}/0001.vcf.gz"
bgzip "${tmp_dir}/0002.vcf" -c > "${tmp_dir}/0002.vcf.gz" && tabix -f "${tmp_dir}/0002.vcf.gz"
bgzip "${tmp_dir}/0003.vcf" -c > "${tmp_dir}/0003.vcf.gz" && tabix -f "${tmp_dir}/0003.vcf.gz"
bcftools merge "${tmp_dir}/0000.vcf.gz" \
${tmp_dir}/0001.vcf.gz "${tmp_dir}/0002.vcf.gz" \
${tmp_dir}/0003.vcf.gz \
--force-samples \
-Oz -o "${wd}/malay_ho1926_panmls_SG10K.GRCh38.chr${i}.vcf.gz"
tabix -f "${wd}/malay_ho1926_panmls_SG10K.GRCh38.chr${i}.vcf.gz"
done &

for i in {1..22};do (bcftools view -S /path/to/merge.sample.info.list  \
${wd}/malay_ho1926_panmls_SG10K.GRCh38.chr${i}.vcf.gz \
-Oz -o "${wd}/malay_ho1926_panmls_SG10K.GRCh38.chr${i}.select1792.vcf.gz" &);done
wait

bcftools concat --threads 16 ${wd}/malay_ho1926_panmls_SG10K.GRCh38.chr*.select1792.vcf.gz \
-Oz -o "${wd}/malay_ho1926_panmls_SG10K.GRCh38.select1792.vcf.gz"

tabix -f "${wd}/malay_ho1926_panmls_SG10K.GRCh38.select1792.vcf.gz"
cd "${wd}" || exit 1
plink --vcf malay_ho1926_panmls_SG10K.GRCh38.select1792.vcf.gz \
      --vcf-half-call missing \
      --make-bed \
      --const-fid 0 \
      --out malay_ho1926_panmls_SG10K.GRCh38.select1792

plink --bfile malay_ho1926_panmls_SG10K.GRCh38.select1792 \
--geno 0.05 --maf 0.05 \
--bp-space 50000 --make-bed \
--out malay_ho1926_panmls_SG10K.GRCh38.select1792.snp.geno005.maf005.ld50k

plink --bfile malay_ho1926_panmls_SG10K.GRCh38.select1792.snp.geno005.maf005.ld50k \
--pca --out malay_ho1926_panmls_SG10K.GRCh38.select1792.snp.geno005.maf005.ld50k.PCA_out

# Then process short-read and long-read data Malaysia Inspect separately
bcftools view -S /path/to/merge.Malaysia_sample.info.list \
malay_ho1926_panmls_SG10K.GRCh38.select1792.vcf.gz \
-Oz -o malay_ho1926_panmls_SG10K.GRCh38.select1792.Malaysia.vcf.gz
tabix -f malay_ho1926_panmls_SG10K.GRCh38.select1792.Malaysia.vcf.gz

plink --vcf malay_ho1926_panmls_SG10K.GRCh38.select1792.Malaysia.vcf.gz \
      --vcf-half-call missing \
      --make-bed \
      --const-fid 0 \
      --out malay_ho1926_panmls_SG10K.GRCh38.select1792.Malaysia

plink --bfile malay_ho1926_panmls_SG10K.GRCh38.select1792.Malaysia \
--geno 0.05 --maf 0.05 \
--bp-space 50000 --make-bed \
--out malay_ho1926_panmls_SG10K.GRCh38.select1792.Malaysia.snp.geno005.maf005.ld50k

plink --bfile malay_ho1926_panmls_SG10K.GRCh38.select1792.Malaysia.snp.geno005.maf005.ld50k \
--pca --out malay_ho1926_panmls_SG10K.GRCh38.select1792.Malaysia.snp.geno005.maf005.ld50k.PCA_out

