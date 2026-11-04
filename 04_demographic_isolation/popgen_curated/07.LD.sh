#!/bin/bash

# Analyze Bajau under selection in SV
# Because SV Cannot readily calculate iHS, XPEHH, Fst and related selection statistics
# First use SNP Calculate selection statistics first, Then inspect SV And SNP linkage disequilibrium
# For variants linked to selected SNP in strong linkage disequilibrium SV Annotate, Found casual Of SV

# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/LD_SV"
# 1. Based on ID.txt Use vcf.gz within SV Select
bcftools view -i 'ID=@filtered_high_r2_results_del.ID.txt' 102Malay_GRCh38_sawfish.jc.sv.del.ld.vcf.gz \
-Oz -o 102Malay_GRCh38_sawfish.jc.sv.del.top001_CMS.ld.vcf.gz

bcftools view -i 'ID=@filtered_high_r2_results_ins.ID.txt' 102Malay_GRCh38_sawfish.jc.sv.ins.ld.vcf.gz \
-Oz -o 102Malay_GRCh38_sawfish.jc.sv.ins.top001_CMS.ld.vcf.gz

# 2. CADD annotated results
bcftools query -f '%CHROM\t%POS\t%SVLEN\n' 102Malay_GRCh38_sawfish.jc.sv.del.top001_CMS.ld.vcf.gz \
> "${wd}/102Malay_GRCh38_sawfish.jc.sv.del.top001_CMS.ld.temp"
awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$2-$3}' "${wd}/102Malay_GRCh38_sawfish.jc.sv.del.top001_CMS.ld.temp" \
> "${wd}/102Malay_GRCh38_sawfish.jc.sv.del.top001_CMS.ld.bed"

bcftools query -f '%CHROM\t%POS\t%SVLEN\n' 102Malay_GRCh38_sawfish.jc.sv.ins.top001_CMS.ld.vcf.gz \
> "${wd}/102Malay_GRCh38_sawfish.jc.sv.ins.top001_CMS.ld.temp"
awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$2+$3}' "${wd}/102Malay_GRCh38_sawfish.jc.sv.ins.top001_CMS.ld.temp" \
> "${wd}/102Malay_GRCh38_sawfish.jc.sv.ins.top001_CMS.ld.bed"
rm ${wd}/*.temp

sed -i "s/chr//g" "${wd}/102Malay_GRCh38_sawfish.jc.sv.del.top001_CMS.ld.bed"
sed -i "s/chr//g" "${wd}/102Malay_GRCh38_sawfish.jc.sv.ins.top001_CMS.ld.bed"


awk 'FNR==NR { key[$1 FS $2 FS $3]; next } ($1 FS $2 FS $3) in key' \
${wd}/102Malay_GRCh38_sawfish.jc.sv.del.top001_CMS.ld.bed \
BajaudelchrA.prescore.bed > 102Malay_GRCh38_sawfish.jc.sv.del.top001_CMS.ld.CADD.bed


awk 'FNR==NR { key[$1 FS $2 FS $3]; next } ($1 FS $2 FS $3) in key' \
${wd}/102Malay_GRCh38_sawfish.jc.sv.ins.top001_CMS.ld.bed \
BajauinschrA.prescore.bed > 102Malay_GRCh38_sawfish.jc.sv.ins.top001_CMS.ld.CADD.bed
