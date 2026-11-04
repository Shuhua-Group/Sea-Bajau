#!/bin/bash

# Keep only GujaratiA
# Process variant calls and summarize variant catalogs.
# Set input, reference and tool paths for your environment before execution.

cd /path/to/admixture || exit 1

# Revise sample selection
# Run a first version level2 version
plink --vcf ../plink/7panel.chrA.level2.vcf.gz \
--vcf-half-call missing \
--geno 0.01 --maf 0.05 --mind 0.1 \
--indep-pairwise 200kb 50 0.4 \
--make-bed --out 7panel.chrA.level2.admixture

for K in {1..10}; do
( /path/to/admixture \
--cv 7panel.chrA.level2.admixture.bed "$K" | \
tee "7panel.chrA.level2.admixture_log${K}.out" &); done

grep "CV error" 7panel.chrA.level2.admixture_log*out

# level2 Add YRI And French Run another version
bcftools view ../7panel.chrA.vcf.gz \
-S ../plink/sample_info/102Malay.ho_1926.PanMLS.SG10K.EGDP.Pacific.2018Cell.for_admixture.list \
-Oz -o 7panel.chrA.level2_AFR_EUR.vcf.gz

plink --vcf 7panel.chrA.level2_AFR_EUR.vcf.gz \
--vcf-half-call missing \
--geno 0.01 --maf 0.05 --mind 0.1 \
--indep-pairwise 200kb 50 0.4 \
--make-bed --out 7panel.chrA.level2_AFR_EUR.admixture

for K in {16..20}; do
(/path/to/admixture \
--cv 7panel.chrA.level2_AFR_EUR.admixture.bed "$K" | \
tee "7panel.chrA.level2_AFR_EUR.admixture_log${K}.out" &); done

grep "CV error" 7panel.chrA.level2_AFR_EUR.admixture_log*out

