#!/bin/bash

# Process variant calls and summarize variant catalogs.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/pop_enrich"
cd "$wd" || exit 1

for chr in  {1..22}  X Y
# for chr in {21..22}
do
# bcftools view -r chr${chr} ../CPC2_uniq.norm2.chm13v2.allchr.rmGRCh38.vcf.gz \
# -Oz -o temp.CPC2.chr${chr}.vcf.gz &
# tabix -p vcf temp.CPC2.chr${chr}.vcf.gz &

python fast_calculate_population_AF.py \
--vcf "temp.CPC2.chr${chr}.vcf.gz" \
--meta 231CPC_info.lan.txt \
--out "CPC2.chr${chr}.AF.miss0.2.tsv" \
--missing 0.2 &
done

cat CPC2.chr*.AF.miss0.2.tsv > CPC2_uniq.allchr.AF.miss0.2.tsv

python calculate_individual_burden.py \
    --vcf ../CPC2_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.vcf.gz \
    --meta 231CPC_info.lan.txt \
    --vclass ../CPC2_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.tsv \
    --out CPC2_uniq.norm2.chm13v2.allchr.individual_burden.tsv
