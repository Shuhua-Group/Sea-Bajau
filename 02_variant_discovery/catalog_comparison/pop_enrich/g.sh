#!/bin/bash

# Process variant calls and summarize variant catalogs.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/pop_enrich"
cd "$wd" || exit 1

for chr in  {1..22}  X Y
# for chr in {21..22}
do
# bcftools view -r chr${chr} ../APC_uniq.norm2.chm13v2.allchr.rmGRCh38.vcf.gz \
# -Oz -o temp.chr${chr}.vcf.gz &
# tabix -p vcf temp.chr${chr}.vcf.gz &

python fast_calculate_population_AF.considersex.py \
--vcf "temp.chr${chr}.vcf.gz" \
--meta ../130APC.v3.info.txt \
--out "chr${chr}.AF.considersex.miss0.2.tsv" \
--missing 0.2 &
done

cat chr*.AF.considersex.tsv > APC_uniq.allchr.AF.considersex.tsv
cat chr*.AF.considersex.miss0.2.tsv > APC_uniq.allchr.AF.considersex.miss0.2.tsv

python calculate_individual_burden.py \
    --vcf ../APC_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.vcf.gz \
    --meta ../130APC.v3.info.txt \
    --vclass ../APC_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.tsv \
    --out APC_uniq.norm2.chm13v2.allchr.individual_burden.tsv
