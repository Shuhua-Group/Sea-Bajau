#!/bin/bash

# Process variant calls and summarize variant catalogs.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/pop_enrich"
cd "$wd" || exit 1

for chr in  {1..22}  X Y
do
# bcftools view -r chr${chr} ../HPRC2_uniq.norm2.chm13v2.allchr.rmGRCh38.vcf.gz \
# -Oz -o temp.HPRC2.chr${chr}.vcf.gz &
# tabix -p vcf temp.HPRC2.chr${chr}.vcf.gz &

python fast_calculate_population_AF.py \
--vcf "temp.HPRC2.chr${chr}.vcf.gz" \
--meta HPRC2.info.txt \
--out "HPRC2.chr${chr}.AF.miss0.2.tsv" \
--missing 0.2 &
done


cat HPRC2.chr*.AF.miss0.2.tsv > HPRC2_uniq.allchr.AF.miss0.2.tsv

python calculate_individual_burden.py \
    --vcf ../HPRC2_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.vcf.gz \
    --meta HPRC2.info.txt \
    --vclass ../HPRC2_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.tsv \
    --out HPRC2_uniq.norm2.chm13v2.allchr.individual_burden.tsv
