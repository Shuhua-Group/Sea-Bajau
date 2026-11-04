#!/bin/bash

# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

conda activate shapeit5_env

for i in {1..22};do
SHAPEIT5_phase_common --input "../../malay_ho1926_SG10K_PanMLS_BajauSaluan.chr${i}.vcf.gz" \
--region "chr${i}" \
--map "/path/to/chr${i}.genetic_map_hg38.txt2" \
--output "malay_ho1926_SG10K_PanMLS_BajauSaluan.chr${i}.shapeited.bcf" \
--thread 32
done

for i in {1..22};do
SHAPEIT5_phase_common --input 69Malay_20Negrito.chr${i}.QC.vcf.gz \
--region chr${i} \
--map /path/to/chr${i}.genetic_map_hg38.txt2 \
--output 69Malay_20Negrito.chr${i}.QC.shapeited.bcf \
--thread 32
bcftools view "69Malay_20Negrito.chr${i}.QC.shapeited.bcf" -Oz -o "69Malay_20Negrito.chr${i}.QC.shapeited.vcf.gz"
done

SHAPEIT5_phase_common --input ../../../../102Malay.GRCh38.dv161.chr1_22XY.norm.QC.SNP.vcf.gz \
--region chr${i} \
--map /path/to/chr${i}.genetic_map_hg38.txt2 \
--output 102Malay.GRCh38.dv161.chr1_22XY.norm.QC.SNP.shapeited.bcf \
--thread 32
