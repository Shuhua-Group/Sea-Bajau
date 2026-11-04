#!/bin/bash

# Run RELATE Estimate effective population size
# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/RELATE"
cd "$wd" || exit 1
# 20Negrito.ngs 33Malay 29pure.Bajau 7unpure.Bajau 29pure.Bajau
for ind in 69Malay_20Negrito;do
for i in 1;do
# # Start from vcf Obtain from the file example.haps And example.sample File
/path/to/RelateFileFormats \
                --mode ConvertFromVcf \
                --haps "${ind}.chr${i}.haps" \
                --sample "${ind}.chr${i}.sample" \
                -i "../ROH/33Malay_36Bajau_20Negrito/${ind}.chr${i}.QC.shapeited"
# Prepare input files
/path/to/PrepareInputFiles.sh \
                 --haps "${ind}.chr${i}.haps" \
                 --sample "${ind}.chr${i}.sample" \
                 --ancestor "/path/to/homo_sapiens_ancestor_${i}.fa" \
                 --mask "/path/to/20160622.chr${i}.pilot_mask.fasta.gz" \
# -o "${ind}.chr${i}_input"
# This N is the effective haplogroup size, Start with defaults 30000 Run a first pass
/path/to/Relate --mode All \
    --haps "${ind}.chr${i}_input.haps.gz" \
    --sample "${ind}.chr${i}_input.sample.gz" \
    --map "/path/to/plink.chr${i}.GRCh38.map" \
    -N 30000 \
    -m 1.25e-8 \
    -o "${ind}.N30000_chr${i}" \
    --memory 0.1 \
    --seed 1
done
done



conda activate r4
# for ind in 33Malay 29pure.Bajau 20Negrito.ngs 29pure.Bajau 69Malay_20Negrito 7unpure.Bajau;do
for ind in 69Malay_20Negrito;do
# /path/to/EstimatePopulationSize.sh \
# -i ${ind}.N30000 \
# --first_chr 1 --last_chr 22 \
# -m 1.25e-8 \
# --poplabels ${ind}.poplabel \
# --seed 1 \
# -o ${ind}.malaylike.N30000.chr1_22.popsize \
# --threads 8


/path/to/RelateCoalescentRate \
--mode FinalizePopulationSize \
--poplabels 69Malay_20Negrito.borneo.poplabel \
--first_chr 1 --last_chr 22 \
-i "${ind}.N30000.chr1_22.popsize" \
-o "${ind}.N30000.RelateCoalescentRate"
done



# for ind in 69Malay_20Negrito;do
# /path/to/EstimatePopulationSize.sh \
# -i ${ind}.N30000 \
# --first_chr 1 --last_chr 22 \
# -m 1.25e-8 \
# --poplabels ${ind}.borneo.poplabel \
# --seed 1 \
# -o ${ind}.borneo.N30000.chr1_22.popsize \
# --threads 16
# done



for ind in 69Malay_20Negrito;do
/path/to/EstimatePopulationSize.sh \
    -i "${ind}.N30000" \
    --first_chr 1 --last_chr 22 \
    -m 1.25e-8 \
    --poplabels 25malay_20negrito_8bornean_29purebajau_7mixedbajau.poplabel \
    --seed 1 \
    -o "${ind}.borneo.N30000.chr1_22.popsize" \
    --threads 8
done
