#!/bin/bash

# Run RELATE Estimate effective population size
# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/RELATE"
cd "$wd" || exit 1
for i in {1..22};do
# Start from vcf Obtain from the file example.haps And example.sample File
/path/to/RelateFileFormats \
                --mode ConvertFromVcf \
                --haps "36Bajau_like.GRCh38.chr${i}.haps" \
                --sample "36Bajau_like.GRCh38.chr${i}.sample" \
                -i "../selection/iHS/36Bajau_like.GRCh38.dv161.chr${i}.norm.QC.SNP.shapeited.reID"
# Prepare input files
/path/to/PrepareInputFiles.sh \
                 --haps "36Bajau_like.GRCh38.chr${i}.haps" \
                 --sample "36Bajau_like.GRCh38.chr${i}.sample" \
                 --ancestor "/path/to/homo_sapiens_ancestor_${i}.fa" \
                 --mask "/path/to/20160622.chr${i}.pilot_mask.fasta.gz" \
                 -o "36Bajau_like.GRCh38.chr${i}_input"
# This N is the effective haplogroup size, Start with defaults 30000 Run a first pass
/path/to/Relate --mode All \
    --haps "36Bajau_like.GRCh38.chr${i}_input.haps.gz" \
    --sample "36Bajau_like.GRCh38.chr${i}_input.sample.gz" \
    --map "/path/to/plink.chr${i}.GRCh38.map" \
    -N 30000 \
    -m 1.25e-8 \
    -o "36Bajau_like.GRCh38.N30000_chr${i}" \
    --seed 1
done &

/path/to/EstimatePopulationSize.sh \
    -i 36Bajau_like.GRCh38.N30000 \
    --first_chr 1 --last_chr 22 \
    -m 1.25e-8 \
    --poplabels 36Bajau_like.GRCh38.poplabels \
    --seed 1 \
    -o 36Bajau_like.GRCh38.N30000.chr1_22.popsize \
    --threads 64

# Then, Use admixture low-quality entries in 9 exclude these samples, Run again
for i in {1..22};do
bcftools view -S /path/to/43Bajau_like.rm7_1stkin.rm_unpure9.list \
../selection/iHS/36Bajau_like.GRCh38.dv161.chr${i}.norm.QC.SNP.shapeited.reID.vcf.gz \
-Oz -o "../selection/iHS/27Bajau.GRCh38.dv161.chr${i}.norm.QC.SNP.shapeited.reID2.vcf.gz"
done
for i in {1..22};do
# Start from vcf Obtain from the file example.haps And example.sample File
/path/to/RelateFileFormats \
                --mode ConvertFromVcf \
                --haps "27Bajau.GRCh38.chr${i}.haps" \
                --sample "27Bajau.GRCh38.chr${i}.sample" \
                -i "../selection/iHS/27Bajau.GRCh38.dv161.chr${i}.norm.QC.SNP.shapeited.reID2"
# Prepare input files
/path/to/PrepareInputFiles.sh \
                 --haps "27Bajau.GRCh38.chr${i}.haps" \
                 --sample "27Bajau.GRCh38.chr${i}.sample" \
                 --ancestor "/path/to/homo_sapiens_ancestor_${i}.fa" \
                 --mask "/path/to/20160622.chr${i}.pilot_mask.fasta.gz" \
                 -o "27Bajau.GRCh38.chr${i}_input"
# This N is the effective haplogroup size, Start with defaults 30000 Run a first pass
/path/to/Relate --mode All \
    --haps "27Bajau.GRCh38.chr${i}_input.haps.gz" \
    --sample "27Bajau.GRCh38.chr${i}_input.sample.gz" \
    --map "/path/to/plink.chr${i}.GRCh38.map" \
    -N 30000 \
    -m 1.25e-8 \
    -o "27Bajau.GRCh38.N30000_chr${i}" \
    --seed 1
done &

/path/to/EstimatePopulationSize.sh \
    -i 27Bajau.GRCh38.N30000 \
    --first_chr 1 --last_chr 22 \
    -m 1.25e-8 \
    --poplabels 27Bajau.GRCh38.poplabels \
    --seed 1 \
    -o 27Bajau.GRCh38.N30000.chr1_22.popsize \
    --threads 30 &
