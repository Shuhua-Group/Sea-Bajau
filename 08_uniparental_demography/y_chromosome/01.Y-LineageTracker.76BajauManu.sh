#!/usr/bin/env bash
# Analyze mitochondrial and Y-chromosome lineage variation.
# Set input, reference and tool paths for your environment before execution.

conda activate Y_LineageTracker_env


# NRY
# 67West_51East.GRCh38.dv161.chrY.Male.snp.norm.vcf.gz
LineageTracker classify \
--vcf 25malay_8bornean_29purebajau_7mixedbajau.GRCh38.dv161.chrY.snp.norm.vcf.gz \
-b 38 \
--chip \
--snp-only \
--mut-info \
--approximate \
--output 25malay_8bornean_29purebajau_7mixedbajau.GRCh38.dv161.chrY.YLT


LineageTracker cluster \
--hg 25malay_8bornean_29purebajau_7mixedbajau.GRCh38.dv161.chrY.YLT.hapresult.hg \
-p male.list \
--method pca \
--level auto \
--freq \
--output 25malay_8bornean_29purebajau_7mixedbajau.GRCh38.dv161.chrY.YLT.clusterresult


LineageTracker cluster \
--hg 67West_51East.GRCh38.dv161.chrY.Male.snp.norm.YLT.hapresult.hg \
-p male.list \
--method pca \
--level 5 \
--freq \
--output 67West_51East.GRCh38.dv161.chrY.Male.snp.norm.YLT.clusterresult.level5

# 67West_51East.GRCh38.dv161.chrY.Male.snp.norm.YLT.hapresult.hg
LineageTracker phylo \
--hg 23bajau.GRCh38.dv161.chrY.YLT.clusterresult.level.hg \
-p male.list \
--output 23bajau.GRCh38.dv161.chrY.YLT.clusterresult.YLT.phylotree \
--seq 23bajau.GRCh38.dv161.chrY.snp.norm.vcf.gz \
--seq-format vcf \
# The two parameters above caused an error in the original run

LineageTracker phylo \
--hg 67West_51East.GRCh38.dv161.chrY.Male.snp.norm.YLT.hapresult.hg \
-p male.list \
--output 67West_51East.GRCh38.dv161.chrY.Male.snp.norm.YLT.hapresult.phylotree \
# The two parameters below caused an error in the original run
--seq 67West_51East.GRCh38.dv161.chrY.norm.vcf.gz \
--seq-format vcf \

# NRY Haplogroup dating
LineageTracker time \
--seq 67West_51East.GRCh38.dv161.chrY.norm.vcf.gz \
--seq-format vcf \
--tree 67West_51East.GRCh38.dv161.chrY.Male.snp.norm.YLT.hapresult.phylotree.tree.nwk \
--cal CF:69000-73000 \
--output 67West_51East.GRCh38.dv161.chrY.Male.snp.norm.YLT.time


# Y-STR
bcftools norm -m -any --threads 32 \
  -Oz -o 67West_51East.GRCh38.dv161.chr1-22XY.Indel.norm.vcf.gz \
  67West_51East.GRCh38.dv161.chr1-22XY.Indel.vcf.gz

tabix -f 67West_51East.GRCh38.dv161.chr1-22XY.Indel.norm.vcf.gz

bcftools index 67West_51East.GRCh38.dv161.chr1-22XY.Indel.norm.vcf.gz
bcftools view -r chrY -S /path/to/118Sample.male.txt \
67West_51East.GRCh38.dv161.chr1-22XY.Indel.norm.vcf.gz \
-Oz -o 67West_51East.GRCh38.dv161.chrY.Male.Indel.norm.vcf.gz


LineageTracker genostr \
--vcf 67West_51East.GRCh38.dv161.chrY.Male.Indel.norm.vcf.gz \
-b 38 \
--output 67West_51East.GRCh38.dv161.chrY.Male.Indel.norm.YLT.genostr

LineageTracker stat \
--matrix 67West_51East.GRCh38.dv161.chrY.Male.Indel.norm.YLT.genostr.str.matrix \
-p male.list \
--gd fst \
--mpd \
--hg \
--predict \
--amova fst \
--output 67West_51East.GRCh38.dv161.chrY.Male.snp.norm.YLT.stat


LineageTracker net \
--matrix 67West_51East.GRCh38.dv161.chrY.Male.Indel.norm.YLT.genostr.str.matrix \
-p male.list \
--fdi \
--output 67West_51East.GRCh38.dv161.chrY.Male.snp.norm.YLT.net

