#!/bin/bash

# Merge first APC/CPC2/HPRC2 three vcf.gz File
# Check the following issue APC And CPC2,APC And HPRC2 Some samples are duplicated, Then update sample list when, Update CPC2, HPRC2 sample names

# Process variant calls and summarize variant catalogs.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/out"
mkdir -p "$wd"
cd "$wd" || exit 1

APC_vcf="/path/to/130APC_CHM13ref_GRCh38plus.wave.vcf.gz"
CPC2_vcf="/path/to/cpc2-20250707-mc-chm13.wave.vcf.gz"
HPRC2_vcf="/path/to/hprc-v2.0-mc-chm13.wave.vcf.gz"
reffa="/path/to/chm13v2.0.fa"

# HPRC2 With mtDNA, Absent from the other two
# for vcf in $APC_vcf $HPRC2_vcf $CPC2_vcf
# do
# echo "========== $vcf =========="
# # bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$vcf" | head
# # bcftools view -h "$vcf" | grep '^##reference'
# bcftools view -h "$vcf" | grep '^##contig'
# done


# (1) Take the intersection; Set-specific differences
## Split first 1-22 X Y Chromosome

for i in `seq 1 22` X Y ;do
mkdir -p "chr$i"
# bcftools view -r chr$i $APC_vcf -O z -o chr$i/APC.chm13v2.chr$i.vcf.gz &
# bcftools view -r chr$i $CPC2_vcf -O z -o chr$i/CPC2.chm13v2.chr$i.vcf.gz &
# bcftools view -r chr$i $HPRC2_vcf -O z -o chr$i/HPRC2.chm13v2.chr$i.vcf.gz &
bcftools norm \
    -f "$reffa" \
    -m -any \
    -c ws \
    -Oz \
    -o "chr$i/APC.chm13v2.chr$i.norm2.vcf.gz" \
    "chr$i/APC.chm13v2.chr$i.vcf.gz" &

bcftools norm \
    -f "$reffa" \
    -m -any \
    -c ws \
    -Oz \
    -o "chr$i/CPC2.chm13v2.chr$i.norm2.vcf.gz" \
    "chr$i/CPC2.chm13v2.chr$i.vcf.gz" &

bcftools norm \
    -f "$reffa" \
    -m -any \
    -c ws \
    -Oz \
    -o "chr$i/HPRC2.chm13v2.chr$i.norm2.vcf.gz" \
    "chr$i/HPRC2.chm13v2.chr$i.vcf.gz"
done

## Then take the intersection
{1..22}
for i in `seq 1 22` X Y;do
mkdir -p "tmp/chr$i"
bcftools sort -Oz -o "chr$i/APC.chm13v2.chr$i.norm2.sorted.vcf.gz" "chr$i/APC.chm13v2.chr$i.norm2.vcf.gz"
bcftools sort -Oz -o "chr$i/CPC2.chm13v2.chr$i.norm2.sorted.vcf.gz" "chr$i/CPC2.chm13v2.chr$i.norm2.vcf.gz"
bcftools sort -Oz -o "chr$i/HPRC2.chm13v2.chr$i.norm2.sorted.vcf.gz" "chr$i/HPRC2.chm13v2.chr$i.norm2.vcf.gz"
tabix -p vcf "chr$i/APC.chm13v2.chr$i.norm2.sorted.vcf.gz"
tabix -p vcf "chr$i/CPC2.chm13v2.chr$i.norm2.sorted.vcf.gz"
tabix -p vcf "chr$i/HPRC2.chm13v2.chr$i.norm2.sorted.vcf.gz"
/path/to/bcftools isec --threads 4 \
chr$i/APC.chm13v2.chr$i.norm2.sorted.vcf.gz "chr$i/CPC2.chm13v2.chr$i.norm2.sorted.vcf.gz" "chr$i/HPRC2.chm13v2.chr$i.norm2.sorted.vcf.gz" \
-n=3 -p "tmp/chr$i"
bgzip "tmp/chr$i/0000.vcf" -c > "tmp/chr$i/0000.vcf.gz" && tabix -f "tmp/chr$i/0000.vcf.gz"
bgzip "tmp/chr$i/0001.vcf" -c > "tmp/chr$i/0001.vcf.gz" && tabix -f "tmp/chr$i/0001.vcf.gz"
bgzip "tmp/chr$i/0002.vcf" -c > "tmp/chr$i/0002.vcf.gz" && tabix -f "tmp/chr$i/0002.vcf.gz"
/path/to/bcftools merge "tmp/chr$i/0000.vcf.gz" "tmp/chr$i/0001.vcf.gz" \
tmp/chr$i/0002.vcf.gz \
--force-samples \
-Oz -o "$wd/APC_CPC2_HPRC2.norm2.isec.chr${i}.vcf.gz"
tabix -p vcf "$wd/APC_CPC2_HPRC2.norm2.isec.chr${i}.vcf.gz" &
done &

## Then calculate set-specific differences
# for i in `seq 1 22`;do
for i in 10;do
mkdir -p "tmp_uniq/chr$i"
/path/to/bcftools isec --threads 4 \
chr$i/APC.chm13v2.chr$i.norm2.sorted.vcf.gz "chr$i/CPC2.chm13v2.chr$i.norm2.sorted.vcf.gz" "chr$i/HPRC2.chm13v2.chr$i.norm2.sorted.vcf.gz" \
-n=1 -p "tmp_uniq/chr$i"
bgzip "tmp_uniq/chr$i/0000.vcf" -c > "tmp_uniq/chr$i/0000.vcf.gz" && tabix -f "tmp_uniq/chr$i/0000.vcf.gz"
bgzip "tmp_uniq/chr$i/0001.vcf" -c > "tmp_uniq/chr$i/0001.vcf.gz" && tabix -f "tmp_uniq/chr$i/0001.vcf.gz"
bgzip "tmp_uniq/chr$i/0002.vcf" -c > "tmp_uniq/chr$i/0002.vcf.gz" && tabix -f "tmp_uniq/chr$i/0002.vcf.gz"
mv "tmp_uniq/chr$i/0000.vcf.gz" "$wd/APC_uniq.norm2.chm13v2.chr$i.vcf.gz"
mv "tmp_uniq/chr$i/0000.vcf.gz.tbi" "$wd/APC_uniq.norm2.chm13v2.chr$i.vcf.gz.tbi"
mv "tmp_uniq/chr$i/0001.vcf.gz" "$wd/CPC2_uniq.norm2.chm13v2.chr$i.vcf.gz"
mv "tmp_uniq/chr$i/0001.vcf.gz.tbi" "$wd/CPC2_uniq.norm2.chm13v2.chr$i.vcf.gz.tbi"
mv "tmp_uniq/chr$i/0002.vcf.gz" "$wd/HPRC2_uniq.norm2.chm13v2.chr$i.vcf.gz"
mv "tmp_uniq/chr$i/0002.vcf.gz.tbi" "$wd/HPRC2_uniq.norm2.chm13v2.chr$i.vcf.gz.tbi"
done &




for type in APC_uniq CPC2_uniq HPRC2_uniq
do

    bcftools concat \
        --threads 8 \
        -Oz \
        -o "${type}.norm2.chm13v2.allchr.vcf.gz" \
        ${type}.norm2.chm13v2.chr{1..22}.vcf.gz \
        "${type}.norm2.chm13v2.chrX.vcf.gz" \
        "${type}.norm2.chm13v2.chrY.vcf.gz"

    tabix -p vcf "${type}.norm2.chm13v2.allchr.vcf.gz"

done


bcftools concat \
    --threads 8 \
    -Oz \
    -o APC_CPC2_HPRC2.norm2.isec.allchr.vcf.gz \
    APC_CPC2_HPRC2.norm2.isec.chr{1..22}.vcf.gz \
    APC_CPC2_HPRC2.norm2.isec.chrX.vcf.gz \
    APC_CPC2_HPRC2.norm2.isec.chrY.vcf.gz

tabix -p vcf APC_CPC2_HPRC2.norm2.isec.allchr.vcf.gz



# bcftools query \
# -f '%CHROM\t%POS\t%REF\t%ALT\n' \
# APC_uniq.norm2.chm13v2.allchr.vcf.gz \
# > APC_uniq.norm2.chm13v2.allchr.alleles.tsv

# python classify_variant.py \
# APC_uniq.norm2.chm13v2.allchr.alleles.tsv \
# APC_uniq.norm2.chm13v2.allchr.class.tsv

# $ cut -f5 APC_uniq.norm2.chm13v2.allchr.class.tsv |
# > tail -n +2 |
# > sort |
# > uniq -c
# 138005 DEL
# 3399469 INDEL
# 382352 INS
# 4945615 SNP

# Compare all
# for i in {1..22} X Y;do
# bcftools query \
# -f '%CHROM\t%POS\t%REF\t%ALT\n' \
# chr${i}/APC.chm13v2.chr${i}.norm2.sorted.vcf.gz \
# >> APC.chm13v2.allchr.norm2.alleles.tsv
# done

# python classify_variant.py \
# APC.chm13v2.allchr.norm2.alleles.tsv \
# APC.chm13v2.allchr.norm2.class.tsv

# cut -f5 APC.chm13v2.allchr.norm2.class.tsv | tail -n +2 | sort | uniq -c
# $ cut -f5 APC.chm13v2.allchr.norm2.class.tsv |
# > tail -n +2 |
# > sort |
# > uniq -c
# 212907 DEL
# 8245797 INDEL
# 510982 INS
# 20663612 SNP

# Extract APC-only genotype, Use GRCh38 Remove

for type in APC_uniq HPRC2_uniq CPC2_uniq ;do
bcftools view \
-S "${type}.list" \
-Oz \
-o "${type}.norm2.chm13v2.allchr.rmGRCh38.vcf.gz" \
${type}.norm2.chm13v2.allchr.vcf.gz

tabix -p vcf "${type}.norm2.chm13v2.allchr.rmGRCh38.vcf.gz"

bcftools +fill-tags \
${type}.norm2.chm13v2.allchr.rmGRCh38.vcf.gz \
-Oz \
-o "${type}.norm2.chm13v2.allchr.rmGRCh38.AF.vcf.gz" \
-- \
-t AC,AN,AF
tabix -p vcf "${type}.norm2.chm13v2.allchr.rmGRCh38.AF.vcf.gz"

done


# Include in the comparison
APC

for type in HPRC2 CPC2;do
for i in {1..22} X Y;do
bcftools sort -Oz -o "chr$i/${type}.chm13v2.chr${i}.norm2.sorted.vcf.gz" \
chr$i/${type}.chm13v2.chr${i}.vcf.gz &
done
bcftools concat \
    --threads 8 \
    -Oz \
    -o "${type}.chm13v2.norm2.allchr.vcf.gz" \
    chr{1..22}/${type}.chm13v2.chr*.norm2.sorted.vcf.gz \
    "chrX/${type}.chm13v2.chrX.norm2.sorted.vcf.gz" \
    "chrY/${type}.chm13v2.chrY.norm2.sorted.vcf.gz"
tabix -p vcf "${type}.chm13v2.norm2.allchr.vcf.gz" &
done

for type in  HPRC2;do
# bcftools view \
# -S ${type}.list \
# -Oz \
# -o ${type}.chm13v2.norm2.allchr.rmGRCh38.vcf.gz \
# ${type}.chm13v2.norm2.allchr.vcf.gz &
# tabix -p vcf ${type}.chm13v2.norm2.allchr.rmGRCh38.vcf.gz &

bcftools +fill-tags \
${type}.chm13v2.norm2.allchr.rmGRCh38.vcf.gz \
-Oz \
-o "${type}.chm13v2.norm2.allchr.rmGRCh38.AF.vcf.gz" \
-- \
-t AC,AN,AF

tabix -p vcf "${type}.chm13v2.norm2.allchr.rmGRCh38.AF.vcf.gz"
done


for type in APC CPC2 HPRC2 ;do
bcftools query \
-f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' \
${type}_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.vcf.gz \
> "${type}_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.raw.tsv"

python classify_variant.py \
${type}_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.raw.tsv \
> "${type}_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.tsv"
done


for type in APC CPC2 HPRC2;do
bcftools query \
-f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' \
${type}.chm13v2.norm2.allchr.rmGRCh38.AF.vcf.gz \
> "${type}.chm13v2.norm2.allchr.rmGRCh38.AF.raw.tsv"

python classify_variant.py \
${type}.chm13v2.norm2.allchr.rmGRCh38.AF.raw.tsv \
> "${type}.chm13v2.norm2.allchr.rmGRCh38.AF.tsv"
done



# APC CPC2
for type in HPRC2;do
python variant_class_composition.py \
-i \
${type}_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.tsv \
${type}.chm13v2.norm2.allchr.rmGRCh38.AF.tsv \
-n \
${type}_uniq \
${type}_all \
-o "${type}_variant_class.rmGRCh38.tsv"
done


awk '$10>=208' APC_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.tsv > APC_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.filter208.tsv
awk '$10>=208' APC.chm13v2.norm2.allchr.rmGRCh38.AF.tsv > APC.chm13v2.norm2.allchr.rmGRCh38.AF.filter208.tsv

awk '$10>=370' CPC2_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.tsv > CPC2_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.filter370.tsv
awk '$10>=370' CPC2.chm13v2.norm2.allchr.rmGRCh38.AF.tsv > CPC2.chm13v2.norm2.allchr.rmGRCh38.AF.filter370.tsv

awk '$10>=370' HPRC2_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.tsv > HPRC2_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.filter370.tsv
awk '$10>=370' HPRC2.chm13v2.norm2.allchr.rmGRCh38.AF.tsv > HPRC2.chm13v2.norm2.allchr.rmGRCh38.AF.filter370.tsv




# Then use R Plot
# library(ggplot2)
# library(dplyr)



# df <- read.table(
# "APC_variant_class.rmGRCh38.tsv.fraction.tsv",
# header=T,
# sep="\t"
#)


# # Set variant class Order
# class_order <- c(
# "COMPLEX",
# "DEL",
# "INS",
# "INDEL",
# "MNP",
# "SNP"
#)


# df$CLASS <- factor(
# df$CLASS,
# levels = class_order
#)


# pdf("APC_variant_class.rmGRCh38.pdf")

# ggplot(
# df,
# aes(
# x = DATASET,
# y = FRACTION,
# fill = CLASS
#)
#) +
# geom_bar(
# stat = "identity",
# position = "stack"
#) +
# theme_classic()
# dev.off()



# Then inspect AFS
# Inspect first AN distribution
# library(dplyr)
# library(ggplot2)


# df <- read.table(
# "APC_unique.norm2.chm13v2.allchr.rmGRCh38.AF.raw.tsv",
# header=FALSE,
# sep="\t"
#)

# colnames(df)=c(
# "CHROM",
# "POS",
# "REF",
# "ALT",
# "AC",
# "AN",
# "AF"
#)


# an <- df %>%
# count(AN)%>%
# mutate(
# fraction=n/sum(n)
#)

# pdf("APC_unique.norm2.chm13v2.allchr.rmGRCh38.AF.raw.AN.pdf")
# ggplot(
# an,
# aes(
# x=AN,
# y=fraction
#)
#)+
# geom_col()+
# theme_classic()+
# xlab("Number of callable alleles (AN)")+
# ylab("Fraction of variants")
# dev.off()

# an_cdf <- df %>%
# count(AN)%>%
# arrange(AN)%>%
# mutate(
# cum_fraction=cumsum(n)/sum(n)
#)

# pdf("APC_unique.norm2.chm13v2.allchr.rmGRCh38.AF.raw.AN.cum.pdf")
# ggplot(
# an_cdf,
# aes(
# AN,
# cum_fraction
#)
#)+
# geom_line()+
# theme_classic()+
# xlab("AN")+
# ylab("Cumulative fraction")
# dev.off()

# Observed pattern:, Most AN remain large . Apply additional filtering, Keep only AN>=0.8*260 Of

# awk '$6>=208' APC_unique.norm2.chm13v2.allchr.rmGRCh38.AF.raw.tsv > APC_unique.norm2.chm13v2.allchr.rmGRCh38.AF.raw.filter208.tsv
# awk '$6>=208' APC.chm13v2.norm2.allchr.rmGRCh38.AF.raw.tsv > APC.chm13v2.norm2.allchr.rmGRCh38.AF.raw.filter208.tsv


awk '$10>=208' APC_unique.norm2.chm13v2.allchr.rmGRCh38.AF.tsv > APC_unique.norm2.chm13v2.allchr.rmGRCh38.AF.filter208.tsv
awk '$10>=208' APC.chm13v2.norm2.allchr.rmGRCh38.AF.tsv > APC.chm13v2.norm2.allchr.rmGRCh38.AF.filter208.tsv

# # Then summarize
# library(data.table)
# library(ggplot2)

# unique <- fread(
# "APC_unique.norm2.chm13v2.allchr.rmGRCh38.AF.filter208.tsv"
#)

# all <- fread(
# "APC.chm13v2.norm2.allchr.rmGRCh38.AF.filter208.tsv"
#)

# unique[,group:="APC-P1-specific"]
# all[,group:="All APC"]


# df <- rbind(unique,all)

# class_order <- c(
# "COMPLEX",
# "DEL",
# "INS",
# "INDEL",
# "MNP",
# "SNP"
#)


# df$CLASS <- factor(
# df$CLASS,
# levels=class_order
#)


# p <- ggplot(
# df,
# aes(
# x=AF,
# color=group
#)
#)+
# geom_density(
# linewidth=1
#)+
# facet_wrap(~ CLASS,scales="free_y")+
# theme_classic()+
# labs(
# x="Allele frequency",
# y="Density"
#)


# ggsave(
# "APC_P1_AF_spectrum_by_class.pdf",
# p,
# width=10,
# height=6
#)

