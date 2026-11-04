#!/bin/bash

# Pacific data are array genotypes !
# Requires additional preprocessing ......

# Process variant calls and summarize variant catalogs.
# Set input, reference and tool paths for your environment before execution.

for pop in IDN BOR MMW MNB;do
zcat "Pacific/201805_Oceania_Stoneking/${pop}.uni.b37.inp.gz" | awk '
BEGIN {
    FS=OFS="\t";
    print "##fileformat=VCFv4.2";
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
}
NR == 1 {
    # Extract sample names (Number 6 column onward)
    split($0, header, "\t");
    samples = "";
    for (i=6; i<=length(header); i++) {
        samples = samples header[i] "\t";
    }
    gsub(/\t$/, "", samples);  #  Remove trailing tabs
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" samples;
    next;
}
{
    chrom = $2;
    pos = $3;
    id = $1;
    split($5, alleles, "/");
    ref = alleles[1];
    alt = alleles[2];

    # Generate genotype fields (A→0/0, B→1/1, H→0/1, U→./.)
    format = "GT";
    printf "%s\t%s\t%s\t%s\t%s\t.\t.\t.\t%s", chrom, pos, id, ref, alt, format;
    for (i=6; i<=NF; i++) {
        val = $i;
        if (val == "A") {
            gt = "0/0";
        } else if (val == "B") {
            gt = "1/1";
        } else if (val == "H") {
            gt = "0/1";
        } else {
            gt = "./.";
        }
        printf "\t%s", gt;
    }
    printf "\n";
}' | bgzip -c > "${pop}.uni.b37.vcf.gz"
done
awk '{print "##contig=<ID=" $1 ",length=" $2 ">"}' \
/path/to/human_g1k_v37.fasta.fai > b37_contigs.txt

for pop in IDN BOR MMW MNB;do
#bcftools annotate -h b37_contigs.txt ${pop}.uni.b37.vcf.gz -Oz -o ${pop}.uni.b37.h.vcf.gz
tabix -f "${pop}.uni.b37.h.vcf.gz"
done



# Convert to hg38 coordinates
# Unexpected output !
# 2025-04-09 01:33:50 [INFO] Total entries: 871188
# 2025-04-09 01:33:50 [INFO] Failed to map: 435092
# Only half the sites were converted
# Because inp.gz For some sites REF/ALT orientation is reversed, No orientation flag is available !!!!
# Update unmap site REF And ALT Reverse orientation, corresponding GT Also convert, Then chain. Successful, A few hundred sites remain unconverted


for pop in IDN BOR MMW MNB
do
CrossMap vcf /path/to/hg19ToHg38.over.chain.gz \
${pop}.uni.b37.h.vcf.gz \
/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
${pop}.uni.hg38.vcf
# Adjust REF/ALT and genotypes
awk '
BEGIN {OFS="\t";}
/^##/ {print; next;}
/^#CHROM/ {print; next;}
{
    ref = $4;
    $4 = $5;
    $5 = ref;
    for (i=10; i<=NF; i++) {
        if ($i ~ /0\/0/) $i = "1/1";
        else if ($i ~ /1\/1/) $i = "0/0";
        else if ($i ~ /0\/1/) $i = "1/0";
        else if ($i ~ /1\/0/) $i = "0/1";
    }
    print;
}' "${pop}.uni.hg38.vcf.unmap" | bgzip -c > "${pop}.uni.b37.trans.vcf.gz"
CrossMap vcf /path/to/hg19ToHg38.over.chain.gz \
${pop}.uni.b37.trans.vcf.gz \
/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
${pop}.uni.hg38.trans.vcf
# Combine sites converted on the first pass ${pop}.uni.hg38.vcf with those converted on the second pass ${pop}.uni.hg38.trans.vcf
bgzip "${pop}.uni.hg38.vcf"
bgzip "${pop}.uni.hg38.trans.vcf"
bcftools sort "${pop}.uni.hg38.vcf.gz" \
-Oz -o "${pop}.uni.hg38.sorted.vcf.gz"
bcftools sort "${pop}.uni.hg38.trans.vcf.gz" \
-Oz -o "${pop}.uni.hg38.trans.sorted.vcf.gz"
tabix -f "${pop}.uni.hg38.sorted.vcf.gz"
tabix -f "${pop}.uni.hg38.trans.sorted.vcf.gz"
bcftools concat -a --threads 16 "${pop}.uni.hg38.sorted.vcf.gz" "${pop}.uni.hg38.trans.sorted.vcf.gz" \
-Oz -o "${pop}.uni.hg38.temp.vcf.gz"
#sort
bcftools sort "${pop}.uni.hg38.temp.vcf.gz" \
-Oz -o "${pop}.uni.hg38.temp.sorted.vcf.gz"
tabix -f "${pop}.uni.hg38.temp.sorted.vcf.gz"
# Add chromosome prefixes
bcftools annotate \
  --rename-chrs /path/to/chr_name_mapping.txt \
  --threads 16 \
  -Oz -o "${pop}.uni.hg38.chr.sorted.vcf.gz" \
  "${pop}.uni.hg38.temp.sorted.vcf.gz"
tabix -f "${pop}.uni.hg38.chr.sorted.vcf.gz"
# Remove temporary files
rm "${pop}.uni.hg38.temp.vcf.gz"
rm "${pop}.uni.hg38.temp.sorted.vcf.gz"
rm "${pop}.uni.hg38.temp.sorted.vcf.gz.tbi"
rm "${pop}.uni.hg38.vcf.gz"
rm "${pop}.uni.hg38.trans.vcf.gz"
rm "${pop}.uni.b37.trans.vcf.gz"
rm "${pop}.uni.hg38.trans.vcf.unmap"
rm "${pop}.uni.hg38.vcf.unmap"
rm "${pop}.uni.hg38.sorted.vcf.gz"
rm "${pop}.uni.hg38.sorted.vcf.gz.tbi"
rm "${pop}.uni.hg38.trans.sorted.vcf.gz"
rm "${pop}.uni.hg38.trans.sorted.vcf.gz.tbi"
rm "${pop}.uni.b37.h.vcf.gz"
rm "${pop}.uni.b37.h.vcf.gz.tbi"
done

# Results file ${pop}.uni.hg38.chr.sorted.vcf.gz
# Merge the files

bcftools merge -Oz -o Pacific.uni.hg38.chr.sorted.vcf.gz \
IDN.uni.hg38.chr.sorted.vcf.gz \
BOR.uni.hg38.chr.sorted.vcf.gz \
MMW.uni.hg38.chr.sorted.vcf.gz \
MNB.uni.hg38.chr.sorted.vcf.gz
tabix -f Pacific.uni.hg38.chr.sorted.vcf.gz

bcftools view Pacific.uni.hg38.chr.sorted.vcf.gz \
-S sample.selected.list \
-Oz -o Pacific.uni.hg38.chr.sorted.selected.vcf.gz


