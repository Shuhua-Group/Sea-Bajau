#!/usr/bin/env bash
# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

cd /path/to/33Malay_36Bajau_20Negrito || exit 1

for i in {1..22};do
# Step 1: Data quality control (Filter MAF < 0.05, Missingness rate > 10% sites)
bcftools view \
  -i 'QUAL >= 20 && MAF > 0.05 && F_MISSING < 0.1' \
  --threads 30 \
  -Oz -o "69Malay_20Negrito.chr${i}.QC.vcf.gz" \
  "69Malay_20Negrito.chr${i}.vcf.gz"
tabix -p vcf "69Malay_20Negrito.chr${i}.QC.vcf.gz"
bcftools index "69Malay_20Negrito.chr${i}.QC.vcf.gz"
done

# Step 2: Merge
bcftools concat \
  -Oz -o 69Malay_20Negrito.chrA.QC.vcf.gz \
  69Malay_20Negrito.chr*.QC.vcf.gz
tabix -p vcf 69Malay_20Negrito.chrA.QC.vcf.gz
# Step 3: Calculate ROH
bcftools roh \
  -G 30 --threads 32 \
  -o 69Malay_20Negrito.chrA.QC.roh.txt \
  69Malay_20Negrito.chrA.QC.vcf.gz

grep ^RG 69Malay_20Negrito.chrA.QC.roh.txt > 69Malay_20Negrito.chrA.QC.roh.RG.txt


# Alternatively use plink Calculate ROH
plink --vcf 69Malay_20Negrito.chrA.QC.vcf.gz \
      --make-bed --vcf-half-call m \
      --out 69Malay_20Negrito.chrA.QC

plink --bfile 69Malay_20Negrito.chrA.QC \
      --homozyg \
      --homozyg-window-snp 50 \
      --homozyg-window-het 1 \
      --homozyg-window-missing 5 \
      --homozyg-snp 100 \
      --homozyg-kb 5000 \
      --out 69Malay_20Negrito.chrA.QC.plink_5Mb_roh


# Step 4: Calculate genetic diversity (Use vcftools)
# Analyze each population separately
# bcftools view -S /path/to/43Bajau_like.list \
# -Oz -o 43Bajau_like.GRCh38.dv161.chr1-22.snp.norm.QC.vcf.gz \
# 67West_51East.GRCh38.dv161.chr1-22.snp.norm.QC.vcf.gz
# bcftools view -S /path/to/38Negrito_like.list \
# -Oz -o 38Negrito_like.GRCh38.dv161.chr1-22.snp.norm.QC.vcf.gz \
# 67West_51East.GRCh38.dv161.chr1-22.snp.norm.QC.vcf.gz
# bcftools view -S /path/to/37Malay_like.list \
# -Oz -o 37Malay_like.GRCh38.dv161.chr1-22.snp.norm.QC.vcf.gz \
# 67West_51East.GRCh38.dv161.chr1-22.snp.norm.QC.vcf.gz
# # Calculate theta pi
# vcftools \
# --gzvcf 43Bajau_like.GRCh38.dv161.chr1-22.snp.norm.QC.vcf.gz \
# --window-pi 10000 \
# --out 43Bajau_like.GRCh38.dv161.chr1-22.snp.norm.QC.10kwindow.pi

# vcftools \
# --gzvcf 38Negrito_like.GRCh38.dv161.chr1-22.snp.norm.QC.vcf.gz \
# --window-pi 10000 \
# --out 38Negrito_like.GRCh38.dv161.chr1-22.snp.norm.QC.10kwindow.pi
# vcftools \
# --gzvcf 37Malay_like.GRCh38.dv161.chr1-22.snp.norm.QC.vcf.gz \
# --window-pi 10000 \
# --out 37Malay_like.GRCh38.dv161.chr1-22.snp.norm.QC.10kwindow.pi

awk '
BEGIN {
    # Define human chromosome lengths (GRCh38)
    genome_length["chr1"] = 248956422;
    genome_length["chr2"] = 242193529;
    genome_length["chr3"] = 198295559;
    genome_length["chr4"] = 190214555;
    genome_length["chr5"] = 181538259;
    genome_length["chr6"] = 170805979;
    genome_length["chr7"] = 159345973;
    genome_length["chr8"] = 145138636;
    genome_length["chr9"] = 138394717;
    genome_length["chr10"] = 133797422;
    genome_length["chr11"] = 135086622;
    genome_length["chr12"] = 133275309;
    genome_length["chr13"] = 114364328;
    genome_length["chr14"] = 107043718;
    genome_length["chr15"] = 101991189;
    genome_length["chr16"] = 90338345;
    genome_length["chr17"] = 83257441;
    genome_length["chr18"] = 80373285;
    genome_length["chr19"] = 58617616;
    genome_length["chr20"] = 64444167;
    genome_length["chr21"] = 46709983;
    genome_length["chr22"] = 50818468;
}

# Process ROH Data (RG Row)
$1 == "RG" {
    sample = $2;
    chrom = $3;
    roh_length = $6 + 0;  #  Coerce to numeric ,  Avoid errors in some  awk  parsing

    # Accumulate totals: ROH Length
    roh_total[sample][chrom] += roh_length;

    # According to ROH Length categories
    if (roh_length < 1000000) {
        short_roh[sample][chrom] += roh_length;
    } else if (roh_length >= 1000000 && roh_length <= 5000000) {
        medium_roh[sample][chrom] += roh_length;
    } else {
        long_roh[sample][chrom] += roh_length;
    }

    # Record the genome length (Only calculate chromosomes with ROH segments)
    genome_used[sample] += (chrom in genome_length ? genome_length[chrom] : 0);
}

END {
    print "Sample\tChromosome\tTotal_ROH_Length\tFROH\tShort_ROH\tMedium_ROH\tLong_ROH";
    for (s in roh_total) {
        for (c in roh_total[s]) {
            total_roh = roh_total[s][c];
            froh = (c in genome_length) ? total_roh / genome_length[c] : "NA";
            short = (s in short_roh && c in short_roh[s]) ? short_roh[s][c] : 0;
            medium = (s in medium_roh && c in medium_roh[s]) ? medium_roh[s][c] : 0;
            long = (s in long_roh && c in long_roh[s]) ? long_roh[s][c] : 0;
            print s, c, total_roh, froh, short, medium, long;
        }
    }
}' 69Malay_20Negrito.chrA.QC.roh.RG.txt > 69Malay_20Negrito.chrA.QC.roh.RG.roh_summary.txt

# Summarize across chromosomes
awk '
BEGIN {
    # Define the total human genome length (GRCh38)
    genome_total = 0;
    genome_total += 248956422;  # chr1
    genome_total += 242193529;  # chr2
    genome_total += 198295559;  # chr3
    genome_total += 190214555;  # chr4
    genome_total += 181538259;  # chr5
    genome_total += 170805979;  # chr6
    genome_total += 159345973;  # chr7
    genome_total += 145138636;  # chr8
    genome_total += 138394717;  # chr9
    genome_total += 133797422;  # chr10
    genome_total += 135086622;  # chr11
    genome_total += 133275309;  # chr12
    genome_total += 114364328;  # chr13
    genome_total += 107043718;  # chr14
    genome_total += 101991189;  # chr15
    genome_total += 90338345;   # chr16
    genome_total += 83257441;   # chr17
    genome_total += 80373285;   # chr18
    genome_total += 58617616;   # chr19
    genome_total += 64444167;   # chr20
    genome_total += 46709983;   # chr21
    genome_total += 50818468;   # chr22
}

# Process ROH Data (RG Row)
$1 == "RG" {
    sample = $2;
    roh_length = $6 + 0;  #  Coerce to numeric ,  Avoid errors in some  awk  parsing

    # Accumulate totals: ROH Length
    roh_total[sample] += roh_length;

    # According to ROH Length categories
    if (roh_length < 1000000) {
        short_roh[sample] += roh_length;
    } else if (roh_length >= 1000000 && roh_length <= 5000000) {
        medium_roh[sample] += roh_length;
    } else {
        long_roh[sample] += roh_length;
    }
}

END {
    print "Sample\tTotal_ROH_Length\tFROH\tShort_ROH\tMedium_ROH\tLong_ROH";
    for (s in roh_total) {
        total_roh = roh_total[s];
        froh = total_roh / genome_total;  #  Calculate  FROH
        short = (s in short_roh) ? short_roh[s] : 0;
        medium = (s in medium_roh) ? medium_roh[s] : 0;
        long = (s in long_roh) ? long_roh[s] : 0;
        print s, total_roh, froh, short, medium, long;
    }
}' 69Malay_20Negrito.chrA.QC.roh.RG.txt > 69Malay_20Negrito.chrA.QC.roh.RG.roh_sample_summary.txt
sed -i "s/ /\t/g" 69Malay_20Negrito.chrA.QC.roh.RG.roh_sample_summary.txt
