#!/bin/bash
# -------------------------------
# TreeMix preparation workflow
# -------------------------------

# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/TreeMix"
PREFIX="4panel.chrA.shapeit5.selected221"

cd "$wd" || exit 1

# -------------------------------
# Step 0: Handle sites with all genotypes missing
# -------------------------------
vcftools --gzvcf "$PREFIX.vcf.gz" \
         --max-missing 1 \
         --recode --stdout | bgzip > "$PREFIX.noN.vcf.gz"

# -------------------------------
# Step 1: LD pruning
# -------------------------------
/path/to/plink2 --vcf "$PREFIX.noN.vcf.gz" \
       --allow-extra-chr \
       --indep-pairwise 50 10 0.2 \
       --out "$PREFIX.LDpruned"

# -------------------------------
# Step 2: Extract LD pruned Sites, Keep sample names unchanged
# -------------------------------
/path/to/plink2 --vcf "$PREFIX.noN.vcf.gz" \
       --allow-extra-chr \
       --extract "$PREFIX.LDpruned.prune.in" \
       --const-fid 0 \   #  Ensure FID does not affect IID
       --recode vcf bgz \
       --out "$PREFIX.LDpruned"

# ------------------------
# 3. Convert PLINK1.9 bed Format (Keep IID Unchanged)
# ------------------------
plink --vcf "$PREFIX.LDpruned.vcf.gz" \
--make-bed \
--allow-extra-chr \
--const-fid 0 \
--out "$PREFIX.LDpruned.treemix"

# ------------------------
# 4. Use PLINK1.9 Calculate per population allele counts
# ------------------------
# Ensure $PREFIX.clust Valid format: FID IID POP
# FID Fill all entries 0 Corresponding --const-fid
# PLINK1.9 Calculate per population MAF
plink --bfile "$PREFIX.LDpruned.treemix" \
      --freq \
      --within "$PREFIX.clust.fixed" \
      --allow-extra-chr \
      --out "$PREFIX.treemix.pop"

python frqstrat2treemix.py "$PREFIX.treemix.pop.frq.strat" "$PREFIX.treemix.input"
gzip "$PREFIX.treemix.input"


for i in {0..5}
do
treemix -i "$PREFIX.treemix.input.gz" \
-m "$i" \
-o "$PREFIX.$i" \
-root 36Bajau \
-bootstrap 100 \
-k 500 > "treemix_${i}_log" &
done

for i in {1..22}
do
for rep in 1; do
treemix -i testout.txt.gz \
-m "$i" \
-o "$PREFIX.${i}_rep${rep}" \
-root Bajau \
-bootstrap 100 \
-k 500 > "treemix_${i}_rep${rep}_log" &
done
done

for m in {0..5}; do
  for rep in {1..5}; do   #  At least  3-5  Times
    treemix -i mydata.treemix.gz -m "$m" -o "TreeMix/m${m}_rep${rep}" -k 500
  done
done
