#!/usr/bin/env bash
# Merge Bajau HiFi SNPs with Ho.1926 GRCh38 SNPs and project Bajau onto Ho.1926 PCs.
# Existing Bajau VCF QC, including MIN(FMT/DP) < 10, is retained unchanged.

# Prepare population datasets and summarize genetic structure.
# Set input, reference and tool paths for your environment before execution.

set -euo pipefail

# ---------------------------------------------------------------------------
# 0. Paths and programs
# ---------------------------------------------------------------------------
wd="/path/to/76Bajau"
reffa="/path/to/Homo_sapiens.GRCh38.dna.chromosome.1_22_X_Y.fa"
ho_template="/path/to/Ho.1926.merge.chr%s.vcf.gz"
PLINK2="/path/to/plink2"
PLINK1="/path/to/plink"
export EIGENSOFT="/path/to/EIG-7.2.1"
export PATH="$EIGENSOFT/bin:$PATH"
threads=8

outroot="${wd}/07.Bajau_ho1926_projection"
mkdir -p "$outroot"/{01.ho_prep,02.plink_prep,03.merge_ld,04.smartpca}

command -v bcftools >/dev/null
command -v tabix >/dev/null
command -v "$PLINK1" >/dev/null
command -v convertf >/dev/null
command -v smartpca >/dev/null

[[ -x "$PLINK2" ]] || { echo "ERROR: cannot execute $PLINK2" >&2; exit 1; }

# Existing Bajau binary dataset, produced from the strictly QC-filtered Bajau VCF.
Bajau_bfile="${wd}/02.plink_prep/76Bajau.QC"
[[ -r "${Bajau_bfile}.bed" && -r "${Bajau_bfile}.bim" && -r "${Bajau_bfile}.fam" ]] || {
  echo "ERROR: cannot find ${Bajau_bfile}.{bed,bim,fam}" >&2
  exit 1
}

# ---------------------------------------------------------------------------
# 1. Ho.1926 VCF preparation: concatenate chr1-22 and harmonize SNP notation
# ---------------------------------------------------------------------------
cd "$outroot/01.ho_prep" || exit 1

: > Ho.1926.chr1_22.vcfs.list
for chr in $(seq 1 22); do
  vcf=$(printf "$ho_template" "$chr")
  [[ -r "$vcf" ]] || { echo "ERROR: missing Ho.1926 VCF: $vcf" >&2; exit 1; }
  printf '%s\n' "$vcf" >> Ho.1926.chr1_22.vcfs.list
done


bcftools concat --threads "$threads" -f Ho.1926.chr1_22.vcfs.list -Oz \
  -o Ho.1926.GRCh38.chr1_22.concat.vcf.gz
tabix -f -p vcf Ho.1926.GRCh38.chr1_22.concat.vcf.gz

# Do not impose Bajau HiFi DP/QUAL filters on an external reference panel.
bcftools norm --threads "$threads" -f "$reffa" -m -any -Oz \
  -o Ho.1926.GRCh38.chr1_22.norm.vcf.gz \
  Ho.1926.GRCh38.chr1_22.concat.vcf.gz
tabix -f -p vcf Ho.1926.GRCh38.chr1_22.norm.vcf.gz

# This expects chr1...chr22 contig names. Rename Ho contigs first if needed.
bcftools view --threads "$threads" -v snps -m2 -M2 \
  -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
  -Oz -o Ho.1926.GRCh38.chrA.norm.biallelic.snps.vcf.gz \
  Ho.1926.GRCh38.chr1_22.norm.vcf.gz
tabix -f -p vcf Ho.1926.GRCh38.chrA.norm.biallelic.snps.vcf.gz

# ---------------------------------------------------------------------------
# 2. PLINK conversion and Bajau sample-ID protection
# ---------------------------------------------------------------------------
cd "$outroot/02.plink_prep" || exit 1

"$PLINK2" --vcf ../01.ho_prep/Ho.1926.GRCh38.chrA.norm.biallelic.snps.vcf.gz \
  --double-id --allow-extra-chr --set-all-var-ids '@:#:$r:$a' \
  --vcf-half-call missing --mind 0.10 --geno 0.05 \
  --make-bed --out Ho.1926.QC

# Prefix Bajau IDs only in this Bajau-Ho analysis to prevent cross-panel duplicate
# IDs from being treated as one individual during the merge.
awk '{print $1, $2, "Bajau_" $2, "Bajau_" $2}' "${Bajau_bfile}.fam" > Bajau.prefix_ids.txt
"$PLINK2" --bfile "$Bajau_bfile" --update-ids Bajau.prefix_ids.txt \
  --make-bed --out 76Bajau.QC.prefixed

n_Bajau=$(awk 'END {print NR}' 76Bajau.QC.prefixed.fam)
[[ "$n_Bajau" -eq 76 ]] || {
  echo "ERROR: expected 76 Bajau samples, found $n_Bajau." >&2
  exit 1
}

# ---------------------------------------------------------------------------
# 3. Exact SNP harmonization, merge, MAF filtering, and LD pruning
# ---------------------------------------------------------------------------
cd "$outroot/03.merge_ld" || exit 1

# IDs encode CHROM:POS:REF:ALT in both datasets, so matching IDs require exact
# coordinate and allele agreement.
awk 'FNR == NR {seen[$2] = 1; next} ($2 in seen) {print $2}' \
  ../02.plink_prep/76Bajau.QC.prefixed.bim \
  ../02.plink_prep/Ho.1926.QC.bim | sort -u > common_snps.txt

"$PLINK2" --bfile ../02.plink_prep/76Bajau.QC.prefixed \
  --extract common_snps.txt --make-bed --out 76Bajau.common
"$PLINK2" --bfile ../02.plink_prep/Ho.1926.QC \
  --extract common_snps.txt --make-bed --out Ho.1926.common

"$PLINK1" --bfile Ho.1926.common --bmerge 76Bajau.common \
  --make-bed --out Merge.common

n_ho=$(awk 'END {print NR}' Ho.1926.common.fam)
n_merged=$(awk 'END {print NR}' Merge.common.fam)
[[ "$n_merged" -eq $((n_ho + n_Bajau)) ]] || {
  echo "ERROR: expected $((n_ho + n_Bajau)) merged samples, found $n_merged." >&2
  exit 1
}

for maf in 0.01 0.02 0.05; do
  prefix="Merge.MAF${maf}"
  "$PLINK2" --bfile Merge.common --maf "$maf" --make-bed --out "$prefix"
  "$PLINK2" --bfile "$prefix" --indep-pairwise 50 5 0.2 --out "${prefix}.LD"
  "$PLINK2" --bfile "$prefix" --extract "${prefix}.LD.prune.in" \
    --make-bed --out "${prefix}.LDpruned"
done

# ---------------------------------------------------------------------------
# 4. Ho.1926-reference projection PCA
# ---------------------------------------------------------------------------
cd "$outroot/04.smartpca" || exit 1

for maf in 0.01 0.02 0.05; do
  prefix="Merge.MAF${maf}"
  fam="../03.merge_ld/${prefix}.LDpruned.fam"

  # Bajau IDs were uniquely prefixed. All remaining individuals are Ho references.
  awk 'BEGIN {OFS="\t"} $2 ~ /^Bajau_/ {print $2, "U", "Bajau_Pop"; next} {print $2, "U", "HO_reference"}' \
    "$fam" > "${prefix}.ind"

  actual_Bajau=$(awk '$3 == "Bajau_Pop" {n++} END {print n+0}' "${prefix}.ind")
  [[ "$actual_Bajau" -eq 76 ]] || {
    echo "ERROR: expected 76 Bajau samples in ${prefix}.ind, found $actual_Bajau." >&2
    exit 1
  }

  printf 'HO_reference\n' > "Ho_reference_pops.MAF${maf}.txt"
  ln -sf "../03.merge_ld/${prefix}.LDpruned.bed" "${prefix}.bed"
  ln -sf "../03.merge_ld/${prefix}.LDpruned.bim" "${prefix}.bim"

  cat > "convertf.MAF${maf}.par" <<EOF
genotypename:    ${prefix}.bed
snpname:         ${prefix}.bim
indivname:       ${prefix}.ind
outputformat:    EIGENSTRAT
genotypeoutname: ${prefix}.geno
snpoutname:      ${prefix}.snp
indivoutname:    ${prefix}.ind.eigen
familynames:     NO
EOF
  convertf -p "convertf.MAF${maf}.par" > "convertf.MAF${maf}.log" 2>&1

  cat > "smartpca.MAF${maf}.par" <<EOF
genotypename:    ${prefix}.geno
snpname:         ${prefix}.snp
indivname:       ${prefix}.ind.eigen
evecoutname:     ${prefix}.evec
evaloutname:     ${prefix}.eval
poplistname:     Ho_reference_pops.MAF${maf}.txt
lsqproject:      YES
numoutevec:      20
outliermode:     0
numthreads:      ${threads}
EOF
  smartpca -p "smartpca.MAF${maf}.par" > "smartpca.MAF${maf}.log" 2>&1

  printf 'MAF=%s\nHO_reference_samples\t' "$maf"
  awk '$3 == "HO_reference" {n++} END {print n+0}' "${prefix}.ind.eigen"
  printf 'Bajau_projected_samples\t'
  awk '$3 == "Bajau_Pop" {n++} END {print n+0}' "${prefix}.ind.eigen"
  printf 'LD_pruned_SNPs\t'
  awk 'END {print NR}' "../03.merge_ld/${prefix}.LDpruned.bim"
done
