#!/usr/bin/env bash
# Prepare three CHM13 graph catalogs for catalog-level comparison.
# Deliberately retains GRCh38 and all other graph donors: this is NOT an AF analysis.
# Process variant calls and summarize variant catalogs.
# Set input, reference and tool paths for your environment before execution.

set -euo pipefail

APC_VCF="/path/to/130APC_CHM13ref_GRCh38plus.wave.vcf.gz"
CPC2_VCF="/path/to/cpc2-20250707-mc-chm13.wave.vcf.gz"
HPRC2_VCF="/path/to/hprc-v2.0-mc-chm13.wave.vcf.gz"
CHM13_FASTA="/PATH/TO/chm13v2.0.fa"
# BED must contain the union of centromeric, telomeric/subtelomeric and satellite
# difficult regions in CHM13 coordinates. It must be bgzip-compressed and tabix-indexed.
DIFFICULT_BED="/PATH/TO/CHM13v2_difficult_regions.bed.gz"
OUTDIR="/path/to/16.threeway_catalog"
SCRIPT_DIR="/path/to/revised_pipeline"
THREADS=12
AUTOSOMES="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"

for exe in bcftools bgzip tabix python3; do command -v "$exe" >/dev/null || { echo "ERROR: $exe unavailable" >&2; exit 1; }; done
[[ -r "$CHM13_FASTA" && -r "${CHM13_FASTA}.fai" ]] || { echo "ERROR: set CHM13_FASTA" >&2; exit 1; }
[[ -r "$DIFFICULT_BED" && -r "${DIFFICULT_BED}.tbi" ]] || { echo "ERROR: provide indexed DIFFICULT_BED" >&2; exit 1; }
[[ -r "$SCRIPT_DIR/classify_biallelic_variant.py" ]] || { echo "ERROR: set SCRIPT_DIR" >&2; exit 1; }
mkdir -p "$OUTDIR" "$OUTDIR/01.prepared" "$OUTDIR/02.keys" "$OUTDIR/03.audit"

prepare_catalog() {
    local label="$1" vcf="$2" out="$OUTDIR/01.prepared/${label}.autosomes.primary.vcf.gz"
    echo "Preparing ${label}" >&2
    # -m -any is mandatory before classifying. -T ^ excludes difficult regions.
    bcftools norm --threads "$THREADS" -f "$CHM13_FASTA" -m -any -c ws "$vcf" -Ou |
      bcftools view --threads "$THREADS" -r "$AUTOSOMES" -m2 -M2 -T "^${DIFFICULT_BED}" -Oz -o "$out"
    tabix -f -p vcf "$out"

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$out" |
      python3 "$SCRIPT_DIR/classify_biallelic_variant.py" > "$OUTDIR/02.keys/${label}.keys.tsv"

    # Header sample inventory is audit metadata only; it does not determine catalog membership.
    {
      printf 'catalog\tsample_id\tis_GRCh38\n'
      bcftools query -l "$vcf" | awk -v x="$label" 'BEGIN{OFS="\t"}{print x,$0,($0=="GRCh38" ? "yes" : "no")}'
    } > "$OUTDIR/03.audit/${label}.sample_inventory.tsv"
    awk 'NR>1{n[$5]++} END{for (x in n) print x,n[x]}' "$OUTDIR/02.keys/${label}.keys.tsv" |
      sort > "$OUTDIR/03.audit/${label}.class_counts.txt"
}

prepare_catalog APC "$APC_VCF"
prepare_catalog CPC2 "$CPC2_VCF"
prepare_catalog HPRC2 "$HPRC2_VCF"

cat "$OUTDIR/03.audit/"*.sample_inventory.tsv | awk 'NR==1 || $1!="catalog"' > "$OUTDIR/03.audit/all_sample_inventory.tsv"
# Disk-backed sort keeps this scalable for tens of millions of alleles.
{
  awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$2,$3,$4,$5,"APC"}' "$OUTDIR/02.keys/APC.keys.tsv"
  awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$2,$3,$4,$5,"CPC2"}' "$OUTDIR/02.keys/CPC2.keys.tsv"
  awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$2,$3,$4,$5,"HPRC2"}' "$OUTDIR/02.keys/HPRC2.keys.tsv"
} | LC_ALL=C sort --parallel="$THREADS" -S 50% -t $'\t' -k1,1 -k2,2n -k3,3 -k4,4 -k6,6 -u \
  > "$OUTDIR/04.exact_membership.source_sorted.tsv"
python3 "$SCRIPT_DIR/make_exact_membership.py" \
  --input "$OUTDIR/04.exact_membership.source_sorted.tsv" \
  --out "$OUTDIR/04.exact_membership"

echo "Done: ${OUTDIR}" >&2
