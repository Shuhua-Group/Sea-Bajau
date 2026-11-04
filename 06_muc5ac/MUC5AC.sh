#!/usr/bin/env bash
# Generate MUC5AC alignment, CDS translation and intron extraction jobs.
# Run the generated jobs in order: align_paf.sh, trans_protein.sh, extract_seqs.sh.
set -euo pipefail

sample_list="${SAMPLE_LIST:-/path/to/HPRC2.list}"
haplotype_map="${HAPLOTYPE_MAP:-/path/to/HPRC2.hap.list}"
assembly_dir="${ASSEMBLY_DIR:-/path/to}"
work_dir="${WORK_DIR:-/path/to}"
minimap2_bin="${MINIMAP2_BIN:-/path/to/minimap2}"
reference_fasta="${REFERENCE_FASTA:-/path/to/chm13v2.0.fa}"
cds_bed="${CDS_BED:-/path/to/muc5ac.CDS.chm13v2.bed}"
intron_bed="${INTRON_BED:-/path/to/muc5ac.introns31-48.21kb3endflanking.chm13v2_coordinates.bed}"
env_activate="${ENV_ACTIVATE:-/path/to/activate}"

for input in "$sample_list" "$haplotype_map"; do
    if [[ ! -r "$input" ]]; then
        printf 'Cannot read input file: %s\n' "$input" >&2
        exit 1
    fi
done
mkdir -p "$work_dir"

write_job_header() {
    local job_name=$1 cpus=$2 name
    printf '#!/usr/bin/env bash\n'
    printf '#SBATCH --job-name=%s\n' "$job_name"
    printf '#SBATCH --partition=batch\n'
    printf '#SBATCH --cpus-per-task=%s\n' "$cpus"
    printf '#SBATCH --mem=32G\n'
    printf 'set -euo pipefail\n\n'
    # Shell-escape paths before embedding them in generated scripts.
    for name in sample outdir asm1 asm2 minimap2_bin reference_fasta cds_bed intron_bed env_activate; do
        printf '%s=%q\n' "$name" "${!name}"
    done
    printf 'cd "$outdir"\n\n'
}

# Use a dedicated descriptor so child commands cannot consume the sample list.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    sample=${sample%$'\r'}
    [[ -z "$sample" || "$sample" == \#* ]] && continue
    if [[ "$sample" == */* || "$sample" == "." || "$sample" == ".." ]]; then
        printf 'Invalid sample identifier: %s\n' "$sample" >&2
        exit 1
    fi
    # Map columns: sample ID, assembly name. Keep the first/last haplotype convention.
    ori_hap1=$(awk -v id="$sample" '$1 == id {print $2; exit}' "$haplotype_map")
    ori_hap2=$(awk -v id="$sample" '$1 == id {hap=$2} END {print hap}' "$haplotype_map")
    if [[ -z "$ori_hap1" || -z "$ori_hap2" || "$ori_hap1" == "$ori_hap2" ]]; then
        printf 'Expected two distinct assemblies for sample: %s\n' "$sample" >&2
        exit 1
    fi
    hap1=$(basename "$ori_hap1" .gz)
    hap2=$(basename "$ori_hap2" .gz)
    asm1="$assembly_dir/$hap1"
    asm2="$assembly_dir/$hap2"
    outdir="$work_dir/$sample"
    mkdir -p "$outdir"

    # Stage 1: align both haplotypes to CHM13.
    {
        write_job_header "${sample}.muc5ac.align" 16
        cat <<'ALIGN'
for hap in 1 2; do
    assembly_var="asm${hap}"
    "$minimap2_bin" -cx asm20 --cs -r2k -t 16 \
        "$reference_fasta" "${!assembly_var}" \
        > "$outdir/$sample.hap${hap}.chm13v2.minimap.paf"
done
ALIGN
    } > "$outdir/align_paf.sh"

    # Stage 2: lift CDS coordinates, concatenate exons and translate to protein.
    {
        write_job_header "${sample}.muc5ac.protein" 8
        cat <<'CDS'
for hap in 1 2; do
    assembly_var="asm${hap}"
    paf="$outdir/$sample.hap${hap}.chm13v2.minimap.paf"
    prefix="$outdir/muc5ac.CDS.hap${hap}"
    source "$env_activate" rustybam
    rustybam liftover --bed "$cds_bed" "$paf" > "$prefix.paf"
    awk 'BEGIN {OFS="\t"} {print $1, $3, $4, ".", "0", $5}' \
        "$prefix.paf" > "$prefix.bed"
    sort -k1,1V -k2,2n "$prefix.bed" > "$prefix.sorted.bed"
    source "$env_activate" seqkit
    seqkit subseq --bed "$prefix.sorted.bed" "${!assembly_var}" > "$prefix.fa"
    samtools faidx "$prefix.fa"
    complete="$outdir/muc5ac.complete_CDS.hap${hap}.fa"
    awk -v name="$sample.hap${hap}.muc5ac.CDS" \
        'BEGIN {print ">" name} /^>/ {next} {printf "%s", $0} END {print ""}' \
        "$prefix.fa" > "$complete"
    samtools faidx "$complete"
    source "$env_activate" emboss
    transeq -sequence "$complete" \
        -outseq "$outdir/$sample.muc5ac.protein.hap${hap}.fa"
done
CDS
    } > "$outdir/trans_protein.sh"

    # Stage 3: extract introns 31-48 and the 21 kb downstream flank.
    {
        write_job_header "${sample}.muc5ac.introns" 8
        cat <<'INTRONS'
for hap in 1 2; do
    assembly_var="asm${hap}"
    paf="$outdir/$sample.hap${hap}.chm13v2.minimap.paf"
    prefix="$outdir/$sample.hap${hap}.muc5ac.intron31-48.21kb3endflanking"
    source "$env_activate" rustybam
    rustybam liftover --bed "$intron_bed" "$paf" > "$prefix.paf"
    awk 'BEGIN {OFS="\t"} {print $1, $3, $4, ".", "0", $5}' \
        "$prefix.paf" > "$prefix.bed"
    sort -k1,1V -k2,2n "$prefix.bed" > "$prefix.sorted.bed"
    source "$env_activate" seqkit
    seqkit subseq --bed "$prefix.sorted.bed" "${!assembly_var}" > "$prefix.fa"
    samtools faidx "$prefix.fa"
    awk -v name="$sample.hap${hap}.muc5ac.intron31-48.21kb3endflanking" \
        'BEGIN {print ">" name} /^>/ {next} {printf "%s", $0} END {print ""}' \
        "$prefix.fa" > "$prefix.connected.fa"
done
INTRONS
    } > "$outdir/extract_seqs.sh"
done 9< "$sample_list"
