#!/usr/bin/env bash
# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

set -euo pipefail

# ================================
# Parameters to configure
# ================================
POPLABEL="25malay_20negrito_8bornean_29purebajau_7mixedbajau.poplabel"
VCF_DIR="/path/to/33Malay_36Bajau_20Negrito"           #  Store VCF directory containing the files
OUT_DIR="/path/to/data"                   #  Output directory
CHROM_LENGTHS=(248956422 242193529 198295559 190214555 181538259 170805979 159345973 145138636 138394717 133797422 135086622 133275309 114364328 107043718 101991189 90338345 83257441 80373285 58617616 64444167 46709983 50818468)  # GRCh38  Chromosome 1-22 Length

# Create the output directory if needed
mkdir -p "$OUT_DIR"

# ================================
# Generate dict
# ================================
eval "$(
awk 'NR>1 {pop[$2]=(pop[$2] ? pop[$2]","$1 : $1)}
     END {
         print "declare -A dict"
         for (p in pop) {
             print "dict[\""p"\"]=\""pop[p]"\""
         }
     }' "$POPLABEL"
)"

# Get all population labels
mapfile -t pops < <(awk 'NR>1 {seen[$2]=1} END{for(p in seen) print p}' "$POPLABEL")

mkdir -p /path/to/temp
export TMPDIR=/path/to/temp

# ================================
# Run in a loop vcf2smc
# ================================
for ((i=0; i<${#pops[@]}; i++)); do
    for k in {1..22}; do
    echo ">>> Run vcf2smc: ${pops[$i]}, chr$k"
    singularity run --bind /path/to/share9 /path/to/smcpp_latest.sif vcf2smc \
    "$vcf" "$out" \
    chr$k \
    ${pops[$i]}:${dict[${pops[$i]}]} \
    -length ${CHROM_LENGTHS[$(($k-1))]}
    echo ">>> Run estimate: ${pops[$i]}, chr$k"
    echo ${pops[$i]}
    mkdir -p ${pops[$i]}
    singularity run --bind /path/to/share9 /path/to/smcpp_latest.sif estimate \
    -o ${pops[$i]}/ \
    1.25e-8 \
    ./path/to/${pops[$i]}.chr*.smc.gz
    singularity run --bind /path/to/share9 /path/to/smcpp_latest.sif plot \
        ${pops[$i]}.plot.pdf ${pops[$i]}/model.final.json
        for ((j=i+1; j<${#pops[@]}; j++)); do
        vcf="$VCF_DIR/69Malay_20Negrito.chr${k}.QC.shapeited.vcf.gz"
        out="$OUT_DIR/${pops[$i]}.chr${k}.smc.gz"
        echo ">>> Run vcf2smc: ${pops[$i]} vs ${pops[$j]}, chr$k"
        singularity run --bind /path/to/share9 /path/to/smcpp_latest.sif vcf2smc \
        "$vcf" "$out" \
        chr$k \
        ${pops[$i]}:${dict[${pops[$i]}]} \
        ${pops[$j]}:${dict[${pops[$j]}]} \
        --length ${CHROM_LENGTHS[$(($k-1))]}

        echo ">>> Run estimate: ${pops[$i]}, ${pops[$j]}, chr$k"

        singularity run --bind /path/to/share9 /path/to/smcpp_latest.sif estimate \
        -o analysis/ \
        1.25e-8 \
        "./path/to/${pops[$i]}.${pops[$j]}.chr${k}.smc.gz"

        done
    done
done

singularity exec /path/to/smcpp_latest.sif ls -l /path/to/data


singularity run /path/to/smcpp_latest.sif vcf2smc \
$vcf data/${pops[$i]}.${pops[$j]}.chr$k.smc.gz \
chr$k \
${pops[$i]}:${dict[${pops[$i]}]} \
${pops[$j]}:${dict[${pops[$j]}]} \
--length ${CHROM_LENGTHS[$(($k-1))]}


