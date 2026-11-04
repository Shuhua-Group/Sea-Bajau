#!/bin/bash


# Analyze MUC5AC sequences and structural haplotypes.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/plot"
mkdir -p "$wd"
cd "$wd" || exit 1


# Read all fasta Read files into an array
mapfile -t seq < <(ls *.connected.fa)

for fa in "${seq[@]}"; do
    samtools faidx "$fa"
    # Calculate the fasta sum of sequence lengths (Take .fai Second column)
    len=$(awk '{sum+=$2} END{print sum}' "${fa}.fai")
    echo -e "$fa\t$len" >> seq_length.txt
done

# Take the left part of each filename 50kb Trim, Keep the right part


# Sort by length (Descending order), Regenerate the ordered file list
# mapfile -t seq_sorted < <(sort -k2,2nr left_cut50k/seq_length.txt | cut -f1)
# Set explicitly seq_sorted
seq_sorted=(
    "muc5ac.chm13v2.connected.fa"
    "NB2388.hap1.muc5ac.connected.fa"
)

# Align consecutive pairs
for ((i=0; i<${#seq_sorted[@]}-1; i++)); do
    prefix1=${seq_sorted[i]%.GSTM4up20k_GSTM3down5k.fa}
    prefix2=${seq_sorted[i+1]%.GSTM4up20k_GSTM3down5k.fa}
    blastn -query "${seq_sorted[i]}" \
           -subject "${seq_sorted[i+1]}" \
           -outfmt 6 \
           -out "${prefix1}.${prefix2}.blast.out"
done

# Clean up
rm seq_length.txt

cat *.blast.out |awk '{if($3 > 0.88 && $4 > 1000) print $0}' > MUC5AC.blast.out
