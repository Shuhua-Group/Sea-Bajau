#!/bin/bash



# Analyze digestive-enzyme copy numbers and haplotypes.
# Set input, reference and tool paths for your environment before execution.

for i in `cut -d "," -f1 ../PGA.ctg.csv|uniq`;do
# start=`grep ^${i} ../PGA.ctg.csv|awk -F',' '{if ($3 != 0 && $3 != 1) print $2}'|head -1`
start=`grep ^${i} ../PGA.ctg.csv|awk -F',' '{if ($3 == 1) print $2}'|tail -1`
end=`grep ^${i} ../PGA.ctg.csv|awk -F',' '{if ($3 == 0) print $2}'|head -1`
echo "Processing ${i} ${start} ${end}"
# i format: HG00097_hap1#HG00097#1#CM094067.1_59988956_60179168_0
# Keep the first # before
sample=`echo ${i} | cut -d "#" -f1`
samtools faidx ../PGA_region_seq.fa.gz "${i}:${start}-${end}" > "${sample}.fa"
samtools faidx "${sample}.fa"
done


# According to PGA Group haplotypes
mkdir groups
while read sample hap; do
  hap_mod=$(echo "$hap" \
    | sed -e 's/^1,//' \
          -e 's/,0$//' \
          -e 's/,/_/g' \
          -e 's/PGA//g')
# echo ${sample} ${hap_mod}
  dir="groups/${hap_mod}"
  mkdir -p "$dir"
  cp "${sample}.fa" "$dir/"
done < ../PGA.GT.out


# conda activate msa_env
cd groups || exit 1

# Inspect only group groups containing more than one sample
for i in `ls`;do
  if [ $(ls $i | wc -l) -gt 1 ];then
    echo "$i" $(ls $i | wc -l)
    # Step 1: Check before within-group alignment (If within-group divergence is too high, Further splitting may be needed, Check the grouping)
    mash sketch -o ${i}.msh $i/*.fa
    # mash dist ${i}.msh ${i}.msh > ${i}.mashdist
    Extract highly divergent sample pairs
    file="${i}.mashdist"
    threshold="0.01" # Custom threshold, Flag differences above this threshold
    awk -v thr="$threshold" '
    {
    # Extract group, Assume group name occurs in / first filename element after _ before
    split($1, a, "/"); split(a[2], a1, "_"); g1=a1[1]
    split($2, b, "/"); split(b[2], b1, "_"); g2=b1[1]
    if ($3 > thr) {
    print g1, $1, $2, $3
    }
    }' OFS="\t" "$file" | sort -k1,1 -k4,4nr > ${i}.diff
    Step 2: Align within groups
    cat $i/*.fa > ${i}.merge.temp.fa
    mafft --auto ${i}.merge.temp.fa > ${i}.aln.fasta
    # Step 3: Use trimal Remove all- gap Column
    trimal -automated1 -in ${i}.aln.fasta -out ${i}.aln.trim.fasta
    # Step 4: Build a phylogeny
    iqtree -s ${i}.aln.trim.fasta -m TEST -bb 1000 -alrt 1000 -nt AUTO -pre "$i"_iqtree
  fi
done

# 3 5
# 3_4 14
# 3_4_3_5 4
# 3_4_4 5
# 3_4_4_4 2
# 3_4_4_4_4 2
# 3_4_4_4_4_5 17
# 3_4_4_4_5 42
# 3_4_4_5 116
# 3_4_5 69
# 3_5 78
# 4_5 5
## Optionally inspect the most common class, 3_4_4_5 ? This interpretation is tentative


# i="3_4_4_4_4_5"
# trimal -automated1 -in ${i}.aln.fasta -out ${i}.aln.trim.fasta
# # Step 4: Build a phylogeny
# iqtree -s ${i}.aln.trim.fasta -m TEST -bb 1000 -alrt 1000 -nt AUTO -pre "$i"_iqtree



# Inspect only group groups containing more than one sample
for group_dir in */; do
  group_name=${group_dir%/}  #  Remove the trailing slash
  sample_count=$(ls "$group_dir" | wc -l)

  if [ "$sample_count" -gt 1 ]; then
    echo "Processing group: $group_name ($sample_count samples)"

    # ===== 1. Screen within-group genetic distances =====
    mash sketch -o "${group_name}.msh" "$group_dir"*.fa
    mash dist "${group_name}.msh" "${group_name}.msh" > "${group_name}.mashdist"

    # # Apply a stricter threshold (0.005=0.5% Difference)
    threshold="0.005"
    awk -v thr="$threshold" '
    {
      split($1, a, "/");
      sample1 = a[length(a)];
      split($2, b, "/");
      sample2 = b[length(b)];

      if ($3 > thr && sample1 != sample2) {
        print sample1, sample2, $3
      }
    }' "${group_name}.mashdist" | sort -k3nr > "${group_name}.high_diff_pairs.tsv"

    # ===== 2. Refine the multiple sequence alignment =====
    cat "$group_dir"*.fa > "${group_name}.merged.fa"

    # Use a more accurate L-INS-i Algorithm (Suitable for the conserved core region + Flanking variation)
    mafft --localpair --maxiterate 1000 --thread 8 \
    --leavegappyregion "${group_name}.merged.fa" > "${group_name}.aln.fasta"

    # ===== 3. Refine alignment trimming =====
    # Apply two-stage trimming
    trimal -in "${group_name}.aln.fasta" -out "${group_name}.aln.trim1.fasta" -gt 0.9 # Keep 90% Sites with called bases in the sample
    trimal -in "${group_name}.aln.trim1.fasta" -out "${group_name}.aln.trim2.fasta" -resoverlap 0.75 -seqoverlap 80

    # ===== 4. Build the phylogenetic tree =====
    iqtree2 -s "${group_name}.aln.trim2.fasta" \
    -m MFP+MERGE -B 1000 --alrt 1000 -T AUTO \
    --prefix "${group_name}_tree" \
    --ufboot

    # ===== 5. Ancestral state reconstruction ===== (Optional)
    # Extract copy-number labels (Expected filename format: sample_haplotype_copies.fa)
    for f in "$group_dir"*.fa; do
    sample=$(basename "$f" | cut -d'_' -f1)
    copies=$(basename "$f" | cut -d'_' -f3 | cut -d'.' -f1)
    echo -e "$sample\t$copies"
    done > "${group_name}.copies.tsv"

    iqtree2 -s "${group_name}.aln.trim2.fasta" \
    -m MFP -asr \
    --label-test "${group_name}.copies.tsv" \
    --prefix "${group_name}_ancestral"

  fi
done



group_name="3_4_4_4_5"
group_dir="3_4_4_4_5/"
echo "Processing group: $group_name ($sample_count samples)"

# ===== 1. Screen within-group genetic distances =====
mash sketch -o "${group_name}.msh" "$group_dir"*.fa
mash dist "${group_name}.msh" "${group_name}.msh" > "${group_name}.mashdist"

# # Apply a stricter threshold (0.005=0.5% Difference)
threshold="0.005"
awk -v thr="$threshold" '
{
  split($1, a, "/");
  sample1 = a[length(a)];
  split($2, b, "/");
  sample2 = b[length(b)];

  if ($3 > thr && sample1 != sample2) {
    print sample1, sample2, $3
  }
}' "${group_name}.mashdist" | sort -k3nr > "${group_name}.high_diff_pairs.tsv"

# ===== 2. Refine the multiple sequence alignment =====
cat "$group_dir"*.fa > "${group_name}.merged.fa"

# # Use a more accurate L-INS-i Algorithm (Suitable for the conserved core region + Flanking variation)
mafft --localpair --maxiterate 1000 --thread 8 \
      --leavegappyregion "${group_name}.merged.fa" > "${group_name}.aln.fasta"

# # ===== 3. Refine alignment trimming =====
# # Apply two-stage trimming
trimal -in "${group_name}.aln.fasta" -out "${group_name}.aln.trim1.fasta" -gt 0.9  #  Keep 90% Sites with called bases in the sample
trimal -in "${group_name}.aln.trim1.fasta" -out "${group_name}.aln.trim2.fasta" -resoverlap 0.75 -seqoverlap 80

# # ===== 4. Build the phylogenetic tree =====
iqtree -s "${group_name}.aln.trim2.fasta"  -m TEST -bb 1000 -alrt 1000 -nt AUTO -pre "$group_name"_iqtree

# # ===== 5. Ancestral state reconstruction ===== (Optional)
# Extract copy-number labels (Expected filename format: sample_haplotype_copies.fa)
for f in "$group_dir"*.fa; do
sample=$(basename "$f" | cut -d'_' -f1)
copies=$(basename "$f" | cut -d'_' -f3 | cut -d'.' -f1)
echo -e "$sample\t$copies"
done > "${group_name}.copies.tsv"

iqtree2 -s "${group_name}.aln.trim2.fasta" \
-m MFP -asr \
--label-test "${group_name}.copies.tsv" \
--prefix "${group_name}_ancestral"
