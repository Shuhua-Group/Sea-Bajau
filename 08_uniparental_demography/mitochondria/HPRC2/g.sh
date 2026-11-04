#!/usr/bin/env bash
#@! /bin/bash

# Analyze mitochondrial and Y-chromosome lineage variation.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/HPRC2"

cd "$wd" || exit 1
for file in `ls /path/to/*.chm13v2.minimap.paf`
do
sample=$(basename "$file" .chm13v2.minimap.paf)
awk '{if ($6 == "chrM") print $1}' "$file" >> "${sample}.contigs.list"
sort "${sample}.contigs.list"|uniq > "${sample}.contigs.uniq.list"
done

rm *.contigs.list
# Then select by size =0 remove these files
for file in `ls *.contigs.uniq.list`
do
if [ -s "$file" ];then
    echo "$file"
else
    rm "$file"
fi
done

# Then extract mtDNA Sequence
asmdir="/path/to/HPRC2"

# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r i <&9 || [[ -n "$i" ]]; do
    [[ -z "$i" ]] && continue
sample=$i

ori_hap1=`grep -w "$i" /path/to/HPRC2.hap.list| head -1|cut -f2`
ori_hap2=`grep -w "$i" /path/to/HPRC2.hap.list| tail -1|cut -f2`
hap1=$(basename $ori_hap1 .gz)
hap2=$(basename $ori_hap2 .gz)
asm1="$asmdir/$hap1"
asm2="$asmdir/$hap2"

if [ -f "${sample}.hap1.contigs.uniq.list" ];then
    # Read one record per line; keep backslashes and a final line without a newline.
    while IFS= read -r contig <&9 || [[ -n "$contig" ]]; do
        [[ -z "$contig" ]] && continue
        samtools faidx "$asm1" "$contig" >> "${sample}.hap1.mtDNA.fasta"
    done 9< "${sample}.hap1.contigs.uniq.list"
fi

if [ -f "${sample}.hap2.contigs.uniq.list" ];then
    # Read one record per line; keep backslashes and a final line without a newline.
    while IFS= read -r contig <&9 || [[ -n "$contig" ]]; do
        [[ -z "$contig" ]] && continue
        samtools faidx "$asm2" "$contig" >> "${sample}.hap2.mtDNA.fasta"
    done 9< "${sample}.hap2.contigs.uniq.list"
fi
done 9< "/path/to/HPRC2.list"

for fasta in `ls *.mtDNA.fasta`
do
    samtools faidx "$fasta"
done

# Found fai with lengths greater than 20000 sequences, report Identify the samples
for fai in `ls *.mtDNA.fasta.fai`
do
    contig=`awk '{if ($2 > 20000) print $1}' $fai`
    if [ -n "$contig" ];then
        echo -e "$fai" >> mtDNA.len20k.list
    else
        continue
    fi
done

# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r fai <&9 || [[ -n "$fai" ]]; do
    [[ -z "$fai" ]] && continue
# awk '{if ($2 < 20000) print $1}' $fai > ${fai%.mtDNA.fasta.fai}.mtDNA.len20k.list
# for contig in `cat ${fai%.mtDNA.fasta.fai}.mtDNA.len20k.list`
# do
# samtools faidx ${fai%.mtDNA.fasta.fai}.mtDNA.fasta $contig >> ${fai%.mtDNA.fasta.fai}.re.mtDNA.fasta
# done
rm ${fai%.mtDNA.fasta.fai}.mtDNA.fasta
mv ${fai%.mtDNA.fasta.fai}.re.mtDNA.fasta ${fai%.mtDNA.fasta.fai}.mtDNA.fasta
done 9< "mtDNA.len20k.list"

# Confirm coverage for all mtDNA.fasta are present
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r i <&9 || [[ -n "$i" ]]; do
    [[ -z "$i" ]] && continue
sample=$i
fasta=$(ls ${sample}.*mtDNA.fasta)
if [ ! -f "$fasta" ];then
echo "$sample mtDNA fasta missing!"
fi
done 9< "/path/to/HPRC2.list"
