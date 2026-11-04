#!/bin/bash
# Process HiFi reads and summarize assembly quality.
# Set input, reference and tool paths for your environment before execution.

samplelist="/path/to/sample.list"
femalelist="/path/to/female.list"
malelist="/path/to/male.list"
prefix="chm13v2"
ref="/path/to/chm13v2.0.fa"
minimap2_cmd="/path/to/minimap2"
thread=16
biotools_env="/path/to/biotools_env"
replace_fasta_cmd="/path/to/02.replace_fasta.py"
wd="/path/to/norDNA_virus.align"
## run alignment


# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
mkdir -p "${wd}/${prefix}/${sample}/log"
asm_align_script="#!/bin/bash
#SBATCH --job-name=m_$sample
#SBATCH --partition=batch
#SBATCH --output=${wd}/${prefix}/${sample}/log/${sample}.minimap2.stdout
#SBATCH --error=${wd}/${prefix}/${sample}/log/${sample}.minimap2.stderr
#SBATCH --mincpus=${thread}
#SBATCH --mem=48g
/usr/bin/time -v ${minimap2_cmd} -ax asm20 --cs -r2k -t ${thread} \
${ref} \
/path/to/${sample}.1.9606.clean.fasta > ${wd}/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sam
samtools view -bS -@${thread} ${wd}/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sam -o ${wd}/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.bam
samtools sort -@${thread} ${wd}/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.bam -o ${wd}/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sorted.bam
samtools index ${wd}/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sorted.bam
/usr/bin/time -v ${minimap2_cmd} -ax asm20 --cs -r2k -t ${thread} \
${ref} \
/path/to/${sample}.2.9606.clean.fasta > ${wd}/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sam
samtools view -bS -@${thread} ${wd}/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sam -o ${wd}/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.bam
samtools sort -@${thread} ${wd}/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.bam -o ${wd}/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sorted.bam
samtools index ${wd}/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sorted.bam
rm ${wd}/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sam
rm ${wd}/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sam
"
echo -e "$asm_align_script" > "${wd}/${prefix}/${sample}/run.align.cmd.sh"
echo -e "sh ${wd}/${prefix}/${sample}/run.align.cmd.sh" >> "$wd/asm_align.${prefix}.run.sh"
echo -e "sbatch ${wd}/${prefix}/${sample}/run.align.cmd.sh" >> "$wd/asm_align.${prefix}.sbatch.sh"
done 9< "${samplelist}"

# Female assignment
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
# 1.1 run assign
#The script is used to assign assembled contigs longer than 100kb to chromosomes (excluding chrM).
# For male samples, all chrX contigs will be assigned to hap1, and all chrY contigs will be assigned to hap2.
mkdir -p "${wd}/../v1_assign/${prefix}/${sample}/log"
asm1=/path/to/${sample}.1.9606.clean.fasta
asm2=/path/to/${sample}.2.9606.clean.fasta
bam1=${wd}/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sorted.bam
bam2=${wd}/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sorted.bam
assign_script="#!bin/bash
#SBATCH --job-name=${sample}_asm_chrom_assigned
#SBATCH --partition=batch
#SBATCH --output=${wd}/../v1_assign/${prefix}/${sample}/log/asm_chrom_assigned.log
#SBATCH --error=${wd}/../v1_assign/${prefix}/${sample}/log/asm_chrom_assigned.err
#SBATCH --mincpus=4
#SBATCH --mem=16g
cd ${wd}/../v1_assign/${prefix}/${sample}
source activate $biotools_env
## Hap1
# seqkit seq -M 100000 $asm1 -j 8 > ${sample}.asm.hap1.unassigned.fa
# seqkit seq -m 100000 $asm1 -j 8 > ${sample}.asm.hap1.m100kb.fa
ln -s $asm1 ${sample}.asm.hap1.m100kb.fa
grep '^>' ${sample}.asm.hap1.m100kb.fa | sed 's/^>//' > ${sample}.hap1.m100kb_contigs.ids
> ${sample}.asm.hap1.fa
> ${sample}.asm.hap1.id_mapping.txt
for i in {1..22} X; do
samtools view -b -F 256 -@8 $bam1 chr\${i} > ${sample}.hap1.chr.bam
samtools fasta ${sample}.hap1.chr.bam > ${sample}.hap1.chr.fa
# Select sequences longer than 100kb Of contigs
# seqkit seq -m 100000 ${sample}.hap1.chr.fa -j 8 > ${sample}.hap1.chr.m100kb.fa
ln -sf ${sample}.hap1.chr.fa ${sample}.hap1.chr.m100kb.fa
# From the original assembly: contigs ID Replace with h1-chr1-1 pattern
grep '^>' ${sample}.hap1.chr.m100kb.fa | sed 's/^>//' > ids.tmp
awk -v i=\${i} 'BEGIN {OFS=\"\t\"} {print \$1, \"h1-chr\" i \"-\" NR}' ids.tmp > ids_mapping.tmp
# Merge temporary results into the final results file
seqkit replace -p '^(\S+)' -r '{kv}\$2' -k ids_mapping.tmp ${sample}.hap1.chr.m100kb.fa -j 8 > ${sample}.hap1.chr.renamed.fa
cat ${sample}.asm.hap1.fa ${sample}.hap1.chr.renamed.fa > ${sample}.asm.hap1.tmp.fa
mv ${sample}.asm.hap1.tmp.fa ${sample}.asm.hap1.fa
cat ${sample}.asm.hap1.id_mapping.txt ids_mapping.tmp > ${sample}.asm.hap1.id_mapping.tmp.txt
mv ${sample}.asm.hap1.id_mapping.tmp.txt ${sample}.asm.hap1.id_mapping.txt
done
# Check unaligned sequences (flag=4) Of contigs
samtools view -b -f 4 -@8 $bam1 > ${sample}.asm.hap1.unaligned.bam
samtools fasta -f 4 ${sample}.asm.hap1.unaligned.bam > ${sample}.asm.hap1.unaligned.fa
if [ ! -s ${sample}.hap1.unaligned.fa ]; then
rm ${sample}.asm.hap1.unaligned.fa ${sample}.asm.hap1.unaligned.bam
fi
rm ids.tmp ids_mapping.tmp ${sample}.hap1.chr.bam ${sample}.hap1.chr.fa ${sample}.hap1.chr.renamed.fa ${sample}.hap1.chr.m100kb.fa ${sample}.asm.hap1.m100kb.fa
samtools faidx ${sample}.asm.hap1.fa

## Hap2
# seqkit seq -M 100000 $asm2 -j 8 > ${sample}.asm.hap2.unassigned.fa
# seqkit seq -m 100000 $asm2 -j 8 > ${sample}.asm.hap2.m100kb.fa
ln -s $asm2 ${sample}.asm.hap2.m100kb.fa
grep '^>' ${sample}.asm.hap2.m100kb.fa | sed 's/^>//' > ${sample}.hap2.m100kb_contigs.ids
> ${sample}.asm.hap2.fa
> ${sample}.asm.hap2.id_mapping.txt
for i in {1..22} X; do
samtools view -b -F 256 -@8 $bam2 chr\${i} > ${sample}.hap2.chr.bam
samtools fasta ${sample}.hap2.chr.bam > ${sample}.hap2.chr.fa
# Select sequences longer than 100kb Of contigs
# seqkit seq -m 100000 ${sample}.hap2.chr.fa -j 8 > ${sample}.hap2.chr.m100kb.fa
ln -sf ${sample}.hap2.chr.fa ${sample}.hap2.chr.m100kb.fa
# From the original assembly: contigs ID Replace with h1-chr1-1 pattern
grep '^>' ${sample}.hap2.chr.m100kb.fa | sed 's/^>//' > ids.tmp
awk -v i=\${i} 'BEGIN {OFS=\"\t\"} {print \$1, \"h2-chr\" i \"-\" NR}' ids.tmp > ids_mapping.tmp
# Merge temporary results into the final results file
seqkit replace -p '^(\S+)' -r '{kv}\$2' -k ids_mapping.tmp ${sample}.hap2.chr.m100kb.fa -j 8 > ${sample}.hap2.chr.renamed.fa
cat ${sample}.asm.hap2.fa ${sample}.hap2.chr.renamed.fa > ${sample}.asm.hap2.tmp.fa
mv ${sample}.asm.hap2.tmp.fa ${sample}.asm.hap2.fa
cat ${sample}.asm.hap2.id_mapping.txt ids_mapping.tmp > ${sample}.asm.hap2.id_mapping.tmp.txt
mv ${sample}.asm.hap2.id_mapping.tmp.txt ${sample}.asm.hap2.id_mapping.txt
done
# Check unaligned sequences (flag=4) Of contigs
samtools view -b -f 4 -@8 $bam2 > ${sample}.asm.hap2.unaligned.bam
samtools fasta -f 4 ${sample}.asm.hap2.unaligned.bam > ${sample}.asm.hap2.unaligned.fa
if [ ! -s ${sample}.hap2.unaligned.fa ]; then
rm ${sample}.asm.hap2.unaligned.fa ${sample}.asm.hap2.unaligned.bam
fi
rm ids.tmp ids_mapping.tmp ${sample}.hap2.chr.bam ${sample}.hap2.chr.fa ${sample}.hap2.chr.renamed.fa ${sample}.hap2.chr.m100kb.fa ${sample}.asm.hap2.m100kb.fa
samtools faidx ${sample}.asm.hap2.fa
cat ${sample}.asm.hap1.id_mapping.txt ${sample}.asm.hap2.id_mapping.txt > ${sample}.asm.id_mapping.txt
"
echo -e "${assign_script}" > "${wd}/../v1_assign/${prefix}/${sample}/run.assign.${prefix}.cmd.sh"
echo -e "sh ${wd}/../v1_assign/${prefix}/${sample}/run.assign.${prefix}.cmd.sh" >> /path/to/assign.run.sh
echo -e "sbatch ${wd}/../v1_assign/${prefix}/${sample}/run.assign.${prefix}.cmd.sh" >> /path/to/assign.sbatch.sh
done 9< "${femalelist}"

# Male assignment
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
mkdir -p "${wd}/../v1_assign/${prefix}/${sample}/log"
# Original and aligned assemblies: bam
asm1=/path/to/${sample}.1.9606.clean.fasta
asm2=/path/to/${sample}.2.9606.clean.fasta
bam1=${wd}/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sorted.bam
bam2=${wd}/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sorted.bam

assign_script="#!/bin/bash
#SBATCH --job-name=${sample}_asm_chrom_assigned
#SBATCH --partition=batch
#SBATCH --output=${wd}/../v1_assign/${prefix}/${sample}/log/asm_chrom_assigned.log
#SBATCH --error=${wd}/../v1_assign/${prefix}/${sample}/log/asm_chrom_assigned.err
#SBATCH --mincpus=4
#SBATCH --mem=16g

cd ${wd}/../v1_assign/${prefix}/${sample}
source activate $biotools_env
# chrXY
chrX_prefix=h1
chrY_prefix=h2
for i in X Y; do
# Merge two sequences hap in chrX/Y Of contigs
samtools view -b -F 256 -@8 $bam1 chr\${i} > ${sample}.hap1.chr\${i}.bam
samtools view -b -F 256 -@8 $bam2 chr\${i} > ${sample}.hap2.chr\${i}.bam
samtools merge -@8 ${sample}.asm.chr\${i}.bam ${sample}.hap1.chr\${i}.bam ${sample}.hap2.chr\${i}.bam
samtools sort -@8 -o ${sample}.asm.chr\${i}.sorted.bam ${sample}.asm.chr\${i}.bam
samtools fasta ${sample}.asm.chr\${i}.sorted.bam > ${sample}.asm.chr\${i}.fa
# Select sequences longer than 100kb Of contigs
# seqkit seq -m 100000 ${sample}.asm.chr\${i}.fa -j 8 > ${sample}.asm.chr\${i}.m100kb.fa
ln -s ${sample}.asm.chr\${i}.fa ${sample}.asm.chr\${i}.m100kb.fa
# From the original assembly: contigs ID Replace with h1-chr1-1 pattern
grep '^>' ${sample}.asm.chr\${i}.m100kb.fa | sed 's/^>//' > ids.tmp
hap_prefix_var=\"chr\${i}_prefix\"
hap_prefix=\${!hap_prefix_var}
awk -v i=\${i} -v hap_prefix=\"\${hap_prefix}\" 'BEGIN {OFS=\"\t\"} {print \$1, hap_prefix \"-chr\" i \"-\" NR}' ids.tmp > ${sample}.asm.chr\${i}.ids_mapping.txt
seqkit replace -p '^(\S+)' -r '{kv}\$2' -k ${sample}.asm.chr\${i}.ids_mapping.txt ${sample}.asm.chr\${i}.m100kb.fa -j 8 > ${sample}.asm.chr\${i}.renamed.fa
done
rm ${sample}.asm.chrX.bam ${sample}.asm.chrX.fa ${sample}.asm.chrX.m100kb.fa ${sample}.asm.chrX.sorted.bam ${sample}.asm.chrY.bam ${sample}.asm.chrY.fa ${sample}.asm.chrY.m100kb.fa ${sample}.asm.chrY.sorted.bam ${sample}.hap1.chrX.bam ${sample}.hap2.chrX.bam ${sample}.hap1.chrY.bam ${sample}.hap2.chrY.bam
# asm1
# seqkit seq -M 100000 $asm1 -j 8 > ${sample}.asm.hap1.unassigned.fa
# seqkit seq -m 100000 $asm1 -j 8 > ${sample}.asm.hap1.m100kb.fa
ln -s $asm1 ${sample}.asm.hap1.m100kb.fa
grep '^>' ${sample}.asm.hap1.m100kb.fa | sed 's/^>//' > ${sample}.hap1.m100kb_contigs.ids
> ${sample}.asm.hap1.fa
> ${sample}.asm.hap1.id_mapping.txt
for i in {1..22}; do
samtools view -b -F 256 -@8 $bam1 chr\${i} > ${sample}.hap1.chr.bam
samtools fasta ${sample}.hap1.chr.bam > ${sample}.hap1.chr.fa
# Select sequences longer than 100kb Of contigs
# seqkit seq -m 100000 ${sample}.hap1.chr.fa -j 8 > ${sample}.hap1.chr.m100kb.fa
ln -s ${sample}.hap1.chr.fa ${sample}.hap1.chr.m100kb.fa
# From the original assembly: contigs ID Replace with h1-chr1-1 pattern
grep '^>' ${sample}.hap1.chr.m100kb.fa | sed 's/^>//' > ids.tmp
awk -v i=\${i} 'BEGIN {OFS=\"\t\"} {print \$1, \"h1-chr\" i \"-\" NR}' ids.tmp > ids_mapping.tmp
seqkit replace -p '^(\S+)' -r '{kv}\$2' -k ids_mapping.tmp ${sample}.hap1.chr.m100kb.fa -j 8 > ${sample}.hap1.chr.renamed.fa
# Merge temporary results into the final results file
cat ${sample}.asm.hap1.fa ${sample}.hap1.chr.renamed.fa > ${sample}.asm.hap1.tmp.fa
mv ${sample}.asm.hap1.tmp.fa ${sample}.asm.hap1.fa
cat ${sample}.asm.hap1.id_mapping.txt ids_mapping.tmp > ${sample}.asm.hap1.id_mapping.tmp.txt
mv ${sample}.asm.hap1.id_mapping.tmp.txt ${sample}.asm.hap1.id_mapping.txt
done
# Check unaligned sequences (flag=4) Of contigs
samtools view -b -f 4 -@8 $bam1 > ${sample}.asm.hap1.unaligned.bam
samtools fasta -f 4 ${sample}.asm.hap1.unaligned.bam > ${sample}.asm.hap1.unaligned.fa
if [ ! -s ${sample}.hap1.unaligned.fa ]; then
rm ${sample}.asm.hap1.unaligned.fa ${sample}.asm.hap1.unaligned.bam
fi
rm ids.tmp ids_mapping.tmp ${sample}.hap1.chr.bam ${sample}.hap1.chr.fa ${sample}.hap1.chr.fa ${sample}.hap1.chr.renamed.fa ${sample}.hap1.chr.m100kb.fa ${sample}.asm.hap1.m100kb.fa
# Merge chrX contigs To hap1 Assembly
cat ${sample}.asm.hap1.fa ${sample}.asm.chrX.renamed.fa > ${sample}.asm.hap1.final.fa
mv ${sample}.asm.hap1.final.fa ${sample}.asm.hap1.fa
samtools faidx ${sample}.asm.hap1.fa
# asm2
# seqkit seq -M 100000 $asm2 -j 8 > ${sample}.asm.hap2.unassigned.fa
# seqkit seq -m 100000 $asm2 -j 8 > ${sample}.asm.hap2.m100kb.fa
ln -s $asm2 ${sample}.asm.hap2.m100kb.fa
grep '^>' ${sample}.asm.hap2.m100kb.fa | sed 's/^>//' > ${sample}.hap2.m100kb_contigs.ids
> ${sample}.asm.hap2.fa
> ${sample}.asm.hap2.id_mapping.txt
for i in {1..22}; do
samtools view -b -F 256 -@8 $bam2 chr\${i} > ${sample}.hap2.chr.bam
samtools fasta ${sample}.hap2.chr.bam > ${sample}.hap2.chr.fa
# Select sequences longer than 100kb Of contigs
# seqkit seq -m 100000 ${sample}.hap2.chr.fa -j 8 > ${sample}.hap2.chr.m100kb.fa
ln -s ${sample}.hap2.chr.fa ${sample}.hap2.chr.m100kb.fa
# From the original assembly: contigs ID Replace with h1-chr1-1 pattern
grep '^>' ${sample}.hap2.chr.m100kb.fa | sed 's/^>//' > ids.tmp
awk -v i=\${i} 'BEGIN {OFS=\"\t\"} {print \$1, \"h2-chr\" i \"-\" NR}' ids.tmp > ids_mapping.tmp
seqkit replace -p '^(\S+)' -r '{kv}\$2' -k ids_mapping.tmp ${sample}.hap2.chr.m100kb.fa -j 8 > ${sample}.hap2.chr.renamed.fa
# Merge temporary results into the final results file
cat ${sample}.asm.hap2.fa ${sample}.hap2.chr.renamed.fa > ${sample}.asm.hap2.tmp.fa
mv ${sample}.asm.hap2.tmp.fa ${sample}.asm.hap2.fa
cat ${sample}.asm.hap2.id_mapping.txt ids_mapping.tmp > ${sample}.asm.hap2.id_mapping.tmp.txt
mv ${sample}.asm.hap2.id_mapping.tmp.txt ${sample}.asm.hap2.id_mapping.txt
done
# Check unaligned sequences (flag=4) Of contigs
samtools view -b -f 4 -@8 $bam2 > ${sample}.asm.hap2.unaligned.bam
samtools fasta -f 4 ${sample}.asm.hap2.unaligned.bam > ${sample}.asm.hap2.unaligned.fa
if [ ! -s ${sample}.hap2.unaligned.fa ]; then
rm ${sample}.asm.hap2.unaligned.fa ${sample}.asm.hap2.unaligned.bam
fi
rm ids.tmp ids_mapping.tmp ${sample}.hap2.chr.bam ${sample}.hap2.chr.fa ${sample}.hap2.chr.renamed.fa ${sample}.hap2.chr.m100kb.fa ${sample}.asm.hap2.m100kb.fa
# Merge chrY contigs To hap2 Assembly
cat ${sample}.asm.hap2.fa ${sample}.asm.chrY.renamed.fa > ${sample}.asm.hap2.final.fa
mv ${sample}.asm.hap2.final.fa ${sample}.asm.hap2.fa
samtools faidx ${sample}.asm.hap2.fa
rm ${sample}.asm.chrX.renamed.fa ${sample}.asm.chrY.renamed.fa
cat ${sample}.asm.hap1.id_mapping.txt ${sample}.asm.hap2.id_mapping.txt ${sample}.asm.chrX.ids_mapping.txt ${sample}.asm.chrY.ids_mapping.txt > ${sample}.asm.id_mapping.txt
"
echo -e "${assign_script}" > "${wd}/../v1_assign/${prefix}/${sample}/run.assign.${prefix}.cmd.sh"
echo -e "sh ${wd}/../v1_assign/${prefix}/${sample}/run.assign.${prefix}.cmd.sh" >> /path/to/assign.run.sh
echo -e "sbatch ${wd}/../v1_assign/${prefix}/${sample}/run.assign.${prefix}.cmd.sh" >> /path/to/assign.sbatch.sh
done 9< "$malelist"

#1.2 run correct strand

# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
# Align original assemblies bam
bam1=${wd}/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sorted.bam
bam2=${wd}/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sorted.bam
outdir=${wd}/../v2_corstrand/${prefix}/${sample}
mkdir -p "${outdir}/log"
corstrand_script="#!/bin/bash
#SBATCH --job-name=${sample}.negative_strand_correction
#SBATCH --partition=batch
#SBATCH --output=${outdir}/log/negative_strand_correction.log
#SBATCH --error=${outdir}/log/negative_strand_correction.err
#SBATCH --mincpus=8
#SBATCH --mem=16g
# Extract two sequences from the original assembly hap All reverse-strand alignments: contigs Merge into a single fa In
samtools view -b -f 16 -@8 $bam1 > $outdir/${sample}.hap1.${prefix}.minimap.f16.sorted.bam
samtools fasta $outdir/${sample}.hap1.${prefix}.minimap.f16.sorted.bam > $outdir/${sample}.hap1.${prefix}.minimap.f16.fa
samtools view -b -f 16 -@8 $bam2 > $outdir/${sample}.hap2.${prefix}.minimap.f16.sorted.bam
samtools fasta $outdir/${sample}.hap2.${prefix}.minimap.f16.sorted.bam > $outdir/${sample}.hap2.${prefix}.minimap.f16.fa
cat $outdir/${sample}.hap1.${prefix}.minimap.f16.fa $outdir/${sample}.hap2.${prefix}.minimap.f16.fa >$outdir/${sample}.asm.${prefix}.minimap.f16.fa
# seqtk Generate reverse-complement sequences
seqtk seq -r $outdir/${sample}.asm.${prefix}.minimap.f16.fa > $outdir/${sample}.asm.${prefix}.minimap.f16.rc.fa
# Based on 01.chrom_assign Generated id_mapping.txt, Replace seq ID As h1-chr1-1 pattern, Less than 100kb Of contigs Keep the original ID
source activate $biotools_env
seqkit replace -p '^(\S+)' -r '{kv}\$2' -k ${wd}/../v1_assign/${prefix}/${sample}/${sample}.asm.id_mapping.txt $outdir/${sample}.asm.chm13v2.minimap.f16.rc.fa --keep-key -j 8 > $outdir/${sample}.asm.chm13v2.minimap.f16.rc.renamed.fa
# Use chrom_assign Reverse-strand sequences in the updated assembly: contigs, Replace with the reverse complement
python $replace_fasta_cmd $outdir/${sample}.asm.${prefix}.minimap.f16.rc.renamed.fa ${wd}/../v1_assign/${prefix}/${sample}/${sample}.asm.hap1.fa $outdir/${sample}.hap1.${prefix}.negative_strand_correction.fa
python $replace_fasta_cmd $outdir/${sample}.asm.${prefix}.minimap.f16.rc.renamed.fa ${wd}/../v1_assign/${prefix}/${sample}/${sample}.asm.hap2.fa $outdir/${sample}.hap2.${prefix}.negative_strand_correction.fa
samtools faidx $outdir/${sample}.hap1.${prefix}.negative_strand_correction.fa
samtools faidx $outdir/${sample}.hap2.${prefix}.negative_strand_correction.fa
rm $outdir/${sample}.asm.${prefix}.minimap.f16.rc.renamed.fa
rm $outdir/${sample}.asm.${prefix}.minimap.f16.rc.fa
rm $outdir/${sample}.asm.${prefix}.minimap.f16.fa
rm $outdir/${sample}.hap1.${prefix}.minimap.f16.fa $outdir/${sample}.hap2.${prefix}.minimap.f16.fa
rm $outdir/${sample}.hap1.${prefix}.minimap.f16.sorted.bam $outdir/${sample}.hap2.${prefix}.minimap.f16.sorted.bam
"
echo -e "${corstrand_script}" > "${outdir}/run.negative_strand_correction.cmd.sh"
echo -e "sh ${outdir}/run.negative_strand_correction.cmd.sh" >> /path/to/negative_strand_correction.run.sh
echo -e "sbatch ${outdir}/run.negative_strand_correction.cmd.sh" >> /path/to/negative_strand_correction.sbatch.sh
done 9< "$samplelist"
