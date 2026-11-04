#!/bin/bash
# Process HiFi reads and summarize assembly quality.
# Set input, reference and tool paths for your environment before execution.

samplelist="/path/to/sample.list"
prefix="chm13v2"
ref="/path/to/chm13v2.0.fa"
minimap2_cmd="/path/to/minimap2"
thread=16
biotools_env="/path/to/biotools_env"
replace_fasta_cmd="/path/to/02.replace_fasta.py"
wd="/path/to/v2_align"
## run alignment



# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
outdir=${wd}/../v2_corstrand/${prefix}/${sample}
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
$outdir/${sample}.hap1.${prefix}.negative_strand_correction.fa > ${wd}/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sam
samtools view -bS -@${thread} ${wd}/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sam -o ${wd}/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.bam
samtools sort -@${thread} ${wd}/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.bam -o ${wd}/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sorted.bam
samtools index ${wd}/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sorted.bam

/usr/bin/time -v ${minimap2_cmd} -ax asm20 --cs -r2k -t ${thread} \
${ref} \
$outdir/${sample}.hap2.${prefix}.negative_strand_correction.fa > ${wd}/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sam
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
