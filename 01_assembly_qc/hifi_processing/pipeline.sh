#!/bin/bash


############################ Input settings to update ################################
# Reduced original HiFi bam
# Process HiFi reads and summarize assembly quality.
# Set input, reference and tool paths for your environment before execution.

ori_dir="/path/to/hifibam"
# Output directory for all analyses
wd="/path/to/analysis"
mkdir -p "$wd/run_scripts"
# Sample list
samplefile="/path/to/sample.list"
# Male / Separate female sample list, For DeepVariant And TRGT Analyze each sex separately call variant
malelist="/path/to/sample.male.list"
femalelist="/path/to/sample.female.list"
# Set the thread count, Shared across this script
thread=32
############################END Input settings to update ################################


############################ Software and reference genomes required for analysis ################################
## Software installed in the original environment, Check availability in your environment
N50_py="/path/to/fastq_N50_cal.py"
reads_length_distribution_plot_py="/path/to/fastq_reads_length_distribution_plot.py"
quality_sh="/path/to/bam_quality_cal.sh"
pb_cpg_tools="/path/to/aligned_bam_to_cpg_scores"
pileup_calling_model="/path/to/pileup_calling_model.v1.tflite"
busco_odb="/path/to/primates_odb10"
meryl="/path/to/meryl"
hifiasm_cmd="/path/to/hifiasm"
trgt_cmd="/path/to/trgt"
inspector_cmd="/path/to/inspector.py"
plot_process_bed_cmd="/path/to/process_bed.py"
plot_process_bed_identifier_cmd="/path/to/process_bed_identifiers.py"
minimap2_cmd="/path/to/minimap2"
replace_fasta_cmd="/path/to/02.replace_fasta.py"
plot_karyotype_cmd="/path/to/karyoPlotR.v2_corstrand.R"
connect_cmd="/path/to/connect_parallel.py"
plot_update_bed_identifiers="/path/to/plot_update_bed_identifiers.py"
screen_assembly_pl="/path/to/screen-assembly.pl"
rdna_fa="/path/to/human-rdna-KY962518.1.fasta.gz"
gxdb="/dev/shm/gxdb"
fcs_py="/path/to/fcs.py"
fcs_gx="/path/to/fcs-gx.sif"

##END Software installed in the original environment, Check availability in your environment

## Check the paths in your environment bind path
dv_cmd="singularity run -B /path/to/share -B /path/to/share9 /path/to/deepvariant-1.6.1.sif /opt/deepvariant/bin/run_deepvariant --model_type PACBIO "
pav_cmd="singularity run -B /path/to/share -B /path/to/share9 /path/to/pav.sif"
##END Check the paths in your environment bind path

## Resources used: conda Environment, Configure the required resources locally
pbsv_env="/path/to/pbsv_env"
quast_env="/path/to/quast"
busco_env="/path/to/busco"
methy_env="/path/to/methy_env"
clair3_env="/path/to/clair3_env"
sniffles2_env="/path/to/sniffles2_env"
cuteSV_env="/path/to/cuteSV_env"
svim_env="/path/to/svim_env"
snakemake_env="/path/to/snakemake_env2"
merqury_env="/path/to/merqury"
biotools_env="/path/to/biotools_env"
r4_env="/path/to/r4"
sawfish_env="/path/to/sawfish_env"
hificnv_env="/path/to/hificnv"
mashmap_env="/path/to/mashmap_env"
##END Resources used: conda Environment, Configure the required resources locally

## Reference genomes used
prefix1="chm13v2"
prefix2="t2tcn1"
prefix3="GRCh38"
ref1="/path/to/chm13v2.0.fa"
ref2="/path/to/CN1_combine.v0.8.fasta"
ref3="/path/to/Homo_sapiens.GRCh38.dna.chromosome.1_22_X_Y.fa"
tr1="/path/to/human_chm13v2.0_maskedY_rCRS.trf.bed"
tr2="/path/to/human_t2tCN1_v08_maskedY_rCRS.trf.bed"
tr3="/path/to/Homo_sapiens.GRCh38.dna.chromosome.1_22_X_Y.trf.bed"
gff1="/path/to/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz"
gff2="/path/to/v0.8_Annotation_latest_CN1_v0.8_Annotation_v202306.gff3.gz"
gff3="/path/to/GCF_000001405.26_GRCh38_genomic.chr.md.gff.gz"
trgt1="/path/to/chm13v2.0_maskedY_rCRS.platinumTRs-v1.0.trgt.bed"
trgt3="/path/to/adotto_TRregions_TRGTFormatWithFlankingSeq_v1.0.bed"
virus_db="/path/to/extracted_sequences_db"
cyto1="/path/to/chm13_cytobands.bed"
cyto2="/path/to/cn1_cytobands.bed"
cyto3="/path/to/GRCh38_cytobands.bed"
gs1="/path/to/chm13_length.txt"
gs2="/path/to/cn1_length.txt"
gs3="/path/to/GRCh38_length.txt"
genome1="/path/to/chm13_genomesize.txt"
genome2="/path/to/cn1_genomesize.txt"
genome3="/path/to/GRCh38_genomesize.txt"
female_expected1="/path/to/expected_cn.chm13v2.XX.bed"
female_expected2="/path/to/expected_cn.cn1.XX.bed"
female_expected3="/path/to/expected_cn.hg38.XX.bed"
male_expected1="/path/to/expected_cn.chm13v2.XY.bed"
male_expected2="/path/to/expected_cn.cn1.XY.bed"
male_expected3="/path/to/expected_cn.hg38.XY.bed"
##END Reference genomes used

############################END Software and reference genomes required for analysis ################################


######################### Step.1 Transfer ccs.bam to ccs.fastq #########################
## This script was used to generate ccs.fastq from the ccs bam
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
mkdir -p "${wd}/01.ccs_fastq/${sample}/log"
script="#!/bin/bash
#SBATCH --job-name=b2q_${sample}
#SBATCH --partition=batch
#SBATCH --output=${wd}/01.ccs_fastq/${sample}/log/bam2ccsfq.stdout
#SBATCH --error=${wd}/01.ccs_fastq/${sample}/log/bam2ccsfq.stderr
#SBATCH --mincpus=${thread}
#SBATCH --mem=64g
samtools fastq -n ${ori_dir}/${sample}.hifi_reads.bam|bgzip -@ ${thread} -c > ${wd}/01.ccs_fastq/${sample}/${sample}.ccs.fastq.gz
"
echo -e "${script}" > "${wd}/01.ccs_fastq/${sample}/run.cmd.sh"
echo -e "sh ${wd}/01.ccs_fastq/${sample}/run.cmd.sh" >> "$wd/run_scripts/01.bam2fq.run.sh"
echo -e "sbatch ${wd}/01.ccs_fastq/${sample}/run.cmd.sh" >> "$wd/run_scripts/01.bam2fq.sbatch.sh"
done 9< "${samplefile}"
#########################END Step.1 Transfer ccs.bam to ccs.fastq #########################

######################### Step.2 QC for HiFi fastq reads #########################
## The script is used for HiFi data quality control.
#It can compute Total_base(G), Read_Number, Read_length(max), and Read_length(mean) for BAM files.
#It can calculate the reads N50 for fastq files and plot the distribution of read lengths.
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sampleID <&9 || [[ -n "$sampleID" ]]; do
    [[ -z "$sampleID" ]] && continue
bam=${ori_dir}/${sampleID}.hifi_reads.bam
fastq=${wd}/01.ccs_fastq/${sampleID}/${sampleID}.ccs.fastq.gz
mkdir -p "${wd}/02.reads.QC/$sampleID"
outdir=${wd}/02.reads.QC/$sampleID
cd "$outdir" || exit 1
script="#!/bin/bash
#SBATCH --job-name=QC_${sampleID}
#SBATCH --partition=batch
#SBATCH --output=${outdir}/QC.out
#SBATCH --error=${outdir}/QC.err
#SBATCH --mincpus=8
#SBATCH --mem=64g
cd $outdir
source activate $biotools_env
samtools stats -@8 $bam > $outdir/${sampleID}_stat.txt
python $N50_py $fastq
python $reads_length_distribution_plot_py $fastq
bash $quality_sh $bam
"
echo -e "${script}" > "$outdir/QC.sh"
echo -e "sh $outdir/QC.sh" >> "$wd/run_scripts/02.reads_qc.run.sh"
echo -e "sbatch $outdir/QC.sh" >> "$wd/run_scripts/02.reads_qc.sbatch.sh"
done 9< "${samplefile}"
######################### END Step.2 QC for HiFi fastq reads #########################

######################### Step.3 Using hifiasm to assembly the HiFi reads #########################
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
mkdir -p "${wd}/03.hifiasm/${sample}/log"
hifiasm_script="#!/bin/bash
#SBATCH --job-name=${sample}_asm
#SBATCH --partition=batch
#SBATCH --output=${wd}/03.hifiasm/${sample}/log/hifiasm.stdout
#SBATCH --error=${wd}/03.hifiasm/${sample}/log/hifiasm.stderr
#SBATCH --mincpus=64
#SBATCH --mem=256g
/usr/bin/time -v $hifiasm_cmd -o ${wd}/03.hifiasm/${sample}/$sample.asm -t 64 ${wd}/01.ccs_fastq/${sample}/${sample}.ccs.fastq.gz
awk '/^S/{print \">\"\$2\"\\\\n\"\$3}' ${wd}/03.hifiasm/${sample}/${sample}.asm.bp.hap1.p_ctg.gfa > ${wd}/03.hifiasm/${sample}/${sample}.asm.bp.hap1.p_ctg.fa
awk '/^S/{print \">\"\$2\"\\\\n\"\$3}' ${wd}/03.hifiasm/${sample}/${sample}.asm.bp.hap2.p_ctg.gfa > ${wd}/03.hifiasm/${sample}/${sample}.asm.bp.hap2.p_ctg.fa
samtools faidx ${wd}/03.hifiasm/${sample}/${sample}.asm.bp.hap1.p_ctg.fa
samtools faidx ${wd}/03.hifiasm/${sample}/${sample}.asm.bp.hap2.p_ctg.fa
"
echo -e "${hifiasm_script}" > "${wd}/03.hifiasm/${sample}/run.hifiasm.cmd.sh"
echo -e "sh ${wd}/03.hifiasm/${sample}/run.hifiasm.cmd.sh" >> "$wd/run_scripts/03.hifiasm.run.sh"
echo -e "sbatch ${wd}/03.hifiasm/${sample}/run.hifiasm.cmd.sh" >> "$wd/run_scripts/03.hifiasm.sbatch.sh"
done 9< "${samplefile}"
#########################END Step.3 Using hifiasm to assembly the HiFi reads #########################

######################### Step.4 Map the HiFi reads to reference genome #########################
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
for i in {1..3}
do
prefix_var="prefix${i}"
ref_var="ref${i}"
prefix=${!prefix_var}
ref=${!ref_var}
mkdir -p "${wd}/04.reads_align/${prefix}/${sample}/log"
align_script="#!/bin/bash
#SBATCH --job-name=${sample}_${prefix}
#SBATCH --partition=batch
#SBATCH --output=${wd}/04.reads_align/${prefix}/${sample}/log/reads_align.stdout
#SBATCH --error=${wd}/04.reads_align/${prefix}/${sample}/log/reads_align.stderr
#SBATCH --mincpus=${thread}
#SBATCH --mem=128g
source activate $pbsv_env
pbmm2 align ${ref} ${ori_dir}/${sample}.hifi_reads.bam ${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam --preset CCS --sample ${sample} --sort
samtools index ${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam
"
echo -e "${align_script}" > "${wd}/04.reads_align/${prefix}/${sample}/run.align.cmd.sh"
echo -e "sh ${wd}/04.reads_align/${prefix}/${sample}/run.align.cmd.sh" >> "$wd/run_scripts/04.align_pbmm2.${prefix}.run.sh"
echo -e "sbatch ${wd}/04.reads_align/${prefix}/${sample}/run.align.cmd.sh" >> "$wd/run_scripts/04.align_pbmm2.${prefix}.sbatch.sh"
done
done 9< "${samplefile}"
#########################END Step.4 Map the HiFi reads to reference genome #########################


######################### Step.5 DeepVariant #########################
## DeepVariant is a deep learning-based variant caller that uses a deep neural network to predict the genotype of a variant in a single read.
## The script is used to call variants from HiFi data using DeepVariant.
# Run male and female samples separately
## Female
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
for i in {1..3}
do
prefix_var="prefix${i}"
ref_var="ref${i}"
tr_var="tr${i}"
prefix=${!prefix_var}
ref=${!ref_var}
mkdir -p "${wd}/05.SNV_calling/${prefix}/dv/${sample}/log"
dv_script="#!/bin/bash
#SBATCH --job-name=${sample}dv
#SBATCH --partition=batch
#SBATCH --output=${wd}/05.SNV_calling/${prefix}/dv/${sample}/log/deepvariant_1.6.1.stdout
#SBATCH --error=${wd}/05.SNV_calling/${prefix}/dv/${sample}/log/deepvariant_1.6.1.stderr
#SBATCH --mincpus=${thread}
#SBATCH --mem=128g

/usr/bin/time -v ${dv_cmd} --ref ${ref} \
--reads ${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam \
--output_vcf=${wd}/05.SNV_calling/${prefix}/dv/${sample}/${sample}.${prefix}.dv161.vcf \
--output_gvcf=${wd}/05.SNV_calling/${prefix}/dv/${sample}/${sample}.${prefix}.dv161.g.vcf \
--num_shards=${thread} \
--sample_name ${sample}
bcftools view -f PASS ${wd}/05.SNV_calling/${prefix}/dv/${sample}/${sample}.${prefix}.dv161.vcf -o ${wd}/05.SNV_calling/${prefix}/dv/${sample}/${sample}.${prefix}.dv161.PASS.vcf
bgzip ${wd}/05.SNV_calling/${prefix}/dv/${sample}/${sample}.${prefix}.dv161.PASS.vcf
tabix -f ${wd}/05.SNV_calling/${prefix}/dv/${sample}/${sample}.${prefix}.dv161.PASS.vcf.gz
"
echo -e "${dv_script}" > "${wd}/05.SNV_calling/${prefix}/dv/${sample}/run.dv161.cmd.sh"
echo -e "sh ${wd}/05.SNV_calling/${prefix}/dv/${sample}/run.dv161.cmd.sh" >> "$wd/run_scripts/05.dv161.${prefix}.female.run.sh"
echo -e "sbatch ${wd}/05.SNV_calling/${prefix}/dv/${sample}/run.dv161.cmd.sh" >> "$wd/run_scripts/05.dv161.${prefix}.female.sbatch.sh"
done
done 9< "${femalelist}"
## Male
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
for i in {1..3}
do
prefix_var="prefix${i}"
ref_var="ref${i}"
tr_var="tr${i}"
prefix=${!prefix_var}
ref=${!ref_var}
mkdir -p "${wd}/05.SNV_calling/${prefix}/dv/${sample}/log"
dv_script="#!/bin/bash
#SBATCH --job-name=${sample}dv
#SBATCH --partition=batch
#SBATCH --output=${wd}/05.SNV_calling/${prefix}/dv/${sample}/log/deepvariant_1.6.1.stdout
#SBATCH --error=${wd}/05.SNV_calling/${prefix}/dv/${sample}/log/deepvariant_1.6.1.stderr
#SBATCH --mincpus=${thread}
#SBATCH --mem=32g

/usr/bin/time -v ${dv_cmd} --ref ${ref} \
--reads ${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam \
--output_vcf=${wd}/05.SNV_calling/${prefix}/dv/${sample}/${sample}.${prefix}.dv161.vcf \
--output_gvcf=${wd}/05.SNV_calling/${prefix}/dv/${sample}/${sample}.${prefix}.dv161.g.vcf \
--num_shards=${thread} \
--sample_name ${sample} \
--haploid_contigs=\"chrX,chrY\"
bcftools view -f PASS ${wd}/05.SNV_calling/${prefix}/dv/${sample}/${sample}.${prefix}.dv161.vcf -o ${wd}/05.SNV_calling/${prefix}/dv/${sample}/${sample}.${prefix}.dv161.PASS.vcf
bgzip ${wd}/05.SNV_calling/${prefix}/dv/${sample}/${sample}.${prefix}.dv161.PASS.vcf
tabix -f ${wd}/05.SNV_calling/${prefix}/dv/${sample}/${sample}.${prefix}.dv161.PASS.vcf.gz
"
echo -e "${dv_script}" > "${wd}/05.SNV_calling/${prefix}/dv/${sample}/run.dv161.cmd.sh"
echo -e "sh ${wd}/05.SNV_calling/${prefix}/dv/${sample}/run.dv161.cmd.sh" >> "$wd/run_scripts/05.dv161.${prefix}.male.run.sh"
echo -e "sbatch ${wd}/05.SNV_calling/${prefix}/dv/${sample}/run.dv161.cmd.sh" >> "$wd/run_scripts/05.dv161.${prefix}.male.sbatch.sh"
done
done 9< "${malelist}"

#########################END Step.5 DeepVariant #########################


######################### Step.6 5mC signal to CpG signal #########################
## 5mC signal to CpG signal is used to quantify the methylation level of CpGs in the genome.
## The script is used to calculate 5mC signal using pb-CpG-tools.
mkdir -p "${wd}/06.5mC/${prefix}"
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
for i in {1..3}
do
prefix_var="prefix${i}"
ref_var="ref${i}"
prefix=${!prefix_var}
ref=${!ref_var}
mkdir -p "${wd}/06.5mC/${prefix}/${sample}/log"
## Without separating haplotype
meth_script="#!/bin/bash
#SBATCH --job-name=5mC_${sample}
#SBATCH --partition=batch
#SBATCH --output=${wd}/06.5mC/${prefix}/${sample}/log/5mC.stdout
#SBATCH --error=${wd}/06.5mC/${prefix}/${sample}/log/5mC.stderr
#SBATCH --mincpus=${thread}
#SBATCH --mem=64g
source activate $methy_env
#/usr/bin/time -v ccsmeth align_hifi --hifireads ${ori_dir}/${sample}.hifi_reads.bam --ref ${ref} -o ${wd}/5mC/${prefix}/${sample}/${sample}.pbmm2.5mc.ccs.alignment.bam --threads ${thread}
${pb_cpg_tools} --min-coverage 4 \
--modsites-mode denovo \
--pileup-mode model \
--min-mapq 10 --threads ${thread} \
--model ${pileup_calling_model} \
--bam ${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam \
--ref ${ref} \
--output-prefix ${wd}/06.5mC/${prefix}/${sample}/${sample}.${prefix}.denovo
${pb_cpg_tools} --min-coverage 4 \
--modsites-mode reference \
--pileup-mode model \
--min-mapq 10 --threads ${thread} \
--model ${pileup_calling_model} \
--bam ${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam \
--ref ${ref} \
--output-prefix ${wd}/06.5mC/${prefix}/${sample}/${sample}.${prefix}.reference
"
echo -e "${meth_script}" > "${wd}/06.5mC/${prefix}/${sample}/run.5mC.cmd.sh"
echo -e "sh ${wd}/06.5mC/${prefix}/${sample}/run.5mC.cmd.sh" >> "$wd/run_scripts/06.5mC.${prefix}.run.sh"
echo -e "sbatch ${wd}/06.5mC/${prefix}/${sample}/run.5mC.cmd.sh" >> "$wd/run_scripts/06.5mC.${prefix}.sbatch.sh"
# Split haplotype
meth_hap_script="#!/bin/bash
#SBATCH --job-name=5mC_${sample}
#SBATCH --partition=batch
#SBATCH --output=${wd}/06.5mC/${prefix}/${sample}/log/5mC.hap.stdout
#SBATCH --error=${wd}/06.5mC/${prefix}/${sample}/log/5mC.hap.stderr
#SBATCH --mincpus=${thread}
source activate $clair3_env
module load samtools/1.15.1
/usr/bin/time -v whatshap phase -o ${wd}/06.5mC/${prefix}/${sample}/${sample}.${prefix}.dv161.PASS.phase.vcf --tag=HP --reference ${ref} ${wd}/05.SNV_calling/${prefix}/dv/${sample}/${sample}.${prefix}.dv161.PASS.vcf.gz ${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam
bgzip ${wd}/06.5mC/${prefix}/${sample}/${sample}.${prefix}.dv161.PASS.phase.vcf
tabix -f ${wd}/06.5mC/${prefix}/${sample}/${sample}.${prefix}.dv161.PASS.phase.vcf.gz
/usr/bin/time -v whatshap haplotag -o ${wd}/06.5mC/${prefix}/${sample}/${sample}.${prefix}.pbmm2.5mC.phase.bam --reference ${ref} ${wd}/06.5mC/${prefix}/${sample}/${sample}.${prefix}.dv161.PASS.phase.vcf.gz ${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam --output-threads=${thread}
samtools index -@${thread} ${wd}/06.5mC/${prefix}/${sample}/${sample}.${prefix}.pbmm2.5mC.phase.bam
${pb_cpg_tools} --min-coverage 4 \
--hap-tag HP \
--modsites-mode denovo \
--pileup-mode model \
--min-mapq 10 --threads ${thread} \
--model ${pileup_calling_model} \
--bam ${wd}/06.5mC/${prefix}/${sample}/${sample}.${prefix}.pbmm2.5mC.phase.bam \
--ref ${ref} \
--output-prefix ${wd}/06.5mC/${prefix}/${sample}/${sample}.${prefix}.denovo

${pb_cpg_tools} --min-coverage 4 \
--hap-tag HP \
--modsites-mode reference \
--pileup-mode model \
--min-mapq 10 --threads ${thread} \
--model ${pileup_calling_model} \
--bam ${wd}/06.5mC/${prefix}/${sample}/${sample}.${prefix}.pbmm2.5mC.phase.bam \
--ref ${ref} \
--output-prefix ${wd}/06.5mC/${prefix}/${sample}/${sample}.${prefix}.reference
"
echo -e "$meth_hap_script" > "${wd}/06.5mC/${prefix}/${sample}/run.5mC.hap.cmd.sh"
echo -e "sh ${wd}/06.5mC/${prefix}/${sample}/run.5mC.hap.cmd.sh" >> "$wd/run_scripts/06.5mC.hap.${prefix}.run.sh"
echo -e "sbatch ${wd}/06.5mC/${prefix}/${sample}/run.5mC.hap.cmd.sh" >> "$wd/run_scripts/06.5mC.hap.${prefix}.sbatch.sh"
done
done 9< "${samplefile}"
#########################END Step.6 5mC signal to CpG signal #########################

######################## Step.7 Evaluate the assembly ########################
# 1. QUAST
# The script is used to calculate assembly metrics such as N50, NA50, and Genome Fraction using QUAST.
while IFS= read -r sample; do
outdir=${wd}/07.hifiasm_QC/${sample}/quast
mkdir -p "$outdir"
cd "$outdir" || exit 1
asm1=${wd}/03.hifiasm/${sample}/${sample}.asm.bp.hap1.p_ctg.fa
asm2=${wd}/03.hifiasm/${sample}/${sample}.asm.bp.hap2.p_ctg.fa
script="#!/bin/bash
#SBATCH --job-name=quast_${sample}
#SBATCH --partition=batch
#SBATCH --output=${outdir}/run_quast.out
#SBATCH --error=${outdir}/run_quast.err
#SBATCH --mincpus=2
#SBATCH --mem=40g
source activate $quast_env
quast.py $asm1 -r $ref1 -g $gff1 -o ${outdir}/${sample}_hap1 -t 10
quast.py $asm2 -r $ref1 -g $gff1 -o ${outdir}/${sample}_hap2 -t 10
"

echo -e "$script" > "${wd}/07.hifiasm_QC/${sample}/quast/run_quast.sh"
echo -e "sh ${wd}/07.hifiasm_QC/${sample}/quast/run_quast.sh" >> "$wd/run_scripts/07.asmQC.quast.run.sh"
echo -e "sbatch ${wd}/07.hifiasm_QC/${sample}/quast/run_quast.sh" >> "$wd/run_scripts/07.asmQC.quast.sbatch.sh"
done < "$samplefile"

# 2. BUSCO
#The script is used to assess assembly gene region completeness using BUSCO.
##To rerun the busco script, add the parameter -f. Must be used when output files with the provided name already exist.
##It is recommended to submit the BUSCO tasks to the SLURM system, as running multiple samples' tasks with multi-threading on a single node is prone to failure.
while IFS= read -r sample; do
outdir1=${wd}/07.hifiasm_QC/${sample}/busco/${sample}_hap1
mkdir -p "$outdir1"
cd "$outdir1" || exit 1
asm1=${wd}/03.hifiasm/${sample}/${sample}.asm.bp.hap1.p_ctg.fa
script1="#!/bin/bash
#SBATCH --job-name=busco_${sample}_hap1
#SBATCH --partition=batch
#SBATCH --output=${outdir1}/run_busco.out
#SBATCH --error=${outdir1}/run_busco.err
#SBATCH --mincpus=8
#SBATCH --mem=32g
source activate $busco_env
busco -i $asm1 -l $busco_odb -o busco -m geno -c 8 --offline #-f
"
echo -e "$script1" > "$outdir1/run_busco.sh"
echo -e "sh $outdir1/run_busco.sh" >> "$wd/run_scripts/07.asmQC.busco.run.sh"
echo -e "sbatch $outdir1/run_busco.sh" >> "$wd/run_scripts/07.asmQC.busco.sbatch.sh"

outdir2=${wd}/07.hifiasm_QC/${sample}/busco/${sample}_hap2
mkdir -p "$outdir2"
cd "$outdir2" || exit 1
asm2=${wd}/03.hifiasm/${sample}/${sample}.asm.bp.hap2.p_ctg.fa
script2="#!/bin/bash
#SBATCH --job-name=busco_${sample}_hap2
#SBATCH --partition=batch
#SBATCH --output=${outdir2}/run_busco.out
#SBATCH --error=${outdir2}/run_busco.err
#SBATCH --mincpus=8
#SBATCH --mem=32g
source activate $busco_env
busco -i $asm2 -l $busco_odb -o busco -m geno -c 8 --offline #-f
"
echo -e "$script2" > "$outdir2/run_busco.sh"
echo -e "sh $outdir2/run_busco.sh" >> "$wd/run_scripts/07.asmQC.busco.run.sh"
echo -e "sbatch $outdir2/run_busco.sh" >> "$wd/run_scripts/07.asmQC.busco.sbatch.sh"
done < "$samplefile"

# 3. merqury
#!/bin/bash
#The script is used to calculate k-mer completeness and QV of the assembly using Merqury.
##To rerun the merqury script, first remove the existing results folder.

while IFS= read -r sample; do
mkdir -p "${wd}/07.hifiasm_QC/${sample}/merqury/log"
mkdir -p "${wd}/07.hifiasm_QC/${sample}/merqury/meryl"

reads=${wd}/01.ccs_fastq/${sample}/${sample}.ccs.fastq.gz
kmers=${wd}/07.hifiasm_QC/${sample}/merqury/meryl/${sample}.meryl
asm1=${wd}/03.hifiasm/${sample}/${sample}.asm.bp.hap1.p_ctg.fa
asm2=${wd}/03.hifiasm/${sample}/${sample}.asm.bp.hap2.p_ctg.fa
script="#!/bin/bash
#SBATCH --job-name=merqury_$sample
#SBATCH --partition=batch
#SBATCH --output=${wd}/07.hifiasm_QC/${sample}/merqury/log/run_merqury.out
#SBATCH --error=${wd}/07.hifiasm_QC/${sample}/merqury/log/run_merqury.err
#SBATCH --mincpus=32
#SBATCH --mem=128g
if [ ! -e \"$kmers\" ]; then
$meryl count k=21 $reads output $kmers
fi
ln -s $asm1 ${sample}.asm.hap1.fasta
ln -s $asm2 ${sample}.asm.hap2.fasta
source activate $merqury_env
cd ${wd}/07.hifiasm_QC/${sample}/merqury/output
merqury.sh $kmers ${sample}.asm.hap1.fasta ${sample}.asm.hap2.fasta merqury
"
mkdir -p "${wd}/07.hifiasm_QC/${sample}/merqury/output"
cd "${wd}/07.hifiasm_QC/${sample}/merqury/output" || exit 1
echo -e "$script" > "${wd}/07.hifiasm_QC/${sample}/merqury/output/run_merqury.sh"
echo -e "sh ${wd}/07.hifiasm_QC/${sample}/merqury/output/run_merqury.sh" >> "$wd/run_scripts/07.asmQC.merqury.run.sh"
echo -e "sbatch ${wd}/07.hifiasm_QC/${sample}/merqury/output/run_merqury.sh" >> "$wd/run_scripts/07.asmQC.merqury.sbatch.sh"
done < "$samplefile"

########################END Step.7 Evaluate the assembly ########################

######################## Step.8 Align the assembly to reference genome #################
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
for i in {1..3}
do
prefix_var="prefix${i}"
ref_var="ref${i}"
prefix=${!prefix_var}
ref=${!ref_var}
mkdir -p "${wd}/08.asm_align/${prefix}/${sample}/log"
asm_align_script="#!/bin/bash
#SBATCH --job-name=m_$sample
#SBATCH --partition=batch
#SBATCH --output=${wd}/08.asm_align/${prefix}/${sample}/log/${sample}.minimap2.stdout
#SBATCH --error=${wd}/08.asm_align/${prefix}/${sample}/log/${sample}.minimap2.stderr
#SBATCH --mincpus=${thread}
#SBATCH --mem=48g
/usr/bin/time -v ${minimap2_cmd} -ax asm5 --cs -r2k -t ${thread} \
${ref} \
${wd}/03.hifiasm/${sample}/${sample}.asm.bp.hap1.p_ctg.fa > ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sam
samtools view -bS -@${thread} ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sam -o ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.bam
samtools sort -@${thread} ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.bam -o ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sorted.bam
samtools index ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sorted.bam
/usr/bin/time -v ${minimap2_cmd} -ax asm5 --cs -r2k -t ${thread} \
${ref} \
${wd}/03.hifiasm/${sample}/${sample}.asm.bp.hap2.p_ctg.fa > ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sam
samtools view -bS -@${thread} ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sam -o ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.bam
samtools sort -@${thread} ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.bam -o ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sorted.bam
samtools index ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sorted.bam
rm ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sam
rm ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sam
"
echo -e "$asm_align_script" > "${wd}/08.asm_align/${prefix}/${sample}/run.align.cmd.sh"
echo -e "sh ${wd}/08.asm_align/${prefix}/${sample}/run.align.cmd.sh" >> "$wd/run_scripts/08.asm_align.${prefix}.run.sh"
echo -e "sbatch ${wd}/08.asm_align/${prefix}/${sample}/run.align.cmd.sh" >> "$wd/run_scripts/08.asm_align.${prefix}.sbatch.sh"
done
done 9< "${samplefile}"

########################END Step.8 Align the assembly to reference genome #################



######################## Step.9 Call SV ########################
#### pbsv/sniffles2/cuteSV/PAV/SVIM-asm/TRGT

# Step.9-1 Run pbsv
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
for i in {1..3}
do
prefix_var="prefix${i}"
ref_var="ref${i}"
tr_var="tr${i}"
prefix=${!prefix_var}
ref=${!ref_var}
tr=${!tr_var}
mkdir -p "${wd}/09.SV_calling/${prefix}/pbsv/${sample}/log"
pbsv_script="#!/bin/bash
#SBATCH --job-name=${sample}pbsv
#SBATCH --partition=batch
#SBATCH --output=${wd}/09.SV_calling/${prefix}/pbsv/${sample}/log/pbsv.stdout
#SBATCH --error=${wd}/09.SV_calling/${prefix}/pbsv/${sample}/log/pbsv.stderr
#SBATCH --mincpus=${thread}
#SBATCH --mem=64g
source activate $pbsv_env
/usr/bin/time -v pbsv discover --hifi -s ${sample} \
${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam \
${wd}/09.SV_calling/${prefix}/pbsv/${sample}/${sample}.${prefix}.pbmm2.svsig.gz \
--tandem-repeats ${tr}
tabix -c '#' -s 3 -b 4 -e 4 ${wd}/09.SV_calling/${prefix}/pbsv/${sample}/${sample}.${prefix}.pbmm2.svsig.gz
/usr/bin/time -v pbsv call ${ref} ${wd}/09.SV_calling/${prefix}/pbsv/${sample}/${sample}.${prefix}.pbmm2.svsig.gz ${wd}/09.SV_calling/${prefix}/pbsv/${sample}/${sample}.${prefix}.pbsv.vcf
"
echo -e "${pbsv_script}" > "${wd}/09.SV_calling/${prefix}/pbsv/${sample}/run.pbsv.cmd.sh"
echo -e "sh ${wd}/09.SV_calling/${prefix}/pbsv/${sample}/run.pbsv.cmd.sh" >>  "$wd/run_scripts/09.pbsv.${prefix}.run.sh"
echo -e "sbatch ${wd}/09.SV_calling/${prefix}/pbsv/${sample}/run.pbsv.cmd.sh" >> "$wd/run_scripts/09.pbsv.${prefix}.sbatch.sh"
done
done 9< "${samplefile}"

# Step.9-2 Run sniffles2
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
for i in {1..3}
do
prefix_var="prefix${i}"
ref_var="ref${i}"
tr_var="tr${i}"
prefix=${!prefix_var}
ref=${!ref_var}
tr=${!tr_var}
mkdir -p "${wd}/09.SV_calling/${prefix}/sniffles2/${sample}/log"
sniffles2_script="#!/bin/bash
#SBATCH --job-name=${sample}sni
#SBATCH --partition=batch
#SBATCH --output=${wd}/09.SV_calling/${prefix}/sniffles2/${sample}/log/sniffles.stdout
#SBATCH --error=${wd}/09.SV_calling/${prefix}/sniffles2/${sample}/log/sniffles.stderr
#SBATCH --mincpus=${thread}
#SBATCH --mem=16g
source activate $sniffles2_env
/usr/bin/time -v sniffles --input ${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam \
--vcf ${wd}/09.SV_calling/${prefix}/sniffles2/${sample}/${sample}.${prefix}.sniffles2.vcf.gz \
--snf ${wd}/09.SV_calling/${prefix}/sniffles2/${sample}/${sample}.${prefix}.sniffles2.snf \
--tandem-repeats ${tr} \
--sample-id ${sample}
"
echo -e "${sniffles2_script}" > "${wd}/09.SV_calling/${prefix}/sniffles2/${sample}/run.sniffles2.cmd.sh"
echo -e "sh ${wd}/09.SV_calling/${prefix}/sniffles2/${sample}/run.sniffles2.cmd.sh" >> "$wd/run_scripts/09.sniffles.${prefix}.run.sh"
echo -e "sbatch ${wd}/09.SV_calling/${prefix}/sniffles2/${sample}/run.sniffles2.cmd.sh" >> "$wd/run_scripts/09.sniffles.${prefix}.sbatch.sh"
done
done 9< "${samplefile}"

# Step.9-3 Run cuteSV
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
for i in {1..3}
do
prefix_var="prefix${i}"
ref_var="ref${i}"
tr_var="tr${i}"
prefix=${!prefix_var}
ref=${!ref_var}
tr=${!tr_var}
mkdir -p "${wd}/09.SV_calling/${prefix}/cuteSV/${sample}/log"
cuteSV_script="#!/bin/bash
#SBATCH --job-name=${sample}cuteSV
#SBATCH --partition=batch
#SBATCH --output=${wd}/09.SV_calling/${prefix}/cuteSV/${sample}/log/cuteSV.stdout
#SBATCH --error=${wd}/09.SV_calling/${prefix}/cuteSV/${sample}/log/cuteSV.stderr
#SBATCH --mincpus=${thread}
#SBATCH --mem=16g
source activate $cuteSV_env
/usr/bin/time -v cuteSV ${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam \
$ref \
${wd}/09.SV_calling/${prefix}/cuteSV/${sample}/${sample}.${prefix}.cuteSV.vcf \
${wd}/09.SV_calling/${prefix}/cuteSV/${sample} \
--max_cluster_bias_INS 1000 \
--diff_ratio_merging_INS 0.9 \
--max_cluster_bias_DEL 1000 \
--diff_ratio_merging_DEL 0.5 \
--threads ${thread} \
--sample ${sample} \
--genotype
"
echo -e "${cuteSV_script}" > "${wd}/09.SV_calling/${prefix}/cuteSV/${sample}/run.cuteSV.cmd.sh"
echo -e "sh ${wd}/09.SV_calling/${prefix}/cuteSV/${sample}/run.cuteSV.cmd.sh" >> "$wd/run_scripts/09.cuteSV.${prefix}.run.sh"
echo -e "sbatch ${wd}/09.SV_calling/${prefix}/cuteSV/${sample}/run.cuteSV.cmd.sh" >> "$wd/run_scripts/09.cuteSV.${prefix}.sbatch.sh"
#cuteSV Run time compared with sniffles2 Nearly identical, Fast
done
done 9< "${samplefile}"

# Step.9-4 Run sawfish
# Male
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
for i in {1..3}
do
prefix_var="prefix${i}"
ref_var="ref${i}"
tr_var="tr${i}"
cn_var="male_expected${i}"
prefix=${!prefix_var}
ref=${!ref_var}
tr=${!tr_var}
cn=${!cn_var}

mkdir -p "${wd}/09.SV_calling/${prefix}/sawfish/${sample}/log"
sawfish_script="#!/bin/bash
#SBATCH --job-name=${sample}sawfish
#SBATCH --partition=batch
#SBATCH --output=${wd}/09.SV_calling/${prefix}/sawfish/${sample}/log/sawfish.stdout
#SBATCH --error=${wd}/09.SV_calling/${prefix}/sawfish/${sample}/log/sawfish.stderr
#SBATCH --mincpus=${thread}
#SBATCH --mem=16g
source activate ${sawfish_env}
cd ${wd}/09.SV_calling/${prefix}/sawfish/${sample}
/usr/bin/time -v  sawfish discover \
--threads ${thread} \
--ref ${ref} \
--bam ${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam \
--output-dir ${sample}_sawfish_discover_dir \
--expected-cn ${cn}

# the --output-dir should not exist
# the --sample should be the output_dir of the discover step
/usr/bin/time -v sawfish joint-call --threads ${thread} \
--sample ${sample}_sawfish_discover_dir \
--output-dir ${sample}_sawfish_call_dir \
--report-supporting-reads
"
echo -e "${sawfish_script}" > "${wd}/09.SV_calling/${prefix}/sawfish/${sample}/run.sawfish.cmd.sh"
echo -e "sh ${wd}/09.SV_calling/${prefix}/sawfish/${sample}/run.sawfish.cmd.sh" >> "$wd/run_scripts/09.sawfish.${prefix}.run.sh"
echo -e "sbatch ${wd}/09.SV_calling/${prefix}/sawfish/${sample}/run.sawfish.cmd.sh" >> "$wd/run_scripts/09.sawfish.${prefix}.sbatch.sh"
done
done 9< "${malelist}"

# Female
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
for i in {1..3}
do
prefix_var="prefix${i}"
ref_var="ref${i}"
tr_var="tr${i}"
cn_var="female_expected${i}"
prefix=${!prefix_var}
ref=${!ref_var}
tr=${!tr_var}
cn=${!cn_var}
mkdir -p "${wd}/09.SV_calling/${prefix}/sawfish/${sample}/log"
sawfish_script="#!/bin/bash
#SBATCH --job-name=${sample}sawfish
#SBATCH --partition=batch
#SBATCH --output=${wd}/09.SV_calling/${prefix}/sawfish/${sample}/log/sawfish.stdout
#SBATCH --error=${wd}/09.SV_calling/${prefix}/sawfish/${sample}/log/sawfish.stderr
#SBATCH --mincpus=${thread}
#SBATCH --mem=16g
source activate ${sawfish_env}
cd ${wd}/09.SV_calling/${prefix}/sawfish/${sample}
/usr/bin/time -v  sawfish discover \
--threads ${thread} \
--ref ${ref} \
--bam ${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam \
--output-dir ${sample}_sawfish_discover_dir \
--expected-cn ${cn}

# the --output-dir should not exist
# the --sample should be the output_dir of the discover step
/usr/bin/time -v sawfish joint-call --threads ${thread} \
--sample ${sample}_sawfish_discover_dir \
--output-dir ${sample}_sawfish_call_dir \
--report-supporting-reads
"
echo -e "${sawfish_script}" > "${wd}/09.SV_calling/${prefix}/sawfish/${sample}/run.sawfish.cmd.sh"
echo -e "sh ${wd}/09.SV_calling/${prefix}/sawfish/${sample}/run.sawfish.cmd.sh" >> "$wd/run_scripts/09.sawfish.${prefix}.run.sh"
echo -e "sbatch ${wd}/09.SV_calling/${prefix}/sawfish/${sample}/run.sawfish.cmd.sh" >> "$wd/run_scripts/09.sawfish.${prefix}.sbatch.sh"
done
done 9< "${femalelist}"

# Step.9-5 Run hificnv
# Male
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
for i in {1..3}
do
prefix_var="prefix${i}"
ref_var="ref${i}"
cn_var="male_expected${i}"
prefix=${!prefix_var}
ref=${!ref_var}
cn=${!cn_var}

mkdir -p "${wd}/09.SV_calling/${prefix}/hificnv/${sample}/log"
hificnv_script="#!/bin/bash
#SBATCH --job-name=${sample}sawfish
#SBATCH --partition=batch
#SBATCH --output=${wd}/09.SV_calling/${prefix}/hificnv/${sample}/log/hificnv.stdout
#SBATCH --error=${wd}/09.SV_calling/${prefix}/hificnv/${sample}/log/hificnv.stderr
#SBATCH --mincpus=${thread}
#SBATCH --mem=16g
source activate ${hificnv_env}
hificnv --ref ${ref} \
--bam ${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam \
--output-prefix ${wd}/09.SV_calling/${prefix}/hificnv/${sample}/${sample}.${prefix}.hificnv \
--expected-cn ${cn}
"
echo -e "${hificnv_script}" > "${wd}/09.SV_calling/${prefix}/hificnv/${sample}/run.hificnv.cmd.sh"
echo -e "sh ${wd}/09.SV_calling/${prefix}/hificnv/${sample}/run.hificnv.cmd.sh" >> "$wd/run_scripts/09.hificnv.${prefix}.run.sh"
echo -e "sbatch ${wd}/09.SV_calling/${prefix}/hificnv/${sample}/run.hificnv.cmd.sh" >> "$wd/run_scripts/09.hificnv.${prefix}.sbatch.sh"
done
done 9< "${malelist}"

# Female
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
for i in {1..3}
do
prefix_var="prefix${i}"
ref_var="ref${i}"
cn_var="female_expected${i}"
prefix=${!prefix_var}
ref=${!ref_var}
cn=${!cn_var}

mkdir -p "${wd}/09.SV_calling/${prefix}/hificnv/${sample}/log"
hificnv_script="#!/bin/bash
#SBATCH --job-name=${sample}sawfish
#SBATCH --partition=batch
#SBATCH --output=${wd}/09.SV_calling/${prefix}/hificnv/${sample}/log/hificnv.stdout
#SBATCH --error=${wd}/09.SV_calling/${prefix}/hificnv/${sample}/log/hificnv.stderr
#SBATCH --mincpus=${thread}
#SBATCH --mem=16g
source activate ${hificnv_env}
hificnv --ref ${ref} \
--bam ${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam \
--output-prefix ${wd}/09.SV_calling/${prefix}/hificnv/${sample}/${sample}.${prefix}.hificnv \
--expected-cn ${cn}
"
echo -e "${hificnv_script}" > "${wd}/09.SV_calling/${prefix}/hificnv/${sample}/run.hificnv.cmd.sh"
echo -e "sh ${wd}/09.SV_calling/${prefix}/hificnv/${sample}/run.hificnv.cmd.sh" >> "$wd/run_scripts/09.hificnv.${prefix}.run.sh"
echo -e "sbatch ${wd}/09.SV_calling/${prefix}/hificnv/${sample}/run.hificnv.cmd.sh" >> "$wd/run_scripts/09.hificnv.${prefix}.sbatch.sh"
done
done 9< "${femalelist}"

# Step.9-6 Run SVIM-asm
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
for i in {1..3}
do
prefix_var="prefix${i}"
ref_var="ref${i}"
prefix=${!prefix_var}
ref=${!ref_var}

mkdir -p "${wd}/09.SV_calling/${prefix}/SVIM-asm/${sample}/SVIM-asm-dip/"
mkdir -p "${wd}/09.SV_calling/${prefix}/SVIM-asm/${sample}/SVIM-asm-hap1/"
mkdir -p "${wd}/09.SV_calling/${prefix}/SVIM-asm/${sample}/SVIM-asm-hap2/"
mkdir -p "${wd}/09.SV_calling/${prefix}/SVIM-asm/${sample}/log"
svim_asm_script="#!/bin/bash
#SBATCH --job-name=${sample}_svim
#SBATCH --partition=batch
#SBATCH --output=${wd}/09.SV_calling/${prefix}/SVIM-asm/${sample}/log/${sample}.svim.out
#SBATCH --error=${wd}/09.SV_calling/${prefix}/SVIM-asm/${sample}/log/${sample}.svim.err
#SBATCH --mincpus=4
#SBATCH --mem=32g
source activate $svim_env
svim-asm haploid ${wd}/09.SV_calling/${prefix}/SVIM-asm/${sample}/SVIM-asm-hap1 ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sorted.bam $ref
svim-asm haploid ${wd}/09.SV_calling/${prefix}/SVIM-asm/${sample}/SVIM-asm-hap2 ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sorted.bam $ref
svim-asm diploid ${wd}/09.SV_calling/${prefix}/SVIM-asm/${sample}/SVIM-asm-dip ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sorted.bam ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sorted.bam $ref
"
echo -e "$svim_asm_script" > "${wd}/09.SV_calling/${prefix}/SVIM-asm/${sample}/run.svim-asm.sh"
echo -e "sh ${wd}/09.SV_calling/${prefix}/SVIM-asm/${sample}/run.svim-asm.sh" >> "$wd/run_scripts/09.svimasm.${prefix}.run.sh"
echo -e "sbatch ${wd}/09.SV_calling/${prefix}/SVIM-asm/${sample}/run.svim-asm.sh" >> "$wd/run_scripts/09.svimasm.${prefix}.sbatch.sh"
done
done 9< "${samplefile}"

# Step.9-7 Run PAV

# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
for i in {1..3}
do
prefix_var="prefix${i}"
ref_var="ref${i}"
prefix=${!prefix_var}
ref=${!ref_var}
mkdir -p "${wd}/09.SV_calling/${prefix}/PAV/${sample}/log"
mkdir -p "${wd}/09.SV_calling/${prefix}/PAV/${sample}/config"
mkdir -p "${wd}/09.SV_calling/${prefix}/PAV/${sample}/input"
ln -sf "${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/${sample}.hap1.${prefix}.negative_strand_correction.fa" "${wd}/09.SV_calling/${prefix}/PAV/${sample}/input/${sample}.h1.fa"
ln -sf "${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/${sample}.hap2.${prefix}.negative_strand_correction.fa" "${wd}/09.SV_calling/${prefix}/PAV/${sample}/input/${sample}.h2.fa"
echo "{" > "${wd}/09.SV_calling/${prefix}/PAV/${sample}/config.json"
echo -e "\"reference\": \"${ref}\"," >> "${wd}/09.SV_calling/${prefix}/PAV/${sample}/config.json"
echo -e "\"asm_pattern\":\"${wd}/09.SV_calling/${prefix}/PAV/${sample}/input/{asm_name}.{hap}.fa\"" >> "${wd}/09.SV_calling/${prefix}/PAV/${sample}/config.json"
echo '}' >> "${wd}/09.SV_calling/${prefix}/PAV/${sample}/config.json"
echo -e "NAME\tHAP1\tHAP2" > "${wd}/09.SV_calling/${prefix}/PAV/${sample}/assemblies.tsv"
echo -e "${sample}\t${wd}/09.SV_calling/${prefix}/PAV/${sample}/input/${sample}.h1.fa\t${wd}/09.SV_calling/${prefix}/PAV/${sample}/input/${sample}.h2.fa" >> "${wd}/09.SV_calling/${prefix}/PAV/${sample}/assemblies.tsv"
pav_content="#!/bin/bash
#SBATCH --job-name=${sample}_pav
#SBATCH --partition=batch
#SBATCH --output=${wd}/09.SV_calling/${prefix}/PAV/${sample}/log/${sample}.pav.out
#SBATCH --error=${wd}/09.SV_calling/${prefix}/PAV/${sample}/log/${sample}.pav.err
#SBATCH --mincpus=${thread}
#SBATCH --mem=64g
source activate $snakemake_env
cd ${wd}/09.SV_calling/${prefix}/PAV/${sample}
$pav_cmd -c ${thread}
"
echo -e "$pav_content" > "${wd}/09.SV_calling/${prefix}/PAV/${sample}/run.pav.sh"
echo -e "sh ${wd}/09.SV_calling/${prefix}/PAV/${sample}/run.pav.sh" >> "$wd/run_scripts/09.pav.${prefix}.run.sh"
echo -e "sbatch ${wd}/09.SV_calling/${prefix}/PAV/${sample}/run.pav.sh" >> "$wd/run_scripts/09.pav.${prefix}.sbatch.sh"
done
done 9< "${samplefile}"

# Step.9-6 Run TRGT
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
for i in 1 3
do
prefix_var="prefix${i}"
ref_var="ref${i}"
trgt_var="trgt${i}"
prefix=${!prefix_var}
ref=${!ref_var}
trgt=${!trgt_var}
mkdir -p "${wd}/09.SV_calling/${prefix}/TRGT/${sample}/log"
trgt_content="#!/bin/bash
#SBATCH --job-name=${sample}_pav
#SBATCH --partition=batch
#SBATCH --output=${wd}/09.SV_calling/${prefix}/TRGT/${sample}/log/${sample}.pav.out
#SBATCH --error=${wd}/09.SV_calling/${prefix}/TRGT/${sample}/log/${sample}.pav.err
#SBATCH --mincpus=${thread}
#SBATCH --mem=64g
$trgt_cmd \
--genome ${ref} \
--repeats ${trgt} \
--reads ${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam \
--output-prefix ${wd}/09.SV_calling/${prefix}/TRGT/${sample}/${sample}
--threads ${thread} \
--karyotype XX \
--sample-name $sample
"
echo -e "$trgt_content" > "${wd}/09.SV_calling/${prefix}/TRGT/${sample}/run.trgt.sh"
echo -e "sh ${wd}/09.SV_calling/${prefix}/TRGT/${sample}/run.trgt.sh" >> "$wd/run_scripts/09.trgt.female.${prefix}.run.sh"
echo -e "sbatch ${wd}/09.SV_calling/${prefix}/TRGT/${sample}/run.trgt.sh" >> "$wd/run_scripts/09.trgt.female.${prefix}.sbatch.sh"
done
done 9< "${femalelist}"

# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
for i in 1 3
do
prefix_var="prefix${i}"
ref_var="ref${i}"
trgt_var="trgt${i}"
prefix=${!prefix_var}
ref=${!ref_var}
trgt=${!trgt_var}
mkdir -p "${wd}/09.SV_calling/${prefix}/TRGT/${sample}/log"
trgt_content="#!/bin/bash
#SBATCH --job-name=${sample}_pav
#SBATCH --partition=batch
#SBATCH --output=${wd}/09.SV_calling/${prefix}/TRGT/${sample}/log/${sample}.pav.out
#SBATCH --error=${wd}/09.SV_calling/${prefix}/TRGT/${sample}/log/${sample}.pav.err
#SBATCH --mincpus=${thread}
#SBATCH --mem=64g
$trgt_cmd \
--genome ${ref} \
--repeats ${trgt} \
--reads ${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam \
--output-prefix ${wd}/09.SV_calling/${prefix}/TRGT/${sample}/${sample}
--threads ${thread} \
--karyotype XY \
--sample-name $sample
"
echo -e "$trgt_content" > "${wd}/09.SV_calling/${prefix}/TRGT/${sample}/run.trgt.sh"
echo -e "sh ${wd}/09.SV_calling/${prefix}/TRGT/${sample}/run.trgt.sh" >> "$wd/run_scripts/09.trgt.male.${prefix}.run.sh"
echo -e "sbatch ${wd}/09.SV_calling/${prefix}/TRGT/${sample}/run.trgt.sh" >> "$wd/run_scripts/09.trgt.male.${prefix}.sbatch.sh"
done
done 9< "${malelist}"

########################END Step.9 Call SV ########################



######################## Step.10 Remap the unmappable contigs to virus genome(EBV & phix174) ########################
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
for i in {1..3}
do
prefix_var="prefix${i}"
ref_var="ref${i}"
prefix=${!prefix_var}
ref=${!ref_var}
mkdir -p "${wd}/10.EBV_contigs/${prefix}/${sample}/log"
ebv_script="#!/bin/bash
#SBATCH --job-name=${sample}_ebv
#SBATCH --partition=batch
#SBATCH --output=${wd}/10.EBV_contigs/${prefix}/${sample}/log/virus_scan.contig.stdout
#SBATCH --error=${wd}/10.EBV_contigs/${prefix}/${sample}/log/virus_scan.contig.stderr
#SBATCH --mincpus=8
#SBATCH --mem=32g
samtools view -f 4 ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sorted.bam > ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.unmap.sam
samtools fasta ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.unmap.sam > ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.unmap.fasta
rm ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.unmap.sam
samtools view -f 4 ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sorted.bam > ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.unmap.sam
samtools fasta ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.unmap.sam > ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.unmap.fasta
rm ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.unmap.sam

blastn -query ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.unmap.fasta -db ${virus_db} -out ${wd}/10.EBV_contigs/${prefix}/${sample}/${sample}.${prefix}.hap1.unmap.map2virus.blast.out -outfmt 6 -max_target_seqs 5 -max_hsps 5
blastn -query ${wd}/08.asm_align/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.unmap.fasta -db ${virus_db} -out ${wd}/10.EBV_contigs/${prefix}/${sample}/${sample}.${prefix}.hap2.unmap.map2virus.blast.out -outfmt 6 -max_target_seqs 5 -max_hsps 5
"
echo -e "${ebv_script}" > "${wd}/10.EBV_contigs/${prefix}/${sample}/run.align_ebv.cmd.sh"
echo -e "sh ${wd}/10.EBV_contigs/${prefix}/${sample}/run.align_ebv.cmd.sh" >> "$wd/run_scripts/10.ebv_contigs.run.sh"
echo -e "sbatch ${wd}/10.EBV_contigs/${prefix}/${sample}/run.align_ebv.cmd.sh" >> "$wd/run_scripts/10.ebv_contigs.sbatch.sh"
done
done 9< "${samplefile}"
########################END Step.10 Remap the unmappable contigs to virus genome(EBV & phix174) ########################

######################## Step.11 run assmebly improvement (assign and connect)########################################
# v0 remove rRNA and virus contigs
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
mkdir -p "${wd}/11.asm_update/assembly_rmrDNA_virus/rm_rDNA/${sample}/input"
indir="${wd}/03.hifiasm/"
asm1=${indir}/${sample}/${sample}.asm.bp.hap1.p_ctg.fa
asm2=${indir}/${sample}/${sample}.asm.bp.hap2.p_ctg.fa
fai1=${indir}/${sample}/${sample}.asm.bp.hap1.p_ctg.fa.fai
fai2=${indir}/${sample}/${sample}.asm.bp.hap2.p_ctg.fa.fai
gfa1=${indir}/${sample}/${sample}.asm.bp.hap1.p_ctg.gfa
gfa2=${indir}/${sample}/${sample}.asm.bp.hap2.p_ctg.gfa
awk -F'\t' '{print "path" "\t" $1 "\t" $1}' "$fai1" > "${wd}/11.asm_update/assembly_rmrDNA_virus/rm_rDNA/${sample}/input/${sample}.hap1.scfmap"
awk -F'\t' '{print "path" "\t" $1 "\t" $1}' "$fai2" > "${wd}/11.asm_update/assembly_rmrDNA_virus/rm_rDNA/${sample}/input/${sample}.hap2.scfmap"
>"${wd}/11.asm_update/assembly_rmrDNA_virus/rm_rDNA/${sample}/input/${sample}.hap1.hifi.csv"
>"${wd}/11.asm_update/assembly_rmrDNA_virus/rm_rDNA/${sample}/input/${sample}.hap2.hifi.csv"
## script
script="#!/bin/bash
#SBATCH --job-name=${sample}.asm_rm_rDNA
#SBATCH --partition=batch
#SBATCH --output=${wd}/11.asm_update/assembly_rmrDNA_virus/rm_rDNA/${sample}/rm_rDNA.log
#SBATCH --error=${wd}/11.asm_update/assembly_rmrDNA_virus/rm_rDNA/${sample}/rm_rDNA.err
#SBATCH --mincpus=4
#SBATCH --mem=16g
source activate $mashmap_env
cd ${wd}/11.asm_update/assembly_rmrDNA_virus/rm_rDNA/${sample}
perl $screen_assembly_pl \
    --assembly $asm1 \
    --graph $gfa1 \
    --graphmap ${wd}/11.asm_update/assembly_rmrDNA_virus/rm_rDNA/${sample}/input/${sample}.hap1.scfmap \
    --hifi-coverage ${wd}/11.asm_update/assembly_rmrDNA_virus/rm_rDNA/${sample}/input/${sample}.hap1.hifi.csv\
    --contaminant rdna $rdna_fa \
    --output ${sample}.hap1.rDNAremoved \
    --threads 4
perl $screen_assembly_pl \
    --assembly $asm2 \
    --graph $gfa2 \
    --graphmap ${wd}/11.asm_update/assembly_rmrDNA_virus/rm_rDNA/${sample}/input/${sample}.hap2.scfmap \
    --hifi-coverage ${wd}/11.asm_update/assembly_rmrDNA_virus/rm_rDNA/${sample}/input/${sample}.hap2.hifi.csv\
    --contaminant rdna $rdna_fa \
    --output ${sample}.hap2.rDNAremoved \
    --threads 4
"
echo -e "$script" > "${wd}/11.asm_update/assembly_rmrDNA_virus/rm_rDNA/${sample}/rm_rDNA.sh"
echo -e "sh ${wd}/11.asm_update/assembly_rmrDNA_virus/rm_rDNA/${sample}/rm_rDNA.sh" >> "$wd/run_scripts/11.rmrDNA.run_rDNA.sh"
echo -e "sbatch ${wd}/11.asm_update/assembly_rmrDNA_virus/rm_rDNA/rm_rDNA.sh" >> "$wd/run_scripts/11.rmrDNA.sbatch_rDNA.sh"


# remove virus contigs
hap1_fa="${wd}/11.asm_update/assembly_rmrDNA_virus/rm_rDNA/${sample}/${sample}.hap1.rDNAremoved.fasta"
hap2_fa="${wd}/11.asm_update/assembly_rmrDNA_virus/rm_rDNA/${sample}/${sample}.hap2.rDNAremoved.fasta"
outdir="${wd}/11.asm_update/assembly_rmrDNA_virus/rm_virus/${sample}"
mkdir -p "${outdir}/hap1"
mkdir -p "${outdir}/hap2"
# screen tax-id 9606 homo sapien
script="#!/bin/bash
#SBATCH --job-name=${sample}.virus_rm
#SBATCH --partition=batch
#SBATCH --output= ${outdir}/rm_virus.log
#SBATCH --error= ${outdir}/rm_virus.err
#SBATCH --mincpus=32
#SBATCH --mem=200g
python3 $fcs_py --image $fcs_gx \
screen genome --fasta ${hap1_fa} \
--out-dir ${outdir}/hap1 \
--gx-db $gxdb --tax-id 9606
python3 $fcs_py --image $fcs_gx \
screen genome --fasta ${hap2_fa} \
--out-dir ${outdir}/hap2 \
--gx-db $gxdb --tax-id 9606
# clean the genome
cd ${outdir}/hap1
cat ${hap1_fa} \
| python3 $fcs_py --image $fcs_gx \
clean genome --action-report ${outdir}/hap1/${sample}.hap1.rDNAremoved.9606.fcs_gx_report.txt \
--output ${sample}.1.9606.clean.fasta --contam-fasta-out ${sample}.1.9606.contam.fasta
cd ${outdir}/hap2
cat ${hap2_fa} \
| python3 $fcs_py --image $fcs_gx \
clean genome --action-report ${outdir}/hap2/${sample}.hap2.rDNAremoved.9606.fcs_gx_report.txt \
--output ${sample}.2.9606.clean.fasta --contam-fasta-out ${sample}.2.9606.contam.fasta

cd ${outdir}/hap1
/usr/bin/time -v ${minimap2_cmd} -ax asm20 --cs -r2k -t ${thread} \
${ref1} \
${sample}.1.9606.clean.fasta > ${sample}.hap1.chm13v2.minimap.sam
samtools view -bS -@${thread} ${sample}.hap1.chm13v2.minimap.sam -o ${sample}.hap1.chm13v2.minimap.bam
samtools sort -@${thread} ${sample}.hap1.chm13v2.minimap.bam -o ${sample}.hap1.chm13v2.minimap.sorted.bam
samtools index ${sample}.hap1.chm13v2.minimap.sorted.bam
rm ${sample}.hap1.chm13v2.minimap.sam
cd ${outdir}/hap2
/usr/bin/time -v ${minimap2_cmd} -ax asm20 --cs -r2k -t ${thread} \
${ref1} \
${sample}.2.9606.clean.fasta > ${sample}.hap2.chm13v2.minimap.sam
samtools view -bS -@${thread} ${sample}.hap2.chm13v2.minimap.sam -o ${sample}.hap2.chm13v2.minimap.bam
samtools sort -@${thread} ${sample}.hap2.chm13v2.minimap.bam -o ${sample}.hap2.chm13v2.minimap.sorted.bam
samtools index ${sample}.hap2.chm13v2.minimap.sorted.bam
rm ${sample}.hap2.chm13v2.minimap.sam
"
echo -e "$script" > "${outdir}/rm_virus.sh"
echo -e "sh ${outdir}/rm_virus.sh" >> "$wd/run_scripts/11.rmvirus.run.sh"
echo -e "sbatch ${outdir}/rm_virus.sh" >> "$wd/run_scripts/11.rmvirus.sbatch.sh"
done 9< "${samplefile}"


for i in {1..3}
do
prefix_var="prefix${i}"
ref_var="ref${i}"
cyto_var="cyto${i}"
gs_var="gs${i}"
gl_var="genome${i}"
prefix=${!prefix_var}
ref=${!ref_var}
cyto=${!cyto_var}
gs=${!gs_var}
gl=${!gl_var}

# Female assignment
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
# 1.1 run assign
#The script is used to assign assembled contigs longer than 100kb to chromosomes (excluding chrM).
# For male samples, all chrX contigs will be assigned to hap1, and all chrY contigs will be assigned to hap2.
mkdir -p "${wd}/11.asm_update/v1_assign/${prefix}/${sample}/log"
asm1=${wd}/11.asm_update/assembly_rmrDNA_virus/rm_virus/${sample}/hap1/${sample}.1.9606.clean.fasta
asm2=${wd}/11.asm_update/assembly_rmrDNA_virus/rm_virus/${sample}/hap2/${sample}.2.9606.clean.fasta
bam1=${wd}/11.asm_update/assembly_rmrDNA_virus/rm_virus/${sample}/hap1/${sample}.hap1.chm13v2.minimap.sorted.bam
bam2=${wd}/11.asm_update/assembly_rmrDNA_virus/rm_virus/${sample}/hap2/${sample}.hap2.chm13v2.minimap.sorted.bam
assign_script="#!bin/bash
#SBATCH --job-name=${sample}_asm_chrom_assigned
#SBATCH --partition=batch
#SBATCH --output=${wd}/11.asm_update/v1_assign/${prefix}/${sample}/log/asm_chrom_assigned.log
#SBATCH --error=${wd}/11.asm_update/v1_assign/${prefix}/${sample}/log/asm_chrom_assigned.err
#SBATCH --mincpus=4
#SBATCH --mem=16g
cd ${wd}/11.asm_update/v1_assign/${prefix}/${sample}
source activate $biotools_env
## Hap1
seqkit seq -M 100000 $asm1 -j 8 > ${sample}.asm.hap1.unassigned.fa
seqkit seq -m 100000 $asm1 -j 8 > ${sample}.asm.hap1.m100kb.fa
grep '^>' ${sample}.asm.hap1.m100kb.fa | sed 's/^>//' > ${sample}.hap1.m100kb_contigs.ids
> ${sample}.asm.hap1.fa
> ${sample}.asm.hap1.id_mapping.txt
for i in {1..22} X; do
samtools view -b -F 256 -@8 $bam1 chr\${i} > ${sample}.hap1.chr.bam
samtools fasta ${sample}.hap1.chr.bam > ${sample}.hap1.chr.fa
# Select sequences longer than 100kb Of contigs
seqkit seq -m 100000 ${sample}.hap1.chr.fa -j 8 > ${sample}.hap1.chr.m100kb.fa
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
seqkit seq -M 100000 $asm2 -j 8 > ${sample}.asm.hap2.unassigned.fa
seqkit seq -m 100000 $asm2 -j 8 > ${sample}.asm.hap2.m100kb.fa
grep '^>' ${sample}.asm.hap2.m100kb.fa | sed 's/^>//' > ${sample}.hap2.m100kb_contigs.ids
> ${sample}.asm.hap2.fa
> ${sample}.asm.hap2.id_mapping.txt
for i in {1..22} X; do
samtools view -b -F 256 -@8 $bam2 chr\${i} > ${sample}.hap2.chr.bam
samtools fasta ${sample}.hap2.chr.bam > ${sample}.hap2.chr.fa
# Select sequences longer than 100kb Of contigs
seqkit seq -m 100000 ${sample}.hap2.chr.fa -j 8 > ${sample}.hap2.chr.m100kb.fa
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
echo -e "${assign_script}" > "${wd}/11.asm_update/v1_assign/${prefix}/${sample}/run.assign.${prefix}.cmd.sh"
echo -e "sh ${wd}/11.asm_update/v1_assign/${prefix}/${sample}/run.assign.${prefix}.cmd.sh" >> "$wd/run_scripts/11.assign.${prefix}.run.sh"
echo -e "sbatch ${wd}/11.asm_update/v1_assign/${prefix}/${sample}/run.assign.${prefix}.cmd.sh" >> "$wd/run_scripts/11.assign.${prefix}.sbatch.sh"
done 9< "${femalelist}"

# Male assignment
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
mkdir -p "${wd}/11.asm_update/v1_assign/${prefix}/${sample}/log"
# Original and aligned assemblies: bam
asm1=${wd}/11.asm_update/assembly_rmrDNA_virus/rm_virus/${sample}/hap1/${sample}.1.9606.clean.fasta
asm2=${wd}/11.asm_update/assembly_rmrDNA_virus/rm_virus/${sample}/hap2/${sample}.2.9606.clean.fasta
bam1=${wd}/11.asm_update/assembly_rmrDNA_virus/rm_virus/${sample}/hap1/${sample}.hap1.chm13v2.minimap.sorted.bam
bam2=${wd}/11.asm_update/assembly_rmrDNA_virus/rm_virus/${sample}/hap2/${sample}.hap2.chm13v2.minimap.sorted.bam

assign_script="#!/bin/bash
#SBATCH --job-name=${sample}_asm_chrom_assigned
#SBATCH --partition=batch
#SBATCH --output=${wd}/11.asm_update/v1_assign/${prefix}/${sample}/log/asm_chrom_assigned.log
#SBATCH --error=${wd}/11.asm_update/v1_assign/${prefix}/${sample}/log/asm_chrom_assigned.err
#SBATCH --mincpus=4
#SBATCH --mem=16g

cd ${wd}/11.asm_update/v1_assign/${prefix}/${sample}
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
seqkit seq -m 100000 ${sample}.asm.chr\${i}.fa -j 8 > ${sample}.asm.chr\${i}.m100kb.fa
# From the original assembly: contigs ID Replace with h1-chr1-1 pattern
grep '^>' ${sample}.asm.chr\${i}.m100kb.fa | sed 's/^>//' > ids.tmp
hap_prefix_var=\"chr\${i}_prefix\"
hap_prefix=\${!hap_prefix_var}
awk -v i=\${i} -v hap_prefix=\"\${hap_prefix}\" 'BEGIN {OFS=\"\t\"} {print \$1, hap_prefix \"-chr\" i \"-\" NR}' ids.tmp > ${sample}.asm.chr\${i}.ids_mapping.txt
seqkit replace -p '^(\S+)' -r '{kv}\$2' -k ${sample}.asm.chr\${i}.ids_mapping.txt ${sample}.asm.chr\${i}.m100kb.fa -j 8 > ${sample}.asm.chr\${i}.renamed.fa
done
rm ${sample}.asm.chrX.bam ${sample}.asm.chrX.fa ${sample}.asm.chrX.m100kb.fa ${sample}.asm.chrX.sorted.bam ${sample}.asm.chrY.bam ${sample}.asm.chrY.fa ${sample}.asm.chrY.m100kb.fa ${sample}.asm.chrY.sorted.bam ${sample}.hap1.chrX.bam ${sample}.hap2.chrX.bam ${sample}.hap1.chrY.bam ${sample}.hap2.chrY.bam
# asm1
seqkit seq -M 100000 $asm1 -j 8 > ${sample}.asm.hap1.unassigned.fa
seqkit seq -m 100000 $asm1 -j 8 > ${sample}.asm.hap1.m100kb.fa
grep '^>' ${sample}.asm.hap1.m100kb.fa | sed 's/^>//' > ${sample}.hap1.m100kb_contigs.ids
> ${sample}.asm.hap1.fa
> ${sample}.asm.hap1.id_mapping.txt
for i in {1..22}; do
samtools view -b -F 256 -@8 $bam1 chr\${i} > ${sample}.hap1.chr.bam
samtools fasta ${sample}.hap1.chr.bam > ${sample}.hap1.chr.fa
# Select sequences longer than 100kb Of contigs
seqkit seq -m 100000 ${sample}.hap1.chr.fa -j 8 > ${sample}.hap1.chr.m100kb.fa
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
seqkit seq -M 100000 $asm2 -j 8 > ${sample}.asm.hap2.unassigned.fa
seqkit seq -m 100000 $asm2 -j 8 > ${sample}.asm.hap2.m100kb.fa
grep '^>' ${sample}.asm.hap2.m100kb.fa | sed 's/^>//' > ${sample}.hap2.m100kb_contigs.ids
> ${sample}.asm.hap2.fa
> ${sample}.asm.hap2.id_mapping.txt
for i in {1..22}; do
samtools view -b -F 256 -@8 $bam2 chr\${i} > ${sample}.hap2.chr.bam
samtools fasta ${sample}.hap2.chr.bam > ${sample}.hap2.chr.fa
# Select sequences longer than 100kb Of contigs
seqkit seq -m 100000 ${sample}.hap2.chr.fa -j 8 > ${sample}.hap2.chr.m100kb.fa
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
echo -e "${assign_script}" > "${wd}/11.asm_update/v1_assign/${prefix}/${sample}/run.assign.${prefix}.cmd.sh"
echo -e "sh ${wd}/11.asm_update/v1_assign/${prefix}/${sample}/run.assign.${prefix}.cmd.sh" >> "$wd/run_scripts/11.assign.${prefix}.run.sh"
echo -e "sbatch ${wd}/11.asm_update/v1_assign/${prefix}/${sample}/run.assign.${prefix}.cmd.sh" >> "$wd/run_scripts/11.assign.${prefix}.sbatch.sh"
done 9< "$malelist"

#1.2 run correct strand

# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
# Align original assemblies bam
bam1=${wd}/11.asm_update/assembly_rmrDNA_virus/rm_virus/${sample}/hap1/${sample}.hap1.chm13v2.minimap.sorted.bam
bam2=${wd}/11.asm_update/assembly_rmrDNA_virus/rm_virus/${sample}/hap2/${sample}.hap2.chm13v2.minimap.sorted.bam
outdir=${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}
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
seqkit replace -p '^(\S+)' -r '{kv}\$2' -k ${wd}/11.asm_update/v1_assign/${prefix}/${sample}/${sample}.asm.id_mapping.txt $outdir/${sample}.asm.chm13v2.minimap.f16.rc.fa --keep-key -j 8 > $outdir/${sample}.asm.chm13v2.minimap.f16.rc.renamed.fa
# Use chrom_assign Reverse-strand sequences in the updated assembly: contigs, Replace with the reverse complement
python $replace_fasta_cmd $outdir/${sample}.asm.${prefix}.minimap.f16.rc.renamed.fa ${wd}/11.asm_update/v1_assign/${prefix}/${sample}/${sample}.asm.hap1.fa $outdir/${sample}.hap1.${prefix}.negative_strand_correction.fa
python $replace_fasta_cmd $outdir/${sample}.asm.${prefix}.minimap.f16.rc.renamed.fa ${wd}/11.asm_update/v1_assign/${prefix}/${sample}/${sample}.asm.hap2.fa $outdir/${sample}.hap2.${prefix}.negative_strand_correction.fa
samtools faidx $outdir/${sample}.hap1.${prefix}.negative_strand_correction.fa
samtools faidx $outdir/${sample}.hap2.${prefix}.negative_strand_correction.fa
rm $outdir/${sample}.asm.${prefix}.minimap.f16.rc.renamed.fa
rm $outdir/${sample}.asm.${prefix}.minimap.f16.rc.fa
rm $outdir/${sample}.asm.${prefix}.minimap.f16.fa
rm $outdir/${sample}.hap1.${prefix}.minimap.f16.fa $outdir/${sample}.hap2.${prefix}.minimap.f16.fa
rm $outdir/${sample}.hap1.${prefix}.minimap.f16.sorted.bam $outdir/${sample}.hap2.${prefix}.minimap.f16.sorted.bam
"
echo -e "${corstrand_script}" > "${outdir}/run.negative_strand_correction.cmd.sh"
echo -e "sh ${outdir}/run.negative_strand_correction.cmd.sh" >> "$wd/run_scripts/11.negative_strand_correction.run.sh"
echo -e "sbatch ${outdir}/run.negative_strand_correction.cmd.sh" >> "$wd/run_scripts/11.negative_strand_correction.sbatch.sh"
done 9< "$samplefile"

#1.3 run plot
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
mkdir -p "${wd}/11.asm_update/plot_pdf/${prefix}/${sample}/log"
mkdir -p "${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed"
readsbam=${wd}/04.reads_align/${prefix}/${sample}/${sample}.${prefix}.pbmm2.sorted.bam
asm1bam=${wd}/11.asm_update/assembly_rmrDNA_virus/rm_virus/${sample}/hap1/${sample}.hap1.chm13v2.minimap.sorted.bam
asm2bam=${wd}/11.asm_update/assembly_rmrDNA_virus/rm_virus/${sample}/hap2/${sample}.hap2.chm13v2.minimap.sorted.bam

plot_script="#!bin/bash
#SBATCH --job-name=${sample}_asm_chrom_assigned
#SBATCH --partition=batch
#SBATCH --output=${wd}/11.asm_update/v1_assign/${prefix}/${sample}/log/asm_chrom_assigned.log
#SBATCH --error=${wd}/11.asm_update/v1_assign/${prefix}/${sample}/log/asm_chrom_assigned.err
#SBATCH --mincpus=8
#SBATCH --mem=16g

#######
samtools view -b -F 256 -@8 $asm1bam | bedtools bamtobed -i - > ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.bed
samtools view -b -F 256 -@8 $asm2bam | bedtools bamtobed -i - > ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.bed
# Based on 01.chrom_assign Generated id_mapping.txt, Replace bed Fourth column contig ID As h1-chr1-1 pattern
python ${plot_update_bed_identifiers} ${wd}/11.asm_update/v1_assign/${prefix}/${sample}/${sample}.asm.id_mapping.txt ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.bed ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.identifiers_updated.bed
python ${plot_update_bed_identifiers} ${wd}/11.asm_update/v1_assign/${prefix}/${sample}/${sample}.asm.id_mapping.txt ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.bed ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.identifiers_updated.bed
# Remove bed unchanged entries contig ID alignment records
awk '\$4 ~ /-/' ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.identifiers_updated.bed > ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.identifiers.bed
awk '\$4 ~ /-/' ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.identifiers_updated.bed > ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.identifiers.bed
# Reassign h1-chrX/h2-chrY alignment records
awk '\$4 ~ /h1/' ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.identifiers.bed > ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.temp_h1_chrX_lines.bed
awk '\$4 ~ /h2/' ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.identifiers.bed > ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.temp_h2_chrY_lines.bed
cat ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.temp_h1_chrX_lines.bed >> ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.identifiers.bed
cat ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.temp_h2_chrY_lines.bed >> ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.identifiers.bed
awk '\$4 !~ /h2/' ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.identifiers.bed > ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.final.bed
awk '\$4 !~ /h1/' ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.identifiers.bed > ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.final.bed
rm ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.bed
rm ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.bed
rm ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.identifiers_updated.bed
rm ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.identifiers_updated.bed
rm ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.identifiers.bed
rm ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.identifiers.bed
rm ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.temp_h1_chrX_lines.bed
rm ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.temp_h2_chrY_lines.bed
#######

python ${plot_process_bed_cmd} ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.final.bed ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.processed.bed
python ${plot_process_bed_identifier_cmd} ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.processed.bed ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.identifiers_filled.bed ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.multi_identifier.bed
python ${plot_process_bed_cmd} ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.final.bed ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.processed.bed
python ${plot_process_bed_identifier_cmd} ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.processed.bed ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.identifiers_filled.bed ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.multi_identifier.bed
# Remove identifiers Column, Sort bed File, Use bedtools Merge overlapping intervals
awk -F'\\t' -v OFS='\\t' '{print \$1, \$2, \$3}' ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.identifiers_filled.bed | sort -k1,1V -k2,2n -k3,3n > ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.identifiers_filled.merged.temp.bed
awk -F'\\t' -v OFS='\\t' '{print \$1, \$2, \$3}' ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.identifiers_filled.bed | sort -k1,1V -k2,2n -k3,3n > ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.identifiers_filled.merged.temp.bed
bedtools merge -i ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.identifiers_filled.merged.temp.bed > ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.identifiers_filled.merged.bed
bedtools merge -i ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.identifiers_filled.merged.temp.bed > ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.identifiers_filled.merged.bed
rm ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.identifiers_filled.merged.temp.bed
rm ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.identifiers_filled.merged.temp.bed

# Generate assemblies gap bed File
bedtools complement -i ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.identifiers_filled.merged.bed -g ${gl} > ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap1.${prefix}.gap.bed
bedtools complement -i ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.identifiers_filled.merged.bed -g ${gl} > ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/bed/${sample}.hap2.${prefix}.gap.bed
source activate ${r4_env}
Rscript ${plot_karyotype_cmd} --wd ${wd} --prefix ${prefix} --sample ${sample} --cyto ${cyto} --gs ${gs} --readsbam ${readsbam}
"
echo -e "$plot_script" > "${wd}/11.asm_update/plot_pdf/${prefix}/${sample}/run.plot.${sample}.${prefix}.cmd.sh"
echo -e "sh ${wd}/11.asm_update/plot_pdf/${prefix}/${sample}/run.plot.${sample}.${prefix}.cmd.sh" >> "$wd/run_scripts/11.plot_assign.${prefix}.run.sh"
echo -e "sbatch ${wd}/11.asm_update/plot_pdf/${prefix}/${sample}/run.plot.${sample}.${prefix}.cmd.sh" >> "$wd/run_scripts/11.plot_assign.${prefix}.sbatch.sh"
done 9< "$samplefile"


#1.4 run connect
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
mkdir -p "${wd}/11.asm_update/v3_connect/${prefix}/${sample}/log"
connect_script="#!bin/bash
#SBATCH --job-name=${sample}_asm_chrom_assigned
#SBATCH --partition=batch
#SBATCH --output=${wd}/11.asm_update/v3_connect/${prefix}/${sample}/log/asm_chrom_assigned.log
#SBATCH --error=${wd}/11.asm_update/v3_connect/${prefix}/${sample}/log/asm_chrom_assigned.err
#SBATCH --mincpus=48
#SBATCH --mem=48g
source activate $biotools_env
python $connect_cmd --fasta_file ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/${sample}.hap1.${prefix}.negative_strand_correction.fa \
                        --fai_file ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/${sample}.hap1.${prefix}.negative_strand_correction.fa.fai \
                        --total_dir ${wd}/11.asm_update/v3_connect/${prefix}/${sample}/hap1

python $connect_cmd --fasta_file ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/${sample}.hap2.${prefix}.negative_strand_correction.fa \
                        --fai_file ${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/${sample}.hap2.${prefix}.negative_strand_correction.fa.fai \
                        --total_dir ${wd}/11.asm_update/v3_connect/${prefix}/${sample}/hap2

# Merge separately hap1 And hap2 Of contigs
cat ${wd}/11.asm_update/v3_connect/${prefix}/${sample}/hap1/chr*/*_merged.final.fa > ${wd}/11.asm_update/v3_connect/${prefix}/${sample}/${sample}.hap1.connect.fa
samtools faidx ${wd}/11.asm_update/v3_connect/${prefix}/${sample}/${sample}.hap1.connect.fa
cat ${wd}/11.asm_update/v3_connect/${prefix}/${sample}/hap2/chr*/*_merged.final.fa > ${wd}/11.asm_update/v3_connect/${prefix}/${sample}/${sample}.hap2.connect.fa
samtools faidx ${wd}/11.asm_update/v3_connect/${prefix}/${sample}/${sample}.hap2.connect.fa
"
echo -e "${connect_script}" > "${wd}/11.asm_update/v3_connect/${prefix}/${sample}/run.connect.${sample}.${prefix}.cmd.sh"
echo -e "sh ${wd}/11.asm_update/v3_connect/${prefix}/${sample}/run.connect.${sample}.${prefix}.cmd.sh" >> "$wd/run_scripts/11.connect.${prefix}.run.sh"
echo -e "sbatch ${wd}/11.asm_update/v3_connect/${prefix}/${sample}/run.connect.${sample}.${prefix}.cmd.sh" >> "$wd/run_scripts/11.connect.${prefix}.sbatch.sh"
done 9< "$samplefile"

#1.5 run align
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
mkdir -p "${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/log"
asm_align_script="#!/bin/bash
#SBATCH --job-name=m_$sample
#SBATCH --partition=batch
#SBATCH --output=${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/log/${sample}.minimap2.stdout
#SBATCH --error=${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/log/${sample}.minimap2.stderr
#SBATCH --mincpus=${thread}
#SBATCH --mem=48g
/usr/bin/time -v ${minimap2_cmd} -ax asm20 --cs -r2k -t ${thread} \
${ref} \
${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/${sample}.hap1.${prefix}.negative_strand_correction.fa > ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sam
samtools view -bS -@${thread} ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sam -o ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.bam
samtools sort -@${thread} ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.bam -o ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sorted.bam
samtools index ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sorted.bam
/usr/bin/time -v ${minimap2_cmd} -ax asm20 --cs -r2k -t ${thread} \
${ref} \
${wd}/11.asm_update/v2_corstrand/${prefix}/${sample}/${sample}.hap2.${prefix}.negative_strand_correction.fa > ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sam
samtools view -bS -@${thread} ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sam -o ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.bam
samtools sort -@${thread} ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.bam -o ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sorted.bam
samtools index ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sorted.bam
rm ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.sam
rm ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.sam
rm ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap1.${prefix}.minimap.bam
rm ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/${sample}.hap2.${prefix}.minimap.bam
"
echo -e "$asm_align_script" > "${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/run.align.cmd.sh"
echo -e "sh ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/run.align.cmd.sh" >> "$wd/run_scripts/11.asm_align.${prefix}.run.sh"
echo -e "sbatch ${wd}/11.asm_update/align_to_ref/${prefix}/${sample}/run.align.cmd.sh" >> "$wd/run_scripts/11.asm_align.${prefix}.sbatch.sh"
done 9< "$samplefile"
done
########################END Step.11 run assmebly improvement (inspector and assign)########################################
