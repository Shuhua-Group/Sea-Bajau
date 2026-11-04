#!/bin/bash

# Analyze digestive-enzyme copy numbers and haplotypes.
# Set input, reference and tool paths for your environment before execution.

cd /path/to/PGA_60_48_8_36 || exit 1

conda activate pggb_env
pggb \
  -i PGA_region_seq.fa.gz \
  -o PGA_60_48_8_36.pggb \
  -p 90 -s 5k -k 419 -t 16 -n 367 \
  -Z

pggb \
  -i PGA_region_seq.fa.gz \
  -o PGA_60_48_8_36.k19.pggb \
  -p 90 -s 5k -k 19 -t 16 -n 367 \
  -Z




# Test one sample

wd="/path/to/2018cell.k19"
cd "$wd" || exit 1
ln -s ../PGA_60_48_8_36.k19* .
# ln -s ../PGAvgindex* .

# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r i <&9 || [[ -n "$i" ]]; do
    [[ -z "$i" ]] && continue
sample=$i
sample_bam=$(ls /path/to/${sample}.sort.rmdup.realign.md.bam)
threads=12
ref="/path/to/human_g1k_v37.1_22XYMT.fa"
# hg38="/path/to/Homo_sapiens.GRCh38.dna.chromosome.1_22_X_Y.fa"
region="11:60930401-61099001"
samtools="/path/to/samtools"
graph="/path/to/PGA_region_seq.fa.gz.bf3285f.37c1a4f.8635317.smooth.final.gfa"
odgi="/path/to/odgi"
gfapack="/path/to/gafpack"
cosigt="/path/to/cosigt"
python="/path/to/python"
tprpy="/path/to/tpr.py"
gfainject="/path/to/gfainject"

# Note: GRCh37 coordinates chr11:60930401-61099001 [ Original GRCh38 coordinates: chr11:61162929-61331529]

# Step 1: Extract from the specified region: BAM File
"${samtools}" view -O bam \
-o "${sample}.PGA.bam" \
-T "$ref" \
-@ "${threads}" \
${sample_bam} \
${region}

# Step 2: From BAM Extract from the file fastq/fastq

"${samtools}" fastq \
-@ "${threads}" \
${sample}.PGA.bam | bgzip > "${sample}.PGA.fastq.gz"


"${samtools}" fasta \
-@ "${threads}" \
${sample}.PGA.bam | bgzip > "${sample}.PGA.fasta.gz"


# # Step 3: odgi #[ Run once ]
# ${odgi} chop \
# -i ${graph} \
# -c 32 -t ${threads} \
# -o - | odgi view \
# -i - \
# -g > PGA_60_48_8_36.k19.pggb.odgi.z.gfa

# ${odgi} build \
# -g PGA_60_48_8_36.k19.pggb.odgi.z.gfa \
# -o - | odgi paths \
# -i - \
# -H | cut -f 1,4- | gzip > PGA_60_48_8_36.k19.pggb.odgi.z.paths.tsv.gz

# ${odgi} paths \
# -i ${graph} \
# -f > PGA_60_48_8_36.k19.pggb.odgi.z.fa


# Step 4: bwa
# # bwa index PGA_60_48_8_36.k19.pggb.odgi.z.fa #[ Run once ]
bwa mem -t "${threads}" PGA_60_48_8_36.k19.pggb.odgi.z.fa "${sample}.PGA.fasta.gz" | samtools sort > "${sample}.PGA.realigned.bam"

# Step 5: vg
# # mkdir PGAvgindex #[ Run once ]
# # vg \
# # autoindex \
# # -w giraffe \
# # -g PGA_60_48_8_36.k19.pggb.odgi.z.gfa \
# # -t ${threads} \
# # -p PGAvgindex \
# # -T PGAvgindex


# vg \
# giraffe \
# -Z PGAvgindex.giraffe.gbz \
# -m PGAvgindex.min \
# -d PGAvgindex.dist \
# -f ${sample}.PGA.fastq.gz \
# -o gaf > ${sample}.x.gaf

# Step 6: gfainject [ Use either this step or step 5 ]

"$gfainject" --gfa PGA_60_48_8_36.k19.pggb.odgi.z.gfa --bam "${sample}.PGA.realigned.bam" > "${sample}.x.gaf"

# # Step 7: gafpack
"${gfapack}" \
--gfa PGA_60_48_8_36.k19.pggb.odgi.z.gfa \
--gaf "${sample}.x.gaf" | bgzip > "${sample}.x.gafpack.gz"

# # Step 8: cosigt

"${cosigt}" -p PGA_60_48_8_36.k19.pggb.odgi.z.paths.tsv.gz \
-g "${sample}.x.gafpack.gz" \
-o "${sample}.PGA_GT" \
-i "${sample}"

# ## Evaluate
zcat "${sample}.PGA_GT/${sample}.PGA_GT.sorted_combos.tsv.gz" |\
sort -k 2 -n -r  | \
awk -v var="${sample}" -F '${sample}' '{print NR "\t" NF-1 "\t" var "\t" $0}' | \
sed 's/haplotype1-/haplotype1_/g' | sed 's/haplotype2-/haplotype2_/g' | \
sed -E 's/([0-9])\-([0-9])/\1_\2/g' | \
tr '-' '\t' | \
sed 's/haplotype1_/haplotype1-/g' | \
sed 's/haplotype2_/haplotype2-/g' \
> "${sample}.PGA_GT/${sample}.PGA_GT.sorted_combos.evaluation.tsv"
# # Plot

# sed -i "s/hap1#h1\t/hap1#h1_/g" ${sample}.PGA_GT/${sample}.PGA_GT.sorted_combos.evaluation.tsv
# sed -i "s/hap2#h2\t/hap2#h2_/g" ${sample}.PGA_GT/${sample}.PGA_GT.sorted_combos.evaluation.tsv
done 9< <(cat -- "../2018cell.list" | head -1)
# $python ${tprpy} ${sample}.PGA_GT/ \
# "" \
# ${sample}.PGA_GT/${sample}.PGA_GT.sorted_combos.evaluation.pdf


# Check whether cell Review source data quality
# Try using HGDP samples

wd="/path/to/HGDP"
mkdir -p "$wd"
cd "$wd" || exit 1
ln -s ../PGA_60_48_8_36.k19* .
ln -s ../PGAvgindex* .


for i in HGDP00547;do
sample=$i
sample_bam=$(ls /path/to/${sample}.dedup.bam)
threads=12
# ref="/path/to/human_g1k_v37.1_22XYMT.fa"
hg38="/path/to/Homo_sapiens.GRCh38.dna.chromosome.1_22_X_Y.fa"
region="chr11:61162929-61331529"
samtools="/path/to/samtools"
graph="/path/to/PGA_region_seq.fa.gz.bf3285f.e219c8e.8635317.smooth.final.gfa"
odgi="/path/to/odgi"
gfapack="/path/to/gafpack"
cosigt="/path/to/cosigt"
python="/path/to/python"
tprpy="/path/to/tpr.py"
gfainject="/path/to/gfainject"

# Note: GRCh37 coordinates chr11:60930401-61099001 [ Original GRCh38 coordinates: chr11:61162929-61331529]

# Step 1: Extract from the specified region: BAM File
# ${samtools} view -O bam \
# -o ${sample}.PGA.bam \
# -T $ref \
# -@ ${threads} \
# ${sample_bam} \
# ${region}

# Step 2: From BAM Extract from the file fastq/fastq

# ${samtools} fastq \
# -@ ${threads} \
# ${sample}.PGA.bam | bgzip > ${sample}.PGA.fastq.gz


# ${samtools} fasta \
# -@ ${threads} \
# ${sample}.PGA.bam | bgzip > ${sample}.PGA.fasta.gz


# # Step 3: odgi #[ Run once ]
# # ${odgi} chop \
# # -i ${graph} \
# # -c 32 -t ${threads} \
# # -o - | odgi view \
# # -i - \
# # -g > PGA_60_48_8_36.k19.pggb.odgi.z.gfa

# # ${odgi} build \
# # -g PGA_60_48_8_36.k19.pggb.odgi.z.gfa \
# # -o - | odgi paths \
# # -i - \
# # -H | cut -f 1,4- | gzip > PGA_60_48_8_36.k19.pggb.odgi.z.paths.tsv.gz

# # ${odgi} paths \
# # -i ${graph} \
# # -f > PGA_60_48_8_36.k19.pggb.odgi.z.fa


# Step 4: bwa
# # bwa index PGA_60_48_8_36.k19.pggb.odgi.z.fa #[ Run once ]
# bwa mem -t ${threads} PGA_60_48_8_36.k19.pggb.odgi.z.fa ${sample}.PGA.fasta.gz | samtools sort > ${sample}.PGA.realigned.bam

# Step 5: vg
# # mkdir PGAvgindex #[ Run once ]
# # vg \
# # autoindex \
# # -w giraffe \
# # -g PGA_60_48_8_36.k19.pggb.odgi.z.gfa \
# # -t ${threads} \
# # -p PGAvgindex \
# # -T PGAvgindex


# vg \
# giraffe \
# -Z PGAvgindex.giraffe.gbz \
# -m PGAvgindex.min \
# -d PGAvgindex.dist \
# -f ${sample}.PGA.fastq.gz \
# -o gaf > ${sample}.x.gaf

# Step 6: gfainject [ Use either this step or step 5 ]

"$gfainject" --gfa PGA_60_48_8_36.k19.pggb.odgi.z.gfa --bam "${sample}.PGA.realigned.bam" > "${sample}.x.gaf"

# # Step 7: gafpack
"${gfapack}" \
--gfa PGA_60_48_8_36.k19.pggb.odgi.z.gfa \
--gaf "${sample}.x.gaf" | bgzip > "${sample}.x.gafpack.gz"

# # Step 8: cosigt

"${cosigt}" -p PGA_60_48_8_36.k19.pggb.odgi.z.paths.tsv.gz \
-g "${sample}.x.gafpack.gz" \
-o "${sample}.PGA_GT" \
-i "${sample}"

# ## Evaluate
zcat "${sample}.PGA_GT/${sample}.PGA_GT.sorted_combos.tsv.gz" |\
sort -k 2 -n -r  | \
awk -v var="${sample}" -F '${sample}' '{print NR "\t" NF-1 "\t" var "\t" $0}' | \
sed 's/haplotype1-/haplotype1_/g' | sed 's/haplotype2-/haplotype2_/g' | \
sed -E 's/([0-9])\-([0-9])/\1_\2/g' | \
tr '-' '\t' | \
sed 's/haplotype1_/haplotype1-/g' | \
sed 's/haplotype2_/haplotype2-/g' \
> "${sample}.PGA_GT/${sample}.PGA_GT.sorted_combos.evaluation.tsv"
# # Plot

sed -i "s/hap1#h1\t/hap1#h1_/g" "${sample}.PGA_GT/${sample}.PGA_GT.sorted_combos.evaluation.tsv"
sed -i "s/hap2#h2\t/hap2#h2_/g" "${sample}.PGA_GT/${sample}.PGA_GT.sorted_combos.evaluation.tsv"

done

