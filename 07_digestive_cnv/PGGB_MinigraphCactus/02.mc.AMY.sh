#!/bin/bash

# Analyze digestive-enzyme copy numbers and haplotypes.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/mc"

cd "$wd" || exit 1

MYJOBSTORE=$wd/js_cactus
MYBUCKET=$wd/MC_graph/
mkdir -p "$MYBUCKET"

PREFIX=AMY_80_56_4_64.MC

seqfile=$wd/seqfile.txt

######## seqfile Prepare again #########
touch temp.list
mkdir asm
for contig in `cut -f1 ../AMY_region_seq.fa.fai`;do
# NA19338_hap2#NA19338#2#CM087780.1_60597798_60769095_0
# Extract sample names and hap Name
sample=`echo $contig | cut -d "#" -f1|cut -d "p" -f1`
hap=`echo $contig | cut -d "#" -f1| cut -d "p" -f2`
# echo $sample $hap
echo "$sample" "$hap" >> temp.list
if [ ! -f "asm/${sample}p${hap}.fa" ]; then
# echo $sample $hap
samtools faidx ../AMY_region_seq.fa.gz "$contig" > "asm/${sample}p${hap}.fa"
fi
done

sort temp.list |uniq > temp2.list
mv temp2.list temp.list


cd asm || exit 1

for i in `ls *.fa`;do
# Replace the first row
name=$(basename $i .fa)
echo "$name" > temp.id
sed -i "s/.hap/_hap/g" temp.id
rename=$(cat temp.id)
sed "1s/.*/>${rename}/" "$i" > "${rename}.rename.fa"
done



for sample in `cut -d " " -f1 temp.list | sort |uniq`;do
if [ "$(grep -w "$sample" temp.list | wc -l)" -eq 2 ]; then
echo "$sample" >> diploid.list
elif [ "$(grep -w "$sample" temp.list | wc -l)" -eq 1 ]; then
echo "$sample" >> haploid.list
fi
done

rm temp.list
# Use haploid.list within chm13v2 And GRCh38 Remove
sed -i '/chm13v2/d' haploid.list
sed -i '/Homo/d' haploid.list


sed -i "s/_ha//g" haploid.list
sed -i "s/.ha//g" haploid.list
sed -i "s/_ha//g" diploid.list
sed -i "s/.ha//g" diploid.list

sort diploid.list|uniq > diploid.list2
sort haploid.list|uniq > haploid.list2
mv diploid.list2 diploid.list
mv haploid.list2 haploid.list


echo -e "GRCh38AMY\t${wd}/asm/GRCh38.AMY.rename.fa" > seqfile.txt
echo -e "CHM13AMY\t${wd}/asm/chm13v2.AMY.rename.fa" >> seqfile.txt

# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
hap1=$(ls ${wd}/asm/${sample}*hap1*.rename.fa)
echo -e "${sample}.1\t${hap1}" >> seqfile.txt
hap2=$(ls ${wd}/asm/${sample}*hap2*.rename.fa)
echo -e "${sample}.2\t${hap2}" >> seqfile.txt
done 9< "diploid.list"


# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
hap=$(ls ${wd}/asm/${sample}*hap*rename.fa)
echo -e "${sample}\t${hap}" >> seqfile.txt
done 9< "haploid.list"



logfile=$wd/cactus_v293.log

source /path/to/activate
mkdir -p "$wd/temp"
export TMPDIR=$wd/temp

cactus-pangenome "${MYJOBSTORE}" "${seqfile}" --outName "${PREFIX}" --outDir "${MYBUCKET}" \
--reference GRCh38AMY CHM13AMY --filter 35 \
--lastTrain --snarlStats --haplo --maxMemory 200G \
--giraffe full clip filter \
--vcf clip filter full \
--vcfReference GRCh38AMY CHM13AMY \
--gbz clip filter full \
--xg clip filter full \
--chrom-vg clip filter full \
--viz --odgi --chrom-vg clip filter \
--chrom-og --gbz clip filter full \
--gfa clip full filter \
--vcfwave \
--logFile "$logfile" \
--mgCores 32 --mapCores 32 --consCores 32 --indexCores 32





# Test one sample

wd="/path/to/2018cell"
cd "$wd" || exit 1
ln -s ../AMY_80_56_4_64* .
ln -s ../AMYvgindex* .

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
graph="/path/to/AMY_80_56_4_64.MC.full.gfa"
odgi="/path/to/odgi"
gfapack="/path/to/gafpack"
cosigt="/path/to/cosigt"
python="/path/to/python"
tprpy="/path/to/tpr.py"
gfainject="/path/to/gfainject"

# Note: GRCh37 coordinates chr11:60930401-61099001 [ Original GRCh38 coordinates: chr11:61162929-61331529]

# Step 1: Extract from the specified region: BAM File
"${samtools}" view -O bam \
-o "${sample}.AMY.bam" \
-T "$ref" \
-@ "${threads}" \
${sample_bam} \
${region}

# Step 2: From BAM Extract from the file fastq/fastq

"${samtools}" fastq \
-@ "${threads}" \
${sample}.AMY.bam | pigz > "${sample}.AMY.fastq.gz"


"${samtools}" fasta \
-@ "${threads}" \
${sample}.AMY.bam | pigz > "${sample}.AMY.fasta.gz"


# # Step 3: odgi #[ Run once ]
"${odgi}" chop  \
-i "${graph}" \
-c 32 -t "${threads}" \
-o - | "${odgi}" view \
-i - \
-g > AMY_80_56_4_64.mc.odgi.z.gfa

"${odgi}" build \
-g AMY_80_56_4_64.mc.odgi.z.gfa \
-o - | "${odgi}" paths \
-i - \
-H | cut -f 1,4- | gzip > AMY_80_56_4_64.mc.odgi.z.paths.tsv.gz

"${odgi}" paths \
-i "${graph}" \
-f > AMY_80_56_4_64.mc.odgi.z.fa


# Step 4: bwa
bwa index AMY_80_56_4_64.mc.odgi.z.fa #[ Run once ]
# bwa mem -t ${threads} AMY_80_56_4_64.mc.odgi.z.fa ${sample}.AMY.fasta.gz | samtools sort > ${sample}.AMY.realigned.bam

# Step 5: vg
# # mkdir AMYvgindex #[ Run once ]
# # vg \
# # autoindex \
# # -w giraffe \
# # -g AMY_80_56_4_64.mc.odgi.z.gfa \
# # -t ${threads} \
# # -p AMYvgindex \
# # -T AMYvgindex


# vg \
# giraffe \
# -Z AMYvgindex.giraffe.gbz \
# -m AMYvgindex.min \
# -d AMYvgindex.dist \
# -f ${sample}.AMY.fastq.gz \
# -o gaf > ${sample}.x.gaf

# Step 6: gfainject [ Use either this step or step 5 ]

"$gfainject" --gfa AMY_80_56_4_64.mc.odgi.z.gfa --bam "${sample}.AMY.realigned.bam" > "${sample}.x.gaf"

# # Step 7: gafpack
"${gfapack}" \
--gfa AMY_80_56_4_64.mc.odgi.z.gfa \
--gaf "${sample}.x.gaf" | pigz > "${sample}.x.gafpack.gz"

# # Step 8: cosigt

"${cosigt}" -p AMY_80_56_4_64.mc.odgi.z.paths.tsv.gz \
-g "${sample}.x.gafpack.gz" \
-o "${sample}.AMY_GT" \
-i "${sample}"

# ## Evaluate
zcat "${sample}.AMY_GT/${sample}.AMY_GT.sorted_combos.tsv.gz" |\
sort -k 2 -n -r  | \
awk -v var="${sample}" -F '${sample}' '{print NR "\t" NF-1 "\t" var "\t" $0}' | \
sed 's/haplotype1-/haplotype1_/g' | sed 's/haplotype2-/haplotype2_/g' | \
sed -E 's/([0-9])\-([0-9])/\1_\2/g' | \
tr '-' '\t' | \
sed 's/haplotype1_/haplotype1-/g' | \
sed 's/haplotype2_/haplotype2-/g' \
> "${sample}.AMY_GT/${sample}.AMY_GT.sorted_combos.evaluation.tsv"
# # Plot

sed -i "s/hap1#h1\t/hap1#h1_/g" "${sample}.AMY_GT/${sample}.AMY_GT.sorted_combos.evaluation.tsv"
sed -i "s/hap2#h2\t/hap2#h2_/g" "${sample}.AMY_GT/${sample}.AMY_GT.sorted_combos.evaluation.tsv"

done 9< "../2018cell.list"
# $python ${tprpy} ${sample}.AMY_GT/ \
# "" \
# ${sample}.AMY_GT/${sample}.AMY_GT.sorted_combos.evaluation.pdf


# Check whether cell Review source data quality
# Try using HGDP samples

wd="/path/to/HGDP"
mkdir -p "$wd"
cd "$wd" || exit 1
ln -s ../AMY_80_56_4_64* .
ln -s ../AMYvgindex* .


for i in HGDP00547;do
sample=$i
sample_bam=$(ls /path/to/${sample}.dedup.bam)
threads=12
# ref="/path/to/human_g1k_v37.1_22XYMT.fa"
hg38="/path/to/Homo_sapiens.GRCh38.dna.chromosome.1_22_X_Y.fa"
region="chr11:61162929-61331529"
samtools="/path/to/samtools"
graph="/path/to/AMY_region_seq.fa.gz.bf3285f.e219c8e.8635317.smooth.final.gfa"
odgi="/path/to/odgi"
gfapack="/path/to/gafpack"
cosigt="/path/to/cosigt"
python="/path/to/python"
tprpy="/path/to/tpr.py"
gfainject="/path/to/gfainject"

# Note: GRCh37 coordinates chr11:60930401-61099001 [ Original GRCh38 coordinates: chr11:61162929-61331529]

# Step 1: Extract from the specified region: BAM File
# ${samtools} view -O bam \
# -o ${sample}.AMY.bam \
# -T $ref \
# -@ ${threads} \
# ${sample_bam} \
# ${region}

# Step 2: From BAM Extract from the file fastq/fastq

# ${samtools} fastq \
# -@ ${threads} \
# ${sample}.AMY.bam | pigz > ${sample}.AMY.fastq.gz


# ${samtools} fasta \
# -@ ${threads} \
# ${sample}.AMY.bam | pigz > ${sample}.AMY.fasta.gz


# # Step 3: odgi #[ Run once ]
"${odgi}" chop  \
-i "${graph}" \
-c 32 -t "${threads}" \
-o - | odgi view \
-i - \
-g > AMY_80_56_4_64.mc.odgi.z.gfa

"${odgi}" build \
-g AMY_80_56_4_64.mc.odgi.z.gfa \
-o - | odgi paths \
-i - \
-H | cut -f 1,4- | gzip > AMY_80_56_4_64.mc.odgi.z.paths.tsv.gz

"${odgi}" paths \
-i "${graph}" \
-f > AMY_80_56_4_64.mc.odgi.z.fa


# Step 4: bwa
# # bwa index AMY_80_56_4_64.mc.odgi.z.fa #[ Run once ]
# bwa mem -t ${threads} AMY_80_56_4_64.mc.odgi.z.fa ${sample}.AMY.fasta.gz | samtools sort > ${sample}.AMY.realigned.bam

# Step 5: vg
# # mkdir AMYvgindex #[ Run once ]
# # vg \
# # autoindex \
# # -w giraffe \
# # -g AMY_80_56_4_64.mc.odgi.z.gfa \
# # -t ${threads} \
# # -p AMYvgindex \
# # -T AMYvgindex


# vg \
# giraffe \
# -Z AMYvgindex.giraffe.gbz \
# -m AMYvgindex.min \
# -d AMYvgindex.dist \
# -f ${sample}.AMY.fastq.gz \
# -o gaf > ${sample}.x.gaf

# Step 6: gfainject [ Use either this step or step 5 ]

"$gfainject" --gfa AMY_80_56_4_64.mc.odgi.z.gfa --bam "${sample}.AMY.realigned.bam" > "${sample}.x.gaf"

# # Step 7: gafpack
"${gfapack}" \
--gfa AMY_80_56_4_64.mc.odgi.z.gfa \
--gaf "${sample}.x.gaf" | pigz > "${sample}.x.gafpack.gz"

# # Step 8: cosigt

"${cosigt}" -p AMY_80_56_4_64.mc.odgi.z.paths.tsv.gz \
-g "${sample}.x.gafpack.gz" \
-o "${sample}.AMY_GT" \
-i "${sample}"

# ## Evaluate
zcat "${sample}.AMY_GT/${sample}.AMY_GT.sorted_combos.tsv.gz" |\
sort -k 2 -n -r  | \
awk -v var="${sample}" -F '${sample}' '{print NR "\t" NF-1 "\t" var "\t" $0}' | \
sed 's/haplotype1-/haplotype1_/g' | sed 's/haplotype2-/haplotype2_/g' | \
sed -E 's/([0-9])\-([0-9])/\1_\2/g' | \
tr '-' '\t' | \
sed 's/haplotype1_/haplotype1-/g' | \
sed 's/haplotype2_/haplotype2-/g' \
> "${sample}.AMY_GT/${sample}.AMY_GT.sorted_combos.evaluation.tsv"
# # Plot

sed -i "s/hap1#h1\t/hap1#h1_/g" "${sample}.AMY_GT/${sample}.AMY_GT.sorted_combos.evaluation.tsv"
sed -i "s/hap2#h2\t/hap2#h2_/g" "${sample}.AMY_GT/${sample}.AMY_GT.sorted_combos.evaluation.tsv"

done

