#!/bin/bash

# Analyze digestive-enzyme copy numbers and haplotypes.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/mc_retry"

cd "$wd" || exit 1

MYJOBSTORE=$wd/js_cactus
MYBUCKET=$wd/MC_graph/
mkdir -p "$MYBUCKET"

PREFIX=PRSS_64_56_6_28.MC

seqfile=$wd/seqfile.txt

######## seqfile Prepare again #########
touch temp.list
for contig in `cut -f1 ../PRSS_region_seq.fa.fai|grep chm13`;do
# NA19338_hap2#NA19338#2#CM087780.1_60597798_60769095_0
# Extract sample names and hap Name
sample=`echo $contig | cut -d "#" -f1`
sample2=`echo $contig | cut -d "#" -f1|cut -d "p" -f1`
hap=`echo $contig | cut -d "#" -f1| cut -d "p" -f2`
echo "$sample2" "$hap" >> temp.list
# echo $contig
var1=$(echo "$contig" | cut -d'#' -f2,3,4 | cut -d'_' -f1)
var2=$(echo "$contig" | awk -F'_' '{print $(NF-2)}')
var3=$(echo "$contig" | awk -F'_' '{print $(NF-1)}')
# start Extend leftward 40000bp, If smaller than 0 Use 1
if [ "$var2" -lt 40000 ]; then
start=1
else
start=$((var2 - 40000))

fi
echo "$start"
samtools faidx "../../asm/${sample}.fa" $var1":"$start"-"$var3 > "asm/${sample2}p${hap}.fa"
done

cd asm || exit 1

for i in `ls *.fa|grep chm13`;do
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


echo -e "GRCh38PRSS\t${wd}/asm/GRCh38.PRSS.rename.fa" > seqfile.txt
echo -e "CHM13PRSS\t${wd}/asm/chm13v2.PRSS.rename.fa" >> seqfile.txt

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
hap=$(ls ${wd}/asm/${sample}*hap*.rename.fa)
echo -e "${sample}\t${hap}" >> seqfile.txt
done 9< "haploid.list"



logfile=$wd/cactus_v293.log

source /path/to/activate
mkdir -p "$wd/temp"
export TMPDIR=$wd/temp

cactus-pangenome "${MYJOBSTORE}" "${seqfile}" --outName "${PREFIX}" --outDir "${MYBUCKET}" \
--reference GRCh38PRSS CHM13PRSS --filter 38 \
--lastTrain --snarlStats --haplo --maxMemory 200G \
--giraffe full clip filter \
--vcf clip filter full \
--vcfReference GRCh38PRSS CHM13PRSS \
--gbz clip filter full \
--xg clip filter full \
--chrom-vg clip filter full \
--viz --odgi --chrom-vg clip filter \
--chrom-og --gbz clip filter full \
--gfa clip full filter \
--vcfwave \
--logFile "$logfile" \
--mgCores 8 --mapCores 8 --consCores 8 --indexCores 8
# --restart


# Then run python Script

mkdir -p "$wd/LDblockshow" && cd "$wd/LDblockshow"

bcftools query -l ../PRSS_64_56_6_28.MC.full.wave.norm.haploid.vcf.gz |grep -v _MINIGRAPH_ > ./newlist

bcftools view -S ./newlist -Oz -o test.vcf.gz \
../PRSS_64_56_6_28.MC.full.wave.norm.haploid.vcf.gz

-Region GRCh38.PRSS:1-59000 \

/path/to/LDBlockShow \
-InVCF test.vcf.gz \
-OutPut PRSS_64_56_6_28.flanking.LDblock \
-Region GRCh38.PRSS:1-89723 \
-SeleVar 2 \
-BlockCut 0.9 \
-MAF 0.001

/path/to/ShowLDSVG \
-InPreFix PRSS_64_56_6_28.flanking.LDblock \
-OutPut PRSS_64_56_6_28.flanking.LDblock.plot \
-crBegin 77,9,95 \
-crMiddle 32,153,137 \
-crEnd 252,230,18 \
-OutPdf
