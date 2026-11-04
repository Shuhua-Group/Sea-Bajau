#!/bin/bash

# IBDNe Estimate effective population size

# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

cd /path/to/33Malay_36Bajau_20Negrito || exit 1
for ind in 9unpure.Bajau 27pure.Bajau 20Negrito.ngs 33Malay;do
for i in {1..22};do
list="${ind}.list"
bcftools view -S "${list}" \
69Malay_20Negrito.chr${i}.QC.shapeited.bcf \
-Oz -o "${ind}.chr${i}.QC.shapeited.vcf.gz"
plink --vcf "${ind}.chr${i}.QC.shapeited.vcf.gz" \
  --cm-map "/path/to/chr${i}.genetic_map_hg38.txt2" "chr${i}" \
  --double-id \
  --out "${ind}.chr${i}.QC.shapeited.map" \
  --recode \
  --set-missing-var-ids @:#
sed -i "s/^/chr/g" "${ind}.chr${i}.QC.shapeited.map.map"
/path/to/java -jar \
/path/to/hap-ibd.jar \
gt=${ind}.chr${i}.QC.shapeited.vcf.gz \
map=${ind}.chr${i}.QC.shapeited.map.map \
out=${ind}.chr${i}.QC.shapeited.map.hapIBD
done
cat ${ind}.chr*.QC.shapeited.map.map > "${ind}.chrA.QC.shapeited.map.map"
zcat ${ind}.chr*.QC.shapeited.map.hapIBD.ibd.gz > "${ind}.chrA.QC.shapeited.map.hapIBD.ibd"
done



for ind in 9unpure.Bajau 27pure.Bajau 20Negrito.ngs 33Malay;do
cat "${ind}.chrA.QC.shapeited.map.hapIBD.ibd"| \
/path/to/java -jar \
/path/to/ibdne.23Apr20.ae9.jar \
map=${ind}.chrA.QC.shapeited.map.map \
out=${ind}.chrA.QC.shapeited.map.hapIBD.IBDNe.6000ge \
gmax=6000
done
