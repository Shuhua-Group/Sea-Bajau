#!/bin/bash

# IBDNe Estimate effective population size

# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

cd /path/to/IBDNe || exit 1
for i in {1..22} ;do
/path/to/java -jar \
/path/to/hap-ibd.jar \
gt=36Bajau_like.GRCh38.dv161.chr${i}.norm.QC.SNP.shapeited.reID2.nomult.vcf.gz \
map=36Bajau_like.GRCh38.dv161.chr${i}.map.map \
out=36Bajau_like.GRCh38.dv161.chr${i}.hapIBD
done

cat 36Bajau_like.GRCh38.dv161.chr*.map.map > 36Bajau_like.GRCh38.dv161.1_22.map.map
zcat 36Bajau_like.GRCh38.dv161.chr*.hapIBD.ibd.gz > 36Bajau_like.GRCh38.dv161.1_22.hapIBD.ibd

cat 36Bajau_like.GRCh38.dv161.1_22.hapIBD.ibd| \
/path/to/java -jar \
/path/to/ibdne.23Apr20.ae9.jar \
map=36Bajau_like.GRCh38.dv161.1_22.map.map \
out=36Bajau_like.GRCh38.dv161.1_22.IBDNe \
gmax=6000

# Results stored in /path/to/IBDNe
# See results in *.ne
# Adjust gmax Inspect results at different distance thresholds

# Then, Use admixture low-quality entries in 9 exclude these samples, Run again
for i in {1..22};do
bcftools view -S /path/to/43Bajau_like.rm7_1stkin.rm_unpure9.list \
36Bajau_like.GRCh38.dv161.chr${i}.norm.QC.SNP.shapeited.reID2.nomult.vcf.gz \
-Oz -o "27Bajau.GRCh38.dv161.chr${i}.norm.QC.SNP.shapeited.reID2.nomult.vcf.gz"
plink --vcf "27Bajau.GRCh38.dv161.chr${i}.norm.QC.SNP.shapeited.reID2.nomult.vcf.gz" \
  --cm-map /path/to/chr7.genetic_map_hg38.txt2 "chr${i}" \
  --double-id \
  --out "27Bajau.GRCh38.dv161.chr${i}.map" \
  --recode \
  --set-missing-var-ids @:#
sed -i "s/^/chr/g" "27Bajau.GRCh38.dv161.chr${i}.map.map"
done

for i in {7..22} ;do
/path/to/java -jar \
/path/to/hap-ibd.jar \
gt=27Bajau.GRCh38.dv161.chr${i}.norm.QC.SNP.shapeited.reID2.nomult.vcf.gz \
map=27Bajau.GRCh38.dv161.chr${i}.map.map \
out=27Bajau.GRCh38.dv161.chr${i}.hapIBD &
done

cat 27Bajau.GRCh38.dv161.chr*.map.map > 27Bajau.GRCh38.dv161.1_22.map.map
zcat 27Bajau.GRCh38.dv161.chr*.hapIBD.ibd.gz > 27Bajau.GRCh38.dv161.1_22.hapIBD.ibd

cat 27Bajau.GRCh38.dv161.1_22.hapIBD.ibd| \
/path/to/java -jar \
/path/to/ibdne.23Apr20.ae9.jar \
map=27Bajau.GRCh38.dv161.1_22.map.map \
out=27Bajau.GRCh38.dv161.1_22.gmax2400.IBDNe \
gmax=2400
