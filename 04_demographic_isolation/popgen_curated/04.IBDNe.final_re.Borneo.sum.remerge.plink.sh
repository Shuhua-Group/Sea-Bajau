#!/bin/bash

# IBDNe Estimate effective population size

# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

cd /path/to/IBDNe_76BajauManu_remerge || exit 1
# for ind in 7unpure.Bajau 29pure.Bajau;do
for list in `ls pop/*.list|grep -v ho1926|grep -v 2018Cell`;do
ind=$(basename "${list}" .list)
rm "${ind}.chrA.shapeit5.map.hapIBD.ibd"
for i in {1..22};do
bcftools view -S "${list}" \
4panel.chr${i}.shapeit5.vcf.gz \
-Oz -o "${ind}.chr${i}.shapeit5.vcf.gz"
awk 'NF>0{gsub("\r",""); print $1"\t"$1}' ${list} > temp.${ind}.keep
plink --vcf ${ind}.chr${i}.shapeit5.vcf.gz \
--cm-map /path/to/chr${i}.genetic_map_hg38.txt2 chr${i} \
--double-id \
--out ${ind}.chr${i}.shapeit5.map \
--recode \
--set-missing-var-ids @:#
sed -i "s/^/chr/g" ${ind}.chr${i}.shapeit5.map.map
/path/to/java -jar \
/path/to/hap-ibd.jar \
gt=${ind}.chr${i}.shapeit5.vcf.gz \
map=/path/to/plink.chr${i}.wchr.GRCh38.map \
out=${ind}.chr${i}.shapeit5.map.hapIBD
cat ${ind}.chr${i}.shapeit5.map.map >> ${ind}.chrA.shapeit5.map.map
zcat "${ind}.chr${i}.shapeit5.map.hapIBD.ibd.gz" >> "${ind}.chrA.shapeit5.map.hapIBD.ibd"
done &
done



for list in `ls pop/*.list|grep -v ho1926|grep -v 2018Cell|tail -2`;do
ind=$(basename "${list}" .list)
cat "${ind}.chrA.shapeit5.map.hapIBD.ibd"| \
/path/to/java -jar \
/path/to/ibdne.23Apr20.ae9.jar \
map=/path/to/plink.chrA.wchr.GRCh38.map \
out=${ind}.chrA.shapeit5.map.hapIBD.IBDNe.300ge \
gmax=300 &
done

# conda activate py39
for list in `ls pop/*.list|grep -v ho1926|grep -v 2018Cell`;do
ind=$(basename "${list}" .list)
for cutoff in 5 10 20;do
python IBDNe_plot.py \
-i "${ind}.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne" \
-o "${ind}.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne.cut${cutoff}ge.pdf" \
-c "${cutoff}"
done
done

for cutoff in 5 10 20;do
python IBDNe_plot.multi.v2.py \
-i 25Malay.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
36Bajau.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
EGDP.Borneo.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
EGDP.Luzon.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
EGDP.Myanmar.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
EGDP.Vietnam.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
Pacific.Barito.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
Pacific.Besemah.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
Pacific.Kalimantan.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
Pacific.Mamanwa.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
Pacific.Manobo.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
Pacific.Semende.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
PanMLS.Jakun.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
PanMLS.Temiar.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
-o "4panel.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne.cut${cutoff}ge.max100ge.pdf" \
--min_gen "${cutoff}" \
--max_gen 100 \
--figsize 10 7
done

for cutoff in 5 10 20;do
python IBDNe_plot.multi.py \
-i 36Bajau.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
Pacific.Barito.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
-o "4panel.Pacific.Barito_36Bajau.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne.cut${cutoff}ge.pdf" \
--min_gen "${cutoff}"
done


for cutoff in 5 10 20;do
python IBDNe_plot.multi.v2.py \
-i 36Bajau.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
Pacific.Barito.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne \
-o "4panel.Pacific.Barito_36Bajau.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne.cut${cutoff}ge.max100ge.pdf" \
--min_gen "${cutoff}" \
--max_gen 100 \
--figsize 10 7
done
