#!/bin/bash

# IBDNe Estimate effective population size

# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

cd /path/to/IBDNe_76BajauManu_remerge || exit 1
# for ind in 7unpure.Bajau 29pure.Bajau;do
for list in `ls pop/*.list|grep -v ho1926|grep -v 2018Cell|grep -v 36Bajau|tail -3`;do
ind=$(basename "${list}" .list)
rm "${ind}_36Bajau.combine.chrA.shapeit5.map.hapIBD.ibd"
for i in {1..22};do
cat pop/36Bajau.list "$list" > "temp.${ind}_36Bajau.combine.list"
bcftools view -S "temp.${ind}_36Bajau.combine.list" \
4panel.chr${i}.shapeit5.vcf.gz \
-Oz -o "${ind}_36Bajau.combine.chr${i}.shapeit5.vcf.gz"
# awk 'NF>0{gsub("\r",""); print $1"\t"$1}' ${list} > temp.${ind}.keep
plink --vcf ${ind}_36Bajau.combine.chr${i}.shapeit5.vcf.gz \
--cm-map /path/to/chr${i}.genetic_map_hg38.txt2 chr${i} \
--double-id \
--out ${ind}_36Bajau.combine.chr${i}.shapeit5.map \
--recode \
--set-missing-var-ids @:#
sed -i "s/^/chr/g" ${ind}_36Bajau.combine.chr${i}.shapeit5.map.map
/path/to/java -jar \
/path/to/hap-ibd.jar \
gt=${ind}_36Bajau.combine.chr${i}.shapeit5.vcf.gz \
map=/path/to/plink.chr${i}.wchr.GRCh38.map \
out=${ind}_36Bajau.combine.chr${i}.shapeit5.map.hapIBD
cat ${ind}_36Bajau.combine.chr${i}.shapeit5.map.map >> ${ind}_36Bajau.combine.chrA.shapeit5.map.map
zcat "${ind}_36Bajau.combine.chr${i}.shapeit5.map.hapIBD.ibd.gz" >> "${ind}_36Bajau.combine.chrA.shapeit5.map.hapIBD.ibd"
done &
done

# Then select ibd Segments shared across populations

# Iterate over pop All files in the directory: list File, Exclude unused
for list in pop/*.list; do
    [[ "$list" =~ ho1926|2018Cell|36Bajau ]] && continue  #  Exclude specified files
    ind=$(basename "$list" .list)
    popA="pop/36Bajau.list"
    popB="$list"
    ibd_file="${ind}_36Bajau.combine.chrA.shapeit5.map.hapIBD.ibd"
    cross_file="${ind}_36Bajau.combine.chrA.shapeit5.map.hapIBD.ibd.cross"
    # Truncate the output file
    > "$cross_file"

    # Read popA into an array a
    awk -v popA="$popA" -v popB="$popB" '
    BEGIN{
        while ((getline lineA < popA) > 0) a[lineA]=1
        while ((getline lineB < popB) > 0) b[lineB]=1
    }
    {
        if (($1 in a && $3 in b) || ($1 in b && $3 in a)) print $0
    }' "$ibd_file" > "$cross_file"

done


#
for list in `ls pop/*.list|grep -v ho1926|grep -v 2018Cell|grep -v 36Bajau|head -10|tail -4`;do
ind=$(basename "${list}" .list)
cat "${ind}_36Bajau.combine.chrA.shapeit5.map.hapIBD.ibd.cross"| \
/path/to/java -jar \
/path/to/ibdne.23Apr20.ae9.jar \
map=/path/to/plink.chrA.wchr.GRCh38.map \
out=${ind}_36Bajau.combine.chrA.shapeit5.map.hapIBD.IBDNe.300ge \
gmax=300 &
done

for list in pop/*.list; do
    # Exclude specified files
    [[ "$list" =~ ho1926|2018Cell|36Bajau ]] && continue
    ind=$(basename "$list" .list)
    # Input file
    file_pop="${ind}.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne"
    file_bajau="36Bajau.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne"
    file_cross="${ind}_36Bajau.combine.chrA.shapeit5.map.hapIBD.IBDNe.300ge.ne"
    # Output file
    out_file="${ind}_36Bajau.RCC"
    # Check that the file exists
    if [[ ! -f "$file_pop" || ! -f "$file_bajau" || ! -f "$file_cross" ]]; then
        echo "Warning: Missing file for $ind, skipping."
        continue
    fi
    echo -e "Ne_${ind}cross\tNe_36Bajau\tNe_${ind}" > "$out_file"
    # Calculate RCC
    # Start from the second row, Assume the first row is a header
    paste "$file_cross" "$file_bajau" "$file_pop" |tail -n +2| \
    awk '{
        Ne12 = $2;          #  Admixed population  Ne
        Ne1 = $6;          #  Population 1 Ne
        Ne2 = $10;         #  Population 2 Ne
        print Ne12"\t"Ne1"\t"Ne2
    }'>> "$out_file"
    echo "RCC file saved: $out_file"
done

