#!/bin/bash

# Use germline2 Identify IBD

# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

wd=/path/to/GERMLINE2
cd "$wd" || exit 1

for list in `ls pop/*.list|grep -v 2018Cell|grep -v ho1926`;do
ind=$(basename "${list}" .list)
mkdir -p input
for i in {1..22};do
bcftools view -S ${list} 4panel.chr${i}.shapeit5.vcf.gz \
-Oz -o input/4panel.${ind}.chr${i}.vcf.gz
/path/to/RelateFileFormats \
--mode ConvertFromVcf \
--haps input/4panel.${ind}.chr${i}.haps \
--sample input/4panel.${ind}.chr${i}.sample \
-i input/4panel.${ind}.chr${i}
# Run g2, This step is fast
mkdir -p 4panel_g2out
## -h is strict mode, hap Do not mix between trans
/path/to/g2 \
-h \
-m 2.0 \
input/4panel.${ind}.chr${i}.haps \
input/4panel.${ind}.chr${i}.sample \
/path/to/chr${i}.genetic_map_hg38.txt2 \
4panel_g2out/4panel.${ind}.chr${i}.m2

## Run an additional comparison without -h Of
/path/to/g2 \
input/4panel.${ind}.chr${i}.haps \
input/4panel.${ind}.chr${i}.sample \
/path/to/chr${i}.genetic_map_hg38.txt2 \
4panel_g2out/4panel.${ind}.chr${i}.nohap
done &
done

# g2 results compared with HapIBD differ
# First convert to ibd Format
# Then run ibdne
# Compare with hapibd Check agreement between results

for list in `ls pop/*.list|grep -v 2018Cell|grep -v ho1926`;do
ind=$(basename "${list}" .list)
rm -f "4panel_g2out/4panel.${ind}.chrA.m2.ibd"
for i in {1..22};do
python germline2ibd.py "4panel_g2out/4panel.${ind}.chr${i}.m2"
ibd="4panel_g2out/4panel.${ind}.chr${i}.m2.ibd.gz"
zcat "${ibd}" >> "4panel_g2out/4panel.${ind}.chrA.m2.ibd"
done
mkdir -p IBDNe_output
cat "4panel_g2out/4panel.${ind}.chrA.m2.ibd"| \
/path/to/java -jar \
/path/to/ibdne.23Apr20.ae9.jar \
map=/path/to/plink.chrA.wchr.GRCh38.map \
out=IBDNe_output/${ind}.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge \
gmax=300
done

# Also calculate cross-population pairs
for list in `ls pop/*.list|grep -v ho1926|grep -v 2018Cell|grep -v 36Bajau`;do
ind=$(basename "${list}" .list)
rm -f "4panel_g2out/4panel.${ind}_36Bajau.combine.chrA.m2.ibd"
for i in {1..22};do
cat pop/36Bajau.list $list > temp.${ind}_36Bajau.combine.list
bcftools view -S temp.${ind}_36Bajau.combine.list \
4panel.chr${i}.shapeit5.vcf.gz \
-Oz -o ${ind}_36Bajau.combine.chr${i}.shapeit5.vcf.gz
Convert first Haps/sample Format
/path/to/RelateFileFormats \
--mode ConvertFromVcf \
--haps input/4panel.${ind}_36Bajau.combine.chr${i}.haps \
--sample input/4panel.${ind}_36Bajau.combine.chr${i}.sample \
-i ${ind}_36Bajau.combine.chr${i}.shapeit5
# g2
/path/to/g2 \
-h \
input/4panel.${ind}_36Bajau.combine.chr${i}.haps \
input/4panel.${ind}_36Bajau.combine.chr${i}.sample \
/path/to/chr${i}.genetic_map_hg38.txt2 \
4panel_g2out/4panel.${ind}_36Bajau.combine.chr${i}.m2
python germline2ibd.py "4panel_g2out/4panel.${ind}_36Bajau.combine.chr${i}.m2"
zcat "4panel_g2out/4panel.${ind}_36Bajau.combine.chr${i}.m2.ibd.gz" >> "4panel_g2out/4panel.${ind}_36Bajau.combine.chrA.m2.ibd"
done
done


# Then select between-population pairs
for list in `ls pop/*.list|grep -v ho1926|grep -v 2018Cell|grep -v 36Bajau`;do
    ind=$(basename "$list" .list)
    popA="pop/36Bajau.list"
    popB="$list"
    ibd_file="4panel_g2out/4panel.${ind}_36Bajau.combine.chrA.m2.ibd"
    cross_file="4panel_g2out/4panel.${ind}_36Bajau.combine.chrA.m2.ibd.cross"
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


# Then run ibdne
for list in `ls pop/*.list|grep -v ho1926|grep -v 2018Cell|grep -v 36Bajau`;do
ind=$(basename "${list}" .list)
mkdir -p IBDNe_output
cat "4panel_g2out/4panel.${ind}_36Bajau.combine.chrA.m2.ibd.cross"| \
/path/to/java -jar \
/path/to/ibdne.23Apr20.ae9.jar \
map=/path/to/plink.chrA.wchr.GRCh38.map \
out=IBDNe_output/${ind}_36Bajau.combine.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge \
gmax=300
done

# Plot the results for inspection
python IBDNe_plot.multi.py \
-i IBDNe_output/25Malay.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne \
IBDNe_output/36Bajau.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne \
IBDNe_output/EGDP.Borneo.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne \
IBDNe_output/EGDP.Luzon.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne \
IBDNe_output/EGDP.Myanmar.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne \
IBDNe_output/EGDP.Vietnam.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne \
IBDNe_output/Pacific.Barito.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne \
IBDNe_output/Pacific.Besemah.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne \
IBDNe_output/Pacific.Kalimantan.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne \
IBDNe_output/Pacific.Mamanwa.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne \
IBDNe_output/Pacific.Manobo.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne \
IBDNe_output/Pacific.Semende.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne \
IBDNe_output/PanMLS.Jakun.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne \
IBDNe_output/PanMLS.Temiar.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne \
-o g2_m2_IBDNe.4panel.300ge.cutoff5ge.pdf \
--min_gen 5

# Run both
python IBDNe_plot.multi.py \
-i IBDNe_output/36Bajau.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne \
IBDNe_output/Pacific.Barito.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne \
-o Pacific.Barito_36Bajau.g2_m2_IBDNe.4panel.300ge.cutoff5ge.pdf \
--min_gen 5

# Then calculate rCCR
for list in `ls pop/*.list|grep -v ho1926|grep -v 2018Cell|grep -v 36Bajau`;do
    mkdir -p rCCR
    ind=$(basename "$list" .list)
    # Input file
    file_pop="IBDNe_output/${ind}.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne"
    file_bajau="IBDNe_output/36Bajau.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne"
    file_cross="IBDNe_output/${ind}_36Bajau.combine.chrA.4panel.shapeit5.g2_m2.IBDNe.300ge.ne"
    # Output file
    out_file="rCCR/${ind}_36Bajau.4panel.g2_m2_IBDNe.RCC"
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
