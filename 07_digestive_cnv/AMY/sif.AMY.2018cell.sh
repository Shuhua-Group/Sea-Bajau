#!/bin/bash




# Analyze digestive-enzyme copy numbers and haplotypes.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/2018cell"
mkdir -p "$wd" && cd "$wd"


# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r i <&9 || [[ -n "$i" ]]; do
    [[ -z "$i" ]] && continue
sample=$i
sample_bam=$(ls /path/to/${sample}.sort.rmdup.realign.md.bam)
ref="/path/to/human_g1k_v37.1_22XYMT.fa"
graph="/path/to/AMY_region_seq.fa.gz.bf3285f.e219c8e.78f50d0.smooth.final.gfa"
mkdir -p "$wd/$sample" && cd "$wd/$sample"
mkdir -p resources/cram
mkdir -p resources/ref
ln -s "$sample_bam" "resources/cram/$sample.bam"
ln -s "$sample_bam.bai" "resources/cram/$sample.bam.bai"
ln -s "$ref" resources/ref/
ln -s "${ref}.fai" resources/ref/

ln -s /path/to/workflow workflow
mkdir -p config
echo -e "samples: config/samples.tsv" > config/config.yaml
echo -e "reference: resources/ref/human_g1k_v37.1_22XYMT.fa" >> config/config.yaml
echo -e "graph: $graph" >> config/config.yaml
echo -e "region: \"1:103948518-104455850\"" >> config/config.yaml

echo -e "sample_id\tcram" > config/samples.tsv
echo -e "$sample\tresources/cram/$sample.bam" >> config/samples.tsv

snakemake -s /path/to/Snakefile \
  --configfile config/config.yaml \
  --use-singularity --singularity-args "-B /path/to/bajau_project" -j 5 cosigt

done 9< <(cat -- "$wd/../2018cell.list" | head -1)
