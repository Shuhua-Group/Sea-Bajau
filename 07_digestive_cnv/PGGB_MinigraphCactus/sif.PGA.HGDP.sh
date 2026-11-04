#!/bin/bash




# Analyze digestive-enzyme copy numbers and haplotypes.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/HGDP"
mkdir -p "$wd" && cd "$wd"
conda activate odgi_v083


# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r i <&9 || [[ -n "$i" ]]; do
    [[ -z "$i" ]] && continue
sample=$(basename $i .dedup.bam)
sample_bam=$i
ref="/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa"
graph="/path/to/PGA_region_seq.fa.gz.bf3285f.37c1a4f.8635317.smooth.final.gfa"

mkdir -p "$wd/$sample" && cd "$wd/$sample"
mkdir -p resources/cram
mkdir -p resources/ref
ln -s "$sample_bam" "resources/cram/$sample.bam"
sample_bai=`echo $sample_bam | sed 's/\.bam$/.bai/'`
ln -s "$sample_bai" "resources/cram/$sample.bam.bai"
ln -s "$ref" resources/ref/
ln -s "${ref}.fai" resources/ref/

ln -s /path/to/workflow workflow
mkdir -p config
echo -e "samples: config/samples.tsv" > config/config.yaml
echo -e "reference: resources/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa" >> config/config.yaml
echo -e "graph: $graph" >> config/config.yaml
echo -e "region: \"chr11:61162929-61331529\"" >> config/config.yaml


echo -e "sample_id\tcram" > config/samples.tsv
echo -e "$sample\tresources/cram/$sample.bam" >> config/samples.tsv

# conda activate odgi_v083
snakemake -s /path/to/Snakefile \
  --configfile config/config.yaml \
  --use-singularity --singularity-args "-B /path/to/share9" -j 5 cosigt
done 9< <(cat -- "../HGDP.ngs.list" | grep bam$|head -264|tail -24) &
