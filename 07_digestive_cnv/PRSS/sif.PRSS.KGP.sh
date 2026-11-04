#!/bin/bash




# Analyze digestive-enzyme copy numbers and haplotypes.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/KGP"
mkdir -p "$wd" && cd "$wd"
conda activate odgi_v083


for i in `ls /path/to/*.dedup.bam /path/to/*.dedup.bam|head -1250`;do
sample=$(basename $i .dedup.bam)
sample_bam=$i
ref="/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa"
graph="/path/to/PRSS_region_seq.fa.gz.bf3285f.e219c8e.3d8788a.smooth.final.gfa"
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
echo -e "region: \"chr7:142737469-142787652\"" >> config/config.yaml


echo -e "sample_id\tcram" > config/samples.tsv
echo -e "$sample\tresources/cram/$sample.bam" >> config/samples.tsv

# conda activate odgi_v083
snakemake -s /path/to/Snakefile \
  --configfile config/config.yaml \
  --use-singularity --singularity-args "-B /path/to/bajau_project" -j 5 cosigt
done &
