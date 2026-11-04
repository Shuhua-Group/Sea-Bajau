#!/bin/bash




# Analyze digestive-enzyme copy numbers and haplotypes.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/HGDP"
mkdir -p "$wd" && cd "$wd"
conda activate odgi_v083

# for i in {1..264};do
# bam=`cat ../HGDP.ngs.list|grep bam$|head -n $i|tail -n 1`
# sample=`basename $bam .dedup.bam`
# if [ ! -f $sample/results/cosigt_results/$sample/best_genotype.tsv ];then
# echo $i
# fi
# done


# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r i <&9 || [[ -n "$i" ]]; do
    [[ -z "$i" ]] && continue
sample=$(basename $i .dedup.bam)
sample_bam=$i
ref="/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa"
graph="/path/to/AMY_region_seq.fa.gz.bf3285f.e219c8e.78f50d0.smooth.final.gfa"
mkdir -p "$wd/$sample" && cd "$wd/$sample"
if [ ! -f "$sample/results/cosigt_results/$sample/best_genotype.tsv" ];then
mkdir -p resources/cram
mkdir -p resources/ref
ln -sf "$sample_bam" "resources/cram/$sample.bam"
sample_bai=`echo $sample_bam | sed 's/\.bam$/.bai/'`
ln -sf "$sample_bai" "resources/cram/$sample.bam.bai"
ln -sf "$ref" resources/ref/
ln -sf "${ref}.fai" resources/ref/
ln -sf /path/to/workflow workflow
mkdir -p config
echo -e "samples: config/samples.tsv" > config/config.yaml
echo -e "reference: resources/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa" >> config/config.yaml
echo -e "graph: $graph" >> config/config.yaml
echo -e "region: \"chr1:103405896-103913228\"" >> config/config.yaml

echo -e "sample_id\tcram" > config/samples.tsv
echo -e "$sample\tresources/cram/$sample.bam" >> config/samples.tsv

snakemake -s /path/to/Snakefile \
  --configfile config/config.yaml \
  --use-singularity --singularity-args "-B /path/to/share9" -j 5 cosigt &
fi
done 9< <(cat -- "../HGDP.ngs.list" | grep bam$)
