#!/bin/bash

# Analyze digestive-enzyme copy numbers and haplotypes.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/AMY_80_56_4_64"



for pop in AFR AMR EUR EAS SAN Bajau Malay;do
for hap in `cut -f1 ${pop}.tsv2|grep -v Haplotype`;do
grep -w "^${hap}" "${wd}/AMY.GT.haplogroup.out" >> "${wd}/AMY.GT.haplogroup.out.${pop}.txt"
done
awk -F'\t' -v pop=${pop} 'BEGIN {OFS="\t"} {print $0,pop}' "${wd}/AMY.GT.haplogroup.out.${pop}.txt" >> "${wd}/AMY.GT.haplogroup.out.bypop.txt"
done


for hap in `cut -f2 AMY.GT.haplogroup.out|sort|uniq -c|awk '{if ($1 > 9) print $2}'`;do
grep -w "${hap}" AMY.GT.haplogroup.out >> "${wd}/AMY.GT.haplogroup.out.MAF005.txt"
done


for pop in AFR AMR EUR EAS SAN Bajau Malay;do
for hap in `cut -f1 ${pop}.tsv2|grep -v Haplotype`;do
grep -w "^${hap}" "${wd}/AMY.GT.haplogroup.out.MAF005.txt" >> "${wd}/AMY.GT.haplogroup.out.MAF005.${pop}.txt"
done
awk -F'\t' -v pop=${pop} 'BEGIN {OFS="\t"} {print $0,pop}' "${wd}/AMY.GT.haplogroup.out.MAF005.${pop}.txt" >> "${wd}/AMY.GT.haplogroup.out.MAF005.bypop.txt"
done
