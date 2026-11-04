#!/bin/bash

# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

wd="/path/to/76BajauManu"
mkdir -p "$wd"
cd "$wd" || exit 1

popA="../pop_info/76BajauManu/36Bajau.list"

mkdir -p output
for i in `ls ../pop_info/76BajauManu/*list | grep 25Malay`; do
    (pop=$(basename "$i" .list) && \
    vcftools --gzvcf ../re_merge/7panel.chrA.vcf.gz \
        --weir-fst-pop "$i" \
        --weir-fst-pop "$popA" \
        --out "output/36Bajau_vs_${pop}.perSNP.fst" && \
    awk -v p=$pop '{if($3 != "-nan" && NR > 1 && $3 >= 0) sum += $3; if($3 != "-nan" && NR > 1) count++} END {print p, sum/count}' \
    "output/36Bajau_vs_${pop}.perSNP.fst.weir.fst" >> output/Fst_mean.txt &)
done

for i in `ls ../pop_info/76BajauManu/*list | grep -v 36Bajau.list`; do
    pop=$(basename "$i" .list)

    awk -F'\t' -v p=$pop '
    NR == 1 { next }  #  Skip  header
    $3 == "-nan" { next }  #  Skip  -nan
    {
        if ($3 >= 0) {
            sum += $3
            sumsq += ($3)^2
        }
        count++
    }
    END {
        if (count > 1) {
            mean = sum / count
            variance = (sumsq / count) - (mean)^2
            if (variance < 0) variance = 0
            stddev = sqrt(variance)
            print p, mean, stddev
        } else {
            print p, "NA", "NA"
        }
    }' "output/36Bajau_vs_${pop}.perSNP.fst.weir.fst" >> output/Fst_mean_std.txt

done



sed -i "s/ /\t/g" output/Fst_mean.txt
sed -i "s/ /\t/g" output/Fst_mean_std.txt
