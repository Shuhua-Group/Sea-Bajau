#!/bin/bash

# Find iHS, Fst, XPEHH intersection
# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

cd /path/to/intersect || exit 1
# 1. iHS And XPEHH intersection
# Malay_XPEHH And iHS No overlap
for i in {1..22};
do grep ^"chr${i}_" XPEHH/36Bajau_vs_37Malay.xp_ehh.fdr0.05.txt \
|cut -f1 |sort|uniq > "./xpehh.chr${i}.pos.txt" && \
cut -f1 "iHS/36Bajau_like.GRCh38.dv161.chr${i}.ihs.abs2.out" |sort|uniq \
> "./iHS.chr${i}.pos.txt" && \
cat "./xpehh.chr${i}.pos.txt" "./iHS.chr${i}.pos.txt" \
> "./xpehh_iHS.chr${i}.txt" && \
xpehh=$(wc -l < "./xpehh.chr${i}.pos.txt") && \
ihs=$(wc -l < "./iHS.chr${i}.pos.txt") && \
sum_count=$(wc -l < "./xpehh_iHS.chr${i}.txt") && \
echo -e "chr${i}\t${xpehh}\t${ihs}\t${sum_count}" >> ./xpehh_ihs_sum.count;
done

# IndoBajau_XPEHH And iHS No overlap
for i in {1..22};
do grep ^"chr${i}_" XPEHH/36Bajau_vs_IndoBajau.xp_ehh.fdr0.05.txt \
|cut -f1 |sort|uniq > "./xpehh.chr${i}.pos.txt" && \
cut -f1 "iHS/36Bajau_like.GRCh38.dv161.chr${i}.ihs.abs2.out" |sort|uniq \
> "./iHS.chr${i}.pos.txt" && \
cat "./xpehh.chr${i}.pos.txt" "./iHS.chr${i}.pos.txt" \
> "./xpehh_iHS.chr${i}.txt" && \
xpehh=$(wc -l < "./xpehh.chr${i}.pos.txt") && \
ihs=$(wc -l < "./iHS.chr${i}.pos.txt") && \
sum_count=$(wc -l < "./xpehh_iHS.chr${i}.txt") && \
echo -e "chr${i}\t${xpehh}\t${ihs}\t${sum_count}" >> ./xpehh_ihs_sum.count;
done

# 2. iHS And Fst intersection
for i in {1..22};
do cut -f1,2 Fst/Bajau_vs_IndoBajau.perSNP.fst.top0.01_FST.txt|sed "s/\t/_/g"|grep ^"chr${i}_" \
|sort|uniq > "./Fst.chr${i}.pos.txt" && \
cut -f1 "iHS/36Bajau_like.GRCh38.dv161.chr${i}.ihs.abs2.out" |sort|uniq \
> "./iHS.chr${i}.pos.txt" && \
cat "./Fst.chr${i}.pos.txt" "./iHS.chr${i}.pos.txt" \
> "./Fst_iHS.chr${i}.txt" && \
Fst=$(wc -l < "./Fst.chr${i}.pos.txt") && \
ihs=$(wc -l < "./iHS.chr${i}.pos.txt") && \
sum_count=$(wc -l < "./Fst_iHS.chr${i}.txt") && \
echo -e "chr${i}\t${Fst}\t${ihs}\t${sum_count}" >> ./Fst_ihs_sum.IndoBajau.count;
done

for i in {1..22};
do cut -f1,2 Fst/Bajau_vs_Malay.perSNP.fst.top0.01_FST.txt|sed "s/\t/_/g"|grep ^"chr${i}_" \
|sort|uniq > "./Fst.chr${i}.pos.txt" && \
cut -f1 "iHS/36Bajau_like.GRCh38.dv161.chr${i}.ihs.abs2.out" |sort|uniq \
> "./iHS.chr${i}.pos.txt" && \
cat "./Fst.chr${i}.pos.txt" "./iHS.chr${i}.pos.txt" \
> "./Fst_iHS.chr${i}.txt" && \
Fst=$(wc -l < "./Fst.chr${i}.pos.txt") && \
ihs=$(wc -l < "./iHS.chr${i}.pos.txt") && \
sum_count=$(wc -l < "./Fst_iHS.chr${i}.txt") && \
echo -e "chr${i}\t${Fst}\t${ihs}\t${sum_count}" >> ./Fst_ihs_sum.Malay.count;
done
