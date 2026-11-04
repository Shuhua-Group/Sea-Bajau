#!/usr/bin/env bash


# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

/path/to/PopLDdecay
	# Usage: PopLDdecay -InVCF <in.vcf.gz> -OutStat <out.stat>

	# -InVCF <str> Input SNP VCF Format
	# -InGenotype <str> Input SNP Genotype Format
	# -OutStat <str> OutPut Stat Dist ~ r^2/D' File

	# -SubPop <str> SubGroup Sample File List[ALLsample]
	# -MaxDist <int> Max Distance (kb) between two SNP [300]
	# -MAF <float> Min minor allele frequency filter [0.005]
	# -Het <float> Max ratio of het allele filter [0.88]
	# -Miss <float> Max ratio of miss allele filter [0.25]
	# -EHH <str> To Run EHH Region decay set StartSite [NA]
	# -OutFilterSNP OutPut the final SNP to calculate
	# -OutType <int> 1: R^2 result 2: R^2 & D' result 3:PairWise LD Out[1]
	# See the Help for more OutType [1-8] details

for list in `ls pop/*.list|grep -v 2018Cell|grep -v ho1926`;do
ind=$(basename "${list}" .list)
/path/to/PopLDdecay \
-InVCF 4panel.chrA.shapeit5.selected221.vcf.gz \
-OutStat "4panel.chrA.shapeit5.selected221.$ind" \
-SubPop "$list"
echo -e "4panel.chrA.shapeit5.selected221.$ind.stat.gz\t${ind}" >> 4panel.14pops.stat.gz.list
# echo -e "4panel.chrA.shapeit5.selected221.$ind\t${ind}" >> 4panel.14pops.list
done


perl /path/to/Plot_MultiPop.pl \
-inList  4panel.14pops.stat.gz.list  \
-output  4panel.14pops.LDdecay.stat.gz

