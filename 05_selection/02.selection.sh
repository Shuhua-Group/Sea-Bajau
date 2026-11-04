#!/usr/bin/env bash
# cd /path/to/selection
# Analyze population differentiation and regional selection statistics.
# Set input, reference and tool paths for your environment before execution.

cd /path/to/Fst || exit 1
# popA="/path/to/43Bajau_like.list"
popA="/path/to/43Bajau_like.rm7_1stkin.list"
# popB="/path/to/38Negrito_like.list"
popB="/path/to/38Negrito_like.rm9_1stkin.list"
popC="/path/to/37Malay_like.list"


# Calculate Fst
# 67West_51East.GRCh38.dv161.chr1-22.snp.norm.QC.vcf.gz
vcftools \
  --gzvcf ../../102Malay.GRCh38.dv161.chr1_22.norm.QC.SNP.vcf.gz \
  --weir-fst-pop "$popA" \
  --weir-fst-pop "$popB" \
  --out Bajau_vs_Negrito.fst \
  --fst-window-size 50000 \
  --fst-window-step 25000 &


vcftools \
  --gzvcf ../../102Malay.GRCh38.dv161.chr1_22.norm.QC.SNP.vcf.gz \
  --weir-fst-pop "$popA" \
  --weir-fst-pop "$popC" \
  --out Bajau_vs_Malay.fst \
  --fst-window-size 50000 \
  --fst-window-step 25000 &

vcftools \
  --gzvcf ../../102Malay.GRCh38.dv161.chr1_22.norm.QC.SNP.vcf.gz \
  --weir-fst-pop "$popB" \
  --weir-fst-pop "$popC" \
  --out Negrito_vs_Malay.fst \
  --fst-window-size 50000 \
  --fst-window-step 25000 &

# Calculate Tajima's D
# 67West_51East.GRCh38.dv161.chr1-22.snp.norm.QC.vcf.gz
vcftools \
  --gzvcf ../../102Malay.GRCh38.dv161.chr1_22.norm.QC.SNP.vcf.gz \
  --keep "$popA" \
  --TajimaD 50000 \
  --out Bajau_50kwindow_TajimaD &

vcftools \
  --gzvcf ../../102Malay.GRCh38.dv161.chr1_22.norm.QC.SNP.vcf.gz \
  --keep "$popB" \
  --TajimaD 50000 \
  --out Negrito_50kwindow_TajimaD &

vcftools \
  --gzvcf ../../102Malay.GRCh38.dv161.chr1_22.norm.QC.SNP.vcf.gz \
  --keep "$popC" \
  --TajimaD 50000 \
  --out Malay_50kwindow_TajimaD &

# Also calculate 2018CellBajau with this study's Bajau Of Fst
cd /path/to/Fst || exit 1
popD="/path/to/90_2018Cell.info.57Bajau.list"
vcftools \
  --gzvcf ../../malay_ho1926_SG10K_PanMLS_BajauSaluan.chr1_22.vcf.gz \
  --weir-fst-pop "$popA" \
  --weir-fst-pop "$popD" \
  --out Bajau_vs_IndoBajau.fst \
  --fst-window-size 50000 \
  --fst-window-step 25000 &

