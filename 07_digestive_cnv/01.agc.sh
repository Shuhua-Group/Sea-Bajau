#!/bin/bash

# 1. Step 1: Create agc File
# Analyze digestive-enzyme copy numbers and haplotypes.
# Set input, reference and tool paths for your environment before execution.

conda activate agc_env
cd /path/to/118malay_cpc1_hprc_cpc2_Viet_Jasa_chm13_cn1.ref_GRCh38 || exit 1

agc create /path/to/Homo_sapiens.GRCh38.dna.chromosome.1_22_X_Y.fa \
/path/to/chm13v2.0.rmchrM.fa \
/path/to/CN1_combine.v0.8.rmchrM.fasta \
../../00.sum/asm_fa/*fa \
/path/to/CPC*fa \
/path/to/NA*fa \
/path/to/ksa*fa \
/path/to/KSA*fa \
/path/to/VN*fa \
/path/to/*fa \
/path/to/*fa \
> ./118malay_cpc1_hprc_cpc2_Viet_Jasa_chm13_cn1.ref_GRCh38.agc

# 2. Step 2: Create an index
# Note: pgr-mdb Do not pass directly agc File, Write the filenames into a file Within, Pass this file
conda activate pgr-tk
# $ which pgr-mdb
# /path/to/pgr-mdb
# $ pwd
# /path/to/07.MAPgraph
echo 118malay_cpc1_hprc_cpc2_Viet_Jasa_chm13_cn1.ref_GRCh38.agc > 118malay_cpc1_hprc_cpc2_Viet_Jasa_chm13_cn1.ref_GRCh38.txt
mkdir w48_k56_r12_m12
cd w48_k56_r12_m12 || exit 1
pgr-mdb -w 48 -k 56 -r 12 -m 12 ../118malay_cpc1_hprc_cpc2_Viet_Jasa_chm13_cn1.ref_GRCh38.txt \
118malay_cpc1_hprc_cpc2_Viet_Jasa_chm13_cn1.ref_GRCh38
