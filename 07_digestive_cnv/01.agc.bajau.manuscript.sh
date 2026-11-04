#!/bin/bash

# Analyze digestive-enzyme copy numbers and haplotypes.
# Set input, reference and tool paths for your environment before execution.

conda activate agc_env

cd /path/to/69Malaysia_125HPRC2 || exit 1

agc create /path/to/Homo_sapiens.GRCh38.dna.chromosome.1_22_X_Y.fa \
/path/to/chm13v2.0.rmchrM.fa \
./asm/*.fa > ./69Malaysia_125HPRC2.agc

conda activate pgr-tk
echo 69Malaysia_125HPRC2.agc > 69Malaysia_125HPRC2.txt
pgr-mdb -w 48 -k 56 -r 12 -m 12 69Malaysia_125HPRC2.txt 69Malaysia_125HPRC2

