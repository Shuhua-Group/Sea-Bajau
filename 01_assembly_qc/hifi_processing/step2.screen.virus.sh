#!/bin/bash

# Process HiFi reads and summarize assembly quality.
# Set input, reference and tool paths for your environment before execution.

gxdb="/dev/shm/gxdb"
mkdir "$gxdb"

python3 /path/to/fcs.py --image /path/to/fcs-gx.sif \
db get --mft "/path/to/all.manifest" --dir "$gxdb"

gxdb="/dev/shm/gxdb"
# Read one record per line; keep backslashes and a final line without a newline.
while IFS= read -r sample <&9 || [[ -n "$sample" ]]; do
    [[ -z "$sample" ]] && continue
hap1_fa="/path/to/${sample}.hap1.rDNAremoved.fasta"
hap2_fa="/path/to/${sample}.hap2.rDNAremoved.fasta"
outdir="/path/to/${sample}"
mkdir -p "${outdir}"
# screen tax-id 9606 homo sapien
script="#!/bin/bash
#SBATCH --job-name=${sample}.virus_rm
#SBATCH --partition=batch
#SBATCH --output= ${outdir}/rm_virus.log
#SBATCH --error= ${outdir}/rm_virus.err
#SBATCH --mincpus=32
#SBATCH --mem=200g

python3 /path/to/fcs.py --image /path/to/fcs-gx.sif \
screen genome --fasta ${hap1_fa} \
--out-dir ${outdir} \
--gx-db "$gxdb" --tax-id 9606
python3 /path/to/fcs.py --image /path/to/fcs-gx.sif \
screen genome --fasta ${hap2_fa} \
--out-dir ${outdir} \
--gx-db "$gxdb" --tax-id 9606
# clean the genome
cd ${outdir}
cat ${hap1_fa} \
| python3 /path/to/fcs.py --image /path/to/fcs-gx.sif \
clean genome --action-report ${outdir}/${sample}.1.9606.fcs_gx_report.txt \
--output ${sample}.1.9606.clean.fasta --contam-fasta-out ${sample}.1.9606.contam.fasta
cat ${hap2_fa} \
| python3 /path/to/fcs.py --image /path/to/fcs-gx.sif \
clean genome --action-report ${outdir}/${sample}.2.9606.fcs_gx_report.txt \
--output ${sample}.2.9606.clean.fasta --contam-fasta-out ${sample}.2.9606.contam.fasta
"
echo -e "$script" > "${outdir}/rm_virus.sh"
echo -e "sbatch ${outdir}/rm_virus.sh" >> /path/to/run.sbatch_virus.sh
done 9< "/path/to/sample.list"

