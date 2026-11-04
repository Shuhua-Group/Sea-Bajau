"""Process variant calls and summarize variant catalogs."""

import gzip
import concurrent.futures
import pandas as pd
from collections import defaultdict


############################################################
# Input VCF
############################################################

# vcf_path = "/path/to/130APC.SVPOP_sv_del.telo.vcf.gz"
# out_path = "/path/to/130APC.SVPOP_sv_del.telo.haplotype_variant_stats.v2.txt"

# vcf_path = "/path/to/130APC.SVPOP_sv_ins.telo.vcf.gz"
# out_path = "/path/to/130APC.SVPOP_sv_ins.telo.haplotype_variant_stats.v2.txt"

# vcf_path = "/path/to/130APC.SVPOP_sv_del.no_telo.cenSat.vcf.gz"
# out_path = "/path/to/130APC.SVPOP_sv_del.no_telo.cenSat.haplotype_variant_stats.v2.txt"

# vcf_path = "/path/to/130APC.SVPOP_sv_ins.no_telo.cenSat.vcf.gz"
# out_path = "/path/to/130APC.SVPOP_sv_ins.no_telo.cenSat.haplotype_variant_stats.v2.txt"

# vcf_path = "/path/to/130APC.SVPOP_sv_del.no_telo-cenSat.vcf.gz"
# out_path = "/path/to/130APC.SVPOP_sv_del.no_telo-cenSat.haplotype_variant_stats.v2.txt"

vcf_path = "/path/to/130APC.SVPOP_sv_ins.no_telo-cenSat.vcf.gz"
out_path = "/path/to/130APC.SVPOP_sv_ins.no_telo-cenSat.haplotype_variant_stats.v2.txt"


############################################################
# Parse VCF header
############################################################

def parse_vcf_header(vcf_path):
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                return fields[9:]
    return []


sample_names = parse_vcf_header(vcf_path)


print(
    f"Detected {len(sample_names)} samples"
)


############################################################
# Build haplotype names
############################################################

# Each sample has two haplotypes:
#
# sample.hap1
# sample.hap2

haplotype_names = []

for sample in sample_names:
    haplotype_names.append(
        f"{sample}.hap1"
    )
    haplotype_names.append(
        f"{sample}.hap2"
    )


############################################################
# Initialize statistics
############################################################

allele_counts = defaultdict(int)

total_lengths = defaultdict(int)


############################################################
# Process one chunk
############################################################

def process_chunk(lines):
    local_allele_counts = defaultdict(int)
    local_total_lengths = defaultdict(int)
    for line in lines:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        ####################################################
        # REF / ALT
        ####################################################
        ref = fields[3]
        alt = fields[4]
        ####################################################
        # SV length
        ####################################################
        L = abs(
            len(ref) - len(alt)
        )
        ####################################################
        # Genotypes
        ####################################################
        genotypes = fields[9:]
        for i, gt_field in enumerate(genotypes):
            sample = sample_names[i]
            ################################################
            # extract GT
            ################################################
            gt = gt_field.split(':')[0]
            ################################################
            # Missing genotype
            ################################################
            if gt in {'.', './.', '.|.'}:
                continue
            ################################################
            # Require phased genotype
            ################################################
            if '|' not in gt:
                # Unphased genotype cannot be assigned
                # to a specific haplotype.
                #
                # For example:
                #
                # 0/1
                #
                # We know one haplotype carries ALT,
                # but cannot determine whether it is
                # hap1 or hap2.
                #
                # Therefore skip this genotype for
                # haplotype-level statistics.
                continue
            ################################################
            # split phased genotype
            ################################################
            alleles = gt.split('|')
            ################################################
            # Make sure diploid genotype
            ################################################
            if len(alleles) != 2:
                continue
            hap1 = f"{sample}.hap1"
            hap2 = f"{sample}.hap2"
            ################################################
            # haplotype 1
            ################################################
            allele1 = alleles[0]
            if allele1 not in {'0', '.'}:
                local_allele_counts[hap1] += 1
                local_total_lengths[hap1] += L
            ################################################
            # haplotype 2
            ################################################
            allele2 = alleles[1]
            if allele2 not in {'0', '.'}:
                local_allele_counts[hap2] += 1
                local_total_lengths[hap2] += L
    return (
        local_allele_counts,
        local_total_lengths
    )


############################################################
# Merge chunk results
############################################################

def merge_counts(
    global_dict,
    local_dict
):
    for k, v in local_dict.items():
        global_dict[k] += v


############################################################
# Parallel processing
############################################################

with gzip.open(vcf_path, 'rt') as f:
    buffer = []
    chunk_size = 10000
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=8
    ) as executor:
        futures = []
        for line in f:
            if line.startswith('#'):
                continue
            buffer.append(line)
            if len(buffer) >= chunk_size:
                futures.append(
                    executor.submit(
                        process_chunk,
                        buffer
                    )
                )
                buffer = []
        ####################################################
        # Remaining variants
        ####################################################
        if buffer:
            futures.append(
                executor.submit(
                    process_chunk,
                    buffer
                )
            )
        ####################################################
        # Collect results
        ####################################################
        for future in concurrent.futures.as_completed(
            futures
        ):
            ac, tl = future.result()
            merge_counts(
                allele_counts,
                ac
            )
            merge_counts(
                total_lengths,
                tl
            )


############################################################
# Generate output
############################################################

df = pd.DataFrame({
    'Haplotype': haplotype_names,
    'Sample': [
        h.rsplit(
            '.hap',
            1
        )[0]
        for h in haplotype_names
    ],
    'AlleleCount': [
        allele_counts[h]
        for h in haplotype_names
    ],
    'TotalLength': [
        total_lengths[h]
        for h in haplotype_names
    ]
})


############################################################
# Save
############################################################



df.to_csv(
    out_path,
    sep='\t',
    index=False
)


print(
    f"Saved {len(df)} haplotypes"
)
