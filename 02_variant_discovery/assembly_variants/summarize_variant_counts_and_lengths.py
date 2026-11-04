"""Process variant calls and summarize variant catalogs."""

import gzip
import concurrent.futures
import pandas as pd
from collections import defaultdict

# vcf_path = "/path/to/130APC.SVPOP_sv_del.telo.vcf.gz" # Input VCF File
# vcf_path = "/path/to/130APC.SVPOP_sv_ins.telo.vcf.gz" # Input VCF File
# vcf_path = "/path/to/130APC.SVPOP_sv_del.no_telo.cenSat.vcf.gz" # Input VCF File
vcf_path = "/path/to/130APC.SVPOP_sv_ins.no_telo.cenSat.vcf.gz"  #  Input  VCF  File




sample_names = []

def parse_vcf_header(vcf_path):
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                return fields[9:]
    return []

sample_names = parse_vcf_header(vcf_path)

# Initialize count accumulators
allele_counts = defaultdict(int)
total_lengths = defaultdict(int)

def process_chunk(lines):
    local_allele_counts = defaultdict(int)
    local_total_lengths = defaultdict(int)
    for line in lines:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        ref, alt = fields[3], fields[4]
        L = abs(len(ref) - len(alt))
        genotypes = fields[9:]
        for i, gt_field in enumerate(genotypes):
            gt = gt_field.split(':')[0]
            alleles = gt.replace('|', '/').split('/')
            for allele in alleles:
                if allele not in {'0', '.'}:
                    local_allele_counts[sample_names[i]] += 1
                    local_total_lengths[sample_names[i]] += L
    return local_allele_counts, local_total_lengths

def merge_counts(global_dict, local_dict):
    for k, v in local_dict.items():
        global_dict[k] += v

# Process chunks with a thread pool
with gzip.open(vcf_path, 'rt') as f:
    buffer = []
    chunk_size = 10000  #  Configurable
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        futures = []
        for line in f:
            if line.startswith('#'):
                continue
            buffer.append(line)
            if len(buffer) >= chunk_size:
                futures.append(executor.submit(process_chunk, buffer))
                buffer = []
        if buffer:
            futures.append(executor.submit(process_chunk, buffer))
        for future in concurrent.futures.as_completed(futures):
            ac, tl = future.result()
            merge_counts(allele_counts, ac)
            merge_counts(total_lengths, tl)

# Save results
df = pd.DataFrame({
    'Sample': sample_names,
    'AlleleCount': [allele_counts[s] for s in sample_names],
    'TotalLength': [total_lengths[s] for s in sample_names]
})

# df.to_csv("/path/to/130APC.SVPOP_sv_del.telo.sample_variant_stats.python.txt", sep='\t', index=False)
# df.to_csv("/path/to/130APC.SVPOP_sv_ins.telo.sample_variant_stats.python.txt", sep='\t', index=False)
# df.to_csv("/path/to/130APC.SVPOP_sv_del.no_telo.cenSat.sample_variant_stats.python.txt", sep='\t', index=False)
df.to_csv("/path/to/130APC.SVPOP_sv_ins.no_telo.cenSat.sample_variant_stats.python.txt", sep='\t', index=False)
