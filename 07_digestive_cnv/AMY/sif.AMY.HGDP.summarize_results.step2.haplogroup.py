"""Analyze digestive-enzyme copy numbers and haplotypes."""

from collections import Counter

def load_sequence_map(file_a):
    seq_map = {}
    with open(file_a) as fin:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            sample_id = parts[0]
            sequence = parts[1]
            seq_map[sample_id] = sequence
    return seq_map

def map_contigs_to_sequences(file_b, seq_map, output_file):
    with open(file_b) as fin, open(output_file, 'w') as fout:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            sample_group = parts[0]
            contig1 = parts[1]
            contig2 = parts[2]
            seq1 = seq_map.get(contig1, "NA")
            seq2 = seq_map.get(contig2, "NA")
            fout.write(f"{sample_group}\t{seq1}\t{seq2}\n")

def count_PGA_per_sample(input_file, output_count_file):
    with open(input_file) as fin, open(output_count_file, 'w') as fout:
        fout.write("sample_hap\thaplogroup\n")
        for line in fin:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            sample = parts[0]
            hap1_seq = parts[1]
            hap2_seq = parts[2]
            # hap1 Count
            fout.write(f"{sample}_hap1\t{hap1_seq}\n")
            # hap2 Count
            fout.write(f"{sample}_hap2\t{hap2_seq}\n")


if __name__ == "__main__":
    # ... Previous code ...
    file_a = "/path/to/AMY.GT.haplogroup.out"        #  First input file
    file_b = "/path/to/HGDP_statistic.merged_output.txt"       #  Second input file
    output_file = "/path/to/HGDP_statistic.merged_output.haplogroup.pattern.txt"  #  Output file
    output_count_file = "/path/to/HGDP_statistic.merged_output.haplogroup.pattern.hap_AMY_counts.txt"
    seq_map = load_sequence_map(file_a)
    map_contigs_to_sequences(file_b, seq_map, output_file)
    print(f"Step1 done! Output saved to {output_file}")
    count_PGA_per_sample(output_file, output_count_file)
    print(f"Step2 done! PGA counts saved to {output_count_file}")
