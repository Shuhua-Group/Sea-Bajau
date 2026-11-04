"""Analyze digestive-enzyme copy numbers and haplotypes."""

import os

def extract_third_col(right_part):
    if not right_part:
        return ""
    first_hash_idx = right_part.find('#')
    if first_hash_idx == -1:
        substr = right_part
    else:
        substr = right_part[:first_hash_idx]
    if '-' in substr:
        substr = substr.split('-')[-1]
    return substr

def parse_line(path, second_col):
    group = os.path.basename(path).split('.')[0]
    if '-' in second_col:
        left_part, right_part = second_col.split('-', 1)
    else:
        left_part = second_col
        right_part = ""
    part1 = left_part.split('#')[0]
    part2 = extract_third_col(right_part)
    return group, part1, part2

def process_files(input_dir, output_file):
    with open(output_file, "w") as fout:
        for root, dirs, files in os.walk(input_dir):
            for filename in files:
                if not filename.endswith(".best_genotype.tsv"):
                    continue
                filepath = os.path.join(root, filename)
                with open(filepath) as fin:
                    for line in fin:
                        line = line.strip()
                        if line.startswith("#"):
                            continue
                        cols = line.split()
                        if len(cols) < 2:
                            continue
                        if float(cols[2]) < 0.8:
                            continue
                        path = cols[0]
                        second_col = cols[1]
                        parsed = parse_line(path, second_col)
                        if parsed:
                            fout.write("\t".join(parsed) + "\n")

if __name__ == "__main__":
    input_dir = "/path/to/2018Cell_statistic"  #  Root directory ,  Set the actual path
    output_file = "/path/to/2018Cell_statistic.merged_output.txt"     #  Results file
    process_files(input_dir, output_file)
    print(f"All done! Results saved to {output_file}")
