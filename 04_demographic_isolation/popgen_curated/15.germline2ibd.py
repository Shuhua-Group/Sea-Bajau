"""Summarize population variation and demographic history."""

import gzip
import sys
import re
import os

# Input file, For example: 4panel.PanMLS.Temiar.chr1
infile = sys.argv[1]

# Extract from the filename chr Name, For example chr1
m = re.search(r'(chr[\w\d]+)', infile)
if not m:
    raise ValueError("Filename must contain 'chr*', e.g. chr1, chrX")
chrom = m.group(1)

outfile = infile + ".ibd.gz"

def split_sample(s):
    """Split AAGC031016D.0 → (AAGC031016D, 1 or 2)"""
    base, hap = s.split(".")
    hap = "1" if hap == "0" else "2"
    return base, hap

with gzip.open(outfile, "wt") as out:
    with open(infile) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            cols = line.split()

            # Original column: 0,1,2,3,4,5,6
            # Use only 0,1,2,3,4
            s1, h1 = split_sample(cols[0])
            s2, h2 = split_sample(cols[1])

            start = cols[2]
            end   = cols[3]
            score = cols[4]

            out.write(f"{s1}\t{h1}\t{s2}\t{h2}\t{chrom}\t{start}\t{end}\t{score}\n")

print(f"Done! Output written to: {outfile}")
