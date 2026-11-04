#!/usr/bin/env python3
"""Stream exact sequence-allele membership from a source-sorted input."""
import argparse
import csv
from collections import Counter


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="CHROM POS REF ALT CLASS SOURCE; sorted by first four fields")
    p.add_argument("--out", required=True)
    a = p.parse_args()
    def emit(w, key, klass, members, counts):
        if key is None: return
        membership = "_".join(x for x in ("APC", "CPC2", "HPRC2") if x in members)
        counts[(membership, klass)] += 1
        chrom, pos, ref, alt = key
        w.writerow({"CHROM":chrom, "POS":pos, "REF":ref, "ALT":alt, "CLASS":klass,
                    "APC":int("APC" in members), "CPC2":int("CPC2" in members), "HPRC2":int("HPRC2" in members),
                    "MEMBERSHIP":membership, "BLOCK_1MB":f"{chrom}:{(int(pos)-1)//1000000}"})

    with open(a.out + ".detailed.tsv", "w", newline="") as out:
        fields = ["CHROM", "POS", "REF", "ALT", "CLASS", "APC", "CPC2", "HPRC2", "MEMBERSHIP", "BLOCK_1MB"]
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t"); w.writeheader()
        counts = Counter(); current_key = None; current_class = None; members = set()
        with open(a.input, newline="") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if len(row) != 6: raise ValueError("expected six fields in source-sorted input")
                key, klass, source = tuple(row[:4]), row[4], row[5]
                if current_key is not None and key != current_key:
                    emit(w, current_key, current_class, members, counts)
                    members = set()
                if current_key == key and klass != current_class:
                    raise ValueError(f"inconsistent class for {key}")
                current_key, current_class = key, klass
                members.add(source)
        emit(w, current_key, current_class, members, counts)
    with open(a.out + ".counts.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter="\t"); w.writerow(["MEMBERSHIP", "CLASS", "N"])
        for (membership, klass), n in sorted(counts.items()): w.writerow([membership, klass, n])


if __name__ == "__main__": main()
