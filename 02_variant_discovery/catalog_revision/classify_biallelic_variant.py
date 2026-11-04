#!/usr/bin/env python3
"""Classify one normalized, explicit, biallelic CHM13 allele per row.

Input:  tab-separated CHROM POS REF ALT [optional extra columns]
Output: input's first four fields plus CLASS, REF_LEN, ALT_LEN, LEN_DIFF.
"""
import sys


def classify(ref: str, alt: str) -> str:
    if "," in alt:
        raise ValueError("multiallelic ALT encountered; run bcftools norm -m -any first")
    if alt.startswith("<") or "[" in alt or "]" in alt:
        raise ValueError("symbolic/BND ALT encountered; this sequence-resolved classifier is not applicable")
    if ref == alt:
        return "REF"
    ref_len, alt_len = len(ref), len(alt)
    if ref_len == alt_len:
        return "SNP" if ref_len == 1 else "MNP"
    diff = alt_len - ref_len
    stem = "INS" if diff > 0 else "DEL"
    return f"INDEL_{stem}" if abs(diff) < 50 else f"SV_{stem}"


def main() -> None:
    print("\t".join(["CHROM", "POS", "REF", "ALT", "CLASS", "REF_LEN", "ALT_LEN", "LEN_DIFF"]))
    for line in sys.stdin:
        if not line.strip() or line.startswith("#"):
            continue
        chrom, pos, ref, alt, *_ = line.rstrip("\n").split("\t")
        print("\t".join([chrom, pos, ref, alt, classify(ref, alt),
                          str(len(ref)), str(len(alt)), str(len(alt) - len(ref))]))


if __name__ == "__main__":
    main()
