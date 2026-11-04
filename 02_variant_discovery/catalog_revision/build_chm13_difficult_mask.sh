#!/usr/bin/env bash
# Build a transparent CHM13v2 exclusion mask from source annotations.
# Do not substitute GRCh38 coordinates or hand-enter coordinates.
# Process variant calls and summarize variant catalogs.
# Set input, reference and tool paths for your environment before execution.

set -euo pipefail

CENTROMERE_BED="/PATH/TO/CHM13v2_centromeres.bed"
TELOMERE_SUBTELOMERE_BED="/PATH/TO/CHM13v2_telomere_subtelomere.bed"
SATELLITE_BED="/PATH/TO/CHM13v2_satellite_or_repeat_difficult.bed"
OUT_BED_GZ="/PATH/TO/CHM13v2_difficult_regions.bed.gz"

for exe in bedtools bgzip tabix; do command -v "$exe" >/dev/null || { echo "ERROR: $exe unavailable" >&2; exit 1; }; done
for x in "$CENTROMERE_BED" "$TELOMERE_SUBTELOMERE_BED" "$SATELLITE_BED"; do
  [[ -r "$x" ]] || { echo "ERROR: cannot read $x" >&2; exit 1; }
done

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
cat "$CENTROMERE_BED" "$TELOMERE_SUBTELOMERE_BED" "$SATELLITE_BED" |
  awk 'BEGIN{OFS="\t"} $0 !~ /^#/ && $0 !~ /^track/ && $0 !~ /^browser/ && NF>=3 {print $1,$2,$3}' |
  bedtools sort -i - | bedtools merge -i - | bgzip -c > "$OUT_BED_GZ"
tabix -f -p bed "$OUT_BED_GZ"
echo "Wrote $(tabix -l "$OUT_BED_GZ" | wc -l | tr -d ' ') contigs to $OUT_BED_GZ" >&2
