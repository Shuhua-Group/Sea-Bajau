#!/usr/bin/env bash

# Prepare population datasets and summarize genetic structure.
# Set input, reference and tool paths for your environment before execution.

set -euo pipefail

wd="/path/to/wd"
ho_info="/path/to/ho_1926.kp.info.txt"
PLINK2="/path/to/plink2"
PLINK1="plink"
ADMIXTURE="admixture"  # Set an absolute path here if ADMIXTURE v1.3.0 is not on PATH.
threads_total=8
parallel_jobs=4
threads_per_job=2
maf=0.01
n_reps=50
cv_folds=10
seed_base=20260905


source_root="${wd}/02.plink_prep"
apc_bfile="${source_root}/QC.prefixed"
ho_bfile="${source_root}/Ho.1926.QC"
outroot="${wd}/admixture"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
manifest_dir="${outroot}/00.manifests"

[[ -x "$PLINK2" ]] || { echo "ERROR: cannot execute PLINK2: $PLINK2" >&2; exit 1; }
command -v "$PLINK1" >/dev/null || { echo "ERROR: PLINK 1.9 unavailable: $PLINK1" >&2; exit 1; }
command -v "$ADMIXTURE" >/dev/null || { echo "ERROR: ADMIXTURE unavailable: $ADMIXTURE" >&2; exit 1; }
[[ -r "$ho_info" && -r "${apc_bfile}.bed" && -r "${ho_bfile}.bed" ]] || {
  echo "ERROR: missing Ho metadata or prepared APC/Ho binary input." >&2; exit 1;
}
[[ "$parallel_jobs" -ge 1 && "$threads_per_job" -ge 1 ]] || exit 1
[[ $((parallel_jobs * threads_per_job)) -le "$threads_total" ]] || {
  echo "ERROR: parallel_jobs * threads_per_job exceeds threads_total." >&2; exit 1;
}

mkdir -p "$outroot"/{00.manifests,01.panel_prep,02.admixture_runs,03.summary}
printf 'Panel\tHo_selected_samples\tHo_reference_samples_after_joint_QC\tAPC_samples_after_joint_QC\tLD_pruned_SNPs\n' \
  > "${outroot}/03.summary/panel_sizes.tsv"

# 1. Reference selection is metadata-only and precedes all ADMIXTURE runs.
Rscript "${script_dir}/make_ho1926_admixture_manifests.R" "$ho_info" "$manifest_dir"

# 2. Build a joint, autosomal, biallelic shared-SNP panel for each manifest.
for spec in "A_asia_focused" "C_broader_reference"; do
  keep="${manifest_dir}/Ho1926_${spec}.keep"
  panel_dir="${outroot}/01.panel_prep/${spec}"
  mkdir -p "$panel_dir"
  cd "$panel_dir" || exit 1

  "$PLINK2" --bfile "$ho_bfile" --keep "$keep" --make-bed --out Ho.selected
  n_ho=$(awk 'END {print NR}' Ho.selected.fam)
  [[ "$n_ho" -gt 0 ]] || { echo "ERROR: no selected Ho samples for $spec." >&2; exit 1; }

  awk 'FNR == NR {id[$2] = 1; next} ($2 in id) {print $2}' \
    "${apc_bfile}.bim" Ho.selected.bim | sort -u > common_snps.txt
  [[ -s common_snps.txt ]] || { echo "ERROR: no shared SNPs for $spec." >&2; exit 1; }

  "$PLINK2" --bfile "$apc_bfile" --extract common_snps.txt --make-bed --out APC.common
  "$PLINK2" --bfile Ho.selected --extract common_snps.txt --make-bed --out Ho.common
  "$PLINK1" --bfile Ho.common --bmerge APC.common --make-bed --out Merge.common

  # Reapply joint-call-rate QC after merging; MAF is evaluated in the joint panel.
  "$PLINK2" --bfile Merge.common --mind 0.10 --geno 0.05 --maf "$maf" \
    --make-bed --out "Merge.MAF${maf}.QC"
  "$PLINK2" --bfile "Merge.MAF${maf}.QC" --indep-pairwise 50 5 0.2 \
    --out "Merge.MAF${maf}.LD"
  "$PLINK2" --bfile "Merge.MAF${maf}.QC" --extract "Merge.MAF${maf}.LD.prune.in" \
    --make-bed --out "Merge.MAF${maf}.LDpruned"

  n_apc=$(awk '$2 ~ /^APC_/ {n++} END {print n+0}' "Merge.MAF${maf}.LDpruned.fam")
  n_ho_final=$(awk '$2 !~ /^APC_/ {n++} END {print n+0}' "Merge.MAF${maf}.LDpruned.fam")
  [[ "$n_apc" -eq 130 ]] || {
    echo "ERROR: expected 130 APC samples after joint QC in $spec, found $n_apc." >&2; exit 1;
  }
  printf '%s\t%s\t%s\t%s\t%s\n' "$spec" "$n_ho" "$n_ho_final" "$n_apc" \
    "$(awk 'END {print NR}' Merge.MAF${maf}.LDpruned.bim)" \
    >> "${outroot}/03.summary/panel_sizes.tsv"
done

# 3. ADMIXTURE CV runs. Each run has an isolated directory, fixed seed, and log.
run_one() {
  local spec="$1" k="$2" rep="$3" seed="$4"
  local input_prefix="${outroot}/01.panel_prep/${spec}/Merge.MAF${maf}.LDpruned"
  local run_dir="${outroot}/02.admixture_runs/${spec}/K${k}/rep$(printf '%03d' "$rep")"
  mkdir -p "$run_dir"
  ln -sf "$input_prefix.bed" "$run_dir/input.bed"
  ln -sf "$input_prefix.bim" "$run_dir/input.bim"
  ln -sf "$input_prefix.fam" "$run_dir/input.fam"
  (
    cd "$run_dir" || exit 1
    "$ADMIXTURE" -j"$threads_per_job" --cv="$cv_folds" -s "$seed" input.bed "$k" \
      > "K${k}.rep$(printf '%03d' "$rep").log" 2>&1
    mv "input.${k}.Q" "K${k}.rep$(printf '%03d' "$rep").Q"
    mv "input.${k}.P" "K${k}.rep$(printf '%03d' "$rep").P"
  )
}

export -f run_one
export outroot maf ADMIXTURE threads_per_job cv_folds

for spec in "A_asia_focused" "C_broader_reference"; do
  for k in $(seq 2 12); do
    for rep in $(seq 1 "$n_reps"); do
      seed=$((seed_base + k * 1000 + rep))
      if [[ "$spec" == "C_broader_reference" ]]; then seed=$((seed + 100000)); fi
      while [[ "$(jobs -pr | wc -l)" -ge "$parallel_jobs" ]]; do wait -n; done
      run_one "$spec" "$k" "$rep" "$seed" &
    done
  done
done
wait

# Collect CV lines only after every expected run has completed successfully.
find "${outroot}/02.admixture_runs" -name '*.log' -print0 | sort -z | \
  xargs -0 grep -H 'CV error' > "${outroot}/03.summary/admixture_cv_raw.txt"

Rscript "${script_dir}/summarize_admixture_cv.R" "$outroot"
